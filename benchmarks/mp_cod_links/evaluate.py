"""Evaluate first-principles parent reconstruction on the MP <-> COD links (links.csv).

For every linked MP entry (DISORDERED_REALIZATION links by default):
  * reconstruct() -> ranked candidate parents;
  * recovery: DisorderMatcher().match_disordered(candidate, cod_parent) for every candidate
    (PRISTINE = same disordered phase incl. composition, DISORDERED_FAMILY = same framework and
    site-species pattern);
  * transfer: for every OTHER MP entry linked to the same COD entry, DisorderMatcher().match(other,
    top-1 candidate) -- does the reconstructed parent catch the other orderings?
False-positive probe (--fp): MP entries of chemical systems that have NO disordered COD entry
(but ordered COD entries): how often does reconstruction propose a disordered parent?

Outputs eval_links.jsonl / eval_fp.jsonl and a summary on stdout.
usage: python evaluate.py [--procs 4] [--max-links 1500] [--fp 400]
"""
from __future__ import annotations

import argparse, collections, csv, gzip, json, os, pickle, random, signal, sys, time, warnings

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import data  # noqa: E402

PARENT_DIR = os.path.join(HERE, "parents")


class _Timeout(Exception):
    pass


def _alarm(*_):
    raise _Timeout


_STRICT = None


def strict_match(a, b, occ_tol=0.15):
    """Same disordered structure: identical element set on every matched site, occupancies within
    occ_tol (stricter than match_disordered, which accepts element-set containment)."""
    global _STRICT
    if _STRICT is None:
        from pymatgen.analysis.structure_matcher import AbstractComparator, StructureMatcher

        class StrictSiteComparator(AbstractComparator):
            def are_equal(self, sp1, sp2):
                if {e.symbol for e in sp1.elements} != {e.symbol for e in sp2.elements}:
                    return False
                return all(abs(sp1[e] - sp2[e]) <= occ_tol for e in sp1.elements)

            def get_hash(self, composition):
                return "-".join(sorted(e.symbol for e in composition.elements))

        _STRICT = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True, scale=True,
                                   attempt_supercell=True, comparator=StrictSiteComparator())
    try:
        return bool(_STRICT.fit(a, b))
    except Exception:
        return False


def cod_features(s):
    """How disordered is the COD entry? min over sites of the majority occupancy."""
    maj = [max(sp.values()) if len(sp) else 0 for sp in (dict(x.species.items()) for x in s)]
    return dict(min_major_occ=round(min(maj), 3))


def _eval(job):
    from opendisorder.matching import DisorderMatcher, MatchKind
    from reconstruct import reconstruct

    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(job["timeout"])
    t = time.time()
    out = dict(mp_id=job["mp_id"], cod_id=job.get("cod_id"), cands=[], transfer=[])
    try:
        cands = reconstruct(job["mp"])
        dm = DisorderMatcher()
        for i, c in enumerate(cands[:8]):
            rec = dict(rank=i, kind=c.kind, groups=c.groups, vac=c.vacancy, index=round(c.index, 3),
                       sg=c.spacegroup, symprec=c.symprec, formula=c.parent.composition.reduced_formula,
                       nsites=len(c.parent))
            if job.get("cod") is not None:
                try:
                    rec["vs_cod"] = dm.match_disordered(c.parent, job["cod"]).kind.value
                except Exception:
                    rec["vs_cod"] = "error"
                rec["strict"] = rec["vs_cod"] == "pristine" and strict_match(c.parent, job["cod"])
            out["cands"].append(rec)
            if i == 0 and job.get("save"):
                os.makedirs(PARENT_DIR, exist_ok=True)
                c.parent.to(filename=os.path.join(PARENT_DIR, f"{job['mp_id']}__{job.get('cod_id')}.cif"))
        if cands and job.get("others"):
            # which candidate to test transfer with: the first that recovers the COD phase, else top-1
            best = next((c for c, r in zip(cands, out["cands"]) if r.get("strict")), cands[0])
            for oid, os_ in job["others"]:
                try:
                    k = dm.match(os_, best.parent).kind.value
                except Exception:
                    k = "error"
                out["transfer"].append((oid, k))
    except _Timeout:
        out["timeout"] = True
    except Exception as exc:
        out["error"] = str(exc)[:200]
    finally:
        signal.alarm(0)
    out["dt"] = round(time.time() - t, 2)
    return out


def load_links(path, kinds):
    rows = [r for r in csv.DictReader(open(path)) if r["kind"] in kinds]
    return rows


def summarize(res, label):
    n = len(res)
    ok = [r for r in res if "timeout" not in r and "error" not in r]
    anyc = sum(bool(r["cands"]) for r in ok)
    top1 = sum(bool(r["cands"]) and r["cands"][0].get("vs_cod") == "pristine" for r in ok)
    top3 = sum(any(c.get("vs_cod") == "pristine" for c in r["cands"][:3]) for r in ok)
    anyp = sum(any(c.get("vs_cod") == "pristine" for c in r["cands"]) for r in ok)
    anyf = sum(any(c.get("vs_cod") in ("pristine", "disordered_family") for c in r["cands"]) for r in ok)
    s1 = sum(bool(r["cands"]) and bool(r["cands"][0].get("strict")) for r in ok)
    sany = sum(any(c.get("strict") for c in r["cands"]) for r in ok)
    tr = collections.Counter(k for r in ok for _, k in r["transfer"])
    print(f"[{label}] strict: top1={s1} any={sany}")
    print(f"[{label}] n={n} ok={len(ok)} with_candidate={anyc} top1_pristine={top1} top3_pristine={top3} "
          f"any_pristine={anyp} any_pristine_or_family={anyf} transfer={dict(tr)}")
    return dict(n=n, ok=len(ok), with_candidate=anyc, top1=top1, top3=top3, any_pristine=anyp,
                any_family=anyf, strict_top1=s1, strict_any=sany, transfer=dict(tr))


def main():
    from multiprocessing import Pool

    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--links", default=os.path.join(HERE, "links.csv"))
    ap.add_argument("--max-links", type=int, default=1500)
    ap.add_argument("--kinds", default="disordered_realization")
    ap.add_argument("--fp", type=int, default=0, help="number of MP entries for the false-positive probe")
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--fp-only", action="store_true")
    a = ap.parse_args()

    if a.fp_only:
        print("[false-positive probe]", false_positive_probe(a, rng := random.Random(1)))
        return
    with gzip.open(os.path.join(HERE, "cod_parsed.pkl.gz"), "rb") as fh:
        cod = {c["cod_id"]: c for c in pickle.load(fh)}
    rows = load_links(a.links, set(a.kinds.split(",")))
    by_cod = collections.defaultdict(list)
    for r in rows:
        by_cod[r["cod_id"]].append(r["mp_id"])
    rng = random.Random(1)
    rng.shuffle(rows)
    # one link per MP entry (its first COD partner) keeps the sample diverse
    seen, sample = set(), []
    for r in rows:
        if r["mp_id"] not in seen:
            seen.add(r["mp_id"])
            sample.append(r)
    sample = sample[: a.max_links]

    mp_cache = {}

    from pymatgen.core import Composition

    mp_formula = {r["mp_id"]: r["mp_formula"] for r in rows}

    def mp_struct(_cs, mid):
        cs = "-".join(sorted(e.symbol for e in Composition(mp_formula[mid]).elements))
        if cs not in mp_cache:
            mp_cache[cs] = {d["material_id"]: data.mp_structure(d) for d in data.mp_docs(cs)}
        return mp_cache[cs][mid]

    jobs = []
    for r in sample:
        c = cod[r["cod_id"]]
        others = [(o, mp_struct(r["chemsys"], o)) for o in by_cod[r["cod_id"]] if o != r["mp_id"]][:5]
        jobs.append(dict(mp_id=r["mp_id"], cod_id=r["cod_id"], mp=mp_struct(r["chemsys"], r["mp_id"]),
                         cod=c["structure"], others=others, timeout=a.timeout, save=True,
                         link=r, cod_feat=cod_features(c["structure"])))
    print("jobs", len(jobs), flush=True)
    res = []
    t0 = time.time()
    with Pool(a.procs) as pool, open(os.path.join(HERE, "eval_links.jsonl"), "w") as fh:
        it = pool.imap_unordered(_eval, [{k: v for k, v in j.items() if k not in ("link", "cod_feat")} for j in jobs])
        for i in range(len(jobs)):
            try:
                out = it.next(timeout=600)  # C-level hangs (spglib / supercell search) ignore SIGALRM
            except Exception:
                print("stuck workers, giving up on", len(jobs) - i, "jobs", flush=True)
                pool.terminate()
                break
            j = next(x for x in jobs if x["mp_id"] == out["mp_id"])
            out.update(link={k: j["link"][k] for k in ("kind", "chemsys", "mp_formula", "cod_formula", "mp_sg",
                                                         "cod_sg", "mp_nsites", "cod_nsites", "theoretical")},
                       **j["cod_feat"])
            res.append(out)
            fh.write(json.dumps(out) + "\n")
            if i % 100 == 0:
                print(i, f"{time.time() - t0:.0f}s", flush=True)
    summ = {"all": summarize(res, "all links")}
    subst = [r for r in res if r["min_major_occ"] <= 0.9]
    summ["substantial"] = summarize(subst, "COD min majority occ <= 0.9")
    nontriv = [r for r in subst if int(r["link"]["mp_nsites"]) > int(r["link"]["cod_nsites"])]
    summ["nontrivial"] = summarize(nontriv, "... and MP cell larger than COD primitive cell")

    if a.fp:
        fp = false_positive_probe(a, rng)
        print("[false-positive probe]", fp)
        summ["fp"] = fp
    json.dump(summ, open(os.path.join(HERE, "eval_summary.json"), "w"), indent=1)


def fp_chemsys():
    conn = data.cod_connect()
    dis_cs = {r[0] for r in conn.execute(
        "SELECT DISTINCT k.value FROM entries e JOIN keys k ON k.cod_id=e.cod_id AND k.kind='cs' "
        "WHERE e.has_disorder=1 OR e.fractional=1")}
    ord_cs = [r[0] for r in conn.execute(
        "SELECT k.value, count(*) FROM entries e JOIN keys k ON k.cod_id=e.cod_id AND k.kind='cs' "
        "WHERE e.has_disorder=0 AND e.fractional=0 AND e.natoms_cell <= 60 GROUP BY k.value")]
    cands = sorted(cs for cs in ord_cs if cs not in dis_cs and "H" not in cs.split("-")
                   and "C" not in cs.split("-") and 2 <= len(cs.split("-")) <= 5)
    random.Random(2).shuffle(cands)
    return cands


def false_positive_probe(a, rng):
    """MP entries whose chemical system has ordered COD entries but no disordered one."""
    from multiprocessing import Pool

    cands = fp_chemsys()
    data.fetch_mp_chemsys(cands[: a.fp // 2])
    jobs = []
    for cs in cands[: a.fp // 2]:
        for d in data.mp_docs(cs)[:4]:
            jobs.append(dict(mp_id=d["material_id"], mp=data.mp_structure(d), timeout=a.timeout,
                             theo=d.get("theoretical"), cs=cs))
    jobs = jobs[: a.fp]
    print("fp jobs", len(jobs), flush=True)
    res = []
    with Pool(a.procs) as pool, open(os.path.join(HERE, "eval_fp.jsonl"), "w") as fh:
        it = pool.imap_unordered(_eval, [{k: v for k, v in j.items() if k not in ("theo", "cs")} for j in jobs])
        for _ in range(len(jobs)):
            try:
                out = it.next(timeout=600)
            except Exception:
                pool.terminate()
                break
            j = next(x for x in jobs if x["mp_id"] == out["mp_id"])
            out.update(theoretical=j["theo"], chemsys=j["cs"])
            res.append(out)
            fh.write(json.dumps(out) + "\n")
    return _fp_summary(res)


def _fp_summary(res):
    ok = [r for r in res if "timeout" not in r and "error" not in r]
    kinds = collections.Counter(r["cands"][0]["kind"] for r in ok if r["cands"])
    sim = sum(any(c["kind"].startswith("similar") for c in r["cands"]) for r in ok)
    out = dict(n=len(res), ok=len(ok), any_candidate=sum(bool(r["cands"]) for r in ok),
               any_similar=sim, top_kind=dict(kinds),
               any_candidate_exp=sum(bool(r["cands"]) for r in ok if not r["theoretical"]),
               n_exp=sum(1 for r in ok if not r["theoretical"]))
    print("[false-positive probe]", out)
    return out


if __name__ == "__main__":
    main()
