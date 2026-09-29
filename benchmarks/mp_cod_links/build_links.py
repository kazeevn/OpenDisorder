"""Ground-truth links "MP ordered structure <-> COD disordered entry" (no ICSD, no MP provenance).

1. COD candidates: inorganic entries flagged disordered/fractional, <= 60 atoms/cell, 2-6 elements;
   parsed structures cached in cod_parsed.pkl.gz (only genuinely disordered ones kept).
2. MP: all summary docs of the same chemical system (data.fetch_mp_chemsys).
3. Cheap prefilter (loose composition, volume/atom, supercell ratio), then
   DisorderMatcher().match(mp_ordered, cod_disordered) with a per-pair timeout.
4. links.csv: every DISORDERED_REALIZATION / DISORDERED_FAMILY pair.

usage: python build_links.py [--procs 4] [--max-per-chemsys 8] [--limit-chemsys N]
"""
from __future__ import annotations

import argparse, collections, csv, gzip, os, pickle, random, signal, sys, time, warnings
from multiprocessing import Pool

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import data  # noqa: E402

COD_PKL = os.path.join(HERE, "cod_parsed.pkl.gz")


# ------------------------------------------------------------------ COD parsing
_conn = None


def _parse(c):
    global _conn
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    _conn = _conn or data.cod_connect()
    try:
        s = data.cod_structure(c["cod_id"], _conn)
    except Exception:
        return None
    if s is None or s.is_ordered or len(s) > 400:
        return None
    try:
        s = s.get_primitive_structure()  # merged-site CIFs stay as they are
    except Exception:
        pass
    return dict(c, structure=s, nsites=len(s))


def load_cod(procs=4):
    if os.path.exists(COD_PKL):
        with gzip.open(COD_PKL, "rb") as fh:
            return pickle.load(fh)
    cands = data.cod_candidates()
    with Pool(procs) as pool:
        out = [r for r in pool.imap_unordered(_parse, cands, chunksize=20) if r]
    with gzip.open(COD_PKL, "wb") as fh:
        pickle.dump(out, fh)
    return out


# ------------------------------------------------------------------ prefilter
def comp_close(c1, c2, atol=0.05, rtol=0.3):
    f1, f2 = c1.fractional_composition, c2.fractional_composition
    return all(abs(f1[e] - f2[e]) <= atol + rtol * (f1[e] + f2[e]) / 2 for e in set(f1) | set(f2))


def prefilter(mp_s, cod_s):
    if not comp_close(mp_s.composition, cod_s.composition):
        return False
    n_cod = cod_s.composition.num_atoms
    v1, v2 = mp_s.volume / len(mp_s), cod_s.volume / max(n_cod, 1e-6)
    if not 0.7 < v1 / v2 < 1.45:
        return False
    return len(mp_s) / max(n_cod, 1) <= 16.5


# ------------------------------------------------------------------ matching
class _Timeout(Exception):
    pass


def _alarm(*_):
    raise _Timeout


def _match(job):
    from opendisorder.matching import DisorderMatcher
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    mp_id, mp_s, cod_id, cod_s, timeout = job
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(timeout)
    t = time.time()
    try:
        res = DisorderMatcher().match(mp_s, cod_s)
        kind, det = res.kind.value, res.details
    except _Timeout:
        kind, det = "timeout", {}
    except Exception as exc:
        kind, det = "error", {"reason": str(exc)[:100]}
    finally:
        signal.alarm(0)
    return mp_id, cod_id, kind, det, time.time() - t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--max-per-chemsys", type=int, default=8)
    ap.add_argument("--limit-chemsys", type=int, default=0)
    ap.add_argument("--timeout", type=int, default=30)
    ap.add_argument("--out", default=os.path.join(HERE, "links.csv"))
    ap.add_argument("--cached-only", action="store_true", help="skip chemsys not yet in mp_cache/")
    a = ap.parse_args()

    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    cod = load_cod(a.procs)
    print("COD disordered parsed:", len(cod), flush=True)
    by_cs = collections.defaultdict(list)
    for c in cod:
        by_cs[c["chemsys"]].append(c)
    rng = random.Random(0)
    css = sorted(by_cs)
    rng.shuffle(css)
    if a.cached_only:
        css = [cs for cs in css if os.path.exists(data._cache_file(cs))]
    if a.limit_chemsys:
        css = css[: a.limit_chemsys]

    jobs, meta = [], {}
    for cs in css:
        entries = by_cs[cs]
        rng.shuffle(entries)
        try:
            docs = data.mp_docs(cs)
        except Exception as exc:
            print("mp fail", cs, exc)
            continue
        if not docs:
            continue
        mps = [(d["material_id"], data.mp_structure(d), d) for d in docs]
        for c in entries[: a.max_per_chemsys]:
            for mid, s, d in mps:
                if prefilter(s, c["structure"]):
                    jobs.append((mid, s, c["cod_id"], c["structure"], a.timeout))
                    meta[mid] = d
            meta[c["cod_id"]] = c
    print("chemsys", len(css), "pairs after prefilter", len(jobs), flush=True)

    t0 = time.time()
    rows, stats = [], collections.Counter()
    with Pool(a.procs) as pool:
        for i, (mid, cid, kind, det, dt) in enumerate(pool.imap_unordered(_match, jobs, chunksize=4)):
            stats[kind] += 1
            if kind in ("disordered_realization", "disordered_family"):
                d, c = meta[mid], meta[cid]
                rows.append(dict(
                    mp_id=mid, cod_id=cid, kind=kind, chemsys=c["chemsys"],
                    mp_formula=d["formula_pretty"], cod_formula=c["structure"].composition.reduced_formula,
                    mp_sg=(d.get("symmetry") or {}).get("number"), cod_sg=c["sg"],
                    mp_nsites=len(d["structure"]["sites"]), cod_nsites=c["nsites"],
                    supercell_sites=det.get("supercell_sites"),
                    max_occ_dev=det.get("max_site_occupancy_deviation"),
                    theoretical=d.get("theoretical"), e_above_hull=d.get("energy_above_hull"),
                ))
            if i % 500 == 0:
                print(i, len(jobs), dict(stats), f"{time.time() - t0:.0f}s", flush=True)
    rows.sort(key=lambda r: (r["chemsys"], r["cod_id"], r["mp_id"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["mp_id"])
        w.writeheader()
        w.writerows(rows)
    print("done", dict(stats), "links", len(rows), f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
