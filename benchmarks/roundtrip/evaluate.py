#!/usr/bin/env python
"""Step 4: match every ordered structure back to COD and score recall / false positives.

For each structure (unrelaxed orderings, MLIP-relaxed orderings, controls):

1. ``DisorderMatcher.match(structure, parent)`` -> recall of the parent COD entry.
2. ``DisorderMatcher.match(structure, decoy)`` for every *decoy*: the other COD entries of the
   same chemical system that ``NoveltyChecker`` would use as layer-2 references (COD
   ``has disorder`` flag or fractional formula).  Each decoy is classified once against the
   parent with ``DisorderMatcher.match_disordered``:
     equivalent  same disordered phase (another refinement / temperature / sample): a
                 DISORDERED_REALIZATION match to it is correct;
     family      same framework and site pattern, different stoichiometry: a realization
                 match is a "soft" false positive (within tolerances of a neighbouring
                 member of a solid-solution series);
     unrelated   different framework: a realization match is a hard false positive.
   Decoys whose composition is outside twice MatterGen's global tolerance from the query
   cannot give a realization and are skipped (counted as negatives) unless --all-decoys.
3. Optionally (``--novelty``) the full ``NoveltyChecker`` on a private working copy of the
   index (entries/keys copied from the read-only COD index, CIFs only for the candidates).

Outputs in --out-dir: results.jsonl (one record per structure and source), decoys.jsonl
(decoy classification), summary.md and summary_*.csv (metrics per bucket).

usage:
  python evaluate.py --selected selected.csv --structures orderings.jsonl relaxed_*.jsonl \
      --cod-db cod.sqlite --out-dir results/ --procs 8 [--novelty] [--shard i/N]
  python evaluate.py --out-dir results/ --summarize-only
"""

from __future__ import annotations

import argparse
import collections
import csv
import glob
import logging
import sqlite3
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

log = logging.getLogger("evaluate")

_G: dict = {}  # per-process state


def _init_worker(cfg: dict):
    warnings.filterwarnings("ignore")
    from opendisorder.matching import DisorderMatcher

    _G.clear()
    _G["cfg"] = cfg
    _G["conn"] = common.cod_connect(cfg["cod_db"])
    _G["matcher"] = DisorderMatcher()
    _G["structs"] = {}
    if cfg.get("work_index"):
        from opendisorder import CODIndex, NoveltyChecker
        from opendisorder.checker import CheckOptions

        _G["index"] = CODIndex(cfg["work_index"])
        _G["checker"] = NoveltyChecker(
            _G["index"], options=CheckOptions(time_limit=cfg["timeout"], exhaustive=cfg["exhaustive"])
        )


def _cod(cid: str):
    st = _G["structs"]
    if cid not in st:
        try:
            st[cid] = common.cod_structure(_G["conn"], cid)
        except Exception:
            st[cid] = None
    return st[cid]


def _match(q, ref) -> tuple[str, dict, float]:
    t0 = time.time()
    try:
        with common.time_limit(_G["cfg"]["timeout"]):
            res = _G["matcher"].match(q, ref)
        return res.kind.value, res.details, time.time() - t0
    except common.Timeout:
        return "timeout", {}, time.time() - t0
    except Exception as exc:
        return "error", {"reason": f"{type(exc).__name__}: {str(exc)[:120]}"}, time.time() - t0


# ----------------------------------------------------------------------------- decoys


def classify_decoys(job) -> list[dict]:
    """Relation of every same-chemsys disordered COD entry to the parent."""
    parent_id, chemsys, max_cell_atoms, max_decoys = job
    parent = _cod(parent_id)
    out = []
    ids = [c for c in common.same_chemsys_disordered_ids(_G["conn"], chemsys, max_cell_atoms) if c != parent_id]
    ids = ids[:max_decoys] if max_decoys else ids
    for cid in ids:
        s = _cod(cid)
        rec = {"parent": parent_id, "decoy": cid}
        if s is None:
            rec["relation"] = "unparsable"
        elif s.is_ordered:
            rec["relation"] = "ordered"  # not a layer-2 reference once parsed
        else:
            rec["decoy_formula"] = s.composition.reduced_formula
            try:
                with common.time_limit(_G["cfg"]["timeout"]):
                    r = _G["matcher"].match_disordered(parent, s)
                    raw = {"pristine": "equivalent", "disordered_family": "family"}.get(r.kind.value, "unrelated")
                    rec["relation_raw"] = raw
                    if raw == "unrelated" and _G["cfg"]["trace"] > 0:
                        # trace occupancies (e.g. 0.75% Ti on an otherwise empty site) add
                        # sites that break match_disordered between two refinements of the
                        # same phase; compare again with them removed
                        p2, s2 = drop_trace(parent, _G["cfg"]["trace"]), drop_trace(s, _G["cfg"]["trace"])
                        if p2 is not None and s2 is not None:
                            r = _G["matcher"].match_disordered(p2, s2)
                            raw = {"pristine": "equivalent", "disordered_family": "family"}.get(r.kind.value, "unrelated")
                rec["relation"] = raw
            except common.Timeout:
                rec["relation"] = "unrelated_timeout"
            except Exception as exc:
                rec["relation"] = "unrelated_error"
                rec["error"] = str(exc)[:120]
        out.append(rec)
    return out


def drop_trace(structure, tol: float):
    """Copy without species of occupancy < tol (sites left empty are removed)."""
    from pymatgen.core import Structure

    species, coords = [], []
    for site in structure:
        sp = {el: x for el, x in site.species.items() if x >= tol}
        if sp:
            species.append(sp)
            coords.append(site.frac_coords)
    if not species:
        return None
    return Structure(structure.lattice, species, coords)


# ----------------------------------------------------------------------------- structures


def evaluate_group(job) -> list[dict]:
    """All structures of one parent: match against the parent and the decoys."""
    from opendisorder.matching import fractional_composition_close
    from pymatgen.core import Structure

    parent_id, records, decoys = job
    cfg = _G["cfg"]
    parent = _cod(parent_id)
    out = []
    for rec in records:
        row = {k: v for k, v in rec.items() if k != "structure"}
        row["source"] = source_of(rec)
        if rec.get("error"):
            row["kind"] = "relax_failed"
            out.append(row)
            continue
        q = Structure.from_dict(rec["structure"])
        if parent is None:
            row["kind"] = "parent_unavailable"
            out.append(row)
            continue
        kind, det, dt = _match(q, parent)
        row.update(
            kind=kind,
            match_time_s=round(dt, 3),
            max_sublattice_dev=det.get("max_sublattice_occupancy_deviation"),
            max_site_dev=det.get("max_site_occupancy_deviation"),
            idealized_sublattices=det.get("idealized_sublattices"),
            reason=det.get("reason"),
        )
        decoy_kinds = {}
        n_skipped = 0
        t0 = time.time()
        for d in decoys:
            if d["relation"] in ("ordered", "unparsable"):
                continue
            ds = _cod(d["decoy"])
            if ds is None:
                continue
            if not cfg["all_decoys"] and not fractional_composition_close(q.composition, ds.composition, 0.02, 0.2):
                n_skipped += 1
                decoy_kinds[d["decoy"]] = "skipped_composition"
                continue
            decoy_kinds[d["decoy"]] = _match(q, ds)[0]
        row["decoy_time_s"] = round(time.time() - t0, 3)
        row["decoy_kinds"] = decoy_kinds
        rel = {d["decoy"]: d["relation"].split("_")[0] for d in decoys}
        for relation in ("equivalent", "family", "unrelated"):
            ks = [k for c, k in decoy_kinds.items() if rel.get(c) == relation]
            row[f"n_decoys_{relation}"] = len(ks)
            row[f"n_real_{relation}"] = sum(k == "disordered_realization" for k in ks)
            row[f"n_fam_{relation}"] = sum(k == "disordered_family" for k in ks)
        row["n_decoy_timeouts"] = sum(k in ("timeout", "error") for k in decoy_kinds.values())
        if "checker" in _G:
            row.update(run_novelty(q, parent_id, rel))
        out.append(row)
    return out


def run_novelty(q, parent_id: str, rel: dict) -> dict:
    t0 = time.time()
    try:
        rep = _G["checker"].check(q, fetch=True)
    except Exception as exc:
        return {"novelty_kind": "error", "novelty_error": str(exc)[:200]}
    ids = rep.cod_ids
    real_ids = [m.reference_id for m in rep.matches if m.kind.value in ("pristine", "disordered_realization")]
    top = rep.best()
    top_rel = None
    if top is not None and top.reference_id:
        top_rel = "parent" if top.reference_id == parent_id else rel.get(top.reference_id, "outside_decoys").split("_")[0]
    return {
        "novelty_kind": rep.kind.value,
        "novelty_is_novel": rep.is_novel,
        "novelty_top_id": top.reference_id if top else None,
        "novelty_top_relation": top_rel,
        "novelty_parent_found": parent_id in ids,
        "novelty_known_ids": real_ids[:20],
        "novelty_errors": len(rep.errors),
        "novelty_time_s": round(time.time() - t0, 2),
    }


def source_of(rec: dict) -> str:
    if rec.get("relax_backend"):
        return f"{rec['relax_backend']}:{rec.get('relax_model', '')}"
    return "unrelaxed"


# ----------------------------------------------------------------------------- work index for NoveltyChecker


def build_work_index(cod_db: str, work: str, structures) -> None:
    """Private, writable copy of the index: all metadata, CIFs/parses only for the candidates."""
    from opendisorder import CODIndex, NoveltyChecker
    from opendisorder.index import SCHEMA
    from opendisorder.matching import clean_structure

    src = str(Path(cod_db).expanduser().resolve())
    Path(work).parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(work)
    db.executescript(SCHEMA)
    db.execute(f"ATTACH DATABASE 'file:{src}?mode=ro' AS src")
    if db.execute("SELECT COUNT(*) FROM entries").fetchone()[0] == 0:
        log.info("copying COD metadata into %s", work)
        with db:
            db.execute("INSERT OR IGNORE INTO entries SELECT * FROM src.entries")
            db.execute("INSERT OR IGNORE INTO keys SELECT * FROM src.keys")
    db.commit()
    index = CODIndex(work)
    checker = NoveltyChecker(index)
    ids: set[str] = set()
    for s in structures:
        ids.update(checker.candidate_ids(clean_structure(s)).all())
    have = {r[0] for r in db.execute("SELECT cod_id FROM cifs")}
    need = sorted(ids - have)
    log.info("work index: %d candidate COD entries, copying %d CIFs", len(ids), len(need))
    for i in range(0, len(need), 500):
        chunk = need[i : i + 500]
        marks = ",".join("?" * len(chunk))
        with db:
            db.execute(f"INSERT OR IGNORE INTO cifs SELECT * FROM src.cifs WHERE cod_id IN ({marks})", chunk)
            db.execute(f"INSERT OR IGNORE INTO structures SELECT * FROM src.structures WHERE cod_id IN ({marks})", chunk)
    db.execute("DETACH DATABASE src")
    db.close()
    index.parse_missing(sorted(ids))  # parse in this process, before the workers only read
    index.close()


# ----------------------------------------------------------------------------- summary


def _rate(num, den):
    return round(num / den, 4) if den else None


def summarize(out_dir: Path):
    import pandas as pd

    rows = list(common.iter_jsonl(sorted(out_dir.glob("results*.jsonl"))))
    if not rows:
        log.warning("no results in %s", out_dir)
        return
    df = pd.DataFrame(rows)
    for c in ("n_real_equivalent", "n_real_family", "n_real_unrelated", "n_decoys_equivalent", "n_decoys_family",
              "n_decoys_unrelated", "n_fam_unrelated"):
        if c not in df:
            df[c] = 0
        df[c] = df[c].fillna(0)
    df["is_real"] = df["kind"] == "disordered_realization"
    df["is_any"] = df["kind"].isin(["disordered_realization", "disordered_family"])
    df["fp_hard"] = df["n_real_unrelated"] > 0
    df["fp_soft"] = df["n_real_family"] > 0
    df["n_atoms_bin"] = pd.cut(df["n_atoms"], [0, 20, 40, 80, 160, 10**6], labels=["<=20", "21-40", "41-80", "81-160", ">160"])
    evaluated = ~df["kind"].isin(["relax_failed", "parent_unavailable"])

    def metrics(g):
        pos = g[(g["expected"] == "realization")]
        neg = g[(g["expected"] == "not_realization")]
        tp = pos["is_real"].sum() + pos["n_real_equivalent"].sum()
        fp_strict = pos["n_real_family"].sum() + pos["n_real_unrelated"].sum()
        d = {
            "n_pos": len(pos),
            "recall_realization": _rate(pos["is_real"].sum(), len(pos)),
            "recall_any": _rate(pos["is_any"].sum(), len(pos)),
            "family_rate": _rate((pos["kind"] == "disordered_family").sum(), len(pos)),
            "none_rate": _rate((pos["kind"] == "none").sum(), len(pos)),
            "timeout_rate": _rate(pos["kind"].isin(["timeout", "error"]).sum(), len(pos)),
            "struct_fp_hard": _rate(pos["fp_hard"].sum(), len(pos)),
            "struct_fp_soft": _rate(pos["fp_soft"].sum(), len(pos)),
            "pair_fp_unrelated": _rate(pos["n_real_unrelated"].sum(), pos["n_decoys_unrelated"].sum()),
            "precision_pairs_strict": _rate(tp, tp + fp_strict),
            "precision_pairs_lenient": _rate(tp + pos["n_real_family"].sum(), tp + fp_strict),
            "n_neg": len(neg),
            "specificity_controls": _rate((~neg["is_real"]).sum(), len(neg)),
        }
        if "novelty_kind" in g:
            nv = pos[pos["novelty_kind"].notna()]
            d["novelty_known"] = _rate((nv["novelty_is_novel"] == False).sum(), len(nv))  # noqa: E712
            d["novelty_top_parent_or_equiv"] = _rate(nv["novelty_top_relation"].isin(["parent", "equivalent"]).sum(), len(nv))
            nn = neg[neg["novelty_kind"].notna()]
            d["novelty_neg_known"] = _rate((nn["novelty_is_novel"] == False).sum(), len(nn))  # noqa: E712
        return pd.Series(d)

    buckets = ["source", "method", "size", "n_atoms_bin", "disorder_type", "structure_class", "kind_ctrl", "exact"]
    # the generator's record kind (positive / ctrl_*) is kept as 'record_kind'; 'kind' is the match kind
    df["kind_ctrl"] = df["record_kind"].fillna("positive") if "record_kind" in df else "positive"
    lines = ["# Round-trip benchmark summary", "", f"{len(df)} evaluated (structure, source) pairs; "
             f"{df['parent'].nunique()} parents; sources: {', '.join(sorted(df['source'].unique()))}", ""]
    lines += ["## Match outcome vs parent, by record kind and source", "",
              pd.crosstab([df["kind_ctrl"], df["source"]], df["kind"]).to_markdown(), ""]
    base = df[evaluated]
    tables = {"overall": base.groupby(lambda _: "all").apply(metrics)}
    for b in buckets:
        if b in base and base[b].notna().any():
            tables[b] = base.groupby(["source", b], observed=True).apply(metrics) if b != "source" else base.groupby(b).apply(metrics)
    key = ["n_pos", "recall_realization", "recall_any", "struct_fp_hard", "struct_fp_soft", "precision_pairs_strict",
           "n_neg", "specificity_controls", "novelty_known", "novelty_top_parent_or_equiv", "novelty_neg_known"]
    for name, t in tables.items():
        t.to_csv(out_dir / f"summary_{name}.csv")  # all columns
        lines += [f"## By {name}", "", t[[c for c in key if c in t]].to_markdown(floatfmt=".3f"), ""]
    # paired robustness: realization before relaxation -> after relaxation
    if "source_id" in base:
        unrel = base[base["source"] == "unrelaxed"].set_index("id")["is_real"]
        rel = base[base["source"] != "unrelaxed"]
        if len(rel):
            rel = rel.assign(before=rel["source_id"].map(unrel))
            paired = rel[rel["before"] == True]  # noqa: E712
            t = paired.groupby(["source", "disorder_type"]).agg(n=("is_real", "size"), kept_realization=("is_real", "mean"))
            t.to_csv(out_dir / "summary_relax_robustness.csv")
            lines += ["## Realization kept after relaxation (paired, structures matched before)", "", t.to_markdown(), ""]
            if "converged" in rel:
                t = rel.groupby("source").agg(n=("id", "size"), converged=("converged", "mean"),
                                              median_steps=("n_steps", "median"), median_time_s=("time_s", "median"),
                                              median_volume_change=("volume_change", "median"))
                lines += ["## Relaxation statistics", "", t.to_markdown(), ""]
    t = base.groupby("source").agg(median_match_s=("match_time_s", "median"), median_decoy_s=("decoy_time_s", "median"))
    lines += ["## Timing", "", t.to_markdown(), ""]
    (out_dir / "summary.md").write_text("\n".join(lines))
    log.info("summary written to %s", out_dir / "summary.md")


# ----------------------------------------------------------------------------- main


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selected")
    ap.add_argument("--structures", nargs="+", help="orderings / relaxed JSON-lines files (globs allowed)")
    ap.add_argument("--cod-db", default=common.DEFAULT_COD_DB)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=60.0, help="seconds per structure comparison")
    ap.add_argument("--max-decoys", type=int, default=100, help="per parent (0 = all)")
    ap.add_argument("--decoy-max-cell-atoms", type=float, default=400.0)
    ap.add_argument("--all-decoys", action="store_true", help="no composition prefilter on decoys")
    ap.add_argument("--trace", type=float, default=0.05, help="occupancies below this are ignored when relating decoys to the parent")
    ap.add_argument("--kinds", help="only these record kinds (positive,ctrl_rattle,...)")
    ap.add_argument("--novelty", action="store_true", help="also run the full NoveltyChecker")
    ap.add_argument("--exhaustive", action="store_true", help="NoveltyChecker: collect every match")
    ap.add_argument("--work-index", help="writable working index for --novelty (default: <out-dir>/work_index.sqlite)")
    ap.add_argument("--shard", help="i/N or slurm/N: subset of parents")
    ap.add_argument("--summarize-only", action="store_true")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    warnings.filterwarnings("ignore")
    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if a.summarize_only:
        summarize(out_dir)
        return
    if not (a.selected and a.structures):
        ap.error("--selected and --structures are required")

    shard, nshards = common.parse_shard(a.shard)
    suffix = f".shard{shard}of{nshards}" if nshards > 1 else ""
    res_path = out_dir / f"results{suffix}.jsonl"
    dec_path = out_dir / f"decoys{suffix}.jsonl"
    with open(a.selected, newline="") as fh:
        parents = {r["cod_id"]: r for r in csv.DictReader(fh) if common.stable_shard(r["cod_id"], nshards) == shard}
    files = [f for pat in a.structures for f in (sorted(glob.glob(pat)) or [pat])]
    done = {(r["id"], r["source"]) for r in common.iter_jsonl(res_path)}
    kinds = set(a.kinds.split(",")) if a.kinds else None
    by_parent = collections.defaultdict(list)
    n = 0
    for rec in common.iter_jsonl(files):
        if rec["parent"] not in parents or (rec["id"], source_of(rec)) in done:
            continue
        if kinds and rec.get("kind") not in kinds:
            continue
        rec["record_kind"] = rec.get("kind", "positive")
        by_parent[rec["parent"]].append(rec)
        n += 1
    log.info("%d structures of %d parents to evaluate (%d done)", n, len(by_parent), len(done))

    cfg = {"cod_db": a.cod_db, "timeout": a.timeout, "all_decoys": a.all_decoys, "exhaustive": a.exhaustive, "trace": a.trace}
    if a.novelty:
        cfg["work_index"] = a.work_index or str(out_dir / f"work_index{suffix}.sqlite")
        from pymatgen.core import Structure

        build_work_index(a.cod_db, cfg["work_index"], (Structure.from_dict(r["structure"]) for rs in by_parent.values() for r in rs if not r.get("error")))

    # ---- decoys (classified once per parent, cached)
    decoys = collections.defaultdict(list)
    for d in common.iter_jsonl(dec_path):
        decoys[d["parent"]].append(d)
    todo = [(p, parents[p].get("chemsys", ""), a.decoy_max_cell_atoms, a.max_decoys) for p in by_parent if p not in decoys]
    pool = Pool(a.procs, initializer=_init_worker, initargs=(cfg,)) if a.procs > 1 else None
    if pool is None:
        _init_worker(cfg)
    mapper = (lambda f, it: pool.imap_unordered(f, it)) if pool else map
    with common.JsonlWriter(dec_path) as w:
        for recs in mapper(classify_decoys, todo):
            for d in recs:
                w.write(d)
                decoys[d["parent"]].append(d)
    rel_counts = collections.Counter(d["relation"] for ds in decoys.values() for d in ds)
    log.info("decoys: %s", dict(rel_counts))

    # ---- structures (one job per parent, largest first)
    jobs = sorted(((p, rs, decoys.get(p, [])) for p, rs in by_parent.items()), key=lambda j: -len(j[1]))
    t0 = time.time()
    with common.JsonlWriter(res_path) as w:
        for i, rows in enumerate(mapper(evaluate_group, jobs)):
            for r in rows:
                w.write(r)
            log.info("[%d/%d] parent %s: %d structures, %s (%.0fs)", i + 1, len(jobs), rows[0]["parent"] if rows else "?",
                     len(rows), dict(collections.Counter(r["kind"] for r in rows)), time.time() - t0)
    if pool:
        pool.close()
    summarize(out_dir)


if __name__ == "__main__":
    main()
