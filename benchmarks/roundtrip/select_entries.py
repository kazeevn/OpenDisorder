#!/usr/bin/env python
"""Step 1: select high-quality, compositionally disordered COD entries as round-trip parents.

Two stages (see README, "Quality criteria"):

A. metadata (COD CSV dump, streamed): has coordinates; COD 'has disorder' flag or a
   fractional formula; not retracted / not flagged with errors or warnings; not theoretical;
   no H/D and only elements the MLIPs cover; 2-6 elements; R1 (Robs, else Rall) <= --max-r;
   not a powder refinement; <= --max-cell-atoms atoms in the conventional cell; ambient
   temperature and pressure (or unspecified).
B. CIF (read-only from the opendisorder index): no powder (_pd_*) tags; no
   _atom_site_disorder_group/assembly (split positions); parses with pymatgen without
   warnings at occupancy tolerance --occupancy-tolerance; genuinely disordered after
   parsing; no site overfull; no two partially occupied sites closer than
   --split-factor x (r_i + r_j) (positional / exclusion disorder); every occupancy a
   fraction with denominator <= --max-denominator (within --rational-tol); the smallest
   exact supercell has <= --max-supercell primitive cells and <= --max-supercell-atoms sites;
   minority occupancy >= --min-occupancy.

Then at most --max-per-chemsys entries per chemical system (lowest R first) and, optionally,
a stratified sample (--limit) over (disorder type x structure class).

usage:
  python select_entries.py --csv all.csv --cod-db cod.sqlite --out selected.csv [--limit 2000]
  python select_entries.py --cod-db cod.sqlite --ids 1000014,9008565 --out smoke.csv   # skip stage A
"""

from __future__ import annotations

import argparse
import ast
import collections
import csv
import json
import logging
import random
import re
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

log = logging.getLogger("select")

# Elements beyond Pu are not covered by the MPtrj/OMat24-trained MLIPs; noble gases are rare in
# COD inorganic entries and poorly represented in training data.  H/D: hydrogen positions from
# X-ray are unreliable and H-containing (often molecular) crystals are not the target here.
DEFAULT_EXCLUDE = "H,D,T,He,Ne,Ar,Kr,Xe,Rn,Po,At,Fr,Ra,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr"


def _float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def stage_a(row: dict, a) -> str | None:
    """Return a rejection reason, or None if the metadata row passes."""
    flags = row.get("flags") or ""
    if "has coordinates" not in flags:
        return "no_coordinates"
    formula = common.cod_formula_elements(row.get("formula")) or common.cod_formula_elements(row.get("calcformula"))
    if not formula:
        return "no_formula"
    fractional = any(abs(v - round(v)) > 1e-6 for v in formula.values())
    if "has disorder" not in flags and not fractional:
        return "not_disordered"
    status = (row.get("status") or "").strip().lower()
    if status in ("retracted", "errors", "warnings"):
        return f"status_{status}"
    if (row.get("method") or "").strip().lower() == "theoretical":
        return "theoretical"
    els = set(formula)
    if els & a.exclude:
        return "excluded_element"
    if not a.min_elements <= len(els) <= a.max_elements:
        return "n_elements"
    r = _float(row.get("Robs"))
    if r is None:
        r = _float(row.get("Rall"))
    if r is None:
        if not a.allow_missing_r:
            return "no_R"
    elif r > a.max_r:
        return "R_high"
    if (row.get("method") or "").strip().lower().startswith("powder") and not a.allow_powder:
        return "powder"
    cell = common.cod_formula_elements(row.get("cellformula"))
    natoms = sum(cell.values()) if cell else None
    if natoms is None or natoms > a.max_cell_atoms:
        return "cell_too_large"
    for key in ("celltemp", "diffrtemp"):
        t = _float(row.get(key))
        if t is not None and not a.tmin <= t <= a.tmax:
            return "non_ambient_T"
    for key in ("cellpressure", "diffrpressure"):
        p = _float(row.get(key))  # kPa in COD
        if p is not None and p > a.max_pressure_kpa:
            return "high_pressure"
    return None


def iter_csv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        yield from csv.DictReader(line for line in fh if not line.startswith("#"))


def meta_record(row: dict) -> dict:
    r = _float(row.get("Robs"))
    return {
        "cod_id": str(row["file"]),
        "cod_formula": (row.get("formula") or "").strip("- ").strip(),
        "sg_number": row.get("sgNumber"),
        "R": r if r is not None else _float(row.get("Rall")),
        "R_kind": "Robs" if r is not None else "Rall",
        "temperature": _float(row.get("diffrtemp")) or _float(row.get("celltemp")),
        "natoms_cell": sum(common.cod_formula_elements(row.get("cellformula")).values()) or None,
        "mineral": row.get("mineral") or "",
        "year": row.get("year") or "",
        "journal": row.get("journal") or "",
        "doi": row.get("doi") or "",
    }


# ----------------------------------------------------------------------------- stage B

_DICT = re.compile(r"(CIF|PMG)=(\{[^}]*\})")
_OCC = re.compile(r"Some occupancies \(\[([^\]]*)\]\) sum to > 1")


def benign_warning(msg: str, stoich_tol: float) -> bool:
    """pymatgen warnings that do not indicate a defective CIF:
    - deprecation notices;
    - 'Incorrect stoichiometry' when the CIF's _chemical_formula_sum and the parsed
      composition agree within ``stoich_tol`` in atomic fraction (COD formulas of
      disordered phases are often rounded);
    - 'occupancies sum to > 1' by floating-point noise (<= 1.005)."""
    if msg.startswith(common.BENIGN_WARNINGS):
        return True
    if msg.startswith("Incorrect stoichiometry"):
        found = dict(_DICT.findall(msg))
        try:
            cif, pmg = ast.literal_eval(found["CIF"]), ast.literal_eval(found["PMG"])
        except Exception:
            return False
        tc, tp = sum(cif.values()), sum(pmg.values())
        if not tc or not tp:
            return False
        return all(abs(cif.get(e, 0) / tc - pmg.get(e, 0) / tp) <= stoich_tol for e in set(cif) | set(pmg))
    m = _OCC.search(msg)
    if m:
        try:
            return max(float(x) for x in m.group(1).split(",")) <= 1.005
        except ValueError:
            return False
    return False


_conn = None
_args = None


def _init(cod_db, args):
    global _conn, _args
    warnings.filterwarnings("ignore")
    _conn = common.cod_connect(cod_db)
    _args = args


def stage_b(meta: dict) -> dict:
    """CIF-level checks; returns meta updated with 'reject' (reason or '') and analysis columns."""
    from opendisorder.index import parse_cif

    a = _args
    cid = meta["cod_id"]
    out = dict(meta, reject="")
    text = common.cod_cif(_conn, cid)
    if text is None:
        out["reject"] = "no_cif"
        return out
    flags = common.cif_flags(text)
    out.update(flags)
    if flags["powder_tags"] and not a.allow_powder:
        out["reject"] = "powder_cif"
        return out
    if flags["disorder_groups"] and not a.allow_disorder_groups:
        out["reject"] = "positional_disorder_tags"
        return out
    try:
        with common.time_limit(a.timeout):
            parsed = parse_cif(cid, text, occupancy_tolerance=a.occupancy_tolerance)
    except common.Timeout:
        out["reject"] = "parse_timeout"
        return out
    msgs = [m for m in parsed.warnings if not benign_warning(m, a.stoich_tol)]
    s = parsed.structure
    if s is None:
        out["reject"] = "parse_error"
        out["parse_error"] = (parsed.error or "")[:200]
        return out
    if msgs and not a.allow_warnings:
        out["reject"] = "parse_warnings"
        out["parse_warnings"] = " | ".join(m.replace("\n", " ")[:120] for m in msgs)[:400]
        return out
    if s.is_ordered:
        out["reject"] = "parsed_ordered"
        return out
    els = sorted({el.symbol for el in s.composition.elements})
    if set(els) & a.exclude:
        out["reject"] = "excluded_element"
        return out
    try:
        with common.time_limit(a.timeout):
            info = common.analyze_disorder(
                s, max_den=a.max_denominator, occ_tol=a.rational_tol, split_factor=a.split_factor
            )
    except common.Timeout:
        out["reject"] = "analysis_timeout"
        return out
    except Exception as exc:  # pragma: no cover - pymatgen corner cases
        out["reject"] = "analysis_error"
        out["parse_error"] = str(exc)[:200]
        return out
    comp = s.composition
    n_occ_prim = comp.num_atoms * info["nsites_prim"] / len(s)  # occupied atoms per primitive cell
    out.update(
        {
            "reduced_formula": comp.reduced_formula,
            "chemsys": "-".join(els),
            "nelements": len(els),
            "structure_class": common.structure_class(els),
            "nsites_conv": len(s),
            "nsites_prim": info["nsites_prim"],
            "atoms_prim": round(n_occ_prim, 3),
            "n_orbits": info["n_orbits"],
            "n_disordered_orbits": info["n_disordered_orbits"],
            "disorder_type": info["disorder_type"],
            "min_supercell": info["min_supercell"],
            "max_rational_error": info["max_rational_error"],
            "min_minority_occupancy": info["min_minority_occupancy"],
            "min_partial_distance": info["min_partial_distance"],
            "min_partial_distance_ratio": info["min_partial_distance_ratio"],
            "orbits": json.dumps(info["orbits"]),
        }
    )
    try:
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        out["crystal_system"] = SpacegroupAnalyzer(s, symprec=0.01).get_crystal_system()
    except Exception:
        out["crystal_system"] = ""
    if info["disorder_type"] == "overfull":
        out["reject"] = "overfull_site"
    elif info["disorder_type"] == "ordered":
        out["reject"] = "parsed_ordered"
    elif info["split_site"] and not a.allow_split_sites:
        out["reject"] = "split_site"
    elif not info["rational_ok"]:
        out["reject"] = "irrational_occupancy"
    elif info["min_minority_occupancy"] < a.min_occupancy:
        out["reject"] = "dilute"
    elif info["min_supercell"] > a.max_supercell or info["min_supercell"] * info["nsites_prim"] > a.max_supercell_atoms:
        out["reject"] = "supercell_too_large"
    if not out["reject"]:
        out["n_same_chemsys_disordered"] = len(common.same_chemsys_disordered_ids(_conn, out["chemsys"])) - 1
    return out


# ----------------------------------------------------------------------------- main


def stratified_sample(rows: list[dict], n: int, seed: int) -> list[dict]:
    rng = random.Random(seed)
    buckets = collections.defaultdict(list)
    for r in rows:
        buckets[(r["disorder_type"], r["structure_class"])].append(r)
    for b in buckets.values():
        rng.shuffle(b)
    keys = sorted(buckets)
    out = []
    while len(out) < n and any(buckets.values()):
        for k in keys:
            if buckets[k] and len(out) < n:
                out.append(buckets[k].pop())
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", help="COD CSV metadata dump (stage A); omit with --ids")
    ap.add_argument("--cod-db", default=common.DEFAULT_COD_DB, help="opendisorder COD index (opened read-only)")
    ap.add_argument("--out", required=True, help="selected entries (CSV)")
    ap.add_argument("--rejected", help="optional CSV with every stage-B rejection and its reason")
    ap.add_argument("--ids", help="comma-separated COD ids: skip stage A, run stage B on these only")
    ap.add_argument("--max-r", type=float, default=0.05)
    ap.add_argument("--allow-missing-r", action="store_true")
    ap.add_argument("--allow-powder", action="store_true")
    ap.add_argument("--allow-warnings", action="store_true")
    ap.add_argument("--allow-disorder-groups", action="store_true")
    ap.add_argument("--allow-split-sites", action="store_true")
    ap.add_argument("--exclude", default=DEFAULT_EXCLUDE)
    ap.add_argument("--min-elements", type=int, default=2)
    ap.add_argument("--max-elements", type=int, default=6)
    ap.add_argument("--max-cell-atoms", type=float, default=80, help="atoms in the conventional cell (COD cellformula)")
    ap.add_argument("--tmin", type=float, default=273.0)
    ap.add_argument("--tmax", type=float, default=313.0)
    ap.add_argument("--max-pressure-kpa", type=float, default=1000.0)
    ap.add_argument("--occupancy-tolerance", type=float, default=1.02)
    ap.add_argument("--stoich-tol", type=float, default=0.01, help="atomic-fraction tolerance, CIF formula vs parse")
    ap.add_argument("--max-denominator", type=int, default=12)
    ap.add_argument("--rational-tol", type=float, default=0.011)
    ap.add_argument("--min-occupancy", type=float, default=0.08)
    ap.add_argument("--split-factor", type=float, default=0.7)
    ap.add_argument("--max-supercell", type=int, default=8, help="primitive cells in the smallest exact supercell")
    ap.add_argument("--max-supercell-atoms", type=int, default=160)
    ap.add_argument("--max-per-chemsys", type=int, default=2)
    ap.add_argument("--limit", type=int, default=0, help="stratified sample of this size (0 = all)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-candidates", type=int, default=0, help="run stage B on a random subset (quick tests)")
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=60.0, help="seconds per CIF for parsing/analysis")
    a = ap.parse_args(argv)
    a.exclude = set(filter(None, a.exclude.split(",")))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    # ---- stage A
    stats_a = collections.Counter()
    metas = []
    if a.ids:
        metas = [{"cod_id": i.strip()} for i in a.ids.split(",") if i.strip()]
    else:
        if not a.csv:
            ap.error("--csv or --ids is required")
        for row in iter_csv(a.csv):
            reason = stage_a(row, a)
            stats_a[reason or "pass"] += 1
            if reason is None:
                metas.append(meta_record(row))
        log.info("stage A: %s", dict(stats_a.most_common()))
        if a.max_candidates and len(metas) > a.max_candidates:
            metas = random.Random(a.seed).sample(metas, a.max_candidates)
            log.info("stage B restricted to a random %d of the stage-A survivors", len(metas))

    # ---- stage B
    stats_b = collections.Counter()
    rows, rejected = [], []
    jobs = metas
    if a.procs > 1:
        pool = Pool(a.procs, initializer=_init, initargs=(a.cod_db, a))
        it = pool.imap_unordered(stage_b, jobs, chunksize=8)
    else:
        _init(a.cod_db, a)
        pool = None
        it = map(stage_b, jobs)
    for i, rec in enumerate(it):
        stats_b[rec["reject"] or "pass"] += 1
        (rejected if rec["reject"] else rows).append(rec)
        if i % 1000 == 999:
            log.info("stage B %d/%d %s", i + 1, len(jobs), dict(stats_b.most_common(8)))
    if pool:
        pool.close()
    log.info("stage B: %s", dict(stats_b.most_common()))

    # ---- diversity: at most N per chemical system, lowest R first
    rows.sort(key=lambda r: (r["chemsys"], r.get("R") if r.get("R") is not None else 1.0, r["cod_id"]))
    per_cs = collections.Counter()
    kept = []
    for r in rows:
        if a.max_per_chemsys and per_cs[r["chemsys"]] >= a.max_per_chemsys:
            continue
        per_cs[r["chemsys"]] += 1
        kept.append(r)
    if a.limit and len(kept) > a.limit:
        kept = stratified_sample(kept, a.limit, a.seed)
    kept.sort(key=lambda r: r["cod_id"])
    log.info(
        "selected %d entries (%d before per-chemsys cap); by type %s; by class %s",
        len(kept),
        len(rows),
        dict(collections.Counter(r["disorder_type"] for r in kept)),
        dict(collections.Counter(r["structure_class"] for r in kept)),
    )
    _write_csv(a.out, kept)
    if a.rejected:
        _write_csv(a.rejected, rejected)
    stats_path = Path(a.out).with_suffix(".stats.json")
    stats_path.write_text(json.dumps({"stage_a": stats_a, "stage_b": stats_b, "selected": len(kept)}, indent=1))


def _write_csv(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields or ["cod_id"])
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
