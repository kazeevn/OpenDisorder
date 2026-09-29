"""Loader for figshare 29145101 (WyFormer / DiffCSP / DiffCSP++ / WyCryst / SymmCD /
CrystalFormer / MiAD / MatterGen / Crys-JEPA / Compositional generated crystals).

Layout: data/<train_dataset>/<model>/[<post-processing>/...]/data.csv.gz
Each data.csv.gz is a pandas CSV; first column `material_id` is the row index
(int, preserved across relaxation stages unless datasets.yaml says
index_preserved: false). Column `structure` is a monty/pymatgen Structure JSON
(as_dict). Wyckoff-only files (WyckoffTransformer, WyckoffLLM-*, WyCryst,
UN-WyckoffTransformer, WyckoffReduction) have no `structure` column and are skipped.

Usage:
    from load_generated import iter_structures, list_sets
    for set_name, sample_id, s in iter_structures(raw_only=True):
        ...   # set_name e.g. "mp_20/DiffCSP", sample_id = int material_id
    # CLI: python load_generated.py [--raw-only] [--include-splits] [--stats]
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Iterator

import pandas as pd
from pymatgen.core import Structure

DATA_ROOT = Path(os.environ.get("FIGSHARE_DATA", Path(__file__).resolve().parent / "data"))

# Directory names that denote a post-processing step (relaxation / DFT) rather
# than a generative model.
_POST = ("CHGNet", "DFT", "eq-V2", "CrySPR")

# Canonical model name for each "raw" generated set (structures as produced by
# the generator, before any ML/DFT relaxation).
RAW_SETS = {
    "mp_20/DiffCSP": "DiffCSP",
    "mp_20/DiffCSP++": "DiffCSP++",
    "mp_20/CrystalFormer": "CrystalFormer",
    "mp_20/SymmCD": "SymmCD",
    "mp_20/MiAD": "MiAD",
    "mp_20/MatterGen/MatterGen_10k": "MatterGen",
    "mp_20/Compositional": "Compositional",
    "mp_20/Crys-JEPA": "Crys-JEPA",
    "mp_20/WyckoffTransformer/DiffCSP++": "WyFormer+DiffCSP++",
    "mp_20/WyckoffTransformer/DiffCSP++10k": "WyFormer+DiffCSP++ (10k)",
    "mp_20/WyckoffTransformer/CrySPR/CHGNet_fix": "WyFormer+CrySPR",
    "mp_20/WyckoffTransformer-letters/DiffCSP++": "WyFormer-letters+DiffCSP++",
    "mp_20/WyckoffTransformer-harmonic/DiffCSP++": "WyFormer-harmonic+DiffCSP++",
    "mp_20/WyckoffLLM-naive/DiffCSP++": "WyckoffLLM-naive+DiffCSP++",
    "mp_20/WyckoffLLM-site-symmetry/DiffCSP++": "WyckoffLLM-site-symmetry+DiffCSP++",
    "mp_20/WyckoffLLM-vanilla/DiffCSP++": "WyckoffLLM-vanilla+DiffCSP++",
    "mp_20_biternary/WyCryst/CrySPR/CHGNet_fix": "WyCryst+CrySPR",
    "mp_20_biternary/WyckoffTransformer/CrySPR/CHGNet_fix": "WyFormer+CrySPR",
    "mpts_52/SymmCD": "SymmCD",
    "mpts_52/WyckoffTransformer/CrySPR/CHGNet_fix": "WyFormer+CrySPR",
}


def list_sets(include_splits: bool = False, raw_only: bool = False,
              include_un: bool = False) -> list[str]:
    """Relative set names (directory under data/) that contain structures."""
    out = []
    for f in sorted(DATA_ROOT.rglob("data.csv.gz")):
        name = str(f.parent.relative_to(DATA_ROOT))
        is_split = "/split/" in name
        if is_split and not include_splits:
            continue
        if "/UN-" in name and not include_un:  # DiffCSP++ on *train* templates (baseline)
            continue
        if raw_only and name not in RAW_SETS and not is_split:
            continue
        cols = pd.read_csv(f, nrows=0).columns
        if "structure" in cols:
            out.append(name)
    return out


def iter_structures(sets: list[str] | None = None, *, include_splits: bool = False,
                    raw_only: bool = False, include_un: bool = False,
                    with_row: bool = False, errors: list | None = None
                    ) -> Iterator[tuple]:
    """Yield (set_name, sample_id, Structure[, row: pd.Series]).

    Parse failures are appended to `errors` as (set_name, sample_id, repr(exc))
    and skipped; the `cif` column is used as a fallback if present.
    """
    if sets is None:
        sets = list_sets(include_splits=include_splits, raw_only=raw_only,
                         include_un=include_un)
    for name in sets:
        df = pd.read_csv(DATA_ROOT / name / "data.csv.gz")
        for _, row in df.iterrows():
            sid = row["material_id"]
            try:
                s = Structure.from_dict(json.loads(row["structure"]))
            except Exception as exc:  # noqa: BLE001
                s = None
                if isinstance(row.get("cif"), str):
                    try:
                        s = Structure.from_str(row["cif"], fmt="cif")
                    except Exception:  # noqa: BLE001
                        pass
                if s is None:
                    if errors is not None:
                        errors.append((name, sid, repr(exc)[:200]))
                    continue
            yield (name, sid, s, row) if with_row else (name, sid, s)


if __name__ == "__main__":
    import argparse
    from collections import Counter, defaultdict

    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-only", action="store_true")
    ap.add_argument("--include-splits", action="store_true")
    ap.add_argument("--include-un", action="store_true")
    a = ap.parse_args()
    errs: list = []
    n = Counter()
    formulas, chemsys, elements = set(), set(), set()
    per_set_formulas = defaultdict(set)
    for name, sid, s in iter_structures(raw_only=a.raw_only, include_splits=a.include_splits,
                                         include_un=a.include_un, errors=errs):
        n[name] += 1
        comp = s.composition
        rf = comp.reduced_formula
        formulas.add(rf)
        per_set_formulas[name].add(rf)
        chemsys.add(comp.chemical_system)
        elements.update(e.symbol for e in comp.elements)
    for k in n:
        print(f"{k:60s} {n[k]:7d} structs {len(per_set_formulas[k]):6d} formulas")
    print("total structures", sum(n.values()), "parse failures", len(errs))
    for e in errs[:20]:
        print("  FAIL", e)
    print("distinct elements", len(elements))
    print("distinct reduced formulas", len(formulas))
    print("distinct chemical systems", len(chemsys))
    sys.exit(0)
