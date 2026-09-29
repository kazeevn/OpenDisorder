"""Screen structures that the literature disputes as "new" against the COD.

Cases (see manifest.csv for provenance of every CIF):

* MatterGen TaCr2O6 and related structures -- Juelsholt, Mater. Horiz. 13, 5672 (2026).
* GNoME predictions -- Cheetham & Seshadri, Chem. Mater. 36, 3490 (2024).
* A-Lab targets (Materials Project entries) -- Leeman et al., PRX Energy 3, 011002 (2024).

Usage::

    python run.py --db ~/.cache/opendisorder/cod.sqlite
"""

from __future__ import annotations

import argparse
import logging
import warnings
from pathlib import Path

import pandas as pd
from pymatgen.core import Structure

from opendisorder import CheckOptions, CODClient, CODIndex, NoveltyChecker

HERE = Path(__file__).resolve().parent


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--processes", type=int, default=4)
    ap.add_argument("--out", type=Path, default=HERE / "results.csv")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.WARNING)
    warnings.filterwarnings("ignore")

    manifest = pd.read_csv(HERE / "manifest.csv")
    structures = [Structure.from_file(HERE / "cifs" / f) for f in manifest["file"]]
    with CODIndex(args.db, client=CODClient()) as index:
        checker = NoveltyChecker(index, options=CheckOptions(exhaustive=True))
        reports = checker.check_many(structures, processes=args.processes)
        rows = []
        for (_, m), s, rep in zip(manifest.iterrows(), structures, reports):
            best = rep.best()
            rows.append(
                {
                    **m.to_dict(),
                    "formula": rep.formula,
                    "kind": rep.kind.value,
                    "verdict": "novel" if rep.is_novel_strict else ("known" if not rep.is_novel else "suspicious"),
                    "best_cod_id": best.reference_id if best else None,
                    "best_cod_formula": (index.entry(best.reference_id) or {}).get("formula") if best else None,
                    "cod_matches": ";".join(f"{x.reference_id}:{x.kind.value}" for x in rep.matches),
                    "n_candidates": rep.n_pristine_candidates + rep.n_disorder_candidates + rep.n_substitution_candidates,
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print(df.groupby(["paper", "verdict"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    main()
