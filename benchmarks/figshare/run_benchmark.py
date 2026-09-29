"""Check the generated crystals of figshare 29145101 against the COD.

Dataset: "Generated crystals (+DFT for some) for WyFormer, DiffCSP, DiffCSP++,
WyCryst, SymmCD, CrystalFormer, MiAD, MatterGen, Crys-JEPA, and Composition",
https://figshare.com/articles/dataset/29145101.  Download it with
``download_figshare.py`` and point ``FIGSHARE_DATA`` at its ``data`` folder.

Usage::

    python run_benchmark.py --db ~/.cache/opendisorder/cod.sqlite --out results/ --processes 4

Per-set results are written to ``<out>/<set>.csv.gz`` (sets already present
are skipped), and a summary to ``<out>/summary.csv`` / ``summary.md``.
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
import warnings
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from load_generated import RAW_SETS, iter_structures  # noqa: E402

from opendisorder import CheckOptions, CODClient, CODIndex, MatchKind, NoveltyChecker  # noqa: E402

CONTROL_SETS = {"mp_20/split/test": "MP-20 test split (control)"}
KINDS = [k.value for k in MatchKind]


def set_label(name: str) -> str:
    return RAW_SETS.get(name) or CONTROL_SETS.get(name) or name


def result_path(out: Path, name: str) -> Path:
    return out / (name.replace("/", "__") + ".csv.gz")


def run_set(checker: NoveltyChecker, name: str, processes: int, fetch: bool) -> pd.DataFrame:
    ids, structures = [], []
    for _, sid, s in iter_structures(sets=[name]):
        ids.append(sid)
        structures.append(s)
    t0 = time.time()
    reports = checker.check_many(structures, processes=processes, fetch=fetch, progress=True)
    logging.info("%s: %d structures in %.0f s", name, len(structures), time.time() - t0)
    rows = []
    for sid, s, rep in zip(ids, structures, reports):
        best = rep.best()
        rows.append(
            {
                "set": name,
                "model": set_label(name),
                "sample_id": sid,
                "formula": rep.formula,
                "chemsys": rep.chemsys,
                "nsites": len(s),
                "kind": rep.kind.value,
                "is_novel": rep.is_novel,
                "is_novel_strict": rep.is_novel_strict,
                "best_cod_id": best.reference_id if best else None,
                "cod_ids": ";".join(f"{m.reference_id}:{m.kind.value}" for m in rep.matches),
                "max_site_occupancy_deviation": (best.details.get("max_site_occupancy_deviation") if best else None),
                "n_pristine_candidates": rep.n_pristine_candidates,
                "n_disorder_candidates": rep.n_disorder_candidates,
                "errors": "; ".join(rep.errors),
            }
        )
    return pd.DataFrame(rows)


def summarize(frames: list[pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for df in frames:
        n = len(df)
        counts = df["kind"].value_counts()
        row = {"set": df["set"].iloc[0], "model": df["model"].iloc[0], "n": n}
        for k in KINDS:
            row[f"{k} %"] = 100.0 * counts.get(k, 0) / n
        row["known (pristine + realization) %"] = 100.0 * (~df["is_novel"]).sum() / n
        row["novel %"] = 100.0 * df["is_novel"].sum() / n
        row["novel (strict) %"] = 100.0 * df["is_novel_strict"].sum() / n
        row["with errors/timeouts"] = int((df["errors"].fillna("") != "").sum())
        rows.append(row)
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("results"))
    ap.add_argument("--sets", nargs="*", help="set names (default: all raw generated sets + MP-20 test control)")
    ap.add_argument("--processes", type=int, default=4)
    ap.add_argument("--offline", action="store_true", help="do not download missing CIFs")
    ap.add_argument("--time-limit", type=float, default=30.0, help="seconds per structure comparison")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    warnings.filterwarnings("ignore")
    args.out.mkdir(parents=True, exist_ok=True)
    sets = args.sets or [*RAW_SETS, *CONTROL_SETS]

    with CODIndex(args.db, client=None if args.offline else CODClient()) as index:
        checker = NoveltyChecker(index, options=CheckOptions(time_limit=args.time_limit))
        frames = []
        for name in sets:
            path = result_path(args.out, name)
            if path.exists() and not args.overwrite:
                logging.info("%s: reusing %s", name, path)
            else:
                run_set(checker, name, args.processes, fetch=not args.offline).to_csv(path, index=False)
            frames.append(pd.read_csv(path, dtype={"best_cod_id": str, "cod_ids": str, "errors": str}))

    summary = summarize(frames)
    summary.to_csv(args.out / "summary.csv", index=False)
    (args.out / "summary.md").write_text(summary.to_markdown(index=False, floatfmt=".1f"))
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.1f}"))


if __name__ == "__main__":
    main()
