"""Fetch the LeMat-Bulk rows needed for the literature parent search.

LeMat-Bulk (LeMaterial / Entalpic, https://huggingface.co/datasets/LeMaterial/LeMat-Bulk, CC-BY-4.0)
is a unified, deduplicated set of Materials Project + Alexandria + OQMD relaxed structures.

Strategy (disk-light): the parquet shards (~3.8 GB in total for the four configs) are downloaded
one at a time into the scratch directory, filtered row group by row group on the ``elements``
column, and deleted immediately.  Peak disk use is one shard (<= 190 MB) plus the filtered output.

Kept: every row whose element set (>= 2 elements, <= MAX_SITES sites) is a subset of
  * the chemical system of an unmatched literature structure (``targets()``), or
  * the chemical system of the known phase the paper names, when that phase has elements the
    literature structure does not (structure-type analogues, e.g. K3Nd(PO4)2 for K3Nd(AsO4)2).

Usage::

    python fetch_lemat.py [--out SCRATCH/lemat_subset.parquet]
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import requests

HERE = Path(__file__).resolve().parent
LIT = HERE.parent / "literature"
SCRATCH = Path(os.environ.get(
    "LEMAT_SCRATCH",
    str(Path(__file__).resolve().parent / "_data")))
REPO = "LeMaterial/LeMat-Bulk"
REVISION = "0dc17eea904b860ad7288141e9870f67f8e6bb2c"  # dataset commit used (HF API sha, 2025-04-15)
CONFIGS = ["compatible_pbe", "compatible_pbesol", "compatible_scan", "non_compatible"]
COLUMNS = ["immutable_id", "elements", "chemical_formula_reduced", "nsites", "lattice_vectors",
           "species_at_sites", "cartesian_site_positions", "energy", "functional", "cross_compatibility"]
MAX_SITES = 250

# chemical systems of the paper's named phases that contain elements absent from the query
# (structure-type analogues); used only to look for the analogue's framework in LeMat.
ANALOGUE_SYSTEMS = {
    "K3Nd(PO4)2": {"K", "Nd", "P", "O"},
    "PbPtGd": {"Pb", "Pt", "Gd"},
    "FeTa4Te4": {"Fe", "Ta", "Te"},
    "Ce2Pt1.5Ga6.5": {"Ce", "Pt", "Ga"},
    "CaCu4P2": {"Ca", "Cu", "P"},
    "RbLa2(P2S6)(PS4)2": {"Rb", "La", "P", "S"},
    "Ge3.9Nb20.4Rh5.9": {"Ge", "Nb", "Rh"},
    "LaPO3/LaPO4 (Ac analogue)": {"La", "P", "O"},
    "LaB6 (Ac analogue)": {"La", "Er", "B"},
}


def unmatched():
    import pandas as pd
    r = pd.read_csv(LIT / "results.csv")
    return r[(r["kind"] == "none") & r["paper"].isin(["Leeman2024", "Cheetham&Seshadri2024"])]


def targets() -> list[frozenset]:
    from pymatgen.core import Structure
    out = set()
    for f in unmatched()["file"]:
        s = Structure.from_file(LIT / "cifs" / f)
        out.add(frozenset(el.symbol for el in s.composition.elements))
    out |= {frozenset(v) for v in ANALOGUE_SYSTEMS.values()}
    return sorted(out, key=lambda x: sorted(x))


def shard_list(config):
    url = f"https://huggingface.co/api/datasets/{REPO}/tree/{REVISION}/{config}"
    return [(x["path"], x["size"]) for x in requests.get(url, timeout=60).json() if x["path"].endswith(".parquet")]


def keep_mask(elements_col, tsets):
    return [bool(els) and len(els) >= 2 and any(set(els) <= t for t in tsets) for els in elements_col]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=SCRATCH / "lemat_subset.parquet")
    args = ap.parse_args(argv)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    tsets = targets()
    print(f"{len(tsets)} target chemical systems", flush=True)
    tables, stats = [], []
    for config in CONFIGS:
        for path, size in shard_list(config):
            local = SCRATCH / Path(path).name
            url = f"https://huggingface.co/datasets/{REPO}/resolve/{REVISION}/{path}"
            t0 = time.time()
            subprocess.run(["curl", "-sSfL", "--retry", "5", "-o", str(local), url], check=True)
            pf = pq.ParquetFile(local)
            n_rows, n_keep = pf.metadata.num_rows, 0
            for rg in range(pf.metadata.num_row_groups):
                meta = pf.read_row_group(rg, columns=["elements", "nsites"])
                mask = keep_mask(meta.column("elements").to_pylist(), tsets)
                mask = pa.array([m and n <= MAX_SITES for m, n in zip(mask, meta.column("nsites").to_pylist())])
                if not pa.compute.any(mask).as_py():
                    continue
                t = pf.read_row_group(rg, columns=COLUMNS).filter(mask)
                t = t.append_column("config", pa.array([config] * len(t)))
                tables.append(t)
                n_keep += len(t)
            local.unlink()
            stats.append((path, size, n_rows, n_keep))
            print(f"{path}: {size/1e6:.0f} MB, {n_rows} rows, kept {n_keep} ({time.time()-t0:.0f}s)", flush=True)
    table = pa.concat_tables(tables)
    pq.write_table(table, args.out, compression="zstd")
    print(f"wrote {args.out}: {len(table)} rows, {args.out.stat().st_size/1e6:.1f} MB; "
          f"scanned {sum(s[2] for s in stats)} rows / {sum(s[1] for s in stats)/1e9:.2f} GB", flush=True)


if __name__ == "__main__":
    sys.exit(main())
