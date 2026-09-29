"""Command line interface: ``opendisorder build-index`` and ``opendisorder check``."""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import warnings
from pathlib import Path

from pymatgen.core import Structure

from .checker import CheckOptions, NoveltyChecker
from .cod import CODClient
from .index import CODIndex

DEFAULT_DB = Path(os.environ.get("OPENDISORDER_DB", Path.home() / ".cache" / "opendisorder" / "cod.sqlite"))


def _build_index(args) -> int:
    client = CODClient()
    csv_path = args.metadata
    if csv_path is None:
        csv_path = args.db.with_name("cod_metadata.csv")
        if not csv_path.exists() or args.refresh:
            logging.info("downloading COD metadata to %s", csv_path)
            client.download_metadata(csv_path)
    with CODIndex(args.db, client=client) as index:
        n = index.build_from_metadata(csv_path)
    print(f"indexed {n} COD entries into {args.db}")
    return 0


def _check(args) -> int:
    options = CheckOptions(
        check_disorder=not args.no_disorder,
        check_siblings=not args.no_siblings,
        exhaustive=args.exhaustive,
        disorder_candidates_flagged_only=args.flagged_only,
    )
    structures = [Structure.from_file(f) for f in args.files]
    with CODIndex(args.db, client=None if args.offline else CODClient(), cif_dir=args.cif_dir) as index:
        if index.n_entries() == 0:
            print(f"index {args.db} is empty; run `opendisorder build-index` first", file=sys.stderr)
            return 1
        checker = NoveltyChecker(index, options=options)
        reports = checker.check_many(structures, processes=args.processes, fetch=not args.offline)
    for f, rep in zip(args.files, reports):
        if args.json:
            print(json.dumps({"file": str(f), **rep.to_dict()}))
        else:
            ids = ", ".join(f"{m.reference_id} ({m.kind.value})" for m in rep.matches) or "-"
            print(f"{f}\t{rep.formula}\t{'novel' if rep.is_novel else 'known'}\t{rep.kind.value}\t{ids}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="opendisorder", description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB, help=f"index database (default {DEFAULT_DB})")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build-index", help="build the COD index from the metadata dump")
    b.add_argument("--metadata", type=Path, help="existing COD CSV dump (result?id=%%&format=csv)")
    b.add_argument("--refresh", action="store_true", help="re-download the metadata dump")
    b.set_defaults(func=_build_index)

    c = sub.add_parser("check", help="check structure files (CIF, POSCAR, ...) against the COD")
    c.add_argument("files", nargs="+", type=Path)
    c.add_argument("--json", action="store_true", help="one JSON report per line")
    c.add_argument("--processes", type=int, default=1)
    c.add_argument("--offline", action="store_true", help="only use CIFs already in the index / cif-dir")
    c.add_argument("--cif-dir", type=Path, help="local mirror of COD CIFs")
    c.add_argument("--no-disorder", action="store_true", help="skip layer 2")
    c.add_argument("--no-siblings", action="store_true", help="skip the sibling-ordering check")
    c.add_argument("--flagged-only", action="store_true", help="layer 2: only entries flagged as disordered")
    c.add_argument("--exhaustive", action="store_true", help="report all matching entries")
    c.set_defaults(func=_check)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING, format="%(levelname)s %(message)s")
    if not args.verbose:
        warnings.filterwarnings("ignore")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
