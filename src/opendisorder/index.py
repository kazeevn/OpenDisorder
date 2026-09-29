"""Local SQLite index of the Crystallography Open Database.

The index is built from COD's CSV metadata dump (one row per entry) and keys
every entry by

* reduced formula -- candidates for the pristine (layer 1) check, and
* chemical system (sorted element set) -- candidates for the disorder
  (layer 2) check, since an ordered realization of a disordered phase
  generally has a different stoichiometry from the one written in the CIF.

Both COD's ``formula`` and ``calcformula`` fields are indexed, because they
disagree for ~15% of entries.  CIFs are fetched lazily (zip batches of 300),
parsed with pymatgen, and the parsed structures are cached in the database, so
each CIF is downloaded and parsed at most once.
"""

from __future__ import annotations

import csv
import io
import json
import logging
import sqlite3
import warnings
import zlib
from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from pymatgen.core import Composition, Structure

from .cod import CODClient

logger = logging.getLogger(__name__)

SCHEMA = """
CREATE TABLE IF NOT EXISTS entries (
    cod_id TEXT PRIMARY KEY,
    formula TEXT,
    calcformula TEXT,
    sg_number INTEGER,
    volume REAL,
    flags TEXT,
    method TEXT,
    svnrevision TEXT,
    fractional INTEGER,
    has_coordinates INTEGER,
    has_disorder INTEGER,
    natoms_cell REAL    -- atoms in the conventional cell (COD cellformula)
);
CREATE TABLE IF NOT EXISTS keys (
    cod_id TEXT,
    kind TEXT,          -- 'rf' (reduced formula) or 'cs' (chemical system)
    value TEXT,
    PRIMARY KEY (kind, value, cod_id)
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS cifs (
    cod_id TEXT PRIMARY KEY,
    data BLOB
);
CREATE TABLE IF NOT EXISTS structures (
    cod_id TEXT PRIMARY KEY,
    ok INTEGER,
    is_ordered INTEGER,
    nsites INTEGER,
    reduced_formula TEXT,
    chemsys TEXT,
    data BLOB,
    warnings TEXT,
    error TEXT
);
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
"""


def cod_formula_to_composition(formula: str | None) -> Composition | None:
    """Parse a COD formula string such as ``"- Fe Li O4 P -"``."""
    if not formula:
        return None
    text = formula.strip().strip("-").strip()
    if not text:
        return None
    try:
        comp = Composition(text)
    except Exception:
        return None
    if comp.num_atoms <= 0:
        return None
    return comp


def composition_keys(comp: Composition) -> tuple[str, str]:
    return comp.reduced_formula, "-".join(sorted({el.symbol for el in comp.elements}))


@dataclass
class ParsedEntry:
    cod_id: str
    structure: Structure | None
    warnings: list[str]
    error: str | None = None


def parse_cif(cod_id: str, cif_text: str, occupancy_tolerance: float = 1.05) -> ParsedEntry:
    """Parse the first structure of a CIF, recording pymatgen warnings."""
    from pymatgen.io.cif import CifParser

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            parser = CifParser(io.StringIO(cif_text), occupancy_tolerance=occupancy_tolerance)
            structures = parser.parse_structures(primitive=False, on_error="ignore")
        except Exception as exc:
            return ParsedEntry(cod_id, None, [str(w.message) for w in caught], f"{type(exc).__name__}: {exc}")
    msgs = sorted({str(w.message)[:300] for w in caught})
    if not structures:
        return ParsedEntry(cod_id, None, msgs, "no structure parsed")
    s = structures[0]
    try:
        s.remove_oxidation_states()
        for prop in list(s.site_properties):
            s.remove_site_property(prop)
    except Exception as exc:
        return ParsedEntry(cod_id, None, msgs, f"cleanup failed: {exc}")
    return ParsedEntry(cod_id, s, msgs)


def _parse_job(args) -> ParsedEntry:
    return parse_cif(*args)


def _pack(obj) -> bytes:
    return zlib.compress(json.dumps(obj).encode())


def _unpack(blob: bytes):
    return json.loads(zlib.decompress(blob).decode())


class CODIndex:
    """SQLite-backed COD index with lazy CIF download and structure caching.

    Args:
        path: SQLite database file.
        client: :class:`CODClient` used to download CIFs (``None`` = offline).
        cif_dir: optional local mirror of COD CIFs (the ``cif/`` tree of the
            COD archive, ``cif/1/00/00/1000041.cif``, or a flat directory);
            used in preference to downloading.
        max_sites: CIFs whose parsed structure has more sites are stored as
            errors (they are molecular crystals, irrelevant for inorganic
            matching and slow to compare).
    """

    def __init__(
        self,
        path: str | Path,
        client: CODClient | None = None,
        cif_dir: str | Path | None = None,
        max_sites: int = 500,
        occupancy_tolerance: float = 1.05,
    ):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.client = client
        self.cif_dir = Path(cif_dir) if cif_dir else None
        self.max_sites = max_sites
        self.occupancy_tolerance = occupancy_tolerance
        self.db = sqlite3.connect(self.path)
        self.db.executescript(SCHEMA)
        self.db.execute("PRAGMA journal_mode=WAL")

    def close(self):
        self.db.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # ------------------------------------------------------------------ build

    def build_from_metadata(self, csv_path: str | Path, include_no_coordinates: bool = False) -> int:
        """(Re)populate ``entries`` and ``keys`` from a COD CSV metadata dump."""
        n = 0
        with open(csv_path, newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(line for line in fh if not line.startswith("#"))
            # metadata tables are rebuilt from scratch (schema may have changed);
            # the CIF and parsed-structure caches are kept.
            self.db.execute("DROP TABLE IF EXISTS entries")
            self.db.execute("DROP TABLE IF EXISTS keys")
            self.db.executescript(SCHEMA)
            with self.db:
                for row in reader:
                    n += self._add_metadata_row(row, include_no_coordinates)
                self.db.execute(
                    "INSERT OR REPLACE INTO meta VALUES ('metadata_source', ?)", (str(Path(csv_path).resolve()),)
                )
        logger.info("indexed %d COD entries", n)
        return n

    def add_entries(self, records: Iterable[dict], include_no_coordinates: bool = False) -> int:
        """Add metadata records (as returned by :meth:`CODClient.search`)."""
        with self.db:
            return sum(self._add_metadata_row(r, include_no_coordinates) for r in records)

    def _add_metadata_row(self, row: dict, include_no_coordinates: bool) -> int:
        flags = row.get("flags") or ""
        has_coords = "has coordinates" in flags
        if not has_coords and not include_no_coordinates:
            return 0
        cod_id = str(row["file"])
        keys = set()
        fractional = False
        for field in ("formula", "calcformula"):
            comp = cod_formula_to_composition(row.get(field))
            if comp is None:
                continue
            fractional |= any(abs(amt - round(amt)) > 1e-6 for amt in comp.values())
            rf, cs = composition_keys(comp)
            keys.add(("rf", rf))
            keys.add(("cs", cs))
        if not keys:
            return 0
        cell_comp = cod_formula_to_composition(row.get("cellformula"))

        def num(x, typ):
            try:
                return typ(x)
            except (TypeError, ValueError):
                return None

        self.db.execute(
            "INSERT OR REPLACE INTO entries VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (
                cod_id,
                row.get("formula"),
                row.get("calcformula"),
                num(row.get("sgNumber"), int),
                num(row.get("vol"), float),
                flags,
                row.get("method"),
                row.get("svnrevision"),
                int(fractional),
                int(has_coords),
                int("has disorder" in flags),
                cell_comp.num_atoms if cell_comp is not None else None,
            ),
        )
        self.db.executemany("INSERT OR IGNORE INTO keys VALUES (?,?,?)", [(cod_id, k, v) for k, v in keys])
        return 1

    # ------------------------------------------------------------------ lookup

    def n_entries(self) -> int:
        return self.db.execute("SELECT COUNT(*) FROM entries").fetchone()[0]

    def _ids(self, kind: str, value: str, disorder_only: bool, max_cell_atoms: float | None) -> list[str]:
        sql = "SELECT k.cod_id FROM keys k JOIN entries e USING (cod_id) WHERE k.kind=? AND k.value=?"
        params: list = [kind, value]
        if disorder_only:
            sql += " AND (e.has_disorder=1 OR e.fractional=1)"
        if max_cell_atoms is not None:
            sql += " AND (e.natoms_cell IS NULL OR e.natoms_cell <= ?)"
            params.append(max_cell_atoms)
        return [r[0] for r in self.db.execute(sql + " ORDER BY k.cod_id", params)]

    def ids_by_reduced_formula(self, reduced_formula: str, max_cell_atoms: float | None = None) -> list[str]:
        """COD ids whose ``formula`` or ``calcformula`` reduces to ``reduced_formula``.

        ``max_cell_atoms`` drops entries with more atoms in the conventional cell.
        """
        return self._ids("rf", reduced_formula, False, max_cell_atoms)

    def ids_by_chemsys(
        self, chemsys: str, disorder_only: bool = False, max_cell_atoms: float | None = None
    ) -> list[str]:
        """COD ids with exactly the elements of ``chemsys`` (e.g. ``"Fe-Li-O"``).

        ``disorder_only`` keeps entries flagged ``has disorder`` by COD or with
        a fractional formula.
        """
        return self._ids("cs", chemsys, disorder_only, max_cell_atoms)

    def entry(self, cod_id: str) -> dict | None:
        cur = self.db.execute("SELECT * FROM entries WHERE cod_id=?", (str(cod_id),))
        row = cur.fetchone()
        return dict(zip([d[0] for d in cur.description], row)) if row else None

    # ------------------------------------------------------------------ CIFs & structures

    def _local_cif(self, cod_id: str) -> str | None:
        if self.cif_dir is None:
            return None
        cid = str(cod_id)
        for p in (
            self.cif_dir / cid[0] / cid[1:3] / cid[3:5] / f"{cid}.cif",
            self.cif_dir / "cif" / cid[0] / cid[1:3] / cid[3:5] / f"{cid}.cif",
            self.cif_dir / f"{cid}.cif",
        ):
            if p.exists():
                return p.read_text(encoding="utf-8", errors="replace")
        return None

    def add_cif(self, cod_id: str, cif_text: str):
        self.db.execute("INSERT OR REPLACE INTO cifs VALUES (?,?)", (str(cod_id), zlib.compress(cif_text.encode())))

    def get_cif(self, cod_id: str) -> str | None:
        row = self.db.execute("SELECT data FROM cifs WHERE cod_id=?", (str(cod_id),)).fetchone()
        if row:
            return zlib.decompress(row[0]).decode()
        return self._local_cif(cod_id)

    def fetch_cifs(self, cod_ids: Iterable[str]) -> int:
        """Make sure the CIFs of ``cod_ids`` are available locally; returns the number downloaded."""
        missing = []
        for cid in dict.fromkeys(str(i) for i in cod_ids):
            if self.db.execute("SELECT 1 FROM structures WHERE cod_id=?", (cid,)).fetchone():
                continue
            if self.db.execute("SELECT 1 FROM cifs WHERE cod_id=?", (cid,)).fetchone() or self._local_cif(cid):
                continue
            missing.append(cid)
        if not missing:
            return 0
        if self.client is None:
            logger.warning("%d CIFs missing and no COD client configured", len(missing))
            return 0
        n = 0
        for i, (cid, text) in enumerate(self.client.iter_cifs(missing)):
            self.add_cif(cid, text)
            n += 1
            if i % 300 == 299:
                self.db.commit()
                logger.info("downloaded %d/%d CIFs", n, len(missing))
        self.db.commit()
        return n

    def parse_missing(self, cod_ids: Iterable[str], processes: int = 1) -> int:
        """Parse (and cache) structures for ``cod_ids`` whose CIFs are available."""
        jobs = []
        for cid in dict.fromkeys(str(i) for i in cod_ids):
            if self.db.execute("SELECT 1 FROM structures WHERE cod_id=?", (cid,)).fetchone():
                continue
            text = self.get_cif(cid)
            if text is not None:
                jobs.append((cid, text, self.occupancy_tolerance))
        if not jobs:
            return 0
        if processes > 1 and len(jobs) > 20:
            with ProcessPoolExecutor(processes) as pool:
                results: Iterator[ParsedEntry] = pool.map(_parse_job, jobs, chunksize=16)
                self._store_parsed(results)
        else:
            self._store_parsed(map(_parse_job, jobs))
        return len(jobs)

    def _store_parsed(self, results: Iterable[ParsedEntry]):
        for i, res in enumerate(results):
            s = res.structure
            if s is not None and len(s) > self.max_sites:
                res = ParsedEntry(res.cod_id, None, res.warnings, f"too many sites ({len(s)})")
                s = None
            if s is None:
                row = (res.cod_id, 0, None, None, None, None, None, json.dumps(res.warnings), res.error)
            else:
                comp = s.composition
                row = (
                    res.cod_id,
                    1,
                    int(s.is_ordered),
                    len(s),
                    comp.reduced_formula,
                    "-".join(sorted({el.symbol for el in comp.elements})),
                    _pack(s.as_dict()),
                    json.dumps(res.warnings),
                    None,
                )
            self.db.execute("INSERT OR REPLACE INTO structures VALUES (?,?,?,?,?,?,?,?,?)", row)
            if i % 500 == 499:
                self.db.commit()
        self.db.commit()

    def prepare(self, cod_ids: Iterable[str], processes: int = 1):
        """Download and parse everything needed for ``cod_ids``."""
        ids = list(cod_ids)
        self.fetch_cifs(ids)
        self.parse_missing(ids, processes=processes)

    def structure_info(self, cod_id: str) -> dict | None:
        cur = self.db.execute(
            "SELECT cod_id, ok, is_ordered, nsites, reduced_formula, chemsys, warnings, error FROM structures WHERE cod_id=?",
            (str(cod_id),),
        )
        row = cur.fetchone()
        return dict(zip([d[0] for d in cur.description], row)) if row else None

    def get_structure(self, cod_id: str, fetch: bool = True) -> Structure | None:
        cid = str(cod_id)
        row = self.db.execute("SELECT ok, data FROM structures WHERE cod_id=?", (cid,)).fetchone()
        if row is None and fetch:
            self.prepare([cid])
            row = self.db.execute("SELECT ok, data FROM structures WHERE cod_id=?", (cid,)).fetchone()
        if row is None or not row[0]:
            return None
        return Structure.from_dict(_unpack(row[1]))

    def get_structures(self, cod_ids: Iterable[str], fetch: bool = True, processes: int = 1) -> dict[str, Structure]:
        ids = [str(i) for i in cod_ids]
        if fetch:
            self.prepare(ids, processes=processes)
        out = {}
        for cid in ids:
            s = self.get_structure(cid, fetch=False)
            if s is not None:
                out[cid] = s
        return out
