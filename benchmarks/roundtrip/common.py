"""Shared helpers for the round-trip benchmark (COD access, I/O, disorder analysis).

Nothing in here writes to the COD index: it is always opened read-only
(``sqlite3`` URI ``mode=ro``), so the benchmark can run next to other jobs
that use the same database.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import re
import signal
import sqlite3
import zlib
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from fractions import Fraction
from pathlib import Path

import numpy as np

DEFAULT_COD_DB = os.environ.get("OD_COD_DB", str(Path("~/.cache/opendisorder/cod.sqlite").expanduser()))

#: pymatgen warnings that do not indicate a problem with the CIF itself
BENIGN_WARNINGS = ("Initializing CifParser from StringIO is deprecated",)

ANIONS = {
    "O": "oxide",
    "S": "chalcogenide",
    "Se": "chalcogenide",
    "Te": "chalcogenide",
    "F": "halide",
    "Cl": "halide",
    "Br": "halide",
    "I": "halide",
    "N": "nitride",
}


# ----------------------------------------------------------------------------- COD access


def cod_connect(path: str | os.PathLike | None = None) -> sqlite3.Connection:
    """Read-only connection to an opendisorder COD index."""
    path = Path(path or DEFAULT_COD_DB).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"COD index not found: {path} (set --cod-db or $OD_COD_DB)")
    return sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=120)


def cod_cif(conn: sqlite3.Connection, cod_id: str) -> str | None:
    row = conn.execute("SELECT data FROM cifs WHERE cod_id=?", (str(cod_id),)).fetchone()
    if row is None or row[0] is None:
        return None
    return zlib.decompress(row[0]).decode("utf-8", "replace")


def cod_structure(conn: sqlite3.Connection, cod_id: str, use_cache: bool = True):
    """Parsed COD structure: the index's cached parse if present, else parse the stored CIF
    with :func:`opendisorder.index.parse_cif` (exactly what ``NoveltyChecker`` would see)."""
    from pymatgen.core import Structure

    if use_cache:
        row = conn.execute("SELECT ok, data FROM structures WHERE cod_id=?", (str(cod_id),)).fetchone()
        if row is not None and row[0] and row[1]:
            return Structure.from_dict(json.loads(zlib.decompress(row[1]).decode()))
    text = cod_cif(conn, cod_id)
    if text is None:
        return None
    from opendisorder.index import parse_cif

    return parse_cif(str(cod_id), text).structure


def same_chemsys_disordered_ids(
    conn: sqlite3.Connection, chemsys: str, max_cell_atoms: float | None = None
) -> list[str]:
    """COD ids of the chemical system that COD flags ``has disorder`` or whose formula is
    fractional (the default layer-2 candidate set of ``NoveltyChecker``)."""
    sql = (
        "SELECT DISTINCT k.cod_id FROM keys k JOIN entries e USING (cod_id) "
        "WHERE k.kind='cs' AND k.value=? AND (e.has_disorder=1 OR e.fractional=1)"
    )
    params: list = [chemsys]
    if max_cell_atoms is not None:
        sql += " AND (e.natoms_cell IS NULL OR e.natoms_cell <= ?)"
        params.append(max_cell_atoms)
    return [r[0] for r in conn.execute(sql + " ORDER BY k.cod_id", params)]


# ----------------------------------------------------------------------------- I/O


def open_text(path: str | os.PathLike, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8") if "b" not in mode else gzip.open(path, mode)
    return open(path, mode.replace("t", ""), encoding="utf-8") if "b" not in mode else open(path, mode)


def iter_jsonl(paths: str | os.PathLike | Iterable) -> Iterator[dict]:
    if isinstance(paths, (str, os.PathLike)):
        paths = [paths]
    for p in paths:
        if not Path(p).exists():
            continue
        with open_text(p, "rt") as fh:
            try:
                for line in fh:
                    line = line.strip()
                    if line:
                        try:
                            yield json.loads(line)
                        except json.JSONDecodeError:  # truncated last line of an interrupted run
                            continue
            except (EOFError, OSError, zlib.error):
                # truncated gzip member of a killed run; plain .jsonl is recommended for
                # outputs that are appended to (resumable runs)
                continue


def done_ids(path: str | os.PathLike, key: str = "id") -> set:
    return {rec.get(key) for rec in iter_jsonl(path)} if Path(path).exists() else set()


class JsonlWriter:
    """Append-mode JSON-lines writer (flushes every record, so runs are resumable)."""

    def __init__(self, path: str | os.PathLike):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        self.path = str(path)
        self.fh = gzip.open(self.path, "at", encoding="utf-8") if self.path.endswith(".gz") else open(
            self.path, "a", encoding="utf-8"
        )

    def write(self, rec: dict):
        self.fh.write(json.dumps(rec, default=_json_default) + "\n")
        self.fh.flush()

    def close(self):
        self.fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if hasattr(o, "value"):
        return o.value
    return str(o)


def stable_shard(key: str, nshards: int) -> int:
    return zlib.crc32(key.encode()) % nshards


def parse_shard(text: str | None) -> tuple[int, int]:
    """``"3/16"`` -> (3, 16); ``None`` -> (0, 1).  Also reads SLURM_ARRAY_TASK_ID if text is 'slurm/N'."""
    if not text:
        return 0, 1
    i, n = text.split("/")
    if i == "slurm":
        i = os.environ["SLURM_ARRAY_TASK_ID"]
    i, n = int(i), int(n)
    if not 0 <= i < n:
        raise ValueError(f"bad shard {text}")
    return i, n


# ----------------------------------------------------------------------------- timeouts


class Timeout(Exception):
    pass


@contextmanager
def time_limit(seconds: float | None):
    if not seconds or not hasattr(signal, "SIGALRM"):
        yield
        return

    def handler(signum, frame):
        raise Timeout()

    old = signal.signal(signal.SIGALRM, handler)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def _limited_child(conn, fn, args, mem_bytes):
    import resource
    import traceback

    try:
        os.setsid()  # own process group, so enum.x & co. are killed with us on timeout
        if mem_bytes:
            resource.setrlimit(resource.RLIMIT_AS, (mem_bytes, mem_bytes))
        conn.send(("ok", fn(*args)))
    except BaseException as exc:  # noqa: BLE001 - report everything to the parent
        conn.send(("error", f"{type(exc).__name__}: {exc}\n{traceback.format_exc()[-500:]}"))
    finally:
        conn.close()


def run_limited(fn, args: tuple, timeout: float | None, mem_gb: float | None):
    """Run ``fn(*args)`` in a forked child with an address-space limit and a wall-clock limit.

    Used for enumlib, whose output (all symmetry-distinct structures of a cell size) can
    explode combinatorially; pymatgen reads all of it into memory."""
    import multiprocessing as mp

    ctx = mp.get_context("fork")
    parent, child = ctx.Pipe(duplex=False)
    p = ctx.Process(target=_limited_child, args=(child, fn, args, int(mem_gb * 2**30) if mem_gb else 0))
    p.start()
    child.close()
    if not parent.poll(timeout):
        try:
            os.killpg(p.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        p.join()
        raise Timeout()
    try:
        status, value = parent.recv()
    except EOFError:
        p.join()
        raise RuntimeError(f"child died (exit code {p.exitcode}), probably out of memory") from None
    p.join()
    if status == "error":
        raise RuntimeError(value)
    return value


# ----------------------------------------------------------------------------- formulas

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")


def cod_formula_elements(formula: str | None) -> dict[str, float]:
    """``"- Fe Li O4 P -"`` -> {'Fe': 1, 'Li': 1, 'O': 4, 'P': 1} (tolerant, no pymatgen)."""
    out: dict[str, float] = {}
    if not formula:
        return out
    for el, amt in _TOKEN.findall(formula.strip().strip("-")):
        try:
            out[el] = out.get(el, 0.0) + (float(amt) if amt else 1.0)
        except ValueError:
            continue
    return out


def structure_class(elements: Iterable[str]) -> str:
    """Coarse chemical class used to bucket the results."""
    els = set(elements)
    kinds = {ANIONS[e] for e in els if e in ANIONS}
    if not kinds:
        return "intermetallic"
    if len(kinds) > 1:
        # e.g. oxyhalides, oxynitrides, oxysulfides
        return "mixed-anion"
    return kinds.pop()


# ----------------------------------------------------------------------------- disorder analysis


def rational(x: float, max_den: int, tol: float) -> Fraction | None:
    """Closest fraction with denominator <= max_den, or None if further than tol."""
    f = Fraction(x).limit_denominator(max_den)
    return f if abs(float(f) - x) <= tol else None


def orbit_labels(structure, symprec: float = 0.01) -> np.ndarray:
    """Symmetry-orbit index of every site (disordered sites are typed by species+occupancy)."""
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    try:
        ds = SpacegroupAnalyzer(structure, symprec=symprec).get_symmetry_dataset()
        eq = np.asarray(ds.equivalent_atoms)
    except Exception:
        eq = np.arange(len(structure))
    _, labels = np.unique(eq, return_inverse=True)
    return labels


def analyze_disorder(
    structure,
    max_den: int = 12,
    occ_tol: float = 0.011,
    ordered_tol: float = 0.02,
    symprec: float = 0.01,
    split_factor: float = 0.7,
    assume_primitive: bool = False,
) -> dict:
    """Describe the compositional disorder of a (conventional or primitive) structure.

    Works on the primitive cell. For every symmetry orbit it records the species, the
    occupancies and their best rational approximation (denominator <= ``max_den``),
    and derives the smallest supercell (in primitive cells) in which every orbit can
    be realized with exactly the rational occupancies.

    ``split_factor``: two sites, at least one of them partially occupied, closer than
    ``split_factor * (r_i + r_j)`` (atomic radii) are taken as positional (split-site)
    disorder or mutually exclusive interstitial sites, which random colorings cannot
    represent physically.
    """
    from pymatgen.core import Element

    prim = structure if assume_primitive else structure.get_primitive_structure()
    labels = orbit_labels(prim, symprec)
    orbits = []
    min_super = 1
    max_err = 0.0
    irrational = False
    for lab in np.unique(labels):
        idx = np.where(labels == lab)[0]
        site = prim[int(idx[0])]
        occ = {el.symbol: float(x) for el, x in site.species.items()}
        total = sum(occ.values())
        vac = max(0.0, 1.0 - total)
        if total > 1.0 + ordered_tol:
            kind = "overfull"
        elif len(occ) == 1 and vac <= ordered_tol:
            kind = "ordered"
        elif len(occ) >= 2 and vac <= ordered_tol:
            kind = "substitutional"
        elif len(occ) == 1:
            kind = "vacancy"
        else:
            kind = "substitutional+vacancy"
        m = len(idx)
        fracs = {}
        if kind != "ordered":
            parts = dict(occ)
            if vac > ordered_tol:
                parts["Vac"] = vac
            for sp, x in parts.items():
                f = rational(x, max_den, occ_tol)
                if f is None:
                    irrational = True
                    f = Fraction(x).limit_denominator(max_den)
                max_err = max(max_err, abs(float(f) - x))
                fracs[sp] = f
                # k * m * f must be an integer
                min_super = math.lcm(min_super, (Fraction(m) * f).denominator)
        orbits.append(
            {
                "orbit": int(lab),
                "multiplicity": m,
                "species": occ,
                "vacancy": round(vac, 4),
                "kind": kind,
                "rational": {k: str(v) for k, v in fracs.items()},
            }
        )
    kinds = {o["kind"] for o in orbits if o["kind"] != "ordered"}
    has_sub = any("substitutional" in k for k in kinds)
    has_vac = any("vacancy" in k for k in kinds)
    dtype = "both" if has_sub and has_vac else "substitutional" if has_sub else "vacancy" if has_vac else "ordered"
    if "overfull" in kinds:
        dtype = "overfull"

    # close partially occupied sites (positional disorder / exclusion)
    min_ratio, min_d = np.inf, np.inf
    partial = [i for i, s in enumerate(prim) if not s.is_ordered]
    if partial:
        dm = prim.distance_matrix
        radius = {}
        for s in prim:
            for el in s.species.elements:
                try:
                    radius[el.symbol] = float(Element(el.symbol).atomic_radius or 1.5)
                except Exception:
                    radius[el.symbol] = 1.5
        rmin = [min(radius[el.symbol] for el in s.species.elements) for s in prim]
        for i in partial:
            for j in range(len(prim)):
                if i == j:
                    continue
                d = dm[i, j]
                ratio = d / (rmin[i] + rmin[j])
                if ratio < min_ratio:
                    min_ratio, min_d = ratio, d
    return {
        "nsites_prim": len(prim),
        "n_orbits": len(orbits),
        "n_disordered_orbits": sum(o["kind"] != "ordered" for o in orbits),
        "disorder_type": dtype,
        "min_supercell": int(min_super),
        "rational_ok": not irrational,
        "max_rational_error": round(max_err, 4),
        "min_partial_distance": None if not np.isfinite(min_d) else round(float(min_d), 3),
        "min_partial_distance_ratio": None if not np.isfinite(min_ratio) else round(float(min_ratio), 3),
        "split_site": bool(np.isfinite(min_ratio) and min_ratio < split_factor),
        "min_minority_occupancy": round(
            min(
                [x for o in orbits if o["kind"] != "ordered" for x in [*o["species"].values(), o["vacancy"]] if x > 0.0]
                or [1.0]
            ),
            4,
        ),
        "orbits": orbits,
    }


CIF_POSITIONAL_TAGS = ("_atom_site_disorder_assembly", "_atom_site_disorder_group")
CIF_POWDER_TAGS = ("_pd_proc_ls_prof_r_factor", "_pd_proc_ls_prof_wr_factor", "_pd_meas_", "_pd_proc_", "_pd_spec_")


def cif_flags(cif_text: str) -> dict:
    """Cheap textual checks on a CIF: powder refinement tags, disorder-group tags."""
    low = cif_text.lower()
    powder = any(t in low for t in CIF_POWDER_TAGS)
    disorder_groups = False
    # the tags are loop columns; any value other than '.'/'?' indicates split positions
    lines = cif_text.splitlines()
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("loop_"):
            cols, j = [], i + 1
            while j < len(lines) and lines[j].strip().startswith("_"):
                cols.append(lines[j].strip().split()[0].lower())
                j += 1
            hits = [cols.index(t) for t in CIF_POSITIONAL_TAGS if t in cols]
            if not hits:
                continue
            while j < len(lines) and lines[j].strip() and not lines[j].strip().startswith(("_", "loop_", "#")):
                vals = lines[j].split()
                if len(vals) == len(cols) and any(vals[h] not in (".", "?") for h in hits):
                    disorder_groups = True
                    break
                j += 1
    single_crystal = "_diffrn_reflns_number" in low or "_refine_ls_r_factor_gt" in low
    return {"powder_tags": powder, "disorder_groups": disorder_groups, "single_crystal_tags": single_crystal}


def clean(structure):
    s = structure.copy()
    try:
        s.remove_oxidation_states()
    except Exception:
        pass
    for p in list(s.site_properties):
        s.remove_site_property(p)
    return s


def structure_to_record(structure) -> dict:
    return clean(structure).as_dict(verbosity=0)


def record_structure(rec: dict):
    from pymatgen.core import Structure

    return Structure.from_dict(rec["structure"])
