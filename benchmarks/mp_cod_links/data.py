"""Data access for the MP <-> COD disorder study.

* COD: read-only access to the opendisorder SQLite index (never writes).
* MP : summary documents (structure, theoretical, e_above_hull) fetched per chemical
  system with the API key in $MP_API_KEY (sent only as the X-API-KEY header) and cached
  as gzipped JSON under ./mp_cache/.
"""
from __future__ import annotations

import gzip, json, os, sqlite3, time, zlib
from functools import lru_cache

from pymatgen.core import Structure

HERE = os.path.dirname(os.path.abspath(__file__))
COD_DB = os.environ.get(
    "OD_COD_DB",
    "/tmp/claude-0/-home-user-OpenDisorder/97bbdad5-5fc7-5d06-ab32-765ce4111df4/scratchpad/db/cod.sqlite",
)
MP_CACHE = os.path.join(HERE, "mp_cache")


# ----------------------------------------------------------------------------- COD
def cod_connect():
    return sqlite3.connect(f"file:{COD_DB}?mode=ro", uri=True, timeout=60)


def cod_candidates(max_atoms=60, min_el=2, max_el=6):
    """Inorganic (no H/D, no C+H) COD entries flagged disordered / fractional."""
    c = cod_connect()
    rows = c.execute(
        "SELECT e.cod_id, k.value, e.natoms_cell, e.sg_number, e.formula FROM entries e "
        "JOIN keys k ON k.cod_id=e.cod_id AND k.kind='cs' "
        "WHERE (e.has_disorder=1 OR e.fractional=1) AND e.natoms_cell <= ?",
        (max_atoms,),
    ).fetchall()
    out = []
    for cid, cs, nat, sg, formula in rows:
        els = cs.split("-")
        if "H" in els or "D" in els or not (min_el <= len(els) <= max_el):
            continue
        out.append(dict(cod_id=cid, chemsys=cs, natoms_cell=nat, sg=sg, formula=formula))
    return out


def cod_structure(cod_id: str, conn=None) -> Structure | None:
    """Parsed COD structure: the index's cached parse if present, else parse the stored CIF."""
    from opendisorder.index import parse_cif

    conn = conn or cod_connect()
    row = conn.execute("SELECT ok, data FROM structures WHERE cod_id=?", (cod_id,)).fetchone()
    if row is not None and row[0] and row[1]:
        return Structure.from_dict(json.loads(zlib.decompress(row[1]).decode()))
    row = conn.execute("SELECT data FROM cifs WHERE cod_id=?", (cod_id,)).fetchone()
    if row is None or row[0] is None:
        return None
    return parse_cif(cod_id, zlib.decompress(row[0]).decode("utf-8", "replace")).structure


# ----------------------------------------------------------------------------- MP
FIELDS = "material_id,formula_pretty,theoretical,energy_above_hull,structure,symmetry"


def _get(path, **params):
    import urllib.parse, urllib.request, urllib.error

    url = "https://api.materialsproject.org" + path + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(
        url,
        headers={"X-API-KEY": os.environ["MP_API_KEY"], "accept": "application/json",
                 "User-Agent": "opendisorder-research/0.1"},
    )
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(3 * (attempt + 1))
                continue
            raise RuntimeError(f"MP API HTTP {e.code}: {e.read()[:200]!r}") from None
        except Exception:
            time.sleep(3 * (attempt + 1))
    raise RuntimeError("MP API unreachable")


def _cache_file(cs):
    return os.path.join(MP_CACHE, cs + ".json.gz")


def fetch_mp_chemsys(chemsys_list, batch=40):
    """Download (and cache) all non-deprecated MP summary docs for the given chemical systems."""
    os.makedirs(MP_CACHE, exist_ok=True)
    todo = [cs for cs in sorted(set(chemsys_list)) if not os.path.exists(_cache_file(cs))]
    for i in range(0, len(todo), batch):
        chunk = todo[i : i + batch]
        docs, skip = [], 0
        while True:
            r = _get("/materials/summary/", chemsys=",".join(chunk), _fields=FIELDS, _limit=1000,
                     _skip=skip, deprecated="false", id_format="legacy")
            docs += r["data"]
            skip += len(r["data"])
            if not r["data"] or skip >= r["meta"]["total_doc"]:
                break
        by = {cs: [] for cs in chunk}
        for d in docs:
            cs = "-".join(sorted({sp["element"] for s in d["structure"]["sites"] for sp in s["species"]}))
            by.setdefault(cs, []).append(d)
        for cs in chunk:
            with gzip.open(_cache_file(cs), "wt") as fh:
                json.dump(by.get(cs, []), fh)
        print(f"MP fetch {i + len(chunk)}/{len(todo)} chemsys, {len(docs)} docs", flush=True)


@lru_cache(maxsize=4096)
def mp_docs(chemsys: str):
    fn = _cache_file(chemsys)
    if not os.path.exists(fn):
        fetch_mp_chemsys([chemsys])
    with gzip.open(fn, "rt") as fh:
        return json.load(fh)


def mp_structure(doc) -> Structure:
    s = Structure.from_dict(doc["structure"])
    s.remove_oxidation_states()
    for p in list(s.site_properties):
        s.remove_site_property(p)
    return s
