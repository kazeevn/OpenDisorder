#!/usr/bin/env python
"""Step 2: generate ordered realizations ("colorings") of the selected disordered COD entries.

For every parent, supercell size and method, ordered structures are written as JSON lines
(one ``Structure.as_dict()`` per record, plus metadata) to ``--out`` (``.jsonl.gz``).

Supercells are built on the primitive cell of the parsed CIF. ``size`` counts primitive
cells; for each size the ``--shapes`` most compact Hermite-normal-form supercells are used.
Per symmetry orbit, the number of atoms of each element (and of vacancies) is rounded
with the largest-remainder method; ``exact`` records whether no rounding was needed.

Methods (``--methods``, comma-separated):
  random      uniform random coloring of every disordered orbit (pure pymatgen/numpy)
  ewald       "MP-like": oxidation states guessed -> OrderDisorderedStructureTransformation,
              lowest-Ewald-energy orderings first (the recipe MP applied to disordered ICSD
              entries, see README); falls back to random if no oxidation states can be found
  enum        symmetry-distinct derivative structures from enumlib (Hart & Forcade) via
              pymatgen's EnumerateStructureTransformation, a random subset of n (needs enum.x
              and makestr.x on PATH, e.g. ``conda install -c conda-forge enumlib``); the
              MatterGen benchmark used this
  enum_ewald  as enum, but the n lowest-Ewald structures
  sqs         special quasirandom structure from sqsgenerator (MIT, pip wheels): one per shape,
              pair correlations of the first shells closest to the random alloy

Negative/robustness controls (``--controls N`` positives per parent):
  ctrl_rattle       Gaussian displacements (--rattle A) and a random +-2% strain; expected to
                    stay a DISORDERED_REALIZATION (a cheap stand-in for relaxation noise)
  ctrl_composition  atoms of the mixed sublattice converted to the majority element (or
                    vacancies filled) until the global composition leaves MatterGen's tolerance;
                    expected DISORDERED_FAMILY or NONE
  ctrl_antisite     one atom of a mixed/vacancy sublattice swapped with an atom of an ordered
                    sublattice that cannot host it (composition unchanged); expected not a
                    realization (species mismatch)

usage:
  python generate_orderings.py --selected selected.csv --cod-db cod.sqlite --out orderings.jsonl.gz \
      --methods random,ewald,sqs --sizes auto --n-per-size 5 --controls 2 [--shard 0/8]
"""

from __future__ import annotations

import argparse
import collections
import csv
import logging
import math
import os
import shutil
import sys
import tempfile
import time
import warnings
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

log = logging.getLogger("generate")

PLACEHOLDERS = ["He", "Ne", "Ar", "Kr", "Xe", "Rn", "Fr", "Ra", "Po", "At", "Es", "Fm", "Md", "No", "Lr"]


# ----------------------------------------------------------------------------- supercells


def hnf_matrices(n: int):
    """All lower-triangular Hermite normal forms with determinant n."""
    for a in range(1, n + 1):
        if n % a:
            continue
        for c in range(1, n // a + 1):
            if (n // a) % c:
                continue
            f = n // (a * c)
            for b in range(c):
                for d in range(f):
                    for e in range(f):
                        yield np.array([[a, 0, 0], [b, c, 0], [d, e, f]])


def compactness(matrix: np.ndarray) -> float:
    """min perpendicular width / V^(1/3) of the lattice (1 for a cube)."""
    a, b, c = matrix
    v = abs(np.linalg.det(matrix))
    widths = [v / np.linalg.norm(np.cross(b, c)), v / np.linalg.norm(np.cross(c, a)), v / np.linalg.norm(np.cross(a, b))]
    return min(widths) / v ** (1 / 3)


def compact_supercells(lattice, size: int, n_shapes: int) -> list[np.ndarray]:
    """The ``n_shapes`` most compact, mutually distinct supercell matrices of ``size`` cells."""
    from pymatgen.core import Lattice

    scored = []
    seen = set()
    for h in hnf_matrices(size):
        m = h @ lattice.matrix
        red = Lattice(m).get_niggli_reduced_lattice()
        key = tuple(np.round(sorted(red.abc), 2)) + tuple(np.round(sorted(red.angles), 0))
        if key in seen:
            continue
        seen.add(key)
        scored.append((compactness(red.matrix), h))
    scored.sort(key=lambda t: -t[0])
    return [h for _, h in scored[:n_shapes]]


# ----------------------------------------------------------------------------- occupancy bookkeeping


def largest_remainder(targets: dict, total: int) -> dict:
    floors = {k: math.floor(v + 1e-9) for k, v in targets.items()}
    left = total - sum(floors.values())
    order = sorted(targets, key=lambda k: -(targets[k] - floors[k]))
    for k in order[: max(0, left)]:
        floors[k] += 1
    return floors


class Supercell:
    """A disordered supercell with per-orbit integer species counts."""

    def __init__(self, prim, orbit_of_site: np.ndarray, matrix: np.ndarray):
        p = prim.copy()
        p.add_site_property("orbit", [int(x) for x in orbit_of_site])
        self.sc = p * matrix
        self.matrix = matrix
        self.orbits = collections.defaultdict(list)
        for i, site in enumerate(self.sc):
            self.orbits[site.properties["orbit"]].append(i)
        self.counts: dict[int, dict] = {}
        self.max_error = 0.0
        for o, idx in self.orbits.items():
            site = self.sc[idx[0]]
            if site.is_ordered:
                continue
            occ = {el.symbol: float(x) for el, x in site.species.items()}
            vac = max(0.0, 1.0 - sum(occ.values()))
            targets = {**occ, None: vac}
            n = len(idx)
            cnt = largest_remainder({k: v * n for k, v in targets.items()}, n)
            self.counts[o] = cnt
            for k, v in targets.items():
                self.max_error = max(self.max_error, abs(cnt[k] / n - v))

    @property
    def exact(self) -> bool:
        return self.max_error < 1e-3

    def disordered_orbits(self):
        return list(self.counts)

    def rounded_disordered(self):
        """Copy of the supercell whose disordered orbits carry the rounded occupancies."""
        s = self.sc.copy()
        for o, cnt in self.counts.items():
            n = len(self.orbits[o])
            sp = {k: v / n for k, v in cnt.items() if k is not None and v > 0}
            for i in self.orbits[o]:
                s.replace(i, sp, properties={"orbit": o})
        return s

    def realize(self, assignment: dict[int, list]):
        """Ordered structure from {orbit: [species or None per site]}; returns (structure, vacant coords)."""
        from pymatgen.core import Structure

        species, coords, orbit, vacant = [], [], [], []
        for i, site in enumerate(self.sc):
            o = site.properties["orbit"]
            if o in assignment:
                el = assignment[o][self.orbits[o].index(i)]
                if el is None:
                    vacant.append((o, site.frac_coords))
                    continue
            else:
                el = site.specie.symbol
            species.append(el)
            coords.append(site.frac_coords)
            orbit.append(o)
        s = Structure(self.sc.lattice, species, coords, site_properties={"orbit": orbit})
        return s, vacant


# ----------------------------------------------------------------------------- oxidation states


def integer_composition(structure, max_den: int = 48):
    """Smallest integer multiple of the (rational) composition, e.g. Er2/3Cu2S2 -> Er2Cu6S6."""
    from fractions import Fraction

    from pymatgen.core import Composition

    fr = {el: Fraction(float(x)).limit_denominator(max_den) for el, x in structure.composition.items()}
    lcm = math.lcm(*[f.denominator for f in fr.values()])
    ints = {el: int(f * lcm) for el, f in fr.items()}
    g = math.gcd(*ints.values())
    return Composition({el: v // g for el, v in ints.items()})


def guess_oxidation_states(structure, timeout: float = 20.0) -> dict | None:
    """Element -> oxidation state for the parent composition (None if no charge-balanced guess).

    Uses pymatgen's ``oxi_state_guesses`` on the integer multiple of the rational composition
    (it fails on fractional formulas such as Er0.6667Cu2S2), then bond valence as a fallback."""
    try:
        with common.time_limit(timeout):
            guesses = integer_composition(structure).oxi_state_guesses(max_sites=-1)
    except Exception:
        guesses = []
    if guesses:
        return {el: round(v) for el, v in guesses[0].items()}
    try:
        from pymatgen.analysis.bond_valence import BVAnalyzer

        with common.time_limit(timeout):
            dec = BVAnalyzer().get_oxi_state_decorated_structure(structure)
        out = {}
        for site in dec:
            for sp in site.species:
                out[sp.symbol] = round(sp.oxi_state)
        return out
    except Exception:
        return None


def ewald_energy(structure, oxi: dict | None) -> float | None:
    if not oxi:
        return None
    from pymatgen.analysis.ewald import EwaldSummation

    try:
        s = structure.copy()
        for p in list(s.site_properties):
            s.remove_site_property(p)
        s.add_oxidation_state_by_element(oxi)
        return EwaldSummation(s).total_energy / len(s)
    except Exception:
        return None


# ----------------------------------------------------------------------------- methods


def gen_random(sc: Supercell, n: int, rng: np.random.Generator, ctx: dict):
    out = []
    for _ in range(n):
        assign = {}
        for o, cnt in sc.counts.items():
            lst = [k for k, v in cnt.items() for _ in range(v)]
            rng.shuffle(lst)
            assign[o] = lst
        s, vac = sc.realize(assign)
        out.append((s, vac, {}))
    return out


def _assignment_from_ordered(sc: Supercell, ordered) -> dict:
    """Map an ordered structure back onto the supercell sites (by fractional coordinates)."""
    from pymatgen.util.coord import pbc_shortest_vectors

    assign = {o: [None] * len(idx) for o, idx in sc.orbits.items() if o in sc.counts}
    fc_sc = sc.sc.frac_coords
    for site in ordered:
        d = np.linalg.norm(pbc_shortest_vectors(sc.sc.lattice, site.frac_coords, fc_sc)[0], axis=-1)
        i = int(np.argmin(d))
        o = sc.sc[i].properties["orbit"]
        if o in assign:
            assign[o][sc.orbits[o].index(i)] = site.specie.symbol
    return assign


def gen_ewald(sc: Supercell, n: int, rng, ctx: dict):
    from pymatgen.transformations.standard_transformations import OrderDisorderedStructureTransformation

    oxi = ctx.get("oxi")
    if not oxi:
        return [(s, v, {"note": "no_oxidation_states_random_fallback"}) for s, v, _ in gen_random(sc, n, rng, ctx)]
    s = sc.rounded_disordered()
    s.remove_site_property("orbit")
    s.add_oxidation_state_by_element(oxi)
    algo = {"fast": 0, "complete": 1, "best_first": 2}[ctx["args"].ewald_algo]
    t = OrderDisorderedStructureTransformation(algo=algo)
    with common.time_limit(ctx["args"].method_timeout):
        ranked = t.apply_transformation(s, return_ranked_list=n)
    out = []
    for r in ranked[:n]:
        o = r["structure"].copy()
        o.remove_oxidation_states()
        assign = _assignment_from_ordered(sc, o)
        real, vac = sc.realize(assign)
        out.append((real, vac, {"ewald_energy_total": r.get("energy")}))
    return out


def _enum_worker(prim_dict: dict, size: int, oxi, lowest: bool, timeout_min: float, workdir: str):
    """Runs in a memory-limited child: all symmetry-distinct derivative structures of ``size``."""
    import warnings as _w

    os.chdir(workdir)  # enumlib/pymatgen write scratch files into the working directory

    from pymatgen.core import Structure
    from pymatgen.transformations.advanced_transformations import EnumerateStructureTransformation

    _w.filterwarnings("ignore")
    prim = Structure.from_dict(prim_dict)
    if oxi:
        prim.add_oxidation_state_by_element(oxi)
    t = EnumerateStructureTransformation(
        min_cell_size=size,
        max_cell_size=size,
        symm_prec=0.1,
        sort_criteria="ewald" if (oxi and lowest) else "nsites",
        timeout=timeout_min,
        n_jobs=1,
    )
    ranked = t.apply_transformation(prim, return_ranked_list=10**9)
    out = []
    for r in ranked:
        s = r["structure"].copy()
        s.remove_oxidation_states()
        out.append((s.as_dict(verbosity=0), r.get("energy")))
    return out


def gen_enum(sc: Supercell, n: int, rng, ctx: dict, lowest: bool = False):
    from pymatgen.core import Structure

    if not (shutil.which("enum.x") or shutil.which("multienum.x")) or not shutil.which("makestr.x"):
        raise RuntimeError("enumlib executables (enum.x/multienum.x, makestr.x) not on PATH")
    a = ctx["args"]
    size = int(round(abs(np.linalg.det(sc.matrix))))
    key = (size, lowest)
    cache = ctx.setdefault("enum_cache", {})
    if key not in cache:
        workdir = tempfile.mkdtemp(prefix="rt_enum_")
        try:
            cache[key] = common.run_limited(
                _enum_worker,
                (ctx["prim_rational"].as_dict(), size, ctx.get("oxi"), lowest, a.method_timeout / 60, workdir),
                timeout=a.method_timeout + 30,
                mem_gb=a.enum_mem_gb,
            )
        finally:
            shutil.rmtree(workdir, ignore_errors=True)
    ranked = cache[key]
    total = len(ranked)
    if lowest:
        picks = list(range(min(n, total)))
    else:
        picks = sorted(rng.choice(total, size=min(n, total), replace=False).tolist())
    out = []
    for i in picks:
        d, energy = ranked[i]
        # enumlib picks its own supercell shape (vacant coordinates unknown -> no controls)
        out.append((Structure.from_dict(d), None, {"enum_index": int(i), "enum_total": total, "ewald_energy_total": energy}))
    return out


def gen_sqs(sc: Supercell, n: int, rng, ctx: dict):
    import sqsgenerator as sg

    if len(sc.counts) > len(PLACEHOLDERS):
        raise RuntimeError("too many disordered orbits for placeholder species")
    ph = {o: PLACEHOLDERS[i] for i, o in enumerate(sc.counts)}
    species = [ph.get(site.properties["orbit"], site.specie.symbol if site.is_ordered else None) for site in sc.sc]
    composition = []
    for o, cnt in sc.counts.items():
        entry = {"sites": ph[o]}
        for k, v in cnt.items():
            if v > 0:
                entry["0" if k is None else k] = int(v)
        composition.append(entry)
    cfg = {
        "iterations": int(ctx["args"].sqs_iterations),
        "sublattice_mode": "split" if len(composition) > 1 else "interact",
        "structure": {
            "supercell": [1, 1, 1],
            "species": species,
            "lattice": sc.sc.lattice.matrix.tolist(),
            "coords": sc.sc.frac_coords.tolist(),
        },
        "composition": composition,
    }
    with common.time_limit(ctx["args"].method_timeout):
        pack = sg.optimize(cfg, level=sg.core.LogLevel.error)
    best = pack.best()
    s = sg.to_pymatgen(best.structure())
    assign = _assignment_from_ordered(sc, s)
    real, vac = sc.realize(assign)
    return [(real, vac, {"sqs_objective": float(np.sum(best.objective))})]


METHODS = {
    "random": gen_random,
    "ewald": gen_ewald,
    "mp_like": gen_ewald,
    "enum": lambda sc, n, rng, ctx: gen_enum(sc, n, rng, ctx, lowest=False),
    "enum_ewald": lambda sc, n, rng, ctx: gen_enum(sc, n, rng, ctx, lowest=True),
    "sqs": gen_sqs,
}


# ----------------------------------------------------------------------------- controls


def ctrl_rattle(s, rng, sigma: float, strain: float):
    t = s.copy()
    eps = rng.uniform(-strain, strain, size=(3, 3))
    eps = (eps + eps.T) / 2
    t.apply_strain(eps, inplace=True) if hasattr(t, "apply_strain") else None
    disp = rng.normal(0, sigma, size=(len(t), 3))
    for i, d in enumerate(disp):
        t.translate_sites([i], d, frac_coords=False, to_unit_cell=True)
    return t


def ctrl_composition(s, vacant, parent_comp, counts: dict, rng):
    """Push the composition of the disordered sublattices out of MatterGen's tolerance."""
    from opendisorder.matching import fractional_composition_close

    t = s.copy()
    orbits = np.array(t.site_properties["orbit"])
    vacant = list(vacant or [])
    changed = 0
    for o, cnt in counts.items():
        occ = {k: v for k, v in cnt.items() if v > 0}
        els = [k for k in occ if k is not None]
        if not els:
            continue
        major = max(els, key=lambda k: occ[k])
        idx = [i for i in np.where(orbits == o)[0] if t[int(i)].specie.symbol != major]
        rng.shuffle(idx)
        for i in idx:  # minority element -> majority element
            t.replace(int(i), major, properties={"orbit": o})
            changed += 1
            if not fractional_composition_close(t.composition, parent_comp):
                return t, changed
        for oo, fc in [v for v in vacant if v[0] == o]:  # fill vacancies
            t.append(major, fc, properties={"orbit": oo})
            changed += 1
            if not fractional_composition_close(t.composition, parent_comp):
                return t, changed
    # other direction: remove atoms of the disordered sublattices (more vacancies)
    t = s.copy()
    changed = 0
    idx = [i for i in range(len(t)) if orbits[i] in counts]
    rng.shuffle(idx)
    remove = []
    for i in idx:
        remove.append(int(i))
        changed += 1
        trial = t.copy()
        trial.remove_sites(remove)
        if len(trial) and not fractional_composition_close(trial.composition, parent_comp):
            return trial, changed
    return None, changed


def ctrl_antisite(s, counts: dict, prim_species: dict, rng):
    """Swap a mixed-sublattice atom with an ordered-sublattice atom that the other sublattice cannot host."""
    t = s.copy()
    orbits = np.array(t.site_properties["orbit"])
    dis = [i for i in range(len(t)) if orbits[i] in counts]
    ordd = [i for i in range(len(t)) if orbits[i] not in counts]
    rng.shuffle(dis)
    rng.shuffle(ordd)
    for i in dis:
        a = t[i].specie.symbol
        for j in ordd:
            b = t[j].specie.symbol
            if a == b or b in prim_species[orbits[i]] or a in prim_species[orbits[j]]:
                continue
            fi, fj = t[i].frac_coords.copy(), t[j].frac_coords.copy()
            t.replace(i, b, coords=fi, properties={"orbit": int(orbits[i])})
            t.replace(j, a, coords=fj, properties={"orbit": int(orbits[j])})
            return t
    return None


# ----------------------------------------------------------------------------- driver


def parse_sizes(text: str, min_super: int, nsites_prim: int, max_atoms: int, n_sizes: int) -> list[int]:
    if text == "auto":
        sizes = [min_super * k for k in range(1, 64) if min_super * k * nsites_prim <= max_atoms]
        return sizes[:n_sizes]
    return [int(x) for x in text.split(",") if int(x) * nsites_prim <= max_atoms]


def make_record(rid, parent_row, s, extra, **meta) -> dict:
    from pymatgen.core import Structure

    clean = Structure(s.lattice, s.species, s.frac_coords)
    rec = {
        "id": rid,
        "parent": parent_row["cod_id"],
        "chemsys": parent_row.get("chemsys", ""),
        "disorder_type": parent_row.get("disorder_type", ""),
        "structure_class": parent_row.get("structure_class", ""),
        "n_atoms": len(clean),
        "formula": clean.composition.reduced_formula,
        **meta,
        **{k: v for k, v in extra.items() if v is not None},
        "structure": clean.as_dict(verbosity=0),
    }
    return rec


def process_parent(row: dict, conn, a, writer) -> collections.Counter:
    from pymatgen.analysis.structure_matcher import StructureMatcher

    stats = collections.Counter()
    cid = row["cod_id"]
    parent = common.cod_structure(conn, cid)
    if parent is None or parent.is_ordered:
        stats["parent_unusable"] += 1
        return stats
    prim = parent.get_primitive_structure()
    orbits = common.orbit_labels(prim)
    rng = np.random.default_rng([a.seed, int(cid)])
    # primitive cell with the rational occupancies (e.g. 0.66 -> 2/3): supercell counts are
    # derived from these; enumlib also needs exact fractions.  The analysis must use the same
    # primitive cell / orbit labels as below.
    info = common.analyze_disorder(prim, assume_primitive=True)
    prim_r = prim.copy()
    for o in info["orbits"]:
        if o["kind"] == "ordered":
            continue
        sp = {k: float(Fraction(v)) for k, v in o["rational"].items() if k != "Vac"}
        for i in np.where(orbits == o["orbit"])[0]:
            prim_r.replace(int(i), sp)
    ctx = {"args": a, "prim_rational": prim_r}
    ctx["oxi"] = guess_oxidation_states(prim_r) if a.oxidation_states else None
    prim_species = {}
    for i, site in enumerate(prim):
        prim_species[int(orbits[i])] = {el.symbol for el in site.species.elements}
    min_super = int(row.get("min_supercell") or info["min_supercell"])
    sizes = parse_sizes(a.sizes, min_super, len(prim), a.max_atoms, a.n_sizes)
    if not sizes:
        stats["no_size_fits"] += 1
        return stats
    sm = StructureMatcher(primitive_cell=False, attempt_supercell=False, scale=False)
    positives = []  # (record, structure with orbit labels, vacant, counts)
    for size in sizes:
        for h_i, h in enumerate(compact_supercells(prim.lattice, size, a.shapes)):
            sc = Supercell(prim_r, orbits, h)
            for method in a.methods:
                if method.startswith("enum") and h_i > 0:
                    continue  # enumlib already covers all supercell shapes of this size
                t0 = time.time()
                try:
                    results = METHODS[method](sc, a.n_per_size, rng, ctx)
                except common.Timeout:
                    stats[f"{method}_timeout"] += 1
                    continue
                except Exception as exc:
                    stats[f"{method}_error"] += 1
                    log.warning("%s %s size %d: %s: %s", cid, method, size, type(exc).__name__, str(exc)[:200])
                    continue
                dt = (time.time() - t0) / max(1, len(results))
                kept: list = []
                for s, vac, extra in results:
                    if a.dedupe and any(sm.fit(s, k[0]) for k in kept):
                        stats[f"{method}_duplicate"] += 1
                        continue
                    kept.append((s, vac, extra))
                energies = [ewald_energy(s, ctx["oxi"]) for s, _, _ in kept]
                order = np.argsort([e if e is not None else np.inf for e in energies])
                rank = {int(j): r for r, j in enumerate(order)}
                for j, (s, vac, extra) in enumerate(kept):
                    rid = f"{cid}-{method}-s{size}-h{h_i}-{j}"
                    rec = make_record(
                        rid,
                        row,
                        s,
                        extra,
                        kind="positive",
                        expected="realization",
                        method=method,
                        size=size,
                        shape=h_i,
                        supercell_matrix=h.tolist(),
                        exact=sc.exact,
                        composition_error=round(sc.max_error, 4),
                        rational_error=info["max_rational_error"],
                        ewald_energy=energies[j],
                        energy_rank=rank[j] if energies[j] is not None else None,
                        gen_time_s=round(dt, 3),
                    )
                    writer.write(rec)
                    stats[f"{method}"] += 1
                    if vac is not None:
                        positives.append((rec, s, vac, sc.counts))

    # ---- controls
    if a.controls and positives:
        picks = rng.choice(len(positives), size=min(a.controls, len(positives)), replace=False)
        for p in picks:
            rec, s, vac, counts = positives[int(p)]
            base = {k: rec[k] for k in ("method", "size", "shape", "exact")}
            t = ctrl_rattle(s, rng, a.rattle, a.strain)
            writer.write(make_record(rec["id"] + "-rattle", row, t, {}, kind="ctrl_rattle", expected="realization", source_id=rec["id"], **base))
            stats["ctrl_rattle"] += 1
            t, nchg = ctrl_composition(s, vac, parent.composition, counts, rng)
            if t is not None:
                writer.write(
                    make_record(rec["id"] + "-comp", row, t, {"n_changed": nchg}, kind="ctrl_composition", expected="not_realization", source_id=rec["id"], **base)
                )
                stats["ctrl_composition"] += 1
            t = ctrl_antisite(s, counts, prim_species, rng)
            if t is not None:
                writer.write(make_record(rec["id"] + "-anti", row, t, {}, kind="ctrl_antisite", expected="not_realization", source_id=rec["id"], **base))
                stats["ctrl_antisite"] += 1
    return stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selected", required=True)
    ap.add_argument("--cod-db", default=common.DEFAULT_COD_DB)
    ap.add_argument("--out", required=True, help="orderings (.jsonl or .jsonl.gz); appended, resumable per parent")
    ap.add_argument("--methods", default="random,ewald,sqs")
    ap.add_argument("--sizes", default="auto", help="'auto' (multiples of the smallest exact supercell) or e.g. 1,2,3,4")
    ap.add_argument("--n-sizes", type=int, default=3, help="with --sizes auto: how many multiples")
    ap.add_argument("--max-atoms", type=int, default=120, help="largest supercell (sites, before vacancies)")
    ap.add_argument("--shapes", type=int, default=1, help="supercell shapes per size")
    ap.add_argument("--n-per-size", type=int, default=5, help="orderings per (size, shape, method)")
    ap.add_argument("--controls", type=int, default=2, help="positives per parent that get negative/robustness controls")
    ap.add_argument("--rattle", type=float, default=0.05, help="Angstrom, std of ctrl_rattle displacements")
    ap.add_argument("--strain", type=float, default=0.02)
    ap.add_argument("--ewald-algo", default="fast", choices=["fast", "complete", "best_first"])
    ap.add_argument("--no-oxidation-states", dest="oxidation_states", action="store_false")
    ap.add_argument("--enum-mem-gb", type=float, default=4.0, help="address-space limit of the enumlib child process")
    ap.add_argument("--sqs-iterations", type=float, default=2e5)
    ap.add_argument("--method-timeout", type=float, default=300.0, help="seconds per (parent, size, shape, method)")
    ap.add_argument("--no-dedupe", dest="dedupe", action="store_false")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--shard", help="i/N (or slurm/N): only parents with crc32(cod_id) %% N == i")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)
    a.methods = [m.strip() for m in a.methods.split(",") if m.strip()]
    for m in a.methods:
        if m not in METHODS:
            ap.error(f"unknown method {m}; choose from {sorted(METHODS)}")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    warnings.filterwarnings("ignore")

    shard, nshards = common.parse_shard(a.shard)
    with open(a.selected, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if common.stable_shard(r["cod_id"], nshards) == shard]
    if a.limit:
        rows = rows[: a.limit]
    done = {rec["parent"] for rec in common.iter_jsonl(a.out)} if Path(a.out).exists() else set()
    conn = common.cod_connect(a.cod_db)
    total = collections.Counter()
    with common.JsonlWriter(a.out) as w:
        for i, row in enumerate(rows):
            if row["cod_id"] in done:
                continue
            t0 = time.time()
            try:
                st = process_parent(row, conn, a, w)
            except Exception as exc:
                log.exception("parent %s failed: %s", row["cod_id"], exc)
                st = collections.Counter(parent_error=1)
            total.update(st)
            log.info("[%d/%d] %s %s: %s (%.1fs)", i + 1, len(rows), row["cod_id"], row.get("reduced_formula", ""), dict(st), time.time() - t0)
    log.info("done: %s", dict(total))


if __name__ == "__main__":
    main()
