"""First-principles reconstruction of a compositionally disordered parent from ONE ordered structure.

No provenance, no external database: only the geometry and chemistry of the ordered structure.

Idea (orbit averaging).  An ordered realization S of a disordered parent P is a supercell of P in
which every mixed site of P is occupied by one of its species (or left empty).  If the species that
share a site in P are made indistinguishable, S recovers the symmetry (and the smaller primitive
cell) of P.  So, for each hypothesis H about which elements share sites:

1. collapse the elements of each group in H onto one representative ("colour-blind" structure);
   for a vacancy hypothesis, first complete the partially filled sublattice with dummy atoms at
   the images of its sites under the symmetry of the rest of the structure (the "framework");
2. find the space group of the collapsed structure (several symprec values: DFT-relaxed orderings
   are distorted away from the parent);
3. every crystallographic orbit gets the average occupation of the real species on it in S
   (dummies = vacancies);
4. symmetrise (spglib refined cell) and reduce to the primitive cell -> parent P_H.

A hypothesis is kept only if it produces at least one mixed / partially occupied orbit.
Hypotheses:
  * "similar": opendisorder.matching.similarity_groups (isovalent ions within 15 % Shannon radius,
    or Pettifor-close metals) -- the chemically most plausible substitutional disorder;
  * "pair": any two elements with the same sign of oxidation state (cation/cation, anion/anion),
    any two elements in intermetallics;
  * "vac": vacancies on the sublattice of one element (optionally together with a "pair"/"similar"
    group), e.g. oxygen-deficient perovskites, Li-deficient cathodes, cation-deficient spinels.

Ranking: chemical plausibility of the hypothesis first (similar < pair < vac < pair+vac), then the
supercell index |S| / |P| (larger = more of S is explained as an ordering), then the symprec needed
(smaller = S is less distorted from P).  `reconstruct()` returns all candidates, best first.
"""
from __future__ import annotations

import itertools, os, sys, warnings
from collections import Counter
from dataclasses import dataclass, field

import numpy as np
from pymatgen.core import Element, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

warnings.filterwarnings("ignore")
VAC = "He"  # dummy for vacancies (never in the inorganic structures we treat; checked below)
SYMPRECS = (0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
KIND_RANK = {"similar": 0, "pair": 1, "vac": 2, "similar+vac": 2, "pair+vac": 3}


@dataclass
class Candidate:
    kind: str
    groups: list
    vacancy: list | None
    symprec: float
    index: float
    spacegroup: int
    parent: Structure = field(repr=False)
    n_mixed_orbits: int = 0

    @property
    def key(self):
        return (KIND_RANK[self.kind], -round(self.index, 2), self.symprec, -self.spacegroup)

    def summary(self) -> str:
        return (f"{self.kind:11s} groups={self.groups} vac={self.vacancy} index={self.index:.2f} "
                f"symprec={self.symprec} -> {describe(self.parent)}")


# ----------------------------------------------------------------------------- chemistry
ANIONS = {"O", "S", "Se", "Te", "F", "Cl", "Br", "I", "N"}


def oxidation_roles(s: Structure) -> dict[str, int]:
    """Sign of the most probable oxidation state per element; {} for intermetallics.

    Falls back to "anion = the most electronegative chalcogen/halogen/N" when pymatgen's
    oxidation-state guess fails (e.g. mixed-valence Mn in KMn3O6)."""
    from opendisorder.matching import guess_ionic_radii

    ionic = guess_ionic_radii(s.composition)
    if len(ionic) == len(s.composition.elements):
        return {el.symbol: int(np.sign(ox)) for el, (ox, _) in ionic.items()}
    els = [el for el in s.composition.elements]
    an = [el for el in els if el.symbol in ANIONS]
    if not an:
        return {}
    xmax = max(el.X for el in an)
    return {el.symbol: (-1 if el.symbol in ANIONS and el.X >= xmax - 0.6 else 1) for el in els}


def hypotheses(s: Structure, allow_vacancies=True, max_vac_pairs=True):
    from opendisorder.matching import guess_ionic_radii, similarity_groups

    comp = s.composition
    els = sorted(el.symbol for el in comp.elements)
    roles = oxidation_roles(s)
    out = []
    sim = [sorted(e.symbol for e in g) for g in similarity_groups(comp.elements, ionic=guess_ionic_radii(comp))
           if len(g) > 1]
    if sim:
        out.append(("similar", sim, None))
    pairs = []
    for a, b in itertools.combinations(els, 2):
        if roles and roles.get(a) != roles.get(b):
            continue
        g = [[a, b]]
        if any(set(g[0]) <= set(x) for x in sim) and len(sim) == 1 and len(sim[0]) == 2:
            continue  # identical to the "similar" hypothesis
        pairs.append(g)
        out.append(("pair", g, None))
    # multi-site-group hypotheses: three same-role elements on one site, or two independent pairs
    for a, b, c in itertools.combinations(els, 3):
        if roles and len({roles.get(a), roles.get(b), roles.get(c)}) > 1:
            continue
        out.append(("pair", [[a, b, c]], None))
    for p1, p2 in itertools.combinations(pairs, 2):
        if not set(p1[0]) & set(p2[0]):
            out.append(("pair", p1 + p2, None))
    if allow_vacancies:
        for e in els:
            out.append(("vac", [], [e]))
        if max_vac_pairs:
            for g in ([("similar", sim)] if sim else []) + [("pair", p) for p in pairs]:
                for grp in g[1]:  # vacancies on the merged sublattice ...
                    out.append((g[0] + "+vac", g[1], list(grp)))
                merged = {e for grp in g[1] for e in grp}
                for e in els:  # ... or on another element's sublattice (e.g. pyrochlore O')
                    if e not in merged:
                        out.append((g[0] + "+vac", g[1], [e]))
    return out


# ----------------------------------------------------------------------------- geometry
def collapse(s: Structure, groups) -> Structure:
    rep = {el: sorted(g)[0] for g in groups for el in g}
    return Structure(s.lattice, [rep.get(x.specie.symbol, x.specie.symbol) for x in s], s.frac_coords)


def complete_vacancies(s: Structure, groups, vac, symprec, rel_dist=0.8):
    """Add dummy VAC atoms at the images of the `vac` sublattice under the framework's symmetry.

    Images closer than ``rel_dist`` x (shortest sublattice-sublattice / sublattice-framework
    distance in S) are rejected: they would be split positions or clashes, not vacancies.
    """
    vac = set(vac)
    frame_sites = [x for x in s if x.specie.symbol not in vac]
    occ = [x.frac_coords for x in s if x.specie.symbol in vac]
    if not frame_sites or not occ:
        return None
    frame = collapse(Structure.from_sites(frame_sites), groups)
    try:
        ops = SpacegroupAnalyzer(frame, symprec=symprec).get_symmetry_operations(cartesian=False)
    except Exception:
        return None
    lat = s.lattice
    d_ff = lat.get_all_distances(occ, occ)
    d_ff = d_ff[d_ff > 0.1].min(initial=9.0)
    d_fr = lat.get_all_distances(occ, frame.frac_coords).min()
    new = []
    for fc0 in occ:
        for op in ops:
            fc = np.mod(op.operate(fc0), 1.0)
            # images of distorted sites scatter by ~symprec around the ideal position
            if lat.get_all_distances([fc], occ + new).min() > max(0.5, 2.5 * symprec):
                new.append(fc)
    if not new or len(new) > 3 * len(occ):
        return None
    d_new = lat.get_all_distances(new, occ + new)
    d_new = d_new[d_new > 0.1].min(initial=9.0)
    if d_new < rel_dist * d_ff or lat.get_all_distances(new, frame.frac_coords).min() < rel_dist * d_fr:
        return None
    out = s.copy()
    for fc in new:
        out.append(VAC, fc)
    return out


def orbit_average(s: Structure, groups, symprec):
    c = collapse(s, groups)
    sga = SpacegroupAnalyzer(c, symprec=symprec)
    ss = sga.get_symmetrized_structure()
    species, n_mixed = [None] * len(s), 0
    for idxs in ss.equivalent_indices:
        cnt = Counter(s[i].specie.symbol for i in idxs)
        occ = {el: n / len(idxs) for el, n in cnt.items() if el != VAC}
        if len(occ) > 1 or (occ and sum(occ.values()) < 1 - 1e-6):
            n_mixed += 1
        for i in idxs:
            species[i] = occ
    keep = [i for i, sp in enumerate(species) if sp]
    d = Structure(s.lattice, [species[i] for i in keep], s.frac_coords[keep])
    return d, n_mixed


def idealise(d: Structure, symprec) -> Structure:
    try:
        ref = SpacegroupAnalyzer(d, symprec=symprec).get_refined_structure()
    except Exception:
        ref = d
    try:
        return ref.get_primitive_structure(tolerance=0.25)
    except Exception:
        return ref


def reconstruct(structure: Structure, allow_vacancies=True, symprecs=SYMPRECS) -> list[Candidate]:
    s = structure.copy()
    s.remove_oxidation_states()
    for p in list(s.site_properties):
        s.remove_site_property(p)
    if not s.is_ordered or Element(VAC) in s.composition:
        return []
    v0 = s.get_primitive_structure().volume  # index = volume ratio (vacancy sites add sites, not volume)
    out = []
    for kind, groups, vac in hypotheses(s, allow_vacancies):
        best = None
        for sp in symprecs:
            base, g = s, [list(x) for x in groups]
            if vac:
                base = complete_vacancies(s, g, vac, sp)
                if base is None:
                    continue
                if any(set(x) == set(vac) for x in g):
                    g = [x + [VAC] if set(x) == set(vac) else x for x in g]
                else:
                    g = g + [list(vac) + [VAC]]
            try:
                d, n_mixed = orbit_average(base, g, sp)
            except Exception:
                continue
            if n_mixed == 0:
                continue
            p = idealise(d, sp)
            if abs(p.composition.num_atoms / len(p) - 1) < 1e-6 and all(x.is_ordered for x in p):
                continue
            index = v0 / p.volume
            if index < 1 - 0.02:
                continue  # a parent must not be larger than the ordering it explains
            try:
                sgn = SpacegroupAnalyzer(p, 0.1).get_space_group_number()
            except Exception:
                sgn = 1
            cand = Candidate(kind, [sorted(x) for x in groups], sorted(vac) if vac else None, sp, index, sgn, p,
                             n_mixed)
            if best is None or (round(index, 2), sgn) > (round(best.index, 2), best.spacegroup):
                best = cand
        if best:
            out.append(best)
    out.sort(key=lambda c: c.key)
    # drop duplicates (same parent reached through different hypotheses)
    uniq, seen = [], set()
    for c in out:
        k = (c.parent.composition.reduced_formula, c.spacegroup, len(c.parent))
        if k not in seen:
            seen.add(k)
            uniq.append(c)
    return uniq


def describe(p: Structure) -> str:
    try:
        sga = SpacegroupAnalyzer(p, 0.1)
        ss = sga.get_symmetrized_structure()
        parts = []
        for sites, w in zip(ss.equivalent_sites, ss.wyckoff_symbols):
            parts.append(f"{w}:" + ",".join(f"{el.symbol}{amt:.3g}" for el, amt in sites[0].species.items()))
        return f"{sga.get_space_group_symbol()} ({len(p)} sites) | " + " ".join(parts)
    except Exception:
        return f"{p.composition.reduced_formula} ({len(p)} sites)"


if __name__ == "__main__":
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "parents")
    os.makedirs(out_dir, exist_ok=True)
    for fn in sys.argv[1:]:
        s = Structure.from_file(fn)
        tag = os.path.basename(fn).rsplit(".", 1)[0]
        res = reconstruct(s)
        print(f"== {tag}  ({len(s.get_primitive_structure())} sites prim, "
              f"{SpacegroupAnalyzer(s, 0.1).get_space_group_symbol()}), {len(res)} candidates")
        for i, c in enumerate(res[:6]):
            print(f"  [{i}] {c.summary()}")
            c.parent.to(filename=os.path.join(out_dir, f"{tag}__parent{i}.cif"))
