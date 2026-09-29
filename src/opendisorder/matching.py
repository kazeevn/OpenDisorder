"""Structure comparison primitives.

Two kinds of match are implemented:

* **Pristine match** -- the query and the reference are the same crystal
  structure (same reduced formula, pymatgen :class:`StructureMatcher` fit).
* **Disorder match** -- the (ordered) query is an ordered realization, i.e. a
  "coloring" / ordered approximant, of a reference that has partial
  occupancies.  Following Zeni et al. (MatterGen, Nature 2025, SI D.3) the
  screen is a :class:`StructureMatcher` with ``attempt_supercell``,
  ``allow_subset`` and :class:`OrderDisorderElementComparator`.  Because that
  screen alone produces false positives (it lets a whole sublattice go
  missing, or lets e.g. Li3FeO4 fit (Li0.5Fe0.5)O), every screen hit is
  verified on the explicit site mapping:

  a. each query site sits on a reference site whose species set contains the
     query element;
  b. reference sites left empty by the query are partially occupied
     (i.e. they have vacancy capacity);
  c. for each class of equivalent reference sites the element fractions of
     the query agree with the reference occupancies within a tolerance.

  If (a) and (b) hold but (c) fails, the query shares the framework and the
  site-species pattern of a known disordered phase but has a different
  stoichiometry (e.g. TaCr2O6 vs. Ta0.5Cr0.5O2, Juelsholt, Mater. Horiz.
  2026).  This is reported as a separate, weaker grade.

* **Sibling ordering** (optional) -- two *ordered* structures with the same
  reduced formula that become identical when chemically similar elements
  (Hume-Rothery-like rules: |dX| <= 1.0 and relative radius difference
  <= 0.3, as in MatterGen) are allowed to share sites.  The query is then a
  different ordering of the same parent structure.
"""

from __future__ import annotations

import itertools
import logging
from dataclasses import dataclass, field
from enum import Enum

import numpy as np
from pymatgen.analysis.structure_matcher import (
    ElementComparator,
    OrderDisorderElementComparator,
    StructureMatcher,
)
from pymatgen.core import Composition, Element, Structure

logger = logging.getLogger(__name__)


class MatchKind(str, Enum):
    """Graded outcome of comparing a query to a reference, strongest first."""

    PRISTINE = "pristine"
    """Same structure (and composition) as an ordered reference entry."""
    DISORDERED_REALIZATION = "disordered_realization"
    """Ordered realization of a disordered reference with consistent composition."""
    DISORDERED_FAMILY = "disordered_family"
    """Same framework and site-species pattern as a disordered reference, but a
    different stoichiometry (e.g. another member of a solid-solution series)."""
    SUBSTITUTIONAL_ORDERING = "substitutional_ordering"
    """Ordering of chemically similar elements over the sites of a known
    structure with fewer elements (e.g. GNoME's TbHo2Tm9Co4 = ordered
    (Tb,Ho,Tm)3Co on the Ho3Co structure, Cheetham & Seshadri 2024)."""
    SIBLING_ORDERING = "sibling_ordering"
    """Different ordering of chemically similar elements on the same parent
    structure as an ordered reference with the same formula."""
    NONE = "none"


#: Order in which match kinds are considered stronger (lower = stronger).
MATCH_STRENGTH = {k: i for i, k in enumerate(MatchKind)}


@dataclass
class MatchResult:
    kind: MatchKind
    reference_id: str | None = None
    details: dict = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.kind is not MatchKind.NONE


# --------------------------------------------------------------------------- helpers


def clean_structure(structure: Structure) -> Structure:
    """Copy of ``structure`` without oxidation states, magmoms or site properties."""
    s = structure.copy()
    s.remove_oxidation_states()
    for prop in list(s.site_properties):
        s.remove_site_property(prop)
    return s


def is_ordered(structure: Structure) -> bool:
    return structure.is_ordered


def chemsys(structure_or_comp: Structure | Composition) -> str:
    comp = structure_or_comp.composition if isinstance(structure_or_comp, Structure) else structure_or_comp
    return "-".join(sorted({el.symbol for el in comp.elements}))


def reduced_formula(structure_or_comp: Structure | Composition) -> str:
    comp = structure_or_comp.composition if isinstance(structure_or_comp, Structure) else structure_or_comp
    return comp.reduced_formula


def fractional_composition_close(c1: Composition, c2: Composition, atol: float = 0.01, rtol: float = 0.1) -> bool:
    """MatterGen's global composition criterion on fractional compositions."""
    f1, f2 = c1.fractional_composition, c2.fractional_composition
    for el in set(f1.elements) | set(f2.elements):
        a, b = f1[el], f2[el]
        if abs(a - b) > atol + rtol * (a + b) / 2:
            return False
    return True


def _site_key(site, ndigits: int = 3) -> tuple:
    """Hashable label of a (possibly disordered) site's species and occupancies."""
    return tuple(sorted((el.symbol, round(float(occ), ndigits)) for el, occ in site.species.items()))


class ElementSetComparator(OrderDisorderElementComparator):
    """Sites match if one's element set contains the other's; structures are
    pre-screened on their element sets rather than on fractional composition,
    so disordered structures with different occupancies can be compared."""

    def get_hash(self, composition):
        return "-".join(sorted({el.symbol for el in composition.elements}))


# --------------------------------------------------------------------------- pristine


@dataclass
class PristineMatcher:
    """Plain :class:`StructureMatcher` match with pymatgen/MatterGen defaults."""

    ltol: float = 0.2
    stol: float = 0.3
    angle_tol: float = 5.0

    def __post_init__(self):
        self._sm = StructureMatcher(
            ltol=self.ltol,
            stol=self.stol,
            angle_tol=self.angle_tol,
            primitive_cell=True,
            scale=True,
            attempt_supercell=False,
            comparator=ElementComparator(),
        )

    def fit(self, query: Structure, reference: Structure) -> bool:
        if query.composition.reduced_composition != reference.composition.reduced_composition:
            # ElementComparator ignores occupancies, so compare compositions explicitly;
            # reduced_composition handles fractional (disordered) compositions as well.
            if not fractional_composition_close(query.composition, reference.composition, atol=1e-3, rtol=0):
                return False
        try:
            return bool(self._sm.fit(query, reference))
        except Exception as exc:  # pragma: no cover - pymatgen corner cases
            logger.debug("PristineMatcher failed: %s", exc)
            return False

    def rms(self, query: Structure, reference: Structure):
        try:
            return self._sm.get_rms_dist(query, reference)
        except Exception:  # pragma: no cover
            return None


# --------------------------------------------------------------------------- disorder


@dataclass
class DisorderMatcher:
    """Match an ordered query against a disordered reference.

    Args:
        ltol, stol, angle_tol: StructureMatcher tolerances.
        occupancy_tol: maximum absolute deviation, per class of equivalent
            reference sites and per element, between the query's site fraction
            and the reference occupancy for a ``DISORDERED_REALIZATION``.
            Also used as the vacancy threshold (a reference site with total
            occupancy above ``1 - occupancy_tol`` must be filled).
        composition_atol, composition_rtol: global composition tolerance
            (MatterGen's criterion); both this and the per-site check must pass.
        max_supercell: largest ratio of query volume per formula to reference
            primitive-cell volume that is attempted; larger ratios are skipped.
    """

    ltol: float = 0.2
    stol: float = 0.3
    angle_tol: float = 5.0
    occupancy_tol: float = 0.1
    composition_atol: float = 0.01
    composition_rtol: float = 0.1
    max_supercell: int = 16

    def __post_init__(self):
        common = dict(
            ltol=self.ltol,
            stol=self.stol,
            angle_tol=self.angle_tol,
            scale=True,
            attempt_supercell=True,
            allow_subset=True,
            comparator=OrderDisorderElementComparator(),
            supercell_size="volume",
        )
        self._screen = StructureMatcher(primitive_cell=True, **common)
        self._mapper = StructureMatcher(primitive_cell=False, **common)
        self._screen_disordered = StructureMatcher(
            ltol=self.ltol,
            stol=self.stol,
            angle_tol=self.angle_tol,
            primitive_cell=True,
            scale=True,
            attempt_supercell=True,
            comparator=ElementSetComparator(),
        )

    # -- public -------------------------------------------------------------

    def match(self, query: Structure, reference: Structure) -> MatchResult:
        """Grade how ``query`` (ordered) relates to ``reference`` (disordered)."""
        if not set(query.composition.elements) <= set(reference.composition.elements):
            return MatchResult(MatchKind.NONE, details={"reason": "elements"})

        q = query.get_primitive_structure()
        r = reference.get_primitive_structure()
        ratio = (q.volume / r.volume) if r.volume > 0 else np.inf
        if ratio > self.max_supercell + 0.5:
            return MatchResult(MatchKind.NONE, details={"reason": "supercell too large", "ratio": ratio})

        try:
            if not self._screen.fit(q, r):
                return MatchResult(MatchKind.NONE, details={"reason": "no framework fit"})
            mapped = self._mapper.get_s2_like_s1(q, r)
        except Exception as exc:
            logger.debug("DisorderMatcher failed: %s", exc)
            return MatchResult(MatchKind.NONE, details={"reason": f"matcher error: {exc}"})
        if mapped is None:
            return MatchResult(MatchKind.NONE, details={"reason": "no mapping"})
        return self._verify(q, mapped)

    def match_disordered(self, query: Structure, reference: Structure) -> MatchResult:
        """Grade a *disordered* query against a disordered reference.

        Both must share the framework and the site-species pattern (each
        query site's element set contained in, or containing, the matched
        reference site's).  Compositions within the global tolerance give
        ``PRISTINE`` (the same disordered phase, occupancies refined
        differently); otherwise ``DISORDERED_FAMILY`` (e.g. Ta1/3Cr2/3O2 vs.
        Ta1/2Cr1/2O2).
        """
        if {el.symbol for el in query.composition.elements} != {el.symbol for el in reference.composition.elements}:
            return MatchResult(MatchKind.NONE, details={"reason": "elements"})
        try:
            fit = self._screen_disordered.fit(query, reference)
        except Exception as exc:
            logger.debug("match_disordered failed: %s", exc)
            return MatchResult(MatchKind.NONE, details={"reason": f"matcher error: {exc}"})
        if not fit:
            return MatchResult(MatchKind.NONE, details={"reason": "no framework fit"})
        comp_ok = fractional_composition_close(
            query.composition, reference.composition, self.composition_atol, self.composition_rtol
        )
        details = {"global_composition_ok": comp_ok, "reference_composition": reference.composition.formula}
        return MatchResult(MatchKind.PRISTINE if comp_ok else MatchKind.DISORDERED_FAMILY, details=details)

    # -- internals ----------------------------------------------------------

    def _verify(self, q: Structure, mapped: Structure) -> MatchResult:
        n = len(q)
        matched, unmatched = list(mapped)[:n], list(mapped)[n:]

        # (a) species compatibility
        for qs, rs in zip(q, matched):
            if qs.specie not in rs.species.elements:
                return MatchResult(MatchKind.NONE, details={"reason": "species mismatch"})
        # (b) empty reference sites must be able to be vacant
        for rs in unmatched:
            if rs.species.num_atoms > 1 - self.occupancy_tol:
                return MatchResult(MatchKind.NONE, details={"reason": "fully occupied site left empty"})

        # (c) per-site-class composition
        expected: dict[tuple, dict[str, float]] = {}
        actual: dict[tuple, dict[str, float]] = {}
        counts: dict[tuple, int] = {}
        for rs in mapped:
            key = _site_key(rs)
            counts[key] = counts.get(key, 0) + 1
            exp = expected.setdefault(key, {})
            for el, occ in rs.species.items():
                exp[el.symbol] = exp.get(el.symbol, 0.0) + float(occ)
        for qs, rs in zip(q, matched):
            act = actual.setdefault(_site_key(rs), {})
            act[qs.specie.symbol] = act.get(qs.specie.symbol, 0.0) + 1.0

        max_dev = 0.0
        for key, exp in expected.items():
            act = actual.get(key, {})
            for el in set(exp) | set(act):
                dev = abs(exp.get(el, 0.0) - act.get(el, 0.0)) / counts[key]
                max_dev = max(max_dev, dev)

        ref_comp = Composition({el: amt for el, amt in mapped.composition.items()})
        comp_ok = fractional_composition_close(q.composition, ref_comp, self.composition_atol, self.composition_rtol)
        details = {
            "max_site_occupancy_deviation": round(max_dev, 4),
            "supercell_sites": len(mapped),
            "reference_supercell_composition": ref_comp.formula,
            "global_composition_ok": comp_ok,
        }
        if max_dev <= self.occupancy_tol + 1e-8 and comp_ok:
            return MatchResult(MatchKind.DISORDERED_REALIZATION, details=details)
        return MatchResult(MatchKind.DISORDERED_FAMILY, details=details)


# --------------------------------------------------------------------------- substitutions


def _pettifor_similar(a: Element, b: Element, max_dm: float, max_dr: float) -> bool:
    try:
        ma, mb = a.mendeleev_no, b.mendeleev_no
        ra, rb = float(a.atomic_radius), float(b.atomic_radius)
    except (TypeError, ValueError, AttributeError):
        return False
    if ma is None or mb is None or not ra or not rb:
        return False
    return abs(ma - mb) <= max_dm and abs(ra - rb) / ((ra + rb) / 2) <= max_dr + 1e-9


def guess_ionic_radii(composition: Composition) -> dict[Element, tuple[int, float]]:
    """``{element: (oxidation state, Shannon radius for CN VI)}`` for the most
    probable charge-balanced oxidation states; empty for intermetallics or
    when no guess is possible."""
    from pymatgen.core import Species

    try:
        guesses = composition.oxi_state_guesses(max_sites=-1)
    except Exception:
        return {}
    if not guesses:
        return {}
    out = {}
    for el, ox in guesses[0].items():
        if abs(ox - round(ox)) > 1e-6 or round(ox) == 0:
            continue
        sp = Species(el, round(ox))
        try:
            r = sp.get_shannon_radius(cn="VI")
        except Exception:
            r = sp.ionic_radius
        if r:
            out[Element(el)] = (round(ox), float(r))
    return out


def _similar_pair(a: Element, b: Element, ionic: dict, max_dm: float, max_dr: float, max_ionic_dr: float) -> bool:
    if a in ionic and b in ionic:
        # ionic compound: isovalent, and either size-matched or chemically close
        (oa, ra), (ob, rb) = ionic[a], ionic[b]
        if oa != ob:
            return False
        return abs(ra - rb) / ((ra + rb) / 2) <= max_ionic_dr + 1e-9 or _pettifor_similar(a, b, max_dm, max_dr)
    return _pettifor_similar(a, b, max_dm, max_dr)


def similarity_groups(
    elements,
    max_dm: float = 6,
    max_dr: float = 0.2,
    ionic: dict | None = None,
    max_ionic_dr: float = 0.15,
) -> list[list[Element]]:
    """Partition elements into groups that can share crystallographic sites.

    Two elements are *Pettifor-similar* if their Pettifor (Mendeleev) numbers
    differ by at most ``max_dm`` and their atomic radii by at most ``max_dr``
    (relative).  If ``ionic`` (from :func:`guess_ionic_radii`) assigns
    oxidation states to both, they are similar when they are isovalent and
    either their ionic radii are within ``max_ionic_dr`` (Goldschmidt's 15%
    rule, e.g. Zr4+/Sn4+, Mg2+/Ni2+) or they are Pettifor-similar (e.g.
    K+/Na+); otherwise Pettifor similarity decides (intermetallics).

    Elements are visited in Pettifor order and join the first group in which
    they are similar to *every* member (complete linkage, no chaining).  The
    first element of a group is its representative.
    """
    ionic = ionic or {}
    els = sorted(set(elements), key=lambda e: (e.mendeleev_no or 0, e.Z))
    groups: list[list[Element]] = []
    for el in els:
        for g in groups:
            if all(_similar_pair(el, o, ionic, max_dm, max_dr, max_ionic_dr) for o in g):
                g.append(el)
                break
        else:
            groups.append([el])
    return groups


def collapse_elements(structure: Structure, groups: list[list[Element]]) -> Structure:
    """Replace every element of each group by the group's representative
    (occupancies of merged elements on one site are summed)."""
    mapping = {el: g[0] for g in groups for el in g[1:]}
    present = {el: rep for el, rep in mapping.items() if el in structure.composition}
    s = structure.copy()
    if present:
        s.replace_species(present)
    return s


# --------------------------------------------------------------------------- sibling orderings


def _similar(a: Element, b: Element, max_dx: float = 1.0, max_dr: float = 0.3) -> bool:
    try:
        xa, xb = a.X, b.X
        ra, rb = a.atomic_radius, b.atomic_radius
    except Exception:
        return False
    if xa is None or xb is None or ra is None or rb is None or np.isnan(xa) or np.isnan(xb):
        return False
    ra, rb = float(ra), float(rb)
    return abs(xa - xb) <= max_dx and abs(ra - rb) / ((ra + rb) / 2) <= max_dr


def substitution_cliques(elements, max_dx: float = 1.0, max_dr: float = 0.3) -> list[set[Element]]:
    """Group elements into transitive classes of Hume-Rothery-similar elements."""
    elements = sorted(set(elements))
    parent = {el: el for el in elements}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in itertools.combinations(elements, 2):
        if _similar(a, b, max_dx, max_dr):
            parent[find(a)] = find(b)
    groups: dict[Element, set[Element]] = {}
    for el in elements:
        groups.setdefault(find(el), set()).add(el)
    return [g for g in groups.values() if len(g) > 1]


def disorder_similar_elements(structure: Structure, max_dx: float = 1.0, max_dr: float = 0.3) -> Structure | None:
    """Replace every site of each similar-element clique by the clique's average occupancy.

    Returns ``None`` if the structure has no pair of similar elements.
    """
    cliques = substitution_cliques(structure.composition.elements, max_dx, max_dr)
    if not cliques:
        return None
    comp = structure.composition
    s = structure.copy()
    for clique in cliques:
        total = sum(comp[el] for el in clique)
        mix = {el: comp[el] / total for el in clique}
        for i, site in enumerate(s):
            if site.specie in clique:
                s.replace(i, mix, coords=site.frac_coords)
    return s


@dataclass
class SiblingOrderingMatcher:
    """MatterGen's ordered-vs-ordered route through Hume-Rothery disordering."""

    ltol: float = 0.2
    stol: float = 0.3
    angle_tol: float = 5.0
    max_dx: float = 1.0
    max_dr: float = 0.3

    def __post_init__(self):
        self._sm = StructureMatcher(
            ltol=self.ltol,
            stol=self.stol,
            angle_tol=self.angle_tol,
            primitive_cell=True,
            scale=True,
            attempt_supercell=True,
            allow_subset=True,
            comparator=OrderDisorderElementComparator(),
        )

    def prepare(self, query: Structure) -> Structure | None:
        return disorder_similar_elements(query, self.max_dx, self.max_dr)

    def fit(self, disordered_query: Structure | None, query: Structure, reference: Structure) -> bool:
        if disordered_query is None or query.composition.reduced_formula != reference.composition.reduced_formula:
            return False
        try:
            return bool(self._sm.fit(disordered_query, reference))
        except Exception as exc:  # pragma: no cover
            logger.debug("SiblingOrderingMatcher failed: %s", exc)
            return False
