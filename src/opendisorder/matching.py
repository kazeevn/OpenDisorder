"""Structure comparison primitives.

Two kinds of match are implemented:

* **Pristine match** -- the query and the reference are the same crystal
  structure (same reduced formula, pymatgen :class:`StructureMatcher` fit).
* **Disorder match** -- the (ordered) query is an ordered realization, i.e. a
  "coloring" / ordered approximant, of a reference that has partial
  occupancies.  Following Zeni et al. (MatterGen, Nature 2025, SI D.4) the
  screen is a :class:`StructureMatcher` with ``attempt_supercell``,
  ``allow_subset`` and :class:`OrderDisorderElementComparator`.  Because that
  screen alone produces false positives (it lets a whole sublattice go
  missing, or lets e.g. Li3FeO4 fit (Li0.5Fe0.5)O), every screen hit is
  verified on the explicit site mapping:

  a. each query site sits on a reference site whose species set contains the
     query element;
  b. reference sites left empty by the query are partially occupied
     (i.e. they have vacancy capacity);
  c. the global composition agrees (MatterGen's tolerance) and every
     sublattice (reference sites hosting the same set of elements) either
     - has query site fractions within ``occupancy_tol`` of the reference
       occupancies, or
     - is realized by its majority occupants: each query atom on a site where
       its element has the largest occupancy, each empty site one where the
       vacancy fraction is largest (the idealized ordering of minor anti-site
       disorder, e.g. the Na0.82Mn0.18 sites of alluaudites).

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


def drop_minor_species(structure: Structure, min_occupancy: float) -> Structure:
    """Remove species whose occupancy on a site is below ``min_occupancy``
    (sites left with nothing are removed).  Tiny occupancies -- trace dopants,
    refinement artefacts such as 0.75% Ti on an otherwise empty spinel site --
    cannot be realized by any small ordered supercell and would otherwise
    break matches."""
    if min_occupancy <= 0 or structure.is_ordered:
        return structure
    s = structure.copy()
    remove = []
    for i, site in enumerate(s):
        kept = {sp: occ for sp, occ in site.species.items() if occ >= min_occupancy}
        if not kept:
            remove.append(i)
        elif len(kept) != len(site.species):
            s.replace(i, kept, coords=site.frac_coords)
    if remove:
        s.remove_sites(remove)
    return s


#: Elements treated as anions when deciding whether a composition is ionic.
ANION_ELEMENTS = frozenset({"H", "N", "O", "F", "S", "Cl", "Se", "Br", "Te", "I"})


def has_anion(composition: Composition) -> bool:
    return any(el.symbol in ANION_ELEMENTS for el in composition.elements)


def volume_per_atom(structure: Structure) -> float:
    """Volume per atom, counting partial occupancies (not sites)."""
    n = structure.composition.num_atoms
    return structure.volume / n if n > 0 else float("nan")


def volume_ratio(query: Structure, reference: Structure) -> float:
    """Density ratio of ``query`` over ``reference``: volume per atom
    (occupancy-weighted) or, if closer to 1, volume per site.

    StructureMatcher's ``scale=True`` normalizes volumes away, so a match onto a
    reference with a very different density (a high-pressure phase, a
    misparsed CIF) is otherwise invisible.  Per atom is right for ordered
    realizations (atoms = occupancies); per site for relations that fill or
    empty vacancies (framework density).  The value closer to 1 is returned."""
    vr = volume_per_atom(reference)
    per_atom = volume_per_atom(query) / vr if vr > 0 else float("nan")
    per_site = (query.volume / len(query)) / (reference.volume / len(reference)) if len(reference) else float("nan")
    return min((per_atom, per_site), key=lambda x: abs(np.log(x)) if x > 0 else np.inf)


def is_generic(structure: Structure) -> bool:
    """Elemental or at most 2 sites per primitive cell (fcc/bcc/hcp, rock salt,
    CsCl ...): matches onto such references carry little information."""
    if len(structure.composition.elements) == 1:
        return True
    try:
        return len(structure.get_primitive_structure()) <= 2
    except Exception:  # pragma: no cover
        return len(structure) <= 2


def minority_fraction(structure: Structure) -> float:
    """Largest non-majority fraction on any site (vacancy counts as a species);
    below ~0.1 the structure is a dilute (terminal) solid solution."""
    worst = 0.0
    for site in structure:
        occ = [float(x) for x in site.species.values()]
        vacancy = max(0.0, 1.0 - sum(occ))
        worst = max(worst, 1.0 - max([*occ, vacancy]))
    return worst


def _element_set_key(site) -> tuple:
    return tuple(sorted(el.symbol for el in site.species.elements))


def _occupancy_deviations(q: Structure, matched, mapped: Structure, key) -> dict[tuple, tuple[float, bool]]:
    """Per class of reference sites (defined by ``key``): the largest per-element
    difference between the query's site fractions and the reference
    occupancies, and whether the query realizes the class by its *majority*
    occupants (each query atom on a site where its element has the largest
    occupancy, each empty site one where vacancy is the largest fraction)."""
    n = len(q)
    expected: dict[tuple, dict[str, float]] = {}
    actual: dict[tuple, dict[str, float]] = {}
    counts: dict[tuple, int] = {}
    majority: dict[tuple, bool] = {}
    for i, rs in enumerate(mapped):
        k = key(rs)
        counts[k] = counts.get(k, 0) + 1
        exp = expected.setdefault(k, {})
        occ = {el.symbol: float(x) for el, x in rs.species.items()}
        for el, x in occ.items():
            exp[el] = exp.get(el, 0.0) + x
        vacancy = max(0.0, 1.0 - sum(occ.values()))
        top = max([*occ.values(), vacancy])
        if i < n:
            el = q[i].specie.symbol
            actual.setdefault(k, {})[el] = actual.get(k, {}).get(el, 0.0) + 1.0
            ok = occ.get(el, 0.0) >= top - 1e-6
        else:
            ok = vacancy >= top - 1e-6
        majority[k] = majority.get(k, True) and ok
    out = {}
    for k, exp in expected.items():
        act = actual.get(k, {})
        dev = max(abs(exp.get(el, 0.0) - act.get(el, 0.0)) / counts[k] for el in set(exp) | set(act))
        out[k] = (dev, majority[k])
    return out


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
        occupancy_tol: maximum absolute deviation, per sublattice (reference
            sites hosting the same set of elements) and per element, between
            the query's site fraction and the reference occupancy for a
            ``DISORDERED_REALIZATION``.  0.2 accommodates the coarse fractions
            a small supercell can realize (multiples of 1/k on a k-site
            sublattice) and redistribution between sublattices; the global
            composition check below is the stoichiometry guard.
        min_occupancy: species with a smaller site occupancy in a disordered
            structure are ignored (see :func:`drop_minor_species`).
        min_volume_ratio, max_volume_ratio: allowed range of the query's
            volume per atom over the reference's (occupancy-weighted); outside
            it there is no match.
        family_max_deviation: a ``DISORDERED_FAMILY`` needs every sublattice
            within this occupancy deviation, and a reference that is not a
            dilute solid solution (some site with a minority fraction of at
            least ``dilute_threshold``).
        vacancy_tol: a reference site with total occupancy above
            ``1 - vacancy_tol`` must not be left empty by the query.
        composition_atol, composition_rtol: global composition tolerance
            (MatterGen's criterion); both this and the sublattice check must pass.
        max_supercell: largest ratio of query volume per formula to reference
            primitive-cell volume that is attempted; larger ratios are skipped.
    """

    ltol: float = 0.2
    stol: float = 0.3
    angle_tol: float = 5.0
    occupancy_tol: float = 0.2
    vacancy_tol: float = 0.1
    min_occupancy: float = 0.05
    min_volume_ratio: float = 0.75
    max_volume_ratio: float = 1.35
    family_max_deviation: float = 0.5
    dilute_threshold: float = 0.1
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
        reference = drop_minor_species(reference, self.min_occupancy)
        if not set(query.composition.elements) <= set(reference.composition.elements):
            return MatchResult(MatchKind.NONE, details={"reason": "elements"})

        vr = volume_ratio(query, reference)
        if not self.min_volume_ratio <= vr <= self.max_volume_ratio:
            return MatchResult(MatchKind.NONE, details={"reason": "volume per atom", "volume_ratio": round(vr, 3)})

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
        res = self._verify(q, mapped)
        if res.kind is MatchKind.DISORDERED_FAMILY and minority_fraction(r) < self.dilute_threshold:
            return MatchResult(MatchKind.NONE, details={"reason": "dilute reference", **res.details})
        res.details["volume_ratio"] = round(vr, 3)
        return res

    def match_disordered(self, query: Structure, reference: Structure) -> MatchResult:
        """Grade a *disordered* query against a disordered reference.

        Both must share the framework and the site-species pattern (each
        query site's element set contained in, or containing, the matched
        reference site's).  Compositions within the global tolerance give
        ``PRISTINE`` (the same disordered phase, occupancies refined
        differently); otherwise ``DISORDERED_FAMILY`` (e.g. Ta1/3Cr2/3O2 vs.
        Ta1/2Cr1/2O2).
        """
        query = drop_minor_species(query, self.min_occupancy)
        reference = drop_minor_species(reference, self.min_occupancy)
        vr = volume_ratio(query, reference)
        if not self.min_volume_ratio <= vr <= self.max_volume_ratio:
            return MatchResult(MatchKind.NONE, details={"reason": "volume per atom", "volume_ratio": round(vr, 3)})
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
        details = {
            "global_composition_ok": comp_ok,
            "reference_composition": reference.composition.formula,
            "volume_ratio": round(vr, 3),
        }
        if comp_ok:
            return MatchResult(MatchKind.PRISTINE, details=details)
        if minority_fraction(reference) < self.dilute_threshold:
            return MatchResult(MatchKind.NONE, details={"reason": "dilute reference", **details})
        return MatchResult(MatchKind.DISORDERED_FAMILY, details=details)

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
            if rs.species.num_atoms > 1 - self.vacancy_tol:
                return MatchResult(MatchKind.NONE, details={"reason": "fully occupied site left empty"})

        # (c) composition per sublattice.  Reference sites are pooled by the set
        # of elements they host (e.g. all {Ti, Cr} sites): the site mapping
        # StructureMatcher returns is only one of several equivalent ones, so
        # a partially ordered reference (two {Ti, Cr} orbits with different
        # ratios) would otherwise look like a mismatch depending on that
        # arbitrary choice.  The strict per-orbit deviation is reported too.
        per_sublattice = _occupancy_deviations(q, matched, mapped, _element_set_key)
        max_dev = max(dev for dev, _ in per_sublattice.values())
        max_dev_orbit = max(dev for dev, _ in _occupancy_deviations(q, matched, mapped, _site_key).values())
        # each sublattice: occupancies reproduced within tolerance, or realized by
        # its majority occupants (idealized ordering of minor anti-site disorder)
        sublattices_ok = all(dev <= self.occupancy_tol + 1e-8 or maj for dev, maj in per_sublattice.values())
        n_majority = sum(1 for dev, maj in per_sublattice.values() if dev > self.occupancy_tol + 1e-8 and maj)

        ref_comp = Composition({el: amt for el, amt in mapped.composition.items()})
        comp_ok = fractional_composition_close(q.composition, ref_comp, self.composition_atol, self.composition_rtol)
        details = {
            "max_sublattice_occupancy_deviation": round(max_dev, 4),
            "max_site_occupancy_deviation": round(max_dev_orbit, 4),
            "supercell_sites": len(mapped),
            "reference_supercell_composition": ref_comp.formula,
            "global_composition_ok": comp_ok,
            "idealized_sublattices": n_majority,
        }
        if comp_ok and sublattices_ok:
            return MatchResult(MatchKind.DISORDERED_REALIZATION, details=details)
        if max_dev > self.family_max_deviation + 1e-8:
            return MatchResult(MatchKind.NONE, details={"reason": "occupancies too far for a family", **details})
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
        # read species from the original (ordered) structure: sites of an
        # earlier clique are already mixed in ``s``
        for i, site in enumerate(structure):
            if site.specie in clique:
                s.replace(i, mix, coords=site.frac_coords)
    return s


@dataclass
class SiblingOrderingMatcher:
    """Different ordering of chemically similar elements, same formula.

    MatterGen merges elements by Hume-Rothery-like rules (|dX| <= 1.0,
    relative radius <= 0.3, transitive cliques), which in practice merges
    cations with anions or every element of a compound (the audit of the
    figshare benchmark found 239 of 374 such hits comparing bare atomic
    packings).  Here the stricter :func:`similarity_groups` are used (ionic
    rule for compounds with anions, Pettifor rule otherwise); structures in
    which every element would be merged, or with anions but no oxidation-state
    guess, are skipped.
    """

    ltol: float = 0.2
    stol: float = 0.3
    angle_tol: float = 5.0
    max_dm: float = 6
    max_dr: float = 0.2
    max_ionic_dr: float = 0.15

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

    def groups(self, composition: Composition) -> list[list[Element]]:
        """Multi-element similarity groups usable for this composition ([] if none)."""
        ionic = guess_ionic_radii(composition)
        if not ionic and has_anion(composition):
            return []
        groups = similarity_groups(composition.elements, self.max_dm, self.max_dr, ionic, self.max_ionic_dr)
        if len(groups) == 1:
            return []  # everything merged: only the bare packing would be compared
        return [g for g in groups if len(g) > 1]

    def prepare(self, query: Structure) -> Structure | None:
        groups = self.groups(query.composition)
        if not groups:
            return None
        comp = query.composition
        s = query.copy()
        for group in groups:
            total = sum(comp[el] for el in group)
            mix = {el: comp[el] / total for el in group}
            for i, site in enumerate(query):
                if site.specie in group:
                    s.replace(i, mix, coords=site.frac_coords)
        return s

    def fit(self, disordered_query: Structure | None, query: Structure, reference: Structure) -> bool:
        if disordered_query is None or query.composition.reduced_formula != reference.composition.reduced_formula:
            return False
        try:
            return bool(self._sm.fit(disordered_query, reference))
        except Exception as exc:  # pragma: no cover
            logger.debug("SiblingOrderingMatcher failed: %s", exc)
            return False
