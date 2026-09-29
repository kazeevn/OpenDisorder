import warnings

import pytest
from pymatgen.core import Lattice, Structure

from opendisorder.matching import (
    DisorderMatcher,
    MatchKind,
    PristineMatcher,
    SiblingOrderingMatcher,
    substitution_cliques,
)

warnings.filterwarnings("ignore")

CAT = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
AN = [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]


def rocksalt_cell(cations, anions=("O",) * 4, abc=(4.25, 4.25, 4.1)):
    species = list(cations) + list(anions)
    coords = CAT[: len(cations)] + AN[: len(anions)]
    return Structure(Lattice.orthorhombic(*abc), species, coords)


def disordered_rocksalt(occ, a=4.2, anion="O"):
    return Structure.from_spacegroup("Fm-3m", Lattice.cubic(a), [occ, anion], [[0, 0, 0], [0.5, 0, 0]])


@pytest.fixture(scope="module")
def dm():
    return DisorderMatcher()


def test_pristine_same_structure_rescaled():
    s = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.64), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    s2 = s.copy()
    s2.scale_lattice(s.volume * 1.1)
    s2.make_supercell([1, 1, 2])
    assert PristineMatcher().fit(s, s2)


def test_pristine_different_structure():
    rs = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.64), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    cscl = Structure.from_spacegroup("Pm-3m", Lattice.cubic(3.3), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    assert not PristineMatcher().fit(rs, cscl)


def test_ordered_realization(dm):
    res = dm.match(rocksalt_cell(["Li", "Li", "Fe", "Fe"]), disordered_rocksalt({"Li": 0.5, "Fe": 0.5}))
    assert res.kind is MatchKind.DISORDERED_REALIZATION


def test_ordered_vacancies(dm):
    res = dm.match(rocksalt_cell(["Li", "Li"]), disordered_rocksalt({"Li": 0.5}))
    assert res.kind is MatchKind.DISORDERED_REALIZATION


def test_off_stoichiometry_is_family_only(dm):
    res = dm.match(rocksalt_cell(["Li", "Li", "Li", "Fe"]), disordered_rocksalt({"Li": 0.5, "Fe": 0.5}))
    assert res.kind is MatchKind.DISORDERED_FAMILY


def test_filled_vacancies_is_family_only(dm):
    res = dm.match(rocksalt_cell(["Li"] * 4), disordered_rocksalt({"Li": 0.5}))
    assert res.kind is MatchKind.DISORDERED_FAMILY


def test_missing_sublattice_rejected(dm):
    res = dm.match(rocksalt_cell(["Li", "Li", "Fe", "Fe"], anions=()), disordered_rocksalt({"Li": 0.5, "Fe": 0.5}))
    assert res.kind is MatchKind.NONE


def test_empty_fully_occupied_site_rejected(dm):
    res = dm.match(rocksalt_cell(["Li", "Li"]), disordered_rocksalt({"Li": 0.5, "Fe": 0.5}))
    assert res.kind is MatchKind.NONE


def test_foreign_element_rejected(dm):
    res = dm.match(rocksalt_cell(["Li", "Li", "Co", "Co"]), disordered_rocksalt({"Li": 0.5, "Fe": 0.5}))
    assert res.kind is MatchKind.NONE


def test_different_framework_rejected(dm):
    zb = Structure.from_spacegroup(
        "F-43m", Lattice.cubic(4.2), [{"Li": 0.5, "Fe": 0.5}, "O"], [[0, 0, 0], [0.25, 0.25, 0.25]]
    )
    res = dm.match(rocksalt_cell(["Li", "Li", "Fe", "Fe"]), zb)
    assert res.kind is MatchKind.NONE


def test_substitution_cliques():
    from pymatgen.core import Element

    cliques = substitution_cliques([Element("Na"), Element("K"), Element("Cl")])
    assert {Element("Na"), Element("K")} in cliques
    assert all(Element("Cl") not in c for c in cliques)


def _nacl_type_ordering(rule, a=6.0):
    """2x2x2 conventional rocksalt cell (in units of a/2) with Na/K on cations chosen by ``rule``."""
    species, coords = [], []
    for n1 in range(4):
        for n2 in range(4):
            for n3 in range(4):
                frac = [n1 / 4, n2 / 4, n3 / 4]
                if (n1 + n2 + n3) % 2 == 0:
                    species.append("Na" if rule(n1, n2, n3) else "K")
                else:
                    species.append("Cl")
                coords.append(frac)
    return Structure(Lattice.cubic(2 * a), species, coords)


def test_sibling_ordering():
    l10 = _nacl_type_ordering(lambda n1, n2, n3: n3 % 2 == 0)
    l11 = _nacl_type_ordering(lambda n1, n2, n3: ((n1 + n2 + n3) // 2) % 2 == 0)
    assert l10.composition == l11.composition
    assert not PristineMatcher().fit(l10, l11)
    sm = SiblingOrderingMatcher()
    assert sm.fit(sm.prepare(l10), l10, l11)
