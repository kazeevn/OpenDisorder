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


def test_disordered_vs_disordered(dm):
    rutile = dict(sg="P4_2/mnm", coords=[[0, 0, 0], [0.3, 0.3, 0]])

    def rut(occ, c=3.0):
        return Structure.from_spacegroup(rutile["sg"], Lattice.tetragonal(4.64, c), [occ, "O"], rutile["coords"])

    half = rut({"Ta": 0.5, "Cr": 0.5})
    assert dm.match_disordered(rut({"Ta": 1 / 3, "Cr": 2 / 3}), half).kind is MatchKind.DISORDERED_FAMILY
    assert dm.match_disordered(rut({"Ta": 0.48, "Cr": 0.52}, c=3.05), half).kind is MatchKind.PRISTINE
    rocksalt = disordered_rocksalt({"Ta": 0.5, "Cr": 0.5})
    assert dm.match_disordered(rocksalt, half).kind is MatchKind.NONE


def test_mattergen_tacr2o6_regression():
    """MatterGen's TaCr2O6 is Ta1/3Cr2/3O2, same disordered rutile as Ta1/2Cr1/2O2
    (COD 9011856); ordered TaCrO4 from Materials Project is a realization of it."""
    from pathlib import Path

    data = Path(__file__).parent / "data"
    cod = Structure.from_file(data / "cod_9011856.cif")
    assert not cod.is_ordered
    dm = DisorderMatcher()
    assert dm.match(Structure.from_file(data / "mattergen_TaCr2O6.cif"), cod).kind is MatchKind.DISORDERED_FAMILY
    assert dm.match(Structure.from_file(data / "mp-753467_TaCrO4.cif"), cod).kind is MatchKind.DISORDERED_REALIZATION


def test_similarity_groups():
    from pymatgen.core import Composition, Element

    from opendisorder.matching import guess_ionic_radii, similarity_groups

    def groups(formula):
        comp = Composition(formula)
        return [{e.symbol for e in g} for g in similarity_groups(comp.elements, ionic=guess_ionic_radii(comp))]

    assert {"Tb", "Ho", "Tm"} in groups("TbHo2Tm9Co4")  # rare earths, not Co
    assert {"Co"} in groups("TbHo2Tm9Co4")
    assert {"Zr", "Sn"} in groups("Ba2ZrSnO6")  # isovalent, similar ionic radii
    assert {"Mg", "Ni"} in groups("Mg3NiO4")
    assert {"Na"} in groups("Na3Ca18Fe(PO4)14")  # heterovalent Na+/Ca2+ stay apart
    assert all(len(g) == 1 for g in similarity_groups([Element("Li"), Element("O")]))


def test_disorder_similar_elements_two_cliques():
    from opendisorder.matching import disorder_similar_elements

    # Na/K and Cl/Br form two separate cliques
    s = rocksalt_cell(["Na", "Na", "K", "K"], anions=("Cl", "Cl", "Br", "Br"), abc=(6.0, 6.0, 6.2))
    d = disorder_similar_elements(s)
    assert d is not None and not d.is_ordered
    assert d.composition.almost_equals(s.composition)


def test_partially_ordered_reference(dm):
    """A reference with two mixed Li/Fe orbits (0.7/0.3 and 0.3/0.7): an ordered
    layering is a realization whichever orbit the matcher maps Li onto."""
    a = {"Li": 0.7, "Fe": 0.3}
    b = {"Li": 0.3, "Fe": 0.7}
    ref = rocksalt_cell([a, a, b, b], abc=(4.2, 4.2, 4.2))
    for order in (["Li", "Li", "Fe", "Fe"], ["Fe", "Fe", "Li", "Li"]):
        res = dm.match(rocksalt_cell(order), ref)
        assert res.kind is MatchKind.DISORDERED_REALIZATION, res.details


def test_majority_ordering_of_antisite_disorder():
    """Minor anti-site disorder: the idealized ordering is a realization even
    when the sublattice occupancies are off by more than the tolerance."""
    dm = DisorderMatcher(occupancy_tol=0.1)
    a = {"Li": 0.88, "Mg": 0.12}  # minor Mg on the Li sites; Mg sites stay pure
    ref = rocksalt_cell([a, a, "Mg", "Mg"], abc=(4.2, 4.2, 4.2))
    res = dm.match(rocksalt_cell(["Li", "Li", "Mg", "Mg"]), ref)
    assert res.details["max_sublattice_occupancy_deviation"] > dm.occupancy_tol
    assert res.kind is MatchKind.DISORDERED_REALIZATION and res.details["idealized_sublattices"] == 1
