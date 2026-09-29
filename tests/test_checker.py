import csv
import warnings

import pytest
from pymatgen.core import Lattice, Structure

from opendisorder import CheckOptions, CODIndex, MatchKind, NoveltyChecker

warnings.filterwarnings("ignore")

FIELDS = ["file", "formula", "calcformula", "sgNumber", "vol", "flags", "method", "svnrevision"]


def _cif(structure):
    from pymatgen.io.cif import CifWriter

    return str(CifWriter(structure))


@pytest.fixture()
def index(tmp_path):
    nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.64), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    lifeo2 = Structure.from_spacegroup(
        "Fm-3m", Lattice.cubic(4.16), [{"Li": 0.5, "Fe": 0.5}, "O"], [[0, 0, 0], [0.5, 0, 0]]
    )
    rows = [
        {"file": "1000001", "formula": "- Cl Na -", "calcformula": "- Cl Na -", "flags": "has coordinates"},
        {"file": "1000002", "formula": "- Fe0.5 Li0.5 O -", "calcformula": "- Fe0.5 Li0.5 O -",
         "flags": "has coordinates,has disorder"},
        {"file": "1000003", "formula": "- C H4 -", "calcformula": "- C H4 -", "flags": ""},
    ]
    csv_path = tmp_path / "meta.csv"
    with open(csv_path, "w", newline="") as fh:
        fh.write("# comment line\n")
        w = csv.DictWriter(fh, FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    idx = CODIndex(tmp_path / "cod.sqlite")
    assert idx.build_from_metadata(csv_path) == 2
    idx.add_cif("1000001", _cif(nacl))
    idx.add_cif("1000002", _cif(lifeo2))
    yield idx
    idx.close()


def test_index_lookup(index):
    assert index.ids_by_reduced_formula("NaCl") == ["1000001"]
    assert index.ids_by_chemsys("Fe-Li-O") == ["1000002"]
    assert index.ids_by_chemsys("Fe-Li-O", disorder_only=True) == ["1000002"]
    assert not index.get_structure("1000002").is_ordered


def test_checker_layers(index):
    checker = NoveltyChecker(index, options=CheckOptions(time_limit=None))
    nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.5), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    rep = checker.check(nacl)
    assert rep.kind is MatchKind.PRISTINE and rep.cod_ids == ["1000001"] and not rep.is_novel

    cat = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    an = [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
    ordered = Structure(Lattice.orthorhombic(4.2, 4.2, 4.1), ["Li", "Li", "Fe", "Fe"] + ["O"] * 4, cat + an)
    rep = checker.check(ordered)
    assert rep.kind is MatchKind.DISORDERED_REALIZATION and rep.cod_ids == ["1000002"]

    novel = Structure.from_spacegroup("Pm-3m", Lattice.cubic(4.0), ["Cs", "Xe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    rep = checker.check(novel)
    assert rep.kind is MatchKind.NONE and rep.is_novel

    reports = checker.check_many([nacl, ordered, novel], processes=1)
    assert [r.kind for r in reports] == [MatchKind.PRISTINE, MatchKind.DISORDERED_REALIZATION, MatchKind.NONE]


def test_check_many_parallel(index):
    checker = NoveltyChecker(index)
    nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.5), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    reports = checker.check_many([nacl, nacl.copy()], processes=2)
    assert all(r.kind is MatchKind.PRISTINE for r in reports)


def test_import_archive(index, tmp_path):
    import io
    import tarfile

    cif = index.get_cif("1000001")
    archive = tmp_path / "cod.tgz"
    with tarfile.open(archive, "w:gz") as tf:
        for name in ("cif/1/00/00/1000001.cif", "cif/1/00/00/1000009.cif", "mysql/data.sql"):
            data = cif.encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    index.db.execute("DELETE FROM cifs")
    # 1000009 is not in the metadata index and the SQL dump is not a CIF
    assert index.import_archive(archive) == 1
    assert index.get_cif("1000001") == cif
    assert index.get_cif("1000009") is None


def test_strip_cif_bulk():
    from opendisorder.index import strip_cif_bulk

    cif = "data_x\n_cell_length_a 5.0\n_shelx_hkl_file\n;\n 1 0 0 12.3 0.4\n 0 1 0 1.2 0.1\n;\n_cell_length_b 5.0\n"
    stripped = strip_cif_bulk(cif)
    assert "_shelx_hkl_file" not in stripped and "12.3" not in stripped
    assert "_cell_length_a 5.0" in stripped and "_cell_length_b 5.0" in stripped
    assert strip_cif_bulk("data_y\n_cell_length_a 1\n") == "data_y\n_cell_length_a 1\n"


def test_substitutional_ordering(index):
    """Ordered (Na,K)Cl on the rock-salt lattice is an ordering of similar
    elements over the known NaCl structure."""
    cat = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    an = [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]
    q = Structure(Lattice.orthorhombic(6.0, 6.0, 6.2), ["Na", "Na", "K", "K"] + ["Cl"] * 4, cat + an)
    checker = NoveltyChecker(index, options=CheckOptions(time_limit=None))
    rep = checker.check(q)
    assert rep.kind is MatchKind.SUBSTITUTIONAL_ORDERING and rep.cod_ids == ["1000001"]
    assert rep.is_novel and not rep.is_novel_strict
    rep = NoveltyChecker(index, options=CheckOptions(check_substitutions=False)).check(q)
    assert rep.kind is MatchKind.NONE


def test_check_many_survives_worker_crash(index, monkeypatch):
    """A worker that dies (e.g. OOM-killed) only costs the offending structure."""
    import os

    import opendisorder.checker as checker_mod

    original = checker_mod._check_group_job

    def crashy(job):
        if any(s.composition.reduced_formula == "CsXe" for _, s in job[2]):
            os._exit(9)
        return original(job)

    monkeypatch.setattr(checker_mod, "_check_group_job", crashy)
    nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.5), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    csxe = Structure.from_spacegroup("Pm-3m", Lattice.cubic(4.0), ["Cs", "Xe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    reports = NoveltyChecker(index).check_many([nacl, csxe, nacl.copy()], processes=2)
    assert reports[0].kind is MatchKind.PRISTINE and reports[2].kind is MatchKind.PRISTINE
    assert reports[1].kind is MatchKind.NONE and "worker process died" in reports[1].errors[0]
