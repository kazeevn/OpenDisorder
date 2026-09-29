"""Compare opendisorder's matcher with MatterGen's ``DisorderedStructureMatcher``.

MatterGen (Zeni et al., Nature 2025; ``mattergen/evaluation/utils/dataset_matcher.py``)
groups reference entries by reduced formula (ordered) or chemical system
(disordered) and calls ``DisorderedStructureMatcher.fit``; we mirror that with
the same COD candidates opendisorder uses (same formula + same chemical
system, disordered only).

Parts:

1. synthetic cases with known answers (as in ``tests/test_matching.py``);
2. the literature cases of ``benchmarks/literature`` (needs an index);
3. optionally, a sample of a figshare results file (needs an index).

Usage::

    python compare.py --db <index.sqlite> [--figshare-results results/mp_20__MatterGen__MatterGen_10k.csv.gz \\
        --figshare-set mp_20/MatterGen/MatterGen_10k --sample 1000]
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import pandas as pd
from pymatgen.core import Lattice, Structure

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mattergen_structure_matcher import DefaultDisorderedStructureMatcher  # noqa: E402

from opendisorder import CheckOptions, CODIndex, NoveltyChecker  # noqa: E402
from opendisorder.matching import DisorderMatcher  # noqa: E402

warnings.filterwarnings("ignore")

CAT = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
AN = [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]


def rocksalt_cell(cations, anions=("O",) * 4, abc=(4.25, 4.25, 4.1)):
    return Structure(Lattice.orthorhombic(*abc), list(cations) + list(anions), CAT[: len(cations)] + AN[: len(anions)])


def disordered_rocksalt(occ, a=4.2):
    return Structure.from_spacegroup("Fm-3m", Lattice.cubic(a), [occ, "O"], [[0, 0, 0], [0.5, 0, 0]])


def synthetic_cases():
    lifeo = disordered_rocksalt({"Li": 0.5, "Fe": 0.5})
    li_vac = disordered_rocksalt({"Li": 0.5})
    zb = Structure.from_spacegroup("F-43m", Lattice.cubic(4.2), [{"Li": 0.5, "Fe": 0.5}, "O"], [[0, 0, 0], [0.25] * 3])
    return [
        # (label, ordered, disordered, should be an ordered realization?)
        ("LiFeO2 ordering of (Li,Fe)O", rocksalt_cell(["Li", "Li", "Fe", "Fe"]), lifeo, True),
        ("ordered vacancies LiO2 of Li0.5O", rocksalt_cell(["Li", "Li"]), li_vac, True),
        ("Li3FeO4 (off-stoichiometry)", rocksalt_cell(["Li", "Li", "Li", "Fe"]), lifeo, False),
        ("vacancies filled: LiO", rocksalt_cell(["Li"] * 4), li_vac, False),
        ("missing O sublattice: LiFe", rocksalt_cell(["Li", "Li", "Fe", "Fe"], anions=()), lifeo, False),
        ("half-empty cation lattice: LiO2 vs (Li,Fe)O", rocksalt_cell(["Li", "Li"]), lifeo, False),
        ("foreign element Co", rocksalt_cell(["Li", "Li", "Co", "Co"]), lifeo, False),
        ("different framework (zinc blende)", rocksalt_cell(["Li", "Li", "Fe", "Fe"]), zb, False),
    ]


def part_synthetic() -> pd.DataFrame:
    mg, ours = DefaultDisorderedStructureMatcher(), DisorderMatcher()
    rows = []
    for label, ordered, disordered, truth in synthetic_cases():
        rows.append(
            {
                "case": label,
                "expected": truth,
                "MatterGen": bool(mg.fit(ordered, disordered)),
                "opendisorder": ours.match(ordered, disordered).kind.value,
            }
        )
    return pd.DataFrame(rows)


def mattergen_matches(index: CODIndex, checker: NoveltyChecker, s: Structure) -> list[str]:
    """COD ids matched by MatterGen's matcher among formula + chemsys candidates."""
    mg = DefaultDisorderedStructureMatcher()
    cand = checker.candidate_ids(s)
    refs = index.get_structures(list(dict.fromkeys(cand.pristine + cand.disorder)), fetch=False)
    hits = []
    for cid, ref in refs.items():
        same_formula = cid in cand.pristine and ref.composition.reduced_formula == s.composition.reduced_formula
        if not (same_formula or not ref.is_ordered):
            continue
        try:
            if mg.fit(s, ref):
                hits.append(cid)
        except Exception:
            pass
    return hits


def part_literature(index: CODIndex) -> pd.DataFrame:
    lit = HERE.parent / "literature"
    results = pd.read_csv(lit / "results.csv", dtype=str)
    checker = NoveltyChecker(index, options=CheckOptions(check_substitutions=False))
    rows = []
    for _, r in results.iterrows():
        s = Structure.from_file(lit / "cifs" / r["file"])
        hits = mattergen_matches(index, checker, s)
        rows.append(
            {
                "paper": r["paper"],
                "formula": r["formula"],
                "opendisorder": r["kind"],
                "MatterGen": "match" if hits else "none",
                "MatterGen COD": ";".join(hits[:3]),
            }
        )
    return pd.DataFrame(rows)


def part_figshare(index: CODIndex, results_file: Path, set_name: str, sample: int) -> pd.DataFrame:
    sys.path.insert(0, str(HERE.parent / "figshare"))
    from load_generated import iter_structures

    res = pd.read_csv(results_file, dtype={"sample_id": str}).set_index("sample_id")
    structures = {str(sid): s for _, sid, s in iter_structures(sets=[set_name])}
    ids = res.sample(min(sample, len(res)), random_state=0).index
    checker = NoveltyChecker(index, options=CheckOptions(check_substitutions=False))
    rows = []
    for sid in ids:
        hits = mattergen_matches(index, checker, structures[sid])
        rows.append({"sample_id": sid, "opendisorder": res.loc[sid, "kind"], "MatterGen": "match" if hits else "none"})
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path)
    ap.add_argument("--figshare-results", type=Path)
    ap.add_argument("--figshare-set")
    ap.add_argument("--sample", type=int, default=1000)
    args = ap.parse_args(argv)

    syn = part_synthetic()
    print("## Synthetic cases\n")
    print(syn.to_markdown(index=False))
    out = ["## Synthetic cases", "", syn.to_markdown(index=False), ""]
    if args.db:
        with CODIndex(args.db) as index:
            lit = part_literature(index)
            tab = pd.crosstab(lit["opendisorder"], lit["MatterGen"])
            print("\n## Literature cases (rows: opendisorder, columns: MatterGen)\n")
            print(tab.to_markdown())
            lit.to_csv(HERE / "literature_comparison.csv", index=False)
            out += ["## Literature cases", "", tab.to_markdown(), ""]
            if args.figshare_results:
                fig = part_figshare(index, args.figshare_results, args.figshare_set, args.sample)
                tab = pd.crosstab(fig["opendisorder"], fig["MatterGen"])
                print(f"\n## {args.figshare_set}: random sample of {len(fig)}\n")
                print(tab.to_markdown())
                fig.to_csv(HERE / "figshare_comparison.csv", index=False)
                out += [f"## {args.figshare_set} (random sample of {len(fig)})", "", tab.to_markdown(), ""]
    (HERE / "results.md").write_text("\n".join(out))


if __name__ == "__main__":
    main()
