# opendisorder

Is this crystal structure already in the [Crystallography Open Database](https://www.crystallography.net/)?
It is, if it matches a COD entry directly, **or** if it is merely an ordered realization of a
COD entry with partial occupancies.

`opendisorder` answers with a graded verdict:

| kind | meaning | counts as known? |
|---|---|---|
| `pristine` | same structure as an ordered COD entry (same reduced formula, `StructureMatcher` fit) | yes |
| `disordered_realization` | ordered realization ("coloring", ordered approximant) of a disordered COD entry with a consistent composition | yes |
| `disordered_family` | same framework and site-species pattern as a disordered COD entry, but a different stoichiometry (another member of a solid-solution series) | no, but suspicious |
| `substitutional_ordering` | ordering of chemically similar elements over the sites of a known COD structure with fewer elements, e.g. GNoME's TbHo₂Tm₉Co₄ as an ordered (Tb,Ho,Tm)₃Co on the Ho₃Co structure | no, but suspicious |
| `sibling_ordering` | a different ordering of chemically similar elements on the same parent structure as an ordered COD entry with the same formula | no, but suspicious |
| `none` | no match | no |

`NoveltyReport.is_novel` is `False` for `pristine` and `disordered_realization`.
`NoveltyReport.is_novel_strict` is `True` only for `none`. The CLI prints `known`, `suspicious`
or `novel`.

A disordered query, such as a Rietveld refinement with mixed sites, is compared with disordered
COD entries of the same chemical system. If the framework and site-species pattern match, the
result is `pristine` when the composition agrees within tolerance and `disordered_family`
otherwise.

## What counts as a disorder match

Recent critiques of generative and high-throughput materials discovery pointed out that many
"new" predicted compounds are ordered versions of known compositionally disordered phases:

- Leeman et al., *PRX Energy* 3, 011002 (2024), on the A-Lab
- Cheetham & Seshadri, *Chem. Mater.* 36, 3490 (2024), on GNoME
- Juelsholt, *Mater. Horiz.* (2026), on MatterGen's TaCr₂O₆

The only operational definition we found is MatterGen's (Zeni et al., *Nature* 2025, SI D.3). An
ordered structure matches a disordered reference if two things hold:

- it is a supercell of the reference framework and each of its sites carries one of the elements
  of the corresponding reference site;
- the compositions agree within `|x₁ − x₂| ≤ 0.01 + 0.1·(x₁ + x₂)/2`.

MatterGen tests this with `StructureMatcher(attempt_supercell=True, allow_subset=True,
comparator=OrderDisorderElementComparator())`.

That screen on its own gives false positives. `allow_subset` lets a whole sublattice go missing,
and the global composition test is diluted by the anions. For example, Li₃FeO₄ passes against
(Li₀.₅Fe₀.₅)O. So `opendisorder` puts the screen (with `supercell_size="volume"`, needed for
ordered vacancies) before an explicit check of the site mapping, in
`opendisorder.matching.DisorderMatcher`:

1. every query atom sits on a reference site whose species include that element;
2. reference sites left empty by the query are partially occupied (they have vacancies);
3. for every class of equivalent reference sites and every element, the query's site fraction is
   within `occupancy_tol = 0.1` of the reference occupancy, and MatterGen's global composition
   test passes as well.

If checks 1 and 2 pass but check 3 fails, the result is `disordered_family`. Juelsholt's TaCr₂O₆
against (Ta₀.₅Cr₀.₅)O₂ falls in this class.

MatterGen's ordered-vs-ordered route becomes `sibling_ordering`. Elements that are similar enough
by Hume-Rothery-like rules (|ΔX| ≤ 1.0, relative atomic radius difference ≤ 0.3, transitive
cliques) are merged into mixed sites, and the structure is compared again.

### Substitutions of similar elements

Many disputed predictions are orderings of chemically similar elements on a known structure that
has *fewer* elements. Examples:

- GNoME's TbHo₂Tm₉Co₄ as an ordered (Tb,Ho,Tm)₃Co (Cheetham & Seshadri 2024);
- the A-Lab target Ba₂ZrSnO₆ as BaZr₀.₅Sn₀.₅O₃ (Leeman et al. 2024).

Their parent phases sit in a different chemical system, so the layers above cannot find them. For
these, `opendisorder` does the following:

1. It groups the query's elements into similarity groups
   (`opendisorder.matching.similarity_groups`):
   - in ionic compounds (oxidation states from pymatgen's `oxi_state_guesses`), elements are
     similar if they are isovalent and their Shannon radii are within 15% (Goldschmidt), e.g.
     Zr⁴⁺/Sn⁴⁺ and Mg²⁺/Ni²⁺, or if they are chemically close, e.g. K⁺/Na⁺;
   - otherwise elements are similar if their Pettifor (Mendeleev) numbers differ by at most 6 and
     their atomic radii by at most 20%, e.g. the rare earths.
   Groups are built by complete linkage, so there is no chaining. MatterGen's Hume-Rothery
   thresholds are far too loose for this step: they would put Co in the same group as the rare
   earths.
2. It looks up COD entries whose element set is a strict subset of the query's, keeping all
   non-grouped elements and at least one element of each group.
3. It maps every group onto one element in both structures and compares them. A pristine match, or
   an ordered realization when the reference is disordered, gives `substitutional_ordering`.

## How COD is indexed

- `opendisorder build-index` downloads COD's CSV metadata dump (all ~530k entries, ~400 MB, about
  a minute) and stores it in SQLite. Both `formula` and `calcformula` are indexed by
  - reduced formula, which supplies layer 1 candidates;
  - chemical system, which supplies layer 2 candidates, since an ordered realization usually has a
    different stoichiometry from the CIF formula.
- CIFs are bulk-loaded from COD's CIF archive (`opendisorder import-archive`). The ~17 GB
  `cod-cifs-mysql.txz` is streamed and decompressed on the fly, not stored, and the CIFs are kept
  zlib-compressed in the database (about 5 GB for all of COD). `--max-cell-atoms` keeps only a
  part of COD when disk is short.
- Entries missing from the archive, such as those added after its snapshot, are downloaded
  lazily through COD's REST API in zip batches of up to 300. A local mirror of COD CIFs can also
  be used (`cif_dir=`).
- CIFs are parsed with pymatgen (`occupancy_tolerance=1.05`), and the parsed structures are cached
  in the same database.
- Layer 2 references are, by default, entries that COD flags `has disorder` or whose formula is
  fractional, plus same-formula entries. `--all-chemsys` considers every entry of the chemical
  system instead.
- COD entries whose conventional cell has more than 4.4× the query's atom count are skipped: they
  cannot contain an ordered realization.

## Usage

```bash
pip install -e .
opendisorder build-index                    # metadata, ~2 min; index in ~/.cache/opendisorder/cod.sqlite
opendisorder import-archive                 # optional bulk load of all CIFs (streams ~17 GB)
opendisorder check my.cif POSCAR --json     # downloads any COD CIFs still missing
```

```python
from pymatgen.core import Structure
from opendisorder import CODClient, CODIndex, NoveltyChecker

with CODIndex("~/.cache/opendisorder/cod.sqlite", client=CODClient()) as index:
    checker = NoveltyChecker(index)
    report = checker.check(Structure.from_file("my.cif"))
    print(report.kind, report.is_novel, report.cod_ids)
    reports = checker.check_many(list_of_structures, processes=4)  # grouped by chemical system
```

The building blocks can be used on their own:

- `PristineMatcher().fit(a, b)`
- `DisorderMatcher().match(ordered, disordered)`, which returns a `MatchResult` with a grade and
  details
- `SiblingOrderingMatcher`

## Case study: MatterGen's TaCr₂O₆

MatterGen (Zeni et al., *Nature* 2025) generated and synthesized "new" TaCr₂O₆. Juelsholt
(*Mater. Horiz.* 13, 5672, 2026) argued that it is the rutile-type solid solution Ta₁/₂Cr₁/₂O₂,
known since 1971 (COD 9011856). Running `opendisorder check --exhaustive` on the structures
involved gives:

| structure | result | COD match |
|---|---|---|
| MatterGen TaCr₂O₆, generated (`data-release/cifs.zip`) | `disordered_family` | 9011856 (Ta₀.₅Cr₀.₅O₂, P4₂/mnm) |
| MatterGen's Rietveld refinement, Ta₁/₃Cr₂/₃O₂ (disordered) | `disordered_family` | 2002668, 9011856 |
| TaCrO₄ in MatterGen's training data: mp-753467, mp-756340, Alexandria agm003743920 | `disordered_realization` | 9011856, 2002578 |
| mp-21512 Ta₂CrO₆ | `pristine` | 1001091 |
| MatterGen Mn₅Fe₃O₈, which Juelsholt calls an ordered rock salt | `disordered_family` | 4344257 ((Mn₀.₅₆Fe₀.₄₄)O, Fm-3m) |
| MatterGen Cr₂MoO₆ | `none` | not in COD |

TaCr₂O₆ is reported as `disordered_family` rather than `disordered_realization`. It has the
framework and the mixed Ta/Cr site of Ta₁/₂Cr₁/₂O₂, but the occupancies differ by 1/6. By
Juelsholt's argument that still makes it not novel, and `NoveltyReport.is_novel_strict` follows
that reading. `tests/data` holds these files as a regression test.

## Case studies: GNoME and the A-Lab

`benchmarks/literature` screens 68 structures that the literature disputes as new:

- the GNoME predictions in Cheetham & Seshadri (2024);
- all A-Lab targets from Leeman et al. (2024), taken from Materials Project;
- the MatterGen cases above.

The per-compound table is in [`benchmarks/literature/README.md`](benchmarks/literature/README.md).

| paper | known | suspicious | no COD match |
|---|---|---|---|
| Juelsholt 2026 (MatterGen) | 4 | 3 | 0 |
| Leeman et al. 2024 (A-Lab) | 6 | 13 | 24 |
| Cheetham & Seshadri 2024 (GNoME) | 0 | 2 | 15 |

What the flags cover and what they miss:

- **A-Lab.** The flags cover most of the disorder cases Leeman et al. name:
  - (Ni,Mg)TiO₃ ilmenite, the pyroxene, (Cu,Mg)₂V₂O₇ and the Zn–Mn–Ni spinel, as ordered
    realizations;
  - the Mg–Ni rock salt, BaZr₀.₅Sn₀.₅O₃, the Zn(Cr,Fe)₂O₄ and Zn–Ni–Sb spinels, the M₂(PO₃)₄
    phosphates and the (Na,K)Pb apatite, as substitutional orderings;
  - the Na-doped KTP phases and the alluaudites, as the same disordered family.
- **A-Lab misses.** Most are cases where the known phase is only in ICSD. Examples are the
  Pb₂(M,Sb)₂O₆.₅ pyrochlores, which have no Pb–Sb–M–O entries in COD.
- **GNoME.** Cheetham & Seshadri's Table 2 compares GNoME compounds with *isostructural analogues
  containing different elements*, such as Pb₂Pd₃SnY₃ against PbPtGd. That is structure-type
  novelty, not disorder, and is out of scope here. Several rare-earth intermetallic parents, such
  as Ln₅Ir₃ and (Dy,Ho)Rh₃, are not in COD.

## Benchmark: generated crystals

`benchmarks/figshare` runs the generated structures of figshare
[29145101](https://figshare.com/articles/dataset/29145101) through the checker. These come from
WyFormer, DiffCSP, DiffCSP++, WyCryst, SymmCD, CrystalFormer, MiAD, MatterGen, Crys-JEPA and
others. See [`benchmarks/figshare/README.md`](benchmarks/figshare/README.md) for results.

## Caveats

- COD is mostly experimental and mostly organic/metal-organic. A "novel" verdict here means novel
  with respect to COD, not ICSD or Materials Project.
- The default tolerances are StructureMatcher's (`ltol=0.2, stol=0.3, angle_tol=5`). Strongly
  distorted relaxed orderings (Jahn–Teller distortions, octahedral tilts) may fall outside them.
- COD's disorder flags and the CIF-derived occupancies are taken at face value. X-ray structures
  often lack H atoms, and split positions (positional disorder) are treated as partially occupied
  sites.

## Tests

```bash
pip install -e '.[test]' && pytest
```
