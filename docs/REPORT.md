# opendisorder: checking crystal structures against the COD, disorder included

*Write-up of the work done on 2026-09-29. The code is in `src/opendisorder`; benchmarks and
their raw results are in `benchmarks/`.*

## Summary

- **What it is.** `opendisorder` is a library and CLI that checks whether a crystal structure
  is already in the [Crystallography Open Database](https://www.crystallography.net/) (COD).
  It returns a graded verdict that goes beyond "same structure": it also recognises an ordered
  structure that is merely one ordering of a known compositionally disordered phase.
- **Open data only.** The reference data is COD (CC0) alone. ICSD and other proprietary
  databases are excluded by design.
- **Literature cases.** On the structures the literature disputes as "new", the library flags
  all the MatterGen cases Juelsholt names, including TaCr₂O₆ against Ta₁/₂Cr₁/₂O₂. It also
  catches 19 of the 43 A-Lab targets, and GNoME's orderings of rare-earth intermetallics.
- **Against MatterGen's matcher.** On 72 literature structures, MatterGen's
  `DisorderedStructureMatcher` finds nothing that `opendisorder` misses. `opendisorder` adds
  ordered-vacancy matches, stoichiometry "families" and substitutional orderings that MatterGen
  cannot find, because it only searches the query's own chemical system.
- **Generated crystals.** On 99.6k generated crystals from 10 generative models (figshare
  29145101), 5–15% match COD. For most models about half of the matches are
  disorder-related rather than exact.
- **COD coverage is the main limit.** COD recognises only about half of the experimentally
  observed MP-20 materials, so "no COD match" is weak evidence of novelty.
- **MP-20.** About half of MP-20 is flagged `theoretical` by Materials Project (MP). A large
  part of those are MP's own ordered approximants of disordered experimental phases, which MP
  flags as theoretical without exception.

## 1. What counts as a match

The literature leaves two notions of "not novel" implicit:

- **Same structure.** Same composition, and a `StructureMatcher` fit.
- **Ordered realization of a known disordered phase.** An ordered structure can be one
  "coloring" of a known compositionally disordered phase. Leeman et al. (*PRX Energy* 2024)
  and Cheetham & Seshadri (*Chem. Mater.* 2024) make this criticism of A-Lab and GNoME;
  Juelsholt (*Mater. Horiz.* 2026) makes it of MatterGen.

The only operational definition we found is MatterGen's (Zeni et al., *Nature* 2025, SI D.4):

- `StructureMatcher(attempt_supercell, allow_subset, OrderDisorderElementComparator)`;
- a global composition tolerance `|x₁−x₂| ≤ 0.01 + 0.1·x̄`;
- Hume-Rothery-like merging of similar elements for ordered-vs-ordered comparisons.

`opendisorder` keeps that structure screen and adds an explicit check of the site mapping.
It reports graded results:

| grade | meaning | known? |
|---|---|---|
| `pristine` | same structure as a COD entry; retried at loose tolerances with a density guard | yes |
| `disordered_realization` | ordered realization of a disordered COD entry | yes |
| `disordered_family` | same framework and site-species pattern, different stoichiometry (TaCr₂O₆ vs Ta₁/₂Cr₁/₂O₂) | suspicious |
| `substitutional_ordering` | ordering of chemically similar elements over a known structure with fewer elements (GNoME TbHo₂Tm₉Co₄ = ordered (Tb,Ho,Tm)₃Co) | suspicious |
| `sibling_ordering` | different ordering of similar elements, same formula | suspicious |
| `none` | no COD evidence | — |

A `disordered_realization` must pass all of the following:

1. Every query atom sits on a reference site that hosts its element.
2. No reference site with occupancy above 0.9 is left empty.
3. The global composition agrees, using MatterGen's tolerance.
4. Every *sublattice* (all sites hosting the same set of elements) is realized. That means
   either the occupancies are reproduced within 0.2, or the sublattice is filled by its
   majority occupants, which is the idealized ordering of minor anti-site disorder.
5. The volume per atom (or per site) is within 0.75–1.35 of the reference.

Two further rules apply:

- occupancies below 0.05 are ignored;
- a `disordered_family` is not assigned against dilute references (all sites ≥ 90% one
  species), nor beyond a 0.5 sublattice deviation.

**Similarity groups.** Used for substitutions and siblings, and built by complete linkage, so
groups never chain:

- *Ionic compounds:* elements must be isovalent, with Shannon radii within 15% or chemically
  close (e.g. Zr⁴⁺/Sn⁴⁺, Mg²⁺/Ni²⁺, K⁺/Na⁺).
- *Intermetallics:* Pettifor-number difference ≤ 6 and atomic radii within 20%.
- *Anion compounds without a charge-balanced oxidation-state guess* are never grouped.

MatterGen's Hume-Rothery rule was dropped for grouping: in the audit it merged cations with
anions, or all elements of a compound.

## 2. Data: an open COD index

- **Metadata.** `opendisorder build-index` loads COD's CSV metadata dump (527k entries,
  ~1.5 min). It indexes every entry by reduced formula and by chemical system, from both
  `formula` and `calcformula`, and records cell atom counts and disorder flags.
- **CIFs.** `opendisorder import-archive` streams COD's 17 GB CIF archive without storing it,
  resuming interrupted downloads with HTTP Range requests. It strips embedded hkl/fcf/res
  data, which is about 80% of the bytes, and stores all 514k CIFs in about 6 GB of SQLite.
  Entries newer than the archive are fetched through COD's REST API.
- **Parsing.** CIFs are parsed lazily with pymatgen and cached. Modulated (superspace) CIFs
  are rejected.
- **Batch checking.** `check_many` groups queries by chemical system and runs each group in a
  forked worker with a memory cap. A crashed worker costs only the offending structure.
- **Sharing.** A read-only mode lets cluster jobs share one index.

## 3. Validation

### Literature cases (`benchmarks/literature`, 72 structures)

| paper | cases | known | suspicious | no COD evidence |
|---|---|---|---|---|
| Juelsholt 2026 (MatterGen TaCr₂O₆ and related) | 7 | 4 | 3 | 0 |
| Leeman et al. 2024 (all A-Lab targets) | 43 | 10 | 9 | 24 |
| Cheetham & Seshadri 2024 (GNoME) | 21 | 0 | 3 | 18 |
| Zeni et al. 2025 SI (Cr₂MoO₆) | 1 | 0 | 0 | 1 |

- **Juelsholt (MatterGen).** All three disputed MatterGen structures are flagged as
  `disordered_family`: TaCr₂O₆, its Rietveld refinement, and Mn₅Fe₃O₈ against (Mn,Fe)O rock
  salt. The TaCrO₄ entries in MatterGen's training data are `disordered_realization`s of COD
  9011856.
- **Leeman et al. (A-Lab).** Most disorder cases the authors name are caught: the Ni/Mg
  ilmenite, pyroxene, alluaudites, spinels, BaZr₀.₅Sn₀.₅O₃, the Mg–Ni rock salt and the
  M₂(PO₃)₄ phosphates.
- **Cheetham & Seshadri (GNoME).** Catches TbHo₂Tm₉Co₄, DyPSi₅Tb₉ and TbSmSeO₂.

**Misses** are phases that exist only in ICSD, such as the Pb₂(M,Sb)₂O₆.₅ pyrochlores and
several Ln–Ir/Rh intermetallics. The table-2 cases of Cheetham & Seshadri are isostructural
analogues with *different* elements; that is structure-type novelty, which is out of scope.

### Comparison with MatterGen's matcher (`benchmarks/mattergen_baseline`)

The comparison uses MatterGen's `DisorderedStructureMatcher`, vendored under its MIT licence,
on the same candidates. Rows are `opendisorder` grades; columns are MatterGen's verdict.

| opendisorder | MatterGen match | MatterGen none |
|---|---|---|
| pristine / disordered_realization | 14 | 0 |
| disordered_family | 0 | 3 |
| substitutional_ordering | 0 | 12 |
| none | 0 | 43 |

On synthetic cases MatterGen misses ordered-vacancy superstructures (LiO₂ as an ordering of
Li₀.₅O). Both reject off-stoichiometric and sublattice-missing structures.

### Audit of match precision (110 sampled hits, before the precision fixes)

| grade | true | debatable | false |
|---|---|---|---|
| pristine | 19 | 1 | 0 |
| disordered_realization | 17 | 2 | 1 |
| substitutional_ordering | 13 | 6 | 1 |
| disordered_family | 8 | 5 | 7 |
| sibling_ordering | 7 | 11 | 2 |

Every failure mode found was fixed:

- superspace CIFs;
- dilute references;
- the loose Hume-Rothery grouping;
- Pettifor fallback for anion compounds without an oxidation-state guess;
- missing density check.

The fixes cut `sibling_ordering` hits from 374 to 9 and `disordered_family` hits by about 25%
in the figshare run.

### MP ↔ COD independent check

MP's task tags record which entries MP generated by ordering disordered ICSD structures. On the
MP-20 test split, `opendisorder` grades those 1,385 entries as:

- 0.5% `pristine`;
- 26.4% `disordered_realization`;
- 19.4% `substitutional_ordering`;
- 48% none, which reflects COD coverage.

The 1,003 entries built from *ordered* ICSD structures come out 15.2% `pristine` and only 4.4%
in the disorder grades.

## 4. Generated crystals (`benchmarks/figshare`)

Values are % of structures; the run uses code at HEAD, including the precision fixes.

| model (training set) | n | pristine | realization | family | substitutional | none |
|---|---|---|---|---|---|---|
| DiffCSP (MP-20) | 10,000 | 6.2 | 0.6 | 0.8 | 3.2 | 89.1 |
| DiffCSP++ | 10,000 | 4.8 | 0.7 | 0.8 | 3.2 | 90.5 |
| CrystalFormer | 999 | 6.2 | 1.2 | 0.6 | 3.9 | 88.1 |
| SymmCD | 9,475 | 5.0 | 0.5 | 0.5 | 3.5 | 90.5 |
| MiAD | 9,996 | 9.3 | 2.4 | 1.8 | 5.4 | 81.1 |
| MatterGen | 10,000 | 4.5 | 0.8 | 1.2 | 4.3 | 89.2 |
| Crys-JEPA | 10,000 | 6.3 | 0.6 | 1.6 | 5.1 | 86.5 |
| Compositional | 10,000 | 3.1 | 0.0 | 0.0 | 0.0 | 96.8 |
| WyFormer + DiffCSP++ (10k) | 9,999 | 3.4 | 0.3 | 0.6 | 2.7 | 93.0 |
| WyFormer + CrySPR | 1,000 | 4.1 | 0.4 | 0.9 | 3.1 | 91.5 |
| WyckoffLLM-vanilla + DiffCSP++ | 1,000 | 13.5 | 0.1 | 2.4 | 1.2 | 82.8 |
| WyCryst + CrySPR (MP-20 biternary) | 994 | 11.1 | 1.9 | 1.3 | 4.5 | 81.2 |
| SymmCD (MPTS-52), WyFormer (MPTS-52), MP-20 test control | | *rerun in progress* | | | | |

`sibling_ordering` is 0.0–0.1% everywhere. The full per-set table and per-structure results are
in `benchmarks/figshare/results/`.

How to read these numbers:

- **Disorder is a large part of the matches.** Realization, family and substitution together
  are comparable to the exact-match rate. MP-20-trained models reproduce the orderings of
  disordered phases that MP-20 itself contains; see section 5.
- **The loose-tolerance retry.** Many generated structures are unrelaxed, so exact matching is
  retried at loose tolerances (0.3 / 0.5 / 10°) with a density guard. It raises every set's
  `pristine` rate, including the `Compositional` baseline (0.4 → 3.1%), whose geometries are
  unphysical. The split between default and loose tolerance will be reported per set.

## 5. COD coverage and what MP-20 really contains

- **COD coverage.** Of the MP-20 test structures MP flags as experimental
  (`theoretical=False`, 3,929), COD recognises 49%. For 1,824 of the 1,999 unmatched ones the
  formula appears nowhere in COD. They are mostly intermetallics, pnictides, borides,
  silicides and oxides, typically only in ICSD or the Pauling File.
  `benchmarks/figshare/mp20_test_experimental_not_in_cod.csv` lists them, as a work-list for
  curating structures from the original papers and depositing them into COD.
- **MP-20's "theoretical" half.** MP flags 53% of MP-20 test entries as `theoretical`,
  although CDVAE describes MP-20 as "originally from ICSD". The flags are the same in MP's
  2022 release. The explanation is that MP flags **every** ordering of a disordered ICSD
  structure as theoretical: `new_disordered_icsd_2017` covers 14,000 materials, 100% flagged
  theoretical, 0% ICSD-linked. An ordered approximant matches no ICSD record. At least 1,379 of
  the 4,833 "theoretical" MP-20 test entries carry such tags. `opendisorder` grades most of the
  theoretical entries it can match as realizations or substitutions (1,141), and only 104 as
  pristine.
- **Consequence.** "Novel with respect to COD" overstates novelty for MP-20-trained models. Two
  better options:
  1. normalise by the experimental control (about 47% known);
  2. use COD *plus* the training set as the reference, so that COD contributes the
     disorder-aware check against experimental records.

## 6. Other work

- **Open alternatives to ICSD** (research only). No open bulk source has the ICSD-only
  *disordered* phases.
  - Materials Cloud MC3D, OQMD and MP are CC-BY but ordered and partly ICSD-derived; whether
    to use them is a policy decision.
  - The durable open route is to curate phases from the original papers and deposit them in
    COD; COD's policy allows this, since coordinates are facts.
- **Reconstructing disordered parents from MP ordered structures, from first principles**
  (`benchmarks/mp_cod_links`, not in the library). This uses 12.9k MP↔COD links built with our
  matcher, and orbit averaging over similar-element or vacancy hypotheses. The top-ranked
  parent is the COD phase in 75% of non-trivial cases (51% exactly), and 82% of other
  orderings of the same phase are recognised. With all hypothesis tiers, the false-positive
  probe is 27%; with high-precision grouping only, 2%. It can tell that a parent is plausible,
  not that it was observed.
- **Round-trip benchmark** (`benchmarks/roundtrip`, designed and smoke-tested, to be run on a
  GPU cluster).
  - Pipeline: 345 high-quality COD disordered parents → orderings (enumlib; an MP-like Ewald
    recipe; random; SQS; plus controls) → relaxation with a permissively licensed potential
    (MACE-MPA-0, ORB-v3) → matched back against the parent and same-chemical-system decoys.
  - Smoke test: 71/71 orderings recovered, 0 false positives, and all controls graded
    correctly.

## 7. Limitations

- **COD coverage** (section 5).
- **Unrelaxed generated structures.** Distorted structures need loose tolerances, and those
  also let unphysical geometries match.
- **Sublattice pooling** cannot tell which of two same-element sublattices the atoms occupy
  (e.g. Li on spinel 8a vs 16c).
- **Structure-type analogues** with different elements are not matched.
- **Similarity rules** are heuristics. Ionic similarity depends on pymatgen's oxidation-state
  guesses.
- **COD records are taken at face value.** Split-site positional disorder is treated as
  partial occupancy.

## Reproduce

```bash
pip install -e '.[test]' && pytest                       # 33 tests
opendisorder --db cod.sqlite build-index                 # metadata
opendisorder --db cod.sqlite import-archive              # all CIFs, streamed
opendisorder --db cod.sqlite check my.cif --exhaustive   # a single structure
python benchmarks/literature/run.py --db cod.sqlite
python benchmarks/mattergen_baseline/compare.py --db cod.sqlite
FIGSHARE_DATA=... python benchmarks/figshare/run_benchmark.py --db cod.sqlite --out results --processes N
```
