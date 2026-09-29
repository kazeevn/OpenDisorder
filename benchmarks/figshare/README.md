# Generated crystals (figshare 29145101) against COD

This benchmark checks the generated structures of the figshare dataset
["Generated crystals for WyFormer, DiffCSP, DiffCSP++, WyCryst, SymmCD, CrystalFormer, MiAD, MatterGen, Crys-JEPA, Composition"](https://figshare.com/articles/dataset/29145101)
against the Crystallography Open Database with `opendisorder`.

It uses the raw generated structures, as they come out of each model. WyFormer and WyCryst
produce only Wyckoff templates, so for those the structures built by CrySPR and relaxed with
CHGNet (fixed symmetry) are used instead. The MP-20 test split serves as a control of real,
experimentally derived materials.

## How it was run

- COD: 513,893 CIFs from the COD archive, plus metadata for 527k entries, as of 2026-09-29.
- Matcher: the code at the commit that added these results, with default options.
- Machine: 4 CPUs, about 1 hour in total. No structure hit a timeout or a worker crash.

```bash
python download_figshare.py <dir>                       # dataset
opendisorder --db cod.sqlite build-index && opendisorder --db cod.sqlite import-archive
FIGSHARE_DATA=<dir>/data python run_benchmark.py --db cod.sqlite --out results --processes 4
```

`results/<set>.csv.gz` holds per-structure verdicts: the grade, the best COD id, every matching
COD id and the occupancy deviations. `results/summary.csv` holds the table below.

## Results (% of structures)

These were computed at the commit "Precision fixes from the match audit", with default options.

**known** =  at default tolerances + . **pristine (loose only)**
counts structures that match only at the loose retry tolerances (0.3 / 0.5 / 10°, with a density
guard). Those are mostly distorted or unphysical geometries: 2.7% of  and 4.3% of
, but only 0.5% of the control. Leave them out when comparing models. The
per-structure column  records which tolerance a match needed.

| model                              |     n |   pristine (default tol.) |   pristine (loose only) |   realization |   family |   substitutional |   sibling |   known (default tol.) |   none |
|:-----------------------------------|------:|--------------------------:|------------------------:|--------------:|---------:|-----------------:|----------:|-----------------------:|-------:|
| DiffCSP                            | 10000 |                       5.8 |                     0.4 |           0.6 |      0.8 |              3.2 |       0.0 |                    6.4 |   89.1 |
| DiffCSP++                          | 10000 |                       4.5 |                     0.3 |           0.7 |      0.8 |              3.2 |       0.0 |                    5.2 |   90.5 |
| CrystalFormer                      |   999 |                       6.2 |                     0.0 |           1.2 |      0.6 |              3.9 |       0.0 |                    7.4 |   88.1 |
| SymmCD                             |  9475 |                       4.3 |                     0.6 |           0.5 |      0.5 |              3.5 |       0.0 |                    4.8 |   90.5 |
| MiAD                               |  9996 |                       9.0 |                     0.3 |           2.4 |      1.8 |              5.4 |       0.0 |                   11.4 |   81.1 |
| MatterGen                          | 10000 |                       4.3 |                     0.2 |           0.8 |      1.2 |              4.3 |       0.0 |                    5.2 |   89.2 |
| Compositional                      | 10000 |                       0.4 |                     2.7 |           0.0 |      0.0 |              0.0 |       0.0 |                    0.4 |   96.8 |
| Crys-JEPA                          | 10000 |                       2.0 |                     4.3 |           0.6 |      1.6 |              5.1 |       0.0 |                    2.6 |   86.5 |
| WyFormer+DiffCSP++                 |  1000 |                       3.7 |                     0.2 |           0.4 |      1.1 |              3.3 |       0.0 |                    4.1 |   91.3 |
| WyFormer+DiffCSP++ (10k)           |  9999 |                       3.2 |                     0.2 |           0.3 |      0.6 |              2.7 |       0.0 |                    3.5 |   93.0 |
| WyFormer+CrySPR                    |  1000 |                       3.6 |                     0.5 |           0.4 |      0.9 |              3.1 |       0.0 |                    4.0 |   91.5 |
| WyFormer-letters+DiffCSP++         |   996 |                       3.0 |                     0.3 |           0.3 |      0.5 |              2.9 |       0.0 |                    3.3 |   93.0 |
| WyFormer-harmonic+DiffCSP++        |  1000 |                       2.7 |                     0.3 |           0.6 |      0.5 |              2.4 |       0.0 |                    3.3 |   93.5 |
| WyckoffLLM-naive+DiffCSP++         |   999 |                       2.2 |                     0.0 |           0.4 |      0.3 |              3.7 |       0.0 |                    2.6 |   93.4 |
| WyckoffLLM-site-symmetry+DiffCSP++ |   999 |                       3.2 |                     0.1 |           0.9 |      0.6 |              4.3 |       0.0 |                    4.1 |   90.9 |
| WyckoffLLM-vanilla+DiffCSP++       |  1000 |                      13.0 |                     0.5 |           0.1 |      2.4 |              1.2 |       0.0 |                   13.1 |   82.8 |
| WyCryst+CrySPR                     |   994 |                      10.6 |                     0.5 |           1.9 |      1.3 |              4.5 |       0.0 |                   12.5 |   81.2 |
| WyFormer+CrySPR                    |   999 |                       2.8 |                     0.5 |           0.5 |      0.8 |              3.1 |       0.0 |                    3.3 |   92.3 |
| SymmCD                             |  9170 |                       3.6 |                     0.5 |           0.2 |      0.3 |              1.1 |       0.0 |                    3.8 |   94.3 |
| WyFormer+CrySPR                    |   995 |                       0.8 |                     0.2 |           0.0 |      0.1 |              0.1 |       0.0 |                    0.8 |   98.8 |
| MP-20 test split (control)         |  9046 |                      21.2 |                     0.5 |           6.5 |      2.2 |              7.1 |       0.0 |                   27.7 |   62.5 |

## Reading the numbers

- **COD coverage limits every number.** Even the MP-20 test split, which consists of real
  materials, is only 28% "known" to COD, and 62% have no COD evidence at all. So "no match" here
  means *no COD evidence*, not *novel*. To estimate a model's rate of rediscovering known
  materials, compare its "known" rate with the control's 28%.
- **Unrelaxed structures count against the model.** Most structures are raw generator output
  with `StructureMatcher` default tolerances. Distorted versions of known structures can
  therefore fall through to `none`.
- **Disorder matters.** Among the generated sets, 5–10% of structures fall into a disorder grade
  (realization, family or substitutional ordering), comparable to or larger than the pristine
  rate. For MatterGen, 4.3% are pristine matches, 0.8% are ordered realizations of a disordered
  COD entry, and a further 6.5% are family or substitutional or sibling orderings.
- **Wyckoff-token LLM.** `WyckoffLLM-vanilla` has an unusually high pristine rate (13%) with
  little disorder. This fits a model that reproduces known templates.

## How much of known chemistry COD covers

CDVAE describes MP-20 as materials "originally from ICSD". Materials Project (MP), however, flags
4,833 of the 9,046 test structures as `theoretical`, meaning no matching ICSD entry, and 284 as
deprecated. The flags were fetched through the MP API on 2026-09-29. MP-20's old ids convert to
MP's new ids as base-26 letters; the formulas agree for 299 of 300 sampled entries. The flags
are the same in MP's 2022-10-28 release on AWS Open Data (4,830 of the 4,833), so this is not a
recent change.

The disorder layers reconcile the two descriptions. Of the MP-20 test structures that MP
flags as theoretical:

- 523 are `disordered_realization`s of disordered COD entries;
- 618 are `substitutional_ordering`s;
- 206 are `disordered_family`;
- 104 are `pristine` matches;
- 100 carry MP's "Ordering of Disordered Crystal" tag.

A large part of this "theoretical" half are MP's own **ordered approximants of compositionally
disordered experimental phases**. An ordered cell matches no ICSD record, so MP flags it
theoretical even though its source was an ICSD entry. This is exactly the situation this
library is built to detect. MP-20-trained generative models learn these orderings and produce
more of them.

**MP flags every ordering of a disordered ICSD entry as theoretical.** MP's task tags record
how an entry was made; they were read through the MP API, and no ICSD ids were used.

| MP task tag | materials | in the 2022 release | `theoretical` | ICSD link |
|---|---|---|---|---|
| `new_ordered_icsd_2017` (ordered ICSD structures) | 43,790 | 17,387 | 12% | 93% |
| `new_disordered_icsd_2017` (MP orderings of disordered ICSD structures) | 14,173 | 14,000 | **100%** | **0%** |
| `ordering of disordered crystal` / `Ordering of Disordered Crystal` | 1,355 | 709 | ~100% | ~0% |
| `solid solution oxide` (enumerated rock-salt solid solutions) | 875 | 152 | 100% | 0% |

These tags give an independent check of `opendisorder`. On the MP-20 test split, the grades
against COD by MP origin are:

| MP origin | n | pristine | realization | substitutional | family | none |
|---|---|---|---|---|---|---|
| ordered ICSD | 1,003 | 15.8% | 1.8% | 0.9% | 1.4% | 80.1% |
| ordering of a disordered ICSD entry | 1,385 | **0.6%** | **26.2%** | **18.6%** | 4.8% | 49.7% |
| other / untagged | 6,658 | 27.1% | 3.0% | 5.6% | 1.7% | 62.5% |

MP's orderings of disordered phases almost never match COD as `pristine`. They land in the
disorder grades instead, and the "none" share reflects COD's coverage.

| MP flag | n | known % | any match % | no same-formula COD entry % |
|---|---|---|---|---|
| experimental (`theoretical=False`) | 3,929 | 47.1 | 49.1 | 47.1 |
| theoretical | 4,833 | 13.0 | 30.4 | 65.1 |

**COD covers about half of the experimental entries.** Of the 2,003 experimental entries
without any match, 1,854 have no same-formula COD entry at all; for 1,824 the formula appears
nowhere in COD's metadata. The remaining 150 are mostly different polymorphs, according to the
audit sample.

The gap is dominated by pnictides, carbides, borides and silicides (633), intermetallics (566)
and oxides (275). The list is in `mp20_test_experimental_not_in_cod.csv`. These structures
typically exist only in ICSD or the Pauling File. They are candidates for curation from the
original publications and deposition into COD.

**Consequence.** "No match in COD" is a weak novelty signal for MP-20-trained models: even for
experimentally known materials, COD recognises only about half. Reasonable ways to report
novelty:

1. **Normalise by the control.** Report the known rate over the control's known rate.
2. **Use COD together with the training set.** Match against COD plus the model's own training
   data (MP-20 train), which is standard practice. COD then contributes what the training set
   lacks: disorder-aware matching against experimental records.

## Audit of match precision

The audit sampled 110 hits across the grades, from a run before the precision fixes in commit
"Precision fixes from the match audit".

| grade | true | debatable | false |
|---|---|---|---|
| pristine | 19 | 1 | 0 |
| disordered_realization | 17 | 2 | 1 |
| substitutional_ordering | 13 | 6 | 1 |
| disordered_family | 8 | 5 | 7 |
| sibling_ordering | 7 | 11 | 2 |
| control (MP-20 test) realization + substitution | 10 | 0 | 0 |

The failure modes found, and fixed since, were:

- modulated (superspace) COD CIFs parsed into wrong averages;
- `disordered_family` assigned against dilute terminal solid solutions;
- MatterGen-style Hume-Rothery cliques in `sibling_ordering` that merged cations with anions;
- a Pettifor fallback that grouped heterovalent ions in compositions without a charge-balanced
  oxidation-state guess;
- no density check.

Many "none" results with same-formula COD entries turned out to be genuinely different
polymorphs; a few are distorted, unrelaxed versions of known structures, which the new
loose-tolerance retry addresses. The `Compositional` set has unphysical geometries: the median
shortest interatomic distance is 1.3 Å, and 78% have contacts below 0.6 × the sum of the atomic
radii. That explains its 99% "no match".

The results table above was computed with these fixes applied.
