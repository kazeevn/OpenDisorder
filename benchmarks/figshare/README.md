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

**known** = `pristine` + `disordered_realization`; **no match** = `none`.
The other grades are described in the top-level README.

| model                              | set                                                        |     n |   pristine |   realization |   family |   substitutional |   sibling |   known |   no match |
|:-----------------------------------|:-----------------------------------------------------------|------:|-----------:|--------------:|---------:|-----------------:|----------:|--------:|-----------:|
| DiffCSP                            | mp_20 / DiffCSP                                            | 10000 |        5.8 |           0.7 |      1.2 |              3.4 |       0.4 |     6.4 |       88.6 |
| DiffCSP++                          | mp_20 / DiffCSP++                                          | 10000 |        4.5 |           0.7 |      1.1 |              3.3 |       0.3 |     5.2 |       90.1 |
| CrystalFormer                      | mp_20 / CrystalFormer                                      |   999 |        6.2 |           1.2 |      0.8 |              4.0 |       0.1 |     7.4 |       87.7 |
| SymmCD                             | mp_20 / SymmCD                                             |  9475 |        4.3 |           0.5 |      0.7 |              3.7 |       0.3 |     4.8 |       90.5 |
| MiAD                               | mp_20 / MiAD                                               |  9996 |        9.0 |           2.4 |      2.4 |              5.5 |       1.0 |    11.4 |       79.7 |
| MatterGen                          | mp_20 / MatterGen / MatterGen_10k                          | 10000 |        4.3 |           0.8 |      1.7 |              4.4 |       0.4 |     5.2 |       88.4 |
| Compositional                      | mp_20 / Compositional                                      | 10000 |        0.4 |           0.0 |      0.1 |              0.1 |       0.1 |     0.4 |       99.4 |
| Crys-JEPA                          | mp_20 / Crys-JEPA                                          | 10000 |        2.0 |           0.6 |      2.7 |              5.5 |       0.3 |     2.6 |       88.9 |
| WyFormer+DiffCSP++                 | mp_20 / WyckoffTransformer / DiffCSP++                     |  1000 |        3.7 |           0.4 |      1.3 |              3.4 |       0.1 |     4.1 |       91.1 |
| WyFormer+DiffCSP++ (10k)           | mp_20 / WyckoffTransformer / DiffCSP++10k                  |  9999 |        3.2 |           0.3 |      0.8 |              2.9 |       0.3 |     3.5 |       92.6 |
| WyFormer+CrySPR                    | mp_20 / WyckoffTransformer / CrySPR / CHGNet_fix           |  1000 |        3.6 |           0.4 |      1.2 |              3.2 |       0.0 |     4.0 |       91.6 |
| WyFormer-letters+DiffCSP++         | mp_20 / WyckoffTransformer-letters / DiffCSP++             |   996 |        3.0 |           0.2 |      1.1 |              2.8 |       0.2 |     3.2 |       92.7 |
| WyFormer-harmonic+DiffCSP++        | mp_20 / WyckoffTransformer-harmonic / DiffCSP++            |  1000 |        2.7 |           0.6 |      0.7 |              2.6 |       0.2 |     3.3 |       93.2 |
| WyckoffLLM-naive+DiffCSP++         | mp_20 / WyckoffLLM-naive / DiffCSP++                       |   999 |        2.2 |           0.4 |      0.4 |              3.7 |       0.1 |     2.6 |       93.2 |
| WyckoffLLM-site-symmetry+DiffCSP++ | mp_20 / WyckoffLLM-site-symmetry / DiffCSP++               |   999 |        3.2 |           0.9 |      0.8 |              4.6 |       0.4 |     4.1 |       90.1 |
| WyckoffLLM-vanilla+DiffCSP++       | mp_20 / WyckoffLLM-vanilla / DiffCSP++                     |  1000 |       13.0 |           0.1 |      2.4 |              1.3 |       0.2 |    13.1 |       83.0 |
| WyCryst+CrySPR                     | mp_20_biternary / WyCryst / CrySPR / CHGNet_fix            |   994 |       10.6 |           1.9 |      1.7 |              4.6 |       0.3 |    12.5 |       80.9 |
| WyFormer+CrySPR                    | mp_20_biternary / WyckoffTransformer / CrySPR / CHGNet_fix |   999 |        2.8 |           0.5 |      1.2 |              3.2 |       0.4 |     3.3 |       91.9 |
| SymmCD                             | mpts_52 / SymmCD                                           |  9170 |        3.6 |           0.2 |      0.4 |              1.1 |       0.2 |     3.8 |       94.4 |
| WyFormer+CrySPR                    | mpts_52 / WyckoffTransformer / CrySPR / CHGNet_fix         |   995 |        0.8 |           0.0 |      0.2 |              0.1 |       0.2 |     0.8 |       98.7 |
| MP-20 test split (control)         | mp_20 / split / test                                       |  9046 |       21.2 |           6.5 |      2.6 |              7.2 |       0.4 |    27.8 |       61.9 |
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

A manual audit of a sample of the matches (precision per grade) is pending and will be added
here.
