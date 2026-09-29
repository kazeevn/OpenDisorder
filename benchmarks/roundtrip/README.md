# Round-trip benchmark: disordered COD entry → ordered realizations → MLIP relaxation → match back

This benchmark asks whether `opendisorder` recognises ordered realizations of a disordered
COD entry as that entry, and how often it wrongly assigns them to a different one. It starts
from experimentally determined structures with compositional disorder, so the correct answer
is known for every structure.

1. **Select** COD entries with a high-quality refinement and purely compositional disorder
   (substitutional mixing and/or vacancies) that small supercells can represent exactly
   (`select_entries.py`).
2. **Order** each entry: several orderings per supercell size, from several generators
   (random, Ewald-ranked "MP-like", enumlib symmetry-distinct, SQS), plus negative and
   robustness controls (`generate_orderings.py`).
3. **Relax** positions and cell with one or more machine-learned interatomic potentials
   (`relax.py`, pluggable backends, resumable and sharded).
4. **Match back**:
   - `DisorderMatcher.match` against the parent, which measures recall;
   - `DisorderMatcher.match` against every other disordered COD entry of the same chemical
     system (the "decoys"), which measures false positives;
   - optionally the full `NoveltyChecker` pipeline.

   Metrics are reported by supercell size, disorder type, relaxation backend, generation
   method and structure class (`evaluate.py`).

Nothing is run by default at scale. `run_all.sh` runs everything on one machine, and
`slurm/pipeline.sh` submits the same steps as dependent array jobs.

## Why a round trip, and what it mirrors

- **MatterGen** (Zeni et al., *Nature* 639, 624 (2025), SI D.4) is the only published
  evaluation of an ordered–disordered matcher we found:
  - 24,256 ordered approximations of 1,196 disordered alloys, made with pymatgen's
    `EnumerateStructureTransformation` (enumlib);
  - 100 % recall and 100 % precision;
  - the authors note that precision "could be inflated since only 32 out of 1,196 disordered
    structures share the same reduced formula".

  We mirror the enumlib arm (`--methods enum`) and address the caveat in three ways:
  - decoys come from the **same chemical system**, not the same formula, and include every
    COD entry that `NoveltyChecker` would consider (typically 1–100 per parent);
  - explicit negative controls sit close to the decision boundary;
  - MLIP relaxation is added, since realistic queries (generative-model outputs, MP/GNoME
    entries) are relaxed structures, not ideal colorings.
- **Order–disorder family trees** (Yamazaki et al., arXiv:2604.21386, 2026) map ordered
  structures onto disordered ICSD parents by group–subgroup relations and report the fraction
  of generated or database structures with a disordered parent. That is a population
  statistic without a ground truth per structure. The round trip supplies that ground truth.
- **DMFlow** (arXiv:2602.04734) built a COD benchmark of substitutionally and positionally
  disordered structures (3–50 atoms per cell). It scores generation with StructureMatcher
  match rates and does not check matching back after ordering and relaxation.
- **The Materials Project** built many of its entries from disordered ICSD structures. We
  reproduce its recipe as the `ewald` / `mp_like` method; see
  [How MP orders disordered structures](#how-the-materials-project-orders-disordered-structures).

## Step 1: quality criteria for COD parents (`select_entries.py`)

Stage A reads the COD CSV dump (528,535 rows in the 2026-09-29 dump used here). Stage B
reads the CIF from the local index, opened read-only.

| criterion | default | why |
|---|---|---|
| `has coordinates`; COD `has disorder` flag **or** fractional formula | required | same candidate definition as `NoveltyChecker` layer 2 |
| not retracted / `errors` / `warnings` status, not `theoretical` | required | the dump query already drops retracted, errors, duplicates and theoretical entries; the 35 `warnings` entries are dropped too |
| no H/D; no noble gases or elements beyond Pu | required | X-ray H positions are unreliable and H-bearing entries are mostly molecular; MPtrj/OMat24-trained MLIPs cover H–Pu except noble gases |
| 2–6 elements | required | unary "disorder" is vacancies in elements; more than 6 elements is rare and slow |
| R1 = `Robs` (else `Rall`) ≤ **0.05**; entries without R rejected | `--max-r 0.05` | occupancies from poor refinements are unreliable; 0.05 is a common threshold for well-refined inorganic single-crystal data |
| single crystal: not `method = powder diffraction`, no `_pd_*` tags | required | Rietveld occupancies are strongly correlated with displacement parameters. COD's `method` field is empty for 99 % of rows, so the CIF tags are what actually decides |
| ≤ 80 atoms in the conventional cell (COD `cellformula`) | `--max-cell-atoms` | keeps supercells MLIP-affordable |
| 273–313 K and ≤ 1 MPa (when stated) | `--tmin/--tmax/--max-pressure-kpa` | MLIP relaxations are at 0 K / 0 GPa; high-T/high-P phases would legitimately transform |
| pymatgen parses without warnings; CIF formula vs parse within 0.01 atomic fraction | `--stoich-tol` | the parse is exactly what `NoveltyChecker` sees; stoichiometry warnings flag broken occupancies or missing atoms |
| site occupancy sums ≤ 1.02, no overfull orbit | `--occupancy-tolerance` | an overfull site is a data error |
| **no positional disorder**: no `_atom_site_disorder_group/assembly`, and no two sites (at least one partial) closer than 0.7 × (r_i + r_j) | `--split-factor 0.7` | split atoms (YBCO chain O, off-centre cations) and mutually exclusive interstitials (Li garnets, Cu₂₋ₓSe) are not compositional disorder; a coloring would put atoms 0.2–1.5 Å apart |
| every occupancy p/q with q ≤ 12 (±0.011) | `--max-denominator` | lets small supercells reproduce the occupancies exactly |
| smallest exact supercell ≤ 8 primitive cells and ≤ 160 sites | `--max-supercell(-atoms)` | the size scan below stays affordable |
| minority occupancy ≥ 0.08 | `--min-occupancy` | dilute dopants need huge cells and fall inside the matcher's tolerance anyway |
| ≤ 2 entries per chemical system (lowest R first) | `--max-per-chemsys` | stops a few well-studied systems (YBCO, Li–Mn–O) from dominating |

Result with the defaults:

- Stage A passes 2,096 entries. Stage B passes 379. Stage B rejections:
  - split sites: 525
  - irrational occupancies: 449
  - supercell too large: 369
  - parse warnings: 126
  - disorder-group tags: 105
  - ordered after parsing: 49
  - dilute: 45
  - no CIF: 33
  - parse error: 16
- After the cap per chemical system, **345 parents** remain:
  - disorder type: 221 substitutional, 103 vacancy, 21 both;
  - structure class: 126 oxides, 117 intermetallics, 67 chalcogenides, 25 mixed-anion,
    6 halides, 4 nitrides;
  - 197 have at least one decoy.
- The smallest exact supercell is 1–8 primitive cells.

An **extended preset** gives 3,686 parents (2,860 substitutional, 597 vacancy, 229 both):

```
--max-r 0.06 --allow-missing-r --max-denominator 16 --rational-tol 0.02 --max-supercell 12
--max-supercell-atoms 240 --max-cell-atoms 120 --max-per-chemsys 3
```

`--limit N` then draws a sample stratified over disorder type × structure class. Most of the
extra entries have no R value, mainly older literature data, so report the two presets
separately.

Each parent's row records the columns used for bucketing:

- disorder type;
- structure class (anion-based: oxide, chalcogenide, halide, nitride, mixed-anion,
  intermetallic);
- crystal system;
- smallest exact supercell;
- per-orbit occupancies and their rational approximants (JSON);
- the number of same-chemsys disordered COD entries.

`--rejected` writes every stage-B rejection with its reason.

## Step 2: ordered representations (`generate_orderings.py`)

Supercells are built on the primitive cell of the parsed CIF, with occupancies replaced by
their rational approximants.

- `size` counts primitive cells. `--sizes auto` uses 1, 2 and 3 times the smallest exact
  supercell, up to `--max-atoms 120`. Explicit sizes such as `1,2,3,4` also allow inexact
  sizes; `exact=false` and `composition_error` record the rounding.
- For each size, the `--shapes` most compact Hermite-normal-form supercells are used. Per
  symmetry orbit, element and vacancy counts are fixed by largest-remainder rounding.

| method | tool (licence) | what it produces | notes |
|---|---|---|---|
| `random` | built in (pymatgen/numpy, MIT/BSD) | uniform random colorings per orbit | always available; the baseline |
| `ewald` = `mp_like` | pymatgen `OrderDisorderedStructureTransformation` (MIT) | lowest-Ewald-energy orderings, ranked | oxidation states from `oxi_state_guesses` on the integer composition, BVAnalyzer as fallback; falls back to random, flagged, if neither works. `--ewald-algo fast/complete/best_first` |
| `enum` / `enum_ewald` | pymatgen `EnumerateStructureTransformation` + **enumlib** (MIT; Hart & Forcade, *PRB* 77, 224115 (2008); *PRB* 80, 014120 (2009)) | all symmetry-distinct derivative structures of a size, then a random subset (`enum`) or the N lowest Ewald (`enum_ewald`) | the MatterGen protocol. Needs `enum.x`/`makestr.x` (`conda install -c conda-forge enumlib`, tested). Enumeration can explode combinatorially, so it runs in a child process capped by `--enum-mem-gb` and `--method-timeout` |
| `sqs` | **sqsgenerator** (MIT, pip wheels; Gehringer et al., *Comput. Phys. Commun.* 2023) | one special quasirandom structure per shape (Warren–Cowley SRO closest to random) | each disordered orbit is its own sublattice (`split` mode); vacancies are handled |

**Primary and fallback.**

- The primary symmetry-distinct generator is enumlib through pymatgen, because the MatterGen
  and MP recipes use it. It is not pip-installable, and for cells with many disordered sites
  its output is huge.
- The always-available fallbacks are `random` + `ewald` + `sqs`, which need only pip wheels.

**Considered, not integrated.**

- `icet` (MPL-2.0; Ångqvist et al., *Adv. Theory Simul.* 2019) offers SQS and enumeration,
  but v4.0 ships only as source (needs a C++17 build). It is the tool to add for
  cluster-expansion ranking.
- `SHRY` (MIT; Prayogo et al., *J. Chem. Inf. Model.* 62, 2909 (2022)) does canonical
  augmentation, which scales better than enumlib for many substitutions. It pins
  pymatgen ≤ 2023.10.4, numpy < 2 and spglib ≤ 2.4, so it needs its own environment.
- `supercell` (Okhotnikov et al., *J. Cheminform.* 8, 17 (2016)) is GPL and has no Python
  API.

**Controls**, generated from `--controls N` positives per parent:

| kind | construction | expected |
|---|---|---|
| `ctrl_rattle` | Gaussian displacements (σ = 0.05 Å) plus a random ±2 % strain | still `disordered_realization`: a cheap stand-in for relaxation noise |
| `ctrl_composition` | minority atoms of the mixed sublattices converted to the majority element, vacancies filled, or atoms removed, until the global composition leaves MatterGen's tolerance | `disordered_family` or `none` |
| `ctrl_antisite` | one atom of a disordered sublattice swapped with an atom of an ordered sublattice that cannot host it; composition unchanged | `none` or `family`. This targets the site-level checks that MatterGen's global composition test lacks |

Ordered structures from *other* parents of the same chemical system act as cross-parent
negatives automatically: they are decoy matches in step 4.

Every record carries `id`, `parent`, `method`, `size`, `shape`, `supercell_matrix`, `exact`,
`composition_error`, `rational_error`, `ewald_energy` (eV/atom, when oxidation states exist),
`energy_rank` (within parent × size × method), `kind`, `expected`, `disorder_type`,
`structure_class` and `structure` (`Structure.as_dict()`). Output is JSON lines; plain
`.jsonl` is recommended because it is appended to and resumable per parent.

## Step 3: MLIP relaxation (`relax.py`)

The protocol is ASE `FrechetCellFilter` + FIRE with `fmax` = 0.05 eV/Å and at most 500
steps, positions and cell together. These are the Matbench Discovery settings (Riebesell
et al., *Nat. Mach. Intell.* 2025) used for MACE-MP-0 and most leaderboard models.

- `--fmax 0.02` tightens convergence.
- `--optimizer bfgs|lbfgs` and `--filter exp|unit|none` change the optimizer and cell filter.
- MACE uses float64 by default for relaxations.

Each record logs:

- initial and final energy per atom;
- number of steps and the converged flag;
- final fmax, volume change and largest displacement;
- wall time, and the error if any. A failed relaxation keeps the input structure and is
  counted as `relax_failed`.

| backend | model default | code licence | weights licence | notes |
|---|---|---|---|---|
| `mace` | `medium-mpa-0` (MACE-MPA-0; also `medium` = MACE-MP-0, `medium-omat-0`) | MIT | MIT for MP-0/MPA-0; OMAT variants carry OMat24 terms | recommended default: permissive and well benchmarked |
| `orb` | `orb_v3_conservative_inf_mpa` | Apache-2.0 | Apache-2.0 | fast, conservative forces; the permissive second opinion |
| `mattersim` | `MatterSim-v1.0.0-5M.pth` | MIT | MIT | |
| `chgnet` | `0.3.0` | BSD-3-Clause-LBNL | same | MPtrj-only, older; a useful weaker baseline |
| `sevennet` | `7net-0` | MIT in current releases (GPL-3.0 in older ones) | per checkpoint | `7net-mf-ompa` uses `modal=mpa` |
| `fairchem` | `uma-s-1p1` (eSEN/UMA) | MIT | FAIR Chemistry licence with acceptable-use policy; gated Hugging Face download | not permissive; optional only |
| `emt` | ASE EMT | LGPL | n/a | Al/Cu/Ag/Au/Ni/Pd/Pt(+H,C,N,O) only; tests the optimizer path |
| `noop` | none | n/a | n/a | tests the pipeline; equals "unrelaxed" |

For a permissive-only pipeline use **MACE-MPA-0 + ORB-v3-mpa**, optionally with MatterSim.
Two architectures show whether matching failures after relaxation come from the potential or
from the structure.

Install one MLIP per environment (`requirements-mlip.txt`): they pin different torch
versions.

## Step 4: evaluation (`evaluate.py`) and metrics

**Decoys.** For each parent, the other COD entries of its chemical system flagged
`has disorder` or with a fractional formula. This is `NoveltyChecker`'s layer-2 candidate
set, capped by `--max-decoys` and `--decoy-max-cell-atoms`. Each decoy is related to the
parent once with `DisorderMatcher.match_disordered`, after removing occupancies below
`--trace 0.05`:

| relation | test | meaning |
|---|---|---|
| `equivalent` | `pristine` | the same phase: another refinement, temperature or sample |
| `family` | `disordered_family` | same framework and site pattern, different stoichiometry |
| `unrelated` | `none` | a different phase |

Removing trace occupancies was necessary. In the smoke test, three refinements of the same
Li-deficient spinel (COD 1008595–1008597) carry 0.75 % Ti on an otherwise empty 16c site.
That extra site made `match_disordered` report them as unrelated to their sister entry
1008793, which inflated false positives from 0 to 27 %.

**Prefilter.** A decoy whose composition differs from the query's by more than twice
MatterGen's tolerance (atol 0.02, rtol 0.2) cannot yield a realization. It is skipped and
counted as a negative. `--all-decoys` disables the prefilter.

Per structure and source (unrelaxed / backend:model), results record:

- the match kind against the parent, with sublattice and site deviations and timings;
- the kind against every decoy;
- counts by relation;
- with `--novelty`, `NoveltyChecker.check` on a private working copy of the index: all
  metadata, plus CIFs for the candidates only. Only this copy is writable; the source index
  is opened read-only.

| metric | definition |
|---|---|
| `recall_realization` | P(kind vs parent = `disordered_realization`), over positives (`expected=realization`: orderings and `ctrl_rattle`) |
| `recall_any` | the same, also counting `disordered_family` |
| `struct_fp_hard` / `struct_fp_soft` | fraction of positives with ≥ 1 realization match to an *unrelated* / *family* decoy |
| `pair_fp_unrelated` | realization matches to unrelated decoys ÷ unrelated decoy pairs evaluated |
| `precision_pairs_strict` | TP / (TP + FP) over (structure, reference) realization pairs. TP = parent or equivalent decoy; FP = family or unrelated decoy. `_lenient` counts family as TP |
| `specificity_controls` | P(kind vs parent ≠ `disordered_realization`) for `ctrl_composition` and `ctrl_antisite` |
| `novelty_known`, `novelty_top_parent_or_equiv` | `NoveltyChecker` end to end: fraction reported known; fraction whose strongest match is the parent or an equivalent |
| `novelty_neg_known` | the same "known" rate for the negative controls. A hit on a family decoy can be correct here: extra Li in a spinel is a realization of the Li-rich spinel |
| relax robustness | P(realization after relaxation \| realization before), paired by `source_id` |

Tables in `summary.md` and `summary_*.csv` break these down by source (none vs each MLIP),
method, supercell size, atom-count bin, disorder type (substitutional / vacancy / both),
structure class, record kind and `exact`. Relaxation statistics (converged fraction, steps,
volume change) and timings are included.

## Smoke test (this repository's container, 2026-09-29)

**Inputs.**

- Four parents, chosen by `--ids`:
  - 1005047 Li₀.₆₆Mn₁.₃₄N, substitutional;
  - 1008793 Li₀.₇₅Ti₂O₄ spinel, vacancy;
  - 1008967 Er₂/₃Cu₂S₂, vacancy;
  - 2002237 Pd₇O₈, vacancy.
- Orderings from `random`, `ewald`, `sqs` and `enum`/`enum_ewald` (enumlib 2.0.4 from
  conda-forge through micromamba) at two sizes each, capped at 60 atoms, plus controls.

**Results.**

- 71 positives (63 orderings + 8 rattled): `disordered_realization` against the parent in
  **71/71** (unrelaxed and `noop`).
- No hard or soft false positives against 25 decoys: 4 equivalent, 7 family, 13 unrelated,
  1 ordered.
- Controls:
  - `ctrl_composition`: `disordered_family` in 6/6;
  - `ctrl_antisite`: `none` in 8/8;
- `NoveltyChecker` reported every positive as known, with the parent or an equivalent entry
  as its top match.

**Relaxation.**

- EMT applies only to Pd₇O₈. As expected for a potential with no Pd–O physics, it destroyed
  the structures: 7 % volume change, 0 of 10 still matched. This shows why relaxed-structure
  recall has to be measured per potential.
- Real MLIPs were not run here: the container had about 7 GB of free disk. The code paths
  import lazily and exit with an installation hint.

**Cost.**

- Median time per match: 0.02–0.04 s against the parent, 0.04 s over all decoys of a
  structure.
- `ewald` with `algo=fast`: about 70 s for three orderings of a 36-atom Li–Mn–N cell; under
  1 s for the other parents.
- enumlib at size 6 for Li–Mn–N hit the 120 s limit.

## Scale, compute and how to run

```bash
pip install -e .  &&  pip install -r benchmarks/roundtrip/requirements.txt
conda install -c conda-forge enumlib            # optional, for --methods enum
pip install mace-torch                          # or another backend, one env each
opendisorder build-index && opendisorder import-archive   # the COD index with CIFs

COD_CSV=~/cod/all.csv COD_DB=~/.cache/opendisorder/cod.sqlite OUT=~/rt BACKENDS="mace orb" \
  DEVICE=cuda PROCS=16 benchmarks/roundtrip/run_all.sh
SMOKE=1 COD_DB=... OUT=/tmp/rt benchmarks/roundtrip/run_all.sh      # 3 parents, noop relaxation
bash benchmarks/roundtrip/slurm/pipeline.sh                          # cluster: array jobs
```

Every step is resumable: finished parents, structure ids and (id, source) pairs are skipped.

- Generation and evaluation shard by parent (`--shard i/N`, or `slurm/N` to read
  `$SLURM_ARRAY_TASK_ID`).
- Relaxation shards by structure id.
- `evaluate.py --summarize-only` merges all `results*.jsonl` shards.

**Expected size.** The default preset is 345 parents × 3 sizes × (2 shapes × (5 random +
5 ewald + 1 sqs) + 5 enum), minus duplicates, plus about 9 controls per parent. That is
**about 25–30k structures**, comparable to MatterGen's 24k. The extended preset with
`--limit 1500` gives about 120k.

**Compute.** These are rough estimates; measure with `relax.py --limit 50` first.

- Relaxation dominates: 20–100-atom cells take typically 100–300 FIRE steps.
  - MACE-MPA-0 medium or ORB-v3 on one modern GPU: about 20–50 ms per step, so roughly
    500–1,500 structures per hour per GPU.
  - On 8 CPU cores, per step is about 10× slower: roughly 30–150 structures per hour.
  - About 30k structures × 2 backends is therefore about 40–120 GPU-hours, or about
    1–2 days on a 16-GPU allocation.
- Evaluation takes 0.05–1 s per structure, including decoys and `NoveltyChecker`, and more
  for the few cells above 100 atoms. That is about 10–30 CPU-hours for 90k (structure,
  source) pairs.
- Generation is minutes per parent, except enumlib on big cells, which is capped.

## How the Materials Project orders disordered structures

The Materials Project is a design reference only. No ICSD data or identifiers are used here.

**Sources.**

- MP database v2019.02 "added over 47,000 new materials from orderings of disordered ICSD
  [structures] as well as compounds from the Pauling File"
  ([database versions](https://docs.materialsproject.org/changes/database-versions)).
- Horton et al., *Nat. Mater.* 24, 1522 (2025) say ordered approximations are built
  "continuously" but "selecting an ordered approximation is not unique and may overlook
  crucial details such as solid solubility". MP is "establishing methods to classify these
  ordered approximations": Woods-Robinson, Horton & Persson, *Patterns* 2023, and cluster
  expansions with smol.
- The MPWorks/emmet scripts that produced the 2019 set do not appear to be public. We did
  not find the parameters documented.

**Provenance survey.** We read the transformation history in MP's provenance API
(`/materials/provenance`) for 200,487 materials. Only history, tags and remarks were read;
no database IDs.

- 2,268 materials carry the tag **"Ordering of Disordered Crystal"**. Their history is CCTBX
  symmetry expansion followed by the pre-pymatgen MIT Java code: "Ordering using
  com.cmcweb.highthroughput.OrderedCandidates (2 March 2009)". Some of these first had Li
  removed, with "OrderedCandidates used to pick optimal structure".
- The pymatgen-era histories show this chain:
  1. `OxidationStateDecorationTransformation` (explicit `oxidation_states`);
  2. `SupercellTransformation` (e.g. `[[1,0,0],[0,2,0],[0,0,1]]`, `[[2,0,0],[0,2,0],[0,0,2]]`);
  3. `PartialRemoveSpecieTransformation` with `algo = 2` (`ALGO_BEST_FIRST`, Ewald
     best-first) and a `fraction_to_remove`, used for vacancies;
  4. `OrderDisorderedStructureTransformation`.
- The 2019 orderings carry only "Materials Project Optimized Structure" in today's API. We
  did not find a `new_disordered_icsd_2017` tag among current provenance tags.

**Resulting recipe.**

1. Decorate oxidation states with fixed values or BVAnalyzer / `AutoOxiStateDecoration`.
2. Take a small supercell, usually the primitive or conventional cell, made commensurate
   with the occupancies. pymatgen's `OrderDisorderedStructureTransformation` rounds
   occupancies within `occ_tol = 0.25`.
3. Model vacancies with `PartialRemoveSpecieTransformation`, or as vacancy "species" in the
   ordering transformation.
4. Rank orderings by Ewald energy (`ALGO_FAST` or `ALGO_BEST_FIRST`), or enumerate with
   enumlib (`EnumerateStructureTransformation`, `symm_prec = 0.1`, default
   `max_cell_size = 1`; newer pymatgen can sort by an M3GNet energy).
5. Keep the lowest one or a few.
6. Relax with DFT: `MPRelaxSet`, GGA(+U); today r2SCAN. Distinct results become separate
   mp-ids after structure matching.

Deng, Zhu, Chu & Ong (*Chem. Mater.* 29, 281 (2017)) describe the same procedure for
Li₆PS₅Cl:

- occupancies adjusted to a charge-balanced stoichiometry;
- all 48 symmetry-distinct orderings of the primitive cell enumerated with enumlib and
  relaxed with DFT;
- the lowest one kept.

They note that larger cells trade cost against finding the ground state, and that random
sampling is used when full enumeration is intractable.

**Our `mp_like` (= `ewald`) mode and `enum_ewald`** reproduce these recipes, with these
deviations:

- the parents are COD entries, not ICSD entries;
- occupancies are rounded per symmetry orbit to exact rational values, not charge-balanced;
- several orderings, sizes and shapes are kept, not only the lowest;
- relaxation uses MLIPs, not DFT.

Reporting `mp_like` separately estimates how well `opendisorder` would recognise MP's own
ordered approximations. `benchmarks/mp_cod_links` measures the same thing directly on MP
entries.

## Known limitations

- **Trace occupancies** (< 5 %) are ignored only when relating decoys. The library itself
  still treats a 0.75 % site as a real site in `match_disordered`; see the suggestions below.
- **Composition controls** cannot be built for dilute vacancy disorder. In Pd₇O₈ one missing
  Pd per 15 atoms stays within MatterGen's global tolerance whatever is done to the
  sublattice. Such parents get only antisite and rattle controls.
- **Oxidation states.** `ewald` falls back to random ordering when no charge-balanced guess
  exists; records note this. Ewald ranking is meaningless for intermetallics.
- **enumlib output size.** Enumeration can exceed memory for cells with many disordered
  sites; these are skipped, and the log records it. Symmetry-distinct counts are therefore
  incomplete for large sizes.
- **Relaxation at 0 K** may legitimately transform an ordering into a different structure
  type. Then `none` is physically right, not a matcher failure. Inspect `volume_change`,
  `max_displacement` and the symmetry of such cases before reading relaxed-recall as matcher
  quality.
- **Decoys** are limited to COD entries flagged disordered or with fractional formulas, the
  same default as `NoveltyChecker`. Unflagged disordered entries are not decoys.
- **Parent selection** favours small, well-refined inorganic cells. Recall on large complex
  structures (zeolites, > 160-site supercells) is not measured.

## Files

- `select_entries.py`: step 1.
- `generate_orderings.py`: step 2.
- `relax.py`: step 3.
- `evaluate.py`: step 4 and summaries.
- `common.py`: read-only COD access, JSONL I/O, disorder analysis, sharding, time and memory
  limits.
- `run_all.sh`: one-machine pipeline.
- `slurm/`: `pipeline.sh` plus one sbatch file per step.
- `requirements.txt`, `requirements-mlip.txt`: extra dependencies with licences.
