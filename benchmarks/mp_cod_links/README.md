# Materials Project ordered structures ↔ COD disordered entries

Research code and data asking whether disordered parents can be reconstructed from Materials
Project (MP) ordered structures **from first principles**. It uses no ICSD identifiers or
ICSD-derived provenance. COD is the ground truth.

## Link set: `links.csv`

`build_links.py` builds 12,928 links: 4,986 `disordered_realization` and 7,942
`disordered_family`.

- **COD side:** inorganic COD entries with partial occupancy, 2–6 elements and at most 60
  atoms per cell. 15,348 parse as genuinely disordered.
- **MP side:** MP summary structures from the same 7,218 chemical systems (36.6k entries).
- **Matching:** each pair is run through `opendisorder.matching.DisorderMatcher.match`.

The non-trivial realization links (the COD entry has a site with majority occupancy ≤ 0.9, and
the MP cell is larger than the COD primitive cell) number 3,229. They cover 2,449 MP entries,
2,092 COD entries and 1,347 chemical systems, including 1,036 vacancy-type links. 85% of the
linked MP entries are marked `theoretical`.

The links were computed with the matcher criteria of 2026-09-29, before sublattice pooling.

## Reconstruction from the ordered structure alone: `reconstruct.py`

`reconstruct.py` uses orbit averaging:

1. Take a hypothesis about which elements share sites (tiers below).
2. Collapse those elements onto one species.
3. Find the higher symmetry with spglib, trying symprec 0.01–0.5.
4. Give each orbit the average occupancy of the real species on it.
5. Symmetrise and reduce to the primitive cell.

Hypothesis tiers, in order:

- `similar`: `opendisorder.matching.similarity_groups`;
- `pair`: same-sign oxidation-state or metal pairs, triples and two independent pairs;
- `vac`: an element's sublattice completed with vacancy sites.

## Evaluation: `evaluate.py` → `eval_summary.json`

The sample is 1,000 MP entries from realization links. On the non-trivial subset (n = 640):

| metric | result |
|---|---|
| any candidate proposed | 84% |
| top-1 parent is the COD phase (`match_disordered` → `pristine`) | 75% |
| top-1 matches exactly (same element set per site, occupancies ±0.15) | 51% |
| other MP orderings of the same COD phase recognised as realizations of the reconstructed parent | 82% |

By disorder type, exact recovery is 68% for substitutional disorder and 32% for vacancy
disorder.

**False-positive probe:** 300 MP entries from chemical systems with no disordered COD entry.

- With all tiers, 27% get some candidate. This is an upper bound, since some of these systems
  may genuinely be disorder-prone.
- With the `similar` tier alone, 2% do (0 of the 181 experimental entries).

**Verdict:** reconstructed parents say that a disordered parent is *plausible*, not that one was
observed. Treat them as a hypothesis layer: a parent is confirmed only when it matches a COD
disordered entry.

## Licensing

MP data is CC BY 4.0 (Jain et al., *APL Mater.* 1, 011002, 2013). Parents computed from MP
structures are adaptations of it. Redistribute them under CC BY 4.0 with attribution to MP and
the material ids, noting that they were modified. COD data is CC0.

## Running

The scripts need `MP_API_KEY` in the environment, sent only as the `X-API-KEY` header. They also
need a COD index built with `opendisorder build-index` and `import-archive`, opened read-only.
Paths inside the scripts point to the original scratch locations and need adjusting.
