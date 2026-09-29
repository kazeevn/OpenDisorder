# Plausible disordered parents from LeMat-Bulk for literature structures with no COD match

This directory looks at the 42 literature structures that `benchmarks/literature` finds **no COD
match** for (`kind = none`): 24 A-Lab targets from Leeman et al. (2024) and 18 GNoME entries from
Cheetham & Seshadri (2024). For each one, it tries to reconstruct a disordered parent and looks
for supporting *computed* orderings in [LeMat-Bulk](https://huggingface.co/datasets/LeMaterial/LeMat-Bulk).

> **What this can and cannot show.** LeMat-Bulk is computational: MP, Alexandria and OQMD relaxed
> DFT structures. A parent reconstructed here, even with several LeMat orderings behind it, shows
> that a disordered parent is **plausible**. It does **not** show that the parent was ever
> observed. Only an experimental record, such as a COD entry, can do that. No ICSD identifiers or
> ICSD data were used. The papers' ICSD-based claims enter only as the formula and space group
> of the phase each paper names.

## Data: LeMat-Bulk

- **Source.** LeMaterial / Entalpic, `LeMaterial/LeMat-Bulk`, dataset commit
  `0dc17eea904b860ad7288141e9870f67f8e6bb2c`, CC-BY-4.0. It is a unified, deduplicated set of
  relaxed structures from Materials Project, Alexandria and OQMD, all three CC-BY-4.0.
- **Access** (`fetch_lemat.py`). The script downloads all 21 parquet shards of the four configs
  (`compatible_pbe`, `compatible_pbesol`, `compatible_scan`, `non_compatible`): 3.76 GB,
  6,725,590 rows. It downloads them one at a time into the scratch directory, filters each
  row group on `elements`, and deletes the shard straight away, so peak disk use is under 200 MB.
- **Rows kept.** A row is kept if its element set (at least 2 elements, at most 250 sites) is a
  subset of the chemical system of a literature structure, or of the system of the phase the
  paper names when that phase is a structure-type analogue with other elements. That leaves
  **56,702 rows**, stored as a 16.6 MB zstd parquet in the scratch directory: 39.6k PBE, 5.8k
  PBEsol, 4.9k SCAN and 6.4k non-compatible.
- **Columns.** `immutable_id, elements, chemical_formula_reduced, nsites, lattice_vectors,
  species_at_sites, cartesian_site_positions, energy, functional, cross_compatibility`, plus the
  config.
- **Runtime.** The fetch takes about 3 minutes and needs no API key.

## Method (`analyze.py`, about 4 minutes on 2 processes)

For each literature structure Q:

1. **Parents of Q.**
   - The orbit-averaging reconstruction of `benchmarks/mp_cod_links/reconstruct.py` is imported
     unchanged, with tiers `similar`, `pair` and `vac`.
   - `similar_full` is added. Each `opendisorder.matching.similarity_groups` group is mixed over
     all of its sites that have the same CrystalNN coordination number (in intermetallics, over
     all of its sites), at the ratio found in Q. This catches orderings that neither lower the
     symmetry nor enlarge the cell, such as Dy/Ho on the two RE sites of CeNi₃-type DyHo₂Rh₉.
     Sites of different coordination, such as dodecahedral Y and octahedral In in a garnet, are
     not mixed.
   - Similarity groups that merge an anion with a non-anion are discarded. Otherwise P and S
     would be merged in AcTl₂P₂S₇, where the oxidation-state guess is incomplete.
2. **LeMat candidate sets.**
   - `same`: element set equal to chemsys(Q).
   - `collapsed`: a strict subset that still covers every similarity group, so the systems are
     equal after `collapse_elements`.
   - `sub`: other subsets that contain every element fully occupying an orbit of a parent. These
     are end members such as Pb₂Sb₂O₇.
   - `analogue`: the named phase's own chemical system.
3. **Evidence.**
   - Q itself in LeMat: `StructureMatcher` fit with the same reduced formula.
   - `DisorderMatcher().match(L, P)` for every candidate L and each parent P of Q (up to 8 per
     Q), after a loose collapsed-composition pre-filter. `DISORDERED_REALIZATION` means another
     ordering at Q's composition. `DISORDERED_FAMILY` is split into:
     - an **ordering at another composition**, when L itself mixes two or more elements of one
       mixed parent site (e.g. KSmB₁₂ for (K,Sm,Lu)B₆);
     - an **end member**, otherwise (e.g. Pb₂Sb₂O₇ for Pb₂(Zr,Sb)₂O₆.₅).
   - Parents are reconstructed from up to 40 `same` candidates, those closest in composition,
     and compared with Q's parents using `DisorderMatcher().match_disordered`. `PRISTINE` counts
     as the same parent and `FAMILY` as another composition.
   - Collapsed siblings: `collapsed`/`same` entries identical to Q once similar elements are
     collapsed, e.g. Tb₅Ir₃, Ho₅Ir₃ and Er₅Ir₃ for Tb₁₆Ho(ErIr₄)₃.
   - Analogues: `fit_anonymous`.
   - Energies: E/atom of other orderings with Q's formula, compared with Q's own LeMat entry
     within the same functional and config.
4. **Consistency with the paper's phase.** A parent is consistent with the named phase if
   either of these holds:
   - the paper claims disorder or the same composition, and the parent's space group
     (spglib, symprec 0.1) equals the named phase's space group;
   - a LeMat entry with the named phase's formula (e.g. Zr₄Ir₂N, Er₂OsB₆, Ba₂GdWO₆, Ni(PO₃)₂)
     realizes or belongs to the family of that parent, or, for a `similar` parent, is a
     collapsed sibling of Q.
5. **Grade.** Every parent is graded and the best is reported. Ties prefer symprec ≤ 0.2, then
   more supporting orderings.

   | grade | criteria |
   |---|---|
   | **STRONG** | `similar` tier, **and** at least 1 other LeMat ordering of the parent (at Q's or another composition), **and** consistent with the paper's phase |
   | **PLAUSIBLE** | consistent with the paper's phase, with a `similar` parent or any LeMat support (orderings or end members); **or** a `similar` parent with at least 1 other LeMat ordering |
   | **WEAK** | a parent without the above, or no parent but LeMat context (Q itself, siblings, an analogue) |
   | **NOTHING** | no parent and no related LeMat entry |

The parent CIFs for STRONG and PLAUSIBLE cases are in `parents/`. Each carries a header with the
hypothesis, grade and supporting LeMat ids. `results.csv` has one row per literature structure,
with all parents (`parents_all`, `parents_graded`) and the evidence columns.

## Results

| paper | STRONG | PLAUSIBLE | WEAK | NOTHING |
|---|---|---|---|---|
| Leeman et al. 2024 (24 A-Lab targets) | 0 | 12 | 12 | 0 |
| Cheetham & Seshadri 2024 (18 GNoME entries) | 0 | 8 | 7 | 3 |

**No STRONG case.** The chemistry is often convincing, but the missing piece is always the same:
LeMat has **no other ordering at Q's composition** for any of the 42 structures, except KBaGdWO₆
and KBaPrWO₆.

- For all 24 A-Lab targets, Q itself is in LeMat, usually as the MP entry plus its Alexandria
  copy. Nobody enumerated alternative orderings.
- For the multi-component GNoME compositions, LeMat usually has nothing in Q's chemical system,
  or only other stoichiometries. Its support therefore comes from end members and collapsed
  siblings, which confirm the *framework* but not the disorder.

### Leeman et al. 2024

| formula | grade | parent (tier, SG, mixed sites) | LeMat support | paper's phase / consistent |
|---|---|---|---|---|
| Ba6Na2Ta2V2O17 | WEAK | — | Q in LeMat | BaNaVO4 + NaVO3 / no |
| Ba6Na2V2Sb2O17 | WEAK | — | Q in LeMat | BaNaVO4 + NaVO3 / no |
| CaFe2P2O9 | WEAK | — | Q in LeMat (Pnma) | crocobelonite CaFe2(PO4)2O [Pnma] / no parent (Q is the known ordered phase) |
| CaGd2Zr(GaO3)4 | PLAUSIBLE | pair, Ia-3d, 12c:Ca⅓Gd⅔ 8a:Zr½Ga½ | end member Gd3Ga5O12 (Ia-3d); Q in LeMat | disordered garnet [Ia-3d] / yes (SG) |
| CaNi(PO3)4 | WEAK | — | Q in LeMat | mixture / no |
| FeSb3Pb4O13 | PLAUSIBLE | pair+vac, Fd-3m, 16d:Fe¼Sb¾, O′½ | end member Pb2Sb2O7 pyrochlore (mp-647520 …) | Pb2(Fe0.5Sb1.5)O6.5 [Fd-3m] / yes (SG) |
| Hf2Sb2Pb4O13 | PLAUSIBLE | pair+vac, Fd-3m, Hf½Sb½, O′½ | Pb2Sb2O7, HfPbO3 | Pb2(HfSb)O6.5 [Fd-3m] / yes |
| InSb3(PO4)6 | WEAK | pair, P2_1/c, In½Sb½ | Sb2(PO4)3 | P2_1/n; Pna2_1 / no |
| InSb3Pb4O13 | PLAUSIBLE | pair+vac, Fd-3m, In¼Sb¾, O′½ | Pb2Sb2O7 | Pb2(In0.5Sb1.5)O6.5 [Fd-3m] / yes (SG) |
| K4MgFe3(PO4)5 | PLAUSIBLE | pair, P-42_1c, 8e:Mg¼Fe¾ | end member K4Fe4(PO4)5 (P-42_1c) | K4MgFe3(PO4)5 [P-42_1c] / yes (SG) |
| KBaGdWO6 | PLAUSIBLE | pair, Fm-3m, 8c:K½Ba½ (rock-salt Gd/W) | **other ordering mp-1523221 (Pnnn), +0.6 meV/atom vs Q**; Ba2GdWO6 | Ba2GdWO6 double perovskite / yes (framework) |
| KBaPrWO6 | WEAK | pair, Fm-3m, K½Ba½ | other ordering mp-1523099, +0.6 meV/atom; Ba2PrWO6 | Ba11W4O23 + Pr6O11 / no (paper names another phase) |
| KMn3O6 | WEAK | pair+vac, C2/m, K½ | Q in LeMat | K0.3MnO2 [P6_3/mmc] / no |
| KNa2Ga3(SiO4)3 | WEAK | — | Q in LeMat | NaGaSiO4 / no |
| KPr9(Si3O13)2 | PLAUSIBLE | pair, P6_3/m, 4f:K¼Pr¾ | end member Pr5Si3O13 (P6_3/m) | Pr9.33Si6O26 apatite [P6_3/m] / yes (SG) |
| MgCuP2O7 | WEAK | — (similar_full rejected: Mg/Cu sites differ in CN) | Q in LeMat, 11 entries with Q's formula | (Mg,Cu)2P2O7 [C2/m] / no |
| MgNi(PO3)4 | PLAUSIBLE | similar_full, C2/c, Mg½Ni½ on both M sites | end members Mg(PO3)2, Ni(PO3)2 (same framework) | M2(PO3)4 / yes (SG + framework) |
| MgTi4(PO4)6 | PLAUSIBLE | pair+vac, R-3c, 6b:Mg½ | MgTi2(PO4)3, Ti2(PO4)3 (R-3c) | Mg0.5Ti2(PO4)3 NASICON [R-3c] / yes (SG) |
| Mn7(P2O7)4 | WEAK | — | Q in LeMat | Mn2(PO3)4 / no |
| MnAgO2 | WEAK | — | Q in LeMat | MnAgO2 [C2/m] / no parent (Q is the known ordered phase) |
| Na3Ca18Fe(PO4)14 | PLAUSIBLE | pair, R3c, Ca½Fe½ and Na⅙Ca⅚ | orderings NaCa10(PO4)7 (mp-1221168 …) | whitlockite [R3c] / yes (SG) |
| Sn2Sb2Pb4O13 | PLAUSIBLE | pair+vac, Fd-3m, Sn½Sb½, O′½ | Pb2Sb2O7, Pb2Sn2O7 (both Fd-3m), PbSnO3 | Pb2SnSbO6.5 [Fd-3m] / yes (SG) |
| Y3In2(GaO4)3 | WEAK | — (Y CN8 / In CN6 not mixed) | Q in LeMat | Y3In2Ga3O12 [Ia-3d] / no parent (Q is the known ordered phase) |
| Zr2Sb2Pb4O13 | PLAUSIBLE | pair+vac, Fd-3m, Zr½Sb½, O′½ | Pb2Sb2O7 | Pb2(ZrSb)O6.5 [Fd-3m] / yes (SG) |

### Cheetham & Seshadri 2024

| formula | grade | parent (tier, SG, mixed sites) | LeMat support | paper's phase / consistent |
|---|---|---|---|---|
| K2SmLuB24 | PLAUSIBLE | pair, Pm-3m, (K½Sm¼Lu¼)B6 | ordering KSmB12 (agm005763225); KB6 | LaB6-type [Pm-3m] / yes (SG) |
| K3Nd(AsO4)2 ×2 | WEAK | — | analogue K3Nd(PO4)2 P2_1/m (mp-17187) fits Q anonymously | phosphate analogue / no |
| DyHo2Rh9 | PLAUSIBLE | similar_full, P6_3/mmc, Dy⅓Ho⅔ on both RE sites | none. LeMat has DyHo2Rh9 only as the PuNi3-type R-3m polytype (agm005945547) | CeNi3-type (Dy,Ho)Rh3 / yes (SG) |
| Tb16Ho(ErIr4)3 | PLAUSIBLE | similar, P4/ncc, 16g:Tb, 4b:Ho¼Er¾ (similar_full: all Ln mixed) | Tb5Ir3, Ho5Ir3, Er5Ir3 all P4/ncc, identical to Q after collapsing | Ln5Ir3 [P4/ncc] / yes (SG + framework) |
| Dy6Y2Ho11Lu(Cd3Ru)2 | WEAK | similar, P6_3/mmc, (Dy,Y,Ho,Lu) | Y10Cd3Ru identical after collapsing | none named / – |
| Dy4Tc3B8Ru | PLAUSIBLE | similar, Pnma, Dy(Tc¾Ru¼)B2 | DyRuB2, DyTcB2 (Pnma), identical after collapsing | B2DyRu [Pnma] / yes (SG + framework) |
| Y3SnPd3Pb2 | WEAK | pair, P-62m, 3f:Sn⅓Pb⅔ | YPdPb end members | PbPtGd analogue / – |
| Hf4Zr11Nb(Ir2N)4 | PLAUSIBLE | similar, Fd-3m, (Zr,Hf) and (Hf,Nb) sites | Zr4Ir2N, Nb4Ir2N (Fd-3m), identical after collapsing | Zr4Ir2N [Fd-3m] / yes (SG + framework) |
| Tb4Er4Tc(B8Os)3 | PLAUSIBLE | similar, Pbam, (Os¾Tc¼) | Er2OsB6, Tb2OsB6 (Pbam), identical after collapsing | B6Er2Os / yes (framework; the paper's "Pban" vs LeMat's Pbam is probably a setting or transcription difference) |
| Ta16AlV3Te16 | WEAK | pair, Pbam, Al¼V¾ | — | FeTa4Te4 analogue [Pbam] / – |
| La2Ga4Pt3Rh | PLAUSIBLE | similar_full, P4/nmm (CaBe2Ge2-type), Pt¾Rh¼ | Q in LeMat (Alexandria agm004517780); orderings La2Ga4PtRh3, LaGa2PtRh (other compositions); La(GaPt)2, La(GaRh)2 | Ce2Pt1.5Ga6.5 analogue [I4/mmm] / no |
| SrTm(Cu2P)4 | WEAK | similar, P4_2/mnm, Sr½Tm½ | Sr(Cu2P)2, Tm(Cu2P)2 end members | CaCu4P2 analogue / – |
| AcTl2P2S7 | NOTHING | — | — | analogue / – |
| Nb10GaOs4 | NOTHING | — | — | sigma phase analogue / – |
| Ac(ErB8)3 | PLAUSIBLE | pair, Pm-3m, (Ac¼Er¾)B6 | ordering AcErB12 (agm005761842); ErB6 | LaB6-type [Pm-3m] / yes (SG + framework) |
| TbSmF30 | WEAK | similar, P6_3/m, Tb½Sm½ | — | not credible / – |
| AcPO3 | NOTHING | — | — | La analogue / – |

### Most convincing cases

- **Rare-earth intermetallics** (Cheetham & Seshadri):
  - Tb₁₆Ho(ErIr₄)₃ reconstructs to Tb₅Ir₃-type P4/ncc, the paper's claim, and collapses onto
    Tb₅Ir₃, Ho₅Ir₃ and Er₅Ir₃ in MP, Alexandria and OQMD.
  - The same holds for Dy₄Tc₃B₈Ru → Dy(Ru,Tc)B₂ (Pnma), Hf₄Zr₁₁Nb(Ir₂N)₄ → (Zr,Hf,Nb)₄Ir₂N
    (Fd-3m) and Tb₄Er₄Tc(B₈Os)₃ → (Tb,Er)₂(Os,Tc)B₆.

  They miss STRONG only because LeMat has no second ordering at a mixed composition.
- **Pb₂(M,Sb)₂O₆.₅ pyrochlores** (Leeman; M = Fe, In, Hf, Sn, Zr). The M/Sb pair plus
  O′-vacancy hypothesis gives Fd-3m with half-occupied O′ in all five cases, exactly the paper's
  claim. The Pb₂Sb₂O₇ (and Pb₂Sn₂O₇) pyrochlores in LeMat are end members of that framework.
  These are `pair`-tier parents (M³⁺/⁴⁺ vs Sb⁵⁺ are not isovalent), so they stay PLAUSIBLE.
- **KBaGdWO₆.** The (K,Ba)₂GdWO₆ double perovskite (Fm-3m) is the only case with a genuine
  second ordering at Q's composition: mp-1523221, Pnnn, 0.6 meV/atom above Q, same PBE-compatible
  set. Two A-site orderings that are practically degenerate are consistent with A-site disorder.
  The paper, however, attributes the product to Ba₂GdWO₆ + Gd₂O₃.
- **Hexaborides.** (K,Sm,Lu)B₆ and (Ac,Er)B₆ are LaB₆-type Pm-3m parents with KSmB₁₂ and AcErB₁₂
  orderings in Alexandria.

## Failure modes and caveats

- **No alternative orderings in LeMat.** MP, Alexandria and OQMD store one ordering, rarely
  several, per composition. Alexandria mostly copies MP or GNoME entries: La₂Ga₄Pt₃Rh is present
  as agm004517780. So the STRONG criterion is almost never testable; we found none. Support comes
  mainly from end members, which confirm the framework, not the disorder.
- **Q is the known ordered phase.** For CaFe₂P₂O₉, MnAgO₂ and Y₃In₂Ga₃O₁₂, the paper's point is prior experimental reports, not disorder. A missing parent is then
  the correct outcome.
- **Mixture claims.** Where the paper attributes the product to a mixture of known phases (the
  Ba₆Na₂M₂V₂O₁₇ pair, CaNi(PO₃)₄, KNa₂Ga₃(SiO₄)₃, Mn₇(P₂O₇)₄, KBaPrWO₆), a parent of Q cannot
  address the claim.
- **Tier reliability.** In `benchmarks/mp_cod_links`, the `pair` and `vac` tiers produced a
  candidate in 27% of MP entries from systems with no disordered COD entry, against 2% for
  `similar`. Pair/vac parents (pyrochlores, garnet, whitlockite, NASICON, perovskites) therefore
  rest on their agreement with the paper's space group, and that agreement is circumstantial.
- **`similar_full` is a formal parent.** It keeps Q's cell and only declares similar elements
  mixed. It is right for (Dy,Ho)Rh₃ but can over-mix inequivalent sites. The coordination-number
  guard catches the Y/In garnet and the Mg/Cu pyrophosphate, but not every case.
- **Polytypes are not recognised as the same parent.** LeMat's DyHo₂Rh₉ is the PuNi₃-type
  polytype, while GNoME's is CeNi₃-type.
- **Space groups.** The paper's space group is compared as an exact number. A parent in a
  subgroup, such as KBaGdWO₆'s other parents, or a setting mismatch (Pbam/Pban) fails the SG
  check. The framework-of-named-phase check covers only formulas present in LeMat.
- **Energies.** Energies are uncorrected VASP totals and are compared only between entries of
  the same formula, functional and config. Q's own GNoME energy is not in LeMat, so the GNoME
  entries have no ΔE.
- **Analogue claims** (K₃Nd(AsO₄)₂ vs K₃Nd(PO₄)₂, CaCu₄P₂, FeTa₄Te₄, PbPtGd, the sigma phase,
  RbLa₂(P₂S₆)(PS₄)₂) are structure-type novelty and fall outside disorder reconstruction. Only
  K₃Nd(AsO₄)₂ has a LeMat analogue that fits anonymously.

## Licence and attribution

The derived files are `results.csv` and the `parents/*.cif` files. They are derived from:

- **LeMat-Bulk**: LeMaterial / Entalpic, CC-BY-4.0, <https://huggingface.co/datasets/LeMaterial/LeMat-Bulk>;
- **Materials Project**: Jain et al., *APL Mater.* 1, 011002 (2013), CC-BY-4.0;
- **Alexandria**: Schmidt et al., CC-BY-4.0;
- **OQMD**: Saal et al., *JOM* 65, 1501 (2013); Kirklin et al., *npj Comput. Mater.* 1, 15010
  (2015), CC-BY-4.0;
- the literature structures in `benchmarks/literature` (MP, CC-BY-4.0; GNoME, CC-BY-4.0).

The parent CIFs are **modified** data: orbit-averaged occupancies and symmetrised cells. Each
lists the LeMat ids it rests on. Redistribute them under CC-BY-4.0 with this attribution.

## Reproduce

```bash
. .venv/bin/activate
python benchmarks/lemat_parents/fetch_lemat.py      # ~3 min, streams 3.8 GB, keeps a 17 MB subset in $LEMAT_SCRATCH
python benchmarks/lemat_parents/analyze.py --processes 2   # ~4 min -> results.csv, parents/
```

`LEMAT_SCRATCH` defaults to the session scratch directory used for this run.
