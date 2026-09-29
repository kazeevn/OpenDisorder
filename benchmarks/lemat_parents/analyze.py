"""Reconstruct plausible disordered parents for the literature structures without a COD match,
using LeMat-Bulk (MP + Alexandria + OQMD) as a source of *other computed orderings*.

For every literature structure Q with kind ``none`` (Leeman et al. 2024, Cheetham & Seshadri 2024):

1. Parents of Q: orbit-averaging reconstruction (``benchmarks/mp_cod_links/reconstruct.py``,
   tiers ``similar`` / ``pair`` / ``vac``), plus ``similar_full`` -- every site of a similarity
   group given the group's global composition, in the symmetrised cell of the collapsed structure
   (catches orderings that neither lower the symmetry nor enlarge the cell, e.g. Dy/Ho on the two
   RE sites of CeNi3-type DyHo2Rh9).
2. LeMat candidates L:
   * ``same``      -- element set equal to chemsys(Q);
   * ``collapsed`` -- element set a strict subset of chemsys(Q) that still contains every
                      similarity group (after collapsing similar elements the systems coincide);
   * ``sub``       -- any other strict subset containing every element that fully occupies an
                      orbit of a parent of Q (end members such as Pb2Sb2O7 for Pb2(Zr,Sb)2O6.5);
   * ``analogue``  -- the chemical system of the paper's named phase if it contains elements
                      absent from Q (structure-type analogues).
3. Evidence:
   * Q in LeMat: StructureMatcher fit with the same reduced formula (or immutable_id = Q's MP id);
   * L is an ordered realization / family member of a parent of Q
     (``DisorderMatcher().match(L, P)``);
   * parents reconstructed from ``same`` candidates vs parents of Q
     (``DisorderMatcher().match_disordered``);
   * ``collapsed`` candidates identical to Q after collapsing similar elements (sibling orderings
     / end members of one substitutional solid solution);
   * analogue: ``StructureMatcher.fit_anonymous(Q, L)``;
   * parent space group vs the space group of the phase the paper names.
4. Grade (see README): STRONG / PLAUSIBLE / WEAK / NOTHING.

Usage::

    python analyze.py [--subset SCRATCH/lemat_subset.parquet] [--processes 2]
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import sys
import time
import warnings
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
LIT = HERE.parent / "literature"
sys.path.insert(0, str(HERE.parent / "mp_cod_links"))
import reconstruct as R  # noqa: E402  (orbit-averaging reconstruction, unmodified)

from pymatgen.analysis.structure_matcher import StructureMatcher  # noqa: E402
from pymatgen.core import Composition, Lattice, Structure  # noqa: E402
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer  # noqa: E402
from pymatgen.symmetry.groups import SpaceGroup  # noqa: E402

from opendisorder.matching import (  # noqa: E402
    ANION_ELEMENTS, DisorderMatcher, MatchKind, collapse_elements, guess_ionic_radii, has_anion, similarity_groups)

SCRATCH = Path(os.environ.get(
    "LEMAT_SCRATCH",
    str(Path(__file__).resolve().parent / "_data")))
MAX_PARENTS = 8          # parents of Q that are tested
MAX_RECON_L = 40         # 'same' candidates whose own parents are reconstructed
MAX_RECON_SITES = 80
MAX_MATCH_SITES = 200
TIER = {"similar": "similar", "similar_full": "similar", "pair": "pair", "vac": "vac",
        "similar+vac": "vac", "pair+vac": "vac"}

# What the paper says Q is (claim type) and the named phase(s) formula / space groups.
#   disorder  : Q is an ordering of a known (compositionally) disordered phase of ~Q's composition
#   same_comp : the same composition is known in another (higher-symmetry) structure
#   mixture   : the synthesis product is a mixture of other known phases
#   analogue  : Q is isostructural with a known phase of other elements
#   none      : nothing named / implausible
CLAIMS = {
    "Ba6Na2Ta2V2O17": ("mixture", "BaNaVO4 + NaVO3", ["P-3m1", "C2/c"]),
    "Ba6Na2V2Sb2O17": ("mixture", "BaNaVO4 + NaVO3", ["P-3m1", "C2/c"]),
    "CaFe2P2O9": ("same_comp", "CaFe2(PO4)2O crocobelonite", ["Pnma"]),
    "CaGd2Zr(GaO3)4": ("disorder", "Ca0.95Zr0.95Gd2.05Ga4.05O12 garnet", ["Ia-3d"]),
    "CaNi(PO3)4": ("mixture", "CaNi3(P2O7)2 + Ca2P2O7 / Ni2(PO3)4", ["P2_1/c", "P4_1"]),
    "FeSb3Pb4O13": ("disorder", "Pb2(Fe0.5Sb1.5)O6.5 pyrochlore", ["Fd-3m"]),
    "Hf2Sb2Pb4O13": ("disorder", "Pb2(HfSb)O6.5 pyrochlore (+PbHfO3)", ["Fd-3m"]),
    "InSb3(PO4)6": ("disorder", "Sb0.5In0.5P2O7; Sb(Sb0.5In0.5)(PO4)3", ["P2_1/n", "Pna2_1"]),
    "InSb3Pb4O13": ("disorder", "Pb2(In0.5Sb1.5)O6.5 pyrochlore (+In2O3)", ["Fd-3m"]),
    "K4MgFe3(PO4)5": ("same_comp", "K4MgFe3(PO4)5", ["P-42_1c"]),
    "KBaGdWO6": ("mixture", "Ba2GdWO6 double perovskite + Gd2O3", ["Fm-3m"]),
    "KBaPrWO6": ("mixture", "Ba11W4O23 + Pr6O11", ["Fd-3m"]),
    "KMn3O6": ("disorder", "K0.3MnO2 / K0.48Mn1.94O5.18 layered", ["P6_3/mmc"]),
    "KNa2Ga3(SiO4)3": ("mixture", "NaGaSiO4", ["P2_1/n"]),
    "KPr9(Si3O13)2": ("disorder", "Pr9.33Si6O26 apatite (K/Pr disorder)", ["P6_3/m"]),
    "MgCuP2O7": ("disorder", "(Mg,Cu)2P2O7 thortveitite-like", ["C2/m"]),
    "MgNi(PO3)4": ("disorder", "M2(PO3)4 (Mg3(PO4)2, Ni3(PO4)2)", ["P2_1/n", "C2/c"]),
    "MgTi4(PO4)6": ("disorder", "Mg0.5Ti2(PO4)3 NASICON", ["R-3c"]),
    "Mn7(P2O7)4": ("mixture", "Mn2(PO3)4", ["C2/c"]),
    "MnAgO2": ("same_comp", "MnAgO2 crednerite", ["C2/m"]),
    "Na3Ca18Fe(PO4)14": ("same_comp", "Na3Ca18Fe(PO4)14 whitlockite", ["R3c"]),
    "Sn2Sb2Pb4O13": ("disorder", "Pb2SnSbO6.5 pyrochlore (+SnO2)", ["Fd-3m"]),
    "Y3In2(GaO4)3": ("same_comp", "Y3In2Ga3O12 garnet", ["Ia-3d"]),
    "Zr2Sb2Pb4O13": ("disorder", "Pb2(ZrSb)O6.5 pyrochlore", ["Fd-3m"]),
    "K2SmLuB24": ("disorder", "LaB6-type hexaboride (K,Sm,Lu)B6", ["Pm-3m"]),
    "K3Nd(AsO4)2": ("analogue", "K3Nd(PO4)2", ["P2_1/m"]),
    "DyHo2Rh9": ("disorder", "(Dy,Ho)Rh3 CeNi3/PuNi3-type", ["P6_3/mmc", "R-3m"]),
    "Tb16Ho(ErIr4)3": ("disorder", "Ln5Ir3 (Tb5Ir3, Pu5Rh3-type)", ["P4/ncc"]),
    "Dy6Y2Ho11Lu(Cd3Ru)2": ("none", "none named (HEA-like)", []),
    "Dy4Tc3B8Ru": ("disorder", "B2DyRu-type Dy(Ru,Tc)B2", ["Pnma"]),
    "Y3SnPd3Pb2": ("analogue", "PbPtGd (ZrNiAl-type)", ["P-62m"]),
    "Hf4Zr11Nb(Ir2N)4": ("disorder", "Zr4Ir2N (filled Ti2Ni), Zr/Hf/Nb disordered", ["Fd-3m"]),
    "Tb4Er4Tc(B8Os)3": ("disorder", "B6Er2Os-type (Er,Tb)2(Os,Tc)B6", ["Pban"]),
    "Ta16AlV3Te16": ("analogue", "FeTa4Te4", ["Pbam"]),
    "La2Ga4Pt3Rh": ("analogue", "Ce2Pt1.5Ga6.5", ["I4/mmm"]),
    "SrTm(Cu2P)4": ("analogue", "CaCu4P2", ["P4_2/mnm"]),
    "AcTl2P2S7": ("analogue", "RbLa2(P2S6)(PS4)2", ["P2_1/c"]),
    "Nb10GaOs4": ("analogue", "sigma phase (Ge3.9Nb20.4Rh5.9)", ["P4_2/mnm"]),
    "Ac(ErB8)3": ("disorder", "LaB6/CaB6-type (Ac,Er)B6", ["Pm-3m"]),
    "TbSmF30": ("none", "not credible (F2 molecules)", []),
    "AcPO3": ("analogue", "La analogue (LaPO3 / LaPO4)", []),
}
# Ordered formulas of the phases the paper names (only those inside chemsys(Q)); a LeMat entry with
# such a formula that is a family member / realization of a parent (or, for a 'similar' parent, a
# collapsed sibling of Q) shows that the parent has the named phase's framework.
NAMED = {
    "K2SmLuB24": ["SmB6", "LuB6"], "Ac(ErB8)3": ["ErB6", "AcB6"], "DyHo2Rh9": ["DyRh3", "HoRh3"],
    "Tb16Ho(ErIr4)3": ["Tb5Ir3", "Ho5Ir3", "Er5Ir3"], "Dy4Tc3B8Ru": ["DyRuB2"],
    "Hf4Zr11Nb(Ir2N)4": ["Zr4Ir2N"], "Tb4Er4Tc(B8Os)3": ["Er2OsB6"],
    "Ba6Na2Ta2V2O17": ["BaNaVO4", "NaVO3"], "Ba6Na2V2Sb2O17": ["BaNaVO4", "NaVO3"],
    "CaFe2P2O9": ["CaFe2P2O9"], "CaNi(PO3)4": ["CaNi3(P2O7)2", "Ca2P2O7", "Ni(PO3)2"],
    "Hf2Sb2Pb4O13": ["PbHfO3"], "K4MgFe3(PO4)5": ["K4MgFe3(PO4)5"], "KBaGdWO6": ["Ba2GdWO6"],
    "KBaPrWO6": ["Ba11W4O23"], "KNa2Ga3(SiO4)3": ["NaGaSiO4"], "MgCuP2O7": ["Mg2P2O7"],
    "MgNi(PO3)4": ["Mg(PO3)2", "Ni(PO3)2"], "Mn7(P2O7)4": ["Mn(PO3)2"], "MnAgO2": ["MnAgO2"],
    "Na3Ca18Fe(PO4)14": ["Na3Ca18Fe(PO4)14"], "Y3In2(GaO4)3": ["Y3In2(GaO4)3"],
}
ANALOGUE_SYSTEMS = {
    "K3Nd(AsO4)2": {"K", "Nd", "P", "O"}, "Y3SnPd3Pb2": {"Pb", "Pt", "Gd"},
    "Ta16AlV3Te16": {"Fe", "Ta", "Te"}, "La2Ga4Pt3Rh": {"Ce", "Pt", "Ga"},
    "SrTm(Cu2P)4": {"Ca", "Cu", "P"}, "AcTl2P2S7": {"Rb", "La", "P", "S"},
    "Nb10GaOs4": {"Ge", "Nb", "Rh"}, "AcPO3": {"La", "P", "O"}, "Ac(ErB8)3": {"La", "Er", "B"},
}


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


def with_timeout(fn, seconds, default):
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(seconds)
    try:
        return fn()
    except Timeout:
        return default
    except Exception:
        return default
    finally:
        signal.alarm(0)


def sg_numbers(symbols):
    out = set()
    for s in symbols:
        s = s.split("(")[0].strip()
        for cand in (s, s.replace("_", ""), s.replace("/", "/")):
            try:
                out.add(SpaceGroup(cand).int_number)
                break
            except Exception:
                continue
    return out


def sg_of(s, symprec=0.1):
    try:
        sga = SpacegroupAnalyzer(s, symprec)
        return sga.get_space_group_number(), sga.get_space_group_symbol()
    except Exception:
        return 1, "P1"


def row_structure(r) -> Structure:
    return Structure(Lattice(np.array(list(r["lattice_vectors"]), dtype=float)), list(r["species_at_sites"]),
                     np.array(list(r["cartesian_site_positions"]), dtype=float), coords_are_cartesian=True)


def clean(s):
    s = s.copy()
    s.remove_oxidation_states()
    for p in list(s.site_properties):
        s.remove_site_property(p)
    return s


def sim_groups(s):
    """opendisorder similarity groups, minus groups that mix anion and non-anion elements; metal
    groups are kept even without an oxidation-state guess (interstitial nitrides such as Zr4Ir2N)."""
    comp = s.composition
    ionic = guess_ionic_radii(comp)
    groups = [g for g in similarity_groups(comp.elements, ionic=ionic) if len(g) > 1]
    # never merge an anion with a non-anion (with an incomplete oxidation-state guess, Pettifor
    # proximity would otherwise merge e.g. P with S in AcTl2P2S7)
    return [g for g in groups if len({el.symbol in ANION_ELEMENTS for el in g}) == 1]


PLACEHOLDERS = ["Ne", "Ar", "Kr", "Xe", "Rn", "He"]


def coordination_numbers(q: Structure) -> list[int] | None:
    from pymatgen.analysis.local_env import CrystalNN
    try:
        nn = CrystalNN()
        return [nn.get_cn(q, i) for i in range(len(q))]
    except Exception:
        return None


def similar_full(q: Structure, groups):
    """Parent in which the elements of each similarity group are mixed over all their sites that
    have the same coordination number (CrystalNN; all sites in intermetallics), at the ratio found
    on those sites in Q.

    Catches orderings that neither lower the symmetry nor enlarge the cell (e.g. Dy/Ho on the two
    RE sites of CeNi3-type DyHo2Rh9) while refusing to mix sites of different coordination
    (e.g. dodecahedral Y and octahedral In in garnet Y3In2Ga3O12)."""
    if not groups:
        return None
    gsym = [{el.symbol for el in g} for g in groups]
    if has_anion(q.composition):
        cns = coordination_numbers(q)
        if cns is None:
            return None
    else:  # intermetallics: similar metals may share sites of different coordination
        cns = [0] * len(q)
    classes = {}  # (group index, cn) -> Counter of elements
    for site, cn in zip(q, cns):
        for gi, g in enumerate(gsym):
            if site.specie.symbol in g:
                classes.setdefault((gi, cn), Counter())[site.specie.symbol] += 1
    mixed = {k: v for k, v in classes.items() if len(v) > 1}
    if not mixed or len(classes) > len(PLACEHOLDERS):
        return None
    label = {k: PLACEHOLDERS[i] for i, k in enumerate(sorted(classes))}
    mixes = {label[k]: {el: n / sum(v.values()) for el, n in v.items()} for k, v in classes.items()}
    species = []
    for site, cn in zip(q, cns):
        sym = site.specie.symbol
        gi = next((i for i, g in enumerate(gsym) if sym in g), None)
        species.append(label[(gi, cn)] if gi is not None else sym)
    c = Structure(q.lattice, species, q.frac_coords)
    best = None
    for sp in R.SYMPRECS:
        try:
            p = SpacegroupAnalyzer(c, sp).get_refined_structure().get_primitive_structure(tolerance=0.25)
        except Exception:
            continue
        if best is None or len(p) < len(best[1]) - 1e-9:
            best = (sp, p)
    if best is None:
        return None
    sp, p = best
    d = Structure(p.lattice, [mixes.get(x.specie.symbol, {x.specie.symbol: 1.0}) for x in p], p.frac_coords)
    return d, sp


def q_parents(q: Structure):
    """[(label, tier, kind, groups, vac, symprec, index, sg, parent)] best first."""
    groups = sim_groups(q)
    out = []
    try:
        cands = R.reconstruct(q)
    except Exception:
        cands = []
    for c in cands:
        out.append(dict(kind=c.kind, tier=TIER[c.kind], groups=c.groups, vac=c.vacancy, symprec=c.symprec,
                        index=round(c.index, 2), parent=c.parent))
    allowed = {frozenset(el.symbol for el in g) for g in groups}
    # drop 'similar' hypotheses of reconstruct.py that merge an anion (see sim_groups)
    out = [x for x in out if not x["kind"].startswith("similar") or all(frozenset(g) in allowed for g in x["groups"])]
    sf = similar_full(q, groups)
    if sf is not None:
        p, sp = sf
        v0 = q.get_primitive_structure().volume
        # skip if an orbit-averaging candidate already produced the same parent
        if not any(R.describe(x["parent"]) == R.describe(p) for x in out):
            out.append(dict(kind="similar_full", tier="similar", groups=[sorted(e.symbol for e in g) for g in groups],
                            vac=None, symprec=sp, index=round(v0 / p.volume, 2), parent=p))
    order = {"similar": 0, "similar_full": 1, "pair": 2, "similar+vac": 3, "vac": 3, "pair+vac": 4}
    out.sort(key=lambda x: (order[x["kind"]], -x["index"]))
    for i, x in enumerate(out):
        x["label"] = f"P{i}"
        x["sg"], x["sg_symbol"] = sg_of(x["parent"])
    return out[:MAX_PARENTS]


def pure_elements(p: Structure):
    """Elements that fully occupy at least one site of p (a realization must contain them)."""
    out = set()
    for site in p:
        items = list(site.species.items())
        if len(items) == 1 and items[0][1] > 0.9:
            out.add(items[0][0].symbol)
    return out


def collapsed_fractions(comp: Composition, groups) -> dict:
    rep = {e: sorted(g)[0] for g in groups for e in g}
    out = Counter()
    for el, amt in comp.items():
        out[rep.get(el.symbol, el.symbol)] += amt
    tot = sum(v for k, v in out.items() if k != R.VAC)
    return {k: v / tot for k, v in out.items() if k != R.VAC}


def composition_compatible(fl: dict, fp: dict, atol=0.04, rtol=0.25) -> bool:
    """Cheap pre-filter: collapsed fractional compositions of L and parent P roughly agree
    (loose enough for family members / end members and vacancy-filled variants)."""
    for k in set(fl) | set(fp):
        a, b = fl.get(k, 0.0), fp.get(k, 0.0)
        if abs(a - b) > atol + rtol * (a + b) / 2:
            return False
    return True


def energy_pa(r):
    return float(r["energy"]) / int(r["nsites"]) if r["energy"] is not None and not pd.isna(r["energy"]) else np.nan


def analyse(job):
    qrow, rows = job
    t0 = time.time()
    dm = DisorderMatcher()
    sm = StructureMatcher()  # default tolerances for identity
    sm_anon = StructureMatcher(ltol=0.3, stol=0.4, angle_tol=8)
    q = clean(Structure.from_file(LIT / "cifs" / qrow["file"]))
    formula = qrow["formula"]
    claim_type, claim_phase, claim_sgs = CLAIMS.get(formula, ("none", "", []))
    claim_nums = sg_numbers(claim_sgs)
    E = {el.symbol for el in q.composition.elements}
    groups = sim_groups(q)
    rep = {el.symbol: g[0].symbol for g in groups for el in g}
    col = lambda els: frozenset(rep.get(e, e) for e in els)  # noqa: E731
    q_sg, q_sgsym = sg_of(q)
    parents = q_parents(q)
    res = dict(file=qrow["file"], paper=qrow["paper"], formula=formula, source_id=qrow["source_id"],
               q_spacegroup=q_sgsym, q_nsites=len(q), similarity_groups=";".join(
                   "/".join(sorted(e.symbol for e in g)) for g in groups),
               claim_type=claim_type, claimed_phase=claim_phase, claimed_sg=";".join(claim_sgs))

    # ---- candidate sets
    L = []
    for r in rows:
        els = set(r["elements"])
        if els == E:
            kind = "same"
        elif els < E and groups and col(els) == col(E):
            kind = "collapsed"
        elif els < E:
            kind = "sub"
        elif formula in ANALOGUE_SYSTEMS and els <= ANALOGUE_SYSTEMS[formula] and not els <= E:
            kind = "analogue"
        else:
            continue
        L.append((kind, r))
    res["n_lemat_same"] = sum(k == "same" for k, _ in L)
    res["n_lemat_collapsed"] = sum(k == "collapsed" for k, _ in L)
    res["n_lemat_sub"] = sum(k == "sub" for k, _ in L)
    res["n_lemat_analogue"] = sum(k == "analogue" for k, _ in L)
    structs = {}

    def get_s(r):
        k = r["immutable_id"] + "|" + r["config"]
        if k not in structs:
            structs[k] = row_structure(r)
        return structs[k]

    # ---- Q itself in LeMat
    q_rf = q.composition.reduced_formula
    q_in, q_e = [], {}
    for kind, r in L:
        if kind != "same":
            continue
        s = get_s(r)
        same_id = r["immutable_id"] == qrow["source_id"]
        if s.composition.reduced_formula == q_rf and (same_id or with_timeout(lambda: sm.fit(q, s), 60, False)):
            q_in.append(f"{r['immutable_id']}[{r['config']}]")
            q_e[(r["functional"], r["config"])] = min(q_e.get((r["functional"], r["config"]), 9e9), energy_pa(r))
    res["q_in_lemat"] = ";".join(dict.fromkeys(q_in))
    # all LeMat polymorphs/orderings with Q's reduced formula (space group, E/atom per functional+config)
    poly = []
    for kind, r in L:
        if kind == "same" and Composition(r["chemical_formula_reduced"]).reduced_formula == q_rf:
            poly.append((r["functional"], r["config"], energy_pa(r), r["immutable_id"], sg_of(get_s(r))[1]))
    poly.sort()
    res["lemat_same_formula"] = ";".join(
        f"{i}({sg},{fu}/{cf[:6]},{e:.4f})" for fu, cf, e, i, sg in poly[:15])

    # ---- (ii) LeMat candidates as realizations / family members of Q's parents
    rel = {p["label"]: dict(real=[], fam=[], pp=[], pfam=[]) for p in parents}
    q_in_ids = {x.split("[")[0] for x in q_in}
    for p in parents:
        P = p["parent"]
        needed = pure_elements(P)
        pels = {el.symbol for el in P.composition.elements}
        pgroups = [list(g) for g in (p["groups"] or [])]
        fp = collapsed_fractions(P.composition, pgroups)
        for kind, r in L:
            if kind == "analogue" or r["nsites"] > MAX_MATCH_SITES:
                continue
            els = set(r["elements"])
            if not (needed <= els <= pels):
                continue
            if r["immutable_id"] in q_in_ids:
                continue
            if not composition_compatible(collapsed_fractions(Composition(r["chemical_formula_reduced"]), pgroups), fp):
                continue
            s = get_s(r)
            m = with_timeout(lambda: dm.match(s, P), 90, None)
            if m is None:
                continue
            tag = (r["immutable_id"], r["config"], r["functional"], s.composition.reduced_formula, energy_pa(r), kind)
            if m.kind is MatchKind.DISORDERED_REALIZATION:
                rel[p["label"]]["real"].append(tag)
            elif m.kind is MatchKind.DISORDERED_FAMILY:
                rel[p["label"]]["fam"].append(tag)

    # ---- (i) parents reconstructed from 'same' candidates vs parents of Q
    same = [r for k, r in L if k == "same" and r["nsites"] <= MAX_RECON_SITES and r["immutable_id"] not in q_in_ids]
    qf = q.composition.fractional_composition

    def cdist(r):
        f = Composition(r["chemical_formula_reduced"]).fractional_composition
        return sum(abs(f[e] - qf[e]) for e in set(f.elements) | set(qf.elements))

    same.sort(key=lambda r: (round(cdist(r), 3), r["nsites"]))
    n_recon = 0
    for r in same[:MAX_RECON_L]:
        s = get_s(r)
        lp = with_timeout(lambda: R.reconstruct(s), 120, [])
        n_recon += 1
        for c in lp[:MAX_PARENTS]:
            for p in parents:
                m = with_timeout(lambda: dm.match_disordered(c.parent, p["parent"]), 60, None)
                if m is None or m.kind is MatchKind.NONE:
                    continue
                tag = (r["immutable_id"], r["config"], r["functional"], s.composition.reduced_formula, energy_pa(r),
                       c.kind)
                key = "pp" if m.kind is MatchKind.PRISTINE else "pfam"
                if tag[0] not in {t[0] for t in rel[p["label"]][key]}:
                    rel[p["label"]][key].append(tag)
    res["n_same_reconstructed"] = n_recon

    # ---- collapsed siblings (identical to Q after collapsing similar elements)
    sib = []
    if groups:
        gl = [list(g) for g in groups]
        qc = collapse_elements(q, gl)
        for kind, r in L:
            if kind not in ("same", "collapsed") or r["nsites"] > MAX_MATCH_SITES or r["immutable_id"] in q_in_ids:
                continue
            s = get_s(r)
            sc = collapse_elements(s, gl)
            if sc.composition.reduced_formula != qc.composition.reduced_formula:
                continue
            if with_timeout(lambda: sm.fit(qc, sc), 60, False):
                sib.append((r["immutable_id"], r["config"], r["functional"], s.composition.reduced_formula,
                            energy_pa(r), kind))
    res["collapsed_siblings"] = ";".join(list(dict.fromkeys(f"{t[0]}({t[3]})" for t in sib))[:12])
    res["n_collapsed_siblings"] = len({t[0] for t in sib})

    # ---- analogue framework
    ana = []
    for kind, r in L:
        if kind != "analogue":
            continue
        s = get_s(r)
        if Composition(s.composition.reduced_formula).anonymized_formula != q.composition.anonymized_formula:
            continue
        if with_timeout(lambda: sm_anon.fit_anonymous(q, s), 60, False):
            ana.append(f"{r['immutable_id']}({s.composition.reduced_formula},{sg_of(s)[1]})")
    res["analogue_matches"] = ";".join(ana[:8])

    # ---- grade every parent
    named = {Composition(f).reduced_formula for f in NAMED.get(formula, [])}
    sib_formulas = {t[3] for t in sib}
    graded = []
    for p in parents:
        rr = rel[p["label"]]
        mixed_sets = [{el.symbol for el in site.species.elements} for site in p["parent"]
                      if len(site.species.elements) > 1]
        fam_ids = {t[0]: t for t in rr["fam"] + rr["pfam"]}

        def is_ordering(t):  # the family member itself mixes >= 2 elements of one mixed parent site
            els = set(Composition(t[3]).as_dict())
            return any(len(els & m) >= 2 for m in mixed_sets)

        n_same_parent = len({t[0] for t in rr["real"]} | {t[0] for t in rr["pp"]})
        n_fam_order = sum(1 for t in fam_ids.values() if is_ordering(t))
        n_end = len(fam_ids) - n_fam_order
        sim = p["tier"] == "similar"
        if sim:
            n_end += len({t[0] for t in sib} - set(fam_ids))
        sg_ok = claim_type in ("disorder", "same_comp") and p["sg"] in claim_nums
        fw_formulas = {t[3] for t in rr["real"] + rr["pp"] + rr["fam"] + rr["pfam"]} | (sib_formulas if sim else set())
        fw_ok = bool(named & fw_formulas)
        claim_ok = sg_ok or fw_ok
        basis = "+".join(x for x, ok in (("space group", sg_ok), ("framework of named phase in LeMat", fw_ok)) if ok)
        n_order = n_same_parent + n_fam_order
        if sim and n_order >= 1 and claim_ok:
            g = "STRONG"
        elif (claim_ok and (sim or n_order + n_end >= 1)) or (sim and n_order >= 1):
            g = "PLAUSIBLE"
        else:
            g = "WEAK"
        graded.append((["STRONG", "PLAUSIBLE", "WEAK"].index(g), p["symprec"] > 0.2, -n_order, -n_end, p["label"], g, p,
                       n_same_parent, n_fam_order, n_end, claim_ok, basis))
    if graded:
        graded.sort(key=lambda x: x[:5])
        _, _, _, _, lab, g, p, nsp, nfo, nf, cok, basis = graded[0]
        res["parents_graded"] = " | ".join(
            f"{x[4]}:{x[5]}(same={x[7]},ord={x[8]},end={x[9]},claim={'y' if x[10] else 'n'})"
            for x in sorted(graded, key=lambda x: x[4]))
    else:
        g, p = ("WEAK" if (q_in or sib or ana) else "NOTHING"), None
    res["grade"] = g
    res["n_parents"] = len(parents)
    res["parents_all"] = " | ".join(f"{x['label']}:{x['kind']}:{x['sg_symbol']}" for x in parents)
    if p is not None:
        rr = rel[p["label"]]
        res.update(best_parent=p["label"], parent_kind=p["kind"], parent_tier=p["tier"],
                   parent_groups=json.dumps(p["groups"]), parent_vac=json.dumps(p["vac"]),
                   parent_spacegroup=p["sg_symbol"], parent_sg_number=p["sg"], parent_index=p["index"],
                   parent_symprec=p["symprec"], parent_formula=p["parent"].composition.reduced_formula,
                   parent_desc=R.describe(p["parent"]), claim_consistent=cok,
                   claim_basis=basis, n_lemat_same_parent=nsp, n_lemat_family_orderings=nfo,
                   n_lemat_end_members=nf,
                   lemat_realizations=";".join(list(dict.fromkeys(f"{t[0]}({t[3]})" for t in rr["real"]))[:12]),
                   lemat_parent_matches=";".join(list(dict.fromkeys(f"{t[0]}({t[3]},{t[5]})" for t in rr["pp"]))[:12]),
                   lemat_family=";".join(list(dict.fromkeys(f"{t[0]}({t[3]})" for t in rr["fam"] + rr["pfam"]))[:12]))
        # energies: other orderings with Q's reduced formula vs Q's own LeMat entry, same functional+config
        diffs, seen_d = [], set()
        for t in rr["real"] + rr["pp"]:
            key = (t[0], t[2], t[1])
            if t[3] == q_rf and (t[2], t[1]) in q_e and not np.isnan(t[4]) and key not in seen_d:
                seen_d.add(key)
                diffs.append((f"{t[2]}/{t[1]}", t[0], (t[4] - q_e[(t[2], t[1])]) * 1000))
        if diffs:
            res["dE_other_orderings_meV_atom"] = ";".join(f"{i}:{d:+.1f} ({fc})" for fc, i, d in sorted(diffs))
        # lowest-energy same-formula ordering among supporters (for GNoME Q, no Q energy in LeMat)
        same_f = [t for t in rr["real"] + rr["pp"] if t[3] == q_rf and not np.isnan(t[4])]
        if same_f:
            t = min(same_f, key=lambda t: t[4])
            res["lowest_same_formula_ordering"] = f"{t[0]} ({t[2]},{t[1]}) {t[4]:.4f} eV/atom"
        res["_parent_structure"] = p["parent"].as_dict()
        res["_supporting_ids"] = sorted({t[0] for t in rr["real"] + rr["pp"] + rr["fam"] + rr["pfam"]})
    res["seconds"] = round(time.time() - t0, 1)
    return res


def rationale(r):
    if r["grade"] == "NOTHING":
        return "no parent reconstructed from Q and no related LeMat entry"
    bits = []
    if isinstance(r.get("parent_spacegroup"), str):
        bits.append(f"{r['parent_kind']} parent {r['parent_spacegroup']} ({r['parent_formula']})")
        bits.append(f"LeMat: {int(r['n_lemat_same_parent'])} other ordering(s) at Q's composition, "
                    f"{int(r['n_lemat_family_orderings'])} ordering(s) at other compositions, "
                    f"{int(r['n_lemat_end_members'])} end-member(s)")
        bits.append(f"consistent with paper's phase ({r['claim_basis']})" if r["claim_consistent"] else
                    ("paper names a mixture/analogue" if r["claim_type"] in ("mixture", "analogue", "none")
                     else "SG differs from paper's phase"))
    else:
        bits.append("no parent reconstructed")
    if isinstance(r.get("q_in_lemat"), str) and r["q_in_lemat"]:
        bits.append("Q itself in LeMat")
    if r.get("n_collapsed_siblings"):
        bits.append(f"{int(r['n_collapsed_siblings'])} collapsed sibling(s)")
    if isinstance(r.get("analogue_matches"), str) and r["analogue_matches"]:
        bits.append("analogue framework in LeMat")
    return "; ".join(bits)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", type=Path, default=SCRATCH / "lemat_subset.parquet")
    ap.add_argument("--processes", type=int, default=2)
    ap.add_argument("--only", nargs="*", default=None, help="formulas to run (debug)")
    ap.add_argument("--out", type=Path, default=HERE / "results.csv")
    args = ap.parse_args(argv)
    lit = pd.read_csv(LIT / "results.csv")
    lit = lit[(lit["kind"] == "none") & lit["paper"].isin(["Leeman2024", "Cheetham&Seshadri2024"])]
    if args.only:
        lit = lit[lit["formula"].isin(args.only)]
    df = pd.read_parquet(args.subset)
    rows = df.to_dict("records")
    for r in rows:
        r["elements"] = list(r["elements"])
    jobs = []
    for _, q in lit.iterrows():
        s = Structure.from_file(LIT / "cifs" / q["file"])
        E = {el.symbol for el in s.composition.elements}
        extra = ANALOGUE_SYSTEMS.get(q["formula"], set())
        sub = [r for r in rows if set(r["elements"]) <= E or (extra and set(r["elements"]) <= extra)]
        jobs.append((q.to_dict(), sub))
    jobs.sort(key=lambda j: -len(j[1]))
    print(f"{len(jobs)} queries; LeMat subset {len(rows)} rows", flush=True)
    out = []
    with Pool(args.processes, maxtasksperchild=1) as pool:
        for res in pool.imap_unordered(analyse, jobs):
            print(f"{res['formula']:22s} {res['grade']:9s} {res.get('parent_spacegroup', '')!s:10s} "
                  f"same={res['n_lemat_same']} real/pp={res.get('n_lemat_same_parent', 0)} "
                  f"ord={res.get('n_lemat_family_orderings', 0)} end={res.get('n_lemat_end_members', 0)} ({res['seconds']}s)", flush=True)
            out.append(res)
    pdir = HERE / "parents"
    pdir.mkdir(exist_ok=True)
    order = {f: i for i, f in enumerate(lit["file"])}
    out.sort(key=lambda r: order[r["file"]])
    for r in out:
        r["rationale"] = rationale(r)
        ps = r.pop("_parent_structure", None)
        ids = r.pop("_supporting_ids", [])
        r["parent_cif"] = ""
        if ps is not None and r["grade"] in ("STRONG", "PLAUSIBLE"):
            tag = Path(r["file"]).stem
            fn = pdir / f"{tag}__parent.cif"
            P = Structure.from_dict(ps)
            from pymatgen.io.cif import CifWriter
            text = str(CifWriter(P, symprec=0.1))
            header = (
                f"# Reconstructed disordered parent (hypothesis, NOT an observed structure) for {r['formula']}\n"
                f"# Query: {r['file']} ({r['paper']}); reconstruction: {r['parent_kind']} "
                f"groups={r['parent_groups']} vac={r['parent_vac']} symprec={r['parent_symprec']}\n"
                f"# Grade {r['grade']}: {r['rationale']}\n"
                f"# Supporting LeMat-Bulk entries (orderings / end members) (LeMaterial/Entalpic, CC-BY-4.0; derived from\n"
                f"# Materials Project (CC-BY-4.0), Alexandria (CC-BY-4.0) and OQMD (CC-BY-4.0)): "
                f"{', '.join(ids[:20]) or 'none'}\n"
                f"# Derived data, modified (orbit-averaged occupancies); redistribute under CC-BY-4.0 with attribution.\n")
            fn.write_text(header + text)
            r["parent_cif"] = f"parents/{fn.name}"
    res = pd.DataFrame(out)
    first = ["paper", "formula", "source_id", "file", "grade", "rationale", "claim_type", "claimed_phase",
             "claimed_sg", "q_spacegroup", "similarity_groups", "parent_kind", "parent_spacegroup", "parent_formula",
             "parent_desc", "claim_consistent", "claim_basis", "n_lemat_same_parent", "n_lemat_family_orderings",
             "n_lemat_end_members", "lemat_realizations",
             "lemat_parent_matches", "lemat_family", "q_in_lemat", "dE_other_orderings_meV_atom",
             "lowest_same_formula_ordering", "lemat_same_formula", "collapsed_siblings", "n_collapsed_siblings", "analogue_matches"]
    cols = [c for c in first if c in res.columns] + [c for c in res.columns if c not in first]
    res[cols].to_csv(args.out, index=False)
    print(res.groupby(["paper", "grade"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    main()
