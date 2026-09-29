"""High-level novelty check of a structure against the COD."""

from __future__ import annotations

import functools
import itertools
import logging
import signal
from collections import defaultdict
from collections.abc import Iterable, Sequence
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field

from pymatgen.core import Composition, Structure

from .index import CODIndex
from .matching import (
    MATCH_STRENGTH,
    DisorderMatcher,
    MatchKind,
    MatchResult,
    PristineMatcher,
    SiblingOrderingMatcher,
    chemsys,
    clean_structure,
    collapse_elements,
    guess_ionic_radii,
    similarity_groups,
)

logger = logging.getLogger(__name__)


@dataclass
class NoveltyReport:
    """Result of checking one structure.

    ``kind`` is the strongest match found; ``matches`` lists every matching
    COD entry that was found (only the first one per layer unless the check
    was exhaustive).
    """

    formula: str
    chemsys: str
    kind: MatchKind = MatchKind.NONE
    matches: list[MatchResult] = field(default_factory=list)
    n_pristine_candidates: int = 0
    n_disorder_candidates: int = 0
    n_substitution_candidates: int = 0
    errors: list[str] = field(default_factory=list)

    @property
    def is_novel(self) -> bool:
        """Novel = no pristine match and no ordered realization of a disordered entry."""
        return self.kind not in (MatchKind.PRISTINE, MatchKind.DISORDERED_REALIZATION)

    @property
    def is_novel_strict(self) -> bool:
        """Novel under a strict reading: no match of any kind, i.e. also not a
        different stoichiometry of a known disordered phase (the case of
        MatterGen's TaCr2O6 vs. Ta1/2Cr1/2O2, Juelsholt 2026) nor a different
        ordering of similar elements on a known structure."""
        return self.kind is MatchKind.NONE

    @property
    def cod_ids(self) -> list[str]:
        return [m.reference_id for m in self.matches if m.reference_id]

    def best(self) -> MatchResult | None:
        return min(self.matches, key=lambda m: MATCH_STRENGTH[m.kind]) if self.matches else None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["kind"] = self.kind.value
        d["is_novel"] = self.is_novel
        d["is_novel_strict"] = self.is_novel_strict
        d["matches"] = [{"cod_id": m.reference_id, "kind": m.kind.value, **m.details} for m in self.matches]
        return d


@functools.lru_cache(maxsize=65536)
def _ionic_radii(reduced_formula: str) -> dict:
    return guess_ionic_radii(Composition(reduced_formula))


class _Timeout(Exception):
    pass


@contextmanager
def _time_limit(seconds: float | None):
    """SIGALRM-based time limit (main thread of a Unix process only)."""
    if not seconds or not hasattr(signal, "SIGALRM"):
        yield
        return

    def handler(signum, frame):
        raise _Timeout()

    try:
        old = signal.signal(signal.SIGALRM, handler)
    except ValueError:  # not in main thread
        yield
        return
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


@dataclass
class CheckOptions:
    check_disorder: bool = True
    check_siblings: bool = True
    exhaustive: bool = False
    disorder_candidates_flagged_only: bool = True
    """Layer 2 references are COD entries of the same chemical system that COD
    flags ``has disorder`` or whose formula is fractional, plus same-formula
    entries.  ``False`` downloads every entry of the chemical system and
    decides from the parsed CIF (slower, catches unflagged disorder)."""
    max_cell_atoms_factor: float | None = 4.4
    """Skip COD entries whose conventional cell has more than this many times
    the number of atoms of the query.  An ordered realization contains at
    least one primitive cell of the reference, and a conventional cell holds
    at most 4 primitive cells; the extra 10% covers composition tolerance."""
    time_limit: float | None = 30.0
    """Seconds allowed per structure comparison (``None`` = no limit)."""
    check_substitutions: bool = True
    """Layer 2c: compare with COD entries whose elements are a subset of the
    query's after merging chemically similar elements (see
    :func:`opendisorder.matching.similarity_groups`)."""
    substitution_max_mendeleev: float = 6
    substitution_max_radius: float = 0.2
    substitution_ionic: bool = True
    """Also treat ions with the same guessed oxidation state and Shannon radii
    within 15% as similar (e.g. Zr4+/Sn4+ in BaZr0.5Sn0.5O3)."""
    substitution_max_ionic_radius: float = 0.15
    max_substitution_chemsys: int = 64


@dataclass
class Candidates:
    pristine: list[str] = field(default_factory=list)
    """same reduced formula"""
    disorder: list[str] = field(default_factory=list)
    """same chemical system"""
    substitution: list[str] = field(default_factory=list)
    """chemical system a strict subset of the query's modulo similar elements"""
    groups: list = field(default_factory=list)

    def all(self) -> list[str]:
        return list(dict.fromkeys([*self.pristine, *self.disorder, *self.substitution]))


class NoveltyChecker:
    """Check whether structures are known in the COD.

    Layer 1 (pristine): COD entries with the same reduced formula whose parsed
    structure is ordered are compared with :class:`PristineMatcher`.

    Layer 2 (disorder): COD entries with the same chemical system whose parsed
    structure has partial occupancies are compared with
    :class:`DisorderMatcher`; same-formula ordered entries are optionally
    compared with :class:`SiblingOrderingMatcher`.

    Layer 2c (substitution): COD entries whose element set is a strict subset
    of the query's once chemically similar elements are merged; the query is
    flagged if, with each similar-element group collapsed onto one element,
    it matches such an entry (or is an ordered realization of it).
    """

    def __init__(
        self,
        index: CODIndex,
        pristine: PristineMatcher | None = None,
        disorder: DisorderMatcher | None = None,
        sibling: SiblingOrderingMatcher | None = None,
        options: CheckOptions | None = None,
    ):
        self.index = index
        self.pristine = pristine or PristineMatcher()
        self.disorder = disorder or DisorderMatcher()
        self.sibling = sibling or SiblingOrderingMatcher()
        self.options = options or CheckOptions()

    # ------------------------------------------------------------------ candidates

    def candidate_ids(self, structure: Structure) -> Candidates:
        """COD ids for layer 1 (same reduced formula), layer 2 (same chemical
        system) and layer 2c (similar-element substitution)."""
        opts = self.options
        rf = structure.composition.reduced_formula
        cs = chemsys(structure)
        max_atoms = opts.max_cell_atoms_factor * structure.composition.num_atoms if opts.max_cell_atoms_factor else None
        cand = Candidates(pristine=self.index.ids_by_reduced_formula(rf, max_cell_atoms=max_atoms))
        if opts.check_disorder:
            cand.disorder = self.index.ids_by_chemsys(
                cs, disorder_only=opts.disorder_candidates_flagged_only, max_cell_atoms=max_atoms
            )
        if opts.check_substitutions and structure.is_ordered:
            cand.groups, cand.substitution = self._substitution_ids(structure, max_atoms)
        return cand

    def _substitution_ids(self, structure: Structure, max_atoms) -> tuple[list, list[str]]:
        opts = self.options
        ionic = {}
        if opts.substitution_ionic:
            try:
                with _time_limit(5.0):
                    ionic = _ionic_radii(structure.composition.reduced_formula)
            except _Timeout:
                ionic = {}
        groups = similarity_groups(
            structure.composition.elements,
            opts.substitution_max_mendeleev,
            opts.substitution_max_radius,
            ionic=ionic,
            max_ionic_dr=opts.substitution_max_ionic_radius,
        )
        if all(len(g) == 1 for g in groups):
            return groups, []
        own = chemsys(structure)
        target = collapse_elements(structure, groups).composition.reduced_formula
        choices = [
            [c for r in range(1, len(g) + 1) for c in itertools.combinations(g, r)] if len(g) > 1 else [tuple(g)]
            for g in groups
        ]
        ids: list[str] = []
        for n, combo in enumerate(itertools.product(*choices)):
            if n >= opts.max_substitution_chemsys:
                break
            cs = "-".join(sorted(el.symbol for part in combo for el in part))
            if cs == own:
                continue
            for cid in self.index.ids_by_chemsys(cs, max_cell_atoms=max_atoms):
                e = self.index.entry(cid)
                # keep entries that could match after collapsing: same collapsed
                # formula, or (possibly) disordered
                if e and (e["has_disorder"] or e["fractional"] or self._collapsed_formula(e, groups) == target):
                    ids.append(cid)
        return groups, ids

    @staticmethod
    def _collapsed_formula(entry: dict, groups) -> str | None:
        from .index import cod_formula_to_composition

        for field_name in ("formula", "calcformula"):
            comp = cod_formula_to_composition(entry.get(field_name))
            if comp is not None:
                rep = {el: g[0] for g in groups for el in g}
                merged: dict = {}
                for el, amt in comp.items():
                    key = rep.get(el, el)
                    merged[key] = merged.get(key, 0) + amt
                from pymatgen.core import Composition

                return Composition(merged).reduced_formula
        return None

    def prefetch(self, structures: Iterable[Structure], processes: int = 1) -> int:
        """Download and parse all COD entries that ``structures`` need."""
        ids: set[str] = set()
        for s in structures:
            ids.update(self.candidate_ids(s).all())
        self.index.prepare(sorted(ids), processes=processes)
        return len(ids)

    # ------------------------------------------------------------------ check

    def check(self, structure: Structure, fetch: bool = True) -> NoveltyReport:
        q = clean_structure(structure)
        cand = self.candidate_ids(q)
        refs = self.index.get_structures(cand.all(), fetch=fetch)
        return self._check_with_refs(q, cand, refs)

    def _compare(self, report: NoveltyReport, fn, *args) -> MatchResult | None:
        try:
            with _time_limit(self.options.time_limit):
                return fn(*args)
        except _Timeout:
            report.errors.append(f"timeout: {getattr(fn, '__qualname__', fn)}")
        except Exception as exc:  # pragma: no cover - defensive
            report.errors.append(f"{type(exc).__name__}: {exc}")
        return None

    def _check_with_refs(self, q: Structure, cand: Candidates, refs: dict[str, Structure]) -> NoveltyReport:
        opts = self.options
        pristine_ids, disorder_ids = cand.pristine, cand.disorder
        report = NoveltyReport(formula=q.composition.reduced_formula, chemsys=chemsys(q))
        rf = q.composition.reduced_formula

        # ---- layer 1: pristine
        same_formula = [
            cid
            for cid in pristine_ids
            if cid in refs
            and refs[cid].is_ordered == q.is_ordered
            and refs[cid].composition.reduced_formula == rf
        ]
        report.n_pristine_candidates = len(same_formula)
        unmatched_ordered = []
        for cid in same_formula:
            if self._compare(report, self.pristine.fit, q, refs[cid]):
                report.matches.append(MatchResult(MatchKind.PRISTINE, cid))
                if not opts.exhaustive:
                    break
            elif q.is_ordered:
                unmatched_ordered.append(cid)
        if report.matches and not opts.exhaustive:
            return self._finish(report)

        if not opts.check_disorder:
            return self._finish(report)

        if not q.is_ordered:
            # disordered query: same disordered phase with other occupancies / family
            others = [cid for cid in dict.fromkeys(disorder_ids) if cid in refs and not refs[cid].is_ordered]
            others = [cid for cid in others if cid not in same_formula]
            report.n_disorder_candidates = len(others)
            for cid in others:
                res = self._compare(report, self.disorder.match_disordered, q, refs[cid])
                if res is not None and res.kind is not MatchKind.NONE:
                    res.reference_id = cid
                    report.matches.append(res)
                    if res.kind is MatchKind.PRISTINE and not opts.exhaustive:
                        break
            return self._finish(report)

        # ---- layer 2a: ordered realization of a disordered entry
        # same-formula disordered entries are layer-2 references too
        disordered = [
            cid for cid in dict.fromkeys([*disorder_ids, *pristine_ids]) if cid in refs and not refs[cid].is_ordered
        ]
        report.n_disorder_candidates = len(disordered)
        found_realization = False
        for cid in disordered:
            res = self._compare(report, self.disorder.match, q, refs[cid])
            if res is not None and res.kind is not MatchKind.NONE:
                res.reference_id = cid
                report.matches.append(res)
                if res.kind is MatchKind.DISORDERED_REALIZATION:
                    found_realization = True
                    if not opts.exhaustive:
                        break
        if found_realization and not opts.exhaustive:
            return self._finish(report)

        # ---- layer 2c: ordering of similar elements on a known structure with fewer elements
        subst = [cid for cid in cand.substitution if cid in refs]
        report.n_substitution_candidates = len(subst)
        if subst:
            qc = collapse_elements(q, cand.groups)
            for cid in subst:
                res = self._compare(report, self._substitution_match, qc, collapse_elements(refs[cid], cand.groups))
                if res is not None:
                    res.reference_id = cid
                    report.matches.append(res)
                    if not opts.exhaustive:
                        break

        # ---- layer 2b: different ordering of similar elements (same formula)
        if opts.check_siblings and unmatched_ordered:
            dq = self._compare(report, self.sibling.prepare, q)
            if dq is not None:
                for cid in unmatched_ordered:
                    if self._compare(report, self.sibling.fit, dq, q, refs[cid]):
                        report.matches.append(MatchResult(MatchKind.SIBLING_ORDERING, cid))
                        if not opts.exhaustive:
                            break
        return self._finish(report)

    def _substitution_match(self, qc: Structure, rc: Structure) -> MatchResult | None:
        if rc.is_ordered:
            if self.pristine.fit(qc, rc):
                return MatchResult(MatchKind.SUBSTITUTIONAL_ORDERING, details={"via": "pristine"})
            return None
        res = self.disorder.match(qc, rc)
        if res.kind is MatchKind.DISORDERED_REALIZATION:
            return MatchResult(MatchKind.SUBSTITUTIONAL_ORDERING, details={"via": "disordered_realization", **res.details})
        return None

    @staticmethod
    def _finish(report: NoveltyReport) -> NoveltyReport:
        best = report.best()
        report.kind = best.kind if best else MatchKind.NONE
        report.matches.sort(key=lambda m: MATCH_STRENGTH[m.kind])
        return report

    # ------------------------------------------------------------------ batch

    def check_many(
        self,
        structures: Sequence[Structure],
        processes: int = 1,
        fetch: bool = True,
        progress: bool = False,
    ) -> list[NoveltyReport]:
        """Check many structures; structures are grouped by chemical system so
        that each COD reference is loaded once per group.

        Downloads happen in the calling process; matching is parallelised over
        chemical-system groups with ``processes`` workers.
        """
        cleaned = [clean_structure(s) for s in structures]
        if fetch:
            self.prefetch(cleaned, processes=processes)
        groups: dict[str, list[int]] = defaultdict(list)
        for i, s in enumerate(cleaned):
            groups[chemsys(s)].append(i)
        jobs = [(self.index.path, self._config(), [(i, cleaned[i]) for i in idx]) for idx in groups.values()]
        # largest groups first for better load balancing
        jobs.sort(key=lambda j: -len(j[2]))
        reports: list[NoveltyReport | None] = [None] * len(cleaned)
        if processes > 1:
            with ProcessPoolExecutor(processes) as pool:
                results = pool.map(_check_group_job, jobs, chunksize=1)
                self._collect(results, reports, len(jobs), progress)
        else:
            self._collect((self._check_group(items) for _, _, items in jobs), reports, len(jobs), progress)
        return reports  # type: ignore[return-value]

    @staticmethod
    def _collect(results, reports, n_jobs, progress):
        for j, group in enumerate(results):
            for i, rep in group:
                reports[i] = rep
            if progress and (j % 100 == 0 or j == n_jobs - 1):
                logger.info("checked %d/%d chemical systems", j + 1, n_jobs)

    def _config(self) -> dict:
        return {
            "pristine": asdict(self.pristine),
            "disorder": asdict(self.disorder),
            "sibling": asdict(self.sibling),
            "options": asdict(self.options),
        }

    def _check_group(self, items: list[tuple[int, Structure]]) -> list[tuple[int, NoveltyReport]]:
        out = []
        ref_cache: dict[str, Structure] = {}
        for i, q in items:
            cand = self.candidate_ids(q)
            need = [cid for cid in cand.all() if cid not in ref_cache]
            ref_cache.update(self.index.get_structures(need, fetch=False))
            out.append((i, self._check_with_refs(q, cand, ref_cache)))
        return out


def _check_group_job(job) -> list[tuple[int, NoveltyReport]]:
    import warnings

    warnings.filterwarnings("ignore")
    db_path, config, items = job
    index = CODIndex(db_path)
    try:
        checker = NoveltyChecker(
            index,
            pristine=PristineMatcher(**config["pristine"]),
            disorder=DisorderMatcher(**config["disorder"]),
            sibling=SiblingOrderingMatcher(**config["sibling"]),
            options=CheckOptions(**config["options"]),
        )
        return checker._check_group(items)
    finally:
        index.close()
