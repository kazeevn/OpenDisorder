#!/usr/bin/env python
"""Step 3: relax the ordered structures with a machine-learned interatomic potential (MLIP).

Positions and cell are relaxed together (ASE ``FrechetCellFilter`` + FIRE by default, the
Matbench Discovery protocol: fmax 0.05 eV/A, <= 500 steps).  One JSON record per structure
is appended to ``--out``; records already present are skipped, so an interrupted job can
simply be restarted.  ``--shard i/N`` (or ``slurm/N`` to read ``$SLURM_ARRAY_TASK_ID``)
processes a disjoint subset for array jobs.

Backends (all imported lazily; install only the one you use, see requirements-mlip.txt):
  mace       MACE-MP-0 / MACE-MPA-0 / MACE-OMAT-0 foundation models   (code MIT; weights MIT
             for mp-0/mpa-0; *omat* weights inherit the OMat24 licence terms)
  chgnet     CHGNet 0.3.0                                           (BSD-3-Clause-LBNL)
  mattersim  MatterSim v1 1M / 5M                                   (MIT)
  sevennet   SevenNet 7net-0 / 7net-l3i5 / 7net-mf-ompa             (code GPL-3.0 in older
             releases, MIT in current ones; check the checkpoint's own terms)
  orb        ORB v3 conservative (orb_v3_conservative_inf_mpa ...)  (Apache-2.0)
  fairchem   eSEN / UMA via fairchem-core                           (code MIT; weights under
             the FAIR Chemistry / OMat24 licence with an acceptable-use policy, gated download)
  emt        ASE EMT (Al, Cu, Ag, Au, Ni, Pd, Pt + H, C, N, O only) - for testing the optimizer
  noop       no relaxation (copies the structure) - for testing the pipeline

usage:
  python relax.py --in orderings.jsonl --backend mace --model medium-mpa-0 --device cuda \
      --out relaxed_mace.jsonl [--shard 3/16] [--kinds positive,ctrl_rattle]
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

log = logging.getLogger("relax")

DEFAULT_MODELS = {
    "mace": "medium-mpa-0",
    "chgnet": "0.3.0",
    "mattersim": "MatterSim-v1.0.0-5M.pth",
    "sevennet": "7net-0",
    "orb": "orb_v3_conservative_inf_mpa",
    "fairchem": "uma-s-1p1",
    "emt": "emt",
    "noop": "noop",
}
EMT_ELEMENTS = {"Al", "Cu", "Ag", "Au", "Ni", "Pd", "Pt", "H", "C", "N", "O"}


def make_calculator(backend: str, model: str, device: str, dtype: str):
    """Instantiate an ASE calculator; heavy dependencies are imported only here."""
    if backend == "mace":
        from mace.calculators import mace_mp

        return mace_mp(model=model, device=device, default_dtype=dtype)
    if backend == "chgnet":
        from chgnet.model import CHGNet
        from chgnet.model.dynamics import CHGNetCalculator

        return CHGNetCalculator(model=CHGNet.load(model_name=model, verbose=False), use_device=device)
    if backend == "mattersim":
        from mattersim.forcefield import MatterSimCalculator

        return MatterSimCalculator(load_path=model, device=device)
    if backend == "sevennet":
        from sevenn.calculator import SevenNetCalculator

        kwargs = {"modal": "mpa"} if "mf" in model else {}
        return SevenNetCalculator(model, device=device, **kwargs)
    if backend == "orb":
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.calculator import ORBCalculator

        orbff = getattr(pretrained, model)(device=device, precision="float32-high")
        return ORBCalculator(orbff, device=device)
    if backend == "fairchem":
        from fairchem.core import FAIRChemCalculator, pretrained_mlip

        predictor = pretrained_mlip.get_predict_unit(model, device=device)
        return FAIRChemCalculator(predictor, task_name="omat")
    if backend == "emt":
        from ase.calculators.emt import EMT

        return EMT()
    if backend == "noop":
        return None
    raise ValueError(f"unknown backend {backend}")


def relax_one(atoms, calc, a) -> dict:
    """Relax positions + cell; returns the result record (without the structure)."""
    import ase.filters
    from ase.filters import FrechetCellFilter, UnitCellFilter
    from ase.optimize import BFGS, FIRE, LBFGS

    ExpCellFilter = getattr(ase.filters, "ExpCellFilter", UnitCellFilter)  # deprecated in newer ASE

    atoms.calc = calc
    t0 = time.time()
    e0 = atoms.get_potential_energy() / len(atoms)
    v0 = atoms.get_volume()
    pos0 = atoms.get_positions().copy()
    filt = {"frechet": FrechetCellFilter, "exp": ExpCellFilter, "unit": UnitCellFilter, "none": None}[a.filter]
    obj = filt(atoms) if filt else atoms
    opt = {"fire": FIRE, "bfgs": BFGS, "lbfgs": LBFGS}[a.optimizer](obj, logfile=None)
    converged = bool(opt.run(fmax=a.fmax, steps=a.steps))
    forces = obj.get_forces()
    fmax = float(np.sqrt((forces**2).sum(axis=1)).max())
    return {
        "e_initial": e0,
        "e_final": atoms.get_potential_energy() / len(atoms),
        "n_steps": int(opt.nsteps),
        "converged": converged,
        "fmax_final": fmax,
        "volume_change": atoms.get_volume() / v0 - 1.0,
        "max_displacement": float(np.linalg.norm(atoms.get_positions() - pos0, axis=1).max()),
        "time_s": round(time.time() - t0, 2),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", nargs="+", required=True, help="orderings JSON-lines file(s)")
    ap.add_argument("--out", required=True, help="relaxed JSON lines (plain .jsonl recommended: appended)")
    ap.add_argument("--backend", required=True, choices=sorted(DEFAULT_MODELS))
    ap.add_argument("--model", help="model name/path (backend default: see DEFAULT_MODELS)")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--dtype", default="float64", help="mace only: float64 recommended for relaxations")
    ap.add_argument("--fmax", type=float, default=0.05, help="eV/A")
    ap.add_argument("--steps", type=int, default=500)
    ap.add_argument("--optimizer", default="fire", choices=["fire", "bfgs", "lbfgs"])
    ap.add_argument("--filter", default="frechet", choices=["frechet", "exp", "unit", "none"])
    ap.add_argument("--kinds", default="positive,ctrl_rattle,ctrl_composition,ctrl_antisite")
    ap.add_argument("--methods", help="only these generation methods (comma-separated)")
    ap.add_argument("--max-atoms", type=int, default=0)
    ap.add_argument("--shard", help="i/N or slurm/N")
    ap.add_argument("--threads", type=int, default=0, help="torch/OpenMP threads (0 = library default)")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    warnings.filterwarnings("ignore")
    if a.threads:
        os.environ.setdefault("OMP_NUM_THREADS", str(a.threads))
        try:
            import torch

            torch.set_num_threads(a.threads)
        except ImportError:
            pass
    model = a.model or DEFAULT_MODELS[a.backend]
    shard, nshards = common.parse_shard(a.shard)
    kinds = set(a.kinds.split(","))
    methods = set(a.methods.split(",")) if a.methods else None
    done = common.done_ids(a.out)

    todo = []
    for rec in common.iter_jsonl(a.inp):
        if rec["id"] in done or rec.get("kind", "positive") not in kinds:
            continue
        if methods and rec.get("method") not in methods:
            continue
        if a.max_atoms and rec["n_atoms"] > a.max_atoms:
            continue
        if common.stable_shard(rec["id"], nshards) != shard:
            continue
        todo.append(rec)
        if a.limit and len(todo) >= a.limit:
            break
    log.info("%d structures to relax with %s/%s on %s (%d already done)", len(todo), a.backend, model, a.device, len(done))
    if not todo:
        return

    try:
        calc = make_calculator(a.backend, model, a.device, a.dtype)
    except ImportError as exc:
        sys.exit(f"backend '{a.backend}' is not installed ({exc}); see requirements-mlip.txt")

    from pymatgen.io.ase import AseAtomsAdaptor

    n_ok = n_conv = 0
    with common.JsonlWriter(a.out) as w:
        for i, rec in enumerate(todo):
            s = common.record_structure(rec)
            meta = {k: v for k, v in rec.items() if k != "structure"}
            out = {**meta, "source_id": rec["id"], "relax_backend": a.backend, "relax_model": model}
            try:
                if a.backend == "noop":
                    res = {"e_initial": None, "e_final": None, "n_steps": 0, "converged": True, "fmax_final": None,
                           "volume_change": 0.0, "max_displacement": 0.0, "time_s": 0.0}
                    relaxed = s
                else:
                    if a.backend == "emt" and not {el.symbol for el in s.composition.elements} <= EMT_ELEMENTS:
                        raise ValueError("EMT does not support these elements")
                    atoms = AseAtomsAdaptor.get_atoms(s)
                    res = relax_one(atoms, calc, a)
                    relaxed = AseAtomsAdaptor.get_structure(atoms)
                out.update(res)
                out["structure"] = common.structure_to_record(relaxed)
                n_ok += 1
                n_conv += bool(res["converged"])
            except Exception as exc:
                out.update({"error": f"{type(exc).__name__}: {str(exc)[:300]}", "structure": rec["structure"], "converged": False})
            w.write(out)
            if i % 20 == 0 or i == len(todo) - 1:
                log.info("%d/%d relaxed (%d ok, %d converged); last %s %s", i + 1, len(todo), n_ok, n_conv, rec["id"],
                         out.get("error") or f"{out.get('n_steps')} steps {out.get('time_s')} s")


if __name__ == "__main__":
    main()
