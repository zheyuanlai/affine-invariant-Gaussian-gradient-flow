"""Run the covariance-bootstrap experiments and write the final CSVs / metadata.

Four experiments share one output directory:

1. ``covariance_bootstrap_scaling``   -> covariance_bootstrap_summary.csv
2. ``dynamic_contraction_benchmark``  -> contraction_benchmark.csv
3. ``wasserstein_bootstrap_then_fr``  -> wasserstein_bootstrap_summary.csv
4. ``stl_noise_floor``                -> stl_floor_summary.csv

Experiments 1-3 are deterministic and CPU/NumPy only. Experiment 4 (stochastic
STL) runs on the chosen backend/device (NumPy on CPU by default; a single CUDA
torch device with ``--device cuda --backend torch``). The deterministic per-step
rows of experiments 1-3 are concatenated into ``results_long.csv``.

Resume: the deterministic block and the stochastic block are each skipped if
their outputs already exist, unless ``--overwrite`` is passed.

Usage::

    python scripts/natural_gradient_covariance_bootstrap/run_experiments.py \
        --config configs/natural_gradient_covariance_bootstrap/covariance_bootstrap.yaml \
        --outdir outputs/natural_gradient_covariance_bootstrap --overwrite --backend auto
"""
from __future__ import annotations

import argparse
import copy
import os
import platform
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.io_utils import load_yaml, ensure_dir, save_dataframe, save_json
from src.common.torch_utils import torch_device_info, resolve_backend
from src.natural_gradient_covariance_bootstrap.targets import build_target
from src.natural_gradient_covariance_bootstrap.runner import (
    run_scaling, run_contraction, run_wasserstein,
    LONG_COLS, SCALING_SUMMARY_COLS, CONTRACTION_COLS, WBOOT_SUMMARY_COLS,
)
from src.natural_gradient_covariance_bootstrap.stochastic import (
    run_stl_floor, STL_FLOOR_COLS,
)

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "..", "configs",
    "natural_gradient_covariance_bootstrap", "covariance_bootstrap.yaml")

DET_FILES = ["results_long.csv", "covariance_bootstrap_summary.csv",
             "contraction_benchmark.csv", "wasserstein_bootstrap_summary.csv"]
STL_FILES = ["stl_floor_summary.csv"]
OWN_FILES = DET_FILES + STL_FILES + ["target_metadata.json", "run_metadata.json"]


def _progress(stage, i, n, msg):
    print(f"  [{stage:12s} {i:3d}/{n}] {msg}", flush=True)


def _apply_smoke(cfg):
    """Shrink every grid to a fast smoke size (in place on a copy)."""
    cfg = copy.deepcopy(cfg)
    sc = cfg["scaling"]
    sc["d"] = [10]
    sc["kappa"] = [10.0, 100.0]
    sc["lambda0"] = [1e-8, 1e-4, 1e-2]
    sc["n_steps"] = 80
    ct = cfg["contraction"]
    ct["d"] = [10]
    ct["kappa"] = [10.0]
    ct["lambda0"] = [1e-6]
    ct["n_steps"] = 800
    wb = cfg["wasserstein"]
    for spec in wb["targets"]:
        spec["d"] = 6
        spec["gh_nodes"] = 24
    wb["lambda0"] = [1e-8, 1e-4]
    wb["c"] = [0.5]
    wb["n_steps"] = 600
    sf = cfg["stl_floor"]
    for spec in sf["targets"]:
        spec["d"] = 4
        spec["gh_nodes"] = 24
    sf["dt_list"] = [0.0625, 0.015625]
    sf["methods"] = ["kl_raw", "kl_stl"]
    sf["n_seeds"] = 3
    sf["horizon"] = 10.0
    sf["max_steps"] = 400
    sf["gh_nodes"] = 20
    return cfg


def _write_target_metadata(cfg, outdir):
    md = {"targets": {}}
    # scaling / contraction Gaussians (one per kappa at max d).
    for kappa in sorted(set(cfg["scaling"]["kappa"]) | set(cfg["contraction"]["kappa"])):
        d = max(list(cfg["scaling"]["d"]) + list(cfg["contraction"]["d"]))
        t = build_target("gaussian", int(d), float(kappa))
        md["targets"][f"gaussian__d{int(d)}__k{kappa:g}"] = t.metadata()
    for spec in cfg["wasserstein"]["targets"]:
        t = build_target(spec["name"], int(spec["d"]), float(spec["kappa"]),
                         gamma=float(spec.get("gamma", 0.0)),
                         n_nodes=int(spec.get("gh_nodes", 80)))
        md["targets"][f"wboot__{spec['name']}__d{int(spec['d'])}__k{spec['kappa']:g}"] = \
            t.metadata()
    for spec in cfg["stl_floor"]["targets"]:
        t = build_target(spec["name"], int(spec["d"]), float(spec["kappa"]),
                         gamma=float(spec.get("gamma", 0.0)),
                         n_nodes=int(spec.get("gh_nodes", 80)))
        md["targets"][f"stl__{spec['name']}__d{int(spec['d'])}__k{spec['kappa']:g}"] = \
            t.metadata()
    save_json(os.path.join(outdir, "target_metadata.json"), md)


def _exists_all(outdir, files):
    return all(os.path.exists(os.path.join(outdir, f)) for f in files)


def run_deterministic(cfg, outdir):
    long1, scaling = run_scaling(cfg, progress=_progress)
    long2, bench = run_contraction(cfg, progress=_progress)
    long3, wboot = run_wasserstein(cfg, progress=_progress)
    save_dataframe(os.path.join(outdir, "results_long.csv"),
                   pd.DataFrame(long1 + long2 + long3), columns=LONG_COLS)
    save_dataframe(os.path.join(outdir, "covariance_bootstrap_summary.csv"),
                   pd.DataFrame(scaling), columns=SCALING_SUMMARY_COLS)
    save_dataframe(os.path.join(outdir, "contraction_benchmark.csv"),
                   pd.DataFrame(bench), columns=CONTRACTION_COLS)
    save_dataframe(os.path.join(outdir, "wasserstein_bootstrap_summary.csv"),
                   pd.DataFrame(wboot), columns=WBOOT_SUMMARY_COLS)


def run_stochastic(cfg, outdir, device, backend_name):
    floor_rows, _ = run_stl_floor(cfg, device=device, backend_name=backend_name,
                                  progress=_progress)
    save_dataframe(os.path.join(outdir, "stl_floor_summary.csv"),
                   pd.DataFrame(floor_rows), columns=STL_FLOOR_COLS)


def parse_args():
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=DEFAULT_CONFIG)
    known, _ = pre.parse_known_args()
    cfg = load_yaml(known.config)
    p = argparse.ArgumentParser(parents=[pre], description=__doc__)
    p.add_argument("--outdir", default=cfg.get("output_dir",
                   "outputs/natural_gradient_covariance_bootstrap"))
    p.add_argument("--device", default="cpu",
                   help="cpu | cuda | cuda:N (stochastic STL only; cuda needs torch)")
    p.add_argument("--backend", default="auto", choices=["auto", "numpy", "torch"])
    p.add_argument("--smoke", action="store_true", help="fast reduced grids")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args(), cfg, known.config


if __name__ == "__main__":
    args, cfg, cfg_path = parse_args()
    if args.smoke:
        cfg = _apply_smoke(cfg)
    # Stochastic backend resolution (deterministic experiments are always NumPy/CPU).
    backend_name = resolve_backend(args.backend, args.device)
    cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES")

    print("=" * 68)
    print("Covariance bootstrap: enhanced global rates for the FR flow")
    print("=" * 68)
    print(f"  config : {cfg_path}{'  (smoke)' if args.smoke else ''}")
    print(f"  stochastic backend/device : {backend_name} / {args.device}")
    print(f"  CUDA_VISIBLE_DEVICES : {cuda_visible}")
    print(f"  outdir : {args.outdir}\n")

    ensure_dir(args.outdir)
    if args.overwrite:
        for name in OWN_FILES:
            pth = os.path.join(args.outdir, name)
            if os.path.exists(pth):
                os.remove(pth)

    t_start = time.time()
    _write_target_metadata(cfg, args.outdir)

    if args.overwrite or not _exists_all(args.outdir, DET_FILES):
        print("[deterministic] experiments 1-3 (CPU/NumPy)")
        run_deterministic(cfg, args.outdir)
    else:
        print("[deterministic] outputs present -> skip (use --overwrite to rerun)")

    if args.overwrite or not _exists_all(args.outdir, STL_FILES):
        print("[stochastic] experiment 4 (STL noise floor)")
        run_stochastic(cfg, args.outdir, args.device, backend_name)
    else:
        print("[stochastic] outputs present -> skip (use --overwrite to rerun)")

    save_json(os.path.join(args.outdir, "run_metadata.json"), {
        "experiment_group": "natural_gradient_covariance_bootstrap",
        "config_path": cfg_path, "smoke": bool(args.smoke), "config": cfg,
        "deterministic_backend": "numpy", "deterministic_device": "cpu",
        "stochastic_backend": backend_name, "stochastic_device": args.device,
        "cuda_visible_devices": cuda_visible,
        "torch_info": torch_device_info(args.device if backend_name == "torch" else "cpu"),
        "wall_time_total": float(time.time() - t_start),
        "python": platform.python_version(), "platform": platform.platform(),
        "numpy": np.__version__,
    })
    print(f"\nWrote {', '.join(OWN_FILES)} -> {args.outdir}")
