"""Experiment B -- algorithm-level STL trajectories.

Runs the four methods ``{riemannian, riemannian_stl, kl, kl_stl}`` from a fixed
initialization on every target config and mini-batch size, with all seeds of a
cell advanced in one batched pass on the chosen backend. The baseline and STL
methods of a ``(config, batch_size)`` cell share their common-random-number
stream, so the comparison is paired.

Writes into ``--outdir``::

    algorithm_results_long.csv   one row per saved step per seed
    algorithm_summary.csv        one row per (config, method, batch size, seed)
    tail_noise_floor.csv         tail / noise-floor statistics, per run
    target_metadata.json         target params, a_star, initialization
    run_metadata.json            config echo, timing, environment

By default the grid runs on a single device (``--device cpu`` or ``--device
cuda``). ``--max-gpus 2`` partitions the cell list across two CUDA devices (the
hard cap is two); it never launches more than two GPU jobs.

Usage::

    python scripts/natural_gradient_stl_variance/run_algorithm_grid.py \
        --config configs/natural_gradient_stl_variance/stl_variance.yaml \
        --outdir outputs/natural_gradient_stl_variance --device cuda --max-gpus 2 --overwrite
"""
from __future__ import annotations

import argparse
import os
import platform
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.io_utils import load_yaml, ensure_dir, save_dataframe, save_json
from src.common.torch_utils import torch_device_info
from src.natural_gradient_stl_variance.targets import build_target
from src.natural_gradient_stl_variance.metrics import GAP_THRESHOLDS, _tol_key
from src.natural_gradient_stl_variance.grid import (
    enumerate_target_configs, build_algorithm_cells, run_cells, gpu_shard_worker,
)

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "..", "configs",
    "natural_gradient_stl_variance", "stl_variance.yaml")

_HIT_COLS = [f"iter_to_{_tol_key(t)}" for t in GAP_THRESHOLDS]

LONG_COLS = [
    "target_name", "kind", "method", "scheme", "stl", "seed", "d", "kappa", "tau",
    "dt", "batch_size", "step", "energy_gap", "sq_mean_error", "rel_cov_fro_error",
    "w2_sq", "min_eig_C", "max_eig_C", "spd_fail",
]
SUMMARY_COLS = [
    "target_name", "kind", "method", "scheme", "stl", "seed", "d", "kappa", "tau",
    "dt", "batch_size", "n_steps", "gap0", "gap_final", "gap_min", "min_eig_C_min",
    "spd_fail", "wall_time_cell", "n_seeds_in_cell",
] + _HIT_COLS + [
    "tail_mean_gap", "tail_median_gap", "tail_std_gap", "final_gap", "tail_steps",
]
TAIL_COLS = [
    "target_name", "kind", "method", "scheme", "stl", "seed", "d", "kappa", "tau",
    "dt", "batch_size", "n_steps",
    "tail_mean_gap", "tail_median_gap", "tail_std_gap", "final_gap", "tail_steps",
    "spd_fail",
]

# Files this script owns; --overwrite clears only these so the algorithm and the
# estimator runs can share an output directory (the README uses one outdir).
OWN_FILES = [
    "algorithm_results_long.csv", "algorithm_summary.csv", "tail_noise_floor.csv",
    "target_metadata.json", "run_metadata.json",
]


def _write_metadata(cfg, outdir):
    """Target metadata (a_star and how it was computed) + the initialization."""
    alg = cfg["algorithm"]
    metadata = {"targets": {}, "initialization": {
        "init_mean_rho": float(alg["init_mean_rho"]),
        "init_cov_scale": float(alg["init_cov_scale"]),
        "dt": float(alg["dt"]),
        "note": ("m0 = m_star + init_mean_rho * sqrt(diag(C_star)); "
                 "C0 = init_cov_scale * C_star"),
    }}
    for kind, d, kappa, tau, gh in enumerate_target_configs(cfg):
        target = build_target(kind, d, kappa, tau=tau, n_nodes=gh or 80)
        metadata["targets"][f"{kind}__d{d}__k{kappa:g}__t{tau:g}"] = target.metadata()
    save_json(os.path.join(outdir, "target_metadata.json"), metadata)


def _write_frames(outdir, long_rows, summary_rows, tail_rows):
    save_dataframe(os.path.join(outdir, "algorithm_results_long.csv"),
                   pd.DataFrame(long_rows), columns=LONG_COLS)
    save_dataframe(os.path.join(outdir, "algorithm_summary.csv"),
                   pd.DataFrame(summary_rows), columns=SUMMARY_COLS)
    save_dataframe(os.path.join(outdir, "tail_noise_floor.csv"),
                   pd.DataFrame(tail_rows), columns=TAIL_COLS)


def _progress(i, n, cell, diag):
    print(f"  [{i:3d}/{n}] {cell['kind']:8s} d={cell['d']:<2d} "
          f"kappa={cell['kappa']:<8g} tau={cell['tau']:<4g} {cell['method']:14s} "
          f"B={cell['batch_size']:<2d} clips={diag['n_clips']} "
          f"({diag['wall_time_cell']:.1f}s)", flush=True)


def run_single(cfg, outdir, device, backend_name):
    cells = build_algorithm_cells(cfg)
    seeds = [int(s) for s in cfg["seeds"]]
    print(f"  {len(cells)} cells x {len(seeds)} seeds on device={device} "
          f"backend={backend_name}\n")
    lr, sr, tr, _ = run_cells(cells, seeds, device=device,
                              backend_name=backend_name, progress=_progress)
    _write_frames(outdir, lr, sr, tr)


def run_two_gpu(cfg, outdir, n_workers):
    """Partition the cell list across ``n_workers`` (<=2) CUDA devices."""
    import multiprocessing as mp

    cells = build_algorithm_cells(cfg)
    seeds = [int(s) for s in cfg["seeds"]]
    shards = [cells[w::n_workers] for w in range(n_workers)]
    print(f"  {len(cells)} cells x {len(seeds)} seeds across {n_workers} GPUs\n")

    ctx = mp.get_context("spawn")
    procs, shard_files = [], []
    for w in range(n_workers):
        paths = {k: os.path.join(outdir, f"_shard{w}_{k}.csv")
                 for k in ("long", "summary", "tail")}
        shard_files.append(paths)
        p = ctx.Process(target=gpu_shard_worker,
                        args=(shards[w], seeds, w, paths))
        p.start()
        procs.append(p)
    for p in procs:
        p.join()
        if p.exitcode != 0:
            raise RuntimeError(f"GPU shard worker exited with code {p.exitcode}")

    # Merge shard CSVs into the final outputs, then remove the shards.
    merged = {k: [] for k in ("long", "summary", "tail")}
    for paths in shard_files:
        for k, path in paths.items():
            if os.path.exists(path):
                merged[k].append(pd.read_csv(path))
    long_df = pd.concat(merged["long"], ignore_index=True)
    summ_df = pd.concat(merged["summary"], ignore_index=True)
    tail_df = pd.concat(merged["tail"], ignore_index=True)
    save_dataframe(os.path.join(outdir, "algorithm_results_long.csv"),
                   long_df, columns=LONG_COLS)
    save_dataframe(os.path.join(outdir, "algorithm_summary.csv"),
                   summ_df, columns=SUMMARY_COLS)
    save_dataframe(os.path.join(outdir, "tail_noise_floor.csv"),
                   tail_df, columns=TAIL_COLS)
    for paths in shard_files:
        for path in paths.values():
            if os.path.exists(path):
                os.remove(path)


def _resolve_backend(device, backend):
    if backend != "auto":
        return backend
    return "torch" if str(device).startswith("cuda") else "numpy"


def _n_gpu_workers(device, max_gpus):
    """Number of GPU worker processes: capped at 2 and at the device count."""
    if not str(device).startswith("cuda"):
        return 1
    want = max(1, min(int(max_gpus), 2))      # hard cap at two GPUs
    if want <= 1:
        return 1
    try:
        from src.common.torch_utils import get_torch
        ndev = get_torch().cuda.device_count()
    except Exception:
        ndev = 1
    return min(want, ndev) if ndev >= 1 else 1


def parse_args():
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=DEFAULT_CONFIG)
    known, _ = pre.parse_known_args()
    cfg = load_yaml(known.config)
    p = argparse.ArgumentParser(parents=[pre], description=__doc__)
    p.add_argument("--outdir", default=cfg["output_dir"])
    p.add_argument("--device", default="cpu",
                   help="cpu | cuda | cuda:N (cuda requires torch + a CUDA device)")
    p.add_argument("--backend", default="auto", choices=["auto", "numpy", "torch"])
    p.add_argument("--max-gpus", type=int, default=1,
                   help="max concurrent CUDA devices for the grid (hard cap 2)")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args(), cfg, known.config


if __name__ == "__main__":
    args, cfg, cfg_path = parse_args()
    backend_name = _resolve_backend(args.device, args.backend)
    n_workers = _n_gpu_workers(args.device, args.max_gpus)
    print("=" * 64)
    print("STL algorithm-level trajectories (Experiment B)")
    print("=" * 64)
    print(f"  config : {cfg_path}")
    print(f"  device : {args.device}   backend: {backend_name}   gpu_workers: {n_workers}")
    print(f"  outdir : {args.outdir}\n")

    ensure_dir(args.outdir)
    if args.overwrite:
        for name in OWN_FILES:
            pth = os.path.join(args.outdir, name)
            if os.path.exists(pth):
                os.remove(pth)

    t_start = time.time()
    _write_metadata(cfg, args.outdir)
    if n_workers >= 2:
        run_two_gpu(cfg, args.outdir, n_workers)
    else:
        run_single(cfg, args.outdir, args.device, backend_name)

    save_json(os.path.join(args.outdir, "run_metadata.json"), {
        "experiment": "algorithm_grid", "config_path": cfg_path, "config": cfg,
        "device": args.device, "backend": backend_name, "gpu_workers": n_workers,
        "torch_info": torch_device_info(args.device if backend_name == "torch" else "cpu"),
        "wall_time_total": float(time.time() - t_start),
        "python": platform.python_version(), "platform": platform.platform(),
        "numpy": np.__version__,
    })
    print(f"\nWrote algorithm_results_long.csv, algorithm_summary.csv, "
          f"tail_noise_floor.csv, target_metadata.json, run_metadata.json "
          f"-> {args.outdir}")
