"""Grid enumeration and cell execution for the STL algorithm experiment.

Kept in the importable package (not in a script's ``__main__``) so that the
optional two-GPU path can re-import the worker under the multiprocessing ``spawn``
start method. ``run_cells`` builds one backend and a target cache and executes a
list of algorithm cells; the two-GPU driver simply partitions the cell list
across two worker processes, each pinned to one CUDA device.

Cell seeds are derived deterministically from the cell identity and the config
``base_seed`` only (not from the device or the method), so a baseline/STL pair in
the same ``(target, batch_size)`` cell shares its common-random-number stream and
results are identical whether the grid runs on one device or two.
"""
from __future__ import annotations

import hashlib

from src.natural_gradient_stl_variance.targets import build_target
from src.natural_gradient_stl_variance.estimators import METHODS
from src.natural_gradient_stl_variance.algorithm import simulate_cell
from src.natural_gradient_stl_variance.linalg import ArrayBackend


def derive_seed(*parts):
    """Stable 60-bit integer seed from ``parts`` (reproducible across processes)."""
    s = "|".join(str(p) for p in parts)
    return int(hashlib.sha256(s.encode()).hexdigest()[:15], 16)


def enumerate_target_configs(cfg):
    """List ``(kind, d, kappa, tau, gh_nodes)`` target configs from a config dict."""
    out = []
    g = cfg.get("gaussian", {})
    for d in g.get("d", []):
        for kappa in g.get("kappa", []):
            out.append(("gaussian", int(d), float(kappa), 0.0, 0))
    lc = cfg.get("log_cosh", {})
    gh = int(lc.get("gh_nodes", 80))
    for d in lc.get("d", []):
        for kappa in lc.get("kappa", []):
            for tau in lc.get("tau", []):
                out.append(("log_cosh", int(d), float(kappa), float(tau), gh))
    return out


def build_algorithm_cells(cfg):
    """Build the full list of algorithm cells (one per config/method/batch size)."""
    alg = cfg["algorithm"]
    base_seed = int(cfg.get("base_seed", 0))
    batch_sizes = [int(b) for b in alg["batch_sizes"]]
    cells = []
    for kind, d, kappa, tau, gh in enumerate_target_configs(cfg):
        n_steps = int(alg["gaussian_n_steps"] if kind == "gaussian"
                      else alg["log_cosh_n_steps"])
        for bs in batch_sizes:
            cell_seed = derive_seed(base_seed, kind, d, kappa, tau, bs)
            for method in METHODS:
                cells.append({
                    "kind": kind, "d": d, "kappa": kappa, "tau": tau,
                    "target_gh_nodes": gh, "method": method, "batch_size": bs,
                    "n_steps": n_steps, "dt": float(alg["dt"]),
                    "init_mean_rho": float(alg["init_mean_rho"]),
                    "init_cov_scale": float(alg["init_cov_scale"]),
                    "n_saved": int(alg["n_saved"]),
                    "tail_frac": float(alg["tail_frac"]),
                    "alg_gh_nodes": int(alg.get("gh_nodes", 40)),
                    "cell_seed": int(cell_seed),
                })
    return cells


def run_cells(cells, seeds, device="cpu", backend_name="numpy", progress=None):
    """Execute ``cells`` on one backend/device; return accumulated record lists.

    ``progress`` (optional) is a callable ``(i, n, cell, diag)`` invoked after each
    cell for logging. Targets are cached across methods/batch sizes.
    """
    bk = ArrayBackend(backend=backend_name, device=device, dtype="float64")
    seeds = [int(s) for s in seeds]
    target_cache = {}
    long_rows, summary_rows, tail_rows, diag_rows = [], [], [], []
    n = len(cells)
    for i, cell in enumerate(cells):
        key = (cell["kind"], cell["d"], cell["kappa"], cell["tau"])
        target = target_cache.get(key)
        if target is None:
            target = build_target(cell["kind"], cell["d"], cell["kappa"],
                                  tau=cell["tau"], n_nodes=cell["target_gh_nodes"] or 80)
            target_cache[key] = target
        lr, sr, tr, diag = simulate_cell(
            target, cell["method"], cell["dt"], cell["n_steps"],
            cell["batch_size"], seeds, bk, cell["cell_seed"],
            n_saved=cell["n_saved"], tail_frac=cell["tail_frac"],
            init_mean_rho=cell["init_mean_rho"],
            init_cov_scale=cell["init_cov_scale"], gh_nodes=cell["alg_gh_nodes"])
        long_rows.extend(lr)
        summary_rows.extend(sr)
        tail_rows.extend(tr)
        diag_row = {"kind": cell["kind"], "d": cell["d"], "kappa": cell["kappa"],
                    "tau": cell["tau"], "batch_size": cell["batch_size"]}
        diag_row.update(diag)
        diag_rows.append(diag_row)
        if progress is not None:
            progress(i + 1, n, cell, diag)
    return long_rows, summary_rows, tail_rows, diag_rows


def gpu_shard_worker(cells, seeds, gpu_id, shard_paths):
    """Worker entry point for the two-GPU path: run a shard and write its CSVs.

    Pins the process to one physical CUDA device via ``CUDA_VISIBLE_DEVICES``
    *before* torch is imported (the backend imports torch lazily), then writes the
    shard's long/summary/tail CSVs to ``shard_paths``.
    """
    import os
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    import pandas as pd
    from src.common.io_utils import save_dataframe

    def _progress(i, n, cell, diag):
        print(f"  [gpu{gpu_id}] [{i:3d}/{n}] {cell['kind']:8s} d={cell['d']:<2d} "
              f"kappa={cell['kappa']:<8g} tau={cell['tau']:<4g} "
              f"{cell['method']:14s} B={cell['batch_size']:<2d} "
              f"clips={diag['n_clips']} ({diag['wall_time_cell']:.1f}s)", flush=True)

    lr, sr, tr, _ = run_cells(cells, seeds, device="cuda",
                              backend_name="torch", progress=_progress)
    save_dataframe(shard_paths["long"], pd.DataFrame(lr))
    save_dataframe(shard_paths["summary"], pd.DataFrame(sr))
    save_dataframe(shard_paths["tail"], pd.DataFrame(tr))
