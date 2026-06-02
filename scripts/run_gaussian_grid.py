"""
Run affine-invariant Gaussian gradient flow on a full parameter grid.

Target: N(0, I_n).  All expectations are exact (closed-form update steps).

Outputs written to <outdir>/:
    results_long.csv   — one row per saved step per run
    summary.csv        — one row per run (convergence statistics)

Usage:
    python scripts/run_gaussian_grid.py [options]

Options:
    --dt DT          time step          (default: 0.02)
    --T  T           total time         (default: 20.0)
    --n  N [N ...]   dimension list     (default: 2 5 10)
    --save-every K   save every K steps (default: 5)
    --outdir DIR     output directory   (default: outputs/gaussian_grid)
"""
import argparse
import csv
import math
import os
import sys
import time

import numpy as np
import yaml

# Allow running from project root without installing the package
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.dynamics import gaussian_step
from src.metrics import compute_all_metrics, kl_energy
from src.initializations import get_initialization, INIT_NAMES
from src.utils import validate_params, make_q_vector


# ---------------------------------------------------------------------------
# Grid specification
# ---------------------------------------------------------------------------

OMEGA_GRID = [1 / 8, 1 / 4, 1 / 2, 1.0, 2.0]


def make_tau_grid(omega, n):
    """Return list of (tau_type, tau_value) for the given omega and n."""
    return [
        ("negative", -omega / (2.0 * n)),
        ("zero",      0.0),
        ("positive", +omega / (2.0 * n)),
    ]


# ---------------------------------------------------------------------------
# Summary statistics from a list of per-step metric dicts
# ---------------------------------------------------------------------------

TOLERANCES = [1e-2, 1e-4, 1e-6]


def _time_to_tol(records, tol):
    """First time (float) at which norm_energy <= tol, or inf if never."""
    for r in records:
        if r["norm_energy"] <= tol:
            return r["time"]
    return math.inf


def compute_summary(run_meta, records):
    """Aggregate per-step records into one summary dict."""
    norm_energies = [r["norm_energy"] for r in records]
    eig_mins = [r["eig_min"] for r in records]
    eig_maxs = [r["eig_max"] for r in records]
    last = records[-1]

    # Monotone energy: every consecutive pair is non-increasing
    monotone = all(
        norm_energies[i + 1] <= norm_energies[i] + 1e-12
        for i in range(len(norm_energies) - 1)
    )

    row = dict(run_meta)
    row.update({
        "final_energy":            last["kl_energy"],
        "final_normalized_energy": last["norm_energy"],
        "time_to_1e_minus_2":      _time_to_tol(records, 1e-2),
        "time_to_1e_minus_4":      _time_to_tol(records, 1e-4),
        "time_to_1e_minus_6":      _time_to_tol(records, 1e-6),
        "monotone_energy_bool":    monotone,
        "min_eig_min_over_time":   min(eig_mins),
        "max_eig_max_over_time":   max(eig_maxs),
    })
    return row


# ---------------------------------------------------------------------------
# Single-run experiment
# ---------------------------------------------------------------------------

def run_single(n, omega, tau_type, tau, init_name, dt, T, save_every):
    """Run one experiment; return (run_meta, list_of_step_dicts)."""
    if not validate_params(omega, tau, n):
        return None, []

    m, C = get_initialization(init_name, n)
    q = make_q_vector(n)
    b = 0.5

    # Initial energy (used to normalise; avoid division by zero)
    E0 = kl_energy(m, C)
    if E0 < 1e-15:
        E0 = 1.0

    run_meta = {
        "n": n, "omega": omega,
        "tau_type": tau_type, "tau_value": tau,
        "init_name": init_name,
        "dt": dt, "T": T,
    }

    num_steps = int(round(T / dt))
    records = []

    # Step 0
    metrics = compute_all_metrics(m, C, E0, q, b)
    metrics["step"] = 0
    metrics["time"] = 0.0
    records.append(metrics)

    for step in range(1, num_steps + 1):
        m, C = gaussian_step(m, C, dt, omega, tau)

        if step % save_every == 0 or step == num_steps:
            metrics = compute_all_metrics(m, C, E0, q, b)
            metrics["step"] = step
            metrics["time"] = step * dt
            records.append(metrics)

    return run_meta, records


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

LONG_COLS = [
    "n", "omega", "tau_type", "tau_value", "init_name", "dt", "T",
    "step", "time",
    "kl_energy", "norm_energy", "mean_error", "cov_error",
    "volume_error", "shape_error", "cosine_error",
    "eig_min", "eig_max", "chi",
]

SUMMARY_COLS = [
    "n", "omega", "tau_type", "tau_value", "init_name", "dt", "T",
    "final_energy", "final_normalized_energy",
    "time_to_1e_minus_2", "time_to_1e_minus_4", "time_to_1e_minus_6",
    "monotone_energy_bool",
    "min_eig_min_over_time", "max_eig_max_over_time",
]


def _open_csv(path, cols):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fh = open(path, "w", newline="")
    writer = csv.DictWriter(fh, fieldnames=cols)
    writer.writeheader()
    return fh, writer


# ---------------------------------------------------------------------------
# Main grid runner
# ---------------------------------------------------------------------------

def run_grid(n_list, omega_grid, dt, T, save_every, outdir):
    os.makedirs(outdir, exist_ok=True)

    long_path    = os.path.join(outdir, "results_long.csv")
    summary_path = os.path.join(outdir, "summary.csv")

    long_fh,    long_writer    = _open_csv(long_path,    LONG_COLS)
    summary_fh, summary_writer = _open_csv(summary_path, SUMMARY_COLS)

    # Count total runs for progress display
    total = len(n_list) * len(omega_grid) * 3 * len(INIT_NAMES)
    done  = 0
    t0    = time.time()

    for n in n_list:
        for omega in omega_grid:
            for tau_type, tau in make_tau_grid(omega, n):
                for init_name in INIT_NAMES:
                    done += 1
                    run_meta, records = run_single(
                        n, omega, tau_type, tau,
                        init_name, dt, T, save_every,
                    )
                    if not records:
                        print(f"  [SKIP] n={n} omega={omega} tau={tau_type}: invalid params")
                        continue

                    # Write long-format rows
                    for rec in records:
                        row = dict(run_meta)
                        row.update(rec)
                        long_writer.writerow({k: row[k] for k in LONG_COLS})

                    # Write summary row
                    summ = compute_summary(run_meta, records)
                    summary_writer.writerow({k: summ[k] for k in SUMMARY_COLS})

                    elapsed = time.time() - t0
                    print(
                        f"  [{done:3d}/{total}] n={n} omega={omega:.4g} "
                        f"tau={tau_type:8s} init={init_name:12s}  "
                        f"({elapsed:.1f}s)"
                    )

    long_fh.close()
    summary_fh.close()
    print(f"\nLong CSV   -> {long_path}")
    print(f"Summary    -> {summary_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    cfg_path = os.path.join(os.path.dirname(__file__), "..", "configs", "gaussian_target.yaml")
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    parser = argparse.ArgumentParser(description="Run Gaussian gradient flow grid.")
    parser.add_argument("--dt",         type=float, default=cfg["defaults"]["dt"])
    parser.add_argument("--T",          type=float, default=cfg["defaults"]["T"])
    parser.add_argument("--n",          type=int,   nargs="+",
                        default=cfg["dimensions"])
    parser.add_argument("--save-every", type=int,   default=cfg["defaults"]["save_every"],
                        dest="save_every")
    parser.add_argument("--outdir",     type=str,   default=cfg["output_dir"])
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    print("=" * 60)
    print("Affine-Invariant Gaussian Gradient Flow — Grid Experiment")
    print("=" * 60)
    print(f"  n list    : {args.n}")
    print(f"  omega grid: {OMEGA_GRID}")
    print(f"  dt={args.dt}  T={args.T}  save_every={args.save_every}")
    print(f"  outdir    : {args.outdir}")
    print()

    run_grid(
        n_list=args.n,
        omega_grid=OMEGA_GRID,
        dt=args.dt,
        T=args.T,
        save_every=args.save_every,
        outdir=args.outdir,
    )
