"""
Run affine-invariant gradient flow on the log-cosh non-Gaussian target.

Target: pi(x) ∝ exp(-V_rho(x))
  V_rho(x) = 0.5||x||^2 + (rho/m) sum_ell log cosh(a_ell^T x)

Outputs written to <outdir>/:
    results_long.csv          — one row per saved step per run
    summary.csv               — one row per run (convergence statistics)
    reference_optimum.npz     — cached reference Gaussian VI optimum
    reference_optimum_meta.json
    target_metadata.json      — target parameters + sample seed

Usage:
    python scripts/run_logconcave_grid.py [options]

Options:
    --n             dimension                 (default: 5)
    --rho           coupling strength         (default: 5.0)
    --K             MC sample count           (default: 4096)
    --K-ref         reference optimum samples (default: 8192)
    --dt            time step                 (default: 0.005)
    --T             total time                (default: 40.0)
    --omega         omega values              (default: 0.25 0.5 1.0)
    --target-seed   seed for A matrix         (default: 123)
    --sample-seed   seed for QMC samples      (default: 2026)
    --outdir        output directory
    --force-optimize  recompute reference optimum even if cached
    --save-every    save every N steps        (default: 10)
"""
import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np
import scipy.linalg
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.targets          import LogCoshTarget
from src.qmc_samples      import make_samples, push_forward
from src.reference_optimum import load_or_compute
from src.lc_dynamics      import logconcave_step
from src.lc_initializations import get_logconcave_initialization, LC_INIT_NAMES
from src.lc_metrics       import compute_lc_metrics, compute_lc_objective
from src.utils            import spd_invsqrt, spd_sqrt, validate_params, make_q_vector


# ---------------------------------------------------------------------------
# Tau grid
# ---------------------------------------------------------------------------

def make_tau_grid(omega, n):
    return [
        ("negative", -omega / (2.0 * n)),
        ("zero",      0.0),
        ("positive", +omega / (2.0 * n)),
    ]


# ---------------------------------------------------------------------------
# Expectation helpers (batch MC)
# ---------------------------------------------------------------------------

def mc_expectations(m, C, Z, target):
    """Return (g, S, obj_val) using fixed samples Z.

    g        = mean_j grad V(theta_j)        shape (n,)
    S        = mean_j Hess V(theta_j)        shape (n, n)
    obj_val  = mean_j V(theta_j) - 0.5 logdet C

    Uses C = L L^T via Cholesky for push-forward.
    """
    n = C.shape[0]
    # Cholesky factor (lower-triangular)
    try:
        L = scipy.linalg.cholesky(C, lower=True)
    except scipy.linalg.LinAlgError:
        # Fallback: symmetrize and add small jitter
        C_safe = 0.5 * (C + C.T) + 1e-10 * np.eye(n)
        L = scipy.linalg.cholesky(C_safe, lower=True)

    Theta = push_forward(m, L, Z)                     # (K, n)
    g     = np.mean(target.batch_grad(Theta), axis=0) # (n,)
    S     = np.mean(target.batch_hess(Theta), axis=0) # (n, n)
    S     = 0.5 * (S + S.T)

    V_mean = float(np.mean(target.batch_value(Theta)))
    log_det_L = float(np.sum(np.log(np.diag(L))))
    obj_val = V_mean - log_det_L                      # = V_mean - 0.5 logdet C

    return g, S, obj_val


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

LONG_COLS = [
    "n", "rho", "m_features", "target_seed", "sample_seed", "K",
    "omega", "tau_type", "tau_value", "init_name", "dt", "T",
    "step", "time",
    "objective", "objective_gap", "normalized_objective_gap",
    "whitened_mean_error", "cov_error", "volume_error", "shape_error",
    "mean_residual", "cov_residual", "trace_residual", "traceless_residual",
    "chi", "eig_min", "eig_max", "cosine_error_to_star",
]

SUMMARY_COLS = [
    "n", "rho", "m_features", "target_seed", "sample_seed", "K",
    "omega", "tau_type", "tau_value", "init_name", "dt", "T",
    "final_objective_gap", "final_normalized_objective_gap",
    "time_to_1e_minus_2", "time_to_1e_minus_4", "time_to_1e_minus_6",
    "monotone_objective_bool",
    "min_eig_min_over_time", "max_eig_max_over_time",
    "initial_chi", "final_chi",
    "initial_volume_error", "initial_shape_error",
    "final_volume_error",   "final_shape_error",
]


def _time_to_tol(records, tol):
    for r in records:
        if r["normalized_objective_gap"] <= tol:
            return r["time"]
    return math.inf


def _open_csv(path, cols):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fh = open(path, "w", newline="")
    w  = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    return fh, w


# ---------------------------------------------------------------------------
# Single run
# ---------------------------------------------------------------------------

def run_single(run_meta, m0, C0, ref, Z, target, dt, T, save_every, q, b, F_star, gap0):
    """Run one gradient flow experiment.

    Returns list of per-step metric dicts (with run_meta merged in).
    """
    n = target.n
    m, C = m0.copy(), C0.copy()

    # Pre-compute C_star_invsqrt once
    C_star_invsqrt = spd_invsqrt(ref["C_star"])
    m_star = ref["m_star"]
    C_star = ref["C_star"]

    omega = run_meta["omega"]
    tau   = run_meta["tau_value"]

    num_steps = int(round(T / dt))
    records   = []

    # Step 0 ----------------------------------------------------------------
    g0, S0, obj0 = mc_expectations(m, C, Z, target)
    met0 = compute_lc_metrics(m, C, g0, S0, m_star, C_star, C_star_invsqrt,
                               obj0, F_star, gap0, q, b)
    met0["step"] = 0
    met0["time"] = 0.0
    records.append({**run_meta, **met0})

    for step in range(1, num_steps + 1):
        # MC expectations at current state
        g, S, obj_val = mc_expectations(m, C, Z, target)
        # Dynamics step
        m, C = logconcave_step(m, C, g, S, dt, omega, tau)

        if step % save_every == 0 or step == num_steps:
            g_s, S_s, obj_s = mc_expectations(m, C, Z, target)
            met = compute_lc_metrics(m, C, g_s, S_s, m_star, C_star,
                                     C_star_invsqrt, obj_s, F_star, gap0, q, b)
            met["step"] = step
            met["time"] = step * dt
            records.append({**run_meta, **met})

    return records


def compute_summary(records):
    """Summarise a list of per-step records into one summary row."""
    norms = [r["normalized_objective_gap"] for r in records]
    monotone = all(norms[i+1] <= norms[i] + 1e-10 for i in range(len(norms)-1))

    first = records[0]
    last  = records[-1]

    row = {k: first[k] for k in [
        "n","rho","m_features","target_seed","sample_seed","K",
        "omega","tau_type","tau_value","init_name","dt","T"
    ]}
    row.update({
        "final_objective_gap":           last["objective_gap"],
        "final_normalized_objective_gap":last["normalized_objective_gap"],
        "time_to_1e_minus_2":            _time_to_tol(records, 1e-2),
        "time_to_1e_minus_4":            _time_to_tol(records, 1e-4),
        "time_to_1e_minus_6":            _time_to_tol(records, 1e-6),
        "monotone_objective_bool":       monotone,
        "min_eig_min_over_time":         min(r["eig_min"] for r in records),
        "max_eig_max_over_time":         max(r["eig_max"] for r in records),
        "initial_chi":                   first["chi"],
        "final_chi":                     last["chi"],
        "initial_volume_error":          first["volume_error"],
        "initial_shape_error":           first["shape_error"],
        "final_volume_error":            last["volume_error"],
        "final_shape_error":             last["shape_error"],
    })
    return row


# ---------------------------------------------------------------------------
# Grid runner
# ---------------------------------------------------------------------------

def run_grid(n, rho, omega_grid, dt, T, K, K_ref, target_seed, sample_seed,
             save_every, outdir, force_optimize):
    os.makedirs(outdir, exist_ok=True)

    # Build target and samples -----------------------------------------------
    target = LogCoshTarget(n=n, rho=rho, seed=target_seed)
    Z_dyn  = make_samples(n, K,     seed=sample_seed)       # dynamics samples
    Z_ref  = make_samples(n, K_ref, seed=sample_seed + 1)   # reference samples

    print(f"  Target: n={n}, rho={rho}, m={target.m_features}, target_seed={target_seed}")
    print(f"  Samples: K={K} (dynamics), K_ref={K_ref} (reference), sample_seed={sample_seed}")

    # Reference optimum -------------------------------------------------------
    ref_path = os.path.join(outdir, "reference_optimum.npz")
    ref = load_or_compute(ref_path, target, Z_ref, force=force_optimize)
    F_star = float(ref["F_star"])
    m_star = ref["m_star"]
    C_star = ref["C_star"]

    print(f"  F_star={F_star:.6f}  ||m_star||={np.linalg.norm(m_star):.3e}"
          f"  grad_m_norm={ref.get('grad_m_norm', float('nan')):.3e}")

    # Save target metadata
    meta = {
        "n": n, "rho": rho,
        "m_features": target.m_features,
        "target_seed": target_seed,
        "sample_seed": sample_seed,
        "K": K, "K_ref": K_ref,
        "dt": dt, "T": T,
        "F_star": F_star,
        "m_star_norm": float(np.linalg.norm(m_star)),
    }
    with open(os.path.join(outdir, "target_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)

    # Test-function vector
    q = make_q_vector(n)
    b = 0.5

    # Open CSVs
    long_path    = os.path.join(outdir, "results_long.csv")
    summary_path = os.path.join(outdir, "summary.csv")
    long_fh,    long_w    = _open_csv(long_path,    LONG_COLS)
    summary_fh, summary_w = _open_csv(summary_path, SUMMARY_COLS)

    total = len(omega_grid) * 3 * len(LC_INIT_NAMES)
    done  = 0
    t0    = time.time()

    for omega in omega_grid:
        for tau_type, tau in make_tau_grid(omega, n):
            if not validate_params(omega, tau, n):
                print(f"  [SKIP] invalid params omega={omega} tau={tau}")
                continue

            for init_name in LC_INIT_NAMES:
                done += 1

                m0, C0 = get_logconcave_initialization(init_name, n, m_star, C_star)

                # Initial objective gap for normalisation
                try:
                    L0 = scipy.linalg.cholesky(C0, lower=True)
                except scipy.linalg.LinAlgError:
                    L0 = np.eye(n)
                obj0 = compute_lc_objective(m0, L0, Z_dyn, target)
                gap0 = obj0 - F_star
                if abs(gap0) < 1e-15:
                    gap0 = 1.0  # avoid division by zero

                run_meta = {
                    "n": n, "rho": rho,
                    "m_features": target.m_features,
                    "target_seed": target_seed,
                    "sample_seed": sample_seed,
                    "K": K,
                    "omega": omega, "tau_type": tau_type,
                    "tau_value": tau,
                    "init_name": init_name,
                    "dt": dt, "T": T,
                }

                records = run_single(
                    run_meta, m0, C0, ref,
                    Z_dyn, target, dt, T, save_every, q, b, F_star, gap0,
                )

                for rec in records:
                    long_w.writerow({k: rec.get(k, "") for k in LONG_COLS})

                summ = compute_summary(records)
                summary_w.writerow({k: summ.get(k, "") for k in SUMMARY_COLS})

                elapsed = time.time() - t0
                print(
                    f"  [{done:3d}/{total}] omega={omega:.3g} tau={tau_type:8s} "
                    f"init={init_name:12s}  ({elapsed:.1f}s)"
                )

    long_fh.close()
    summary_fh.close()
    print(f"\nLong CSV   -> {long_path}")
    print(f"Summary    -> {summary_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    cfg_path = os.path.join(os.path.dirname(__file__), "..", "configs", "logconcave_target.yaml")
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)
    d = cfg["defaults"]

    p = argparse.ArgumentParser(description="Run log-concave gradient flow grid.")
    p.add_argument("--n",             type=int,   default=d["n"])
    p.add_argument("--rho",           type=float, default=d["rho"])
    p.add_argument("--K",             type=int,   default=d["K"])
    p.add_argument("--K-ref",         type=int,   default=d["K_ref"],  dest="K_ref")
    p.add_argument("--dt",            type=float, default=d["dt"])
    p.add_argument("--T",             type=float, default=d["T"])
    p.add_argument("--omega",         type=float, nargs="+",
                   default=cfg["omega_grid"])
    p.add_argument("--target-seed",   type=int,   default=d["target_seed"], dest="target_seed")
    p.add_argument("--sample-seed",   type=int,   default=d["sample_seed"], dest="sample_seed")
    p.add_argument("--save-every",    type=int,   default=d["save_every"],  dest="save_every")
    p.add_argument("--outdir",        type=str,   default=cfg["output_dir"])
    p.add_argument("--force-optimize",action="store_true", dest="force_optimize")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    print("=" * 60)
    print("Affine-Invariant Gradient Flow — Log-Concave Target")
    print("=" * 60)
    print(f"  n={args.n}  rho={args.rho}  K={args.K}  K_ref={args.K_ref}")
    print(f"  dt={args.dt}  T={args.T}  save_every={args.save_every}")
    print(f"  omega grid: {args.omega}")
    print(f"  outdir: {args.outdir}")
    print()

    run_grid(
        n=args.n, rho=args.rho,
        omega_grid=args.omega,
        dt=args.dt, T=args.T,
        K=args.K, K_ref=args.K_ref,
        target_seed=args.target_seed,
        sample_seed=args.sample_seed,
        save_every=args.save_every,
        outdir=args.outdir,
        force_optimize=args.force_optimize,
    )
