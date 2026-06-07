"""Run the Riemannian-vs-KL discretization stepsize grid (CPU, deterministic).

For every ``(target, lambda, method, dt)`` the runner simulates ``N = ceil(T/dt)``
steps from the target's initial condition, classifies the run (SPD-feasible /
stable / monotone / accurate), and writes three CSVs plus ``target_metadata.json``:

    results_long.csv      one row per saved step per run
    summary.csv           one row per run
    stepsize_summary.csv  one row per (target, lambda, method)
    target_metadata.json  target params, a_star, theory bounds, ODE diagnostics

Usage::

    python scripts/natural_gradient_discretization_stepsize/run_stepsize_grid.py \
        --config configs/natural_gradient_discretization_stepsize/stepsize_grid.yaml \
        --outdir outputs/natural_gradient_discretization_stepsize --overwrite

Add ``--smoke`` for the fast reduced grid defined in the config.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.io_utils import load_yaml, ensure_dir, save_dataframe, save_json
from src.natural_gradient_discretization_stepsize.targets import build_target
from src.natural_gradient_discretization_stepsize.optimize_star import compute_star
from src.natural_gradient_discretization_stepsize.ode_reference import integrate_reference
from src.natural_gradient_discretization_stepsize.metrics import (
    theory_stepsize_bounds, dt_theory_for_method,
)
from src.natural_gradient_discretization_stepsize.runner import (
    simulate_run, simulate_scalar_diagnostic, stepsize_summary_rows,
)

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "..", "configs",
    "natural_gradient_discretization_stepsize", "stepsize_grid.yaml")

LONG_COLS = [
    "target_name", "lambda", "method", "dt", "n", "t",
    "energy", "energy_gap", "mean_error_to_star", "cov_error_to_star",
    "ode_mean_error", "ode_cov_error", "min_eig_C", "max_eig_C",
    "is_spd", "energy_increase", "wall_time_cumulative",
]

SUMMARY_COLS = [
    "target_name", "lambda", "method", "dt", "T", "N", "F0", "F_final",
    "gap_final", "gap_min", "spd_feasible", "stable", "monotone", "accurate",
    "num_energy_increases", "max_energy_increase", "max_gap_ratio",
    "min_eig_C_min", "max_eig_C_max", "wall_time_total",
    "time_to_gap_1e_minus_4", "time_to_gap_1e_minus_6",
    "iter_to_gap_1e_minus_4", "iter_to_gap_1e_minus_6",
    "terminal_accuracy_error", "failed_at_step",
]

STEPSIZE_COLS = [
    "target_name", "lambda", "method", "theory_bound_available",
    "dt_theory_riem", "dt_theory_kl", "dt_theory_for_method",
    "dt_max_spd", "dt_max_stable", "dt_max_monotone", "dt_max_accurate",
    "stable_over_theory_ratio", "monotone_over_theory_ratio",
    "accurate_over_theory_ratio",
]

SCALAR_COLS = ["experiment", "method", "C0", "dt", "n", "t", "C", "m"]


def _apply_smoke(cfg):
    """Shallow-merge the optional ``smoke:`` block over the base config."""
    for k, v in cfg.get("smoke", {}).items():
        cfg[k] = v
    return cfg


def _build_target(name, lam, cfg):
    if name == "smooth_logconcave":
        st = cfg.get("smooth_target", {})
        return build_target(name, lam, delta=float(st.get("delta", 0.05)),
                            gamma=float(st.get("gamma", 1.0)),
                            n_nodes=int(st.get("gh_nodes", 80)))
    return build_target(name, lam)


# ---------------------------------------------------------------------------
# Main grid
# ---------------------------------------------------------------------------

def run_grid(cfg, outdir):
    ensure_dir(outdir)
    T = float(cfg["T"])
    dt_grid = [float(x) for x in cfg["dt_grid"]]
    ode_cfg = cfg.get("ode_reference", {})

    long_rows, summary_rows = [], []
    theory_by_key = {}
    metadata = {"config": {k: cfg[k] for k in
                           ("targets", "lambdas", "methods", "dt_grid", "T")},
                "targets": {}}

    total = len(cfg["targets"]) * len(cfg["lambdas"]) * len(cfg["methods"]) * len(dt_grid)
    done = 0
    t0 = time.time()

    for target_name in cfg["targets"]:
        for lam in cfg["lambdas"]:
            target = _build_target(target_name, lam, cfg)
            # Reference optimum and theory bounds (per target/lambda).
            m_star, C_star, F_star, opt_diag = compute_star(target)
            bounds = theory_stepsize_bounds(
                getattr(target, "alpha", None), getattr(target, "beta", None), target.C0)
            ode = integrate_reference(
                target, T,
                rtol=float(ode_cfg.get("rtol", 1e-10)),
                atol=float(ode_cfg.get("atol", 1e-12)),
                n_eval=int(ode_cfg.get("n_eval", 400)))

            meta_key = f"{target_name}__lam{lam:g}"
            md = target.metadata()
            md.update({
                "m_star": np.asarray(m_star).tolist(),
                "C_star": np.asarray(C_star).tolist(),
                "F_star": float(F_star),
                "optimizer": opt_diag,
                "theory_bounds": bounds,
                "ode_solver": {k: ode[k] for k in
                               ("success", "status", "message", "nfev", "njev",
                                "rtol", "atol", "method", "T")},
                "ode_terminal_m": np.asarray(ode["m_T"]).tolist(),
                "ode_terminal_C": np.asarray(ode["C_T"]).tolist(),
            })
            metadata["targets"][meta_key] = md

            for method in cfg["methods"]:
                key = (target_name, lam, method)
                tb = dict(bounds)
                tb["dt_theory_for_method"] = dt_theory_for_method(bounds, method)
                theory_by_key[key] = tb
                for dt in dt_grid:
                    records, summ = simulate_run(
                        method, target, dt, T, F_star, m_star, C_star, ode,
                        max_saved_rows=int(cfg.get("max_saved_rows", 400)))
                    for rec in records:
                        row = {"target_name": target_name, "lambda": lam,
                               "method": method, "dt": dt}
                        row.update(rec)
                        long_rows.append(row)
                    summary_rows.append(summ)
                    done += 1
                    flags = "".join(c for c, k in (
                        ("S", "spd_feasible"), ("T", "stable"),
                        ("M", "monotone"), ("A", "accurate")) if summ[k])
                    print(f"  [{done:4d}/{total}] {target_name:22s} "
                          f"lam={lam:<5g} {method:10s} dt={dt:<6g} "
                          f"[{flags or '-':4s}] ({time.time() - t0:5.1f}s)")

    step_rows = stepsize_summary_rows(summary_rows, theory_by_key)

    long_df = pd.DataFrame(long_rows).reindex(columns=LONG_COLS)
    summary_df = pd.DataFrame(summary_rows).reindex(columns=SUMMARY_COLS)
    step_df = pd.DataFrame(step_rows).reindex(columns=STEPSIZE_COLS)
    save_dataframe(os.path.join(outdir, "results_long.csv"), long_df)
    save_dataframe(os.path.join(outdir, "summary.csv"), summary_df)
    save_dataframe(os.path.join(outdir, "stepsize_summary.csv"), step_df)
    save_json(os.path.join(outdir, "target_metadata.json"), metadata)
    return long_df, summary_df, step_df


def run_scalar_diagnostic(cfg, outdir):
    """Scalar covariance diagnostic + mean-overshoot probe (target D)."""
    sd = cfg.get("scalar_diagnostic", {})
    T = float(sd.get("T", 15.0))
    rows = []
    for method in cfg["methods"]:
        for C0 in [float(x) for x in sd.get("C0_values", [0.01, 0.1, 10.0, 100.0])]:
            for dt in [float(x) for x in sd.get("dt_values", [0.5, 1.0, 2.0, 5.0])]:
                for rec in simulate_scalar_diagnostic(method, C0, dt, T, m0=0.0):
                    rows.append({"experiment": "covariance", "method": method,
                                 "C0": C0, "dt": dt, **rec})
        mo = sd.get("mean_overshoot", {"m0": 1.0, "C0": 100.0})
        C0 = float(mo.get("C0", 100.0))
        for dt in [float(x) for x in sd.get("dt_values", [0.5, 1.0, 2.0, 5.0])]:
            for rec in simulate_scalar_diagnostic(method, C0, dt, T, m0=float(mo.get("m0", 1.0))):
                rows.append({"experiment": "mean_overshoot", "method": method,
                             "C0": C0, "dt": dt, **rec})
    df = pd.DataFrame(rows).reindex(columns=SCALAR_COLS)
    save_dataframe(os.path.join(outdir, "scalar_diagnostic.csv"), df)
    return df


def parse_args():
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=DEFAULT_CONFIG)
    pre.add_argument("--smoke", action="store_true")
    known, _ = pre.parse_known_args()
    cfg = load_yaml(known.config)
    if known.smoke:
        cfg = _apply_smoke(cfg)
    parser = argparse.ArgumentParser(parents=[pre], description=__doc__)
    parser.add_argument("--outdir", default=cfg["output_dir"])
    parser.add_argument("--overwrite", action="store_true",
                        help="Clear the output directory before writing.")
    args = parser.parse_args()
    return args, cfg


if __name__ == "__main__":
    args, cfg = parse_args()
    print("=" * 64)
    print("Riemannian vs KL discretization — stepsize grid")
    print("=" * 64)
    print(f"  config : {args.config}{'  [smoke]' if args.smoke else ''}")
    print(f"  targets: {cfg['targets']}")
    print(f"  lambdas: {cfg['lambdas']}   methods: {cfg['methods']}")
    print(f"  dt grid: {cfg['dt_grid']}   T={cfg['T']}")
    print(f"  outdir : {args.outdir}")
    print()
    if args.overwrite and os.path.isdir(args.outdir):
        for entry in os.listdir(args.outdir):
            p = os.path.join(args.outdir, entry)
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    ensure_dir(args.outdir)

    run_grid(cfg, args.outdir)
    run_scalar_diagnostic(cfg, args.outdir)
    print(f"\nWrote results_long.csv, summary.csv, stepsize_summary.csv, "
          f"scalar_diagnostic.csv, target_metadata.json -> {args.outdir}")
