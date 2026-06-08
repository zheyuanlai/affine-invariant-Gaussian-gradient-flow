"""Run the supplementary theoretical-rate benchmark (CPU, deterministic).

For every globally-smooth ``(target, lambda, method, c)`` the runner steps the
scheme for ``N = ceil(T/dt)`` steps at ``dt = c * dt_ref`` on the common
Riemannian-scale reference ``dt_ref = 1/(beta * lambda_max)``, then compares the
observed contraction against the method-specific theorem contraction factor. It
writes four files:

    rate_results_long.csv      one row per saved step per run
    rate_summary.csv           one row per (target, lambda, method, c)
    rate_tolerance_summary.csv one row per (target, lambda, method, c, eps)
    rate_metadata.json         config, per-target alpha/beta/dt_ref, theory grid

Usage::

    python scripts/natural_gradient_discretization_stepsize/run_rate_benchmark.py \
        --config configs/natural_gradient_discretization_stepsize/rate_benchmark.yaml \
        --outdir outputs/natural_gradient_discretization_stepsize --overwrite

Add ``--smoke`` for the fast reduced grid defined in the config. The benchmark
never overwrites the stepsize-study CSVs: all files are prefixed ``rate_``.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.io_utils import load_yaml, ensure_dir, save_dataframe, save_json
from src.natural_gradient_discretization_stepsize.targets import build_target
from src.natural_gradient_discretization_stepsize.optimize_star import compute_star
from src.natural_gradient_discretization_stepsize import rate_metrics as rm
from src.natural_gradient_discretization_stepsize.rate_runner import simulate_rate_run

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "..", "configs",
    "natural_gradient_discretization_stepsize", "rate_benchmark.yaml")

LONG_COLS = [
    "target", "lambda", "method", "alpha", "beta", "lambda_min", "lambda_max",
    "dt_ref", "c", "dt", "n", "time", "gap", "gap_raw", "is_floor_limited",
    "q_riem_theory", "q_kl_formula", "r_riem_theory", "r_kl_formula",
]

SUMMARY_COLS = [
    "target", "lambda", "method", "alpha", "beta", "lambda_min", "lambda_max",
    "dt_ref", "c", "dt", "N", "T_actual", "initial_gap", "final_gap_raw",
    "final_gap_for_logs", "q_theory", "r_theory", "q_hat_terminal",
    "r_hat_terminal", "q_hat_fit", "r_hat_fit", "fit_num_points",
    "fit_window_status", "theory_terminal_bound", "terminal_slack",
    "log10_terminal_slack", "floor_limited_final", "r_continuous", "status",
]

TOL_COLS = [
    "target", "lambda", "method", "c", "dt", "eps", "N_obs", "T_obs",
    "obs_status", "N_theory", "T_theory", "theory_status", "N_theory_over_N_obs",
]


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

def run_rate_grid(cfg, outdir):
    ensure_dir(outdir)
    T = float(cfg["T"])
    c_grid = [float(x) for x in cfg["c_grid"]]
    eps_list = [float(x) for x in cfg["eps_values"]]
    gap_floor = float(cfg.get("gap_floor", rm.GAP_FLOOR))

    long_rows, summary_rows, tol_rows = [], [], []
    metadata = {
        "config": {k: cfg[k] for k in
                   ("targets", "lambdas", "methods", "c_grid", "T",
                    "eps_values", "gap_floor")},
        "note": ("Common Riemannian-scale reference dt_ref=1/(beta*lambda_max) "
                 "for BOTH schemes; the KL contraction formula is a formal "
                 "benchmark evaluated on this grid, not a proof that the KL "
                 "theorem holds at these stepsizes."),
        "targets": {},
    }

    total = len(cfg["targets"]) * len(cfg["lambdas"]) * len(cfg["methods"]) * len(c_grid)
    done = 0
    t0 = time.time()

    for target_name in cfg["targets"]:
        for lam in cfg["lambdas"]:
            target = _build_target(target_name, lam, cfg)
            alpha = float(target.alpha)
            beta = float(target.beta)
            lam_min, lam_max = rm.spectral_bounds(alpha, beta, target.C0)
            dtr = rm.dt_ref(beta, lam_max)
            _, _, F_star, opt_diag = compute_star(target)

            meta_key = f"{target_name}__lam{lam:g}"
            md = target.metadata()
            md.update({
                "alpha": alpha, "beta": beta,
                "lambda_min": lam_min, "lambda_max": lam_max,
                "dt_ref": dtr, "F_star": float(F_star),
                "r_continuous": rm.r_continuous(alpha, lam_min),
                "optimizer": opt_diag,
                "c_grid": c_grid,
                "dt_grid": [c * dtr for c in c_grid],
            })
            metadata["targets"][meta_key] = md

            for method in cfg["methods"]:
                for c in c_grid:
                    records, summ, trows = simulate_rate_run(
                        method, target, c, dtr, alpha, beta, lam_min, lam_max,
                        T, F_star, eps_list, gap_floor=gap_floor,
                        max_saved_rows=int(cfg.get("max_saved_rows", 400)))
                    for rec in records:
                        row = {"target": target_name, "lambda": lam,
                               "method": method, "alpha": alpha, "beta": beta,
                               "lambda_min": lam_min, "lambda_max": lam_max,
                               "dt_ref": dtr, "c": c}
                        row.update(rec)
                        long_rows.append(row)
                    summary_rows.append(summ)
                    tol_rows.extend(trows)
                    done += 1
                    print(f"  [{done:3d}/{total}] {target_name:20s} "
                          f"lam={lam:<5g} {method:10s} c={c:<5g} "
                          f"dt={c * dtr:.3e} slack10={summ['log10_terminal_slack']:+.1f} "
                          f"({time.time() - t0:5.1f}s)")

    long_df = pd.DataFrame(long_rows).reindex(columns=LONG_COLS)
    summary_df = pd.DataFrame(summary_rows).reindex(columns=SUMMARY_COLS)
    tol_df = pd.DataFrame(tol_rows).reindex(columns=TOL_COLS)
    save_dataframe(os.path.join(outdir, "rate_results_long.csv"), long_df)
    save_dataframe(os.path.join(outdir, "rate_summary.csv"), summary_df)
    save_dataframe(os.path.join(outdir, "rate_tolerance_summary.csv"), tol_df)
    save_json(os.path.join(outdir, "rate_metadata.json"), metadata)
    return long_df, summary_df, tol_df


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
                        help="Remove the rate_* files in the output dir before writing.")
    args = parser.parse_args()
    return args, cfg


_RATE_FILES = ["rate_results_long.csv", "rate_summary.csv",
               "rate_tolerance_summary.csv", "rate_metadata.json"]


if __name__ == "__main__":
    args, cfg = parse_args()
    print("=" * 64)
    print("Theoretical-rate benchmark — Riemannian-scale stepsize grid")
    print("=" * 64)
    print(f"  config : {args.config}{'  [smoke]' if args.smoke else ''}")
    print(f"  targets: {cfg['targets']}")
    print(f"  lambdas: {cfg['lambdas']}   methods: {cfg['methods']}")
    print(f"  c grid : {cfg['c_grid']}   T={cfg['T']}")
    print(f"  outdir : {args.outdir}")
    print()
    ensure_dir(args.outdir)
    if args.overwrite:
        # Only remove this benchmark's own files; never touch the stepsize CSVs.
        for fn in _RATE_FILES:
            p = os.path.join(args.outdir, fn)
            if os.path.exists(p):
                os.remove(p)

    run_rate_grid(cfg, args.outdir)
    print(f"\nWrote {', '.join(_RATE_FILES)} -> {args.outdir}")

