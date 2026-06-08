"""Run the WFR Gaussian gradient-flow grid (CPU, deterministic).

For every ``(target, Lambda, epsilon, method)`` the runner simulates the named
method from the far underdispersed initialization with a shared Fisher--Rao step
``main_dt``, stepping to a fixed expectation-batch budget, classifies the run, and
writes:

    results_long.csv     one row per saved step per run
    summary.csv          one row per run (hitting times in iters and batches)
    hitting_times.csv     compact main hitting-time table (one row per run)
    schedule_sweep.csv   fixed-schedule c-sweep + theory marker (hitting batches)
    target_metadata.json target params, alpha/beta, a_star/E_star, schedules
    run_metadata.json    config echo, timing, environment

Usage::

    python scripts/wfr_gradient_flow/run_wfr_grid.py \
        --config configs/wfr_gradient_flow/wfr_grid.yaml \
        --outdir outputs/wfr_gradient_flow --overwrite

Add ``--smoke`` to load the reduced grid from the smoke config instead.
"""
from __future__ import annotations

import argparse
import os
import platform
import shutil
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.io_utils import load_yaml, ensure_dir, save_dataframe, save_json
from src.natural_gradient_discretization_stepsize.optimize_star import compute_star
from src.wfr_gradient_flow.targets import build_target
from src.wfr_gradient_flow.schedules import build_schedule, theory_mu_min
from src.wfr_gradient_flow.runner import simulate_run
from src.wfr_gradient_flow.metrics import GAP_THRESHOLDS, _tol_key

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "..", "configs",
    "wfr_gradient_flow", "wfr_grid.yaml")

_HIT_COLS = []
for _tol in GAP_THRESHOLDS:
    _HIT_COLS += [f"iter_to_{_tol_key(_tol)}", f"batches_to_{_tol_key(_tol)}"]

LONG_COLS = [
    "target_name", "Lambda", "epsilon", "method", "schedule_name", "dt",
    "h_n", "lambda_n", "s_n", "iteration", "expectation_batches",
    "objective", "objective_gap", "mean_error_norm", "covariance_error_fro",
    "min_eig_C", "max_eig_C", "spd_feasible", "monotone_so_far", "status",
]

SUMMARY_COLS = [
    "target_name", "Lambda", "epsilon", "method", "schedule_name", "dt",
    "alpha", "beta", "h_first", "h_final", "lambda_first",
    "n_iters", "n_batches", "max_batches", "gap0", "gap_final", "gap_min",
    "min_eig_C_min", "spd_feasible", "monotone", "status", "failed_at_iter",
    "wall_time_total",
] + _HIT_COLS

HIT_COLS = [
    "target_name", "Lambda", "epsilon", "method", "schedule_name", "dt",
    "batches_to_1e_minus_1", "batches_to_1e_minus_3", "batches_to_1e_minus_6",
    "iter_to_1e_minus_1", "iter_to_1e_minus_3", "iter_to_1e_minus_6",
    "gap_final",
]

SWEEP_COLS = [
    "target_name", "Lambda", "epsilon", "schedule_kind", "c", "h", "dt",
    "batches_to_1e_minus_3", "iter_to_1e_minus_3", "gap_final",
    "spd_feasible", "monotone", "status",
]

DTSWEEP_COLS = [
    "target_name", "Lambda", "epsilon", "dt", "c", "h",
    "batches_to_1e_minus_3", "gap_final", "spd_feasible", "monotone", "status",
]


def _build_target(name, Lambda, epsilon, cfg):
    if name == "smooth_log_cosh":
        st = cfg.get("smooth_target", {})
        return build_target(name, Lambda, epsilon, r=float(cfg.get("r", 8.0)),
                            delta=float(st.get("delta", 0.05)),
                            gamma=float(st.get("gamma", 1.0)),
                            n_nodes=int(st.get("gh_nodes", 80)))
    return build_target(name, Lambda, epsilon, r=float(cfg.get("r", 8.0)))


def _star(target):
    """Reference optimum ``(m_star, C_star, F_star, diag)`` (closed form or L-BFGS-B)."""
    return compute_star(target)


# ---------------------------------------------------------------------------
# Main method comparison (one run per target/Lambda/epsilon/method)
# ---------------------------------------------------------------------------

def run_main_grid(cfg, outdir):
    """Five-method comparison at the shared Fisher--Rao step ``main_dt``."""
    dt = float(cfg["main_dt"])
    max_batches = int(cfg["max_batches"])
    c = float(cfg["fixed_c"])
    adp = cfg.get("adaptive", {})
    h_max_frac = float(adp.get("h_max_frac", 0.9))
    s0 = float(adp.get("s0", 0.5))
    max_saved = int(cfg.get("max_saved_rows", 600))

    long_rows, summary_rows, hit_rows = [], [], []
    metadata = {"targets": {}}
    combos = [(t, L, e) for t in cfg["targets"]
              for L in cfg["Lambdas"] for e in cfg["epsilons"]]
    total = len(combos) * len(cfg["methods"])
    done = 0
    t0 = time.time()

    for tname, Lambda, eps in combos:
        target = _build_target(tname, Lambda, eps, cfg)
        m_star, C_star, F_star, opt_diag = _star(target)
        mu_min = theory_mu_min(target.C0, target.beta)
        meta_key = f"{tname}__L{Lambda:g}__eps{eps:g}"
        md = target.metadata()
        md.update({
            "m_star": np.asarray(m_star).tolist(),
            "C_star": np.asarray(C_star).tolist(),
            "F_star": float(F_star), "optimizer": opt_diag,
            "mu_min_theory": float(mu_min),
            "gap0": float(_gap0(target, F_star)),
            "schedules": {
                "fr_only": {"h": 0.0},
                "w_only": {"h": c / target.beta, "c": c},
                "wfr_fixed": {"h": c / target.beta, "c": c},
                "wfr_theory": {"h": mu_min,
                               "note": "h=mu_min=min(lambda_min(C0),1/beta)"},
                "wfr_adaptive": {"h_max": h_max_frac / target.beta,
                                 "h_max_frac": h_max_frac, "s0": s0},
            },
        })
        metadata["targets"][meta_key] = md

        for method in cfg["methods"]:
            schedule = build_schedule(method, beta=target.beta, C0=target.C0,
                                      c=c, h_max_frac=h_max_frac, s0=s0)
            records, summ = simulate_run(
                method, target, dt, max_batches, schedule,
                F_star, m_star, C_star, max_saved_rows=max_saved)
            long_rows.extend(records)
            summary_rows.append(summ)
            hit_rows.append({k: summ.get(k) for k in HIT_COLS})
            done += 1
            print(f"  [{done:3d}/{total}] {tname:16s} L={Lambda:<7g} "
                  f"eps={eps:<6g} {method:13s} "
                  f"b3={summ.get('batches_to_1e_minus_3'):>6} "
                  f"gapf={summ['gap_final']:.2e} "
                  f"[{summ['status']}] ({time.time() - t0:5.1f}s)")

    save_dataframe(os.path.join(outdir, "results_long.csv"),
                   pd.DataFrame(long_rows).reindex(columns=LONG_COLS))
    save_dataframe(os.path.join(outdir, "summary.csv"),
                   pd.DataFrame(summary_rows).reindex(columns=SUMMARY_COLS))
    save_dataframe(os.path.join(outdir, "hitting_times.csv"),
                   pd.DataFrame(hit_rows).reindex(columns=HIT_COLS))
    save_json(os.path.join(outdir, "target_metadata.json"), metadata)
    return metadata


def _gap0(target, F_star):
    return (target.energy_gap(target.m0, target.C0)
            if hasattr(target, "energy_gap")
            else target.objective(target.m0, target.C0) - F_star)


# ---------------------------------------------------------------------------
# Schedule sweep: fixed WFR schedule h = c/beta vs the theory marker h = mu_min
# ---------------------------------------------------------------------------

def _fixed_run(target, h, dt, max_batches, F_star, m_star, C_star,
               schedule_name="wfr_fixed", max_saved=300):
    """One full-WFR run at a constant Wasserstein step ``h`` (c=0 => fr_only)."""
    from src.wfr_gradient_flow.schedules import _ConstantSchedule
    method = "fr_only" if h == 0.0 else "wfr_fixed"
    sched = _ConstantSchedule(schedule_name, h)
    return simulate_run(method, target, dt, max_batches, sched,
                        F_star, m_star, C_star, max_saved_rows=max_saved)


def run_schedule_sweep(cfg, outdir):
    """Fixed-c sweep of the WFR schedule + the theorem-bound h=mu_min marker."""
    dt = float(cfg["main_dt"])
    max_batches = int(cfg["max_batches"])
    c_grid = [float(x) for x in cfg.get("schedule_sweep", {}).get("c_grid", [])]
    rows = []
    combos = [(t, L, e) for t in cfg["targets"]
              for L in cfg["Lambdas"] for e in cfg["epsilons"]]
    for tname, Lambda, eps in combos:
        target = _build_target(tname, Lambda, eps, cfg)
        m_star, C_star, F_star, _ = _star(target)
        beta = target.beta
        for c in c_grid:
            h = c / beta
            _, s = _fixed_run(target, h, dt, max_batches, F_star, m_star, C_star)
            rows.append(_sweep_row(tname, Lambda, eps, "fixed_c", c, h, dt, s))
        # Theory marker: h = mu_min (reported as a separate schedule_kind).
        mu_min = theory_mu_min(target.C0, beta)
        _, s = _fixed_run(target, mu_min, dt, max_batches, F_star, m_star, C_star,
                          schedule_name="wfr_theory")
        rows.append(_sweep_row(tname, Lambda, eps, "theory_mu_min",
                               mu_min * beta, mu_min, dt, s))
        print(f"  sweep {tname:16s} L={Lambda:<7g} eps={eps:<6g} done")
    save_dataframe(os.path.join(outdir, "schedule_sweep.csv"),
                   pd.DataFrame(rows).reindex(columns=SWEEP_COLS))


def _sweep_row(tname, Lambda, eps, kind, c, h, dt, s):
    return {
        "target_name": tname, "Lambda": Lambda, "epsilon": eps,
        "schedule_kind": kind, "c": float(c), "h": float(h), "dt": float(dt),
        "batches_to_1e_minus_3": s.get("batches_to_1e_minus_3"),
        "iter_to_1e_minus_3": s.get("iter_to_1e_minus_3"),
        "gap_final": s["gap_final"], "spd_feasible": s["spd_feasible"],
        "monotone": s["monotone"], "status": s["status"],
    }


# ---------------------------------------------------------------------------
# Optional Delta t / c stability heatmap sweep for wfr_fixed (Figure 5)
# ---------------------------------------------------------------------------

def run_dt_sweep(cfg, outdir):
    """Delta t x c grid for wfr_fixed: hitting batches + stability flags."""
    ds = cfg.get("dt_sweep", {})
    dt_grid = [float(x) for x in ds.get("dt_grid", [])]
    c_grid = [float(x) for x in ds.get("c_grid", [])]
    if not dt_grid or not c_grid:
        return
    max_batches = int(cfg["max_batches"])
    rows = []
    combos = [(t, L, e) for t in cfg["targets"]
              for L in cfg["Lambdas"] for e in cfg["epsilons"]]
    for tname, Lambda, eps in combos:
        target = _build_target(tname, Lambda, eps, cfg)
        m_star, C_star, F_star, _ = _star(target)
        beta = target.beta
        for dt in dt_grid:
            for c in c_grid:
                h = c / beta
                _, s = _fixed_run(target, h, dt, max_batches,
                                  F_star, m_star, C_star)
                rows.append({
                    "target_name": tname, "Lambda": Lambda, "epsilon": eps,
                    "dt": dt, "c": c, "h": h,
                    "batches_to_1e_minus_3": s.get("batches_to_1e_minus_3"),
                    "gap_final": s["gap_final"],
                    "spd_feasible": s["spd_feasible"],
                    "monotone": s["monotone"], "status": s["status"],
                })
        print(f"  dt-sweep {tname:16s} L={Lambda:<7g} eps={eps:<6g} done")
    save_dataframe(os.path.join(outdir, "dt_sweep.csv"),
                   pd.DataFrame(rows).reindex(columns=DTSWEEP_COLS))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=DEFAULT_CONFIG)
    pre.add_argument("--smoke", action="store_true",
                     help="Load the reduced smoke config instead of --config.")
    known, _ = pre.parse_known_args()
    cfg_path = known.config
    if known.smoke:
        cfg_path = os.path.join(os.path.dirname(__file__), "..", "..", "configs",
                                "wfr_gradient_flow", "wfr_smoke.yaml")
    cfg = load_yaml(cfg_path)
    parser = argparse.ArgumentParser(parents=[pre], description=__doc__)
    parser.add_argument("--outdir", default=cfg["output_dir"])
    parser.add_argument("--overwrite", action="store_true",
                        help="Clear the output directory before writing.")
    args = parser.parse_args()
    return args, cfg, cfg_path


if __name__ == "__main__":
    args, cfg, cfg_path = parse_args()
    print("=" * 64)
    print("Wasserstein--Fisher--Rao Gaussian gradient flow -- grid")
    print("=" * 64)
    print(f"  config : {cfg_path}{'  [smoke]' if args.smoke else ''}")
    print(f"  targets: {cfg['targets']}")
    print(f"  Lambdas: {cfg['Lambdas']}   epsilons: {cfg['epsilons']}")
    print(f"  methods: {cfg['methods']}")
    print(f"  main_dt: {cfg['main_dt']}   max_batches: {cfg['max_batches']}")
    print(f"  outdir : {args.outdir}\n")

    if args.overwrite and os.path.isdir(args.outdir):
        for entry in os.listdir(args.outdir):
            p = os.path.join(args.outdir, entry)
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    ensure_dir(args.outdir)

    t_start = time.time()
    print("main method comparison:")
    run_main_grid(cfg, args.outdir)
    print("schedule sweep:")
    run_schedule_sweep(cfg, args.outdir)
    print("dt sweep:")
    run_dt_sweep(cfg, args.outdir)

    save_json(os.path.join(args.outdir, "run_metadata.json"), {
        "config_path": cfg_path, "config": cfg,
        "wall_time_total": float(time.time() - t_start),
        "python": platform.python_version(), "platform": platform.platform(),
        "numpy": np.__version__, "smoke": bool(args.smoke),
    })
    print(f"\nWrote results_long.csv, summary.csv, hitting_times.csv, "
          f"schedule_sweep.csv, dt_sweep.csv, target_metadata.json, "
          f"run_metadata.json -> {args.outdir}")
