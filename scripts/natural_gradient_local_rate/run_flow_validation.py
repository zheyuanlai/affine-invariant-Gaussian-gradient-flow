"""Riemannian flow validation: check the simulated decay rate against gamma_loc.

For each grid point we obtain the slow eigenvector of L_star (reused from the
linearized-rate stage if available, otherwise re-estimated), then run the
exponential-map flow from a small perturbation along it for each
(epsilon, Delta_t) and fit log R^2 vs t. ``fit_gamma_flow`` should match
``gamma_loc`` (ratio ~ 1). This validates the linearized rate; it is not the
primary evidence.

Usage:
    python scripts/natural_gradient_local_rate/run_flow_validation.py \
        --config configs/natural_gradient_local_rate/smoke.yaml [--outdir DIR] [--overwrite]

Output: <base_dir>/flow_validation/{trajectories.csv, summary.csv, config.json}
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))
sys.path.insert(0, _HERE)

import _common  # noqa: E402
from src.common.io_utils import save_dataframe, save_json  # noqa: E402
from src.natural_gradient_local_rate.linearized_rate import estimate_gamma_loc  # noqa: E402
from src.natural_gradient_local_rate.riemannian_flow import run_flow_validation  # noqa: E402

DEFAULT_CONFIG = os.path.join(_HERE, "..", "..", "configs",
                              "natural_gradient_local_rate", "smoke.yaml")


def parse_args():
    p = argparse.ArgumentParser(description="Validate the Riemannian flow against gamma_loc.")
    p.add_argument("--config", default=DEFAULT_CONFIG)
    p.add_argument("--smoke", action="store_true",
                   help="Shortcut for --config <.../smoke.yaml>")
    p.add_argument("--outdir", default=None, help="Override outputs.base_dir")
    p.add_argument("--eigenvector-dir", default=None,
                   help="Where to look for cached eigenvectors "
                        "(default: <base_dir>/linearized_rate_grid/eigenvectors)")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def _load_slow_mode(point, cfg, args, Z, eig, cs):
    """Return (gamma_loc, u_star, X_star), reusing a cached eigenvector if present."""
    base = args.outdir if args.outdir else cfg["outputs"]["base_dir"]
    edir = args.eigenvector_dir or os.path.join(base, "linearized_rate_grid", "eigenvectors")
    path = os.path.join(edir, _common.point_key(point) + ".npz")
    if os.path.exists(path):
        d = np.load(path)
        return float(d["gamma_loc"]), d["u_star"], d["X_star"]
    pot = _common.make_potential(cfg, point, Z)
    gamma, (u_star, X_star) = estimate_gamma_loc(pot, Z, chunk_size=cs,
                                                 return_eigenvector=True, **eig)
    return gamma, u_star, X_star


def main():
    args = parse_args()
    if args.smoke:
        args.config = DEFAULT_CONFIG
    cfg = _common.load_config(args.config)
    outdir = _common.stage_dir(cfg, args, "flow_validation")
    summary_path = os.path.join(outdir, "summary.csv")
    if os.path.exists(summary_path) and not args.overwrite:
        print(f"[exists] {summary_path} (use --overwrite to regenerate)")
        return

    eig = _common.eigsh_opts(cfg)
    cs = _common.chunk_size(cfg)
    fv = cfg["flow_validation"]
    dts = [float(x) for x in fv["Delta_t"]]
    epsilons = [float(x) for x in fv["epsilon"]]
    n_steps = int(fv["n_steps"])
    window = (float(fv.get("fit_start_fraction", 0.05)),
              float(fv.get("fit_end_fraction", 0.5)))

    points = list(_common.grid_points(cfg))
    summary_rows, traj_frames = [], []
    t0 = time.time()
    for i, point in enumerate(points, 1):
        Z = _common.make_bank(cfg, point)
        pot = _common.make_potential(cfg, point, Z)
        gamma, u_star, X_star = _load_slow_mode(point, cfg, args, Z, eig, cs)

        for dt in dts:
            for eps in epsilons:
                out = run_flow_validation(pot, Z, u_star, X_star, epsilon=eps,
                                          Delta_t=dt, n_steps=n_steps,
                                          fit_window=window, chunk_size=cs)
                s = out["summary"]
                fit_gamma = s["fit_gamma_flow"]
                ratio = (fit_gamma / gamma) if gamma else float("nan")
                row = {
                    "family": point["family"], "N_theta": point["N_theta"],
                    "kappa_target": point["kappa_target"], "seed": point["seed"],
                    "gamma_loc": gamma, "Delta_t": dt, "epsilon": eps,
                    "fit_gamma_flow": fit_gamma,
                    "fit_slope_log_R2": s["fit_slope_log_R2"],
                    "fit_r2_flow": s["fit_r2_flow"],
                    "ratio_fit_over_gamma": ratio,
                    "epsilon_used": s["epsilon_used"],
                    "warning": s["warning"],
                }
                summary_rows.append(row)

                tdf = pd.DataFrame(out["trajectory"])
                tdf.insert(0, "key", _common.point_key(point))
                tdf["Delta_t"] = dt
                tdf["epsilon"] = eps
                traj_frames.append(tdf)

        print(f"  [{i:3d}/{len(points)}] {_common.point_key(point):32s} "
              f"gamma_loc={gamma:.4f}  ({time.time()-t0:.1f}s)")

    save_dataframe(summary_path, pd.DataFrame(summary_rows))
    save_dataframe(os.path.join(outdir, "trajectories.csv"),
                   pd.concat(traj_frames, ignore_index=True))
    save_json(os.path.join(outdir, "config.json"),
              {"config_path": os.path.abspath(args.config), "config": cfg})

    print(f"\nSummary      -> {summary_path}")
    print(f"Trajectories -> {os.path.join(outdir, 'trajectories.csv')}")


if __name__ == "__main__":
    main()
