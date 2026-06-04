"""Local-rate experiment: estimate gamma_loc = lambda_min(L_star).

gamma_loc is the primary numerical estimate of the local convergence rate; it is
always computed from the self-adjoint H_sym in Fisher--Rao-whitened coordinates.
The full operator-norm suite (raw / symmetrized / diagonal / separable-exact) is
estimated alongside, and the slow eigenvector (u_star, X_star) is saved so the
flow-validation stage can reuse it.

Usage:
    python scripts/natural_gradient_local_rate/run_linearized_rate_grid.py \
        --config configs/natural_gradient_local_rate/smoke.yaml [--outdir DIR] [--overwrite]

Output: <base_dir>/linearized_rate_grid/{results_long.csv, summary.csv,
        eigenvectors/<key>.npz, config.json}
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
from src.common.io_utils import save_dataframe, save_json, ensure_dir  # noqa: E402
from src.natural_gradient_local_rate.estimator_suite import compute_row  # noqa: E402

DEFAULT_CONFIG = os.path.join(_HERE, "..", "..", "configs",
                              "natural_gradient_local_rate", "smoke.yaml")


def parse_args():
    p = argparse.ArgumentParser(description="Estimate gamma_loc (and Lambda_hat) over a grid.")
    p.add_argument("--config", default=DEFAULT_CONFIG)
    p.add_argument("--smoke", action="store_true",
                   help="Shortcut for --config <.../smoke.yaml>")
    p.add_argument("--outdir", default=None, help="Override outputs.base_dir")
    p.add_argument("--overwrite", action="store_true")
    _common.add_backend_cli_args(p)
    return p.parse_args()


def main():
    args = parse_args()
    if args.smoke:
        args.config = DEFAULT_CONFIG
    cfg = _common.load_config(args.config)
    outdir = _common.stage_dir(cfg, args, "linearized_rate_grid")
    long_path = os.path.join(outdir, "results_long.csv")
    if os.path.exists(long_path) and not args.overwrite:
        print(f"[exists] {long_path} (use --overwrite to regenerate)")
        return

    opts = _common.operator_opts(cfg)
    opts["compute_gamma_loc"] = True   # gamma_loc is the point of this stage
    _common.apply_cli_overrides(opts, args)
    run_id, group = _common.run_context()
    lr_cfg = cfg.get("linearized_rate", {})
    save_eigs = bool(lr_cfg.get("save_eigenvectors", True))
    eigdir = ensure_dir(os.path.join(outdir, "eigenvectors")) if save_eigs else None

    points = list(_common.grid_points(cfg))
    rows = []
    t0 = time.time()
    for i, point in enumerate(points, 1):
        point["M_mc"] = _common.grid_M_mc(cfg, point)
        Z = _common.make_bank(cfg, point)
        pot = _common.make_potential(cfg, point, Z)
        row = compute_row(pot, Z, point, opts, run_id=run_id, experiment_group=group)

        if save_eigs and "_u_star" in row:
            np.savez(os.path.join(eigdir, _common.point_key(point) + ".npz"),
                     u_star=row["_u_star"], X_star=row["_X_star"],
                     gamma_loc=row["gamma_loc"], N_theta=point["N_theta"])
        rows.append(row)

        print(f"  [{i:3d}/{len(points)}] {_common.point_key(point):32s} "
              f"gamma_loc={row['gamma_loc']:.4f} sym={row['Lambda_hat_full_sym']:.4f} "
              f"diag={row['Lambda_hat_diag']:.4f} [{row['status']}]  "
              f"({time.time()-t0:.1f}s)")

    df = _common.order_columns(pd.DataFrame(rows))
    save_dataframe(long_path, df)

    keys = ["potential_family", "N_theta", "kappa_target"]
    agg = df.groupby(keys).agg(
        gamma_loc_mean=("gamma_loc", "mean"),
        gamma_loc_std=("gamma_loc", "std"),
        gamma_loc_min=("gamma_loc", "min"),
        inverse_gamma_loc_mean=("inverse_gamma_loc", "mean"),
        inverse_gamma_over_logkappa_mean=("inverse_gamma_over_logkappa", "mean"),
        Lambda_hat_full_sym_mean=("Lambda_hat_full_sym", "mean"),
        Lambda_hat_diag_mean=("Lambda_hat_diag", "mean"),
        Lambda_hat_separable_exact_mean=("Lambda_hat_separable_exact", "mean"),
        current_bound_rate=("current_bound_rate", "first"),
        conjecture_bound_rate=("conjecture_bound_rate", "first"),
        self_adjoint_error_L_star_max=("self_adjoint_error_L_star", "max"),
        n_seeds=("seed", "count"),
    ).reset_index()
    save_dataframe(os.path.join(outdir, "summary.csv"), agg)
    save_json(os.path.join(outdir, "config.json"),
              {"config_path": os.path.abspath(args.config), "config": cfg,
               "run_id": run_id})

    print(f"\nLong CSV -> {long_path}")
    print(f"Summary  -> {os.path.join(outdir, 'summary.csv')}")
    if save_eigs:
        print(f"Eigvecs  -> {eigdir}")


if __name__ == "__main__":
    main()
