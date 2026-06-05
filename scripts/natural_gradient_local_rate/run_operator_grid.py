"""Operator-norm experiment: estimate the operator norm Lambda_hat.

For each grid point we build a fixed sample bank, construct the (centered)
potential on that same bank, and estimate every operator-norm mode side by side:

* Lambda_hat_full_sym   -- lambda_max of the self-adjoint H_sym (the default);
* Lambda_hat_raw_forward -- largest real eigenvalue of the uncorrected forward
  H_lin (diagnostic only);
* Lambda_hat_diag        -- diagonal-restricted operator norm (lambda_max of
  A = G - 1 1^T), a separable sanity check;
* Lambda_hat_separable_exact -- Gauss--Hermite ground truth (separable only).

Self-adjointness errors and eigen-residuals are recorded alongside. The
symmetrized estimator is the headline; raw_forward is kept for comparison.

Usage:
    python scripts/natural_gradient_local_rate/run_operator_grid.py \
        --config configs/natural_gradient_local_rate/smoke.yaml [--outdir DIR] [--overwrite]

Output: <base_dir>/operator_grid/{results_long.csv, summary.csv, config.json}
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
from src.natural_gradient_local_rate.estimator_suite import compute_row  # noqa: E402

DEFAULT_CONFIG = os.path.join(_HERE, "..", "..", "configs",
                              "natural_gradient_local_rate", "smoke.yaml")


def parse_args():
    p = argparse.ArgumentParser(description="Estimate Lambda_hat over a grid.")
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
    outdir = _common.stage_dir(cfg, args, "operator_grid")
    long_path = os.path.join(outdir, "results_long.csv")
    if os.path.exists(long_path) and not args.overwrite:
        print(f"[exists] {long_path} (use --overwrite to regenerate)")
        return
    opts = _common.operator_opts(cfg)
    _common.apply_cli_overrides(opts, args)
    run_id, group = _common.run_context()

    points = list(_common.grid_points(cfg))
    rows = []
    t0 = time.time()
    for i, point in enumerate(points, 1):
        point["M_mc"] = _common.grid_M_mc(cfg, point)
        Z = _common.make_bank(cfg, point)
        pot = _common.make_potential_for_opts(cfg, point, Z, opts)
        row = compute_row(pot, Z, point, opts, run_id=run_id, experiment_group=group)
        rows.append(row)
        print(f"  [{i:3d}/{len(points)}] {_common.point_key(point):32s} "
              f"sym={row['Lambda_hat_full_sym']:.4f} diag={row['Lambda_hat_diag']:.4f} "
              f"exact={row['Lambda_hat_separable_exact']:.4f} [{row['status']}] "
              f"({time.time()-t0:.1f}s)")

    df = _common.order_columns(pd.DataFrame(rows))
    save_dataframe(long_path, df)

    # summary: mean/std over seeds per (family, N_theta, kappa_target)
    keys = ["potential_family", "N_theta", "kappa_target"]
    agg = df.groupby(keys).agg(
        Lambda_hat_full_sym_mean=("Lambda_hat_full_sym", "mean"),
        Lambda_hat_full_sym_std=("Lambda_hat_full_sym", "std"),
        Lambda_hat_diag_mean=("Lambda_hat_diag", "mean"),
        Lambda_hat_separable_exact_mean=("Lambda_hat_separable_exact", "mean"),
        full_minus_diag_mean=("full_minus_diag", "mean"),
        lambda_over_logkappa_mean=("lambda_over_logkappa", "mean"),
        self_adjoint_error_H_sym_max=("self_adjoint_error_H_sym", "max"),
        n_seeds=("seed", "count"),
    ).reset_index()
    save_dataframe(os.path.join(outdir, "summary.csv"), agg)
    save_json(os.path.join(outdir, "config.json"),
              {"config_path": os.path.abspath(args.config), "config": cfg,
               "run_id": run_id})

    print(f"\nLong CSV -> {long_path}")
    print(f"Summary  -> {os.path.join(outdir, 'summary.csv')}")


if __name__ == "__main__":
    main()
