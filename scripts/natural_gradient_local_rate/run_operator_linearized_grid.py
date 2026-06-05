"""Joint operator-norm + linearized-rate grid runner.

This runner is the production entry point when both ``operator_grid`` and
``linearized_rate_grid`` outputs are needed for the same config. It builds each
sample bank and potential once, calls the shared estimator suite once, then
writes both stage outputs:

* ``<base>/operator_grid/{results_long.csv, summary.csv, config.json}``
* ``<base>/linearized_rate_grid/{results_long.csv, summary.csv, eigenvectors, config.json}``

On the torch backend this avoids the expensive duplicate dense accumulation that
would happen if ``run_operator_grid.py`` and ``run_linearized_rate_grid.py`` were
run separately.

Usage:
    python scripts/natural_gradient_local_rate/run_operator_linearized_grid.py \
        --config configs/natural_gradient_local_rate/production_all.yaml \
        --backend torch --device cuda --overwrite
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
from src.common.io_utils import ensure_dir, save_dataframe, save_json  # noqa: E402
from src.natural_gradient_local_rate.estimator_suite import compute_row  # noqa: E402

DEFAULT_CONFIG = os.path.join(_HERE, "..", "..", "configs",
                              "natural_gradient_local_rate", "smoke.yaml")


def parse_args():
    p = argparse.ArgumentParser(
        description="Estimate Lambda_hat and gamma_loc over one shared grid pass.")
    p.add_argument("--config", default=DEFAULT_CONFIG)
    p.add_argument("--smoke", action="store_true",
                   help="Shortcut for --config <.../smoke.yaml>")
    p.add_argument("--outdir", default=None, help="Override outputs.base_dir")
    p.add_argument("--overwrite", action="store_true")
    _common.add_backend_cli_args(p)
    return p.parse_args()


def _stage_dir(base, stage):
    return ensure_dir(os.path.join(base, stage))


def _operator_summary(df):
    keys = ["potential_family", "N_theta", "kappa_target"]
    return df.groupby(keys).agg(
        Lambda_hat_full_sym_mean=("Lambda_hat_full_sym", "mean"),
        Lambda_hat_full_sym_std=("Lambda_hat_full_sym", "std"),
        Lambda_hat_diag_mean=("Lambda_hat_diag", "mean"),
        Lambda_hat_separable_exact_mean=("Lambda_hat_separable_exact", "mean"),
        full_minus_diag_mean=("full_minus_diag", "mean"),
        lambda_over_logkappa_mean=("lambda_over_logkappa", "mean"),
        self_adjoint_error_H_sym_max=("self_adjoint_error_H_sym", "max"),
        n_seeds=("seed", "count"),
    ).reset_index()


def _linearized_summary(df):
    keys = ["potential_family", "N_theta", "kappa_target"]
    return df.groupby(keys).agg(
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


def main():
    args = parse_args()
    if args.smoke:
        args.config = DEFAULT_CONFIG
    cfg = _common.load_config(args.config)
    base = args.outdir if args.outdir else cfg["outputs"]["base_dir"]
    opdir = _stage_dir(base, "operator_grid")
    lrdir = _stage_dir(base, "linearized_rate_grid")
    op_long = os.path.join(opdir, "results_long.csv")
    lr_long = os.path.join(lrdir, "results_long.csv")
    existing = [p for p in (op_long, lr_long) if os.path.exists(p)]
    if existing and not args.overwrite:
        print("[exists] " + ", ".join(existing) + " (use --overwrite to regenerate)")
        return

    opts = _common.operator_opts(cfg)
    opts["compute_gamma_loc"] = True
    _common.apply_cli_overrides(opts, args)
    run_id, group = _common.run_context()
    print(f"[joint backend={opts['backend']} device={opts['device']} dtype={opts['dtype']} "
          f"eigensolver={opts['eigensolver']}]")

    lr_cfg = cfg.get("linearized_rate", {})
    save_eigs = bool(lr_cfg.get("save_eigenvectors", True))
    eigdir = ensure_dir(os.path.join(lrdir, "eigenvectors")) if save_eigs else None

    points = list(_common.grid_points(cfg))
    rows = []
    t0 = time.time()
    for i, point in enumerate(points, 1):
        point["M_mc"] = _common.grid_M_mc(cfg, point)
        Z = _common.make_bank(cfg, point)
        pot = _common.make_potential_for_opts(cfg, point, Z, opts)
        row = compute_row(pot, Z, point, opts, run_id=run_id, experiment_group=group)

        if save_eigs and "_u_star" in row:
            np.savez(os.path.join(eigdir, _common.point_key(point) + ".npz"),
                     u_star=row["_u_star"], X_star=row["_X_star"],
                     gamma_loc=row["gamma_loc"], N_theta=point["N_theta"])
        rows.append(row)

        print(f"  [{i:3d}/{len(points)}] {_common.point_key(point):32s} "
              f"sym={row['Lambda_hat_full_sym']:.4f} "
              f"gamma={row['gamma_loc']:.4f} "
              f"diag={row['Lambda_hat_diag']:.4f} [{row['status']}] "
              f"({time.time()-t0:.1f}s)")

    df = _common.order_columns(pd.DataFrame(rows))
    save_dataframe(op_long, df)
    save_dataframe(lr_long, df)
    save_dataframe(os.path.join(opdir, "summary.csv"), _operator_summary(df))
    save_dataframe(os.path.join(lrdir, "summary.csv"), _linearized_summary(df))

    meta = {"config_path": os.path.abspath(args.config), "config": cfg,
            "run_id": run_id, "runner": "run_operator_linearized_grid.py"}
    save_json(os.path.join(opdir, "config.json"), meta)
    save_json(os.path.join(lrdir, "config.json"), meta)

    print(f"\nOperator long -> {op_long}")
    print(f"Operator summary -> {os.path.join(opdir, 'summary.csv')}")
    print(f"Linearized long -> {lr_long}")
    print(f"Linearized summary -> {os.path.join(lrdir, 'summary.csv')}")
    if save_eigs:
        print(f"Eigvecs -> {eigdir}")


if __name__ == "__main__":
    main()
