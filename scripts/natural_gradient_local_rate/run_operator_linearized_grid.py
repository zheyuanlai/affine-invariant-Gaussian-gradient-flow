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
        --config configs/natural_gradient_local_rate/gpu_lowdim_operator_full.yaml \
        --backend torch --device cuda --dtype float64 \
        --chunk-size 1048576 --overwrite
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
        tau_H_sq_mean=("tau_H_sq", "mean"),
        coupling_bound_rate_mean=("coupling_bound_rate", "mean"),
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
        tau_H_sq_mean=("tau_H_sq", "mean"),
        coupling_bound_rate_mean=("coupling_bound_rate", "mean"),
        gamma_over_coupling_bound_min=("gamma_over_coupling_bound", "min"),
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
    suffix = _common.shard_suffix(args)
    opdir = _stage_dir(base, "operator_grid")
    lrdir = _stage_dir(base, "linearized_rate_grid")
    op_long = os.path.join(opdir, f"results_long{suffix}.csv")
    lr_long = os.path.join(lrdir, f"results_long{suffix}.csv")
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

    all_points = list(_common.grid_points(cfg))
    points = _common.apply_grid_shard(all_points, args)
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
              f"tau2={row['tau_H_sq']:.4f} "
              f"diag={row['Lambda_hat_diag']:.4f} [{row['status']}] "
              f"({time.time()-t0:.1f}s)")

    df = _common.order_columns(pd.DataFrame(rows))
    save_dataframe(op_long, df)
    save_dataframe(lr_long, df)
    op_summary = os.path.join(opdir, f"summary{suffix}.csv")
    lr_summary = os.path.join(lrdir, f"summary{suffix}.csv")
    save_dataframe(op_summary, _operator_summary(df))
    save_dataframe(lr_summary, _linearized_summary(df))

    meta = {"config_path": os.path.abspath(args.config), "config": cfg,
            "run_id": run_id, "runner": "run_operator_linearized_grid.py",
            "num_shards": int(args.num_shards), "shard_index": int(args.shard_index),
            "n_points_total": len(all_points), "n_points_this_shard": len(points)}
    save_json(os.path.join(opdir, f"config{suffix}.json"), meta)
    save_json(os.path.join(lrdir, f"config{suffix}.json"), meta)

    print(f"\nOperator long -> {op_long}")
    print(f"Operator summary -> {op_summary}")
    print(f"Linearized long -> {lr_long}")
    print(f"Linearized summary -> {lr_summary}")
    if save_eigs:
        print(f"Eigvecs -> {eigdir}")


if __name__ == "__main__":
    main()
