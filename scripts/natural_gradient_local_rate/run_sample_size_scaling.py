"""Sample-size scaling: how the operator estimates depend on M_mc.

For each base grid point (fixed dimension / kappa / potential family / seed) the
script sweeps several Monte-Carlo sample sizes ``M_mc`` and records the full
estimator suite at each. This is the decisive diagnostic for telling apart a
*real* dimension-dependence trend from *finite-sample spectral noise*: for a
separable control the dimension-free exact benchmark is M-independent, so if
Lambda_hat_full_sym shrinks toward the diagonal/exact value as M grows, the
high-dimensional inflation was an estimator artifact.

Usage:
    python scripts/natural_gradient_local_rate/run_sample_size_scaling.py \
        --config configs/natural_gradient_local_rate/sample_size_scaling.yaml \
        [--outdir DIR] [--overwrite]

Output: <base_dir>/sample_size_scaling/{results_long.csv, summary.csv, config.json}
"""
import argparse
import os
import sys
import time

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))
sys.path.insert(0, _HERE)

import _common  # noqa: E402
from src.common.io_utils import save_dataframe, save_json  # noqa: E402
from src.natural_gradient_local_rate.estimator_suite import compute_row  # noqa: E402

DEFAULT_CONFIG = os.path.join(_HERE, "..", "..", "configs",
                              "natural_gradient_local_rate", "sample_size_scaling.yaml")


def parse_args():
    p = argparse.ArgumentParser(description="Sample-size scaling of the operator estimators.")
    p.add_argument("--config", default=DEFAULT_CONFIG)
    p.add_argument("--outdir", default=None, help="Override outputs.base_dir")
    p.add_argument("--overwrite", action="store_true")
    _common.add_backend_cli_args(p)
    return p.parse_args()


def main():
    args = parse_args()
    cfg = _common.load_config(args.config)
    outdir = _common.stage_dir(cfg, args, "sample_size_scaling")
    long_path = os.path.join(outdir, "results_long.csv")
    if os.path.exists(long_path) and not args.overwrite:
        print(f"[exists] {long_path} (use --overwrite to regenerate)")
        return

    opts = _common.operator_opts(cfg)
    opts["compute_gamma_loc"] = True
    _common.apply_cli_overrides(opts, args)
    run_id, group = _common.run_context()
    print(f"[backend={opts['backend']} device={opts['device']} dtype={opts['dtype']} "
          f"eigensolver={opts['eigensolver']}]")

    points = list(_common.scaling_points(cfg))
    rows = []
    t0 = time.time()
    for i, point in enumerate(points, 1):
        Z = _common.make_bank(cfg, point)
        pot = _common.make_potential(cfg, point, Z)
        row = compute_row(pot, Z, point, opts, run_id=run_id, experiment_group=group)
        rows.append(row)
        print(f"  [{i:3d}/{len(points)}] {_common.point_key(point):28s} "
              f"M={point['M_mc']:<7d} sym={row['Lambda_hat_full_sym']:.4f} "
              f"diag={row['Lambda_hat_diag']:.4f} exact={row['Lambda_hat_separable_exact']:.4f} "
              f"gamma={row['gamma_loc']:.4f} [{row['status']}]  ({time.time()-t0:.1f}s)")

    df = _common.order_columns(pd.DataFrame(rows))
    save_dataframe(long_path, df)

    # summary: mean over seeds per (family, N_theta, kappa, M_mc)
    keys = ["potential_family", "N_theta", "kappa_target", "M_mc"]
    agg = df.groupby(keys).agg(
        Lambda_hat_full_sym_mean=("Lambda_hat_full_sym", "mean"),
        Lambda_hat_raw_forward_mean=("Lambda_hat_raw_forward", "mean"),
        Lambda_hat_diag_mean=("Lambda_hat_diag", "mean"),
        Lambda_hat_separable_exact_mean=("Lambda_hat_separable_exact", "mean"),
        gamma_loc_mean=("gamma_loc", "mean"),
        full_sym_minus_exact_mean=("full_sym_minus_exact", "mean"),
        self_adjoint_error_H_sym_max=("self_adjoint_error_H_sym", "max"),
        self_adjoint_error_L_star_max=("self_adjoint_error_L_star", "max"),
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
