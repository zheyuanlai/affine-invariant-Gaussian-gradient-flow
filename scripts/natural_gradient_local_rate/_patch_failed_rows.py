"""Patch the failed rows of an existing joint-runner output in place.

Re-runs only the grid points whose ``status != "ok"`` in an existing
operator_grid/linearized_rate_grid results CSV (e.g. CUDA-OOM rows from a
contended GPU), on the current device, then merges the fresh rows back, rewrites
the eigenvector .npz files for the patched points, and regenerates both stage
summaries. The good rows are left untouched (same run_id preserved).

Usage:
    CUDA_VISIBLE_DEVICES=2 python scripts/natural_gradient_local_rate/_patch_failed_rows.py \
        --config configs/natural_gradient_local_rate/gpu_lowdim_operator_full.yaml \
        --outdir outputs/natural_gradient_local_rate \
        --device cuda --dtype float64 --chunk-size 1048576
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))
sys.path.insert(0, _HERE)

import _common  # noqa: E402
from src.common.io_utils import save_dataframe  # noqa: E402
from src.natural_gradient_local_rate.estimator_suite import compute_row  # noqa: E402
from run_operator_linearized_grid import _operator_summary, _linearized_summary  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--device", default="cuda")
    p.add_argument("--dtype", default="float64")
    p.add_argument("--chunk-size", type=int, default=None)
    args = p.parse_args()

    cfg = _common.load_config(args.config)
    opdir = os.path.join(args.outdir, "operator_grid")
    lrdir = os.path.join(args.outdir, "linearized_rate_grid")
    op_long = os.path.join(opdir, "results_long.csv")
    lr_long = os.path.join(lrdir, "results_long.csv")
    eigdir = os.path.join(lrdir, "eigenvectors")

    df = pd.read_csv(op_long)
    with open(os.path.join(opdir, "config.json")) as f:
        run_id = json.load(f).get("run_id", "")
    group = "natural_gradient_local_rate"

    bad_mask = df["status"] != "ok"
    bad = df[bad_mask]
    print(f"[patch] {len(bad)}/{len(df)} rows need re-running (run_id={run_id})")

    opts = _common.operator_opts(cfg)
    opts["compute_gamma_loc"] = True
    opts["backend"] = "torch"
    opts["device"] = args.device
    opts["dtype"] = args.dtype
    if args.chunk_size is not None:
        opts["chunk_size"] = args.chunk_size

    key_cols = ["potential_family", "N_theta", "kappa_target", "seed"]

    def _row_key(s):
        return (str(s["potential_family"]), int(s["N_theta"]),
                float(s["kappa_target"]), int(s["seed"]))

    # original column order (without the leading-underscore eigenvector blobs)
    orig_cols = [c for c in df.columns if not c.startswith("_")]
    good_rows = [r.to_dict() for _, r in df[~bad_mask].iterrows()]
    patched_rows = []

    fixed = 0
    t0 = time.time()
    for _, r in bad.iterrows():
        point = {
            "family": r["potential_family"],
            "N_theta": int(r["N_theta"]),
            "kappa_target": float(r["kappa_target"]),
            "seed": int(r["seed"]),
        }
        point["M_mc"] = _common.grid_M_mc(cfg, point)
        Z = _common.make_bank(cfg, point)
        pot = _common.make_potential_for_opts(cfg, point, Z, opts)
        row = compute_row(pot, Z, point, opts, run_id=run_id, experiment_group=group)
        if row["status"] != "ok":
            print(f"  STILL FAILED {_common.point_key(point)}: {row['error_message']}")
            patched_rows.append(r.to_dict())   # keep the old (failed) row as-is
            continue
        if "_u_star" in row:
            np.savez(os.path.join(eigdir, _common.point_key(point) + ".npz"),
                     u_star=row["_u_star"], X_star=row["_X_star"],
                     gamma_loc=row["gamma_loc"], N_theta=point["N_theta"])
        patched_rows.append({k: v for k, v in row.items() if not k.startswith("_")})
        fixed += 1
        print(f"  [{fixed}/{len(bad)}] {_common.point_key(point):30s} "
              f"sym={row['Lambda_hat_full_sym']:.4f} gamma={row['gamma_loc']:.4f} "
              f"[{row['status']}] ({time.time()-t0:.1f}s)")

    # rebuild the frame: untouched good rows + freshly computed rows, then
    # restore the canonical row order by sorting on the grid key.
    merged = pd.DataFrame(good_rows + patched_rows)
    merged["_k"] = merged.apply(_row_key, axis=1)
    orig_order = [(_row_key(r)) for _, r in df.iterrows()]
    order_index = {k: i for i, k in enumerate(orig_order)}
    merged["_ord"] = merged["_k"].map(order_index)
    merged = merged.sort_values("_ord").drop(columns=["_k", "_ord"]).reset_index(drop=True)
    df = _common.order_columns(merged)
    save_dataframe(op_long, df)
    save_dataframe(lr_long, df)
    save_dataframe(os.path.join(opdir, "summary.csv"), _operator_summary(df))
    save_dataframe(os.path.join(lrdir, "summary.csv"), _linearized_summary(df))
    remaining = int((df["status"] != "ok").sum())
    print(f"\n[patch] fixed {fixed} rows; {remaining} still not ok; total {len(df)}")


if __name__ == "__main__":
    main()
