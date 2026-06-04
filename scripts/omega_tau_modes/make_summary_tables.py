"""
Recompute summary.csv from an existing results_long.csv.

Useful if you want to add new summary columns or change tolerance thresholds
without re-running the dynamics.

Usage:
    python scripts/make_summary_tables.py [--indir DIR] [--outdir DIR]

Defaults to outputs/omega_tau_modes/gaussian_grid/ for both in and out.
"""
import argparse
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))


TOLERANCES = {
    "time_to_1e_minus_2": 1e-2,
    "time_to_1e_minus_4": 1e-4,
    "time_to_1e_minus_6": 1e-6,
}

RUN_KEYS = ["n", "omega", "tau_type", "tau_value", "init_name", "dt", "T"]


def time_to_tol(group_df, tol):
    """First time in group where norm_energy <= tol, or inf."""
    hit = group_df[group_df["norm_energy"] <= tol]
    if hit.empty:
        return math.inf
    return float(hit["time"].iloc[0])


def summarise_group(grp):
    """Compute summary statistics for one (run) group."""
    grp = grp.sort_values("time")
    ne = grp["norm_energy"].values

    monotone = bool(np.all(np.diff(ne) <= 1e-12))

    row = {}
    for k in RUN_KEYS:
        row[k] = grp[k].iloc[0]
    row["final_energy"]            = float(grp["kl_energy"].iloc[-1])
    row["final_normalized_energy"] = float(ne[-1])
    for col, tol in TOLERANCES.items():
        row[col] = time_to_tol(grp, tol)
    row["monotone_energy_bool"]    = monotone
    row["min_eig_min_over_time"]   = float(grp["eig_min"].min())
    row["max_eig_max_over_time"]   = float(grp["eig_max"].max())
    return pd.Series(row)


def make_summary(long_csv, summary_csv):
    print(f"Reading {long_csv} ...")
    df = pd.read_csv(long_csv)

    print("Computing summary statistics ...")
    summary = df.groupby(RUN_KEYS, sort=False).apply(
        summarise_group, include_groups=False
    ).reset_index(drop=True)

    os.makedirs(os.path.dirname(summary_csv), exist_ok=True)
    summary.to_csv(summary_csv, index=False)
    print(f"Summary written -> {summary_csv}  ({len(summary)} rows)")


def parse_args():
    parser = argparse.ArgumentParser(description="Recompute summary.csv from results_long.csv.")
    parser.add_argument("--indir",  default="outputs/omega_tau_modes/gaussian_grid")
    parser.add_argument("--outdir", default=None,
                        help="Output directory (default: same as indir)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    outdir   = args.outdir or args.indir
    long_csv    = os.path.join(args.indir,  "results_long.csv")
    summary_csv = os.path.join(outdir,      "summary.csv")
    make_summary(long_csv, summary_csv)
