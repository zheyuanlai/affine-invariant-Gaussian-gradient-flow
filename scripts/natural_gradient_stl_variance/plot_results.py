"""Generate the STL-variance figures from the committed CSVs.

Reads ``estimator_variance.csv``, ``algorithm_results_long.csv`` and
``tail_noise_floor.csv`` from ``--outdir`` and writes paired PDF/PNG figures into
``<outdir>/figures``. No dynamics are re-run; this is pure post-processing that
reuses the figure builders in
:mod:`src.natural_gradient_stl_variance.plotting` (the same builders the report
asset generator uses).

Usage::

    python scripts/natural_gradient_stl_variance/plot_results.py \
        --outdir outputs/natural_gradient_stl_variance
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.plotting_style import apply_style, save_figure
from src.natural_gradient_stl_variance import plotting as P


def _savefig(fig, figs_dir, name):
    if fig is None:
        print(f"  [skip] {name} (no data)")
        return
    save_figure(fig, os.path.join(figs_dir, name))
    plt.close(fig)
    print(f"  fig  {name}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
    figs_dir = os.path.join(args.outdir, "figures")
    os.makedirs(figs_dir, exist_ok=True)
    apply_style()

    est_path = os.path.join(args.outdir, "estimator_variance.csv")
    long_path = os.path.join(args.outdir, "algorithm_results_long.csv")
    tail_path = os.path.join(args.outdir, "tail_noise_floor.csv")

    if os.path.exists(est_path):
        est_df = pd.read_csv(est_path)
        for kind in sorted(est_df.kind.unique()):
            _savefig(P.fig_estimator_variance_ratio(est_df, kind), figs_dir,
                     f"stl_estimator_variance_ratio_{kind}")
            _savefig(P.fig_variance_by_distance(est_df, kind), figs_dir,
                     f"stl_variance_by_distance_{kind}")
    else:
        print(f"  [skip] estimator figures: {est_path} not found")

    if os.path.exists(long_path):
        long_df = pd.read_csv(long_path)
        for kind in sorted(long_df.kind.unique()):
            _savefig(P.fig_algorithm_gap(long_df, kind), figs_dir,
                     f"stl_algorithm_gap_{kind}")
    else:
        print(f"  [skip] algorithm-gap figures: {long_path} not found")

    if os.path.exists(tail_path):
        tail_df = pd.read_csv(tail_path)
        _savefig(P.fig_noise_floor(tail_df), figs_dir, "stl_noise_floor")
        _savefig(P.fig_batchsize_effect(tail_df), figs_dir, "stl_batchsize_effect")
    else:
        print(f"  [skip] noise-floor figures: {tail_path} not found")

    print(f"\nWrote figures -> {figs_dir}")


if __name__ == "__main__":
    main()
