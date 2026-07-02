"""Render the covariance-bootstrap figures from the final CSVs.

Reads the CSVs in ``--outdir`` and writes paired PNG/PDF figures into
``--outdir/figures``. No dynamics are re-run. The figure builders are the shared
single source of truth in :mod:`src.natural_gradient_covariance_bootstrap.plotting`.
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.plotting_style import apply_style, save_figure
import matplotlib.pyplot as plt
from src.natural_gradient_covariance_bootstrap import plotting as P

DEFAULT_OUTDIR = "outputs/natural_gradient_covariance_bootstrap"


def _save(fig, figdir, name):
    if fig is None:
        return
    paths = save_figure(fig, os.path.join(figdir, name))
    plt.close(fig)
    print(f"  fig  {name}  ->  {os.path.basename(paths[0])}, {os.path.basename(paths[1])}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = p.parse_args()
    apply_style()
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(figdir, exist_ok=True)

    long_df = pd.read_csv(os.path.join(args.outdir, "results_long.csv"))
    scaling = pd.read_csv(os.path.join(args.outdir, "covariance_bootstrap_summary.csv"))
    bench = pd.read_csv(os.path.join(args.outdir, "contraction_benchmark.csv"))
    wboot = pd.read_csv(os.path.join(args.outdir, "wasserstein_bootstrap_summary.csv"))

    _save(P.fig_covariance_envelope(long_df, scheme="kl"), figdir,
          "fig_covariance_envelope")
    _save(P.fig_warmup_scaling(scaling), figdir, "fig_warmup_scaling")
    _save(P.fig_dynamic_contraction(bench), figdir, "fig_dynamic_contraction")
    _save(P.fig_contraction_factors(bench), figdir, "fig_contraction_factors")
    _save(P.fig_wasserstein_bootstrap(wboot), figdir, "fig_wasserstein_bootstrap")
    _save(P.fig_three_stage(long_df, scheme="kl"), figdir, "fig_three_stage")

    floor_path = os.path.join(args.outdir, "stl_floor_summary.csv")
    if os.path.exists(floor_path):
        floor = pd.read_csv(floor_path)
        _save(P.fig_stl_noise_floor(floor), figdir, "fig_stl_noise_floor")
    else:
        print("  [skip] fig_stl_noise_floor: stl_floor_summary.csv not found")
    print(f"\nWrote figures -> {figdir}")


if __name__ == "__main__":
    main()
