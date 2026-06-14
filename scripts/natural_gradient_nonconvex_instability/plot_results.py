"""Generate figures and a LaTeX table for the nonconvex-instability outputs.

Usage:

    python scripts/natural_gradient_nonconvex_instability/plot_results.py \
        --outdir outputs/natural_gradient_nonconvex_instability
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.plotting_style import apply_style  # noqa: E402
from src.natural_gradient_nonconvex_instability.plotting import (  # noqa: E402
    build_all_figures,
    summary_table_tex,
)


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="outputs/natural_gradient_nonconvex_instability")
    args = parser.parse_args()
    apply_style()
    figs_dir = os.path.join(args.outdir, "figures")
    long_df = pd.read_csv(os.path.join(args.outdir, "results_long.csv"))
    summary_df = pd.read_csv(os.path.join(args.outdir, "summary.csv"))
    kl_df = pd.read_csv(os.path.join(args.outdir, "kl_pole_summary.csv"))
    bw_df = pd.read_csv(os.path.join(args.outdir, "wasserstein_bound_summary.csv"))
    build_all_figures(long_df, kl_df, bw_df, figs_dir)
    _write(os.path.join(args.outdir, "tab_nonconvex_summary.tex"),
           summary_table_tex(summary_df, kl_df, bw_df))
    print(f"Wrote figures and tab_nonconvex_summary.tex under {args.outdir}")


if __name__ == "__main__":
    main()

