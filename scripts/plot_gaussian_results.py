"""
Generate all four manuscript figures from the experiment CSVs.

Figures are written to <figdir>/n{N}/ for each dimension N present in the data.

Usage:
    python scripts/plot_gaussian_results.py [--indir DIR] [--figdir DIR] [--n N [N ...]]

Defaults:
    --indir  outputs/gaussian_grid
    --figdir outputs/gaussian_grid/figures
    --n      all dimensions found in results_long.csv
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.plotting import generate_all_figures


def parse_args():
    parser = argparse.ArgumentParser(description="Generate Gaussian gradient flow figures.")
    parser.add_argument("--indir",  default="outputs/gaussian_grid",
                        help="Directory containing results_long.csv and summary.csv")
    parser.add_argument("--figdir", default=None,
                        help="Output directory for figures (default: <indir>/figures)")
    parser.add_argument("--n", type=int, nargs="+", default=None,
                        help="Dimensions to plot (default: all in CSV)")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    figdir   = args.figdir or os.path.join(args.indir, "figures")
    long_csv = os.path.join(args.indir, "results_long.csv")
    summ_csv = os.path.join(args.indir, "summary.csv")

    for path in [long_csv, summ_csv]:
        if not os.path.exists(path):
            print(f"ERROR: File not found: {path}")
            print("Run 'python scripts/run_gaussian_grid.py' first.")
            sys.exit(1)

    # Determine which n values to plot
    if args.n:
        n_list = sorted(args.n)
    else:
        n_list = sorted(pd.read_csv(long_csv, usecols=["n"])["n"].unique())

    print(f"Generating figures for n = {n_list}")
    print(f"Output root: {figdir}")
    generate_all_figures(long_csv, summ_csv, figdir, n_list=n_list)

    # Print final directory tree
    print("\nFigures written:")
    for n in n_list:
        n_dir = os.path.join(figdir, f"n{n}")
        files = sorted(f for f in os.listdir(n_dir) if f.endswith((".png", ".pdf")))
        for f in files:
            print(f"  {os.path.join(n_dir, f)}")
    print("Done.")
