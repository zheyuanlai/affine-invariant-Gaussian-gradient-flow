"""
Generate all five log-concave experiment figures.

Usage:
    python scripts/plot_logconcave_results.py [options]

Options:
    --results   path to results_long.csv  (default: outputs/logconcave_grid/results_long.csv)
    --summary   path to summary.csv       (default: outputs/logconcave_grid/summary.csv)
    --outdir    figure output directory   (default: <results_dir>/figures)
    --n         dimension to plot         (default: from CSV)
    --rho       rho value to plot         (default: from CSV)
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.lc_plotting import generate_lc_figures


def parse_args():
    p = argparse.ArgumentParser(description="Plot log-concave gradient flow results.")
    p.add_argument("--results", default="outputs/logconcave_grid/results_long.csv")
    p.add_argument("--summary", default="outputs/logconcave_grid/summary.csv")
    p.add_argument("--outdir",  default=None,
                   help="Figure output directory (default: same dir as --results + /figures)")
    p.add_argument("--n",   type=int,   default=None)
    p.add_argument("--rho", type=float, default=None)
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    for path in [args.results, args.summary]:
        if not os.path.exists(path):
            print(f"ERROR: File not found: {path}")
            print("Run 'python scripts/run_logconcave_grid.py' first.")
            sys.exit(1)

    figdir = args.outdir or os.path.join(os.path.dirname(args.results), "figures")

    # Determine n and rho from CLI or from CSV
    df = pd.read_csv(args.summary)
    n   = args.n   if args.n   is not None else int(df["n"].iloc[0])
    rho = args.rho if args.rho is not None else float(df["rho"].iloc[0])

    print(f"Generating log-concave figures for n={n}, rho={rho}")
    print(f"Output: {figdir}")

    generate_lc_figures(args.results, args.summary, figdir, n=n, rho=rho)

    print("\nFigures written:")
    for f in sorted(os.listdir(figdir)):
        if f.endswith((".png", ".pdf")):
            print(f"  {os.path.join(figdir, f)}")
    print("Done.")
