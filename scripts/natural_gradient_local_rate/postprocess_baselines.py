"""Add Gaussian/separable baseline corrections to a natural-gradient CSV.

Usage:
    python scripts/natural_gradient_local_rate/postprocess_baselines.py \
        --input outputs/natural_gradient_local_rate/sample_size_scaling/results_long.csv \
        --out outputs/natural_gradient_local_rate/sample_size_scaling/results_long_with_baselines.csv
"""
import argparse
import os
import sys

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))

from src.common.io_utils import save_dataframe  # noqa: E402
from src.natural_gradient_local_rate.baseline_postprocessing import (  # noqa: E402
    add_baseline_corrections, aggregate_seed_summary,
)
from src.natural_gradient_local_rate.scaling_diagnostics import (  # noqa: E402
    add_noise_scale_columns, fit_scaling_diagnostics,
)


def parse_args():
    p = argparse.ArgumentParser(description="Postprocess baseline-corrected metrics.")
    p.add_argument("--input", required=True, help="Input results_long.csv")
    p.add_argument("--out", required=True, help="Output results_long_with_baselines.csv")
    p.add_argument("--summary-out", default=None,
                   help="Optional seed-aggregate summary CSV")
    p.add_argument("--convergence-out", default=None,
                   help="Optional per-seed convergence-fit summary CSV")
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input)
    corrected = add_baseline_corrections(add_noise_scale_columns(df))
    save_dataframe(args.out, corrected)
    if args.summary_out:
        save_dataframe(args.summary_out, aggregate_seed_summary(corrected))
    if args.convergence_out:
        save_dataframe(args.convergence_out, fit_scaling_diagnostics(corrected))
    print(f"Wrote baseline-corrected CSV -> {args.out}")


if __name__ == "__main__":
    main()
