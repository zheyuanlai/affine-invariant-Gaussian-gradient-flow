#!/usr/bin/env python3
"""Run all deterministic sharp-rate experiments."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.natural_gradient_sharpness.runner import run_all


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/natural_gradient_sharpness/sharp_rates.yaml")
    parser.add_argument("--outdir", default="outputs/natural_gradient_sharpness")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()) and not args.overwrite:
        parser.error(f"{outdir} is nonempty; pass --overwrite")
    with open(args.config, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    _, metadata = run_all(config, outdir)
    print(f"wrote sharp-rate outputs to {outdir}")
    for key, value in metadata["rate_summary"].items():
        print(f"{key}: {value:.6g}")


if __name__ == "__main__":
    main()
