"""Run the nonconvex Gaussian-gradient instability experiments.

The runner writes:

    results_long.csv
    summary.csv
    kl_pole_summary.csv
    wasserstein_bound_summary.csv
    target_metadata.json
    run_metadata.json

Usage:

    python scripts/natural_gradient_nonconvex_instability/run_experiments.py \
        --config configs/natural_gradient_nonconvex_instability/nonconvex_instability.yaml \
        --outdir outputs/natural_gradient_nonconvex_instability --overwrite
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.io_utils import ensure_dir, load_yaml  # noqa: E402
from src.natural_gradient_nonconvex_instability.runner import (  # noqa: E402
    run_all,
    shallow_apply_smoke,
)


DEFAULT_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "..", "configs",
    "natural_gradient_nonconvex_instability", "nonconvex_instability.yaml")


def parse_args():
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=DEFAULT_CONFIG)
    pre.add_argument("--smoke", action="store_true")
    known, _ = pre.parse_known_args()
    cfg = load_yaml(known.config)
    if known.smoke:
        cfg = shallow_apply_smoke(cfg)
    parser = argparse.ArgumentParser(parents=[pre], description=__doc__)
    parser.add_argument("--outdir", default=cfg["output_dir"])
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    return args, cfg


def main():
    args, cfg = parse_args()
    print("=" * 72)
    print("Natural-gradient nonconvex instability experiments")
    print("=" * 72)
    print(f"  config : {args.config}{'  [smoke]' if args.smoke else ''}")
    print(f"  outdir : {args.outdir}")
    print(f"  GH nodes: {cfg.get('gh_nodes', 160)}")
    print()
    if args.overwrite and os.path.isdir(args.outdir):
        for entry in os.listdir(args.outdir):
            p = os.path.join(args.outdir, entry)
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    ensure_dir(args.outdir)
    run_all(cfg, args.outdir)
    print(f"\nWrote nonconvex instability outputs -> {args.outdir}")


if __name__ == "__main__":
    main()

