"""Run the operator-grid, linearized-rate-grid and flow-validation stages.

By default (or with --smoke) this uses the fast smoke config. A larger grid is
run only if you explicitly pass --config <your.yaml>; nothing expensive runs
automatically.

Usage:
    python scripts/natural_gradient_local_rate/run_all.py --smoke
    python scripts/natural_gradient_local_rate/run_all.py \
        --config configs/natural_gradient_local_rate/smoke.yaml
"""
import argparse
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))
sys.path.insert(0, _HERE)

import run_operator_grid  # noqa: E402
import run_linearized_rate_grid  # noqa: E402
import run_flow_validation  # noqa: E402

DEFAULT_CONFIG = os.path.join(_HERE, "..", "..", "configs",
                              "natural_gradient_local_rate", "smoke.yaml")


def parse_args():
    p = argparse.ArgumentParser(description="Run all natural-gradient local-rate stages.")
    p.add_argument("--config", default=DEFAULT_CONFIG)
    p.add_argument("--smoke", action="store_true",
                   help="Shortcut for --config <.../smoke.yaml>")
    p.add_argument("--outdir", default=None)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def _run(module, config, outdir, overwrite):
    argv = [module.__name__, "--config", config]
    if outdir:
        argv += ["--outdir", outdir]
    if overwrite:
        argv += ["--overwrite"]
    old = sys.argv
    sys.argv = argv
    try:
        module.main()
    finally:
        sys.argv = old


def main():
    args = parse_args()
    config = DEFAULT_CONFIG if args.smoke else args.config
    for stage, module in [
        ("operator grid", run_operator_grid),
        ("linearized rate grid", run_linearized_rate_grid),
        ("flow validation", run_flow_validation),
    ]:
        print("=" * 70)
        print(f"STAGE: {stage}   (config: {os.path.basename(config)})")
        print("=" * 70)
        _run(module, config, args.outdir, args.overwrite)
        print()


if __name__ == "__main__":
    main()
