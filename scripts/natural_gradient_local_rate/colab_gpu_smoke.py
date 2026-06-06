"""Colab/A100 GPU smoke helper for the natural-gradient local-rate experiment.

Prints the torch/CUDA environment, then runs the GPU smoke config on the
resolved device. Intended to be the first thing you run in a Colab session after
installing torch, to confirm the GPU path works before launching a production
sweep.

Usage (in Colab, repo root on sys.path):
    python scripts/natural_gradient_local_rate/colab_gpu_smoke.py            # device auto
    python scripts/natural_gradient_local_rate/colab_gpu_smoke.py --device cuda
"""
import argparse
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))
sys.path.insert(0, _HERE)

import _common  # noqa: E402
import run_sample_size_scaling  # noqa: E402
from src.common.torch_utils import torch_device_info  # noqa: E402

GPU_SMOKE_CONFIG = os.path.join(_HERE, "..", "..", "configs",
                                "natural_gradient_local_rate", "gpu_smoke.yaml")


def main():
    p = argparse.ArgumentParser(description="Colab GPU smoke for natural_gradient_local_rate.")
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--dtype", choices=["float64", "float32"], default="float64")
    p.add_argument("--config", default=GPU_SMOKE_CONFIG)
    p.add_argument("--outdir", default="outputs/natural_gradient_local_rate/gpu_smoke")
    args = p.parse_args()

    info = torch_device_info(args.device, args.dtype)
    print("=== torch device info ===")
    for k, v in info.items():
        print(f"  {k}: {v}")
    if not info["torch_available"]:
        print("\nPyTorch is not installed. See requirements-gpu.txt and the "
              "README 'GPU backend' section for setup.")
        sys.exit(1)
    print()

    # Delegate to the real runner with backend=torch and the requested device.
    sys.argv = [
        "run_sample_size_scaling.py",
        "--config", args.config,
        "--backend", "torch",
        "--device", args.device,
        "--dtype", args.dtype,
        "--outdir", args.outdir,
        "--overwrite",
    ]
    run_sample_size_scaling.main()


if __name__ == "__main__":
    main()
