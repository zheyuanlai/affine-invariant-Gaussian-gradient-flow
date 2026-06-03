"""
Print and validate a cached reference optimum.

Usage:
    python scripts/check_logconcave_reference.py [--path PATH] [--n N] [--rho RHO]

If --path is not given, defaults to outputs/logconcave_grid/reference_optimum.npz.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.reference_optimum import load_reference_optimum
from src.targets import LogCoshTarget
from src.qmc_samples import make_samples


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--path", default="outputs/logconcave_grid/reference_optimum.npz")
    p.add_argument("--n",           type=int,   default=None)
    p.add_argument("--rho",         type=float, default=None)
    p.add_argument("--target-seed", type=int,   default=123, dest="target_seed")
    p.add_argument("--sample-seed", type=int,   default=2026, dest="sample_seed")
    p.add_argument("--K-ref",       type=int,   default=8192, dest="K_ref")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if not os.path.exists(args.path):
        print(f"ERROR: {args.path} not found.")
        print("Run 'python scripts/run_logconcave_grid.py' first.")
        sys.exit(1)

    ref = load_reference_optimum(args.path)

    print("=" * 50)
    print("Reference Gaussian VI Optimum")
    print("=" * 50)
    for k in ["n", "rho", "m_features", "target_seed", "K_ref",
              "F_star", "grad_m_norm", "cov_residual_norm", "converged", "optim_message"]:
        if k in ref:
            print(f"  {k:25s}: {ref[k]}")

    m_star = ref["m_star"]
    C_star = ref["C_star"]
    print(f"\n  ||m_star||              : {np.linalg.norm(m_star):.4e}")
    print(f"  C_star eigenvalues (min, max): "
          f"{np.linalg.eigvalsh(C_star).min():.4e}, {np.linalg.eigvalsh(C_star).max():.4e}")

    # Revalidate gradient if target params available
    n   = ref.get("n",   args.n)
    rho = ref.get("rho", args.rho)
    if n is not None and rho is not None:
        print(f"\n  Re-validating gradient at a_star (n={n}, rho={rho}) ...")
        tgt = LogCoshTarget(n=n, rho=rho, seed=ref.get("target_seed", args.target_seed))
        K_ref = ref.get("K_ref", args.K_ref)
        Z_ref = make_samples(n, K_ref, seed=ref.get("sample_seed", args.sample_seed) + 1
                             if "sample_seed" in ref else args.sample_seed + 1)
        from src.qmc_samples import push_forward
        import scipy.linalg
        L_star = ref["L_star"]
        Theta_ref = push_forward(m_star, L_star, Z_ref)
        g_ref = np.mean(tgt.batch_grad(Theta_ref), axis=0)
        print(f"  ||g_ref||  (should be ~0): {np.linalg.norm(g_ref):.4e}")
