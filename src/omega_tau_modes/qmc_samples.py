"""
Fixed sample generation for Monte Carlo / quasi-Monte Carlo expectations.

Strategy
--------
Use a Sobol low-discrepancy sequence (scipy.stats.qmc.Sobol) transformed to
standard normal via the inverse CDF (scipy.stats.norm.ppf).  Sobol requires
K to be a power of 2 and n <= 21201; for other sizes we fall back to a seeded
NumPy normal generator.

The returned Z matrix is fixed for the lifetime of an experiment so that all
(omega, tau, init) runs share the same randomness — common random numbers for
fair comparison.

Usage
-----
    Z = make_samples(n, K, seed)          # shape (K, n), dtype float64
    Theta = push_forward(m, L, Z)         # shape (K, n)
"""
import numpy as np


def make_samples(n: int, K: int, seed: int = 0) -> np.ndarray:
    """Generate K standard-normal samples in R^n.

    Tries Sobol QMC first (better uniformity, requires K = power-of-two and
    scipy >= 1.7).  Falls back to seeded NumPy if Sobol is unavailable or n
    exceeds the Sobol dimension limit.

    Args:
        n    : dimension
        K    : number of samples
        seed : base seed for reproducibility

    Returns:
        Z : ndarray, shape (K, n), dtype float64 — i.i.d. N(0, I_n) samples
    """
    # --- Try Sobol QMC ---
    try:
        from scipy.stats.qmc import Sobol
        from scipy.stats import norm as _norm

        # Sobol supports up to ~21201 dimensions; K must be 2^k
        max_sobol_dim = 21201
        K_is_pow2 = (K & (K - 1)) == 0
        if n <= max_sobol_dim and K_is_pow2 and K >= 4:
            sampler = Sobol(d=n, scramble=True, seed=seed)
            u = sampler.random(K)                           # uniform (K, n)
            # Clip away exact 0 and 1 before norm.ppf to avoid ±inf
            u = np.clip(u, 1e-12, 1.0 - 1e-12)
            Z = _norm.ppf(u).astype(np.float64)
            return Z
    except Exception:
        pass

    # --- Fallback: seeded NumPy ---
    rng = np.random.default_rng(seed)
    return rng.standard_normal((K, n)).astype(np.float64)


def push_forward(m: np.ndarray, L: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """Map standard-normal samples Z to N(m, C) where C = L L^T.

    Args:
        m : mean, shape (n,)
        L : lower-triangular Cholesky factor, shape (n, n)
        Z : standard-normal samples, shape (K, n)

    Returns:
        Theta : shape (K, n), each row is m + L z_j
    """
    return m[np.newaxis, :] + Z @ L.T
