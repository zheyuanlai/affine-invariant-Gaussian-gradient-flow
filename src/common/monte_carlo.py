"""Reproducible Gaussian Monte Carlo sample banks.

A single fixed sample bank ``Z ~ N(0, I_{N_theta})`` is reused throughout a run
(common random numbers), so that every operator application and every grid point
that shares a seed sees identical noise. Antithetic sampling (``z`` and ``-z``)
is available to reduce variance of odd moments. All arrays are float64.
"""
from __future__ import annotations

import numpy as np

from src.common.spd import symmetric_sqrt


def gaussian_samples(N_theta, M, seed, antithetic=True):
    """Return a fixed bank of ``M`` standard-Gaussian samples, shape ``(M, N_theta)``.

    With ``antithetic=True`` the bank is ``[z_1..z_h, -z_1..-z_h]`` using
    ``h = M // 2`` independent base draws (one extra independent draw is appended
    if ``M`` is odd). The result is deterministic given ``seed``.
    """
    if M <= 0:
        raise ValueError(f"M must be positive, got {M}")
    rng = np.random.default_rng(seed)
    if antithetic:
        half = M // 2
        base = rng.standard_normal((half, N_theta))
        Z = np.concatenate([base, -base], axis=0)
        if 2 * half < M:  # odd M -> one extra independent sample
            Z = np.concatenate([Z, rng.standard_normal((M - 2 * half, N_theta))], axis=0)
    else:
        Z = rng.standard_normal((M, N_theta))
    return np.ascontiguousarray(Z, dtype=np.float64)


def transform_gaussian_samples(Z, m, C):
    """Map standard-Gaussian samples to ``N(m, C)``: ``theta_j = m + C^{1/2} z_j``.

    ``Z`` has shape ``(M, N_theta)``; returns ``Theta`` of the same shape.
    """
    Z = np.asarray(Z, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).reshape(1, -1)
    C_sqrt = symmetric_sqrt(C)
    # rows of (Z @ C_sqrt) are (C_sqrt z_j) because C_sqrt is symmetric
    return m + Z @ C_sqrt
