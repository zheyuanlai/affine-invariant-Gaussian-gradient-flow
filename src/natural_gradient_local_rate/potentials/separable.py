"""Separable raw feature: ``Phi(theta) = sum_i phi(theta_i + c_i)``.

A coordinate-separable control family with a bounded-curvature nonlinearity
``phi`` (default ``log cosh``, so ``phi'' = sech^2 in (0, 1]``). Because
``Hess Phi`` and ``M = E[Hess Phi]`` are diagonal, there is a clean rigorous
bound ``L_A = max_i max(M_ii, 1 - M_ii)`` on ``||Hess Phi - M||_op``.
"""
from __future__ import annotations

import numpy as np

from src.natural_gradient_local_rate.potentials.base import (
    RawFeaturePotential, get_nonlinearity,
)


class SeparablePotential(RawFeaturePotential):
    def __init__(self, N_theta, seed=0, feature_scale=1.0, phi="log_cosh"):
        self.N_theta = int(N_theta)
        self.seed = int(seed)
        self.feature_scale = float(feature_scale)
        self.phi_name = phi
        self._f, self._fp, self._fpp = get_nonlinearity(phi)

        rng = np.random.default_rng(seed)
        # per-coordinate offsets break the symmetry so M is non-trivial
        self.c = self.feature_scale * rng.standard_normal(self.N_theta)

    def phi_value(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        return np.sum(self._f(T + self.c[None, :]), axis=1)

    def phi_grad(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        return self._fp(T + self.c[None, :])

    def phi_hess(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        d = self._fpp(T + self.c[None, :])           # (M, N)
        M = T.shape[0]
        H = np.zeros((M, self.N_theta, self.N_theta), dtype=np.float64)
        idx = np.arange(self.N_theta)
        H[:, idx, idx] = d
        return H

    def deterministic_LA(self, M):
        # Hess Phi - M is diagonal with entries phi''(.) - M_ii, phi'' in (0, 1].
        mu = np.clip(np.diag(M), 0.0, 1.0)
        return float(np.max(np.maximum(mu, 1.0 - mu)))

    def raw_metadata(self):
        return {
            "family": "separable",
            "phi": self.phi_name,
            "feature_scale": self.feature_scale,
            "r": self.N_theta,
        }
