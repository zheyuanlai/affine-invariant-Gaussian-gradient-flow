"""Radial-tail raw feature: ``Phi(theta) = h(||theta||^2)``.

With ``h(R) = log cosh( scale * (R / N_theta - shift) )`` and ``R = ||theta||^2``:

    h'(R)    = (scale / N) * tanh(u),     u = scale * (R/N - shift)
    h''(R)   = (scale / N)^2 * sech^2(u)
    grad Phi = 2 h'(R) theta
    Hess Phi = 2 h'(R) I + 4 h''(R) theta theta^T

The shift defaults near 1 so the activation sits in the bulk of ``R/N ~ 1``
under ``theta ~ N(0, I)``; the heavy quadratic-in-``theta`` Hessian term stresses
the Gaussian tail. No clean deterministic ``L_A`` (the curvature grows with
``||theta||``), so the centering uses an empirical estimate.
"""
from __future__ import annotations

import numpy as np

from src.natural_gradient_local_rate.potentials.base import (
    RawFeaturePotential, get_nonlinearity,
)


class RadialTailPotential(RawFeaturePotential):
    def __init__(self, N_theta, seed=0, scale=1.0, shift=1.0, phi="log_cosh"):
        self.N_theta = int(N_theta)
        self.seed = int(seed)
        self.scale = float(scale)
        self.shift = float(shift)
        self.phi_name = phi
        self._f, self._fp, self._fpp = get_nonlinearity(phi)
        # radial features are isotropic; the seed is recorded for bookkeeping.

    def _u(self, R):
        return self.scale * (R / self.N_theta - self.shift)

    def phi_value(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        R = np.sum(T * T, axis=1)
        return self._f(self._u(R))

    def phi_grad(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        R = np.sum(T * T, axis=1)
        hp = (self.scale / self.N_theta) * self._fp(self._u(R))   # h'(R), (M,)
        return 2.0 * hp[:, None] * T

    def phi_hess(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        R = np.sum(T * T, axis=1)
        u = self._u(R)
        hp = (self.scale / self.N_theta) * self._fp(u)            # h'(R)
        hpp = (self.scale / self.N_theta) ** 2 * self._fpp(u)     # h''(R)
        I = np.eye(self.N_theta)
        outer = np.einsum("mi,mj->mij", T, T, optimize=True)
        return 2.0 * hp[:, None, None] * I[None, :, :] + 4.0 * hpp[:, None, None] * outer

    def raw_metadata(self):
        return {
            "family": "radial_tail",
            "phi": self.phi_name,
            "scale": self.scale,
            "shift": self.shift,
            "r": 1,
        }
