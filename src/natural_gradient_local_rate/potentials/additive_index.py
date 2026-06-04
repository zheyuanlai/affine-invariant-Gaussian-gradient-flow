"""Additive-index (ridge-sum) raw feature.

    Phi(theta) = (1/sqrt(r)) sum_{l=1}^r a_l phi(w_l^T theta + c_l)

with random unit directions ``w_l``, sign/normal coefficients ``a_l`` and random
offsets ``c_l``. Gradient and Hessian are computed in batch:

    grad Phi = W^T ( coeff * phi'(s) ),          s = W theta + c,  coeff = a/sqrt(r)
    Hess Phi = sum_l coeff_l phi''(s_l) w_l w_l^T = W^T diag(coeff * phi''(s)) W

No clean deterministic ``L_A`` is available for the coupled Hessian, so the
centering uses an empirical estimate with a safety factor.
"""
from __future__ import annotations

import numpy as np

from src.natural_gradient_local_rate.potentials.base import (
    RawFeaturePotential, get_nonlinearity,
)


class RidgeSumFeature(RawFeaturePotential):
    """Shared implementation for additive-index and dense random-feature maps."""

    family = "ridge_sum"

    def __init__(self, N_theta, r, seed=0, feature_scale=1.0, phi="log_cosh",
                 coeff_dist="signs", offset_scale=1.0):
        self.N_theta = int(N_theta)
        self.r = int(r)
        self.seed = int(seed)
        self.feature_scale = float(feature_scale)
        self.offset_scale = float(offset_scale)
        self.phi_name = phi
        self.coeff_dist = coeff_dist
        self._f, self._fp, self._fpp = get_nonlinearity(phi)

        rng = np.random.default_rng(seed)
        W = rng.standard_normal((self.r, self.N_theta))
        W /= np.linalg.norm(W, axis=1, keepdims=True)  # unit rows
        self.W = self.feature_scale * W                 # (r, N)

        if coeff_dist == "signs":
            self.a = rng.choice([-1.0, 1.0], size=self.r)
        elif coeff_dist == "normal":
            self.a = rng.standard_normal(self.r)
        else:
            raise ValueError(f"unknown coeff_dist '{coeff_dist}'")
        self.c = self.offset_scale * rng.standard_normal(self.r)
        self.coeff = self.a / np.sqrt(self.r)           # (r,)

    def _pre(self, T):
        return T @ self.W.T + self.c[None, :]           # (M, r)

    def phi_value(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        return self._f(self._pre(T)) @ self.coeff

    def phi_grad(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        s = self._pre(T)
        return (self._fp(s) * self.coeff[None, :]) @ self.W   # (M, N)

    def phi_hess(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        s = self._pre(T)
        d = self._fpp(s) * self.coeff[None, :]                # (M, r)
        # Hess[m] = W^T diag(d[m]) W
        return np.einsum("ml,li,lj->mij", d, self.W, self.W, optimize=True)

    def raw_metadata(self):
        return {
            "family": self.family,
            "phi": self.phi_name,
            "r": self.r,
            "feature_scale": self.feature_scale,
            "offset_scale": self.offset_scale,
            "coeff_dist": self.coeff_dist,
        }


class AdditiveIndexPotential(RidgeSumFeature):
    """Additive-index model: a sum of ``r`` random single-index ridge features."""

    family = "additive_index"

    def __init__(self, N_theta, r, seed=0, feature_scale=1.0, phi="log_cosh",
                 coeff_dist="signs", offset_scale=1.0):
        super().__init__(N_theta, r, seed=seed, feature_scale=feature_scale,
                         phi=phi, coeff_dist=coeff_dist, offset_scale=offset_scale)
