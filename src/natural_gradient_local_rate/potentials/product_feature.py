"""Product-feature raw potentials for adversarial coupling searches.

This family is designed to stress the transverse mean-covariance coupling while
remaining Hessian-compatible: every Hessian comes from an explicit smooth scalar
feature

    Phi(theta) = (1/sqrt(r)) sum_l a_l tanh(w_l^T theta + c_l)
                                      tanh(v_l^T theta + d_l).

The paired directions ``w_l`` and ``v_l`` are orthonormalized pairwise. The
mixed Hessian term

    tanh'(s_l) tanh'(t_l) (w_l v_l^T + v_l w_l^T)

is the intended stressor for the first-Hermite coupling ``tau_H``; the centered
wrapper still enforces ``E[grad V]=0`` and ``E[Hess V]=I`` on the sample bank.
"""
from __future__ import annotations

import numpy as np

from src.natural_gradient_local_rate.potentials.base import RawFeaturePotential


def _tanh_derivatives(x):
    """Return ``tanh(x)``, first derivative and second derivative."""
    f = np.tanh(x)
    fp = 1.0 - f * f
    fpp = -2.0 * f * fp
    return f, fp, fpp


class ProductFeaturePotential(RawFeaturePotential):
    """Smooth product-feature stress test with pairwise orthogonal directions."""

    family = "product_feature"

    def __init__(self, N_theta, r=None, seed=0, feature_multiplier=4,
                 feature_scale=1.0, offset_scale=1.0):
        self.N_theta = int(N_theta)
        if self.N_theta < 2:
            raise ValueError("product_feature requires N_theta >= 2")
        self.r = int(r) if r is not None else int(feature_multiplier) * int(N_theta)
        self.seed = int(seed)
        self.feature_multiplier = int(feature_multiplier)
        self.feature_scale = float(feature_scale)
        self.offset_scale = float(offset_scale)

        rng = np.random.default_rng(seed)
        W = rng.standard_normal((self.r, self.N_theta))
        W /= np.linalg.norm(W, axis=1, keepdims=True)

        V = rng.standard_normal((self.r, self.N_theta))
        # Pairwise Gram-Schmidt so each product feature has a genuine transverse
        # direction. If a random draw is nearly parallel, resample that row.
        V -= np.sum(V * W, axis=1, keepdims=True) * W
        tiny = np.linalg.norm(V, axis=1) < 1e-12
        while np.any(tiny):
            V[tiny] = rng.standard_normal((int(np.sum(tiny)), self.N_theta))
            V[tiny] -= np.sum(V[tiny] * W[tiny], axis=1, keepdims=True) * W[tiny]
            tiny = np.linalg.norm(V, axis=1) < 1e-12
        V /= np.linalg.norm(V, axis=1, keepdims=True)

        self.W = self.feature_scale * W
        self.V = self.feature_scale * V
        self.a = rng.choice([-1.0, 1.0], size=self.r)
        self.c = self.offset_scale * rng.standard_normal(self.r)
        self.d = self.offset_scale * rng.standard_normal(self.r)
        self.coeff = self.a / np.sqrt(self.r)

    def _pre(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        s = T @ self.W.T + self.c[None, :]
        t = T @ self.V.T + self.d[None, :]
        return s, t

    def phi_value(self, theta_batch):
        s, t = self._pre(theta_batch)
        fs = np.tanh(s)
        ft = np.tanh(t)
        return (fs * ft) @ self.coeff

    def phi_grad(self, theta_batch):
        s, t = self._pre(theta_batch)
        fs, fps, _ = _tanh_derivatives(s)
        ft, fpt, _ = _tanh_derivatives(t)
        return ((fps * ft * self.coeff[None, :]) @ self.W
                + (fs * fpt * self.coeff[None, :]) @ self.V)

    def phi_hess(self, theta_batch):
        s, t = self._pre(theta_batch)
        fs, fps, fpps = _tanh_derivatives(s)
        ft, fpt, fppt = _tanh_derivatives(t)
        c = self.coeff[None, :]
        ww = np.einsum("ml,li,lj->mij", c * fpps * ft, self.W, self.W,
                       optimize=True)
        vv = np.einsum("ml,li,lj->mij", c * fs * fppt, self.V, self.V,
                       optimize=True)
        wv = np.einsum("ml,li,lj->mij", c * fps * fpt, self.W, self.V,
                       optimize=True)
        vw = np.einsum("ml,li,lj->mij", c * fps * fpt, self.V, self.W,
                       optimize=True)
        H = ww + vv + wv + vw
        return 0.5 * (H + np.transpose(H, (0, 2, 1)))

    def raw_metadata(self):
        return {
            "family": self.family,
            "phi": "tanh_product",
            "r": self.r,
            "feature_multiplier": self.feature_multiplier,
            "feature_scale": self.feature_scale,
            "offset_scale": self.offset_scale,
        }
