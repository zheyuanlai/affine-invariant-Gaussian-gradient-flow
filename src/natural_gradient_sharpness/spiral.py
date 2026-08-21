"""Gaussian-flow simulation inside the logarithmic spiral's flat valley."""
from __future__ import annotations

import math
import numpy as np

from .quadrature import normal_rule


class SpiralValleyTarget:
    """The exact flat-core dynamics of Appendix B.1's spiral target.

    The construction chooses a radius so large that every quadrature point
    remains in ``chi=1`` and ``|phi|<delta_K/2``.  In that certified tube the
    full potential, gradient, and Hessian reduce to the explicit formulas used
    below.  ``containment_ratio`` is recorded so a run cannot silently use the
    local formula outside its stated domain.
    """

    def __init__(self, K: int, alpha=1.0, radius_factor=200.0, gh_order=10):
        self.K = int(K)
        self.alpha = float(alpha)
        self.h = 2.0 * self.alpha
        self.c = math.sqrt(self.K - 2.0) / 3.0
        self.delta = 1e-4 / self.K
        self.D = 18.0 * self.h * (self.K - 1.0) / (self.K + 7.0)
        self.r0 = float(radius_factor) * self.K ** 3 / math.sqrt(self.alpha)
        self.radius_initial = math.e ** 2 * self.r0
        self.theta0 = -2.0 * self.c
        self.m0 = self.radius_initial * np.array([math.cos(self.theta0), math.sin(self.theta0)])
        sK = math.sqrt(self.K - 2.0)
        Q = np.array([
            [(self.K + 1.0) ** 2, 2.0 * sK * (self.K + 1.0)],
            [2.0 * sK * (self.K + 1.0), 6.0 * (self.K + 1.0)],
        ]) / (self.K + 7.0)
        R = np.array([
            [math.cos(self.theta0), -math.sin(self.theta0)],
            [math.sin(self.theta0), math.cos(self.theta0)],
        ])
        P0 = self.h * R @ Q @ R.T
        self.C0 = np.linalg.inv(P0)
        z, w = normal_rule(gh_order)
        zz1, zz2 = np.meshgrid(z, z, indexing="ij")
        ww1, ww2 = np.meshgrid(w, w, indexing="ij")
        self.z = np.column_stack([zz1.ravel(), zz2.ravel()])
        self.weights = (ww1 * ww2).ravel()
        self.max_normalized_phase = 0.0
        self.min_radius_over_r0 = math.inf

    @staticmethod
    def _sym(C):
        return 0.5 * (C + C.T)

    def _point_grad_hess(self, points):
        x, y = points[:, 0], points[:, 1]
        r = np.sqrt(x * x + y * y)
        theta = np.arctan2(y, x)
        phase = theta + self.c * np.log(r / self.r0)
        phase = (phase + 0.5 * math.pi) % math.pi - 0.5 * math.pi
        self.max_normalized_phase = max(
            self.max_normalized_phase, float(np.max(np.abs(phase)) / (0.5 * self.delta)))
        self.min_radius_over_r0 = min(self.min_radius_over_r0, float(np.min(r / self.r0)))

        grad_r = r * (self.h + self.D * (0.5 * phase ** 2 + 0.5 * self.c * phase))
        grad_th = 0.5 * r * self.D * phase
        er = np.column_stack([x / r, y / r])
        et = np.column_stack([-y / r, x / r])
        grad = grad_r[:, None] * er + grad_th[:, None] * et

        # Exact polar-frame Hessian on chi=1, f(phi)=phi^2/2.  Retaining the
        # phase-dependent terms is essential for transverse stability; only at
        # phi=0 does this reduce to the constant rank-one matrix in the proof.
        b11 = self.h + self.D * (
            0.5 * phase ** 2 + 1.5 * self.c * phase + 0.5 * self.c ** 2)
        b12 = self.D * (0.5 * phase + 0.5 * self.c)
        b22 = self.h + self.D * (
            0.5 * phase ** 2 + 0.5 * self.c * phase + 0.5)
        hess = np.empty((points.shape[0], 2, 2), dtype=np.float64)
        for i in range(points.shape[0]):
            R = np.column_stack([er[i], et[i]])
            hess[i] = R @ np.array([[b11[i], b12[i]], [b12[i], b22[i]]]) @ R.T
        return grad, hess

    def expectations(self, m, C):
        L = np.linalg.cholesky(self._sym(C))
        points = np.asarray(m)[None, :] + self.z @ L.T
        grad, hess = self._point_grad_hess(points)
        return np.einsum("i,ij->j", self.weights, grad), np.einsum("i,ijk->jk", self.weights, hess)

    def rhs(self, _time, state):
        m = np.asarray(state[:2], dtype=np.float64)
        C = np.array([[state[2], state[3]], [state[3], state[4]]], dtype=np.float64)
        g, H = self.expectations(m, C)
        dm = -C @ g
        dC = self._sym(C - C @ H @ C)
        return np.array([dm[0], dm[1], dC[0, 0], dC[0, 1], dC[1, 1]])

    @property
    def state0(self):
        return np.array([self.m0[0], self.m0[1], self.C0[0, 0], self.C0[0, 1], self.C0[1, 1]])
