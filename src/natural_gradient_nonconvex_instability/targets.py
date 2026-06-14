"""One-dimensional non-log-concave target and deterministic expectations.

The target is

    V_R(x) = 0.5 x^2 - 2 R^2 log cosh(x / R),

with derivatives

    V'_R(x)  = x - 2 R tanh(x / R),
    V''_R(x) = 1 - 2 sech^2(x / R) = -1 + 2 tanh(x / R)^2.

All Gaussian expectations are computed by fixed Gauss-Hermite quadrature.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.polynomial.hermite import hermgauss


def stable_log_cosh(z):
    """Stable elementwise ``log(cosh(z))``."""
    z = np.asarray(z, dtype=np.float64)
    return np.logaddexp(z, -z) - math.log(2.0)


def sech2_from_tanh(z):
    """Compute ``sech(z)^2`` through ``tanh`` to avoid overflow."""
    t = np.tanh(np.asarray(z, dtype=np.float64))
    return 1.0 - t * t


def mode_location_unit():
    """Positive solution of ``y = 2 tanh(y)`` for deterministic starts."""
    # The fixed value is enough for starts and avoids a scipy import at module
    # import time. The absolute error is below 1e-12.
    return 1.9150080481545375


@dataclass(frozen=True)
class ExpectationResult:
    V: float
    Vp: float
    Vpp: float


class NonconvexLogCoshTarget:
    """Deterministic Gauss-Hermite expectations for ``V_R``."""

    def __init__(self, R: float, n_nodes: int = 160):
        if R <= 0:
            raise ValueError("R must be positive")
        if n_nodes < 8:
            raise ValueError("n_nodes must be at least 8")
        self.R = float(R)
        self.n_nodes = int(n_nodes)
        nodes, weights = hermgauss(self.n_nodes)
        self.nodes = nodes.astype(np.float64)
        self.weights = (weights / math.sqrt(math.pi)).astype(np.float64)

    def points(self, m: float, c: float):
        """Quadrature support for ``X ~ N(m, c)``."""
        if c <= 0:
            raise ValueError("c must be positive")
        return float(m) + math.sqrt(2.0 * float(c)) * self.nodes

    def potential(self, x):
        x = np.asarray(x, dtype=np.float64)
        z = x / self.R
        return 0.5 * x * x - 2.0 * self.R * self.R * stable_log_cosh(z)

    def grad(self, x):
        x = np.asarray(x, dtype=np.float64)
        return x - 2.0 * self.R * np.tanh(x / self.R)

    def hess(self, x):
        z = np.asarray(x, dtype=np.float64) / self.R
        t = np.tanh(z)
        return -1.0 + 2.0 * t * t

    def expectations(self, m: float, c: float) -> ExpectationResult:
        x = self.points(m, c)
        w = self.weights
        return ExpectationResult(
            V=float(np.dot(w, self.potential(x))),
            Vp=float(np.dot(w, self.grad(x))),
            Vpp=float(np.dot(w, self.hess(x))),
        )

    def A(self, m: float, c: float) -> float:
        """``A_R(m,c) = E[V''_R(X)]``."""
        return self.expectations(m, c).Vpp

    def A0(self, c: float) -> float:
        """Symmetric-core Hessian expectation ``A_R(0,c)``."""
        return self.A(0.0, c)

    def energy(self, m: float, c: float) -> float:
        """Unnormalized Gaussian variational objective ``E[V] - 0.5 log c``."""
        if c <= 0:
            return float("inf")
        ev = self.expectations(m, c).V
        return float(ev - 0.5 * math.log(c))

    def objective_and_grad_m_logc(self, z):
        """Objective and gradient in variables ``(m, log c)``.

        Bonnet and Price identities give
        ``d_m E[V] = E[V']`` and ``d_c E[V] = 0.5 E[V'']``. Therefore
        ``d_{log c} F = 0.5 * (c A - 1)``.
        """
        m = float(z[0])
        log_c = float(z[1])
        c = math.exp(log_c)
        ex = self.expectations(m, c)
        f = ex.V - 0.5 * log_c
        grad = np.array([ex.Vp, 0.5 * (c * ex.Vpp - 1.0)], dtype=np.float64)
        return float(f), grad

    def metadata(self):
        return {
            "target_name": "nonconvex_logcosh",
            "R": self.R,
            "gh_nodes": self.n_nodes,
            "formula": "0.5*x^2 - 2*R^2*log(cosh(x/R))",
            "hessian_lower_bound": -1.0,
            "hessian_upper_bound": 1.0,
            "mode_location_unit": mode_location_unit(),
        }


def fr_gradient_sq(c: float, A: float) -> float:
    return float((1.0 - float(c) * float(A)) ** 2)


def bw_gradient_sq(c: float, A: float) -> float:
    return float(fr_gradient_sq(c, A) / float(c))

