"""Deterministic quadrature used by the sharpness constructions.

The high-dimensional ridge and shell targets are rotationally symmetric in
``m`` transverse coordinates.  Their Gaussian expectations therefore reduce
to two scalar variables: a standard normal ``T`` and ``U ~ chi2(m)``.  The
chi-square integral is performed in survival-probability coordinates.  A
tail-adapted composite Gauss--Legendre rule resolves the narrow shell without
ever constructing an ``m``-dimensional vector.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def normal_rule(order: int = 28):
    """Nodes and weights for expectation under a standard normal."""
    x, w = np.polynomial.hermite.hermgauss(int(order))
    return np.sqrt(2.0) * x, w / np.sqrt(np.pi)


def _legendre_interval(a: float, b: float, order: int):
    x, w = np.polynomial.legendre.leggauss(int(order))
    nodes = 0.5 * ((b - a) * x + (b + a))
    weights = 0.5 * (b - a) * w
    return nodes, weights


def chi2_survival_rule(m: int, order: int = 14, focus_probability: float | None = None):
    """Nodes/weights for ``E[f(U)]``, ``U ~ chi2(m)``.

    In survival coordinates ``q = P(U >= u)``, expectation is simply
    ``int_0^1 f(chi2.isf(q)) dq``.  Composite intervals are logarithmically
    refined near both endpoints and, when supplied, around the shell tail
    probability.  The Gauss--Legendre nodes never include 0 or 1.
    """
    cuts = {
        0.0, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3,
        1e-2, 5e-2, 0.2, 0.5, 0.8, 0.95, 0.99, 0.999, 0.9999,
        0.999999, 1.0 - 1e-10, 1.0,
    }
    if focus_probability is not None:
        p = float(focus_probability)
        for factor in (1e-3, 1e-2, 0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, 1000.0):
            cuts.add(float(np.clip(p * factor, 0.0, 1.0)))
    cuts = sorted(c for c in cuts if 0.0 <= c <= 1.0)
    q_parts, w_parts = [], []
    for a, b in zip(cuts[:-1], cuts[1:]):
        if b <= a:
            continue
        q, w = _legendre_interval(a, b, order)
        q_parts.append(q)
        w_parts.append(w)
    q = np.concatenate(q_parts)
    weights = np.concatenate(w_parts)
    u = stats.chi2.isf(q, int(m))
    return u, weights


class RadialGaussianRule:
    """Product quadrature for ``T ~ N(0,1)``, ``U ~ chi2(m)``."""

    def __init__(self, m: int, normal_order: int = 28, radial_order: int = 14,
                 focus_probability: float | None = None):
        self.m = int(m)
        t, wt = normal_rule(normal_order)
        u, wu = chi2_survival_rule(m, radial_order, focus_probability)
        self.t = t[None, :]
        self.u = u[:, None]
        self.weights = wu[:, None] * wt[None, :]

    def mean(self, values):
        values = np.asarray(values, dtype=np.float64)
        return float(np.sum(self.weights * values))

    @property
    def size(self):
        return int(self.weights.size)
