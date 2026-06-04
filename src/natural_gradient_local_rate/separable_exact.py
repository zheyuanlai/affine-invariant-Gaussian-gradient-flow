"""Exact / Gauss--Hermite quadrature benchmark for *separable* control potentials.

For a separable potential ``V(theta) = sum_i V_i(theta_i)`` the Hessian is
diagonal and ``Hess V_ii`` depends only on ``theta_i``. The diagonal-mode
operator (see :func:`src.natural_gradient_local_rate.operators.diagonal_A_matrix`)
is then exactly diagonal with entries

    A_ii = E_{Z ~ N(0,1)}[ V_i''(Z) (Z^2 - 1) ],

a one-dimensional expectation that has no Monte-Carlo error when evaluated by
Gauss--Hermite quadrature. The largest such entry,

    Lambda_hat_separable_exact = max_i E[ V_i''(Z) (Z^2 - 1) ],

is the dimension-free ground truth against which the Monte-Carlo diagonal and
full-operator estimates must be compared. For non-separable potentials this
benchmark is undefined and the caller should record ``NaN`` with a status
message.

The expectation is taken against the *centered* potential's own second
derivative (evaluated through ``batch_hess`` on per-coordinate quadrature
grids), so no assumption about the centering algebra is baked in here.
"""
from __future__ import annotations

import numpy as np

from src.natural_gradient_local_rate.potentials import (
    GaussianPotential, CenteredPotential, SeparablePotential,
)


def is_separable(potential):
    """True iff ``potential`` is a coordinate-separable control potential.

    Covers the analytic Gaussian (``V = 0.5||theta||^2``) and any
    :class:`CenteredPotential` wrapping a :class:`SeparablePotential`.
    """
    if isinstance(potential, GaussianPotential):
        return True
    if isinstance(potential, CenteredPotential):
        return isinstance(potential.raw, SeparablePotential)
    return False


def _gauss_hermite_probabilists(n_nodes):
    """Nodes/weights for ``E_{Z~N(0,1)}[f(Z)] = sum_k w_k f(x_k)`` (sum w_k = 1)."""
    x, w = np.polynomial.hermite.hermgauss(int(n_nodes))   # weight exp(-x^2)
    return np.sqrt(2.0) * x, w / np.sqrt(np.pi)


def separable_exact_diagonal(potential, n_nodes=80):
    """Per-coordinate exact diagonal coefficients ``E[V_i''(Z)(Z^2 - 1)]``.

    Returns an array of length ``N_theta``. Raises ``ValueError`` if the
    potential is not separable.
    """
    if not is_separable(potential):
        raise ValueError("separable_exact_diagonal requires a separable potential")
    N = potential.N_theta
    nodes, weights = _gauss_hermite_probabilists(n_nodes)
    centered = (nodes * nodes - 1.0) * weights            # w_k (x_k^2 - 1)
    out = np.empty(N, dtype=np.float64)
    for i in range(N):
        # Evaluate Hess V at theta = node_k * e_i (other coords 0); the (i, i)
        # entry equals V_i''(node_k) because the potential is separable.
        Theta = np.zeros((nodes.size, N), dtype=np.float64)
        Theta[:, i] = nodes
        Hii = potential.batch_hess(Theta)[:, i, i]        # (n_nodes,)
        out[i] = float(np.dot(Hii, centered))
    return out


def separable_exact_lambda(potential, n_nodes=80):
    """``(value, status)`` for ``Lambda_hat_separable_exact = max_i E[V_i''(Z)(Z^2-1)]``.

    ``value`` is ``NaN`` and ``status`` explains why when the potential is not
    separable; otherwise ``status == "ok"``.
    """
    if not is_separable(potential):
        return float("nan"), "not_separable"
    diag = separable_exact_diagonal(potential, n_nodes=n_nodes)
    return float(np.max(diag)), "ok"
