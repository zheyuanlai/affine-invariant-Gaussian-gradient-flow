"""Mean estimators and covariance updates for the STL experiment (NumPy path).

This module is the single-matrix, NumPy reference used by the tests and by the
CPU smoke runs. The batched (NumPy/torch) implementations used by the production
estimator-variance and algorithm grids live in
:mod:`src.natural_gradient_stl_variance.estimator_variance` and
:mod:`src.natural_gradient_stl_variance.algorithm`; they implement the *same*
formulas through the :mod:`src.natural_gradient_stl_variance.linalg` backend.

Mean estimators
---------------
For a Gaussian state ``a = (m, C)`` and a sample ``theta ~ N(m, C)``:

* baseline:  ``b_base(theta) = score_post(theta)``
* STL:       ``b_stl(theta)  = score_post(theta) + C^{-1}(theta - m)``.

``E_q[C^{-1}(theta - m)] = 0`` so the STL estimator is unbiased for
``E_q[score_post]``; at a well-specified Gaussian optimum it is pointwise zero.

Covariance updates
------------------
With the one-sample Hessian estimator ``S = Hess_log_post(theta)`` the covariance
update is held *fixed* across baseline and STL (the experiment isolates the mean
block). Two schemes:

* Riemannian:  ``C' = e^{dt} C^{1/2} exp(dt C^{1/2} S C^{1/2}) C^{1/2}``
* KL/Bregman:  ``C' = (1 + dt) (C^{-1} - dt S)^{-1}``.

Both have an exactly equivalent *Hessian-residual* form with ``K = S + C^{-1}``:

* Riemannian:  ``C' = C^{1/2} exp(dt C^{1/2} K C^{1/2}) C^{1/2}``
* KL/Bregman:  ``C' = (1 + dt) ((1 + dt) C^{-1} - dt K)^{-1}``,

verified numerically in the tests. This is why the Hessian-residual rewrite is
not a genuinely different estimator for the direct Hessian sampler used here.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg

from src.common.spd import (
    symmetrize, symmetric_sqrt, symmetric_expm,
)

METHODS = ["riemannian", "riemannian_stl", "kl", "kl_stl"]
SCHEMES = ["riemannian", "kl"]


def method_scheme_stl(method):
    """Split a method name into ``(scheme, uses_stl)``."""
    if method == "riemannian":
        return "riemannian", False
    if method == "riemannian_stl":
        return "riemannian", True
    if method == "kl":
        return "kl", False
    if method == "kl_stl":
        return "kl", True
    raise ValueError(f"unknown method '{method}' (known: {METHODS})")


# ---------------------------------------------------------------------------
# Mean estimators (batched over a leading sample axis)
# ---------------------------------------------------------------------------

def stl_correction(theta, m, C):
    """STL score-residual correction ``C^{-1}(theta - m)`` (batched over rows).

    ``theta`` has shape ``(..., d)``; ``m`` shape ``(d,)``; ``C`` shape ``(d, d)``.
    Returns the same shape as ``theta``.
    """
    theta = np.asarray(theta, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64)
    C = symmetrize(C)
    d = C.shape[0]
    L = np.linalg.cholesky(C)
    C_inv = scipy.linalg.cho_solve((L, True), np.eye(d))
    # (theta - m) C^{-1}^T == C^{-1}(theta - m) row-wise since C^{-1} is symmetric.
    return (theta - m) @ C_inv


def baseline_mean_estimator(score_vals):
    """Baseline mean estimator ``b_base = score_post(theta)`` (identity wrapper)."""
    return np.asarray(score_vals, dtype=np.float64)


def stl_mean_estimator(score_vals, theta, m, C):
    """STL mean estimator ``b_stl = score_post(theta) + C^{-1}(theta - m)``."""
    return np.asarray(score_vals, dtype=np.float64) + stl_correction(theta, m, C)


# ---------------------------------------------------------------------------
# Covariance updates -- direct form (sampled Hessian S)
# ---------------------------------------------------------------------------

def riemannian_cov_step(C, S, dt):
    """Riemannian covariance update ``e^{dt} C^{1/2} exp(dt C^{1/2} S C^{1/2}) C^{1/2}``."""
    C = symmetrize(C)
    S = symmetrize(S)
    C_sqrt = symmetric_sqrt(C)
    inner = symmetrize(C_sqrt @ S @ C_sqrt)
    expA = symmetric_expm(dt * inner)
    return symmetrize(np.exp(dt) * (C_sqrt @ expA @ C_sqrt))


def kl_cov_step(C, S, dt):
    """KL/Bregman covariance update ``(1 + dt)(C^{-1} - dt S)^{-1}``."""
    C = symmetrize(C)
    S = symmetrize(S)
    d = C.shape[0]
    L = np.linalg.cholesky(C)
    C_inv = scipy.linalg.cho_solve((L, True), np.eye(d))
    P = symmetrize(C_inv - dt * S)
    LP = np.linalg.cholesky(P)
    P_inv = scipy.linalg.cho_solve((LP, True), np.eye(d))
    return symmetrize((1.0 + dt) * P_inv)


# ---------------------------------------------------------------------------
# Covariance updates -- Hessian-residual form (K = S + C^{-1})
# ---------------------------------------------------------------------------

def riemannian_cov_step_residual(C, K, dt):
    """Riemannian update in residual form ``C^{1/2} exp(dt C^{1/2} K C^{1/2}) C^{1/2}``."""
    C = symmetrize(C)
    K = symmetrize(K)
    C_sqrt = symmetric_sqrt(C)
    inner = symmetrize(C_sqrt @ K @ C_sqrt)
    expA = symmetric_expm(dt * inner)
    return symmetrize(C_sqrt @ expA @ C_sqrt)


def kl_cov_step_residual(C, K, dt):
    """KL update in residual form ``(1 + dt)((1 + dt) C^{-1} - dt K)^{-1}``."""
    C = symmetrize(C)
    K = symmetrize(K)
    d = C.shape[0]
    L = np.linalg.cholesky(C)
    C_inv = scipy.linalg.cho_solve((L, True), np.eye(d))
    P = symmetrize((1.0 + dt) * C_inv - dt * K)
    LP = np.linalg.cholesky(P)
    P_inv = scipy.linalg.cho_solve((LP, True), np.eye(d))
    return symmetrize((1.0 + dt) * P_inv)


def cov_step(scheme, C, S, dt):
    """Direct covariance update for the named ``scheme`` (``riemannian`` / ``kl``)."""
    if scheme == "riemannian":
        return riemannian_cov_step(C, S, dt)
    if scheme == "kl":
        return kl_cov_step(C, S, dt)
    raise ValueError(f"unknown scheme '{scheme}' (known: {SCHEMES})")


def residual_from_sampled(C, S):
    """Hessian-residual sample ``K = S + C^{-1}`` (the STL covariance rewrite)."""
    C = symmetrize(C)
    d = C.shape[0]
    L = np.linalg.cholesky(C)
    C_inv = scipy.linalg.cho_solve((L, True), np.eye(d))
    return symmetrize(np.asarray(S, dtype=np.float64) + C_inv)
