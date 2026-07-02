"""Fisher--Rao covariance steps and the one-step Wasserstein/Bures bootstrap.

Curvature convention: ``G = E[grad V]`` and ``A = E[grad^2 V]`` (SPD for a
log-concave target). The Gaussian natural-gradient flow is ``dm/dt = -C G``,
``dC/dt = C - C A C`` and one step of size ``dt`` is:

Mean (shared)::

    m_next = m - dt C G

Fisher--Rao Riemannian/geodesic covariance update::

    C_next = C^{1/2} exp( dt (I - C^{1/2} A C^{1/2}) ) C^{1/2}

Fisher--Rao KL/Bregman covariance update::

    C_next = (1 + dt) (C^{-1} + dt A)^{-1}.

These equal the ``H_disc = -A`` discretization/WFR/STL updates
(``exp(dt(I - .)) = e^{dt} exp(dt C^{1/2}(-A)C^{1/2})`` and ``C^{-1} + dt A =
C^{-1} - dt H_disc``); the equivalence is asserted in the tests.

One-step Wasserstein/Bures bootstrap (``eta = c/beta``)::

    m_b = m - eta G
    M   = I - eta A
    C_tilde = M C M
    C_b = 1/2 ( C_tilde + 2 eta I + [C_tilde (C_tilde + 4 eta I)]^{1/2} ),

computed spectrally from ``C_tilde = Q diag(w) Q^T`` (the two factors commute):
``C_b = Q diag( 1/2 (w + 2 eta + sqrt(w (w + 4 eta))) ) Q^T``. Since ``w >= 0``
this satisfies ``lambda_min(C_b) >= eta = c/beta`` (each mapped eigenvalue is
``>= 1/2 * 2 eta``). This is the WFR Wasserstein half-step at transport size
``eta`` in the curvature convention; we implement it here so the group is
self-contained and the ``A``-convention is explicit.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg

from src.common.spd import symmetrize, symmetric_sqrt, symmetric_expm, eigh_spd

SCHEMES = ["riemannian", "kl"]


def mean_step(m, C, G, dt):
    """Shared explicit mean update ``m - dt C G``."""
    m = np.asarray(m, dtype=np.float64)
    C = symmetrize(C)
    return m - dt * (C @ np.asarray(G, dtype=np.float64))


def riemannian_cov_step(C, A, dt):
    """Riemannian/geodesic covariance update in the curvature convention.

    ``C_next = C^{1/2} exp(dt (I - C^{1/2} A C^{1/2})) C^{1/2}``.
    """
    C = symmetrize(C)
    A = symmetrize(A)
    d = C.shape[0]
    C_sqrt = symmetric_sqrt(C)
    inner = np.eye(d) - symmetrize(C_sqrt @ A @ C_sqrt)
    expA = symmetric_expm(dt * inner)
    return symmetrize(C_sqrt @ expA @ C_sqrt)


def kl_cov_step(C, A, dt):
    """KL/Bregman covariance update ``(1 + dt)(C^{-1} + dt A)^{-1}``.

    The precision ``P = C^{-1} + dt A`` is SPD for every ``dt > 0`` when ``A`` is
    PSD (log-concave target), so the update is unconditionally well posed; a
    :class:`numpy.linalg.LinAlgError` is raised only if ``P`` fails to be SPD.
    """
    C = symmetrize(C)
    A = symmetrize(A)
    d = C.shape[0]
    L = np.linalg.cholesky(C)
    C_inv = scipy.linalg.cho_solve((L, True), np.eye(d))
    P = symmetrize(C_inv + dt * A)
    LP = np.linalg.cholesky(P)
    P_inv = scipy.linalg.cho_solve((LP, True), np.eye(d))
    return symmetrize((1.0 + dt) * P_inv)


def cov_step(scheme, C, A, dt):
    """Fisher--Rao covariance update for the named ``scheme`` (``riemannian``/``kl``)."""
    if scheme == "riemannian":
        return riemannian_cov_step(C, A, dt)
    if scheme == "kl":
        return kl_cov_step(C, A, dt)
    raise ValueError(f"unknown scheme '{scheme}' (known: {SCHEMES})")


def bures_bootstrap_cov_step(C, A, eta, spd_tol=1e-12):
    """Wasserstein/Bures bootstrap covariance step ``C_tilde = (I-eta A) C (I-eta A)``.

    Returns ``C_b = 1/2 (C_tilde + 2 eta I + [C_tilde (C_tilde + 4 eta I)]^{1/2})``,
    computed spectrally from ``C_tilde``. ``lambda_min(C_b) >= eta``.
    """
    C = symmetrize(C)
    A = symmetrize(A)
    eta = float(eta)
    d = C.shape[0]
    M = np.eye(d) - eta * A
    C_tilde = symmetrize(M @ C @ M)
    w, Q = np.linalg.eigh(C_tilde)
    if w[0] < -abs(spd_tol):
        raise ValueError(
            f"bures_bootstrap_cov_step: C_tilde not PSD (min eigenvalue {w[0]:.3e})")
    w = np.clip(w, 0.0, None)
    lam_next = 0.5 * (w + 2.0 * eta + np.sqrt(w * (w + 4.0 * eta)))
    return symmetrize((Q * lam_next) @ Q.T)


def bures_bootstrap_step(m, C, G, A, eta):
    """One Wasserstein/Bures bootstrap step ``(m, C) -> (m_b, C_b)``.

    ``m_b = m - eta G`` and ``C_b`` from :func:`bures_bootstrap_cov_step`. ``G``,
    ``A`` are the expectations already evaluated at ``(m, C)``.
    """
    m_b = np.asarray(m, dtype=np.float64) - float(eta) * np.asarray(G, dtype=np.float64)
    C_b = bures_bootstrap_cov_step(C, A, eta)
    return m_b, C_b


def eig_extremes(C):
    """``(lambda_min, lambda_max, finite_ok)`` for a covariance matrix."""
    finite_ok = bool(np.all(np.isfinite(C)))
    if not finite_ok:
        return float("nan"), float("nan"), False
    w = eigh_spd(C)[0]
    return float(w[0]), float(w[-1]), True
