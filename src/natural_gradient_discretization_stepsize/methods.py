"""The two time discretizations of the Gaussian natural gradient flow.

Both schemes share the explicit mean update and differ only in the covariance
update. With ``g = E[grad log rho_post]`` and ``H = E[grad^2 log rho_post]`` the
continuous flow is ``dm/dt = C g``, ``dC/dt = C + C H C``, and one step of size
``dt`` is:

Mean (shared)::

    m_{n+1} = m_n + dt * C_n g_n

Riemannian-distance covariance update::

    C_{n+1} = e^{dt} C_n^{1/2} exp( dt C_n^{1/2} H_n C_n^{1/2} ) C_n^{1/2}

KL / Bregman covariance update::

    C_{n+1} = (1 + dt) ( C_n^{-1} - dt H_n )^{-1}.

Every covariance is symmetrized after the update. The Riemannian step uses a
symmetric eigendecomposition for the matrix square root and exponential; the KL
step builds the precision ``C^{-1} - dt H`` and inverts it through a Cholesky
factorization (clear and stable in the 1-/2-dimensional regime studied here).
"""
from __future__ import annotations

import numpy as np
import scipy.linalg

from src.common.spd import symmetrize, symmetric_sqrt, symmetric_expm, eigh_spd

METHOD_NAMES = ["riemannian", "kl"]


def mean_step(m, C, g, dt):
    """Shared explicit mean update ``m + dt * C g``."""
    m = np.asarray(m, dtype=np.float64)
    C = symmetrize(C)
    return m + dt * (C @ np.asarray(g, dtype=np.float64))


def riemannian_cov_step(C, H, dt):
    """Riemannian-distance covariance update.

    ``C_{n+1} = e^{dt} C^{1/2} exp(dt C^{1/2} H C^{1/2}) C^{1/2}``. Returns an
    exactly symmetric SPD matrix (in exact arithmetic) built from the symmetric
    eigendecomposition.
    """
    C = symmetrize(C)
    H = symmetrize(H)
    C_sqrt = symmetric_sqrt(C)
    inner = symmetrize(C_sqrt @ H @ C_sqrt)
    expA = symmetric_expm(dt * inner)
    C_next = np.exp(dt) * (C_sqrt @ expA @ C_sqrt)
    return symmetrize(C_next)


def kl_cov_step(C, H, dt):
    """KL / Bregman covariance update ``(1 + dt)(C^{-1} - dt H)^{-1}``.

    The precision ``P = C^{-1} - dt H`` is formed explicitly (``C^{-1}`` via a
    Cholesky solve) and inverted through its Cholesky factor. For a log-concave
    target ``H`` is negative semidefinite, so ``P`` is SPD for every ``dt > 0``
    and the update is unconditionally well posed; a :class:`numpy.linalg.LinAlgError`
    is raised only if ``P`` fails to be SPD (caught by the runner as a failure).
    """
    C = symmetrize(C)
    H = symmetrize(H)
    d = C.shape[0]
    # C^{-1} via Cholesky solve (stable; avoids an explicit generic inverse).
    L = np.linalg.cholesky(C)
    C_inv = scipy.linalg.cho_solve((L, True), np.eye(d))
    P = symmetrize(C_inv - dt * H)
    # (1 + dt) P^{-1} via Cholesky of the precision.
    LP = np.linalg.cholesky(P)
    P_inv = scipy.linalg.cho_solve((LP, True), np.eye(d))
    return symmetrize((1.0 + dt) * P_inv)


def _cov_step(method, C, H, dt):
    if method == "riemannian":
        return riemannian_cov_step(C, H, dt)
    if method == "kl":
        return kl_cov_step(C, H, dt)
    raise ValueError(f"unknown method '{method}' (known: {METHOD_NAMES})")


def discretization_step(method, target, m, C, dt):
    """One full step ``(m, C) -> (m_next, C_next)`` for the named scheme.

    Returns ``(m_next, C_next, diag)`` where ``diag`` carries the eigenvalue
    extremes of ``C_next`` and an ``spd_ok`` / ``finite_ok`` flag. ``g`` and
    ``H`` are evaluated once at the current state and the *same* ``g`` drives the
    shared mean step for both schemes.
    """
    g, H = target.g_H(m, C)
    m_next = mean_step(m, C, g, dt)
    C_next = _cov_step(method, C, H, dt)

    finite_ok = bool(np.all(np.isfinite(m_next)) and np.all(np.isfinite(C_next)))
    if finite_ok:
        w = eigh_spd(C_next)[0]
        min_eig, max_eig = float(w[0]), float(w[-1])
    else:
        min_eig, max_eig = float("nan"), float("nan")
    diag = {
        "min_eig_C": min_eig,
        "max_eig_C": max_eig,
        "spd_ok": bool(finite_ok and min_eig > 0.0),
        "finite_ok": finite_ok,
    }
    return m_next, C_next, diag
