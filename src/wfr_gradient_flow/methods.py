"""The WFR Gaussian forward--backward splitting and its two half-steps.

The Wasserstein--Fisher--Rao (WFR) Gaussian flow with transport strength
``lambda_t`` is, writing ``g = E[grad log rho_post]`` and ``H = E[grad^2 log
rho_post]``,

    dm/dt = (C + lambda I) g,
    dC/dt = C + C H C + lambda (2 I + C H + H C).

Its forward--backward discretization, with Wasserstein step size
``h_n = lambda_n Delta t`` and Fisher--Rao step size ``Delta t``, interleaves a
Wasserstein (Bures--transport) half-step and a Fisher--Rao (KL/natural-gradient)
half-step. With ``g_n, H_n`` evaluated at ``a_n = (m_n, C_n)``:

Wasserstein half-step (``h = h_n``)::

    m_{n+1/2} = m_n + h g_n
    M_n       = I + h H_n
    C_tilde   = M_n C_n M_n
    C_{n+1/2} = 1/2 ( C_tilde + 2 h I + [ C_tilde (C_tilde + 4 h I) ]^{1/2} )

Fisher--Rao half-step (``dt = Delta t``), with ``g_{n+1/2}, H_{n+1/2}`` evaluated
at ``a_{n+1/2}``::

    m_{n+1} = m_{n+1/2} + dt C_{n+1/2} g_{n+1/2}
    C_{n+1} = (1 + dt) ( C_{n+1/2}^{-1} - dt H_{n+1/2} )^{-1}.

The Fisher--Rao half-step is exactly the KL/Bregman covariance update and shared
explicit mean update of the natural-gradient discretization; we reuse that
single source of truth (:mod:`...natural_gradient_discretization_stepsize.methods`).

Cost accounting (expectation batches per iteration): ``fr_only`` and ``w_only``
evaluate the expectations once; the full WFR splitting evaluates them twice (once
at ``a_n``, once at ``a_{n+1/2}``).

Numerics. The Wasserstein covariance square root is computed spectrally: in exact
arithmetic ``C_tilde`` and ``C_tilde + 4hI`` commute, so with the symmetric
eigendecomposition ``C_tilde = V diag(w) V^T`` the closed form

    C_{n+1/2} = V diag( 1/2 ( w + 2h + sqrt(w^2 + 4 h w) ) ) V^T

is exact, symmetric, and SPD whenever ``C_tilde`` is PSD and ``h >= 0``. We never
silently repair a non-SPD ``C_tilde``; a negative ``w`` raises and the runner
records the step as failed.
"""
from __future__ import annotations

import numpy as np

from src.common.spd import symmetrize, eigh_spd
# Fisher--Rao half-step = KL covariance step + shared explicit mean step.
from src.natural_gradient_discretization_stepsize.methods import (
    kl_cov_step, mean_step,
)

METHOD_NAMES = ["fr_only", "w_only", "wfr_fixed", "wfr_theory", "wfr_adaptive"]
# Methods that take the full two-half-step WFR splitting (2 expectation batches).
WFR_METHODS = ["wfr_fixed", "wfr_theory", "wfr_adaptive"]


def _eig_diag(C):
    """Eigenvalue extremes + SPD/finite flags for a covariance matrix."""
    finite_ok = bool(np.all(np.isfinite(C)))
    if finite_ok:
        w = eigh_spd(C)[0]
        min_eig, max_eig = float(w[0]), float(w[-1])
    else:
        min_eig, max_eig = float("nan"), float("nan")
    return {
        "min_eig_C": min_eig, "max_eig_C": max_eig,
        "spd_ok": bool(finite_ok and min_eig > 0.0), "finite_ok": finite_ok,
    }


# ---------------------------------------------------------------------------
# Wasserstein (Bures-transport) half-step
# ---------------------------------------------------------------------------

def wasserstein_cov_step(C, H, h, spd_tol=1e-12):
    """Wasserstein covariance forward--backward update.

    ``C_{n+1/2} = 1/2 ( C_tilde + 2 h I + [C_tilde (C_tilde + 4hI)]^{1/2} )`` with
    ``C_tilde = (I + h H) C (I + h H)``. Computed spectrally from the symmetric
    eigendecomposition of ``C_tilde`` (the two factors commute). ``h = 0`` returns
    ``C`` unchanged (the W step is the identity). Raises ``ValueError`` if
    ``C_tilde`` fails to be PSD beyond rounding.
    """
    C = symmetrize(C)
    H = symmetrize(H)
    h = float(h)
    if h == 0.0:
        return C
    d = C.shape[0]
    M = np.eye(d, dtype=np.float64) + h * H
    C_tilde = symmetrize(M @ C @ M)
    w, V = np.linalg.eigh(C_tilde)
    if w[0] < -abs(spd_tol):
        raise ValueError(
            f"wasserstein_cov_step: C_tilde not PSD (min eigenvalue {w[0]:.3e})")
    w = np.clip(w, 0.0, None)
    # 1/2 ( w + 2h + sqrt(w (w + 4h)) ), the spectral form of the matrix update.
    lam_next = 0.5 * (w + 2.0 * h + np.sqrt(w * (w + 4.0 * h)))
    return symmetrize((V * lam_next) @ V.T)


def wasserstein_mean_step(m, g, h):
    """Wasserstein explicit mean update ``m + h g``."""
    m = np.asarray(m, dtype=np.float64)
    return m + float(h) * np.asarray(g, dtype=np.float64)


def wasserstein_step(target, m, C, g, H, h):
    """One Wasserstein half-step ``(m, C) -> (m_half, C_half)`` using given g, H.

    ``g, H`` are the expectations already evaluated at ``(m, C)`` (no new batch).
    """
    m_half = wasserstein_mean_step(m, g, h)
    C_half = wasserstein_cov_step(C, H, h)
    return m_half, symmetrize(C_half)


# ---------------------------------------------------------------------------
# Fisher--Rao (KL / natural-gradient) half-step -- reused single source of truth
# ---------------------------------------------------------------------------

def fisher_rao_step(target, m, C, g, H, dt):
    """One Fisher--Rao half-step using the KL covariance + shared mean update.

    ``m_next = m + dt C g`` and ``C_next = (1 + dt)(C^{-1} - dt H)^{-1}`` with
    ``g, H`` already evaluated at ``(m, C)``.
    """
    m_next = mean_step(m, C, g, dt)
    C_next = kl_cov_step(C, H, dt)
    return m_next, symmetrize(C_next)


# ---------------------------------------------------------------------------
# Composed one-iteration step (dispatches on the method family)
# ---------------------------------------------------------------------------

def wfr_step(method, target, m, C, g, H, h, dt):
    """One full iteration ``(m, C) -> (m_next, C_next)`` for the named method.

    ``g, H`` are the expectations *already* evaluated at ``a_n = (m, C)`` (one
    batch, shared with the schedule that picked ``h``). Returns
    ``(m_next, C_next, diag)`` where ``diag`` carries the eigenvalue extremes /
    SPD flag of ``C_next`` and ``n_batches`` (expectation batches used this
    iteration). The count is exact: ``fr_only`` / ``w_only`` use the one batch
    already spent; the full WFR splitting spends a second at the intermediate
    state ``a_{n+1/2}`` (two total), except in the degenerate ``h = 0`` case where
    the Wasserstein step is the identity and the same ``a_n`` batch drives the
    Fisher--Rao step (one batch, recovering ``fr_only``).
    """
    if method == "fr_only":
        m_next, C_next = fisher_rao_step(target, m, C, g, H, dt)
        n_batches = 1
    elif method == "w_only":
        m_next, C_next = wasserstein_step(target, m, C, g, H, h)
        n_batches = 1
    elif method in WFR_METHODS:
        m_half, C_half = wasserstein_step(target, m, C, g, H, h)
        if float(h) == 0.0:
            g2, H2, n_batches = g, H, 1   # W step is identity: reuse a_n batch
        else:
            g2, H2 = target.g_H(m_half, C_half)   # batch at a_{n+1/2}
            n_batches = 2
        m_next, C_next = fisher_rao_step(target, m_half, C_half, g2, H2, dt)
    else:
        raise ValueError(f"unknown method '{method}' (known: {METHOD_NAMES})")

    m_next = np.asarray(m_next, dtype=np.float64)
    C_next = symmetrize(C_next)
    diag = _eig_diag(C_next)
    diag["n_batches"] = int(n_batches)
    return m_next, C_next, diag
