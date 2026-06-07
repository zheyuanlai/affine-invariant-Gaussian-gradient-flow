"""Reference variational optimum ``a_star = (m_star, C_star)`` for a target.

Minimises ``F(m, C) = -1/2 log det C + E_{N(m,C)}[V(theta)]`` over ``m`` and
``C in SPD(d)`` with the Cholesky parameterization ``C = L L^T``, ``L`` lower
triangular with positive diagonal (stored as off-diagonal reals plus log-diagonal
``L_ii = exp(eta_i)``). The optimization is unconstrained in those coordinates.

The gradient is analytic through the Bonnet/Price identities
``d/dm E[V] = E[grad V] = -g`` and ``d/dC E[V] = 1/2 E[grad^2 V] = -H/2`` (with
``g, H`` the target's drift fields), so

    grad_m F = -g,
    G_C      = -1/2 (H + C^{-1}),
    grad_L F = 2 G_C L.

For the analytic-Gaussian target the optimum is closed-form (``m=0``,
``C=Q^{-1}``) and no optimization is run. Several deterministic restarts guard
against a bad local basin; the lowest-objective converged result wins.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg
from scipy.optimize import minimize

from src.common.spd import symmetrize


def _pack(m, L):
    d = len(m)
    rows, cols = np.tril_indices(d, k=-1)
    return np.concatenate([m, L[rows, cols], np.log(np.diag(L))])


def _unpack(params, d):
    m = params[:d].copy()
    n_off = d * (d - 1) // 2
    off = params[d:d + n_off]
    eta = params[d + n_off:]
    L = np.zeros((d, d), dtype=np.float64)
    rows, cols = np.tril_indices(d, k=-1)
    L[rows, cols] = off
    np.fill_diagonal(L, np.exp(eta))
    return m, L


def _objective_and_grad(params, target):
    d = len(target.m0)
    m, L = _unpack(params, d)
    C = symmetrize(L @ L.T)
    F = target.objective(m, C)

    g, H = target.g_H(m, C)
    C_inv = scipy.linalg.cho_solve((L, True), np.eye(d))
    grad_m = -np.asarray(g, dtype=np.float64)            # d/dm E[V] = -g
    G_C = -0.5 * (symmetrize(H) + symmetrize(C_inv))      # d/dC F
    grad_L = 2.0 * (G_C @ L)

    rows, cols = np.tril_indices(d, k=-1)
    grad_off = grad_L[rows, cols]
    grad_eta = np.diag(grad_L) * np.diag(L)               # chain rule eta -> L_ii
    return float(F), np.concatenate([grad_m, grad_off, grad_eta])


def _starts(target):
    """Deterministic restart initial points ``(m, L)``."""
    d = len(target.m0)
    m0 = np.asarray(target.m0, dtype=np.float64)
    yield np.zeros(d), np.eye(d)                           # canonical
    yield m0.copy(), np.eye(d)                             # at the IC mean
    yield np.zeros(d), 0.5 * np.eye(d)                     # tighter
    yield np.zeros(d), 2.0 * np.eye(d)                     # looser


def compute_star(target, maxiter=5000, gtol=1e-10):
    """Return ``(m_star, C_star, F_star, diagnostics)`` for ``target``.

    Uses the closed form when the target exposes :meth:`star`; otherwise runs the
    Cholesky-parameterized L-BFGS-B with deterministic restarts.
    """
    if hasattr(target, "star"):
        m_star, C_star, F_star = target.star()
        g, H = target.g_H(m_star, C_star)
        C_inv = np.linalg.inv(symmetrize(C_star))
        diag = {
            "method": "analytic",
            "converged": True,
            "grad_norm": 0.0,
            "stationarity_g_norm": float(np.linalg.norm(g)),
            "stationarity_cov_residual": float(
                np.linalg.norm(C_inv + symmetrize(H), "fro")),
            "n_restarts": 0,
        }
        return m_star, symmetrize(C_star), float(F_star), diag

    d = len(target.m0)
    best = None
    n_ok = 0
    for m_init, L_init in _starts(target):
        res = minimize(
            _objective_and_grad, _pack(m_init, L_init), args=(target,),
            method="L-BFGS-B", jac=True,
            options={"maxiter": maxiter, "gtol": gtol, "ftol": 1e-16},
        )
        n_ok += int(bool(res.success))
        if best is None or res.fun < best.fun:
            best = res

    m_star, L_star = _unpack(best.x, d)
    C_star = symmetrize(L_star @ L_star.T)
    F_star = float(best.fun)

    _, grad = _objective_and_grad(best.x, target)
    g, H = target.g_H(m_star, C_star)
    C_inv = np.linalg.inv(C_star)
    diag = {
        "method": "lbfgsb_cholesky",
        "converged": bool(best.success),
        "message": str(best.message),
        "grad_norm": float(np.linalg.norm(grad)),
        "stationarity_g_norm": float(np.linalg.norm(g)),
        "stationarity_cov_residual": float(np.linalg.norm(C_inv + symmetrize(H), "fro")),
        "n_restarts": len(list(_starts(target))),
        "n_converged_restarts": int(n_ok),
        "nit": int(best.nit),
    }
    return m_star, C_star, F_star, diag
