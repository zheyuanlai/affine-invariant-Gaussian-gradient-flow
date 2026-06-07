"""High-accuracy continuous-time reference for the natural gradient flow.

Integrates ``dm/dt = C g(a)``, ``dC/dt = C + C H(a) C`` from the target's initial
condition with a tight adaptive solver (``scipy.integrate.solve_ivp``). The flat
state is ``[m (d entries), vech(C)]`` storing the upper triangle of the symmetric
``C``; ``C`` is re-symmetrized whenever the vector field is evaluated.

The ODE reference is used only for the discretization-*accuracy* metric and the
overlaid reference curves in the figures. The stability conclusions never depend
on it.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from src.common.spd import symmetrize


def _triu_idx(d):
    return np.triu_indices(d)


def pack_state(m, C):
    """Flatten ``(m, C)`` to ``[m, upper-triangle(C)]``."""
    m = np.asarray(m, dtype=np.float64).ravel()
    C = symmetrize(C)
    r, c = _triu_idx(C.shape[0])
    return np.concatenate([m, C[r, c]])


def unpack_state(y, d):
    """Inverse of :func:`pack_state`; returns ``(m, C)`` with ``C`` symmetric."""
    y = np.asarray(y, dtype=np.float64)
    m = y[:d].copy()
    r, c = _triu_idx(d)
    C = np.zeros((d, d), dtype=np.float64)
    C[r, c] = y[d:]
    C[c, r] = y[d:]
    return m, C


def integrate_reference(target, T, rtol=1e-10, atol=1e-12, n_eval=400, max_step=None):
    """Integrate the flow to time ``T`` from ``(target.m0, target.C0)``.

    Returns a dict with the dense-output sample arrays (``t``, ``m``, ``C``), the
    terminal state ``(m_T, C_T)``, a dense interpolant ``eval`` (a callable
    ``t -> (m, C)`` valid on ``[0, T]``), and the solver diagnostics.
    """
    m0 = np.asarray(target.m0, dtype=np.float64).ravel()
    C0 = symmetrize(target.C0)
    d = m0.size

    def rhs(t, y):
        m, C = unpack_state(y, d)
        C = symmetrize(C)
        g, H = target.g_H(m, C)
        dm = C @ g
        dC = symmetrize(C + C @ H @ C)
        r, c = _triu_idx(d)
        return np.concatenate([dm, dC[r, c]])

    t_eval = np.linspace(0.0, T, int(n_eval))
    kwargs = dict(method="Radau", rtol=rtol, atol=atol, t_eval=t_eval, dense_output=True)
    if max_step is not None:
        kwargs["max_step"] = float(max_step)
    sol = solve_ivp(rhs, (0.0, T), pack_state(m0, C0), **kwargs)

    ms = np.empty((sol.t.size, d), dtype=np.float64)
    Cs = np.empty((sol.t.size, d, d), dtype=np.float64)
    for k in range(sol.t.size):
        mk, Ck = unpack_state(sol.y[:, k], d)
        ms[k] = mk
        Cs[k] = Ck
    m_T, C_T = unpack_state(sol.y[:, -1], d)

    T_end = float(T)

    def eval_state(t):
        """Dense interpolant of the reference state, clamped to ``[0, T]``."""
        tc = float(min(max(t, 0.0), T_end))
        return unpack_state(sol.sol(tc), d)

    return {
        "t": sol.t.copy(),
        "m": ms,
        "C": Cs,
        "m_T": m_T,
        "C_T": symmetrize(C_T),
        "eval": eval_state,
        "success": bool(sol.success),
        "status": int(sol.status),
        "message": str(sol.message),
        "nfev": int(sol.nfev),
        "njev": int(getattr(sol, "njev", 0)),
        "rtol": float(rtol),
        "atol": float(atol),
        "method": "Radau",
        "T": float(T),
    }
