"""Riemannian exponential-map discretization of the Gaussian natural gradient flow.

With ``g(a) = -E_{rho_a}[grad V]`` and ``H(a) = -E_{rho_a}[Hess V]`` the flow is

    dm/dt = C g(a),     dC/dt = C + C H(a) C,

and the SPD-preserving exponential-map step is

    m_{n+1} = m_n + Delta_t * C_n g_n
    C_{n+1} = C_n^{1/2} exp( Delta_t * [ I + C_n^{1/2} H_n C_n^{1/2} ] ) C_n^{1/2}.

This is a *validation* of the linearized local rate ``gamma_loc``: started as a
small perturbation along the slow eigenvector of ``L_star``, the squared distance
``R^2 = ||m||^2 + 0.5 ||C - I||_F^2`` decays like ``exp(-2 gamma_loc t)``, so a
linear fit of ``log R^2`` vs ``t`` has slope ``-2 gamma_loc``.
"""
from __future__ import annotations

import numpy as np

from src.common.spd import symmetrize, symmetric_sqrt, symmetric_expm, eigh_spd
from src.common.monte_carlo import transform_gaussian_samples


def _expectations(potential, m, C, Z, chunk_size=None):
    """Return ``g = -E[grad V]`` and ``H = -E[Hess V]`` under ``N(m, C)``.

    Accumulated in chunks so the ``(M, N, N)`` Hessian batch is never fully
    materialized for large problems.
    """
    Theta = transform_gaussian_samples(Z, m, C)
    M, N = Theta.shape
    cs = int(chunk_size) if chunk_size else M
    g_acc = np.zeros(N)
    H_acc = np.zeros((N, N))
    for s in range(0, M, cs):
        Tc = Theta[s:s + cs]
        g_acc += potential.batch_grad(Tc).sum(axis=0)
        H_acc += potential.batch_hess(Tc).sum(axis=0)
    g = -g_acc / M
    H = -symmetrize(H_acc / M)
    return g, H


def natural_gradient_vector_field(potential, m, C, Z, chunk_size=None):
    """Return ``(dm, dC, g, H)`` for the natural gradient flow at ``(m, C)``."""
    g, H = _expectations(potential, m, C, Z, chunk_size=chunk_size)
    dm = C @ g
    dC = C + C @ H @ C
    return dm, symmetrize(dC), g, H


def riemannian_step(potential, m, C, Z, Delta_t, chunk_size=None):
    """One exponential-map natural-gradient step; returns ``(m_next, C_next, diag)``."""
    g, H = _expectations(potential, m, C, Z, chunk_size=chunk_size)
    N = C.shape[0]
    C_sqrt = symmetric_sqrt(C)
    m_next = m + Delta_t * (C @ g)
    A = np.eye(N) + C_sqrt @ H @ C_sqrt
    expA = symmetric_expm(Delta_t * symmetrize(A))
    C_next = symmetrize(C_sqrt @ expA @ C_sqrt)
    w = eigh_spd(C_next)[0]
    diag = {
        "min_eig_C": float(w[0]),
        "max_eig_C": float(w[-1]),
        "spd_ok": bool(w[0] > 0.0),
    }
    return m_next, C_next, diag


def _r_squared(m, C):
    """``R^2 = ||m||^2 + 0.5 ||C - I||_F^2`` (Fisher--Rao squared distance to a_star)."""
    N = C.shape[0]
    D = C - np.eye(N)
    return float(m @ m + 0.5 * np.sum(D * D))


def _linear_fit(t, y):
    """Least-squares fit ``y = a + b t``; return ``(slope, intercept, r2)``."""
    t = np.asarray(t, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    A = np.column_stack([np.ones_like(t), t])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    intercept, slope = float(coef[0]), float(coef[1])
    resid = y - (intercept + slope * t)
    ss_res = float(resid @ resid)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return slope, intercept, r2


def _make_spd_initial(X_star, epsilon, floor=0.25):
    """Return ``(C0, epsilon_used, warning)`` with ``C0 = I + epsilon X_star`` SPD."""
    N = X_star.shape[0]
    I = np.eye(N)
    lam_min = float(np.linalg.eigvalsh(symmetrize(X_star))[0])
    eps = float(epsilon)
    warning = None
    if 1.0 + eps * lam_min < floor:
        # shrink epsilon so that min eig(C0) >= floor
        eps_new = (1.0 - floor) / abs(lam_min) if lam_min < 0 else eps
        eps_new = min(eps, max(eps_new, 0.0))
        warning = (f"reduced epsilon {eps:.3e} -> {eps_new:.3e} to keep C0 SPD "
                   f"(min eig X_star={lam_min:.3e})")
        eps = eps_new
    C0 = symmetrize(I + eps * X_star)
    return C0, eps, warning


def run_flow_validation(potential, Z, u_star, X_star, epsilon, Delta_t, n_steps,
                        fit_window=(0.05, 0.5), chunk_size=None):
    """Simulate the flow from a small slow-mode perturbation and fit the decay rate.

    Returns a dict with ``trajectory`` (list of per-step records) and ``summary``
    (the fitted slope / rate and bookkeeping). ``fit_window`` is the
    ``(start_fraction, end_fraction)`` of the time horizon used for the fit.
    """
    N = potential.N_theta
    u_star = np.asarray(u_star, dtype=np.float64)
    X_star = symmetrize(X_star)

    C0, eps_used, warn = _make_spd_initial(X_star, epsilon)
    m = eps_used * u_star
    C = C0

    traj = []
    for n in range(n_steps + 1):
        t = n * Delta_t
        w = np.linalg.eigvalsh(symmetrize(C))
        R2 = _r_squared(m, C)
        traj.append({
            "step": n,
            "t": t,
            "R2": R2,
            "log_R2": float(np.log(R2)) if R2 > 0 else float("-inf"),
            "min_eig_C": float(w[0]),
            "max_eig_C": float(w[-1]),
            "norm_m": float(np.linalg.norm(m)),
            "norm_C_minus_I": float(np.linalg.norm(C - np.eye(N))),
        })
        if n < n_steps:
            m, C, _ = riemannian_step(potential, m, C, Z, Delta_t, chunk_size=chunk_size)

    # Fit log R^2 vs t over the requested window.
    t_all = np.array([r["t"] for r in traj])
    logR_all = np.array([r["log_R2"] for r in traj])
    t_max = t_all[-1] if t_all[-1] > 0 else 1.0
    lo, hi = fit_window
    mask = (t_all >= lo * t_max) & (t_all <= hi * t_max) & np.isfinite(logR_all)
    if mask.sum() >= 2:
        slope, intercept, r2 = _linear_fit(t_all[mask], logR_all[mask])
    else:
        slope, intercept, r2 = float("nan"), float("nan"), float("nan")

    summary = {
        "epsilon_requested": float(epsilon),
        "epsilon_used": float(eps_used),
        "Delta_t": float(Delta_t),
        "n_steps": int(n_steps),
        "fit_start_fraction": float(lo),
        "fit_end_fraction": float(hi),
        "fit_slope_log_R2": slope,
        "fit_gamma_flow": -0.5 * slope,
        "fit_r2_flow": r2,
        "warning": warn,
    }
    return {"trajectory": traj, "summary": summary}
