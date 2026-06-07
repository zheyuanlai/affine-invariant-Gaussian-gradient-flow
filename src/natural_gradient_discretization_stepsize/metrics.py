"""Stepsize feasibility classification and theoretical stepsize bounds.

A single run produces a trajectory of states ``(m_n, C_n)`` together with the
objective ``F_n = F(a_n)`` and the eigenvalue extremes of ``C_n``. The four
nested feasibility criteria below classify the run; the maxima of feasible
``dt`` over the stepsize grid (computed in :func:`stepsize_summary`) are the
quantities compared against the theoretical sufficient bounds.

Criteria (with ``F_star = F(a_star)``)
--------------------------------------
* **SPD-feasible**: no NaN/Inf and ``min_eig(C_n) > spd_tol`` for all ``n``.
* **Stable**: SPD-feasible, ``F_N < F_0``, and no explosion,
  ``F_n - F_star <= explosion_factor (F_0 - F_star)`` for all ``n``.
* **Monotone**: SPD-feasible and ``F_{n+1} <= F_n + monotone_tol`` for all ``n``.
* **Accurate**: SPD-feasible and the terminal discretization error against the
  ODE reference is ``<= accurate_tol``,
  ``||m_N - m(T)||_2 + ||C_N - C(T)||_F / ||C(T)||_F``.
"""
from __future__ import annotations

import math

import numpy as np

SPD_TOL = 1e-12
EXPLOSION_FACTOR = 1e3
MONOTONE_TOL = 1e-10
ACCURATE_TOL = 1e-2


def classify_run(F, energy_gap, min_eig_C, *, F0, F_star,
                 terminal_accuracy_error,
                 spd_tol=SPD_TOL, explosion_factor=EXPLOSION_FACTOR,
                 monotone_tol=MONOTONE_TOL, accurate_tol=ACCURATE_TOL):
    """Classify one full trajectory.

    Parameters
    ----------
    F, energy_gap, min_eig_C : 1-D arrays over all steps ``n = 0..N`` (full
        trajectory, not decimated).
    F0, F_star : float
        Initial objective and reference optimum objective.
    terminal_accuracy_error : float
        ``||m_N - m(T)|| + ||C_N - C(T)||_F / ||C(T)||_F`` (NaN if unavailable).

    Returns
    -------
    dict with the four booleans and supporting diagnostics.
    """
    F = np.asarray(F, dtype=np.float64)
    min_eig_C = np.asarray(min_eig_C, dtype=np.float64)

    finite_ok = bool(np.all(np.isfinite(F)) and np.all(np.isfinite(min_eig_C)))
    spd_feasible = bool(finite_ok and np.all(min_eig_C > spd_tol))

    gap0 = F0 - F_star
    gaps = F - F_star
    # Explosion guard relative to the initial gap (gap0 > 0 in every run here).
    if gap0 > 0:
        max_gap_ratio = float(np.max(gaps) / gap0) if finite_ok else math.inf
    else:
        max_gap_ratio = float(np.max(np.abs(gaps))) if finite_ok else math.inf
    no_explosion = bool(finite_ok and max_gap_ratio <= explosion_factor)
    improved = bool(finite_ok and F[-1] < F0)
    stable = bool(spd_feasible and improved and no_explosion)

    if finite_ok and F.size >= 2:
        diffs = np.diff(F)
        num_increases = int(np.sum(diffs > monotone_tol))
        max_increase = float(np.max(diffs)) if diffs.size else 0.0
    else:
        num_increases = (F.size - 1) if F.size >= 1 else 0
        max_increase = math.inf
    monotone = bool(spd_feasible and num_increases == 0)

    acc_err = float(terminal_accuracy_error)
    accurate = bool(spd_feasible and np.isfinite(acc_err) and acc_err <= accurate_tol)

    return {
        "spd_feasible": spd_feasible,
        "stable": stable,
        "monotone": monotone,
        "accurate": accurate,
        "num_energy_increases": num_increases,
        "max_energy_increase": max_increase,
        "max_gap_ratio": max_gap_ratio,
        "terminal_accuracy_error": acc_err,
        "finite_ok": finite_ok,
    }


# ---------------------------------------------------------------------------
# Theoretical sufficient stepsize bounds
# ---------------------------------------------------------------------------

def theory_stepsize_bounds(alpha, beta, C0):
    """Theoretical sufficient stepsize bounds for the two schemes.

    With ``lambda_0_min/max`` the extreme eigenvalues of ``C0`` and

        lambda_min = min(lambda_0_min, 1/beta),
        lambda_max = max(lambda_0_max, 1/alpha),

    the proof Lipschitz constants and stepsizes are

        L_Riem = beta lambda_max,                       dt_Riem = 1 / L_Riem,
        L_KL   = beta lambda_max max{1, lambda_max^3 / (2 lambda_min^3)},
                                                        dt_KL   = 1 / L_KL.

    Returns ``None`` for both stepsizes when ``alpha``/``beta`` are unavailable
    (non-globally-smooth target).
    """
    if alpha is None or beta is None:
        return {
            "theory_bound_available": False,
            "lambda_min": None, "lambda_max": None,
            "L_riem": None, "L_kl": None,
            "dt_theory_riem": None, "dt_theory_kl": None,
        }
    w0 = np.linalg.eigvalsh(0.5 * (np.asarray(C0) + np.asarray(C0).T))
    lam0_min, lam0_max = float(w0[0]), float(w0[-1])
    lam_min = min(lam0_min, 1.0 / beta)
    lam_max = max(lam0_max, 1.0 / alpha)
    L_riem = beta * lam_max
    kl_factor = max(1.0, lam_max ** 3 / (2.0 * lam_min ** 3))
    L_kl = beta * lam_max * kl_factor
    return {
        "theory_bound_available": True,
        "lambda_min": lam_min, "lambda_max": lam_max,
        "L_riem": float(L_riem), "L_kl": float(L_kl),
        "dt_theory_riem": float(1.0 / L_riem),
        "dt_theory_kl": float(1.0 / L_kl),
    }


def dt_theory_for_method(bounds, method):
    """Pick the method's theoretical stepsize from a :func:`theory_stepsize_bounds` dict."""
    if not bounds.get("theory_bound_available"):
        return None
    return bounds["dt_theory_riem"] if method == "riemannian" else bounds["dt_theory_kl"]


def max_feasible_dt(dts, flags):
    """Largest ``dt`` whose flag is True; ``nan`` if none pass.

    ``dts`` and ``flags`` are paired sequences (one entry per stepsize for a
    fixed target/lambda/method). This is a plain max over passing stepsizes — it
    does not assume monotonicity in ``dt`` (the runner records every stepsize, so
    a non-monotone feasibility boundary is captured faithfully).
    """
    passing = [float(dt) for dt, ok in zip(dts, flags) if bool(ok)]
    return max(passing) if passing else float("nan")
