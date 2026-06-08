"""Theoretical-rate benchmark utilities (Riemannian-scale stepsize grid).

This module supports the *supplementary* rate-benchmark experiment of the
discretization-stepsize group. It compares the theorem-predicted contraction
factors of the two schemes against the observed numerical contraction, on a
common *Riemannian-scale* reference stepsize

    dt_ref = 1 / (beta * lambda_max),

with the same spectral bounds the stepsize study already uses,

    lambda_min = min(lambda_0_min, 1 / beta),
    lambda_max = max(lambda_0_max, 1 / alpha).

We deliberately do NOT use the tiny KL proof stepsize to size the grid: the
KL contraction formula is *evaluated* on the common Riemannian-scale grid
``dt = c * dt_ref`` for ``c in {0.05, ..., 1.0}``. The point is to ask how
conservative the theoretical contraction factors are at the practically
meaningful (Riemannian-scale) stepsizes identified by the stability study.

All quantities are deterministic and closed-form; there is no simulation here
(see :mod:`...rate_runner` for the trajectory driver that consumes these).

Notation
--------
``q`` is a per-step contraction factor (gap multiplied by ``q`` each step);
``r = -log(q)/dt`` is the corresponding per-unit-time rate, so a trajectory
obeying ``gap_n = q^n gap_0`` decays like ``exp(-r * n dt)``.
"""
from __future__ import annotations

import math

import numpy as np

# Default floor applied to the *raw* energy gap before any log-rate computation
# (machine-precision plateaus would otherwise give -inf log-rates). The raw gap
# is always kept alongside in the CSVs.
GAP_FLOOR = 1e-16

# Fit window (as a fraction of the initial gap) for the empirical mid-regime
# log-linear rate fit; points outside this band (early transient, machine-
# precision plateau) are excluded from the fit.
FIT_REL_HI = 1e-2
FIT_REL_LO = 1e-10
FIT_MIN_POINTS = 4


# ---------------------------------------------------------------------------
# Spectral bounds and the common Riemannian-scale reference stepsize
# ---------------------------------------------------------------------------

def spectral_bounds(alpha, beta, C0):
    """Return ``(lambda_min, lambda_max)`` used by the rate benchmark.

    Same convention as the stepsize study:
    ``lambda_min = min(lambda_0_min, 1/beta)``,
    ``lambda_max = max(lambda_0_max, 1/alpha)`` with ``lambda_0_*`` the extreme
    eigenvalues of ``C0``. Requires globally-smooth constants ``(alpha, beta)``.
    """
    if alpha is None or beta is None:
        raise ValueError("spectral_bounds requires globally-smooth (alpha, beta)")
    w0 = np.linalg.eigvalsh(0.5 * (np.asarray(C0) + np.asarray(C0).T))
    lam0_min, lam0_max = float(w0[0]), float(w0[-1])
    lam_min = min(lam0_min, 1.0 / beta)
    lam_max = max(lam0_max, 1.0 / alpha)
    return float(lam_min), float(lam_max)


def dt_ref(beta, lambda_max):
    """Common Riemannian-scale reference stepsize ``1/(beta * lambda_max)``."""
    return 1.0 / (float(beta) * float(lambda_max))


# ---------------------------------------------------------------------------
# Theoretical contraction factors and per-unit-time rates
# ---------------------------------------------------------------------------

def q_riem_theory(dt, alpha, beta, lambda_min, lambda_max):
    """Riemannian theorem per-step contraction factor.

        q = 1 - alpha * lambda_min * dt * (2 - beta * lambda_max * dt).

    For ``dt in (0, 1/(beta lambda_max)]`` (i.e. ``c in (0, 1]`` of ``dt_ref``)
    the bracket ``(2 - beta lambda_max dt) >= 1 > 0``, so ``q < 1``; and the
    decrement is bounded so ``q > 0`` over the tested grid (verified in tests).
    """
    blm = beta * lambda_max * dt
    return 1.0 - alpha * lambda_min * dt * (2.0 - blm)


def kl_kappa(dt, alpha, beta, lambda_min, lambda_max):
    """KL-analysis ``kappa(dt)`` factor (the min of the two competing terms).

        kappa = lambda_min * min{ (1 + dt alpha lambda_min)/(2(1+dt)),
                                  1/(4 (1 + dt beta lambda_max)^2) }.
    """
    term1 = (1.0 + dt * alpha * lambda_min) / (2.0 * (1.0 + dt))
    term2 = 1.0 / (4.0 * (1.0 + dt * beta * lambda_max) ** 2)
    return lambda_min * min(term1, term2)


def q_kl_formula(dt, alpha, beta, lambda_min, lambda_max):
    """KL-analysis per-step contraction factor ``1 - 2 alpha kappa(dt) dt``.

    Evaluated on the *common Riemannian-scale* grid (NOT the tiny KL proof
    stepsize). This is a formal/benchmark contraction factor: it is what the
    current KL analysis would predict at these stepsizes, not a claim that the
    KL theorem has been proved for them.
    """
    kappa = kl_kappa(dt, alpha, beta, lambda_min, lambda_max)
    return 1.0 - 2.0 * alpha * kappa * dt


def per_unit_rate(q, dt):
    """Per-unit-time rate ``r = -log(q)/dt`` (NaN if ``q`` is not in ``(0, 1]``)."""
    if not (math.isfinite(q) and q > 0.0 and q <= 1.0 and dt > 0.0):
        return float("nan")
    return -math.log(q) / dt


def r_continuous(alpha, lambda_min):
    """Continuous-time lower-bound rate ``2 alpha lambda_min``."""
    return 2.0 * float(alpha) * float(lambda_min)


def q_theory_for_method(method, dt, alpha, beta, lambda_min, lambda_max):
    """Method-specific per-step contraction factor."""
    if method == "riemannian":
        return q_riem_theory(dt, alpha, beta, lambda_min, lambda_max)
    if method == "kl":
        return q_kl_formula(dt, alpha, beta, lambda_min, lambda_max)
    raise ValueError(f"unknown method '{method}'")


# ---------------------------------------------------------------------------
# Observed-rate diagnostics from a trajectory of energy gaps
# ---------------------------------------------------------------------------

def terminal_rate(initial_gap, final_gap_raw, N, dt, gap_floor=GAP_FLOOR):
    """Terminal contraction diagnostics from the first/last gap.

    Returns a dict with the floored final gap, the per-step terminal factor
    ``q_hat_terminal = (final/initial)^(1/N)`` and the per-unit-time rate
    ``r_hat_terminal = -log(final/initial)/(N dt)``, plus a ``floor_limited``
    flag set when the raw final gap is at or below ``gap_floor``.
    """
    final_for_logs = max(float(final_gap_raw), gap_floor)
    floor_limited = bool(float(final_gap_raw) <= gap_floor)
    out = {
        "final_gap_raw": float(final_gap_raw),
        "final_gap_for_logs": float(final_for_logs),
        "floor_limited_final": int(floor_limited),
        "q_hat_terminal": float("nan"),
        "r_hat_terminal": float("nan"),
    }
    if N >= 1 and initial_gap > 0.0 and final_for_logs > 0.0 and math.isfinite(final_for_logs):
        ratio = final_for_logs / initial_gap
        if ratio > 0.0:
            out["q_hat_terminal"] = float(ratio ** (1.0 / N))
            out["r_hat_terminal"] = float(-math.log(ratio) / (N * dt))
    return out


def fitted_rate(times, gaps, initial_gap, dt,
                rel_lo=FIT_REL_LO, rel_hi=FIT_REL_HI, min_points=FIT_MIN_POINTS,
                gap_floor=GAP_FLOOR):
    """Log-linear empirical rate over a robust mid-regime window.

    Fits ``log(gap_n) = a - r t`` for steps whose *raw* gap lies in
    ``[rel_lo, rel_hi] * initial_gap`` (excluding the early transient and the
    machine-precision plateau). Returns ``r_hat_fit``, ``q_hat_fit =
    exp(-r_hat_fit dt)``, ``fit_num_points`` and ``fit_window_status`` (``ok``
    or ``insufficient_points``). Falls back to NaN rates when too few points.
    """
    times = np.asarray(times, dtype=np.float64)
    gaps = np.asarray(gaps, dtype=np.float64)
    out = {"r_hat_fit": float("nan"), "q_hat_fit": float("nan"),
           "fit_num_points": 0, "fit_window_status": "insufficient_points"}
    if initial_gap <= 0.0 or times.size == 0:
        return out
    lo = rel_lo * initial_gap
    hi = rel_hi * initial_gap
    mask = (np.isfinite(gaps) & (gaps > max(lo, gap_floor)) & (gaps <= hi)
            & np.isfinite(times))
    n_pts = int(np.sum(mask))
    out["fit_num_points"] = n_pts
    if n_pts < min_points:
        return out
    t_fit = times[mask]
    log_gap = np.log(gaps[mask])
    # Least-squares slope of log_gap vs t; rate r = -slope.
    slope, _ = np.polyfit(t_fit, log_gap, 1)
    r_fit = float(-slope)
    out["r_hat_fit"] = r_fit
    out["q_hat_fit"] = float(math.exp(-r_fit * dt))
    out["fit_window_status"] = "ok"
    return out


# ---------------------------------------------------------------------------
# Iteration-to-tolerance comparison (observed vs theory)
# ---------------------------------------------------------------------------

def n_obs_to_tol(gaps, eps):
    """First step index ``n`` with ``gap_n <= eps`` (raw gaps), else ``None``.

    ``gaps`` is the full per-step raw-gap trajectory (index 0 = initial state).
    Returns ``None`` ("not reached") if the tolerance is never met.
    """
    gaps = np.asarray(gaps, dtype=np.float64)
    for n in range(gaps.size):
        if math.isfinite(gaps[n]) and gaps[n] <= eps:
            return int(n)
    return None


def n_theory_to_tol(eps, initial_gap, q_theory):
    """Theoretical iteration count ``ceil(log(eps/gap0)/log(q_theory))``.

    Returns ``None`` when the bound is undefined (``q_theory`` not a strict
    contraction in ``(0, 1)``, or ``eps >= initial_gap`` so zero steps suffice
    we still report the ceiling, which may be <= 0 -> clamp to 0).
    """
    if not (math.isfinite(q_theory) and 0.0 < q_theory < 1.0):
        return None
    if initial_gap <= 0.0 or eps <= 0.0:
        return None
    target = eps / initial_gap
    if target >= 1.0:
        return 0
    n = math.ceil(math.log(target) / math.log(q_theory))
    return int(max(n, 0))


