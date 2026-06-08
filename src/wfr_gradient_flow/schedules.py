"""Wasserstein transport-strength schedules ``h_n`` for the WFR splitting.

The Wasserstein step size is ``h_n = lambda_n Delta t``; we schedule ``h_n``
directly and recover ``lambda_n = h_n / Delta t`` only for reporting. Each
schedule is a small stateful object exposing ``h(n, m, C, g, H)`` so that the
adaptive schedule can read the local curvature ``H`` already evaluated at the
current state (no extra expectation batch). The five methods of the study are:

* ``fr_only``     -- ``h_n = 0`` (the Fisher--Rao / natural-gradient scheme).
* ``w_only``      -- Wasserstein-only forward--backward step, ``h_n = c / beta``.
* ``wfr_fixed``   -- full WFR splitting with a fixed ``h_n = c / beta``.
* ``wfr_theory``  -- the theorem-bound-optimal constant ``h_n = mu_min`` with
  ``mu_min = min(lambda_min(C_0), 1/beta)``, the maximizer of the proven discrete
  W contribution ``alpha h / (1 + h/mu_min)^2``.
* ``wfr_adaptive`` -- the practical curvature-adaptive schedule
  ``h_n = h_max / (1 + (s_n / s0)^2)`` with
  ``s_n = lambda_min(C_n^{1/2}(-H_n)C_n^{1/2})``, ``h_max = 0.9/beta``, ``s0=0.5``.

``s_n`` measures how well the covariance is calibrated to the local curvature:
``s_n << 1`` means the state is underdispersed relative to the curvature (the
Fisher--Rao step is throttled and large transport helps warmup), while
``s_n ~ 1`` means the covariance is locally calibrated (the Fisher--Rao step is
reliable and transport decays toward zero).
"""
from __future__ import annotations

import numpy as np

from src.common.spd import symmetrize, symmetric_sqrt, eigh_spd

SCHEDULE_NAMES = ["fr_only", "w_only", "wfr_fixed", "wfr_theory", "wfr_adaptive"]


def theory_mu_min(C0, beta):
    """Theorem-bound optimal constant ``mu_min = min(lambda_min(C0), 1/beta)``.

    This maximizes the proven discrete W contribution
    ``f(h) = alpha h / (1 + h/mu_min)^2`` whose maximizer is ``h = mu_min``
    (``f'(h) = alpha (1 - h/mu_min)/(1 + h/mu_min)^3``).
    """
    lam0_min = float(eigh_spd(C0)[0][0])
    return float(min(lam0_min, 1.0 / float(beta)))


def calibration_ratio(C, H):
    """``s = lambda_min(C^{1/2} (-H) C^{1/2})``, the curvature calibration ratio.

    ``-H = E[Hess V]`` is SPD for a strongly log-concave target, so the symmetric
    product ``C^{1/2}(-H)C^{1/2}`` is SPD and ``s > 0``. ``s`` is the smallest
    eigenvalue of the dimensionless curvature-in-covariance-metric operator;
    ``s ~ 1`` is local calibration, ``s << 1`` is underdispersion.
    """
    C = symmetrize(C)
    Hbar = symmetrize(-np.asarray(H, dtype=np.float64))
    C_sqrt = symmetric_sqrt(C)
    prod = symmetrize(C_sqrt @ Hbar @ C_sqrt)
    return float(eigh_spd(prod)[0][0])


# ---------------------------------------------------------------------------
# Schedule objects
# ---------------------------------------------------------------------------

class _ConstantSchedule:
    """Constant Wasserstein step ``h_n = h_const`` (covers fr_only/fixed/theory)."""

    def __init__(self, name, h_const):
        self.name = name
        self.h_const = float(h_const)

    def h(self, n, m, C, g, H):
        return self.h_const, {}


class _AdaptiveSchedule:
    """Curvature-adaptive ``h_n = h_max / (1 + (s_n/s0)^2)``."""

    name = "wfr_adaptive"

    def __init__(self, h_max, s0):
        self.h_max = float(h_max)
        self.s0 = float(s0)

    def h(self, n, m, C, g, H):
        s = calibration_ratio(C, H)
        h = self.h_max / (1.0 + (s / self.s0) ** 2)
        return float(h), {"s": float(s)}


def build_schedule(method, *, beta, C0, c=0.5, h_max_frac=0.9, s0=0.5):
    """Construct the schedule object for ``method``.

    Parameters
    ----------
    method : one of :data:`SCHEDULE_NAMES`.
    beta : float        global smoothness constant of the target.
    C0 : array          initial covariance (for the theory schedule's mu_min).
    c : float           fixed-schedule fraction (``h = c/beta``) for w_only / wfr_fixed.
    h_max_frac : float  adaptive ceiling fraction (``h_max = h_max_frac / beta``).
    s0 : float          adaptive calibration scale.

    Returns a schedule object with a ``h(n, m, C, g, H) -> (h, info)`` method and
    a ``name`` attribute. ``info`` may carry the per-step ``s`` for the adaptive
    schedule (logged along the trajectory).
    """
    beta = float(beta)
    if method == "fr_only":
        return _ConstantSchedule("fr_only", 0.0)
    if method == "w_only":
        return _ConstantSchedule("w_only", c / beta)
    if method == "wfr_fixed":
        return _ConstantSchedule("wfr_fixed", c / beta)
    if method == "wfr_theory":
        return _ConstantSchedule("wfr_theory", theory_mu_min(C0, beta))
    if method == "wfr_adaptive":
        return _AdaptiveSchedule(h_max=h_max_frac / beta, s0=s0)
    raise ValueError(f"unknown method '{method}' (known: {SCHEDULE_NAMES})")
