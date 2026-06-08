"""Per-trajectory metrics for the WFR runs: hitting times and feasibility.

A *run* produces, at every iteration ``n = 0, 1, ...``, the objective gap
``E(a_n) - E_star``, the cumulative number of expectation batches spent, and the
eigenvalue extremes of ``C_n``. Because the WFR splitting costs two expectation
batches per iteration while the single-step methods cost one, the *expectation
batch count* is the fair cost axis; we record both iteration count and batch
count and report hitting times in both units.
"""
from __future__ import annotations

import math

import numpy as np

# Objective-gap thresholds for hitting-time reporting.
GAP_THRESHOLDS = [1e-1, 1e-3, 1e-6]
SPD_TOL = 1e-12


def _hitting(gaps, iters, batches, tol):
    """First ``(iteration, batches)`` at which ``gap <= tol``; inf / -1 if never."""
    for k in range(len(gaps)):
        g = gaps[k]
        if math.isfinite(g) and g <= tol:
            return int(iters[k]), int(batches[k])
    return -1, -1


def hitting_table(gaps, iters, batches, thresholds=GAP_THRESHOLDS):
    """Map each threshold to ``(iter_to, batches_to)`` (``-1`` if never reached)."""
    out = {}
    for tol in thresholds:
        it, ba = _hitting(gaps, iters, batches, tol)
        key = _tol_key(tol)
        out[f"iter_to_{key}"] = it
        out[f"batches_to_{key}"] = ba
    return out


def _tol_key(tol):
    """Stable column-name suffix for a tolerance, e.g. 1e-3 -> '1e_minus_3'."""
    exp = int(round(math.log10(tol)))
    return f"1e_minus_{-exp}" if exp < 0 else f"1e{exp}"


def classify_trajectory(gaps, min_eigs, *, spd_tol=SPD_TOL, monotone_tol=1e-9):
    """SPD-feasibility and monotonicity over a full trajectory.

    Returns ``spd_feasible`` (no NaN/Inf, ``min_eig(C_n) > spd_tol`` throughout)
    and ``monotone`` (the objective gap is non-increasing to within
    ``monotone_tol``). Both are computed from the full trajectory.
    """
    gaps = np.asarray(gaps, dtype=np.float64)
    min_eigs = np.asarray(min_eigs, dtype=np.float64)
    finite_ok = bool(np.all(np.isfinite(gaps)) and np.all(np.isfinite(min_eigs)))
    spd_feasible = bool(finite_ok and np.all(min_eigs > spd_tol))
    if finite_ok and gaps.size >= 2:
        monotone = bool(np.all(np.diff(gaps) <= monotone_tol))
    else:
        monotone = False
    return {"spd_feasible": spd_feasible, "monotone": monotone}


def monotone_prefix(gaps, monotone_tol=1e-9):
    """Boolean per-step flag: gap non-increasing on every step up to ``n``."""
    gaps = np.asarray(gaps, dtype=np.float64)
    flags = np.ones(gaps.size, dtype=bool)
    ok = True
    for n in range(1, gaps.size):
        ok = ok and math.isfinite(gaps[n]) and (gaps[n] - gaps[n - 1] <= monotone_tol)
        flags[n] = ok
    return flags
