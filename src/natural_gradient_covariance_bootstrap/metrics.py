"""Trajectory metrics: first-crossing (N_cov), hitting times, tail floor.

A deterministic trajectory records at every iteration ``n = 0, 1, ...`` the
covariance eigenvalue extremes and the energy gap ``Delta_n``. ``N_cov`` is the
first iteration at which ``lambda_min(C_n)`` reaches the curvature scale
``1/(2 beta)``; hitting times are the first iterations at which the energy gap
drops below the reporting thresholds. All classifications use the *full*
trajectory (before any decimation of the saved rows).
"""
from __future__ import annotations

import math

import numpy as np

# Energy-gap thresholds for hitting-time reporting.
GAP_THRESHOLDS = [1e-1, 1e-3, 1e-6]


def tol_key(tol):
    """Stable column-name suffix for a tolerance, e.g. 1e-3 -> '1e_minus_3'."""
    exp = int(round(math.log10(tol)))
    return f"1e_minus_{-exp}" if exp < 0 else f"1e{exp}"


def first_crossing_ge(values, threshold):
    """First index ``n`` with ``values[n] >= threshold`` (``-1`` if never)."""
    values = np.asarray(values, dtype=np.float64)
    hits = np.where(np.isfinite(values) & (values >= threshold))[0]
    return int(hits[0]) if hits.size else -1


def first_crossing_le(values, threshold):
    """First index ``n`` with ``values[n] <= threshold`` (``-1`` if never)."""
    values = np.asarray(values, dtype=np.float64)
    hits = np.where(np.isfinite(values) & (values <= threshold))[0]
    return int(hits[0]) if hits.size else -1


def n_cov(min_eigs, beta):
    """``N_cov`` = first iteration with ``lambda_min(C_n) >= 1/(2 beta)``."""
    return first_crossing_ge(min_eigs, 1.0 / (2.0 * float(beta)))


def hitting_times(gaps, thresholds=GAP_THRESHOLDS):
    """Map each gap threshold to the first iteration below it (``-1`` if never)."""
    out = {}
    for tol in thresholds:
        out[f"iter_to_{tol_key(tol)}"] = first_crossing_le(gaps, tol)
    return out


def tail_floor(gaps, tail_frac=0.25):
    """Tail statistics of the energy gap over the last ``tail_frac`` of the run."""
    gaps = np.asarray(gaps, dtype=np.float64)
    n = gaps.size
    if n == 0:
        return {"tail_mean_gap": float("nan"), "tail_median_gap": float("nan"),
                "tail_std_gap": float("nan"), "final_gap": float("nan"),
                "tail_steps": 0}
    k = max(1, int(round(tail_frac * n)))
    tail = gaps[-k:]
    finite = tail[np.isfinite(tail)]
    final = float(gaps[-1]) if math.isfinite(gaps[-1]) else float("nan")
    if finite.size == 0:
        return {"tail_mean_gap": float("nan"), "tail_median_gap": float("nan"),
                "tail_std_gap": float("nan"), "final_gap": final, "tail_steps": int(k)}
    return {
        "tail_mean_gap": float(np.mean(finite)),
        "tail_median_gap": float(np.median(finite)),
        "tail_std_gap": float(np.std(finite)),
        "final_gap": final, "tail_steps": int(k),
    }


def decimate_indices(n_rows, max_rows, keep=()):
    """Indices to persist: ``max_rows`` evenly spaced, plus first/last and ``keep``."""
    if n_rows <= max_rows:
        return list(range(n_rows))
    idx = set(np.linspace(0, n_rows - 1, num=max_rows).round().astype(int).tolist())
    idx.add(0)
    idx.add(n_rows - 1)
    for k in keep:
        if 0 <= k < n_rows:
            idx.add(int(k))
    return sorted(idx)


HIT_KEYS = [f"iter_to_{tol_key(t)}" for t in GAP_THRESHOLDS]
