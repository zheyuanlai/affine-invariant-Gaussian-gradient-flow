"""Trajectory metrics for the algorithm-level STL study.

A *trajectory* is one stochastic natural-gradient run (fixed target, scheme, STL
flag, batch size, seed). Per step we record the energy gap and the covariance
eigenvalue extremes; richer per-step metrics (squared mean error, relative
covariance Frobenius error, W2^2 to the optimum) are computed at the saved
(decimated) steps. The *tail / noise-floor* metrics summarize the last fraction
of the run -- the stochastic floor STL is meant to lower.
"""
from __future__ import annotations

import math

import numpy as np

from src.common.spd import symmetrize, symmetric_sqrt

# Energy-gap thresholds for hitting-time reporting.
GAP_THRESHOLDS = [1e-1, 1e-2, 1e-3]


def w2_squared_gaussian(m, C, m_star, C_star):
    """Squared 2-Wasserstein distance between ``N(m, C)`` and ``N(m_star, C_star)``.

    ``W2^2 = ||m - m_star||^2 + Tr(C + C_star - 2 (C_star^{1/2} C C_star^{1/2})^{1/2})``.
    """
    m = np.asarray(m, dtype=np.float64)
    m_star = np.asarray(m_star, dtype=np.float64)
    C = symmetrize(C)
    C_star = symmetrize(C_star)
    Cs_sqrt = symmetric_sqrt(C_star)
    cross = symmetric_sqrt(symmetrize(Cs_sqrt @ C @ Cs_sqrt))
    bures = float(np.trace(C) + np.trace(C_star) - 2.0 * np.trace(cross))
    return float(np.dot(m - m_star, m - m_star)) + max(bures, 0.0)


def relative_cov_fro_error(C, C_star):
    """``||C - C_star||_F / ||C_star||_F``."""
    C = symmetrize(C)
    C_star = symmetrize(C_star)
    denom = np.linalg.norm(C_star, "fro")
    return float(np.linalg.norm(C - C_star, "fro") / denom) if denom > 0 else float("nan")


def _tol_key(tol):
    """Stable column-name suffix for a tolerance, e.g. 1e-3 -> '1e_minus_3'."""
    exp = int(round(math.log10(tol)))
    return f"1e_minus_{-exp}" if exp < 0 else f"1e{exp}"


def hitting_times(gaps, thresholds=GAP_THRESHOLDS):
    """First iteration index at which the energy gap drops to each threshold.

    Returns ``{iter_to_<tol>: k or -1}`` (``-1`` if never reached within the run).
    """
    gaps = np.asarray(gaps, dtype=np.float64)
    out = {}
    for tol in thresholds:
        hit = -1
        for k, g in enumerate(gaps):
            if math.isfinite(g) and g <= tol:
                hit = int(k)
                break
        out[f"iter_to_{_tol_key(tol)}"] = hit
    return out


def tail_noise_floor(gaps, tail_frac=0.2):
    """Tail statistics of the energy gap over the last ``tail_frac`` of the run.

    Returns the tail mean / median / std of the energy gap, the final gap, and the
    number of tail steps used. Non-finite gaps are dropped from the tail stats.
    """
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
                "tail_std_gap": float("nan"), "final_gap": final,
                "tail_steps": int(k)}
    return {
        "tail_mean_gap": float(np.mean(finite)),
        "tail_median_gap": float(np.median(finite)),
        "tail_std_gap": float(np.std(finite)),
        "final_gap": final,
        "tail_steps": int(k),
    }
