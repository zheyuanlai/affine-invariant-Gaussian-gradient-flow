"""Run-level driver for the WFR experiment: simulate one run and summarize it.

A *run* is a fixed ``(target, Lambda, epsilon, method, dt, schedule params)``. The
run is stepped until it has spent a fixed *expectation-batch budget* (the fair
cost axis: the single-step methods spend one batch per iteration, the full WFR
splitting two), keeping the full trajectory in memory. The same ``g_n, H_n``
evaluated at ``a_n`` drive both the schedule (which may read the local curvature)
and the step, so the batch count is exact.

Per-step records carry the objective, the objective gap ``E(a_n) - E_star``, the
errors to ``a_star``, the eigenvalue extremes of ``C_n``, the scheduled Wasserstein
step ``h_n`` (and the recovered ``lambda_n = h_n / Delta t`` and adaptive ``s_n``),
the cumulative iteration and expectation-batch counts, and a monotone-so-far flag.
The summary aggregates the full trajectory and records hitting times to the gap
thresholds in both iterations and expectation batches.
"""
from __future__ import annotations

import math
import time

import numpy as np

from src.common.spd import symmetrize, eigh_spd
from src.wfr_gradient_flow.methods import wfr_step, WFR_METHODS
from src.wfr_gradient_flow.metrics import (
    hitting_table, classify_trajectory, monotone_prefix, GAP_THRESHOLDS,
)


def _gap(target, m, C, F_star):
    """Objective gap ``E(a) - E_star`` (analytic when the target provides it)."""
    if hasattr(target, "energy_gap"):
        return float(target.energy_gap(m, C))
    return float(target.objective(m, C) - F_star)


def simulate_run(method, target, dt, max_batches, schedule,
                 F_star, m_star, C_star, *, hard_iter_cap=200_000,
                 max_saved_rows=600):
    """Simulate one WFR run and return ``(records, summary)``.

    Steps the named ``method`` from the target's initial condition with
    Fisher--Rao step ``dt`` and the Wasserstein step supplied by ``schedule``,
    until the cumulative expectation-batch count reaches ``max_batches`` (or a
    failure / the ``hard_iter_cap`` safety bound). ``records`` is decimated to at
    most ``max_saved_rows`` rows (always including the first and last step);
    classification and hitting times always use the full trajectory.
    """
    d = len(target.m0)
    m = np.asarray(target.m0, dtype=np.float64).copy()
    C = symmetrize(target.C0)

    # Full-trajectory arrays for classification / hitting times.
    gaps_full, iters_full, batches_full, mineig_full = [], [], [], []
    records = []
    status = "ok"
    failed_iter = -1

    def push_full(gap, it, ba, mineig):
        gaps_full.append(gap); iters_full.append(it)
        batches_full.append(ba); mineig_full.append(mineig)

    # Step 0 (initial state, zero batches spent).
    w0 = eigh_spd(C)[0]
    gap0 = _gap(target, m, C, F_star)
    push_full(gap0, 0, 0, float(w0[0]))
    rec0 = _record(target, method, schedule, dt, 0, 0, gap0,
                   target.objective(m, C), m, C, m_star, C_star,
                   h_n=0.0, lam_n=0.0, s_n=float("nan"),
                   min_eig=float(w0[0]), max_eig=float(w0[-1]),
                   spd=bool(w0[0] > 0.0), monotone=True, status="ok")
    records.append(rec0)

    batches = 0
    wall = 0.0
    n = 0
    while batches < max_batches and n < hard_iter_cap:
        n += 1
        t0 = time.perf_counter()
        try:
            g, H = target.g_H(m, C)                       # batch at a_n
            h_n, info = schedule.h(n - 1, m, C, g, H)
            m_next, C_next, diag = wfr_step(method, target, m, C, g, H, h_n, dt)
        except Exception as exc:  # noqa: BLE001 -- record, never crash the grid
            status = f"error:{type(exc).__name__}"
            failed_iter = n
            wall += time.perf_counter() - t0
            break
        wall += time.perf_counter() - t0

        batches += int(diag["n_batches"])
        if diag["finite_ok"]:
            obj = target.objective(m_next, C_next)
            gap = _gap(target, m_next, C_next, F_star)
        else:
            obj, gap = float("nan"), float("nan")

        push_full(gap, n, batches, diag["min_eig_C"])
        s_n = float(info.get("s", float("nan")))
        lam_n = h_n / dt if dt > 0 else float("nan")
        rec = _record(target, method, schedule, dt, n, batches, gap, obj,
                      m_next, C_next, m_star, C_star, h_n=h_n, lam_n=lam_n,
                      s_n=s_n, min_eig=diag["min_eig_C"], max_eig=diag["max_eig_C"],
                      spd=bool(diag["spd_ok"]), monotone=True, status="ok")
        records.append(rec)

        m, C = m_next, C_next
        if (not diag["finite_ok"]) or diag["min_eig_C"] <= 1e-12:
            status = "spd_loss" if status == "ok" else status
            failed_iter = n
            break

    # Monotone-so-far flag from the full gap trajectory, mapped onto saved rows.
    mono_flags = monotone_prefix(gaps_full)
    iter_to_pos = {it: k for k, it in enumerate(iters_full)}
    for rec in records:
        k = iter_to_pos.get(rec["iteration"])
        rec["monotone_so_far"] = int(bool(mono_flags[k])) if k is not None else 0
        rec["status"] = status if rec is records[-1] else rec.get("status", "ok")

    # Decimate saved rows (keep first + last) without touching full arrays.
    records = _decimate(records, max_saved_rows)

    summary = _summarize(target, method, schedule, dt, F_star, m_star, C_star,
                         gaps_full, iters_full, batches_full, mineig_full,
                         records, wall, status, failed_iter, max_batches)
    return records, summary


def _record(target, method, schedule, dt, it, batches, gap, obj, m, C,
            m_star, C_star, *, h_n, lam_n, s_n, min_eig, max_eig, spd,
            monotone, status):
    m = np.asarray(m, dtype=np.float64)
    mean_err = float(np.linalg.norm(m - np.asarray(m_star)))
    cov_err = float(np.linalg.norm(symmetrize(C) - symmetrize(C_star), "fro"))
    return {
        "target_name": target.name, "Lambda": target.Lambda,
        "epsilon": target.epsilon, "method": method,
        "schedule_name": schedule.name, "dt": float(dt),
        "h_n": float(h_n), "lambda_n": float(lam_n), "s_n": float(s_n),
        "iteration": int(it), "expectation_batches": int(batches),
        "objective": float(obj), "objective_gap": float(gap),
        "mean_error_norm": mean_err, "covariance_error_fro": cov_err,
        "min_eig_C": float(min_eig), "max_eig_C": float(max_eig),
        "spd_feasible": int(bool(spd)), "monotone_so_far": int(bool(monotone)),
        "status": status,
    }


def _decimate(records, max_rows):
    if len(records) <= max_rows:
        return records
    n = len(records)
    step = max(1, n // max_rows)
    keep = set(range(0, n, step))
    keep.add(0)
    keep.add(n - 1)
    return [r for i, r in enumerate(records) if i in keep]


def _summarize(target, method, schedule, dt, F_star, m_star, C_star,
               gaps, iters, batches, mineigs, records, wall, status,
               failed_iter, max_batches):
    gaps_arr = np.asarray(gaps, dtype=np.float64)
    finite = gaps_arr[np.isfinite(gaps_arr)]
    cls = classify_trajectory(gaps, mineigs)
    hits = hitting_table(gaps, iters, batches)
    # The Wasserstein step actually used (constant for all but adaptive).
    h_used = [r["h_n"] for r in records if r["iteration"] > 0]
    h_final = float(records[-1]["h_n"]) if records else float("nan")
    summary = {
        "target_name": target.name, "Lambda": target.Lambda,
        "epsilon": target.epsilon, "method": method,
        "schedule_name": schedule.name, "dt": float(dt),
        "beta": float(getattr(target, "beta", float("nan"))),
        "alpha": float(getattr(target, "alpha", float("nan"))),
        "h_first": float(h_used[0]) if h_used else 0.0,
        "h_final": h_final,
        "lambda_first": (float(h_used[0]) / dt) if (h_used and dt > 0) else 0.0,
        "n_iters": int(iters[-1]) if iters else 0,
        "n_batches": int(batches[-1]) if batches else 0,
        "max_batches": int(max_batches),
        "gap0": float(gaps_arr[0]) if gaps_arr.size else float("nan"),
        "gap_final": float(gaps_arr[-1]) if gaps_arr.size else float("nan"),
        "gap_min": float(np.min(finite)) if finite.size else float("nan"),
        "min_eig_C_min": float(np.min([e for e in mineigs if math.isfinite(e)]))
        if any(math.isfinite(e) for e in mineigs) else float("nan"),
        "spd_feasible": int(cls["spd_feasible"]),
        "monotone": int(cls["monotone"]),
        "status": status, "failed_at_iter": int(failed_iter),
        "wall_time_total": float(wall),
    }
    summary.update(hits)
    return summary
