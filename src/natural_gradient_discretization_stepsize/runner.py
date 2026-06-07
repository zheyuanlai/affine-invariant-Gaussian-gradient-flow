"""Run-level driver: simulate one discretization run and summarize it.

A *run* is a fixed ``(target, lambda, method, dt)``. :func:`simulate_run` steps
the scheme for ``N = ceil(T / dt)`` steps from the target's initial condition,
keeping the full trajectory in memory (objective, energy gap, eigenvalue extremes,
errors to ``a_star`` and to the ODE reference). Classification
(:mod:`...metrics`) is always computed from the *full* trajectory; the long-format
CSV may be decimated for size without affecting any classification.
"""
from __future__ import annotations

import math
import time

import numpy as np

from src.common.spd import symmetrize, eigh_spd
from src.natural_gradient_discretization_stepsize.methods import discretization_step
from src.natural_gradient_discretization_stepsize.metrics import classify_run


def _cov_fro(C):
    return float(np.linalg.norm(symmetrize(C), "fro"))


def _accuracy_error(m, C, m_ref, C_ref):
    """``||m - m_ref||_2 + ||C - C_ref||_F / ||C_ref||_F``."""
    num = float(np.linalg.norm(np.asarray(m) - np.asarray(m_ref)))
    cref = _cov_fro(C_ref)
    rel = float(np.linalg.norm(symmetrize(C) - symmetrize(C_ref), "fro"))
    rel = rel / cref if cref > 0 else rel
    return num + rel


def simulate_run(method, target, dt, T, F_star, m_star, C_star, ode_ref,
                 max_saved_rows=400):
    """Simulate one run and return ``(records, summary)``.

    ``records`` is a list of per-saved-step dicts (decimated to at most
    ``max_saved_rows`` rows, always including the first and last step).
    ``summary`` aggregates the *full* trajectory.
    """
    d = len(target.m0)
    N = int(math.ceil(T / dt))
    save_every = max(1, N // max_saved_rows)

    m = np.asarray(target.m0, dtype=np.float64).copy()
    C = symmetrize(target.C0)
    F0 = target.objective(m, C)

    # Full-trajectory arrays (for classification).
    F_full = np.empty(N + 1, dtype=np.float64)
    gap_full = np.empty(N + 1, dtype=np.float64)
    min_eig_full = np.empty(N + 1, dtype=np.float64)

    records = []
    ode_eval = ode_ref["eval"]

    def record_step(n, t, F_n, gap_n, mn, Cn, diag, wall):
        m_err = float(np.linalg.norm(mn - m_star))
        cov_err = float(np.linalg.norm(symmetrize(Cn) - symmetrize(C_star), "fro"))
        m_ref, C_ref = ode_eval(t)
        ode_m_err = float(np.linalg.norm(mn - m_ref))
        cref = _cov_fro(C_ref)
        ode_cov_err = float(np.linalg.norm(symmetrize(Cn) - symmetrize(C_ref), "fro"))
        ode_cov_err = ode_cov_err / cref if cref > 0 else ode_cov_err
        records.append({
            "n": n, "t": t,
            "energy": F_n, "energy_gap": gap_n,
            "mean_error_to_star": m_err, "cov_error_to_star": cov_err,
            "ode_mean_error": ode_m_err, "ode_cov_error": ode_cov_err,
            "min_eig_C": diag["min_eig_C"], "max_eig_C": diag["max_eig_C"],
            "is_spd": int(bool(diag["spd_ok"])),
            "energy_increase": int(0),  # filled below for n>=1
            "wall_time_cumulative": wall,
        })

    # Step 0.
    w0 = eigh_spd(C)[0]
    diag0 = {"min_eig_C": float(w0[0]), "max_eig_C": float(w0[-1]),
             "spd_ok": bool(w0[0] > 0.0), "finite_ok": True}
    F_full[0] = F0
    gap_full[0] = F0 - F_star
    min_eig_full[0] = float(w0[0])
    record_step(0, 0.0, F0, F_full[0], m, C, diag0, 0.0)

    failed_n = None
    wall = 0.0
    last_n_filled = 0
    for n in range(1, N + 1):
        t = n * dt
        t0 = time.perf_counter()
        try:
            m_next, C_next, diag = discretization_step(method, target, m, C, dt)
        except Exception:
            diag = {"min_eig_C": float("nan"), "max_eig_C": float("nan"),
                    "spd_ok": False, "finite_ok": False}
            m_next = np.full(d, np.nan)
            C_next = np.full((d, d), np.nan)
        wall += time.perf_counter() - t0

        if diag["finite_ok"]:
            F_n = target.objective(m_next, C_next)
        else:
            F_n = float("nan")
        gap_n = F_n - F_star

        F_full[n] = F_n
        gap_full[n] = gap_n
        min_eig_full[n] = diag["min_eig_C"]
        last_n_filled = n

        if n % save_every == 0 or n == N or not diag["finite_ok"] or diag["min_eig_C"] <= 0:
            record_step(n, t, F_n, gap_n, m_next, C_next, diag, wall)
            if len(records) >= 2 and math.isfinite(F_n):
                prev_F = records[-2]["energy"]
                records[-1]["energy_increase"] = int(
                    math.isfinite(prev_F) and (F_n - prev_F) > 0.0)

        m, C = m_next, C_next
        if (not diag["finite_ok"]) or diag["min_eig_C"] <= 1e-12:
            failed_n = n
            break

    # Truncate full arrays if the run failed early.
    n_last = last_n_filled
    F_full = F_full[:n_last + 1]
    gap_full = gap_full[:n_last + 1]
    min_eig_full = min_eig_full[:n_last + 1]

    # Terminal accuracy: discrete final state vs ODE reference at t_final.
    if failed_n is None:
        m_ref_T, C_ref_T = ode_eval(N * dt)
        acc_err = _accuracy_error(m, C, m_ref_T, C_ref_T)
    else:
        acc_err = float("inf")

    cls = classify_run(F_full, gap_full, min_eig_full,
                       F0=F0, F_star=F_star, terminal_accuracy_error=acc_err)

    summary = _summarize(target, method, dt, T, N, F0, F_star,
                         F_full, gap_full, min_eig_full, records, cls, wall, failed_n)
    return records, summary


def _first_time_to_gap(records, gap_full, dt, tol):
    """First ``(time, iter)`` at which the full-trajectory gap drops to ``<= tol``."""
    for n in range(gap_full.size):
        if math.isfinite(gap_full[n]) and gap_full[n] <= tol:
            return float(n * dt), int(n)
    return math.inf, -1


def _summarize(target, method, dt, T, N, F0, F_star,
               F_full, gap_full, min_eig_full, records, cls, wall, failed_n):
    finite_gaps = gap_full[np.isfinite(gap_full)]
    gap_min = float(np.min(finite_gaps)) if finite_gaps.size else float("nan")
    finite_eig = min_eig_full[np.isfinite(min_eig_full)]
    max_eig_vals = [r["max_eig_C"] for r in records if math.isfinite(r["max_eig_C"])]
    t_1e4, i_1e4 = _first_time_to_gap(records, gap_full, dt, 1e-4)
    t_1e6, i_1e6 = _first_time_to_gap(records, gap_full, dt, 1e-6)
    F_final = float(F_full[-1]) if F_full.size else float("nan")
    return {
        "target_name": target.name, "lambda": getattr(target, "lam", float("nan")),
        "method": method, "dt": dt, "T": T, "N": N,
        "F0": float(F0), "F_final": F_final,
        "gap_final": float(gap_full[-1]) if gap_full.size else float("nan"),
        "gap_min": gap_min,
        "spd_feasible": int(cls["spd_feasible"]), "stable": int(cls["stable"]),
        "monotone": int(cls["monotone"]), "accurate": int(cls["accurate"]),
        "num_energy_increases": cls["num_energy_increases"],
        "max_energy_increase": cls["max_energy_increase"],
        "max_gap_ratio": cls["max_gap_ratio"],
        "min_eig_C_min": float(np.min(finite_eig)) if finite_eig.size else float("nan"),
        "max_eig_C_max": float(max(max_eig_vals)) if max_eig_vals else float("nan"),
        "terminal_accuracy_error": cls["terminal_accuracy_error"],
        "failed_at_step": -1 if failed_n is None else int(failed_n),
        "wall_time_total": float(wall),
        "time_to_gap_1e_minus_4": t_1e4, "time_to_gap_1e_minus_6": t_1e6,
        "iter_to_gap_1e_minus_4": i_1e4, "iter_to_gap_1e_minus_6": i_1e6,
    }


# ---------------------------------------------------------------------------
# Scalar diagnostic (target D)
# ---------------------------------------------------------------------------

def simulate_scalar_diagnostic(method, C0, dt, T, m0=0.0):
    """Scalar ``N(0,1)`` diagnostic trajectory for one ``(method, C0, dt)``.

    Steps the scalar covariance recursion (closed form, matching the dense
    implementation) and the shared mean recursion ``m_{n+1} = (1 - dt C_n) m_n``.
    Returns a list of per-step dicts ``{n, t, C, m}``.
    """
    N = int(math.ceil(T / dt))
    C = float(C0)
    m = float(m0)
    out = [{"n": 0, "t": 0.0, "C": C, "m": m}]
    for n in range(1, N + 1):
        m = (1.0 - dt * C) * m
        if method == "kl":
            C = (1.0 + dt) * C / (1.0 + dt * C)
        elif method == "riemannian":
            C = C * math.exp(dt * (1.0 - C))
        else:
            raise ValueError(f"unknown method '{method}'")
        if not (math.isfinite(C) and math.isfinite(m)):
            out.append({"n": n, "t": n * dt, "C": float("nan"), "m": float("nan")})
            break
        out.append({"n": n, "t": n * dt, "C": C, "m": m})
    return out


# ---------------------------------------------------------------------------
# Stepsize summary (one row per target/lambda/method)
# ---------------------------------------------------------------------------

def stepsize_summary_rows(summaries, theory_by_key):
    """Aggregate per-run summaries into one row per ``(target, lambda, method)``.

    ``summaries`` is the list of per-run summary dicts; ``theory_by_key`` maps
    ``(target_name, lambda, method)`` -> theory-bounds dict (see
    :func:`...metrics.theory_stepsize_bounds`) augmented with the per-method
    ``dt_theory_for_method``.
    """
    from collections import defaultdict
    from src.natural_gradient_discretization_stepsize.metrics import max_feasible_dt

    groups = defaultdict(list)
    for s in summaries:
        groups[(s["target_name"], s["lambda"], s["method"])].append(s)

    rows = []
    for key, runs in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])):
        target_name, lam, method = key
        runs = sorted(runs, key=lambda r: r["dt"])
        dts = [r["dt"] for r in runs]
        th = theory_by_key.get(key, {})
        dt_spd = max_feasible_dt(dts, [r["spd_feasible"] for r in runs])
        dt_stable = max_feasible_dt(dts, [r["stable"] for r in runs])
        dt_mono = max_feasible_dt(dts, [r["monotone"] for r in runs])
        dt_acc = max_feasible_dt(dts, [r["accurate"] for r in runs])
        dt_theory = th.get("dt_theory_for_method")

        def ratio(x):
            if dt_theory and dt_theory > 0 and math.isfinite(x):
                return float(x / dt_theory)
            return float("nan")

        rows.append({
            "target_name": target_name, "lambda": lam, "method": method,
            "theory_bound_available": int(bool(th.get("theory_bound_available"))),
            "dt_theory_riem": th.get("dt_theory_riem"),
            "dt_theory_kl": th.get("dt_theory_kl"),
            "dt_theory_for_method": dt_theory,
            "dt_max_spd": dt_spd, "dt_max_stable": dt_stable,
            "dt_max_monotone": dt_mono, "dt_max_accurate": dt_acc,
            "stable_over_theory_ratio": ratio(dt_stable),
            "monotone_over_theory_ratio": ratio(dt_mono),
            "accurate_over_theory_ratio": ratio(dt_acc),
        })
    return rows
