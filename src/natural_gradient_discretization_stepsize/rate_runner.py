"""Run-level driver for the theoretical-rate benchmark (Riemannian-scale grid).

A *rate run* is a fixed ``(target, lambda, method, c)`` where the stepsize is
``dt = c * dt_ref`` on the common Riemannian-scale reference
``dt_ref = 1/(beta * lambda_max)`` (see :mod:`...rate_metrics`). Each run steps
the named scheme for ``N = ceil(T / dt)`` steps from the target's initial
condition, keeping the full per-step energy gap in memory for the observed-rate
and tolerance diagnostics; the long-format CSV may be decimated for size, but
every classification / rate / tolerance value uses the full trajectory.

Only globally-smooth targets (Gaussian posterior, smooth log-concave) are run
here -- the non-smooth quartic has no global beta and hence no clean theorem
rate to compare against.
"""
from __future__ import annotations

import math

import numpy as np

from src.common.spd import symmetrize
from src.natural_gradient_discretization_stepsize.methods import discretization_step
from src.natural_gradient_discretization_stepsize import rate_metrics as rm


def _energy_gap(target, m, C, F_star):
    """Energy gap ``E(a) - E(a_star)`` (analytic for Gaussian, objective else)."""
    if hasattr(target, "energy_gap"):
        return float(target.energy_gap(m, C))
    return float(target.objective(m, C) - F_star)


def simulate_rate_run(method, target, c, dt_ref, alpha, beta,
                      lambda_min, lambda_max, T, F_star,
                      eps_list, gap_floor=rm.GAP_FLOOR, max_saved_rows=400):
    """Simulate one rate run; return ``(records, summary, tol_rows)``.

    ``records`` is the decimated per-step long-format list (first and last step
    always included). ``summary`` is the one-row-per-run dict. ``tol_rows`` is
    one row per tolerance in ``eps_list``.
    """
    dt = c * dt_ref
    N = int(math.ceil(T / dt))
    save_every = max(1, N // max_saved_rows)

    # Method-specific theoretical contraction factor at this dt (state-free).
    q_riem = rm.q_riem_theory(dt, alpha, beta, lambda_min, lambda_max)
    q_kl = rm.q_kl_formula(dt, alpha, beta, lambda_min, lambda_max)
    r_riem = rm.per_unit_rate(q_riem, dt)
    r_kl = rm.per_unit_rate(q_kl, dt)
    q_theory = q_riem if method == "riemannian" else q_kl
    r_theory = r_riem if method == "riemannian" else r_kl

    m = np.asarray(target.m0, dtype=np.float64).copy()
    C = symmetrize(target.C0)
    gap0 = _energy_gap(target, m, C, F_star)

    gap_full = np.empty(N + 1, dtype=np.float64)
    time_full = np.empty(N + 1, dtype=np.float64)
    gap_full[0] = gap0
    time_full[0] = 0.0

    records = []

    def record_step(n, t, gap_raw):
        gap_floored = max(gap_raw, gap_floor) if math.isfinite(gap_raw) else float("nan")
        records.append({
            "n": n, "time": t,
            "gap": gap_floored, "gap_raw": gap_raw,
            "is_floor_limited": int(math.isfinite(gap_raw) and gap_raw <= gap_floor),
            "q_riem_theory": q_riem, "q_kl_formula": q_kl,
            "r_riem_theory": r_riem, "r_kl_formula": r_kl,
        })

    record_step(0, 0.0, gap0)

    failed_n = None
    last_n = 0
    for n in range(1, N + 1):
        t = n * dt
        try:
            m_next, C_next, diag = discretization_step(method, target, m, C, dt)
        except Exception:
            diag = {"finite_ok": False, "min_eig_C": float("nan")}
            m_next = np.full(len(target.m0), np.nan)
            C_next = np.full((len(target.m0),) * 2, np.nan)
        if diag["finite_ok"]:
            gap_n = _energy_gap(target, m_next, C_next, F_star)
        else:
            gap_n = float("nan")
        gap_full[n] = gap_n
        time_full[n] = t
        last_n = n
        if n % save_every == 0 or n == N or not diag["finite_ok"]:
            record_step(n, t, gap_n)
        m, C = m_next, C_next
        if (not diag["finite_ok"]) or diag.get("min_eig_C", 0.0) <= 1e-12:
            failed_n = n
            break

    gap_full = gap_full[:last_n + 1]
    time_full = time_full[:last_n + 1]
    N_actual = last_n

    summary = _summarize_rate_run(
        target, method, c, dt, dt_ref, alpha, beta, lambda_min, lambda_max,
        N_actual, T, gap0, gap_full, time_full, q_theory, r_theory,
        gap_floor, failed_n)
    tol_rows = _tolerance_rows(
        target, method, c, dt, gap_full, gap0, q_theory, eps_list, failed_n)
    return records, summary, tol_rows


# ---------------------------------------------------------------------------
# Per-run summary and tolerance rows
# ---------------------------------------------------------------------------

def _summarize_rate_run(target, method, c, dt, dt_ref, alpha, beta,
                        lambda_min, lambda_max, N, T, gap0, gap_full, time_full,
                        q_theory, r_theory, gap_floor, failed_n):
    final_gap_raw = float(gap_full[-1]) if gap_full.size else float("nan")
    term = rm.terminal_rate(gap0, final_gap_raw, N, dt, gap_floor=gap_floor)
    fit = rm.fitted_rate(time_full, gap_full, gap0, dt, gap_floor=gap_floor)

    # Method-specific theory terminal bound and terminal slack.
    if math.isfinite(q_theory) and 0.0 < q_theory <= 1.0 and N >= 1:
        theory_terminal_bound = (q_theory ** N) * gap0
    else:
        theory_terminal_bound = float("nan")
    final_for_logs = term["final_gap_for_logs"]
    if math.isfinite(theory_terminal_bound) and final_for_logs > 0.0:
        terminal_slack = theory_terminal_bound / final_for_logs
        log10_slack = math.log10(terminal_slack) if terminal_slack > 0.0 else float("nan")
    else:
        terminal_slack = float("nan")
        log10_slack = float("nan")

    status = "ok" if failed_n is None else "failed"
    return {
        "target": target.name, "lambda": getattr(target, "lam", float("nan")),
        "method": method, "alpha": alpha, "beta": beta,
        "lambda_min": lambda_min, "lambda_max": lambda_max,
        "dt_ref": dt_ref, "c": c, "dt": dt, "N": N, "T_actual": N * dt,
        "initial_gap": gap0,
        "final_gap_raw": final_gap_raw,
        "final_gap_for_logs": final_for_logs,
        "q_theory": float(q_theory), "r_theory": float(r_theory),
        "q_hat_terminal": term["q_hat_terminal"],
        "r_hat_terminal": term["r_hat_terminal"],
        "q_hat_fit": fit["q_hat_fit"], "r_hat_fit": fit["r_hat_fit"],
        "fit_num_points": fit["fit_num_points"],
        "fit_window_status": fit["fit_window_status"],
        "theory_terminal_bound": float(theory_terminal_bound),
        "terminal_slack": float(terminal_slack),
        "log10_terminal_slack": float(log10_slack),
        "floor_limited_final": term["floor_limited_final"],
        "r_continuous": rm.r_continuous(alpha, lambda_min),
        "status": status,
    }


def _tolerance_rows(target, method, c, dt, gap_full, gap0, q_theory,
                    eps_list, failed_n):
    rows = []
    for eps in eps_list:
        n_obs = rm.n_obs_to_tol(gap_full, eps)
        n_theory = rm.n_theory_to_tol(eps, gap0, q_theory)
        if n_obs is None:
            obs_status, t_obs, n_obs_v = "not_reached", float("nan"), float("nan")
        else:
            obs_status, t_obs, n_obs_v = "reached", n_obs * dt, n_obs
        if n_theory is None:
            theory_status, t_theory, n_theory_v = "undefined", float("nan"), float("nan")
        else:
            theory_status, t_theory, n_theory_v = "ok", n_theory * dt, n_theory
        if (n_obs is not None) and (n_theory is not None) and n_obs > 0:
            ratio = float(n_theory) / float(n_obs)
        else:
            ratio = float("nan")
        rows.append({
            "target": target.name, "lambda": getattr(target, "lam", float("nan")),
            "method": method, "c": c, "dt": dt, "eps": eps,
            "N_obs": n_obs_v, "T_obs": t_obs, "obs_status": obs_status,
            "N_theory": n_theory_v, "T_theory": t_theory, "theory_status": theory_status,
            "N_theory_over_N_obs": ratio,
        })
    return rows

