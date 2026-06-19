"""Deterministic runners for the nonconvex-instability experiment group."""
from __future__ import annotations

import math
import os
import platform
import sys
import time

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from src.common.io_utils import ensure_dir, save_dataframe, save_json
from src.natural_gradient_nonconvex_instability.methods import (
    clip_smoothness_constant,
    clipped_kl_next,
    gaussian_kl_divergence,
    kl_next,
    riemannian_next_log_c,
    theorem_safe_dt,
    wasserstein_fb_next,
)
from src.natural_gradient_nonconvex_instability.targets import (
    NonconvexLogCoshTarget,
    bw_gradient_sq,
    fr_gradient_sq,
    mode_location_unit,
)


RESULTS_LONG_COLS = [
    "experiment", "method", "R", "epsilon", "dt", "eta", "n", "N",
    "c", "log_c", "A", "energy", "fr_grad_sq", "bw_grad_sq",
    "running_min_bw_grad_sq", "theorem_envelope", "status", "stop_reason",
    "denom", "ctilde",
]

SUMMARY_COLS = [
    "experiment", "method", "R", "epsilon", "dt", "eta", "num_steps",
    "status", "stop_reason", "final_c", "final_log_c", "max_fr_grad_sq",
    "max_fr_grad_norm", "max_bw_grad_sq", "running_min_bw_grad_sq",
    "theorem_envelope_final", "bw_bound_satisfied_all_N", "F0", "F_star_num",
    "Delta", "wall_time_seconds",
]

KL_POLE_COLS = [
    "R", "dt", "epsilon", "c0", "denom0", "c1", "status",
    "fr_grad_sq_after_first", "bw_grad_sq_after_first",
    "expected_c1_2_over_epsilon", "expected_fr_grad_sq_4_over_epsilon_sq",
    "c1_ratio_to_reference", "fr_grad_sq_ratio_to_reference",
    "in_pole_fit_window", "fit_slope_log_fr_vs_log_epsilon",
    "fit_slope_pole_window_log_fr_vs_log_epsilon",
]

BW_BOUND_COLS = [
    "R", "eta", "N", "running_min_bw_grad_sq", "theorem_envelope",
    "bound_satisfied", "F0", "F_star_num", "Delta",
]

CLIPPED_KL_COLS = [
    "R", "lambda_minus", "lambda_plus", "beta", "L_clip", "dt_safety", "dt",
    "dt_rule", "step", "c", "c_tilde", "c_next", "A", "F",
    "D_kl_current_to_next", "running_min_D", "prefix_bound", "bound_margin",
    "clipped_lower_active", "clipped_upper_active", "denominator", "status",
]

CLIPPED_SUMMARY_COLS = [
    "R", "lambda_minus", "lambda_plus", "dt", "dt_rule", "dt_times_L_clip",
    "num_steps", "final_c", "min_c", "max_c", "final_running_min_D",
    "final_bound", "min_bound_margin", "max_violation", "first_upper_clip_step",
    "first_lower_clip_step", "energy_monotone", "denominator_positive",
    "theorem_check_pass",
]

# Absolute tolerance on the theorem-check residual (running_min_D - prefix_bound).
CLIPPED_THEOREM_TOL = 1.0e-9


def shallow_apply_smoke(cfg):
    """Shallow-merge ``smoke:`` over a config dict."""
    cfg = dict(cfg)
    for key, value in cfg.get("smoke", {}).items():
        cfg[key] = value
    return cfg


def _finite_max(values):
    arr = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    return float(np.max(arr)) if arr.size else float("nan")


def _terminal_summary(rows, *, experiment, method, R, epsilon=np.nan,
                      dt=np.nan, eta=np.nan, start_time=None,
                      F0=np.nan, F_star=np.nan, Delta=np.nan,
                      bw_bound_ok=np.nan):
    last = rows[-1] if rows else {}
    max_fr = _finite_max([r.get("fr_grad_sq", np.nan) for r in rows])
    max_bw = _finite_max([r.get("bw_grad_sq", np.nan) for r in rows])
    return {
        "experiment": experiment,
        "method": method,
        "R": R,
        "epsilon": epsilon,
        "dt": dt,
        "eta": eta,
        "num_steps": int(last.get("n", 0)),
        "status": last.get("status", "empty"),
        "stop_reason": last.get("stop_reason", ""),
        "final_c": last.get("c", np.nan),
        "final_log_c": last.get("log_c", np.nan),
        "max_fr_grad_sq": max_fr,
        "max_fr_grad_norm": math.sqrt(max_fr) if np.isfinite(max_fr) else np.nan,
        "max_bw_grad_sq": max_bw,
        "running_min_bw_grad_sq": last.get("running_min_bw_grad_sq", np.nan),
        "theorem_envelope_final": last.get("theorem_envelope", np.nan),
        "bw_bound_satisfied_all_N": bw_bound_ok,
        "F0": F0,
        "F_star_num": F_star,
        "Delta": Delta,
        "wall_time_seconds": time.perf_counter() - start_time if start_time else np.nan,
    }


def _state_row(experiment, method, target, c, log_c, n, *,
               epsilon=np.nan, dt=np.nan, eta=np.nan, N=np.nan,
               status="ok", stop_reason="", running_min=np.nan,
               envelope=np.nan, denom=np.nan, ctilde=np.nan):
    A = target.A0(c)
    energy = target.energy(0.0, c)
    fr = fr_gradient_sq(c, A)
    bw = bw_gradient_sq(c, A)
    return {
        "experiment": experiment,
        "method": method,
        "R": target.R,
        "epsilon": epsilon,
        "dt": dt,
        "eta": eta,
        "n": int(n),
        "N": N,
        "c": float(c),
        "log_c": float(log_c),
        "A": float(A),
        "energy": float(energy),
        "fr_grad_sq": float(fr),
        "bw_grad_sq": float(bw),
        "running_min_bw_grad_sq": running_min,
        "theorem_envelope": envelope,
        "status": status,
        "stop_reason": stop_reason,
        "denom": denom,
        "ctilde": ctilde,
    }


def estimate_F_star(target, *, m_starts=None, logc_starts=None,
                    maxiter=300, gtol=1e-10):
    """Robust deterministic multistart minimization over ``(m, log c)``."""
    R = target.R
    ystar = mode_location_unit()
    if m_starts is None:
        m_starts = [0.0, -ystar * R, ystar * R]
    if logc_starts is None:
        logc_starts = [-2.0, 0.0, 2.0]
    upper_logc = math.log(max(10.0, 8.0 * R * R))
    bounds = [(-4.0 * R, 4.0 * R), (-8.0, upper_logc)]

    best = None
    starts = []
    for m0 in m_starts:
        for ell0 in logc_starts:
            x0 = np.array([float(m0), float(ell0)], dtype=np.float64)
            starts.append(x0.tolist())

            def fun(z):
                f, g = target.objective_and_grad_m_logc(z)
                return f, g

            res = minimize(
                fun,
                x0,
                method="L-BFGS-B",
                jac=True,
                bounds=bounds,
                options={"maxiter": int(maxiter), "gtol": float(gtol), "ftol": 1e-14},
            )
            if best is None or float(res.fun) < float(best.fun):
                best = res

    if best is None:
        raise RuntimeError("no optimizer runs were attempted")
    m_star = float(best.x[0])
    logc_star = float(best.x[1])
    c_star = math.exp(logc_star)
    _, grad = target.objective_and_grad_m_logc(best.x)
    return {
        "F_star_num": float(best.fun),
        "m_star": m_star,
        "c_star": c_star,
        "log_c_star": logc_star,
        "grad_norm_m_logc": float(np.linalg.norm(grad)),
        "success": bool(best.success),
        "status": int(best.status),
        "message": str(best.message),
        "nit": int(best.nit),
        "nfev": int(best.nfev),
        "starts": starts,
        "bounds": bounds,
    }


def simulate_riemannian_cascade(target, *, dt=1.0, c0=1.0, max_iter=20,
                                stop_fraction=0.25, max_log_c=700.0):
    """Run the Riemannian Fisher-Rao covariance cascade at ``m=0``."""
    start = time.perf_counter()
    rows = []
    c = float(c0)
    log_c = math.log(c)
    threshold = stop_fraction * target.R * target.R
    for n in range(int(max_iter) + 1):
        row = _state_row(
            "riemannian_cascade", "riemannian", target, c, log_c, n, dt=dt)
        if c > threshold:
            row["status"] = "stopped"
            row["stop_reason"] = "c_threshold"
            rows.append(row)
            break
        rows.append(row)
        A = row["A"]
        if n == int(max_iter):
            rows[-1]["status"] = "stopped"
            rows[-1]["stop_reason"] = "max_iter"
            break
        next_log = riemannian_next_log_c(log_c, c, A, dt)
        if (not math.isfinite(next_log)) or next_log > float(max_log_c):
            rows[-1]["status"] = "stopped"
            rows[-1]["stop_reason"] = "log_c_overflow_guard"
            break
        log_c = next_log
        c = math.exp(log_c)
        if not (math.isfinite(c) and c > 0.0):
            rows[-1]["status"] = "stopped"
            rows[-1]["stop_reason"] = "overflow"
            break

    summary = _terminal_summary(
        rows, experiment="riemannian_cascade", method="riemannian",
        R=target.R, dt=dt, start_time=start)
    summary["c_threshold"] = threshold
    return rows, summary


def simulate_kl_pole(target, *, epsilon, dt=1.0, max_iter=2):
    """Run one/two KL steps from ``c0 = 1 - epsilon``."""
    start = time.perf_counter()
    eps = float(epsilon)
    c = 1.0 - eps
    log_c = math.log(c)
    rows = []
    for n in range(int(max_iter) + 1):
        row = _state_row("kl_pole", "kl", target, c, log_c, n, epsilon=eps, dt=dt)
        rows.append(row)
        if n == int(max_iter):
            rows[-1]["status"] = "stopped"
            rows[-1]["stop_reason"] = "max_iter"
            break
        step = kl_next(c, row["A"], dt)
        rows[-1]["denom"] = step.denom
        if step.status != "ok":
            rows[-1]["status"] = "stopped"
            rows[-1]["stop_reason"] = step.status
            break
        c = step.c_next
        if not (math.isfinite(c) and c > 0.0):
            rows[-1]["status"] = "stopped"
            rows[-1]["stop_reason"] = "overflow"
            break
        log_c = math.log(c)
    summary = _terminal_summary(
        rows, experiment="kl_pole", method="kl", R=target.R,
        epsilon=eps, dt=dt, start_time=start)
    return rows, summary


def simulate_wasserstein_fb(target, star, *, eta=0.9, c0=1.0, max_iter=20000):
    """Run Bures-Wasserstein forward-backward at ``m=0``."""
    start = time.perf_counter()
    F0 = target.energy(0.0, c0)
    F_star = float(star["F_star_num"])
    Delta = F0 - F_star
    rows = []
    bound_rows = []
    running_min = float("inf")
    c = float(c0)
    for n in range(int(max_iter)):
        log_c = math.log(c)
        A = target.A0(c)
        fr = fr_gradient_sq(c, A)
        bw = fr / c
        running_min = min(running_min, bw)
        N = n + 1
        envelope = 150.0 * Delta / (float(eta) * N)
        row = {
            "experiment": "wasserstein_fb",
            "method": "wasserstein_fb",
            "R": target.R,
            "epsilon": np.nan,
            "dt": np.nan,
            "eta": eta,
            "n": n,
            "N": N,
            "c": c,
            "log_c": log_c,
            "A": A,
            "energy": target.energy(0.0, c),
            "fr_grad_sq": fr,
            "bw_grad_sq": bw,
            "running_min_bw_grad_sq": running_min,
            "theorem_envelope": envelope,
            "status": "ok",
            "stop_reason": "",
            "denom": np.nan,
            "ctilde": np.nan,
        }
        step = wasserstein_fb_next(c, A, eta)
        row["ctilde"] = step.ctilde
        rows.append(row)
        bound_rows.append({
            "R": target.R,
            "eta": eta,
            "N": N,
            "running_min_bw_grad_sq": running_min,
            "theorem_envelope": envelope,
            "bound_satisfied": int(running_min <= envelope * (1.0 + 1e-12)),
            "F0": F0,
            "F_star_num": F_star,
            "Delta": Delta,
        })
        if step.status != "ok":
            rows[-1]["status"] = "stopped"
            rows[-1]["stop_reason"] = step.status
            break
        c = step.c_next
    else:
        rows[-1]["status"] = "stopped"
        rows[-1]["stop_reason"] = "max_iter"

    bw_ok = int(all(r["bound_satisfied"] for r in bound_rows))
    summary = _terminal_summary(
        rows, experiment="wasserstein_fb", method="wasserstein_fb",
        R=target.R, eta=eta, start_time=start, F0=F0, F_star=F_star,
        Delta=Delta, bw_bound_ok=bw_ok)
    return rows, summary, bound_rows


def simulate_clipped_kl(target, *, lambda_minus=0.5, lambda_plus=2.0, beta=1.0,
                        dt_safety=0.9, c0=1.0, num_steps=300, dt=None,
                        dt_rule="theorem_safe"):
    """Projected KL covariance scheme with eigenvalue clipping (Theorem 2.18).

    The covariance is evolved by the unprojected KL step and then clipped back
    into ``[lambda_minus, lambda_plus]``. For each prefix horizon ``N`` we record
    the running-minimum Bregman displacement ``D_min(N) = min_{0<=n<N} D_n`` with
    ``D_n = KL(rho_{a_n} || rho_{a_{n+1}})`` and the energy-drop envelope
    ``B_N = (dt / N) {F(c_0) - F(c_N)}``. The theorem check is ``D_min(N) <= B_N``.

    By default the theorem-safe stepsize ``dt = dt_safety / L_clip`` is used, with
    ``L_clip = beta * max(lambda_plus, lambda_plus^4 / lambda_minus^3)`` and
    ``dt_safety < 1`` (``dt_rule="theorem_safe"``). Passing an explicit ``dt``
    overrides this and the ``dt_rule`` label records which rule produced it; this
    is used to probe the deliberately non-theorem-safe stepsize
    ``dt = 1 / (beta * lambda_plus)`` (``dt * L_clip > 1``), which lies outside
    the Theorem 2.18 condition.
    """
    L_clip = clip_smoothness_constant(beta, lambda_minus, lambda_plus)
    if dt is None:
        dt = theorem_safe_dt(beta, lambda_minus, lambda_plus, dt_safety)
    dt = float(dt)
    F0 = target.energy(0.0, float(c0))

    rows = []
    c = float(c0)
    running_min = float("inf")
    first_upper = -1
    first_lower = -1
    energy_monotone = True
    denom_positive = True
    max_violation = 0.0
    min_margin = float("inf")

    for n in range(int(num_steps)):
        A = target.A0(c)
        F_cur = target.energy(0.0, c)
        step = clipped_kl_next(c, A, dt, lambda_minus, lambda_plus)
        denom = step.denom if step.denom is not None else float("nan")
        if not (np.isfinite(denom) and denom > 0.0):
            denom_positive = False
        if step.status != "ok":
            rows.append({
                "R": target.R, "lambda_minus": float(lambda_minus),
                "lambda_plus": float(lambda_plus), "beta": float(beta),
                "L_clip": L_clip, "dt_safety": float(dt_safety), "dt": dt,
                "dt_rule": str(dt_rule),
                "step": int(n), "c": float(c), "c_tilde": step.ctilde,
                "c_next": np.nan, "A": float(A), "F": float(F_cur),
                "D_kl_current_to_next": np.nan, "running_min_D": running_min,
                "prefix_bound": np.nan, "bound_margin": np.nan,
                "clipped_lower_active": 0, "clipped_upper_active": 0,
                "denominator": denom, "status": step.status,
            })
            break

        c_tilde = float(step.ctilde)
        c_next = float(step.c_next)
        lower_active = int(c_tilde < float(lambda_minus))
        upper_active = int(c_tilde > float(lambda_plus))
        if upper_active and first_upper < 0:
            first_upper = int(n)
        if lower_active and first_lower < 0:
            first_lower = int(n)

        D_n = gaussian_kl_divergence(c, c_next)
        running_min = min(running_min, D_n)
        F_next = target.energy(0.0, c_next)
        N = n + 1
        prefix_bound = (dt / N) * (F0 - F_next)
        bound_margin = prefix_bound - running_min
        min_margin = min(min_margin, bound_margin)
        max_violation = max(max_violation, max(0.0, running_min - prefix_bound))
        if F_next > F_cur + 1e-12:
            energy_monotone = False

        rows.append({
            "R": target.R, "lambda_minus": float(lambda_minus),
            "lambda_plus": float(lambda_plus), "beta": float(beta),
            "L_clip": L_clip, "dt_safety": float(dt_safety), "dt": dt,
            "dt_rule": str(dt_rule),
            "step": int(n), "c": float(c), "c_tilde": c_tilde,
            "c_next": c_next, "A": float(A), "F": float(F_cur),
            "D_kl_current_to_next": float(D_n), "running_min_D": float(running_min),
            "prefix_bound": float(prefix_bound), "bound_margin": float(bound_margin),
            "clipped_lower_active": lower_active,
            "clipped_upper_active": upper_active,
            "denominator": float(denom), "status": "ok",
        })
        c = c_next

    cs = [r["c"] for r in rows if np.isfinite(r["c"])]
    summary = {
        "R": target.R,
        "lambda_minus": float(lambda_minus),
        "lambda_plus": float(lambda_plus),
        "dt": dt,
        "dt_rule": str(dt_rule),
        "dt_times_L_clip": float(dt * L_clip),
        "num_steps": int(len(rows)),
        "final_c": rows[-1]["c_next"] if rows else np.nan,
        "min_c": float(min(cs)) if cs else np.nan,
        "max_c": float(max(cs)) if cs else np.nan,
        "final_running_min_D": float(running_min) if np.isfinite(running_min) else np.nan,
        "final_bound": rows[-1]["prefix_bound"] if rows else np.nan,
        "min_bound_margin": float(min_margin) if np.isfinite(min_margin) else np.nan,
        "max_violation": float(max_violation),
        "first_upper_clip_step": int(first_upper),
        "first_lower_clip_step": int(first_lower),
        "energy_monotone": bool(energy_monotone),
        "denominator_positive": bool(denom_positive),
        "theorem_check_pass": bool(denom_positive
                                   and max_violation <= CLIPPED_THEOREM_TOL),
    }
    return rows, summary


def _kl_pole_summary(kl_rows):
    first_rows = []
    for eps, group in kl_rows.groupby("epsilon", dropna=True):
        g = group.sort_values("n")
        r0 = g[g["n"] == 0].iloc[0]
        r1s = g[g["n"] == 1]
        status = str(g.iloc[-1]["stop_reason"] or g.iloc[-1]["status"])
        if len(r1s):
            r1 = r1s.iloc[0]
            c1 = float(r1["c"])
            fr1 = float(r1["fr_grad_sq"])
            bw1 = float(r1["bw_grad_sq"])
        else:
            c1 = fr1 = bw1 = np.nan
        ref_c = 2.0 / float(eps)
        ref_fr = 4.0 / (float(eps) ** 2)
        first_rows.append({
            "R": float(r0["R"]),
            "dt": float(r0["dt"]),
            "epsilon": float(eps),
            "c0": float(r0["c"]),
            "denom0": float(r0["denom"]),
            "c1": c1,
            "status": status,
            "fr_grad_sq_after_first": fr1,
            "bw_grad_sq_after_first": bw1,
            "expected_c1_2_over_epsilon": ref_c,
            "expected_fr_grad_sq_4_over_epsilon_sq": ref_fr,
            "c1_ratio_to_reference": c1 / ref_c if np.isfinite(c1) else np.nan,
            "fr_grad_sq_ratio_to_reference": fr1 / ref_fr if np.isfinite(fr1) else np.nan,
            "in_pole_fit_window": int(ref_c <= 0.02 * float(r0["R"]) ** 2),
            "fit_slope_log_fr_vs_log_epsilon": np.nan,
            "fit_slope_pole_window_log_fr_vs_log_epsilon": np.nan,
        })
    out = pd.DataFrame(first_rows).sort_values("epsilon", ascending=False)
    finite = out[np.isfinite(out.fr_grad_sq_after_first) & (out.fr_grad_sq_after_first > 0)]
    if len(finite) >= 2:
        slope, _ = np.polyfit(
            np.log(finite.epsilon.values),
            np.log(finite.fr_grad_sq_after_first.values),
            1,
        )
        out["fit_slope_log_fr_vs_log_epsilon"] = float(slope)
    pole = finite[finite.in_pole_fit_window == 1]
    if len(pole) >= 2:
        slope, _ = np.polyfit(
            np.log(pole.epsilon.values),
            np.log(pole.fr_grad_sq_after_first.values),
            1,
        )
        out["fit_slope_pole_window_log_fr_vs_log_epsilon"] = float(slope)
    return out.reindex(columns=KL_POLE_COLS)


def run_all(cfg, outdir):
    """Run all experiments and write the requested output files."""
    ensure_dir(outdir)
    gh_nodes = int(cfg.get("gh_nodes", 160))
    targets = {}
    target_meta = {
        "description": (
            "V_R(x)=0.5*x^2 - 2*R^2*log(cosh(x/R)); "
            "V''_R(x)=1-2*sech(x/R)^2 in [-1,1]."
        ),
        "quadrature": {"type": "Gauss-Hermite", "nodes": gh_nodes},
        "targets": {},
    }

    def get_target(R):
        R = float(R)
        if R not in targets:
            targets[R] = NonconvexLogCoshTarget(R, n_nodes=gh_nodes)
            target_meta["targets"][f"R={R:g}"] = targets[R].metadata()
        return targets[R]

    all_rows = []
    summaries = []

    rcfg = cfg.get("riemannian_cascade", {})
    for R in [float(x) for x in rcfg.get("R_values", [20, 50, 100, 300, 1000])]:
        rows, summary = simulate_riemannian_cascade(
            get_target(R),
            dt=float(rcfg.get("dt", 1.0)),
            c0=float(rcfg.get("c0", 1.0)),
            max_iter=int(rcfg.get("max_iter", 20)),
            stop_fraction=float(rcfg.get("stop_fraction", 0.25)),
            max_log_c=float(rcfg.get("max_log_c", 700.0)),
        )
        all_rows.extend(rows)
        summaries.append(summary)

    kcfg = cfg.get("kl_pole", {})
    ktarget = get_target(float(kcfg.get("R", 1000.0)))
    for eps in [float(x) for x in kcfg.get("epsilons", [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6])]:
        rows, summary = simulate_kl_pole(
            ktarget, epsilon=eps, dt=float(kcfg.get("dt", 1.0)),
            max_iter=int(kcfg.get("max_iter", 2)))
        all_rows.extend(rows)
        summaries.append(summary)

    wcfg = cfg.get("wasserstein_fb", {})
    opt_cfg = cfg.get("optimizer", {})
    star_by_R = {}
    bw_bound_rows = []
    for R in [float(x) for x in wcfg.get("R_values", [20, 50, 100, 300, 1000])]:
        target = get_target(R)
        star = estimate_F_star(
            target,
            maxiter=int(opt_cfg.get("maxiter", 300)),
            gtol=float(opt_cfg.get("gtol", 1e-10)),
        )
        star_by_R[R] = star
        target_meta["targets"][f"R={R:g}"]["F_star_num"] = star["F_star_num"]
        target_meta["targets"][f"R={R:g}"]["optimizer"] = star
        rows, summary, bnd = simulate_wasserstein_fb(
            target,
            star,
            eta=float(wcfg.get("eta", 0.9)),
            c0=float(wcfg.get("c0", 1.0)),
            max_iter=int(wcfg.get("max_iter", 20000)),
        )
        all_rows.extend(rows)
        summaries.append(summary)
        bw_bound_rows.extend(bnd)

    ccfg = cfg.get("clipped_kl_stationarity", {})
    clipped_rows = []
    clipped_summaries = []
    c_lam_m = float(ccfg.get("lambda_minus", 0.5))
    c_lam_p = float(ccfg.get("lambda_plus", 2.0))
    c_beta = float(ccfg.get("beta", 1.0))
    c_dt_safety = float(ccfg.get("dt_safety", 0.9))
    c_c0 = float(ccfg.get("c0", 1.0))
    c_steps = int(ccfg.get("num_steps", 300))
    # Optional non-theorem-safe stepsize study: dt = 1 / (beta * lambda_plus),
    # which gives dt * L_clip > 1 and so violates the Theorem 2.18 condition.
    compare_rs = bool(ccfg.get("compare_riemannian_scale_dt", False))
    for R in [float(x) for x in ccfg.get("R_values", [20, 50, 100, 300, 1000])]:
        target = get_target(R)
        rows, summary = simulate_clipped_kl(
            target, lambda_minus=c_lam_m, lambda_plus=c_lam_p, beta=c_beta,
            dt_safety=c_dt_safety, c0=c_c0, num_steps=c_steps,
            dt_rule="theorem_safe",
        )
        clipped_rows.extend(rows)
        clipped_summaries.append(summary)
        if compare_rs:
            rs_rows, rs_summary = simulate_clipped_kl(
                target, lambda_minus=c_lam_m, lambda_plus=c_lam_p, beta=c_beta,
                dt_safety=c_dt_safety, c0=c_c0, num_steps=c_steps,
                dt=1.0 / (c_beta * c_lam_p), dt_rule="riemannian_scale",
            )
            clipped_rows.extend(rs_rows)
            clipped_summaries.append(rs_summary)

    long_df = pd.DataFrame(all_rows).reindex(columns=RESULTS_LONG_COLS)
    summary_df = pd.DataFrame(summaries).reindex(columns=SUMMARY_COLS)
    kl_df = _kl_pole_summary(long_df[long_df.experiment == "kl_pole"])
    bw_df = pd.DataFrame(bw_bound_rows).reindex(columns=BW_BOUND_COLS)
    clipped_df = pd.DataFrame(clipped_rows).reindex(columns=CLIPPED_KL_COLS)
    clipped_summary_df = pd.DataFrame(clipped_summaries).reindex(
        columns=CLIPPED_SUMMARY_COLS)

    save_dataframe(os.path.join(outdir, "results_long.csv"), long_df)
    save_dataframe(os.path.join(outdir, "summary.csv"), summary_df)
    save_dataframe(os.path.join(outdir, "kl_pole_summary.csv"), kl_df)
    save_dataframe(os.path.join(outdir, "wasserstein_bound_summary.csv"), bw_df)
    save_dataframe(os.path.join(outdir, "clipped_kl_stationarity.csv"), clipped_df)
    save_dataframe(os.path.join(outdir, "clipped_kl_summary.csv"), clipped_summary_df)
    save_json(os.path.join(outdir, "target_metadata.json"), target_meta)
    save_json(os.path.join(outdir, "run_metadata.json"), {
        "config": cfg,
        "python": sys.version,
        "platform": platform.platform(),
        "created_at_unix": time.time(),
        "output_files": [
            "results_long.csv",
            "summary.csv",
            "kl_pole_summary.csv",
            "wasserstein_bound_summary.csv",
            "clipped_kl_stationarity.csv",
            "clipped_kl_summary.csv",
            "target_metadata.json",
            "run_metadata.json",
        ],
    })
    return (long_df, summary_df, kl_df, bw_df, clipped_df,
            clipped_summary_df, target_meta)
