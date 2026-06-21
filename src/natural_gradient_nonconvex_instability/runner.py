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
from src.common.theory_constants import (
    THEORY_VERSION_PROJECTED_KL,
    projected_kl_theory_constants,
)
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
    "R", "beta", "lambda_minus", "lambda_plus", "L_clip",
    "dt_projected_KL_theory", "safety", "dt_used", "dt_rule", "step",
    "c", "c_tilde", "c_next", "A", "F", "D_kl_current_to_next", "running_min_D",
    "prefix_bound", "bound_margin", "clipped_lower_active", "clipped_upper_active",
    "denominator", "status", "theory_version",
]

CLIPPED_SUMMARY_COLS = [
    "R", "beta", "lambda_minus", "lambda_plus", "L_clip",
    "dt_projected_KL_theory", "dt_used", "dt_rule", "safety", "dt_times_L_clip",
    "num_steps", "final_c", "min_c", "max_c", "final_D_min", "final_B_N",
    "min_bound_margin", "max_violation", "theorem_check_pass",
    "boundary_active_fraction", "first_upper_clip_step", "first_lower_clip_step",
    "energy_monotone", "denominator_positive", "theory_version",
]

# Theory sweep CSV: shows dt_theory is independent of lambda_minus and scales
# like 1/lambda_plus (projected-KL theorem with L_clip = 2 beta lambda_plus).
CLIPPED_SWEEP_COLS = [
    "sweep_kind", "beta", "lambda_minus", "lambda_plus", "L_clip",
    "dt_projected_KL_theory", "theory_version",
]

# Absolute tolerance on the theorem-check residual (running_min_D - prefix_bound).
CLIPPED_THEOREM_TOL = 1.0e-9

# dt = 1/(beta*lambda_plus) is exactly 2x the projected-KL theorem edge
# 1/(2*beta*lambda_plus), so dt*L_clip = 2: outside the theorem, by a factor 2.
OUTSIDE_THEOREM_2X_RULE = "outside_theorem_2x"


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
                        safety=0.5, c0=1.0, num_steps=300, dt=None,
                        dt_rule="theorem_safe"):
    """Projected KL covariance scheme with eigenvalue clipping.

    The covariance is evolved by the unprojected KL step and then clipped back
    into ``[lambda_minus, lambda_plus]``. For each prefix horizon ``N`` we record
    the running-minimum Bregman displacement ``D_min(N) = min_{0<=n<N} D_n`` with
    ``D_n = KL(rho_{a_n} || rho_{a_{n+1}})`` and the energy-drop envelope
    ``B_N = (dt / N) {F(c_0) - F(c_N)}``. The theorem check is ``D_min(N) <= B_N``.

    Current projected-KL theorem (global Hessian bound ``||grad^2 V|| <= beta``):

        L_clip = 2 * beta * lambda_plus,
        dt_projected_KL_theory = 1 / L_clip = 1 / (2 * beta * lambda_plus).

    The theorem-safe scale depends on ``lambda_plus`` and ``beta`` only, not on
    ``lambda_minus``. By default ``dt = safety / L_clip`` with ``safety <= 1``
    (``dt_rule="theorem_safe"``). Passing an explicit ``dt`` overrides this; it is
    used to probe the deliberately non-theorem-safe stepsize
    ``dt = 1 / (beta * lambda_plus)``, which is exactly twice the theorem edge so
    ``dt * L_clip = 2`` (``dt_rule="outside_theorem_2x"``); the theorem gives no
    guarantee there.

    The clipped certificate is a constrained Bregman stationarity guarantee over
    the covariance-truncated feasible set ``[lambda_minus, lambda_plus]``. It does
    NOT certify that the unconstrained Fisher-Rao gradient is small.
    """
    consts = projected_kl_theory_constants(beta, lambda_plus)
    L_clip = consts["L_clip"]
    dt_theory = consts["dt_projected_KL_theory"]
    if dt is None:
        dt = theorem_safe_dt(beta, lambda_plus, safety)
    dt = float(dt)
    F0 = target.energy(0.0, float(c0))

    rows = []
    c = float(c0)
    running_min = float("inf")
    first_upper = -1
    first_lower = -1
    n_boundary_active = 0
    energy_monotone = True
    denom_positive = True
    max_violation = 0.0
    min_margin = float("inf")

    def _base(n):
        return {
            "R": target.R, "beta": float(beta),
            "lambda_minus": float(lambda_minus), "lambda_plus": float(lambda_plus),
            "L_clip": float(L_clip), "dt_projected_KL_theory": float(dt_theory),
            "safety": float(safety), "dt_used": dt, "dt_rule": str(dt_rule),
            "step": int(n), "theory_version": THEORY_VERSION_PROJECTED_KL,
        }

    for n in range(int(num_steps)):
        A = target.A0(c)
        F_cur = target.energy(0.0, c)
        step = clipped_kl_next(c, A, dt, lambda_minus, lambda_plus)
        denom = step.denom if step.denom is not None else float("nan")
        if not (np.isfinite(denom) and denom > 0.0):
            denom_positive = False
        if step.status != "ok":
            row = _base(n)
            row.update({
                "c": float(c), "c_tilde": step.ctilde, "c_next": np.nan,
                "A": float(A), "F": float(F_cur),
                "D_kl_current_to_next": np.nan, "running_min_D": running_min,
                "prefix_bound": np.nan, "bound_margin": np.nan,
                "clipped_lower_active": 0, "clipped_upper_active": 0,
                "denominator": denom, "status": step.status,
            })
            rows.append(row)
            break

        c_tilde = float(step.ctilde)
        c_next = float(step.c_next)
        lower_active = int(c_tilde < float(lambda_minus))
        upper_active = int(c_tilde > float(lambda_plus))
        if upper_active and first_upper < 0:
            first_upper = int(n)
        if lower_active and first_lower < 0:
            first_lower = int(n)
        if lower_active or upper_active:
            n_boundary_active += 1

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

        row = _base(n)
        row.update({
            "c": float(c), "c_tilde": c_tilde, "c_next": c_next,
            "A": float(A), "F": float(F_cur),
            "D_kl_current_to_next": float(D_n), "running_min_D": float(running_min),
            "prefix_bound": float(prefix_bound), "bound_margin": float(bound_margin),
            "clipped_lower_active": lower_active,
            "clipped_upper_active": upper_active,
            "denominator": float(denom), "status": "ok",
        })
        rows.append(row)
        c = c_next

    cs = [r["c"] for r in rows if np.isfinite(r["c"])]
    n_recorded = len(rows)
    summary = {
        "R": target.R,
        "beta": float(beta),
        "lambda_minus": float(lambda_minus),
        "lambda_plus": float(lambda_plus),
        "L_clip": float(L_clip),
        "dt_projected_KL_theory": float(dt_theory),
        "dt_used": dt,
        "dt_rule": str(dt_rule),
        "safety": float(safety),
        "dt_times_L_clip": float(dt * L_clip),
        "num_steps": int(n_recorded),
        "final_c": rows[-1]["c_next"] if rows else np.nan,
        "min_c": float(min(cs)) if cs else np.nan,
        "max_c": float(max(cs)) if cs else np.nan,
        "final_D_min": float(running_min) if np.isfinite(running_min) else np.nan,
        "final_B_N": rows[-1]["prefix_bound"] if rows else np.nan,
        "min_bound_margin": float(min_margin) if np.isfinite(min_margin) else np.nan,
        "max_violation": float(max_violation),
        "theorem_check_pass": bool(denom_positive
                                   and max_violation <= CLIPPED_THEOREM_TOL),
        "boundary_active_fraction": (float(n_boundary_active) / n_recorded
                                     if n_recorded else np.nan),
        "first_upper_clip_step": int(first_upper),
        "first_lower_clip_step": int(first_lower),
        "energy_monotone": bool(energy_monotone),
        "denominator_positive": bool(denom_positive),
        "theory_version": THEORY_VERSION_PROJECTED_KL,
    }
    return rows, summary


def clipped_kl_theory_sweep(*, beta=1.0, lambda_minus_fixed=0.5,
                            lambda_plus_fixed=2.0,
                            lambda_minus_values=(0.1, 0.25, 0.5, 1.0),
                            lambda_plus_values=(1.0, 2.0, 4.0, 8.0)):
    """Demonstrate the projected-KL theorem-safe stepsize scaling.

    Two sweeps over the new constant ``L_clip = 2 beta lambda_plus``:

    * ``vary_lambda_minus`` (``lambda_plus`` fixed): ``dt_projected_KL_theory`` is
      unchanged -- the theorem-safe scale does not depend on ``lambda_minus``.
    * ``vary_lambda_plus`` (``lambda_minus`` fixed): ``dt_projected_KL_theory``
      scales like ``1 / lambda_plus``.
    """
    rows = []
    for lm in lambda_minus_values:
        consts = projected_kl_theory_constants(beta, lambda_plus_fixed)
        rows.append({
            "sweep_kind": "vary_lambda_minus", "beta": float(beta),
            "lambda_minus": float(lm), "lambda_plus": float(lambda_plus_fixed),
            "L_clip": consts["L_clip"],
            "dt_projected_KL_theory": consts["dt_projected_KL_theory"],
            "theory_version": THEORY_VERSION_PROJECTED_KL,
        })
    for lp in lambda_plus_values:
        consts = projected_kl_theory_constants(beta, lp)
        rows.append({
            "sweep_kind": "vary_lambda_plus", "beta": float(beta),
            "lambda_minus": float(lambda_minus_fixed), "lambda_plus": float(lp),
            "L_clip": consts["L_clip"],
            "dt_projected_KL_theory": consts["dt_projected_KL_theory"],
            "theory_version": THEORY_VERSION_PROJECTED_KL,
        })
    return rows


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
    # Theorem-safe safety fractions for dt = safety / L_clip (L_clip = 2 beta lambda_+).
    safety_values = [float(s) for s in ccfg.get("safety_values", [0.25, 0.5, 0.9])]
    primary_safety = float(ccfg.get("dt_safety", 0.5))
    c_c0 = float(ccfg.get("c0", 1.0))
    c_steps = int(ccfg.get("num_steps", 300))
    # Outside-theorem comparison: dt = 1/(beta*lambda_plus) = 2x the theorem edge
    # 1/(2 beta lambda_plus), so dt*L_clip = 2 (outside_theorem_2x). The theorem
    # gives no guarantee here; we record whether the envelope still holds.
    compare_outside = bool(ccfg.get("compare_outside_theorem_2x", True))
    for R in [float(x) for x in ccfg.get("R_values", [20, 50, 100, 300, 1000])]:
        target = get_target(R)
        for safety in safety_values:
            # The primary safety is labelled "theorem_safe" (used by the figures);
            # the others carry a safety-tagged rule for the multi-safety table.
            rule = ("theorem_safe" if safety == primary_safety
                    else f"theorem_safe_s{safety:g}")
            rows, summary = simulate_clipped_kl(
                target, lambda_minus=c_lam_m, lambda_plus=c_lam_p, beta=c_beta,
                safety=safety, c0=c_c0, num_steps=c_steps, dt_rule=rule,
            )
            clipped_rows.extend(rows)
            clipped_summaries.append(summary)
        if compare_outside:
            os_rows, os_summary = simulate_clipped_kl(
                target, lambda_minus=c_lam_m, lambda_plus=c_lam_p, beta=c_beta,
                c0=c_c0, num_steps=c_steps,
                dt=1.0 / (c_beta * c_lam_p), dt_rule=OUTSIDE_THEOREM_2X_RULE,
            )
            clipped_rows.extend(os_rows)
            clipped_summaries.append(os_summary)

    # Theory-only sweep: dt_theory is independent of lambda_minus and scales like
    # 1/lambda_plus (projected-KL theorem, L_clip = 2 beta lambda_plus).
    scfg = ccfg.get("theory_sweep", {})
    sweep_rows = clipped_kl_theory_sweep(
        beta=c_beta,
        lambda_minus_fixed=c_lam_m, lambda_plus_fixed=c_lam_p,
        lambda_minus_values=[float(x) for x in
                             scfg.get("lambda_minus_values", [0.1, 0.25, 0.5, 1.0])],
        lambda_plus_values=[float(x) for x in
                            scfg.get("lambda_plus_values", [1.0, 2.0, 4.0, 8.0])],
    )

    long_df = pd.DataFrame(all_rows).reindex(columns=RESULTS_LONG_COLS)
    summary_df = pd.DataFrame(summaries).reindex(columns=SUMMARY_COLS)
    kl_df = _kl_pole_summary(long_df[long_df.experiment == "kl_pole"])
    bw_df = pd.DataFrame(bw_bound_rows).reindex(columns=BW_BOUND_COLS)
    clipped_df = pd.DataFrame(clipped_rows).reindex(columns=CLIPPED_KL_COLS)
    clipped_summary_df = pd.DataFrame(clipped_summaries).reindex(
        columns=CLIPPED_SUMMARY_COLS)
    clipped_sweep_df = pd.DataFrame(sweep_rows).reindex(columns=CLIPPED_SWEEP_COLS)

    save_dataframe(os.path.join(outdir, "results_long.csv"), long_df)
    save_dataframe(os.path.join(outdir, "summary.csv"), summary_df)
    save_dataframe(os.path.join(outdir, "kl_pole_summary.csv"), kl_df)
    save_dataframe(os.path.join(outdir, "wasserstein_bound_summary.csv"), bw_df)
    save_dataframe(os.path.join(outdir, "clipped_kl_stationarity.csv"), clipped_df)
    save_dataframe(os.path.join(outdir, "clipped_kl_summary.csv"), clipped_summary_df)
    save_dataframe(os.path.join(outdir, "clipped_kl_sweep.csv"), clipped_sweep_df)
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
            "clipped_kl_sweep.csv",
            "target_metadata.json",
            "run_metadata.json",
        ],
    })
    return (long_df, summary_df, kl_df, bw_df, clipped_df,
            clipped_summary_df, clipped_sweep_df, target_meta)
