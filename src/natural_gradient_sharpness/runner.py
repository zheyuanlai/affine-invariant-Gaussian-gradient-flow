"""Experiment orchestration and rate diagnostics for sharpness constructions."""
from __future__ import annotations

import json
import math
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
from scipy import integrate

from .bump_train import BumpTrainTarget
from .local_targets import (
    ChenLogBump,
    RidgeTarget,
    ShellTarget,
    kl_jacobian,
)
from .profiles import IntegratedFlatTop
from .spiral import SpiralValleyTarget


def loglog_slope(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if np.sum(mask) < 2:
        return float("nan"), float("nan")
    lx, ly = np.log(x[mask]), np.log(y[mask])
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = slope * lx + intercept
    denom = np.sum((ly - np.mean(ly)) ** 2)
    r2 = 1.0 - np.sum((ly - pred) ** 2) / denom if denom > 0 else 1.0
    return float(slope), float(r2)


def _first_crossing(x, y, threshold):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    idx = np.flatnonzero(y <= threshold)
    if idx.size == 0:
        return float("nan")
    j = int(idx[0])
    if j == 0:
        return float(x[0])
    x0, x1, y0, y1 = x[j - 1], x[j], y[j - 1], y[j]
    if y1 == y0:
        return float(x1)
    return float(x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))


def run_spiral_grid(config):
    summaries, long_rows = [], []
    for K in config["K"]:
        target = SpiralValleyTarget(
            K, radius_factor=config.get("radius_factor", 200.0),
            gh_order=config.get("gh_order", 10))
        t_end = config.get("time_fraction", 0.012) * K
        t_eval = np.linspace(0.0, t_end, config.get("num_times", 241))
        sol = integrate.solve_ivp(
            target.rhs, (0.0, t_end), target.state0, t_eval=t_eval,
            rtol=config.get("rtol", 2e-9), atol=config.get("atol", 1e-10),
            max_step=K / config.get("max_steps_per_K", 1800.0))
        if not sol.success:
            raise RuntimeError(sol.message)
        # ``solve_ivp`` may evaluate and reject trial stages outside the
        # accepted trajectory. Recompute containment only on accepted states.
        target.max_normalized_phase = 0.0
        target.min_radius_over_r0 = math.inf
        for accepted in sol.y.T:
            C_acc = np.array([
                [accepted[2], accepted[3]],
                [accepted[3], accepted[4]],
            ])
            target.expectations(accepted[:2], C_acc)
        radius = np.linalg.norm(sol.y[:2].T, axis=1)
        ratio = radius / radius[0]
        fit = np.polyfit(sol.t, np.log(radius), 1)
        time_98 = _first_crossing(sol.t, ratio, config.get("radius_threshold", 0.98))
        eig_min = np.array([
            np.linalg.eigvalsh(np.array([[row[2], row[3]], [row[3], row[4]]]))[0]
            for row in sol.y.T
        ])
        summaries.append({
            "construction": "logarithmic_spiral",
            "K": K,
            "kappa_star_lower": float(K),
            "kappa_star_upper": 1.03 * K,
            "t_end": t_end,
            "radius_ratio_end": float(ratio[-1]),
            "certified_mean_gap_ratio_end": float(ratio[-1] ** 2),
            "fitted_radial_rate": float(-fit[0]),
            "K_times_fitted_rate": float(-fit[0] * K),
            "time_to_radius_threshold": time_98,
            "time_over_K": time_98 / K,
            "max_normalized_phase": target.max_normalized_phase,
            "min_radius_over_r0": target.min_radius_over_r0,
            "min_covariance_eigenvalue": float(np.min(eig_min)),
            "containment_pass": bool(target.max_normalized_phase < 1.0 and target.min_radius_over_r0 > 1.0),
        })
        save_every = max(1, len(sol.t) // config.get("max_saved_rows", 160))
        for j in range(0, len(sol.t), save_every):
            long_rows.append({
                "K": K, "time": sol.t[j], "radius": radius[j],
                "radius_ratio": ratio[j], "certified_mean_gap_ratio": ratio[j] ** 2,
            })
        if long_rows[-1]["time"] != sol.t[-1]:
            long_rows.append({
                "K": K, "time": sol.t[-1], "radius": radius[-1],
                "radius_ratio": ratio[-1], "certified_mean_gap_ratio": ratio[-1] ** 2,
            })
    return pd.DataFrame(summaries), pd.DataFrame(long_rows)


def run_bump_grid(config, profile):
    summaries, long_rows = [], []
    relative_threshold = float(config.get("gap_threshold", 0.8))
    for kappa in config["kappa"]:
        target = BumpTrainTarget(
            kappa, step_fraction=config.get("step_fraction", 0.5),
            train_time=config.get("train_time", 0.15),
            initial_scale=config.get("initial_scale", 100.0),
            lipschitz=config.get("lipschitz", 1.0),
            gh_order=config.get("gh_order", 32), profile=profile)
        for method in ("riemannian", "kl"):
            state = np.array([target.mean0, target.cov0])
            gap0 = target.gap(*state)
            gaps = [gap0]
            mean_error, cov_error = [0.0], [0.0]
            states = [state.copy()]
            for j in range(target.num_train_steps):
                state = target.step(method, state)
                states.append(state.copy())
                gaps.append(target.gap(*state))
                mean_error.append(abs(state[0] - target.centers_path[j + 1]))
                cov_error.append(abs(kappa * state[1] - 1.0))
            gaps = np.asarray(gaps)
            ratios = gaps / gap0
            hit = _first_crossing(np.arange(len(ratios)), ratios, relative_threshold)
            summaries.append({
                "construction": "flat_top_bump_train",
                "kappa": kappa,
                "method": method,
                "dt": target.dt,
                "kappa_over_dt": kappa / target.dt,
                "num_train_steps": target.num_train_steps,
                "train_steps_over_kappa2": target.num_train_steps / kappa ** 2,
                "final_gap_ratio": float(ratios[-1]),
                "gap_threshold": relative_threshold,
                "steps_to_gap_threshold": hit,
                "steps_to_threshold_over_kappa2": hit / kappa ** 2,
                "max_center_shadow_error_over_width": float(np.max(mean_error) / target.width),
                "max_relative_covariance_error": float(np.max(cov_error)),
                "min_spacing_over_width": target.min_spacing_over_width,
                "hessian_lipschitz": target.lipschitz,
                "shadowing_pass": bool(np.max(mean_error) < target.width / 8.0 and np.max(cov_error) < 0.125),
            })
            save_every = max(1, target.num_train_steps // config.get("max_saved_rows", 180))
            for j in range(0, target.num_train_steps + 1, save_every):
                long_rows.append({
                    "kappa": kappa, "method": method, "step": j,
                    "time": j * target.dt, "mean": states[j][0],
                    "covariance": states[j][1], "gap_ratio": ratios[j],
                })
    return pd.DataFrame(summaries), pd.DataFrame(long_rows)


def _local_row(target, regime, tuning_parameter):
    if regime == "fixed_dim_log":
        predicted = math.log(math.e * target.kappa)
        n = 1
        slow_block = "scalar"
        alpha, beta = target.alpha, target.beta
    elif regime == "dimension_branch":
        predicted = math.sqrt(target.n * math.log(math.e * target.kappa))
        n, slow_block = target.n, target.slow_block
        alpha, beta = target.alpha, target.beta
    else:
        predicted = math.sqrt(target.kappa)
        n, slow_block = target.n, target.slow_block
        alpha, beta = target.alpha, target.beta
    return {
        "construction": target.__class__.__name__,
        "regime": regime,
        "tuning_parameter": tuning_parameter,
        "n": n,
        "alpha": alpha,
        "beta": beta,
        "kappa": target.kappa,
        "gamma": target.gamma,
        "inverse_gap": 1.0 / target.gamma,
        "Lambda": target.Lambda,
        "predicted_inverse_gap_scale": predicted,
        "inverse_gap_over_predicted": 1.0 / (target.gamma * predicted),
        "slow_block": slow_block,
        "isotropy_error": getattr(target, "isotropy_error", 0.0),
        "hessian_lipschitz_certificate": getattr(target, "hessian_lipschitz_certificate", float("nan")),
        "Delta": getattr(target, "Delta", float("nan")),
        "tail_probability": getattr(target, "tail_probability", float("nan")),
    }


def run_local_spectral(config):
    rows, objects = [], []
    for location in config["chen_locations"]:
        target = ChenLogBump(float(location))
        rows.append(_local_row(target, "fixed_dim_log", location))
    qcfg = config.get("quadrature", {})
    for K in config["ridge_K"]:
        target = RidgeTarget(K, **qcfg)
        rows.append(_local_row(target, "uniform_saturation", K))
        objects.append(target)
    for m in config["shell_subcritical_m"]:
        depth = max(1.0, float(config.get("shell_depth_factor", 1.0)) * math.log(m))
        target = ShellTarget(m, depth, **qcfg)
        rows.append(_local_row(target, "dimension_branch", depth))
        objects.append(target)
    for m in config["shell_saturation_m"]:
        target = ShellTarget(m, 1.0, **qcfg)
        rows.append(_local_row(target, "uniform_saturation", 1.0))
        objects.append(target)
    return pd.DataFrame(rows), objects


def _norm_from_state(target, state, dt=0.0, kl=False):
    x = target.scalar_coordinates(state)
    if kl:
        return float(math.sqrt(x[0] ** 2 + (1.0 + dt) * np.dot(x[1:], x[1:])))
    return float(np.linalg.norm(x))


def run_local_dynamics(objects, config):
    rows = []
    allowed = set(config.get("constructions", ["ridge", "shell"]))
    max_targets = int(config.get("max_targets_per_construction", 2))
    max_n = int(config.get("max_n", 10**18))
    selected = []
    for name in allowed:
        candidates = [obj for obj in objects if obj.construction == name and obj.n <= max_n]
        if candidates:
            indices = np.unique(np.linspace(0, len(candidates) - 1, min(max_targets, len(candidates))).astype(int))
            selected.extend(candidates[i] for i in indices)
    eps = float(config.get("perturbation", 2e-5))
    horizon = float(config.get("decay_exponents", 2.0))
    safety = float(config.get("stepsize_safety", 0.5))
    for target in selected:
        L = target.blocks["scalar"]
        eigvals, eigvecs = np.linalg.eigh(L)
        v = eigvecs[:, 0]
        state0 = target.state_from_scalar_coordinates(v, eps)
        t_end = horizon / target.gamma
        t_eval = np.linspace(0.0, t_end, int(config.get("continuous_points", 80)))
        sol = integrate.solve_ivp(
            target.reduced_rhs, (0.0, t_end), state0, t_eval=t_eval,
            rtol=2e-9, atol=1e-12, max_step=t_end / 250.0)
        norms = np.array([_norm_from_state(target, state) for state in sol.y.T])
        rate = -math.log(norms[-1] / norms[0]) / t_end
        rows.append({
            "construction": target.construction, "n": target.n, "kappa": target.kappa,
            "method": "continuous", "dt": 0.0, "steps": 0,
            "spectral_gamma": target.gamma, "spectral_Lambda": target.Lambda,
            "predicted_contraction": math.exp(-target.gamma * t_end),
            "observed_contraction": norms[-1] / norms[0],
            "observed_rate": rate, "observed_over_gamma": rate / target.gamma,
            "effective_iterations": float("nan"),
        })
        dt = safety / target.Lambda
        for method in ("riemannian", "kl"):
            if method == "riemannian":
                J = np.eye(3) - dt * L
                values, vectors = np.linalg.eigh(J)
                idx = int(np.argmax(np.abs(values)))
                q_pred, v_method = float(abs(values[idx])), vectors[:, idx]
            else:
                J = kl_jacobian(L, dt, mean_coordinates=1)
                values, vectors = np.linalg.eig(J)
                idx = int(np.argmax(np.abs(values)))
                q_pred = float(abs(values[idx]))
                v_method = np.real(vectors[:, idx])
                v_method /= np.linalg.norm(v_method)
            state = target.state_from_scalar_coordinates(v_method, eps)
            norm0 = _norm_from_state(target, state, dt, method == "kl")
            n_steps = int(math.ceil(horizon / max(1e-15, 1.0 - q_pred)))
            for _ in range(n_steps):
                state = target.reduced_step(method, state, dt)
            normN = _norm_from_state(target, state, dt, method == "kl")
            q_obs = (normN / norm0) ** (1.0 / n_steps)
            rows.append({
                "construction": target.construction, "n": target.n, "kappa": target.kappa,
                "method": method, "dt": dt, "steps": n_steps,
                "spectral_gamma": target.gamma, "spectral_Lambda": target.Lambda,
                "predicted_contraction": q_pred,
                "observed_contraction": q_obs,
                "observed_rate": -math.log(q_obs) / dt,
                "observed_over_gamma": (-math.log(q_obs) / dt) / target.gamma,
                "effective_iterations": 1.0 / max(1e-15, 1.0 - q_obs),
            })
    return pd.DataFrame(rows)


def _simulate_discrete_to_entry(target, method, entry_radius, max_steps, save_every=1):
    state = np.array([target.mean0, target.cov0])
    records = [(0, state.copy(), target.error_norm(state))]
    for n in range(1, max_steps + 1):
        state = target.step(method, state)
        error = target.error_norm(state)
        if n % save_every == 0 or error <= entry_radius:
            records.append((n, state.copy(), error))
        if error <= entry_radius:
            return n, state, records
    raise RuntimeError(f"{method} did not enter the local ball")


def run_global_to_local(config, profile):
    kappa = int(config.get("kappa", 16))
    target = BumpTrainTarget(
        kappa, step_fraction=config.get("step_fraction", 0.5),
        train_time=config.get("train_time", 0.15),
        initial_scale=config.get("initial_scale", 100.0),
        lipschitz=config.get("lipschitz", 1.0),
        gh_order=config.get("gh_order", 32), profile=profile)
    entry = float(config.get("entry_radius", 0.1))
    tol = float(config.get("tolerance", 1e-6))
    local_dt = float(config.get("local_dt", 0.5))
    summaries, rows = [], []

    # Continuous flow: the same target exhibits an O(kappa) far-field phase and
    # a unit-rate quadratic local phase.
    t_end = target.train_time * kappa + math.log(target.mean0) + math.log(kappa) + 15.0
    t_eval = np.linspace(0.0, t_end, int(config.get("continuous_points", 700)))
    sol = integrate.solve_ivp(
        target.rhs, (0.0, t_end), [target.mean0, target.cov0], t_eval=t_eval,
        rtol=2e-8, atol=1e-10, max_step=config.get("continuous_max_step", 0.04))
    errors = np.array([target.error_norm(s) for s in sol.y.T])
    entry_time = _first_crossing(sol.t, errors, entry)
    tol_time = _first_crossing(sol.t, errors, tol)
    summaries.append({
        "method": "continuous", "kappa": kappa, "global_dt": 0.0,
        "entry_time_or_steps": entry_time, "fixed_total": tol_time,
        "two_stage_total": tol_time, "tail_savings": 0.0,
    })
    for t, state, error in zip(sol.t, sol.y.T, errors):
        rows.append({
            "method": "continuous", "variant": "continuous", "progress": t,
            "elapsed_time": t, "error": error,
            "phase": "global" if t < entry_time else "local",
        })

    max_entry_steps = target.num_train_steps + int((math.log(target.mean0) + math.log(kappa) + 20.0) / target.dt)
    for method in ("riemannian", "kl"):
        n_entry, entry_state, prefix = _simulate_discrete_to_entry(
            target, method, entry, max_entry_steps,
            save_every=max(1, max_entry_steps // config.get("max_saved_rows", 350)))
        for n, _state, error in prefix:
            rows.append({
                "method": method, "variant": "shared_global", "progress": n,
                "elapsed_time": n * target.dt, "error": error, "phase": "global",
            })

        totals = {}
        for variant, dt in (("fixed", target.dt), ("two_stage", local_dt)):
            state = entry_state.copy()
            n_tail = 0
            max_tail = int(math.ceil((math.log(entry / tol) + 10.0) / min(dt, 0.9))) + 10
            while target.error_norm(state) > tol and n_tail < max_tail:
                state = target.step(method, state, dt=dt)
                n_tail += 1
                if n_tail <= 5 or n_tail % max(1, max_tail // 160) == 0 or target.error_norm(state) <= tol:
                    rows.append({
                        "method": method, "variant": variant,
                        "progress": n_entry + n_tail,
                        "elapsed_time": n_entry * target.dt + n_tail * dt,
                        "error": target.error_norm(state), "phase": "local",
                    })
            if target.error_norm(state) > tol:
                raise RuntimeError(f"{method}/{variant} failed to reach tolerance")
            totals[variant] = n_entry + n_tail
        summaries.append({
            "method": method, "kappa": kappa, "global_dt": target.dt,
            "entry_time_or_steps": n_entry, "fixed_total": totals["fixed"],
            "two_stage_total": totals["two_stage"],
            "tail_savings": totals["fixed"] - totals["two_stage"],
        })
    return pd.DataFrame(summaries), pd.DataFrame(rows)


def summarize_rates(spiral, bump, local):
    out = {}
    out["continuous_global_kappa_slope"], out["continuous_global_kappa_r2"] = loglog_slope(
        spiral["kappa_star_lower"], spiral["time_to_radius_threshold"])
    for method in ("riemannian", "kl"):
        frame = bump[bump.method == method]
        out[f"{method}_global_kappa2_slope"], out[f"{method}_global_kappa2_r2"] = loglog_slope(
            frame["kappa"], frame["steps_to_gap_threshold"])
    for regime in ("fixed_dim_log", "dimension_branch"):
        frame = local[local.regime == regime].sort_values("predicted_inverse_gap_scale")
        out[f"{regime}_predicted_scale_slope"], out[f"{regime}_predicted_scale_r2"] = loglog_slope(
            frame["predicted_inverse_gap_scale"], frame["inverse_gap"])
        tail = frame.iloc[len(frame) // 2:]
        out[f"{regime}_tail_slope"], out[f"{regime}_tail_r2"] = loglog_slope(
            tail["predicted_inverse_gap_scale"], tail["inverse_gap"])
    uniform = local[local.regime == "uniform_saturation"]
    for construction, frame in uniform.groupby("construction"):
        frame = frame.sort_values("predicted_inverse_gap_scale")
        key = construction.replace("Target", "").lower()
        out[f"uniform_saturation_{key}_predicted_scale_slope"], out[f"uniform_saturation_{key}_predicted_scale_r2"] = loglog_slope(
            frame["predicted_inverse_gap_scale"], frame["inverse_gap"])
        tail = frame.iloc[len(frame) // 2:]
        out[f"uniform_saturation_{key}_tail_slope"], out[f"uniform_saturation_{key}_tail_r2"] = loglog_slope(
            tail["predicted_inverse_gap_scale"], tail["inverse_gap"])
    return out


def run_all(config, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    profile = IntegratedFlatTop(config.get("profile_grid_size", 20001))
    spiral, spiral_long = run_spiral_grid(config["global_continuous"])
    bump, bump_long = run_bump_grid(config["global_discrete"], profile)
    local, objects = run_local_spectral(config["local_spectral"])
    local_dyn = run_local_dynamics(objects, config["local_dynamics"])
    gtl, gtl_long = run_global_to_local(config["global_to_local"], profile)
    rate_summary = summarize_rates(spiral, bump, local)

    outputs = {
        "global_continuous_summary.csv": spiral,
        "global_continuous_trajectories.csv": spiral_long,
        "global_discrete_summary.csv": bump,
        "global_discrete_trajectories.csv": bump_long,
        "local_spectral_summary.csv": local,
        "local_dynamics_summary.csv": local_dyn,
        "global_to_local_summary.csv": gtl,
        "global_to_local_trajectories.csv": gtl_long,
    }
    for name, frame in outputs.items():
        frame.to_csv(outdir / name, index=False)
    metadata = {
        "python": sys.version,
        "platform": platform.platform(),
        "quadrature": "deterministic Gauss-Hermite x tail-adapted chi-square survival quadrature",
        "profile": {
            "I_phi": profile.I_phi,
            "M_phi": profile.M_phi,
            "grid_size": config.get("profile_grid_size", 20001),
        },
        "rate_summary": rate_summary,
        "config": config,
    }
    (outdir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return outputs, metadata
