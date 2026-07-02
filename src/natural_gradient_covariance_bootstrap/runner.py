"""Deterministic trajectory simulation and the three deterministic experiment
drivers (covariance-bootstrap scaling, dynamic contraction, Wasserstein bootstrap).

A *trajectory* is a fixed ``(target, scheme, dt, C0)`` run of the Fisher--Rao
covariance flow, optionally preceded by one Wasserstein/Bures bootstrap step. All
expectations are exact; everything is CPU/NumPy. The full trajectory drives the
classification (``N_cov``, hitting times, contraction); the saved long rows are
decimated.
"""
from __future__ import annotations

import time

import numpy as np

from src.common.spd import symmetrize
from src.natural_gradient_covariance_bootstrap.methods import (
    mean_step, cov_step, bures_bootstrap_step, eig_extremes,
)
from src.natural_gradient_covariance_bootstrap.envelopes import (
    envelope_sequence, frozen_lower_bound, lambda_max_bound, contraction_factor,
    envelope_L0,
)
from src.natural_gradient_covariance_bootstrap import metrics as M


def _init_state(target, lambda0, init_mean_rho):
    """Initial ``(m0, C0)``: ``C0 = lambda0 I`` and a balanced far mean offset."""
    d = target.d
    m_star, C_star = target.a_star()
    s = np.sqrt(np.clip(np.diag(C_star), 0.0, None))
    m0 = np.asarray(m_star, dtype=np.float64) + float(init_mean_rho) * s
    C0 = float(lambda0) * np.eye(d, dtype=np.float64)
    return m0, C0, m_star, C_star


def deterministic_trajectory(target, scheme, dt, lambda0, n_steps, *,
                             init_mean_rho=3.0, bootstrap_c=None,
                             track_envelope=True):
    """Simulate one deterministic trajectory and return full-length arrays.

    Returns a dict with the initial state, per-iteration arrays (``iteration``,
    ``gap``, ``min_eig``, ``max_eig``, ``mean_err``, ``local_norm`` and -- when
    ``track_envelope`` -- ``L``, ``q_dynamic``, ``q_frozen``), the scalar
    constants, the bootstrap info, wall time, and status.
    """
    alpha, beta = float(target.alpha), float(target.beta)
    m0, C0, m_star, C_star = _init_state(target, lambda0, init_mean_rho)
    m, C = m0.copy(), symmetrize(C0)

    lam0_min, lam0_max, _ = eig_extremes(C0)
    lam_frozen = frozen_lower_bound(lam0_min, beta)
    lam_max_bd = lambda_max_bound(lam0_max, alpha)

    do_boot = bootstrap_c is not None
    eta = (float(bootstrap_c) / beta) if do_boot else float("nan")

    it, gap, mineig, maxeig, mean_err, local_norm = [], [], [], [], [], []
    status = "ok"
    wall = 0.0
    bootstrapped = 0
    boot_info = {"bootstrap_c": bootstrap_c, "eta": eta, "bootstrapped": 0,
                 "lambda_min_C0": lam0_min}

    def record(n, m, C):
        lmn, lmx, finite = eig_extremes(C)
        it.append(int(n))
        gap.append(float(target.energy_gap(m, C)) if finite else float("nan"))
        mineig.append(lmn)
        maxeig.append(lmx)
        mean_err.append(float(np.linalg.norm(np.asarray(m) - m_star)))
        # Fisher--Rao mean-gradient (drift) norm ||C^{1/2} G|| = sqrt(G^T C G).
        G = target.G(m, C)
        local_norm.append(float(np.sqrt(max(float(G @ (symmetrize(C) @ G)), 0.0))))
        return finite

    record(0, m, C)

    n = 0
    # One Wasserstein/Bures bootstrap step (iteration 1), only if it lifts C.
    if do_boot and lam0_min < eta:
        t0 = time.perf_counter()
        G0 = target.G(m, C)
        A0 = target.A_matrix(m, C)
        m, C = bures_bootstrap_step(m, C, G0, A0, eta)
        wall += time.perf_counter() - t0
        n = 1
        bootstrapped = 1
        boot_info["bootstrapped"] = 1
        finite = record(n, m, C)
        if not finite:
            status = "boot_nonfinite"

    # Fisher--Rao tail.
    while n < n_steps and status == "ok":
        n += 1
        t0 = time.perf_counter()
        try:
            G = target.G(m, C)
            A = target.A_matrix(m, C)
            m = mean_step(m, C, G, dt)
            C = cov_step(scheme, C, A, dt)
        except Exception as exc:  # noqa: BLE001 -- record, never crash the grid
            status = f"error:{type(exc).__name__}"
            wall += time.perf_counter() - t0
            break
        wall += time.perf_counter() - t0
        finite = record(n, m, C)
        if not finite or mineig[-1] <= 0.0:
            status = "spd_loss" if status == "ok" else status
            break

    it = np.asarray(it, dtype=np.int64)
    gap = np.asarray(gap, dtype=np.float64)
    mineig = np.asarray(mineig, dtype=np.float64)
    maxeig = np.asarray(maxeig, dtype=np.float64)
    mean_err = np.asarray(mean_err, dtype=np.float64)
    local_norm = np.asarray(local_norm, dtype=np.float64)

    out = {
        "iteration": it, "gap": gap, "min_eig": mineig, "max_eig": maxeig,
        "mean_err": mean_err, "local_norm": local_norm,
        "alpha": alpha, "beta": beta, "dt": float(dt), "lambda0": float(lambda0),
        "lambda_frozen": lam_frozen, "lambda_max_bound": lam_max_bd,
        "n_iters": int(it[-1]) if it.size else 0, "wall_time": float(wall),
        "status": status, "boot": boot_info,
        "m_star": m_star, "C_star": C_star,
    }

    if track_envelope and not do_boot:
        n_env = int(it[-1]) if it.size else 0
        L_seq = envelope_sequence(scheme, n_env, lam0_min, dt, beta)
        L = L_seq[it]                      # align envelope to recorded iterations
        q_dyn = contraction_factor(scheme, L, dt, alpha, beta, lam_max_bd)
        q_frz = contraction_factor(scheme, lam_frozen, dt, alpha, beta, lam_max_bd)
        out["L"] = np.asarray(L, dtype=np.float64)
        out["q_dynamic"] = np.asarray(q_dyn, dtype=np.float64)
        out["q_frozen"] = float(q_frz)
    return out


def _gap_ratio(gap):
    """Per-step empirical contraction ``Delta_n / Delta_{n-1}`` (NaN at n=0)."""
    ratio = np.full(gap.shape, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio[1:] = gap[1:] / gap[:-1]
    return ratio


def _base_row(experiment, target, scheme, method, traj, **extra):
    row = {
        "experiment": experiment, "target": target.name, "scheme": scheme,
        "method": method, "d": int(target.d), "kappa": float(target.kappa),
        "alpha": float(target.alpha), "beta": float(target.beta),
        "gamma": float(getattr(target, "gamma", 0.0)),
        "lambda0": float(traj["lambda0"]), "dt": float(traj["dt"]),
    }
    row.update(extra)
    return row


# ---------------------------------------------------------------------------
# Long-row emission (shared schema for results_long.csv)
# ---------------------------------------------------------------------------

LONG_COLS = [
    "experiment", "target", "scheme", "method", "d", "kappa", "alpha", "beta",
    "gamma", "lambda0", "c", "bootstrap", "dt", "iteration", "energy_gap",
    "gap_ratio", "lambda_min_C", "lambda_max_C", "L_n", "lambda_frozen",
    "q_dynamic", "q_frozen", "mean_err_norm", "local_norm", "status",
]


def _emit_long_rows(experiment, target, scheme, method, traj, c, max_saved_rows,
                    keep_iters=()):
    gap = traj["gap"]
    ratio = _gap_ratio(gap)
    it = traj["iteration"]
    has_env = "L" in traj
    keep = [int(np.where(it == k)[0][0]) for k in keep_iters if k in set(it.tolist())]
    idx = M.decimate_indices(it.size, max_saved_rows, keep=keep)
    rows = []
    for i in idx:
        row = _base_row(experiment, target, scheme, method, traj,
                        c=(float(c) if c is not None else float("nan")),
                        bootstrap=int(traj["boot"]["bootstrapped"]))
        row.update({
            "iteration": int(it[i]), "energy_gap": float(gap[i]),
            "gap_ratio": float(ratio[i]),
            "lambda_min_C": float(traj["min_eig"][i]),
            "lambda_max_C": float(traj["max_eig"][i]),
            "L_n": float(traj["L"][i]) if has_env else float("nan"),
            "lambda_frozen": float(traj["lambda_frozen"]),
            "q_dynamic": float(traj["q_dynamic"][i]) if has_env else float("nan"),
            "q_frozen": float(traj["q_frozen"]) if has_env else float("nan"),
            "mean_err_norm": float(traj["mean_err"][i]),
            "local_norm": float(traj["local_norm"][i]),
            "status": traj["status"],
        })
        rows.append(row)
    return rows


# ===========================================================================
# Experiment 1: covariance_bootstrap_scaling
# ===========================================================================

SCALING_SUMMARY_COLS = [
    "target", "scheme", "d", "kappa", "alpha", "beta", "dt", "lambda0",
    "log_inv_beta_lambda0", "inv_lambda0", "N_cov", "N_cov_envelope",
    "N_cov_theory_slope", "gap0", "gap_final", "min_eig_final", "n_iters",
    "wall_time", "status",
]


def run_scaling(cfg, progress=None):
    """Experiment 1: covariance burn-in ``N_cov`` vs ``log(1/(beta lambda0))``."""
    sc = cfg["scaling"]
    dt = float(sc["dt"])
    n_steps = int(sc["n_steps"])
    init_rho = float(cfg.get("init_mean_rho", 3.0))
    from src.natural_gradient_covariance_bootstrap.targets import build_target

    long_rows, summary_rows = [], []
    combos = [(int(d), float(k), float(l0), s)
              for d in sc["d"] for k in sc["kappa"]
              for l0 in sc["lambda0"] for s in sc["schemes"]]
    for i, (d, kappa, lam0, scheme) in enumerate(combos):
        target = build_target("gaussian", d, kappa)
        traj = deterministic_trajectory(target, scheme, dt, lam0, n_steps,
                                        init_mean_rho=init_rho, track_envelope=True)
        beta = target.beta
        ncov = M.n_cov(traj["min_eig"], beta)
        # Envelope-predicted N_cov (first envelope index reaching 1/(2 beta)).
        L_full = envelope_sequence(scheme, n_steps, min(lam0, traj["max_eig"][0]),
                                   dt, beta)
        ncov_env = M.first_crossing_ge(L_full, 1.0 / (2.0 * beta))
        ncov_iter = ncov if ncov >= 0 else traj["iteration"][
            int(np.argmax(traj["min_eig"]))]
        long_rows.extend(_emit_long_rows(
            "covariance_bootstrap_scaling", target, scheme, scheme, traj, None,
            int(cfg.get("max_saved_rows", 250)), keep_iters=(ncov_iter,)))
        summary_rows.append({
            "target": "gaussian", "scheme": scheme, "d": d, "kappa": kappa,
            "alpha": target.alpha, "beta": beta, "dt": dt, "lambda0": lam0,
            "log_inv_beta_lambda0": float(np.log(1.0 / (beta * lam0))),
            "inv_lambda0": float(1.0 / lam0),
            "N_cov": int(ncov), "N_cov_envelope": int(ncov_env),
            "N_cov_theory_slope": float(1.0 / np.log(1.0 + dt)),
            "gap0": float(traj["gap"][0]), "gap_final": float(traj["gap"][-1]),
            "min_eig_final": float(traj["min_eig"][-1]),
            "n_iters": int(traj["n_iters"]), "wall_time": float(traj["wall_time"]),
            "status": traj["status"],
        })
        if progress:
            progress("scaling", i + 1, len(combos),
                     f"d={d} kappa={kappa:g} l0={lam0:g} {scheme} N_cov={ncov}")
    return long_rows, summary_rows


# ===========================================================================
# Experiment 2: dynamic_contraction_benchmark
# ===========================================================================

CONTRACTION_COLS = [
    "target", "scheme", "d", "kappa", "alpha", "beta", "dt", "lambda0",
    "iteration", "energy_gap", "gap_ratio_emp", "q_dynamic", "q_frozen", "L_n",
    "lambda_frozen", "gap_dynamic_envelope", "gap_frozen_envelope",
    "lambda_min_C", "status",
]


def run_contraction(cfg, progress=None):
    """Experiment 2: dynamic vs frozen energy-gap contraction (theorem-safe dt)."""
    ct = cfg["contraction"]
    safety = float(ct.get("dt_safety", 0.5))
    n_steps = int(ct["n_steps"])
    init_rho = float(cfg.get("init_mean_rho", 3.0))
    from src.natural_gradient_covariance_bootstrap.targets import build_target

    long_rows, bench_rows = [], []
    combos = [(int(d), float(k), float(l0), s)
              for d in ct["d"] for k in ct["kappa"]
              for l0 in ct["lambda0"] for s in ct["schemes"]]
    for i, (d, kappa, lam0, scheme) in enumerate(combos):
        target = build_target("gaussian", d, kappa)
        alpha, beta = target.alpha, target.beta
        lam_max_bd = lambda_max_bound(lam0, alpha)
        dt = safety / (beta * lam_max_bd)          # theorem-safe stepsize
        traj = deterministic_trajectory(target, scheme, dt, lam0, n_steps,
                                        init_mean_rho=init_rho, track_envelope=True)
        it = traj["iteration"]
        gap = traj["gap"]
        ratio = _gap_ratio(gap)
        q_dyn = traj["q_dynamic"]
        q_frz = traj["q_frozen"]
        # Cumulative-product envelopes from the observed gap0.
        gap0 = float(gap[0])
        # q at step k predicts Delta_{k+1}; envelope at iteration n = gap0 * prod_{k<n} q_k.
        cum_dyn = np.concatenate([[1.0], np.cumprod(np.clip(q_dyn[:-1], 0.0, None))])
        gap_dyn_env = gap0 * cum_dyn
        gap_frz_env = gap0 * (max(q_frz, 0.0) ** it.astype(np.float64))
        long_rows.extend(_emit_long_rows(
            "dynamic_contraction_benchmark", target, scheme, scheme, traj, None,
            int(cfg.get("max_saved_rows", 250))))
        idx = M.decimate_indices(it.size, int(cfg.get("max_saved_rows", 250)))
        for j in idx:
            bench_rows.append({
                "target": "gaussian", "scheme": scheme, "d": d, "kappa": kappa,
                "alpha": alpha, "beta": beta, "dt": dt, "lambda0": lam0,
                "iteration": int(it[j]), "energy_gap": float(gap[j]),
                "gap_ratio_emp": float(ratio[j]),
                "q_dynamic": float(q_dyn[j]), "q_frozen": float(q_frz),
                "L_n": float(traj["L"][j]), "lambda_frozen": float(traj["lambda_frozen"]),
                "gap_dynamic_envelope": float(gap_dyn_env[j]),
                "gap_frozen_envelope": float(gap_frz_env[j]),
                "lambda_min_C": float(traj["min_eig"][j]), "status": traj["status"],
            })
        if progress:
            progress("contraction", i + 1, len(combos),
                     f"d={d} kappa={kappa:g} l0={lam0:g} {scheme} dt={dt:.2e}")
    return long_rows, bench_rows


# ===========================================================================
# Experiment 3: wasserstein_bootstrap_then_fr
# ===========================================================================

WBOOT_SUMMARY_COLS = [
    "target", "method", "scheme", "d", "kappa", "alpha", "beta", "gamma",
    "lambda0", "c", "eta", "dt", "bootstrapped", "lambda_min_C0",
    "lambda_min_after", "lambda_max_after", "gap0", "gap_after",
] + M.HIT_KEYS + ["gap_final", "n_iters", "wall_time", "status"]


def run_wasserstein(cfg, progress=None):
    """Experiment 3: pure FR vs one Wasserstein bootstrap then FR tail."""
    wb = cfg["wasserstein"]
    dt = float(wb["dt"])
    n_steps = int(wb["n_steps"])
    init_rho = float(cfg.get("init_mean_rho", 3.0))
    from src.natural_gradient_covariance_bootstrap.targets import build_target

    long_rows, summary_rows = [], []
    tgt_specs = wb["targets"]          # list of {name, kappa, gamma?}
    combos = []
    for spec in tgt_specs:
        for lam0 in wb["lambda0"]:
            for scheme in wb["schemes"]:
                combos.append((spec, float(lam0), scheme, None))          # pure
                for c in wb["c"]:
                    combos.append((spec, float(lam0), scheme, float(c)))  # wboot

    for i, (spec, lam0, scheme, c) in enumerate(combos):
        target = build_target(spec["name"], int(spec["d"]), float(spec["kappa"]),
                              gamma=float(spec.get("gamma", 0.0)),
                              n_nodes=int(spec.get("gh_nodes", 80)))
        method = scheme if c is None else f"wboot_{scheme}"
        traj = deterministic_trajectory(target, scheme, dt, lam0, n_steps,
                                        init_mean_rho=init_rho, bootstrap_c=c,
                                        track_envelope=False)
        it, gap, mineig, maxeig = (traj["iteration"], traj["gap"],
                                   traj["min_eig"], traj["max_eig"])
        hits = M.hitting_times(gap)
        # "after warmup" = state at iteration 1 (the bootstrap for wboot, else the
        # first FR step); for a skipped bootstrap this is the first FR step.
        after_idx = 1 if it.size > 1 else 0
        eta = (c / target.beta) if c is not None else float("nan")
        long_rows.extend(_emit_long_rows(
            "wasserstein_bootstrap_then_fr", target, scheme, method, traj, c,
            int(cfg.get("max_saved_rows", 250))))
        summ = {
            "target": target.name, "method": method, "scheme": scheme,
            "d": int(target.d), "kappa": float(target.kappa),
            "alpha": float(target.alpha), "beta": float(target.beta),
            "gamma": float(getattr(target, "gamma", 0.0)),
            "lambda0": lam0, "c": (float(c) if c is not None else float("nan")),
            "eta": float(eta), "dt": dt,
            "bootstrapped": int(traj["boot"]["bootstrapped"]),
            "lambda_min_C0": float(mineig[0]),
            "lambda_min_after": float(mineig[after_idx]),
            "lambda_max_after": float(maxeig[after_idx]),
            "gap0": float(gap[0]), "gap_after": float(gap[after_idx]),
            "gap_final": float(gap[-1]), "n_iters": int(traj["n_iters"]),
            "wall_time": float(traj["wall_time"]), "status": traj["status"],
        }
        summ.update(hits)
        summary_rows.append(summ)
        if progress:
            progress("wasserstein", i + 1, len(combos),
                     f"{target.name} l0={lam0:g} {method} boot={summ['bootstrapped']}")
    return long_rows, summary_rows
