"""Experiment 4 -- stochastic STL trajectories and the noise floor.

We run the stochastic Gaussian natural-gradient scheme with the raw and the STL
mean estimators, on the KL and (optionally) the Riemannian covariance update, and
measure the asymptotic energy-gap noise floor. The covariance update uses the
*sampled* Hessian for both raw and STL (STL changes only the mean block).

One full step (per seed) with a mini-batch of size ``B``::

    X_{n,b} = m_n + C_n^{1/2} xi_{n,b},        xi_{n,b} ~ N(0, I)
    b_raw   = -grad V(X)                        (raw)
    b_stl   = -grad V(X) + C_n^{-1}(X - m_n)    (STL)
    b_bar   = mean_b b
    S_bar   = mean_b grad^2 V(X)                (diagonal for these targets)
    m_{n+1} = m_n + dt C_n b_bar
    C_{n+1} = C^{1/2} exp(dt (I - C^{1/2} S_bar C^{1/2})) C^{1/2}   (Riemannian)
              or (1 + dt)(C_n^{-1} + dt S_bar)^{-1}                 (KL).

All seeds of a cell are advanced in one batched pass on the chosen backend
(NumPy, or a single CUDA torch device). The evaluation energy gap is the *exact*
``Delta(m_n, C_n)`` of the state (closed form / Gauss--Hermite), never an estimate
that reuses the update sample. The intrinsic Hessian fluctuation ``Psi`` is
evaluated deterministically for reference (it is exactly ``0`` for the Gaussian
target and positive for ``log cosh``).
"""
from __future__ import annotations

import hashlib
import math
import time

import numpy as np

from src.common.spd import symmetrize
from src.natural_gradient_discretization_stepsize.targets import gauss_hermite_nodes
from src.natural_gradient_stl_variance.linalg import ArrayBackend
from src.natural_gradient_covariance_bootstrap import metrics as M
from src.natural_gradient_covariance_bootstrap.targets import build_target

SPD_FLOOR = 1e-14
METHODS = ["kl_raw", "kl_stl", "riemannian_raw", "riemannian_stl"]


def _derive_seed(*parts):
    """Stable 60-bit integer seed from ``parts`` (reproducible across processes)."""
    s = "|".join(str(p) for p in parts)
    return int(hashlib.sha256(s.encode()).hexdigest()[:15], 16)


def method_scheme_stl(method):
    """Split a method name into ``(scheme, uses_stl)``."""
    mapping = {
        "kl_raw": ("kl", False), "kl_stl": ("kl", True),
        "riemannian_raw": ("riemannian", False),
        "riemannian_stl": ("riemannian", True),
    }
    if method not in mapping:
        raise ValueError(f"unknown method '{method}' (known: {METHODS})")
    return mapping[method]


def _energy_gap_batched(bk, target, m_bk, C_bk, H_bk, gh_nodes, gh_weights,
                        log_const, F_star):
    """Batched exact energy gap ``Delta(m, C)`` for all seeds (shape ``(S,)``)."""
    d = target.d
    diagC = bk.diagonal(C_bk)
    _, logdetC = bk.slogdet(C_bk)
    if target.kind == "gaussian":
        mHm = bk.sum(H_bk * (m_bk * m_bk), axis=1)
        TrHC = bk.sum(H_bk * diagC, axis=1)
        return 0.5 * (mHm + TrHC - (log_const + logdetC) - d)
    # log-cosh: objective(m,C) - F_star, objective = -0.5 logdet C + E_q[V].
    quad = 0.5 * bk.sum(H_bk * (m_bk * m_bk + diagC), axis=1)
    sd = bk.sqrt(bk.clip_min(diagC, 0.0))
    y = m_bk[..., None] + sd[..., None] * gh_nodes        # (S, d, n_nodes)
    ay = bk.abs(y)
    lc = ay + bk.log1p(bk.exp(-2.0 * ay)) - np.log(2.0)
    e_logcosh = bk.sum(lc * gh_weights, axis=-1)          # (S, d)
    eV = quad + target.gamma * bk.sum(e_logcosh, axis=1)
    return (-0.5 * logdetC + eV) - F_star


def simulate_stl_cell(target, method, dt, n_steps, batch_size, seeds, backend,
                      cell_seed, *, tail_frac=0.25, init_mean_rho=2.0,
                      init_cov_scale=0.5, gh_nodes=40):
    """Simulate every seed of one ``(target, method, dt, batch_size)`` cell.

    Returns ``(summary_rows, diag)``: one summary dict per seed (tail floor,
    status) and a cell-level diagnostics dict (wall time, SPD clip count, the
    deterministic ``Psi`` references, and the seed-median tail floor).
    """
    bk = backend
    scheme, use_stl = method_scheme_stl(method)
    S = len(seeds)
    d = target.d
    need_inv = use_stl or scheme == "kl"

    m_star, C_star = target.a_star()
    F_star = float(target.objective(m_star, C_star))
    H_bk = bk.asarray(target.H_diag)
    gamma = float(getattr(target, "gamma", 0.0))
    log_const = float(np.sum(np.log(target.H_diag)))       # log det H

    s = np.sqrt(np.clip(np.diag(C_star), 0.0, None))
    m0 = np.asarray(m_star, dtype=np.float64) + init_mean_rho * s
    C0 = symmetrize(init_cov_scale * np.asarray(C_star, dtype=np.float64))
    m_bk = bk.asarray(np.broadcast_to(m0, (S, d)).copy())
    C_bk = bk.asarray(np.broadcast_to(C0, (S, d, d)).copy())

    nodes_np, weights_np = gauss_hermite_nodes(gh_nodes)
    gh_nodes_bk = bk.asarray(nodes_np)
    gh_weights_bk = bk.asarray(weights_np)
    eye = bk.eye(d, batch=(S,))

    gen = bk.generator(cell_seed)
    gaps_full = np.empty((n_steps + 1, S), dtype=np.float64)
    spd_fail = np.zeros(S, dtype=bool)
    n_clips = 0

    def record_gap(n, m_bk, C_bk):
        nonlocal n_clips
        gap = bk.to_numpy(_energy_gap_batched(
            bk, target, m_bk, C_bk, H_bk, gh_nodes_bk, gh_weights_bk,
            log_const, F_star))
        w = bk.to_numpy(bk.eigvalsh(C_bk))
        min_eig = w[:, 0]
        fail = (~np.isfinite(gap)) | (min_eig <= SPD_FLOOR)
        spd_fail[fail] = True
        gaps_full[n] = gap

    record_gap(0, m_bk, C_bk)

    t0 = time.perf_counter()
    for n in range(1, n_steps + 1):
        Z = bk.randn((S, batch_size, d), gen)
        C_sqrt = bk.sqrtm_sym(C_bk)                        # (S, d, d)
        C_inv = bk.inv(C_bk) if need_inv else None
        theta = m_bk[:, None, :] + Z @ C_sqrt              # (S, B, d)

        score = -(theta * H_bk)
        if gamma != 0.0:
            score = score - gamma * bk.tanh(theta)
        if use_stl:
            b = score + (theta - m_bk[:, None, :]) @ C_inv
        else:
            b = score
        b_bar = bk.mean(b, axis=1)                         # (S, d)

        # Sampled Hessian (curvature): H + gamma sech^2(theta).
        if gamma != 0.0:
            t = bk.tanh(theta)
            hd = H_bk + gamma * (1.0 - t * t)
        else:
            hd = bk.asarray(np.broadcast_to(target.H_diag, theta.shape).copy())
        S_bar = bk.diag_embed(bk.mean(hd, axis=1))         # (S, d, d)

        m_bk = m_bk + dt * bk.matvec(C_bk, b_bar)
        if scheme == "riemannian":
            inner = eye - bk.sym(C_sqrt @ S_bar @ C_sqrt)
            expA = bk.expm_sym(dt * inner)
            C_bk = bk.sym(C_sqrt @ expA @ C_sqrt)
        else:  # kl
            P = bk.sym(C_inv + dt * S_bar)
            C_bk = bk.sym((1.0 + dt) * bk.inv(P))

        w = bk.eigvalsh(C_bk)
        wmin = bk.to_numpy(w)[:, 0]
        if np.any(wmin < SPD_FLOOR):
            n_clips += int(np.sum(wmin < SPD_FLOOR))
            ew, eV = bk.eigh(C_bk)
            ew = bk.clip_min(ew, SPD_FLOOR)
            C_bk = bk.sym((eV * ew[..., None, :]) @ bk.transpose(eV))

        record_gap(n, m_bk, C_bk)
    wall = time.perf_counter() - t0

    # Deterministic Psi references (exact; 0 for the Gaussian target).
    psi_star = float(target.Psi(m_star, C_star))
    m_final_mean = bk.to_numpy(m_bk).mean(axis=0)
    C_final_mean = symmetrize(bk.to_numpy(C_bk).mean(axis=0))
    psi_tail = float(target.Psi(m_final_mean, C_final_mean))

    summary_rows = []
    for si, seed in enumerate(seeds):
        tail = M.tail_floor(gaps_full[:, si], tail_frac=tail_frac)
        row = {
            "target": target.name, "kind": target.kind, "method": method,
            "scheme": scheme, "stl": int(use_stl), "seed": int(seed),
            "d": int(d), "kappa": float(target.kappa), "gamma": gamma,
            "dt": float(dt), "batch_size": int(batch_size), "n_steps": int(n_steps),
            "spd_fail": int(bool(spd_fail[si])),
        }
        row.update(tail)
        summary_rows.append(row)

    diag = {
        "method": method, "scheme": scheme, "stl": int(use_stl), "dt": float(dt),
        "wall_time_cell": float(wall), "n_clips": int(n_clips), "n_seeds": int(S),
        "psi_star": psi_star, "psi_tail": psi_tail,
        "tail_median_floor": float(np.median(
            [r["tail_median_gap"] for r in summary_rows])),
        "final_gap_median": float(np.median([r["final_gap"] for r in summary_rows])),
    }
    return summary_rows, diag


# ===========================================================================
# Experiment 4 driver: stl_noise_floor
# ===========================================================================

STL_FLOOR_COLS = [
    "target", "kind", "method", "scheme", "stl", "d", "kappa", "gamma", "dt",
    "batch_size", "n_seeds", "n_steps", "tail_median_gap", "tail_mean_gap",
    "tail_median_gap_std", "final_gap", "psi_star", "psi_tail", "dt_psi_star",
    "floor_over_dt_psi", "n_clips", "spd_fail_any", "wall_time_cell",
]


def _n_steps_for_dt(dt, horizon, min_steps, max_steps):
    """Number of steps so the flow horizon ``dt * n_steps`` is ~ ``horizon``."""
    return int(min(max_steps, max(min_steps, math.ceil(horizon / dt))))


def run_stl_floor(cfg, device="cpu", backend_name="numpy", progress=None):
    """Experiment 4: stochastic STL noise floor vs ``dt`` for both targets.

    Returns ``(floor_rows, diag_rows)``: one aggregated (over seeds) row per
    ``(target, method, dt)`` cell for ``stl_floor_summary.csv``, and the cell
    diagnostics.
    """
    sf = cfg["stl_floor"]
    base_seed = int(cfg.get("base_seed", 0))
    n_seeds = int(sf["n_seeds"])
    seeds = list(range(n_seeds))
    batch_size = int(sf.get("batch_size", 1))
    tail_frac = float(sf.get("tail_frac", 0.25))
    horizon = float(sf.get("horizon", 30.0))
    min_steps = int(sf.get("min_steps", 400))
    max_steps = int(sf.get("max_steps", 8000))
    methods = list(sf.get("methods", METHODS))
    gh_nodes = int(sf.get("gh_nodes", 40))

    bk = ArrayBackend(backend=backend_name, device=device, dtype="float64")
    target_cache = {}
    floor_rows, diag_rows = [], []
    cells = [(spec, float(dt), method)
             for spec in sf["targets"] for dt in sf["dt_list"]
             for method in methods]
    for i, (spec, dt, method) in enumerate(cells):
        key = (spec["name"], int(spec["d"]), float(spec["kappa"]),
               float(spec.get("gamma", 0.0)))
        target = target_cache.get(key)
        if target is None:
            target = build_target(spec["name"], int(spec["d"]), float(spec["kappa"]),
                                  gamma=float(spec.get("gamma", 0.0)),
                                  n_nodes=int(spec.get("gh_nodes", 80)))
            target_cache[key] = target
        n_steps = _n_steps_for_dt(dt, horizon, min_steps, max_steps)
        cell_seed = _derive_seed(base_seed, spec["name"], spec["d"], spec["kappa"],
                                 spec.get("gamma", 0.0), method, dt, batch_size)
        summ, diag = simulate_stl_cell(
            target, method, dt, n_steps, batch_size, seeds, bk, cell_seed,
            tail_frac=tail_frac, init_mean_rho=float(sf.get("init_mean_rho", 2.0)),
            init_cov_scale=float(sf.get("init_cov_scale", 0.5)), gh_nodes=gh_nodes)

        tail_med = np.array([r["tail_median_gap"] for r in summ], dtype=np.float64)
        tail_mean = np.array([r["tail_mean_gap"] for r in summ], dtype=np.float64)
        final = np.array([r["final_gap"] for r in summ], dtype=np.float64)
        floor = float(np.median(tail_med[np.isfinite(tail_med)])) \
            if np.any(np.isfinite(tail_med)) else float("nan")
        dt_psi = dt * diag["psi_star"]
        floor_rows.append({
            "target": target.name, "kind": target.kind, "method": method,
            "scheme": diag["scheme"], "stl": diag["stl"], "d": int(target.d),
            "kappa": float(target.kappa), "gamma": float(getattr(target, "gamma", 0.0)),
            "dt": dt, "batch_size": batch_size, "n_seeds": n_seeds,
            "n_steps": n_steps, "tail_median_gap": floor,
            "tail_mean_gap": float(np.median(tail_mean[np.isfinite(tail_mean)]))
            if np.any(np.isfinite(tail_mean)) else float("nan"),
            "tail_median_gap_std": float(np.std(tail_med[np.isfinite(tail_med)]))
            if np.any(np.isfinite(tail_med)) else float("nan"),
            "final_gap": float(np.median(final[np.isfinite(final)]))
            if np.any(np.isfinite(final)) else float("nan"),
            "psi_star": diag["psi_star"], "psi_tail": diag["psi_tail"],
            "dt_psi_star": dt_psi,
            "floor_over_dt_psi": (floor / dt_psi) if dt_psi > 0 else float("nan"),
            "n_clips": diag["n_clips"],
            "spd_fail_any": int(any(r["spd_fail"] for r in summ)),
            "wall_time_cell": diag["wall_time_cell"],
        })
        diag_rows.append({"target": target.name, "dt": dt, **diag})
        if progress:
            progress("stl_floor", i + 1, len(cells),
                     f"{target.name} {method} dt={dt:.4g} floor={floor:.2e} "
                     f"({diag['wall_time_cell']:.1f}s)")
    return floor_rows, diag_rows
