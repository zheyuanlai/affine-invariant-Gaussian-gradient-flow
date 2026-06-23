"""Experiment B -- algorithm-level stochastic natural-gradient trajectories.

We run the stochastic Gaussian natural-gradient scheme with and without the
mean-block STL estimator. Four methods are compared:

* ``riemannian``      -- baseline mean ``b_base``, Riemannian covariance update;
* ``riemannian_stl``  -- STL mean ``b_stl``,  Riemannian covariance update;
* ``kl``              -- baseline mean ``b_base``, KL covariance update;
* ``kl_stl``          -- STL mean ``b_stl``,  KL covariance update.

Within a scheme the covariance update uses the *same* sampled Hessian ``S_n`` for
baseline and STL, so the experiment isolates the mean estimator.

One full step (per seed) with a mini-batch of size ``B``::

    theta_{n,b} = m_n + C_n^{1/2} z_{n,b},        z_{n,b} ~ N(0, I)
    b_bar       = mean_b [ score_post(theta) (+ C_n^{-1}(theta - m_n) if STL) ]
    S_bar       = mean_b Hess_log_post(theta)     (diagonal for these targets)
    m_{n+1}     = m_n + dt C_n b_bar
    C_{n+1}     = riemannian / KL covariance update with S_bar.

All seeds in a cell are advanced in a single batched pass (leading axis ``S``)
on the chosen backend, and the four methods sharing a ``(config, batch_size)``
cell are driven from the *same* common-random-number stream (the generator is
re-seeded identically per method), so the baseline/STL comparison is paired.
"""
from __future__ import annotations

import time

import numpy as np

from src.common.spd import symmetrize
from src.natural_gradient_discretization_stepsize.targets import gauss_hermite_nodes
from src.natural_gradient_stl_variance.estimators import method_scheme_stl
from src.natural_gradient_stl_variance.metrics import (
    hitting_times, tail_noise_floor, GAP_THRESHOLDS, _tol_key,
)

SPD_FLOOR = 1e-12          # emergency eigenvalue floor (recorded, never silent)


def _save_steps(n_steps, n_saved):
    """Decimated step indices to persist (always includes 0 and ``n_steps``)."""
    n_saved = min(int(n_saved), n_steps + 1)
    idx = np.unique(np.linspace(0, n_steps, num=n_saved).round().astype(int))
    return set(int(i) for i in idx)


def _init_state(bk, target, S, m_star, C_star, init_mean_rho, init_cov_scale):
    """Batched initial ``(m, C)`` shared by all seeds (they diverge via noise)."""
    d = target.d
    s = np.sqrt(np.clip(np.diag(C_star), 0.0, None))
    m0 = np.asarray(m_star, dtype=np.float64) + init_mean_rho * s
    C0 = symmetrize(init_cov_scale * np.asarray(C_star, dtype=np.float64))
    m_bk = bk.asarray(np.broadcast_to(m0, (S, d)).copy())
    C_bk = bk.asarray(np.broadcast_to(C0, (S, d, d)).copy())
    return m_bk, C_bk, m0, C0


def _energy_gap_batched(bk, target, m_bk, C_bk, a_bk, gh_nodes, gh_weights,
                        log_const, F_star, m_star_bk, C_star_bk):
    """Batched energy gap ``E(a) - E(a_star)`` for all seeds (shape ``(S,)``)."""
    d = target.d
    diagC = bk.diagonal(C_bk)                       # (S, d)
    _, logdetC = bk.slogdet(C_bk)                   # (S,)
    if target.kind == "gaussian":
        mAm = bk.sum(a_bk * (m_bk * m_bk), axis=1)
        TrAC = bk.sum(a_bk * diagC, axis=1)
        gap = 0.5 * (mAm + TrAC - (log_const + logdetC) - d)
        return gap
    # log-cosh: objective(m,C) - F_star, objective = -0.5 logdet C + E_q[V].
    quad = 0.5 * bk.sum(a_bk * (m_bk * m_bk + diagC), axis=1)
    # E[log cosh(theta_i)] for theta_i ~ N(m_i, C_ii) via Gauss--Hermite.
    sd = bk.sqrt(bk.clip_min(diagC, 0.0))           # (S, d)
    y = m_bk[..., None] + sd[..., None] * gh_nodes  # (S, d, n_nodes)
    ay = bk.abs(y)
    lc = ay + bk.log1p(bk.exp(-2.0 * ay)) - np.log(2.0)
    e_logcosh = bk.sum(lc * gh_weights, axis=-1)    # (S, d)
    eV = quad + target.tau * bk.sum(e_logcosh, axis=1)
    objective = -0.5 * logdetC + eV
    return objective - F_star


def _w2sq_batched(bk, m_bk, C_bk, m_star_bk, Cs_sqrt_bk, tr_Cstar):
    """Batched squared 2-Wasserstein distance to ``a_star`` (shape ``(S,)``)."""
    inner = bk.sym(Cs_sqrt_bk @ C_bk @ Cs_sqrt_bk)
    cross = bk.sqrtm_sym(inner)
    bures = bk.trace(C_bk) + tr_Cstar - 2.0 * bk.trace(cross)
    dm = m_bk - m_star_bk
    return bk.sum(dm * dm, axis=1) + bk.clip_min(bures, 0.0)


def simulate_cell(target, method, dt, n_steps, batch_size, seeds, backend,
                  cell_seed, *, n_saved=160, tail_frac=0.2,
                  init_mean_rho=3.0, init_cov_scale=0.5, gh_nodes=40):
    """Simulate every seed of one ``(target, method, batch_size)`` cell (batched).

    Returns ``(long_rows, summary_rows, tail_rows, diag)``:

    * ``long_rows``    -- one dict per saved step per seed (the trajectory CSV);
    * ``summary_rows`` -- one dict per seed (final metrics, hitting times, status);
    * ``tail_rows``    -- one dict per seed (tail / noise-floor statistics);
    * ``diag``         -- cell-level diagnostics (wall time, SPD clip count).
    """
    bk = backend
    scheme, use_stl = method_scheme_stl(method)
    S = len(seeds)
    d = target.d
    need_inv = use_stl or scheme == "kl"

    m_star, C_star = target.a_star()
    F_star = float(target.objective(m_star, C_star))
    a_bk = bk.asarray(target.a_diag)
    tau = float(target.tau)
    log_const = float(np.sum(np.log(target.a_diag)))   # log det A
    m_bk, C_bk, m0, C0 = _init_state(bk, target, S, m_star, C_star,
                                     init_mean_rho, init_cov_scale)

    # Constants for batched metrics.
    m_star_bk = bk.asarray(np.broadcast_to(m_star, (S, d)).copy())
    C_star_bk = bk.asarray(np.broadcast_to(C_star, (S, d, d)).copy())
    Cs_sqrt_bk = bk.sqrtm_sym(bk.asarray(C_star))      # (d, d)
    tr_Cstar = float(np.trace(C_star))
    cstar_fro = float(np.linalg.norm(symmetrize(C_star), "fro"))
    nodes_np, weights_np = gauss_hermite_nodes(gh_nodes)
    gh_nodes_bk = bk.asarray(nodes_np)
    gh_weights_bk = bk.asarray(weights_np)

    save_at = _save_steps(n_steps, n_saved)
    gen = bk.generator(cell_seed)

    gaps_full = np.empty((n_steps + 1, S), dtype=np.float64)
    mineig_full = np.empty((n_steps + 1, S), dtype=np.float64)
    spd_fail = np.zeros(S, dtype=bool)
    n_clips = 0
    long_rows = []

    def _record_step(n, m_bk, C_bk):
        """Compute batched metrics at step ``n`` and append saved long rows."""
        nonlocal n_clips
        w = bk.eigvalsh(C_bk)                         # (S, d), ascending
        w_np = bk.to_numpy(w)
        min_eig = w_np[:, 0]
        max_eig = w_np[:, -1]
        gap = bk.to_numpy(_energy_gap_batched(
            bk, target, m_bk, C_bk, a_bk, gh_nodes_bk, gh_weights_bk,
            log_const, F_star, m_star_bk, C_star_bk))
        finite = np.isfinite(gap) & np.isfinite(min_eig)
        fail = (~finite) | (min_eig <= SPD_FLOOR)
        spd_fail[fail] = True
        gaps_full[n] = gap
        mineig_full[n] = min_eig
        if n in save_at:
            m_np = bk.to_numpy(m_bk)
            dm = m_np - np.asarray(m_star)
            sq_mean_err = np.sum(dm * dm, axis=1)
            cov_fro = bk.to_numpy(bk.norm_fro(C_bk - C_star_bk))
            rel_cov_fro = cov_fro / cstar_fro if cstar_fro > 0 else np.full(S, np.nan)
            w2sq = bk.to_numpy(_w2sq_batched(
                bk, m_bk, C_bk, m_star_bk, Cs_sqrt_bk, tr_Cstar))
            for si, seed in enumerate(seeds):
                long_rows.append({
                    "target_name": target.name, "kind": target.kind,
                    "method": method, "scheme": scheme,
                    "stl": int(use_stl), "seed": int(seed),
                    "d": int(d), "kappa": float(target.kappa), "tau": tau,
                    "dt": float(dt), "batch_size": int(batch_size),
                    "step": int(n),
                    "energy_gap": float(gap[si]),
                    "sq_mean_error": float(sq_mean_err[si]),
                    "rel_cov_fro_error": float(rel_cov_fro[si]),
                    "w2_sq": float(w2sq[si]),
                    "min_eig_C": float(min_eig[si]),
                    "max_eig_C": float(max_eig[si]),
                    "spd_fail": int(bool(fail[si])),
                })

    _record_step(0, m_bk, C_bk)

    t0 = time.perf_counter()
    for n in range(1, n_steps + 1):
        Z = bk.randn((S, batch_size, d), gen)
        C_sqrt = bk.sqrtm_sym(C_bk)                   # (S, d, d)
        C_inv = bk.inv(C_bk) if need_inv else None
        theta = m_bk[:, None, :] + Z @ C_sqrt         # (S, B, d)

        score = -(theta * a_bk)
        if tau != 0.0:
            score = score - tau * bk.tanh(theta)
        if use_stl:
            corr = (theta - m_bk[:, None, :]) @ C_inv
            b = score + corr
        else:
            b = score
        b_bar = bk.mean(b, axis=1)                    # (S, d)

        # Sampled Hessian (diagonal): -(a + tau sech^2(theta)).
        if tau != 0.0:
            t = bk.tanh(theta)
            hd = -(a_bk + tau * (1.0 - t * t))
        else:
            hd = -bk.asarray(np.broadcast_to(target.a_diag, theta.shape).copy())
        S_diag_bar = bk.mean(hd, axis=1)              # (S, d)
        S_bar = bk.diag_embed(S_diag_bar)             # (S, d, d)

        m_bk = m_bk + dt * bk.matvec(C_bk, b_bar)
        if scheme == "riemannian":
            inner = bk.sym(C_sqrt @ S_bar @ C_sqrt)
            expA = bk.expm_sym(dt * inner)
            edt = float(np.exp(dt))   # python float multiplies both backends
            C_bk = bk.sym(edt * (C_sqrt @ expA @ C_sqrt))
        else:  # kl
            P = bk.sym(C_inv - dt * S_bar)
            C_bk = bk.sym((1.0 + dt) * bk.inv(P))

        # Emergency SPD guard (recorded, not hidden).
        w = bk.eigvalsh(C_bk)
        wmin = bk.to_numpy(w)[:, 0]
        if np.any(wmin < SPD_FLOOR):
            n_clips += int(np.sum(wmin < SPD_FLOOR))
            ew, eV = bk.eigh(C_bk)
            ew = bk.clip_min(ew, SPD_FLOOR)
            C_bk = bk.sym((eV * ew[..., None, :]) @ bk.transpose(eV))

        _record_step(n, m_bk, C_bk)
    wall = time.perf_counter() - t0

    summary_rows, tail_rows = [], []
    for si, seed in enumerate(seeds):
        gseq = gaps_full[:, si]
        hits = hitting_times(gseq)
        tail = tail_noise_floor(gseq, tail_frac=tail_frac)
        common = {
            "target_name": target.name, "kind": target.kind,
            "method": method, "scheme": scheme, "stl": int(use_stl),
            "seed": int(seed), "d": int(d), "kappa": float(target.kappa),
            "tau": tau, "dt": float(dt), "batch_size": int(batch_size),
            "n_steps": int(n_steps),
        }
        summ = dict(common)
        summ.update({
            "gap0": float(gseq[0]), "gap_final": float(gseq[-1]),
            "gap_min": float(np.nanmin(gseq[np.isfinite(gseq)]))
            if np.any(np.isfinite(gseq)) else float("nan"),
            "min_eig_C_min": float(np.nanmin(mineig_full[:, si])),
            "spd_fail": int(bool(spd_fail[si])),
            "wall_time_cell": float(wall), "n_seeds_in_cell": int(S),
        })
        summ.update(hits)
        summ.update({k: tail[k] for k in
                     ("tail_mean_gap", "tail_median_gap", "tail_std_gap",
                      "final_gap", "tail_steps")})
        summary_rows.append(summ)

        trow = dict(common)
        trow.update(tail)
        trow["spd_fail"] = int(bool(spd_fail[si]))
        tail_rows.append(trow)

    diag = {
        "method": method, "wall_time_cell": float(wall),
        "n_clips": int(n_clips), "n_seeds": int(S),
        "m0": np.asarray(m0).tolist(), "C0_diag": np.diag(C0).tolist(),
    }
    return long_rows, summary_rows, tail_rows, diag


HIT_KEYS = [f"iter_to_{_tol_key(t)}" for t in GAP_THRESHOLDS]
