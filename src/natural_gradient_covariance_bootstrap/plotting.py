"""Figure builders for the covariance-bootstrap report.

Each builder takes the already-loaded final DataFrame(s) and returns a matplotlib
figure; they are the single source of truth shared by the per-group
``plot_results.py`` and by ``reports/make_report_assets.py``. Matplotlib only
(no seaborn); the shared rcParams are applied by the caller.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

GAP_FLOOR = 1e-16          # display floor for log-scale gaps
FLOOR_DISP = 1e-14         # display floor for the STL noise floor

SCHEME_STYLE = {
    "riemannian": {"color": "#1f77b4", "label": "Riemannian"},
    "kl": {"color": "#d62728", "label": "KL"},
}
TARGET_TEX = {"gaussian": "Gaussian", "log_cosh": "log-cosh"}
STL_METHOD_STYLE = {
    "kl_raw": {"color": "#d62728", "ls": "-", "label": "KL raw"},
    "kl_stl": {"color": "#d62728", "ls": "--", "label": "KL STL"},
    "riemannian_raw": {"color": "#1f77b4", "ls": "-", "label": "Riem. raw"},
    "riemannian_stl": {"color": "#1f77b4", "ls": "--", "label": "Riem. STL"},
}


def _clip(y, floor=GAP_FLOOR):
    return np.clip(np.asarray(y, dtype=np.float64), floor, None)


def _match(series, value, rel=1e-6):
    """Relative-tolerance mask (``np.isclose``'s ``atol=1e-8`` conflates tiny
    ``lambda0`` values like ``1e-10`` and ``1e-8``)."""
    return np.abs(np.asarray(series, dtype=np.float64) - value) <= rel * abs(value)


# ---------------------------------------------------------------------------
# fig_covariance_envelope: lambda_min(C_n) vs theoretical L_n
# ---------------------------------------------------------------------------

def fig_covariance_envelope(long_df, scheme="kl", kappa=None):
    """Actual ``lambda_min(C_n)`` (solid) vs envelope ``L_n`` (dashed) by ``lambda0``."""
    df = long_df[(long_df.experiment == "covariance_bootstrap_scaling")
                 & (long_df.scheme == scheme)].copy()
    if kappa is None:
        kappa = sorted(df.kappa.unique())[len(sorted(df.kappa.unique())) // 2]
    df = df[np.isclose(df.kappa, kappa)]
    d_val = sorted(df.d.unique())[-1]
    df = df[df.d == d_val]
    beta = float(df.beta.iloc[0])
    lam0s = sorted(df.lambda0.unique())
    cmap = plt.cm.viridis(np.linspace(0.1, 0.88, len(lam0s)))
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    for k, lam0 in enumerate(lam0s):
        s = df[_match(df.lambda0, lam0)].sort_values("iteration")
        if s.empty:
            continue
        ax.plot(s.iteration, _clip(s.lambda_min_C), color=cmap[k], lw=1.7,
                label=rf"$\lambda_0=10^{{{int(round(np.log10(lam0)))}}}$")
        ax.plot(s.iteration, _clip(s.L_n), color=cmap[k], lw=1.2, ls="--", alpha=0.9)
    ax.axhline(1.0 / (2.0 * beta), color="k", ls=":", lw=1.1,
               label=r"$1/(2\beta)$")
    ax.set_yscale("log")
    ax.set_xlabel("iteration $n$")
    ax.set_ylabel(r"$\lambda_{\min}(C_n)$  (solid),  $L_n$ (dashed)")
    ax.set_title(rf"{SCHEME_STYLE[scheme]['label']} covariance envelope "
                 rf"($\kappa={kappa:g}$, $d={int(d_val)}$)")
    ax.legend(fontsize=7.5, loc="lower right", ncol=2)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# fig_warmup_scaling: N_cov vs log(1/(beta lambda0)) and vs 1/lambda0
# ---------------------------------------------------------------------------

def _fit_line(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 2:
        return None
    a, b = np.polyfit(x[ok], y[ok], 1)
    yhat = a * x[ok] + b
    ss_res = float(np.sum((y[ok] - yhat) ** 2))
    ss_tot = float(np.sum((y[ok] - y[ok].mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(a), float(b), float(r2)


def fig_warmup_scaling(summary_df):
    """Left: ``N_cov`` vs ``log(1/(beta lambda0))`` (linear). Right: vs ``1/lambda0``."""
    df = summary_df[summary_df.N_cov >= 0].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    schemes = sorted(df.scheme.unique())
    markers = {"kl": "o", "riemannian": "s"}
    kappas = sorted(df.kappa.unique())
    cmap = plt.cm.plasma(np.linspace(0.1, 0.78, len(kappas)))
    ax = axes[0]
    for scheme in schemes:
        for ki, kappa in enumerate(kappas):
            s = df[(df.scheme == scheme) & np.isclose(df.kappa, kappa)]
            if s.empty:
                continue
            ax.scatter(s.log_inv_beta_lambda0, s.N_cov, color=cmap[ki],
                       marker=markers.get(scheme, "o"), s=42,
                       edgecolor="white", linewidth=0.5,
                       label=rf"{scheme}, $\kappa={kappa:g}$")
    # fit and theory slope
    fit = _fit_line(df.log_inv_beta_lambda0, df.N_cov)
    xs = np.linspace(0, float(df.log_inv_beta_lambda0.max()) * 1.05, 50)
    if fit:
        a, b, r2 = fit
        ax.plot(xs, a * xs + b, "k-", lw=1.4,
                label=rf"fit slope $={a:.2f}$ ($R^2={r2:.3f}$)")
    slope_th = float(df.N_cov_theory_slope.iloc[0])
    ax.plot(xs, slope_th * xs, "k--", lw=1.2,
            label=rf"theory $1/\log(1+\Delta t)={slope_th:.2f}$")
    ax.set_xlabel(r"$\log(1/(\beta\lambda_0))$")
    ax.set_ylabel(r"burn-in $N_{\mathrm{cov}}$")
    ax.set_title("logarithmic burn-in")
    ax.legend(fontsize=7, loc="upper left", ncol=1)

    ax = axes[1]
    for scheme in schemes:
        for ki, kappa in enumerate(kappas):
            s = df[(df.scheme == scheme) & np.isclose(df.kappa, kappa)]
            if s.empty:
                continue
            ax.scatter(s.inv_lambda0, s.N_cov, color=cmap[ki],
                       marker=markers.get(scheme, "o"), s=42,
                       edgecolor="white", linewidth=0.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"$1/\lambda_0$")
    ax.set_ylabel(r"burn-in $N_{\mathrm{cov}}$")
    ax.set_title(r"$N_{\mathrm{cov}}$ is flat in $1/\lambda_0$ (log $x$): "
                 r"not a $1/\lambda_0$ law")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# fig_dynamic_contraction: observed gap, dynamic envelope, frozen envelope
# ---------------------------------------------------------------------------

def fig_dynamic_contraction(bench_df):
    """Observed gap vs dynamic and frozen envelopes; per-panel (kappa, scheme)."""
    keys = (bench_df[["kappa", "scheme"]].drop_duplicates()
            .sort_values(["kappa", "scheme"]).itertuples(index=False))
    keys = list(keys)
    ncol = 2
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.4 * ncol, 3.6 * nrow),
                             squeeze=False)
    for idx, (kappa, scheme) in enumerate(keys):
        ax = axes[idx // ncol][idx % ncol]
        s = bench_df[(np.isclose(bench_df.kappa, kappa))
                     & (bench_df.scheme == scheme)].sort_values("iteration")
        if s.empty:
            continue
        ax.semilogy(s.iteration, _clip(s.energy_gap), color="k", lw=1.9,
                    label="observed gap")
        ax.semilogy(s.iteration, _clip(s.gap_dynamic_envelope),
                    color="#2ca02c", ls="--", lw=1.6, label="dynamic envelope")
        ax.semilogy(s.iteration, _clip(s.gap_frozen_envelope),
                    color="#9467bd", ls=":", lw=1.6, label="frozen envelope")
        ax.set_title(rf"{SCHEME_STYLE[scheme]['label']}, $\kappa={kappa:g}$")
        ax.grid(True, which="both", alpha=0.25)
        if idx // ncol == nrow - 1:
            ax.set_xlabel("iteration $n$")
        if idx % ncol == 0:
            ax.set_ylabel("energy gap")
    for idx in range(len(keys), nrow * ncol):
        axes[idx // ncol][idx % ncol].axis("off")
    axes[0][0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Dynamic (growing-$L_n$) vs frozen contraction envelope", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


def fig_contraction_factors(bench_df):
    """Empirical per-step contraction vs dynamic and frozen predictions."""
    keys = list((bench_df[["kappa", "scheme"]].drop_duplicates()
                 .sort_values(["kappa", "scheme"]).itertuples(index=False)))
    ncol = 2
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.4 * ncol, 3.4 * nrow),
                             squeeze=False)
    for idx, (kappa, scheme) in enumerate(keys):
        ax = axes[idx // ncol][idx % ncol]
        s = bench_df[(np.isclose(bench_df.kappa, kappa))
                     & (bench_df.scheme == scheme)].sort_values("iteration")
        s = s[s.iteration > 0]
        if s.empty:
            continue
        ax.plot(s.iteration, s.gap_ratio_emp, color="k", lw=1.7, label="observed")
        ax.plot(s.iteration, s.q_dynamic, color="#2ca02c", ls="--", lw=1.5,
                label=r"$q_{\mathrm{dyn}}(L_n)$")
        ax.plot(s.iteration, s.q_frozen, color="#9467bd", ls=":", lw=1.5,
                label=r"$q_{\mathrm{frozen}}$")
        ax.set_title(rf"{SCHEME_STYLE[scheme]['label']}, $\kappa={kappa:g}$")
        ax.set_ylim(top=1.001)
        ax.grid(True, which="both", alpha=0.25)
        if idx // ncol == nrow - 1:
            ax.set_xlabel("iteration $n$")
        if idx % ncol == 0:
            ax.set_ylabel(r"$\Delta_{n+1}/\Delta_n$")
    for idx in range(len(keys), nrow * ncol):
        axes[idx // ncol][idx % ncol].axis("off")
    axes[0][0].legend(fontsize=8, loc="lower right")
    fig.suptitle("Per-step contraction: observed vs dynamic/frozen theory", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


# ---------------------------------------------------------------------------
# fig_wasserstein_bootstrap: hitting times and covariance lift
# ---------------------------------------------------------------------------

def _wboot_int(x):
    x = np.asarray(x, dtype=np.float64)
    return np.where(x < 0, np.nan, x)


def fig_wasserstein_bootstrap(wboot_df, c=0.5):
    """Pure FR vs W-boot: warm-up hitting iterations and post-warmup covariance."""
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.4), squeeze=False)
    targets = [t for t in ["gaussian", "log_cosh"] if t in set(wboot_df.target.unique())]
    schemes = sorted(wboot_df.scheme.unique())
    marker = {"kl": "o", "riemannian": "s"}
    for col, target in enumerate(targets):
        sub = wboot_df[wboot_df.target == target]
        beta = float(sub.beta.iloc[0])
        ax_hit = axes[0][col]
        ax_cov = axes[1][col]
        for scheme in schemes:
            pure = sub[(sub.scheme == scheme) & (sub.method == scheme)].sort_values("lambda0")
            wb = sub[(sub.scheme == scheme) & (sub.method == f"wboot_{scheme}")
                     & np.isclose(sub.c, c)].sort_values("lambda0")
            col0 = SCHEME_STYLE[scheme]["color"]
            if not pure.empty:
                ax_hit.plot(pure.lambda0, _wboot_int(pure.iter_to_1e_minus_1),
                            color=col0, ls="-", marker=marker[scheme],
                            label=f"pure {scheme}")
                ax_cov.plot(pure.lambda0, pure.lambda_min_after, color=col0, ls="-",
                            marker=marker[scheme], label=f"pure {scheme}")
            if not wb.empty:
                ax_hit.plot(wb.lambda0, _wboot_int(wb.iter_to_1e_minus_1),
                            color=col0, ls="--", marker=marker[scheme],
                            label=f"W-boot$\\to${scheme}")
                ax_cov.plot(wb.lambda0, wb.lambda_min_after, color=col0, ls="--",
                            marker=marker[scheme], label=f"W-boot$\\to${scheme}")
        ax_cov.axhline(c / beta, color="k", ls=":", lw=1.1, label=r"$c/\beta$")
        for ax in (ax_hit, ax_cov):
            ax.set_xscale("log")
        ax_hit.set_yscale("log")
        ax_cov.set_yscale("log")
        ax_hit.set_title(rf"{TARGET_TEX[target]}: iterations to gap $<10^{{-1}}$")
        ax_cov.set_title(rf"{TARGET_TEX[target]}: $\lambda_{{\min}}(C)$ after warm-up")
        ax_cov.set_xlabel(r"initial $\lambda_0$")
        ax_hit.set_ylabel("warm-up iterations")
        ax_cov.set_ylabel(r"$\lambda_{\min}$ after step 1")
        ax_hit.grid(True, which="both", alpha=0.25)
        ax_cov.grid(True, which="both", alpha=0.25)
    axes[0][0].legend(fontsize=7.5, loc="upper right")
    axes[1][0].legend(fontsize=7.5, loc="lower right")
    fig.suptitle(rf"One Wasserstein/Bures bootstrap ($c={c:g}$) removes the "
                 r"$\lambda_0$-dependent warm-up", y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    return fig


# ---------------------------------------------------------------------------
# fig_stl_noise_floor: final stochastic floor vs dt (Gaussian | log-cosh)
# ---------------------------------------------------------------------------

def fig_stl_noise_floor(floor_df):
    """Final STL/raw noise floor vs ``dt`` for both targets, with the ``dt*Psi`` line."""
    targets = [t for t in ["gaussian", "log_cosh"] if t in set(floor_df.target.unique())]
    fig, axes = plt.subplots(1, len(targets), figsize=(5.4 * len(targets), 4.2),
                             squeeze=False)
    for col, target in enumerate(targets):
        ax = axes[0][col]
        sub = floor_df[floor_df.target == target]
        for method in ["kl_raw", "kl_stl", "riemannian_raw", "riemannian_stl"]:
            s = sub[sub.method == method].sort_values("dt")
            if s.empty:
                continue
            st = STL_METHOD_STYLE[method]
            ax.plot(s.dt, _clip(s.tail_median_gap, FLOOR_DISP), color=st["color"],
                    ls=st["ls"], marker="o", ms=4, label=st["label"])
        # dt * Psi_star reference (informative only for the non-Gaussian target).
        s = sub.sort_values("dt")
        psi_star = float(s.psi_star.iloc[0]) if not s.empty else 0.0
        if psi_star > 0:
            ax.plot(s.dt, s.dt * psi_star, color="k", ls="-.", lw=1.3,
                    label=r"$\Delta t\,\Psi_\star$")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.axhline(FLOOR_DISP, color="0.6", ls=":", lw=1.0)
        ax.set_xlabel(r"step size $\Delta t$")
        if col == 0:
            ax.set_ylabel("tail energy-gap floor")
        ax.set_title(rf"{TARGET_TEX[target]} ($\Psi_\star={psi_star:.2g}$)")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("STL noise floor: zero for the Gaussian target, "
                 r"$O(\Delta t\,\Psi)$ for log-cosh", y=0.998)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


# ---------------------------------------------------------------------------
# fig_three_stage: covariance, energy gap, local norm with stage markers
# ---------------------------------------------------------------------------

def fig_three_stage(long_df, scheme="kl", kappa=None, lambda0=None):
    """Covariance, energy gap, and FR mean-drift norm with burn-in stage markers."""
    df = long_df[(long_df.experiment == "covariance_bootstrap_scaling")
                 & (long_df.scheme == scheme)].copy()
    if df.empty:
        return None
    if kappa is None:
        kappa = sorted(df.kappa.unique())[-1]
    df = df[np.isclose(df.kappa, kappa)]
    df = df[df.d == sorted(df.d.unique())[-1]]
    if lambda0 is None:
        lambda0 = sorted(df.lambda0.unique())[0]        # smallest lambda0
    s = df[_match(df.lambda0, lambda0)].sort_values("iteration")
    if s.empty:
        return None
    beta = float(s.beta.iloc[0])
    thr = 1.0 / (2.0 * beta)
    lmn = s.lambda_min_C.values
    it = s.iteration.values
    ncov = int(it[np.argmax(lmn >= thr)]) if np.any(lmn >= thr) else None
    gap = _clip(s.energy_gap.values)
    hit_idx = np.where(gap <= 1e-3)[0]
    n_hit = int(it[hit_idx[0]]) if hit_idx.size else None

    fig, axes = plt.subplots(3, 1, figsize=(6.4, 7.4), sharex=True)
    axes[0].semilogy(it, _clip(lmn), color="#1f77b4", lw=1.8)
    axes[0].semilogy(it, _clip(s.L_n.values), color="#1f77b4", ls="--", lw=1.2,
                     label=r"$L_n$")
    axes[0].axhline(thr, color="k", ls=":", lw=1.0, label=r"$1/(2\beta)$")
    axes[0].set_ylabel(r"$\lambda_{\min}(C_n)$")
    axes[0].legend(fontsize=8, loc="lower right")
    axes[1].semilogy(it, gap, color="k", lw=1.8)
    axes[1].set_ylabel("energy gap")
    axes[2].semilogy(it, _clip(s.local_norm.values), color="#2ca02c", lw=1.8)
    axes[2].set_ylabel(r"FR drift $\|C^{1/2}G\|$")
    axes[2].set_xlabel("iteration $n$")
    for ax in axes:
        if ncov is not None:
            ax.axvline(ncov, color="#ff7f0e", ls="-", lw=1.1, alpha=0.8)
        if n_hit is not None:
            ax.axvline(n_hit, color="#8c564b", ls="-", lw=1.1, alpha=0.8)
        ax.grid(True, which="both", alpha=0.25)
    lbl = []
    if ncov is not None:
        lbl.append(rf"$N_{{\mathrm{{cov}}}}={ncov}$")
    if n_hit is not None:
        lbl.append(rf"gap$<10^{{-3}}$ at $n={n_hit}$")
    axes[0].set_title(rf"Three stages ({SCHEME_STYLE[scheme]['label']}, "
                      rf"$\kappa={kappa:g}$, $\lambda_0={lambda0:g}$): "
                      + ", ".join(lbl))
    fig.tight_layout()
    return fig
