"""Regenerate every figure and table fragment used by the LaTeX reports.

Reads only the final experiment CSVs and writes publication-ready PDF/PNG
figures and ``booktabs`` LaTeX table fragments into ``reports/assets/``. No
dynamics are re-run; this is pure post-processing of committed outputs.

Final inputs
------------
* omega/tau flows :  ``outputs/gaussian_grid/{summary,results_long}.csv``
                     ``outputs/logconcave_grid/{summary,results_long}.csv``
* local rate      :  ``outputs/natural_gradient_local_rate/operator_grid/results_long.csv``

Usage
-----
    python reports/make_report_assets.py [--repo-root .]

Outputs
-------
    reports/assets/figs/*.pdf, *.png
    reports/assets/tab_*.tex
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
sys.path.insert(0, _ROOT)

from src.common.plotting_style import apply_style, save_figure  # noqa: E402
from src.omega_tau_modes.plotting_utils import (  # noqa: E402
    plot_heatmap_with_not_reached_mask, speedup_diverging_norm,
)
from src.wfr_gradient_flow.plotting import (  # noqa: E402
    METHOD_STYLE, METHOD_ORDER, PHASE_METHODS, TARGET_TITLE, clip_gap, GAP_FLOOR,
)

# ---------------------------------------------------------------------------
# Paths and shared constants
# ---------------------------------------------------------------------------
ASSETS = os.path.join(_HERE, "assets")
FIGS = os.path.join(ASSETS, "figs")

GAUSS_DIR = os.path.join(_ROOT, "outputs", "gaussian_grid")
LOGC_DIR = os.path.join(_ROOT, "outputs", "logconcave_grid")
LR_DIR = os.path.join(_ROOT, "outputs", "natural_gradient_local_rate")
WFR_DIR = os.path.join(_ROOT, "outputs", "wfr_gradient_flow")

# WFR display order / labels.
WFR_TARGET_ORDER = ["gaussian", "smooth_log_cosh"]
WFR_METHOD_TEX = {
    "fr_only": "FR-only", "w_only": "W-only", "wfr_fixed": "WFR-fixed",
    "wfr_theory": "WFR-theory", "wfr_adaptive": "WFR-adaptive",
}
WFR_TARGET_TEX = {"gaussian": "Gaussian", "smooth_log_cosh": "smooth log-cosh"}

# Display order / labels for the five mode initializations.
INIT_ORDER = ["mean_only", "volume_high", "volume_low", "shape_only", "mixed"]
INIT_LABEL = {
    "mean_only": "mean", "volume_high": "vol-high", "volume_low": "vol-low",
    "shape_only": "shape", "mixed": "mixed",
}
# Potential families (local rate) display order / labels.
FAM_ORDER = ["gaussian", "separable", "additive_index", "random_feature", "radial_tail"]
FAM_LABEL = {
    "gaussian": "Gaussian", "separable": "separable",
    "additive_index": "additive-index", "random_feature": "random-feature",
    "radial_tail": "radial-tail",
}
FAM_COLOR = {
    "gaussian": "#444444", "separable": "#1f77b4", "additive_index": "#2ca02c",
    "random_feature": "#d62728", "radial_tail": "#9467bd",
}
FAM_MARKER = {
    "gaussian": "o", "separable": "s", "additive_index": "^",
    "random_feature": "D", "radial_tail": "v",
}


def _savefig(fig, name):
    paths = save_figure(fig, os.path.join(FIGS, name))
    plt.close(fig)
    print(f"  fig  {name}  ->  {os.path.basename(paths[0])}, {os.path.basename(paths[1])}")


def _write_table(name, body):
    path = os.path.join(ASSETS, name)
    with open(path, "w") as fh:
        fh.write(body if body.endswith("\n") else body + "\n")
    print(f"  tab  {name}")


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    return f"{x:.{nd}f}"


# ===========================================================================
# omega / tau experiments
# ===========================================================================

def _ot_pivot(summary, tau_type, tol_col):
    """Time-to-tolerance matrix (rows = INIT_ORDER, cols = sorted omega)."""
    omegas = sorted(summary["omega"].unique())
    sub = summary[summary["tau_type"] == tau_type]
    mat = np.full((len(INIT_ORDER), len(omegas)), np.nan)
    for i, init in enumerate(INIT_ORDER):
        for j, om in enumerate(omegas):
            r = sub[(sub["init_name"] == init) & (sub["omega"] == om)]
            if len(r):
                v = float(r[tol_col].iloc[0])
                mat[i, j] = v if np.isfinite(v) else np.nan
    return mat, omegas


def _ot_speedup(summary, tol_col):
    """Speedup matrices T(tau_-)/T(0) and T(tau_+)/T(0): rows=init, cols=omega."""
    omegas = sorted(summary["omega"].unique())
    neg = np.full((len(INIT_ORDER), len(omegas)), np.nan)
    pos = np.full((len(INIT_ORDER), len(omegas)), np.nan)
    for i, init in enumerate(INIT_ORDER):
        for j, om in enumerate(omegas):
            base = summary[(summary.init_name == init) & (summary.omega == om)
                           & (summary.tau_type == "zero")][tol_col]
            if not len(base) or not np.isfinite(base.iloc[0]) or base.iloc[0] == 0:
                continue
            t0 = float(base.iloc[0])
            for mat, tt in ((neg, "negative"), (pos, "positive")):
                r = summary[(summary.init_name == init) & (summary.omega == om)
                            & (summary.tau_type == tt)][tol_col]
                if len(r) and np.isfinite(r.iloc[0]):
                    mat[i, j] = float(r.iloc[0]) / t0
    return neg, pos, omegas


def _fig_time_to_tol(summary, tol_col, T, name, title):
    mat, omegas = _ot_pivot(summary, "zero", tol_col)
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    plot_heatmap_with_not_reached_mask(
        ax, mat, T=T,
        xticklabels=[f"{o:g}" for o in omegas],
        yticklabels=[INIT_LABEL[i] for i in INIT_ORDER],
        cmap="viridis_r", cbar_label="time to tol.", fig=fig)
    ax.set_xlabel(r"$\omega$")
    ax.set_title(title)
    fig.tight_layout()
    _savefig(fig, name)


def _fig_tau_speedup(summary, tol_col, name, title):
    neg, pos, omegas = _ot_speedup(summary, tol_col)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0), sharey=True)
    norm = speedup_diverging_norm(0.5, 1.0, 1.5)
    labels = [INIT_LABEL[i] for i in INIT_ORDER]
    for ax, mat, sub in ((axes[0], neg, r"$\tau_-=-\omega/2n$"),
                         (axes[1], pos, r"$\tau_+=+\omega/2n$")):
        im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", norm=norm, origin="upper")
        ax.set_xticks(range(len(omegas)))
        ax.set_xticklabels([f"{o:g}" for o in omegas])
        ax.set_yticks(range(len(INIT_ORDER)))
        ax.set_yticklabels(labels)
        ax.set_xlabel(r"$\omega$")
        ax.set_title(sub)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if np.isfinite(mat[i, j]):
                    ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                            fontsize=7.5,
                            color="white" if abs(mat[i, j] - 1) > 0.32 else "black")
    fig.colorbar(im, ax=axes, label=r"$T(\tau)/T(0)$", fraction=0.046, pad=0.04)
    fig.suptitle(title, y=1.02, fontsize=11)
    _savefig(fig, name)


def _fig_speedup_vs_chi(summary, tol_col, name):
    """Realized tau_- speedup vs initial trace-dominance chi_0 (log-concave)."""
    omegas = sorted(summary["omega"].unique())
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    cmap = plt.get_cmap("viridis")
    for k, om in enumerate(omegas):
        xs, ys = [], []
        for init in INIT_ORDER:
            base = summary[(summary.init_name == init) & (summary.omega == om)
                           & (summary.tau_type == "zero")]
            neg = summary[(summary.init_name == init) & (summary.omega == om)
                          & (summary.tau_type == "negative")]
            if not len(base) or not len(neg):
                continue
            t0, tn = float(base[tol_col].iloc[0]), float(neg[tol_col].iloc[0])
            chi0 = float(base["initial_chi"].iloc[0])
            if np.isfinite(t0) and t0 > 0 and np.isfinite(tn):
                xs.append(chi0)
                ys.append(tn / t0)
        ax.scatter(xs, ys, color=cmap(0.15 + 0.7 * k / max(1, len(omegas) - 1)),
                   s=42, label=rf"$\omega={om:g}$", edgecolor="white", linewidth=0.5, zorder=3)
    ax.axhline(1.0, color="0.5", lw=1.0, ls="--", zorder=1)
    ax.axhline(0.5, color="0.7", lw=0.9, ls=":", zorder=1)
    ax.set_xlabel(r"initial trace dominance $\chi_0$")
    ax.set_ylabel(r"$T(\tau_-)/T(0)$")
    ax.set_ylim(0.35, 1.25)
    ax.legend(loc="lower left")
    fig.tight_layout()
    _savefig(fig, name)


def _tab_ot_speedup(gauss, logc):
    """Speedup factors T(tau)/T(0) at omega=0.5 for both targets, by init."""
    def col(summary, tol_col, om=0.5):
        out = {}
        for init in INIT_ORDER:
            base = summary[(summary.init_name == init) & (summary.omega == om)
                           & (summary.tau_type == "zero")][tol_col]
            t0 = float(base.iloc[0]) if len(base) else np.nan
            row = [t0]
            for tt in ("negative", "positive"):
                r = summary[(summary.init_name == init) & (summary.omega == om)
                            & (summary.tau_type == tt)][tol_col]
                rr = float(r.iloc[0]) if len(r) else np.nan
                row.append(rr / t0 if np.isfinite(t0) and t0 else np.nan)
            out[init] = row
        return out
    g = col(gauss, "time_to_1e_minus_4")
    l = col(logc, "time_to_1e_minus_2")
    lines = [
        r"\begin{tabular}{lcccccc}",
        r"\toprule",
        r" & \multicolumn{3}{c}{Gaussian ($10^{-4}$)} & \multicolumn{3}{c}{log-concave ($10^{-2}$)} \\",
        r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
        r"initialization & $T(0)$ & $T(\tau_-)/T_0$ & $T(\tau_+)/T_0$ & $T(0)$ & $T(\tau_-)/T_0$ & $T(\tau_+)/T_0$ \\",
        r"\midrule",
    ]
    for init in INIT_ORDER:
        gv, lv = g[init], l[init]
        lines.append(
            f"{INIT_LABEL[init]} & {_fmt(gv[0],1)} & {_fmt(gv[1],2)} & {_fmt(gv[2],2)} "
            f"& {_fmt(lv[0],1)} & {_fmt(lv[1],2)} & {_fmt(lv[2],2)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _write_table("tab_ot_speedup.tex", "\n".join(lines))


def _tab_ot_metadata(gauss_meta, logc_meta):
    g, l = gauss_meta, logc_meta
    lines = [
        r"\begin{tabular}{lll}",
        r"\toprule",
        r"quantity & Gaussian target & log-concave target \\",
        r"\midrule",
        rf"dimension $n$ & $\{{2,5,10\}}$ & $5$ \\",
        rf"$\omega$ grid & $\{{0.125,0.25,0.5,1,2\}}$ & $\{{0.25,0.5,1\}}$ \\",
        rf"$\tau$ per $\omega$ & $\{{-\omega/2n,\,0,\,+\omega/2n\}}$ & $\{{-\omega/2n,\,0,\,+\omega/2n\}}$ \\",
        rf"step $\Delta t$ & $0.02$ & ${l['dt']:g}$ \\",
        rf"horizon $T$ & $20$ & ${l['T']:g}$ \\",
        rf"coupling $\rho$ & -- & ${l['rho']:g}$ \\",
        rf"features $m=4n$ & -- & ${l['m_features']}$ \\",
        rf"dynamics samples $K$ & exact & ${l['K']}$ \\",
        rf"reference objective $F_\star$ & $0$ (exact) & ${l['F_star']:.4f}$ \\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    _write_table("tab_ot_metadata.tex", "\n".join(lines))


def build_omega_tau_assets():
    print("omega/tau assets:")
    gauss = pd.read_csv(os.path.join(GAUSS_DIR, "summary.csv"))
    logc = pd.read_csv(os.path.join(LOGC_DIR, "summary.csv"))
    gauss5 = gauss[gauss.n == 5]
    import json
    with open(os.path.join(LOGC_DIR, "target_metadata.json")) as fh:
        lc_meta = json.load(fh)

    _fig_time_to_tol(gauss5, "time_to_1e_minus_4", 20.0,
                     "fig_ot_gaussian_time_to_tol",
                     r"Gaussian, $n=5$: time to $F/F_0\leq10^{-4}$ ($\tau=0$)")
    _fig_tau_speedup(gauss5, "time_to_1e_minus_4",
                     "fig_ot_gaussian_tau_speedup", r"Gaussian target ($n=5$)")
    _fig_time_to_tol(logc, "time_to_1e_minus_2", 40.0,
                     "fig_ot_logconcave_time_to_tol",
                     r"log-concave, $n=5$: time to gap $\leq10^{-2}$ ($\tau=0$)")
    _fig_tau_speedup(logc, "time_to_1e_minus_2",
                     "fig_ot_logconcave_tau_speedup", r"log-concave target ($n=5$)")
    _fig_speedup_vs_chi(logc, "time_to_1e_minus_2", "fig_ot_speedup_vs_chi")
    _tab_ot_speedup(gauss5, logc)
    _tab_ot_metadata({}, lc_meta)


# ===========================================================================
# natural-gradient local rate (final gpu_lowdim_operator_full grid only)
# ===========================================================================

def _load_lr():
    path = os.path.join(LR_DIR, "operator_grid", "results_long.csv")
    df = pd.read_csv(path)
    df = df[df["status"] == "ok"].copy()
    return df


def _seed_mean(df, value):
    """Mean and std over seeds, indexed by (potential_family, N_theta, kappa_target)."""
    g = df.groupby(["potential_family", "N_theta", "kappa_target"])[value]
    return g.mean(), g.std()


def _fig_gamma_vs_Ntheta(df, kappa, name):
    """gamma_loc vs N_theta at fixed kappa, one line per family (mean +/- std)."""
    mean, std = _seed_mean(df, "gamma_loc")
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    for fam in FAM_ORDER:
        Ns, ms, ss = [], [], []
        for N in sorted(df.N_theta.unique()):
            key = (fam, N, float(kappa))
            if key in mean.index:
                Ns.append(N); ms.append(mean[key]); ss.append(std.get(key, 0.0) or 0.0)
        if not Ns:
            continue
        ms, ss = np.array(ms), np.array(ss)
        ax.plot(Ns, ms, marker=FAM_MARKER[fam], color=FAM_COLOR[fam],
                label=FAM_LABEL[fam], markersize=4)
        ax.fill_between(Ns, ms - ss, ms + ss, color=FAM_COLOR[fam], alpha=0.15, lw=0)
    ax.set_xlabel(r"dimension $N_\theta$")
    ax.set_ylabel(r"local rate $\gamma_{\mathrm{loc}}$")
    ax.set_title(rf"$\kappa={kappa:g}$")
    ax.set_ylim(0.55, 1.03)
    ax.legend(loc="lower right", ncol=2)
    fig.tight_layout()
    _savefig(fig, name)


def _fig_gamma_vs_kappa(df, N, name):
    """gamma_loc vs kappa at fixed N_theta, one line per family."""
    mean, std = _seed_mean(df, "gamma_loc")
    kappas = sorted(df.kappa_target.unique())
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    for fam in FAM_ORDER:
        ks, ms, ss = [], [], []
        for k in kappas:
            key = (fam, int(N), float(k))
            if key in mean.index:
                ks.append(k); ms.append(mean[key]); ss.append(std.get(key, 0.0) or 0.0)
        if not ks:
            continue
        ms, ss = np.array(ms), np.array(ss)
        ax.plot(ks, ms, marker=FAM_MARKER[fam], color=FAM_COLOR[fam],
                label=FAM_LABEL[fam], markersize=4)
        ax.fill_between(ks, ms - ss, ms + ss, color=FAM_COLOR[fam], alpha=0.15, lw=0)
    ax.set_xscale("log")
    ax.set_xlabel(r"conditioning $\kappa$")
    ax.set_ylabel(r"local rate $\gamma_{\mathrm{loc}}$")
    ax.set_title(rf"$N_\theta={N}$")
    ax.set_ylim(0.55, 1.03)
    ax.legend(loc="lower left", ncol=2)
    fig.tight_layout()
    _savefig(fig, name)


def _fig_gamma_grid(df, name):
    """Small-multiples: gamma_loc vs N_theta, one panel per kappa, lines per family."""
    mean, _ = _seed_mean(df, "gamma_loc")
    kappas = sorted(df.kappa_target.unique())
    ncol = 3
    nrow = int(np.ceil(len(kappas) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(8.4, 2.6 * nrow),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for a, kappa in enumerate(kappas):
        ax = axes[a]
        for fam in FAM_ORDER:
            Ns, ms = [], []
            for N in sorted(df.N_theta.unique()):
                key = (fam, N, float(kappa))
                if key in mean.index:
                    Ns.append(N); ms.append(mean[key])
            if Ns:
                ax.plot(Ns, ms, marker=FAM_MARKER[fam], color=FAM_COLOR[fam],
                        label=FAM_LABEL[fam], markersize=3, lw=1.3)
        ax.set_title(rf"$\kappa={kappa:g}$", fontsize=10)
        ax.set_ylim(0.55, 1.03)
        ax.grid(True, alpha=0.3)
    for a in range(len(kappas), len(axes)):
        axes[a].axis("off")
    for ax in axes[-ncol:]:
        ax.set_xlabel(r"$N_\theta$")
    for r in range(nrow):
        axes[r * ncol].set_ylabel(r"$\gamma_{\mathrm{loc}}$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    _savefig(fig, name)


def _fig_lambda_artifact(df, name):
    """Separable: full-operator Lambda_hat (noisy) vs dimension-free exact benchmark."""
    sub = df[df.potential_family == "separable"]
    mean_full, _ = _seed_mean(sub, "Lambda_hat_full_sym")
    mean_exact, _ = _seed_mean(sub, "Lambda_hat_separable_exact")
    kappas = [5.0, 100.0]
    fig, axes = plt.subplots(1, len(kappas), figsize=(7.4, 3.1), sharey=True)
    for ax, kappa in zip(np.atleast_1d(axes), kappas):
        Ns = sorted(sub.N_theta.unique())
        full = [mean_full.get(("separable", N, kappa), np.nan) for N in Ns]
        exact = [mean_exact.get(("separable", N, kappa), np.nan) for N in Ns]
        ax.plot(Ns, full, marker="s", color="#d62728",
                label=r"$\widehat\Lambda_{\mathrm{full\;sym}}$ (MC)")
        ax.plot(Ns, exact, marker="o", color="#1f77b4",
                label=r"$\Lambda_{\mathrm{exact}}$ (Gauss--Hermite)")
        ax.set_xlabel(r"$N_\theta$")
        ax.set_title(rf"separable, $\kappa={kappa:g}$")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(r"operator-norm estimate $\widehat\Lambda$")
    axes[0].legend(loc="upper left")
    fig.tight_layout()
    _savefig(fig, name)


def _tab_lr_family_ranges(df):
    """gamma_loc range (min..max of seed-mean over grid) and inverse rate by family."""
    mean, _ = _seed_mean(df, "gamma_loc")
    lines = [
        r"\begin{tabular}{lccc}",
        r"\toprule",
        r"family & $\min\gamma_{\mathrm{loc}}$ & $\max\gamma_{\mathrm{loc}}$ & $\max(1/\gamma_{\mathrm{loc}})$ \\",
        r"\midrule",
    ]
    for fam in FAM_ORDER:
        vals = [v for (f, _, _), v in mean.items() if f == fam and np.isfinite(v)]
        if not vals:
            continue
        lo, hi = min(vals), max(vals)
        lines.append(f"{FAM_LABEL[fam]} & {_fmt(lo,3)} & {_fmt(hi,3)} & {_fmt(1.0/lo,3)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _write_table("tab_lr_family_ranges.tex", "\n".join(lines))


def _tab_lr_metadata(df):
    import json
    with open(os.path.join(LR_DIR, "operator_grid", "config.json")) as fh:
        meta = json.load(fh)
    cfg = meta["config"]
    r0 = df.iloc[0]
    Ns = sorted(df.N_theta.unique())
    ks = sorted(df.kappa_target.unique())
    lines = [
        r"\begin{tabular}{ll}",
        r"\toprule",
        r"quantity & value \\",
        r"\midrule",
        rf"dimensions $N_\theta$ & $1,\dots,{max(Ns)}$ ({len(Ns)} values) \\",
        rf"conditioning $\kappa$ & $\{{{','.join(f'{k:g}' for k in ks)}\}}$ \\",
        rf"families & Gaussian, separable, additive-index, random-feature, radial-tail \\",
        rf"seeds & $\{{0,1,2\}}$ \\",
        rf"rows & ${len(df)}$ (all status ok) \\",
        rf"Monte Carlo samples $M$ & ${int(r0['M_mc'])}=2^{{22}}$ \\",
        rf"estimator & {cfg['operator']['estimator']} ($H_{{\mathrm{{sym}}}}$) \\",
        rf"backend / device / dtype & {r0['backend']} / {r0['device_name']} / {r0['dtype']} \\",
        rf"run id & \texttt{{{meta['run_id'].replace('_', chr(92)+'_')}}} \\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    _write_table("tab_lr_metadata.tex", "\n".join(lines))


def build_local_rate_assets():
    print("local-rate assets:")
    df = _load_lr()
    _fig_gamma_vs_Ntheta(df, 100.0, "fig_lr_gamma_vs_Ntheta_k100")
    _fig_gamma_vs_kappa(df, 16, "fig_lr_gamma_vs_kappa_N16")
    _fig_gamma_grid(df, "fig_lr_gamma_vs_Ntheta_grid")
    _fig_lambda_artifact(df, "fig_lr_lambda_artifact")
    _tab_lr_family_ranges(df)
    _tab_lr_metadata(df)


# ===========================================================================
# discretization stepsize (Riemannian vs KL)
# ===========================================================================

DISC_DIR = os.path.join(_ROOT, "outputs", "natural_gradient_discretization_stepsize")


def _tab_disc_stepsize(step_df):
    """Empirical vs theoretical maximum stepsize, by target/lambda/method.

    Reports the largest stable, monotone, and accurate stepsizes together with
    the monotone/theory and accurate/theory ratios, so the proof-artifact gap is
    visible in both the monotone and the (stricter) accurate class.
    """
    tt = {"gaussian_posterior": "Gaussian", "literature_logconcave": "non-smooth",
          "smooth_logconcave": "smooth"}
    lines = [
        r"\begin{tabular}{lllccccc}",
        r"\toprule",
        r"target & $\lambda$ & method & $\Delta t_{\mathrm{theory}}$ & "
        r"$\Delta t_{\max}^{\mathrm{stab}}$ & $\Delta t_{\max}^{\mathrm{mono}}$ & "
        r"$\Delta t_{\max}^{\mathrm{acc}}$ & "
        r"$\tfrac{\mathrm{mono}}{\mathrm{theory}}$ / "
        r"$\tfrac{\mathrm{acc}}{\mathrm{theory}}$ \\",
        r"\midrule",
    ]
    order = {"gaussian_posterior": 0, "smooth_logconcave": 1, "literature_logconcave": 2}
    step_df = step_df.sort_values(
        by=["target_name", "lambda", "method"],
        key=lambda s: s.map(order) if s.name == "target_name" else s)
    for _, r in step_df.iterrows():
        avail = r["theory_bound_available"]
        theory_s = _fmt_sci(r["dt_theory_for_method"]) if avail else "--"
        mono_ratio_s = _fmt_sci(r["monotone_over_theory_ratio"]) if avail else "--"
        acc_ratio_s = _fmt_sci(r["accurate_over_theory_ratio"]) if avail else "--"
        ratio_s = f"{mono_ratio_s} / {acc_ratio_s}" if avail else "--"
        lines.append(
            f"{tt.get(r['target_name'], r['target_name'])} & {r['lambda']:g} & "
            f"{r['method']} & {theory_s} & {_fmt(r['dt_max_stable'], 3)} & "
            f"{_fmt(r['dt_max_monotone'], 3)} & {_fmt(r['dt_max_accurate'], 3)} & "
            f"{ratio_s} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _write_table("tab_disc_stepsize.tex", "\n".join(lines))


def _fmt_sci(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    if x == 0:
        return "0"
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / 10.0 ** exp
    if -2 <= exp <= 2:
        return f"{x:.3g}"
    return rf"${mant:.1f}\times10^{{{exp}}}$"


def _tab_disc_metadata(meta):
    """Per-target reference optimum and theory constants."""
    tt = {"gaussian_posterior": "Gaussian posterior",
          "smooth_logconcave": "smooth log-concave",
          "literature_logconcave": "non-smooth log-concave"}
    lines = [
        r"\begin{tabular}{llcccc}",
        r"\toprule",
        r"target & $\lambda$ & $\alpha$ & $\beta$ & $\mathcal{E}_\star$ & global smooth \\",
        r"\midrule",
    ]
    keys = sorted(meta["targets"].keys(),
                  key=lambda k: ({"gaussian_posterior": 0, "smooth_logconcave": 1,
                                  "literature_logconcave": 2}[k.split("__")[0]],
                                 meta["targets"][k]["lambda"]))
    for k in keys:
        m = meta["targets"][k]
        a = _fmt(m["alpha"], 3) if m["alpha"] is not None else "--"
        b = _fmt(m["beta"], 3) if m["beta"] is not None else "--"
        lines.append(
            f"{tt.get(m['target_name'], m['target_name'])} & {m['lambda']:g} & "
            f"{a} & {b} & {_fmt(m['F_star'], 4)} & "
            f"{'yes' if m['has_theory'] else 'no'} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _write_table("tab_disc_metadata.tex", "\n".join(lines))


# Stepsizes used by the matched-stepsize convergence comparison (must match the
# per-target sets the appendix-derived figures use).
_CONV_DT = {
    "gaussian_posterior": [0.05, 0.1, 0.5, 1.0],
    "smooth_logconcave": [0.05, 0.1, 0.5, 1.0],
    "literature_logconcave": [0.002, 0.005, 0.01, 0.02],
}
_GAP_EPS = 1e-16


def _tab_disc_convergence(long_df):
    """Head-to-head convergence speed: terminal energy gap, Riemannian vs KL.

    For each (target, dt) we report the terminal energy gap (geometric mean over
    the three lambda) of each scheme at a matched stepsize, the winner, and the
    factor by which the winner's terminal gap is smaller. This isolates the
    per-step convergence-rate comparison from the stepsize-range question.
    """
    tt = {"gaussian_posterior": "Gaussian", "smooth_logconcave": "smooth",
          "literature_logconcave": "non-smooth"}

    def gmean_gap(target, method, dt):
        gaps = []
        for lam in (0.01, 0.1, 1.0):
            r = long_df[(long_df.target_name == target)
                        & (np.isclose(long_df["lambda"], lam))
                        & (long_df.method == method) & (np.isclose(long_df.dt, dt))]
            if len(r):
                r = r.sort_values("n")
                gaps.append(max(float(r.energy_gap.iloc[-1]), _GAP_EPS))
        return float(np.exp(np.mean(np.log(gaps)))) if gaps else np.nan

    lines = [
        r"\begin{tabular}{llcccc}",
        r"\toprule",
        r"target & $\Delta t$ & gap$_{\mathrm{Riem}}$ & gap$_{\mathrm{KL}}$ & "
        r"winner & factor \\",
        r"\midrule",
    ]
    for ti, target in enumerate(["gaussian_posterior", "smooth_logconcave",
                                 "literature_logconcave"]):
        for di, dt in enumerate(_CONV_DT[target]):
            gr, gk = gmean_gap(target, "riemannian", dt), gmean_gap(target, "kl", dt)
            if not (np.isfinite(gr) and np.isfinite(gk)):
                continue
            winner = "Riem." if gr < gk else "KL"
            factor = max(gr, gk) / min(gr, gk)
            tlabel = tt[target] if di == 0 else ""
            lines.append(
                f"{tlabel} & {dt:g} & {_fmt_sci(gr)} & {_fmt_sci(gk)} & "
                f"{winner} & {_fmt_sci(factor)} \\\\")
        if ti < 2:
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _write_table("tab_disc_convergence.tex", "\n".join(lines))


# Rate-benchmark aggregation: geometric means over the runs in each cell.
_RATE_TT = {"gaussian_posterior": "Gaussian", "smooth_logconcave": "smooth"}


def _gmean(vals):
    v = np.asarray([x for x in vals if x is not None and np.isfinite(x) and x > 0],
                   dtype=float)
    return float(np.exp(np.mean(np.log(v)))) if v.size else np.nan


def _tab_discretization_rate_summary(rate_summary):
    """Aggregated rate benchmark: observed vs theoretical contraction by c.

    One row per (target, method, c): geometric-mean observed rate, theoretical
    rate, observed/theory ratio, terminal slack, run count, and floor-limited
    count. Aggregation is over the three lambda values.
    """
    df = rate_summary[rate_summary.status == "ok"].copy()
    lines = [
        r"\begin{tabular}{lllcccccc}",
        r"\toprule",
        r"target & method & $c$ & $\bar r_{\mathrm{obs}}$ & $\bar r_{\mathrm{th}}$ & "
        r"$\overline{r_{\mathrm{obs}}/r_{\mathrm{th}}}$ & "
        r"$\overline{\mathrm{slack}}$ & $n$ & floor \\",
        r"\midrule",
    ]
    targets = [t for t in ("gaussian_posterior", "smooth_logconcave")
               if t in set(df.target.unique())]
    for ti, target in enumerate(targets):
        for method in ("riemannian", "kl"):
            cs = sorted(df[(df.target == target) & (df.method == method)].c.unique())
            for ci, c in enumerate(cs):
                cell = df[(df.target == target) & (df.method == method)
                          & (np.isclose(df.c, c))]
                if cell.empty:
                    continue
                ratio = cell.r_hat_terminal.values / cell.r_theory.values
                tlabel = _RATE_TT[target] if (method == "riemannian" and ci == 0) else ""
                mlabel = method if ci == 0 else ""
                n_floor = int(cell.floor_limited_final.sum())
                lines.append(
                    f"{tlabel} & {mlabel} & {c:g} & "
                    f"{_fmt_sci(_gmean(cell.r_hat_terminal))} & "
                    f"{_fmt_sci(_gmean(cell.r_theory))} & "
                    f"{_fmt_sci(_gmean(ratio))} & "
                    f"{_fmt_sci(_gmean(cell.terminal_slack))} & "
                    f"{len(cell)} & {n_floor} \\\\")
            if not (ti == len(targets) - 1 and method == "kl"):
                lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _write_table("tab_discretization_rate_summary.tex", "\n".join(lines))


def build_discretization_assets():
    """Figures + tables for the Riemannian-vs-KL discretization report.

    Reuses the per-group plot builders (single source of truth) writing into
    ``reports/assets/figs``; adds two booktabs tables read from the final CSVs.
    """
    print("discretization assets:")
    sys.path.insert(0, os.path.join(_ROOT, "scripts", "natural_gradient_discretization_stepsize"))
    import importlib
    pr = importlib.import_module("plot_results")
    import json

    long_df = pd.read_csv(os.path.join(DISC_DIR, "results_long.csv"))
    summary = pd.read_csv(os.path.join(DISC_DIR, "summary.csv"))
    step_df = pd.read_csv(os.path.join(DISC_DIR, "stepsize_summary.csv"))
    scalar_df = pd.read_csv(os.path.join(DISC_DIR, "scalar_diagnostic.csv"))
    with open(os.path.join(DISC_DIR, "target_metadata.json")) as fh:
        meta = json.load(fh)

    for tname in ("gaussian_posterior", "literature_logconcave", "smooth_logconcave"):
        pr.fig_energy_gap(long_df, tname, FIGS, meta, summary=summary)
        pr.fig_stability_heatmap(summary, step_df, tname, FIGS)
    pr.fig_theory_vs_empirical(step_df, FIGS)
    pr.fig_scalar_diagnostic(scalar_df, FIGS)
    pr.fig_time_to_tolerance(summary, FIGS)
    # Convergence-speed study (matched-stepsize Riemannian vs KL, convergent dt).
    pr.fig_convergence_speed_grid(long_df, summary, FIGS)
    pr.fig_iterations_vs_stepsize(summary, FIGS)
    _tab_disc_stepsize(step_df)
    _tab_disc_metadata(meta)
    _tab_disc_convergence(long_df)

    # Supplementary theoretical-rate benchmark (rate_* CSVs). Built only when
    # present so the rest of the discretization assets never block on it.
    rate_summary_path = os.path.join(DISC_DIR, "rate_summary.csv")
    if os.path.exists(rate_summary_path):
        rate_long = pd.read_csv(os.path.join(DISC_DIR, "rate_results_long.csv"))
        rate_summary = pd.read_csv(rate_summary_path)
        rate_tol = pd.read_csv(os.path.join(DISC_DIR, "rate_tolerance_summary.csv"))
        pr.build_rate_figs(rate_long, rate_summary, rate_tol, FIGS)
        _tab_discretization_rate_summary(rate_summary)
    else:
        print("  [skip] rate benchmark assets: rate_summary.csv not found "
              "(run run_rate_benchmark.py)")


# ===========================================================================
# WFR Gaussian gradient flow
# ===========================================================================

def _wfr_fmt_int(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)) or int(x) < 0:
        return r"$\infty$"
    return f"{int(x)}"


def _wfr_fmt_sci(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    if x == 0:
        return "0"
    if x < 1e-12:
        return r"$<\!10^{-12}$"
    exp = int(np.floor(np.log10(abs(x))))
    if -2 <= exp <= 2:
        return f"{x:.2g}"
    mant = x / 10.0 ** exp
    return rf"${mant:.1f}\times10^{{{exp}}}$"


def _wfr_fig_phase_separation(long_df, target_name, warmup_iters=20):
    """Gap vs iteration for one target, all (Lambda, eps): four methods.

    Two columns -- a warmup zoom (linear x, log y up to ``warmup_iters``) and the
    full trajectory (log x, log y) -- one row per ``(Lambda, eps)`` combination.
    The warmup column exposes the Wasserstein early-descent advantage; the full
    column exposes the Fisher--Rao tail advantage and the WFR combination.
    """
    sub = long_df[long_df.target_name == target_name]
    combos = (sub[["Lambda", "epsilon"]].drop_duplicates()
              .sort_values(["Lambda", "epsilon"]).itertuples(index=False))
    combos = list(combos)
    nrow = len(combos)
    fig, axes = plt.subplots(nrow, 2, figsize=(10.5, 3.0 * nrow),
                             squeeze=False)
    for i, (Lambda, eps) in enumerate(combos):
        cell = sub[(np.isclose(sub.Lambda, Lambda)) & (np.isclose(sub.epsilon, eps))]
        ax_w, ax_f = axes[i, 0], axes[i, 1]
        for method in PHASE_METHODS:
            s = cell[cell.method == method].sort_values("iteration")
            if s.empty:
                continue
            st = METHOD_STYLE[method]
            b = s.iteration.values
            g = clip_gap(s.objective_gap.values)
            ax_w.plot(b, g, color=st["color"], ls=st["ls"], label=st["label"], lw=1.6)
            ax_f.plot(b, g, color=st["color"], ls=st["ls"], label=st["label"], lw=1.6)
        ax_w.set_xlim(0, warmup_iters)
        ax_w.set_yscale("log")
        ax_f.set_xscale("log")
        ax_f.set_yscale("log")
        ax_w.set_ylabel(rf"$\Lambda={Lambda:g},\ \varepsilon={eps:g}$"
                        "\n" r"energy gap")
        for ax in (ax_w, ax_f):
            ax.grid(True, which="both", alpha=0.25)
            ax.set_ylim(GAP_FLOOR, None)
        if i == 0:
            ax_w.set_title("warmup (linear iterations)")
            ax_f.set_title("full trajectory (log iterations)")
        if i == nrow - 1:
            ax_w.set_xlabel("iteration")
            ax_f.set_xlabel("iteration")
    axes[0, 1].legend(fontsize=8, loc="upper right", framealpha=0.9)
    fig.suptitle(f"{TARGET_TITLE[target_name]}: WFR phase separation "
                 "(gap vs iteration)", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    return fig


def _wfr_fig_adaptive_schedule(long_df, summary_df):
    """The wfr_adaptive transport ``h_n beta`` and calibration ``s_n`` over a run.

    One column per target; the top row shows the normalized transport
    ``h_n beta`` (so the 0.9 ceiling is target-independent) and the bottom row the
    calibration ratio ``s_n`` driving it, both vs expectation batches, for every
    ``(Lambda, eps)``. The schedule keeps transport near the ceiling while
    ``s_n`` is small (underdispersed warmup) and decays it as ``s_n`` approaches
    one (the covariance becomes locally calibrated).
    """
    fig, axes = plt.subplots(2, len(WFR_TARGET_ORDER), figsize=(5.2 * len(WFR_TARGET_ORDER), 6.4),
                             squeeze=False, sharex="col")
    for j, target_name in enumerate(WFR_TARGET_ORDER):
        sub = long_df[(long_df.target_name == target_name)
                      & (long_df.method == "wfr_adaptive")]
        beta_map = (summary_df[summary_df.target_name == target_name]
                    .set_index(["Lambda", "epsilon"]).beta.to_dict())
        combos = (sub[["Lambda", "epsilon"]].drop_duplicates()
                  .sort_values(["Lambda", "epsilon"]).itertuples(index=False))
        cmap = plt.cm.viridis(np.linspace(0.1, 0.85, 4))
        for k, (Lambda, eps) in enumerate(combos):
            s = sub[(np.isclose(sub.Lambda, Lambda)) & (np.isclose(sub.epsilon, eps))]
            s = s[s.iteration > 0].sort_values("iteration")
            if s.empty:
                continue
            beta = beta_map.get((Lambda, eps), 1.0)
            lbl = rf"$\Lambda={Lambda:g},\varepsilon={eps:g}$"
            axes[0, j].plot(s.iteration, s.h_n * beta, color=cmap[k % 4], lw=1.5, label=lbl)
            axes[1, j].plot(s.iteration, s.s_n, color=cmap[k % 4], lw=1.5, label=lbl)
        axes[0, j].set_title(TARGET_TITLE[target_name])
        axes[0, j].axhline(0.9, color="k", ls=":", lw=1, alpha=0.6)
        axes[0, j].set_ylabel(r"normalized transport $h_n\beta$")
        axes[1, j].axhline(0.5, color="k", ls=":", lw=1, alpha=0.6)
        axes[1, j].set_ylabel(r"calibration $s_n=\lambda_{\min}(C^{1/2}(-H)C^{1/2})$")
        axes[1, j].set_xlabel("iteration")
        axes[1, j].set_xscale("log")
        for r in (0, 1):
            axes[r, j].grid(True, which="both", alpha=0.25)
        axes[1, j].set_yscale("log")
    axes[0, 0].legend(fontsize=8, loc="best")
    fig.suptitle("WFR-adaptive transport schedule and curvature calibration", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    return fig


def _wfr_fig_error_decomposition(long_df):
    """Mean-error and covariance-error vs iteration on the hardest start per target.

    Two columns per target (mean-error norm, covariance-error Frobenius), top row
    Gaussian and bottom row smooth, at ``Lambda=1000, eps=1e-3``. This exposes the
    *mechanism* behind the phase separation rather than asserting it: at the
    underdispersed start the covariance error (~10^3) dwarfs the mean error (~10),
    so the covariance is the bottleneck; W-only barely dents the covariance error
    (it inflates volume but cannot calibrate shape), while FR-only and WFR collapse
    it once the mean has arrived.
    """
    hard = {"Lambda": 1000.0, "epsilon": 1e-3}
    fig, axes = plt.subplots(len(WFR_TARGET_ORDER), 2,
                             figsize=(11.0, 3.6 * len(WFR_TARGET_ORDER)), squeeze=False)
    for i, tname in enumerate(WFR_TARGET_ORDER):
        cell = long_df[(long_df.target_name == tname)
                       & np.isclose(long_df.Lambda, hard["Lambda"])
                       & np.isclose(long_df.epsilon, hard["epsilon"])]
        for col, (ycol, ylab) in enumerate([
                ("mean_error_norm", r"mean error $\|m-m_\star\|$"),
                ("covariance_error_fro", r"covariance error $\|C-C_\star\|_F$")]):
            ax = axes[i][col]
            for method in PHASE_METHODS:
                s = cell[cell.method == method].sort_values("iteration")
                if s.empty:
                    continue
                st = METHOD_STYLE[method]
                y = np.clip(s[ycol].values, GAP_FLOOR, None)
                ax.semilogy(s.iteration.values, y, color=st["color"],
                            ls=st["ls"], label=st["label"], lw=1.6)
            ax.set_xlim(0, 60)
            ax.grid(True, which="both", alpha=0.25)
            ax.set_ylabel(ylab)
            if i == len(WFR_TARGET_ORDER) - 1:
                ax.set_xlabel("iteration")
            if col == 0:
                ax.set_title(TARGET_TITLE[tname] + r"  ($\Lambda=1000,\varepsilon=10^{-3}$)",
                             loc="left", fontsize=10)
    axes[0][0].legend(fontsize=8, loc="best")
    fig.suptitle("Phase-separation mechanism: the covariance, not the mean, is the "
                 "bottleneck", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    return fig


def _wfr_fig_schedule_sweep(sweep_df):
    """Hitting iterations to gap < 1e-3 vs the fixed fraction ``c``.

    One column per target; one curve per ``(Lambda, eps)``. ``c = 0`` is FR-only;
    the theorem-bound schedule ``h = mu_min`` is overlaid as a star at its own
    ``c = mu_min beta``. More transport monotonically reduces the iteration count
    -- each WFR iteration blends a transport move and a natural-gradient move --
    until it saturates, while the theorem bound is stranded in the interior.
    """
    fig, axes = plt.subplots(1, len(WFR_TARGET_ORDER), figsize=(5.4 * len(WFR_TARGET_ORDER), 3.8),
                             squeeze=False)
    ycol = "iter_to_1e_minus_3"
    ylab = r"iterations to gap $<10^{-3}$"
    for j, target_name in enumerate(WFR_TARGET_ORDER):
        sub = sweep_df[sweep_df.target_name == target_name]
        combos = list(sub[["Lambda", "epsilon"]].drop_duplicates()
                      .sort_values(["Lambda", "epsilon"]).itertuples(index=False))
        cmap = plt.cm.plasma(np.linspace(0.05, 0.8, 4))
        ax = axes[0, j]
        for k, (Lambda, eps) in enumerate(combos):
            fx = sub[(sub.schedule_kind == "fixed_c")
                     & (np.isclose(sub.Lambda, Lambda)) & (np.isclose(sub.epsilon, eps))]
            fx = fx.sort_values("c")
            yv = fx[ycol].replace(-1, np.nan)
            lbl = rf"$\Lambda={Lambda:g},\varepsilon={eps:g}$"
            ax.plot(fx.c, yv, "o-", color=cmap[k % 4], lw=1.5, ms=4, label=lbl)
            th = sub[(sub.schedule_kind == "theory_mu_min")
                     & (np.isclose(sub.Lambda, Lambda)) & (np.isclose(sub.epsilon, eps))]
            if not th.empty:
                yt = th[ycol].replace(-1, np.nan).iloc[0]
                ax.plot(th.c.iloc[0], yt, "*", color=cmap[k % 4], ms=13,
                        markeredgecolor="k", markeredgewidth=0.5, zorder=5)
        ax.set_ylabel(ylab)
        ax.grid(True, which="both", alpha=0.25)
        ax.set_title(TARGET_TITLE[target_name])
        ax.set_xlabel(r"fixed transport fraction $c$  ($h=c/\beta$)")
    axes[0, 0].legend(fontsize=8, loc="best")
    fig.suptitle(r"WFR schedule sweep: transport cuts iterations to convergence; "
                 r"$\star=$ theory $h=\mu_{\min}$", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def _wfr_fig_dt_heatmap(dt_df):
    """Per-target ``Delta t`` x ``c`` heatmap of hitting iterations to gap<1e-3.

    One panel per ``(target, Lambda, eps)``; cells show log10 iterations-to-1e-3
    (white = never reached within budget / not SPD-feasible). Demonstrates the
    joint robustness of the WFR-fixed scheme across the discretization step
    ``Delta t`` and the transport fraction ``c``.
    """
    keys = (dt_df[["target_name", "Lambda", "epsilon"]].drop_duplicates()
            .sort_values(["target_name", "Lambda", "epsilon"]).itertuples(index=False))
    keys = list(keys)
    n = len(keys)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.5 * ncol, 3.1 * nrow), squeeze=False)
    dts = sorted(dt_df.dt.unique())
    cs = sorted(dt_df.c.unique())
    for idx, (tname, Lambda, eps) in enumerate(keys):
        ax = axes[idx // ncol][idx % ncol]
        cell = dt_df[(dt_df.target_name == tname)
                     & (np.isclose(dt_df.Lambda, Lambda)) & (np.isclose(dt_df.epsilon, eps))]
        Z = np.full((len(cs), len(dts)), np.nan)
        for _, row in cell.iterrows():
            i, k = cs.index(row["c"]), dts.index(row["dt"])
            b = row.batches_to_1e_minus_3
            if b is not None and b >= 0 and row.spd_feasible == 1:
                # c=0 is fr_only (1 batch/iter); c>0 is the WFR splitting (2/iter).
                iters = b if row["c"] == 0 else b / 2.0
                Z[i, k] = np.log10(iters)
        im = ax.imshow(Z, origin="lower", aspect="auto", cmap="viridis_r")
        ax.set_xticks(range(len(dts))); ax.set_xticklabels([f"{d:g}" for d in dts], fontsize=7)
        ax.set_yticks(range(len(cs))); ax.set_yticklabels([f"{c:g}" for c in cs], fontsize=7)
        ax.set_title(rf"{TARGET_TITLE[tname]}" "\n" rf"$\Lambda={Lambda:g},\varepsilon={eps:g}$", fontsize=8)
        if idx % ncol == 0:
            ax.set_ylabel(r"$c$")
        if idx // ncol == nrow - 1:
            ax.set_xlabel(r"$\Delta t$")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=r"$\log_{10}$ iters")
    for idx in range(n, nrow * ncol):
        axes[idx // ncol][idx % ncol].axis("off")
    fig.suptitle(r"WFR-fixed: iterations to gap $<10^{-3}$ over $\Delta t\times c$", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    return fig


def _tab_wfr_hitting(summary_df):
    """Main hitting-iteration table: iterations to gap thresholds + final gap, per run."""
    lines = [
        r"\begin{tabular}{llrrrrrr}", r"\toprule",
        r"target & method & $\Lambda$ & $\varepsilon$ & "
        r"$n_{10^{-1}}$ & $n_{10^{-3}}$ & $n_{10^{-6}}$ & gap$_{\mathrm{final}}$ \\",
        r"\midrule",
    ]
    df = summary_df.sort_values(["target_name", "Lambda", "epsilon", "method"])
    prev = None
    for _, r in df.iterrows():
        tgt = WFR_TARGET_TEX.get(r.target_name, r.target_name)
        key = (r.target_name, r.Lambda, r.epsilon)
        if prev is not None and key != prev:
            lines.append(r"\midrule")
        show_tgt = tgt if key != prev else ""
        lines.append(
            f"{show_tgt} & {WFR_METHOD_TEX.get(r.method, r.method)} & "
            f"{r.Lambda:g} & {r.epsilon:g} & "
            f"{_wfr_fmt_int(r.iter_to_1e_minus_1)} & "
            f"{_wfr_fmt_int(r.iter_to_1e_minus_3)} & "
            f"{_wfr_fmt_int(r.iter_to_1e_minus_6)} & {_wfr_fmt_sci(r.gap_final)} \\\\")
        prev = key
    lines += [r"\bottomrule", r"\end{tabular}"]
    _write_table("tab_wfr_hitting.tex", "\n".join(lines))


def _tab_wfr_iters(summary_df):
    """Iterations to gap < 1e-3, per run, at the harder start eps=1e-3.

    The headline comparison: the WFR splitting reaches the tolerance in the
    fewest iterations (about half of pure Fisher--Rao), W-only is far behind.
    """
    lines = [
        r"\begin{tabular}{ll rr}", r"\toprule",
        r" & & \multicolumn{2}{c}{iterations to gap $<10^{-3}$} \\",
        r"\cmidrule(lr){3-4}",
        r"target & method & $\Lambda{=}100$ & $\Lambda{=}1000$ \\", r"\midrule",
    ]
    # Use eps = 1e-3 (the harder, more underdispersed start) for this summary.
    df = summary_df[np.isclose(summary_df.epsilon, 1e-3)]
    for ti, tname in enumerate(WFR_TARGET_ORDER):
        if ti:
            lines.append(r"\midrule")
        for mi, method in enumerate(METHOD_ORDER):
            sub = df[(df.target_name == tname) & (df.method == method)]
            def cell(col, Lam):
                r = sub[np.isclose(sub.Lambda, Lam)]
                return _wfr_fmt_int(r[col].iloc[0]) if len(r) else "--"
            show = WFR_TARGET_TEX.get(tname, tname) if mi == 0 else ""
            lines.append(
                f"{show} & {WFR_METHOD_TEX.get(method, method)} & "
                f"{cell('iter_to_1e_minus_3', 100)} & {cell('iter_to_1e_minus_3', 1000)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _write_table("tab_wfr_iters.tex", "\n".join(lines))


def _tab_wfr_sweep(sweep_df):
    """Schedule-sweep table: fixed c, h=c/beta, iterations to 1e-3, final gap."""
    lines = [
        r"\begin{tabular}{llrrrr}", r"\toprule",
        r"target & $(\Lambda,\varepsilon)$ & $c$ & $h=c/\beta$ & "
        r"$n_{10^{-3}}$ & gap$_{\mathrm{final}}$ \\", r"\midrule",
    ]
    df = sweep_df[sweep_df.schedule_kind == "fixed_c"].sort_values(
        ["target_name", "Lambda", "epsilon", "c"])
    prev = None
    for _, r in df.iterrows():
        key = (r.target_name, r.Lambda, r.epsilon)
        if prev is not None and key != prev:
            lines.append(r"\midrule")
        show = (WFR_TARGET_TEX.get(r.target_name, r.target_name) if key != prev else "")
        combo = rf"$({r.Lambda:g},{r.epsilon:g})$" if key != prev else ""
        lines.append(
            f"{show} & {combo} & {r.c:g} & {_fmt(r.h, 4)} & "
            f"{_wfr_fmt_int(r.iter_to_1e_minus_3)} & {_wfr_fmt_sci(r.gap_final)} \\\\")
        prev = key
    lines += [r"\bottomrule", r"\end{tabular}"]
    _write_table("tab_wfr_sweep.tex", "\n".join(lines))


def _tab_wfr_metadata(metadata):
    """Target metadata: alpha, beta, condition number beta/alpha, eps, gap0."""
    lines = [
        r"\begin{tabular}{llrrrrr}", r"\toprule",
        r"target & $\Lambda$ & $\varepsilon$ & $\alpha$ & $\beta$ & "
        r"$\beta/\alpha$ & gap$_0$ \\", r"\midrule",
    ]
    rows = []
    for key, md in metadata["targets"].items():
        rows.append((md["target_name"], md["Lambda"], md["epsilon"],
                     md["alpha"], md["beta"], md.get("gap0", float("nan"))))
    for tname, Lambda, eps, alpha, beta, gap0 in sorted(rows):
        cond = beta / alpha if alpha else float("nan")
        lines.append(
            f"{WFR_TARGET_TEX.get(tname, tname)} & {Lambda:g} & {eps:g} & "
            f"{_fmt(alpha, 3)} & {_fmt(beta, 3)} & {_fmt(cond, 1)} & "
            f"{_wfr_fmt_sci(gap0)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _write_table("tab_wfr_metadata.tex", "\n".join(lines))


def build_wfr_assets():
    """Figures + tables for the WFR Gaussian gradient-flow report.

    Reads the final CSVs in ``outputs/wfr_gradient_flow`` and writes figures into
    ``reports/assets/figs`` and table fragments into ``reports/assets``. Raises
    ``FileNotFoundError`` (via pandas) if the required CSVs are absent, so the
    top-level driver can skip the group cleanly on a checkout without outputs.
    """
    print("WFR Gaussian gradient flow:")
    import json
    long_df = pd.read_csv(os.path.join(WFR_DIR, "results_long.csv"))
    summary_df = pd.read_csv(os.path.join(WFR_DIR, "summary.csv"))
    sweep_df = pd.read_csv(os.path.join(WFR_DIR, "schedule_sweep.csv"))
    with open(os.path.join(WFR_DIR, "target_metadata.json")) as fh:
        metadata = json.load(fh)

    for tname in WFR_TARGET_ORDER:
        _savefig(_wfr_fig_phase_separation(long_df, tname),
                 f"wfr_phase_separation_{tname}")
    _savefig(_wfr_fig_adaptive_schedule(long_df, summary_df), "wfr_adaptive_schedule")
    _savefig(_wfr_fig_error_decomposition(long_df), "wfr_error_decomposition")
    _savefig(_wfr_fig_schedule_sweep(sweep_df), "wfr_schedule_sweep")

    dt_path = os.path.join(WFR_DIR, "dt_sweep.csv")
    if os.path.exists(dt_path):
        dt_df = pd.read_csv(dt_path)
        if not dt_df.empty:
            _savefig(_wfr_fig_dt_heatmap(dt_df), "wfr_dt_heatmap")

    _tab_wfr_hitting(summary_df)
    _tab_wfr_iters(summary_df)
    _tab_wfr_sweep(sweep_df)
    _tab_wfr_metadata(metadata)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=_ROOT)
    p.add_argument("--only",
                   choices=["omega_tau", "local_rate", "discretization", "wfr"],
                   default=None, help="Build only one group's assets.")
    args = p.parse_args()
    os.makedirs(FIGS, exist_ok=True)
    apply_style()
    builders = {
        "omega_tau": build_omega_tau_assets,
        "local_rate": build_local_rate_assets,
        "discretization": build_discretization_assets,
        "wfr": build_wfr_assets,
    }
    selected = [args.only] if args.only else list(builders)
    failures = []
    for name in selected:
        try:
            builders[name]()
        except FileNotFoundError as e:
            # A group whose final outputs are absent (e.g. gitignored config.json
            # not regenerated) should not block the other groups' assets.
            print(f"  [skip] {name}: missing input ({e})")
            failures.append(name)
    print(f"\nAll report assets written under {ASSETS}")
    if failures:
        print(f"Skipped (missing inputs): {', '.join(failures)}")


if __name__ == "__main__":
    main()
