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

# ---------------------------------------------------------------------------
# Paths and shared constants
# ---------------------------------------------------------------------------
ASSETS = os.path.join(_HERE, "assets")
FIGS = os.path.join(ASSETS, "figs")

GAUSS_DIR = os.path.join(_ROOT, "outputs", "gaussian_grid")
LOGC_DIR = os.path.join(_ROOT, "outputs", "logconcave_grid")
LR_DIR = os.path.join(_ROOT, "outputs", "natural_gradient_local_rate")

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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=_ROOT)
    p.parse_args()
    os.makedirs(FIGS, exist_ok=True)
    apply_style()
    build_omega_tau_assets()
    build_local_rate_assets()
    print(f"\nAll report assets written under {ASSETS}")


if __name__ == "__main__":
    main()
