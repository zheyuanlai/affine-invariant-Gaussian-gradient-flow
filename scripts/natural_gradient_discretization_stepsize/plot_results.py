"""Generate the discretization-stepsize figures from the committed CSVs.

Reads ``results_long.csv``, ``summary.csv``, ``stepsize_summary.csv``,
``scalar_diagnostic.csv`` and ``target_metadata.json`` from ``--outdir`` and
writes paired PDF/PNG figures into ``<outdir>/figures/``. No dynamics are
re-run; this is pure post-processing.

Figures
-------
1-3  energy gap vs elapsed time n*dt (Gaussian / literature / smooth), 3 lambda panels
4-6  stepsize stability heatmaps (one figure per target, Riemannian + KL panels)
7    theory vs empirical stepsize (Gaussian + smooth; the main proof-artifact figure)
8    scalar covariance diagnostic
9    wall-clock / time-to-tolerance summary

Usage::

    python scripts/natural_gradient_discretization_stepsize/plot_results.py \
        --outdir outputs/natural_gradient_discretization_stepsize
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.plotting_style import apply_style, save_figure
from src.natural_gradient_discretization_stepsize.plotting import (
    METHOD_STYLE, LAMBDA_ORDER, REPRESENTATIVE_DT, REPRESENTATIVE_DT_BY_TARGET,
    CONVERGENCE_DT_BY_TARGET, semilogy_clipped, classify_level,
    stability_heatmap, add_class_legend,
)

TARGET_TITLE = {
    "gaussian_posterior": "Gaussian posterior",
    "literature_logconcave": "literature log-concave (quartic)",
    "smooth_logconcave": "smooth strongly log-concave",
}


def _savefig(fig, figs_dir, name):
    paths = save_figure(fig, os.path.join(figs_dir, name))
    plt.close(fig)
    print(f"  fig  {name}")
    return paths


# ---------------------------------------------------------------------------
# Figures 1-3: energy gap vs elapsed time n*dt
# ---------------------------------------------------------------------------

def fig_energy_gap(long_df, target_name, figs_dir, meta=None, summary=None):
    """Energy gap vs elapsed time n*dt at representative *convergent* stepsizes.

    The representative stepsizes are chosen per target so that they converge
    for every lambda (see ``REPRESENTATIVE_DT_BY_TARGET``); the quartic target
    has unbounded curvature and therefore much smaller stable stepsizes than the
    globally smooth targets. As a safety net, any run flagged non-stable in
    ``summary`` is skipped so a divergent stepsize can never reach the figure.
    """
    sub = long_df[long_df.target_name == target_name]
    lambdas = [l for l in LAMBDA_ORDER if l in set(sub["lambda"].unique())]
    fig, axes = plt.subplots(1, len(lambdas), figsize=(4.4 * len(lambdas), 3.4),
                             sharey=True)
    axes = np.atleast_1d(axes)
    rep_pool = REPRESENTATIVE_DT_BY_TARGET.get(target_name, REPRESENTATIVE_DT)
    rep_dt = [d for d in rep_pool if d in set(sub["dt"].unique())]
    alphas = np.linspace(0.45, 1.0, len(rep_dt))

    def _is_stable(lam, method, dt):
        if summary is None:
            return True
        r = summary[(summary.target_name == target_name)
                    & (np.isclose(summary["lambda"], lam))
                    & (summary.method == method) & (np.isclose(summary.dt, dt))]
        return bool(r.stable.iloc[0]) if len(r) else True

    for ax, lam in zip(axes, lambdas):
        for method in ("riemannian", "kl"):
            st = METHOD_STYLE[method]
            for a, dt in zip(alphas, rep_dt):
                if not _is_stable(lam, method, dt):
                    continue
                run = sub[(sub["lambda"] == lam) & (sub.method == method)
                          & (np.isclose(sub.dt, dt))].sort_values("n")
                if run.empty:
                    continue
                semilogy_clipped(ax, run.t.values, run.energy_gap.values,
                                 color=st["color"], ls=st["ls"], alpha=float(a),
                                 lw=1.5, label=f"{st['label']} dt={dt:g}")
        ax.set_title(rf"$\lambda={lam:g}$")
        ax.set_xlabel(r"$n\Delta t$")
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel(r"energy gap $\mathcal{E}(a_n)-\mathcal{E}(a_\star)$")
    h, l = axes[-1].get_legend_handles_labels()
    axes[-1].legend(h, l, fontsize=7, loc="upper right", ncol=1)
    fig.suptitle(f"{TARGET_TITLE[target_name]}: energy gap vs time "
                 r"(convergent $\Delta t$)", y=1.02)
    fig.tight_layout()
    return _savefig(fig, figs_dir, f"fig_gap_{target_name}")


# ---------------------------------------------------------------------------
# Figures 4-6: stepsize stability heatmaps
# ---------------------------------------------------------------------------

def _class_matrix(summary, target_name, method, dts, lambdas):
    mat = np.zeros((len(lambdas), len(dts)), dtype=float)
    for i, lam in enumerate(lambdas):
        for j, dt in enumerate(dts):
            r = summary[(summary.target_name == target_name) & (summary.method == method)
                        & (np.isclose(summary["lambda"], lam)) & (np.isclose(summary.dt, dt))]
            mat[i, j] = classify_level(int(r.stable.iloc[0]), int(r.monotone.iloc[0])) \
                if len(r) else np.nan
    return mat


def fig_stability_heatmap(summary, step_df, target_name, figs_dir):
    sub = summary[summary.target_name == target_name]
    dts = sorted(sub.dt.unique())
    lambdas = sorted(sub["lambda"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.2), sharey=True)
    for ax, method in zip(axes, ("riemannian", "kl")):
        mat = _class_matrix(sub, target_name, method, dts, lambdas)
        th = step_df[(step_df.target_name == target_name) & (step_df.method == method)]
        dt_theory = None
        if len(th) and int(th.theory_bound_available.iloc[0]) == 1:
            dt_theory = float(th.dt_theory_for_method.iloc[0])
        stability_heatmap(ax, mat, dts, lambdas,
                          title=f"{METHOD_STYLE[method]['label']}", dt_theory=dt_theory)
    fig.suptitle(f"{TARGET_TITLE[target_name]}: stepsize stability"
                 + (r"  (dotted: theoretical $\Delta t$)"
                    if target_name != "literature_logconcave" else ""), y=1.05)
    add_class_legend(fig)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return _savefig(fig, figs_dir, f"fig_heatmap_{target_name}")


# ---------------------------------------------------------------------------
# Figure 7: theory vs empirical stepsize (the proof-artifact figure)
# ---------------------------------------------------------------------------

def fig_theory_vs_empirical(step_df, figs_dir):
    targets = ["gaussian_posterior", "smooth_logconcave"]
    fig, axes = plt.subplots(1, len(targets), figsize=(5.2 * len(targets), 3.8),
                             sharey=True)
    axes = np.atleast_1d(axes)
    for ax, tname in zip(axes, targets):
        sub = step_df[(step_df.target_name == tname)
                      & (step_df.theory_bound_available == 1)]
        lambdas = sorted(sub["lambda"].unique())
        x = np.arange(len(lambdas))
        width = 0.38
        for k, method in enumerate(("riemannian", "kl")):
            st = METHOD_STYLE[method]
            ms = sub[sub.method == method].set_index("lambda")
            theory = [ms.loc[l, "dt_theory_for_method"] for l in lambdas]
            mono = [ms.loc[l, "dt_max_monotone"] for l in lambdas]
            stab = [ms.loc[l, "dt_max_stable"] for l in lambdas]
            off = (k - 0.5) * width
            ax.bar(x + off, theory, width * 0.9, color=st["color"], alpha=0.35,
                   label=f"{st['label']} theory", edgecolor=st["color"])
            ax.scatter(x + off, mono, marker="o", color=st["color"], zorder=5,
                       s=42, label=f"{st['label']} empirical monotone")
            ax.scatter(x + off, stab, marker="^", facecolor="white", zorder=5,
                       edgecolor=st["color"], s=42, label=f"{st['label']} empirical stable")
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{l:g}" for l in lambdas])
        ax.set_xlabel(r"$\lambda$")
        ax.set_title(TARGET_TITLE[tname])
        ax.grid(True, which="both", axis="y", alpha=0.25)
    axes[0].set_ylabel(r"stepsize $\Delta t$")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=7.5, frameon=False,
               bbox_to_anchor=(0.5, -0.08))
    fig.suptitle(r"Theoretical vs empirical maximum stepsize", y=1.03)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return _savefig(fig, figs_dir, "fig_theory_vs_empirical")


# ---------------------------------------------------------------------------
# Figure 8: scalar covariance diagnostic
# ---------------------------------------------------------------------------

def fig_scalar_diagnostic(scalar_df, figs_dir):
    cov = scalar_df[scalar_df.experiment == "covariance"]
    C0s = sorted(cov.C0.unique())
    dts = sorted(cov.dt.unique())
    fig, axes = plt.subplots(1, len(C0s), figsize=(3.2 * len(C0s), 3.2), sharey=False)
    axes = np.atleast_1d(axes)
    cmap = plt.get_cmap("viridis")
    for ax, C0 in zip(axes, C0s):
        for k, dt in enumerate(dts):
            shade = cmap(0.15 + 0.7 * k / max(1, len(dts) - 1))
            for method in ("riemannian", "kl"):
                st = METHOD_STYLE[method]
                run = cov[(np.isclose(cov.C0, C0)) & (np.isclose(cov.dt, dt))
                          & (cov.method == method)].sort_values("n")
                if run.empty:
                    continue
                ax.plot(run.t.values, run.C.values, color=shade, ls=st["ls"], lw=1.4,
                        label=(f"dt={dt:g}" if method == "riemannian" else None))
        ax.axhline(1.0, color="0.5", lw=1.0, ls=":")
        ax.set_yscale("log")
        ax.set_title(rf"$C_0={C0:g}$")
        ax.set_xlabel(r"$n\Delta t$")
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel(r"covariance $C_n$")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=len(dts), fontsize=7.5, frameon=False,
               bbox_to_anchor=(0.5, -0.06))
    fig.suptitle(r"Scalar $\mathcal{N}(0,1)$ covariance diagnostic "
                 r"(solid: Riemannian, dashed: KL)", y=1.04)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return _savefig(fig, figs_dir, "fig_scalar_covariance")


# ---------------------------------------------------------------------------
# Figure 9: time-to-tolerance / wall-clock summary
# ---------------------------------------------------------------------------

def fig_time_to_tolerance(summary, figs_dir):
    """Iterations-to-tolerance for stable runs (dimension is small, so wall-clock
    differences are noisy; iterations-to-tolerance is the meaningful axis)."""
    sub = summary[(summary.stable == 1) & np.isfinite(summary.iter_to_gap_1e_minus_4)
                  & (summary.iter_to_gap_1e_minus_4 >= 0)]
    targets = sorted(sub.target_name.unique())
    fig, axes = plt.subplots(1, len(targets), figsize=(4.0 * len(targets), 3.2),
                             sharey=True)
    axes = np.atleast_1d(axes)
    for ax, tname in zip(axes, targets):
        ts = sub[sub.target_name == tname]
        for method in ("riemannian", "kl"):
            st = METHOD_STYLE[method]
            mm = ts[ts.method == method].sort_values("dt")
            if mm.empty:
                continue
            ax.loglog(mm.dt.values, mm.iter_to_gap_1e_minus_4.values, marker="o",
                      color=st["color"], ls=st["ls"], label=st["label"], markersize=4)
        ax.set_title(TARGET_TITLE.get(tname, tname), fontsize=9)
        ax.set_xlabel(r"$\Delta t$")
        ax.grid(True, which="both", alpha=0.25)
    axes[0].set_ylabel(r"iterations to gap $\leq 10^{-4}$")
    axes[-1].legend(fontsize=8)
    fig.suptitle(r"Iterations to tolerance for stable runs", y=1.03)
    fig.tight_layout()
    return _savefig(fig, figs_dir, "fig_time_to_tolerance")


# ---------------------------------------------------------------------------
# Appendix figure A1: comprehensive convergence-speed grid
# ---------------------------------------------------------------------------

APPENDIX_TARGETS = ["gaussian_posterior", "smooth_logconcave", "literature_logconcave"]


def fig_convergence_speed_grid(long_df, summary, figs_dir):
    """One per-target figure each: rows = lambda, columns = stepsize.

    Every panel shows exactly two curves -- Riemannian (solid) vs KL (dashed) --
    at one (lambda, dt) cell, so the matched-stepsize method comparison is
    immediate and uncluttered. Stepsizes are the curated convergent set for the
    target (``CONVERGENCE_DT_BY_TARGET``), stable for both methods at every
    lambda. Returns the list of saved-path tuples (one per target).
    """
    out = []
    for tname in APPENDIX_TARGETS:
        out.append(_fig_convergence_speed_target(long_df, summary, tname, figs_dir))
    return out


def _fig_convergence_speed_target(long_df, summary, target_name, figs_dir):
    sub = long_df[long_df.target_name == target_name]
    lambdas = [l for l in LAMBDA_ORDER if l in set(sub["lambda"].unique())]
    dts = [d for d in CONVERGENCE_DT_BY_TARGET[target_name]
           if d in set(sub["dt"].unique())]
    nrow, ncol = len(lambdas), len(dts)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.9 * ncol, 2.5 * nrow),
                             sharex=False, sharey="row", squeeze=False)
    for i, lam in enumerate(lambdas):
        for j, dt in enumerate(dts):
            ax = axes[i, j]
            for method in ("riemannian", "kl"):
                st = METHOD_STYLE[method]
                r = summary[(summary.target_name == target_name)
                            & (np.isclose(summary["lambda"], lam))
                            & (summary.method == method) & (np.isclose(summary.dt, dt))]
                if len(r) and not bool(r.stable.iloc[0]):
                    continue
                run = sub[(np.isclose(sub["lambda"], lam)) & (sub.method == method)
                          & (np.isclose(sub.dt, dt))].sort_values("n")
                if run.empty:
                    continue
                semilogy_clipped(ax, run.t.values, run.energy_gap.values,
                                 color=st["color"], ls=st["ls"], lw=1.6,
                                 label=st["label"])
            if i == 0:
                ax.set_title(rf"$\Delta t={dt:g}$")
            if j == 0:
                ax.set_ylabel(rf"$\lambda={lam:g}$" + "\n" + r"gap", fontsize=8.5)
            if i == nrow - 1:
                ax.set_xlabel(r"$n\Delta t$")
            ax.grid(True, which="both", alpha=0.22)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=2, fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(f"{TARGET_TITLE[target_name]}: convergence at matched stepsize "
                 r"(Riemannian vs KL)", y=1.0)
    fig.tight_layout(rect=(0, 0.05, 1, 0.98))
    return _savefig(fig, figs_dir, f"fig_convergence_speed_{target_name}")


# ---------------------------------------------------------------------------
# Cost of a stepsize (iterations & elapsed n*dt to tolerance)
# ---------------------------------------------------------------------------

def fig_iterations_vs_stepsize(summary, figs_dir):
    """Iterations- and physical-time-to-tolerance vs stepsize for stable runs.

    Top row: iterations to gap <= 1e-6 vs dt (log-log), one line per
    (target, method), averaged over lambda. The slope is approximately -1
    (iterations ~ 1/dt). Bottom row: the same in elapsed time n*dt, which is nearly
    flat in dt -- a larger convergent stepsize costs proportionally fewer
    iterations for the same progress, with no penalty in elapsed n*dt.
    """
    sub = summary[(summary.stable == 1)
                  & (summary.iter_to_gap_1e_minus_6 >= 0)
                  & np.isfinite(summary.time_to_gap_1e_minus_6)]
    fig, axes = plt.subplots(2, len(APPENDIX_TARGETS),
                             figsize=(4.0 * len(APPENDIX_TARGETS), 6.0),
                             sharex=True)
    axes = np.atleast_2d(axes)
    for j, tname in enumerate(APPENDIX_TARGETS):
        ts = sub[sub.target_name == tname]
        ax_it, ax_t = axes[0, j], axes[1, j]
        for method in ("riemannian", "kl"):
            st = METHOD_STYLE[method]
            mm = ts[ts.method == method]
            if mm.empty:
                continue
            it = mm.groupby("dt")["iter_to_gap_1e_minus_6"].mean().sort_index()
            tt = mm.groupby("dt")["time_to_gap_1e_minus_6"].mean().sort_index()
            ax_it.loglog(it.index.values, it.values, marker="o", ms=4,
                         color=st["color"], ls=st["ls"], label=st["label"])
            ax_t.semilogx(tt.index.values, tt.values, marker="o", ms=4,
                          color=st["color"], ls=st["ls"], label=st["label"])
        # 1/dt reference slope on the iteration panel.
        if not ts.empty:
            dref = np.array(sorted(ts.dt.unique()), dtype=float)
            base = ts[ts.method == "kl"]["iter_to_gap_1e_minus_6"].max()
            if np.isfinite(base) and base > 0:
                ax_it.loglog(dref, base * dref.min() / dref, color="0.6",
                             ls=":", lw=1.0, label=r"$\propto 1/\Delta t$")
        ax_it.set_title(TARGET_TITLE[tname], fontsize=9)
        ax_it.grid(True, which="both", alpha=0.25)
        ax_t.grid(True, which="both", alpha=0.25)
        ax_t.set_xlabel(r"$\Delta t$")
    axes[0, 0].set_ylabel(r"iterations to gap $\leq 10^{-6}$")
    axes[1, 0].set_ylabel(r"elapsed $n\Delta t$ to gap $\leq 10^{-6}$")
    axes[0, -1].legend(fontsize=7.5, loc="upper right")
    fig.suptitle(r"Cost of a stepsize: iterations $\propto 1/\Delta t$, "
                 r"elapsed $n\Delta t$ roughly flat", y=1.01)
    fig.tight_layout()
    return _savefig(fig, figs_dir, "fig_iterations_vs_stepsize")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def build_all(outdir):
    figs_dir = os.path.join(outdir, "figures")
    os.makedirs(figs_dir, exist_ok=True)
    apply_style()
    long_df = pd.read_csv(os.path.join(outdir, "results_long.csv"))
    summary = pd.read_csv(os.path.join(outdir, "summary.csv"))
    step_df = pd.read_csv(os.path.join(outdir, "stepsize_summary.csv"))
    scalar_df = pd.read_csv(os.path.join(outdir, "scalar_diagnostic.csv"))
    with open(os.path.join(outdir, "target_metadata.json")) as fh:
        meta = json.load(fh)

    for tname in ("gaussian_posterior", "literature_logconcave", "smooth_logconcave"):
        if tname in set(long_df.target_name.unique()):
            fig_energy_gap(long_df, tname, figs_dir, meta, summary=summary)
            fig_stability_heatmap(summary, step_df, tname, figs_dir)
    fig_theory_vs_empirical(step_df, figs_dir)
    fig_scalar_diagnostic(scalar_df, figs_dir)
    fig_time_to_tolerance(summary, figs_dir)
    # Appendix: comprehensive convergence-speed study (convergent dt only).
    fig_convergence_speed_grid(long_df, summary, figs_dir)
    fig_iterations_vs_stepsize(summary, figs_dir)
    print(f"\nFigures written under {figs_dir}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default="outputs/natural_gradient_discretization_stepsize")
    args = p.parse_args()
    build_all(args.outdir)


if __name__ == "__main__":
    main()
