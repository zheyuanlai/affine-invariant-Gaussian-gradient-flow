"""Plotting and table helpers for the nonconvex-instability report."""
from __future__ import annotations

import os

import numpy as np
import matplotlib.pyplot as plt

from src.common.plotting_style import save_figure


R_COLORS = {
    20.0: "#1f77b4",
    50.0: "#2ca02c",
    100.0: "#ff7f0e",
    300.0: "#d62728",
    1000.0: "#9467bd",
}


def _color_for_R(R):
    return R_COLORS.get(float(R), plt.get_cmap("viridis")(0.5))


def _savefig(fig, figs_dir, name):
    paths = save_figure(fig, os.path.join(figs_dir, name))
    plt.close(fig)
    return paths


def _positive_finite(y):
    y = np.asarray(y, dtype=float)
    return np.where(np.isfinite(y) & (y > 0.0), y, np.nan)


def fig_riemannian_cascade(long_df, figs_dir):
    sub = long_df[long_df.experiment == "riemannian_cascade"]
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for R in sorted(sub.R.dropna().unique()):
        run = sub[np.isclose(sub.R, R)].sort_values("n")
        ax.plot(
            run.n.values,
            _positive_finite(run.fr_grad_sq.values),
            marker="o",
            color=_color_for_R(R),
            label=rf"$R={R:g}$",
        )
    ax.set_yscale("log")
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"FR gradient squared $(1-cA_R(c))^2$")
    ax.set_title("Riemannian Fisher--Rao cascade")
    ax.legend(ncol=2)
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    return _savefig(fig, figs_dir, "fig_nonconvex_riemannian_cascade")


def fig_kl_pole(kl_df, figs_dir):
    df = kl_df.sort_values("epsilon", ascending=False)
    eps = df.epsilon.values.astype(float)
    y = _positive_finite(df.fr_grad_sq_after_first.values)
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.loglog(eps, y, "o-", color="#d62728", label="after first KL step")
    ref = 4.0 / (eps ** 2)
    ax.loglog(eps, ref, "k--", lw=1.2, label=r"$4\varepsilon^{-2}$ reference")
    ax.invert_xaxis()
    ax.set_xlabel(r"$\varepsilon$ in $c_0=1-\varepsilon$")
    ax.set_ylabel(r"FR gradient squared at $c_1$")
    ax.set_title("KL/Bregman pole near loss of SPD")
    ax.legend()
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    return _savefig(fig, figs_dir, "fig_nonconvex_kl_pole")


def fig_wasserstein_bound(bw_df, figs_dir):
    fig, ax = plt.subplots(figsize=(5.8, 3.9))
    for R in sorted(bw_df.R.dropna().unique()):
        run = bw_df[np.isclose(bw_df.R, R)].sort_values("N")
        color = _color_for_R(R)
        ax.loglog(
            run.N.values,
            _positive_finite(run.running_min_bw_grad_sq.values),
            color=color,
            lw=1.6,
            label=rf"$R={R:g}$ run. min",
        )
        ax.loglog(
            run.N.values,
            _positive_finite(run.theorem_envelope.values),
            color=color,
            ls="--",
            lw=1.1,
            alpha=0.75,
        )
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"BW gradient squared")
    ax.set_title(r"Running minimum vs $150\Delta/(\eta N)$ envelope")
    ax.legend(ncol=2, fontsize=8)
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    return _savefig(fig, figs_dir, "fig_nonconvex_wasserstein_bound")


def fig_curvature_diagnostic(long_df, figs_dir):
    sub = long_df[long_df.experiment == "riemannian_cascade"]
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    for R in sorted(sub.R.dropna().unique()):
        run = sub[np.isclose(sub.R, R)].sort_values("n")
        ax.plot(
            run.n.values,
            run.A.values,
            marker="o",
            color=_color_for_R(R),
            label=rf"$R={R:g}$",
        )
    ax.axhline(-1.0, color="0.4", ls=":", lw=1.0)
    ax.axhline(0.0, color="0.55", ls="--", lw=1.0)
    ax.axhline(1.0, color="0.4", ls=":", lw=1.0)
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$A_R(c_n)=\mathbb{E}[V_R''(X_n)]$")
    ax.set_title("Curvature sampled by unstable trajectories")
    ax.legend(ncol=2)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    return _savefig(fig, figs_dir, "fig_nonconvex_curvature_diagnostic")


def build_all_figures(long_df, kl_df, bw_df, figs_dir):
    os.makedirs(figs_dir, exist_ok=True)
    fig_riemannian_cascade(long_df, figs_dir)
    fig_kl_pole(kl_df, figs_dir)
    fig_wasserstein_bound(bw_df, figs_dir)
    fig_curvature_diagnostic(long_df, figs_dir)


def _fmt_sci(x):
    if x is None or not np.isfinite(float(x)):
        return "--"
    x = float(x)
    if x == 0.0:
        return "0"
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / (10.0 ** exp)
    if -2 <= exp <= 2:
        return f"{x:.3g}"
    return rf"${mant:.2g}\times10^{{{exp}}}$"


def _fmt_bool(x):
    if isinstance(x, str):
        return x
    return "yes" if bool(x) else "no"


def _tex_text(x):
    text = str(x)
    return text.replace("\\", r"\textbackslash{}").replace("_", r"\_")


def summary_table_tex(summary_df, kl_df, bw_df):
    """Manuscript-ready compact summary table."""
    riem = summary_df[summary_df.experiment == "riemannian_cascade"].sort_values("R")
    bw_last = (bw_df.sort_values("N").groupby("R", as_index=False).tail(1)
               .sort_values("R"))
    slope_col = ("fit_slope_pole_window_log_fr_vs_log_epsilon"
                 if "fit_slope_pole_window_log_fr_vs_log_epsilon" in kl_df
                 else "fit_slope_log_fr_vs_log_epsilon")
    slope = float(kl_df[slope_col].dropna().iloc[0]) \
        if len(kl_df[slope_col].dropna()) else np.nan

    lines = [
        r"\begin{tabular}{lrrrr}",
        r"\toprule",
        r"experiment & $R$ & max / slope & terminal value & check \\",
        r"\midrule",
    ]
    for _, r in riem.iterrows():
        lines.append(
            "Riemannian cascade"
            f" & {r.R:g} & {_fmt_sci(r.max_fr_grad_norm)}"
            f" & {_fmt_sci(r.final_c)} & {_tex_text(r.stop_reason)} \\\\")
    lines.append(r"\midrule")
    lines.append(
        "KL pole"
        f" & {float(kl_df.R.iloc[0]):g} & {slope:.3f}"
        f" & {_fmt_sci(float(kl_df.fr_grad_sq_after_first.max()))}"
        r" & pole-window slope \\")
    lines.append(r"\midrule")
    for _, r in bw_last.iterrows():
        lines.append(
            "BW forward--backward"
            f" & {r.R:g} & {_fmt_sci(r.running_min_bw_grad_sq)}"
            f" & {_fmt_sci(r.theorem_envelope)}"
            f" & {_fmt_bool(r.bound_satisfied)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)
