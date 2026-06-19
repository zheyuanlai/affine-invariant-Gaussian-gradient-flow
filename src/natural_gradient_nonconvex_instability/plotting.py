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


def _by_rule(df, rule):
    """Select one stepsize regime; tolerate a legacy single-rule frame."""
    if "dt_rule" in df.columns:
        return df[df.dt_rule == rule]
    return df


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


def fig_clipped_kl_stationarity(clipped_df, figs_dir):
    """Running-minimum Bregman displacement vs the energy-drop envelope.

    Solid: ``D_min(N) = min_{0<=n<N} KL(rho_{a_n}||rho_{a_{n+1}})``.
    Dashed: ``B_N = (dt / N) {F(c_0) - F(c_N)}``. The theorem check is
    ``D_min(N) <= B_N``; once the upper covariance constraint binds, the
    projected step is stationary at the boundary and ``D_n`` (hence ``D_min``)
    collapses to zero (off the log axis).
    """
    clipped_df = _by_rule(clipped_df, "theorem_safe")
    fig, ax = plt.subplots(figsize=(5.8, 3.9))
    for R in sorted(clipped_df.R.dropna().unique()):
        run = clipped_df[np.isclose(clipped_df.R, R)].sort_values("step")
        N = run.step.values + 1
        color = _color_for_R(R)
        ax.loglog(
            N,
            _positive_finite(run.running_min_D.values),
            color=color,
            lw=1.6,
            label=rf"$R={R:g}$ $D_{{\min}}(N)$",
        )
        ax.loglog(
            N,
            _positive_finite(run.prefix_bound.values),
            color=color,
            ls="--",
            lw=1.1,
            alpha=0.75,
        )
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"Bregman displacement $\mathrm{KL}(\rho_{a_n}\|\rho_{a_{n+1}})$")
    ax.set_title(r"Projected KL: $D_{\min}(N)$ vs envelope $B_N$")
    ax.legend(ncol=2, fontsize=8)
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    return _savefig(fig, figs_dir, "fig_nonconvex_clipped_kl_stationarity")


def fig_clipped_covariance(clipped_df, figs_dir):
    """Clipped covariance trajectory with the feasible interval overlaid."""
    clipped_df = _by_rule(clipped_df, "theorem_safe")
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    lam_m = float(clipped_df.lambda_minus.iloc[0])
    lam_p = float(clipped_df.lambda_plus.iloc[0])
    for R in sorted(clipped_df.R.dropna().unique()):
        run = clipped_df[np.isclose(clipped_df.R, R)].sort_values("step")
        ax.plot(
            run.step.values,
            run.c_next.values,
            color=_color_for_R(R),
            lw=1.6,
            label=rf"$R={R:g}$",
        )
    ax.axhline(lam_p, color="0.35", ls="--", lw=1.0, label=rf"$\lambda_+={lam_p:g}$")
    ax.axhline(lam_m, color="0.6", ls=":", lw=1.0, label=rf"$\lambda_-={lam_m:g}$")
    ax.set_xlabel("step $n$")
    ax.set_ylabel(r"clipped covariance $c_{n+1}$")
    ax.set_title(r"Projected covariance pinned to $[\lambda_-,\lambda_+]$")
    ax.legend(ncol=2, fontsize=8)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    return _savefig(fig, figs_dir, "fig_nonconvex_clipped_covariance")


def fig_clipped_kl_largestep(clipped_df, figs_dir):
    """Non-theorem-safe stepsize ``dt = 1/(beta lambda_+)`` (``dt L_clip > 1``).

    Left: ``D_min(N)`` (solid) vs the envelope ``B_N`` (dashed) for the larger
    stepsize; right: the covariance trajectory, which reaches the upper bound in
    a single step. The Theorem 2.18 condition does not hold here, yet the
    stationarity envelope is still satisfied empirically because the one large
    step makes a big energy drop and lands on the active upper constraint, where
    the projected step is stationary and ``D_n = 0``.
    """
    rs = _by_rule(clipped_df, "riemannian_scale")
    if rs.empty:
        return None
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.6, 3.7))
    for R in sorted(rs.R.dropna().unique()):
        run = rs[np.isclose(rs.R, R)].sort_values("step")
        color = _color_for_R(R)
        N = run.step.values + 1
        axL.loglog(N, _positive_finite(run.running_min_D.values), color=color,
                   lw=1.6, marker="o", markersize=5,
                   label=rf"$R={R:g}$ $D_{{\min}}(N)$")
        axL.loglog(N, _positive_finite(run.prefix_bound.values), color=color,
                   ls="--", lw=1.1, alpha=0.75)
        axR.plot(run.step.values, run.c_next.values, color=color, lw=1.6,
                 marker="o", markersize=3, label=rf"$R={R:g}$")
    axL.set_xlabel(r"$N$")
    axL.set_ylabel(r"Bregman displacement $\mathrm{KL}(\rho_{a_n}\|\rho_{a_{n+1}})$")
    axL.set_title(r"$D_{\min}(N)$ vs envelope $B_N$")
    axL.legend(ncol=2, fontsize=8)
    axL.grid(True, which="both", alpha=0.25)
    lam_p = float(rs.lambda_plus.iloc[0])
    lam_m = float(rs.lambda_minus.iloc[0])
    axR.axhline(lam_p, color="0.35", ls="--", lw=1.0, label=rf"$\lambda_+={lam_p:g}$")
    axR.axhline(lam_m, color="0.6", ls=":", lw=1.0, label=rf"$\lambda_-={lam_m:g}$")
    axR.set_xlim(-0.3, 10)
    axR.set_xlabel("step $n$")
    axR.set_ylabel(r"clipped covariance $c_{n+1}$")
    axR.set_title(r"Covariance reaches $\lambda_+$ in one step")
    axR.legend(ncol=2, fontsize=8)
    axR.grid(True, alpha=0.25)
    fig.suptitle(r"Projected KL beyond the theorem-safe regime: "
                 r"$\Delta t=1/(\beta\lambda_+)$, $\Delta t\,L_{\mathrm{clip}}=64$",
                 y=1.02, fontsize=11)
    fig.tight_layout()
    return _savefig(fig, figs_dir, "fig_nonconvex_clipped_kl_largestep")


def build_all_figures(long_df, kl_df, bw_df, clipped_df, figs_dir):
    os.makedirs(figs_dir, exist_ok=True)
    fig_riemannian_cascade(long_df, figs_dir)
    fig_kl_pole(kl_df, figs_dir)
    fig_wasserstein_bound(bw_df, figs_dir)
    fig_curvature_diagnostic(long_df, figs_dir)
    fig_clipped_kl_stationarity(clipped_df, figs_dir)
    fig_clipped_covariance(clipped_df, figs_dir)
    fig_clipped_kl_largestep(clipped_df, figs_dir)


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


def clipped_kl_summary_table_tex(clipped_summary_df):
    """Per-run projected-KL stationarity summary (Theorem 2.18 check)."""
    df = _by_rule(clipped_summary_df, "theorem_safe").sort_values("R")
    lines = [
        r"\begin{tabular}{lrrrrl}",
        r"\toprule",
        r"$R$ & $\Delta t$ & $D_{\min}(N)$ & $B_N$ & max. viol. & check \\",
        r"\midrule",
    ]
    for _, r in df.iterrows():
        lines.append(
            f"{r.R:g} & {_fmt_sci(r['dt'])} & {_fmt_sci(r.final_running_min_D)}"
            f" & {_fmt_sci(r.final_bound)} & {_fmt_sci(r.max_violation)}"
            f" & {_fmt_bool(r.theorem_check_pass)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def clipped_kl_largestep_table_tex(clipped_summary_df):
    """Compare the theorem-safe and non-theorem-safe (``dt=1/(beta lambda_+)``)
    stepsize rules: stepsize, ``dt*L_clip``, energy monotonicity, max violation,
    and the envelope check."""
    safe = _by_rule(clipped_summary_df, "theorem_safe").sort_values("R")
    rs = _by_rule(clipped_summary_df, "riemannian_scale").sort_values("R")
    lines = [
        r"\begin{tabular}{llrrrrll}",
        r"\toprule",
        r"$R$ & regime & $\Delta t$ & $\Delta t\,L_{\mathrm{clip}}$ & "
        r"$D_{\min}(N)$ & max. viol. & energy mono. & check \\",
        r"\midrule",
    ]

    def _row(r, label, show_R):
        dtl = r["dt_times_L_clip"] if "dt_times_L_clip" in r else float("nan")
        rlab = f"{r.R:g}" if show_R else ""
        return (f"{rlab} & {label} & {_fmt_sci(r['dt'])} & {_fmt_sci(dtl)}"
                f" & {_fmt_sci(r.final_running_min_D)} & {_fmt_sci(r.max_violation)}"
                f" & {_fmt_bool(r.energy_monotone)}"
                f" & {_fmt_bool(r.theorem_check_pass)} \\\\")

    rs_by_R = {float(r.R): r for _, r in rs.iterrows()}
    for _, r in safe.iterrows():
        lines.append(_row(r, "theorem-safe", show_R=True))
        rr = rs_by_R.get(float(r.R))
        if rr is not None:
            lines.append(_row(rr, r"$1/(\beta\lambda_+)$", show_R=False))
        lines.append(r"\midrule")
    if lines[-1] == r"\midrule":
        lines.pop()
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)
