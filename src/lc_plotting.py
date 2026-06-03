"""
Plotting routines for log-concave gradient flow experiments.

Figures produced (all get ``fig_logconcave_`` prefix):
  1. fig_logconcave_tau_effect_n{n}_rho{r}       — clean main tau comparison
     fig_logconcave_tau_effect_full_n{n}_rho{r}  — full appendix diagnostic grid
  2. fig_logconcave_omega_sweep_n{n}_rho{r}      — omega sweep, normalised gap
  3. fig_logconcave_time_to_tol_n{n}_rho{r}      — time-to-1e-4, one panel per tau
  4. fig_logconcave_tau_speedup_n{n}_rho{r}      — speedup ratio T(tau)/T(tau=0)
  5. fig_logconcave_speedup_vs_chi_n{n}_rho{r}   — scatter: initial_chi vs speedup

All figures saved as PNG (dpi=200) and PDF.

Display-only cleanups (residual floor clipping, chi masking when the covariance
residual is tiny, not-reached heatmap masking, diverging speedup normalisation)
are imported from scripts/plotting_utils.py.  The experiment data is untouched.
Reuses style constants from src/plotting.py.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Reuse shared style constants from the Gaussian plotting module
from src.plotting import (
    INIT_LABELS, INIT_ORDER, INIT_ORDER_CLEAN, TAU_COLORS, TAU_LABELS,
    _save,
)

# Shared display-only helpers (scripts/plotting_utils.py)
_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)
from plotting_utils import (                       # noqa: E402
    PLOT_FLOOR_OBJECTIVE_LOGCONCAVE,
    PLOT_FLOOR_RESIDUAL_LOGCONCAVE,
    PLOT_FLOOR_MEAN_LOGCONCAVE,
    CHI_COV_RESIDUAL_MIN,
    get_metric_column,
    clip_for_log_plot,
    semilogy_clipped,
    mask_chi_when_residual_small,
    is_reached,
    format_not_reached_label,
    plot_heatmap_with_not_reached_mask,
    speedup_diverging_norm,
    apply_manuscript_style,
)

# ---------------------------------------------------------------------------
# Log-concave metric labels and display floors
# ---------------------------------------------------------------------------

LC_METRIC_COLS = [
    "normalized_objective_gap",
    "whitened_mean_error",
    "cov_residual",
    "trace_residual",
    "traceless_residual",
    "chi",
]
LC_METRIC_LABELS = {
    "normalized_objective_gap": "Normalised objective gap",
    "whitened_mean_error":      r"Whitened mean error",
    "cov_residual":             r"Cov. residual $\|I - B\|_F$",
    "trace_residual":           r"Trace residual",
    "traceless_residual":       r"Traceless residual",
    "chi":                      r"$\chi$ (trace dominance)",
}

# Per-metric display floors (display only — data is not modified).
LC_METRIC_FLOORS = {
    "normalized_objective_gap": PLOT_FLOOR_OBJECTIVE_LOGCONCAVE,
    "whitened_mean_error":      PLOT_FLOOR_MEAN_LOGCONCAVE,
    "cov_residual":             PLOT_FLOOR_RESIDUAL_LOGCONCAVE,
    "trace_residual":           PLOT_FLOOR_RESIDUAL_LOGCONCAVE,
    "traceless_residual":       PLOT_FLOOR_RESIDUAL_LOGCONCAVE,
}

# Columns that hold the covariance-residual magnitude (for chi masking).
_COV_RESIDUAL_CANDIDATES = ["cov_residual", "cov_error"]


def _rho_tag(n, rho):
    """Return e.g. 'n5_rho5' (integer rho) or 'n5_rho2.5' (fractional rho)."""
    rho_str = str(int(rho)) if float(rho) == int(rho) else str(rho)
    return f"n{n}_rho{rho_str}"


def _plot_metric_series(ax, s_tt, metric, color, label):
    """Plot one (init, tau) curve for a given metric with display cleanups.

    * chi is plotted on a linear axis with masking where the covariance
      residual is too small for the ratio to be meaningful;
    * every other metric is plotted on a clipped semilogy axis.
    """
    x = s_tt["time"].values
    if metric == "chi":
        res_col = get_metric_column(s_tt, _COV_RESIDUAL_CANDIDATES)
        if res_col is not None:
            y = mask_chi_when_residual_small(
                s_tt["chi"].values, s_tt[res_col].values, CHI_COV_RESIDUAL_MIN)
        else:
            y = s_tt["chi"].values
        ax.plot(x, y, color=color, lw=1.7, label=label)
        ax.set_ylim(0, 1.05)
    else:
        floor = LC_METRIC_FLOORS.get(metric, PLOT_FLOOR_RESIDUAL_LOGCONCAVE)
        semilogy_clipped(ax, x, s_tt[metric].values, floor,
                         color=color, lw=1.7, label=label)


# ---------------------------------------------------------------------------
# Figure 1 (appendix): full tau-effect diagnostic grid
# ---------------------------------------------------------------------------

def plot_lc_tau_effect(df, outdir, n=5, rho=5, omega=0.5):
    """Full diagnostic grid: 5 inits x 6 metrics, comparing the three tau types.

    Residual curves are clipped to a display floor; the chi panel is masked
    where the covariance residual is too small for the ratio to be meaningful.
    """
    tag = _rho_tag(n, rho)
    sub = df[(df["n"] == n) & (df["rho"] == rho) & np.isclose(df["omega"], omega)].copy()
    if sub.empty:
        print(f"[plot_lc_tau_effect] No data for {tag}, omega={omega}. Skipping.")
        return

    inits   = INIT_ORDER
    metrics = LC_METRIC_COLS
    ncols, nrows = len(metrics), len(inits)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 2.7 * nrows), sharey="col")

    for ri, init in enumerate(inits):
        for cj, metric in enumerate(metrics):
            ax = axes[ri, cj]
            s_init = sub[sub["init_name"] == init]
            for tt in ["negative", "zero", "positive"]:
                s_tt = s_init[s_init["tau_type"] == tt].sort_values("time")
                if s_tt.empty:
                    continue
                _plot_metric_series(ax, s_tt, metric, TAU_COLORS[tt], TAU_LABELS[tt])

            ax.set_xlim(left=0)
            if ri == nrows - 1:
                ax.set_xlabel("Time")
            if cj == 0:
                ax.set_ylabel(INIT_LABELS[init], fontsize=9)
            if ri == 0:
                ax.set_title(LC_METRIC_LABELS[metric], fontsize=8)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", fontsize=8,
                   bbox_to_anchor=(1.0, 1.0), framealpha=0.9)

    fig.suptitle(
        rf"$\tau$ effect (full grid) — $n={n}$, $\rho={rho}$, $\omega={omega}$ (log-concave)",
        y=1.005, fontsize=11,
    )
    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_logconcave_tau_effect_full_{tag}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Figure 1b (main): clean tau-effect figure
# ---------------------------------------------------------------------------

def plot_lc_tau_effect_clean(df, outdir, n=5, rho=5, omega=0.5):
    """Clean main tau-effect figure focusing on trace vs traceless behaviour.

    Rows  : volume_high, shape_only, mixed.
    Cols  : normalised objective gap, trace residual, traceless residual, chi
            (chi masked where the covariance residual is tiny).

    Saved as ``fig_logconcave_tau_effect_n{n}_rho{rho}`` (clean main).
    """
    tag = _rho_tag(n, rho)
    sub = df[(df["n"] == n) & (df["rho"] == rho) & np.isclose(df["omega"], omega)].copy()
    if sub.empty:
        print(f"[plot_lc_tau_effect_clean] No data for {tag}, omega={omega}. Skipping.")
        return

    inits   = INIT_ORDER_CLEAN
    metrics = ["normalized_objective_gap", "trace_residual", "traceless_residual", "chi"]
    ncols, nrows = len(metrics), len(inits)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 2.8 * nrows),
                             sharey="col", squeeze=False)

    for ri, init in enumerate(inits):
        for cj, metric in enumerate(metrics):
            ax = axes[ri, cj]
            s_init = sub[sub["init_name"] == init]
            for tt in ["negative", "zero", "positive"]:
                s_tt = s_init[s_init["tau_type"] == tt].sort_values("time")
                if s_tt.empty:
                    continue
                _plot_metric_series(ax, s_tt, metric, TAU_COLORS[tt], TAU_LABELS[tt])

            ax.set_xlim(left=0)
            if ri == nrows - 1:
                ax.set_xlabel("Time")
            if cj == 0:
                ax.set_ylabel(INIT_LABELS[init], fontsize=11)
            if ri == 0:
                ax.set_title(LC_METRIC_LABELS[metric], fontsize=10)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", fontsize=9,
                   bbox_to_anchor=(1.0, 1.0), framealpha=0.9)

    fig.suptitle(
        rf"$\tau$ effect — $n={n}$, $\rho={rho}$, $\omega={omega}$ (log-concave target)",
        y=1.01, fontsize=12,
    )
    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_logconcave_tau_effect_{tag}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Figure 2: omega sweep for tau=0
# ---------------------------------------------------------------------------

def plot_lc_omega_sweep(df, outdir, n=5, rho=5):
    """Normalised objective gap vs time for all omega values (tau=0).

    Curves are clipped to ``PLOT_FLOOR_OBJECTIVE_LOGCONCAVE`` so the MC/QMC
    numerical floor on the objective gap is not drawn as a noisy tail.
    """
    tag = _rho_tag(n, rho)
    sub = df[(df["n"] == n) & (df["rho"] == rho) & (df["tau_type"] == "zero")].copy()
    if sub.empty:
        print(f"[plot_lc_omega_sweep] No data for {tag}, tau=0. Skipping.")
        return

    gap_col = get_metric_column(sub, ["normalized_objective_gap", "norm_energy"])
    floor   = PLOT_FLOOR_OBJECTIVE_LOGCONCAVE

    omegas = sorted(sub["omega"].unique())
    cmap   = plt.get_cmap("plasma", len(omegas))
    o_colors = {o: cmap(i) for i, o in enumerate(omegas)}

    inits = INIT_ORDER
    fig, axes = plt.subplots(1, len(inits), figsize=(3.9 * len(inits), 3.5))
    if len(inits) == 1:
        axes = [axes]

    for cj, init in enumerate(inits):
        ax = axes[cj]
        s = sub[sub["init_name"] == init]
        for omega in omegas:
            s_o = s[np.isclose(s["omega"], omega)].sort_values("time")
            if s_o.empty:
                continue
            semilogy_clipped(ax, s_o["time"], s_o[gap_col], floor,
                             color=o_colors[omega], lw=1.9, label=rf"$\omega={omega}$")
        ax.axhline(floor, color="0.75", lw=0.8, ls=":", zorder=0)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=floor * 0.5)
        ax.set_xlabel("Time")
        if cj == 0:
            ax.set_ylabel("Normalised objective gap")
        ax.set_title(INIT_LABELS[init], fontsize=11)
        ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(rf"$\omega$ sweep — $n={n}$, $\rho={rho}$, $\tau=0$ (log-concave target)",
                 fontsize=12)
    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_logconcave_omega_sweep_{tag}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Figure 3: time-to-tolerance heatmap (one panel per tau type)
# ---------------------------------------------------------------------------

def _lc_time_matrix(s_tt, inits, omegas, tol_col, T_max):
    mat = np.full((len(inits), len(omegas)), np.nan)
    for ri, init in enumerate(inits):
        for cj, omega in enumerate(omegas):
            rows = s_tt[(s_tt["init_name"] == init) & np.isclose(s_tt["omega"], omega)]
            if rows.empty:
                continue
            v = float(rows[tol_col].iloc[0])
            mat[ri, cj] = v if is_reached(v, T_max) else np.nan
    return mat


def plot_lc_time_to_tol_heatmap(sdf, outdir, n=5, rho=5, tol_col="time_to_1e_minus_4"):
    """Time-to-1e-4 heatmaps (rows=init, cols=omega), one panel per tau type.

    Not-reached cells are masked (light gray) and annotated ``>T``.  All three
    panels share one colour scale computed from the finite reached times only.

    Saved as ``fig_logconcave_time_to_tol_n{n}_rho{rho}``.
    """
    tag = _rho_tag(n, rho)
    sub = sdf[(sdf["n"] == n) & (sdf["rho"] == rho)].copy()
    if sub.empty:
        print(f"[plot_lc_time_to_tol_heatmap] No data for {tag}. Skipping.")
        return

    omegas    = sorted(sub["omega"].unique())
    inits     = INIT_ORDER
    tau_types = ["negative", "zero", "positive"]
    T_max     = float(sub["T"].max())

    # Build all matrices first to derive a common colour scale.
    mats = {tt: _lc_time_matrix(sub[sub["tau_type"] == tt], inits, omegas, tol_col, T_max)
            for tt in tau_types}
    all_finite = np.concatenate([m[np.isfinite(m)] for m in mats.values()
                                 if np.isfinite(m).any()]) if any(
        np.isfinite(m).any() for m in mats.values()) else np.array([0.0, 1.0])
    vmin, vmax = float(all_finite.min()), float(all_finite.max())

    fig, axes = plt.subplots(1, 3, figsize=(4.6 * 3, 4.5), sharey=True)

    last_im = None
    for ax, tt in zip(axes, tau_types):
        last_im = plot_heatmap_with_not_reached_mask(
            ax, mats[tt], T=T_max,
            xticklabels=[str(o) for o in omegas],
            yticklabels=[INIT_LABELS[i] for i in inits],
            cmap="viridis_r", fig=None, vmin=vmin, vmax=vmax,
        )
        ax.set_xlabel(r"$\omega$")
        ax.set_title(TAU_LABELS[tt], fontsize=9)

    cbar = fig.colorbar(last_im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("Time")

    fig.suptitle(
        rf"Time to normalised gap $\leq 10^{{-4}}$ — $n={n}$, $\rho={rho}$ (log-concave)",
        fontsize=11,
    )

    path = os.path.join(outdir, f"fig_logconcave_time_to_tol_{tag}")
    fig.savefig(path + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(path + ".pdf", bbox_inches="tight")
    print(f"Saved {path}.png / .pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4: tau speedup heatmap (diverging, centred at 1)
# ---------------------------------------------------------------------------

def _lc_speedup_matrix(sub, inits, target_omegas, tau_side, tol_col, T_max):
    mat = np.full((len(inits), len(target_omegas)), np.nan)
    for ri, init in enumerate(inits):
        for cj, omega in enumerate(target_omegas):
            r_zero = sub[(sub["init_name"] == init) &
                         np.isclose(sub["omega"], omega) &
                         (sub["tau_type"] == "zero")]
            r_tau  = sub[(sub["init_name"] == init) &
                         np.isclose(sub["omega"], omega) &
                         (sub["tau_type"] == tau_side)]
            if r_zero.empty or r_tau.empty:
                continue
            t0  = float(r_zero[tol_col].iloc[0])
            t_t = float(r_tau[tol_col].iloc[0])
            z_ok = is_reached(t0, T_max)
            t_ok = is_reached(t_t, T_max)
            if not z_ok and not t_ok:
                mat[ri, cj] = 1.0
            elif not z_ok:
                mat[ri, cj] = 0.4
            elif not t_ok:
                mat[ri, cj] = 2.5
            else:
                mat[ri, cj] = t_t / t0
    return mat


def plot_lc_tau_speedup_heatmap(sdf, outdir, n=5, rho=5, tol_col="time_to_1e_minus_4"):
    """Ratio T(tau)/T(tau=0) for tau_neg and tau_pos, omega in {1/4, 1/2, 1}.

    Diverging colour map centred exactly at ratio = 1 (no speedup).

    Saved as ``fig_logconcave_tau_speedup_n{n}_rho{rho}``.
    """
    tag = _rho_tag(n, rho)
    target_omegas = [0.25, 0.5, 1.0]
    sub = sdf[
        (sdf["n"] == n) & (sdf["rho"] == rho) &
        sdf["omega"].apply(lambda o: any(np.isclose(o, target_omegas)))
    ].copy()
    if sub.empty:
        print(f"[plot_lc_tau_speedup_heatmap] No data for {tag}. Skipping.")
        return

    T_max = float(sub["T"].max())
    inits = INIT_ORDER
    norm  = speedup_diverging_norm(vmin=0.5, vcenter=1.0, vmax=1.5)
    cmap  = plt.get_cmap("RdYlGn_r")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)

    last_im = None
    for ax, tau_side in zip(axes, ["negative", "positive"]):
        mat = _lc_speedup_matrix(sub, inits, target_omegas, tau_side, tol_col, T_max)
        masked = np.ma.masked_invalid(mat)
        cmap_obj = cmap.copy()
        cmap_obj.set_bad(color="0.85")

        last_im = ax.imshow(masked, aspect="auto", cmap=cmap_obj,
                            norm=norm, origin="upper")
        ax.set_xticks(range(len(target_omegas)))
        ax.set_xticklabels([str(o) for o in target_omegas])
        ax.set_yticks(range(len(inits)))
        ax.set_yticklabels([INIT_LABELS[i] for i in inits])
        ax.set_xlabel(r"$\omega$")
        sign = "-" if tau_side == "negative" else "+"
        ax.set_title(rf"$\tau {sign}$ : $T(\tau)/T(\tau{{=}}0)$", fontsize=11)

        for ri in range(len(inits)):
            for cj in range(len(target_omegas)):
                v = mat[ri, cj]
                if not np.isfinite(v):
                    continue
                ax.text(cj, ri, f"{v:.2f}", ha="center", va="center",
                        fontsize=8, color="black")

    cbar = fig.colorbar(last_im, ax=axes, fraction=0.046, pad=0.04)
    cbar.set_label(r"$T(\tau)/T(\tau{=}0)$  (1 = no change)")

    fig.suptitle(
        rf"$\tau$ speedup ratio — $n={n}$, $\rho={rho}$ (log-concave target)",
        fontsize=12,
    )

    path = os.path.join(outdir, f"fig_logconcave_tau_speedup_{tag}")
    fig.savefig(path + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(path + ".pdf", bbox_inches="tight")
    print(f"Saved {path}.png / .pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 5: speedup vs chi scatter
# ---------------------------------------------------------------------------

def plot_lc_speedup_vs_chi(sdf, outdir, n=5, rho=5, tol_col="time_to_1e_minus_4"):
    """Scatter: x = initial_chi, y = T(tau=0)/T(tau_neg) (speedup > 1 is better).

    chi diagnoses the trace dominance of the *covariance* residual, not the
    total objective.  mean_only is dominated by mean error, so chi does not
    predict its speedup — it is therefore drawn as hollow gray markers and
    labelled "mean-dominated (excluded)", separate from the main trend.
    """
    tag = _rho_tag(n, rho)
    sub = sdf[(sdf["n"] == n) & (sdf["rho"] == rho)].copy()
    if sub.empty:
        print(f"[plot_lc_speedup_vs_chi] No data for {tag}. Skipping.")
        return

    omegas = sorted(sub["omega"].unique())
    markers = {"volume_high": "s", "volume_low": "^", "shape_only": "D", "mixed": "P"}
    cmap   = plt.get_cmap("plasma", len(omegas))
    o_col  = {o: cmap(i) for i, o in enumerate(omegas)}

    fig, ax = plt.subplots(figsize=(7.2, 5.2))

    def _speedup(init, omega):
        r_zero = sub[(sub["init_name"] == init) & np.isclose(sub["omega"], omega) &
                     (sub["tau_type"] == "zero")]
        r_neg  = sub[(sub["init_name"] == init) & np.isclose(sub["omega"], omega) &
                     (sub["tau_type"] == "negative")]
        if r_zero.empty or r_neg.empty:
            return None
        chi_val = float(r_zero["initial_chi"].iloc[0])
        t_zero  = float(r_zero[tol_col].iloc[0])
        t_neg   = float(r_neg[tol_col].iloc[0])
        if not is_reached(t_zero, float(r_zero["T"].iloc[0])) or \
           not is_reached(t_neg, float(r_neg["T"].iloc[0])) or t_neg < 1e-15:
            return None
        return chi_val, t_zero / t_neg

    # Main trend: covariance-structured initializations (exclude mean_only).
    for omega in omegas:
        for init in [i for i in INIT_ORDER if i != "mean_only"]:
            res = _speedup(init, omega)
            if res is None:
                continue
            chi_val, speedup = res
            ax.scatter(chi_val, speedup, color=o_col[omega], marker=markers[init],
                       s=90, zorder=3, edgecolors="black", linewidths=0.4,
                       label=f"{INIT_LABELS[init]}|{omega}")

    # mean_only: hollow gray, shown but excluded from the trend.
    for omega in omegas:
        res = _speedup("mean_only", omega)
        if res is None:
            continue
        chi_val, speedup = res
        ax.scatter(chi_val, speedup, facecolors="none", edgecolors="0.55",
                   marker="o", s=70, zorder=2,
                   label="mean-dominated (excluded)|x")

    ax.axhline(1.0, color="gray", lw=1.2, ls="--", zorder=1)
    ax.text(0.02, 1.02, "no speedup", color="gray", fontsize=8,
            transform=ax.get_yaxis_transform())
    ax.set_xlabel(r"Initial $\chi$ (trace dominance of covariance residual)", fontsize=11)
    ax.set_ylabel(r"Speedup $T(\tau{=}0)/T(\tau_-)$", fontsize=11)
    ax.set_title(rf"$\tau_-$ speedup vs trace dominance — $n={n}$, $\rho={rho}$",
                 fontsize=11)
    ax.set_xlim(0, 1.05)
    ax.text(0.5, -0.16,
            r"$\chi$ diagnoses covariance-residual trace dominance, not total objective dominance.",
            ha="center", va="top", fontsize=8, color="0.35", transform=ax.transAxes)

    # Compact legend: dedupe by init label (ignore the |omega suffix).
    seen, handles, labels_leg = set(), [], []
    for h, l in zip(*ax.get_legend_handles_labels()):
        key = l.split("|")[0]
        if key not in seen:
            seen.add(key)
            handles.append(h)
            labels_leg.append(key)
    ax.legend(handles, labels_leg, fontsize=8, loc="upper left",
              framealpha=0.85, ncol=2)

    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_logconcave_speedup_vs_chi_{tag}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

def generate_lc_figures(long_csv, summary_csv, figures_dir,
                         n=5, rho=5, omega_half=0.5):
    """Load CSVs and produce all log-concave figures (clean main + appendix)."""
    apply_manuscript_style()
    os.makedirs(figures_dir, exist_ok=True)
    df  = pd.read_csv(long_csv)
    sdf = pd.read_csv(summary_csv)

    plot_lc_tau_effect_clean(df, figures_dir, n=n, rho=rho, omega=omega_half)
    plot_lc_tau_effect(df, figures_dir, n=n, rho=rho, omega=omega_half)
    plot_lc_omega_sweep(df, figures_dir, n=n, rho=rho)
    plot_lc_time_to_tol_heatmap(sdf, figures_dir, n=n, rho=rho)
    plot_lc_tau_speedup_heatmap(sdf, figures_dir, n=n, rho=rho)
    plot_lc_speedup_vs_chi(sdf, figures_dir, n=n, rho=rho)
