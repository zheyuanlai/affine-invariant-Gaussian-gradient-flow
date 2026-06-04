"""
Plotting routines for affine-invariant Gaussian gradient flow experiments.

All figures are saved as both PNG (dpi=200) and PDF.
Requires pandas for DataFrame filtering; matplotlib for rendering.
No seaborn dependency.

Display-only cleanups (curve clipping, not-reached heatmap masking, diverging
speedup normalisation) live in ``scripts/plotting_utils.py`` and are applied
here.  None of these change the underlying experiment data.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd

# Shared display-only helpers (curve clipping, heatmap masking, diverging
# speedup normalisation). None of these change the underlying experiment data.
from src.omega_tau_modes.plotting_utils import (
    PLOT_FLOOR_OBJECTIVE_GAUSSIAN,
    PLOT_FLOOR_RESIDUAL_GAUSSIAN,
    get_metric_column,
    clip_for_log_plot,
    semilogy_clipped,
    series_is_identically_small,
    is_reached,
    format_not_reached_label,
    plot_heatmap_with_not_reached_mask,
    speedup_diverging_norm,
    apply_manuscript_style,
)


# ---------------------------------------------------------------------------
# Style constants
# ---------------------------------------------------------------------------

METRIC_COLS = [
    "norm_energy", "mean_error", "cov_error",
    "volume_error", "shape_error", "cosine_error",
]
METRIC_LABELS = {
    "norm_energy":  "Normalised KL energy",
    "mean_error":   r"Mean error $\|m\|_2$",
    "cov_error":    r"Rel. cov. error $\|C-I\|_F/\sqrt{n}$",
    "volume_error": r"Volume error $|\log\det C / n|$",
    "shape_error":  r"Shape error (trace-free $\log C$)",
    "cosine_error": "Cosine test-fn error",
}

# Per-metric display floor: the objective (energy) uses a looser floor than the
# residual-type diagnostics.  Curves are pinned to these floors for display.
GAUSSIAN_METRIC_FLOORS = {
    "norm_energy":  PLOT_FLOOR_OBJECTIVE_GAUSSIAN,
    "mean_error":   PLOT_FLOOR_RESIDUAL_GAUSSIAN,
    "cov_error":    PLOT_FLOOR_RESIDUAL_GAUSSIAN,
    "volume_error": PLOT_FLOOR_RESIDUAL_GAUSSIAN,
    "shape_error":  PLOT_FLOOR_RESIDUAL_GAUSSIAN,
    "cosine_error": PLOT_FLOOR_RESIDUAL_GAUSSIAN,
}

INIT_LABELS = {
    "mean_only":    "Mean only",
    "volume_high":  "Volume high",
    "volume_low":   "Volume low",
    "shape_only":   "Shape only",
    "mixed":        "Mixed",
}

INIT_ORDER = ["mean_only", "volume_high", "volume_low", "shape_only", "mixed"]

# Reduced row set for the clean main tau-effect figure (drops modes whose error
# panels are zero by construction, keeping the figure readable).
INIT_ORDER_CLEAN = ["volume_high", "shape_only", "mixed"]

TAU_COLORS = {
    "negative": "#d62728",   # red
    "zero":     "#1f77b4",   # blue
    "positive": "#2ca02c",   # green
}
TAU_LABELS = {
    "negative": r"$\tau < 0$  ($\tau = -\omega/2n$)",
    "zero":     r"$\tau = 0$",
    "positive": r"$\tau > 0$  ($\tau = +\omega/2n$)",
}


def _save(fig, path_no_ext):
    """Save figure as PNG and PDF."""
    fig.savefig(path_no_ext + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(path_no_ext + ".pdf", bbox_inches="tight")
    plt.close(fig)


def _semilogy_safe(ax, x, y, **kwargs):
    """Plot with semilogy, skipping non-positive y values gracefully.

    Retained for backward compatibility; new code should prefer
    ``semilogy_clipped`` which pins a display floor instead of dropping points.
    """
    mask = np.asarray(y, dtype=float) > 0
    xs = np.asarray(x)[mask]
    ys = np.asarray(y)[mask]
    if len(xs):
        ax.semilogy(xs, ys, **kwargs)


# ---------------------------------------------------------------------------
# Figure 1 (appendix): full tau-effect diagnostic grid
# ---------------------------------------------------------------------------

def plot_tau_effect(df, outdir, n=5, omega=0.5):
    """Full diagnostic grid: 5 inits x 6 metrics, comparing the three tau types.

    Display cleanups:
      * each curve is clipped to a per-metric display floor (no machine-precision
        tails);
      * a panel that is identically zero by construction (e.g. covariance error
        for the mean-only init) is annotated "identically zero" instead of being
        shown as an empty log axis.

    Saved as ``fig_gaussian_tau_effect_full_n{n}`` (appendix).
    """
    sub = df[(df["n"] == n) & (np.isclose(df["omega"], omega))].copy()
    if sub.empty:
        print(f"[plot_tau_effect] No data for n={n}, omega={omega}. Skipping.")
        return

    inits = INIT_ORDER
    metrics = METRIC_COLS
    ncols = len(metrics)
    nrows = len(inits)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3.4 * ncols, 2.7 * nrows),
        sharey="col",
    )

    for row_i, init in enumerate(inits):
        for col_j, metric in enumerate(metrics):
            ax = axes[row_i, col_j]
            s_init = sub[sub["init_name"] == init]
            floor = GAUSSIAN_METRIC_FLOORS.get(metric, PLOT_FLOOR_RESIDUAL_GAUSSIAN)

            # Detect a panel that is identically zero by construction.
            all_vals = s_init[metric].to_numpy(dtype=float) if not s_init.empty else np.array([])
            if series_is_identically_small(all_vals, floor):
                ax.set_xticks([])
                ax.set_yticks([])
                ax.text(0.5, 0.5, "identically zero\nby construction",
                        ha="center", va="center", fontsize=8, color="0.45",
                        transform=ax.transAxes)
                if row_i == nrows - 1:
                    ax.set_xlabel("Time")
                if col_j == 0:
                    ax.set_ylabel(INIT_LABELS[init], fontsize=9)
                if row_i == 0:
                    ax.set_title(METRIC_LABELS[metric], fontsize=9)
                continue

            for tt in ["negative", "zero", "positive"]:
                s_tt = s_init[s_init["tau_type"] == tt].sort_values("time")
                if s_tt.empty:
                    continue
                semilogy_clipped(
                    ax, s_tt["time"], s_tt[metric], floor,
                    color=TAU_COLORS[tt], lw=1.6, label=TAU_LABELS[tt],
                )

            ax.set_xlim(left=0)
            if row_i == nrows - 1:
                ax.set_xlabel("Time")
            if col_j == 0:
                ax.set_ylabel(INIT_LABELS[init], fontsize=9)
            if row_i == 0:
                ax.set_title(METRIC_LABELS[metric], fontsize=9)

    handles, labels = axes[1, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles, labels,
            loc="upper right", fontsize=9,
            bbox_to_anchor=(1.0, 1.0),
            framealpha=0.9,
        )

    fig.suptitle(
        rf"$\tau$ effect (full diagnostic grid) — $n={n}$, $\omega={omega}$",
        y=1.005, fontsize=12,
    )
    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_gaussian_tau_effect_full_n{n}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Figure 1b (main): clean tau-effect figure
# ---------------------------------------------------------------------------

def plot_tau_effect_clean(df, outdir, n=5, omega=0.5):
    """Clean main tau-effect figure.

    Rows  : volume_high, shape_only, mixed  (modes with informative covariance
            structure; the empty mean-error / cov-error panels are dropped).
    Cols  : normalised KL energy, covariance error, volume error, shape error,
            cosine error.

    Saved as ``fig_gaussian_tau_effect_n{n}`` (clean main).
    """
    sub = df[(df["n"] == n) & (np.isclose(df["omega"], omega))].copy()
    if sub.empty:
        print(f"[plot_tau_effect_clean] No data for n={n}, omega={omega}. Skipping.")
        return

    inits = INIT_ORDER_CLEAN
    metrics = ["norm_energy", "cov_error", "volume_error", "shape_error", "cosine_error"]
    ncols, nrows = len(metrics), len(inits)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3.4 * ncols, 2.8 * nrows),
        sharey="col", squeeze=False,
    )

    for row_i, init in enumerate(inits):
        for col_j, metric in enumerate(metrics):
            ax = axes[row_i, col_j]
            s_init = sub[sub["init_name"] == init]
            floor = GAUSSIAN_METRIC_FLOORS.get(metric, PLOT_FLOOR_RESIDUAL_GAUSSIAN)

            for tt in ["negative", "zero", "positive"]:
                s_tt = s_init[s_init["tau_type"] == tt].sort_values("time")
                if s_tt.empty:
                    continue
                semilogy_clipped(
                    ax, s_tt["time"], s_tt[metric], floor,
                    color=TAU_COLORS[tt], lw=1.9, label=TAU_LABELS[tt],
                )

            ax.set_xlim(left=0)
            if row_i == nrows - 1:
                ax.set_xlabel("Time")
            if col_j == 0:
                ax.set_ylabel(INIT_LABELS[init], fontsize=11)
            if row_i == 0:
                ax.set_title(METRIC_LABELS[metric], fontsize=10)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", fontsize=9,
                   bbox_to_anchor=(1.0, 1.0), framealpha=0.9)

    fig.suptitle(rf"$\tau$ effect — $n={n}$, $\omega={omega}$ (Gaussian target)",
                 y=1.01, fontsize=12)
    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_gaussian_tau_effect_n{n}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Figure 2: omega sweep for tau=0
# ---------------------------------------------------------------------------

def plot_omega_sweep(df, outdir, n=5):
    """Normalised energy vs time for all omega values (tau=0), one panel per init.

    Display cleanup: curves are clipped to ``PLOT_FLOOR_OBJECTIVE_GAUSSIAN`` so
    they no longer dive into the 1e-15..1e-16 tail; a subtle dashed floor line
    marks where the display floor sits.
    """
    sub = df[(df["n"] == n) & (df["tau_type"] == "zero")].copy()
    if sub.empty:
        print(f"[plot_omega_sweep] No data for n={n}, tau=0. Skipping.")
        return

    energy_col = get_metric_column(sub, ["norm_energy", "normalized_objective_gap"])
    floor = PLOT_FLOOR_OBJECTIVE_GAUSSIAN

    omegas = sorted(sub["omega"].unique())
    cmap = plt.get_cmap("plasma", len(omegas))
    omega_colors = {o: cmap(i) for i, o in enumerate(omegas)}

    inits = INIT_ORDER
    fig, axes = plt.subplots(1, len(inits), figsize=(3.9 * len(inits), 3.5))
    if len(inits) == 1:
        axes = [axes]

    for col_j, init in enumerate(inits):
        ax = axes[col_j]
        s_init = sub[sub["init_name"] == init]

        for omega in omegas:
            s_o = s_init[np.isclose(s_init["omega"], omega)].sort_values("time")
            if s_o.empty:
                continue
            semilogy_clipped(
                ax, s_o["time"], s_o[energy_col], floor,
                color=omega_colors[omega], lw=1.9, label=rf"$\omega={omega}$",
            )

        # Subtle display-floor reference line.
        ax.axhline(floor, color="0.75", lw=0.8, ls=":", zorder=0)

        ax.set_xlim(left=0)
        ax.set_ylim(bottom=floor * 0.5)
        ax.set_xlabel("Time")
        if col_j == 0:
            ax.set_ylabel("Normalised KL energy")
        ax.set_title(INIT_LABELS[init], fontsize=11)
        ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(rf"$\omega$ sweep — $n={n}$, $\tau=0$ (Gaussian target)", fontsize=12)
    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_gaussian_omega_sweep_n{n}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Figure 3: time-to-tolerance heatmap (tau=0)
# ---------------------------------------------------------------------------

def _time_to_tol_matrix(sub, inits, omegas, tol_col, T_max):
    """Build a matrix of reached times; NaN where not reached / missing."""
    mat = np.full((len(inits), len(omegas)), np.nan)
    for r_i, init in enumerate(inits):
        for c_j, omega in enumerate(omegas):
            rows = sub[(sub["init_name"] == init) & np.isclose(sub["omega"], omega)]
            if rows.empty:
                continue
            val = float(rows[tol_col].iloc[0])
            mat[r_i, c_j] = val if is_reached(val, T_max) else np.nan
    return mat


def plot_time_to_tol_heatmap(summary_df, outdir, n=5, tol_col="time_to_1e_minus_4"):
    """Heatmap of time-to-1e-4 (rows=init, cols=omega, tau=0).

    Not-reached cells are masked (light gray) and annotated ``>T`` rather than
    drawn with a fabricated value such as ``T + 2``.  The colour scale uses only
    the finite reached times.

    Saved as ``fig_gaussian_time_to_tol_n{n}``.
    """
    sub = summary_df[(summary_df["n"] == n) & (summary_df["tau_type"] == "zero")].copy()
    if sub.empty:
        print(f"[plot_time_to_tol_heatmap] No data for n={n}, tau=0. Skipping.")
        return

    omegas = sorted(sub["omega"].unique())
    inits  = INIT_ORDER
    T_max  = float(sub["T"].max())

    mat = _time_to_tol_matrix(sub, inits, omegas, tol_col, T_max)

    fig, ax = plt.subplots(figsize=(7, 4))
    plot_heatmap_with_not_reached_mask(
        ax, mat, T=T_max,
        xticklabels=[str(o) for o in omegas],
        yticklabels=[INIT_LABELS[i] for i in inits],
        cmap="viridis_r", cbar_label="Time", fig=fig,
    )
    ax.set_xlabel(r"$\omega$")
    ax.set_title(rf"Time to $E/E_0 \leq 10^{{-4}}$ — $n={n}$, $\tau=0$", fontsize=11)
    fig.tight_layout()

    path = os.path.join(outdir, f"fig_gaussian_time_to_tol_n{n}")
    fig.savefig(path + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(path + ".pdf", bbox_inches="tight")
    print(f"Saved {path}.png / .pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 4: tau speedup heatmap (diverging, centred at 1)
# ---------------------------------------------------------------------------

def _speedup_matrix(sub, inits, target_omegas, tau_side, tol_col, T_max):
    """Ratio T(tau)/T(tau=0); NaN where undefined.

    Stores the raw ratio for finite/finite cells.  Special cases:
      both not reached -> 1.0 (no change),
      only tau reached -> small ratio (tau much faster),
      only zero reached -> large ratio (tau much slower).
    """
    mat = np.full((len(inits), len(target_omegas)), np.nan)
    for r_i, init in enumerate(inits):
        for c_j, omega in enumerate(target_omegas):
            r_zero = sub[(sub["init_name"] == init) &
                         np.isclose(sub["omega"], omega) &
                         (sub["tau_type"] == "zero")]
            r_tau  = sub[(sub["init_name"] == init) &
                         np.isclose(sub["omega"], omega) &
                         (sub["tau_type"] == tau_side)]
            if r_zero.empty or r_tau.empty:
                continue
            t_zero = float(r_zero[tol_col].iloc[0])
            t_tau  = float(r_tau[tol_col].iloc[0])
            z_ok = is_reached(t_zero, T_max)
            t_ok = is_reached(t_tau, T_max)
            if not z_ok and not t_ok:
                mat[r_i, c_j] = 1.0
            elif not z_ok:
                mat[r_i, c_j] = 0.4          # tau reached, zero did not -> fast
            elif not t_ok:
                mat[r_i, c_j] = 2.5          # zero reached, tau did not -> slow
            else:
                mat[r_i, c_j] = t_tau / t_zero
    return mat


def plot_tau_speedup_heatmap(summary_df, outdir, n=5, tol_col="time_to_1e_minus_4"):
    """Ratio T(tau)/T(tau=0) for tau_neg and tau_pos, omega in {1/4, 1/2, 1}.

    Uses a diverging colour map centred exactly at ratio = 1 (no speedup):
        ratio < 1  faster than tau=0   ratio > 1  slower than tau=0.

    Saved as ``fig_gaussian_tau_speedup_n{n}``.
    """
    target_omegas = [0.25, 0.5, 1.0]
    sub = summary_df[
        (summary_df["n"] == n) &
        summary_df["omega"].apply(lambda o: any(np.isclose(o, target_omegas)))
    ].copy()
    if sub.empty:
        print(f"[plot_tau_speedup_heatmap] No data for n={n}. Skipping.")
        return

    T_max = float(sub["T"].max())
    inits = INIT_ORDER
    norm  = speedup_diverging_norm(vmin=0.5, vcenter=1.0, vmax=1.5)
    cmap  = plt.get_cmap("RdYlGn_r")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)

    last_im = None
    for ax, tau_side in zip(axes, ["negative", "positive"]):
        mat = _speedup_matrix(sub, inits, target_omegas, tau_side, tol_col, T_max)
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

        for r_i in range(len(inits)):
            for c_j in range(len(target_omegas)):
                v = mat[r_i, c_j]
                if not np.isfinite(v):
                    continue
                ax.text(c_j, r_i, f"{v:.2f}", ha="center", va="center",
                        fontsize=8, color="black")

    cbar = fig.colorbar(last_im, ax=axes, fraction=0.046, pad=0.04)
    cbar.set_label(r"$T(\tau)/T(\tau{=}0)$  (1 = no change)")

    fig.suptitle(
        rf"$\tau$ speedup ratio — $n={n}$ (Gaussian target)",
        fontsize=12,
    )

    path = os.path.join(outdir, f"fig_gaussian_tau_speedup_n{n}")
    fig.savefig(path + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(path + ".pdf", bbox_inches="tight")
    print(f"Saved {path}.png / .pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

def generate_all_figures(long_csv, summary_csv, figures_base_dir, n_list=None):
    """Load CSVs and produce all Gaussian figures for every n in n_list.

    Each n gets its own subdirectory: <figures_base_dir>/n{N}/.
    Both the clean main figures and the full appendix grid are produced.
    """
    apply_manuscript_style()

    df  = pd.read_csv(long_csv)
    sdf = pd.read_csv(summary_csv)

    if n_list is None:
        n_list = sorted(df["n"].unique())

    for n in n_list:
        n_dir = os.path.join(figures_base_dir, f"n{n}")
        os.makedirs(n_dir, exist_ok=True)
        print(f"\n--- n={n}  ->  {n_dir} ---")
        plot_tau_effect_clean(df, n_dir, n=n)     # clean main
        plot_tau_effect(df, n_dir, n=n)           # full appendix grid
        plot_omega_sweep(df, n_dir, n=n)
        plot_time_to_tol_heatmap(sdf, n_dir, n=n)
        plot_tau_speedup_heatmap(sdf, n_dir, n=n)
