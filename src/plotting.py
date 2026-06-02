"""
Plotting routines for affine-invariant Gaussian gradient flow experiments.

All figures are saved as both PNG (dpi=200) and PDF.
Requires pandas for DataFrame filtering; matplotlib for rendering.
No seaborn dependency.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd


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

INIT_LABELS = {
    "mean_only":    "Mean only",
    "volume_high":  "Volume high",
    "volume_low":   "Volume low",
    "shape_only":   "Shape only",
    "mixed":        "Mixed",
}

INIT_ORDER = ["mean_only", "volume_high", "volume_low", "shape_only", "mixed"]

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
    """Plot with semilogy, skipping non-positive y values gracefully."""
    mask = np.asarray(y, dtype=float) > 0
    xs = np.asarray(x)[mask]
    ys = np.asarray(y)[mask]
    if len(xs):
        ax.semilogy(xs, ys, **kwargs)


# ---------------------------------------------------------------------------
# Figure 1: tau effect for omega=0.5, n=5
# ---------------------------------------------------------------------------

def plot_tau_effect(df, outdir, n=5, omega=0.5):
    """Fig 1: compare tau_negative / tau_zero / tau_positive.

    Layout: 5 rows (initializations) × 6 columns (metrics).
    Each subplot: three semilogy curves, one per tau type.
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
        figsize=(3.5 * ncols, 2.8 * nrows),
        sharey="col",
    )

    for row_i, init in enumerate(inits):
        for col_j, metric in enumerate(metrics):
            ax = axes[row_i, col_j]
            s_init = sub[sub["init_name"] == init]

            for tt in ["negative", "zero", "positive"]:
                s_tt = s_init[s_init["tau_type"] == tt].sort_values("time")
                if s_tt.empty:
                    continue
                _semilogy_safe(
                    ax, s_tt["time"], s_tt[metric],
                    color=TAU_COLORS[tt],
                    lw=1.6,
                    label=TAU_LABELS[tt],
                )

            ax.set_xlim(left=0)
            if row_i == nrows - 1:
                ax.set_xlabel("Time")
            if col_j == 0:
                ax.set_ylabel(INIT_LABELS[init], fontsize=9)
            if row_i == 0:
                ax.set_title(METRIC_LABELS[metric], fontsize=9)

    # Shared legend (top-right subplot)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles, labels,
            loc="upper right", fontsize=9,
            bbox_to_anchor=(1.0, 1.0),
            framealpha=0.9,
        )

    fig.suptitle(
        rf"$\tau$ effect  —  $n={n}$, $\omega={omega}$",
        y=1.01, fontsize=12,
    )
    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_tau_effect_omega_half_n{n}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Figure 2: omega sweep for tau=0, n=5
# ---------------------------------------------------------------------------

def plot_omega_sweep(df, outdir, n=5):
    """Fig 2: normalised energy vs time for all omega values (tau=0).

    Layout: 1 row × 5 columns (one per initialisation).
    Each subplot: one curve per omega value.
    """
    sub = df[(df["n"] == n) & (df["tau_type"] == "zero")].copy()
    if sub.empty:
        print(f"[plot_omega_sweep] No data for n={n}, tau=0. Skipping.")
        return

    omegas = sorted(sub["omega"].unique())
    cmap = plt.cm.get_cmap("plasma", len(omegas))
    omega_colors = {o: cmap(i) for i, o in enumerate(omegas)}

    inits = INIT_ORDER
    fig, axes = plt.subplots(1, len(inits), figsize=(3.8 * len(inits), 3.4))
    if len(inits) == 1:
        axes = [axes]

    for col_j, init in enumerate(inits):
        ax = axes[col_j]
        s_init = sub[sub["init_name"] == init]

        for omega in omegas:
            s_o = s_init[np.isclose(s_init["omega"], omega)].sort_values("time")
            if s_o.empty:
                continue
            _semilogy_safe(
                ax, s_o["time"], s_o["norm_energy"],
                color=omega_colors[omega],
                lw=1.6,
                label=rf"$\omega={omega}$",
            )

        ax.set_xlim(left=0)
        ax.set_xlabel("Time")
        if col_j == 0:
            ax.set_ylabel("Normalised KL energy")
        ax.set_title(INIT_LABELS[init], fontsize=10)
        ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(rf"$\omega$ sweep  —  $n={n}$, $\tau=0$", fontsize=12)
    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_omega_sweep_tau_zero_n{n}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Figure 3: time-to-tolerance heatmap (n=5, tau=0)
# ---------------------------------------------------------------------------

def plot_time_to_tol_heatmap(summary_df, outdir, n=5, tol_col="time_to_1e_minus_4"):
    """Fig 3: heatmap of time-to-1e-4 for rows=init, cols=omega (tau=0)."""
    sub = summary_df[(summary_df["n"] == n) & (summary_df["tau_type"] == "zero")].copy()
    if sub.empty:
        print(f"[plot_time_to_tol_heatmap] No data for n={n}, tau=0. Skipping.")
        return

    omegas = sorted(sub["omega"].unique())
    inits  = INIT_ORDER

    # Build matrix; replace inf with a capped value for display
    T_max = float(sub["T"].max())
    mat = np.full((len(inits), len(omegas)), np.nan)
    for r_i, init in enumerate(inits):
        for c_j, omega in enumerate(omegas):
            rows = sub[(sub["init_name"] == init) & np.isclose(sub["omega"], omega)]
            if not rows.empty:
                val = float(rows[tol_col].iloc[0])
                mat[r_i, c_j] = val if np.isfinite(val) else T_max * 1.1

    fig, ax = plt.subplots(figsize=(7, 4))
    im = ax.imshow(mat, aspect="auto", cmap="viridis_r", origin="upper")
    ax.set_xticks(range(len(omegas)))
    ax.set_xticklabels([str(o) for o in omegas])
    ax.set_yticks(range(len(inits)))
    ax.set_yticklabels([INIT_LABELS[i] for i in inits])
    ax.set_xlabel(r"$\omega$")
    ax.set_title(
        rf"Time to $E/E_0 \leq 10^{{-4}}$  —  $n={n}$, $\tau=0$",
        fontsize=11,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Time")

    # Annotate cells
    for r_i in range(len(inits)):
        for c_j in range(len(omegas)):
            v = mat[r_i, c_j]
            txt = f"{v:.1f}" if np.isfinite(v) else "∞"
            ax.text(c_j, r_i, txt, ha="center", va="center",
                    fontsize=8, color="white")

    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_time_to_tol_heatmap_n{n}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Figure 4: tau speedup heatmap (n=5)
# ---------------------------------------------------------------------------

def plot_tau_speedup_heatmap(summary_df, outdir, n=5, tol_col="time_to_1e_minus_4"):
    """Fig 4: ratio T(tau) / T(tau_zero) for tau_neg and tau_pos.

    omega in {1/4, 1/2, 1}; rows = initialisation; two heatmaps side by side.
    Ratio > 1 => tau slows convergence; ratio < 1 => tau accelerates it.
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

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)

    for ax, tau_side in zip(axes, ["negative", "positive"]):
        mat = np.full((len(inits), len(target_omegas)), np.nan)

        for r_i, init in enumerate(inits):
            for c_j, omega in enumerate(target_omegas):
                r_zero = sub[
                    (sub["init_name"] == init) &
                    np.isclose(sub["omega"], omega) &
                    (sub["tau_type"] == "zero")
                ]
                r_tau = sub[
                    (sub["init_name"] == init) &
                    np.isclose(sub["omega"], omega) &
                    (sub["tau_type"] == tau_side)
                ]
                if r_zero.empty or r_tau.empty:
                    continue
                t_zero = float(r_zero[tol_col].iloc[0])
                t_tau  = float(r_tau[tol_col].iloc[0])
                # Both finite: compute ratio; inf/inf -> nan; inf/finite -> cap
                if not np.isfinite(t_zero) and not np.isfinite(t_tau):
                    mat[r_i, c_j] = 1.0
                elif not np.isfinite(t_zero):
                    mat[r_i, c_j] = 0.0
                elif not np.isfinite(t_tau):
                    mat[r_i, c_j] = 3.0   # display cap
                else:
                    mat[r_i, c_j] = t_tau / t_zero

        vmin, vmax = 0.2, 3.0
        im = ax.imshow(mat, aspect="auto", cmap="RdYlGn_r",
                       vmin=vmin, vmax=vmax, origin="upper")
        ax.set_xticks(range(len(target_omegas)))
        ax.set_xticklabels([str(o) for o in target_omegas])
        ax.set_yticks(range(len(inits)))
        ax.set_yticklabels([INIT_LABELS[i] for i in inits])
        ax.set_xlabel(r"$\omega$")
        tau_sign = "-" if tau_side == "negative" else "+"
        ax.set_title(rf"$\tau {tau_sign}$ speedup ratio")

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(r"$T(\tau)/T(\tau{=}0)$")

        for r_i in range(len(inits)):
            for c_j in range(len(target_omegas)):
                v = mat[r_i, c_j]
                txt = f"{v:.2f}" if np.isfinite(v) else "?"
                color = "black" if 0.8 < v < 2.0 else "white"
                ax.text(c_j, r_i, txt, ha="center", va="center",
                        fontsize=8, color=color)

    fig.suptitle(
        rf"$\tau$ speedup ratio $T(\tau)/T(\tau{{=}}0)$  —  $n={n}$",
        fontsize=11,
    )
    fig.tight_layout()
    stem = os.path.join(outdir, f"fig_tau_speedup_heatmap_n{n}")
    _save(fig, stem)
    print(f"Saved {stem}.png / .pdf")


# ---------------------------------------------------------------------------
# Convenience wrapper: generate all four figures
# ---------------------------------------------------------------------------

def generate_all_figures(long_csv, summary_csv, figures_base_dir, n_list=None):
    """Load CSVs and produce all four manuscript figures for every n in n_list.

    Each n gets its own subdirectory: <figures_base_dir>/n{N}/.

    Args:
        long_csv         : path to results_long.csv
        summary_csv      : path to summary.csv
        figures_base_dir : parent directory; per-n subdirs are created inside it
        n_list           : list of n values to plot (default: all n in the CSV)
    """
    df  = pd.read_csv(long_csv)
    sdf = pd.read_csv(summary_csv)

    if n_list is None:
        n_list = sorted(df["n"].unique())

    for n in n_list:
        n_dir = os.path.join(figures_base_dir, f"n{n}")
        os.makedirs(n_dir, exist_ok=True)
        print(f"\n--- n={n}  ->  {n_dir} ---")
        plot_tau_effect(df, n_dir, n=n)
        plot_omega_sweep(df, n_dir, n=n)
        plot_time_to_tol_heatmap(sdf, n_dir, n=n)
        plot_tau_speedup_heatmap(sdf, n_dir, n=n)
