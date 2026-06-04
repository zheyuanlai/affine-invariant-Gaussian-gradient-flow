"""Plot natural-gradient local-rate results.

Figures (PNG + PDF), faceted by potential family:
  1. Lambda_hat vs N_theta
  2. Lambda_hat / (1 + log kappa) vs N_theta
  3. 1 / gamma_loc vs N_theta
  4. (1 / gamma_loc) / (1 + log kappa) vs N_theta
  5. Lambda_hat vs (1 + log kappa)
  6. Flow validation: log R^2 vs t with the predicted slope -2 gamma_loc overlaid

Usage:
    python scripts/natural_gradient_local_rate/plot_results.py \
        --input outputs/natural_gradient_local_rate \
        --outdir outputs/natural_gradient_local_rate/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))

from src.common.plotting_style import apply_style, save_figure  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def _read_csv(path):
    return pd.read_csv(path) if os.path.exists(path) else None


def _load(input_dir):
    op = _read_csv(os.path.join(input_dir, "operator_grid", "results_long.csv"))
    lr = _read_csv(os.path.join(input_dir, "linearized_rate_grid", "results_long.csv"))
    flow_traj = _read_csv(os.path.join(input_dir, "flow_validation", "trajectories.csv"))
    flow_sum = _read_csv(os.path.join(input_dir, "flow_validation", "summary.csv"))
    return op, lr, flow_traj, flow_sum


def _facet_vs(df, xcol, ycol, linecol, title, xlabel, ylabel, outpath):
    """One subplot per family; one line per value of ``linecol`` (mean over seeds)."""
    if df is None or ycol not in df.columns or xcol not in df.columns:
        print(f"  [skip] {os.path.basename(outpath)} (missing {ycol}/{xcol})")
        return
    families = sorted(df["family"].unique())
    n = len(families)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.6), squeeze=False)
    cmap = plt.get_cmap("viridis")
    line_vals = sorted(df[linecol].unique())
    for ax, fam in zip(axes[0], families):
        sub = df[df["family"] == fam]
        for k, lv in enumerate(line_vals):
            s = sub[sub[linecol] == lv]
            g = s.groupby(xcol)[ycol].mean().reset_index().sort_values(xcol)
            if g.empty:
                continue
            color = cmap(k / max(1, len(line_vals) - 1))
            ax.plot(g[xcol], g[ycol], marker="o", color=color,
                    label=f"{linecol.split('_')[0]}={lv:g}")
        ax.set_title(fam)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(title=linecol.split("_")[0], fontsize=8)
    fig.suptitle(title)
    save_figure(fig, outpath)
    plt.close(fig)
    print(f"  [ok]   {os.path.basename(outpath)}")


def _plot_flow(flow_traj, flow_sum, outpath, max_curves=8):
    if flow_traj is None or flow_sum is None:
        print(f"  [skip] {os.path.basename(outpath)} (no flow data)")
        return
    # one representative (epsilon, Delta_t) per key: smallest epsilon, smallest dt
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    keys = list(flow_traj["key"].unique())[:max_curves]
    cmap = plt.get_cmap("tab10")
    for i, key in enumerate(keys):
        sub = flow_traj[flow_traj["key"] == key]
        eps = sub["epsilon"].min()
        dt = sub["Delta_t"].min()
        curve = sub[(sub["epsilon"] == eps) & (sub["Delta_t"] == dt)].sort_values("t")
        if curve.empty:
            continue
        color = cmap(i % 10)
        ax.plot(curve["t"], curve["log_R2"], color=color, lw=1.4, label=key)
        # predicted slope -2 gamma_loc through the first fitted point
        row = flow_sum[(flow_sum["key"] == key)] if "key" in flow_sum.columns else None
        gl = None
        if row is not None and not row.empty and "gamma_loc" in row.columns:
            gl = float(row["gamma_loc"].iloc[0])
        if gl is not None:
            t0 = curve["t"].iloc[0]
            y0 = curve["log_R2"].iloc[0]
            ax.plot(curve["t"], y0 - 2.0 * gl * (curve["t"] - t0),
                    color=color, ls="--", lw=0.9)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\log R^2$")
    ax.set_title(r"Flow decay (solid) vs predicted slope $-2\,\gamma_{loc}$ (dashed)")
    ax.legend(fontsize=7, ncol=2)
    save_figure(fig, outpath)
    plt.close(fig)
    print(f"  [ok]   {os.path.basename(outpath)}")


def main():
    p = argparse.ArgumentParser(description="Plot natural-gradient local-rate results.")
    p.add_argument("--input", default="outputs/natural_gradient_local_rate")
    p.add_argument("--outdir", default=None)
    args = p.parse_args()
    outdir = args.outdir or os.path.join(args.input, "figures")
    os.makedirs(outdir, exist_ok=True)
    apply_style()

    op, lr, flow_traj, flow_sum = _load(args.input)
    # Prefer the operator grid for Lambda_hat; fall back to the rate grid.
    lam_df = op if (op is not None and "Lambda_hat" in op.columns) else lr

    print(f"Writing figures to {outdir}")
    _facet_vs(lam_df, "N_theta", "Lambda_hat", "kappa_target",
              "Operator norm vs dimension", r"$N_\theta$", r"$\hat\Lambda$",
              os.path.join(outdir, "fig1_lambda_vs_Ntheta"))
    _facet_vs(lam_df, "N_theta", "lambda_over_logkappa", "kappa_target",
              "Normalized operator norm vs dimension", r"$N_\theta$",
              r"$\hat\Lambda / (1+\log\kappa)$",
              os.path.join(outdir, "fig2_lambda_over_logkappa_vs_Ntheta"))
    _facet_vs(lr, "N_theta", "inverse_gamma_loc", "kappa_target",
              "Inverse local rate vs dimension", r"$N_\theta$",
              r"$1/\gamma_{loc}$",
              os.path.join(outdir, "fig3_inv_gamma_vs_Ntheta"))
    _facet_vs(lr, "N_theta", "inverse_gamma_over_logkappa", "kappa_target",
              "Normalized inverse local rate vs dimension", r"$N_\theta$",
              r"$1/\gamma_{loc}/(1+\log\kappa)$",
              os.path.join(outdir, "fig4_inv_gamma_over_logkappa_vs_Ntheta"))
    _facet_vs(lam_df, "log_kappa_factor", "Lambda_hat", "N_theta",
              "Operator norm vs conditioning", r"$1+\log\kappa$", r"$\hat\Lambda$",
              os.path.join(outdir, "fig5_lambda_vs_logkappa"))
    _plot_flow(flow_traj, flow_sum, os.path.join(outdir, "fig6_flow_logR2"))
    print("Done.")


if __name__ == "__main__":
    main()
