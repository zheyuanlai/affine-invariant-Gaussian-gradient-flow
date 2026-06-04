"""Estimator-diagnostic figures for the natural-gradient local-rate experiment.

Reads the operator-grid, linearized-rate-grid and sample-size-scaling results
and writes PNG+PDF figures (matplotlib only, no seaborn):

  1. Raw vs symmetrized Lambda_hat vs N_theta.
  2. H self-adjointness error (raw vs symmetrized) vs N_theta.
  3. Separable: full_sym vs diagonal vs exact Lambda vs N_theta.
  4. Lambda_hat vs M_mc (log x), with the exact benchmark overlaid.
  5. gamma_loc vs M_mc (log x).
  6. self_adjoint_error_L_star vs M_mc (log x).
  7. full_over_diag vs N_theta for separable controls.

Usage:
    python scripts/natural_gradient_local_rate/plot_estimator_diagnostics.py \
        --input outputs/natural_gradient_local_rate \
        --outdir outputs/natural_gradient_local_rate/figures/estimator_diagnostics
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

SEPARABLE_FAMILIES = ("gaussian", "separable")


def _read_csv(path):
    return pd.read_csv(path) if os.path.exists(path) else None


def _grid_df(input_dir):
    """Prefer the linearized-rate grid (full suite incl. gamma), else operator grid."""
    lr = _read_csv(os.path.join(input_dir, "linearized_rate_grid", "results_long.csv"))
    op = _read_csv(os.path.join(input_dir, "operator_grid", "results_long.csv"))
    if lr is not None and "Lambda_hat_full_sym" in lr.columns:
        return lr
    return op


def _has(df, *cols):
    return df is not None and all(c in df.columns for c in cols) and len(df) > 0


def _mean_over(df, xcol, ycol, by=("potential_family",)):
    """Group by ``by + (xcol,)`` and average ``ycol`` (NaNs ignored)."""
    g = df.dropna(subset=[ycol]).groupby(list(by) + [xcol])[ycol].mean().reset_index()
    return g.sort_values(xcol)


def _families(df):
    return sorted(df["potential_family"].unique())


# --- grid figures (vs N_theta) -------------------------------------------------

def plot_raw_vs_sym(df, outpath):
    if not _has(df, "potential_family", "N_theta",
                "Lambda_hat_raw_forward", "Lambda_hat_full_sym"):
        print(f"  [skip] {os.path.basename(outpath)} (missing columns)")
        return
    fams = _families(df)
    fig, axes = plt.subplots(1, len(fams), figsize=(4.2 * len(fams), 3.6), squeeze=False)
    for ax, fam in zip(axes[0], fams):
        sub = df[df["potential_family"] == fam]
        for col, lab, mk in [("Lambda_hat_raw_forward", "raw_forward", "s"),
                             ("Lambda_hat_full_sym", "symmetrized", "o")]:
            g = _mean_over(sub, "N_theta", col, by=())
            if not g.empty:
                ax.plot(g["N_theta"], g[col], marker=mk, label=lab)
        ax.set_title(fam)
        ax.set_xlabel(r"$N_\theta$")
        ax.set_ylabel(r"$\hat\Lambda$")
        ax.legend(fontsize=8)
    fig.suptitle("Raw forward vs symmetrized operator norm")
    save_figure(fig, outpath)
    plt.close(fig)
    print(f"  [ok]   {os.path.basename(outpath)}")


def plot_self_adjoint_H(df, outpath):
    if not _has(df, "self_adjoint_error_H_raw", "self_adjoint_error_H_sym"):
        print(f"  [skip] {os.path.basename(outpath)} (missing columns)")
        return
    fams = _families(df)
    fig, axes = plt.subplots(1, len(fams), figsize=(4.2 * len(fams), 3.6), squeeze=False)
    for ax, fam in zip(axes[0], fams):
        sub = df[df["potential_family"] == fam]
        for col, lab, mk in [("self_adjoint_error_H_raw", "raw_forward", "s"),
                             ("self_adjoint_error_H_sym", "symmetrized", "o")]:
            g = _mean_over(sub, "N_theta", col, by=())
            if not g.empty:
                ax.plot(g["N_theta"], np.maximum(g[col], 1e-18), marker=mk, label=lab)
        ax.set_yscale("log")
        ax.set_title(fam)
        ax.set_xlabel(r"$N_\theta$")
        ax.set_ylabel("relative self-adjointness error")
        ax.legend(fontsize=8)
    fig.suptitle("H self-adjointness error: raw forward vs symmetrized")
    save_figure(fig, outpath)
    plt.close(fig)
    print(f"  [ok]   {os.path.basename(outpath)}")


def plot_separable_full_diag_exact(df, outpath):
    if not _has(df, "Lambda_hat_full_sym", "Lambda_hat_diag"):
        print(f"  [skip] {os.path.basename(outpath)} (missing columns)")
        return
    sub = df[df["potential_family"].isin(SEPARABLE_FAMILIES)]
    if sub.empty:
        print(f"  [skip] {os.path.basename(outpath)} (no separable controls)")
        return
    fams = sorted(sub["potential_family"].unique())
    fig, axes = plt.subplots(1, len(fams), figsize=(4.2 * len(fams), 3.6), squeeze=False)
    for ax, fam in zip(axes[0], fams):
        s = sub[sub["potential_family"] == fam]
        for col, lab, mk in [("Lambda_hat_full_sym", "full_sym", "o"),
                             ("Lambda_hat_diag", "diagonal", "s"),
                             ("Lambda_hat_separable_exact", "exact (quad)", "^")]:
            if col not in s.columns:
                continue
            g = _mean_over(s, "N_theta", col, by=())
            if not g.empty:
                ls = "--" if "exact" in lab else "-"
                ax.plot(g["N_theta"], g[col], marker=mk, ls=ls, label=lab)
        ax.set_title(fam)
        ax.set_xlabel(r"$N_\theta$")
        ax.set_ylabel(r"$\hat\Lambda$")
        ax.legend(fontsize=8)
    fig.suptitle("Separable controls: full vs diagonal vs exact operator norm")
    save_figure(fig, outpath)
    plt.close(fig)
    print(f"  [ok]   {os.path.basename(outpath)}")


def plot_full_over_diag_vs_N(df, outpath):
    if not _has(df, "full_over_diag"):
        print(f"  [skip] {os.path.basename(outpath)} (missing columns)")
        return
    sub = df[df["potential_family"].isin(SEPARABLE_FAMILIES)]
    if sub.empty:
        print(f"  [skip] {os.path.basename(outpath)} (no separable controls)")
        return
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for fam in sorted(sub["potential_family"].unique()):
        g = _mean_over(sub[sub["potential_family"] == fam], "N_theta", "full_over_diag", by=())
        if not g.empty:
            ax.plot(g["N_theta"], g["full_over_diag"], marker="o", label=fam)
    ax.axhline(1.0, color="0.6", lw=0.9, ls=":")
    ax.set_xlabel(r"$N_\theta$")
    ax.set_ylabel(r"$\hat\Lambda_{\mathrm{full}} / \hat\Lambda_{\mathrm{diag}}$")
    ax.set_title("Full / diagonal operator-norm ratio (separable controls)")
    ax.legend(fontsize=8)
    save_figure(fig, outpath)
    plt.close(fig)
    print(f"  [ok]   {os.path.basename(outpath)}")


# --- sample-size-scaling figures (vs M_mc) ------------------------------------

def _scaling_facets(df, ycol, title, ylabel, outpath, logy=False, overlay_exact=False):
    if not _has(df, "M_mc", ycol):
        print(f"  [skip] {os.path.basename(outpath)} (missing {ycol})")
        return
    fams = _families(df)
    fig, axes = plt.subplots(1, len(fams), figsize=(4.4 * len(fams), 3.8), squeeze=False)
    cmap = plt.get_cmap("viridis")
    for ax, fam in zip(axes[0], fams):
        sub = df[df["potential_family"] == fam]
        Ns = sorted(sub["N_theta"].unique())
        for k, N in enumerate(Ns):
            s = sub[sub["N_theta"] == N]
            g = _mean_over(s, "M_mc", ycol, by=())
            if g.empty:
                continue
            color = cmap(k / max(1, len(Ns) - 1))
            ax.plot(g["M_mc"], np.maximum(g[ycol], 1e-18) if logy else g[ycol],
                    marker="o", color=color, label=fr"$N_\theta={int(N)}$")
            if overlay_exact and "Lambda_hat_separable_exact" in s.columns:
                ge = _mean_over(s, "M_mc", "Lambda_hat_separable_exact", by=())
                if not ge.empty and np.isfinite(ge["Lambda_hat_separable_exact"]).any():
                    ax.axhline(float(ge["Lambda_hat_separable_exact"].iloc[-1]),
                               color=color, ls="--", lw=0.9)
        ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
        ax.set_title(fam)
        ax.set_xlabel(r"$M_{\mathrm{mc}}$")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=8)
    fig.suptitle(title)
    save_figure(fig, outpath)
    plt.close(fig)
    print(f"  [ok]   {os.path.basename(outpath)}")


def main():
    p = argparse.ArgumentParser(description="Estimator-diagnostic figures.")
    p.add_argument("--input", default="outputs/natural_gradient_local_rate")
    p.add_argument("--outdir", default=None)
    args = p.parse_args()
    outdir = args.outdir or os.path.join(args.input, "figures", "estimator_diagnostics")
    os.makedirs(outdir, exist_ok=True)
    apply_style()

    grid = _grid_df(args.input)
    scaling = _read_csv(os.path.join(args.input, "sample_size_scaling", "results_long.csv"))

    print(f"Writing estimator-diagnostic figures to {outdir}")
    plot_raw_vs_sym(grid, os.path.join(outdir, "diag1_raw_vs_sym_lambda"))
    plot_self_adjoint_H(grid, os.path.join(outdir, "diag2_self_adjoint_error_H"))
    plot_separable_full_diag_exact(grid, os.path.join(outdir, "diag3_separable_full_diag_exact"))
    _scaling_facets(scaling, "Lambda_hat_full_sym",
                    "Operator norm vs sample size", r"$\hat\Lambda_{\mathrm{full\ sym}}$",
                    os.path.join(outdir, "diag4_lambda_vs_Mmc"), overlay_exact=True)
    _scaling_facets(scaling, "gamma_loc",
                    "Local rate vs sample size", r"$\gamma_{\mathrm{loc}}$",
                    os.path.join(outdir, "diag5_gamma_vs_Mmc"))
    _scaling_facets(scaling, "self_adjoint_error_L_star",
                    "L_star self-adjointness error vs sample size",
                    "relative self-adjointness error",
                    os.path.join(outdir, "diag6_self_adjoint_L_vs_Mmc"), logy=True)
    plot_full_over_diag_vs_N(grid, os.path.join(outdir, "diag7_full_over_diag_vs_Ntheta"))
    print("Done.")


if __name__ == "__main__":
    main()
