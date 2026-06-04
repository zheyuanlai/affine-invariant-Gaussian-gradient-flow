"""Estimator-diagnostic figures for natural-gradient local-rate experiments.

The default plots treat ``Lambda_hat_full_sym`` as a finite-sample diagnostic,
not as direct evidence of dimension dependence. Baseline-corrected columns are
used when available; otherwise they are computed in memory from
``results_long.csv``.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))

from src.common.plotting_style import apply_style, save_figure  # noqa: E402
from src.natural_gradient_local_rate.baseline_postprocessing import (  # noqa: E402
    add_baseline_corrections,
)
from src.natural_gradient_local_rate.scaling_diagnostics import (  # noqa: E402
    add_noise_scale_columns, fit_scaling_diagnostics,
)
import matplotlib.pyplot as plt  # noqa: E402


def _read_csv(path):
    return pd.read_csv(path) if os.path.exists(path) else None


def _ensure_corrected(df):
    if df is None:
        return None
    if "full_sym_minus_gaussian" in df.columns and "sqrt_p_over_M" in df.columns:
        return df
    return add_baseline_corrections(add_noise_scale_columns(df))


def _load_scaling(input_dir):
    base = os.path.join(input_dir, "sample_size_scaling")
    corrected = _read_csv(os.path.join(base, "results_long_with_baselines.csv"))
    if corrected is not None:
        return _ensure_corrected(corrected)
    raw = _read_csv(os.path.join(base, "results_long.csv"))
    return _ensure_corrected(raw)


def _load_grid(input_dir):
    for stage in ("linearized_rate_grid", "operator_grid"):
        df = _read_csv(os.path.join(input_dir, stage, "results_long.csv"))
        if df is not None and len(df):
            return _ensure_corrected(df)
    return None


def _has(df, *cols):
    return df is not None and len(df) > 0 and all(c in df.columns for c in cols)


def _num(df, col):
    return pd.to_numeric(df[col], errors="coerce")


def _families(df):
    return sorted(df["potential_family"].dropna().astype(str).unique())


def _mean(df, by, ycol):
    if not _has(df, ycol):
        return pd.DataFrame()
    cols = list(by) + [ycol]
    g = df[cols].copy()
    g[ycol] = pd.to_numeric(g[ycol], errors="coerce")
    return g.dropna(subset=[ycol]).groupby(list(by), dropna=False)[ycol].mean().reset_index()


def _save(fig, outpath):
    save_figure(fig, outpath)
    plt.close(fig)
    print(f"  [ok] {os.path.basename(outpath)}")


def _line_by_N(ax, df, ycol, *, label_prefix="", linestyle="-", marker="o"):
    if not _has(df, "M_mc", "N_theta", ycol):
        return
    for N in sorted(df["N_theta"].dropna().unique()):
        sub = df[df["N_theta"] == N]
        g = _mean(sub, ["M_mc"], ycol).sort_values("M_mc")
        if not g.empty:
            label = f"{label_prefix}N={int(N)}"
            ax.plot(g["M_mc"], g[ycol], marker=marker, ls=linestyle, label=label)


def plot_lambda_vs_M(df, outpath):
    if not _has(df, "M_mc", "N_theta", "Lambda_hat_full_sym"):
        return
    fams = _families(df)
    fig, axes = plt.subplots(1, len(fams), figsize=(4.2 * len(fams), 3.7), squeeze=False)
    for ax, fam in zip(axes[0], fams):
        _line_by_N(ax, df[df["potential_family"] == fam], "Lambda_hat_full_sym")
        ax.set_xscale("log")
        ax.set_title(fam)
        ax.set_xlabel("M_mc")
        ax.set_ylabel("Lambda full sym")
        ax.legend(fontsize=7)
    fig.suptitle("Full operator estimates vs sample size (finite-sample diagnostic)")
    _save(fig, outpath)


def plot_lambda_with_baselines(df, outpath):
    if not _has(df, "M_mc", "N_theta", "Lambda_hat_full_sym"):
        return
    fams = _families(df)
    fig, axes = plt.subplots(1, len(fams), figsize=(4.5 * len(fams), 3.8), squeeze=False)
    for ax, fam in zip(axes[0], fams):
        sub = df[df["potential_family"] == fam]
        _line_by_N(ax, sub, "Lambda_hat_full_sym", label_prefix="full ")
        if "Lambda_hat_gaussian_baseline" in sub.columns:
            _line_by_N(ax, sub, "Lambda_hat_gaussian_baseline",
                       label_prefix="gaussian ", linestyle="--", marker=".")
        if "Lambda_hat_separable_baseline" in sub.columns:
            _line_by_N(ax, sub, "Lambda_hat_separable_baseline",
                       label_prefix="separable ", linestyle=":", marker=".")
        ax.set_xscale("log")
        ax.set_title(fam)
        ax.set_xlabel("M_mc")
        ax.set_ylabel("Lambda full sym")
        ax.legend(fontsize=6)
    fig.suptitle("Full operator with matched Gaussian/separable baselines")
    _save(fig, outpath)


def plot_separable_full_diag_exact(df, outpath):
    if not _has(df, "M_mc", "Lambda_hat_full_sym", "Lambda_hat_diag"):
        return
    sub = df[df["potential_family"].astype(str) == "separable"]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    for col, label, style in [
        ("Lambda_hat_full_sym", "full sym", "-"),
        ("Lambda_hat_diag", "diagonal", "--"),
        ("Lambda_hat_separable_exact", "exact", ":"),
    ]:
        g = _mean(sub, ["N_theta", "M_mc"], col)
        for N in sorted(g["N_theta"].unique()):
            s = g[g["N_theta"] == N].sort_values("M_mc")
            ax.plot(s["M_mc"], s[col], ls=style, marker="o", label=f"{label} N={int(N)}")
    ax.set_xscale("log")
    ax.set_xlabel("M_mc")
    ax.set_ylabel("Lambda")
    ax.set_title("Separable controls: full, diagonal, exact vs sample size")
    ax.legend(fontsize=6, ncol=2)
    _save(fig, outpath)


def plot_gaussian_true_zero(df, outpath):
    sub = df[df["potential_family"].astype(str) == "gaussian"] if df is not None else pd.DataFrame()
    if sub.empty or not _has(sub, "Lambda_hat_full_sym"):
        return
    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    _line_by_N(ax, sub, "Lambda_hat_full_sym")
    ax.axhline(0.0, color="0.25", lw=1.0, ls="--", label="true Lambda=0")
    ax.set_xscale("log")
    ax.set_xlabel("M_mc")
    ax.set_ylabel("Lambda full sym")
    ax.set_title("Gaussian baseline: finite-sample full operator noise")
    ax.legend(fontsize=7)
    _save(fig, outpath)


def plot_random_gap(df, col, ylabel, title, outpath):
    sub = df[df["potential_family"].astype(str) == "random_feature"] if df is not None else pd.DataFrame()
    if sub.empty or not _has(sub, col):
        return
    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    _line_by_N(ax, sub, col)
    ax.axhline(0.0, color="0.5", lw=0.9, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("M_mc")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=7)
    _save(fig, outpath)


def plot_gamma_gaussian_true(df, outpath):
    if not _has(df, "M_mc", "gamma_loc"):
        return
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for fam in _families(df):
        g = _mean(df[df["potential_family"] == fam], ["M_mc"], "gamma_loc").sort_values("M_mc")
        if not g.empty:
            ax.plot(g["M_mc"], g["gamma_loc"], marker="o", label=fam)
    ax.axhline(1.0, color="0.25", lw=1.0, ls="--", label="Gaussian true gamma=1")
    ax.set_xscale("log")
    ax.set_xlabel("M_mc")
    ax.set_ylabel("gamma_loc")
    ax.set_title("Local rate vs sample size")
    ax.legend(fontsize=7)
    _save(fig, outpath)


def plot_gamma_compare(df, outpath):
    if not _has(df, "M_mc", "N_theta", "gamma_loc"):
        return
    Ns = sorted(df["N_theta"].dropna().unique())
    fig, axes = plt.subplots(1, len(Ns), figsize=(4.2 * len(Ns), 3.7), squeeze=False)
    for ax, N in zip(axes[0], Ns):
        subN = df[df["N_theta"] == N]
        for fam in _families(subN):
            g = _mean(subN[subN["potential_family"] == fam], ["M_mc"], "gamma_loc").sort_values("M_mc")
            if not g.empty:
                ax.plot(g["M_mc"], g["gamma_loc"], marker="o", label=fam)
        ax.axhline(1.0, color="0.6", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_title(f"N={int(N)}")
        ax.set_xlabel("M_mc")
        ax.set_ylabel("gamma_loc")
        ax.legend(fontsize=7)
    fig.suptitle("gamma_loc comparison across baselines and potentials")
    _save(fig, outpath)


def plot_self_adjoint_errors(df, outpath):
    if not _has(df, "M_mc"):
        return
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for col, label in [
        ("self_adjoint_error_H_sym", "H_sym"),
        ("self_adjoint_error_L_star", "L_star"),
    ]:
        if col not in df.columns:
            continue
        g = _mean(df, ["M_mc"], col).sort_values("M_mc")
        if not g.empty:
            ax.plot(g["M_mc"], np.maximum(g[col], 1e-18), marker="o", label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("M_mc")
    ax.set_ylabel("relative self-adjointness error")
    ax.set_title("Self-adjointness checks vs sample size")
    ax.legend(fontsize=7)
    _save(fig, outpath)


def plot_gap_vs_noise_scale(df, outpath):
    if not _has(df, "sqrt_p_over_M"):
        return
    work = df.copy()
    if "full_sym_minus_separable" in work.columns:
        work["gap_for_plot"] = _num(work, "full_sym_minus_separable")
    elif "full_sym_minus_gaussian" in work.columns:
        work["gap_for_plot"] = _num(work, "full_sym_minus_gaussian")
    else:
        return
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for fam in _families(work):
        sub = work[work["potential_family"] == fam]
        g = _mean(sub, ["sqrt_p_over_M"], "gap_for_plot").sort_values("sqrt_p_over_M")
        if not g.empty:
            ax.plot(g["sqrt_p_over_M"], g["gap_for_plot"], marker="o", label=fam)
    ax.axhline(0.0, color="0.5", lw=0.9, ls=":")
    ax.set_xlabel("sqrt(p / M_mc)")
    ax.set_ylabel("full sym minus baseline")
    ax.set_title("Baseline-corrected gap vs finite-sample noise scale")
    ax.legend(fontsize=7)
    _save(fig, outpath)


def plot_convergence_fit(df, outpath):
    if df is None or df.empty:
        return
    fits = fit_scaling_diagnostics(df)
    if fits.empty:
        return
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for fam in sorted(fits["potential_family"].unique()):
        sub = fits[fits["potential_family"] == fam]
        g = sub.groupby("N_theta", dropna=False)["Lambda_inf_fit"].mean().reset_index()
        ax.plot(g["N_theta"], g["Lambda_inf_fit"], marker="o", label=fam)
    ax.axhline(0.0, color="0.5", lw=0.9, ls=":")
    ax.set_xlabel("N_theta")
    ax.set_ylabel("fitted Lambda_inf")
    ax.set_title("Convergence fit: Lambda_inf from Lambda = Lambda_inf + c sqrt(p/M)")
    ax.legend(fontsize=7)
    _save(fig, outpath)


def plot_gap_heatmap(df, outpath):
    if df is None or df.empty:
        return
    candidates = [
        ("separable", "full_sym_minus_separable_exact", "separable full minus exact"),
        ("random_feature", "full_sym_minus_separable", "random_feature full minus separable"),
    ]
    panels = []
    for fam, col, title in candidates:
        sub = df[df["potential_family"].astype(str) == fam]
        if not sub.empty and col in sub.columns:
            pivot = sub.pivot_table(index="N_theta", columns="M_mc", values=col, aggfunc="mean")
            if not pivot.empty:
                panels.append((pivot, title))
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(5.0 * len(panels), 3.8), squeeze=False)
    for ax, (pivot, title) in zip(axes[0], panels):
        im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", origin="lower")
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels([str(int(c)) for c in pivot.columns], rotation=45, ha="right")
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels([str(int(i)) for i in pivot.index])
        ax.set_xlabel("M_mc")
        ax.set_ylabel("N_theta")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Baseline-corrected full-operator gaps")
    _save(fig, outpath)


def plot_full_minus_diag_vs_N(df, outpath):
    if not _has(df, "N_theta", "full_sym_minus_diag"):
        return
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for fam in _families(df):
        g = _mean(df[df["potential_family"] == fam], ["N_theta"], "full_sym_minus_diag").sort_values("N_theta")
        if not g.empty:
            ax.plot(g["N_theta"], g["full_sym_minus_diag"], marker="o", label=fam)
    ax.axhline(0.0, color="0.5", lw=0.9, ls=":")
    ax.set_xlabel("N_theta")
    ax.set_ylabel("full sym minus diagonal")
    ax.set_title("Full minus diagonal operator diagnostic")
    ax.legend(fontsize=7)
    _save(fig, outpath)


def plot_deprecated_ratio(df, outpath):
    if not _has(df, "N_theta", "full_over_diag"):
        return
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for fam in _families(df):
        g = _mean(df[df["potential_family"] == fam], ["N_theta"], "full_over_diag").sort_values("N_theta")
        if not g.empty:
            ax.plot(g["N_theta"], g["full_over_diag"], marker="o", label=fam)
    ax.axhline(1.0, color="0.6", lw=0.9, ls=":")
    ax.set_xlabel("N_theta")
    ax.set_ylabel("full / diagonal")
    ax.set_title("Deprecated diagnostic: ratio unstable when diagonal is near zero")
    ax.legend(fontsize=7)
    _save(fig, outpath)


def main():
    p = argparse.ArgumentParser(description="Estimator-diagnostic figures.")
    p.add_argument("--input", default="outputs/natural_gradient_local_rate")
    p.add_argument("--outdir", default=None)
    p.add_argument("--include-ratio-plot", action="store_true",
                   help="Also write the deprecated full/diagonal ratio plot.")
    args = p.parse_args()
    outdir = args.outdir or os.path.join(args.input, "figures", "estimator_diagnostics")
    os.makedirs(outdir, exist_ok=True)
    apply_style()

    scaling = _load_scaling(args.input)
    grid = _load_grid(args.input)
    print(f"Writing estimator-diagnostic figures to {outdir}")

    plot_lambda_vs_M(scaling, os.path.join(outdir, "diag1_lambda_full_sym_vs_Mmc"))
    plot_lambda_with_baselines(scaling, os.path.join(outdir, "diag2_lambda_with_baseline_overlay"))
    plot_separable_full_diag_exact(scaling, os.path.join(outdir, "diag3_separable_full_diag_exact_vs_Mmc"))
    plot_gaussian_true_zero(scaling, os.path.join(outdir, "diag4_gaussian_full_sym_true_zero"))
    plot_random_gap(scaling, "full_sym_minus_separable",
                    "full sym minus separable baseline",
                    "Random feature excess over separable baseline",
                    os.path.join(outdir, "diag5_random_minus_separable"))
    plot_random_gap(scaling, "full_sym_minus_gaussian",
                    "full sym minus Gaussian baseline",
                    "Random feature excess over Gaussian baseline",
                    os.path.join(outdir, "diag6_random_minus_gaussian"))
    plot_gamma_gaussian_true(scaling, os.path.join(outdir, "diag7_gamma_with_gaussian_true_line"))
    plot_gamma_compare(scaling, os.path.join(outdir, "diag8_gamma_family_comparison"))
    plot_self_adjoint_errors(scaling, os.path.join(outdir, "diag9_self_adjoint_errors_vs_Mmc"))
    plot_gap_vs_noise_scale(scaling, os.path.join(outdir, "diag10_gap_vs_sqrt_p_over_M"))
    plot_convergence_fit(scaling, os.path.join(outdir, "diag11_convergence_fit_Lambda_inf"))
    plot_gap_heatmap(scaling, os.path.join(outdir, "diag12_gap_heatmap"))
    plot_full_minus_diag_vs_N(grid, os.path.join(outdir, "diag13_full_minus_diag_vs_Ntheta"))
    if args.include_ratio_plot:
        plot_deprecated_ratio(grid, os.path.join(outdir, "deprecated_full_over_diag_vs_Ntheta"))
    print("Done.")


if __name__ == "__main__":
    main()
