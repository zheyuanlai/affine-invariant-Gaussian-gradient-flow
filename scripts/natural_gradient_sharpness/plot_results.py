#!/usr/bin/env python3
"""Plot the sharp-rate experiment outputs."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _save(fig, outdir, stem):
    fig.tight_layout()
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_global(outdir):
    spiral = pd.read_csv(outdir / "global_continuous_summary.csv")
    bump = pd.read_csv(outdir / "global_discrete_summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    ax = axes[0]
    ax.loglog(spiral.kappa_star_lower, spiral.time_to_radius_threshold, "o-", label="spiral flow")
    ref = spiral.time_to_radius_threshold.iloc[0] * spiral.kappa_star_lower / spiral.kappa_star_lower.iloc[0]
    ax.loglog(spiral.kappa_star_lower, ref, "--", color="0.35", label=r"$\kappa_\star$")
    ax.set(xlabel=r"$\kappa_\star$", ylabel="time to 2% mean contraction", title="Continuous global localization")
    ax.legend(frameon=False)
    ax = axes[1]
    for method, marker in (("riemannian", "o"), ("kl", "s")):
        frame = bump[bump.method == method]
        ax.loglog(frame.kappa, frame.steps_to_gap_threshold, marker + "-", label=method)
    ref = bump.steps_to_gap_threshold.dropna().iloc[0] * bump.kappa.unique() ** 2 / bump.kappa.unique()[0] ** 2
    ax.loglog(bump.kappa.unique(), ref, "--", color="0.35", label=r"$\kappa^2$")
    ax.set(xlabel=r"$\kappa$", ylabel="iterations to fixed gap reduction", title="Both fixed-step schemes")
    ax.legend(frameon=False)
    _save(fig, outdir, "sharp_global_rates")


def plot_local(outdir):
    local = pd.read_csv(outdir / "local_spectral_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.9))
    regimes = [
        ("fixed_dim_log", r"$\log(e\kappa)$", "1-D bump"),
        ("dimension_branch", r"$\sqrt{n\log(e\kappa)}$", "gamma shell"),
        ("uniform_saturation", r"$\sqrt{\kappa}$", "uniform saturation"),
    ]
    for ax, (regime, xlabel, title) in zip(axes, regimes):
        frame = local[local.regime == regime]
        for construction, group in frame.groupby("construction"):
            ax.loglog(group.predicted_inverse_gap_scale, group.inverse_gap, "o-", label=construction)
        lo = min(frame.predicted_inverse_gap_scale.min(), frame.inverse_gap.min())
        hi = max(frame.predicted_inverse_gap_scale.max(), frame.inverse_gap.max())
        ratio = np.median(frame.inverse_gap / frame.predicted_inverse_gap_scale)
        ax.loglog([lo, hi], ratio * np.array([lo, hi]), "--", color="0.35", label="slope 1")
        ax.set(xlabel=xlabel, ylabel=r"measured $\gamma^{-1}$", title=title)
        ax.legend(frameon=False, fontsize=8)
    _save(fig, outdir, "sharp_local_phase_diagram")

    dyn = pd.read_csv(outdir / "local_dynamics_summary.csv")
    disc = dyn[dyn.method != "continuous"]
    fig, ax = plt.subplots(figsize=(5.4, 4.3))
    for method, marker in (("riemannian", "o"), ("kl", "s")):
        frame = disc[disc.method == method]
        ax.loglog(1.0 / (frame.dt * frame.spectral_gamma), frame.effective_iterations,
                  marker, label=method)
    lim = ax.get_xlim()
    ax.loglog(lim, lim, "--", color="0.35", label="predicted slope 1")
    ax.set(xlabel=r"$(\Delta t\,\gamma)^{-1}$", ylabel="measured e-fold iterations",
           title="Local discrete rates")
    ax.legend(frameon=False)
    _save(fig, outdir, "sharp_local_discrete_rates")


def plot_global_to_local(outdir):
    traj = pd.read_csv(outdir / "global_to_local_trajectories.csv")
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8))
    c = traj[(traj.method == "continuous") & (traj.variant == "continuous")]
    axes[0].semilogy(c.elapsed_time, c.error)
    axes[0].set(xlabel="flow time", ylabel="Fisher--Rao local norm proxy", title="Continuous flow")
    for ax, method in zip(axes[1:], ("riemannian", "kl")):
        frame = traj[traj.method == method]
        global_part = frame[frame.variant == "shared_global"]
        ax.semilogy(global_part.progress, global_part.error, color="0.2", label="global step")
        for variant, style in (("fixed", "--"), ("two_stage", "-")):
            tail = frame[frame.variant == variant]
            ax.semilogy(tail.progress, tail.error, style, label=variant)
        ax.set(xlabel="iteration", title=method.capitalize())
        ax.legend(frameon=False, fontsize=8)
    _save(fig, outdir, "sharp_global_to_local")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="outputs/natural_gradient_sharpness")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    plot_global(outdir)
    plot_local(outdir)
    plot_global_to_local(outdir)
    print(f"wrote figures to {outdir}")


if __name__ == "__main__":
    main()
