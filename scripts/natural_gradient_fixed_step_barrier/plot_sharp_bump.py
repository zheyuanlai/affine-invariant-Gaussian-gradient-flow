"""Plot iterations-to-constant-factor-reduction against kappa on the sharp bump train.

Reads ``sharp_bump_summary.csv`` and writes ``figures/sharp_bump_scaling.{png,pdf}``:
one panel per train family, log-log, with ``kappa`` and ``kappa^2`` guides.

Usage::

    python scripts/natural_gradient_fixed_step_barrier/plot_sharp_bump.py \
        --indir outputs/natural_gradient_fixed_step_barrier
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.io_utils import ensure_dir
from src.common.plotting_style import apply_style, save_figure
import matplotlib.pyplot as plt

ARM_STYLE = {
    "theory": ("o-", "tab:red", r"$\Delta t=\gamma/\kappa$ (certified)"),
    "const": ("s-", "tab:blue", r"$\Delta t=\gamma$ (order one)"),
    "relcurv": ("^--", "tab:green", r"$\Delta t=\gamma/\max(cA,1)$ (adaptive)"),
}
FAMILY_TITLE = {
    "manuscript": "Appendix C train (built for $\\Delta t=\\gamma/\\kappa$)",
    "retuned": "train retuned to each arm's own $\\Delta t$",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", default="outputs/natural_gradient_fixed_step_barrier")
    ap.add_argument("--scheme", default="riemannian", choices=["riemannian", "kl"])
    args = ap.parse_args()

    df = pd.read_csv(os.path.join(args.indir, "sharp_bump_summary.csv"))
    df = df[(df["scheme"] == args.scheme) & (df["n_half"] > 0)]
    families = [f for f in ("manuscript", "retuned") if f in set(df["family"])]

    apply_style()
    fig, axes = plt.subplots(1, len(families), figsize=(5.2 * len(families), 4.0),
                             sharey=True)
    axes = np.atleast_1d(axes)
    for ax, family in zip(axes, families):
        sub = df[df["family"] == family]
        for arm, (fmt, color, label) in ARM_STYLE.items():
            g = sub[sub["arm"] == arm].sort_values("kappa")
            if g.empty:
                continue
            ax.loglog(g["kappa"], g["n_half"], fmt, color=color, label=label)
        k = np.array(sorted(sub["kappa"].unique()), dtype=float)
        anchor = float(sub[sub["arm"] == "theory"]["n_half"].min())
        ax.loglog(k, anchor * (k / k[0]), ":", color="0.5", lw=1.0)
        ax.loglog(k, anchor * (k / k[0]) ** 2, ":", color="0.3", lw=1.0)
        ax.text(k[-1], anchor * (k[-1] / k[0]), r"  $\kappa$", color="0.5", va="center")
        ax.text(k[-1], anchor * (k[-1] / k[0]) ** 2, r"  $\kappa^2$", color="0.3", va="center")
        ax.set_xlabel(r"condition number $\kappa$")
        ax.set_title(FAMILY_TITLE.get(family, family))
    axes[0].set_ylabel(r"iterations to $\Delta\mathcal{E}\leq\Delta\mathcal{E}_0/2$")
    axes[0].legend(loc="upper left")

    outdir = ensure_dir(os.path.join(args.indir, "figures"))
    path = save_figure(fig, os.path.join(outdir, "sharp_bump_scaling"))
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
