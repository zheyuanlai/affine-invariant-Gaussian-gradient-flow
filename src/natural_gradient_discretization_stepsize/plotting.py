"""Display-only plotting helpers for the discretization stepsize experiment.

Matplotlib only (no seaborn). None of these functions mutate or re-save the
experiment data; they only control how committed CSV values are displayed.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Representative stepsizes overlaid on the energy-gap convergence figures.
# Kept as a sane default for globally smooth targets.
REPRESENTATIVE_DT = [0.1, 0.5, 1.0]

# Per-target representative stepsizes for the energy-gap figures. These are
# chosen to be *convergent* for every lambda of that target, so the energy-gap
# curves show genuine convergence rather than blow-up. The quartic
# (literature) target has unbounded curvature, so its stable stepsizes are far
# smaller than the smooth/Gaussian targets'.
REPRESENTATIVE_DT_BY_TARGET = {
    "gaussian_posterior": [0.1, 0.5, 1.0],
    "smooth_logconcave": [0.1, 0.5, 1.0],
    "literature_logconcave": [0.005, 0.01, 0.02],
}

# Per-target stepsizes used by the appendix convergence-speed figures, where each
# stepsize gets its own subplot column (rows = lambda) so a panel shows only the
# two methods. Chosen so that BOTH schemes are stable at every lambda, spanning
# fine -> coarse over the convergent range of that target.
CONVERGENCE_DT_BY_TARGET = {
    "gaussian_posterior": [0.05, 0.1, 0.5, 1.0],
    "smooth_logconcave": [0.05, 0.1, 0.5, 1.0],
    "literature_logconcave": [0.002, 0.005, 0.01, 0.02],
}

# Display floor for the energy-gap log axes (data is never modified).
GAP_FLOOR = 1e-12

METHOD_STYLE = {
    "riemannian": {"color": "#1f77b4", "label": "Riemannian", "ls": "-"},
    "kl": {"color": "#d62728", "label": "KL", "ls": "--"},
}

LAMBDA_ORDER = [0.01, 0.1, 1.0]

# Three-level stability classification used by the heatmaps.
#   2 = stable and monotone
#   1 = stable but non-monotone
#   0 = unstable / diverged (not stable)
CLASS_COLORS = ["#b2182b", "#fdae61", "#1a9850"]   # red / amber / green
CLASS_LABELS = ["unstable / diverged", "stable, non-monotone", "stable & monotone"]


def classify_level(stable, monotone):
    """Map ``(stable, monotone)`` flags to the 3-level heatmap code 0/1/2."""
    if not bool(stable):
        return 0
    return 2 if bool(monotone) else 1


def clip_for_log(y, floor=GAP_FLOOR):
    """Clip from below to ``floor`` and map non-finite to NaN (display only)."""
    y = np.asarray(y, dtype=float)
    return np.where(np.isfinite(y), np.maximum(y, floor), np.nan)


def semilogy_clipped(ax, x, y, floor=GAP_FLOOR, **kwargs):
    yy = clip_for_log(y, floor)
    x = np.asarray(x, dtype=float)
    if np.any(np.isfinite(yy)):
        ax.semilogy(x, yy, **kwargs)


def class_cmap_norm():
    """A discrete colormap + norm for the 0/1/2 stability classification."""
    cmap = mcolors.ListedColormap(CLASS_COLORS)
    norm = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    return cmap, norm


def stability_heatmap(ax, mat, dts, lambdas, *, title=None, dt_theory=None):
    """Render a stability classification heatmap.

    ``mat`` has shape ``(len(lambdas), len(dts))`` with entries in ``{0, 1, 2}``.
    ``dt_theory`` (optional) draws a vertical marker at the theoretical stepsize.
    """
    cmap, norm = class_cmap_norm()
    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm, origin="upper")
    ax.set_xticks(range(len(dts)))
    ax.set_xticklabels([f"{d:g}" for d in dts], rotation=45, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(lambdas)))
    ax.set_yticklabels([f"{l:g}" for l in lambdas])
    ax.set_xlabel(r"$\Delta t$")
    ax.set_ylabel(r"$\lambda$")
    if title:
        ax.set_title(title)
    if dt_theory is not None and np.isfinite(dt_theory):
        # Place the marker between the dt columns straddling dt_theory.
        xs = np.array(dts, dtype=float)
        pos = np.searchsorted(xs, dt_theory) - 0.5
        ax.axvline(pos, color="black", lw=1.4, ls=":")
    return im


def add_class_legend(fig, loc="lower center", ncol=3):
    """Attach a shared legend describing the three stability levels."""
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", markersize=9,
                          markerfacecolor=c, markeredgecolor="0.3", label=lab)
               for c, lab in zip(CLASS_COLORS, CLASS_LABELS)]
    fig.legend(handles=handles, loc=loc, ncol=ncol, frameon=False,
               bbox_to_anchor=(0.5, -0.02))


def dt_shades(dts, cmap_name="viridis", lo=0.12, hi=0.88):
    """Map a sorted list of stepsizes to evenly spaced colormap shades.

    Smaller ``dt`` -> lighter, larger ``dt`` -> darker, so the legend reads as
    "coarser stepsize = darker". Returns a dict ``dt -> rgba``.
    """
    cmap = plt.get_cmap(cmap_name)
    n = max(1, len(dts) - 1)
    return {dt: cmap(lo + (hi - lo) * k / n) for k, dt in enumerate(sorted(dts))}
