"""Shared matplotlib plotting defaults (matplotlib only, no seaborn).

Figures are saved as both PNG and PDF with concise labels and no oversized
titles. Import :func:`apply_style` once at the top of a plotting script and use
:func:`save_figure` to write paired PNG/PDF output.
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")  # headless backend; safe for scripts and CI
import matplotlib.pyplot as plt  # noqa: E402


RC_DEFAULTS = {
    "figure.dpi": 120,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "legend.frameon": False,
    "lines.linewidth": 1.6,
    "lines.markersize": 5,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
}


def apply_style():
    """Apply the shared rcParams in place."""
    plt.rcParams.update(RC_DEFAULTS)


def save_figure(fig, path_no_ext, formats=("png", "pdf")):
    """Save ``fig`` to ``<path_no_ext>.<ext>`` for each ext; return the paths."""
    out_dir = os.path.dirname(path_no_ext)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    paths = []
    for ext in formats:
        p = f"{path_no_ext}.{ext}"
        fig.savefig(p)
        paths.append(p)
    return paths
