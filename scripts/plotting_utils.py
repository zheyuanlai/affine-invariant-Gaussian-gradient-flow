"""
Shared display-only helpers for manuscript-ready figures.

This module is imported by both ``src/plotting.py`` (Gaussian target) and
``src/lc_plotting.py`` (log-concave target).  Every function here changes only
how values are *displayed* — none of them mutate, recompute, or re-save the
underlying experiment data.  The CSV files are never touched.

Contents
--------
Display-floor constants:
    PLOT_FLOOR_OBJECTIVE_GAUSSIAN, PLOT_FLOOR_RESIDUAL_GAUSSIAN,
    PLOT_FLOOR_OBJECTIVE_LOGCONCAVE, PLOT_FLOOR_RESIDUAL_LOGCONCAVE,
    PLOT_FLOOR_MEAN_LOGCONCAVE, CHI_COV_RESIDUAL_MIN, NOT_REACHED_LABEL_MODE

Schema helpers:
    get_metric_column(df, candidates)

Curve helpers:
    clip_for_log_plot(y, floor)
    semilogy_clipped(ax, x, y, floor, **kwargs)
    series_is_identically_small(y, floor)

Chi helpers:
    mask_chi_when_residual_small(chi, cov_residual, threshold)

Heatmap helpers:
    is_reached(value, T)
    format_not_reached_label(T)
    plot_heatmap_with_not_reached_mask(ax, mat, ...)
    speedup_diverging_norm(vmin, vcenter, vmax)

Style:
    apply_manuscript_style()
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# ---------------------------------------------------------------------------
# Display-floor constants
# ---------------------------------------------------------------------------
# Curves below these floors are pinned to the floor for log-scale display so
# the main figures do not dive into machine-precision (1e-15..1e-16) tails.
# The *data* is never modified — only the plotted y-values.

PLOT_FLOOR_OBJECTIVE_GAUSSIAN   = 1e-12   # Gaussian normalised KL energy
PLOT_FLOOR_RESIDUAL_GAUSSIAN    = 1e-14   # Gaussian mean/cov/shape/etc. errors

PLOT_FLOOR_OBJECTIVE_LOGCONCAVE = 1e-5    # log-concave normalised objective gap
PLOT_FLOOR_RESIDUAL_LOGCONCAVE  = 1e-10   # log-concave residual diagnostics
PLOT_FLOOR_MEAN_LOGCONCAVE      = 1e-5    # log-concave whitened mean error

# chi = (trace residual)^2 / (n * ||residual||_F^2) is a small/small ratio that
# becomes numerically meaningless once the covariance residual is tiny.
CHI_COV_RESIDUAL_MIN            = 1e-7

# How to label time-to-tolerance cells that never reached tolerance.
NOT_REACHED_LABEL_MODE          = "greater_than_T"


# ---------------------------------------------------------------------------
# Schema-compatibility helper
# ---------------------------------------------------------------------------

def get_metric_column(df, candidates):
    """Return the first column name in ``candidates`` that exists in ``df``.

    Lets the same plotting code consume both the Gaussian schema (e.g.
    ``norm_energy``, ``mean_error``, ``cosine_error``) and the log-concave
    schema (``normalized_objective_gap``, ``whitened_mean_error``,
    ``cosine_error_to_star``).

    Returns None if none of the candidates are present.
    """
    for col in candidates:
        if col in df.columns:
            return col
    return None


# ---------------------------------------------------------------------------
# Curve display helpers
# ---------------------------------------------------------------------------

def clip_for_log_plot(y, floor):
    """Clip values from below to ``floor`` for log-scale *display* only.

    Non-finite entries (inf / NaN) are mapped to NaN so matplotlib simply skips
    them (breaking the line) instead of drawing spurious points.

    Returns a float ndarray of the same shape as ``y``.
    """
    y = np.asarray(y, dtype=float)
    return np.where(np.isfinite(y), np.maximum(y, floor), np.nan)


def semilogy_clipped(ax, x, y, floor, **kwargs):
    """semilogy with values clipped from below to ``floor`` (display only)."""
    yy = clip_for_log_plot(y, floor)
    x = np.asarray(x, dtype=float)
    if np.any(np.isfinite(yy)):
        ax.semilogy(x, yy, **kwargs)


def series_is_identically_small(y, floor):
    """True if every finite value in ``y`` is at or below ``floor``.

    Used to detect panels that are zero by construction (e.g. covariance error
    for a mean-only initialization) so they can be annotated rather than shown
    as an empty log-scale axis.
    """
    y = np.asarray(y, dtype=float)
    finite = y[np.isfinite(y)]
    if finite.size == 0:
        return True
    return float(np.max(np.abs(finite))) <= floor


# ---------------------------------------------------------------------------
# chi masking
# ---------------------------------------------------------------------------

def mask_chi_when_residual_small(chi, cov_residual, threshold=CHI_COV_RESIDUAL_MIN):
    """Return ``chi`` with entries blanked (NaN) where the residual is tiny.

    chi is a ratio of two small numbers once convergence is reached; below the
    residual threshold it is dominated by numerical noise and should not be
    plotted.  Only the displayed values change; the data is untouched.
    """
    chi = np.asarray(chi, dtype=float).copy()
    res = np.asarray(cov_residual, dtype=float)
    chi[~(res > threshold)] = np.nan
    return chi


# ---------------------------------------------------------------------------
# Heatmap helpers
# ---------------------------------------------------------------------------

def is_reached(value, T):
    """True iff a time-to-tolerance value is finite and within the horizon T."""
    return bool(np.isfinite(value) and value <= T + 1e-9)


def format_not_reached_label(T):
    """Label for cells that never reached tolerance within the horizon T."""
    if NOT_REACHED_LABEL_MODE == "greater_than_T":
        return f">{T:g}"
    return "n/r"


def plot_heatmap_with_not_reached_mask(
    ax, mat, *, T, xticklabels, yticklabels,
    cmap="viridis_r", cbar_label="Time", masked_color="0.85",
    value_fmt="{:.1f}", fig=None, annotate=True, vmin=None, vmax=None,
):
    """Render a heatmap where not-reached / missing cells (NaN) are masked.

    ``mat`` must already contain ``np.nan`` for not-reached or missing cells
    and finite reached times elsewhere.  When ``vmin``/``vmax`` are not given,
    the colour scale is derived ONLY from the finite (reached) cells, so an
    unreachable run never inflates the scale.  Pass explicit ``vmin``/``vmax``
    to share one scale across several panels.  Masked cells are drawn in a
    light neutral colour and annotated ``>T``.

    Returns the AxesImage handle.
    """
    mat = np.asarray(mat, dtype=float)
    finite = mat[np.isfinite(mat)]
    if vmin is None:
        vmin = float(finite.min()) if finite.size else 0.0
    if vmax is None:
        vmax = float(finite.max()) if finite.size else 1.0
    if vmin == vmax:
        vmax = vmin + 1.0

    masked = np.ma.masked_invalid(mat)
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(color=masked_color)

    im = ax.imshow(masked, aspect="auto", cmap=cmap_obj, origin="upper",
                   vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(xticklabels)))
    ax.set_xticklabels(xticklabels)
    ax.set_yticks(range(len(yticklabels)))
    ax.set_yticklabels(yticklabels)
    if fig is not None:
        fig.colorbar(im, ax=ax, label=cbar_label)

    if annotate:
        nr_label = format_not_reached_label(T)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if np.isfinite(v):
                    frac = (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5
                    color = "white" if frac > 0.55 else "black"
                    ax.text(j, i, value_fmt.format(v), ha="center", va="center",
                            fontsize=8, color=color)
                else:
                    ax.text(j, i, nr_label, ha="center", va="center",
                            fontsize=8, color="0.30")
    return im


def speedup_diverging_norm(vmin=0.5, vcenter=1.0, vmax=1.5):
    """Diverging colour norm centred at ratio = 1 (no speedup).

    ratio < 1  -> faster than tau = 0
    ratio = 1  -> no change   (the visual centre)
    ratio > 1  -> slower than tau = 0
    """
    return mcolors.TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)


# ---------------------------------------------------------------------------
# Manuscript style
# ---------------------------------------------------------------------------

def apply_manuscript_style():
    """Apply readable, consistent rcParams for manuscript figures.

    Display-only: affects font sizes and line defaults, never the data.
    Safe to call repeatedly.
    """
    plt.rcParams.update({
        "figure.dpi":        110,
        "savefig.dpi":       200,
        "font.size":         11,
        "axes.titlesize":    11,
        "axes.labelsize":    11,
        "xtick.labelsize":   9,
        "ytick.labelsize":   9,
        "legend.fontsize":   8,
        "legend.framealpha": 0.9,
        "lines.linewidth":   1.8,
        "axes.grid":         False,
    })
