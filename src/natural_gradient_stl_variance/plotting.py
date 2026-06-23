"""Display helpers and figure builders for the STL variance experiment.

Matplotlib only. None of these mutate or re-save experiment data; they read the
committed CSVs and return :class:`matplotlib.figure.Figure` objects. The same
builders are used by the per-group ``plot_results.py`` (writing into
``<outdir>/figures``) and by ``reports/make_report_assets.py`` (writing the
report-named figures into ``reports/assets/figs``).
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from src.natural_gradient_stl_variance.states import STATE_NAMES

# Algorithm method display: paired by scheme (color), baseline vs STL (style).
METHOD_STYLE = {
    "riemannian":     {"color": "#1f77b4", "ls": "-",  "label": "Riemannian (base)"},
    "riemannian_stl": {"color": "#1f77b4", "ls": "--", "label": "Riemannian (STL)"},
    "kl":             {"color": "#d62728", "ls": "-",  "label": "KL (base)"},
    "kl_stl":         {"color": "#d62728", "ls": "--", "label": "KL (STL)"},
}
METHOD_ORDER = ["riemannian", "riemannian_stl", "kl", "kl_stl"]
SCHEME_PAIR = {"riemannian": ("riemannian", "riemannian_stl"),
               "kl": ("kl", "kl_stl")}

KIND_TITLE = {"gaussian": "well-specified Gaussian",
              "log_cosh": "misspecified log-cosh"}

GAP_FLOOR = 1e-16


def _clip(y, floor=GAP_FLOOR):
    return np.clip(np.asarray(y, dtype=np.float64), floor, None)


def _config_label(d, kappa, tau, kind):
    if kind == "log_cosh":
        return rf"$d{{=}}{int(d)},\kappa{{=}}{kappa:g},\tau{{=}}{tau:g}$"
    return rf"$d{{=}}{int(d)},\kappa{{=}}{kappa:g}$"


# ---------------------------------------------------------------------------
# Experiment A: estimator-level variance ratio
# ---------------------------------------------------------------------------

def fig_estimator_variance_ratio(est_df, kind):
    """Median Fisher--Rao variance ratio ``Var_FR(stl)/Var_FR(base)`` by state.

    One panel per dimension ``d``; one line per ``(kappa[, tau])`` config. The
    ratio is the median over seeds (clipped to a tiny positive floor for the log
    axis); the dashed line at ratio ``= 1`` separates variance reduction (below)
    from variance inflation (above).
    """
    sub = est_df[est_df.kind == kind].copy()
    ds = sorted(sub.d.unique())
    fig, axes = plt.subplots(1, len(ds), figsize=(4.4 * len(ds), 3.6),
                             squeeze=False, sharey=True)
    states = [s for s in STATE_NAMES if s in set(sub.state.unique())]
    xpos = np.arange(len(states))
    for j, d in enumerate(ds):
        ax = axes[0, j]
        cell = sub[sub.d == d]
        configs = (cell[["kappa", "tau"]].drop_duplicates()
                   .sort_values(["kappa", "tau"]).itertuples(index=False))
        cmap = plt.cm.viridis(np.linspace(0.1, 0.85, max(1, len(set(
            map(tuple, cell[["kappa", "tau"]].values))))))
        for k, (kappa, tau) in enumerate(configs):
            ys = []
            for st in states:
                r = cell[(cell.state == st) & np.isclose(cell.kappa, kappa)
                         & np.isclose(cell.tau, tau)]
                ys.append(np.median(r.ratio_fr.values) if len(r) else np.nan)
            ax.plot(xpos, _clip(ys), "o-", color=cmap[k % len(cmap)],
                    label=_config_label(d, kappa, tau, kind), lw=1.6, ms=5)
        ax.axhline(1.0, color="0.4", ls="--", lw=1.0)
        ax.set_yscale("log")
        ax.set_xticks(xpos)
        ax.set_xticklabels(states, rotation=35, ha="right", fontsize=8)
        ax.set_title(rf"$d={int(d)}$")
        ax.grid(True, which="both", alpha=0.25)
        if j == 0:
            ax.set_ylabel(r"$\mathrm{Var}_{FR}(b_{\mathrm{stl}})/"
                          r"\mathrm{Var}_{FR}(b_{\mathrm{base}})$")
    axes[0, -1].legend(fontsize=7.5, loc="best")
    fig.suptitle(f"Estimator-level FR variance ratio -- {KIND_TITLE[kind]} target",
                 y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


def fig_variance_by_distance(est_df, kind):
    """FR variance ratio vs mean-distance to the optimum (mean-perturbation states).

    Only the ``optimum/near/medium/far`` states (which perturb the mean) are
    shown, against the offset multiplier ``rho``; one line per ``(d, kappa[, tau])``.
    """
    sub = est_df[(est_df.kind == kind)
                 & est_df.state.isin(["optimum", "near", "medium", "far"])].copy()
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    configs = (sub[["d", "kappa", "tau"]].drop_duplicates()
               .sort_values(["d", "kappa", "tau"]).itertuples(index=False))
    configs = list(configs)
    cmap = plt.cm.plasma(np.linspace(0.05, 0.85, max(1, len(configs))))
    for k, (d, kappa, tau) in enumerate(configs):
        cell = sub[(sub.d == d) & np.isclose(sub.kappa, kappa)
                   & np.isclose(sub.tau, tau)]
        agg = (cell.groupby("distance_to_optimum").ratio_fr.median()
               .sort_index())
        if agg.empty:
            continue
        ax.plot(agg.index.values + 1e-3, _clip(agg.values), "o-",
                color=cmap[k % len(cmap)], lw=1.5, ms=4,
                label=_config_label(d, kappa, tau, kind))
    ax.axhline(1.0, color="0.4", ls="--", lw=1.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"mean offset $\rho$ (in marginal std; $\rho{=}0$ is the optimum)")
    ax.set_ylabel(r"FR variance ratio (stl/base)")
    ax.set_title(f"STL variance ratio vs distance -- {KIND_TITLE[kind]}")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7.5, loc="best")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Experiment B: algorithm-level trajectories
# ---------------------------------------------------------------------------

def _representative_config(df, kind):
    """Pick a representative ``(d, kappa, tau, batch_size)`` cell for a kind."""
    sub = df[df.kind == kind]
    if sub.empty:
        return None
    d = sorted(sub.d.unique())[len(sub.d.unique()) // 2]
    sub = sub[sub.d == d]
    kappa = sorted(sub.kappa.unique())[-1]            # hardest conditioning
    sub = sub[np.isclose(sub.kappa, kappa)]
    tau = sorted(sub.tau.unique())[-1]
    sub = sub[np.isclose(sub.tau, tau)]
    bs = sorted(sub.batch_size.unique())[0]           # smallest batch (noisiest)
    return int(d), float(kappa), float(tau), int(bs)


def fig_algorithm_gap(long_df, kind):
    """Energy gap vs iteration (median over seeds) for the four methods.

    A representative cell (mid dimension, hardest conditioning, smallest batch)
    is shown; the STL curves (dashed) sit below their baseline (solid) once the
    iterate is near the optimum -- a lower stochastic noise floor.
    """
    rep = _representative_config(long_df, kind)
    if rep is None:
        return None
    d, kappa, tau, bs = rep
    cell = long_df[(long_df.kind == kind) & (long_df.d == d)
                   & np.isclose(long_df.kappa, kappa)
                   & np.isclose(long_df.tau, tau)
                   & (long_df.batch_size == bs)]
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    for method in METHOD_ORDER:
        s = cell[cell.method == method]
        if s.empty:
            continue
        agg = s.groupby("step").energy_gap.median().sort_index()
        st = METHOD_STYLE[method]
        ax.plot(agg.index.values, _clip(agg.values), color=st["color"],
                ls=st["ls"], lw=1.7, label=st["label"])
    ax.set_yscale("log")
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"energy gap $\mathcal{E}(a_n)-\mathcal{E}(a_\star)$")
    ax.set_title(rf"{KIND_TITLE[kind]}: gap vs iteration "
                 rf"({_config_label(d, kappa, tau, kind)}, $B{{=}}{bs}$)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    return fig


def fig_noise_floor(tail_df):
    """Tail (noise-floor) energy gap, baseline vs STL, per scheme and target.

    For every ``(kind, scheme)`` we scatter the per-cell median tail gap of the
    baseline against the STL method across all configs and batch sizes; points
    below the diagonal are cells where STL lowers the noise floor.
    """
    kinds = sorted(tail_df.kind.unique())
    fig, axes = plt.subplots(1, len(kinds), figsize=(4.8 * len(kinds), 4.2),
                             squeeze=False)
    keys = ["d", "kappa", "tau", "batch_size"]
    for j, kind in enumerate(kinds):
        ax = axes[0, j]
        sub = tail_df[tail_df.kind == kind]
        med = (sub.groupby(keys + ["method"]).tail_median_gap.median()
               .reset_index())
        for scheme, (base, stl) in SCHEME_PAIR.items():
            b = med[med.method == base].set_index(keys).tail_median_gap
            s = med[med.method == stl].set_index(keys).tail_median_gap
            common = b.index.intersection(s.index)
            if len(common) == 0:
                continue
            xb = _clip(b.loc[common].values)
            ys = _clip(s.loc[common].values)
            ax.scatter(xb, ys, s=42, alpha=0.8,
                       color=METHOD_STYLE[stl]["color"],
                       label=f"{scheme}")
        lim_lo = GAP_FLOOR
        all_vals = _clip(med.tail_median_gap.values)
        lim_hi = float(np.nanmax(all_vals)) * 3 if all_vals.size else 1.0
        ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], color="0.4", ls="--", lw=1.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("baseline tail gap")
        ax.set_ylabel("STL tail gap")
        ax.set_title(KIND_TITLE[kind])
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=8, loc="best", title="scheme")
    fig.suptitle("Stochastic noise floor: STL vs baseline (below diagonal = STL lower)",
                 y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def fig_batchsize_effect(tail_df, kind=None):
    """Tail noise-floor gap vs mini-batch size, per method.

    One panel per kind; lines are the median (over configs and seeds) tail gap of
    each method against the batch size. Both baseline and STL floors fall as the
    batch grows; STL stays below its paired baseline.
    """
    kinds = sorted(tail_df.kind.unique()) if kind is None else [kind]
    fig, axes = plt.subplots(1, len(kinds), figsize=(4.8 * len(kinds), 3.8),
                             squeeze=False, sharey=True)
    for j, kd in enumerate(kinds):
        ax = axes[0, j]
        sub = tail_df[tail_df.kind == kd]
        for method in METHOD_ORDER:
            s = sub[sub.method == method]
            if s.empty:
                continue
            agg = s.groupby("batch_size").tail_median_gap.median().sort_index()
            st = METHOD_STYLE[method]
            ax.plot(agg.index.values, _clip(agg.values), marker="o",
                    color=st["color"], ls=st["ls"], lw=1.6, ms=5, label=st["label"])
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xlabel("mini-batch size $B$")
        ax.set_title(KIND_TITLE[kd])
        ax.grid(True, which="both", alpha=0.25)
        if j == 0:
            ax.set_ylabel("median tail energy gap")
    axes[0, -1].legend(fontsize=8, loc="best")
    fig.suptitle("Noise floor vs mini-batch size", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig
