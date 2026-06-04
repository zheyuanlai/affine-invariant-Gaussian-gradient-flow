"""Baseline corrections for natural-gradient local-rate result tables.

The raw ``Lambda_hat_full_sym`` estimator is intentionally left untouched by the
runner scripts. This module adds matched Gaussian/separable baselines after the
fact, so downstream summaries and plots can distinguish finite-sample spectral
inflation from possible excess behavior of coupled potentials.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

MATCH_KEYS = ["N_theta", "kappa_target", "M_mc"]
SEED_KEY = "seed"

BASELINE_COLUMNS = [
    "Lambda_hat_gaussian_baseline",
    "gamma_loc_gaussian_baseline",
    "full_sym_minus_gaussian",
    "gamma_minus_gaussian",
    "inverse_gamma_minus_gaussian_inverse",
    "Lambda_hat_separable_baseline",
    "gamma_loc_separable_baseline",
    "full_sym_minus_separable",
    "gamma_minus_separable",
    "full_sym_over_separable",
    "full_sym_minus_separable_exact",
    "full_sym_minus_diag",
    "diag_minus_exact",
    "abs_diag_minus_exact",
    "full_sym_excess_over_gaussian",
    "full_sym_excess_over_separable",
    "noise_baseline_available",
    "baseline_match",
    "gaussian_baseline_match",
    "separable_baseline_match",
    "baseline_correction_status",
]


def _nan():
    return float("nan")


def add_true_benchmark_columns(df):
    """Return a copy with analytic true columns for Gaussian rows.

    ``Lambda_true`` and ``gamma_true`` are only asserted for the analytic
    Gaussian target. Other families keep ``NaN`` here; separable exact quadrature
    remains a separate benchmark.
    """
    out = df.copy()
    fam = out.get("potential_family", out.get("family", pd.Series("", index=out.index))).astype(str)
    is_gaussian = fam == "gaussian"
    if "Lambda_true" not in out.columns:
        out["Lambda_true"] = np.nan
    if "gamma_true" not in out.columns:
        out["gamma_true"] = np.nan
    if "baseline_type" not in out.columns:
        out["baseline_type"] = ""
    out.loc[is_gaussian, "Lambda_true"] = 0.0
    out.loc[is_gaussian, "gamma_true"] = 1.0
    out.loc[is_gaussian, "baseline_type"] = "gaussian"
    sep_exact = out.get("Lambda_hat_separable_exact")
    if sep_exact is not None:
        mask = (~is_gaussian) & np.isfinite(pd.to_numeric(sep_exact, errors="coerce"))
        out.loc[mask, "baseline_type"] = out.loc[mask, "baseline_type"].replace("", "separable_exact")
    return out


def _row_value(row, col):
    if col not in row.index:
        return _nan()
    try:
        return float(row[col])
    except Exception:
        return _nan()


def _baseline_lookup(df, family, value_cols):
    """Build exact-seed and group-mean lookup dicts for one baseline family."""
    fam_df = df[df["potential_family"].astype(str) == family].copy()
    if fam_df.empty:
        return {}, {}
    for col in value_cols:
        if col not in fam_df.columns:
            fam_df[col] = np.nan
        fam_df[col] = pd.to_numeric(fam_df[col], errors="coerce")
    exact = {}
    for _, row in fam_df.iterrows():
        key = tuple(row[k] for k in MATCH_KEYS + [SEED_KEY])
        exact[key] = {col: _row_value(row, col) for col in value_cols}
    group = {}
    grouped = fam_df.groupby(MATCH_KEYS, dropna=False)
    for key, sub in grouped:
        if not isinstance(key, tuple):
            key = (key,)
        group[key] = {col: float(sub[col].mean()) for col in value_cols}
    return exact, group


def _get_baseline(row, exact, group, value_col):
    exact_key = tuple(row[k] for k in MATCH_KEYS + [SEED_KEY])
    if exact_key in exact:
        return exact[exact_key].get(value_col, _nan()), "seed_exact"
    group_key = tuple(row[k] for k in MATCH_KEYS)
    if group_key in group:
        return group[group_key].get(value_col, _nan()), "group_mean"
    return _nan(), "missing"


def _safe_diff(a, b):
    return float(a - b) if np.isfinite(a) and np.isfinite(b) else _nan()


def _safe_ratio(a, b):
    return float(a / b) if np.isfinite(a) and np.isfinite(b) and b != 0.0 else _nan()


def _inverse(x):
    return float(1.0 / x) if np.isfinite(x) and x != 0.0 else _nan()


def add_baseline_corrections(df):
    """Return ``df`` with Gaussian/separable baseline correction columns.

    Baselines are matched first by ``(N_theta, kappa_target, M_mc, seed)``. If no
    seed-specific row exists, the mean baseline over seeds at
    ``(N_theta, kappa_target, M_mc)`` is used and labeled as ``group_mean``. All
    missing baselines are explicit ``NaN`` with status text.
    """
    if df is None:
        return None
    out = add_true_benchmark_columns(df)
    required = ["potential_family"] + MATCH_KEYS + [SEED_KEY]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise ValueError(f"cannot baseline-correct results; missing columns {missing}")

    for col in ("Lambda_hat_full_sym", "gamma_loc", "inverse_gamma_loc",
                "Lambda_hat_diag", "Lambda_hat_separable_exact"):
        if col not in out.columns:
            out[col] = np.nan
        out[col] = pd.to_numeric(out[col], errors="coerce")

    g_exact, g_group = _baseline_lookup(
        out, "gaussian", ["Lambda_hat_full_sym", "gamma_loc", "inverse_gamma_loc"])
    s_exact, s_group = _baseline_lookup(
        out, "separable", ["Lambda_hat_full_sym", "gamma_loc", "inverse_gamma_loc"])

    records = []
    for _, row in out.iterrows():
        full = _row_value(row, "Lambda_hat_full_sym")
        gamma = _row_value(row, "gamma_loc")
        inv_gamma = _row_value(row, "inverse_gamma_loc")
        diag = _row_value(row, "Lambda_hat_diag")
        exact = _row_value(row, "Lambda_hat_separable_exact")

        g_lam, g_match = _get_baseline(row, g_exact, g_group, "Lambda_hat_full_sym")
        g_gamma, _ = _get_baseline(row, g_exact, g_group, "gamma_loc")
        g_inv, _ = _get_baseline(row, g_exact, g_group, "inverse_gamma_loc")
        s_lam, s_match = _get_baseline(row, s_exact, s_group, "Lambda_hat_full_sym")
        s_gamma, _ = _get_baseline(row, s_exact, s_group, "gamma_loc")
        s_inv, _ = _get_baseline(row, s_exact, s_group, "inverse_gamma_loc")

        statuses = []
        if g_match == "missing":
            statuses.append("missing_gaussian")
        if s_match == "missing":
            statuses.append("missing_separable")
        status = "ok" if not statuses else ("partial:" + ",".join(statuses)
                                             if len(statuses) == 1 else "missing:all_baselines")
        any_baseline = g_match != "missing" or s_match != "missing"
        if g_match == s_match:
            match = g_match
        elif any_baseline:
            match = "mixed"
        else:
            match = "missing"

        rec = {
            "Lambda_hat_gaussian_baseline": g_lam,
            "gamma_loc_gaussian_baseline": g_gamma,
            "full_sym_minus_gaussian": _safe_diff(full, g_lam),
            "gamma_minus_gaussian": _safe_diff(gamma, g_gamma),
            "inverse_gamma_minus_gaussian_inverse": _safe_diff(inv_gamma, g_inv),
            "Lambda_hat_separable_baseline": s_lam,
            "gamma_loc_separable_baseline": s_gamma,
            "full_sym_minus_separable": _safe_diff(full, s_lam),
            "gamma_minus_separable": _safe_diff(gamma, s_gamma),
            "full_sym_over_separable": _safe_ratio(full, s_lam),
            "full_sym_minus_separable_exact": _safe_diff(full, exact),
            "full_sym_minus_diag": _safe_diff(full, diag),
            "diag_minus_exact": _safe_diff(diag, exact),
            "abs_diag_minus_exact": abs(_safe_diff(diag, exact))
            if np.isfinite(_safe_diff(diag, exact)) else _nan(),
            "full_sym_excess_over_gaussian": _safe_diff(full, g_lam),
            "full_sym_excess_over_separable": _safe_diff(full, s_lam),
            "noise_baseline_available": bool(any_baseline),
            "baseline_match": match,
            "gaussian_baseline_match": g_match,
            "separable_baseline_match": s_match,
            "baseline_correction_status": status,
        }
        records.append(rec)

    baseline_df = pd.DataFrame.from_records(records, index=out.index)
    for col in BASELINE_COLUMNS:
        out[col] = baseline_df[col]
    return out


def aggregate_seed_summary(df):
    """Mean/std summary over seeds for raw and baseline-corrected quantities."""
    if df is None or df.empty:
        return pd.DataFrame()
    keys = ["potential_family", "N_theta", "kappa_target", "M_mc"]
    candidates = [
        "Lambda_hat_full_sym", "Lambda_hat_raw_forward", "Lambda_hat_diag",
        "Lambda_hat_separable_exact", "Lambda_true", "gamma_loc", "gamma_true",
        "inverse_gamma_loc", "full_sym_minus_exact", "full_sym_minus_diag",
        "full_sym_minus_gaussian", "full_sym_minus_separable", "full_sym_over_separable",
        "full_sym_minus_separable_exact", "diag_minus_exact", "abs_diag_minus_exact",
        "gamma_minus_gaussian", "gamma_minus_separable",
        "inverse_gamma_minus_gaussian_inverse", "self_adjoint_error_H_sym",
        "self_adjoint_error_L_star",
    ]
    agg_spec = {}
    for col in candidates:
        if col in df.columns:
            agg_spec[f"{col}_mean"] = (col, "mean")
            agg_spec[f"{col}_std"] = (col, "std")
    if "seed" in df.columns:
        agg_spec["n_seeds"] = ("seed", "count")
    if "baseline_correction_status" in df.columns:
        agg_spec["baseline_correction_status"] = ("baseline_correction_status", "first")
    summary = df.groupby(keys, dropna=False).agg(**agg_spec).reset_index()
    if "M_mc" in summary.columns:
        max_m = summary.groupby(["potential_family", "N_theta", "kappa_target"])["M_mc"].transform("max")
        summary["is_largest_M_row"] = summary["M_mc"] == max_m
    return summary
