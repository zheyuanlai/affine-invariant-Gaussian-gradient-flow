"""Sample-size convergence diagnostics for local-rate estimator tables."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd


def sym_dim(N_theta):
    """Dimension of the symmetric-matrix operator domain."""
    N = int(N_theta)
    return N * (N + 1) // 2


def add_noise_scale_columns(df):
    """Add row-level ``p_sym_dim``, ``M_over_p`` and ``sqrt_p_over_M`` columns."""
    out = df.copy()
    if "N_theta" not in out.columns or "M_mc" not in out.columns:
        return out
    p = out["N_theta"].astype(int).map(sym_dim).astype(float)
    M = pd.to_numeric(out["M_mc"], errors="coerce").astype(float)
    out["p_sym_dim"] = p.astype(int)
    out["M_over_p"] = M / p
    out["sqrt_p_over_M"] = np.sqrt(p / M)
    return out


def _linear_fit(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x[mask], dtype=float)
    y = np.asarray(y[mask], dtype=float)
    if x.size < 2 or np.allclose(x, x[0]):
        return float("nan"), float("nan"), float("nan"), "insufficient_data"
    A = np.column_stack([np.ones_like(x), x])
    intercept, slope = np.linalg.lstsq(A, y, rcond=None)[0]
    pred = intercept + slope * x
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    return float(intercept), float(slope), float(r2), "ok"


def _last_finite(sub, col):
    if col not in sub.columns:
        return float("nan")
    vals = pd.to_numeric(sub[col], errors="coerce")
    vals = vals[np.isfinite(vals)]
    return float(vals.iloc[-1]) if len(vals) else float("nan")


def _first_finite(sub, col):
    if col not in sub.columns:
        return float("nan")
    vals = pd.to_numeric(sub[col], errors="coerce")
    vals = vals[np.isfinite(vals)]
    return float(vals.iloc[0]) if len(vals) else float("nan")


def _close(a, b, tolerance):
    return np.isfinite(a) and np.isfinite(b) and abs(a - b) <= tolerance


def _tail_stable(values, tolerance):
    vals = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    if vals.size < 3:
        return False
    return abs(vals[-1] - vals[-2]) <= max(tolerance, 0.1 * max(1.0, abs(vals[-1])))


def _group_warning(family, status, lam_inf, largest_lam, exact, tail_values,
                   tolerance):
    if status != "ok":
        return False, status
    fam = str(family)
    if fam == "gaussian":
        ok = _close(lam_inf, 0.0, tolerance) or abs(largest_lam) <= tolerance
        return ok, "" if ok else "gaussian_full_not_close_to_zero"
    if fam == "separable":
        if not np.isfinite(exact):
            return False, "missing_separable_exact"
        ok = _close(lam_inf, exact, tolerance) or _close(largest_lam, exact, tolerance)
        return ok, "" if ok else "separable_full_not_close_to_exact"
    stable = _tail_stable(tail_values, tolerance)
    return stable, "" if stable else "coupled_tail_not_stable"


def fit_scaling_diagnostics(df, *, y_col="Lambda_hat_full_sym", tolerance=0.1):
    """Fit ``Lambda(M) = Lambda_inf + c * sqrt(p/M)`` per family/N/kappa/seed.

    The same fit also represents ``Lambda_inf + c / sqrt(M)`` up to a constant
    slope rescaling within each fixed-``N`` group. The returned table is
    deliberately conservative: convergence flags are false unless the largest-M
    or fitted intercept is close to the relevant Gaussian/separable benchmark,
    or a coupled-family tail is visibly stable.
    """
    if df is None or df.empty:
        return pd.DataFrame()
    work = add_noise_scale_columns(df)
    keys = ["potential_family", "N_theta", "kappa_target", "seed"]
    missing = [c for c in keys + ["M_mc", y_col] if c not in work.columns]
    if missing:
        raise ValueError(f"cannot fit scaling diagnostics; missing columns {missing}")

    rows = []
    for key, sub in work.groupby(keys, dropna=False):
        fam, N, kappa, seed = key
        sub = sub.sort_values("M_mc")
        y = pd.to_numeric(sub[y_col], errors="coerce").to_numpy(dtype=float)
        x = pd.to_numeric(sub["sqrt_p_over_M"], errors="coerce").to_numpy(dtype=float)
        lam_inf, slope, r2, status = _linear_fit(x, y)
        largest_M = int(pd.to_numeric(sub["M_mc"], errors="coerce").max())
        p = sym_dim(N)
        largest_lam = _last_finite(sub, y_col)
        exact = _last_finite(sub, "Lambda_hat_separable_exact")
        g_base = _last_finite(sub, "Lambda_hat_gaussian_baseline")
        s_base = _last_finite(sub, "Lambda_hat_separable_baseline")
        conv, warning = _group_warning(fam, status, lam_inf, largest_lam, exact, y, tolerance)

        rows.append({
            "potential_family": fam,
            "N_theta": int(N),
            "kappa_target": float(kappa),
            "seed": int(seed),
            "p_sym_dim": int(p),
            "M_over_p": float(largest_M / p),
            "sqrt_p_over_M": float(math.sqrt(p / largest_M)),
            "Lambda_inf_fit": lam_inf,
            "Lambda_inf_fit_status": status,
            "Lambda_inf_minus_exact": lam_inf - exact if np.isfinite(lam_inf) and np.isfinite(exact) else float("nan"),
            "Lambda_inf_minus_gaussian": lam_inf - g_base if np.isfinite(lam_inf) and np.isfinite(g_base) else (
                lam_inf if str(fam) == "gaussian" and np.isfinite(lam_inf) else float("nan")
            ),
            "Lambda_inf_minus_separable": lam_inf - s_base if np.isfinite(lam_inf) and np.isfinite(s_base) else float("nan"),
            "fit_slope_noise": slope,
            "fit_r2_scaling": r2,
            "n_M_values": int(np.isfinite(y).sum()),
            "largest_M_mc": largest_M,
            "largest_M_Lambda_hat_full_sym": largest_lam,
            "largest_M_gamma_loc": _last_finite(sub, "gamma_loc"),
            "largest_M_full_sym_minus_gaussian": _last_finite(sub, "full_sym_minus_gaussian"),
            "largest_M_full_sym_minus_separable": _last_finite(sub, "full_sym_minus_separable"),
            "largest_M_full_sym_minus_exact": _last_finite(sub, "full_sym_minus_separable_exact"),
            "largest_M_diag_minus_exact": _last_finite(sub, "diag_minus_exact"),
            "converged_flag": bool(conv),
            "convergence_warning": warning,
        })
    return pd.DataFrame(rows)


def aggregate_convergence_summary(fits):
    """Aggregate per-seed convergence fits to one row per family/N/kappa."""
    if fits is None or fits.empty:
        return pd.DataFrame()
    keys = ["potential_family", "N_theta", "kappa_target"]
    numeric = [
        "p_sym_dim", "M_over_p", "sqrt_p_over_M", "Lambda_inf_fit",
        "Lambda_inf_minus_exact", "Lambda_inf_minus_gaussian",
        "Lambda_inf_minus_separable", "fit_slope_noise", "fit_r2_scaling",
        "n_M_values", "largest_M_mc", "largest_M_Lambda_hat_full_sym",
        "largest_M_gamma_loc", "largest_M_full_sym_minus_gaussian",
        "largest_M_full_sym_minus_separable", "largest_M_full_sym_minus_exact",
        "largest_M_diag_minus_exact",
    ]
    agg = {}
    for col in numeric:
        if col in fits.columns:
            agg[f"{col}_mean"] = (col, "mean")
            if col not in ("p_sym_dim", "n_M_values", "largest_M_mc"):
                agg[f"{col}_std"] = (col, "std")
    agg["n_seeds"] = ("seed", "count")
    agg["all_seeds_converged"] = ("converged_flag", "all")
    agg["any_seed_converged"] = ("converged_flag", "any")
    agg["convergence_warning"] = (
        "convergence_warning",
        lambda s: ";".join(sorted({str(x) for x in s if str(x)})),
    )
    return fits.groupby(keys, dropna=False).agg(**agg).reset_index()
