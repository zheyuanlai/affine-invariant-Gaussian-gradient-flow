"""
Extract Markdown result tables from existing experiment summaries.

Read-only utility: loads already-computed CSV/JSON outputs and prints
GitHub-flavoured Markdown tables for inclusion in README.md.
Never runs dynamics or regenerates experiment data.

Usage:
    python scripts/extract_readme_results.py

Reads:
    outputs/gaussian_grid/summary.csv
    outputs/logconcave_grid/summary.csv
    outputs/logconcave_grid/target_metadata.json
    outputs/logconcave_grid/reference_optimum_meta.json
"""
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

INIT_ORDER = ["mean_only", "volume_high", "volume_low", "shape_only", "mixed"]
INIT_LABEL = {
    "mean_only":   "mean-only",
    "volume_high": "volume-high",
    "volume_low":  "volume-low",
    "shape_only":  "shape-only",
    "mixed":       "mixed",
}


def fmt_time(v, T):
    if v is None or not np.isfinite(v):
        return f">{T:g}"
    return f"{v:.1f}"


def fmt_ratio(num, den):
    if not np.isfinite(num) or not np.isfinite(den) or den == 0:
        return "—"
    return f"{num / den:.2f}"


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def time_to_tol_table(df, tol_col, omegas, T):
    header = ["Initialization"] + [f"$\\omega={o:g}$" for o in omegas]
    rows = []
    for init in INIT_ORDER:
        row = [INIT_LABEL[init]]
        for o in omegas:
            sel = df[(df.init_name == init) & np.isclose(df.omega, o)
                     & (df.tau_type == "zero")]
            v = float(sel[tol_col].iloc[0]) if not sel.empty else np.nan
            row.append(fmt_time(v, T))
        rows.append(row)
    return md_table(header, rows)


def speedup_table_fixed_omega(df, tol_col, omega):
    header = ["Initialization",
              "$T(\\tau{=}0)$", "$T(\\tau_-)/T_0$", "$T(\\tau_+)/T_0$"]
    rows = []
    for init in INIT_ORDER:
        sel = df[(df.init_name == init) & np.isclose(df.omega, omega)]
        t0 = sel[sel.tau_type == "zero"][tol_col]
        tn = sel[sel.tau_type == "negative"][tol_col]
        tp = sel[sel.tau_type == "positive"][tol_col]
        t0 = float(t0.iloc[0]) if not t0.empty else np.nan
        tn = float(tn.iloc[0]) if not tn.empty else np.nan
        tp = float(tp.iloc[0]) if not tp.empty else np.nan
        T  = float(sel["T"].iloc[0]) if not sel.empty else np.nan
        rows.append([INIT_LABEL[init], fmt_time(t0, T),
                     fmt_ratio(tn, t0), fmt_ratio(tp, t0)])
    return md_table(header, rows)


def speedup_table_across_omega(df, tol_col, omegas):
    header = ["Initialization"]
    for o in omegas:
        header += [f"$\\tau_-$ ($\\omega={o:g}$)", f"$\\tau_+$ ($\\omega={o:g}$)"]
    rows = []
    for init in INIT_ORDER:
        row = [INIT_LABEL[init]]
        for o in omegas:
            sel = df[(df.init_name == init) & np.isclose(df.omega, o)]
            t0 = sel[sel.tau_type == "zero"][tol_col]
            tn = sel[sel.tau_type == "negative"][tol_col]
            tp = sel[sel.tau_type == "positive"][tol_col]
            t0 = float(t0.iloc[0]) if not t0.empty else np.nan
            tn = float(tn.iloc[0]) if not tn.empty else np.nan
            tp = float(tp.iloc[0]) if not tp.empty else np.nan
            row += [fmt_ratio(tn, t0), fmt_ratio(tp, t0)]
        rows.append(row)
    return md_table(header, rows)


def chi_table(df, omega):
    header = ["Initialization", "$\\chi_0$ (initial trace dominance)"]
    rows = []
    for init in INIT_ORDER:
        sel = df[(df.init_name == init) & np.isclose(df.omega, omega)
                 & (df.tau_type == "zero")]
        v = float(sel["initial_chi"].iloc[0]) if not sel.empty else np.nan
        rows.append([INIT_LABEL[init], f"{v:.3f}"])
    return md_table(header, rows)


def main():
    g  = pd.read_csv(os.path.join(ROOT, "outputs/gaussian_grid/summary.csv"))
    lc = pd.read_csv(os.path.join(ROOT, "outputs/logconcave_grid/summary.csv"))
    with open(os.path.join(ROOT, "outputs/logconcave_grid/target_metadata.json")) as f:
        meta = json.load(f)
    with open(os.path.join(ROOT, "outputs/logconcave_grid/reference_optimum_meta.json")) as f:
        ref = json.load(f)

    g5 = g[g.n == 5]
    g_omegas = sorted(g5.omega.unique())
    T_g = float(g5["T"].iloc[0])
    lc_omegas = sorted(lc.omega.unique())
    T_lc = float(lc["T"].iloc[0])

    print("=" * 70)
    print("GAUSSIAN TARGET  (n = 5)")
    print("=" * 70)
    print("\n### G1. Time to E/E0 <= 1e-4, tau = 0\n")
    print(time_to_tol_table(g5, "time_to_1e_minus_4", g_omegas, T_g))
    print("\n### G2. tau-speedup T(tau)/T(0), tol 1e-4, omega = 0.5\n")
    print(speedup_table_fixed_omega(g5, "time_to_1e_minus_4", 0.5))
    print("\n### G3. tau-speedup T(tau)/T(0), tol 1e-4, omega in {0.25,0.5,1.0}\n")
    print(speedup_table_across_omega(g5, "time_to_1e_minus_4", [0.25, 0.5, 1.0]))

    print("\n" + "=" * 70)
    print("LOG-CONCAVE TARGET  (n = 5, rho = 5)")
    print("=" * 70)
    print(f"\nMetadata: n={meta['n']}, rho={meta['rho']}, m={meta['m_features']}, "
          f"K={meta['K']}, K_ref={meta['K_ref']}, dt={meta['dt']}, T={meta['T']}")
    print(f"F_star={meta['F_star']:.6f}, ||m_star||={meta['m_star_norm']:.3e}")
    print(f"grad_m_norm={ref['grad_m_norm']:.3e}, "
          f"cov_residual_norm={ref['cov_residual_norm']:.6f}, "
          f"converged={ref['converged']}")
    print("\n### L1. Time to normalized gap <= 1e-2, tau = 0\n")
    print(time_to_tol_table(lc, "time_to_1e_minus_2", lc_omegas, T_lc))
    print("\n### L2. Time to normalized gap <= 1e-4, tau = 0 (near objective floor)\n")
    print(time_to_tol_table(lc, "time_to_1e_minus_4", lc_omegas, T_lc))
    print("\n### L3. tau-speedup T(tau)/T(0), tol 1e-2, omega = 0.5\n")
    print(speedup_table_fixed_omega(lc, "time_to_1e_minus_2", 0.5))
    print("\n### L4. tau-speedup T(tau)/T(0), tol 1e-2, omega in {0.25,0.5,1.0}\n")
    print(speedup_table_across_omega(lc, "time_to_1e_minus_2", [0.25, 0.5, 1.0]))
    print("\n### L5. Initial trace dominance chi, omega = 0.5\n")
    print(chi_table(lc, 0.5))


if __name__ == "__main__":
    sys.exit(main())
