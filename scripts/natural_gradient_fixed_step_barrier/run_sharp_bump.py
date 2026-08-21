"""Run the sharp bump-train experiment (CPU, deterministic).

For every ``(family, kappa, scheme, arm)`` the runner builds the Appendix C
potential ``V_kappa``, starts at ``(m, c) = (x_0, 1/kappa)`` and iterates the
scalar scheme until the objective gap has fallen by ``tol_rel``, then fits
``log n_half`` against ``log kappa``. It writes four files:

    sharp_bump_long.csv        one row per saved step per run
    sharp_bump_summary.csv     one row per (family, kappa, scheme, arm)
    sharp_bump_slopes.csv      log-log slope of n_half vs kappa, per (family, scheme, arm)
    sharp_bump_gamma_sweep.csv the fixed order-one step swept over its own constant
    sharp_bump_metadata.json   config, per-train geometry, environment

Usage::

    python scripts/natural_gradient_fixed_step_barrier/run_sharp_bump.py \
        --config configs/natural_gradient_fixed_step_barrier/sharp_bump.yaml \
        --outdir outputs/natural_gradient_fixed_step_barrier --overwrite

Add ``--smoke`` for the reduced grid defined in the config.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.io_utils import load_yaml, ensure_dir, save_dataframe, save_json
from src.natural_gradient_fixed_step_barrier.bump_target import BumpTrain
from src.natural_gradient_fixed_step_barrier.runner import nominal_dt, simulate

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "..", "configs",
    "natural_gradient_fixed_step_barrier", "sharp_bump.yaml")

LONG_COLS = ["family", "kappa", "scheme", "arm", "n", "m", "c", "gap", "rel_gap", "cA", "A"]

SUMMARY_COLS = [
    "family", "kappa", "scheme", "arm", "gamma", "dt_nominal", "dt_mean",
    "s_train", "N_train", "x0", "initial_gap", "final_gap", "final_rel_gap",
    "n_steps", "n_half", "n_tol", "tol_rel", "n_half_over_N_train",
    "shadow_mean_err_over_w", "shadow_cov_err", "cA_min", "cA_max",
    "c_min", "c_max", "max_energy_increase", "monotone", "status",
]

SLOPE_COLS = ["family", "scheme", "arm", "metric", "n_points", "kappa_min",
              "kappa_max", "slope", "intercept", "r2",
              "slope_upper", "kappa_min_upper", "n_points_upper",
              "predicted_exponent"]


def merge_smoke(cfg):
    """Shallow-merge the ``smoke`` block over the base config."""
    smoke = cfg.get("smoke") or {}
    out = dict(cfg)
    out.update(smoke)
    out.pop("smoke", None)
    return out


def train_gain(family, arm, gamma, kappa):
    """Per-step mean gain ``s = dt * c_kappa`` the train is built against."""
    dt = (gamma / kappa) if family == "manuscript" else nominal_dt(arm, gamma, kappa)
    return dt / kappa


def loglog_slope(kappas, values):
    """Least-squares slope/intercept/R^2 of ``log values`` on ``log kappas``."""
    k = np.asarray(kappas, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)
    ok = np.isfinite(k) & np.isfinite(v) & (k > 0) & (v > 0)
    if int(ok.sum()) < 2:
        return float("nan"), float("nan"), float("nan"), int(ok.sum())
    x, y = np.log(k[ok]), np.log(v[ok])
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), float(r2), int(ok.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    cfg = load_yaml(args.config)
    if args.smoke:
        cfg = merge_smoke(cfg)
    else:
        cfg.pop("smoke", None)
    outdir = args.outdir or cfg["output_dir"]
    ensure_dir(outdir)

    summary_path = os.path.join(outdir, "sharp_bump_summary.csv")
    if os.path.exists(summary_path) and not args.overwrite:
        raise SystemExit(f"{summary_path} exists; pass --overwrite")

    gamma = float(cfg["gamma"])
    tcfg = cfg["train"]
    t0 = time.time()
    long_rows, summaries, trains = [], [], {}

    for family in cfg["families"]:
        for kappa in cfg["kappas"]:
            kappa = float(kappa)
            for arm in cfg["arms"]:
                s = train_gain(family, arm, gamma, kappa)
                key = (family, kappa, arm)
                train = BumpTrain(kappa, s, T=tcfg["T"], Y=tcfg["Y"],
                                  LH=tcfg["LH"], gh_nodes=tcfg["gh_nodes"])
                trains[f"{family}|kappa={kappa:g}|{arm}"] = train.metadata()
                max_steps = max(int(cfg["max_steps_floor"]),
                                int(cfg["max_steps_factor"]) * train.N)
                for scheme in cfg["schemes"]:
                    recs, summ = simulate(
                        train, scheme, arm, gamma, max_steps,
                        tol_rel=float(cfg["tol_rel"]),
                        max_saved_rows=int(cfg["max_saved_rows"]))
                    summ["family"] = family
                    summaries.append(summ)
                    for r in recs:
                        r.update(family=family, kappa=kappa, scheme=scheme, arm=arm)
                    long_rows.extend(recs)
                    print(f"[{family:10s}] kappa={kappa:6g} {scheme:10s} {arm:8s} "
                          f"dt={summ['dt_nominal']:.5g} N_train={train.N:6d} "
                          f"n_half={summ['n_half']:7d} n_tol={summ['n_tol']:7d} "
                          f"shadow={summ['shadow_mean_err_over_w']:.2e} "
                          f"{summ['status']}", flush=True)

    df_long = pd.DataFrame(long_rows).reindex(columns=LONG_COLS)
    df_sum = pd.DataFrame(summaries).reindex(columns=SUMMARY_COLS)

    slope_rows = []
    for (family, scheme, arm), grp in df_sum.groupby(["family", "scheme", "arm"]):
        grp = grp.sort_values("kappa")
        # The construction is asymptotic in kappa, so also fit the upper half of the
        # grid: the full-range slope is dragged down by the small-kappa points.
        upper = grp.tail(max(2, len(grp) // 2))
        for metric in ("n_half", "n_tol"):
            slope, intercept, r2, npts = loglog_slope(grp["kappa"], grp[metric])
            slope_u, _, _, npts_u = loglog_slope(upper["kappa"], upper[metric])
            slope_rows.append({
                "family": family, "scheme": scheme, "arm": arm, "metric": metric,
                "n_points": npts,
                "kappa_min": float(grp["kappa"].min()), "kappa_max": float(grp["kappa"].max()),
                "slope": slope, "intercept": intercept, "r2": r2,
                "slope_upper": slope_u, "n_points_upper": npts_u,
                "kappa_min_upper": float(upper["kappa"].min()),
                # thm:sharp-disc certifies Omega(T kappa / dt): 2 at dt ~ 1/kappa, 1 at dt ~ 1.
                "predicted_exponent": 2.0 if arm == "theory" else 1.0,
            })
    df_slope = pd.DataFrame(slope_rows).reindex(columns=SLOPE_COLS)

    # Robustness of the FIXED order-one step in its own constant gamma. Each gamma
    # gets its own retuned train (s = gamma/kappa) and only the const arm is run.
    gamma_rows = []
    for gamma_s in cfg.get("gamma_sweep", []):
        gamma_s = float(gamma_s)
        for kappa in cfg["kappas"]:
            kappa = float(kappa)
            train = BumpTrain(kappa, gamma_s / kappa, T=tcfg["T"], Y=tcfg["Y"],
                              LH=tcfg["LH"], gh_nodes=tcfg["gh_nodes"])
            max_steps = max(int(cfg["max_steps_floor"]),
                            int(cfg["max_steps_factor"]) * train.N)
            for scheme in cfg["schemes"]:
                _, summ = simulate(train, scheme, "const", gamma_s, max_steps,
                                   tol_rel=float(cfg["tol_rel"]),
                                   max_saved_rows=2)
                summ["family"] = "retuned"
                gamma_rows.append(summ)
    df_gamma = pd.DataFrame(gamma_rows).reindex(columns=SUMMARY_COLS)
    if not df_gamma.empty:
        gslopes = []
        for (gamma_s, scheme), grp in df_gamma.groupby(["gamma", "scheme"]):
            grp = grp.sort_values("kappa")
            upper = grp.tail(max(2, len(grp) // 2))
            slope, _, r2, npts = loglog_slope(grp["kappa"], grp["n_half"])
            slope_u, _, _, _ = loglog_slope(upper["kappa"], upper["n_half"])
            gslopes.append({
                "gamma": gamma_s, "scheme": scheme, "n_points": npts,
                "slope": slope, "r2": r2, "slope_upper": slope_u,
                "all_monotone": int(grp["monotone"].min()),
                "all_converged": int((grp["status"] == "ok").all()),
                "statuses": "|".join(sorted(set(grp["status"]))),
            })
        df_gamma = df_gamma.merge(pd.DataFrame(gslopes), on=["gamma", "scheme"],
                                  how="left", suffixes=("", "_fit"))

    save_dataframe(os.path.join(outdir, "sharp_bump_long.csv"), df_long)
    save_dataframe(summary_path, df_sum)
    save_dataframe(os.path.join(outdir, "sharp_bump_slopes.csv"), df_slope)
    save_dataframe(os.path.join(outdir, "sharp_bump_gamma_sweep.csv"), df_gamma)
    save_json(os.path.join(outdir, "sharp_bump_metadata.json"), {
        "config": cfg, "config_path": os.path.abspath(args.config),
        "smoke": bool(args.smoke), "trains": trains,
        "wall_time_seconds": time.time() - t0,
        "numpy_version": np.__version__,
    })

    print("\nlog-log slope of iterations vs kappa (n_half = constant-factor reduction):")
    show = ["family", "scheme", "arm", "slope", "r2", "slope_upper",
            "kappa_min_upper", "predicted_exponent"]
    print(df_slope[df_slope["metric"] == "n_half"][show].to_string(index=False))
    if not df_gamma.empty:
        print("\nfixed order-one step swept over its own constant gamma "
              "(retuned train, const arm):")
        gshow = ["gamma", "scheme", "slope", "slope_upper", "all_monotone",
                 "all_converged", "statuses"]
        print(df_gamma[gshow].drop_duplicates().sort_values(["gamma", "scheme"])
              .to_string(index=False))
    print(f"\nwrote {outdir} in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
