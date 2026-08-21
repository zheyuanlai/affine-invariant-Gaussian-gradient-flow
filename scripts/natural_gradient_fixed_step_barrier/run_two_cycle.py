"""Run the exact two-cycle counterexample (CPU, deterministic).

Three studies, all on the same family of targets:

1. **Cycle.** For every ``(kappa, gamma, scheme)`` with ``scheme`` a Fisher--Rao
   scheme, start at ``(M, 1/p)`` and check that the orbit is an exact nonoptimal
   period-two cycle: the mean flips sign every step, the amplitude and the gap are
   constant, and the gap never falls.
2. **Bures--Wasserstein.** The same targets under the BW forward--backward step at
   ``eta = mult / beta``. Inside its certified range (``mult <= 1``) the cycle cannot
   form; the edge is at ``mult = 2``, where ``|1 - eta r| = 1`` at ``r = beta``.
3. **Basin.** The initial mean is perturbed by a relative amount, to measure how
   large a neighbourhood of the orbit fails to converge -- i.e. whether the
   counterexample is robust or a knife edge.

Writes::

    two_cycle_summary.csv     one row per (kappa, gamma, scheme)  [studies 1 and 2]
    two_cycle_basin.csv       one row per (kappa, gamma, scheme, perturbation)
    two_cycle_targets.csv     realized geometry of each target (curvature range, LH)
    two_cycle_long.csv        decimated trajectories
    two_cycle_metadata.json   config and environment

Usage::

    python scripts/natural_gradient_fixed_step_barrier/run_two_cycle.py \
        --config configs/natural_gradient_fixed_step_barrier/two_cycle.yaml \
        --outdir outputs/natural_gradient_fixed_step_barrier --overwrite
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
from src.natural_gradient_fixed_step_barrier.two_cycle import TwoCycleTarget
from src.natural_gradient_fixed_step_barrier.runner import simulate_two_cycle

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "..", "configs",
    "natural_gradient_fixed_step_barrier", "two_cycle.yaml")

SUMMARY_COLS = [
    "kappa", "gamma", "scheme", "family", "dt", "dt_over_certified", "perturb",
    "p", "r", "q", "c_cycle", "M", "initial_gap", "final_gap", "final_rel_gap",
    "min_rel_gap", "max_rel_gap", "amplitude_min", "amplitude_max",
    "n_steps", "sign_flips", "n_tol", "tol_rel", "converged", "status",
]
TARGET_COLS = [
    "kappa", "gamma", "q", "p", "r", "c_cycle", "M", "x1", "x2", "d",
    "width_factor", "V2_min", "V2_max", "kappa_realized", "hessian_lipschitz",
    "A_at_M", "b_over_M_at_M", "cA_at_M", "mean_multiplier",
    "c_star", "energy_star", "gap_at_cycle",
]
LONG_COLS = ["kappa", "gamma", "scheme", "dt", "perturb", "n", "m", "c", "gap",
             "rel_gap", "amplitude", "cA"]


def merge_smoke(cfg):
    out = dict(cfg)
    out.update(cfg.get("smoke") or {})
    out.pop("smoke", None)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    cfg = load_yaml(args.config)
    cfg = merge_smoke(cfg) if args.smoke else {k: v for k, v in cfg.items() if k != "smoke"}
    outdir = args.outdir or cfg["output_dir"]
    ensure_dir(outdir)
    summary_path = os.path.join(outdir, "two_cycle_summary.csv")
    if os.path.exists(summary_path) and not args.overwrite:
        raise SystemExit(f"{summary_path} exists; pass --overwrite")

    tcfg = cfg["target"]
    tol_rel = float(cfg["tol_rel"])
    t0 = time.time()
    summaries, basin_rows, target_rows, long_rows = [], [], [], []

    for kappa in cfg["kappas"]:
        kappa = float(kappa)
        for gamma in cfg["gammas"]:
            gamma = float(gamma)
            target = TwoCycleTarget(kappa, gamma, LH=tcfg["LH"],
                                    width_factor=tcfg["width_factor"],
                                    gh_nodes=tcfg["gh_nodes"])
            md = target.metadata()
            md["width_factor"] = target.width_factor
            target_rows.append(md)
            # The Fisher-Rao certified step at this state, 1/(beta lambda_max) with
            # lambda_max the a priori bound max(c_0, 1/alpha) = 1.
            dt_cert_fr = 1.0 / kappa

            def record(recs, summ, scheme, dt, family, perturb=0.0):
                summ.update(family=family, dt_over_certified=dt / dt_cert_fr,
                            perturb=perturb)
                for r in recs:
                    r.update(kappa=kappa, gamma=gamma, scheme=scheme, dt=dt,
                             perturb=perturb)
                long_rows.extend(recs)
                return summ

            for scheme in cfg["fisher_rao_schemes"]:
                recs, s = simulate_two_cycle(
                    target, scheme, gamma, int(cfg["max_steps"]), tol_rel=tol_rel,
                    max_saved_rows=int(cfg["max_saved_rows"]))
                summaries.append(record(recs, s, scheme, gamma, "fisher_rao"))
                print(f"[cycle] kappa={kappa:6g} gamma={gamma:4g} {scheme:10s} "
                      f"status={s['status']:9s} flips={s['sign_flips']:5d} "
                      f"amp=[{s['amplitude_min']:.6f},{s['amplitude_max']:.6f}] "
                      f"conv={s['converged']}", flush=True)

                for pert in cfg["perturbations"]:
                    pr, ps = simulate_two_cycle(
                        target, scheme, gamma, int(cfg["max_steps_perturb"]),
                        tol_rel=tol_rel, max_saved_rows=20, perturb=float(pert))
                    basin_rows.append(record(pr, ps, scheme, gamma, "fisher_rao",
                                             perturb=float(pert)))

            for mult in cfg["bw_eta_multipliers"]:
                eta = float(mult) / kappa          # beta = kappa
                recs, s = simulate_two_cycle(
                    target, "bures_wasserstein", eta, int(cfg["max_steps"]),
                    tol_rel=tol_rel, max_saved_rows=int(cfg["max_saved_rows"]))
                s["bw_eta_multiplier"] = float(mult)
                summaries.append(record(recs, s, "bures_wasserstein", eta,
                                        "bures_wasserstein"))
                print(f"[ BW  ] kappa={kappa:6g} gamma={gamma:4g} eta={mult:g}/beta   "
                      f"status={s['status']:9s} n_tol={s['n_tol']:6d} "
                      f"conv={s['converged']}", flush=True)

    # ---- barrier sweep: the two counterexamples combined over a range of fixed dt ----
    # For each fixed dt at one kappa: if dt > 2/kappa the two-cycle exists and the
    # hitting time is infinite; otherwise the retuned bump train is run and its
    # measured blocked-iteration count is reported. The envelope is minimized at
    # dt ~ 1/kappa, which is where the Theta(kappa^2) comes from.
    from src.natural_gradient_fixed_step_barrier.bump_target import BumpTrain
    from src.natural_gradient_fixed_step_barrier.runner import simulate as simulate_bump

    barrier_rows = []
    bk = float(cfg.get("barrier_kappa", 256))
    for dt in cfg.get("barrier_dt_grid", []):
        dt = float(dt)
        cycles = dt > 2.0 / bk
        row = {"kappa": bk, "dt": dt, "dt_times_kappa": dt * bk,
               "two_cycle_exists": int(cycles)}
        if cycles:
            row.update(n_half=float("inf"), bump_N_train=float("nan"),
                       regime="two_cycle", scheme="both")
        else:
            train = BumpTrain(bk, dt / bk)
            _, s = simulate_bump(train, "riemannian", "const", dt,
                                 max_steps=max(20000, 8 * train.N), tol_rel=1e-6)
            row.update(n_half=float(s["n_half"]), bump_N_train=float(train.N),
                       regime="bump_train", scheme="riemannian")
        barrier_rows.append(row)
        print(f"[barrier] dt={dt:.6g} (dt*kappa={dt*bk:8.3f})  regime={row['regime']:10s} "
              f"n_half={row['n_half']}", flush=True)
    df_barrier = pd.DataFrame(barrier_rows)

    df_sum = pd.DataFrame(summaries)
    keep = [c for c in SUMMARY_COLS if c in df_sum.columns] + \
           [c for c in df_sum.columns if c not in SUMMARY_COLS]
    df_sum = df_sum[keep]
    df_basin = pd.DataFrame(basin_rows).reindex(columns=SUMMARY_COLS)
    df_tgt = pd.DataFrame(target_rows).reindex(columns=TARGET_COLS)
    df_long = pd.DataFrame(long_rows).reindex(columns=LONG_COLS)

    save_dataframe(summary_path, df_sum)
    save_dataframe(os.path.join(outdir, "two_cycle_basin.csv"), df_basin)
    save_dataframe(os.path.join(outdir, "two_cycle_barrier.csv"), df_barrier)
    save_dataframe(os.path.join(outdir, "two_cycle_targets.csv"), df_tgt)
    save_dataframe(os.path.join(outdir, "two_cycle_long.csv"), df_long)
    save_json(os.path.join(outdir, "two_cycle_metadata.json"), {
        "config": cfg, "config_path": os.path.abspath(args.config),
        "smoke": bool(args.smoke), "wall_time_seconds": time.time() - t0,
        "numpy_version": np.__version__,
    })

    fr = df_sum[df_sum.family == "fisher_rao"]
    bw = df_sum[df_sum.family == "bures_wasserstein"]
    print(f"\nFisher-Rao runs:        {len(fr):3d}, converged: {int(fr.converged.sum())}")
    print(f"Bures-Wasserstein runs: {len(bw):3d}, converged: {int(bw.converged.sum())}")
    print("\nBW convergence by stepsize multiple of 1/beta:")
    print(bw.groupby("bw_eta_multiplier")["converged"].agg(["sum", "count"]).to_string())
    print(f"\nwrote {outdir} in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
