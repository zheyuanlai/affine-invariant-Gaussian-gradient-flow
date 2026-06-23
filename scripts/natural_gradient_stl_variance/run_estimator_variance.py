"""Experiment A -- estimator-level STL variance comparison.

For every target config and every fixed Gaussian state we draw ``M`` samples and
compare the baseline mean estimator ``b_base = score_post(theta)`` against the STL
estimator ``b_stl = score_post(theta) + C^{-1}(theta - m)``, reporting biases,
Euclidean and Fisher--Rao variances, and the STL/baseline variance ratios. The
sampling/reductions run on the chosen backend (NumPy for ``--device cpu``, a single
CUDA torch device for ``--device cuda``).

Writes into ``--outdir``::

    estimator_variance.csv          one row per (config, state, seed)
    estimator_variance_summary.csv  one row per (config, state): median ratios etc.
    target_metadata.json            target params, a_star, how a_star was computed
    run_metadata.json               config echo, timing, environment

Usage::

    python scripts/natural_gradient_stl_variance/run_estimator_variance.py \
        --config configs/natural_gradient_stl_variance/stl_variance_smoke.yaml \
        --outdir outputs/natural_gradient_stl_variance_smoke --device cpu --overwrite
"""
from __future__ import annotations

import argparse
import os
import platform
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.common.io_utils import load_yaml, ensure_dir, save_dataframe, save_json
from src.common.torch_utils import torch_device_info
from src.natural_gradient_stl_variance.targets import build_target
from src.natural_gradient_stl_variance.states import build_states
from src.natural_gradient_stl_variance.estimator_variance import evaluate_state, make_backend
from src.natural_gradient_stl_variance.grid import enumerate_target_configs

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(__file__), "..", "..", "configs",
    "natural_gradient_stl_variance", "stl_variance.yaml")

LONG_COLS = [
    "target_name", "kind", "d", "kappa", "tau", "state", "distance_to_optimum",
    "M", "seed", "backend", "device",
    "bias_base", "bias_stl",
    "var_eucl_base", "var_eucl_stl", "var_fr_base", "var_fr_stl",
    "ratio_euclidean", "ratio_fr",
    "var_eucl_base_exact", "var_fr_base_exact",
    "var_eucl_stl_exact", "var_fr_stl_exact",
]

SUMMARY_COLS = [
    "target_name", "kind", "d", "kappa", "tau", "state", "distance_to_optimum",
    "n_seeds", "M",
    "ratio_fr_median", "ratio_fr_mean", "ratio_fr_std",
    "ratio_euclidean_median",
    "var_fr_base_median", "var_fr_stl_median",
    "bias_base_median", "bias_stl_median",
    "var_fr_base_exact", "var_fr_stl_exact", "ratio_fr_exact",
]


def _resolve_backend(device, backend):
    if backend != "auto":
        return backend
    return "torch" if str(device).startswith("cuda") else "numpy"


def run(cfg, outdir, device, backend_name):
    bk = make_backend(backend=backend_name, device=device, dtype="float64")
    seeds = [int(s) for s in cfg["seeds"]]
    est = cfg["estimator"]
    M = int(est["M"])
    chunk = int(est.get("chunk_size", 200_000))
    state_params = est.get("states", {})

    long_rows, summary_rows = [], []
    metadata = {"targets": {}, "states": {}}
    configs = enumerate_target_configs(cfg)
    t0 = time.time()
    for ci, (kind, d, kappa, tau, gh) in enumerate(configs):
        target = build_target(kind, d, kappa, tau=tau, n_nodes=gh or 80)
        m_star, C_star = target.a_star()
        states, used_params = build_states(m_star, C_star, state_params)
        metadata["targets"][f"{kind}__d{d}__k{kappa:g}__t{tau:g}"] = target.metadata()
        metadata["states"]["params"] = used_params
        for (state_name, m, C, dist) in states:
            rows = [evaluate_state(target, state_name, m, C, dist, M, seed, bk,
                                   chunk_size=chunk) for seed in seeds]
            long_rows.extend(rows)
            summary_rows.append(_summarize(rows))
        print(f"  [{ci + 1:2d}/{len(configs)}] {kind:8s} d={d:<2d} kappa={kappa:<8g} "
              f"tau={tau:<4g} ({time.time() - t0:5.1f}s)")

    save_dataframe(os.path.join(outdir, "estimator_variance.csv"),
                   pd.DataFrame(long_rows), columns=LONG_COLS)
    save_dataframe(os.path.join(outdir, "estimator_variance_summary.csv"),
                   pd.DataFrame(summary_rows), columns=SUMMARY_COLS)
    # Estimator-specific metadata is prefixed so the estimator and algorithm runs
    # can share an output directory without overwriting each other's provenance.
    save_json(os.path.join(outdir, "estimator_target_metadata.json"), metadata)
    return metadata, time.time() - t0


# Files this script owns; --overwrite clears only these (not a shared directory).
OWN_FILES = [
    "estimator_variance.csv", "estimator_variance_summary.csv",
    "estimator_target_metadata.json", "estimator_run_metadata.json",
]


def _summarize(rows):
    df = pd.DataFrame(rows)
    r0 = rows[0]
    ratio_fr = df["ratio_fr"].to_numpy(dtype=float)
    finite = ratio_fr[np.isfinite(ratio_fr)]
    base_ex = float(r0["var_fr_base_exact"])
    stl_ex = float(r0["var_fr_stl_exact"])
    ratio_ex = (stl_ex / base_ex) if (np.isfinite(base_ex) and base_ex > 0) else float("nan")
    return {
        "target_name": r0["target_name"], "kind": r0["kind"],
        "d": r0["d"], "kappa": r0["kappa"], "tau": r0["tau"],
        "state": r0["state"], "distance_to_optimum": r0["distance_to_optimum"],
        "n_seeds": len(rows), "M": r0["M"],
        "ratio_fr_median": float(np.median(finite)) if finite.size else float("nan"),
        "ratio_fr_mean": float(np.mean(finite)) if finite.size else float("nan"),
        "ratio_fr_std": float(np.std(finite)) if finite.size else float("nan"),
        "ratio_euclidean_median": float(np.nanmedian(df["ratio_euclidean"])),
        "var_fr_base_median": float(np.median(df["var_fr_base"])),
        "var_fr_stl_median": float(np.median(df["var_fr_stl"])),
        "bias_base_median": float(np.median(df["bias_base"])),
        "bias_stl_median": float(np.median(df["bias_stl"])),
        "var_fr_base_exact": base_ex, "var_fr_stl_exact": stl_ex,
        "ratio_fr_exact": ratio_ex,
    }


def parse_args():
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=DEFAULT_CONFIG)
    known, _ = pre.parse_known_args()
    cfg = load_yaml(known.config)
    p = argparse.ArgumentParser(parents=[pre], description=__doc__)
    p.add_argument("--outdir", default=cfg["output_dir"])
    p.add_argument("--device", default="cpu",
                   help="cpu | cuda | cuda:N (cuda requires torch + a CUDA device)")
    p.add_argument("--backend", default="auto", choices=["auto", "numpy", "torch"])
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args(), cfg, known.config


if __name__ == "__main__":
    args, cfg, cfg_path = parse_args()
    backend_name = _resolve_backend(args.device, args.backend)
    print("=" * 64)
    print("STL estimator-level variance (Experiment A)")
    print("=" * 64)
    print(f"  config : {cfg_path}")
    print(f"  device : {args.device}   backend: {backend_name}")
    print(f"  outdir : {args.outdir}\n")

    ensure_dir(args.outdir)
    if args.overwrite:
        for name in OWN_FILES:
            p = os.path.join(args.outdir, name)
            if os.path.exists(p):
                os.remove(p)

    t_start = time.time()
    metadata, run_wall = run(cfg, args.outdir, args.device, backend_name)
    save_json(os.path.join(args.outdir, "estimator_run_metadata.json"), {
        "experiment": "estimator_variance", "config_path": cfg_path, "config": cfg,
        "device": args.device, "backend": backend_name,
        "torch_info": torch_device_info(args.device if backend_name == "torch" else "cpu"),
        "wall_time_total": float(time.time() - t_start),
        "python": platform.python_version(), "platform": platform.platform(),
        "numpy": np.__version__,
    })
    print(f"\nWrote estimator_variance.csv, estimator_variance_summary.csv, "
          f"estimator_target_metadata.json, estimator_run_metadata.json "
          f"-> {args.outdir}")
