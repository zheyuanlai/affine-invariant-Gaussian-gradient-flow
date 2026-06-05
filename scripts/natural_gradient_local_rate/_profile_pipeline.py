"""TEMPORARY profiling harness for the low-dimensional high-M torch path.

Replicates the joint-runner per-row loop (bank generation -> CPU potential
construction -> compute_row) but wraps each phase in a CUDA-synchronized timer so
we can see where wall time actually goes. Prints a per-row table + per-phase
totals. This is a benchmarking tool, not production; safe to delete.

Usage:
    python scripts/natural_gradient_local_rate/_profile_pipeline.py \
        --config configs/natural_gradient_local_rate/_profile_lowdim.yaml \
        --device cuda --dtype float64 [--chunk-size N]
"""
import argparse
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", ".."))
sys.path.insert(0, _HERE)

import _common  # noqa: E402
from src.natural_gradient_local_rate.estimator_suite import compute_row  # noqa: E402


def _sync():
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.synchronize()
    except Exception:
        pass


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--device", default="cuda")
    p.add_argument("--dtype", default="float64")
    p.add_argument("--backend", default="torch")
    p.add_argument("--chunk-size", type=int, default=None)
    args = p.parse_args()

    cfg = _common.load_config(args.config)
    opts = _common.operator_opts(cfg)
    opts["compute_gamma_loc"] = True
    opts["backend"] = args.backend
    opts["device"] = args.device
    opts["dtype"] = args.dtype
    if args.chunk_size is not None:
        opts["chunk_size"] = args.chunk_size
    run_id, group = _common.run_context()

    points = list(_common.grid_points(cfg))
    print(f"[profile backend={opts['backend']} device={opts['device']} "
          f"dtype={opts['dtype']} chunk_size={opts['chunk_size']}] {len(points)} rows\n")

    header = (f"{'family':16s} {'N':>3s} {'bank_s':>8s} {'pot_s':>9s} "
              f"{'row_s':>8s} {'mat_s':>8s} {'eigh_s':>8s} {'rt_int':>8s} "
              f"{'peakMB':>9s} {'total_s':>8s}")
    print(header)
    print("-" * len(header))

    tot = {"bank": 0.0, "pot": 0.0, "row": 0.0, "mat": 0.0, "eigh": 0.0, "total": 0.0}
    rows_meta = []
    for point in points:
        point["M_mc"] = _common.grid_M_mc(cfg, point)

        _sync(); t = time.time()
        Z = _common.make_bank(cfg, point)
        _sync(); bank_s = time.time() - t

        _sync(); t = time.time()
        pot = _common.make_potential_for_opts(cfg, point, Z, opts)
        _sync(); pot_s = time.time() - t

        _sync(); t = time.time()
        row = compute_row(pot, Z, point, opts, run_id=run_id, experiment_group=group)
        _sync(); row_s = time.time() - t

        mat_s = float(row.get("matrix_construction_runtime_seconds", float("nan")))
        eigh_s = float(row.get("eigh_runtime_seconds", float("nan")))
        rt_int = float(row.get("runtime_seconds", float("nan")))
        peak = float(row.get("gpu_peak_memory_mb", float("nan")))
        total_s = bank_s + pot_s + row_s
        if row["status"] != "ok":
            print(f"  ERROR {point['family']} N={point['N_theta']}: {row['error_message']}")
        print(f"{point['family']:16s} {point['N_theta']:3d} {bank_s:8.3f} {pot_s:9.3f} "
              f"{row_s:8.3f} {mat_s:8.3f} {eigh_s:8.3f} {rt_int:8.3f} "
              f"{peak:9.1f} {total_s:8.3f}")
        tot["bank"] += bank_s; tot["pot"] += pot_s; tot["row"] += row_s
        tot["mat"] += mat_s; tot["eigh"] += eigh_s; tot["total"] += total_s
        rows_meta.append((point["family"], point["N_theta"], bank_s, pot_s, row_s))

    print("-" * len(header))
    n = len(points)
    print(f"{'TOTAL':16s} {'':>3s} {tot['bank']:8.3f} {tot['pot']:9.3f} "
          f"{tot['row']:8.3f} {tot['mat']:8.3f} {tot['eigh']:8.3f} {'':>8s} "
          f"{'':>9s} {tot['total']:8.3f}")
    print(f"\nPer-phase share of wall time (total {tot['total']:.2f}s over {n} rows):")
    for k in ("bank", "pot", "row"):
        print(f"  {k:6s}: {tot[k]:8.2f}s  ({100*tot[k]/max(tot['total'],1e-9):5.1f}%)  "
              f"avg {tot[k]/n:.3f}s/row")
    print(f"  (within row: matrix_construction {tot['mat']:.2f}s, eigh {tot['eigh']:.2f}s)")


if __name__ == "__main__":
    main()
