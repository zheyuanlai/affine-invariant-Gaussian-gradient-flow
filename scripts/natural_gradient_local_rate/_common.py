"""Shared glue for the natural-gradient local-rate runner scripts.

Centralizes config loading, grid iteration, deterministic sample-bank seeding,
and potential construction so that the operator / linearized-rate / flow scripts
all build identical banks and potentials for a given grid point.
"""
from __future__ import annotations

import itertools
import os

import numpy as np

from src.common.io_utils import load_yaml, ensure_dir, run_id as _make_run_id
from src.common.monte_carlo import gaussian_samples
from src.natural_gradient_local_rate.potentials import build_potential
from src.natural_gradient_local_rate.estimator_suite import (
    default_opts, EXPERIMENT_GROUP,
)

FAMILY_CODE = {
    "gaussian": 0, "separable": 1, "additive_index": 2,
    "random_feature": 3, "radial_tail": 4, "product_feature": 5,
}

_BANK_BASE = 20260603  # fixed base so banks are reproducible across stages

# Canonical required output schema (Part 8 of the spec), in order. The runner
# scripts reindex to REQUIRED_COLUMNS first, then append any extra columns
# (potential metadata, reference quantities, the legacy ``Lambda_hat`` alias)
# so old columns are preserved rather than dropped.
REQUIRED_COLUMNS = [
    "run_id", "experiment_group", "potential_family", "seed", "N_theta",
    "kappa_target", "rho", "M_mc", "operator_estimator", "backend", "chunk_size",
    "Lambda_hat_raw_forward", "Lambda_hat_full_sym", "Lambda_hat_diag",
    "Lambda_hat_separable_exact", "full_over_diag", "full_minus_diag",
    "diag_minus_exact", "tau_H", "tau_H_sq", "coupling_bound_rate",
    "tau_top_total", "tau_top_longitudinal", "tau_top_mixed",
    "tau_top_transverse", "tau_top_longitudinal_fraction",
    "tau_top_mixed_fraction", "tau_top_transverse_fraction",
    "tau_top_X_longitudinal_norm_sq", "tau_top_X_mixed_norm_sq",
    "tau_top_X_transverse_norm_sq", "tau_top_svd_gap",
    "gamma_over_coupling_bound", "gamma_loc", "inverse_gamma_loc",
    "self_adjoint_error_H_raw", "self_adjoint_error_H_sym",
    "self_adjoint_error_L_star", "eig_residual_H", "eig_residual_L_star",
    "current_bound_rate", "conjecture_bound_rate", "lambda_over_logkappa",
    "inverse_gamma_over_logkappa", "runtime_seconds", "status", "error_message",
]


def load_config(path):
    return load_yaml(path)


def run_context():
    """A ``(run_id, experiment_group)`` pair shared by all rows of one invocation."""
    return _make_run_id("nglr"), EXPERIMENT_GROUP


def order_columns(df):
    """Put the required columns first (preserving extras) for a results frame."""
    extras = [c for c in df.columns if c not in REQUIRED_COLUMNS and not c.startswith("_")]
    cols = [c for c in REQUIRED_COLUMNS if c in df.columns] + extras
    return df.reindex(columns=cols)


def operator_opts(cfg):
    """Estimator options, reading the ``operator:`` block with legacy fallbacks.

    Older configs put ``chunk_size`` under ``monte_carlo`` and the eigensolver
    tolerances under ``eigsh``; those are honored when no ``operator:`` block is
    present so existing configs keep working.
    """
    o = default_opts()
    op = cfg.get("operator", {}) or {}
    mc = cfg.get("monte_carlo", {}) or {}
    eg = cfg.get("eigsh", {}) or {}
    o["estimator"] = str(op.get("estimator", o["estimator"]))
    o["diagonal_benchmark"] = bool(op.get("diagonal_benchmark", o["diagonal_benchmark"]))
    o["separable_exact_benchmark"] = bool(
        op.get("separable_exact_benchmark", o["separable_exact_benchmark"]))
    o["compute_gamma_loc"] = bool(op.get("compute_gamma_loc", o["compute_gamma_loc"]))
    o["backend"] = str(op.get("backend", o["backend"]))
    o["chunk_size"] = op.get("chunk_size", mc.get("chunk_size", o["chunk_size"]))
    o["eigsh_tol"] = float(op.get("eigsh_tol", eg.get("tol", o["eigsh_tol"])))
    o["eigsh_maxiter"] = int(op.get("eigsh_maxiter", eg.get("maxiter", o["eigsh_maxiter"])))
    o["quadrature_nodes"] = int(op.get("quadrature_nodes", o["quadrature_nodes"]))
    o["self_adjoint_probes"] = int(op.get("self_adjoint_probes", o["self_adjoint_probes"]))
    # --- torch / GPU backend keys (ignored by the numpy path) ---
    o["device"] = str(op.get("device", o["device"]))
    o["dtype"] = str(op.get("dtype", o["dtype"]))
    o["eigensolver"] = str(op.get("eigensolver", o["eigensolver"]))
    o["explicit_dense_max_N_theta"] = int(
        op.get("explicit_dense_max_N_theta", o["explicit_dense_max_N_theta"]))
    o["basis_block_size"] = int(op.get("basis_block_size", o["basis_block_size"]))
    return o


def add_backend_cli_args(parser):
    """Add the shared backend/device CLI overrides to an argparse parser."""
    parser.add_argument("--backend", choices=["numpy", "torch", "auto"], default=None,
                        help="Override operator.backend")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default=None,
                        help="Override operator.device (torch backend)")
    parser.add_argument("--dtype", choices=["float64", "float32"], default=None,
                        help="Override operator.dtype (torch backend)")
    parser.add_argument("--chunk-size", type=int, default=None,
                        help="Override operator.chunk_size")
    parser.add_argument("--basis-block-size", type=int, default=None,
                        help="Override operator.basis_block_size (torch dense builder)")
    parser.add_argument("--explicit-dense-max-N-theta", type=int, default=None,
                        dest="explicit_dense_max_N_theta",
                        help="Override operator.explicit_dense_max_N_theta (torch dense path)")
    parser.add_argument("--num-shards", type=int, default=1,
                        help="Split grid points into this many deterministic shards")
    parser.add_argument("--shard-index", type=int, default=0,
                        help="Run only this zero-based shard index")
    return parser


def apply_cli_overrides(opts, args):
    """Apply non-None CLI overrides onto an ``opts`` dict (CLI wins over YAML)."""
    if getattr(args, "backend", None) is not None:
        opts["backend"] = args.backend
    if getattr(args, "device", None) is not None:
        opts["device"] = args.device
    if getattr(args, "dtype", None) is not None:
        opts["dtype"] = args.dtype
    if getattr(args, "chunk_size", None) is not None:
        opts["chunk_size"] = args.chunk_size
    if getattr(args, "basis_block_size", None) is not None:
        opts["basis_block_size"] = args.basis_block_size
    if getattr(args, "explicit_dense_max_N_theta", None) is not None:
        opts["explicit_dense_max_N_theta"] = args.explicit_dense_max_N_theta
    return opts


def grid_points(cfg):
    """Yield grid-point dicts (Cartesian product of the ``grid`` section)."""
    g = cfg["grid"]
    for N, kappa, fam, seed in itertools.product(
        g["N_theta"], g["kappa_target"], g["potential_family"], g["seeds"],
    ):
        yield {
            "N_theta": int(N),
            "kappa_target": float(kappa),
            "family": str(fam),
            "seed": int(seed),
        }


def apply_grid_shard(points, args):
    """Return the deterministic shard of ``points`` selected by CLI args."""
    num = int(getattr(args, "num_shards", 1))
    idx = int(getattr(args, "shard_index", 0))
    if num <= 0:
        raise ValueError(f"--num-shards must be positive, got {num}")
    if idx < 0 or idx >= num:
        raise ValueError(f"--shard-index must be in [0, {num}), got {idx}")
    if num == 1:
        return list(points)
    return [pt for k, pt in enumerate(points) if k % num == idx]


def shard_suffix(args):
    """Filesystem suffix for sharded outputs, or empty string for unsharded runs."""
    num = int(getattr(args, "num_shards", 1))
    idx = int(getattr(args, "shard_index", 0))
    return "" if num == 1 else f"_shard{idx:02d}-of-{num:02d}"


def point_key(point):
    """Filesystem-safe identifier for a grid point."""
    return (f"{point['family']}_N{point['N_theta']}"
            f"_k{point['kappa_target']:g}_s{point['seed']}")


def bank_seed(point):
    """Deterministic sample-bank seed (shared across stages for a point)."""
    s = _BANK_BASE
    s = s * 1000003 + point["seed"]
    s = s * 1000003 + point["N_theta"]
    s = s * 1000003 + int(round(point["kappa_target"] * 1000))
    s = s * 1000003 + FAMILY_CODE.get(point["family"], 99)
    return s % (2 ** 31)


def grid_M_mc(cfg, point):
    """Resolve ``M_mc`` for a point: per-point value wins, else ``monte_carlo.M_mc``."""
    if point.get("M_mc") is not None:
        return int(point["M_mc"])
    return int(cfg["monte_carlo"]["M_mc"])


def make_bank(cfg, point):
    mc = cfg["monte_carlo"]
    return gaussian_samples(point["N_theta"], grid_M_mc(cfg, point),
                            seed=bank_seed(point),
                            antithetic=bool(mc.get("antithetic", True)))


def scaling_points(cfg):
    """Yield grid points for sample-size scaling: base grid x each ``M_mc``.

    ``monte_carlo.M_mc`` is a *list* here; every base grid point is repeated once
    per ``M_mc`` value (with that value stored on the point).
    """
    m_list = cfg["monte_carlo"]["M_mc"]
    m_list = m_list if isinstance(m_list, (list, tuple)) else [m_list]
    for base in grid_points(cfg):
        for m in m_list:
            pt = dict(base)
            pt["M_mc"] = int(m)
            yield pt


def make_potential(cfg, point, Z):
    pot = cfg.get("potential", {})
    return build_potential(
        point["family"], point["N_theta"], point["kappa_target"], point["seed"],
        feature_multiplier=int(pot.get("feature_multiplier", 4)),
        phi=str(pot.get("phi", "log_cosh")),
        safety_factor=float(pot.get("empirical_LA_safety_factor", 2.0)),
        Z_ref=Z,
    )


def make_potential_for_opts(cfg, point, Z, opts):
    """Build the potential, running centering on the GPU when the torch/CUDA path is active.

    For ``backend="torch"`` (or ``auto`` resolving to CUDA) the expensive bank
    reductions (``E[grad Phi]``, ``E[Hess Phi]``, eigenvalue extremes) run on the
    device via :func:`torch_backend.build_centered_potential_gpu`, which returns
    the *identical* NumPy :class:`CenteredPotential` (byte-for-byte the CPU build
    on the same bank). The numpy/CPU path and the Gaussian family always use the
    plain CPU :func:`make_potential`. Any failure falls back to the CPU builder so
    a row is never lost to a centering hiccup.
    """
    backend = str(opts.get("backend", "numpy")).lower()
    family = point["family"]
    if backend == "numpy" or family == "gaussian":
        return make_potential(cfg, point, Z)
    try:
        from src.common.torch_utils import resolve_backend, resolve_device, resolve_dtype
        if resolve_backend(backend, opts.get("device", "auto")) != "torch":
            return make_potential(cfg, point, Z)
        from src.natural_gradient_local_rate import torch_backend as tb
        if not tb.torch_supports_family(family):
            return make_potential(cfg, point, Z)
        device = resolve_device(opts.get("device", "auto"))
        if device.type != "cuda":  # torch-on-CPU centering is no faster than numpy
            return make_potential(cfg, point, Z)
        dtype = resolve_dtype(opts.get("dtype", "float64"))
        pot = cfg.get("potential", {})
        return tb.build_centered_potential_gpu(
            family, point, pot, Z, device, dtype,
            chunk_size=opts.get("chunk_size"),
            safety_factor=float(pot.get("empirical_LA_safety_factor", 2.0)),
            phi=str(pot.get("phi", "log_cosh")),
            feature_multiplier=int(pot.get("feature_multiplier", 4)))
    except Exception:
        # Never lose a row to a centering fast-path issue: fall back to CPU.
        return make_potential(cfg, point, Z)


def chunk_size(cfg):
    return cfg["monte_carlo"].get("chunk_size", None)


def eigsh_opts(cfg):
    e = cfg.get("eigsh", {})
    return dict(eigsh_tol=float(e.get("tol", 1e-6)),
                eigsh_maxiter=int(e.get("maxiter", 1000)))


def stage_dir(cfg, args, stage):
    """Resolve the output directory for a stage (``--outdir`` overrides config)."""
    base = args.outdir if getattr(args, "outdir", None) else cfg["outputs"]["base_dir"]
    out = os.path.join(base, stage)
    ensure_dir(out)
    return out


def potential_meta_row(pot, point):
    """Common metadata columns for a results row."""
    md = pot.metadata()
    keep = ["alpha_target", "beta_target", "rho", "L_A", "L_A_is_empirical",
            "norm_mean_grad", "norm_mean_hess_minus_I",
            "empirical_min_hess_eig", "empirical_max_hess_eig"]
    row = {
        "family": point["family"],
        "N_theta": point["N_theta"],
        "kappa_target": point["kappa_target"],
        "seed": point["seed"],
    }
    for k in keep:
        row[k] = md.get(k)
    return row
