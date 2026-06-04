"""One-stop assembly of all natural-gradient local-rate estimator columns.

:func:`compute_row` runs, for a single ``(potential, sample bank)`` pair, every
estimator and diagnostic the pipeline reports and returns a flat dict ready to
become one CSV row. It is the single entry point shared by the operator-grid,
linearized-rate-grid and sample-size-scaling runner scripts (and by the tests),
so all stages emit an identical, comparable schema.

Estimator modes (Part 1 of the spec):

* ``raw_forward``   -- uncorrected forward ``H_lin`` (diagnostic / backward
  comparison only; not self-adjoint on a finite bank);
* ``symmetrized``   -- self-adjoint ``H_sym`` (the default; used for every
  eigenvalue computation, including the ``L_star`` covariance block);
* ``diagonal_restricted`` -- ``lambda_max`` of the diagonal-mode matrix
  ``A = G - 1 1^T`` (separable sanity check);
* ``separable_exact`` -- Gauss--Hermite quadrature ground truth for separable
  control potentials.

``compute_row`` always records the raw / symmetrized / diagonal / exact numbers
side by side; ``operator_estimator`` names which one is the *headline*
``Lambda_hat`` (default ``symmetrized``). The ``L_star`` covariance block always
uses ``H_sym``, never the raw forward estimator.

Backends: ``backend="numpy"`` (default) runs the NumPy/SciPy CPU path below;
``backend="torch"`` dispatches to the PyTorch GPU backend
(:mod:`src.natural_gradient_local_rate.torch_backend`); ``backend="auto"`` uses
torch only when it is importable *and* the requested device resolves to CUDA,
otherwise NumPy. The torch path emits the same CSV schema plus GPU bookkeeping
columns. torch is imported lazily so the NumPy-only environment is unaffected.
"""
from __future__ import annotations

import time

import numpy as np

from src.common.symspace import sym_norm, pack_tangent_fr
from src.natural_gradient_local_rate import diagnostics
from src.natural_gradient_local_rate.operators import (
    apply_H_sym_banked, make_L_star_operator,
    self_adjoint_error_H, self_adjoint_error_L_star,
)
from src.natural_gradient_local_rate.linearized_rate import (
    estimate_lambda_hat, estimate_gamma_loc, estimate_diagonal_lambda,
    RAW_DENSE_LIMIT,
)
from src.natural_gradient_local_rate.separable_exact import separable_exact_lambda

EXPERIMENT_GROUP = "natural_gradient_local_rate"


def default_opts():
    """Default estimator options (mirrors the ``operator:`` config block)."""
    return dict(
        estimator="symmetrized",
        diagonal_benchmark=True,
        separable_exact_benchmark=True,
        compute_gamma_loc=True,
        backend="numpy",
        chunk_size=None,
        eigsh_tol=1e-6,
        eigsh_maxiter=1000,
        quadrature_nodes=80,
        self_adjoint_probes=8,
        self_adjoint_seed=0,
        raw_dense_limit=RAW_DENSE_LIMIT,
        # --- torch / GPU backend options (ignored by the numpy path) ---
        device="auto",
        dtype="float64",
        eigensolver="auto",
        explicit_dense_max_N_theta=64,
        basis_block_size=32,
    )


def _safe_ratio(a, b):
    if b is None or not np.isfinite(b) or b == 0.0:
        return float("nan")
    return float(a) / float(b)


def _safe_diff(a, b):
    if a is None or b is None or not (np.isfinite(a) and np.isfinite(b)):
        return float("nan")
    return float(a) - float(b)


def compute_row(potential, Z, point, opts=None, *, run_id="",
                experiment_group=EXPERIMENT_GROUP):
    """Compute every estimator/diagnostic column for one grid point.

    ``Z`` is the fixed sample bank the ``potential`` was centered on. ``point``
    is a dict with ``family``, ``N_theta``, ``kappa_target``, ``seed`` and
    ``M_mc``; ``opts`` mirrors :func:`default_opts`. Returns a flat dict
    including ``status`` / ``error_message`` / ``runtime_seconds``; on an
    unexpected failure the numeric columns are ``NaN`` and the run continues.
    """
    o = default_opts()
    if opts:
        o.update(opts)

    # --- backend dispatch (torch is imported lazily; numpy needs no torch) ---
    backend = str(o["backend"]).lower()
    if backend not in ("numpy", "torch", "auto"):
        raise ValueError(f"unknown backend {backend!r} (expected numpy/torch/auto)")
    if backend != "numpy":
        from src.common.torch_utils import resolve_backend
        from src.natural_gradient_local_rate import torch_backend
        resolved = resolve_backend(backend, o.get("device", "auto"))
        if resolved == "torch" and not torch_backend.torch_supports_family(point["family"]):
            # auto -> fall back to numpy for unsupported families; explicit torch
            # is left to fail clearly inside compute_row_torch.
            if backend == "auto":
                resolved = "numpy"
        if resolved == "torch":
            return torch_backend.compute_row_torch(
                potential, Z, point, o, run_id=run_id,
                experiment_group=experiment_group)
        # fell back to numpy (auto without CUDA, or unsupported family): record
        # the backend actually used so the CSV column is accurate.
        o["backend"] = "numpy"
        # else: fall through to the numpy path below

    md = potential.metadata()
    N = int(point["N_theta"])
    kappa = float(point["kappa_target"])
    beta = float(md.get("beta_target", 1.0))
    cs = o["chunk_size"]

    row = {
        "run_id": run_id,
        "experiment_group": experiment_group,
        "potential_family": point["family"],
        "family": point["family"],          # legacy alias for existing plot scripts
        "seed": int(point["seed"]),
        "N_theta": N,
        "kappa_target": kappa,
        "rho": float(md.get("rho", 0.0)),
        "M_mc": int(point["M_mc"]),
        "operator_estimator": o["estimator"],
        "backend": o["backend"],
        "chunk_size": cs if cs is not None else "",
        "quadrature_nodes": int(o["quadrature_nodes"]),
        # numeric columns default to NaN until computed
        "Lambda_hat_raw_forward": float("nan"),
        "Lambda_hat_full_sym": float("nan"),
        "Lambda_hat_diag": float("nan"),
        "Lambda_hat_separable_exact": float("nan"),
        "diag_offdiag_norm": float("nan"),
        "full_over_diag": float("nan"),
        "full_minus_diag": float("nan"),
        "diag_minus_exact": float("nan"),
        "full_sym_minus_exact": float("nan"),
        "gamma_loc": float("nan"),
        "inverse_gamma_loc": float("nan"),
        "self_adjoint_error_H_raw": float("nan"),
        "self_adjoint_error_H_sym": float("nan"),
        "self_adjoint_error_L_star": float("nan"),
        "eig_residual_H": float("nan"),
        "eig_residual_L_star": float("nan"),
        "separable_exact_status": "",
        "status": "ok",
        "error_message": "",
    }
    row.update(diagnostics.true_benchmark_columns(point["family"]))

    t0 = time.time()
    try:
        Zbank = Z

        # --- symmetrized full operator (the default / headline Lambda_hat) ---
        lam_sym, X_eig = estimate_lambda_hat(
            potential, Zbank, estimator="symmetrized", chunk_size=cs,
            eigsh_tol=o["eigsh_tol"], eigsh_maxiter=o["eigsh_maxiter"],
            return_eigenvector=True)
        row["Lambda_hat_full_sym"] = float(lam_sym)
        HX = apply_H_sym_banked(potential, X_eig, Zbank, chunk_size=cs)
        row["eig_residual_H"] = sym_norm(HX - lam_sym * X_eig) / max(1.0, abs(lam_sym))

        # --- raw forward operator (diagnostic / backward comparison) ---
        row["Lambda_hat_raw_forward"] = float(estimate_lambda_hat(
            potential, Zbank, estimator="raw_forward", chunk_size=cs,
            eigsh_tol=o["eigsh_tol"], eigsh_maxiter=o["eigsh_maxiter"],
            raw_dense_limit=o["raw_dense_limit"]))

        # --- self-adjointness errors (the clean algebraic identities) ---
        row["self_adjoint_error_H_raw"] = self_adjoint_error_H(
            potential, Zbank, estimator="raw_forward",
            n_probe=o["self_adjoint_probes"], seed=o["self_adjoint_seed"], chunk_size=cs)
        row["self_adjoint_error_H_sym"] = self_adjoint_error_H(
            potential, Zbank, estimator="symmetrized",
            n_probe=o["self_adjoint_probes"], seed=o["self_adjoint_seed"], chunk_size=cs)

        # --- diagonal-restricted benchmark ---
        if o["diagonal_benchmark"]:
            diag = estimate_diagonal_lambda(potential, Zbank, chunk_size=cs)
            row["Lambda_hat_diag"] = diag["Lambda_hat_diag"]
            row["diag_offdiag_norm"] = diag["diag_offdiag_norm"]

        # --- separable exact (Gauss--Hermite) benchmark ---
        if o["separable_exact_benchmark"]:
            lam_exact, exact_status = separable_exact_lambda(
                potential, n_nodes=o["quadrature_nodes"])
            row["Lambda_hat_separable_exact"] = lam_exact
            row["separable_exact_status"] = exact_status
            if row["baseline_type"] == "" and exact_status == "ok":
                row["baseline_type"] = "separable_exact"

        # --- derived comparisons ---
        row["full_over_diag"] = _safe_ratio(row["Lambda_hat_full_sym"], row["Lambda_hat_diag"])
        row["full_minus_diag"] = _safe_diff(row["Lambda_hat_full_sym"], row["Lambda_hat_diag"])
        row["diag_minus_exact"] = _safe_diff(row["Lambda_hat_diag"], row["Lambda_hat_separable_exact"])
        row["full_sym_minus_exact"] = _safe_diff(row["Lambda_hat_full_sym"], row["Lambda_hat_separable_exact"])

        # --- local rate gamma_loc (always via H_sym + Fisher--Rao metric) ---
        if o["compute_gamma_loc"]:
            gamma, (u_star, X_star) = estimate_gamma_loc(
                potential, Zbank, chunk_size=cs, eigsh_tol=o["eigsh_tol"],
                eigsh_maxiter=o["eigsh_maxiter"], return_eigenvector=True)
            row["gamma_loc"] = float(gamma)
            row["inverse_gamma_loc"] = (1.0 / gamma) if gamma != 0.0 else float("inf")
            row["self_adjoint_error_L_star"] = self_adjoint_error_L_star(
                potential, Zbank, n_probe=o["self_adjoint_probes"],
                seed=o["self_adjoint_seed"], chunk_size=cs)
            op = make_L_star_operator(potential, Zbank, covariance_weight="half",
                                      chunk_size=cs)
            y = pack_tangent_fr(u_star, X_star)
            resid = np.linalg.norm(op.matvec(y) - gamma * y) / max(1.0, abs(gamma))
            row["eig_residual_L_star"] = float(resid)
            row["_u_star"] = u_star
            row["_X_star"] = X_star

        # --- reference / headline bookkeeping ---
        headline = (row["Lambda_hat_raw_forward"] if o["estimator"] == "raw_forward"
                    else row["Lambda_hat_full_sym"])
        row["Lambda_hat"] = float(headline)   # legacy alias for existing plot scripts
        ref = diagnostics.reference_columns(
            N, kappa, beta, Lambda_hat=headline,
            gamma_loc=(row["gamma_loc"] if np.isfinite(row["gamma_loc"]) else None))
        row.update(ref)
        # potential-construction bookkeeping (kept from the meta)
        for k in ("alpha_target", "beta_target", "L_A", "L_A_is_empirical",
                  "norm_mean_grad", "norm_mean_hess_minus_I",
                  "empirical_min_hess_eig", "empirical_max_hess_eig"):
            if k in md:
                row[k] = md[k]

    except Exception as exc:  # one bad point must not kill the grid
        row["status"] = "error"
        row["error_message"] = f"{type(exc).__name__}: {exc}"

    row["runtime_seconds"] = float(time.time() - t0)
    return row
