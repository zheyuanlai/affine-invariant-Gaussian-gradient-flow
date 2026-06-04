"""Separable controls: Monte-Carlo diagonal estimate vs Gauss--Hermite quadrature.

For a separable potential the diagonal-mode coefficients A_ii are 1-D
expectations that Gauss--Hermite quadrature evaluates exactly. The Monte-Carlo
diagonal must converge to that quadrature benchmark as M_mc grows, and the
finite-sample off-diagonal leakage must shrink. The full-operator norm
(Lambda_hat_full_sym) is *not* the same quantity and must not be read as the
exact value without this comparison.
"""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.natural_gradient_local_rate.potentials import build_potential
from src.natural_gradient_local_rate.operators import diagonal_A_matrix
from src.natural_gradient_local_rate.linearized_rate import estimate_lambda_hat
from src.natural_gradient_local_rate.separable_exact import (
    separable_exact_diagonal, separable_exact_lambda, is_separable,
)

N = 6
KAPPA = 10.0


def _diag_error(M):
    Z = gaussian_samples(N, M, seed=0, antithetic=True)
    pot = build_potential("separable", N, kappa_target=KAPPA, seed=0, Z_ref=Z)
    A = diagonal_A_matrix(pot, Z)
    exact = separable_exact_diagonal(pot, n_nodes=80)
    err = float(np.max(np.abs(np.diag(A) - exact)))
    offdiag = float(np.linalg.norm(A - np.diag(np.diag(A))))
    return err, offdiag


def test_mc_diagonal_converges_to_quadrature():
    err_small, off_small = _diag_error(2048)
    err_large, off_large = _diag_error(65536)
    assert err_large < err_small, f"diag error did not shrink: {err_small} -> {err_large}"
    assert err_large < 0.02, f"diag error at large M too big: {err_large}"
    # off-diagonal leakage (zero in expectation for separable) also shrinks
    assert off_large < off_small


def test_quadrature_benchmark_is_finite_and_stable():
    Z = gaussian_samples(N, 4096, seed=0, antithetic=True)
    pot = build_potential("separable", N, kappa_target=KAPPA, seed=0, Z_ref=Z)
    assert is_separable(pot)
    lam40, st40 = separable_exact_lambda(pot, n_nodes=40)
    lam80, st80 = separable_exact_lambda(pot, n_nodes=80)
    assert st40 == "ok" and st80 == "ok"
    assert np.isfinite(lam80)
    assert -2.0 < lam80 < 6.0          # within the manuscript's tail-cut envelope
    assert abs(lam40 - lam80) < 1e-2   # quadrature has converged in node count


def test_full_sym_differs_from_exact():
    # The full-operator norm acts on all symmetric perturbations and is a
    # different (larger) quantity than the diagonal-mode exact value; this
    # documents why full_sym must be compared, not interpreted alone.
    Z = gaussian_samples(N, 16384, seed=0, antithetic=True)
    pot = build_potential("separable", N, kappa_target=KAPPA, seed=0, Z_ref=Z)
    full_sym = estimate_lambda_hat(pot, Z, estimator="symmetrized")
    exact, _ = separable_exact_lambda(pot, n_nodes=80)
    assert full_sym - exact > 0.1


def test_non_separable_quadrature_is_nan():
    Z = gaussian_samples(N, 4096, seed=0, antithetic=True)
    pot = build_potential("random_feature", N, kappa_target=KAPPA, seed=0, Z_ref=Z)
    assert not is_separable(pot)
    val, status = separable_exact_lambda(pot, n_nodes=80)
    assert np.isnan(val)
    assert status == "not_separable"
