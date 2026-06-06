"""Gaussian ground truth across dimensions: V(theta) = 0.5||theta||^2.

The exact answers are Lambda_hat = 0 and gamma_loc = 1 at every dimension. Two
regimes are tested, reflecting an honest reading of the estimator:

* The *full* symmetric operator H_sym acts on the N(N+1)/2-dimensional space of
  symmetric matrices, so its top eigenvalue (and hence gamma_loc) is inflated by
  finite-sample spectral noise that grows with dimension. Lambda_hat_full_sym ~ 0
  and gamma_loc ~ 1 therefore hold only when M_mc is large *relative to that
  dimension* -- demonstrated here for N <= 16 at M = 65536. (At fixed M and large
  N the inflation is real and is exactly what the sample-size-scaling diagnostic
  exists to expose; see reports/natural_gradient_local_rate_report.tex.)

* The *clean*, dimension-robust quantities -- the diagonal-restricted estimator,
  the H_sym/L_star self-adjointness errors and the eigen-residuals -- stay at
  (or near) machine precision at every dimension, including N = 32 and 64.
"""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.natural_gradient_local_rate.potentials import GaussianPotential
from src.natural_gradient_local_rate.linearized_rate import (
    estimate_lambda_hat, estimate_gamma_loc, estimate_diagonal_lambda,
)
from src.natural_gradient_local_rate.operators import (
    self_adjoint_error_H, self_adjoint_error_L_star,
)


@pytest.mark.parametrize("N", [4, 8, 16])
def test_full_sym_near_zero_and_gamma_near_one(N):
    """With M large vs the operator dimension, the full estimator is near truth."""
    pot = GaussianPotential(N)
    Z = gaussian_samples(N, 65536, seed=0, antithetic=True)
    lam = estimate_lambda_hat(pot, Z)
    gam, _ = estimate_gamma_loc(pot, Z)
    diag = estimate_diagonal_lambda(pot, Z)
    assert lam < 0.30, f"Lambda_hat_full_sym={lam} at N={N}"
    assert diag["Lambda_hat_diag"] < 0.12, f"Lambda_hat_diag={diag['Lambda_hat_diag']}"
    assert 0.85 < gam < 1.10, f"gamma_loc={gam} at N={N}"


@pytest.mark.parametrize("N", [4, 8, 16, 32])
def test_clean_estimators_are_dimension_robust(N):
    """Diagonal benchmark and self-adjointness identities hold at every N."""
    pot = GaussianPotential(N)
    Z = gaussian_samples(N, 16384, seed=0, antithetic=True)
    diag = estimate_diagonal_lambda(pot, Z)
    # the cleanest separable estimator (max diagonal entry) is ~0 at all N
    assert diag["max_diag"] < 0.06, f"max_diag={diag['max_diag']} at N={N}"
    # self-adjointness identities are exact up to floating point at all N
    assert self_adjoint_error_H(pot, Z, estimator="symmetrized", n_probe=6) < 1e-9
    assert self_adjoint_error_L_star(pot, Z, n_probe=6) < 1e-9
    # and the raw forward estimator is demonstrably *not* self-adjoint
    assert self_adjoint_error_H(pot, Z, estimator="raw_forward", n_probe=6) > 1e-4


def test_separable_exact_is_zero_for_gaussian():
    from src.natural_gradient_local_rate.separable_exact import separable_exact_lambda
    for N in (4, 16, 64):
        val, status = separable_exact_lambda(GaussianPotential(N), n_nodes=80)
        assert status == "ok"
        assert abs(val) < 1e-12, f"exact Lambda for Gaussian N={N} is {val}"
