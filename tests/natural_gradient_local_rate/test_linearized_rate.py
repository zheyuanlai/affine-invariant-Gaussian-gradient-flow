"""Tests for the linearized-rate estimators and the L_star operator."""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.natural_gradient_local_rate.potentials import build_potential, GaussianPotential
from src.natural_gradient_local_rate.operators import make_L_star_operator
from src.natural_gradient_local_rate.linearized_rate import (
    estimate_lambda_hat, estimate_gamma_loc, _to_dense,
)
from src.natural_gradient_local_rate import diagnostics


@pytest.mark.parametrize("N", [4, 9])
def test_L_star_half_operator_is_symmetric(N):
    # The Fisher--Rao-whitened L_star must be a symmetric matrix (eigsh validity).
    pot = build_potential("random_feature", N, kappa_target=8.0, seed=0,
                          centering_samples=4096)
    Z = gaussian_samples(N, 4096, seed=1, antithetic=True)
    op = make_L_star_operator(pot, Z, covariance_weight="half")
    Mat = _to_dense(op)
    asym = np.max(np.abs(Mat - Mat.T))
    assert asym < 1e-9, f"asymmetry {asym}"


def test_plain_operator_not_required_symmetric():
    # 'plain' coordinates are generally not symmetric; this documents the choice.
    N = 5
    pot = build_potential("random_feature", N, kappa_target=20.0, seed=0,
                          centering_samples=4096)
    Z = gaussian_samples(N, 4096, seed=1, antithetic=True)
    op_plain = make_L_star_operator(pot, Z, covariance_weight="plain")
    dim = op_plain.shape[0]
    Mat = np.column_stack([op_plain.matvec(np.eye(dim)[:, i]) for i in range(dim)])
    assert np.max(np.abs(Mat - Mat.T)) > 1e-6  # asymmetric, as expected


@pytest.mark.parametrize("family", ["separable", "random_feature"])
def test_centered_rate_is_positive_and_bounded(family):
    N = 6
    pot = build_potential(family, N, kappa_target=10.0, seed=0, centering_samples=8192)
    Z = gaussian_samples(N, 16384, seed=2, antithetic=True)
    lam = estimate_lambda_hat(pot, Z)
    gam, eig = estimate_gamma_loc(pot, Z)
    assert lam > 0.0
    assert 0.0 < gam <= 1.05            # PD generator, rate <~ 1
    u_star, X_star = eig
    assert u_star.shape == (N,)
    assert X_star.shape == (N, N)
    assert np.allclose(X_star, X_star.T, atol=1e-10)


def test_eigsh_path_matches_dense_path():
    # N=9 uses eigsh; compare against an explicit dense diagonalization.
    N = 9
    pot = build_potential("random_feature", N, kappa_target=10.0, seed=0,
                          centering_samples=8192)
    Z = gaussian_samples(N, 8192, seed=3, antithetic=True)
    gam_eigsh, _ = estimate_gamma_loc(pot, Z, eigsh_tol=1e-9)
    op = make_L_star_operator(pot, Z, covariance_weight="half")
    w = np.linalg.eigvalsh(_to_dense(op))
    assert gam_eigsh == pytest.approx(float(w[0]), rel=1e-5)


def test_reference_columns():
    cols = diagnostics.reference_columns(N_theta=16, kappa=10.0, beta_target=1.8,
                                         Lambda_hat=2.0, gamma_loc=0.25)
    lk = 1.0 + np.log(10.0)
    assert cols["log_kappa_factor"] == pytest.approx(lk)
    assert cols["lambda_over_logkappa"] == pytest.approx(2.0 / lk)
    assert cols["inverse_gamma_loc"] == pytest.approx(4.0)
    assert cols["inverse_gamma_over_logkappa"] == pytest.approx(4.0 / lk)
    # current bound uses min(beta-1, ...) = min(0.8, ...) = 0.8 here
    assert cols["current_bound_rate"] == pytest.approx(1.0 / (4.0 + 0.8))
    assert 0.0 < cols["conjecture_bound_rate"] < 1.0


def test_lambda_hat_eigenvector_shape():
    N = 5
    pot = build_potential("separable", N, kappa_target=5.0, seed=0, centering_samples=4096)
    Z = gaussian_samples(N, 8192, seed=0, antithetic=True)
    lam, X_eig = estimate_lambda_hat(pot, Z, return_eigenvector=True)
    assert X_eig.shape == (N, N)
    assert np.allclose(X_eig, X_eig.T, atol=1e-10)
