"""Analytic Gaussian ground truth: V(theta) = 0.5 ||theta||^2.

Then grad V = theta, Hess V = I, the centering conditions hold exactly, the
local operators vanish, Lambda_hat ~ 0 and gamma_loc ~ 1 (within MC error).
"""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.common.symspace import random_symmetric, sym_norm
from src.natural_gradient_local_rate.potentials import GaussianPotential
from src.natural_gradient_local_rate.operators import (
    apply_T, apply_T_star, apply_H_lin, apply_H_sym,
)
from src.natural_gradient_local_rate.linearized_rate import (
    estimate_lambda_hat, estimate_gamma_loc,
)


def test_grad_hess_exact():
    N = 5
    pot = GaussianPotential(N)
    rng = np.random.default_rng(0)
    theta = rng.standard_normal(N)
    assert np.allclose(pot.grad(theta), theta)
    assert np.allclose(pot.hess(theta), np.eye(N))


def test_local_operators_vanish():
    N = 5
    pot = GaussianPotential(N)
    Z = gaussian_samples(N, 16384, seed=0, antithetic=True)
    H = pot.batch_hess(Z)
    rng = np.random.default_rng(1)
    X = random_symmetric(N, rng)
    u = rng.standard_normal(N)
    # T[X] and T_star[u] are exactly zero (antithetic), H_lin[X] ~ 0 (MC).
    assert np.linalg.norm(apply_T(X, H, Z)) < 1e-9
    assert sym_norm(apply_T_star(u, H, Z)) < 1e-9
    assert sym_norm(apply_H_lin(X, H, Z)) < 0.15
    assert sym_norm(apply_H_sym(X, H, Z)) < 0.15


@pytest.mark.parametrize("N", [3, 5])
def test_lambda_hat_near_zero_gamma_near_one(N):
    pot = GaussianPotential(N)
    Z = gaussian_samples(N, 32768, seed=0, antithetic=True)
    lam = estimate_lambda_hat(pot, Z)
    gam, _ = estimate_gamma_loc(pot, Z)
    assert lam < 0.12, f"Lambda_hat={lam}"
    assert 0.88 < gam < 1.05, f"gamma_loc={gam}"
