"""Adjoint / self-adjointness identities for the local operators.

On a fixed sample bank these identities are exact (up to floating point), which
guards against indexing or scaling mistakes in the operator code.
"""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.common.symspace import sym_inner, random_symmetric
from src.natural_gradient_local_rate.potentials import build_potential
from src.natural_gradient_local_rate.operators import (
    apply_T, apply_T_star, apply_H_sym, apply_H_lin, apply_H_lin_adjoint,
)


@pytest.mark.parametrize("family", ["separable", "random_feature", "radial_tail"])
def test_T_Tstar_adjoint(family):
    N = 6
    pot = build_potential(family, N, kappa_target=8.0, seed=0, centering_samples=4096)
    Z = gaussian_samples(N, 4096, seed=3, antithetic=True)
    H = pot.batch_hess(Z)
    rng = np.random.default_rng(0)
    for _ in range(6):
        u = rng.standard_normal(N)
        X = random_symmetric(N, rng)
        lhs = float(np.dot(apply_T(X, H, Z), u))   # <T[X], u>
        rhs = sym_inner(X, apply_T_star(u, H, Z))   # <X, T_star[u]>
        assert lhs == pytest.approx(rhs, rel=1e-9, abs=1e-9)


@pytest.mark.parametrize("family", ["separable", "random_feature"])
def test_H_sym_self_adjoint(family):
    N = 6
    pot = build_potential(family, N, kappa_target=8.0, seed=1, centering_samples=4096)
    Z = gaussian_samples(N, 4096, seed=5, antithetic=True)
    H = pot.batch_hess(Z)
    rng = np.random.default_rng(1)
    for _ in range(6):
        X = random_symmetric(N, rng)
        Y = random_symmetric(N, rng)
        lhs = sym_inner(apply_H_sym(X, H, Z), Y)
        rhs = sym_inner(X, apply_H_sym(Y, H, Z))
        assert lhs == pytest.approx(rhs, rel=1e-9, abs=1e-9)


def test_H_lin_adjoint_is_frobenius_adjoint():
    # <H_lin[X], Y> == <X, H_lin*[Y]> exactly on a fixed bank.
    N = 5
    pot = build_potential("random_feature", N, kappa_target=5.0, seed=2, centering_samples=4096)
    Z = gaussian_samples(N, 4096, seed=7, antithetic=True)
    H = pot.batch_hess(Z)
    rng = np.random.default_rng(2)
    X = random_symmetric(N, rng)
    Y = random_symmetric(N, rng)
    lhs = sym_inner(apply_H_lin(X, H, Z), Y)
    rhs = sym_inner(X, apply_H_lin_adjoint(Y, H, Z))
    assert lhs == pytest.approx(rhs, rel=1e-9, abs=1e-9)


def test_outputs_are_symmetric_matrices():
    N = 5
    pot = build_potential("random_feature", N, kappa_target=5.0, seed=0, centering_samples=2048)
    Z = gaussian_samples(N, 2048, seed=1, antithetic=True)
    H = pot.batch_hess(Z)
    rng = np.random.default_rng(0)
    X = random_symmetric(N, rng)
    u = rng.standard_normal(N)
    for S in [apply_T_star(u, H, Z), apply_H_lin(X, H, Z), apply_H_sym(X, H, Z)]:
        assert np.allclose(S, S.T, atol=1e-12)
