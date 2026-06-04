"""Finite-sample self-adjointness of the symmetrized operator H_sym.

On a fixed sample bank the forward estimator H_lin is generally *not* Frobenius
self-adjoint, while H_sym = 0.5 (H_lin + H_lin*) is self-adjoint by construction.
The symmetrized self-adjointness error should sit at machine precision and be
orders of magnitude smaller than the raw forward error.
"""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.common.symspace import sym_inner, random_symmetric
from src.natural_gradient_local_rate.potentials import build_potential
from src.natural_gradient_local_rate.operators import (
    apply_H_lin, apply_H_sym, self_adjoint_error_H,
)


@pytest.mark.parametrize("family", ["separable", "random_feature", "radial_tail"])
def test_H_sym_is_self_adjoint_raw_is_not(family):
    N = 6
    pot = build_potential(family, N, kappa_target=10.0, seed=0, centering_samples=4096)
    Z = gaussian_samples(N, 8192, seed=1, antithetic=True)
    err_sym = self_adjoint_error_H(pot, Z, estimator="symmetrized", n_probe=8)
    err_raw = self_adjoint_error_H(pot, Z, estimator="raw_forward", n_probe=8)
    assert err_sym < 1e-9, f"H_sym self-adjointness error {err_sym}"
    assert err_raw > 1e-4, f"raw forward error unexpectedly tiny: {err_raw}"
    # symmetrization improves self-adjointness by many orders of magnitude
    assert err_sym < err_raw / 1e3


def test_H_sym_inner_product_symmetry_directly():
    # <X, H_sym[Y]> == <H_sym[X], Y> for explicit random symmetric pairs.
    N = 5
    pot = build_potential("random_feature", N, kappa_target=8.0, seed=2,
                          centering_samples=4096)
    Z = gaussian_samples(N, 4096, seed=3, antithetic=True)
    H = pot.batch_hess(Z)
    rng = np.random.default_rng(0)
    for _ in range(6):
        X = random_symmetric(N, rng)
        Y = random_symmetric(N, rng)
        lhs = sym_inner(X, apply_H_sym(Y, H, Z))
        rhs = sym_inner(apply_H_sym(X, H, Z), Y)
        assert lhs == pytest.approx(rhs, rel=1e-9, abs=1e-9)


def test_H_sym_output_is_symmetric_matrix():
    N = 5
    pot = build_potential("separable", N, kappa_target=5.0, seed=0, centering_samples=4096)
    Z = gaussian_samples(N, 4096, seed=1, antithetic=True)
    H = pot.batch_hess(Z)
    rng = np.random.default_rng(1)
    X = random_symmetric(N, rng)
    S = apply_H_sym(X, H, Z)
    assert np.allclose(S, S.T, atol=1e-12)
    # the raw forward output is also a symmetric *matrix* even though the
    # operator is not self-adjoint
    R = apply_H_lin(X, H, Z)
    assert np.allclose(R, R.T, atol=1e-12)
