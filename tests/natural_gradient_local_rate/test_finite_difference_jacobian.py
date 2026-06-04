"""Finite-difference Jacobian of the flow vs the analytic generator.

The flow vector field ``F(m, C) = (C g, C + C H C)`` has Jacobian
``DF(a_star) = J_star = -L_star`` at the equilibrium ``a_star = (0, I)``.

The flow uses a reparameterization (pathwise) estimator while ``apply_L_star``
uses the Stein-form operators; these are different finite-sample estimators of
the same limit, so the match is exact only as ``M -> infinity``. We therefore use
a large bank and an MC-scaled tolerance, and check the sign explicitly.
"""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.common.symspace import pack_tangent, random_symmetric
from src.natural_gradient_local_rate.potentials import build_potential, GaussianPotential
from src.natural_gradient_local_rate.riemannian_flow import natural_gradient_vector_field
from src.natural_gradient_local_rate.operators import apply_L_star


def _fd_jacobian(pot, Z, v, Y, s=1e-5):
    """Finite-difference of F along (v, Y) at a_star = (0, I): ~ J_star (v, Y)."""
    N = v.size
    m = s * v
    C = np.eye(N) + s * Y
    dm, dC, _, _ = natural_gradient_vector_field(pot, m, C, Z)  # F(a_star + s xi) (F(a_star)=0)
    return dm / s, dC / s


def test_sign_convention_gaussian():
    # Gaussian: J_star = -I exactly on the mean block, so fd_dm ~ -v.
    N = 4
    Z = gaussian_samples(N, 16384, seed=0, antithetic=True)
    pot = GaussianPotential(N)
    rng = np.random.default_rng(0)
    v = rng.standard_normal(N)
    Y = random_symmetric(N, rng)
    fd_dm, fd_dC = _fd_jacobian(pot, Z, v, Y)
    # L_star(v, Y) mean-block is exactly v (T[Y] = 0), and J_star = -L_star.
    assert np.allclose(fd_dm, -v, atol=1e-3)


def test_fd_jacobian_matches_minus_L_star():
    N = 4
    M = 65536
    Z = gaussian_samples(N, M, seed=0, antithetic=True)
    pot = build_potential("random_feature", N, 8.0, seed=0, Z_ref=Z)
    rng = np.random.default_rng(1)
    rel_errs, cosines = [], []
    for _ in range(4):
        v = rng.standard_normal(N)
        Y = random_symmetric(N, rng)
        fd_dm, fd_dC = _fd_jacobian(pot, Z, v, Y)
        uL, XL = apply_L_star(pot, v, Y, Z)         # L_star(v, Y)
        fd = pack_tangent(fd_dm, fd_dC)
        analytic = pack_tangent(-uL, -XL)           # J_star = -L_star
        rel_errs.append(np.linalg.norm(fd - analytic) / np.linalg.norm(analytic))
        cosines.append(float(fd @ analytic / (np.linalg.norm(fd) * np.linalg.norm(analytic))))
    assert np.mean(rel_errs) < 0.03, f"rel errors {rel_errs}"
    assert min(cosines) > 0.99, f"cosines {cosines}"
