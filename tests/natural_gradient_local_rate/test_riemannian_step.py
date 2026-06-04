"""First-order consistency of the Riemannian exponential-map step.

For small ``Delta_t``:
    (m_next - m) / Delta_t  == C g            (exact; same g as the field)
    (C_next - C) / Delta_t  ~= C + C H C      (first order in Delta_t)
"""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.common.symspace import random_symmetric
from src.natural_gradient_local_rate.potentials import build_potential, GaussianPotential
from src.natural_gradient_local_rate.riemannian_flow import (
    natural_gradient_vector_field, riemannian_step,
)


def _state(N, rng, scale=0.1):
    m = scale * rng.standard_normal(N)
    C = np.eye(N) + scale * random_symmetric(N, rng, normalize=True)
    return m, C


@pytest.mark.parametrize("family", ["gaussian", "random_feature"])
def test_mean_step_is_exact(family):
    N = 4
    Z = gaussian_samples(N, 8192, seed=0, antithetic=True)
    pot = GaussianPotential(N) if family == "gaussian" else \
        build_potential(family, N, 5.0, seed=0, Z_ref=Z)
    rng = np.random.default_rng(1)
    m, C = _state(N, rng)
    dm, dC, g, H = natural_gradient_vector_field(pot, m, C, Z)
    dt = 1e-3
    m_next, C_next, diag = riemannian_step(pot, m, C, Z, dt)
    # mean update is linear in dt with slope C g == dm
    assert np.allclose((m_next - m) / dt, dm, atol=1e-9)
    assert diag["spd_ok"]


@pytest.mark.parametrize("family", ["gaussian", "random_feature"])
def test_cov_step_first_order(family):
    N = 4
    Z = gaussian_samples(N, 8192, seed=0, antithetic=True)
    pot = GaussianPotential(N) if family == "gaussian" else \
        build_potential(family, N, 5.0, seed=0, Z_ref=Z)
    rng = np.random.default_rng(2)
    m, C = _state(N, rng)
    dm, dC, g, H = natural_gradient_vector_field(pot, m, C, Z)

    def cov_err(dt):
        _, C_next, _ = riemannian_step(pot, m, C, Z, dt)
        return np.linalg.norm((C_next - C) / dt - dC)

    e1, e2 = cov_err(1e-3), cov_err(5e-4)
    # consistent (error -> 0) and first order (halving dt ~ halves error)
    assert e2 < e1
    assert e1 < 5e-2
    assert e2 / e1 == pytest.approx(0.5, abs=0.2)


def test_step_preserves_spd_and_symmetry():
    N = 5
    Z = gaussian_samples(N, 4096, seed=0, antithetic=True)
    pot = build_potential("separable", N, 10.0, seed=0, Z_ref=Z)
    rng = np.random.default_rng(3)
    m, C = _state(N, rng, scale=0.2)
    for _ in range(20):
        m, C, diag = riemannian_step(pot, m, C, Z, 0.05)
        assert diag["spd_ok"]
        assert np.allclose(C, C.T, atol=1e-12)
