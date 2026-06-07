"""Tests for the two discretization schemes.

Covers first-order consistency with the continuous vector field at small ``dt``,
SPD preservation, the unconditional SPD property of the KL step for a
log-concave target, and symmetry of every covariance update.
"""
import numpy as np
import pytest

from src.common.spd import symmetrize
from src.natural_gradient_discretization_stepsize.targets import (
    GaussianPosteriorTarget, SmoothLogconcaveTarget,
)
from src.natural_gradient_discretization_stepsize.methods import (
    riemannian_cov_step, kl_cov_step, mean_step, discretization_step,
)


def _vector_field(target, m, C):
    g, H = target.g_H(m, C)
    dm = C @ g
    dC = symmetrize(C + C @ H @ C)
    return dm, dC


@pytest.mark.parametrize("method", ["riemannian", "kl"])
@pytest.mark.parametrize("lam", [0.01, 0.1, 1.0])
def test_first_order_consistency(method, lam):
    """One small step matches the ODE vector field to first order in dt."""
    T = GaussianPosteriorTarget(lam)
    m = np.array([1.5, -0.8])
    C = np.array([[1.2, 0.25], [0.25, 0.9]])
    dm_true, dC_true = _vector_field(T, m, C)
    dt = 1e-6
    m1, C1, _ = discretization_step(method, T, m, C, dt)
    np.testing.assert_allclose((m1 - m) / dt, dm_true, rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose((C1 - C) / dt, dC_true, rtol=1e-4, atol=1e-6)


@pytest.mark.parametrize("method", ["riemannian", "kl"])
def test_consistency_smooth_target(method):
    """Consistency also holds for the non-separable smooth target."""
    T = SmoothLogconcaveTarget(0.1)
    m = np.array([0.3, 0.6])
    C = np.array([[1.4, -0.2], [-0.2, 1.1]])
    dm_true, dC_true = _vector_field(T, m, C)
    dt = 1e-6
    m1, C1, _ = discretization_step(method, T, m, C, dt)
    np.testing.assert_allclose((m1 - m) / dt, dm_true, rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose((C1 - C) / dt, dC_true, rtol=1e-4, atol=1e-6)


@pytest.mark.parametrize("method", ["riemannian", "kl"])
def test_spd_preserved_small_step(method):
    """A small stable step preserves SPD on the Gaussian target."""
    T = GaussianPosteriorTarget(0.1)
    m, C = T.m0.copy(), T.C0.copy()
    for _ in range(50):
        m, C, diag = discretization_step(method, T, m, C, 0.05)
        assert diag["spd_ok"]
        np.testing.assert_allclose(C, C.T, atol=1e-14)


@pytest.mark.parametrize("dt", [0.1, 1.0, 10.0, 100.0])
def test_kl_unconditionally_spd_for_logconcave(dt):
    """For H <= 0 (log-concave) the KL precision C^{-1} - dt H is SPD for any dt."""
    T = GaussianPosteriorTarget(0.1)        # H = -Q <= 0
    C = np.array([[0.5, 0.0], [0.0, 2.0]])
    _, H = T.g_H(T.m0, C)
    C_next = kl_cov_step(C, H, dt)
    eig = np.linalg.eigvalsh(C_next)
    assert eig[0] > 0.0
    np.testing.assert_allclose(C_next, C_next.T, atol=1e-14)


def test_riemannian_matches_scalar_form():
    """The 1x1 Riemannian step equals the scalar closed form C exp(dt(1-C))."""
    C = np.array([[3.0]])
    H = np.array([[-1.0]])      # scalar N(0,1): H = -1
    dt = 0.7
    out = riemannian_cov_step(C, H, dt)[0, 0]
    expected = 3.0 * np.exp(dt * (1.0 - 3.0))
    np.testing.assert_allclose(out, expected, atol=1e-12)


def test_kl_matches_scalar_form():
    """The 1x1 KL step equals the scalar closed form (1+dt)C/(1+dt C)."""
    C = np.array([[3.0]])
    H = np.array([[-1.0]])
    dt = 0.7
    out = kl_cov_step(C, H, dt)[0, 0]
    expected = (1.0 + dt) * 3.0 / (1.0 + dt * 3.0)
    np.testing.assert_allclose(out, expected, atol=1e-12)


def test_mean_step_shared():
    """Both schemes use the identical explicit mean step m + dt C g."""
    T = GaussianPosteriorTarget(0.1)
    m = np.array([2.0, 3.0])
    C = np.array([[1.1, 0.2], [0.2, 0.8]])
    dt = 0.3
    g, _ = T.g_H(m, C)
    expected = m + dt * (C @ g)
    mr, _, _ = discretization_step("riemannian", T, m, C, dt)
    mk, _, _ = discretization_step("kl", T, m, C, dt)
    np.testing.assert_allclose(mr, expected, atol=1e-13)
    np.testing.assert_allclose(mk, expected, atol=1e-13)
    np.testing.assert_allclose(mr, mk, atol=1e-13)
