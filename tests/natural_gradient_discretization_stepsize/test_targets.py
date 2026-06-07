"""Tests for the discretization-stepsize targets.

Covers the analytic Gaussian-posterior formulas, the literature quartic target
(analytic expected gradient/Hessian/energy vs Monte Carlo), and the smooth
log-concave Gauss--Hermite quadrature (sech^2 expectation bounds and Hessian
spectrum bounds).
"""
import numpy as np
import pytest

from src.natural_gradient_discretization_stepsize.targets import (
    GaussianPosteriorTarget, LiteratureLogconcaveTarget, SmoothLogconcaveTarget,
    gaussian_expectation_1d, _sech2,
)


# ---------------------------------------------------------------------------
# Target A: Gaussian posterior
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("lam", [0.01, 0.1, 1.0])
def test_gaussian_grad_hess_exact(lam):
    """g = -Q m and H = -Q exactly, for any state."""
    T = GaussianPosteriorTarget(lam)
    rng = np.random.default_rng(0)
    m = rng.standard_normal(2)
    C = np.array([[1.3, 0.2], [0.2, 0.7]])
    g, H = T.g_H(m, C)
    np.testing.assert_allclose(g, -(T.Q @ m), atol=1e-14)
    np.testing.assert_allclose(H, -T.Q, atol=1e-14)


@pytest.mark.parametrize("lam", [0.01, 0.1, 1.0])
def test_gaussian_energy_gap_zero_at_star(lam):
    """Energy gap vanishes at m=0, C=Q^{-1}; matches F - F_star elsewhere."""
    T = GaussianPosteriorTarget(lam)
    m_star, C_star, F_star = T.star()
    assert abs(T.energy_gap(m_star, C_star)) < 1e-12
    m = np.array([2.0, -1.0])
    C = np.array([[0.8, 0.1], [0.1, 1.5]])
    np.testing.assert_allclose(T.energy_gap(m, C), T.objective(m, C) - F_star, atol=1e-10)


def test_gaussian_energy_gap_nonnegative():
    """The analytic energy gap is non-negative (KL >= 0)."""
    T = GaussianPosteriorTarget(0.1)
    rng = np.random.default_rng(1)
    for _ in range(20):
        m = rng.standard_normal(2) * 3
        A = rng.standard_normal((2, 2))
        C = A @ A.T + 0.1 * np.eye(2)
        assert T.energy_gap(m, C) >= -1e-10


# ---------------------------------------------------------------------------
# Target B: literature quartic (analytic vs Monte Carlo)
# ---------------------------------------------------------------------------

def _mc_grad_hess_V(lam, m, C, n=4_000_000, seed=0):
    """Monte Carlo E[grad V], E[Hess V] for the literature quartic potential."""
    sl = np.sqrt(lam)
    rng = np.random.default_rng(seed)
    L = np.linalg.cholesky(C)
    Z = rng.standard_normal((n, 2))
    X = m + Z @ L.T
    x, y = X[:, 0], X[:, 1]
    gx = (lam * x - sl * y) / 10.0
    gy = (-sl * x + y) / 10.0 + y ** 3 / 5.0
    grad = np.array([gx.mean(), gy.mean()])
    h22 = (1.0 / 10.0 + 3.0 * y ** 2 / 5.0).mean()
    H = np.array([[lam / 10.0, -sl / 10.0], [-sl / 10.0, h22]])
    return grad, H


@pytest.mark.parametrize("lam", [0.1, 1.0])
def test_literature_grad_hess_vs_mc(lam):
    """Analytic E[grad V]/E[Hess V] match a large Monte Carlo estimate."""
    T = LiteratureLogconcaveTarget(lam)
    m = np.array([0.5, -0.3])
    C = np.array([[1.2, 0.3], [0.3, 0.9]])
    g, H = T.g_H(m, C)               # g = -E[grad V], H = -E[Hess V]
    mc_grad, mc_H = _mc_grad_hess_V(lam, m, C, seed=7)
    np.testing.assert_allclose(-g, mc_grad, rtol=0, atol=5e-3)
    np.testing.assert_allclose(-H, mc_H, rtol=0, atol=5e-3)


def test_literature_expected_V_moments():
    """E[V] uses the exact Gaussian moments E[Y^2], E[Y^4] and E[(sqrt(lam)X-Y)^2]."""
    lam = 0.1
    T = LiteratureLogconcaveTarget(lam)
    m = np.array([0.4, 0.7])
    C = np.array([[1.1, 0.2], [0.2, 0.8]])
    mx, my = m
    C11, C12, C22 = C[0, 0], C[0, 1], C[1, 1]
    sl = np.sqrt(lam)
    E_y4 = my ** 4 + 6 * my ** 2 * C22 + 3 * C22 ** 2
    mean_U = sl * mx - my
    var_U = lam * C11 - 2 * sl * C12 + C22
    expected = (mean_U ** 2 + var_U) / 20.0 + E_y4 / 20.0
    np.testing.assert_allclose(T.expected_V(m, C), expected, atol=1e-12)


def test_literature_no_theory_bound():
    """The quartic target reports no global (alpha, beta)."""
    T = LiteratureLogconcaveTarget(1.0)
    assert T.has_theory is False
    assert T.alpha is None and T.beta is None


# ---------------------------------------------------------------------------
# Target C: smooth log-concave Gauss--Hermite quadrature
# ---------------------------------------------------------------------------

def test_sech2_expectation_in_unit_interval():
    """E[sech^2(Y)] lies in (0, 1] for any Gaussian Y."""
    for mean, var in [(0.0, 1.0), (3.0, 4.0), (-2.0, 0.25), (10.0, 9.0)]:
        val = gaussian_expectation_1d(_sech2, mean, var, n_nodes=80)
        assert 0.0 < val <= 1.0 + 1e-12


@pytest.mark.parametrize("lam", [0.01, 0.1, 1.0])
def test_smooth_hessian_within_alpha_beta(lam):
    """The expected-Hessian spectrum lies in [alpha, beta] for several states."""
    T = SmoothLogconcaveTarget(lam)
    rng = np.random.default_rng(3)
    for _ in range(10):
        m = rng.standard_normal(2) * 2
        A = rng.standard_normal((2, 2))
        C = A @ A.T + 0.5 * np.eye(2)
        _, H = T.g_H(m, C)
        eig = np.linalg.eigvalsh(-H)     # -H = E[Hess V]
        assert eig[0] >= T.alpha - 1e-9
        assert eig[-1] <= T.beta + 1e-9


def test_smooth_alpha_beta_values():
    """alpha = delta, beta = delta + (1+lam)/10 + gamma."""
    T = SmoothLogconcaveTarget(0.1, delta=0.05, gamma=1.0)
    assert abs(T.alpha - 0.05) < 1e-15
    assert abs(T.beta - (0.05 + 1.1 / 10.0 + 1.0)) < 1e-15


def test_gauss_hermite_normalization():
    """Gauss--Hermite weights integrate constants and the variance exactly."""
    assert abs(gaussian_expectation_1d(lambda y: np.ones_like(y), 2.0, 3.0) - 1.0) < 1e-13
    assert abs(gaussian_expectation_1d(lambda y: y ** 2, 0.0, 5.0) - 5.0) < 1e-10
