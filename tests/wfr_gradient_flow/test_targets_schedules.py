"""Tests for the WFR targets and transport schedules.

Covers:
* Gaussian target expectations ``g = -Gamma^{-1} m``, ``H = -Gamma^{-1}`` and the
  exact energy gap;
* smooth log-cosh target Hessian bounds ``alpha I <= Hess V <= beta I`` and the
  eigenvalue beta formula;
* theory schedule ``h = mu_min`` as the maximizer of the proven W contribution;
* adaptive schedule monotonic decay in the calibration ratio ``s``.
"""
import numpy as np
import pytest

from src.wfr_gradient_flow.targets import (
    GaussianTarget, SmoothLogCoshTarget, smooth_beta, flat_valley_init,
)
from src.wfr_gradient_flow.schedules import (
    theory_mu_min, calibration_ratio, build_schedule,
)


# ---------------------------------------------------------------------------
# Gaussian target
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("Lambda", [100.0, 1000.0])
def test_gaussian_grad_hess_exact(Lambda):
    """g = -Gamma^{-1} m and H = -Gamma^{-1} exactly, for any state."""
    T = GaussianTarget(Lambda, epsilon=1e-2)
    rng = np.random.default_rng(0)
    m = rng.standard_normal(2)
    C = np.array([[1.3, 0.2], [0.2, 0.7]])
    g, H = T.g_H(m, C)
    np.testing.assert_allclose(g, -(T.Gamma_inv @ m), atol=1e-13)
    np.testing.assert_allclose(H, -T.Gamma_inv, atol=1e-13)
    # Gamma is the covariance; precision is diag(1, 1/Lambda).
    np.testing.assert_allclose(T.Gamma_inv, np.diag([1.0, 1.0 / Lambda]), atol=1e-13)


@pytest.mark.parametrize("Lambda", [100.0, 1000.0])
def test_gaussian_star_and_gap(Lambda):
    """Optimum is m=0, C=Gamma; the energy gap vanishes there and is exact."""
    T = GaussianTarget(Lambda, epsilon=1e-2)
    m_star, C_star, F_star = T.star()
    np.testing.assert_allclose(m_star, np.zeros(2), atol=1e-13)
    np.testing.assert_allclose(C_star, np.diag([1.0, Lambda]), atol=1e-10)
    assert abs(T.energy_gap(m_star, C_star)) < 1e-10
    m = np.array([2.0, -1.0])
    C = np.array([[0.8, 0.1], [0.1, 1.5]])
    np.testing.assert_allclose(
        T.energy_gap(m, C), T.objective(m, C) - F_star, atol=1e-9)


def test_gaussian_alpha_beta():
    """alpha = 1/Lambda, beta = 1 (spectrum of Gamma^{-1} = diag(1, 1/Lambda))."""
    T = GaussianTarget(100.0, epsilon=1e-2)
    assert abs(T.alpha - 1.0 / 100.0) < 1e-14
    assert abs(T.beta - 1.0) < 1e-14


def test_gaussian_energy_gap_nonnegative():
    """The analytic energy gap is a KL divergence, hence non-negative."""
    T = GaussianTarget(100.0, epsilon=1e-2)
    rng = np.random.default_rng(1)
    for _ in range(20):
        m = rng.standard_normal(2) * 3
        A = rng.standard_normal((2, 2))
        C = A @ A.T + 0.1 * np.eye(2)
        assert T.energy_gap(m, C) >= -1e-9


# ---------------------------------------------------------------------------
# Initialization
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("Lambda", [100.0, 1000.0])
@pytest.mark.parametrize("eps", [1e-2, 1e-3])
def test_flat_valley_init(Lambda, eps):
    """C0 = eps I and m0 has norm r along (1, sqrt(Lambda))/sqrt(1+Lambda)."""
    m0, C0 = flat_valley_init(Lambda, eps, r=8.0)
    np.testing.assert_allclose(C0, eps * np.eye(2), atol=1e-15)
    assert abs(np.linalg.norm(m0) - 8.0) < 1e-10
    direction = np.array([1.0, np.sqrt(Lambda)]) / np.sqrt(1.0 + Lambda)
    np.testing.assert_allclose(m0 / np.linalg.norm(m0), direction, atol=1e-12)


# ---------------------------------------------------------------------------
# Smooth log-cosh target
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("Lambda", [100.0, 1000.0])
def test_smooth_hessian_within_alpha_beta(Lambda):
    """E[Hess V] spectrum lies in [alpha, beta] for several states."""
    T = SmoothLogCoshTarget(Lambda, epsilon=1e-2)
    rng = np.random.default_rng(3)
    for _ in range(10):
        m = rng.standard_normal(2) * 2
        A = rng.standard_normal((2, 2))
        C = A @ A.T + 0.5 * np.eye(2)
        _, H = T.g_H(m, C)
        eig = np.linalg.eigvalsh(-H)        # -H = E[Hess V]
        assert eig[0] >= T.alpha - 1e-9
        assert eig[-1] <= T.beta + 1e-9


@pytest.mark.parametrize("Lambda", [100.0, 1000.0])
def test_smooth_beta_is_max_eig(Lambda):
    """beta equals lambda_max of the worst-case (sech^2=1) Hessian matrix."""
    delta, gamma = 0.05, 1.0
    M = np.array([
        [delta + Lambda / 10.0, -np.sqrt(Lambda) / 10.0],
        [-np.sqrt(Lambda) / 10.0, delta + 1.0 / 10.0 + gamma],
    ])
    expected = float(np.linalg.eigvalsh(M)[-1])
    assert abs(smooth_beta(Lambda, delta, gamma) - expected) < 1e-12
    T = SmoothLogCoshTarget(Lambda, epsilon=1e-2)
    assert abs(T.beta - expected) < 1e-12
    assert abs(T.alpha - delta) < 1e-15


def test_smooth_beta_tighter_than_trace_bound():
    """The eigenvalue beta is no larger than the loose trace-style bound."""
    Lambda, delta, gamma = 1000.0, 0.05, 1.0
    loose = delta + (1.0 + Lambda) / 10.0 + gamma
    assert smooth_beta(Lambda, delta, gamma) <= loose + 1e-9


# ---------------------------------------------------------------------------
# Theory schedule h = mu_min maximizes f(h) = alpha h / (1 + h/mu_min)^2
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("eps", [1e-2, 1e-3])
def test_theory_mu_min_value(eps):
    """mu_min = min(lambda_min(C0), 1/beta); here C0=eps I and beta>=1 so =eps."""
    T = GaussianTarget(100.0, epsilon=eps)     # beta = 1
    mu = theory_mu_min(T.C0, T.beta)
    assert abs(mu - min(eps, 1.0 / T.beta)) < 1e-14


def test_theory_mu_min_maximizes_contribution():
    """h = mu_min numerically maximizes alpha h / (1 + h/mu_min)^2."""
    alpha, mu_min = 0.3, 0.05
    f = lambda h: alpha * h / (1.0 + h / mu_min) ** 2
    grid = np.linspace(1e-4, 20.0 * mu_min, 20001)
    h_best = grid[int(np.argmax(f(grid)))]
    assert abs(h_best - mu_min) < 5e-4
    # The closed-form maximizer beats any other grid point.
    assert f(mu_min) >= f(grid).max() - 1e-12


# ---------------------------------------------------------------------------
# Adaptive schedule: large when s << s0, small when s >> s0
# ---------------------------------------------------------------------------

def test_adaptive_decays_with_calibration():
    """h_n -> h_max as s -> 0 and decays monotonically as s grows."""
    T = GaussianTarget(100.0, epsilon=1e-2)
    sched = build_schedule("wfr_adaptive", beta=T.beta, C0=T.C0,
                           h_max_frac=0.9, s0=0.5)
    h_max = 0.9 / T.beta
    # Underdispersed: tiny C => s << s0 => h near h_max.
    C_small = 1e-4 * np.eye(2)
    _, H = T.g_H(T.m0, C_small)
    h_small, info = sched.h(0, T.m0, C_small, None, H)
    assert info["s"] < 0.5
    assert h_small > 0.8 * h_max
    # Calibrated: C ~ Gamma => s ~ O(1) => h throttled well below h_max.
    C_cal = T.Gamma.copy()
    _, H2 = T.g_H(np.zeros(2), C_cal)
    h_cal, info2 = sched.h(0, np.zeros(2), C_cal, None, H2)
    assert info2["s"] > info["s"]
    assert h_cal < h_small


def test_calibration_ratio_units():
    """s = lambda_min(C^{1/2}(-H)C^{1/2}); for C=Gamma, -H=Gamma^{-1} => s=1."""
    T = GaussianTarget(100.0, epsilon=1e-2)
    _, H = T.g_H(np.zeros(2), T.Gamma)
    s = calibration_ratio(T.Gamma, H)          # C^{1/2} Gamma^{-1} C^{1/2} = I
    assert abs(s - 1.0) < 1e-9
