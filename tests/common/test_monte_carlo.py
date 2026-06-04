"""Tests for the Monte Carlo sample-bank helpers."""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples, transform_gaussian_samples


@pytest.mark.parametrize("N,M", [(3, 8), (5, 100), (1, 4)])
def test_shape(N, M):
    Z = gaussian_samples(N, M, seed=0)
    assert Z.shape == (M, N)
    assert Z.dtype == np.float64


def test_reproducible():
    Z1 = gaussian_samples(4, 64, seed=123)
    Z2 = gaussian_samples(4, 64, seed=123)
    Z3 = gaussian_samples(4, 64, seed=124)
    assert np.array_equal(Z1, Z2)
    assert not np.array_equal(Z1, Z3)


def test_antithetic_mean_is_exactly_zero_even_M():
    Z = gaussian_samples(6, 100, seed=7, antithetic=True)
    # antithetic with even M -> sample mean is exactly zero
    assert np.allclose(Z.mean(axis=0), 0.0, atol=1e-14)
    # second half is the negation of the first half
    half = 50
    assert np.allclose(Z[:half], -Z[half:], atol=1e-14)


def test_non_antithetic_differs():
    Z_anti = gaussian_samples(4, 50, seed=1, antithetic=True)
    Z_plain = gaussian_samples(4, 50, seed=1, antithetic=False)
    assert not np.array_equal(Z_anti, Z_plain)


def test_odd_M_antithetic():
    Z = gaussian_samples(3, 11, seed=2, antithetic=True)
    assert Z.shape == (11, 3)


def test_transform_pushforward():
    N = 4
    Z = gaussian_samples(N, 2, seed=0, antithetic=False)
    m = np.arange(N, dtype=float)
    rng = np.random.default_rng(0)
    B = rng.standard_normal((N, N))
    C = B @ B.T + np.eye(N)
    Theta = transform_gaussian_samples(Z, m, C)
    assert Theta.shape == (2, N)
    # check theta_0 = m + C^{1/2} z_0 explicitly
    from src.common.spd import symmetric_sqrt
    expected0 = m + symmetric_sqrt(C) @ Z[0]
    assert np.allclose(Theta[0], expected0, atol=1e-10)


def test_transform_moments():
    # Large bank -> empirical mean/cov approximate (m, C)
    N = 3
    Z = gaussian_samples(N, 40000, seed=5, antithetic=True)
    m = np.array([1.0, -2.0, 0.5])
    C = np.array([[2.0, 0.3, 0.0], [0.3, 1.0, -0.2], [0.0, -0.2, 1.5]])
    Theta = transform_gaussian_samples(Z, m, C)
    assert np.allclose(Theta.mean(axis=0), m, atol=0.05)
    emp_cov = np.cov(Theta.T, bias=True)
    assert np.allclose(emp_cov, C, atol=0.1)
