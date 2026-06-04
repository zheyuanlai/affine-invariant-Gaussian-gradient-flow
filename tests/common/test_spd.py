"""Tests for SPD matrix utilities."""
import numpy as np
import pytest
import scipy.linalg

from src.common.spd import (
    symmetrize, is_spd, eigh_spd, check_spd,
    symmetric_sqrt, symmetric_invsqrt, symmetric_expm,
)


def _random_spd(N, rng, eps=0.5):
    A = rng.standard_normal((N, N))
    return A @ A.T + eps * np.eye(N)


@pytest.mark.parametrize("N", [1, 2, 4, 6])
def test_sqrt_squares_to_original(N):
    rng = np.random.default_rng(N)
    A = _random_spd(N, rng)
    S = symmetric_sqrt(A)
    assert np.allclose(S, S.T, atol=1e-12)
    assert np.allclose(S @ S, A, atol=1e-9)


@pytest.mark.parametrize("N", [1, 2, 4, 6])
def test_invsqrt_inverse_of_sqrt(N):
    rng = np.random.default_rng(10 + N)
    A = _random_spd(N, rng)
    Si = symmetric_invsqrt(A)
    assert np.allclose(Si @ A @ Si, np.eye(N), atol=1e-8)
    assert np.allclose(Si @ symmetric_sqrt(A), np.eye(N), atol=1e-8)


@pytest.mark.parametrize("N", [1, 2, 3, 5])
def test_expm_matches_scipy(N):
    rng = np.random.default_rng(20 + N)
    A = symmetrize(0.3 * rng.standard_normal((N, N)))  # small symmetric
    assert np.allclose(symmetric_expm(A), scipy.linalg.expm(A), atol=1e-9)


def test_check_spd_catches_non_spd():
    # Negative definite
    with pytest.raises(ValueError):
        check_spd(-np.eye(3), name="C")
    # Indefinite
    with pytest.raises(ValueError):
        check_spd(np.diag([1.0, -1e-3, 2.0]))
    # Singular (zero eigenvalue)
    with pytest.raises(ValueError):
        check_spd(np.diag([1.0, 0.0, 2.0]))
    # Valid SPD returns (min, max) eigenvalues
    lo, hi = check_spd(np.diag([2.0, 5.0]))
    assert lo == pytest.approx(2.0)
    assert hi == pytest.approx(5.0)


def test_is_spd():
    assert is_spd(np.eye(4))
    assert not is_spd(-np.eye(4))
    assert not is_spd(np.diag([1.0, 0.0]))


def test_eigh_floor():
    A = np.diag([-1.0, 2.0, 3.0])
    w, V = eigh_spd(A, floor=0.0)
    assert w[0] == pytest.approx(0.0)
    assert np.allclose(V @ V.T, np.eye(3), atol=1e-12)


def test_symmetrize_defensive():
    A = np.array([[1.0, 2.0], [0.0, 1.0]])
    S = symmetrize(A)
    assert np.allclose(S, [[1.0, 1.0], [1.0, 1.0]])
