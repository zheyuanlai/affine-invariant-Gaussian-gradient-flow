"""
Unit tests for the metrics module.

Each test targets a specific mathematical property of the metric definitions.
"""
import numpy as np
import pytest

from src.metrics import kl_energy, compute_all_metrics
from src.initializations import get_initialization
from src.utils import make_q_vector


# ---------------------------------------------------------------------------
# KL energy tests
# ---------------------------------------------------------------------------

def test_kl_energy_nonnegative():
    """KL energy must be >= 0 for any (m, C)."""
    rng = np.random.default_rng(0)
    for _ in range(20):
        n = rng.integers(2, 8)
        m = rng.standard_normal(n)
        A = rng.standard_normal((n, n))
        C = A @ A.T + np.eye(n) * 0.5
        assert kl_energy(m, C) >= -1e-12, "KL energy negative"


def test_kl_energy_zero_at_target():
    """KL energy must equal zero at the target m=0, C=I."""
    for n in [2, 5, 10]:
        E = kl_energy(np.zeros(n), np.eye(n))
        np.testing.assert_allclose(E, 0.0, atol=1e-14,
                                   err_msg=f"KL energy non-zero at target for n={n}")


def test_kl_energy_positive_away_from_target():
    """KL energy must be strictly positive when (m, C) != (0, I)."""
    n = 5
    assert kl_energy(np.ones(n), np.eye(n)) > 0
    assert kl_energy(np.zeros(n), 2.0 * np.eye(n)) > 0
    assert kl_energy(np.zeros(n), 0.5 * np.eye(n)) > 0


# ---------------------------------------------------------------------------
# Shape error tests
# ---------------------------------------------------------------------------

def test_shape_error_zero_for_scalar_covariance():
    """Shape error must be zero for C = s*I (pure scaling, no anisotropy)."""
    q = make_q_vector(5)
    b = 0.5
    for s in [0.1, 1.0, 2.5, 10.0]:
        m = np.zeros(5)
        C = s * np.eye(5)
        E0 = max(kl_energy(m, C), 1.0)
        mets = compute_all_metrics(m, C, E0, q, b)
        np.testing.assert_allclose(
            mets["shape_error"], 0.0, atol=1e-13,
            err_msg=f"Shape error non-zero for C={s}*I",
        )


# ---------------------------------------------------------------------------
# Volume error tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", [2, 5, 10])
def test_volume_error_zero_for_unit_determinant(n):
    """Volume error must be zero for shape_only initialisation (det(C0) = 1)."""
    q = make_q_vector(n)
    b = 0.5
    m, C = get_initialization("shape_only", n)
    E0 = max(kl_energy(m, C), 1.0)
    mets = compute_all_metrics(m, C, E0, q, b)
    np.testing.assert_allclose(
        mets["volume_error"], 0.0, atol=1e-13,
        err_msg=f"Volume error non-zero for shape_only n={n}",
    )


def test_volume_error_positive_for_volume_initializations():
    """volume_high and volume_low should have strictly positive volume error."""
    n = 5
    q = make_q_vector(n)
    b = 0.5
    for init_name in ["volume_high", "volume_low"]:
        m, C = get_initialization(init_name, n)
        E0 = max(kl_energy(m, C), 1.0)
        mets = compute_all_metrics(m, C, E0, q, b)
        assert mets["volume_error"] > 1e-10, \
            f"Expected positive volume error for {init_name}"


# ---------------------------------------------------------------------------
# Chi (trace dominance ratio) tests
# ---------------------------------------------------------------------------

def test_chi_equals_one_for_equal_residuals():
    """chi = 1 when all lambda_i equal the same value != 1 (max trace dominance)."""
    n = 6
    lam = 3.0   # uniform eigenvalue, residuals all = 1 - 3 = -2
    C = lam * np.eye(n)
    m = np.zeros(n)
    q = make_q_vector(n)
    E0 = max(kl_energy(m, C), 1.0)
    mets = compute_all_metrics(m, C, E0, q, b=0.5)
    np.testing.assert_allclose(mets["chi"], 1.0, atol=1e-14)


def test_chi_equals_one_over_n_for_single_residual():
    """chi = 1/n when only one eigenvalue differs from 1 (shape-dominated)."""
    n = 5
    diag_vals = np.ones(n)
    diag_vals[0] = 3.0   # only one non-unit eigenvalue
    C = np.diag(diag_vals)
    m = np.zeros(n)
    q = make_q_vector(n)
    E0 = max(kl_energy(m, C), 1.0)
    mets = compute_all_metrics(m, C, E0, q, b=0.5)
    np.testing.assert_allclose(mets["chi"], 1.0 / n, atol=1e-13)


def test_chi_is_one_at_target():
    """chi should be 1.0 at target (all residuals zero, denominator handled)."""
    n = 4
    m = np.zeros(n)
    C = np.eye(n)
    q = make_q_vector(n)
    E0 = 1.0
    mets = compute_all_metrics(m, C, E0, q, b=0.5)
    assert mets["chi"] == 1.0


# ---------------------------------------------------------------------------
# Cosine test-function
# ---------------------------------------------------------------------------

def test_cosine_error_zero_at_target():
    """Cosine error must be zero at (m=0, C=I)."""
    for n in [2, 5, 10]:
        m = np.zeros(n)
        C = np.eye(n)
        q = make_q_vector(n)
        E0 = 1.0
        mets = compute_all_metrics(m, C, E0, q, b=0.5)
        np.testing.assert_allclose(
            mets["cosine_error"], 0.0, atol=1e-14,
            err_msg=f"Cosine error non-zero at target for n={n}",
        )


# ---------------------------------------------------------------------------
# Normalised energy monotonicity (spot-check over a short run)
# ---------------------------------------------------------------------------

def test_norm_energy_bounded_by_one_at_step_zero():
    """norm_energy at step 0 should equal 1.0 (by construction)."""
    from src.dynamics import gaussian_step
    n = 5
    m, C = get_initialization("volume_high", n)
    q    = make_q_vector(n)
    E0   = max(kl_energy(m, C), 1e-15)
    mets = compute_all_metrics(m, C, E0, q, b=0.5)
    np.testing.assert_allclose(mets["norm_energy"], 1.0, atol=1e-14)
