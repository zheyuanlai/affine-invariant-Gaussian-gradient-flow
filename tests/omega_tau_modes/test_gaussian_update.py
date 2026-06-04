"""
Unit tests for the Gaussian gradient flow update step.

Tests verify mathematical properties that must hold exactly (up to floating-
point precision) for the closed-form Gaussian-target dynamics.
"""
import numpy as np
import pytest

from src.omega_tau_modes.dynamics import gaussian_step
from src.omega_tau_modes.utils import validate_params


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_steps(m0, C0, dt, omega, tau, n_steps=50):
    m, C = m0.copy(), C0.copy()
    for _ in range(n_steps):
        m, C = gaussian_step(m, C, dt, omega, tau)
    return m, C


# ---------------------------------------------------------------------------
# Test 1: C0 = I with nonzero mean — covariance must stay I
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", [2, 5, 10])
@pytest.mark.parametrize("omega", [0.25, 1.0, 2.0])
@pytest.mark.parametrize("tau_factor", [-0.4, 0.0, 0.4])
def test_identity_covariance_is_invariant(n, omega, tau_factor):
    """Starting from C = I, covariance must remain I for all (omega, tau)."""
    tau = tau_factor * omega / n  # ensures omega + n*tau != 0

    m0 = np.ones(n, dtype=np.float64) * 3.0 / np.sqrt(n)
    C0 = np.eye(n, dtype=np.float64)

    _, C_final = run_steps(m0, C0, dt=0.05, omega=omega, tau=tau, n_steps=100)

    # C must still equal I
    np.testing.assert_allclose(
        C_final, np.eye(n),
        atol=1e-12,
        err_msg=f"Covariance drifted from I for n={n}, omega={omega}, tau_factor={tau_factor}",
    )


# ---------------------------------------------------------------------------
# Test 2: Isotropic covariance C = lambda*I stays isotropic after one step
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", [2, 5])
@pytest.mark.parametrize("lam", [0.25, 0.5, 2.0, 4.0])
@pytest.mark.parametrize("omega", [0.5, 1.0])
@pytest.mark.parametrize("tau_factor", [-0.3, 0.0, 0.3])
def test_isotropic_covariance_stays_isotropic(n, lam, omega, tau_factor):
    """Isotropic C = lambda*I must have all eigenvalues equal after each step."""
    tau = tau_factor * omega / n

    m0 = np.zeros(n, dtype=np.float64)
    C0 = lam * np.eye(n, dtype=np.float64)

    m1, C1 = gaussian_step(m0, C0, dt=0.02, omega=omega, tau=tau)

    eigvals = np.linalg.eigvalsh(C1)
    np.testing.assert_allclose(
        eigvals, np.full(n, eigvals[0]),
        rtol=1e-12,
        err_msg=f"Eigenvalues not equal after step for n={n}, lam={lam}",
    )


# ---------------------------------------------------------------------------
# Test 3: Parameter validation rejects omega + n*tau <= 0
# ---------------------------------------------------------------------------

def test_invalid_params_negative_denominator():
    """gaussian_step must raise ValueError when omega + n*tau <= 0."""
    n, omega = 5, 1.0
    # tau = -omega/n makes omega + n*tau = 0 (boundary, should still raise)
    tau_boundary = -omega / n
    tau_invalid  = -omega / n - 0.01

    m = np.zeros(n, dtype=np.float64)
    C = np.eye(n, dtype=np.float64)

    with pytest.raises(ValueError):
        gaussian_step(m, C, dt=0.02, omega=omega, tau=tau_boundary)

    with pytest.raises(ValueError):
        gaussian_step(m, C, dt=0.02, omega=omega, tau=tau_invalid)


def test_invalid_params_nonpositive_omega():
    """gaussian_step must raise ValueError when omega <= 0."""
    n = 4
    m = np.zeros(n, dtype=np.float64)
    C = np.eye(n, dtype=np.float64)

    with pytest.raises(ValueError):
        gaussian_step(m, C, dt=0.02, omega=0.0, tau=0.0)

    with pytest.raises(ValueError):
        gaussian_step(m, C, dt=0.02, omega=-1.0, tau=0.0)


def test_validate_params_grid():
    """validate_params agrees with the analytical validity conditions."""
    assert validate_params(omega=1.0, tau=0.0,   n=5) is True
    assert validate_params(omega=0.5, tau=0.1,   n=5) is True
    assert validate_params(omega=0.5, tau=-0.09, n=5) is True   # 0.5 - 0.45 > 0
    assert validate_params(omega=0.5, tau=-0.10, n=5) is False  # 0.5 - 0.50 = 0
    assert validate_params(omega=0.5, tau=-0.11, n=5) is False  # 0.5 - 0.55 < 0
    assert validate_params(omega=0.0, tau=0.0,   n=5) is False
    assert validate_params(omega=-1., tau=0.0,   n=5) is False


# ---------------------------------------------------------------------------
# Test 4: Mean converges to zero, covariance converges to I
#          (long run, should reach machine-precision neighbourhood)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n", [2, 5])
def test_convergence_to_target(n):
    """After many steps, (m, C) should be very close to (0, I)."""
    omega, tau = 0.5, 0.0
    m0 = 3.0 * np.ones(n) / np.sqrt(n)
    C0 = 4.0 * np.eye(n)

    m, C = run_steps(m0, C0, dt=0.05, omega=omega, tau=tau, n_steps=500)

    assert np.linalg.norm(m) < 1e-6,  f"Mean did not converge: ||m||={np.linalg.norm(m)}"
    assert np.linalg.norm(C - np.eye(n), 'fro') < 1e-6, \
        f"Covariance did not converge: ||C-I||={np.linalg.norm(C - np.eye(n))}"
