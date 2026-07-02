"""Unit tests for the covariance-bootstrap targets, methods, and envelopes.

Covers: the Gaussian energy gap vanishing at ``(0, H^{-1})``; the KL envelope
closed form vs the recurrence; the actual KL update lower-bounded by ``L_n``; SPD
preservation of the Riemannian update; the Bures bootstrap lift
``lambda_min(C_b) >= c/beta``; ``Psi == 0`` for the Gaussian target and ``> 0`` for
log-cosh; the curvature/discretization sign equivalence; and the log-cosh optimum.
"""
import numpy as np
import pytest

from src.common.spd import symmetrize, eigh_spd
from src.natural_gradient_discretization_stepsize.methods import (
    riemannian_cov_step as disc_riem, kl_cov_step as disc_kl,
)
from src.natural_gradient_covariance_bootstrap.targets import build_target
from src.natural_gradient_covariance_bootstrap.methods import (
    riemannian_cov_step, kl_cov_step, mean_step, bures_bootstrap_cov_step,
    bures_bootstrap_step,
)
from src.natural_gradient_covariance_bootstrap.envelopes import (
    envelope_L0, kl_envelope_step, kl_envelope_closed, kl_envelope_sequence,
    riemannian_envelope_sequence, frozen_lower_bound,
)


# ---------------------------------------------------------------------------
# 1. Gaussian energy gap vanishes at the optimum
# ---------------------------------------------------------------------------

def test_gaussian_energy_gap_zero_at_optimum():
    T = build_target("gaussian", d=8, kappa=1000.0)
    m_star, C_star = T.a_star()
    assert abs(T.energy_gap(m_star, C_star)) < 1e-10
    # and it is strictly positive away from the optimum.
    assert T.energy_gap(m_star, 0.3 * C_star) > 1e-6
    assert T.energy_gap(m_star + 1.0, C_star) > 1e-6


def test_gaussian_oracles():
    T = build_target("gaussian", d=5, kappa=100.0)
    rng = np.random.default_rng(0)
    m = rng.standard_normal(5)
    C = np.diag(rng.uniform(0.5, 2.0, 5))
    np.testing.assert_allclose(T.G(m, C), T.H_diag * m, atol=1e-12)
    np.testing.assert_allclose(np.diag(T.A_matrix(m, C)), T.H_diag, atol=1e-12)


# ---------------------------------------------------------------------------
# 2. KL envelope closed form matches the recurrence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dt,beta,lam0", [(0.5, 100.0, 1e-9), (0.1, 10.0, 1e-4),
                                          (1.0, 1000.0, 1e-6)])
def test_kl_envelope_closed_form_matches_recurrence(dt, beta, lam0):
    seq = kl_envelope_sequence(60, lam0, dt, beta)
    L0 = envelope_L0(lam0, beta)
    n = np.arange(seq.size)
    closed = kl_envelope_closed(n, L0, dt, beta)
    np.testing.assert_allclose(seq, closed, rtol=1e-12, atol=1e-18)
    # and the recurrence step reproduces the sequence.
    for k in range(seq.size - 1):
        assert abs(kl_envelope_step(seq[k], dt, beta) - seq[k + 1]) < 1e-15
    # envelope monotonically increases to 1/beta.
    assert seq[0] <= seq[-1] <= 1.0 / beta + 1e-15


# ---------------------------------------------------------------------------
# 3. Actual KL covariance update satisfies lambda_min(C_n) >= L_n
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("dt", [0.25, 0.5, 1.0])
def test_kl_update_lower_bounded_by_envelope(dt):
    T = build_target("gaussian", d=6, kappa=1000.0)
    beta = T.beta
    lam0 = 1e-9
    C = lam0 * np.eye(T.d)
    A = T.A_matrix(np.zeros(T.d), C)
    Ls = kl_envelope_sequence(40, lam0, dt, beta)
    for n in range(1, 41):
        C = kl_cov_step(C, A, dt)
        lam_min = eigh_spd(C)[0][0]
        assert lam_min >= Ls[n] - 1e-14, f"n={n}: {lam_min} < {Ls[n]}"


def test_riemannian_envelope_lower_bounds_short_run():
    T = build_target("gaussian", d=6, kappa=100.0)
    beta = T.beta
    lam0, dt = 1e-8, 0.5
    C = lam0 * np.eye(T.d)
    A = T.A_matrix(np.zeros(T.d), C)
    Ls = riemannian_envelope_sequence(40, lam0, dt, beta)
    for n in range(1, 41):
        C = riemannian_cov_step(C, A, dt)
        assert eigh_spd(C)[0][0] >= Ls[n] - 1e-14


# ---------------------------------------------------------------------------
# 4. Riemannian covariance update preserves SPD
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", [0, 1, 2])
def test_riemannian_preserves_spd(seed):
    rng = np.random.default_rng(seed)
    d = 5
    B = rng.standard_normal((d, d))
    C = symmetrize(B @ B.T + 0.1 * np.eye(d))
    G = rng.standard_normal((d, d))
    A = symmetrize(G @ G.T + 0.05 * np.eye(d))   # SPD curvature
    C_next = riemannian_cov_step(C, A, 0.4)
    w = np.linalg.eigvalsh(symmetrize(C_next))
    assert w[0] > 0
    np.testing.assert_allclose(C_next, C_next.T, atol=1e-12)


# ---------------------------------------------------------------------------
# Curvature <-> discretization (H_disc = -A) equivalence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", [0, 1])
def test_methods_equal_disc_convention(seed):
    rng = np.random.default_rng(seed)
    d = 4
    B = rng.standard_normal((d, d))
    C = symmetrize(B @ B.T + np.eye(d))
    G = rng.standard_normal((d, d))
    A = symmetrize(G @ G.T + 0.1 * np.eye(d))    # SPD curvature
    dt = 0.3
    np.testing.assert_allclose(riemannian_cov_step(C, A, dt),
                               disc_riem(C, -A, dt), atol=1e-10)
    np.testing.assert_allclose(kl_cov_step(C, A, dt), disc_kl(C, -A, dt), atol=1e-10)


# ---------------------------------------------------------------------------
# 5. Wasserstein/Bures bootstrap lifts lambda_min(C_b) >= c/beta
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("c", [0.25, 0.5, 1.0])
def test_bures_bootstrap_lifts_covariance(c):
    T = build_target("gaussian", d=10, kappa=1000.0)
    beta = T.beta
    eta = c / beta
    C0 = 1e-10 * np.eye(T.d)                     # tiny covariance
    A = T.A_matrix(np.zeros(T.d), C0)
    C_b = bures_bootstrap_cov_step(C0, A, eta)
    lam_min = np.linalg.eigvalsh(symmetrize(C_b))[0]
    assert lam_min >= eta - 1e-16
    # the full step also moves the mean.
    G = T.G(np.full(T.d, 2.0), C0)
    m_b, C_b2 = bures_bootstrap_step(np.full(T.d, 2.0), C0, G, A, eta)
    np.testing.assert_allclose(m_b, np.full(T.d, 2.0) - eta * G, atol=1e-12)
    assert np.linalg.eigvalsh(symmetrize(C_b2))[0] >= eta - 1e-16


# ---------------------------------------------------------------------------
# 6. Psi == 0 for the Gaussian target, > 0 for log-cosh
# ---------------------------------------------------------------------------

def test_gaussian_psi_zero():
    T = build_target("gaussian", d=6, kappa=100.0)
    m_star, C_star = T.a_star()
    assert T.Psi(m_star, C_star) == 0.0
    assert T.Psi(np.ones(6), 0.5 * C_star) == 0.0


def test_logcosh_psi_positive_and_metadata():
    T = build_target("log_cosh", d=5, kappa=50.0, gamma=1.0)
    m_star, C_star = T.a_star()
    psi = T.Psi(m_star, C_star)
    assert psi > 0.0
    md = T.metadata()
    assert md["Psi_star"] > 0.0
    assert md["a_star_diagnostics"]["stationarity_residual_max"] < 1e-9


def test_logcosh_astar_stationarity():
    T = build_target("log_cosh", d=5, kappa=100.0, gamma=1.0)
    m_star, C_star = T.a_star()
    np.testing.assert_allclose(m_star, 0.0, atol=1e-14)
    # C_star^{-1} = A(m_star, C_star) (the natural-gradient stationarity condition).
    np.testing.assert_allclose(1.0 / np.diag(C_star), T.A_diag(m_star, C_star),
                               atol=1e-9)


# ---------------------------------------------------------------------------
# Mean step sign and frozen lower bound
# ---------------------------------------------------------------------------

def test_mean_step_moves_toward_zero():
    T = build_target("gaussian", d=3, kappa=10.0)
    m = np.array([1.0, 1.0, 1.0])
    C = np.diag(1.0 / T.H_diag)                  # C = C_star
    m1 = mean_step(m, C, T.G(m, C), dt=0.5)
    # m - dt C H m = (1 - dt) m at C = H^{-1}.
    np.testing.assert_allclose(m1, 0.5 * m, atol=1e-12)


def test_frozen_lower_bound():
    assert frozen_lower_bound(1e-8, 100.0) == 1e-8
    assert frozen_lower_bound(1.0, 100.0) == 1.0 / 100.0
