"""
Unit tests for the log-concave extension.

Tests cover:
  1.  LogCoshTarget gradient and Hessian shapes
  2.  Hessian is SPD with min eigenvalue >= 1
  3.  batch_grad and batch_hess match single-point versions
  4.  Cholesky packing/unpacking round-trips
  5.  Reference objective gradient has correct shape and finite values
  6.  One logconcave_step preserves SPD covariance
  7.  Parameter validation rejects omega + n*tau <= 0
  8.  Log-concave initializations produce SPD covariance
  9.  At a_star, objective gap is approximately zero (MC noise tolerance)
 10.  volume_error is |log s| for C = s * C_star (with R = sI => log det R / n = log s)
 11.  shape_error is zero for C = s * C_star (R = sI => logR is scalar * I)
"""
import numpy as np
import pytest
import scipy.linalg

from src.targets import LogCoshTarget
from src.lc_dynamics import logconcave_step
from src.lc_initializations import get_logconcave_initialization, LC_INIT_NAMES
from src.lc_metrics import compute_lc_metrics, compute_lc_objective
from src.reference_optimum import _pack, _unpack, _objective_and_grad
from src.qmc_samples import make_samples, push_forward
from src.utils import spd_invsqrt, make_q_vector


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def target5():
    return LogCoshTarget(n=5, rho=3.0, seed=42)


@pytest.fixture(scope="module")
def target2():
    return LogCoshTarget(n=2, rho=1.0, seed=0)


@pytest.fixture(scope="module")
def Z5():
    return make_samples(5, 256, seed=0)


@pytest.fixture(scope="module")
def Z2():
    return make_samples(2, 256, seed=0)


# ---------------------------------------------------------------------------
# 1. Gradient and Hessian shapes
# ---------------------------------------------------------------------------

def test_grad_shape(target5):
    x = np.ones(5)
    g = target5.grad(x)
    assert g.shape == (5,)

def test_hess_shape(target5):
    x = np.ones(5)
    H = target5.hess(x)
    assert H.shape == (5, 5)

def test_batch_shapes(target5):
    X = np.random.default_rng(1).standard_normal((10, 5))
    assert target5.batch_grad(X).shape  == (10, 5)
    assert target5.batch_hess(X).shape  == (10, 5, 5)
    assert target5.batch_value(X).shape == (10,)


# ---------------------------------------------------------------------------
# 2. Hessian is SPD and min eigenvalue >= 1
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n,rho", [(2, 0.5), (5, 3.0), (5, 10.0)])
def test_hessian_min_eigenvalue(n, rho):
    """Hess V(x) >= I for all x, because sech^2 >= 0."""
    tgt = LogCoshTarget(n=n, rho=rho, seed=99)
    rng = np.random.default_rng(7)
    for _ in range(10):
        x = rng.standard_normal(n) * 3.0
        H = tgt.hess(x)
        # Must be symmetric
        np.testing.assert_allclose(H, H.T, atol=1e-14)
        eigs = np.linalg.eigvalsh(H)
        assert eigs.min() >= 1.0 - 1e-12, f"min eig = {eigs.min():.4f} < 1"


# ---------------------------------------------------------------------------
# 3. batch_grad and batch_hess match single-point versions
# ---------------------------------------------------------------------------

def test_batch_grad_matches_single(target5):
    X = np.random.default_rng(2).standard_normal((8, 5))
    batch_g  = target5.batch_grad(X)
    single_g = np.stack([target5.grad(X[i]) for i in range(8)])
    np.testing.assert_allclose(batch_g, single_g, atol=1e-12)

def test_batch_hess_matches_single(target2):
    X = np.random.default_rng(3).standard_normal((6, 2))
    batch_H  = target2.batch_hess(X)
    single_H = np.stack([target2.hess(X[i]) for i in range(6)])
    np.testing.assert_allclose(batch_H, single_H, atol=1e-12)


# ---------------------------------------------------------------------------
# 4. Cholesky packing / unpacking
# ---------------------------------------------------------------------------

def test_cholesky_pack_unpack_identity():
    n = 4
    m = np.array([1.0, -0.5, 2.0, 0.3])
    L = np.eye(n)
    params = _pack(m, L)
    m2, L2 = _unpack(params, n)
    np.testing.assert_allclose(m, m2, atol=1e-15)
    np.testing.assert_allclose(L, L2, atol=1e-15)

def test_cholesky_pack_unpack_random():
    rng = np.random.default_rng(5)
    n = 5
    m = rng.standard_normal(n)
    # Build a random lower-triangular L with positive diagonal
    L = np.tril(rng.standard_normal((n, n)))
    np.fill_diagonal(L, np.abs(np.diag(L)) + 0.5)
    params = _pack(m, L)
    m2, L2 = _unpack(params, n)
    np.testing.assert_allclose(m, m2, atol=1e-14)
    np.testing.assert_allclose(L, L2, atol=1e-14)


# ---------------------------------------------------------------------------
# 5. Reference objective gradient has correct shape and finite values
# ---------------------------------------------------------------------------

def test_ref_grad_shape_and_finite(target2, Z2):
    n = 2
    m0 = np.zeros(n)
    L0 = np.eye(n)
    params = _pack(m0, L0)
    F, grad = _objective_and_grad(params, target2, Z2)
    assert np.isfinite(F)
    assert grad.shape == params.shape
    assert np.all(np.isfinite(grad))


# ---------------------------------------------------------------------------
# 6. One logconcave_step preserves SPD covariance
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("omega,tau_factor", [(0.5, 0.0), (0.25, -0.4), (1.0, 0.3)])
def test_step_preserves_spd(target5, Z5, omega, tau_factor):
    """logconcave_step must return an SPD covariance."""
    n = 5
    tau = tau_factor * omega / n
    rng = np.random.default_rng(11)
    m = rng.standard_normal(n)
    A = rng.standard_normal((n, n))
    C = A @ A.T + 1.5 * np.eye(n)

    L = scipy.linalg.cholesky(C, lower=True)
    Theta = push_forward(m, L, Z5)
    g = np.mean(target5.batch_grad(Theta), axis=0)
    S = np.mean(target5.batch_hess(Theta), axis=0)
    S = 0.5 * (S + S.T)

    _, C_next = logconcave_step(m, C, g, S, dt=0.01, omega=omega, tau=tau)

    # Must be symmetric
    np.testing.assert_allclose(C_next, C_next.T, atol=1e-12)
    # All eigenvalues must be positive
    eigs = np.linalg.eigvalsh(C_next)
    assert eigs.min() > 0, f"Non-positive eigenvalue after step: {eigs.min():.4e}"


# ---------------------------------------------------------------------------
# 7. Parameter validation
# ---------------------------------------------------------------------------

def test_logconcave_step_rejects_bad_params(target5, Z5):
    n = 5
    m = np.zeros(n)
    C = np.eye(n)
    g = np.zeros(n)
    S = np.eye(n)

    with pytest.raises(ValueError):
        logconcave_step(m, C, g, S, dt=0.01, omega=0.0, tau=0.0)   # omega=0

    with pytest.raises(ValueError):
        logconcave_step(m, C, g, S, dt=0.01, omega=-1.0, tau=0.0)  # omega<0

    with pytest.raises(ValueError):
        logconcave_step(m, C, g, S, dt=0.01, omega=1.0, tau=-1.0/n)  # omega+n*tau=0


# ---------------------------------------------------------------------------
# 8. Log-concave initializations produce SPD covariance
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", LC_INIT_NAMES)
@pytest.mark.parametrize("n", [2, 5])
def test_lc_init_spd(name, n):
    """Every initialization must yield an SPD covariance."""
    rng = np.random.default_rng(99)
    A = rng.standard_normal((n, n))
    C_star = A @ A.T + np.eye(n)
    m_star = rng.standard_normal(n)

    _, C0 = get_logconcave_initialization(name, n, m_star, C_star)
    eigs = np.linalg.eigvalsh(C0)
    assert eigs.min() > 0, f"{name}: min eigenvalue = {eigs.min():.4e}"


# ---------------------------------------------------------------------------
# 9. At a_star, objective gap is approximately zero
# ---------------------------------------------------------------------------

def test_objective_gap_near_zero_at_star(target5, Z5):
    """At the reference optimum, objective gap should be small (MC noise)."""
    n = 5
    # Use the same target; find approx reference via short optimisation
    from src.reference_optimum import compute_reference_optimum
    Z_ref = make_samples(n, 1024, seed=77)
    ref = compute_reference_optimum(target5, Z_ref, maxiter=500, gtol=1e-4)
    m_star, L_star = ref["m_star"], ref["L_star"]
    F_star = ref["F_star"]
    obj_at_star = compute_lc_objective(m_star, L_star, Z5, target5)
    # Gap should be small relative to F_star; exact zero only with infinite samples
    gap = abs(obj_at_star - F_star)
    assert gap < 2.0, f"Objective gap at a_star too large: {gap:.4f}"


# ---------------------------------------------------------------------------
# 10. volume_error = |log s| for C = s * C_star
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("s", [0.25, 0.5, 1.0, 2.0, 4.0])
def test_volume_error_scaled_cstar(s):
    """For C = s * C_star,  R = sI,  volume_error = |log(s)|."""
    n = 4
    rng = np.random.default_rng(33)
    A = rng.standard_normal((n, n))
    C_star = A @ A.T + np.eye(n)
    C = s * C_star
    C_star_invsqrt = spd_invsqrt(C_star)
    q = make_q_vector(n)
    # dummy values for unused args
    m = np.zeros(n); g = np.zeros(n); S = np.eye(n)
    mets = compute_lc_metrics(m, C, g, S, np.zeros(n), C_star,
                               C_star_invsqrt, 0.0, 0.0, 1.0, q, 0.5)
    np.testing.assert_allclose(mets["volume_error"], abs(np.log(s)), atol=1e-12)


# ---------------------------------------------------------------------------
# 11. shape_error = 0 for C = s * C_star
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("s", [0.25, 1.0, 3.0])
def test_shape_error_zero_for_scaled_cstar(s):
    """For C = s * C_star,  R = sI  is isotropic,  shape_error = 0."""
    n = 5
    rng = np.random.default_rng(44)
    A = rng.standard_normal((n, n))
    C_star = A @ A.T + np.eye(n)
    C = s * C_star
    C_star_invsqrt = spd_invsqrt(C_star)
    q = make_q_vector(n)
    m = np.zeros(n); g = np.zeros(n); S = np.eye(n)
    mets = compute_lc_metrics(m, C, g, S, np.zeros(n), C_star,
                               C_star_invsqrt, 0.0, 0.0, 1.0, q, 0.5)
    np.testing.assert_allclose(mets["shape_error"], 0.0, atol=1e-11)
