"""Tests for the STL targets, estimators, and the two experiment cores.

Covers the target score/Hessian formulas, STL unbiasedness and the pointwise-zero
property at the Gaussian optimum, the Hessian-residual covariance-update
equivalence for both schemes, the deterministic log-cosh optimum, and a short
algorithm run showing STL lowers the stochastic noise floor.
"""
import numpy as np
import pytest

from src.common.spd import symmetrize
from src.common.torch_utils import torch_available
from src.natural_gradient_stl_variance.targets import build_target
from src.natural_gradient_stl_variance.estimators import (
    stl_mean_estimator, stl_correction,
    riemannian_cov_step, kl_cov_step,
    riemannian_cov_step_residual, kl_cov_step_residual,
    residual_from_sampled, method_scheme_stl, METHODS,
)
from src.natural_gradient_stl_variance.linalg import ArrayBackend
from src.natural_gradient_stl_variance.states import build_states
from src.natural_gradient_stl_variance.estimator_variance import evaluate_state
from src.natural_gradient_stl_variance.algorithm import simulate_cell


# ---------------------------------------------------------------------------
# Target score / Hessian formulas (vs finite differences of V)
# ---------------------------------------------------------------------------

def _V_gaussian(theta, a):
    return 0.5 * float(np.sum(a * theta ** 2))


def _V_logcosh(theta, a, tau):
    return 0.5 * float(np.sum(a * theta ** 2)) + tau * float(np.sum(np.log(np.cosh(theta))))


def _num_grad(fV, theta, eps=1e-6):
    g = np.zeros_like(theta)
    for i in range(theta.size):
        tp = theta.copy(); tp[i] += eps
        tm = theta.copy(); tm[i] -= eps
        g[i] = (fV(tp) - fV(tm)) / (2 * eps)
    return g


def _num_hess_diag(fV, theta, eps=1e-4):
    h = np.zeros_like(theta)
    f0 = fV(theta)
    for i in range(theta.size):
        tp = theta.copy(); tp[i] += eps
        tm = theta.copy(); tm[i] -= eps
        h[i] = (fV(tp) - 2 * f0 + fV(tm)) / eps ** 2
    return h


def test_gaussian_score_hessian():
    T = build_target("gaussian", d=5, kappa=50.0)
    rng = np.random.default_rng(0)
    theta = rng.standard_normal(5)
    # score = -A theta = -grad V; Hess_log_post = -A.
    np.testing.assert_allclose(T.score(theta), -(T.a_diag * theta), atol=1e-12)
    np.testing.assert_allclose(T.score(theta),
                               -_num_grad(lambda t: _V_gaussian(t, T.a_diag), theta),
                               atol=1e-6)
    np.testing.assert_allclose(T.hess_diag(theta), -T.a_diag, atol=1e-12)


def test_logcosh_score_hessian():
    T = build_target("log_cosh", d=5, kappa=50.0, tau=0.7)
    rng = np.random.default_rng(1)
    theta = rng.standard_normal(5)
    expected_score = -(T.a_diag * theta) - T.tau * np.tanh(theta)
    expected_hess = -(T.a_diag + T.tau * (1.0 - np.tanh(theta) ** 2))
    np.testing.assert_allclose(T.score(theta), expected_score, atol=1e-12)
    np.testing.assert_allclose(
        T.score(theta), -_num_grad(lambda t: _V_logcosh(t, T.a_diag, T.tau), theta),
        atol=1e-6)
    np.testing.assert_allclose(T.hess_diag(theta), expected_hess, atol=1e-12)
    np.testing.assert_allclose(
        T.hess_diag(theta), -_num_hess_diag(lambda t: _V_logcosh(t, T.a_diag, T.tau), theta),
        atol=1e-4)


# ---------------------------------------------------------------------------
# STL unbiasedness and pointwise-zero at the Gaussian optimum
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kind,tau", [("gaussian", 0.0), ("log_cosh", 1.0)])
def test_stl_unbiased(kind, tau):
    """The STL mean estimator is unbiased for E_q[score_post] (MC check)."""
    T = build_target(kind, d=3, kappa=10.0, tau=tau)
    m = np.array([0.3, -0.5, 0.8])
    C = np.diag([0.7, 1.2, 0.4])
    rng = np.random.default_rng(2)
    Z = rng.standard_normal((500000, 3))
    theta = m + Z @ np.linalg.cholesky(C).T
    b_stl = stl_mean_estimator(T.score(theta), theta, m, C)
    ref = T.exp_score(m, C)
    np.testing.assert_allclose(b_stl.mean(axis=0), ref, atol=0.03)


def test_stl_pointwise_zero_at_gaussian_optimum():
    """At the exact Gaussian optimum the STL estimator is pointwise zero."""
    T = build_target("gaussian", d=6, kappa=1000.0)
    m_star, C_star = T.a_star()
    rng = np.random.default_rng(3)
    theta = m_star + rng.standard_normal((50, 6)) @ np.linalg.cholesky(C_star).T
    b_stl = stl_mean_estimator(T.score(theta), theta, m_star, C_star)
    assert np.max(np.abs(b_stl)) < 1e-9


def test_stl_correction_zero_mean():
    """E_q[C^{-1}(theta - m)] = 0 numerically (the source of unbiasedness)."""
    m = np.array([1.0, -2.0])
    C = np.array([[2.0, 0.3], [0.3, 1.0]])
    rng = np.random.default_rng(4)
    theta = m + rng.standard_normal((200000, 2)) @ np.linalg.cholesky(C).T
    corr = stl_correction(theta, m, C)
    np.testing.assert_allclose(corr.mean(axis=0), np.zeros(2), atol=5e-3)


# ---------------------------------------------------------------------------
# Hessian-residual covariance-update equivalence (both schemes)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", [0, 1, 2])
def test_residual_equivalence_both_schemes(seed):
    """Direct update with S equals the residual update with K = S + C^{-1}."""
    rng = np.random.default_rng(seed)
    d = 4
    A = rng.standard_normal((d, d))
    C = symmetrize(A @ A.T + np.eye(d))
    B = rng.standard_normal((d, d))
    S = -symmetrize(B @ B.T + 0.1 * np.eye(d))   # negative definite (log-concave)
    dt = 0.37
    K = residual_from_sampled(C, S)
    np.testing.assert_allclose(riemannian_cov_step(C, S, dt),
                               riemannian_cov_step_residual(C, K, dt), atol=1e-10)
    np.testing.assert_allclose(kl_cov_step(C, S, dt),
                               kl_cov_step_residual(C, K, dt), atol=1e-10)


def test_method_scheme_stl_mapping():
    assert method_scheme_stl("riemannian") == ("riemannian", False)
    assert method_scheme_stl("riemannian_stl") == ("riemannian", True)
    assert method_scheme_stl("kl") == ("kl", False)
    assert method_scheme_stl("kl_stl") == ("kl", True)
    assert set(METHODS) == {"riemannian", "riemannian_stl", "kl", "kl_stl"}


# ---------------------------------------------------------------------------
# Log-cosh deterministic optimum
# ---------------------------------------------------------------------------

def test_logcosh_astar_stationarity():
    """The log-cosh optimum is diagonal, m=0, with a tiny stationarity residual."""
    T = build_target("log_cosh", d=5, kappa=100.0, tau=1.0, n_nodes=80)
    m_star, C_star = T.a_star()
    np.testing.assert_allclose(m_star, np.zeros(5), atol=1e-14)
    # C_star is diagonal.
    off = C_star - np.diag(np.diag(C_star))
    assert np.max(np.abs(off)) < 1e-14
    diag = T.a_star_diagnostics()
    assert diag["stationarity_residual_max"] < 1e-9
    # C_star^{-1} = a + tau E[sech^2] >= a (the Gaussian-part lower bound).
    assert np.all(1.0 / np.diag(C_star) >= T.a_diag - 1e-12)


# ---------------------------------------------------------------------------
# Estimator-level core: zero variance at the optimum, exact match
# ---------------------------------------------------------------------------

def test_estimator_variance_zero_at_gaussian_optimum():
    T = build_target("gaussian", d=4, kappa=100.0)
    m_star, C_star = T.a_star()
    bk = ArrayBackend("numpy", "cpu", "float64")
    states, _ = build_states(m_star, C_star)
    opt = [s for s in states if s[0] == "optimum"][0]
    row = evaluate_state(T, *opt, M=50000, seed=11, backend=bk)
    assert row["var_fr_stl"] < 1e-18
    assert row["ratio_fr"] < 1e-12
    # MC baseline variance matches the closed form Tr(C A C A) = d at the optimum.
    np.testing.assert_allclose(row["var_fr_base"], row["var_fr_base_exact"], rtol=0.05)
    np.testing.assert_allclose(row["var_fr_base_exact"], T.d, atol=1e-9)


def test_estimator_variance_far_can_inflate():
    """STL inflates variance when the covariance is under-dispersed (C too small).

    For the Gaussian target the STL fluctuation is ``(C^{-1} - A)(theta - m)``, so
    the FR variance ratio depends only on ``C`` vs ``C_star = A^{-1}``: it is zero
    whenever ``C = C_star`` (any mean), below one when over-dispersed, and *above*
    one when under-dispersed (``C^{-1}`` blows up). This is the regime where STL
    can hurt.
    """
    T = build_target("gaussian", d=4, kappa=100.0)
    m_star, C_star = T.a_star()
    bk = ArrayBackend("numpy", "cpu", "float64")
    states, _ = build_states(m_star, C_star)
    ratios = {}
    for name, m, C, dist in states:
        row = evaluate_state(T, name, m, C, dist, M=40000, seed=5, backend=bk)
        ratios[name] = row["ratio_fr"]
    # Mean-only perturbations keep C = C_star, so STL variance is ~0 there.
    assert ratios["optimum"] < 1e-6
    assert ratios["near"] < 1e-6 and ratios["far"] < 1e-6
    # Covariance misspecification: under-dispersion inflates, over-dispersion reduces.
    assert ratios["underdispersed"] > 1.0
    assert ratios["overdispersed"] < 1.0


# ---------------------------------------------------------------------------
# Algorithm-level core: STL lowers the stochastic noise floor
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("scheme", ["riemannian", "kl"])
def test_algorithm_stl_lowers_noise_floor(scheme):
    T = build_target("gaussian", d=4, kappa=100.0)
    bk = ArrayBackend("numpy", "cpu", "float64")
    seeds = list(range(8))
    floors = {}
    for stl in (False, True):
        method = scheme + ("_stl" if stl else "")
        _, _, tail, _ = simulate_cell(T, method, dt=0.1, n_steps=400,
                                      batch_size=4, seeds=seeds, backend=bk,
                                      cell_seed=2024, n_saved=40)
        floors[stl] = np.median([t["tail_median_gap"] for t in tail])
    assert floors[True] < floors[False]


def test_algorithm_spd_preserved_and_finite():
    T = build_target("log_cosh", d=4, kappa=100.0, tau=1.0)
    bk = ArrayBackend("numpy", "cpu", "float64")
    long_rows, summary, tail, diag = simulate_cell(
        T, "kl_stl", dt=0.1, n_steps=120, batch_size=4, seeds=[0, 1, 2],
        backend=bk, cell_seed=7, n_saved=30)
    assert diag["n_clips"] == 0
    assert all(s["spd_fail"] == 0 for s in summary)
    assert all(np.isfinite(r["energy_gap"]) for r in long_rows)
    assert all(r["min_eig_C"] > 0 for r in long_rows)


@pytest.mark.skipif(not torch_available(), reason="torch not installed")
def test_torch_cpu_backend_matches_exact_gaussian():
    """The torch backend reproduces the Gaussian closed-form variance."""
    T = build_target("gaussian", d=3, kappa=50.0)
    m_star, C_star = T.a_star()
    bk = ArrayBackend("torch", "cpu", "float64")
    states, _ = build_states(m_star, C_star)
    # Under-dispersed state has a nonzero STL variance (C != C_star).
    und = [s for s in states if s[0] == "underdispersed"][0]
    row = evaluate_state(T, *und, M=60000, seed=9, backend=bk)
    np.testing.assert_allclose(row["var_fr_base"], row["var_fr_base_exact"], rtol=0.06)
    np.testing.assert_allclose(row["var_fr_stl"], row["var_fr_stl_exact"], rtol=0.08)
