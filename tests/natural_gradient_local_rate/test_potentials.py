"""Tests for the natural-gradient potential families."""
import numpy as np
import pytest

from src.natural_gradient_local_rate.potentials import (
    build_potential, FEATURE_FAMILIES, ALL_FAMILIES,
    GaussianPotential, SeparablePotential, CenteredPotential,
)


def _fd_grad(pot, theta, eps=1e-6):
    n = theta.size
    g = np.zeros(n)
    for i in range(n):
        e = np.zeros(n); e[i] = eps
        g[i] = (pot.value(theta + e) - pot.value(theta - e)) / (2 * eps)
    return g


def _fd_hess(pot, theta, eps=1e-6):
    n = theta.size
    H = np.zeros((n, n))
    for i in range(n):
        e = np.zeros(n); e[i] = eps
        H[:, i] = (pot.grad(theta + e) - pot.grad(theta - e)) / (2 * eps)
    return 0.5 * (H + H.T)


@pytest.mark.parametrize("family", ALL_FAMILIES)
def test_grad_hess_finite_difference(family):
    pot = build_potential(family, N_theta=5, kappa_target=5.0, seed=1,
                          centering_samples=4096)
    rng = np.random.default_rng(0)
    theta = rng.standard_normal(5)
    assert np.allclose(pot.grad(theta), _fd_grad(pot, theta), atol=1e-6)
    assert np.allclose(pot.hess(theta), _fd_hess(pot, theta), atol=1e-6)
    # Hessian symmetric
    assert np.allclose(pot.hess(theta), pot.hess(theta).T, atol=1e-12)


@pytest.mark.parametrize("family", ALL_FAMILIES)
def test_batch_matches_single(family):
    pot = build_potential(family, N_theta=4, kappa_target=3.0, seed=2,
                          centering_samples=2048)
    rng = np.random.default_rng(3)
    Theta = rng.standard_normal((6, 4))
    bv, bg, bh = pot.batch_value(Theta), pot.batch_grad(Theta), pot.batch_hess(Theta)
    assert bv.shape == (6,)
    assert bg.shape == (6, 4)
    assert bh.shape == (6, 4, 4)
    for j in range(6):
        assert pot.value(Theta[j]) == pytest.approx(bv[j], abs=1e-10)
        assert np.allclose(pot.grad(Theta[j]), bg[j], atol=1e-10)
        assert np.allclose(pot.hess(Theta[j]), bh[j], atol=1e-10)


def test_gaussian_exact():
    pot = GaussianPotential(6)
    rng = np.random.default_rng(0)
    theta = rng.standard_normal(6)
    assert np.allclose(pot.grad(theta), theta)
    assert np.allclose(pot.hess(theta), np.eye(6))
    md = pot.metadata()
    assert md["kappa_target"] == 1.0 and md["rho"] == 0.0


@pytest.mark.parametrize("family", FEATURE_FAMILIES)
def test_centering_conditions(family):
    # E[grad V] ~ 0 and E[Hess V] = I on the centering bank (machine-zero).
    pot = build_potential(family, N_theta=5, kappa_target=4.0, seed=4,
                          centering_samples=8192)
    md = pot.metadata()
    assert md["norm_mean_grad"] < 1e-10
    assert md["norm_mean_hess_minus_I"] < 1e-10


@pytest.mark.parametrize("family", FEATURE_FAMILIES)
@pytest.mark.parametrize("kappa", [2.0, 5.0, 20.0])
def test_nominal_kappa(family, kappa):
    # By construction beta_target / alpha_target == kappa_target exactly.
    pot = build_potential(family, N_theta=4, kappa_target=kappa, seed=0,
                          centering_samples=4096)
    md = pot.metadata()
    assert md["beta_target"] / md["alpha_target"] == pytest.approx(kappa, rel=1e-9)
    assert md["alpha_target"] > 0.0
    assert md["rho"] > 0.0


def test_separable_uses_deterministic_LA():
    pot = build_potential("separable", N_theta=5, kappa_target=5.0, seed=0,
                          centering_samples=4096)
    assert pot.metadata()["L_A_is_empirical"] is False


@pytest.mark.parametrize("family", ["additive_index", "random_feature", "radial_tail",
                                    "product_feature"])
def test_coupled_families_empirical_LA(family):
    pot = build_potential(family, N_theta=5, kappa_target=5.0, seed=0,
                          centering_samples=4096, safety_factor=2.0)
    assert pot.metadata()["L_A_is_empirical"] is True


def test_metadata_required_keys():
    pot = build_potential("random_feature", N_theta=4, kappa_target=10.0, seed=0,
                          centering_samples=2048)
    md = pot.metadata()
    for key in ["kappa_target", "alpha_target", "beta_target", "rho",
                "L_A", "estimated_L_A", "L_A_is_empirical",
                "norm_mean_grad", "norm_mean_hess_minus_I",
                "empirical_min_hess_eig", "empirical_max_hess_eig", "N_theta"]:
        assert key in md, f"missing metadata key {key}"


def test_empirical_hess_eig_within_strong_convexity():
    # Empirical Hessian stays positive (strongly convex) on the centering bank.
    pot = build_potential("separable", N_theta=5, kappa_target=10.0, seed=7,
                          centering_samples=8192)
    md = pot.metadata()
    assert md["empirical_min_hess_eig"] > 0.0
    assert md["empirical_max_hess_eig"] >= md["empirical_min_hess_eig"]
