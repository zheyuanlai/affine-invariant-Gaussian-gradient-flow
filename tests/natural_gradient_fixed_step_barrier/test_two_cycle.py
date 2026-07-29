"""Checks for the exact two-cycle counterexample.

The first four tests verify that the constructed potential is genuinely admissible
--- if any of them fails the counterexample is vacuous, because a potential whose
curvature range is narrower than ``[1, kappa]`` would only be reproducing the
textbook stability limit at a smaller condition number. The last three are the
claims themselves.
"""
import math

import numpy as np
import pytest

from src.natural_gradient_fixed_step_barrier.runner import bw_step, scheme_step, simulate_two_cycle
from src.natural_gradient_fixed_step_barrier.two_cycle import TwoCycleTarget, cycle_constants

KAPPAS = [64.0, 256.0]
GAMMAS = [0.5, 1.0]
SCHEMES = ["riemannian", "kl"]


@pytest.fixture(scope="module")
def targets():
    return {(k, g): TwoCycleTarget(k, g) for k in KAPPAS for g in GAMMAS}


@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_curvature_bounds_are_attained(targets, kappa, gamma):
    """``min V'' = 1`` and ``max V'' = kappa``, both attained: the condition number is real."""
    t = targets[(kappa, gamma)]
    v2min, v2max = t.curvature_range()
    assert v2min == pytest.approx(1.0, rel=1e-9)
    assert v2max == pytest.approx(kappa, rel=1e-9)
    assert v2max / v2min == pytest.approx(kappa, rel=1e-8)


@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_hessian_lipschitz_budget(targets, kappa, gamma):
    """``||V'''||_inf <= LH``, dimension-free and independent of kappa."""
    t = targets[(kappa, gamma)]
    assert t.hessian_lipschitz() <= t.LH * (1.0 + 1e-6)


@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_matched_state_identities(targets, kappa, gamma):
    """``A(M,c) = p`` and ``b(M,c) = r M`` -- the two prescribed Gaussian moments."""
    t = targets[(kappa, gamma)]
    b, A = t.b_A(t.M, t.c)
    assert A == pytest.approx(t.p, rel=1e-10)
    assert b == pytest.approx(t.r * t.M, rel=1e-10)
    assert t.c * A == pytest.approx(1.0, rel=1e-12)      # covariance maps stationary
    assert 1.0 < t.p < kappa and 1.0 < t.r < kappa


@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_potential_is_even_and_optimizer_is_interior(targets, kappa, gamma):
    """``V`` even, ``V'`` odd, and the optimizer is ``(0, 1/kappa)`` with a positive gap at the cycle."""
    t = targets[(kappa, gamma)]
    xs = np.linspace(-2.0 * t.M, 2.0 * t.M, 2001)
    assert np.allclose(t.V0(xs), t.V0(-xs), rtol=1e-12, atol=0)
    assert np.allclose(t.V1(xs), -t.V1(-xs), rtol=1e-12, atol=1e-9)
    assert t.c_star == pytest.approx(1.0 / kappa, rel=1e-12)
    assert t.gap(0.0, t.c_star) == pytest.approx(0.0, abs=1e-9)
    assert t.gap(t.M, t.c) > 0.0                          # the cycle is NOT optimal


@pytest.mark.parametrize("scheme", SCHEMES)
@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_exact_two_cycle_for_both_fisher_rao_schemes(targets, kappa, gamma, scheme):
    """One step maps ``(M, c) -> (-M, c)`` and the next maps back, for both schemes."""
    t = targets[(kappa, gamma)]
    b, A = t.b_A(t.M, t.c)
    m1, c1 = scheme_step(scheme, t.M, t.c, b, A, gamma)
    assert m1 == pytest.approx(-t.M, rel=1e-10)
    assert c1 == pytest.approx(t.c, rel=1e-12)
    b1, A1 = t.b_A(m1, c1)
    m2, c2 = scheme_step(scheme, m1, c1, b1, A1, gamma)
    assert m2 == pytest.approx(t.M, rel=1e-10)
    assert c2 == pytest.approx(t.c, rel=1e-12)

    _, s = simulate_two_cycle(t, scheme, gamma, max_steps=500)
    assert s["converged"] == 0 and s["status"] == "cycling"
    assert s["sign_flips"] == 500                      # flips every single step
    assert s["min_rel_gap"] == pytest.approx(1.0, rel=1e-9)


@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_bures_wasserstein_converges_on_the_same_target(targets, kappa, gamma):
    """At its certified ``eta <= 1/beta`` the BW scheme converges where Fisher--Rao cycles."""
    t = targets[(kappa, gamma)]
    for mult in (0.25, 1.0):
        _, s = simulate_two_cycle(t, "bures_wasserstein", mult / kappa, max_steps=2000)
        assert s["converged"] == 1, (kappa, gamma, mult, s["status"])
        assert s["sign_flips"] == 0                    # monotone, no oscillation
    # The BW mean multiplier is 1 - eta * secant, and the secant is at most beta,
    # so it stays in [0, 1) for every eta <= 1/beta -- no covariance factor appears.
    b, _ = t.b_A(t.M, t.c)
    assert 0.0 <= 1.0 - (1.0 / kappa) * b / t.M < 1.0


def test_cycle_constants_require_gamma_above_two_over_kappa():
    """Below ``2/kappa`` the construction is infeasible: that is where the bump train takes over."""
    q, p, r = cycle_constants(256.0, 0.5)
    assert r / p == pytest.approx(2.0 / 0.5)
    with pytest.raises(ValueError):
        cycle_constants(256.0, 1.0 / 256.0)            # gamma < 2/kappa
    with pytest.raises(ValueError):
        cycle_constants(256.0, 2.5)                    # gamma > 2


def test_bw_step_has_no_covariance_preconditioner():
    """``bw_step`` mean update is ``m - eta b``; the Fisher--Rao one is ``m - dt c b``."""
    m, c, b, A, h = 3.0, 0.01, 5.0, 2.0, 0.1
    m_bw, _ = bw_step(m, c, b, A, h)
    m_fr, _ = scheme_step("riemannian", m, c, b, A, h)
    assert m_bw == pytest.approx(m - h * b)
    assert m_fr == pytest.approx(m - h * c * b)
