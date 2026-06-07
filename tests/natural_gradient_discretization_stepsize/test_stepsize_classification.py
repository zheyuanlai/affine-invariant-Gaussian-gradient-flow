"""Tests for the scalar diagnostic and the end-to-end stepsize classification.

Includes the scalar KL closed-form check (target D), the convexity of the
classification (monotone => stable => SPD-feasible), and the central
proof-artifact observation that the KL scheme is empirically stable for
stepsizes far above its theoretical sufficient bound.
"""
import numpy as np
import pytest

from src.natural_gradient_discretization_stepsize.targets import (
    GaussianPosteriorTarget, ScalarGaussianTarget,
)
from src.natural_gradient_discretization_stepsize.methods import (
    kl_cov_step, riemannian_cov_step,
)
from src.natural_gradient_discretization_stepsize.ode_reference import integrate_reference
from src.natural_gradient_discretization_stepsize.optimize_star import compute_star
from src.natural_gradient_discretization_stepsize.metrics import theory_stepsize_bounds
from src.natural_gradient_discretization_stepsize.runner import (
    simulate_run, simulate_scalar_diagnostic,
)


# ---------------------------------------------------------------------------
# Scalar diagnostic (target D)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("C0", [0.01, 0.1, 10.0, 100.0])
@pytest.mark.parametrize("dt", [0.5, 1.0, 2.0, 5.0])
def test_scalar_kl_closed_form(C0, dt):
    """The scalar KL recursion matches (1+dt) C / (1 + dt C) at every step."""
    traj = simulate_scalar_diagnostic("kl", C0, dt, T=10.0, m0=0.0)
    for k in range(1, len(traj)):
        C_prev = traj[k - 1]["C"]
        expected = (1.0 + dt) * C_prev / (1.0 + dt * C_prev)
        assert np.isclose(traj[k]["C"], expected, atol=1e-12)


@pytest.mark.parametrize("C0", [0.01, 0.1, 10.0, 100.0])
def test_scalar_kl_fixed_point_contraction(C0):
    """KL scalar covariance contracts toward C=1: |C_{n+1}-1| < |C_n-1|."""
    dt = 2.0
    traj = simulate_scalar_diagnostic("kl", C0, dt, T=20.0, m0=0.0)
    err = [abs(r["C"] - 1.0) for r in traj]
    assert err[-1] < err[0]
    # closed form: C_{n+1}-1 = (C_n-1)/(1+dt C_n), strictly contracting for dt>0
    for k in range(1, len(traj)):
        assert err[k] <= err[k - 1] + 1e-14


def test_scalar_dense_matches_closed_form():
    """The 1x1 dense covariance steps match the scalar closed forms used by the diagnostic."""
    C = np.array([[10.0]])
    H = np.array([[-1.0]])
    dt = 1.5
    assert np.isclose(kl_cov_step(C, H, dt)[0, 0],
                      ScalarGaussianTarget.kl_cov_closed_form(10.0, dt), atol=1e-12)
    assert np.isclose(riemannian_cov_step(C, H, dt)[0, 0],
                      ScalarGaussianTarget.riemannian_cov_closed_form(10.0, dt), atol=1e-12)


def test_scalar_mean_overshoot():
    """The shared explicit mean step can overshoot for a large stepsize."""
    # m_{n+1} = (1 - dt C_n) m_n; C0=100, dt=5 gives a huge first overshoot.
    traj = simulate_scalar_diagnostic("kl", C0=100.0, dt=5.0, T=10.0, m0=1.0)
    assert abs(traj[1]["m"]) > abs(traj[0]["m"])


# ---------------------------------------------------------------------------
# Classification convexity and the proof-artifact observation
# ---------------------------------------------------------------------------

def _run(method, target, dt, T=15.0):
    m_star, C_star, F_star, _ = compute_star(target)
    ode = integrate_reference(target, T, rtol=1e-8, atol=1e-10, n_eval=120)
    _, summ = simulate_run(method, target, dt, T, F_star, m_star, C_star, ode,
                           max_saved_rows=200)
    return summ


@pytest.mark.parametrize("method", ["riemannian", "kl"])
@pytest.mark.parametrize("dt", [0.01, 0.1, 1.0, 5.0])
def test_classification_nesting(method, dt):
    """monotone => stable => spd_feasible and accurate => spd_feasible."""
    T = GaussianPosteriorTarget(0.1)
    s = _run(method, T, dt)
    if s["monotone"]:
        assert s["stable"]
    if s["stable"]:
        assert s["spd_feasible"]
    if s["accurate"]:
        assert s["spd_feasible"]


def test_kl_stable_far_above_theory_bound():
    """Central observation: KL is empirically stable at dt >> dt_theory_kl.

    For the Gaussian posterior at lambda=0.1 the theoretical KL stepsize is
    2.5e-5, yet a dt=1.0 KL run is SPD-feasible, stable and monotone.
    """
    T = GaussianPosteriorTarget(0.1)
    bounds = theory_stepsize_bounds(T.alpha, T.beta, T.C0)
    assert bounds["dt_theory_kl"] < 1e-3        # very small theoretical bound
    s = _run("kl", T, 1.0)
    assert s["spd_feasible"] and s["stable"] and s["monotone"]
    # dt=1.0 is at least 1000x the theoretical KL stepsize.
    assert 1.0 / bounds["dt_theory_kl"] > 1e3


def test_gaussian_run_converges_at_small_dt():
    """A small-dt Gaussian run drives the energy gap to ~0 and is accurate."""
    T = GaussianPosteriorTarget(1.0)
    s = _run("riemannian", T, 0.01)
    assert s["spd_feasible"] and s["stable"] and s["monotone"] and s["accurate"]
    assert s["gap_final"] < 1e-2
