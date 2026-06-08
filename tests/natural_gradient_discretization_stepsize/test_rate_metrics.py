"""Tests for the theoretical-rate benchmark utilities and runner.

Basename is globally unique across tests/ (flat default import mode).
"""
import math

import numpy as np
import pytest

from src.natural_gradient_discretization_stepsize.targets import (
    GaussianPosteriorTarget, build_target,
)
from src.natural_gradient_discretization_stepsize import rate_metrics as rm
from src.natural_gradient_discretization_stepsize.rate_runner import simulate_rate_run
from src.natural_gradient_discretization_stepsize.optimize_star import compute_star


# Spectral bounds reused across cases (Gaussian lam=0.1 style: ill-conditioned).
def _bounds(lam=0.1):
    T = GaussianPosteriorTarget(lam)
    lam_min, lam_max = rm.spectral_bounds(T.alpha, T.beta, T.C0)
    return T.alpha, T.beta, lam_min, lam_max


C_GRID = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0]


# ---------------------------------------------------------------------------
# Contraction factors are valid probabilities in (0, 1) on the c-grid
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("lam", [0.01, 0.1, 1.0])
def test_q_riem_theory_in_unit_interval(lam):
    """q_riem_theory in (0, 1) for c in (0, 1] of dt_ref."""
    alpha, beta, lam_min, lam_max = _bounds(lam)
    dtr = rm.dt_ref(beta, lam_max)
    for c in C_GRID:
        q = rm.q_riem_theory(c * dtr, alpha, beta, lam_min, lam_max)
        assert 0.0 < q < 1.0, (lam, c, q)


@pytest.mark.parametrize("lam", [0.01, 0.1, 1.0])
def test_q_kl_formula_in_unit_interval(lam):
    """q_kl_formula in (0, 1) on the same Riemannian-scale grid."""
    alpha, beta, lam_min, lam_max = _bounds(lam)
    dtr = rm.dt_ref(beta, lam_max)
    for c in C_GRID:
        q = rm.q_kl_formula(c * dtr, alpha, beta, lam_min, lam_max)
        assert 0.0 < q < 1.0, (lam, c, q)


def test_q_riem_boundary_at_dt_ref():
    """At c=1 (dt=dt_ref) the Riemannian bracket (2 - beta lam_max dt) = 1."""
    alpha, beta, lam_min, lam_max = _bounds(0.1)
    dtr = rm.dt_ref(beta, lam_max)
    # beta*lam_max*dt_ref = 1 exactly, so bracket = 1 and q = 1 - alpha*lam_min*dt.
    expected = 1.0 - alpha * lam_min * dtr
    assert np.isclose(rm.q_riem_theory(dtr, alpha, beta, lam_min, lam_max), expected)


def test_kl_more_conservative_than_riem_rate():
    """At matched Riemannian-scale dt the KL formula rate is <= Riemannian rate."""
    alpha, beta, lam_min, lam_max = _bounds(0.1)
    dtr = rm.dt_ref(beta, lam_max)
    for c in C_GRID:
        dt = c * dtr
        r_riem = rm.per_unit_rate(rm.q_riem_theory(dt, alpha, beta, lam_min, lam_max), dt)
        r_kl = rm.per_unit_rate(rm.q_kl_formula(dt, alpha, beta, lam_min, lam_max), dt)
        assert r_kl <= r_riem + 1e-12, (c, r_kl, r_riem)


def test_per_unit_rate_invalid():
    """per_unit_rate returns NaN for non-contractive q."""
    assert math.isnan(rm.per_unit_rate(1.5, 0.1))
    assert math.isnan(rm.per_unit_rate(0.0, 0.1))
    assert math.isnan(rm.per_unit_rate(-0.2, 0.1))
