"""Tests for the centralized theory constants (improved KL / projected-KL proof).

Basename is globally unique across tests/ (flat default import mode).
"""
import numpy as np
import pytest

from src.common import theory_constants as tc


# Representative smooth log-concave (alpha, beta, lambda0_min, lambda0_max) cases.
SMOOTH_CASES = [
    (0.1, 1.0, 0.5, 2.0),
    (0.01, 1.0, 0.5, 2.0),
    (1.0, 1.0, 0.5, 2.0),
    (0.05, 1.15, 1.0, 1.0),
    (0.1, 2.0, 1e-3, 5.0),
]


@pytest.mark.parametrize("alpha,beta,l0min,l0max", SMOOTH_CASES)
def test_kl_and_riem_share_theorem_safe_stepsize(alpha, beta, l0min, l0max):
    """For smooth log-concave targets dt_kl_theory == dt_riem_theory."""
    riem = tc.riemannian_theory_constants(alpha, beta, l0min, l0max)
    kl = tc.kl_theory_constants(alpha, beta, l0min, l0max)
    assert riem["L_Riem"] == kl["L_KL"]
    assert riem["dt_Riem_theory"] == kl["dt_KL_theory"]
    # Both equal 1/(beta*lambda_max).
    _, lam_max = tc.natural_gradient_spectral_bounds(l0min, l0max, alpha, beta)
    assert np.isclose(kl["dt_KL_theory"], 1.0 / (beta * lam_max))


def test_kl_stepsize_has_no_lambda_min_cubic_dependence():
    """dt_KL is invariant to lambda0_min (the old lambda_min^3 penalty is gone)."""
    alpha, beta, l0max = 0.1, 1.0, 2.0
    dts = [tc.kl_theory_constants(alpha, beta, l0min, l0max)["dt_KL_theory"]
           for l0min in (1e-6, 1e-3, 0.1, 0.5, 1.0)]
    assert max(dts) - min(dts) == 0.0
    # And it equals the Riemannian stepsize, which also ignores lambda_min.
    riem = tc.riemannian_theory_constants(alpha, beta, 1e-6, l0max)
    assert dts[0] == riem["dt_Riem_theory"]


@pytest.mark.parametrize("alpha,beta,l0min,l0max", SMOOTH_CASES)
def test_q_kl_in_unit_interval_on_theorem_safe_grid(alpha, beta, l0min, l0max):
    """q_KL(dt) in (0, 1) for theorem-safe grid points dt = c/(beta*lambda_max)."""
    lam_min, lam_max = tc.natural_gradient_spectral_bounds(l0min, l0max, alpha, beta)
    dt_theory = 1.0 / (beta * lam_max)
    for c in (0.05, 0.1, 0.25, 0.5, 0.75, 1.0):
        dt = c * dt_theory
        q = tc.q_kl(dt, alpha, beta, lam_min, lam_max)
        assert 0.0 < q < 1.0, (alpha, beta, c, q)
        qr = tc.q_riem(dt, alpha, beta, lam_min, lam_max)
        assert 0.0 < qr < 1.0, (alpha, beta, c, qr)


def test_projected_kl_constant_is_2_beta_lambda_plus_independent_of_lambda_minus():
    """L_clip = 2*beta*lambda_plus and does not depend on lambda_minus."""
    for beta in (0.5, 1.0, 2.0):
        for lam_plus in (1.0, 2.0, 5.0):
            consts = tc.projected_kl_theory_constants(beta, lam_plus)
            assert consts["L_clip"] == 2.0 * beta * lam_plus
            assert np.isclose(consts["dt_projected_KL_theory"],
                              1.0 / (2.0 * beta * lam_plus))
    # projected_kl_theory_constants has no lambda_minus argument at all.
    import inspect
    params = inspect.signature(tc.projected_kl_theory_constants).parameters
    assert "lambda_minus" not in params and set(params) == {"beta", "lambda_plus"}


def test_projected_kl_theorem_safe_dt_scales_like_inverse_lambda_plus():
    """dt_projected_KL_theory scales like 1/lambda_plus at fixed beta."""
    beta = 1.0
    lam_plus_values = [1.0, 2.0, 4.0, 8.0]
    dts = [tc.projected_kl_theory_constants(beta, lp)["dt_projected_KL_theory"]
           for lp in lam_plus_values]
    # dt * lambda_plus is constant (= 1/(2 beta)).
    products = [dt * lp for dt, lp in zip(dts, lam_plus_values)]
    assert np.allclose(products, products[0])
    assert np.isclose(products[0], 1.0 / (2.0 * beta))
    # Doubling lambda_plus halves the theorem-safe stepsize.
    assert np.isclose(dts[1], dts[0] / 2.0)
    assert np.isclose(dts[2], dts[0] / 4.0)


def test_deprecated_formulas_recover_old_values():
    """Deprecated helpers still reproduce the obsolete (historical) constants."""
    # Old KL cubic penalty.
    assert tc.deprecated_old_kl_stepsize_factor(0.5, 10.0) == 10.0 ** 3 / (2.0 * 0.5 ** 3)
    # Old clipped relative-smoothness constant.
    assert tc.deprecated_old_projected_kl_smoothness_constant(1.0, 0.5, 2.0) == \
        1.0 * max(2.0, 2.0 ** 4 / 0.5 ** 3)
