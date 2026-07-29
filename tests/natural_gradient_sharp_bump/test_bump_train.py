"""Checks that the bump train reproduces Appendix C, and that the arms separate.

Every assertion below is a statement the manuscript proves, so a failure means the
implementation has drifted from ``thm:sharp-disc`` rather than that the theory is
wrong. The last two tests are the experiment's actual claims.
"""
import math

import numpy as np
import pytest

from src.natural_gradient_sharp_bump import BumpTrain, simulate, step_size

GAMMA = 0.5
KAPPAS = [8.0, 32.0, 128.0]


def manuscript_train(kappa, **kw):
    """The Appendix C family: per-step mean gain ``s = gamma/kappa^2``."""
    return BumpTrain(kappa, GAMMA / kappa ** 2, **kw)


@pytest.mark.parametrize("kappa", KAPPAS)
def test_regularity_and_condition_number(kappa):
    """lem:bump-regularity: ``1 <= V'' <= kappa`` (both attained), ``|V'''| <= LH``."""
    t = manuscript_train(kappa)
    theta = np.linspace(t.x_N - 2 * t.w, t.x0 + 2 * t.w, 200_001)
    v2 = t.V2(theta)
    assert v2.min() == pytest.approx(1.0, abs=1e-9)
    assert v2.max() == pytest.approx(kappa, rel=1e-9)
    assert float(t.V2(0.0)) == pytest.approx(1.0, abs=1e-12)
    v3 = np.gradient(v2, theta)
    assert np.abs(v3).max() <= t.LH * (1.0 + 1e-6)


@pytest.mark.parametrize("kappa", KAPPAS)
def test_center_geometry(kappa):
    """lem:bump-centers: ``x_N >= 3/4 x_0`` and spacing ``>= 3 Y kappa^2 / 4``."""
    t = manuscript_train(kappa)
    assert t.x_N >= 0.75 * t.x0
    assert t.min_spacing >= 0.75 * t.Y * kappa ** 2
    assert t.min_spacing > 2.0 * t.w        # supports pairwise disjoint


@pytest.mark.parametrize("kappa", KAPPAS)
def test_gaussian_averages_at_centers(kappa):
    """lem:bump-averages: ``A = kappa`` and ``b = x_j + H (N - j + 1/2)`` at every center."""
    t = manuscript_train(kappa)
    for j in (0, t.N // 2, t.N):
        b, A = t.b_A(t.centers[j], t.c_kappa)
        assert A == pytest.approx(kappa, rel=1e-10)
        assert b == pytest.approx(t.centers[j] + t.H * (t.N - j + 0.5), rel=1e-10)


@pytest.mark.parametrize("kappa", KAPPAS)
def test_optimizer_is_standard_gaussian(kappa):
    """``V(theta) = theta^2/2`` near the origin, so ``a_star = (0, 1)`` and ``E_star = 1/2``."""
    t = manuscript_train(kappa)
    theta = np.linspace(-6.0, 6.0, 1001)
    assert np.allclose(t.V0(theta), 0.5 * theta ** 2, rtol=0, atol=1e-9)
    assert t.energy_star == pytest.approx(0.5, abs=1e-10)
    assert t.gap(0.0, 1.0) == pytest.approx(0.0, abs=1e-10)
    # a_star is a minimum: perturbing either coordinate raises the objective.
    assert t.gap(0.3, 1.0) > 0.0 and t.gap(0.0, 1.6) > 0.0 and t.gap(0.0, 0.6) > 0.0


@pytest.mark.parametrize("scheme", ["riemannian", "kl"])
@pytest.mark.parametrize("kappa", KAPPAS)
def test_theory_arm_is_blocked_for_N_train_steps(scheme, kappa):
    """thm:sharp-disc: at ``dt = gamma/kappa`` the gap stays macroscopic for ``N`` steps.

    Also checks lem:bump-shadowing -- the iterates track the ideal centers ``x_j``
    and hold ``c_n = 1/kappa`` -- which is what makes the lower bound bite.
    """
    t = manuscript_train(kappa)
    _, s = simulate(t, scheme, "theory", GAMMA, max_steps=t.N, tol_rel=1e-6)
    assert s["n_steps"] == t.N                  # never reached the tolerance
    assert s["n_half"] == -1                    # gap never even halved
    assert s["final_rel_gap"] >= 0.5
    assert s["shadow_mean_err_over_w"] < 1e-6   # |m_n - x_n| << w
    assert s["shadow_cov_err"] < 1e-6           # |kappa c_n - 1| small
    assert s["cA_max"] == pytest.approx(1.0, rel=1e-9)


@pytest.mark.parametrize("scheme", ["riemannian", "kl"])
def test_order_one_step_beats_the_certified_step(scheme):
    """The experiment's claim: on the retuned train the ``O(1)`` step costs ``O(kappa)``.

    The certified step ``gamma/kappa`` costs ``~kappa^2``, so the ratio must grow with
    ``kappa``; both arms stay monotone, so the ``O(1)`` step is not merely faster by
    being unstable.
    """
    ratios = []
    for kappa in (32.0, 128.0):
        theory = manuscript_train(kappa)
        _, s_theory = simulate(theory, scheme, "theory", GAMMA,
                               max_steps=8 * theory.N, tol_rel=1e-6)
        retuned = BumpTrain(kappa, GAMMA / kappa)        # rebuilt for dt = gamma
        _, s_const = simulate(retuned, scheme, "const", GAMMA,
                              max_steps=8 * max(retuned.N, 1000), tol_rel=1e-6)
        assert s_theory["status"] == s_const["status"] == "ok"
        assert s_theory["monotone"] == 1 and s_const["monotone"] == 1
        assert s_const["n_half"] > 0 and s_theory["n_half"] > 0
        ratios.append(s_theory["n_half"] / s_const["n_half"])
    # kappa 32 -> 128 is a factor 4; an extra factor of kappa gives >= 4x in the ratio.
    assert ratios[1] > 3.5 * ratios[0]


def test_relative_curvature_rule_only_shrinks_the_step():
    """``relcurv`` is a stability test: ``dt <= gamma``, with equality when ``cA <= 1``."""
    assert step_size("relcurv", GAMMA, 100.0, c=1.0, A=1.0) == pytest.approx(GAMMA)
    assert step_size("relcurv", GAMMA, 100.0, c=0.01, A=1.0) == pytest.approx(GAMMA)
    assert step_size("relcurv", GAMMA, 100.0, c=1.0, A=4.0) == pytest.approx(GAMMA / 4.0)


def test_retuned_train_matches_manuscript_at_the_certified_step():
    """``s = dt/kappa`` with ``dt = gamma/kappa`` reproduces the Appendix C family."""
    a = manuscript_train(64.0)
    b = BumpTrain(64.0, (GAMMA / 64.0) / 64.0)
    assert a.N == b.N and a.x0 == pytest.approx(b.x0)
    assert np.allclose(a.centers, b.centers)
    assert b.N == math.floor(b.T / b.s)
