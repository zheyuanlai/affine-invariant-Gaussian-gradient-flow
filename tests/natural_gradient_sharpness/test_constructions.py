import math

import numpy as np

from src.natural_gradient_sharpness.bump_train import BumpTrainTarget
from src.natural_gradient_sharpness.local_targets import (
    ChenLogBump,
    RidgeTarget,
    SQRT2,
    ShellTarget,
    kl_jacobian,
)
from src.natural_gradient_sharpness.profiles import IntegratedFlatTop
from src.natural_gradient_sharpness.spiral import SpiralValleyTarget


def test_flat_top_profile_and_bump_train_shadow_exactly():
    profile = IntegratedFlatTop(5001)
    assert abs(profile.I_phi - 1.5) < 2e-6
    target = BumpTrainTarget(8, gh_order=12, profile=profile)
    assert target.min_spacing_over_width > 2.0
    for method in ("riemannian", "kl"):
        state = np.array([target.mean0, target.cov0])
        for j in range(5):
            state = target.step(method, state)
            assert abs(state[0] - target.centers_path[j + 1]) < 1e-8
            assert abs(target.kappa * state[1] - 1.0) < 1e-12


def _finite_difference_scalar_generator(target, eps=2e-5):
    jac = np.empty((3, 3))
    for j in range(3):
        direction = np.zeros(3)
        direction[j] = 1.0
        plus = target.reduced_rhs(0.0, target.state_from_scalar_coordinates(direction, eps))
        minus = target.reduced_rhs(0.0, target.state_from_scalar_coordinates(direction, -eps))
        physical = (plus - minus) / (2.0 * eps)
        jac[:, j] = [
            physical[0], physical[1] / SQRT2,
            math.sqrt(target.m) * physical[2] / SQRT2,
        ]
    return -jac


def test_ridge_and_shell_are_isotropic_and_match_nonlinear_jacobian():
    targets = [
        RidgeTarget(4, normal_order=24, radial_order=12),
        ShellTarget(16, math.log(16), normal_order=36, radial_order=18),
    ]
    for target in targets:
        assert target.isotropy_error < 2e-10
        assert target.gamma > 0.0
        assert target.slow_block == "scalar"
        fd = _finite_difference_scalar_generator(target)
        assert np.max(np.abs(target.blocks["scalar"] - fd)) < 3e-3


def test_chen_bump_exhibits_logarithmic_inverse_gap_growth():
    coarse = ChenLogBump(0.1)
    fine = ChenLogBump(0.01)
    assert fine.kappa > coarse.kappa
    assert 1.0 / fine.gamma > 1.0 / coarse.gamma
    normalized = [
        1.0 / (coarse.gamma * math.log(math.e * coarse.kappa)),
        1.0 / (fine.gamma * math.log(math.e * fine.kappa)),
    ]
    assert max(normalized) / min(normalized) < 1.15


def test_discrete_linearizations_are_contractions_at_safe_local_step():
    target = RidgeTarget(4, normal_order=20, radial_order=10)
    L = target.blocks["scalar"]
    dt = 0.5 / target.Lambda
    assert np.max(np.abs(np.linalg.eigvalsh(np.eye(3) - dt * L))) < 1.0
    assert np.max(np.abs(np.linalg.eigvals(kl_jacobian(L, dt)))) < 1.0


def test_spiral_initialization_is_in_flat_core_and_has_kappa_rate():
    target = SpiralValleyTarget(16, gh_order=6)
    derivative = target.rhs(0.0, target.state0)
    m, dm = target.m0, derivative[:2]
    radial_rate = -float(m @ dm) / float(m @ m)
    assert abs(radial_rate * target.K - 48.0 / 17.0) < 1e-8
    assert target.max_normalized_phase < 0.2
