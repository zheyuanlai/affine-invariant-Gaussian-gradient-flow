import math

import numpy as np

from src.natural_gradient_nonconvex_instability.methods import (
    clip_smoothness_constant,
    clipped_kl_next,
    gaussian_kl_divergence,
    theorem_safe_dt,
)
from src.natural_gradient_nonconvex_instability.runner import (
    CLIPPED_THEOREM_TOL,
    simulate_clipped_kl,
)
from src.natural_gradient_nonconvex_instability.targets import (
    NonconvexLogCoshTarget,
)


def test_gaussian_kl_diagnostic_nonneg_and_zero_iff_equal():
    # Zero exactly when covariances coincide.
    for c in [0.25, 1.0, 4.0, 100.0]:
        assert gaussian_kl_divergence(c, c) == 0.0
    # Strictly positive and nonnegative otherwise.
    rng = np.random.default_rng(0)
    for _ in range(200):
        c0 = float(rng.uniform(1e-3, 1e3))
        c1 = float(rng.uniform(1e-3, 1e3))
        d = gaussian_kl_divergence(c0, c1)
        assert d >= -1e-15
        if not math.isclose(c0, c1, rel_tol=1e-12):
            assert d > 0.0


def test_clipped_kl_update_always_in_feasible_interval():
    lam_m, lam_p = 0.5, 2.0
    dt = theorem_safe_dt(beta=1.0, lambda_minus=lam_m, lambda_plus=lam_p, safety=0.9)
    rng = np.random.default_rng(1)
    for _ in range(500):
        c = float(rng.uniform(lam_m, lam_p))
        A = float(rng.uniform(-1.0, 1.0))  # |V''| <= 1 => |A| <= 1
        step = clipped_kl_next(c, A, dt, lam_m, lam_p)
        assert step.status == "ok"
        assert lam_m - 1e-12 <= step.c_next <= lam_p + 1e-12


def test_theorem_safe_stepsize_keeps_denominator_positive():
    lam_m, lam_p, beta = 0.5, 2.0, 1.0
    dt = theorem_safe_dt(beta=beta, lambda_minus=lam_m, lambda_plus=lam_p, safety=0.9)
    L_clip = clip_smoothness_constant(beta, lam_m, lam_p)
    assert math.isclose(dt, 0.9 / L_clip, rel_tol=1e-15)
    # denom = 1 + dt c A >= 1 - dt beta lambda_plus, which must stay positive.
    assert 1.0 - dt * beta * lam_p > 0.0
    # Exhaustive over feasible (c, A): denominator stays strictly positive.
    for c in np.linspace(lam_m, lam_p, 50):
        for A in np.linspace(-beta, beta, 50):
            step = clipped_kl_next(float(c), float(A), dt, lam_m, lam_p)
            assert step.denom > 0.0
            assert step.status == "ok"


def test_smoke_clipped_run_satisfies_theorem_envelope():
    target = NonconvexLogCoshTarget(R=1000.0, n_nodes=80)
    rows, summary = simulate_clipped_kl(
        target, lambda_minus=0.5, lambda_plus=2.0, beta=1.0,
        dt_safety=0.9, c0=1.0, num_steps=60)
    assert summary["dt_rule"] == "theorem_safe"
    assert summary["dt_times_L_clip"] < 1.0
    assert summary["denominator_positive"] is True
    assert summary["max_violation"] <= 1e-6
    assert summary["theorem_check_pass"] is True
    # Every recorded prefix obeys D_min(N) <= B_N up to tolerance.
    for r in rows:
        assert r["running_min_D"] <= r["prefix_bound"] + 1e-9
    # The covariance is pinned inside the feasible interval.
    assert 0.5 - 1e-12 <= summary["min_c"] <= summary["max_c"] <= 2.0 + 1e-12
    assert summary["max_violation"] <= CLIPPED_THEOREM_TOL or summary["theorem_check_pass"]


def test_riemannian_scale_stepsize_is_outside_theorem_but_envelope_holds():
    # dt = 1 / (beta * lambda_plus) exceeds 1 / L_clip, so the Theorem 2.18
    # condition fails (dt * L_clip > 1), yet on this target the stationarity
    # envelope is still satisfied empirically and the run stays SPD-feasible.
    lam_m, lam_p, beta = 0.5, 2.0, 1.0
    dt = 1.0 / (beta * lam_p)
    L = clip_smoothness_constant(beta, lam_m, lam_p)
    target = NonconvexLogCoshTarget(R=1000.0, n_nodes=80)
    rows, summary = simulate_clipped_kl(
        target, lambda_minus=lam_m, lambda_plus=lam_p, beta=beta,
        c0=1.0, num_steps=40, dt=dt, dt_rule="riemannian_scale")
    assert summary["dt_rule"] == "riemannian_scale"
    assert math.isclose(summary["dt"], dt, rel_tol=1e-15)
    assert summary["dt_times_L_clip"] == dt * L > 1.0  # outside the theorem
    assert summary["denominator_positive"] is True     # clip masks the pole here
    assert 0.5 - 1e-12 <= summary["min_c"] <= summary["max_c"] <= 2.0 + 1e-12
    # Envelope still holds and the covariance reaches the upper bound at once.
    assert summary["max_violation"] <= 1e-9
    assert summary["theorem_check_pass"] is True
    assert summary["first_upper_clip_step"] == 0
