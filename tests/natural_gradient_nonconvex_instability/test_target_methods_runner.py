import math
import os

import numpy as np

from src.common.io_utils import load_yaml
from src.natural_gradient_nonconvex_instability.methods import (
    kl_next,
    riemannian_next,
    wasserstein_fb_next,
)
from src.natural_gradient_nonconvex_instability.runner import (
    run_all,
    shallow_apply_smoke,
)
from src.natural_gradient_nonconvex_instability.targets import (
    NonconvexLogCoshTarget,
    sech2_from_tanh,
    stable_log_cosh,
)


def test_target_derivatives_and_hessian_bounds():
    target = NonconvexLogCoshTarget(R=7.0, n_nodes=40)
    xs = np.linspace(-70.0, 70.0, 301)
    h = 1e-5
    fd_grad = (target.potential(xs + h) - target.potential(xs - h)) / (2.0 * h)
    np.testing.assert_allclose(target.grad(xs), fd_grad, rtol=1e-6, atol=2e-7)

    fd_hess = (target.grad(xs + h) - target.grad(xs - h)) / (2.0 * h)
    np.testing.assert_allclose(target.hess(xs), fd_hess, rtol=1e-6, atol=2e-7)
    assert np.min(target.hess(xs)) >= -1.0 - 1e-14
    assert np.max(target.hess(xs)) <= 1.0 + 1e-14


def test_stable_logcosh_and_sech2_are_finite():
    z = np.array([-1000.0, -10.0, 0.0, 10.0, 1000.0])
    out = stable_log_cosh(z)
    assert np.all(np.isfinite(out))
    s2 = sech2_from_tanh(z)
    assert np.all(np.isfinite(s2))
    assert np.all((s2 >= 0.0) & (s2 <= 1.0))
    np.testing.assert_allclose(s2[2], 1.0)


def test_quadrature_finite_behavior():
    target = NonconvexLogCoshTarget(R=1000.0, n_nodes=80)
    for c in [1e-4, 1.0, 1e2, 1e5, 1e6]:
        ex = target.expectations(0.0, c)
        assert math.isfinite(ex.V)
        assert math.isfinite(ex.Vp)
        assert math.isfinite(ex.Vpp)
        assert -1.0 - 1e-12 <= ex.Vpp <= 1.0 + 1e-12
        assert math.isfinite(target.energy(0.0, c))


def test_method_formulas_under_constant_negative_curvature():
    c = 1.0
    A = -1.0
    dt = 1.0
    r = riemannian_next(c, A, dt)
    np.testing.assert_allclose(r.c_next, math.exp(2.0), rtol=1e-14)

    eps = 1e-3
    k = kl_next(1.0 - eps, A, dt)
    np.testing.assert_allclose(k.denom, eps, rtol=1e-13)
    np.testing.assert_allclose(k.c_next, 2.0 * (1.0 - eps) / eps, rtol=1e-13)

    eta = 0.9
    w = wasserstein_fb_next(c, A, eta)
    ctilde = (1.0 + eta) ** 2 * c
    expected = 0.5 * (ctilde + 2.0 * eta + math.sqrt(ctilde * (ctilde + 4.0 * eta)))
    np.testing.assert_allclose(w.ctilde, ctilde, rtol=1e-14)
    np.testing.assert_allclose(w.c_next, expected, rtol=1e-14)


def test_kl_non_spd_detection():
    step = kl_next(c=2.0, A=-1.0, dt=1.0)
    assert step.status == "non_spd"
    assert step.denom <= 0.0
    assert math.isnan(step.c_next)


def test_smoke_runner_writes_requested_outputs(tmp_path):
    cfg_path = os.path.join(
        os.path.dirname(__file__), "..", "..", "configs",
        "natural_gradient_nonconvex_instability", "nonconvex_instability.yaml")
    cfg = shallow_apply_smoke(load_yaml(cfg_path))
    run_all(cfg, tmp_path)
    expected = [
        "results_long.csv",
        "summary.csv",
        "kl_pole_summary.csv",
        "wasserstein_bound_summary.csv",
        "clipped_kl_stationarity.csv",
        "clipped_kl_summary.csv",
        "target_metadata.json",
        "run_metadata.json",
    ]
    for name in expected:
        path = tmp_path / name
        assert path.exists()
        assert path.stat().st_size > 0

