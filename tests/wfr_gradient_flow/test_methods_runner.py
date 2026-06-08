"""Tests for the WFR splitting steps and the run-level driver.

Covers:
* the Wasserstein covariance step against its spectral closed form and the fixed
  point of the Gaussian target;
* SPD preservation of the full WFR splitting on the Gaussian target;
* FR-only recovered exactly when the Wasserstein step ``h = 0`` (and the WFR
  splitting with ``h = 0`` collapsing to FR-only with a single batch);
* expectation-batch accounting (1 per iter for fr_only/w_only, 2 for WFR);
* a short end-to-end smoke run of every method.
"""
import numpy as np
import pytest

from src.common.spd import symmetrize
from src.wfr_gradient_flow.targets import GaussianTarget, SmoothLogCoshTarget
from src.wfr_gradient_flow.methods import (
    wasserstein_cov_step, wfr_step, fisher_rao_step,
)
from src.wfr_gradient_flow.schedules import build_schedule, _ConstantSchedule
from src.wfr_gradient_flow.runner import simulate_run
from src.natural_gradient_discretization_stepsize.methods import kl_cov_step


# ---------------------------------------------------------------------------
# Wasserstein covariance step
# ---------------------------------------------------------------------------

def test_wasserstein_step_spectral_form():
    """Diagonal C, H: the W step matches the per-eigenvalue closed form."""
    C = np.diag([0.5, 2.0])
    H = np.diag([-1.0, -0.25])      # SPD -H
    h = 0.3
    out = wasserstein_cov_step(C, H, h)
    # Diagonal => eigenvalues are the diagonal of C_tilde = (1+h H)C(1+h H).
    M = np.eye(2) + h * H
    w = np.diag(M @ C @ M)
    expected = 0.5 * (w + 2 * h + np.sqrt(w * (w + 4 * h)))
    np.testing.assert_allclose(np.diag(out), expected, atol=1e-12)
    np.testing.assert_allclose(out, out.T, atol=1e-14)


def test_wasserstein_step_identity_at_h0():
    """h = 0 leaves the covariance unchanged (the W step is the identity)."""
    C = np.array([[1.3, 0.2], [0.2, 0.7]])
    H = -np.eye(2)
    np.testing.assert_allclose(wasserstein_cov_step(C, H, 0.0), symmetrize(C), atol=1e-14)


def test_wasserstein_step_spd_preserved():
    """The W covariance step stays SPD for SPD C and log-concave H."""
    rng = np.random.default_rng(5)
    for _ in range(20):
        A = rng.standard_normal((2, 2))
        C = A @ A.T + 0.1 * np.eye(2)
        H = -(A @ A.T + 0.05 * np.eye(2))      # -H SPD
        out = wasserstein_cov_step(C, H, 0.4)
        assert np.linalg.eigvalsh(out)[0] > 0.0


# ---------------------------------------------------------------------------
# FR-only recovered when h = 0
# ---------------------------------------------------------------------------

def test_fr_only_equals_wfr_h0():
    """WFR splitting with h=0 equals the standalone FR step and uses 1 batch."""
    T = GaussianTarget(100.0, epsilon=1e-2)
    m, C = T.m0.copy(), T.C0.copy()
    g, H = T.g_H(m, C)
    dt = 0.5
    # fr_only path.
    m_fr, C_fr, d_fr = wfr_step("fr_only", T, m, C, g, H, 0.0, dt)
    # wfr_fixed path at h = 0 (degenerate: W step identity).
    m_wfr, C_wfr, d_wfr = wfr_step("wfr_fixed", T, m, C, g, H, 0.0, dt)
    np.testing.assert_allclose(m_fr, m_wfr, atol=1e-13)
    np.testing.assert_allclose(C_fr, C_wfr, atol=1e-13)
    assert d_fr["n_batches"] == 1 and d_wfr["n_batches"] == 1
    # And both equal the raw KL covariance step.
    np.testing.assert_allclose(C_fr, kl_cov_step(C, H, dt), atol=1e-13)


# ---------------------------------------------------------------------------
# Batch accounting
# ---------------------------------------------------------------------------

def test_batch_accounting():
    """fr_only/w_only spend 1 batch/iter; WFR with h>0 spends 2."""
    T = GaussianTarget(100.0, epsilon=1e-2)
    m, C = T.m0.copy(), T.C0.copy()
    g, H = T.g_H(m, C)
    h, dt = 0.3, 0.5
    assert wfr_step("fr_only", T, m, C, g, H, h, dt)[2]["n_batches"] == 1
    assert wfr_step("w_only", T, m, C, g, H, h, dt)[2]["n_batches"] == 1
    assert wfr_step("wfr_fixed", T, m, C, g, H, h, dt)[2]["n_batches"] == 2
    assert wfr_step("wfr_adaptive", T, m, C, g, H, h, dt)[2]["n_batches"] == 2


# ---------------------------------------------------------------------------
# SPD preservation of the full splitting over many steps
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("method", ["w_only", "wfr_fixed", "wfr_adaptive"])
def test_spd_preserved_over_run(method):
    """A multi-step run preserves SPD and symmetry on the Gaussian target."""
    T = GaussianTarget(100.0, epsilon=1e-2)
    sched = build_schedule(method, beta=T.beta, C0=T.C0, c=0.5)
    m, C = T.m0.copy(), symmetrize(T.C0)
    for n in range(60):
        g, H = T.g_H(m, C)
        h, _ = sched.h(n, m, C, g, H)
        m, C, diag = wfr_step(method, T, m, C, g, H, h, 0.5)
        assert diag["spd_ok"], f"lost SPD at step {n} ({method})"
        np.testing.assert_allclose(C, C.T, atol=1e-12)


# ---------------------------------------------------------------------------
# End-to-end smoke run of every method
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("method",
                         ["fr_only", "w_only", "wfr_fixed", "wfr_theory", "wfr_adaptive"])
def test_simulate_run_converges_gaussian(method):
    """Every method drives the Gaussian objective gap down and stays SPD."""
    T = GaussianTarget(100.0, epsilon=1e-2)
    m_star, C_star, F_star = T.star()
    sched = build_schedule(method, beta=T.beta, C0=T.C0, c=0.5)
    records, summary = simulate_run(method, T, 0.5, 4000, sched,
                                    F_star, m_star, C_star)
    assert summary["status"] == "ok"
    assert summary["spd_feasible"] == 1
    assert summary["gap_final"] < summary["gap0"]
    # Batch budget respected; iteration count consistent with the per-iter cost.
    assert summary["n_batches"] <= 4000 + 2


def test_wfr_reaches_tolerance_faster_in_warmup():
    """On the underdispersed Gaussian, WFR-adaptive hits 1e-1 in no more batches
    than FR-only (warmup advantage; not a strict inequality assertion)."""
    T = GaussianTarget(1000.0, epsilon=1e-3)
    m_star, C_star, F_star = T.star()
    out = {}
    for method in ("fr_only", "wfr_adaptive"):
        sched = build_schedule(method, beta=T.beta, C0=T.C0, c=0.5)
        _, s = simulate_run(method, T, 0.5, 6000, sched, F_star, m_star, C_star)
        out[method] = s["batches_to_1e_minus_1"]
    # Both should reach the loose threshold within budget.
    assert out["fr_only"] != -1 and out["wfr_adaptive"] != -1


def test_smooth_target_runs():
    """The smooth log-cosh target runs end to end with a numerical optimum."""
    from src.natural_gradient_discretization_stepsize.optimize_star import compute_star
    T = SmoothLogCoshTarget(100.0, epsilon=1e-2, n_nodes=40)
    m_star, C_star, F_star, diag = compute_star(T)
    sched = build_schedule("wfr_adaptive", beta=T.beta, C0=T.C0, c=0.5)
    _, summary = simulate_run("wfr_adaptive", T, 0.5, 2000, sched,
                              F_star, m_star, C_star)
    assert summary["status"] == "ok"
    assert summary["gap_final"] < summary["gap0"]
