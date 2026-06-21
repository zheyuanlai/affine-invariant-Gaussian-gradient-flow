"""Tests for the theoretical stepsize bounds and the run classification."""
import numpy as np

from src.natural_gradient_discretization_stepsize.targets import GaussianPosteriorTarget
from src.natural_gradient_discretization_stepsize.metrics import (
    theory_stepsize_bounds, dt_theory_for_method, max_feasible_dt, classify_run,
)


def test_theory_bounds_formula():
    """Improved KL proof: dt_Riem = dt_KL = 1/(beta lam_max), no cubic penalty."""
    alpha, beta = 0.1, 1.0
    C0 = np.diag([0.5, 2.0])
    b = theory_stepsize_bounds(alpha, beta, C0)
    # lam_min = min(0.5, 1/beta=1) = 0.5 ; lam_max = max(2.0, 1/alpha=10) = 10
    assert b["lambda_min"] == 0.5 and b["lambda_max"] == 10.0
    assert np.isclose(b["L_riem"], 1.0 * 10.0)
    # KL Lipschitz constant is now exactly beta*lambda_max (cubic factor removed).
    assert np.isclose(b["L_kl"], 1.0 * 10.0)
    assert np.isclose(b["dt_theory_riem"], 1.0 / 10.0)
    assert np.isclose(b["dt_theory_kl"], 1.0 / 10.0)
    assert b["theory_version"] == "kl_beta_lambda_max_no_cubic_penalty"


def test_kl_and_riemann_share_theorem_safe_stepsize():
    """For smooth log-concave targets the KL and Riemannian stepsizes are equal."""
    for lam in (0.01, 0.1, 1.0):
        T = GaussianPosteriorTarget(lam)
        b = theory_stepsize_bounds(T.alpha, T.beta, T.C0)
        assert np.isclose(b["dt_theory_kl"], b["dt_theory_riem"], rtol=0, atol=1e-15)


def test_kl_stepsize_has_no_lambda_min_cubic_dependence():
    """dt_KL = 1/(beta lam_max) is invariant to lambda0_min (no lam_min^3 term)."""
    alpha, beta = 0.1, 1.0
    # Same lambda0_max=2 but very different lambda0_min: the cubic penalty would
    # have made dt_KL differ by orders of magnitude; the new bound is identical.
    b_wide = theory_stepsize_bounds(alpha, beta, np.diag([1e-4, 2.0]))
    b_narrow = theory_stepsize_bounds(alpha, beta, np.diag([1.0, 2.0]))
    assert np.isclose(b_wide["dt_theory_kl"], b_narrow["dt_theory_kl"])
    assert np.isclose(b_wide["dt_theory_kl"], b_wide["dt_theory_riem"])


def test_theory_unavailable():
    """No (alpha, beta) -> theory_bound_available False and None stepsizes."""
    b = theory_stepsize_bounds(None, None, np.eye(2))
    assert b["theory_bound_available"] is False
    assert dt_theory_for_method(b, "kl") is None
    assert dt_theory_for_method(b, "riemannian") is None


def test_max_feasible_dt():
    """max_feasible_dt returns the largest passing stepsize (not assuming monotone)."""
    dts = [0.01, 0.1, 0.5, 1.0, 5.0]
    flags = [1, 1, 1, 0, 0]
    assert max_feasible_dt(dts, flags) == 0.5
    # non-monotone boundary: a later dt passes again
    assert max_feasible_dt(dts, [1, 0, 1, 0, 0]) == 0.5
    assert np.isnan(max_feasible_dt(dts, [0, 0, 0, 0, 0]))


# ---------------------------------------------------------------------------
# classify_run on synthetic trajectories
# ---------------------------------------------------------------------------

def test_classify_monotone_decreasing():
    """A clean monotone-decreasing trajectory is SPD/stable/monotone."""
    F = np.array([10.0, 5.0, 2.0, 1.0, 0.5])
    gap = F - 0.0
    min_eig = np.full(5, 0.5)
    c = classify_run(F, gap, min_eig, F0=10.0, F_star=0.0,
                     terminal_accuracy_error=1e-3)
    assert c["spd_feasible"] and c["stable"] and c["monotone"] and c["accurate"]
    assert c["num_energy_increases"] == 0


def test_classify_nonmonotone_but_stable():
    """A bump that still ends below F0 is stable but not monotone."""
    F = np.array([10.0, 11.0, 6.0, 2.0])     # one increase, ends well below F0
    c = classify_run(F, F - 0.0, np.full(4, 0.4), F0=10.0, F_star=0.0,
                     terminal_accuracy_error=0.5)
    assert c["spd_feasible"] and c["stable"]
    assert not c["monotone"] and c["num_energy_increases"] == 1
    assert not c["accurate"]                  # accuracy error 0.5 > 1e-2


def test_classify_non_spd():
    """A non-SPD eigenvalue fails every criterion."""
    F = np.array([10.0, 5.0, 2.0])
    min_eig = np.array([0.5, 1e-15, 0.5])     # dips below SPD tol
    c = classify_run(F, F - 0.0, min_eig, F0=10.0, F_star=0.0,
                     terminal_accuracy_error=1e-4)
    assert not c["spd_feasible"]
    assert not c["stable"] and not c["monotone"] and not c["accurate"]


def test_classify_explosion():
    """A blow-up beyond explosion_factor * gap0 is not stable."""
    F = np.array([10.0, 50.0, 1e7])
    c = classify_run(F, F - 0.0, np.full(3, 0.5), F0=10.0, F_star=0.0,
                     terminal_accuracy_error=np.inf)
    assert not c["stable"]
    assert c["max_gap_ratio"] > 1e3


def test_classify_nan_fails():
    """NaN/Inf anywhere fails SPD-feasibility."""
    F = np.array([10.0, np.nan, 2.0])
    c = classify_run(F, F - 0.0, np.array([0.5, np.nan, 0.5]),
                     F0=10.0, F_star=0.0, terminal_accuracy_error=np.inf)
    assert not c["spd_feasible"] and not c["finite_ok"]
