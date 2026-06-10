"""Tests for the Hessian first-Hermite coupling ``tau_H``.

The revised local-rate route uses the deterministic implication

    gamma_loc >= 1 / (tau_H^2 + 3),

where ``tau_H = ||T||_op`` is the mean-covariance coupling. These tests keep the
estimator and the finite-bank inequality tied to the implemented operators.
"""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.natural_gradient_local_rate.estimator_suite import compute_row, default_opts
from src.natural_gradient_local_rate.linearized_rate import estimate_gamma_loc, estimate_tau_H
from src.natural_gradient_local_rate.potentials import build_potential


def test_tau_H_gaussian_is_zero_on_antithetic_bank():
    N = 6
    Z = gaussian_samples(N, 4096, seed=0, antithetic=True)
    pot = build_potential("gaussian", N, kappa_target=1.0, seed=0, Z_ref=Z)
    tau = estimate_tau_H(pot, Z)
    assert tau["tau_H"] == pytest.approx(0.0, abs=1e-12)
    assert tau["coupling_bound_rate"] == pytest.approx(1.0 / 3.0)


@pytest.mark.parametrize("family", ["separable", "random_feature", "product_feature"])
def test_coupling_bound_is_below_gamma_loc(family):
    N = 5
    Z = gaussian_samples(N, 8192, seed=1, antithetic=True)
    pot = build_potential(family, N, kappa_target=10.0, seed=2, Z_ref=Z,
                          feature_multiplier=6)
    tau = estimate_tau_H(pot, Z)
    gamma, _ = estimate_gamma_loc(pot, Z, eigsh_tol=1e-8)
    assert gamma + 1e-10 >= tau["coupling_bound_rate"]


def test_compute_row_reports_tau_columns_for_product_feature():
    point = {
        "family": "product_feature",
        "N_theta": 4,
        "kappa_target": 20.0,
        "seed": 0,
        "M_mc": 4096,
    }
    Z = gaussian_samples(point["N_theta"], point["M_mc"], seed=3, antithetic=True)
    pot = build_potential(point["family"], point["N_theta"], point["kappa_target"],
                          point["seed"], Z_ref=Z, feature_multiplier=6)
    opts = default_opts()
    opts.update({"chunk_size": 1024, "eigsh_tol": 1e-7})
    row = compute_row(pot, Z, point, opts)
    assert row["status"] == "ok"
    assert np.isfinite(row["tau_H"])
    assert np.isfinite(row["tau_H_sq"])
    assert np.isfinite(row["coupling_bound_rate"])
    assert np.isfinite(row["tau_top_total"])
    assert np.isfinite(row["tau_top_longitudinal"])
    assert np.isfinite(row["tau_top_mixed"])
    assert np.isfinite(row["tau_top_transverse"])
    assert row["tau_top_total"] == pytest.approx(row["tau_H"], abs=1e-10)
    assert (
        row["tau_top_longitudinal"]
        + row["tau_top_mixed"]
        + row["tau_top_transverse"]
    ) == pytest.approx(row["tau_top_total"], abs=1e-10)
    assert (
        row["tau_top_X_longitudinal_norm_sq"]
        + row["tau_top_X_mixed_norm_sq"]
        + row["tau_top_X_transverse_norm_sq"]
    ) == pytest.approx(1.0, abs=1e-10)
    assert row["gamma_over_coupling_bound"] >= 1.0
