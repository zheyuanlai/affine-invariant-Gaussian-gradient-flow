import pandas as pd
import pytest

from src.natural_gradient_local_rate.scaling_diagnostics import (
    fit_scaling_diagnostics, sym_dim,
)


def test_scaling_fit_recovers_intercept():
    N = 4
    p = sym_dim(N)
    a = 0.25
    c = 1.3
    rows = []
    for M in [10_000, 40_000, 160_000, 640_000]:
        rows.append({
            "potential_family": "separable",
            "N_theta": N,
            "kappa_target": 5.0,
            "seed": 0,
            "M_mc": M,
            "Lambda_hat_full_sym": a + c * (p / M) ** 0.5,
            "Lambda_hat_separable_exact": a,
            "gamma_loc": 0.8,
        })
    fits = fit_scaling_diagnostics(pd.DataFrame(rows), tolerance=0.02)
    row = fits.iloc[0]
    assert row["Lambda_inf_fit_status"] == "ok"
    assert row["Lambda_inf_fit"] == pytest.approx(a, abs=1e-10)
    assert row["Lambda_inf_minus_exact"] == pytest.approx(0.0, abs=1e-10)
    assert row["fit_r2_scaling"] == pytest.approx(1.0)
    assert bool(row["converged_flag"]) is True


def test_scaling_fit_flags_unconverged_with_insufficient_data():
    df = pd.DataFrame([{
        "potential_family": "random_feature",
        "N_theta": 4,
        "kappa_target": 5.0,
        "seed": 0,
        "M_mc": 1000,
        "Lambda_hat_full_sym": 2.0,
        "gamma_loc": -0.1,
    }])
    row = fit_scaling_diagnostics(df).iloc[0]
    assert row["Lambda_inf_fit_status"] == "insufficient_data"
    assert bool(row["converged_flag"]) is False
    assert row["convergence_warning"] == "insufficient_data"
