import numpy as np
import pandas as pd
import pytest

from src.natural_gradient_local_rate.baseline_postprocessing import add_baseline_corrections


def _row(fam, seed, M, lam, gamma):
    return {
        "potential_family": fam,
        "family": fam,
        "N_theta": 4,
        "kappa_target": 5.0,
        "M_mc": M,
        "seed": seed,
        "Lambda_hat_full_sym": lam,
        "gamma_loc": gamma,
        "inverse_gamma_loc": 1.0 / gamma,
        "Lambda_hat_diag": lam / 2.0,
        "Lambda_hat_separable_exact": 0.25 if fam == "separable" else np.nan,
    }


def test_baseline_postprocessing_seed_exact_and_group_mean():
    df = pd.DataFrame([
        _row("gaussian", 0, 100, 0.20, 0.90),
        _row("separable", 0, 100, 0.30, 0.80),
        _row("random_feature", 0, 100, 0.55, 0.70),
        _row("gaussian", 1, 100, 0.22, 0.91),
        _row("separable", 1, 100, 0.32, 0.81),
        _row("random_feature", 9, 100, 0.62, 0.75),
    ])
    out = add_baseline_corrections(df)

    exact = out[(out["potential_family"] == "random_feature") & (out["seed"] == 0)].iloc[0]
    assert exact["gaussian_baseline_match"] == "seed_exact"
    assert exact["separable_baseline_match"] == "seed_exact"
    assert exact["full_sym_minus_gaussian"] == pytest.approx(0.35)
    assert exact["full_sym_minus_separable"] == pytest.approx(0.25)
    assert exact["full_sym_over_separable"] == pytest.approx(0.55 / 0.30)
    assert bool(exact["noise_baseline_available"]) is True

    mean = out[(out["potential_family"] == "random_feature") & (out["seed"] == 9)].iloc[0]
    assert mean["gaussian_baseline_match"] == "group_mean"
    assert mean["separable_baseline_match"] == "group_mean"
    assert mean["Lambda_hat_gaussian_baseline"] == pytest.approx(0.21)
    assert mean["Lambda_hat_separable_baseline"] == pytest.approx(0.31)


def test_baseline_postprocessing_missing_baseline_is_nan_and_statused():
    df = pd.DataFrame([_row("random_feature", 0, 200, 0.50, 0.70)])
    out = add_baseline_corrections(df)
    row = out.iloc[0]
    assert np.isnan(row["full_sym_minus_gaussian"])
    assert np.isnan(row["full_sym_minus_separable"])
    assert bool(row["noise_baseline_available"]) is False
    assert row["baseline_match"] == "missing"
    assert row["baseline_correction_status"] == "missing:all_baselines"
