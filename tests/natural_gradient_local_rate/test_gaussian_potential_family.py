import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.natural_gradient_local_rate.estimator_suite import compute_row
from src.natural_gradient_local_rate.potentials import build_potential, GaussianPotential


def _row(N, M, *, backend="numpy", device="cpu"):
    Z = gaussian_samples(N, M, seed=10, antithetic=True)
    pot = build_potential("gaussian", N, kappa_target=5.0, seed=0, Z_ref=Z)
    point = {"family": "gaussian", "N_theta": N, "kappa_target": 5.0,
             "seed": 0, "M_mc": M}
    return compute_row(
        pot, Z, point,
        {"backend": backend, "device": device, "dtype": "float64",
         "compute_gamma_loc": True, "eigsh_tol": 1e-6, "eigsh_maxiter": 1000},
    )


def test_gaussian_family_builds_cpu_and_records_true_columns():
    pot = build_potential("gaussian", N_theta=4, kappa_target=5.0, seed=3)
    assert isinstance(pot, GaussianPotential)
    row = _row(3, 2048)
    assert row["status"] == "ok"
    assert row["baseline_type"] == "gaussian"
    assert row["Lambda_true"] == pytest.approx(0.0)
    assert row["gamma_true"] == pytest.approx(1.0)
    assert row["Lambda_hat_separable_exact"] == pytest.approx(0.0, abs=1e-12)


def test_gaussian_cpu_full_noise_decreases_and_gamma_improves_with_M():
    small = _row(3, 512)
    large = _row(3, 8192)
    assert large["Lambda_hat_full_sym"] < small["Lambda_hat_full_sym"]
    assert abs(large["gamma_loc"] - 1.0) < abs(small["gamma_loc"] - 1.0)
    assert abs(large["Lambda_hat_diag"]) < 0.05
    assert large["self_adjoint_error_H_sym"] < 1e-10
    assert large["self_adjoint_error_L_star"] < 1e-10


def test_gaussian_torch_cpu_if_available_records_true_columns():
    torch = pytest.importorskip("torch")
    row = _row(3, 1024, backend="torch", device="cpu")
    assert row["status"] == "ok"
    assert row["backend"] == "torch"
    assert row["device"] == "cpu"
    assert row["Lambda_true"] == pytest.approx(0.0)
    assert row["gamma_true"] == pytest.approx(1.0)
    assert np.isfinite(row["Lambda_hat_full_sym"])


@pytest.mark.gpu
def test_gaussian_torch_cuda_tiny_if_available():
    torch = pytest.importorskip("torch")
    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    row = _row(3, 512, backend="torch", device="cuda")
    assert row["status"] == "ok"
    assert row["device"].startswith("cuda")
    assert row["Lambda_true"] == pytest.approx(0.0)
    assert row["gamma_true"] == pytest.approx(1.0)
