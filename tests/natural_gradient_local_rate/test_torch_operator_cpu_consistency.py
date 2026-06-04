"""Torch CPU backend must reproduce the NumPy estimators on a small case.

These run on torch *CPU* (no CUDA needed) and compare against the corrected
NumPy/SciPy path to float64 precision. A CUDA-only test (CPU vs GPU) is included
and skipped when no GPU is present.
"""
import numpy as np
import pytest

torch = pytest.importorskip("torch")

from src.common.monte_carlo import gaussian_samples
from src.common.symspace import random_symmetric, sym_norm
from src.natural_gradient_local_rate.potentials import build_potential, GaussianPotential
from src.natural_gradient_local_rate.operators import apply_H_sym_banked, HessianBank, _hsym_from_bank
from src.natural_gradient_local_rate.linearized_rate import (
    estimate_lambda_hat, estimate_gamma_loc, estimate_diagonal_lambda,
)
from src.natural_gradient_local_rate import torch_backend as tb
from src.natural_gradient_local_rate.estimator_suite import compute_row

_CPU = "cpu"


def _pair(fam, N=4, M=4096, seed=0):
    Z = gaussian_samples(N, M, seed=seed + 1, antithetic=True)
    cpu = (GaussianPotential(N) if fam == "gaussian"
           else build_potential(fam, N, kappa_target=8.0, seed=seed, Z_ref=Z))
    pot_t = tb.torch_potential_from_cpu(cpu, torch.device(_CPU), torch.float64)
    Zt = torch.as_tensor(Z, dtype=torch.float64)
    return cpu, Z, pot_t, Zt


@pytest.mark.parametrize("fam", ["gaussian", "separable", "random_feature"])
def test_batch_hess_matches_numpy(fam):
    cpu, Z, pot_t, Zt = _pair(fam)
    H_np = cpu.batch_hess(Z[:7])
    H_t = pot_t.batch_hess(Zt[:7]).numpy()
    assert H_t == pytest.approx(H_np, abs=1e-10)


@pytest.mark.parametrize("fam", ["separable", "random_feature"])
def test_H_sym_apply_matches_numpy(fam):
    cpu, Z, pot_t, Zt = _pair(fam)
    rng = np.random.default_rng(0)
    X = random_symmetric(cpu.N_theta, rng)
    H_np = _hsym_from_bank(HessianBank(cpu, Z), X)
    H_t = tb.torch_apply_H_sym(pot_t, Zt, torch.as_tensor(X, dtype=torch.float64)).numpy()
    assert H_t == pytest.approx(H_np, abs=1e-9)


@pytest.mark.parametrize("fam", ["gaussian", "separable", "random_feature"])
def test_lambda_gamma_diag_match_numpy(fam):
    cpu, Z, pot_t, Zt = _pair(fam)
    lam_np = estimate_lambda_hat(cpu, Z)
    gam_np, _ = estimate_gamma_loc(cpu, Z)
    diag_np = estimate_diagonal_lambda(cpu, Z)["Lambda_hat_diag"]
    lam_t, resH = tb.torch_estimate_lambda_hat_dense(pot_t, Zt)
    gam_t, resL, _, _ = tb.torch_estimate_gamma_loc_dense(pot_t, Zt)
    diag_t, _ = tb.torch_diagonal_lambda(pot_t, Zt)
    assert lam_t == pytest.approx(lam_np, abs=1e-8)
    assert gam_t == pytest.approx(gam_np, abs=1e-8)
    assert diag_t == pytest.approx(diag_np, abs=1e-8)
    assert resH < 1e-8 and resL < 1e-8


@pytest.mark.parametrize("fam", ["separable", "random_feature"])
def test_compute_row_torch_matches_numpy_row(fam):
    """End-to-end: the torch row equals the numpy row on the same potential/bank."""
    cpu, Z, _, _ = _pair(fam, N=4, M=8192)
    point = {"family": fam, "N_theta": 4, "kappa_target": 8.0, "seed": 0, "M_mc": Z.shape[0]}
    opts_np = {"backend": "numpy"}
    opts_t = {"backend": "torch", "device": "cpu", "dtype": "float64"}
    r_np = compute_row(cpu, Z, point, opts_np)
    r_t = compute_row(cpu, Z, point, opts_t)
    assert r_np["status"] == "ok" and r_t["status"] == "ok"
    assert r_t["backend"] == "torch" and r_np["backend"] == "numpy"
    for key in ["Lambda_hat_full_sym", "Lambda_hat_diag", "gamma_loc",
                "Lambda_hat_separable_exact"]:
        a, b = r_np[key], r_t[key]
        if np.isnan(a) and np.isnan(b):
            continue
        assert b == pytest.approx(a, abs=1e-7), f"{key}: numpy={a} torch={b}"
    # self-adjointness diagnostics near machine precision on both paths
    assert r_t["self_adjoint_error_H_sym"] < 1e-9
    assert r_t["self_adjoint_error_L_star"] < 1e-9


@pytest.mark.gpu
@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
@pytest.mark.parametrize("fam", ["separable", "random_feature"])
def test_cuda_matches_cpu(fam):
    """torch CUDA dense estimates match torch CPU within float64 tolerance."""
    cpu, Z, _, _ = _pair(fam, N=4, M=8192)
    Zt_cpu = torch.as_tensor(Z, dtype=torch.float64, device="cpu")
    Zt_gpu = torch.as_tensor(Z, dtype=torch.float64, device="cuda")
    pot_cpu = tb.torch_potential_from_cpu(cpu, torch.device("cpu"), torch.float64)
    pot_gpu = tb.torch_potential_from_cpu(cpu, torch.device("cuda"), torch.float64)
    lam_c, _ = tb.torch_estimate_lambda_hat_dense(pot_cpu, Zt_cpu)
    lam_g, _ = tb.torch_estimate_lambda_hat_dense(pot_gpu, Zt_gpu)
    gam_c, _, _, _ = tb.torch_estimate_gamma_loc_dense(pot_cpu, Zt_cpu)
    gam_g, _, _, _ = tb.torch_estimate_gamma_loc_dense(pot_gpu, Zt_gpu)
    assert lam_g == pytest.approx(lam_c, abs=1e-8)
    assert gam_g == pytest.approx(gam_c, abs=1e-8)
