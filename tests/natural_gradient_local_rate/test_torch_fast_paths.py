"""Gaussian-analytic and separable-diagonal fast accumulation paths.

Both fast paths must produce the *same* ``(G, T_mat, A)`` triple -- and hence the
same ``Lambda_hat_full_sym`` / ``Lambda_hat_diag`` / ``gamma_loc`` /
self-adjointness -- as the generic dense accumulator, on the same bank. They run
on torch CPU (and on CUDA when available).
"""
import numpy as np
import pytest

torch = pytest.importorskip("torch")

from src.common.monte_carlo import gaussian_samples
from src.natural_gradient_local_rate.potentials import build_potential, GaussianPotential
from src.natural_gradient_local_rate import torch_backend as tb


def _devices():
    devs = [torch.device("cpu")]
    if torch.cuda.is_available():
        devs.append(torch.device("cuda"))
    return devs


@pytest.mark.parametrize("device", _devices())
@pytest.mark.parametrize("N", [2, 4, 8])
@pytest.mark.parametrize("M", [512, 2048])
@pytest.mark.parametrize("antithetic", [True, False])
def test_gaussian_analytic_matches_generic(N, M, antithetic, device):
    """Gaussian analytic (G, T_mat, A) == generic, and downstream estimates agree."""
    Z = gaussian_samples(N, M, seed=1, antithetic=antithetic)
    Zt = torch.as_tensor(Z, dtype=torch.float64, device=device)
    pot = tb.torch_potential_from_cpu(GaussianPotential(N), device, torch.float64)
    assert tb._resolve_accumulation_mode(pot) == "gaussian_analytic"

    g_fast = tb._accumulate_dense(pot, Zt, mode="gaussian_analytic")
    g_gen = tb._accumulate_dense(pot, Zt, mode="generic_dense")
    for a, b, name in zip(g_fast, g_gen, ("G", "T_mat", "A")):
        assert a.cpu().numpy() == pytest.approx(b.cpu().numpy(), abs=1e-10), name

    # full estimate suite via the auto-dispatched path equals the generic numbers
    est_fast = tb.torch_operator_estimates(pot, Zt, compute_gamma_loc=True)
    assert est_fast["torch_accumulation_mode"] == "gaussian_analytic"
    # build a generic-mode reference H_sym / L_star by hand
    G, T_mat, _ = g_gen
    H = 0.5 * (G + G.T)
    lam_gen = float(torch.linalg.eigvalsh(H)[-1].item())
    assert est_fast["Lambda_hat_full_sym"] == pytest.approx(lam_gen, abs=1e-9)


@pytest.mark.parametrize("device", _devices())
@pytest.mark.parametrize("N", [4, 8])
@pytest.mark.parametrize("M", [1024, 2048])
def test_separable_fast_matches_generic(N, M, device):
    """Separable diagonal fast path == generic on (G, T_mat, A) and all estimates."""
    Z = gaussian_samples(N, M, seed=2, antithetic=True)
    cpu = build_potential("separable", N, kappa_target=8.0, seed=0, Z_ref=Z)
    pot = tb.torch_potential_from_cpu(cpu, device, torch.float64)
    Zt = torch.as_tensor(Z, dtype=torch.float64, device=device)
    assert tb._resolve_accumulation_mode(pot) == "separable_diagonal_fast"

    g_fast = tb._accumulate_dense(pot, Zt, mode="separable_diagonal_fast")
    g_gen = tb._accumulate_dense(pot, Zt, mode="generic_dense")
    for a, b, name in zip(g_fast, g_gen, ("G", "T_mat", "A")):
        assert a.cpu().numpy() == pytest.approx(b.cpu().numpy(), abs=1e-10), name

    est_fast = tb.torch_operator_estimates(pot, Zt, compute_gamma_loc=True)
    # force a generic comparison by temporarily evaluating the generic suite
    G, T_mat, A = g_gen
    H = 0.5 * (G + G.T)
    lam_gen = float(torch.linalg.eigvalsh(H)[-1].item())
    A_sym = 0.5 * (A + A.T)
    diag_gen = float(torch.linalg.eigvalsh(A_sym)[-1].item())
    assert est_fast["torch_accumulation_mode"] == "separable_diagonal_fast"
    assert est_fast["Lambda_hat_full_sym"] == pytest.approx(lam_gen, abs=1e-9)
    assert est_fast["Lambda_hat_diag"] == pytest.approx(diag_gen, abs=1e-9)
    assert est_fast["self_adjoint_error_H_sym"] < 1e-9
    assert est_fast["self_adjoint_error_L_star"] < 1e-9


@pytest.mark.parametrize("N", [2, 4, 8])
@pytest.mark.parametrize("M", [512, 2048])
@pytest.mark.parametrize("antithetic", [True, False])
def test_gaussian_analytic_compute_row_matches_numpy(N, M, antithetic):
    """End-to-end compute_row (torch gaussian-analytic) == numpy row for gaussian."""
    from src.natural_gradient_local_rate.estimator_suite import compute_row
    Z = gaussian_samples(N, M, seed=5, antithetic=antithetic)
    cpu = GaussianPotential(N)
    point = {"family": "gaussian", "N_theta": N, "kappa_target": 1.0, "seed": 0, "M_mc": M}
    r_np = compute_row(cpu, Z, point, {"backend": "numpy"})
    r_t = compute_row(cpu, Z, point, {"backend": "torch", "device": "cpu", "dtype": "float64"})
    assert r_t["torch_accumulation_mode"] == "gaussian_analytic"
    for key in ["Lambda_hat_full_sym", "Lambda_hat_diag", "gamma_loc"]:
        a, b = r_np[key], r_t[key]
        if np.isnan(a) and np.isnan(b):
            continue
        assert b == pytest.approx(a, abs=1e-7), f"{key}: numpy={a} torch={b}"
