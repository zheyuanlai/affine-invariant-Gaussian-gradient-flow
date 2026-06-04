"""Gaussian ground truth on the torch dense path.

For ``V = 0.5||theta||^2`` the exact answers are Lambda_hat = 0 and gamma_loc = 1.
With M large vs the operator dimension (small N) the torch dense estimates should
be near those values, and the clean diagonal / self-adjointness diagnostics hold
at every tested dimension (mirrors the NumPy Gaussian test, on torch CPU).
"""
import pytest

torch = pytest.importorskip("torch")

from src.common.monte_carlo import gaussian_samples
from src.natural_gradient_local_rate.potentials import GaussianPotential
from src.natural_gradient_local_rate import torch_backend as tb


def _setup(N, M, seed=0):
    Z = gaussian_samples(N, M, seed=seed, antithetic=True)
    pot_t = tb.torch_potential_from_cpu(GaussianPotential(N), torch.device("cpu"), torch.float64)
    return pot_t, torch.as_tensor(Z, dtype=torch.float64)


@pytest.mark.parametrize("N", [4, 8, 16])
def test_gaussian_full_sym_near_zero_gamma_near_one(N):
    pot_t, Zt = _setup(N, 65536)
    lam, resH = tb.torch_estimate_lambda_hat_dense(pot_t, Zt)
    gam, resL, _, _ = tb.torch_estimate_gamma_loc_dense(pot_t, Zt)
    assert lam < 0.30, f"Lambda_hat_full_sym={lam} at N={N}"
    assert 0.85 < gam < 1.10, f"gamma_loc={gam} at N={N}"
    assert resH < 1e-8 and resL < 1e-8


@pytest.mark.parametrize("N", [4, 8, 16, 32])
def test_gaussian_clean_diagnostics_dimension_robust(N):
    pot_t, Zt = _setup(N, 16384)
    # The cleanest, dimension-robust separable estimator is max_i A_ii (~0 here);
    # the full N x N diagonal *eigenvalue* inflates with sqrt(N/M) and is not used.
    A = tb.torch_diagonal_A_matrix(pot_t, Zt)
    max_diag = float(A.diagonal().max().item())
    assert max_diag < 0.06, f"max_diag={max_diag} at N={N}"
    est = tb.torch_operator_estimates(pot_t, Zt, self_adjoint_probes=6)
    assert est["self_adjoint_error_H_sym"] < 1e-9
    assert est["self_adjoint_error_L_star"] < 1e-9
    # raw forward is demonstrably not self-adjoint
    assert est["self_adjoint_error_H_raw"] > 1e-4
