"""Affine-invariant Gaussian gradient flow experiments, parameterized by (omega, tau).

Study of parameter effects (omega, tau) in affine-invariant Gaussian gradient
flows, for both an exact Gaussian target N(0, I_n) and a strongly log-concave
non-Gaussian target. See ``reports/affine_invariant_omega_tau_report.tex`` for
the detailed write-up.

A few convenience re-exports are provided for interactive use; plotting modules
are intentionally *not* imported here to avoid importing matplotlib eagerly.
"""
from src.omega_tau_modes.dynamics import gaussian_step
from src.omega_tau_modes.metrics import compute_all_metrics, kl_energy
from src.omega_tau_modes.initializations import get_initialization, INIT_NAMES
from src.omega_tau_modes.utils import validate_params, make_q_vector

__all__ = [
    "gaussian_step",
    "compute_all_metrics",
    "kl_energy",
    "get_initialization",
    "INIT_NAMES",
    "validate_params",
    "make_q_vector",
]
