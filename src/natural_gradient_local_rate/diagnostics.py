"""Theoretical reference bounds for the natural-gradient local rate.

All manuscript constants live here so they are easy to revise in one place.
``kappa`` is ``kappa_target`` and ``beta_target`` is the nominal smoothness.

    Gamma_current        = min( beta - 1, N_theta * (3 + 4/sqrt(pi) (1 + log kappa)) )
    current_bound_rate   = 1 / (4 + Gamma_current)
    D_kappa              = (4 + 4/sqrt(pi)) (1 + log kappa)
    conjecture_bound_rate = 1 / (4 + D_kappa)
"""
from __future__ import annotations

import numpy as np

# --- revisable constants ---
FOUR_OVER_SQRT_PI = 4.0 / np.sqrt(np.pi)
BULK_CONST = 3.0          # the "3 +" in Gamma_current (tail-cut bound on |G_ij - 1|)
DENOM_CONST = 4.0         # the "4 +" in the rate denominators
CONJECTURE_LEAD = 4.0 + FOUR_OVER_SQRT_PI  # leading constant of D_kappa


def log_kappa_factor(kappa):
    """``1 + log(kappa)``."""
    return 1.0 + np.log(kappa)


def gamma_current(N_theta, kappa, beta_target):
    """Manuscript-style ``Gamma = min(beta-1, N_theta (3 + 4/sqrt(pi)(1+log k)))``."""
    return float(min(
        beta_target - 1.0,
        N_theta * (BULK_CONST + FOUR_OVER_SQRT_PI * log_kappa_factor(kappa)),
    ))


def current_bound_rate(N_theta, kappa, beta_target):
    """Current manuscript bound on the local rate: ``1 / (4 + Gamma_current)``."""
    return 1.0 / (DENOM_CONST + gamma_current(N_theta, kappa, beta_target))


def D_kappa(kappa):
    """Conjectural dimension-free quantity ``(4 + 4/sqrt(pi))(1 + log kappa)``."""
    return float(CONJECTURE_LEAD * log_kappa_factor(kappa))


def conjecture_bound_rate(kappa):
    """Conjectural dimension-free rate ``1 / (4 + D_kappa)``."""
    return 1.0 / (DENOM_CONST + D_kappa(kappa))


def reference_columns(N_theta, kappa, beta_target, Lambda_hat=None, gamma_loc=None):
    """Assemble the theoretical reference columns for a results row."""
    lk = float(log_kappa_factor(kappa))
    out = {
        "log_kappa_factor": lk,
        "Gamma_current": gamma_current(N_theta, kappa, beta_target),
        "current_bound_rate": current_bound_rate(N_theta, kappa, beta_target),
        "D_kappa": D_kappa(kappa),
        "conjecture_bound_rate": conjecture_bound_rate(kappa),
    }
    if Lambda_hat is not None:
        out["lambda_over_logkappa"] = float(Lambda_hat) / lk
    if gamma_loc is not None and gamma_loc != 0.0:
        out["inverse_gamma_loc"] = 1.0 / float(gamma_loc)
        out["inverse_gamma_over_logkappa"] = (1.0 / float(gamma_loc)) / lk
    return out
