"""Centralized theory constants for the Gaussian natural-gradient discretizations.

Single source of truth for the stepsize and contraction constants used by the
experiment runners, summary scripts, tests, and report-asset generation, so that
no theory formula is hard-coded in more than one place.

Current theory (the improved KL proof)
--------------------------------------
For the ``alpha``-strongly-log-concave, ``beta``-smooth target the Riemannian /
geodesic and the KL / Bregman discretizations of the Gaussian natural-gradient
flow now share the *same* theorem-safe stepsize scale::

    L_Riem = L_KL = beta * lambda_max,
    dt_Riem_theory = dt_KL_theory = 1 / (beta * lambda_max).

The KL proof no longer carries the obsolete cubic penalty
``max{1, lambda_max^3 / (2 lambda_min^3)}``; that factor is gone. The remaining
theoretical difference between the two schemes is the per-step *contraction*
factor (:func:`q_riem` vs :func:`q_kl`), not the admissible stepsize.

Spectral bounds (log-concave case, ``Lemma cov-bound`` in the manuscript)::

    lambda_min = min(lambda0_min, 1 / beta),
    lambda_max = max(lambda0_max, 1 / alpha),

with ``lambda0_min``/``lambda0_max`` the extreme eigenvalues of the initial
covariance ``C0``.

Projected / clipped KL (non-log-concave) theory
-----------------------------------------------
Under the global Hessian bound ``||grad^2 V|| <= beta`` and covariance clipped to
``[lambda_minus, lambda_plus]`` the current projected-KL theorem uses::

    L_clip = 2 * beta * lambda_plus,
    dt_projected_KL_theory = 1 / (2 * beta * lambda_plus),

which depends on ``lambda_plus`` and ``beta`` only -- *not* on ``lambda_minus``.
The certificate is a constrained Bregman stationarity bound on the
covariance-truncated feasible set,
``min_{0<=n<N} KL(rho_{a_n} || rho_{a_{n+1}}) <= (dt / N) (E(a_0) - E(a_N))``;
it does not certify that the unconstrained Fisher-Rao gradient is small.

Deprecated formulas (cubic KL penalty, ``lambda_plus^4 / lambda_minus^3`` clip
constant) are kept under ``deprecated_old_*`` names for historical comparison
only and must never appear in current summaries or report prose.
"""
from __future__ import annotations

# Tags recorded in output files so a CSV/figure is self-describing about which
# theory generated it.
THEORY_VERSION_LOGCONCAVE = "kl_beta_lambda_max_no_cubic_penalty"
THEORY_VERSION_PROJECTED_KL = "projected_kl_L_clip_2_beta_lambda_plus"


# ---------------------------------------------------------------------------
# Log-concave spectral bounds and theorem-safe stepsizes
# ---------------------------------------------------------------------------

def natural_gradient_spectral_bounds(lambda0_min, lambda0_max, alpha, beta):
    """Return ``(lambda_min, lambda_max)`` for the log-concave covariance bounds.

        lambda_min = min(lambda0_min, 1 / beta),
        lambda_max = max(lambda0_max, 1 / alpha).
    """
    lambda_min = min(float(lambda0_min), 1.0 / float(beta))
    lambda_max = max(float(lambda0_max), 1.0 / float(alpha))
    return float(lambda_min), float(lambda_max)


def riemannian_theory_constants(alpha, beta, lambda0_min, lambda0_max):
    """Riemannian/geodesic theorem constants.

    Returns ``lambda_min, lambda_max, L_Riem, dt_Riem_theory`` with
    ``L_Riem = beta * lambda_max`` and ``dt_Riem_theory = 1 / L_Riem``.
    """
    lam_min, lam_max = natural_gradient_spectral_bounds(
        lambda0_min, lambda0_max, alpha, beta)
    L_riem = float(beta) * lam_max
    return {
        "lambda_min": lam_min,
        "lambda_max": lam_max,
        "L_Riem": float(L_riem),
        "dt_Riem_theory": float(1.0 / L_riem),
    }


def kl_theory_constants(alpha, beta, lambda0_min, lambda0_max):
    """KL/Bregman theorem constants under the improved proof.

    Returns ``lambda_min, lambda_max, L_KL, dt_KL_theory`` with
    ``L_KL = beta * lambda_max`` (no ``lambda_max^3 / lambda_min^3`` penalty) and
    ``dt_KL_theory = 1 / L_KL``. By construction this equals the Riemannian
    theorem-safe stepsize for any smooth log-concave target.
    """
    lam_min, lam_max = natural_gradient_spectral_bounds(
        lambda0_min, lambda0_max, alpha, beta)
    # Improved KL proof: the Bregman relative-smoothness constant is beta*lambda_max,
    # matching the Riemannian L-smoothness constant exactly.
    L_kl = float(beta) * lam_max
    return {
        "lambda_min": lam_min,
        "lambda_max": lam_max,
        "L_KL": float(L_kl),
        "dt_KL_theory": float(1.0 / L_kl),
    }


# ---------------------------------------------------------------------------
# Per-step contraction factors
# ---------------------------------------------------------------------------

def q_riem(dt, alpha, beta, lambda_min, lambda_max):
    """Riemannian theorem per-step contraction factor.

        q_Riem(dt) = 1 - alpha * lambda_min * dt * (2 - beta * lambda_max * dt).
    """
    return 1.0 - alpha * lambda_min * dt * (2.0 - beta * lambda_max * dt)


def q_kl(dt, alpha, beta, lambda_min, lambda_max):
    """KL theorem per-step contraction factor (improved proof).

        q_KL(dt) = 1 - alpha * lambda_min * dt
                       / (2 * (1 + dt) * (1 + dt * beta * lambda_max)).

    This is the proven contraction of ``Theorem (conv-KL)`` and is valid on the
    same theorem-safe stepsize range ``dt <= 1 / (beta * lambda_max)`` as the
    Riemannian factor.
    """
    return 1.0 - alpha * lambda_min * dt / (
        2.0 * (1.0 + dt) * (1.0 + dt * beta * lambda_max))


# ---------------------------------------------------------------------------
# Projected / clipped KL (non-log-concave) theorem constants
# ---------------------------------------------------------------------------

def projected_kl_theory_constants(beta, lambda_plus):
    """Projected-KL theorem constants.

    Returns ``L_clip = 2 * beta * lambda_plus`` and
    ``dt_projected_KL_theory = 1 / L_clip = 1 / (2 * beta * lambda_plus)``. The
    theorem-safe scale depends on ``lambda_plus`` and ``beta`` only, not on
    ``lambda_minus``.
    """
    L_clip = 2.0 * float(beta) * float(lambda_plus)
    return {
        "L_clip": float(L_clip),
        "dt_projected_KL_theory": float(1.0 / L_clip),
    }


# ---------------------------------------------------------------------------
# Deprecated formulas -- historical / regression comparison ONLY.
# Never use these in current summaries, README prose, or report narrative.
# ---------------------------------------------------------------------------

def deprecated_old_kl_stepsize_factor(lambda_min, lambda_max):
    """OBSOLETE cubic KL stepsize penalty ``max{1, lambda_max^3/(2 lambda_min^3)}``.

    Superseded by the improved KL proof (``L_KL = beta * lambda_max``). Retained
    only so historical comparisons can recover the old, conservative value; it is
    never used as current theory.
    """
    return max(1.0, float(lambda_max) ** 3 / (2.0 * float(lambda_min) ** 3))


def deprecated_old_projected_kl_smoothness_constant(beta, lambda_minus, lambda_plus):
    """OBSOLETE clipped relative-smoothness constant
    ``beta * max(lambda_plus, lambda_plus^4 / lambda_minus^3)``.

    Superseded by ``L_clip = 2 * beta * lambda_plus``. Retained for historical
    comparison only; never used as current theory.
    """
    return float(beta) * max(
        float(lambda_plus), float(lambda_plus) ** 4 / float(lambda_minus) ** 3)
