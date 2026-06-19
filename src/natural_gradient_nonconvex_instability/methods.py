"""Scalar covariance updates for the nonconvex-instability experiments."""
from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ScalarStep:
    c_next: float
    status: str
    denom: float | None = None
    ctilde: float | None = None


def riemannian_next(c: float, A: float, dt: float) -> ScalarStep:
    """Riemannian Fisher-Rao scalar covariance step."""
    exponent = float(dt) * (1.0 - float(c) * float(A))
    try:
        c_next = float(c) * math.exp(exponent)
    except OverflowError:
        return ScalarStep(float("inf"), "overflow")
    status = "ok" if math.isfinite(c_next) and c_next > 0.0 else "overflow"
    return ScalarStep(c_next, status)


def riemannian_next_log_c(log_c: float, c: float, A: float, dt: float) -> float:
    """Log-domain Riemannian update used by the runner's overflow guard."""
    return float(log_c) + float(dt) * (1.0 - float(c) * float(A))


def kl_next(c: float, A: float, dt: float) -> ScalarStep:
    """KL/Bregman Fisher-Rao scalar covariance step."""
    denom = 1.0 + float(dt) * float(c) * float(A)
    if denom <= 0.0 or not math.isfinite(denom):
        return ScalarStep(float("nan"), "non_spd", denom=denom)
    c_next = (1.0 + float(dt)) * float(c) / denom
    status = "ok" if math.isfinite(c_next) and c_next > 0.0 else "overflow"
    return ScalarStep(c_next, status, denom=denom)


def wasserstein_fb_next(c: float, A: float, eta: float) -> ScalarStep:
    """Bures-Wasserstein forward-backward scalar covariance step."""
    c = float(c)
    A = float(A)
    eta = float(eta)
    ctilde = (1.0 - eta * A) ** 2 * c
    radicand = ctilde * (ctilde + 4.0 * eta)
    if radicand < 0.0 or not math.isfinite(radicand):
        return ScalarStep(float("nan"), "non_finite", ctilde=ctilde)
    c_next = 0.5 * (ctilde + 2.0 * eta + math.sqrt(radicand))
    status = "ok" if math.isfinite(c_next) and c_next > 0.0 else "non_finite"
    return ScalarStep(c_next, status, ctilde=ctilde)


def clip_scalar(x: float, lo: float, hi: float) -> float:
    """Scalar covariance clip ``min(hi, max(lo, x))`` (eigenvalue clip in 1-D)."""
    return float(min(float(hi), max(float(lo), float(x))))


def clipped_kl_next(c: float, A: float, dt: float,
                    lambda_minus: float, lambda_plus: float) -> ScalarStep:
    """Projected KL/Bregman scalar covariance step.

    The unprojected KL covariance update is
    ``ctilde = (1 + dt) c / (1 + dt c A)``; the projected step clips the result
    back into the feasible interval ``[lambda_minus, lambda_plus]`` (the 1-D
    eigenvalue clip of Theorem 2.18). ``ctilde`` is the unclipped update and
    ``c_next`` the clipped covariance.
    """
    denom = 1.0 + float(dt) * float(c) * float(A)
    if denom <= 0.0 or not math.isfinite(denom):
        return ScalarStep(float("nan"), "non_spd", denom=denom)
    ctilde = (1.0 + float(dt)) * float(c) / denom
    if not (math.isfinite(ctilde) and ctilde > 0.0):
        return ScalarStep(float("nan"), "overflow", denom=denom, ctilde=ctilde)
    c_next = clip_scalar(ctilde, lambda_minus, lambda_plus)
    return ScalarStep(c_next, "ok", denom=denom, ctilde=ctilde)


def gaussian_kl_divergence(c_from: float, c_to: float) -> float:
    """``KL(N(0, c_from) || N(0, c_to))`` for scalar Gaussians.

    Equal to ``0.5 (r - 1 - log r)`` with ``r = c_from / c_to``; nonnegative and
    zero iff ``c_from == c_to``.
    """
    r = float(c_from) / float(c_to)
    return 0.5 * (r - 1.0 - math.log(r))


def clip_smoothness_constant(beta: float, lambda_minus: float,
                             lambda_plus: float) -> float:
    """Theorem 2.18 Bregman smoothness constant ``L_clip``.

    ``L_clip = beta * max(lambda_plus, lambda_plus^4 / lambda_minus^3)``.
    """
    beta = float(beta)
    lam_m = float(lambda_minus)
    lam_p = float(lambda_plus)
    return beta * max(lam_p, lam_p ** 4 / lam_m ** 3)


def theorem_safe_dt(beta: float, lambda_minus: float, lambda_plus: float,
                    safety: float = 0.9) -> float:
    """Theorem-safe stepsize ``dt = safety / L_clip`` for the clipped KL scheme."""
    return float(safety) / clip_smoothness_constant(beta, lambda_minus, lambda_plus)

