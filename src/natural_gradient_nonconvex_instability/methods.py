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

