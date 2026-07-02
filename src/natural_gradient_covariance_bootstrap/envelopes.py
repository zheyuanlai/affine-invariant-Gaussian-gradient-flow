"""Covariance lower envelopes and dynamic/frozen contraction factors.

The manuscript note replaces the *frozen* covariance lower bound
``lambda_frozen = min(lambda0, 1/beta)`` -- a constant valid for all time -- by a
*growing* lower envelope ``L_n`` that tracks the covariance burn-in. Both the KL
and the Riemannian Fisher--Rao covariance updates admit such an envelope with
``L_n <= lambda_min(C_n)``.

KL covariance envelope::

    L_0     = min(lambda0, 1/(2 beta))
    L_{n+1} = (1 + dt) L_n / (1 + dt beta L_n)
    L_n     = 1 / ( beta (1 + (1/(beta L_0) - 1) (1 + dt)^{-n}) )   (closed form)

The closed form solves the recurrence exactly; the bound ``lambda_min(C_n) >=
L_n`` holds for every ``dt > 0`` (per coordinate ``x/(1 + dt beta x)`` is
increasing and ``h_i <= beta``). ``L_n -> 1/beta`` geometrically, so the burn-in
to any fixed fraction of ``1/beta`` is ``O(log(1/(beta lambda0)))`` -- logarithmic
in ``1/lambda0``, not ``1/lambda0``.

Riemannian covariance envelope::

    L_0     = min(lambda0, 1/(2 beta))
    L_{n+1} = e^{dt} L_n / (1 + (e - 1) dt beta L_n)

(the rational lower bound to ``e^{dt} L e^{-dt beta L}``, valid for
``dt beta L <= 1``).

The dynamic and frozen per-step energy-gap contraction factors reuse the theorem
constants ``q_riem`` / ``q_kl`` from :mod:`src.common.theory_constants` with the
covariance lower bound set to the growing ``L_n`` (dynamic) or the constant
``lambda_frozen`` (frozen), and the covariance upper bound set to
``lambda_max_bound = max(lambda0_max, 1/alpha)``.
"""
from __future__ import annotations

import numpy as np

from src.common.theory_constants import q_riem, q_kl

E_CONST = float(np.e)


# ---------------------------------------------------------------------------
# Spectral constants
# ---------------------------------------------------------------------------

def envelope_L0(lambda0, beta):
    """Envelope start ``L_0 = min(lambda0, 1/(2 beta))``."""
    return float(min(float(lambda0), 1.0 / (2.0 * float(beta))))


def frozen_lower_bound(lambda0, beta):
    """Old frozen covariance lower bound ``lambda_frozen = min(lambda0, 1/beta)``."""
    return float(min(float(lambda0), 1.0 / float(beta)))


def lambda_max_bound(lambda0_max, alpha):
    """Covariance upper bound ``max(lambda0_max, 1/alpha)``."""
    return float(max(float(lambda0_max), 1.0 / float(alpha)))


# ---------------------------------------------------------------------------
# KL envelope
# ---------------------------------------------------------------------------

def kl_envelope_step(L, dt, beta):
    """One KL envelope step ``(1 + dt) L / (1 + dt beta L)``."""
    return (1.0 + dt) * L / (1.0 + dt * beta * L)


def kl_envelope_closed(n, L0, dt, beta):
    """Closed-form KL envelope ``L_n`` (accepts scalar or array ``n``)."""
    n = np.asarray(n, dtype=np.float64)
    x0 = 1.0 / (beta * L0)
    return 1.0 / (beta * (1.0 + (x0 - 1.0) * (1.0 + dt) ** (-n)))


def kl_envelope_sequence(n_steps, lambda0, dt, beta):
    """KL envelope ``[L_0, ..., L_{n_steps}]`` from the recurrence."""
    L = envelope_L0(lambda0, beta)
    out = [L]
    for _ in range(int(n_steps)):
        L = kl_envelope_step(L, dt, beta)
        out.append(L)
    return np.asarray(out, dtype=np.float64)


# ---------------------------------------------------------------------------
# Riemannian envelope
# ---------------------------------------------------------------------------

def riemannian_envelope_step(L, dt, beta):
    """One Riemannian envelope step ``e^{dt} L / (1 + (e - 1) dt beta L)``."""
    return np.exp(dt) * L / (1.0 + (E_CONST - 1.0) * dt * beta * L)


def riemannian_envelope_sequence(n_steps, lambda0, dt, beta):
    """Riemannian envelope ``[L_0, ..., L_{n_steps}]`` from the recurrence."""
    L = envelope_L0(lambda0, beta)
    out = [L]
    for _ in range(int(n_steps)):
        L = riemannian_envelope_step(L, dt, beta)
        out.append(L)
    return np.asarray(out, dtype=np.float64)


def envelope_sequence(scheme, n_steps, lambda0, dt, beta):
    """Covariance lower envelope for the named ``scheme``."""
    if scheme == "kl":
        return kl_envelope_sequence(n_steps, lambda0, dt, beta)
    if scheme == "riemannian":
        return riemannian_envelope_sequence(n_steps, lambda0, dt, beta)
    raise ValueError(f"unknown scheme '{scheme}' (known: riemannian, kl)")


# ---------------------------------------------------------------------------
# Dynamic / frozen contraction factors (reusing theory_constants)
# ---------------------------------------------------------------------------

def contraction_factor(scheme, lam_lower, dt, alpha, beta, lam_max):
    """Per-step energy-gap contraction factor for the named scheme.

    ``lam_lower`` is the covariance lower bound: the growing ``L_n`` gives the
    dynamic factor, the constant ``lambda_frozen`` gives the frozen factor.
    ``lam_max`` is ``lambda_max_bound``. Accepts a scalar or an array
    ``lam_lower`` (an array yields the per-iteration dynamic factor).
    """
    if scheme == "riemannian":
        return q_riem(dt, alpha, beta, lam_lower, lam_max)
    if scheme == "kl":
        return q_kl(dt, alpha, beta, lam_lower, lam_max)
    raise ValueError(f"unknown scheme '{scheme}' (known: riemannian, kl)")
