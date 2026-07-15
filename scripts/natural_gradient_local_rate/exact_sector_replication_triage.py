"""Deterministic exact-sector quadrature for replicated radial-softplus modules.

This script is a triage tool, not a proof certificate.  O(d)-symmetry reduces
the Gaussian expectations to T ~ N(0, 1) and S ~ chi-square(d).  Tensor-product
Gauss-Hermite/generalized-Laguerre quadrature then checks:

* the exact shared-longitudinal isotropy equation J delta^2 Q = M;
* the normalized first-Hermite transverse witness;
* the affine Rayleigh witness used in the bracket report.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import brentq
from scipy.special import gamma, roots_genlaguerre


def gaussian_sectors(dimension: int, nodes: int = 64):
    """Return quadrature nodes for T and S with normalized product weights."""
    hermite_x, hermite_w = hermgauss(nodes)
    t = math.sqrt(2.0) * hermite_x
    t_weights = hermite_w / math.sqrt(math.pi)

    shape = dimension / 2.0
    laguerre_x, laguerre_w = roots_genlaguerre(nodes, shape - 1.0)
    s_norm = 2.0 * laguerre_x
    s_weights = laguerre_w / gamma(shape)

    weights = t_weights[:, None] * s_weights[None, :]
    return t[:, None], s_norm[None, :], weights


def sector_moments(
    dimension: int,
    sigma: float,
    delta: float,
    sectors,
):
    """Evaluate Q, M, and the transverse coupling core at fixed delta."""
    t, s_norm, weights = sectors
    radius = np.sqrt(sigma * sigma + s_norm)
    big_r = math.sqrt(sigma * sigma + dimension)
    phase = (radius - big_r - delta * t) / sigma

    # Stable logistic evaluation.
    p = np.where(
        phase >= 0.0,
        1.0 / (1.0 + np.exp(-phase)),
        np.exp(phase) / (1.0 + np.exp(phase)),
    )
    q = p * (1.0 - p) / sigma

    q_mean = float(np.sum(weights * q))
    transverse_mean = float(
        np.sum(
            weights
            * (
                q * s_norm / radius**2
                + p
                * (
                    (dimension - 1.0) / radius
                    + sigma * sigma / radius**3
                )
            )
        )
        / dimension
    )
    coupling_core = float(np.sum(weights * q * s_norm / radius))
    return q_mean, transverse_mean, coupling_core, (
        p,
        q,
        radius,
        s_norm,
        weights,
    )


def evaluate(dimension: int, sigma: float, replicas: int, nodes: int = 64):
    """Solve exact isotropy and return the audited witness quantities."""
    sectors = gaussian_sectors(dimension, nodes)

    def isotropy_residual(delta: float) -> float:
        q_mean, transverse_mean, *_ = sector_moments(
            dimension, sigma, delta, sectors
        )
        return replicas * delta * delta * q_mean - transverse_mean

    delta = brentq(isotropy_residual, 0.0, 1.0)
    q_mean, transverse_mean, coupling_core, auxiliary = sector_moments(
        dimension, sigma, delta, sectors
    )
    p, q, radius, s_norm, weights = auxiliary

    big_r = math.sqrt(sigma * sigma + dimension)
    alpha = sigma / big_r
    coefficient = (1.0 - alpha) / transverse_mean

    coupling = (
        coefficient
        * delta
        * math.sqrt(replicas / dimension)
        * coupling_core
    )

    lam = 2.0 * delta * big_r / dimension
    residual = s_norm / radius - dimension / big_r
    one_module_energy = float(
        np.sum(
            weights
            * (
                q * lam * lam * residual**2
                + p
                * lam
                * lam
                * sigma
                * sigma
                * s_norm
                / radius**3
            )
        )
    )
    x_squared = replicas * dimension * lam * lam
    four_q = (
        alpha * (4.0 + x_squared)
        + coefficient * replicas * one_module_energy
        + x_squared
    )
    fisher_squared = 1.0 + 0.5 * x_squared
    rayleigh = 0.25 * four_q / fisher_squared

    return {
        "d": dimension,
        "sigma": sigma,
        "J": replicas,
        "delta2": delta * delta,
        "J_delta2": replicas * delta * delta,
        "Q": q_mean,
        "M": transverse_mean,
        "tau_witness2": coupling * coupling,
        "rayleigh": rayleigh,
        "sigma_over_R": alpha,
    }


def main() -> None:
    for n in (2, 3):
        dimension = n**4
        for replicas in (1, 4, 16):
            row = evaluate(dimension, float(n), replicas)
            print(" ".join(f"{key}={value:.10g}" for key, value in row.items()))


if __name__ == "__main__":
    main()
