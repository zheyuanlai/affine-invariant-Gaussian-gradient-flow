"""Symmetry-reduced local sharpness constructions.

This module implements the three local mechanisms appearing in the paper:

* the one-dimensional smoothed bump of Chen et al. (inverse gap ``log kappa``);
* the gamma shell (inverse gap ``sqrt(n log kappa)`` and saturation);
* the convex ridge (inverse gap ``sqrt(kappa)`` with uniform Hessian Lipschitzness).

The shell and ridge never allocate an ambient-dimensional vector.  O(m)
symmetry reduces the full Fisher--Rao Hessian to scalar, vector, and traceless
irreducible blocks of sizes 3, 2, and 1.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy import integrate, optimize, stats

from .profiles import smooth_step, smooth_step_derivative
from .quadrature import RadialGaussianRule


SQRT2 = math.sqrt(2.0)


def _ridge_fp_fpp(x, K, center):
    """Triangular ``F''`` and its exact integral ``F'``."""
    d = np.asarray(x, dtype=np.float64) - float(center)
    fpp = np.maximum(float(K) - np.abs(d), 0.0)
    fp = np.zeros_like(d)
    left = (d > -K) & (d <= 0.0)
    fp[left] = 0.5 * (d[left] + K) ** 2
    right = (d > 0.0) & (d < K)
    fp[right] = 0.5 * K ** 2 + K * d[right] - 0.5 * d[right] ** 2
    fp[d >= K] = K ** 2
    return fp, fpp


def _max_hessian_increment(m, sigma, fp_fpp, profile_grid):
    """Numerically maximize the O(m)-symmetric Hessian operator norm.

    The profile variable and ``rho`` are independent because the distinguished
    coordinate can shift the ridge/shell argument for every ``rho``.  A dense
    deterministic grid also includes the ``rho -> infinity`` rank-one limit.
    """
    rho_min = math.sqrt(m)
    rho = rho_min * np.geomspace(1.0, 1e4, 260)
    a = np.sqrt(np.maximum(1.0 - m / rho ** 2, 0.0))
    best = 0.0
    for x in np.asarray(profile_grid, dtype=np.float64):
        fp, fpp = fp_fpp(np.array([x]))
        fp, fpp = float(fp[0]), float(fpp[0])
        aa = fpp * sigma ** 2
        bb = fpp * sigma * a
        dd = fpp * a ** 2 + fp * m / rho ** 3
        lam = 0.5 * (aa + dd + np.sqrt((aa - dd) ** 2 + 4.0 * bb ** 2))
        tangential = fp / rho
        best = max(best, float(np.max(lam)), float(np.max(tangential)))
        best = max(best, fpp * (1.0 + sigma ** 2))
    return best


def invariant_blocks(target):
    """Return metric-orthonormal scalar/vector/traceless generator blocks."""
    m = target.m
    rule = target.rule
    T, U, W = rule.t, rule.u, rule.weights
    rho = np.sqrt(m + U)
    fp, fpp = target.fp_fpp(target.profile_argument(rho, T))
    alpha, scale, sigma = target.alpha, target.scale, target.sigma

    # Scalar sector: mean_t, covariance_tt, normalized transverse trace.
    wt = [2.0 + 0.0 * T + 0.0 * U, SQRT2 * T + 0.0 * U, 0.0 * T + 0.0 * U]
    d = [
        2.0 * sigma + 0.0 * T + 0.0 * U,
        SQRT2 * sigma * T + 0.0 * U,
        SQRT2 * U / (math.sqrt(m) * rho) + 0.0 * T,
    ]
    scalar = np.empty((3, 3), dtype=np.float64)
    for i in range(3):
        for j in range(i, 3):
            euclid = wt[i] * wt[j]
            if i == 2 and j == 2:
                euclid = euclid + 2.0 * U / m
            tang = 0.0
            if i == 2 and j == 2:
                tang = 2.0 * U / rho ** 3
            bform = alpha * euclid + scale * (fpp * d[i] * d[j] + fp * tang)
            extra = 0.5 if i == j and i > 0 else 0.0
            scalar[i, j] = scalar[j, i] = 0.25 * float(np.sum(W * bform)) + extra

    # Vector sector: mean_y and t-y cross covariance (multiplicity m).
    htan = 1.0 / rho - U / (m * rho ** 3)
    # Conditional angular averaging is applied analytically below.
    vector = np.empty((2, 2), dtype=np.float64)
    rank00 = 4.0 * U / (m * rho ** 2)
    rank01 = 2.0 * U / (m * rho) * (sigma + T / rho)
    rank11 = U / m * (sigma + T / rho) ** 2
    e00 = 4.0 + 0.0 * T + 0.0 * U
    e01 = 2.0 * T + 0.0 * U
    e11 = T ** 2 + U / m
    tan00, tan01, tan11 = 4.0 * htan, 2.0 * T * htan, T ** 2 * htan
    ranks = [[rank00, rank01], [rank01, rank11]]
    eus = [[e00, e01], [e01, e11]]
    tans = [[tan00, tan01], [tan01, tan11]]
    for i in range(2):
        for j in range(i, 2):
            bform = alpha * eus[i][j] + scale * (fpp * ranks[i][j] + fp * tans[i][j])
            extra = 0.5 if i == j == 1 else 0.0
            vector[i, j] = vector[j, i] = 0.25 * float(np.sum(W * bform)) + extra

    # Transverse traceless covariance sector.
    angular2 = 4.0 * U ** 2 / (m * (m + 2.0))
    euclid = 2.0 * U / m
    rank = angular2 / rho ** 2
    tang = 2.0 * U / (m * rho) - angular2 / rho ** 3
    bform = alpha * euclid + scale * (fpp * rank + fp * tang)
    traceless = 0.25 * float(np.sum(W * bform)) + 0.5

    return {
        "scalar": 0.5 * (scalar + scalar.T),
        "vector": 0.5 * (vector + vector.T),
        "traceless": np.array([[traceless]], dtype=np.float64),
    }


def block_spectrum(blocks):
    rows = []
    for name, block in blocks.items():
        for value in np.linalg.eigvalsh(block):
            rows.append((float(value), name))
    rows.sort()
    return rows


def kl_jacobian(block, dt, mean_coordinates=1):
    """KL/Bregman derivative from the continuous generator block.

    In star-metric orthonormal coordinates the covariance rows use effective
    step ``dt/(1+dt)`` while mean rows use ``dt``.
    """
    d = np.full(block.shape[0], float(dt) / (1.0 + float(dt)))
    d[:int(mean_coordinates)] = float(dt)
    return np.eye(block.shape[0]) - d[:, None] * block


class RadialLocalTarget:
    """Common reduced dynamics for the ridge and shell targets."""

    def fp_fpp(self, x):  # pragma: no cover - abstract protocol
        raise NotImplementedError

    def profile_argument(self, rho, standard_t):
        return rho + self.sigma * standard_t

    def expected_fields(self, mean_t, cov_t, cov_y):
        T, U, W = self.rule.t, self.rule.u, self.rule.weights
        x_t = float(mean_t) + math.sqrt(float(cov_t)) * T
        r2 = float(cov_y) * U
        rho = np.sqrt(self.m + r2)
        fp, fpp = self.fp_fpp(self.profile_argument(rho, x_t))
        efp = float(np.sum(W * fp))
        h_tt = self.alpha + self.scale * self.sigma ** 2 * float(np.sum(W * fpp))
        h_yy_u = (
            fpp * r2 / rho ** 2
            + fp * (self.m / rho - r2 / rho ** 3)
        ) / self.m
        h_yy = self.alpha + self.scale * float(np.sum(W * h_yy_u))
        grad_t = self.alpha * float(mean_t) + self.scale * self.sigma * efp - self.ell_t
        return grad_t, h_tt, h_yy

    def reduced_rhs(self, _time, state):
        mean_t, cov_t, cov_y = map(float, state)
        grad_t, h_tt, h_yy = self.expected_fields(mean_t, cov_t, cov_y)
        return np.array([
            -cov_t * grad_t,
            cov_t - cov_t ** 2 * h_tt,
            cov_y - cov_y ** 2 * h_yy,
        ])

    def reduced_step(self, method, state, dt):
        mean_t, cov_t, cov_y = map(float, state)
        grad_t, h_tt, h_yy = self.expected_fields(mean_t, cov_t, cov_y)
        mean_next = mean_t - dt * cov_t * grad_t
        if method == "riemannian":
            ct_next = cov_t * math.exp(dt * (1.0 - cov_t * h_tt))
            cy_next = cov_y * math.exp(dt * (1.0 - cov_y * h_yy))
        elif method == "kl":
            ct_next = (1.0 + dt) * cov_t / (1.0 + dt * cov_t * h_tt)
            cy_next = (1.0 + dt) * cov_y / (1.0 + dt * cov_y * h_yy)
        else:
            raise ValueError(f"unknown method {method!r}")
        return np.array([mean_next, ct_next, cy_next])

    def state_from_scalar_coordinates(self, coordinates, amplitude):
        v = np.asarray(coordinates, dtype=np.float64)
        eps = float(amplitude)
        return np.array([
            eps * v[0],
            1.0 + eps * SQRT2 * v[1],
            1.0 + eps * SQRT2 * v[2] / math.sqrt(self.m),
        ])

    def scalar_coordinates(self, state):
        mean_t, cov_t, cov_y = np.asarray(state, dtype=np.float64)
        return np.array([
            mean_t,
            (cov_t - 1.0) / SQRT2,
            math.sqrt(self.m) * (cov_y - 1.0) / SQRT2,
        ])


class RidgeTarget(RadialLocalTarget):
    """Optimizer-whitened convex-ridge family from the sharpness note."""

    construction = "ridge"

    def __init__(self, K: float, normal_order=28, radial_order=14):
        self.K = float(K)
        m_float = self.K ** 4 / 8.0
        self.m = int(round(m_float))
        if abs(self.m - m_float) > 1e-8:
            raise ValueError("K must make K^4/8 an integer")
        self.n = self.m + 1
        self.rho0 = self.K ** 2 / 2.0
        self.center = self.rho0 + 0.5
        self.rule = RadialGaussianRule(self.m, normal_order, radial_order)

        def phi(sigma):
            rho = np.sqrt(self.m + self.rule.u)
            fp, fpp = _ridge_fp_fpp(rho + sigma * self.rule.t, self.K, self.center)
            kbar = self.rule.mean(fpp)
            G = self.rule.mean((
                fpp * self.rule.u / rho ** 2
                + fp * (self.m / rho - self.rule.u / rho ** 3)
            ) / self.m)
            return sigma ** 2 * kbar - G

        hi = 2.0 / math.sqrt(self.K)
        self.sigma = float(optimize.brentq(phi, 0.0, hi, xtol=2e-13, rtol=2e-13))
        rho = np.sqrt(self.m + self.rule.u)
        fp, fpp = self.fp_fpp(rho + self.sigma * self.rule.t)
        self.kbar = self.rule.mean(fpp)
        self.a = self.sigma ** 2 * self.kbar
        self.alpha = 1.0 - self.a
        self.scale = 1.0
        self.ell_t = self.sigma * self.rule.mean(fp)
        g0, htt0, hyy0 = self.expected_fields(0.0, 1.0, 1.0)
        self.isotropy_error = max(abs(g0), abs(htt0 - 1.0), abs(hyy0 - 1.0))
        grid = np.linspace(self.center - self.K, self.center + self.K, 401)
        inc = _max_hessian_increment(
            self.m, self.sigma, self.fp_fpp, grid)
        self.beta = self.alpha + inc
        self.kappa = self.beta / self.alpha
        self.hessian_lipschitz_certificate = 64.0
        self.blocks = invariant_blocks(self)
        spectrum = block_spectrum(self.blocks)
        self.gamma, self.slow_block = spectrum[0]
        self.Lambda = spectrum[-1][0]

    def fp_fpp(self, x):
        return _ridge_fp_fpp(x, self.K, self.center)


class ShellTarget(RadialLocalTarget):
    """Optimizer-whitened gamma-shell family."""

    construction = "shell"

    def __init__(self, m: int, depth: float, a0=0.2,
                 normal_order=28, radial_order=14):
        self.m = int(m)
        self.n = self.m + 1
        self.depth = float(depth)
        if not (1.0 <= self.depth <= self.m):
            raise ValueError("shell depth must satisfy 1 <= depth <= m")
        f = lambda delta: 0.5 * (delta - self.m * math.log1p(delta / self.m)) - self.depth
        self.Delta = float(optimize.brentq(f, 1e-14, 20.0 * (math.sqrt(self.m * self.depth) + self.depth)))
        self.u0 = self.m + self.Delta
        self.rho0 = math.sqrt(2.0 * self.m + self.Delta)
        logp = float(stats.chi2.logsf(self.u0, self.m))
        self.tail_probability = math.exp(logp)
        log_g0 = math.log(2.0 * self.rho0) + float(stats.chi2.logpdf(self.u0, self.m))
        self.eta = math.exp(log_g0 - logp)
        self.rule = RadialGaussianRule(
            self.m, normal_order, radial_order, self.tail_probability)

        def root_value(sigma):
            rho = np.sqrt(self.m + self.rule.u)
            fp, fpp = self.fp_fpp(rho - self.rho0 + sigma * self.rule.t)
            Q = self.rule.mean(fpp)
            M = self.rule.mean((
                fpp * self.rule.u / rho ** 2
                + fp * (self.m / rho - self.rule.u / rho ** 3)
            ) / self.m)
            return sigma ** 2 * Q - M

        center = 1.0 / math.sqrt(self.Delta)
        scan = center * np.geomspace(0.02, 50.0, 80)
        vals = [root_value(s) for s in scan]
        bracket = None
        for left, right, vl, vr in zip(scan[:-1], scan[1:], vals[:-1], vals[1:]):
            if vl == 0.0 or vl * vr < 0.0:
                bracket = (left, right)
                break
        if bracket is None:
            raise RuntimeError("failed to bracket shell isotropization root")
        self.sigma = float(optimize.brentq(root_value, *bracket, xtol=2e-12, rtol=2e-12))
        rho = np.sqrt(self.m + self.rule.u)
        fp, fpp = self.fp_fpp(rho - self.rho0 + self.sigma * self.rule.t)
        self.Q = self.rule.mean(fpp)
        self.M = self.rule.mean((
            fpp * self.rule.u / rho ** 2
            + fp * (self.m / rho - self.rule.u / rho ** 3)
        ) / self.m)
        self.alpha = float(a0) / self.Delta
        self.scale = (1.0 - self.alpha) / self.M
        self.ell_t = self.scale * self.sigma * self.rule.mean(fp)
        g0, htt0, hyy0 = self.expected_fields(0.0, 1.0, 1.0)
        self.isotropy_error = max(abs(g0), abs(htt0 - 1.0), abs(hyy0 - 1.0))
        profile_grid = np.linspace(-1.0 / self.eta, 1.0 / self.eta, 401)
        inc = _max_hessian_increment(self.m, self.sigma, self.fp_fpp, profile_grid)
        self.beta = self.alpha + self.scale * inc
        self.kappa = self.beta / self.alpha
        self.blocks = invariant_blocks(self)
        spectrum = block_spectrum(self.blocks)
        self.gamma, self.slow_block = spectrum[0]
        self.Lambda = spectrum[-1][0]

    def fp_fpp(self, x):
        z = self.eta * np.asarray(x, dtype=np.float64)
        return smooth_step(z), self.eta * smooth_step_derivative(z)

    def profile_argument(self, rho, standard_t):
        return rho - self.rho0 + self.sigma * standard_t


@dataclass
class ChenLogBump:
    """One-dimensional sharp logarithmic example of Chen et al., Thm. 5.7."""

    bump_location: float

    def __post_init__(self):
        g = float(self.bump_location)
        if not (0.0 < g < math.exp(-1.0)):
            raise ValueError("bump_location must lie in (0,e^-1)")
        self.sigma = g ** 1.5
        self.C_star = -g ** 2 / (2.0 * math.log(g)) - g ** 3
        self.alpha = 1.0 / ((-math.log(g)) * self.C_star)
        S = self.C_star + self.sigma ** 2
        lo, hi = g - self.sigma, g + self.sigma
        p = stats.norm.cdf(hi / math.sqrt(S)) - stats.norm.cdf(lo / math.sqrt(S))
        self.beta = self.alpha + (1.0 / self.C_star - self.alpha) / p
        self.kappa = self.beta / self.alpha
        cvar = self.sigma ** 2 * self.C_star / S
        coef = self.C_star / S

        def truncated_moment(power):
            integrand = lambda x: (
                (coef * x if power == 1 else cvar + (coef * x) ** 2)
                * stats.norm.pdf(x, scale=math.sqrt(S))
            )
            return integrate.quad(integrand, lo, hi, epsabs=1e-14, epsrel=1e-12, limit=200)[0]

        A1 = (self.beta - self.alpha) * truncated_moment(1)
        A2 = self.alpha * self.C_star + (self.beta - self.alpha) * truncated_moment(2)
        T_whitened = math.sqrt(self.C_star) * A1
        self.block = np.array([
            [1.0, T_whitened / SQRT2],
            [T_whitened / SQRT2, 0.5 * (1.0 + A2)],
        ])
        eig = np.linalg.eigvalsh(self.block)
        self.gamma = float(eig[0])
        self.Lambda = float(eig[-1])
        self.inverse_gap_scale = math.log(math.e * self.kappa)
