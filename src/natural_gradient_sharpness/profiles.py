"""Smooth cutoffs and their integrated flat-top bump profiles."""
from __future__ import annotations

import numpy as np
from scipy import integrate, optimize, special


def smooth_step(x):
    """A C-infinity step: zero on ``(-inf,-1]`` and one on ``[1,inf)``."""
    x = np.asarray(x, dtype=np.float64)
    out = np.zeros_like(x)
    out[x >= 1.0] = 1.0
    mask = (x > -1.0) & (x < 1.0)
    xm = x[mask]
    logit = -1.0 / (xm + 1.0) + 1.0 / (1.0 - xm)
    out[mask] = special.expit(logit)
    return out


def smooth_step_derivative(x):
    x = np.asarray(x, dtype=np.float64)
    out = np.zeros_like(x)
    mask = (x > -1.0) & (x < 1.0)
    xm = x[mask]
    h = smooth_step(xm)
    dlogit = 1.0 / (xm + 1.0) ** 2 + 1.0 / (1.0 - xm) ** 2
    out[mask] = h * (1.0 - h) * dlogit
    return out


def flat_top_bump(x):
    """Even C-infinity bump equal to one on ``|x|<=1/2`` and zero on ``|x|>=1``."""
    x = np.asarray(x, dtype=np.float64)
    return 1.0 - smooth_step(4.0 * np.abs(x) - 3.0)


def flat_top_bump_prime(x):
    x = np.asarray(x, dtype=np.float64)
    sign = np.sign(x)
    return -4.0 * sign * smooth_step_derivative(4.0 * np.abs(x) - 3.0)


class IntegratedFlatTop:
    """Numerical ``phi``, ``Phi=int phi`` and ``J=int Phi`` on the transition."""

    def __init__(self, grid_size: int = 20001):
        grid_size = int(grid_size)
        if grid_size % 2 == 0:
            grid_size += 1
        self.grid = np.linspace(-1.0, 1.0, grid_size)
        phi = flat_top_bump(self.grid)
        self.phi_grid = phi
        self.Phi_grid = integrate.cumulative_trapezoid(phi, self.grid, initial=0.0)
        self.J_grid = integrate.cumulative_trapezoid(self.Phi_grid, self.grid, initial=0.0)
        self.I_phi = float(self.Phi_grid[-1])
        self.J_at_one = float(self.J_grid[-1])
        result = optimize.minimize_scalar(
            lambda z: -float(abs(flat_top_bump_prime(np.array([z]))[0])),
            bounds=(0.5, 1.0), method="bounded", options={"xatol": 1e-14})
        self.M_phi = float(-result.fun)

    def phi(self, x):
        return flat_top_bump(x)

    def Phi(self, x):
        x = np.asarray(x, dtype=np.float64)
        return np.where(
            x <= -1.0, 0.0,
            np.where(x >= 1.0, self.I_phi,
                     np.interp(x, self.grid, self.Phi_grid)))

    def J(self, x):
        x = np.asarray(x, dtype=np.float64)
        return np.where(
            x <= -1.0, 0.0,
            np.where(x >= 1.0,
                     self.J_at_one + self.I_phi * (x - 1.0),
                     np.interp(x, self.grid, self.J_grid)))
