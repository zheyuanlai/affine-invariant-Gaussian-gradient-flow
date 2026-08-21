"""One-dimensional flat-top bump train for sharp global discrete rates."""
from __future__ import annotations

import math
import numpy as np

from .profiles import IntegratedFlatTop
from .quadrature import normal_rule


class BumpTrainTarget:
    """The fixed-step lower construction from Appendix B.2 of the paper.

    Expectations are exact up to one-dimensional Gauss--Hermite quadrature.
    Evaluation is O(log N): disjoint bump supports allow saturated bumps to be
    accumulated with prefix sums, leaving at most one active transition.
    """

    def __init__(self, kappa: int, step_fraction=0.5, train_time=0.15,
                 initial_scale=100.0, lipschitz=1.0, gh_order=32,
                 profile: IntegratedFlatTop | None = None):
        self.kappa = int(kappa)
        self.step_fraction = float(step_fraction)
        self.dt = self.step_fraction / self.kappa
        self.train_time = float(train_time)
        self.initial_scale = float(initial_scale)
        self.lipschitz = float(lipschitz)
        self.profile = profile or IntegratedFlatTop()
        self.width = (self.kappa - 1.0) * self.profile.M_phi / self.lipschitz
        self.bump_area = self.width * self.profile.I_phi
        self.score_jump = (self.kappa - 1.0) * self.bump_area
        self.num_train_steps = int(math.floor(
            self.train_time * self.kappa ** 2 / self.step_fraction))
        self.mean0 = self.initial_scale * self.kappa ** 4 / self.step_fraction
        centers = [self.mean0]
        for j in range(self.num_train_steps):
            B = centers[-1] + self.score_jump * (self.num_train_steps - j + 0.5)
            centers.append(centers[-1] - self.step_fraction * B / self.kappa ** 2)
        # Construction order is descending; evaluation order is ascending.
        self.centers_path = np.asarray(centers, dtype=np.float64)
        self.centers = self.centers_path[::-1].copy()
        self.prefix_centers = np.concatenate([[0.0], np.cumsum(self.centers)])
        self.cov0 = 1.0 / self.kappa
        self.gh_z, self.gh_w = normal_rule(gh_order)
        spacings = self.centers_path[:-1] - self.centers_path[1:]
        self.min_spacing_over_width = float(np.min(spacings) / self.width)

    def _positive_vp_vpp(self, x):
        x = np.asarray(x, dtype=np.float64)
        flat = x.ravel()
        vp = np.empty_like(flat)
        vpp = np.empty_like(flat)
        for idx, value in enumerate(flat):
            if value < 0.0:
                p, pp = self._positive_vp_vpp(np.array([-value]))
                vp[idx], vpp[idx] = -p[0], pp[0]
                continue
            left = int(np.searchsorted(self.centers, value - self.width, side="right"))
            right = int(np.searchsorted(self.centers, value + self.width, side="left"))
            full_count = left
            score = value + (self.kappa - 1.0) * self.width * self.profile.I_phi * full_count
            curvature = 1.0
            if right > left:
                # Disjointness guarantees at most one active support.
                center = self.centers[left]
                arg = (value - center) / self.width
                score += (self.kappa - 1.0) * self.width * float(self.profile.Phi(arg))
                curvature += (self.kappa - 1.0) * float(self.profile.phi(arg))
            vp[idx], vpp[idx] = score, curvature
        return vp.reshape(x.shape), vpp.reshape(x.shape)

    def potential(self, x):
        """Even potential normalized by ``V(0)=0``."""
        x = np.abs(np.asarray(x, dtype=np.float64))
        flat = x.ravel()
        out = np.empty_like(flat)
        I, J1, w = self.profile.I_phi, self.profile.J_at_one, self.width
        for idx, value in enumerate(flat):
            left = int(np.searchsorted(self.centers, value - w, side="right"))
            right = int(np.searchsorted(self.centers, value + w, side="left"))
            count = left
            # For u>=1: J(u)=J(1)+I(u-1).
            full = count * (J1 - I) + I * (count * value / w - self.prefix_centers[count] / w)
            active = 0.0
            if right > left:
                active = float(self.profile.J((value - self.centers[left]) / w))
            out[idx] = 0.5 * value ** 2 + (self.kappa - 1.0) * w ** 2 * (full + active)
        return out.reshape(x.shape)

    def expectations(self, mean, covariance):
        x = float(mean) + math.sqrt(float(covariance)) * self.gh_z
        vp, vpp = self._positive_vp_vpp(x)
        return float(self.gh_w @ vp), float(self.gh_w @ vpp)

    def objective(self, mean, covariance):
        x = float(mean) + math.sqrt(float(covariance)) * self.gh_z
        return float(self.gh_w @ self.potential(x) - 0.5 * math.log(covariance))

    @property
    def optimum_objective(self):
        # The remote bumps have zero quadrature mass at N(0,1); the target is
        # exactly quadratic throughout the resolved optimizer neighborhood.
        return 0.5

    def gap(self, mean, covariance):
        return self.objective(mean, covariance) - self.optimum_objective

    def rhs(self, _time, state):
        mean, covariance = map(float, state)
        b, A = self.expectations(mean, covariance)
        return np.array([-covariance * b, covariance - covariance ** 2 * A])

    def step(self, method, state, dt=None):
        mean, covariance = map(float, state)
        dt = self.dt if dt is None else float(dt)
        b, A = self.expectations(mean, covariance)
        mean_next = mean - dt * covariance * b
        if method == "riemannian":
            cov_next = covariance * math.exp(dt * (1.0 - covariance * A))
        elif method == "kl":
            cov_next = (1.0 + dt) * covariance / (1.0 + dt * covariance * A)
        else:
            raise ValueError(f"unknown method {method!r}")
        return np.array([mean_next, cov_next])

    @staticmethod
    def error_norm(state):
        mean, covariance = map(float, state)
        return math.sqrt(mean ** 2 + 0.5 * (covariance - 1.0) ** 2)
