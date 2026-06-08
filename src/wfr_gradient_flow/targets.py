"""Deterministic targets for the WFR Gaussian gradient-flow experiment.

Each target exposes, for the current variational state ``a = (m, C)``:

* :meth:`g_H(m, C)` -> ``(g, H)`` with ``g = E[grad log rho_post] = -E[grad V]``
  and ``H = E[grad^2 log rho_post] = -E[grad^2 V]`` (so ``log rho_post = -V`` up
  to a constant);
* :meth:`objective(m, C)` -> ``E(a) = E_{N(m,C)}[V] - 1/2 log det C`` (the KL up
  to an ``a``-independent constant);
* global strong-log-concavity / smoothness constants ``(alpha, beta)``.

Two production targets are implemented:

* :class:`GaussianTarget` -- the exact ill-conditioned Gaussian posterior
  ``rho_post = N(0, Gamma)``, ``Gamma = diag(1, Lambda)``. Here ``Gamma`` is the
  *covariance*; the precision is ``Gamma^{-1} = diag(1, 1/Lambda)`` so
  ``g = -Gamma^{-1} m`` and ``H = -Gamma^{-1}`` are exact and state-independent,
  ``alpha = 1/Lambda``, ``beta = 1``, and the optimum is closed form
  ``m_star = 0``, ``C_star = Gamma`` with an analytic energy gap.
* :class:`SmoothLogCoshTarget` -- the smooth nonseparable strongly log-concave
  posterior reused from the discretization experiments,
  ``V = (sqrt(Lambda) x - y)^2/20 + (delta/2)(x^2+y^2) + gamma log cosh(y)``,
  with ``alpha = delta`` and the *eigenvalue* smoothness constant
  ``beta = lambda_max([[delta+Lambda/10, -sqrt(Lambda)/10],
                       [-sqrt(Lambda)/10, delta+1/10+gamma]])``.

The anisotropy parameter is named ``Lambda`` throughout (the symbol ``lambda``
is reserved for the WFR transport strength ``lambda_n``). All expectations are
exact: closed form for the Gaussian target, deterministic Gauss--Hermite
quadrature (reused from the discretization module) for the ``log cosh`` term.
"""
from __future__ import annotations

import numpy as np

from src.common.spd import symmetrize
# Reuse the deterministic 1-D Gauss--Hermite quadrature from the discretization
# experiment (single source of truth for the smooth-target expectations).
from src.natural_gradient_discretization_stepsize.targets import (
    gaussian_expectation_1d, _sech2,
)

TARGET_NAMES = ["gaussian", "smooth_log_cosh"]


# ---------------------------------------------------------------------------
# Initialization shared by both targets
# ---------------------------------------------------------------------------

def flat_valley_init(Lambda, epsilon, r=8.0):
    """Far, underdispersed initial state ``(m0, C0)``.

    ``C0 = epsilon I`` (underdispersed) and the mean is placed a distance ``r``
    along the weak-curvature / valley direction
    ``(1, sqrt(Lambda)) / sqrt(1 + Lambda)``. For the Gaussian target this is the
    broad direction of ``Gamma = diag(1, Lambda)``; for the smooth target it
    aligns with the valley of the coupled term ``sqrt(Lambda) x - y``.
    """
    Lambda = float(Lambda)
    direction = np.array([1.0, np.sqrt(Lambda)], dtype=np.float64)
    direction /= np.sqrt(1.0 + Lambda)
    m0 = float(r) * direction
    C0 = float(epsilon) * np.eye(2, dtype=np.float64)
    return m0, C0


# ---------------------------------------------------------------------------
# Target 1: exact ill-conditioned Gaussian posterior
# ---------------------------------------------------------------------------

class GaussianTarget:
    """Exact Gaussian posterior ``rho_post = N(0, Gamma)``, ``Gamma=diag(1,Lambda)``.

    ``V(theta) = 1/2 theta^T Gamma^{-1} theta`` so ``grad V = Gamma^{-1} theta``
    and ``grad^2 V = Gamma^{-1}`` are exact and state-independent:

        g = -E[grad V] = -Gamma^{-1} m,    H = -E[grad^2 V] = -Gamma^{-1}.

    The optimum is exact, ``m_star = 0``, ``C_star = Gamma``, and the energy gap
    equals the exact Gaussian KL.
    """

    name = "gaussian"
    d = 2
    has_theory = True

    def __init__(self, Lambda, epsilon, r=8.0):
        self.Lambda = float(Lambda)
        self.epsilon = float(epsilon)
        self.r = float(r)
        self.Gamma = np.diag([1.0, self.Lambda]).astype(np.float64)
        self.Gamma_inv = np.diag([1.0, 1.0 / self.Lambda]).astype(np.float64)
        self.m0, self.C0 = flat_valley_init(self.Lambda, self.epsilon, self.r)
        # grad^2 V = Gamma^{-1} = diag(1, 1/Lambda) globally.
        self.alpha = float(min(1.0, 1.0 / self.Lambda))
        self.beta = float(max(1.0, 1.0 / self.Lambda))

    def g_H(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        g = -(self.Gamma_inv @ m)
        H = -self.Gamma_inv.copy()
        return g, H

    def expected_V(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        return 0.5 * (float(m @ self.Gamma_inv @ m)
                      + float(np.trace(self.Gamma_inv @ C)))

    def objective(self, m, C):
        C = symmetrize(C)
        _, logdet = np.linalg.slogdet(C)
        return -0.5 * float(logdet) + self.expected_V(m, C)

    def star(self):
        m_star = np.zeros(2, dtype=np.float64)
        C_star = self.Gamma.copy()
        return m_star, symmetrize(C_star), self.objective(m_star, C_star)

    def energy_gap(self, m, C):
        """Analytic ``E(a) - E(a_star) = KL(N(m,C) || N(0,Gamma))``."""
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        GiC = self.Gamma_inv @ C
        _, logdet = np.linalg.slogdet(GiC)
        return 0.5 * (float(m @ self.Gamma_inv @ m) + float(np.trace(GiC))
                      - float(logdet) - self.d)

    def metadata(self):
        return {
            "target_name": self.name, "Lambda": self.Lambda, "d": self.d,
            "epsilon": self.epsilon, "r": self.r,
            "Gamma_diag": [1.0, self.Lambda],
            "m0": self.m0.tolist(), "C0": self.C0.tolist(),
            "has_theory": True, "alpha": self.alpha, "beta": self.beta,
            "beta_formula": "max eig(Gamma^{-1}) = max(1, 1/Lambda) = 1",
            "description": "exact Gaussian posterior N(0, Gamma), Gamma=diag(1,Lambda)",
        }


# ---------------------------------------------------------------------------
# Target 2: smooth nonseparable strongly log-concave posterior
# ---------------------------------------------------------------------------

def smooth_beta(Lambda, delta, gamma):
    """Eigenvalue smoothness constant of the smooth log-cosh target.

    ``beta = lambda_max(M)`` with the worst-case Hessian (``sech^2 = 1``)

        M = [[delta + Lambda/10,   -sqrt(Lambda)/10],
             [-sqrt(Lambda)/10,    delta + 1/10 + gamma]].
    """
    a = delta + Lambda / 10.0
    c = delta + 1.0 / 10.0 + gamma
    b = -np.sqrt(Lambda) / 10.0
    return float(0.5 * (a + c) + np.sqrt((0.5 * (a - c)) ** 2 + b * b))


class SmoothLogCoshTarget:
    """Smooth strongly log-concave nonseparable posterior.

        V(x,y) = (sqrt(Lambda) x - y)^2/20 + (delta/2)(x^2+y^2)
                 + gamma log cosh(y),    delta=0.05, gamma=1.

    Globally strongly log-concave and smooth with ``alpha = delta`` and the
    eigenvalue smoothness constant :func:`smooth_beta`. Expectations of
    ``tanh(Y)``, ``sech^2(Y)`` and ``log cosh(Y)`` under ``Y ~ N(m_y, C_22)`` are
    evaluated by deterministic Gauss--Hermite quadrature. The optimum is found
    numerically (no closed form), so this class exposes no ``star`` method.
    """

    name = "smooth_log_cosh"
    d = 2
    has_theory = True

    def __init__(self, Lambda, epsilon, r=8.0, delta=0.05, gamma=1.0, n_nodes=80):
        self.Lambda = float(Lambda)
        self.epsilon = float(epsilon)
        self.r = float(r)
        self.sl = float(np.sqrt(self.Lambda))
        self.delta = float(delta)
        self.gamma = float(gamma)
        self.n_nodes = int(n_nodes)
        self.m0, self.C0 = flat_valley_init(self.Lambda, self.epsilon, self.r)
        self.alpha = self.delta
        self.beta = smooth_beta(self.Lambda, self.delta, self.gamma)

    def g_H(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        mx, my = float(m[0]), float(m[1])
        C22 = float(C[1, 1])
        lam, sl, delta, gamma = self.Lambda, self.sl, self.delta, self.gamma
        E_tanh = gaussian_expectation_1d(np.tanh, my, C22, self.n_nodes)
        E_sech2 = gaussian_expectation_1d(_sech2, my, C22, self.n_nodes)
        gradx = (lam * mx - sl * my) / 10.0 + delta * mx
        grady = (-sl * mx + my) / 10.0 + delta * my + gamma * E_tanh
        g = -np.array([gradx, grady], dtype=np.float64)
        H = -np.array([
            [lam / 10.0 + delta,              -sl / 10.0],
            [-sl / 10.0, 1.0 / 10.0 + delta + gamma * E_sech2],
        ], dtype=np.float64)
        return g, symmetrize(H)

    def expected_V(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        mx, my = float(m[0]), float(m[1])
        C11, C12, C22 = float(C[0, 0]), float(C[0, 1]), float(C[1, 1])
        lam, sl, delta, gamma = self.Lambda, self.sl, self.delta, self.gamma
        mean_U = sl * mx - my
        var_U = lam * C11 - 2.0 * sl * C12 + C22
        E_U2 = mean_U ** 2 + var_U
        E_quad = mx ** 2 + C11 + my ** 2 + C22
        E_logcosh = gaussian_expectation_1d(
            lambda y: np.log(np.cosh(y)), my, C22, self.n_nodes)
        return E_U2 / 20.0 + 0.5 * delta * E_quad + gamma * E_logcosh

    def objective(self, m, C):
        C = symmetrize(C)
        _, logdet = np.linalg.slogdet(C)
        return -0.5 * float(logdet) + self.expected_V(m, C)

    def metadata(self):
        return {
            "target_name": self.name, "Lambda": self.Lambda, "d": self.d,
            "epsilon": self.epsilon, "r": self.r,
            "delta": self.delta, "gamma": self.gamma, "gh_nodes": self.n_nodes,
            "m0": self.m0.tolist(), "C0": self.C0.tolist(),
            "has_theory": True, "alpha": self.alpha, "beta": self.beta,
            "beta_formula": ("lambda_max([[delta+Lambda/10, -sqrt(Lambda)/10], "
                             "[-sqrt(Lambda)/10, delta+1/10+gamma]])"),
            "description": ("smooth strongly log-concave "
                            "V=(sqrt(Lambda)x-y)^2/20 + delta/2 |theta|^2 + gamma logcosh(y)"),
        }


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

_TARGET_CLASSES = {
    "gaussian": GaussianTarget,
    "smooth_log_cosh": SmoothLogCoshTarget,
}


def build_target(name, Lambda, epsilon, **kwargs):
    """Construct a WFR target by ``name``, anisotropy ``Lambda`` and ``epsilon``."""
    if name not in _TARGET_CLASSES:
        raise ValueError(f"unknown target '{name}' (known: {sorted(_TARGET_CLASSES)})")
    return _TARGET_CLASSES[name](Lambda, epsilon, **kwargs)
