"""Deterministic targets for the discretization stepsize experiment.

Each target exposes, for the current variational state ``a = (m, C)``:

* :meth:`g_H(m, C)` -> ``(g, H)`` with ``g = E[grad log rho_post] = -E[grad V]``
  and ``H = E[grad^2 log rho_post] = -E[grad^2 V]``, the drift fields of the
  natural gradient flow ``dm/dt = C g``, ``dC/dt = C + C H C``;
* :meth:`objective(m, C)` -> ``F = -1/2 log det C + E_{N(m,C)}[V]`` (the KL up to
  an additive constant independent of ``a``);
* the initial condition ``(m0, C0)``;
* analytic or numerically-optimized ``(m_star, C_star, F_star)``;
* global strong-convexity / smoothness constants ``(alpha, beta)`` when the
  target is globally smooth (``has_theory``).

All expectations are exact: closed form for the Gaussian and quartic targets,
deterministic Gauss--Hermite quadrature for the ``log cosh`` smooth target. There
is no Monte Carlo anywhere in this module.
"""
from __future__ import annotations

import numpy as np

from src.common.spd import symmetrize

# Lambda values used by every two-dimensional target.
LAMBDA_GRID = [0.01, 0.1, 1.0]

TARGET_NAMES = ["gaussian_posterior", "literature_logconcave", "smooth_logconcave"]


# ---------------------------------------------------------------------------
# Gauss--Hermite quadrature for 1-D Gaussian expectations
# ---------------------------------------------------------------------------

def gauss_hermite_nodes(n_nodes=80):
    """Nodes/weights for ``E_{Z~N(0,1)}[f(Z)] = sum_k w_k f(x_k)`` (``sum w_k = 1``)."""
    x, w = np.polynomial.hermite.hermgauss(int(n_nodes))   # weight exp(-x^2)
    return np.sqrt(2.0) * x, w / np.sqrt(np.pi)


def gaussian_expectation_1d(f, mean, var, n_nodes=80):
    """``E[f(Y)]`` for ``Y ~ N(mean, var)`` via Gauss--Hermite quadrature."""
    nodes, weights = gauss_hermite_nodes(n_nodes)
    y = mean + np.sqrt(max(var, 0.0)) * nodes
    return float(np.dot(weights, f(y)))


def _sech2(z):
    """Numerically stable ``sech^2(z) = 1 - tanh^2(z)``."""
    t = np.tanh(z)
    return 1.0 - t * t


# ---------------------------------------------------------------------------
# Target A: exact Gaussian posterior
# ---------------------------------------------------------------------------

class GaussianPosteriorTarget:
    """Gaussian posterior ``rho_post = N(0, Q^{-1})`` with ``Q = diag(1, lambda)``.

    ``V(theta) = 1/2 theta^T Q theta`` so ``grad V = Q theta`` and
    ``grad^2 V = Q`` are exact and state-independent:

        g = -E[grad V] = -Q m,      H = -E[grad^2 V] = -Q.

    The optimal variational Gaussian is exact, ``m_star = 0``, ``C_star = Q^{-1}``,
    and the energy gap is analytic.
    """

    name = "gaussian_posterior"
    d = 2
    has_theory = True

    def __init__(self, lam):
        self.lam = float(lam)
        self.Q = np.diag([1.0, self.lam]).astype(np.float64)
        self.m0 = np.array([10.0, 10.0], dtype=np.float64)
        self.C0 = np.diag([0.5, 2.0]).astype(np.float64)
        # grad^2 V = Q globally => alpha = min eig(Q), beta = max eig(Q).
        self.alpha = float(min(1.0, self.lam))
        self.beta = float(max(1.0, self.lam))

    def g_H(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        g = -(self.Q @ m)
        H = -self.Q.copy()
        return g, H

    def expected_V(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        return 0.5 * (float(m @ self.Q @ m) + float(np.trace(self.Q @ C)))

    def objective(self, m, C):
        C = symmetrize(C)
        sign, logdet = np.linalg.slogdet(C)
        return -0.5 * float(logdet) + self.expected_V(m, C)

    def star(self):
        m_star = np.zeros(2, dtype=np.float64)
        C_star = np.linalg.inv(self.Q)
        return m_star, symmetrize(C_star), self.objective(m_star, C_star)

    def energy_gap(self, m, C):
        """Analytic ``E(a) - E(a_star) = 1/2[m^T Q m + Tr(QC) - log det(QC) - d]``."""
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        QC = self.Q @ C
        sign, logdet = np.linalg.slogdet(QC)
        return 0.5 * (float(m @ self.Q @ m) + float(np.trace(QC)) - float(logdet) - self.d)

    def metadata(self):
        return {
            "target_name": self.name, "lambda": self.lam, "d": self.d,
            "Q_diag": [1.0, self.lam], "m0": self.m0.tolist(), "C0": self.C0.tolist(),
            "has_theory": True, "alpha": self.alpha, "beta": self.beta,
            "description": "Gaussian posterior N(0, Q^{-1}), Q=diag(1,lambda)",
        }


# ---------------------------------------------------------------------------
# Target B: non-separable log-concave posterior with a quartic tail
# ---------------------------------------------------------------------------

class LiteratureLogconcaveTarget:
    """Non-separable log-concave posterior ``V(x,y) = (sqrt(lam) x - y)^2/20 + y^4/20``.

    Convex but *not* globally smooth (the quartic term makes ``grad^2 V`` grow
    without bound), so no theorem-based global ``(alpha, beta)`` is reported
    (``has_theory == False``). All required expectations under ``N(m, C)`` close
    in terms of the Gaussian moments of ``Y``:

        E[Y^2] = my^2 + C22,  E[Y^3] = my^3 + 3 my C22,
        E[Y^4] = my^4 + 6 my^2 C22 + 3 C22^2.
    """

    name = "literature_logconcave"
    d = 2
    has_theory = False

    def __init__(self, lam):
        self.lam = float(lam)
        self.sl = float(np.sqrt(self.lam))
        self.m0 = np.array([10.0, 10.0], dtype=np.float64)
        self.C0 = 4.0 * np.eye(2, dtype=np.float64)
        self.alpha = None
        self.beta = None

    def g_H(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        mx, my = float(m[0]), float(m[1])
        C22 = float(C[1, 1])
        lam, sl = self.lam, self.sl
        # E[grad V]
        E_y3 = my ** 3 + 3.0 * my * C22
        gradx = (lam * mx - sl * my) / 10.0
        grady = (-sl * mx + my) / 10.0 + E_y3 / 5.0
        g = -np.array([gradx, grady], dtype=np.float64)
        # E[grad^2 V]
        E_y2 = my ** 2 + C22
        H = -np.array([
            [lam / 10.0,        -sl / 10.0],
            [-sl / 10.0, 1.0 / 10.0 + 0.6 * E_y2],
        ], dtype=np.float64)
        return g, symmetrize(H)

    def expected_V(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        mx, my = float(m[0]), float(m[1])
        C11, C12, C22 = float(C[0, 0]), float(C[0, 1]), float(C[1, 1])
        lam, sl = self.lam, self.sl
        # E[(sqrt(lam) X - Y)^2] = mean^2 + var of U = sqrt(lam) X - Y
        mean_U = sl * mx - my
        var_U = lam * C11 - 2.0 * sl * C12 + C22
        E_U2 = mean_U ** 2 + var_U
        E_y4 = my ** 4 + 6.0 * my ** 2 * C22 + 3.0 * C22 ** 2
        return E_U2 / 20.0 + E_y4 / 20.0

    def objective(self, m, C):
        C = symmetrize(C)
        sign, logdet = np.linalg.slogdet(C)
        return -0.5 * float(logdet) + self.expected_V(m, C)

    def metadata(self):
        return {
            "target_name": self.name, "lambda": self.lam, "d": self.d,
            "m0": self.m0.tolist(), "C0": self.C0.tolist(),
            "has_theory": False, "alpha": None, "beta": None,
            "description": "non-separable quartic V=(sqrt(lam)x-y)^2/20 + y^4/20",
        }


# ---------------------------------------------------------------------------
# Target C: smooth strongly log-concave non-separable posterior
# ---------------------------------------------------------------------------

class SmoothLogconcaveTarget:
    """Smooth strongly log-concave non-separable posterior.

        V(x,y) = (sqrt(lam) x - y)^2/20 + (delta/2)(x^2 + y^2) + gamma log cosh(y).

    Globally strongly convex and smooth, with constants

        alpha = delta,    beta = delta + (1 + lam)/10 + gamma,

    derived from the rank-one coupled quadratic (top eigenvalue ``(1+lam)/10``),
    the ``delta I`` ridge, and ``0 <= gamma sech^2 <= gamma``. Expectations
    involving ``tanh(Y)``, ``sech^2(Y)`` and ``log cosh(Y)`` are evaluated by
    deterministic Gauss--Hermite quadrature over ``Y ~ N(my, C22)``.
    """

    name = "smooth_logconcave"
    d = 2
    has_theory = True

    def __init__(self, lam, delta=0.05, gamma=1.0, n_nodes=80):
        self.lam = float(lam)
        self.sl = float(np.sqrt(self.lam))
        self.delta = float(delta)
        self.gamma = float(gamma)
        self.n_nodes = int(n_nodes)
        self.m0 = np.array([10.0, 10.0], dtype=np.float64)
        self.C0 = 4.0 * np.eye(2, dtype=np.float64)
        self.alpha = self.delta
        self.beta = self.delta + (1.0 + self.lam) / 10.0 + self.gamma

    def g_H(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        mx, my = float(m[0]), float(m[1])
        C22 = float(C[1, 1])
        lam, sl, delta, gamma = self.lam, self.sl, self.delta, self.gamma
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
        lam, sl, delta, gamma = self.lam, self.sl, self.delta, self.gamma
        mean_U = sl * mx - my
        var_U = lam * C11 - 2.0 * sl * C12 + C22
        E_U2 = mean_U ** 2 + var_U
        E_quad = mx ** 2 + C11 + my ** 2 + C22
        E_logcosh = gaussian_expectation_1d(
            lambda y: np.log(np.cosh(y)), my, C22, self.n_nodes)
        return E_U2 / 20.0 + 0.5 * delta * E_quad + gamma * E_logcosh

    def objective(self, m, C):
        C = symmetrize(C)
        sign, logdet = np.linalg.slogdet(C)
        return -0.5 * float(logdet) + self.expected_V(m, C)

    def metadata(self):
        return {
            "target_name": self.name, "lambda": self.lam, "d": self.d,
            "delta": self.delta, "gamma": self.gamma, "gh_nodes": self.n_nodes,
            "m0": self.m0.tolist(), "C0": self.C0.tolist(),
            "has_theory": True, "alpha": self.alpha, "beta": self.beta,
            "description": ("smooth strongly log-concave "
                            "V=(sqrt(lam)x-y)^2/20 + delta/2 |theta|^2 + gamma logcosh(y)"),
        }


# ---------------------------------------------------------------------------
# Target D: scalar Gaussian diagnostic
# ---------------------------------------------------------------------------

class ScalarGaussianTarget:
    """Scalar diagnostic ``rho_post = N(0, 1)`` (``d = 1``, ``V = theta^2/2``).

    Here ``g = -m`` and ``H = -1`` so the covariance recursions reduce to closed
    forms used directly by the diagnostic plots:

        KL:        C_{n+1} = (1 + dt) C_n / (1 + dt C_n)
        Riemann:   C_{n+1} = C_n exp(dt (1 - C_n))
        mean:      m_{n+1} = (1 - dt C_n) m_n   (shared by both schemes).

    The optimum is ``m_star = 0``, ``C_star = 1``, ``F_star = -1/2`` (objective
    ``F = -1/2 log C + (m^2 + C)/2``).
    """

    name = "scalar_gaussian"
    d = 1
    has_theory = True

    def __init__(self, m0=0.0, C0=1.0):
        self.m0 = np.array([float(m0)], dtype=np.float64)
        self.C0 = np.array([[float(C0)]], dtype=np.float64)
        self.alpha = 1.0
        self.beta = 1.0

    def g_H(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        return -m.copy(), -np.ones((1, 1), dtype=np.float64)

    def expected_V(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        return 0.5 * (float(m[0]) ** 2 + float(C[0, 0]))

    def objective(self, m, C):
        C = symmetrize(C)
        return -0.5 * float(np.log(C[0, 0])) + self.expected_V(m, C)

    def star(self):
        m_star = np.zeros(1, dtype=np.float64)
        C_star = np.eye(1, dtype=np.float64)
        return m_star, C_star, self.objective(m_star, C_star)

    @staticmethod
    def kl_cov_closed_form(C, dt):
        """Closed-form scalar KL covariance update ``(1 + dt) C / (1 + dt C)``."""
        return (1.0 + dt) * C / (1.0 + dt * C)

    @staticmethod
    def riemannian_cov_closed_form(C, dt):
        """Closed-form scalar Riemannian covariance update ``C exp(dt (1 - C))``."""
        return C * np.exp(dt * (1.0 - C))

    def metadata(self):
        return {
            "target_name": self.name, "d": self.d,
            "m0": self.m0.tolist(), "C0": self.C0.tolist(),
            "has_theory": True, "alpha": self.alpha, "beta": self.beta,
            "description": "scalar Gaussian N(0,1) diagnostic",
        }


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

_TARGET_CLASSES = {
    "gaussian_posterior": GaussianPosteriorTarget,
    "literature_logconcave": LiteratureLogconcaveTarget,
    "smooth_logconcave": SmoothLogconcaveTarget,
}


def build_target(name, lam, **kwargs):
    """Construct a two-dimensional target by ``name`` and ``lambda``."""
    if name not in _TARGET_CLASSES:
        raise ValueError(f"unknown target '{name}' (known: {sorted(_TARGET_CLASSES)})")
    return _TARGET_CLASSES[name](lam, **kwargs)
