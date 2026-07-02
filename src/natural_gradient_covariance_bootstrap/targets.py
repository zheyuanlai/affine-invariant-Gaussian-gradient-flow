"""Targets for the covariance-bootstrap experiment (curvature convention).

Each target exposes, for a diagonal Gaussian state ``a = (m, C)``:

* pointwise fields (batched over leading axes):
  :meth:`grad_V(X)`, :meth:`hess_V_diag(X)`, :meth:`score(X) = -grad_V(X)`;
* deterministic reference expectations
  :meth:`G(m, C) = E_q[grad V]` and :meth:`A_diag(m, C) = diag E_q[grad^2 V]`;
* :meth:`a_star()` -> ``(m_star, C_star)``, :meth:`objective(m, C)`,
  :meth:`energy_gap(m, C) = E(a) - E(a_star)``;
* :meth:`Psi(m, C) = E || C^{1/2}(grad^2 V(X) - A) C^{1/2} ||_F^2`` (the intrinsic
  stochastic Hessian-fluctuation intensity; exactly ``0`` for the Gaussian target);
* the strong-log-concavity / smoothness constants ``(alpha, beta)`` and metadata.

Everything is diagonal: the curvature ``H`` is diagonal and both targets keep a
diagonal state diagonal, so the exact expectations are one-dimensional
Gauss--Hermite quadratures on the marginals. The code nevertheless returns full
matrices for ``A`` and ``C_star`` so the callers use ordinary matrix algebra.

This is the ``A = E[grad^2 V]`` (curvature) convention; it equals the
``H_disc = E[grad^2 log rho_post] = -A`` convention used by the discretization,
WFR, and STL groups. All expectations are exact (closed form / Gauss--Hermite);
there is no Monte Carlo in this module.
"""
from __future__ import annotations

import numpy as np

from src.common.spd import symmetrize
# Single source of truth for the 1-D Gaussian Gauss--Hermite quadrature.
from src.natural_gradient_discretization_stepsize.targets import (
    gauss_hermite_nodes, _sech2,
)

TARGET_NAMES = ["gaussian", "log_cosh"]


def curvature_diagonal(d, alpha, beta):
    """Diagonal curvature ``H`` of length ``d`` log-spaced between ``alpha`` and ``beta``.

    For ``d == 1`` returns ``[beta]``. For ``d >= 2`` the condition number of
    ``H = diag(h)`` is exactly ``beta / alpha``.
    """
    d = int(d)
    alpha = float(alpha)
    beta = float(beta)
    if d <= 1:
        return np.array([beta], dtype=np.float64)
    return np.logspace(np.log10(alpha), np.log10(beta), d).astype(np.float64)


# ---------------------------------------------------------------------------
# Diagonal Gauss--Hermite helpers (E[f(Y_i)] for Y_i ~ N(mean_i, var_i))
# ---------------------------------------------------------------------------

def _e_diag(f, mean, var, nodes, weights):
    """``E[f(Y_i)]`` for a vector of independent ``Y_i ~ N(mean_i, var_i)``."""
    mean = np.asarray(mean, dtype=np.float64).reshape(-1, 1)
    sd = np.sqrt(np.maximum(np.asarray(var, dtype=np.float64), 0.0)).reshape(-1, 1)
    y = mean + sd * nodes.reshape(1, -1)          # (d, n_nodes)
    return f(y) @ weights                          # (d,)


def _log_cosh(y):
    """Numerically stable ``log cosh(y) = |y| + log1p(exp(-2|y|)) - log 2``."""
    ay = np.abs(y)
    return ay + np.log1p(np.exp(-2.0 * ay)) - np.log(2.0)


# ---------------------------------------------------------------------------
# Target A: exact diagonal Gaussian
# ---------------------------------------------------------------------------

class GaussianTarget:
    """Exact diagonal Gaussian target ``V(x) = 1/2 x^T H x``.

    ``grad V = H x``, ``grad^2 V = H`` are exact and state-independent, so
    ``G = H m``, ``A = H``, the optimum is ``(0, H^{-1})`` and the energy gap is
    the analytic Gaussian KL. ``Psi = 0`` (deterministic curvature).
    """

    name = "gaussian"
    kind = "gaussian"
    has_theory = True

    def __init__(self, d, kappa, alpha=1.0):
        self.d = int(d)
        self.kappa = float(kappa)
        self.gamma = 0.0
        self.H_diag = curvature_diagonal(self.d, alpha, alpha * self.kappa)
        self.H = np.diag(self.H_diag)
        self.alpha = float(self.H_diag.min())
        self.beta = float(self.H_diag.max())

    # -- pointwise fields -----------------------------------------------------
    def grad_V(self, X):
        return np.asarray(X, dtype=np.float64) * self.H_diag

    def hess_V_diag(self, X):
        X = np.asarray(X, dtype=np.float64)
        return np.broadcast_to(self.H_diag, X.shape).copy()

    def score(self, X):
        return -self.grad_V(X)

    # -- deterministic reference expectations ---------------------------------
    def G(self, m, C):
        """``E_q[grad V] = H m``."""
        return self.H_diag * np.asarray(m, dtype=np.float64)

    def A_diag(self, m, C):
        """Diagonal of ``E_q[grad^2 V] = H`` (state-independent)."""
        return self.H_diag.copy()

    def A_matrix(self, m, C):
        return np.diag(self.A_diag(m, C))

    def a_star(self):
        m_star = np.zeros(self.d, dtype=np.float64)
        C_star = np.diag(1.0 / self.H_diag)
        return m_star, symmetrize(C_star)

    def objective(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        _, logdet = np.linalg.slogdet(C)
        eV = 0.5 * (float(m @ self.H @ m) + float(np.trace(self.H @ C)))
        return -0.5 * float(logdet) + eV

    def energy_gap(self, m, C):
        """Analytic ``Delta(m,C) = 1/2 (m^T H m + tr(HC) - logdet(HC) - d)``."""
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        HC = self.H @ C
        _, logdet = np.linalg.slogdet(HC)
        return 0.5 * (float(m @ self.H @ m) + float(np.trace(HC))
                      - float(logdet) - self.d)

    def Psi(self, m, C):
        """Intrinsic stochastic Hessian fluctuation; zero for a Gaussian target."""
        return 0.0

    def metadata(self):
        m_star, C_star = self.a_star()
        return {
            "target_name": self.name, "kind": self.kind,
            "d": self.d, "kappa": self.kappa, "gamma": self.gamma,
            "H_diag": self.H_diag.tolist(),
            "alpha": self.alpha, "beta": self.beta,
            "m_star": m_star.tolist(), "C_star_diag": np.diag(C_star).tolist(),
            "a_star_method": "closed_form (m=0, C=H^{-1})",
            "has_theory": True, "Psi_star": 0.0,
            "description": "exact diagonal Gaussian, V=0.5 x^T H x",
        }


# ---------------------------------------------------------------------------
# Target B: smooth strongly log-concave log-cosh
# ---------------------------------------------------------------------------

class LogCoshTarget:
    """Separable smooth strongly log-concave ``V = 1/2 x^T H x + gamma sum log cosh(x_i)``.

    ``grad V = H x + gamma tanh(x)``, ``grad^2 V = H + gamma diag(sech^2 x_i)``, so
    ``alpha = min H`` and ``beta = max H + gamma``. The optimum has ``m_star = 0``
    (symmetry) and diagonal ``C_star`` from a per-coordinate fixed point. Reference
    expectations, the optimum, the energy gap and ``Psi`` use deterministic
    Gauss--Hermite quadrature on the diagonal marginals.
    """

    name = "log_cosh"
    kind = "log_cosh"
    has_theory = True

    def __init__(self, d, kappa, gamma, alpha=1.0, n_nodes=80,
                 fp_tol=1e-13, fp_max_iter=500):
        self.d = int(d)
        self.kappa = float(kappa)
        self.gamma = float(gamma)
        self.n_nodes = int(n_nodes)
        self.fp_tol = float(fp_tol)
        self.fp_max_iter = int(fp_max_iter)
        self.H_diag = curvature_diagonal(self.d, alpha, alpha * self.kappa)
        self.H = np.diag(self.H_diag)
        self._nodes, self._weights = gauss_hermite_nodes(self.n_nodes)
        # grad^2 V = H + gamma sech^2 with sech^2 in (0, 1].
        self.alpha = float(self.H_diag.min())
        self.beta = float(self.H_diag.max() + self.gamma)
        self._c_star_diag = None

    # -- pointwise fields -----------------------------------------------------
    def grad_V(self, X):
        X = np.asarray(X, dtype=np.float64)
        return X * self.H_diag + self.gamma * np.tanh(X)

    def hess_V_diag(self, X):
        X = np.asarray(X, dtype=np.float64)
        return self.H_diag + self.gamma * _sech2(X)

    def score(self, X):
        return -self.grad_V(X)

    # -- deterministic reference expectations ---------------------------------
    def G(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        var = np.diag(symmetrize(C))
        e_tanh = _e_diag(np.tanh, m, var, self._nodes, self._weights)
        return self.H_diag * m + self.gamma * e_tanh

    def A_diag(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        var = np.diag(symmetrize(C))
        e_sech2 = _e_diag(_sech2, m, var, self._nodes, self._weights)
        return self.H_diag + self.gamma * e_sech2

    def A_matrix(self, m, C):
        return np.diag(self.A_diag(m, C))

    def a_star(self):
        if self._c_star_diag is None:
            c = 1.0 / self.H_diag
            for _ in range(self.fp_max_iter):
                q = _e_diag(_sech2, np.zeros(self.d), c, self._nodes, self._weights)
                c_new = 1.0 / (self.H_diag + self.gamma * q)
                if np.max(np.abs(c_new - c)) < self.fp_tol:
                    c = c_new
                    break
                c = c_new
            self._c_star_diag = c
        m_star = np.zeros(self.d, dtype=np.float64)
        return m_star, symmetrize(np.diag(self._c_star_diag))

    def a_star_diagnostics(self):
        m_star, C_star = self.a_star()
        c = self._c_star_diag
        residual = float(np.max(np.abs(1.0 / c - self.A_diag(m_star, C_star))))
        return {"method": "diagonal_fixed_point", "fp_tol": self.fp_tol,
                "gh_nodes": self.n_nodes, "stationarity_residual_max": residual,
                "c_star_diag": c.tolist()}

    def expected_V(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        var = np.diag(symmetrize(C))
        quad = 0.5 * float(self.H_diag @ (m * m + var))
        e_lc = _e_diag(_log_cosh, m, var, self._nodes, self._weights)
        return quad + self.gamma * float(np.sum(e_lc))

    def objective(self, m, C):
        C = symmetrize(C)
        _, logdet = np.linalg.slogdet(C)
        return -0.5 * float(logdet) + self.expected_V(m, C)

    def energy_gap(self, m, C):
        m_star, C_star = self.a_star()
        return self.objective(m, C) - self.objective(m_star, C_star)

    def Psi(self, m, C):
        """``Psi = sum_i c_i^2 gamma^2 Var(sech^2 X_i)`` (diagonal state)."""
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        c = np.diag(C)
        var = c
        e1 = _e_diag(_sech2, m, var, self._nodes, self._weights)
        e2 = _e_diag(lambda y: _sech2(y) ** 2, m, var, self._nodes, self._weights)
        var_sech2 = np.maximum(e2 - e1 ** 2, 0.0)
        return float(np.sum(c ** 2 * self.gamma ** 2 * var_sech2))

    def metadata(self):
        m_star, C_star = self.a_star()
        diag = self.a_star_diagnostics()
        return {
            "target_name": self.name, "kind": self.kind,
            "d": self.d, "kappa": self.kappa, "gamma": self.gamma,
            "H_diag": self.H_diag.tolist(),
            "alpha": self.alpha, "beta": self.beta,
            "m_star": m_star.tolist(), "C_star_diag": np.diag(C_star).tolist(),
            "a_star_method": ("diagonal fixed point c_i = 1/(h_i + gamma "
                              "E_{N(0,c_i)}[sech^2]); m_star=0 by symmetry"),
            "a_star_diagnostics": diag,
            "Psi_star": self.Psi(m_star, C_star),
            "gh_nodes": self.n_nodes, "has_theory": True,
            "description": ("smooth strongly log-concave "
                            "V=0.5 x^T H x + gamma sum log cosh(x_i)"),
        }


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

def build_target(name, d, kappa, gamma=0.0, alpha=1.0, n_nodes=80):
    """Construct a target by ``name`` with dimension ``d`` and conditioning ``kappa``.

    ``gamma`` and ``n_nodes`` are used only by the ``log_cosh`` target.
    """
    if name == "gaussian":
        return GaussianTarget(d, kappa, alpha=alpha)
    if name == "log_cosh":
        return LogCoshTarget(d, kappa, gamma, alpha=alpha, n_nodes=n_nodes)
    raise ValueError(f"unknown target '{name}' (known: {TARGET_NAMES})")
