"""Targets for the STL variance experiment.

Both targets have a *diagonal* curvature, so the Gaussian VI optimum is diagonal
and the whole stochastic flow stays diagonal when started diagonal (the code
nevertheless uses general matrix linear algebra so it is correct for any state).

Each target exposes, for the current variational state ``a = (m, C)``:

* model parameters ``a_diag`` (the diagonal of ``A``), ``tau`` and ``kind``;
* :meth:`score(theta)`     -> ``score_post(theta) = grad log rho_post = -grad V``;
* :meth:`hess_diag(theta)` -> the diagonal of ``Hess_log_post(theta) = -grad^2 V``
  (both targets have a diagonal Hessian);
* :meth:`exp_score(m, C)`  -> the exact / quadrature reference ``E_q[score_post]``;
* :meth:`a_star()`         -> the Gaussian VI optimum ``(m_star, C_star)``;
* :meth:`objective(m, C)`  -> ``E_{N(m,C)}[V] - 1/2 log det C`` (KL up to a const.);
* :meth:`energy_gap(m, C)` -> ``E(a) - E(a_star)`` (analytic for the Gaussian
  target, Gauss--Hermite for ``log cosh``);
* :meth:`metadata()`.

**Target 1 -- well-specified anisotropic Gaussian.**
``V(theta) = 1/2 theta^T A theta`` with ``A = diag(a)``, ``a`` log-spaced in
``[1, kappa]``. Then ``score = -A theta``, ``Hess_log_post = -A`` and the optimum
is exact: ``m_star = 0``, ``C_star = A^{-1}``. At the optimum the STL estimator is
pointwise zero.

**Target 2 -- misspecified smooth log-cosh.**
``V(theta) = 1/2 theta^T A theta + tau sum_i log cosh(theta_i)`` (separable,
strongly log-concave, globally smooth, non-Gaussian). ``score = -A theta - tau
tanh(theta)``, ``Hess_log_post = -A - tau diag(sech^2(theta_i))``. The optimum has
``m_star = 0`` by symmetry; the diagonal ``C_star`` is found by a deterministic
fixed-point solve of the per-coordinate stationarity condition
``c_i = 1/(a_i + tau E_{N(0,c_i)}[sech^2])`` (Gauss--Hermite quadrature).
"""
from __future__ import annotations

import numpy as np

from src.common.spd import symmetrize
# Single source of truth for the 1-D Gaussian Gauss--Hermite quadrature.
from src.natural_gradient_discretization_stepsize.targets import (
    gauss_hermite_nodes, _sech2,
)

TARGET_NAMES = ["gaussian", "log_cosh"]


def anisotropic_diagonal(d, kappa):
    """Diagonal curvature ``a`` of length ``d`` log-spaced between ``1`` and ``kappa``.

    For ``d == 1`` returns ``[1.0]``. The condition number of ``A = diag(a)`` is
    exactly ``kappa`` for ``d >= 2``.
    """
    d = int(d)
    kappa = float(kappa)
    if d <= 1:
        return np.ones(1, dtype=np.float64)
    return np.logspace(0.0, np.log10(kappa), d).astype(np.float64)


def _e_logcosh_diag(mean, var, nodes, weights):
    """``E[log cosh(Y)]`` for ``Y_i ~ N(mean_i, var_i)`` over a vector of coords."""
    mean = np.asarray(mean, dtype=np.float64).reshape(-1, 1)
    sd = np.sqrt(np.maximum(np.asarray(var, dtype=np.float64), 0.0)).reshape(-1, 1)
    y = mean + sd * nodes.reshape(1, -1)          # (d, n_nodes)
    # log cosh(y) computed stably as |y| + log1p(exp(-2|y|)) - log 2.
    ay = np.abs(y)
    lc = ay + np.log1p(np.exp(-2.0 * ay)) - np.log(2.0)
    return lc @ weights                            # (d,)


def _e_sech2_diag(mean, var, nodes, weights):
    """``E[sech^2(Y)]`` for ``Y_i ~ N(mean_i, var_i)`` over a vector of coords."""
    mean = np.asarray(mean, dtype=np.float64).reshape(-1, 1)
    sd = np.sqrt(np.maximum(np.asarray(var, dtype=np.float64), 0.0)).reshape(-1, 1)
    y = mean + sd * nodes.reshape(1, -1)
    return _sech2(y) @ weights


def _e_tanh_diag(mean, var, nodes, weights):
    """``E[tanh(Y)]`` for ``Y_i ~ N(mean_i, var_i)`` over a vector of coords."""
    mean = np.asarray(mean, dtype=np.float64).reshape(-1, 1)
    sd = np.sqrt(np.maximum(np.asarray(var, dtype=np.float64), 0.0)).reshape(-1, 1)
    y = mean + sd * nodes.reshape(1, -1)
    return np.tanh(y) @ weights


# ---------------------------------------------------------------------------
# Target 1: well-specified anisotropic Gaussian
# ---------------------------------------------------------------------------

class GaussianTarget:
    """Exact anisotropic Gaussian target ``rho_post = N(0, A^{-1})``, ``A=diag(a)``."""

    name = "gaussian"
    kind = "gaussian"
    has_theory = True

    def __init__(self, d, kappa):
        self.d = int(d)
        self.kappa = float(kappa)
        self.tau = 0.0
        self.a_diag = anisotropic_diagonal(self.d, self.kappa)
        self.A = np.diag(self.a_diag)
        self.alpha = float(self.a_diag.min())
        self.beta = float(self.a_diag.max())

    # -- pointwise sampling-based fields (batched over leading axes) ----------
    def score(self, theta):
        """``score_post(theta) = -A theta`` (broadcast over leading axes)."""
        return -(np.asarray(theta, dtype=np.float64) * self.a_diag)

    def hess_diag(self, theta):
        """Diagonal of ``Hess_log_post = -A`` (state-independent)."""
        theta = np.asarray(theta, dtype=np.float64)
        return np.broadcast_to(-self.a_diag, theta.shape).copy()

    # -- exact reference expectations ----------------------------------------
    def exp_score(self, m, C):
        """``E_q[score_post] = -A m`` (exact)."""
        return -(self.a_diag * np.asarray(m, dtype=np.float64))

    def exp_hess_diag(self, m, C):
        """``E_q[Hess_log_post] = -A`` (exact, state-independent)."""
        return -self.a_diag.copy()

    def a_star(self):
        """Gaussian VI optimum ``(m_star, C_star) = (0, A^{-1})``."""
        m_star = np.zeros(self.d, dtype=np.float64)
        C_star = np.diag(1.0 / self.a_diag)
        return m_star, symmetrize(C_star)

    def objective(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        _, logdet = np.linalg.slogdet(C)
        eV = 0.5 * (float(m @ self.A @ m) + float(np.trace(self.A @ C)))
        return -0.5 * float(logdet) + eV

    def energy_gap(self, m, C):
        """Analytic ``E(a) - E(a_star) = KL(N(m,C) || N(0, A^{-1}))``."""
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        AC = self.A @ C
        _, logdet = np.linalg.slogdet(AC)
        return 0.5 * (float(m @ self.A @ m) + float(np.trace(AC))
                      - float(logdet) - self.d)

    def metadata(self):
        m_star, C_star = self.a_star()
        return {
            "target_name": self.name, "kind": self.kind,
            "d": self.d, "kappa": self.kappa, "tau": self.tau,
            "a_diag": self.a_diag.tolist(),
            "alpha": self.alpha, "beta": self.beta,
            "m_star": m_star.tolist(),
            "C_star_diag": np.diag(C_star).tolist(),
            "a_star_method": "closed_form (m=0, C=A^{-1})",
            "has_theory": True,
            "description": "well-specified anisotropic Gaussian, V=0.5 theta^T A theta",
        }


# ---------------------------------------------------------------------------
# Target 2: misspecified smooth log-cosh
# ---------------------------------------------------------------------------

class LogCoshTarget:
    """Separable smooth strongly log-concave ``log cosh`` target.

    ``V(theta) = 1/2 theta^T A theta + tau sum_i log cosh(theta_i)`` with
    ``A = diag(a)`` log-spaced in ``[1, kappa]``. The optimum mean is zero by
    symmetry; the diagonal optimal covariance is computed by a deterministic
    fixed-point solve (see :meth:`a_star`).
    """

    name = "log_cosh"
    kind = "log_cosh"
    has_theory = True

    def __init__(self, d, kappa, tau, n_nodes=80, fp_tol=1e-13, fp_max_iter=500):
        self.d = int(d)
        self.kappa = float(kappa)
        self.tau = float(tau)
        self.n_nodes = int(n_nodes)
        self.fp_tol = float(fp_tol)
        self.fp_max_iter = int(fp_max_iter)
        self.a_diag = anisotropic_diagonal(self.d, self.kappa)
        self.A = np.diag(self.a_diag)
        self._nodes, self._weights = gauss_hermite_nodes(self.n_nodes)
        # Global strong log-concavity / smoothness of grad^2 V = A + tau sech^2:
        # sech^2 in (0, 1], so alpha = min a_i, beta = max a_i + tau.
        self.alpha = float(self.a_diag.min())
        self.beta = float(self.a_diag.max() + self.tau)
        self._c_star_diag = None
        self._a_star_diag = None  # alias kept for clarity in metadata

    # -- pointwise sampling-based fields --------------------------------------
    def score(self, theta):
        """``score_post(theta) = -A theta - tau tanh(theta)``."""
        theta = np.asarray(theta, dtype=np.float64)
        return -(theta * self.a_diag) - self.tau * np.tanh(theta)

    def hess_diag(self, theta):
        """Diagonal of ``Hess_log_post = -A - tau diag(sech^2(theta_i))``."""
        theta = np.asarray(theta, dtype=np.float64)
        return -(self.a_diag + self.tau * _sech2(theta))

    # -- reference expectations (Gauss--Hermite) ------------------------------
    def exp_score(self, m, C):
        """``E_q[score_post]_i = -a_i m_i - tau E[tanh(theta_i)]`` (marginal)."""
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        var = np.diag(C)
        e_tanh = _e_tanh_diag(m, var, self._nodes, self._weights)
        return -(self.a_diag * m) - self.tau * e_tanh

    def exp_hess_diag(self, m, C):
        """``E_q[Hess_log_post]_ii = -a_i - tau E[sech^2(theta_i)]`` (marginal)."""
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        var = np.diag(C)
        e_sech2 = _e_sech2_diag(m, var, self._nodes, self._weights)
        return -(self.a_diag + self.tau * e_sech2)

    def a_star(self):
        """Gaussian VI optimum ``(m_star, C_star)`` (deterministic fixed point).

        ``m_star = 0`` by symmetry. The diagonal ``C_star = diag(c)`` solves, per
        coordinate, the natural-gradient stationarity condition
        ``c_i^{-1} = a_i + tau E_{N(0,c_i)}[sech^2]``, i.e. ``C_star^{-1} =
        -E_q[Hess_log_post]`` at ``m=0``. The fixed-point map is a contraction
        (``E[sech^2]`` decreases in ``c``) and converges from ``c = 1/a``.
        """
        if self._c_star_diag is None:
            c = 1.0 / self.a_diag                       # Gaussian-part start
            for _ in range(self.fp_max_iter):
                q = _e_sech2_diag(np.zeros(self.d), c, self._nodes, self._weights)
                c_new = 1.0 / (self.a_diag + self.tau * q)
                if np.max(np.abs(c_new - c)) < self.fp_tol:
                    c = c_new
                    break
                c = c_new
            self._c_star_diag = c
        m_star = np.zeros(self.d, dtype=np.float64)
        C_star = np.diag(self._c_star_diag)
        return m_star, symmetrize(C_star)

    def a_star_diagnostics(self):
        """Stationarity residual ``|| C_star^{-1} + E_q[Hess] ||`` at the optimum."""
        m_star, C_star = self.a_star()
        c = self._c_star_diag
        e_hess = self.exp_hess_diag(m_star, C_star)     # = -(a + tau E[sech2])
        residual = float(np.max(np.abs(1.0 / c + e_hess)))
        return {
            "method": "diagonal_fixed_point",
            "fp_tol": self.fp_tol, "fp_max_iter": self.fp_max_iter,
            "gh_nodes": self.n_nodes,
            "stationarity_residual_max": residual,
            "c_star_diag": c.tolist(),
        }

    def expected_V(self, m, C):
        m = np.asarray(m, dtype=np.float64)
        C = symmetrize(C)
        var = np.diag(C)
        quad = 0.5 * float(self.a_diag @ (m * m + var))
        e_logcosh = _e_logcosh_diag(m, var, self._nodes, self._weights)
        return quad + self.tau * float(np.sum(e_logcosh))

    def objective(self, m, C):
        C = symmetrize(C)
        _, logdet = np.linalg.slogdet(C)
        return -0.5 * float(logdet) + self.expected_V(m, C)

    def energy_gap(self, m, C):
        m_star, C_star = self.a_star()
        return self.objective(m, C) - self.objective(m_star, C_star)

    def metadata(self):
        m_star, C_star = self.a_star()
        diag = self.a_star_diagnostics()
        return {
            "target_name": self.name, "kind": self.kind,
            "d": self.d, "kappa": self.kappa, "tau": self.tau,
            "a_diag": self.a_diag.tolist(),
            "alpha": self.alpha, "beta": self.beta,
            "m_star": m_star.tolist(),
            "C_star_diag": np.diag(C_star).tolist(),
            "a_star_method": ("diagonal fixed point c_i = 1/(a_i + tau "
                              "E_{N(0,c_i)}[sech^2]); m_star=0 by symmetry"),
            "a_star_diagnostics": diag,
            "F_star": float(self.objective(m_star, C_star)),
            "gh_nodes": self.n_nodes,
            "has_theory": True,
            "description": ("misspecified smooth strongly log-concave "
                            "V=0.5 theta^T A theta + tau sum log cosh(theta_i)"),
        }


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

_TARGET_CLASSES = {"gaussian": GaussianTarget, "log_cosh": LogCoshTarget}


def build_target(name, d, kappa, tau=0.0, n_nodes=80):
    """Construct a target by ``name`` with dimension ``d``, conditioning ``kappa``.

    ``tau`` and ``n_nodes`` are used only by the ``log_cosh`` target.
    """
    if name == "gaussian":
        return GaussianTarget(d, kappa)
    if name == "log_cosh":
        return LogCoshTarget(d, kappa, tau, n_nodes=n_nodes)
    raise ValueError(f"unknown target '{name}' (known: {TARGET_NAMES})")
