"""Potential interface and the centered-feature construction.

A *potential* defines the target ``rho_post(theta) ~ exp(-V(theta))`` through its
value/gradient/Hessian of ``V``. Every potential exposes both single-point and
batched evaluators and a ``metadata()`` dict.

The natural-gradient local-rate experiment works in equilibrium-whitened
coordinates where the equilibrium is ``a_star = (0, I)``. The
:class:`CenteredPotential` wrapper turns an arbitrary raw feature ``Phi`` into a
potential ``V`` whose stationary conditions hold at the standard Gaussian:

    V(theta)      = 0.5||theta||^2 + rho * (Phi(theta) - b^T theta - 0.5 theta^T M theta)
    grad V(theta) = theta          + rho * (grad Phi(theta) - b - M theta)
    Hess V(theta) = I              + rho * (Hess Phi(theta) - M)

with ``b = E[grad Phi(Z)]`` and ``M = E[Hess Phi(Z)]`` estimated on a fixed
sample bank ``Z_ref ~ N(0, I)``. By construction

    E_{Z~N(0,I)}[grad V(Z)] ~= 0,    E_{Z~N(0,I)}[Hess V(Z)] ~= I,

(exactly on ``Z_ref`` when antithetic sampling makes ``mean(Z_ref) = 0``).
"""
from __future__ import annotations

import abc

import numpy as np

from src.common.monte_carlo import gaussian_samples


# ---------------------------------------------------------------------------
# Scalar nonlinearities (bounded curvature)
# ---------------------------------------------------------------------------

def log_cosh(s):
    """Numerically stable ``log cosh(s)``."""
    a = np.abs(s)
    return a + np.log1p(np.exp(-2.0 * a)) - np.log(2.0)


def log_cosh_prime(s):
    """d/ds log cosh(s) = tanh(s)."""
    return np.tanh(s)


def log_cosh_double(s):
    """d^2/ds^2 log cosh(s) = sech^2(s) = 1 - tanh^2(s), in (0, 1]."""
    t = np.tanh(s)
    return 1.0 - t * t


NONLINEARITIES = {
    "log_cosh": (log_cosh, log_cosh_prime, log_cosh_double),
}

# Tight bound on |phi''| for each nonlinearity (used for deterministic L_A).
NONLINEARITY_CURVATURE_BOUND = {
    "log_cosh": 1.0,
}


def get_nonlinearity(name):
    if name not in NONLINEARITIES:
        raise ValueError(f"unknown nonlinearity '{name}', have {list(NONLINEARITIES)}")
    return NONLINEARITIES[name]


# ---------------------------------------------------------------------------
# Base potential
# ---------------------------------------------------------------------------

class BasePotential(abc.ABC):
    """Abstract potential ``V`` for ``rho_post ~ exp(-V)``.

    Subclasses must implement the batched evaluators ``batch_value``,
    ``batch_grad``, ``batch_hess`` and ``metadata``. Single-point ``value``,
    ``grad`` and ``hess`` are derived from them.
    """

    N_theta: int
    kappa_target: float
    seed: int

    # -- single point (derived from batch) --
    def value(self, theta):
        return float(self.batch_value(np.atleast_2d(theta))[0])

    def grad(self, theta):
        return self.batch_grad(np.atleast_2d(theta))[0]

    def hess(self, theta):
        return self.batch_hess(np.atleast_2d(theta))[0]

    # -- batched (abstract) --
    @abc.abstractmethod
    def batch_value(self, theta_batch):
        """Return ``V(theta_j)`` for a batch, shape ``(M,)``."""

    @abc.abstractmethod
    def batch_grad(self, theta_batch):
        """Return ``grad V(theta_j)`` for a batch, shape ``(M, N_theta)``."""

    @abc.abstractmethod
    def batch_hess(self, theta_batch):
        """Return ``Hess V(theta_j)`` for a batch, shape ``(M, N_theta, N_theta)``."""

    @abc.abstractmethod
    def metadata(self):
        """Return a JSON-serializable dict describing the potential."""


# ---------------------------------------------------------------------------
# Gaussian potential V = 0.5 ||theta||^2  (target N(0, I), kappa = 1)
# ---------------------------------------------------------------------------

class GaussianPotential(BasePotential):
    """The exact Gaussian case ``V(theta) = 0.5 ||theta||^2``.

    Then ``grad V = theta``, ``Hess V = I``, the equilibrium ``a_star = (0, I)``
    is exact, and all the local operators vanish (``Lambda_hat = 0``,
    ``gamma_loc = 1``). Used as the analytic ground-truth in tests.
    """

    def __init__(self, N_theta, seed=0):
        self.N_theta = int(N_theta)
        self.kappa_target = 1.0
        self.seed = int(seed)

    def batch_value(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        return 0.5 * np.sum(T * T, axis=1)

    def batch_grad(self, theta_batch):
        return np.asarray(theta_batch, dtype=np.float64).copy()

    def batch_hess(self, theta_batch):
        M = np.asarray(theta_batch).shape[0]
        return np.broadcast_to(np.eye(self.N_theta), (M, self.N_theta, self.N_theta)).copy()

    def metadata(self):
        return {
            "family": "gaussian",
            "N_theta": self.N_theta,
            "kappa_target": 1.0,
            "alpha_target": 1.0,
            "beta_target": 1.0,
            "rho": 0.0,
            "L_A": 0.0,
            "L_A_is_empirical": False,
            "seed": self.seed,
        }


# ---------------------------------------------------------------------------
# Raw feature interface
# ---------------------------------------------------------------------------

class RawFeaturePotential(abc.ABC):
    """A raw feature map ``Phi`` (batched value / grad / Hessian).

    Subclasses provide ``phi_value``, ``phi_grad``, ``phi_hess`` for ``Phi`` and
    ``raw_metadata``. ``deterministic_LA`` may return a rigorous bound on
    ``sup_theta || Hess Phi(theta) - M ||_op`` (``M = E[Hess Phi(Z)]``) or
    ``None`` if no clean bound is available (then the centering falls back to an
    empirical estimate).
    """

    N_theta: int
    seed: int

    @abc.abstractmethod
    def phi_value(self, theta_batch):
        """``Phi(theta_j)``, shape ``(M,)``."""

    @abc.abstractmethod
    def phi_grad(self, theta_batch):
        """``grad Phi(theta_j)``, shape ``(M, N_theta)``."""

    @abc.abstractmethod
    def phi_hess(self, theta_batch):
        """``Hess Phi(theta_j)``, shape ``(M, N_theta, N_theta)``."""

    def deterministic_LA(self, M_center):
        """Optional rigorous bound on ``sup_theta ||Hess Phi - M||_op``."""
        return None

    @abc.abstractmethod
    def raw_metadata(self):
        """Family-specific metadata dict."""


# ---------------------------------------------------------------------------
# Centered potential
# ---------------------------------------------------------------------------

def _mean_hess_phi(raw, Z, chunk):
    """E[Hess Phi] over the bank ``Z``, accumulated in chunks (no full storage)."""
    N = Z.shape[1]
    acc = np.zeros((N, N), dtype=np.float64)
    for start in range(0, Z.shape[0], chunk):
        block = raw.phi_hess(Z[start:start + chunk])  # (c, N, N)
        acc += block.sum(axis=0)
    M = acc / Z.shape[0]
    return 0.5 * (M + M.T)


def _hess_dev_eig_extremes(raw, Z, M, chunk):
    """Global min / max eigenvalue of ``Hess Phi(Z_j) - M`` over the bank.

    Returns ``(emin, emax)``. This single pass yields both the empirical
    operator-norm deviation ``max(|emin|, |emax|)`` (for an empirical ``L_A``)
    and, after ``rho`` is known, the Hessian-eigenvalue extremes of ``Hess V``.
    """
    emin, emax = np.inf, -np.inf
    for start in range(0, Z.shape[0], chunk):
        block = raw.phi_hess(Z[start:start + chunk]) - M  # (c, N, N)
        w = np.linalg.eigvalsh(0.5 * (block + np.transpose(block, (0, 2, 1))))
        emin = min(emin, float(w[:, 0].min()))
        emax = max(emax, float(w[:, -1].max()))
    return emin, emax


class CenteredPotential(BasePotential):
    """Wrap a :class:`RawFeaturePotential` into a centered, conditioned ``V``.

    The coupling strength ``rho`` is chosen from ``kappa_target`` via the
    conservative bound: if ``Hess V = I + rho * A(theta)`` with
    ``||A(theta)||_op <= L_A`` then

        rho          = ((kappa - 1) / (kappa + 1)) / L_A
        alpha_target = 1 - rho * L_A = 2 / (kappa + 1)
        beta_target  = 1 + rho * L_A = 2 kappa / (kappa + 1)

    so the *nominal* condition number ``beta_target / alpha_target == kappa``.
    """

    def __init__(self, raw, kappa_target, *, centering_samples=8192,
                 centering_seed=None, safety_factor=2.0, la_mode="auto",
                 centering_chunk=2048, Z_ref=None, precomputed_stats=None):
        self.raw = raw
        self.N_theta = int(raw.N_theta)
        self.kappa_target = float(kappa_target)
        self.seed = int(getattr(raw, "seed", 0))
        self.safety_factor = float(safety_factor)
        self.la_mode = la_mode
        self.centering_chunk = int(centering_chunk)

        if precomputed_stats is not None:
            # Fast path: the three heavy bank reductions (b = E[grad Phi],
            # M = E[Hess Phi], and the (emin, emax) eigenvalue extremes of
            # Hess Phi - M over the bank) were already computed elsewhere (e.g. on
            # the GPU, see torch_backend.centering_stats_on_device). The reductions
            # are the *only* thing that touches the full bank; every line below is
            # the identical cheap centering algebra, so the resulting potential is
            # byte-for-byte the CPU potential it stands in for. ``Z_ref`` is not
            # stored (it is never read after construction) to avoid holding the 4M
            # bank twice.
            self.centering_samples = int(precomputed_stats["centering_samples"])
            self.centering_seed = precomputed_stats.get("centering_seed")
            self._Z_ref = None
            self.b = np.asarray(precomputed_stats["b"], dtype=np.float64).reshape(-1)
            self.M = np.asarray(precomputed_stats["M"], dtype=np.float64)
            self.M = 0.5 * (self.M + self.M.T)
            emin_dev = float(precomputed_stats["emin_dev"])
            emax_dev = float(precomputed_stats["emax_dev"])
            mean_Z = np.asarray(precomputed_stats["mean_Z"], dtype=np.float64).reshape(-1)
        else:
            if Z_ref is not None:
                # Share the operator/flow bank so that (0, I) is the *exact* discrete
                # fixed point (E[grad V] = 0, E[Hess V] = I hold on this bank).
                Z_ref = np.ascontiguousarray(Z_ref, dtype=np.float64)
                self.centering_samples = int(Z_ref.shape[0])
                self.centering_seed = None
            else:
                if centering_seed is None:
                    centering_seed = (self.seed + 1) * 100003 + 7
                self.centering_seed = int(centering_seed)
                self.centering_samples = int(centering_samples)
                Z_ref = gaussian_samples(self.N_theta, self.centering_samples,
                                         seed=self.centering_seed, antithetic=True)
            self._Z_ref = Z_ref

            # b = E[grad Phi] (cheap, full); M = E[Hess Phi] (chunked accumulation).
            grad_phi = raw.phi_grad(Z_ref)            # (M, N)
            self.b = np.mean(grad_phi, axis=0)
            self.M = _mean_hess_phi(raw, Z_ref, self.centering_chunk)

            # Single chunked pass: eigenvalue extremes of (Hess Phi - M).
            emin_dev, emax_dev = _hess_dev_eig_extremes(raw, Z_ref, self.M, self.centering_chunk)
            mean_Z = np.mean(Z_ref, axis=0)

        # --- choose L_A ---
        det_LA = None
        if la_mode in ("auto", "deterministic"):
            det_LA = raw.deterministic_LA(self.M)
        if det_LA is not None and la_mode != "empirical":
            self.L_A = float(det_LA)
            self.L_A_is_empirical = False
        else:
            emp = max(abs(emin_dev), abs(emax_dev))  # empirical ||Hess Phi - M||_op
            self.L_A = float(self.safety_factor * emp)
            self.L_A_is_empirical = True
            if la_mode == "deterministic":
                import warnings
                warnings.warn(
                    f"{type(raw).__name__}: no deterministic L_A available; "
                    "using empirical estimate with safety factor."
                )

        if self.L_A <= 0.0:
            raise ValueError("L_A must be positive; the raw feature appears to be constant.")

        kappa = self.kappa_target
        self.rho = ((kappa - 1.0) / (kappa + 1.0)) / self.L_A
        self.alpha_target = 1.0 - self.rho * self.L_A
        self.beta_target = 1.0 + self.rho * self.L_A

        # --- centering diagnostics ---
        # mean grad Phi == b exactly, so the residual reduces to the (anti-thetic)
        # sample-mean term; mean Hess V == I exactly because M == mean Hess Phi.
        # ``mean_Z`` was set above (from the bank, or supplied in precomputed_stats).
        mean_gradV = mean_Z - self.rho * (mean_Z @ self.M)
        self.norm_mean_grad = float(np.linalg.norm(mean_gradV))
        self.norm_mean_hess_minus_I = 0.0  # exact by construction

        # Hess V eigenvalues = 1 + rho * eig(Hess Phi - M); use the extremes.
        self.empirical_min_hess_eig = float(1.0 + self.rho * emin_dev)
        self.empirical_max_hess_eig = float(1.0 + self.rho * emax_dev)

    # -- batched evaluators --
    def batch_value(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        quad = 0.5 * np.sum(T * T, axis=1)
        phi = self.raw.phi_value(T)
        lin = T @ self.b
        gauss = 0.5 * np.einsum("mi,ij,mj->m", T, self.M, T, optimize=True)
        return quad + self.rho * (phi - lin - gauss)

    def batch_grad(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        gphi = self.raw.phi_grad(T)
        return T + self.rho * (gphi - self.b[None, :] - T @ self.M)

    def batch_hess(self, theta_batch):
        T = np.asarray(theta_batch, dtype=np.float64)
        Mb = T.shape[0]
        hphi = self.raw.phi_hess(T)
        I = np.eye(self.N_theta)
        H = I[None, :, :] + self.rho * (hphi - self.M[None, :, :])
        return 0.5 * (H + np.transpose(H, (0, 2, 1)))

    def metadata(self):
        md = {
            "family": type(self.raw).__name__,
            "N_theta": self.N_theta,
            "seed": self.seed,
            "kappa_target": self.kappa_target,
            "alpha_target": self.alpha_target,
            "beta_target": self.beta_target,
            "rho": self.rho,
            "L_A": self.L_A,
            "estimated_L_A": self.L_A,
            "L_A_is_empirical": self.L_A_is_empirical,
            "safety_factor": self.safety_factor,
            "centering_samples": self.centering_samples,
            "centering_seed": self.centering_seed,
            "norm_mean_grad": self.norm_mean_grad,
            "norm_mean_hess_minus_I": self.norm_mean_hess_minus_I,
            "empirical_min_hess_eig": self.empirical_min_hess_eig,
            "empirical_max_hess_eig": self.empirical_max_hess_eig,
        }
        md.update(self.raw.raw_metadata())
        return md
