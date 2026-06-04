"""
Initial conditions for log-concave gradient flow experiments.

All five initializations are defined **relative to the reference Gaussian
optimum a_star = (m_star, C_star)**, not relative to N(0, I).

This is essential because the true posterior is no longer isotropic.
Initializations in the a_star coordinate frame ensure that the five
conditions probe the same geometric regimes as in the Gaussian case.

Let  C_sqrt = C_star^{1/2}  (symmetric positive-definite square root).

Initializations
---------------
mean_only   : m0 = m_star + r * C_sqrt @ 1/sqrt(n),  C0 = C_star,  r = 3
              Tests convergence from a large mean offset with correct covariance.

volume_high : m0 = m_star,  C0 = 4 * C_star
              Tests convergence when covariance is inflated (volume too high).

volume_low  : m0 = m_star,  C0 = 0.25 * C_star
              Tests convergence when covariance is deflated (volume too low).

shape_only  : m0 = m_star,
              C0 = C_sqrt @ diag(e^r, e^{-r}, 1,...,1) @ C_sqrt,  r = 2
              Volume (det) matches C_star; only shape (anisotropy) differs.

mixed       : m0 = m_star + 2 * C_sqrt @ 1/sqrt(n),
              C0 = C_sqrt @ [2 * diag(e^r, e^{-r}, 1,...,1)] @ C_sqrt,  r=1.5
              Nonzero mean offset + scaled anisotropic covariance.
"""
import numpy as np
from src.omega_tau_modes.utils import spd_sqrt


LC_INIT_NAMES = ["mean_only", "volume_high", "volume_low", "shape_only", "mixed"]


def _aniso_diag(n, r):
    """Return diagonal (exp(r), exp(-r), 1, ..., 1) of length n."""
    d = np.ones(n, dtype=np.float64)
    d[0] = np.exp(r)
    if n >= 2:
        d[1] = np.exp(-r)
    return d


def get_logconcave_initialization(name: str, n: int,
                                   m_star: np.ndarray,
                                   C_star: np.ndarray):
    """Return (m0, C0) for the named initialization relative to a_star.

    Args:
        name   : one of LC_INIT_NAMES
        n      : dimension
        m_star : reference mean, shape (n,)
        C_star : reference covariance, shape (n, n), SPD

    Returns:
        (m0, C0) : copies of initial mean and covariance, both float64 SPD

    Raises:
        ValueError : if name is not recognised
    """
    C_sqrt = spd_sqrt(C_star)      # C_star^{1/2}

    ones_dir = np.ones(n, dtype=np.float64) / np.sqrt(n)

    if name == "mean_only":
        r = 3.0
        m0 = m_star + r * (C_sqrt @ ones_dir)
        C0 = C_star.copy()

    elif name == "volume_high":
        m0 = m_star.copy()
        C0 = 4.0 * C_star

    elif name == "volume_low":
        m0 = m_star.copy()
        C0 = 0.25 * C_star

    elif name == "shape_only":
        # det(D) = exp(r)*exp(-r)*1*...*1 = 1, so det(C0) = det(C_star)
        r = 2.0
        m0 = m_star.copy()
        D = np.diag(_aniso_diag(n, r))
        C0 = C_sqrt @ D @ C_sqrt

    elif name == "mixed":
        r = 1.5
        s = 2.0
        m0 = m_star + 2.0 * (C_sqrt @ ones_dir)
        D = s * np.diag(_aniso_diag(n, r))
        C0 = C_sqrt @ D @ C_sqrt

    else:
        raise ValueError(
            f"Unknown log-concave initialization '{name}'. "
            f"Choose from: {LC_INIT_NAMES}."
        )

    # Symmetrise for numerical cleanliness
    C0 = 0.5 * (C0 + C0.T)
    return m0.copy(), C0.copy()
