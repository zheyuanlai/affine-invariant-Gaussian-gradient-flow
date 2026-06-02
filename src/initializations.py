"""
Initial conditions (m0, C0) for affine-invariant Gaussian gradient flow experiments.

Five initializations designed to probe different convergence regimes:

  mean_only   — pure mean offset, identity covariance
                tests mean-dominated convergence
  volume_high — zero mean, inflated isotropic covariance
                tests volume shrinkage
  volume_low  — zero mean, deflated isotropic covariance
                tests volume expansion
  shape_only  — zero mean, anisotropic covariance with det = 1
                tests shape (eigenvector spread) convergence
  mixed       — nonzero mean + scaled anisotropic covariance
                tests simultaneous mean / volume / shape convergence
"""
import numpy as np

INIT_NAMES = ["mean_only", "volume_high", "volume_low", "shape_only", "mixed"]


def _anisotropic_diag(n, r):
    """Return diagonal (exp(r), exp(-r), 1, ..., 1) of length n."""
    d = np.ones(n, dtype=np.float64)
    d[0] = np.exp(r)
    if n >= 2:
        d[1] = np.exp(-r)
    return d


def get_initialization(name, n):
    """Return (m0, C0) for the named initialization in R^n.

    Args:
        name : one of INIT_NAMES
        n    : dimension (positive integer)

    Returns:
        (m0, C0) : copies of the initial mean (shape (n,)) and
                   covariance (shape (n, n)), both float64

    Raises:
        ValueError : if name is not recognised
    """
    if name == "mean_only":
        # Large mean offset; covariance already at target.
        # Tests whether mean converges without covariance interference.
        r  = 3.0
        m0 = r * np.ones(n, dtype=np.float64) / np.sqrt(n)
        C0 = np.eye(n, dtype=np.float64)

    elif name == "volume_high":
        # Covariance too large (volume 4^n times target).
        m0 = np.zeros(n, dtype=np.float64)
        C0 = 4.0 * np.eye(n, dtype=np.float64)

    elif name == "volume_low":
        # Covariance too small (volume (1/4)^n times target).
        m0 = np.zeros(n, dtype=np.float64)
        C0 = 0.25 * np.eye(n, dtype=np.float64)

    elif name == "shape_only":
        # Anisotropic covariance; det(C0) = exp(r)*exp(-r)*1*...*1 = 1,
        # so volume is already correct; only shape (eigenvalue spread) differs.
        r  = 2.0
        m0 = np.zeros(n, dtype=np.float64)
        C0 = np.diag(_anisotropic_diag(n, r))

    elif name == "mixed":
        # Non-zero mean + scaled anisotropic covariance.
        r  = 1.5
        s  = 2.0
        m0 = 2.0 * np.ones(n, dtype=np.float64) / np.sqrt(n)
        C0 = s * np.diag(_anisotropic_diag(n, r))

    else:
        raise ValueError(
            f"Unknown initialization '{name}'. "
            f"Choose from: {INIT_NAMES}."
        )

    return m0.copy(), C0.copy()
