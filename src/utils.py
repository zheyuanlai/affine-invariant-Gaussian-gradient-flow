"""
Utility functions for SPD matrix operations and parameter validation.

All matrices are float64.  Eigenvalues are clamped to >= 1e-300 to keep
logarithms finite; in practice this only matters for degenerate test states.
"""
import numpy as np


# ---------------------------------------------------------------------------
# Eigendecomposition
# ---------------------------------------------------------------------------

def spd_eigh(C):
    """Eigendecomposition of a symmetric positive-(semi)definite matrix.

    Returns (eigvals, eigvecs) where eigvals >= 1e-300.
    numpy.linalg.eigh is used (exploits symmetry, returns sorted eigvals).
    """
    eigvals, eigvecs = np.linalg.eigh(C)
    eigvals = np.maximum(eigvals, 1e-300)
    return eigvals, eigvecs


# ---------------------------------------------------------------------------
# Parameter validation
# ---------------------------------------------------------------------------

def validate_params(omega, tau, n):
    """Return True iff (omega, tau, n) satisfy the well-posedness constraints.

    Constraints:
        omega > 0
        omega + n * tau > 0   (denominator in the covariance update)
    """
    if omega <= 0.0:
        return False
    if omega + n * tau <= 0.0:
        return False
    return True


# ---------------------------------------------------------------------------
# Cosine expectation (closed form for Gaussian)
# ---------------------------------------------------------------------------

def cosine_expectation(m, C, q, b):
    """Compute E_{theta ~ N(m, C)}[ cos(q^T theta + b) ].

    Closed form: exp(-0.5 * q^T C q) * cos(q^T m + b).
    """
    qCq = float(q @ C @ q)
    return float(np.exp(-0.5 * qCq) * np.cos(float(q @ m) + b))


# ---------------------------------------------------------------------------
# Test-function direction vector
# ---------------------------------------------------------------------------

def make_q_vector(n):
    """Return the normalised test direction q = (1, 2, ..., n) / ||.||_2."""
    q = np.arange(1, n + 1, dtype=np.float64)
    return q / np.linalg.norm(q)
