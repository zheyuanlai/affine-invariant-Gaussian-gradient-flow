"""Symmetric / symmetric-positive-definite (SPD) matrix utilities.

All routines assume real symmetric inputs but symmetrize defensively. Matrix
square roots, inverse square roots, and exponentials of symmetric matrices are
computed through the eigendecomposition (numerically stable and exactly
symmetric in exact arithmetic). All computation is in float64.
"""
from __future__ import annotations

import numpy as np

# Default eigenvalue tolerance used to decide positive-definiteness.
DEFAULT_SPD_TOL = 1e-10


def symmetrize(A):
    """Return the symmetric part ``0.5 (A + A^T)`` as a float64 array."""
    A = np.asarray(A, dtype=np.float64)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"symmetrize expects a square matrix, got shape {A.shape}")
    return 0.5 * (A + A.T)


def is_spd(A, tol=DEFAULT_SPD_TOL):
    """Return True if (the symmetric part of) ``A`` is SPD with min eig > tol."""
    try:
        w = np.linalg.eigvalsh(symmetrize(A))
    except (np.linalg.LinAlgError, ValueError):
        return False
    return bool(w[0] > tol)


def eigh_spd(A, floor=None):
    """Eigendecomposition of the symmetric part of ``A``.

    Returns ``(w, V)`` with ascending eigenvalues ``w`` and orthonormal columns
    ``V`` (so ``A ~= V diag(w) V^T``). If ``floor`` is not None, eigenvalues are
    clipped from below to ``floor`` (useful before a sqrt / log).
    """
    w, V = np.linalg.eigh(symmetrize(A))
    if floor is not None:
        w = np.maximum(w, floor)
    return w, V


def check_spd(A, name="C", tol=DEFAULT_SPD_TOL):
    """Raise ``ValueError`` unless ``A`` is SPD; return ``(min_eig, max_eig)``."""
    A = np.asarray(A, dtype=np.float64)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"{name} must be a square matrix, got shape {A.shape}")
    asym = float(np.max(np.abs(A - A.T))) if A.size else 0.0
    w = np.linalg.eigvalsh(symmetrize(A))
    if w[0] <= tol:
        raise ValueError(
            f"{name} is not SPD: min eigenvalue {w[0]:.3e} <= tol {tol:.3e} "
            f"(max|A - A^T| = {asym:.3e})"
        )
    return float(w[0]), float(w[-1])


def symmetric_sqrt(A, tol=DEFAULT_SPD_TOL):
    """Symmetric (principal) square root of a PSD symmetric matrix."""
    w, V = eigh_spd(A)
    if w[0] < -np.abs(tol):
        raise ValueError(
            f"symmetric_sqrt: matrix is not PSD (min eigenvalue {w[0]:.3e})"
        )
    sqrt_w = np.sqrt(np.clip(w, 0.0, None))
    return (V * sqrt_w) @ V.T


def symmetric_invsqrt(A, tol=DEFAULT_SPD_TOL):
    """Symmetric inverse square root of an SPD symmetric matrix."""
    w, V = eigh_spd(A)
    if w[0] <= np.abs(tol):
        raise ValueError(
            f"symmetric_invsqrt: matrix is not SPD (min eigenvalue {w[0]:.3e})"
        )
    inv_sqrt_w = 1.0 / np.sqrt(w)
    return (V * inv_sqrt_w) @ V.T


def symmetric_expm(A):
    """Matrix exponential of a symmetric matrix via its eigendecomposition.

    For symmetric ``A`` this equals ``scipy.linalg.expm(A)`` up to rounding and
    is guaranteed symmetric (and SPD).
    """
    w, V = np.linalg.eigh(symmetrize(A))
    return (V * np.exp(w)) @ V.T
