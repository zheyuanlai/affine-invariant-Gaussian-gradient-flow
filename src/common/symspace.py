"""Frobenius-isometric vectorization of symmetric matrices.

A symmetric ``N_theta x N_theta`` matrix has ``N_theta (N_theta + 1) / 2``
independent entries. We store the upper triangle (row-major) with the
convention:

* diagonal entries are stored directly,
* off-diagonal entries are stored scaled by ``sqrt(2)``.

This makes the map an isometry between the Frobenius inner product on symmetric
matrices and the Euclidean inner product on vectors::

    dot(sym_to_vec(X), sym_to_vec(Y)) == Tr(X @ Y)
    norm(sym_to_vec(X))               == ||X||_F

which is exactly what the natural-gradient operator code relies on so that
``scipy.sparse.linalg.eigsh`` sees a genuinely symmetric matrix.
"""
from __future__ import annotations

import numpy as np

_SQRT2 = np.sqrt(2.0)


def sym_dim(N_theta):
    """Dimension of the space of symmetric ``N_theta x N_theta`` matrices."""
    return N_theta * (N_theta + 1) // 2


def _triu_offdiag_mask(N_theta):
    rows, cols = np.triu_indices(N_theta)
    return rows, cols, rows != cols


def sym_to_vec(X):
    """Vectorize a symmetric matrix (sqrt(2) scaling on off-diagonals)."""
    X = np.asarray(X, dtype=np.float64)
    N = X.shape[0]
    Xs = 0.5 * (X + X.T)
    rows, cols, offdiag = _triu_offdiag_mask(N)
    v = Xs[rows, cols].astype(np.float64, copy=True)
    v[offdiag] *= _SQRT2
    return v


def vec_to_sym(v, N_theta):
    """Inverse of :func:`sym_to_vec`."""
    v = np.asarray(v, dtype=np.float64)
    rows, cols, offdiag = _triu_offdiag_mask(N_theta)
    vals = v.astype(np.float64, copy=True)
    vals[offdiag] /= _SQRT2
    X = np.zeros((N_theta, N_theta), dtype=np.float64)
    X[rows, cols] = vals
    X[cols, rows] = vals  # diagonal is written twice with the same value
    return X


def sym_inner(X, Y):
    """Frobenius inner product ``Tr(X @ Y)`` for symmetric ``X``, ``Y``."""
    X = symmetrize_(X)
    Y = symmetrize_(Y)
    return float(np.sum(X * Y))


def sym_norm(X):
    """Frobenius norm ``||X||_F`` (equals ``norm(sym_to_vec(X))``)."""
    return float(np.sqrt(sym_inner(X, X)))


def random_symmetric(N_theta, rng, normalize=False):
    """Random symmetric matrix; if ``normalize`` scale to unit Frobenius norm."""
    A = rng.standard_normal((N_theta, N_theta))
    X = 0.5 * (A + A.T)
    if normalize:
        nrm = sym_norm(X)
        if nrm > 0.0:
            X = X / nrm
    return X


def pack_tangent(u, X):
    """Pack a tangent vector ``(u, X)`` into a single 1-D array ``[u, vec(X)]``."""
    u = np.asarray(u, dtype=np.float64).ravel()
    return np.concatenate([u, sym_to_vec(X)])


def unpack_tangent(w, N_theta):
    """Inverse of :func:`pack_tangent`; returns ``(u, X)``."""
    w = np.asarray(w, dtype=np.float64)
    u = w[:N_theta].copy()
    X = vec_to_sym(w[N_theta:], N_theta)
    return u, X


def tangent_dim(N_theta):
    """Dimension of a packed tangent vector ``(u, X)``."""
    return N_theta + sym_dim(N_theta)


# ---------------------------------------------------------------------------
# Fisher--Rao-whitened packing of tangent vectors (u, X)
# ---------------------------------------------------------------------------
#
# The local Fisher--Rao inner product at equilibrium ``a_star = (0, I)`` is
#
#     <(u1, X1), (u2, X2)>_star = u1 . u2 + 0.5 Tr(X1 X2).
#
# Packing a tangent vector as ``y = (u, sym_to_vec(X) / sqrt(2))`` makes the
# *Euclidean* dot product of packed vectors equal this Fisher--Rao inner
# product (because ``sym_to_vec`` is a Frobenius isometry and the ``1/sqrt(2)``
# turns ``Tr(X1 X2)`` into ``0.5 Tr(X1 X2)``). An operator written in these
# packed coordinates is therefore a genuinely symmetric matrix exactly when it
# is Fisher--Rao self-adjoint, which is what ``scipy.sparse.linalg.eigsh``
# requires for the local-rate eigenproblem.

def pack_tangent_fr(u, X):
    """Pack ``(u, X)`` into Fisher--Rao-whitened coords ``y = (u, vec(X)/sqrt2)``."""
    u = np.asarray(u, dtype=np.float64).ravel()
    return np.concatenate([u, sym_to_vec(X) / _SQRT2])


def unpack_tangent_fr(y, N_theta):
    """Inverse of :func:`pack_tangent_fr`; returns ``(u, X)``."""
    y = np.asarray(y, dtype=np.float64)
    u = y[:N_theta].copy()
    X = vec_to_sym(_SQRT2 * y[N_theta:], N_theta)
    return u, X


def fisher_rao_inner(t1, t2):
    """Fisher--Rao inner product of tangents ``t = (u, X)`` at equilibrium."""
    u1, X1 = t1
    u2, X2 = t2
    u1 = np.asarray(u1, dtype=np.float64).ravel()
    u2 = np.asarray(u2, dtype=np.float64).ravel()
    return float(u1 @ u2) + 0.5 * sym_inner(X1, X2)


def fisher_rao_norm(t):
    """Fisher--Rao norm ``sqrt(||u||^2 + 0.5 ||X||_F^2)`` of a tangent ``(u, X)``."""
    return float(np.sqrt(fisher_rao_inner(t, t)))


def symmetrize_(A):
    A = np.asarray(A, dtype=np.float64)
    return 0.5 * (A + A.T)
