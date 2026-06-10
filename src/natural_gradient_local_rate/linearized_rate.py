"""Estimate the local rate ``gamma_loc`` and operator norm ``Lambda_hat``.

* ``Lambda_hat`` = largest eigenvalue of the self-adjoint ``H_sym`` (the
  symmetric-matrix operator behind the dimension-free conjecture).
* ``gamma_loc`` = smallest eigenvalue of the linearized positive generator
  ``L_star`` (in Fisher--Rao-whitened coordinates), the local convergence rate.

For small operators we form the dense matrix and use ``numpy.linalg.eigh``
(robust and exact up to MC error); for larger operators we use
``scipy.sparse.linalg.eigsh`` matrix-free, with a shift so that the *smallest*
eigenvalue is found as the *largest* of a shifted operator (more reliable for
Lanczos than ``which='SA'``).
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

from src.common.symspace import sym_norm, sym_to_vec, vec_to_sym
from src.natural_gradient_local_rate.operators import (
    make_H_linear_operator, make_H_forward_operator, make_L_star_operator,
    diagonal_A_matrix, T_matrix, unpack_weighted_eigenvector,
)

# Operators with dimension <= this are diagonalized densely (exact, robust).
DENSE_THRESHOLD = 32
# Raw-forward operators with dimension <= this are diagonalized densely (a
# non-symmetric dense eig is reliable); larger ones fall back to ARPACK ``eigs``.
RAW_DENSE_LIMIT = 1024


def _to_matrix(linop):
    """Materialize ``linop`` as a dense matrix (no symmetrization)."""
    dim = linop.shape[0]
    cols = [linop.matvec(np.eye(dim)[:, i]) for i in range(dim)]
    return np.column_stack(cols)


def _to_dense(linop):
    Mat = _to_matrix(linop)
    return 0.5 * (Mat + Mat.T)


def _largest_real_eig(linop, tol, maxiter, dense_limit=RAW_DENSE_LIMIT):
    """Largest real part of the spectrum of a (possibly non-symmetric) operator.

    Dense ``numpy.linalg.eigvals`` below ``dense_limit`` (reliable); ARPACK
    ``scipy.sparse.linalg.eigs(which='LR')`` above it, returning ``NaN`` if it
    fails to converge (this is a diagnostic quantity, never a reported rate).
    """
    dim = linop.shape[0]
    if dim <= dense_limit:
        return float(np.max(np.real(np.linalg.eigvals(_to_matrix(linop)))))
    from scipy.sparse.linalg import eigs
    try:
        vals = eigs(linop, k=1, which="LR", tol=tol, maxiter=maxiter,
                    return_eigenvectors=False)
        return float(np.max(np.real(vals)))
    except Exception:
        return float("nan")


def _largest_eigh(linop, tol, maxiter):
    dim = linop.shape[0]
    if dim <= DENSE_THRESHOLD:
        w, V = np.linalg.eigh(_to_dense(linop))
        return float(w[-1]), V[:, -1]
    vals, vecs = eigsh(linop, k=1, which="LA", tol=tol, maxiter=maxiter)
    return float(vals[0]), vecs[:, 0]


def _smallest_eigh(linop, tol, maxiter):
    """Smallest eigenvalue via a shift: largest eig of ``(c I - L)`` -> ``c - .``."""
    dim = linop.shape[0]
    if dim <= DENSE_THRESHOLD:
        w, V = np.linalg.eigh(_to_dense(linop))
        return float(w[0]), V[:, 0]
    # shift past the top of the spectrum, then find the largest of (c I - L).
    top, _ = _largest_eigh(linop, tol, maxiter)
    shift = top * 1.0 + 1.0

    def shifted_matvec(v):
        return shift * v - linop.matvec(v)

    shifted = LinearOperator((dim, dim), matvec=shifted_matvec,
                             rmatvec=shifted_matvec, dtype=np.float64)
    vals, vecs = eigsh(shifted, k=1, which="LA", tol=tol, maxiter=maxiter)
    return float(shift - vals[0]), vecs[:, 0]


def estimate_lambda_hat(potential, Z, estimator="symmetrized", eigsh_tol=1e-6,
                        eigsh_maxiter=1000, chunk_size=None,
                        return_eigenvector=False, raw_dense_limit=RAW_DENSE_LIMIT):
    """Largest-eigenvalue operator norm ``Lambda_hat``.

    ``estimator="symmetrized"`` (default) returns ``lambda_max`` of the
    self-adjoint ``H_sym`` (and optionally the eigen-matrix ``X_eig``).
    ``estimator="raw_forward"`` returns the largest real part of the spectrum of
    the uncorrected, non-self-adjoint forward ``H_lin`` (diagnostic only; no
    eigenvector).
    """
    if estimator == "symmetrized":
        op = make_H_linear_operator(potential, Z, chunk_size=chunk_size)
        lam, vec = _largest_eigh(op, eigsh_tol, eigsh_maxiter)
        if return_eigenvector:
            return lam, vec_to_sym(vec, potential.N_theta)
        return lam
    if estimator == "raw_forward":
        op = make_H_forward_operator(potential, Z, chunk_size=chunk_size)
        lam = _largest_real_eig(op, eigsh_tol, eigsh_maxiter, dense_limit=raw_dense_limit)
        if return_eigenvector:
            return lam, None
        return lam
    raise ValueError(f"unknown estimator {estimator!r} "
                     "(expected 'symmetrized' or 'raw_forward')")


def estimate_diagonal_lambda(potential, Z, chunk_size=None):
    """Diagonal-restricted operator norm and finite-sample off-diagonal leakage.

    Returns a dict with ``Lambda_hat_diag`` (``lambda_max`` of the symmetric
    part of ``A = G - 1 1^T``, matching the manuscript's diagonal-mode quadratic
    form), ``diag_offdiag_norm`` (Frobenius norm of the off-diagonal part of
    ``A``; ~0 for separable potentials up to MC noise), ``max_diag``
    (``max_i A_ii``, the cleanest separable estimator), and the matrix ``A``.
    """
    A = diagonal_A_matrix(potential, Z, chunk_size=chunk_size)
    A_sym = 0.5 * (A + A.T)
    lam = float(np.linalg.eigvalsh(A_sym)[-1])
    offdiag = A - np.diag(np.diag(A))
    return {
        "Lambda_hat_diag": lam,
        "diag_offdiag_norm": float(np.linalg.norm(offdiag)),
        "max_diag": float(np.max(np.diag(A))),
        "A": A,
    }


def _decompose_T_top_mode(Tmat, N):
    """Decompose the top ``T`` singular mode into longitudinal/mixed/transverse blocks.

    For the top singular pair ``T[X] = tau_H w`` with ``||w||=||X||_F=1``, write

        ``X = a ww^T + wb^T + bw^T + B``, with ``b perpendicular to w`` and
        ``Bw = 0``.

    The returned contributions are the three summands of ``w^T T[X]``. They are
    diagnostic only: a large transverse contribution is the numerically relevant
    signal for the hard analytic subproblem in the manuscript.
    """
    if Tmat.size == 0:
        return {}
    U, S, Vt = np.linalg.svd(Tmat, full_matrices=False)
    if S.size == 0:
        return {}
    tau = float(S[0])
    w = U[:, 0]
    X = vec_to_sym(Vt[0, :], N)
    a = float(w @ X @ w)
    b = X @ w - a * w
    B = X - a * np.outer(w, w) - np.outer(w, b) - np.outer(b, w)
    B = 0.5 * (B + B.T)

    X_long = a * np.outer(w, w)
    X_mixed = np.outer(w, b) + np.outer(b, w)

    def contribution(Y):
        return float(w @ (Tmat @ sym_to_vec(Y)))

    total = contribution(X)
    denom = total if abs(total) > 1e-14 else np.nan
    long = contribution(X_long)
    mixed = contribution(X_mixed)
    trans = contribution(B)
    return {
        "tau_top_total": total,
        "tau_top_longitudinal": long,
        "tau_top_mixed": mixed,
        "tau_top_transverse": trans,
        "tau_top_longitudinal_fraction": float(long / denom) if np.isfinite(denom) else float("nan"),
        "tau_top_mixed_fraction": float(mixed / denom) if np.isfinite(denom) else float("nan"),
        "tau_top_transverse_fraction": float(trans / denom) if np.isfinite(denom) else float("nan"),
        "tau_top_X_longitudinal_norm_sq": float(a * a),
        "tau_top_X_mixed_norm_sq": float(2.0 * (b @ b)),
        "tau_top_X_transverse_norm_sq": float(sym_norm(B) ** 2),
        "tau_top_svd_gap": float(tau - (float(S[1]) if S.size > 1 else 0.0)),
    }


def estimate_tau_H(potential, Z, chunk_size=None):
    """Estimate the Hessian first-Hermite coupling ``tau_H = ||T||_op``.

    ``T`` is represented in Frobenius-isometric coordinates on ``Sym(N)``:

        ``T_mat @ sym_to_vec(X) = E[Tr(X Hess V(Z)) Z]``.

    The equivalent adjoint form is

        ``tau_H = sup_{||w||=1} ||E[Hess V(Z) (w^T Z)]||_F``.

    Returns a dict containing ``tau_H``, ``tau_H_sq`` and the deterministic
    coupling-only rate bound ``1 / (tau_H_sq + 3)`` from the revised
    Schur-complement reduction.
    """
    Tmat = T_matrix(potential, Z, chunk_size=chunk_size)
    # The nonzero singular values of T are the square-roots of eig(T T^T);
    # this is only N x N, so dense eigvalsh is cheap and stable.
    gram = Tmat @ Tmat.T
    tau_sq = max(0.0, float(np.linalg.eigvalsh(0.5 * (gram + gram.T))[-1]))
    out = {
        "tau_H": float(np.sqrt(tau_sq)),
        "tau_H_sq": tau_sq,
        "coupling_bound_rate": 1.0 / (tau_sq + 3.0),
        "T_matrix": Tmat,
    }
    out.update(_decompose_T_top_mode(Tmat, potential.N_theta))
    return out


def estimate_gamma_loc(potential, Z, eigsh_tol=1e-6, eigsh_maxiter=1000,
                       chunk_size=None, return_eigenvector=True):
    """Smallest eigenvalue ``gamma_loc`` of ``L_star`` and its slow mode.

    Returns ``(gamma_loc, (u_star, X_star))`` (the eigenpair is ``None`` when
    ``return_eigenvector`` is False).
    """
    op = make_L_star_operator(potential, Z, covariance_weight="half",
                              chunk_size=chunk_size)
    gam, vec = _smallest_eigh(op, eigsh_tol, eigsh_maxiter)
    if return_eigenvector:
        u_star, X_star = unpack_weighted_eigenvector(vec, potential.N_theta)
        return gam, (u_star, X_star)
    return gam, None
