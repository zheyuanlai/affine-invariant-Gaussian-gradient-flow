"""Matrix-free natural-gradient local operators ``T``, ``T*``, ``H_lin`` and the
linearized positive generator ``L_star = -J_star``.

In equilibrium-whitened coordinates with ``Z ~ N(0, I)`` and
``H_j = Hess V(Z_j)`` (the equilibrium ``theta = Z`` since ``a_star = (0, I)``):

    T[X]        = mean_j Tr(X H_j) Z_j                       (vector)
    T_star[u]   = mean_j (Z_j . u) H_j                       (sym matrix)
    H_lin[X]    = mean_j (Z_j^T X Z_j - Tr X) H_j            (sym matrix)
    H_lin*[X]   = mean_j Tr(H_j X) (Z_j Z_j^T - I)           (Frobenius adjoint)

    L_star(u, X) = ( u + 0.5 T[X],  X + T_star[u] + 0.5 H_sym[X] )

with ``H_sym = 0.5 (H_lin + H_lin*)``. On a *fixed* sample bank ``T`` and
``T_star`` are exact Frobenius adjoints, so the only operator whose MC estimator
is not automatically self-adjoint is ``H_lin``; symmetrizing it makes ``L_star``
exactly self-adjoint with respect to the Fisher--Rao inner product at
equilibrium, ``<(u1,X1),(u2,X2)> = u1.u2 + 0.5 Tr(X1 X2)``, which is required for
``scipy.sparse.linalg.eigsh``.

Memory: ``H_j`` has shape ``(M, N, N)``. The :class:`HessianBank` precomputes it
once when it fits (so repeated ``eigsh`` matvecs reuse it) and otherwise
recomputes it chunk-by-chunk; ``chunk_size`` bounds the per-call working set.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import LinearOperator

from src.common.spd import symmetrize
from src.common.symspace import (
    sym_dim, sym_to_vec, vec_to_sym, tangent_dim,
    sym_inner, random_symmetric, pack_tangent_fr, unpack_tangent_fr,
)


# ---------------------------------------------------------------------------
# Low-level applications on an explicit (hess_batch, Z)
# ---------------------------------------------------------------------------

def apply_T(X, hess_batch, Z):
    """``T[X] = mean_j Tr(X H_j) Z_j`` -> vector ``(N,)``."""
    Xs = symmetrize(X)
    d = np.einsum("ij,mij->m", Xs, hess_batch, optimize=True)   # Tr(X H_j)
    return (Z.T @ d) / Z.shape[0]


def apply_T_star(u, hess_batch, Z):
    """``T_star[u] = mean_j (Z_j . u) H_j`` -> symmetric matrix ``(N, N)``."""
    e = Z @ np.asarray(u, dtype=np.float64)                     # Z_j . u
    S = np.tensordot(e, hess_batch, axes=([0], [0])) / Z.shape[0]
    return symmetrize(S)


def apply_H_lin(X, hess_batch, Z):
    """``H_lin[X] = mean_j (Z_j^T X Z_j - Tr X) H_j`` -> symmetric matrix."""
    Xs = symmetrize(X)
    trX = np.trace(Xs)
    q = np.einsum("mi,ij,mj->m", Z, Xs, Z, optimize=True)       # Z_j^T X Z_j
    S = np.tensordot(q - trX, hess_batch, axes=([0], [0])) / Z.shape[0]
    return symmetrize(S)


def apply_H_lin_adjoint(X, hess_batch, Z):
    """Frobenius adjoint ``H_lin*[X] = mean_j Tr(H_j X) (Z_j Z_j^T - I)``."""
    Xs = symmetrize(X)
    d = np.einsum("ij,mij->m", Xs, hess_batch, optimize=True)   # Tr(H_j X)
    M = Z.shape[0]
    zz = np.einsum("m,mi,mj->ij", d, Z, Z, optimize=True) / M
    return symmetrize(zz - (d.mean()) * np.eye(Z.shape[1]))


def apply_H_sym(X, hess_batch, Z):
    """Self-adjoint estimator ``0.5 (H_lin + H_lin*)[X]``."""
    return 0.5 * (apply_H_lin(X, hess_batch, Z) + apply_H_lin_adjoint(X, hess_batch, Z))


# ---------------------------------------------------------------------------
# Hessian bank (precompute-or-chunk)
# ---------------------------------------------------------------------------

class HessianBank:
    """Provides ``(Z_chunk, Hess_chunk)`` pairs, precomputed when memory allows."""

    def __init__(self, potential, Z, chunk_size=None, max_precompute_bytes=2e9):
        self.potential = potential
        self.Z = np.ascontiguousarray(Z, dtype=np.float64)
        self.M, self.N = self.Z.shape
        self.chunk_size = int(chunk_size) if chunk_size else self.M
        nbytes = self.M * self.N * self.N * 8
        self.precomputed = nbytes <= max_precompute_bytes
        self._H = potential.batch_hess(self.Z) if self.precomputed else None

    def chunks(self):
        cs = self.chunk_size
        for s in range(0, self.M, cs):
            Zc = self.Z[s:s + cs]
            Hc = self._H[s:s + cs] if self.precomputed else self.potential.batch_hess(Zc)
            yield Zc, Hc


# ---------------------------------------------------------------------------
# Chunked applications using a HessianBank
# ---------------------------------------------------------------------------

def _hforward_from_bank(bank, X):
    """Raw forward estimator ``H_lin[X]`` accumulated over the bank's chunks.

    Returned *symmetrized for storage* (the output of ``H_lin`` is a symmetric
    matrix even though the operator itself is not Frobenius self-adjoint on a
    finite bank).
    """
    Xs = symmetrize(X)
    trX = np.trace(Xs)
    N, M = bank.N, bank.M
    hlin = np.zeros((N, N))
    for Zc, Hc in bank.chunks():
        q = np.einsum("mi,ij,mj->m", Zc, Xs, Zc, optimize=True)  # Z_j^T X Z_j
        hlin += np.tensordot(q - trX, Hc, axes=([0], [0]))
    return symmetrize(hlin / M)


def _hsym_from_bank(bank, X):
    """``H_sym[X]`` accumulated over the bank's chunks."""
    Xs = symmetrize(X)
    trX = np.trace(Xs)
    N, M = bank.N, bank.M
    hlin = np.zeros((N, N))
    adj_zz = np.zeros((N, N))
    adj_d = 0.0
    for Zc, Hc in bank.chunks():
        d = np.einsum("ij,mij->m", Xs, Hc, optimize=True)        # Tr(X H_j)
        q = np.einsum("mi,ij,mj->m", Zc, Xs, Zc, optimize=True)  # Z_j^T X Z_j
        hlin += np.tensordot(q - trX, Hc, axes=([0], [0]))
        adj_zz += np.einsum("m,mi,mj->ij", d, Zc, Zc, optimize=True)
        adj_d += float(d.sum())
    H_lin = hlin / M
    H_adj = adj_zz / M - (adj_d / M) * np.eye(N)
    return symmetrize(0.5 * (H_lin + H_adj))


def _diagonal_A_from_bank(bank):
    """Diagonal-mode coefficient matrix ``A_ij = mean_j Hess V_ii(Z_j) (Z_j[i]^2 - 1)``.

    Equivalently ``A = G - 1 1^T`` with ``G_ij = E[Hess V_ii(Z) Z_j^2]``. For a
    *diagonal* perturbation ``X = diag(lambda)`` the diagonal of ``H_lin[X]`` is
    ``A lambda``, so ``lambda_max`` of (the symmetric part of) ``A`` is the
    diagonal-restricted operator norm. ``A`` is exactly diagonal in expectation
    for separable potentials; off-diagonal mass measures finite-sample leakage.
    """
    N, M = bank.N, bank.M
    A = np.zeros((N, N))
    for Zc, Hc in bank.chunks():
        diagH = np.einsum("mii->mi", Hc)        # (c, N): Hess V_ii(Z_j)
        Z2m1 = Zc * Zc - 1.0                     # (c, N): Z_j[k]^2 - 1
        A += diagH.T @ Z2m1                      # sum_j diagH[j,i] * Z2m1[j,k]
    return A / M


def _lstar_from_bank(bank, u, X):
    """Return ``L_star(u, X) = (u + 0.5 T[X], X + T_star[u] + 0.5 H_sym[X])``."""
    Xs = symmetrize(X)
    u = np.asarray(u, dtype=np.float64)
    trX = np.trace(Xs)
    N, M = bank.N, bank.M
    t_sum = np.zeros(N)
    tstar = np.zeros((N, N))
    hlin = np.zeros((N, N))
    adj_zz = np.zeros((N, N))
    adj_d = 0.0
    for Zc, Hc in bank.chunks():
        d = np.einsum("ij,mij->m", Xs, Hc, optimize=True)        # Tr(X H_j)
        t_sum += Zc.T @ d
        e = Zc @ u                                               # Z_j . u
        tstar += np.tensordot(e, Hc, axes=([0], [0]))
        q = np.einsum("mi,ij,mj->m", Zc, Xs, Zc, optimize=True)  # Z_j^T X Z_j
        hlin += np.tensordot(q - trX, Hc, axes=([0], [0]))
        adj_zz += np.einsum("m,mi,mj->ij", d, Zc, Zc, optimize=True)
        adj_d += float(d.sum())
    T_X = t_sum / M
    T_star_u = symmetrize(tstar / M)
    H_lin = hlin / M
    H_adj = adj_zz / M - (adj_d / M) * np.eye(N)
    H_sym = symmetrize(0.5 * (H_lin + H_adj))
    u_out = u + 0.5 * T_X
    X_out = symmetrize(Xs + T_star_u + 0.5 * H_sym)
    return u_out, X_out


def apply_L_star(potential, u, X, Z, chunk_size=None, max_precompute_bytes=2e9):
    """Apply ``L_star`` once (builds a one-shot HessianBank)."""
    bank = HessianBank(potential, Z, chunk_size=chunk_size,
                       max_precompute_bytes=max_precompute_bytes)
    return _lstar_from_bank(bank, u, X)


# ---------------------------------------------------------------------------
# LinearOperators for eigsh
# ---------------------------------------------------------------------------

def make_H_linear_operator(potential, Z, chunk_size=None, max_precompute_bytes=2e9):
    """Self-adjoint ``H_sym`` as a ``LinearOperator`` on ``sym_to_vec(X)``."""
    N = potential.N_theta
    d = sym_dim(N)
    bank = HessianBank(potential, Z, chunk_size=chunk_size,
                       max_precompute_bytes=max_precompute_bytes)

    def matvec(v):
        X = vec_to_sym(v, N)
        return sym_to_vec(_hsym_from_bank(bank, X))

    return LinearOperator((d, d), matvec=matvec, rmatvec=matvec, dtype=np.float64)


def make_H_forward_operator(potential, Z, chunk_size=None, max_precompute_bytes=2e9):
    """Raw forward ``H_lin`` as a ``LinearOperator`` on ``sym_to_vec(X)``.

    The matvec is the forward estimator and the rmatvec is its true Frobenius
    adjoint, so this operator is *not* symmetric on a finite bank (that is the
    point: it is the uncorrected estimator kept for diagnostics). Use a
    non-symmetric eigensolver on it.
    """
    N = potential.N_theta
    d = sym_dim(N)
    bank = HessianBank(potential, Z, chunk_size=chunk_size,
                       max_precompute_bytes=max_precompute_bytes)

    def matvec(v):
        return sym_to_vec(_hforward_from_bank(bank, vec_to_sym(v, N)))

    def rmatvec(v):
        # Frobenius adjoint: H_sym = 0.5 (forward + adjoint) -> adjoint = 2 H_sym - forward.
        X = vec_to_sym(v, N)
        return sym_to_vec(symmetrize(2.0 * _hsym_from_bank(bank, X)
                                     - _hforward_from_bank(bank, X)))

    return LinearOperator((d, d), matvec=matvec, rmatvec=rmatvec, dtype=np.float64)


def apply_H_sym_banked(potential, X, Z, chunk_size=None, max_precompute_bytes=2e9):
    """Apply the self-adjoint ``H_sym`` once (chunked); returns a symmetric matrix."""
    bank = HessianBank(potential, Z, chunk_size=chunk_size,
                       max_precompute_bytes=max_precompute_bytes)
    return _hsym_from_bank(bank, X)


def diagonal_A_matrix(potential, Z, chunk_size=None, max_precompute_bytes=2e9):
    """Diagonal-mode coefficient matrix ``A = G - 1 1^T`` (shape ``(N, N)``)."""
    bank = HessianBank(potential, Z, chunk_size=chunk_size,
                       max_precompute_bytes=max_precompute_bytes)
    return _diagonal_A_from_bank(bank)


def make_L_star_operator(potential, Z, covariance_weight="half",
                         chunk_size=None, max_precompute_bytes=2e9):
    """``L_star`` as a ``LinearOperator`` on a packed tangent ``(u, X)``.

    ``covariance_weight="half"`` returns the operator in Fisher--Rao-whitened
    coordinates ``v = (u, (1/sqrt2) sym_to_vec(X))`` so that the Euclidean inner
    product equals the Fisher--Rao inner product and the matrix is *symmetric*
    (use this for ``eigsh``). ``"plain"`` returns the raw operator on
    ``(u, sym_to_vec(X))`` (not symmetric; for diagnostics only).
    """
    N = potential.N_theta
    dim = tangent_dim(N)
    bank = HessianBank(potential, Z, chunk_size=chunk_size,
                       max_precompute_bytes=max_precompute_bytes)

    if covariance_weight == "half":
        def matvec(v):
            u, X = unpack_tangent_fr(v, N)         # X = sqrt2 * vec_to_sym(v[N:])
            u_out, X_out = _lstar_from_bank(bank, u, X)
            return pack_tangent_fr(u_out, X_out)   # (u_out, vec(X_out)/sqrt2)
    elif covariance_weight == "plain":
        def matvec(v):
            u = v[:N]
            X = vec_to_sym(v[N:], N)
            u_out, X_out = _lstar_from_bank(bank, u, X)
            out = np.empty_like(v)
            out[:N] = u_out
            out[N:] = sym_to_vec(X_out)
            return out
    else:
        raise ValueError(f"covariance_weight must be 'half' or 'plain', got {covariance_weight!r}")

    return LinearOperator((dim, dim), matvec=matvec, rmatvec=matvec, dtype=np.float64)


def unpack_weighted_eigenvector(v, N_theta):
    """Recover ``(u, X)`` from a ``covariance_weight='half'`` eigenvector."""
    return unpack_tangent_fr(v, N_theta)


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def adjoint_residual(potential, Z, n_probe=4, seed=0, chunk_size=None):
    """Relative residual of the adjoint identity ``<T[X], u> == <X, T_star[u]>``.

    This identity is exact on a fixed bank, so the residual reflects only
    floating-point error; it is a guard against indexing/scaling mistakes.
    """
    rng = np.random.default_rng(seed)
    N = potential.N_theta
    H = potential.batch_hess(np.ascontiguousarray(Z))
    worst = 0.0
    for _ in range(n_probe):
        u = rng.standard_normal(N)
        A = rng.standard_normal((N, N))
        X = 0.5 * (A + A.T)
        lhs = float(np.dot(apply_T(X, H, Z), u))
        rhs = sym_inner(X, apply_T_star(u, H, Z))
        denom = max(abs(lhs), abs(rhs), 1e-300)
        worst = max(worst, abs(lhs - rhs) / denom)
    return worst


def _relative_self_adjoint_error(apply_op, N, n_probe, seed, reduce):
    """Relative self-adjointness error of ``apply_op`` over random symmetric pairs.

    ``error = |<X, op[Y]> - <op[X], Y>| / max(1, |<X, op[Y]>|, |<op[X], Y>|)``
    reduced over ``n_probe`` random symmetric ``(X, Y)`` pairs (``reduce`` is
    ``"max"`` or ``"mean"``). Near machine precision for a self-adjoint ``op``.
    """
    rng = np.random.default_rng(seed)
    errs = []
    for _ in range(n_probe):
        X = random_symmetric(N, rng)
        Y = random_symmetric(N, rng)
        a = sym_inner(X, apply_op(Y))
        b = sym_inner(apply_op(X), Y)
        errs.append(abs(a - b) / max(1.0, abs(a), abs(b)))
    return float(np.max(errs) if reduce == "max" else np.mean(errs))


def self_adjoint_error_H(potential, Z, estimator="symmetrized", n_probe=8,
                         seed=0, chunk_size=None, reduce="max",
                         max_precompute_bytes=2e9):
    """Relative Frobenius self-adjointness error of the ``H`` estimator.

    ``estimator="raw_forward"`` probes the uncorrected ``H_lin`` (generally
    nonzero on a finite bank); ``"symmetrized"`` probes ``H_sym`` (near machine
    precision by construction).
    """
    N = potential.N_theta
    bank = HessianBank(potential, Z, chunk_size=chunk_size,
                       max_precompute_bytes=max_precompute_bytes)
    if estimator == "raw_forward":
        apply_op = lambda X: _hforward_from_bank(bank, X)
    elif estimator == "symmetrized":
        apply_op = lambda X: _hsym_from_bank(bank, X)
    else:
        raise ValueError(f"unknown estimator {estimator!r}")
    return _relative_self_adjoint_error(apply_op, N, n_probe, seed, reduce)


def self_adjoint_error_L_star(potential, Z, n_probe=8, seed=0, chunk_size=None,
                              reduce="max", max_precompute_bytes=2e9):
    """Relative self-adjointness error of ``L_star`` in Fisher--Rao packed coords.

    Probes ``|dot(y1, L y2) - dot(L y1, y2)| / max(1, ...)`` for random packed
    Euclidean vectors ``y``; near machine precision iff ``L_star`` is genuinely
    symmetric in these coordinates (the precondition for ``eigsh``).
    """
    op = make_L_star_operator(potential, Z, covariance_weight="half",
                              chunk_size=chunk_size,
                              max_precompute_bytes=max_precompute_bytes)
    dim = op.shape[0]
    rng = np.random.default_rng(seed)
    errs = []
    for _ in range(n_probe):
        y1 = rng.standard_normal(dim)
        y2 = rng.standard_normal(dim)
        a = float(y1 @ op.matvec(y2))
        b = float(op.matvec(y1) @ y2)
        errs.append(abs(a - b) / max(1.0, abs(a), abs(b)))
    return float(np.max(errs) if reduce == "max" else np.mean(errs))
