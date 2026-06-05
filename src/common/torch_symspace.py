"""Torch equivalents of the Frobenius-isometric symmetric-matrix vectorization.

Convention matches :mod:`src.common.symspace` *exactly* (upper triangle,
row-major; diagonal stored directly; off-diagonal scaled by sqrt(2)), so that

    dot(torch_sym_to_vec(X), torch_sym_to_vec(Y)) == Tr(X @ Y)

and the Fisher--Rao packing ``y = (u, sym_to_vec(X)/sqrt2)`` agrees with the
NumPy path. The triangle ordering is taken from ``numpy.triu_indices`` so it is
identical to the CPU helpers bit-for-bit (up to floating point).

torch is imported lazily via :mod:`src.common.torch_utils`; importing this module
without torch installed is fine until a function is actually called.
"""
from __future__ import annotations

import math

import numpy as np

from src.common.torch_utils import get_torch

_SQRT2 = math.sqrt(2.0)


def torch_sym_dim(N_theta):
    return N_theta * (N_theta + 1) // 2


# Cache of (rows, cols, offdiag) upper-triangle index tensors, keyed by
# ``(N_theta, device_str, dtype-irrelevant)``. Rebuilding ``np.triu_indices`` and
# moving three small tensors to the device every chunk is pure overhead in the
# tight Monte-Carlo accumulation loop, so the indices are memoized per device.
_TRIU_CACHE = {}


def _triu(N_theta, device):
    """Row-major upper-triangle indices + off-diagonal mask as torch tensors.

    Memoized per ``(N_theta, device)`` so the tight accumulation loop reuses the
    same index tensors instead of rebuilding/transferring them every chunk.
    """
    torch = get_torch()
    key = (int(N_theta), str(device))
    cached = _TRIU_CACHE.get(key)
    if cached is not None:
        return cached
    rows_np, cols_np = np.triu_indices(N_theta)
    offdiag_np = rows_np != cols_np
    out = (torch.as_tensor(rows_np, dtype=torch.long, device=device),
           torch.as_tensor(cols_np, dtype=torch.long, device=device),
           torch.as_tensor(offdiag_np, dtype=torch.bool, device=device))
    _TRIU_CACHE[key] = out
    return out


def _symmetrize(X):
    return 0.5 * (X + X.transpose(-1, -2))


def torch_sym_to_vec(X):
    """Vectorize a symmetric matrix ``(N, N) -> (p,)`` (sqrt2 on off-diagonals)."""
    Xs = _symmetrize(X)
    N = Xs.shape[-1]
    rows, cols, offdiag = _triu(N, Xs.device)
    v = Xs[rows, cols].clone()
    v[offdiag] = v[offdiag] * _SQRT2
    return v


def torch_sym_to_vec_batch(Xb):
    """Batched vectorize ``(B, N, N) -> (B, p)``."""
    Xs = _symmetrize(Xb)
    N = Xs.shape[-1]
    rows, cols, offdiag = _triu(N, Xs.device)
    v = Xs[:, rows, cols].clone()
    v[:, offdiag] = v[:, offdiag] * _SQRT2
    return v


def torch_outer_minus_I_symvec(Zc, N_theta):
    """Direct ``sym_to_vec(z z^T - I)`` for a batch of vectors, ``(B, N) -> (B, p)``.

    Equivalent to ``torch_sym_to_vec_batch(einsum('bi,bj->bij', Zc, Zc) - I)`` but
    avoids ever materializing the ``(B, N, N)`` outer-product tensor: for the
    upper-triangle index pairs ``(i, j)`` the vectorized entry is

    * diagonal (``i == j``):    ``z_i^2 - 1``
    * off-diagonal (``i < j``): ``sqrt(2) * z_i * z_j``

    matching the Frobenius-isometric convention of :func:`torch_sym_to_vec_batch`
    (``z z^T - I`` is already symmetric, so no explicit symmetrization is needed).
    """
    rows, cols, offdiag = _triu(N_theta, Zc.device)
    # Zc[:, rows] and Zc[:, cols] are (B, p); their product gives z_i z_j per pair.
    v = Zc[:, rows] * Zc[:, cols]                 # (B, p)
    v[:, offdiag] = v[:, offdiag] * _SQRT2        # sqrt(2) on off-diagonals
    v[:, ~offdiag] = v[:, ~offdiag] - 1.0         # subtract Tr(I) part on diagonal
    return v


def torch_vec_to_sym(v, N_theta):
    """Inverse of :func:`torch_sym_to_vec` ``(p,) -> (N, N)``."""
    torch = get_torch()
    rows, cols, offdiag = _triu(N_theta, v.device)
    vals = v.clone()
    vals[offdiag] = vals[offdiag] / _SQRT2
    X = torch.zeros((N_theta, N_theta), dtype=v.dtype, device=v.device)
    X[rows, cols] = vals
    X[cols, rows] = vals  # diagonal written twice with the same value
    return X


def torch_sym_inner(X, Y):
    """Frobenius inner product ``Tr(X @ Y)`` for symmetric ``X``, ``Y`` (scalar tensor)."""
    return (_symmetrize(X) * _symmetrize(Y)).sum()


def torch_pack_tangent_fr(u, X):
    """Pack ``(u, X)`` into Fisher--Rao coords ``y = (u, sym_to_vec(X)/sqrt2)``."""
    torch = get_torch()
    return torch.cat([u.reshape(-1), torch_sym_to_vec(X) / _SQRT2])


def torch_unpack_tangent_fr(y, N_theta):
    """Inverse of :func:`torch_pack_tangent_fr`; returns ``(u, X)``."""
    u = y[:N_theta].clone()
    X = torch_vec_to_sym(_SQRT2 * y[N_theta:], N_theta)
    return u, X


def torch_fisher_rao_inner(u1, X1, u2, X2):
    """Fisher--Rao inner product ``u1.u2 + 0.5 Tr(X1 X2)`` (scalar tensor)."""
    return (u1.reshape(-1) @ u2.reshape(-1)) + 0.5 * torch_sym_inner(X1, X2)
