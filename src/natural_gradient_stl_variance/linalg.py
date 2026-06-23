"""A small batched linear-algebra backend (NumPy or a single CUDA torch device).

The estimator-variance and algorithm grids run the *same* formulas on either
backend. All matrix operations are batched over the leading axes (the last two
axes are the ``d x d`` matrix); samples carry an extra leading axis. Everything is
float64 by default.

The backend is deliberately thin: it wraps the handful of operations the
experiment needs (standard-normal draws, symmetric square root / matrix
exponential via eigendecomposition, inverse, eigenvalue extremes, reductions)
and hides the few NumPy/torch API differences (transpose of the last two axes,
``diag_embed``, generator-based RNG).

torch is imported lazily through :mod:`src.common.torch_utils`, so importing this
module without torch installed is fine until a torch backend is requested.
"""
from __future__ import annotations

import numpy as np

from src.common.torch_utils import get_torch, resolve_device, resolve_dtype


class ArrayBackend:
    """Batched float64 linear algebra on NumPy or a single torch device."""

    def __init__(self, backend="numpy", device="cpu", dtype="float64"):
        self.backend = str(backend).lower()
        if self.backend not in ("numpy", "torch"):
            raise ValueError(f"unknown backend {backend!r} (expected numpy/torch)")
        self.is_torch = self.backend == "torch"
        self.dtype_str = str(dtype)
        if self.is_torch:
            self.torch = get_torch()
            self.device = resolve_device(device)
            self.dtype = resolve_dtype(dtype)
        else:
            self.torch = None
            self.device = "cpu"
            self.dtype = np.float64 if "64" in self.dtype_str else np.float32

    # -- info ----------------------------------------------------------------
    def device_str(self):
        return str(self.device)

    # -- creation / conversion -----------------------------------------------
    def asarray(self, x):
        if self.is_torch:
            return self.torch.as_tensor(np.asarray(x, dtype=np.float64),
                                        dtype=self.dtype, device=self.device)
        return np.asarray(x, dtype=self.dtype)

    def to_numpy(self, x):
        if self.is_torch:
            return x.detach().to("cpu", dtype=self.torch.float64).numpy()
        return np.asarray(x, dtype=np.float64)

    def eye(self, d, batch=()):
        if self.is_torch:
            E = self.torch.eye(d, dtype=self.dtype, device=self.device)
        else:
            E = np.eye(d, dtype=self.dtype)
        if batch:
            E = self.broadcast_eye(E, batch)
        return E

    def broadcast_eye(self, E, batch):
        shape = tuple(batch) + E.shape
        if self.is_torch:
            return E.expand(shape).clone()
        return np.broadcast_to(E, shape).copy()

    def diag_embed(self, v):
        """Build a batched diagonal matrix from a vector ``(..., d) -> (..., d, d)``."""
        if self.is_torch:
            return self.torch.diag_embed(v)
        d = v.shape[-1]
        out = np.zeros(v.shape + (d,), dtype=self.dtype)
        idx = np.arange(d)
        out[..., idx, idx] = v
        return out

    # -- random --------------------------------------------------------------
    def generator(self, seed):
        if self.is_torch:
            g = self.torch.Generator(device=self.device)
            g.manual_seed(int(seed))
            return g
        return np.random.default_rng(int(seed))

    def randn(self, shape, generator):
        if self.is_torch:
            return self.torch.randn(*shape, generator=generator,
                                    dtype=self.dtype, device=self.device)
        return generator.standard_normal(size=shape).astype(self.dtype)

    # -- batched matrix ops (last two axes are the matrix) -------------------
    def transpose(self, A):
        if self.is_torch:
            return A.transpose(-1, -2)
        return np.swapaxes(A, -1, -2)

    def sym(self, A):
        return 0.5 * (A + self.transpose(A))

    def matmul(self, A, B):
        return A @ B

    def matvec(self, A, v):
        """Batched matrix-vector ``(..., d, d) x (..., d) -> (..., d)``."""
        if self.is_torch:
            return (A @ v.unsqueeze(-1)).squeeze(-1)
        return np.squeeze(A @ v[..., None], axis=-1)

    def eigh(self, A):
        if self.is_torch:
            return self.torch.linalg.eigh(self.sym(A))
        return np.linalg.eigh(self.sym(A))

    def eigvalsh(self, A):
        if self.is_torch:
            return self.torch.linalg.eigvalsh(self.sym(A))
        return np.linalg.eigvalsh(self.sym(A))

    def _from_eigh(self, w, V, fw):
        """Reconstruct ``V diag(fw(w)) V^T`` for a batched eigendecomposition."""
        scaled = V * fw[..., None, :]
        return self.sym(scaled @ self.transpose(V))

    def sqrtm_sym(self, A):
        w, V = self.eigh(A)
        w = self.clip_min(w, 0.0)
        return self._from_eigh(w, V, self.sqrt(w))

    def expm_sym(self, A):
        w, V = self.eigh(A)
        return self._from_eigh(w, V, self.exp(w))

    def inv(self, A):
        if self.is_torch:
            return self.torch.linalg.inv(self.sym(A))
        return np.linalg.inv(self.sym(A))

    def slogdet(self, A):
        if self.is_torch:
            sign, logabs = self.torch.linalg.slogdet(A)
            return sign, logabs
        sign, logabs = np.linalg.slogdet(A)
        return sign, logabs

    def trace(self, A):
        if self.is_torch:
            return A.diagonal(dim1=-2, dim2=-1).sum(-1)
        return np.trace(A, axis1=-2, axis2=-1)

    def diagonal(self, A):
        if self.is_torch:
            return A.diagonal(dim1=-2, dim2=-1)
        return np.diagonal(A, axis1=-2, axis2=-1)

    # -- elementwise / reductions --------------------------------------------
    def exp(self, x):
        return self.torch.exp(x) if self.is_torch else np.exp(x)

    def sqrt(self, x):
        return self.torch.sqrt(x) if self.is_torch else np.sqrt(x)

    def tanh(self, x):
        return self.torch.tanh(x) if self.is_torch else np.tanh(x)

    def abs(self, x):
        return self.torch.abs(x) if self.is_torch else np.abs(x)

    def log(self, x):
        return self.torch.log(x) if self.is_torch else np.log(x)

    def log1p(self, x):
        return self.torch.log1p(x) if self.is_torch else np.log1p(x)

    def clip_min(self, x, lo):
        if self.is_torch:
            return self.torch.clamp(x, min=lo)
        return np.maximum(x, lo)

    def mean(self, x, axis):
        if self.is_torch:
            return x.mean(dim=axis)
        return np.mean(x, axis=axis)

    def sum(self, x, axis=None):
        if self.is_torch:
            return x.sum() if axis is None else x.sum(dim=axis)
        return np.sum(x, axis=axis)

    def isfinite_all(self, x):
        if self.is_torch:
            return bool(self.torch.isfinite(x).all().item())
        return bool(np.all(np.isfinite(x)))

    def norm_fro(self, A, axis=(-2, -1)):
        if self.is_torch:
            return self.torch.linalg.matrix_norm(A, ord="fro")
        return np.linalg.norm(A, ord="fro", axis=axis)

    def norm_vec(self, v, axis=-1):
        if self.is_torch:
            return self.torch.linalg.vector_norm(v, dim=axis)
        return np.linalg.norm(v, axis=axis)
