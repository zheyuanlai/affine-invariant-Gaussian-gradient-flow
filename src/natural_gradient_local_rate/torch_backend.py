"""PyTorch GPU backend for the natural-gradient local-rate operators.

This backend reproduces the *corrected* CPU estimators (symmetrized ``H_sym``,
Fisher--Rao-scaled ``L_star``, diagonal benchmark) on a torch device, for
production-scale sample-size / operator / local-rate sweeps on Colab A100. It
does **not** replace or weaken the NumPy/SciPy CPU path; it is selected only when
``operator.backend`` resolves to ``"torch"``.

Design
------
* The torch potential is built by *copying parameters* (``rho``, ``M``, the raw
  feature ``W/c/coeff`` or separable offsets ``c``) from an already-constructed
  CPU potential, so a torch run and a CPU run on the same bank agree to floating
  point. Only ``Hess V(Z)`` is needed by the operators, so only ``batch_hess`` is
  implemented in torch.
* The dense ``H_sym`` and ``L_star`` matrices are built as **chunked matmuls**
  rather than per-basis-vector loops. In ``sym_to_vec`` (Frobenius-isometric)
  coordinates with ``Vw[j]=vec(W_j)`` and ``Vq[j]=vec(Z_j Z_j^T - I)``::

      G      = (1/M) Vw^T @ Vq            (p x p)   forward H_lin matrix
      H_sym  = 0.5 (G + G^T)
      T_mat  = (1/M) Z^T @ Vw             (N x p)   T matrix
      A      = (1/M) diag(W)^T @ (Z^2-1)  (N x N)   diagonal benchmark

  and the Fisher--Rao-packed ``L_star`` (D = N + p) is the symmetric block matrix
  ``[[I_N, T_mat/sqrt2], [T_mat^T/sqrt2, I_p + 0.5 H_sym]]``.
* Eigenvalues come from ``torch.linalg.eigh`` / ``eigvalsh`` entirely on device;
  no SciPy, no CPU/GPU transfers inside the solve.
"""
from __future__ import annotations

import math
import time

import numpy as np

from src.common.torch_utils import (
    get_torch, resolve_device, resolve_dtype, torch_device_info,
)
from src.common.torch_symspace import (
    torch_sym_dim, torch_sym_to_vec, torch_sym_to_vec_batch, torch_vec_to_sym,
    torch_unpack_tangent_fr,
)
from src.natural_gradient_local_rate import diagnostics
from src.natural_gradient_local_rate.separable_exact import separable_exact_lambda

_SQRT2 = math.sqrt(2.0)

# Potential families the torch backend can reconstruct analytically.
TORCH_SUPPORTED_FAMILIES = (
    "gaussian", "separable", "additive_index", "ridge_sum",
    "random_feature", "radial_tail",
)
# Raw-forward eigenvalue (non-symmetric) is a CPU-style diagnostic; on GPU we
# only compute it for small operators to keep production sweeps fast.
TORCH_RAW_EIG_MAX_N = 32


def torch_supports_family(family):
    return str(family) in TORCH_SUPPORTED_FAMILIES


# ---------------------------------------------------------------------------
# Sample bank
# ---------------------------------------------------------------------------

def torch_gaussian_samples(N_theta, M, seed, antithetic=True, device="cpu",
                           dtype=None):
    """Reproducible standard-Gaussian bank ``(M, N_theta)`` on ``device``.

    Drawn with a CPU ``torch.Generator`` (deterministic given ``seed``) then moved
    to ``device``. With ``antithetic=True`` the bank is ``[base, -base]`` from
    ``M // 2`` base draws (one extra independent row if ``M`` is odd).
    """
    torch = get_torch()
    if M <= 0:
        raise ValueError(f"M must be positive, got {M}")
    dtype = dtype if dtype is not None else torch.float64
    gen = torch.Generator(device="cpu")
    gen.manual_seed(int(seed))
    if antithetic:
        half = M // 2
        base = torch.randn((half, N_theta), generator=gen, dtype=dtype)
        Z = torch.cat([base, -base], dim=0)
        if 2 * half < M:
            extra = torch.randn((M - 2 * half, N_theta), generator=gen, dtype=dtype)
            Z = torch.cat([Z, extra], dim=0)
    else:
        Z = torch.randn((M, N_theta), generator=gen, dtype=dtype)
    return Z.to(device)


# ---------------------------------------------------------------------------
# Torch potentials (analytic Hess V), reconstructed from a CPU potential
# ---------------------------------------------------------------------------

def _phi_first_derivative(name):
    """Return ``phi'`` for the named nonlinearity (only ``log_cosh`` supported)."""
    torch = get_torch()
    if name == "log_cosh":
        return lambda s: torch.tanh(s)              # d/ds log cosh(s)
    raise NotImplementedError(
        f"torch backend only supports the 'log_cosh' nonlinearity, got {name!r}")


def _phi_second_derivative(name):
    """Return ``phi''`` for the named nonlinearity (only ``log_cosh`` supported)."""
    torch = get_torch()
    if name == "log_cosh":
        return lambda s: 1.0 - torch.tanh(s) ** 2   # sech^2
    raise NotImplementedError(
        f"torch backend only supports the 'log_cosh' nonlinearity, got {name!r}")


class _TorchPotential:
    """Common interface: ``N_theta`` and ``batch_hess(theta_batch) -> (B, N, N)``."""

    N_theta: int

    def batch_hess(self, theta_batch):  # pragma: no cover - abstract
        raise NotImplementedError


class TorchGaussianPotential(_TorchPotential):
    """``V = 0.5||theta||^2`` -> ``Hess V = I`` (independent of theta)."""

    def __init__(self, N_theta, device, dtype):
        self.N_theta = int(N_theta)
        self.device = device
        self.dtype = dtype
        torch = get_torch()
        self._I = torch.eye(self.N_theta, device=device, dtype=dtype)

    def batch_hess(self, theta_batch):
        B = theta_batch.shape[0]
        return self._I.expand(B, self.N_theta, self.N_theta)


class TorchSeparablePotential(_TorchPotential):
    """Centered separable log-cosh: ``Hess V`` diagonal, entry depends on theta_i."""

    def __init__(self, N_theta, rho, M_diag, c, phi_name, device, dtype):
        torch = get_torch()
        self.N_theta = int(N_theta)
        self.device = device
        self.dtype = dtype
        self.rho = float(rho)
        self.M_diag = _to_t(M_diag, device, dtype).reshape(-1)
        self.c = _to_t(c, device, dtype).reshape(-1)
        self._fpp = _phi_second_derivative(phi_name)

    def batch_hess(self, theta_batch):
        torch = get_torch()
        T = theta_batch
        d = self._fpp(T + self.c)                      # (B, N) phi''(theta_i + c_i)
        diagV = 1.0 + self.rho * (d - self.M_diag)     # (B, N)
        return torch.diag_embed(diagV)                 # (B, N, N)


class TorchRidgeSumPotential(_TorchPotential):
    """Centered ridge-sum (additive-index / random-feature) log-cosh potential."""

    def __init__(self, N_theta, rho, M, W, c, coeff, phi_name, device, dtype):
        torch = get_torch()
        self.N_theta = int(N_theta)
        self.device = device
        self.dtype = dtype
        self.rho = float(rho)
        self.M = _to_t(M, device, dtype)                       # (N, N)
        self.W = _to_t(W, device, dtype)                       # (r, N)
        self.c = _to_t(c, device, dtype).reshape(-1)           # (r,)
        self.coeff = _to_t(coeff, device, dtype).reshape(-1)   # (r,)
        self._fpp = _phi_second_derivative(phi_name)
        self._I = torch.eye(self.N_theta, device=device, dtype=dtype)

    def batch_hess(self, theta_batch):
        torch = get_torch()
        T = theta_batch
        s = T @ self.W.T + self.c                       # (B, r)
        d = self._fpp(s) * self.coeff                   # (B, r)
        hphi = torch.einsum("br,ri,rj->bij", d, self.W, self.W)   # (B, N, N)
        H = self._I + self.rho * (hphi - self.M)
        return 0.5 * (H + H.transpose(-1, -2))


class TorchRadialTailPotential(_TorchPotential):
    """Centered radial-tail feature ``Phi(theta) = h(||theta||^2)`` (log-cosh ``h``).

    With ``R = ||theta||^2`` and ``u = scale * (R / N - shift)``::

        h'(R)    = (scale / N) * phi'(u)
        h''(R)   = (scale / N)^2 * phi''(u)
        Hess Phi = 2 h'(R) I + 4 h''(R) theta theta^T

    and the centered Hessian is ``Hess V = I + rho (Hess Phi - M)`` (matching the
    NumPy :class:`RadialTailPotential` / :class:`CenteredPotential` algebra). The
    full mean Hessian ``M`` is copied from the CPU potential so the GPU and CPU
    operators agree to floating point on the same bank.
    """

    def __init__(self, N_theta, rho, M, scale, shift, phi_name, device, dtype):
        torch = get_torch()
        self.N_theta = int(N_theta)
        self.device = device
        self.dtype = dtype
        self.rho = float(rho)
        self.scale = float(scale)
        self.shift = float(shift)
        self.M = _to_t(M, device, dtype)                  # (N, N)
        self._fp = _phi_first_derivative(phi_name)        # phi'  = tanh
        self._fpp = _phi_second_derivative(phi_name)      # phi'' = sech^2
        self._I = torch.eye(self.N_theta, device=device, dtype=dtype)

    def batch_hess(self, theta_batch):
        torch = get_torch()
        T = theta_batch                                     # (B, N)
        R = torch.sum(T * T, dim=1)                         # (B,)
        u = self.scale * (R / self.N_theta - self.shift)    # (B,)
        hp = (self.scale / self.N_theta) * self._fp(u)            # (B,) h'(R)
        hpp = (self.scale / self.N_theta) ** 2 * self._fpp(u)    # (B,) h''(R)
        outer = torch.einsum("bi,bj->bij", T, T)                 # (B, N, N)
        hphi = (2.0 * hp[:, None, None] * self._I
                + 4.0 * hpp[:, None, None] * outer)              # (B, N, N)
        H = self._I + self.rho * (hphi - self.M)
        return 0.5 * (H + H.transpose(-1, -2))


def torch_potential_from_cpu(cpu_pot, device, dtype):
    """Build a torch potential by copying parameters from a CPU potential object."""
    # Imported here to avoid any import-time coupling; these are NumPy-only.
    from src.natural_gradient_local_rate.potentials import (
        GaussianPotential, CenteredPotential, SeparablePotential,
    )
    from src.natural_gradient_local_rate.potentials.additive_index import RidgeSumFeature
    from src.natural_gradient_local_rate.potentials.radial_tail import RadialTailPotential

    if isinstance(cpu_pot, GaussianPotential):
        return TorchGaussianPotential(cpu_pot.N_theta, device, dtype)
    if isinstance(cpu_pot, CenteredPotential):
        raw = cpu_pot.raw
        if isinstance(raw, SeparablePotential):
            return TorchSeparablePotential(
                cpu_pot.N_theta, cpu_pot.rho, np.diag(cpu_pot.M), raw.c,
                raw.phi_name, device, dtype)
        if isinstance(raw, RidgeSumFeature):
            return TorchRidgeSumPotential(
                cpu_pot.N_theta, cpu_pot.rho, cpu_pot.M, raw.W, raw.c, raw.coeff,
                raw.phi_name, device, dtype)
        if isinstance(raw, RadialTailPotential):
            return TorchRadialTailPotential(
                cpu_pot.N_theta, cpu_pot.rho, cpu_pot.M, raw.scale, raw.shift,
                raw.phi_name, device, dtype)
        raise NotImplementedError(
            f"torch backend does not support raw feature {type(raw).__name__!r} "
            "(supported: SeparablePotential, RidgeSumFeature, RadialTailPotential). "
            "Use backend='numpy'.")
    raise NotImplementedError(
        f"torch backend does not support potential {type(cpu_pot).__name__!r}.")


# ---------------------------------------------------------------------------
# Chunking
# ---------------------------------------------------------------------------

def _to_t(x, device, dtype):
    """``torch.as_tensor`` on a writable contiguous copy (avoids read-only views)."""
    torch = get_torch()
    return torch.as_tensor(np.ascontiguousarray(np.asarray(x, dtype=np.float64)),
                           device=device, dtype=dtype)


def _chunk_bounds(M, chunk_size):
    cs = int(chunk_size) if chunk_size else M
    cs = max(1, cs)
    return [(s, min(s + cs, M)) for s in range(0, M, cs)]


def _as_device_tensor(Z, device, dtype):
    torch = get_torch()
    if isinstance(Z, np.ndarray):
        return torch.as_tensor(Z, device=device, dtype=dtype)
    return Z.to(device=device, dtype=dtype)


# ---------------------------------------------------------------------------
# Matrix-free operator applications (mirror the NumPy operators)
# ---------------------------------------------------------------------------

def torch_apply_H_forward(pot_t, Z, X, chunk_size=None):
    """``H_forward[X] = mean_j (Z_j^T X Z_j - Tr X) W_j`` (symmetric matrix)."""
    torch = get_torch()
    N, M = Z.shape[1], Z.shape[0]
    Xs = 0.5 * (X + X.T)
    trX = torch.trace(Xs)
    acc = torch.zeros((N, N), device=Z.device, dtype=Z.dtype)
    for s, e in _chunk_bounds(M, chunk_size):
        Zc = Z[s:e]
        W = pot_t.batch_hess(Zc)
        q = torch.einsum("bi,ij,bj->b", Zc, Xs, Zc) - trX
        acc += torch.einsum("b,bij->ij", q, W)
    acc /= M
    return 0.5 * (acc + acc.T)


def torch_apply_H_adjoint(pot_t, Z, Y, chunk_size=None):
    """``H_adjoint[Y] = mean_j Tr(Y W_j) (Z_j Z_j^T - I)`` (symmetric matrix)."""
    torch = get_torch()
    N, M = Z.shape[1], Z.shape[0]
    Ys = 0.5 * (Y + Y.T)
    I = torch.eye(N, device=Z.device, dtype=Z.dtype)
    acc = torch.zeros((N, N), device=Z.device, dtype=Z.dtype)
    for s, e in _chunk_bounds(M, chunk_size):
        Zc = Z[s:e]
        W = pot_t.batch_hess(Zc)
        d = torch.einsum("ij,bij->b", Ys, W)            # Tr(Y W_j)
        zz = torch.einsum("b,bi,bj->ij", d, Zc, Zc)
        acc += zz - d.sum() * I
    acc /= M
    return 0.5 * (acc + acc.T)


def torch_apply_H_sym(pot_t, Z, X, chunk_size=None):
    """Self-adjoint estimator ``0.5 (H_forward + H_adjoint)[X]``."""
    return 0.5 * (torch_apply_H_forward(pot_t, Z, X, chunk_size)
                  + torch_apply_H_adjoint(pot_t, Z, X, chunk_size))


def torch_apply_T(pot_t, Z, X, chunk_size=None):
    """``T[X] = mean_j Tr(X W_j) Z_j`` -> vector ``(N,)``."""
    torch = get_torch()
    N, M = Z.shape[1], Z.shape[0]
    Xs = 0.5 * (X + X.T)
    acc = torch.zeros(N, device=Z.device, dtype=Z.dtype)
    for s, e in _chunk_bounds(M, chunk_size):
        Zc = Z[s:e]
        W = pot_t.batch_hess(Zc)
        d = torch.einsum("ij,bij->b", Xs, W)            # Tr(X W_j)
        acc += Zc.T @ d
    return acc / M


def torch_apply_T_star(pot_t, Z, u, chunk_size=None):
    """``T_star[u] = mean_j (Z_j . u) W_j`` -> symmetric matrix ``(N, N)``."""
    torch = get_torch()
    N, M = Z.shape[1], Z.shape[0]
    u = u.reshape(-1)
    acc = torch.zeros((N, N), device=Z.device, dtype=Z.dtype)
    for s, e in _chunk_bounds(M, chunk_size):
        Zc = Z[s:e]
        W = pot_t.batch_hess(Zc)
        e_ = Zc @ u
        acc += torch.einsum("b,bij->ij", e_, W)
    acc /= M
    return 0.5 * (acc + acc.T)


# ---------------------------------------------------------------------------
# Dense operator construction (chunked matmuls)
# ---------------------------------------------------------------------------

def _accumulate_dense(pot_t, Z, chunk_size=None):
    """Single chunked pass returning ``(G, T_mat, A)`` on device.

    ``G = (1/M) Vw^T Vq`` (p x p forward H_lin matrix), ``T_mat = (1/M) Z^T Vw``
    (N x p), ``A = (1/M) diag(W)^T (Z^2 - 1)`` (N x N diagonal benchmark).
    """
    torch = get_torch()
    N, M = Z.shape[1], Z.shape[0]
    p = torch_sym_dim(N)
    I = torch.eye(N, device=Z.device, dtype=Z.dtype)
    G = torch.zeros((p, p), device=Z.device, dtype=Z.dtype)
    T_mat = torch.zeros((N, p), device=Z.device, dtype=Z.dtype)
    A = torch.zeros((N, N), device=Z.device, dtype=Z.dtype)
    for s, e in _chunk_bounds(M, chunk_size):
        Zc = Z[s:e]
        W = pot_t.batch_hess(Zc)                          # (B, N, N)
        Vw = torch_sym_to_vec_batch(W)                    # (B, p)
        ZZ = torch.einsum("bi,bj->bij", Zc, Zc) - I       # (B, N, N)
        Vq = torch_sym_to_vec_batch(ZZ)                   # (B, p)
        G += Vw.T @ Vq
        T_mat += Zc.T @ Vw
        diagH = torch.diagonal(W, dim1=-2, dim2=-1)       # (B, N)
        A += diagH.T @ (Zc * Zc - 1.0)
    return G / M, T_mat / M, A / M


def torch_dense_H_sym_matrix(pot_t, Z, chunk_size=None, basis_block_size=None):
    """Dense symmetric ``H_sym`` matrix ``(p, p)`` on device.

    ``basis_block_size`` is accepted for API compatibility; the matmul builder is
    vectorized over the sample bank (via ``chunk_size``) and needs no basis
    blocking.
    """
    G, _, _ = _accumulate_dense(pot_t, Z, chunk_size)
    return 0.5 * (G + G.T)


def torch_dense_L_star_matrix(pot_t, Z, chunk_size=None, basis_block_size=None):
    """Dense symmetric ``L_star`` matrix ``(D, D)`` in Fisher--Rao packed coords."""
    torch = get_torch()
    G, T_mat, _ = _accumulate_dense(pot_t, Z, chunk_size)
    N = Z.shape[1]
    p = G.shape[0]
    H_sym = 0.5 * (G + G.T)
    D = N + p
    L = torch.zeros((D, D), device=Z.device, dtype=Z.dtype)
    L[:N, :N] = torch.eye(N, device=Z.device, dtype=Z.dtype)
    L[:N, N:] = T_mat / _SQRT2
    L[N:, :N] = T_mat.T / _SQRT2
    L[N:, N:] = torch.eye(p, device=Z.device, dtype=Z.dtype) + 0.5 * H_sym
    return 0.5 * (L + L.T)


def torch_diagonal_A_matrix(pot_t, Z, chunk_size=None):
    """Diagonal-mode coefficient matrix ``A = G - 1 1^T`` (shape ``(N, N)``)."""
    _, _, A = _accumulate_dense(pot_t, Z, chunk_size)
    return A


def torch_diagonal_lambda(pot_t, Z, chunk_size=None):
    """Diagonal benchmark: ``(Lambda_hat_diag, diag_offdiag_norm)`` from ``A``.

    ``Lambda_hat_diag`` is ``lambda_max`` of the *symmetric part* of ``A`` (the
    manuscript diagonal-mode quadratic form). The cleanest separable estimator is
    ``max_i A_ii`` -- use :func:`torch_diagonal_A_matrix` for that.
    """
    torch = get_torch()
    A = torch_diagonal_A_matrix(pot_t, Z, chunk_size)
    A_sym = 0.5 * (A + A.T)
    lam = torch.linalg.eigvalsh(A_sym)[-1]
    offdiag = A - torch.diag(torch.diagonal(A))
    return float(lam.item()), float(torch.linalg.norm(offdiag).item())


# ---------------------------------------------------------------------------
# Self-adjointness probes on dense matrices
# ---------------------------------------------------------------------------

def _self_adjoint_error_matrix(Mt, n_probe=8, seed=0):
    """Relative asymmetry of a dense operator matrix via random Euclidean probes."""
    torch = get_torch()
    dim = Mt.shape[0]
    gen = torch.Generator(device="cpu")
    gen.manual_seed(int(seed))
    worst = 0.0
    for _ in range(n_probe):
        a = torch.randn(dim, generator=gen, dtype=Mt.dtype).to(Mt.device)
        b = torch.randn(dim, generator=gen, dtype=Mt.dtype).to(Mt.device)
        lhs = float((a @ (Mt @ b)).item())
        rhs = float(((Mt @ a) @ b).item())
        worst = max(worst, abs(lhs - rhs) / max(1.0, abs(lhs), abs(rhs)))
    return worst


# ---------------------------------------------------------------------------
# Dense GPU eigen-estimators
# ---------------------------------------------------------------------------

def _residual(Mt, w, v):
    torch = get_torch()
    return float((torch.linalg.norm(Mt @ v - w * v) / max(1.0, abs(float(w)))).item())


def torch_estimate_lambda_hat_dense(pot_t, Z, chunk_size=None, basis_block_size=None):
    """Largest eigenvalue of dense ``H_sym`` with its eigen-residual."""
    torch = get_torch()
    H = torch_dense_H_sym_matrix(pot_t, Z, chunk_size, basis_block_size)
    w, V = torch.linalg.eigh(H)
    lam = w[-1]
    return float(lam.item()), _residual(H, lam, V[:, -1])


def torch_estimate_gamma_loc_dense(pot_t, Z, chunk_size=None, basis_block_size=None):
    """Smallest eigenvalue of dense ``L_star`` (gamma_loc) with eigen-residual and slow mode."""
    torch = get_torch()
    L = torch_dense_L_star_matrix(pot_t, Z, chunk_size, basis_block_size)
    w, V = torch.linalg.eigh(L)
    gamma = w[0]
    v0 = V[:, 0]
    u_star, X_star = torch_unpack_tangent_fr(v0, pot_t.N_theta)
    return (float(gamma.item()), _residual(L, gamma, v0),
            u_star.detach().cpu().numpy(), X_star.detach().cpu().numpy())


# ---------------------------------------------------------------------------
# Optional matrix-free subspace iteration (future-proofing for N > dense max)
# ---------------------------------------------------------------------------

def torch_subspace_largest_eig(apply_A, dim, block_size=8, n_iter=100, tol=1e-6,
                               device="cpu", dtype=None, seed=0):
    """Block subspace iteration for the largest eigenvalue of a symmetric operator.

    ``apply_A`` maps a ``(dim, k)`` block to a ``(dim, k)`` block. This is a simple
    (non-production) solver kept for ``N_theta`` beyond the dense limit; the dense
    path is the supported production route in this pass.
    """
    torch = get_torch()
    dtype = dtype if dtype is not None else torch.float64
    gen = torch.Generator(device="cpu")
    gen.manual_seed(int(seed))
    Q = torch.randn((dim, block_size), generator=gen, dtype=dtype).to(device)
    Q, _ = torch.linalg.qr(Q)
    prev = None
    for _ in range(int(n_iter)):
        Q, _ = torch.linalg.qr(apply_A(Q))
        T = Q.T @ apply_A(Q)
        w = torch.linalg.eigvalsh(0.5 * (T + T.T))
        top = float(w[-1].item())
        if prev is not None and abs(top - prev) <= tol * max(1.0, abs(top)):
            break
        prev = top
    return top, Q


# ---------------------------------------------------------------------------
# Full GPU estimate (one chunked pass + dense eigh)
# ---------------------------------------------------------------------------

def torch_operator_estimates(pot_t, Z, *, chunk_size=None, basis_block_size=None,
                             self_adjoint_probes=8, self_adjoint_seed=0,
                             compute_gamma_loc=True, compute_raw_forward="auto"):
    """Run the full GPU operator suite from a single dense accumulation pass.

    Returns a dict of numeric estimates plus timing / matrix-size metadata. All
    heavy linear algebra stays on ``Z.device``; only final scalars are moved off.
    """
    torch = get_torch()
    N, M = Z.shape[1], Z.shape[0]
    p = torch_sym_dim(N)
    D = N + p
    out = {
        "operator_matrix_dim": int(D),
        "dense_matrix_memory_mb": float(D * D * torch.finfo(Z.dtype).bits / 8 / 1e6),
    }

    t_build = time.time()
    G, T_mat, A = _accumulate_dense(pot_t, Z, chunk_size)
    H_sym = 0.5 * (G + G.T)
    L = torch.zeros((D, D), device=Z.device, dtype=Z.dtype)
    L[:N, :N] = torch.eye(N, device=Z.device, dtype=Z.dtype)
    L[:N, N:] = T_mat / _SQRT2
    L[N:, :N] = T_mat.T / _SQRT2
    L[N:, N:] = torch.eye(p, device=Z.device, dtype=Z.dtype) + 0.5 * H_sym
    out["matrix_construction_runtime_seconds"] = float(time.time() - t_build)

    t_eig = time.time()
    wH, VH = torch.linalg.eigh(H_sym)
    lam = wH[-1]
    out["Lambda_hat_full_sym"] = float(lam.item())
    out["eig_residual_H"] = _residual(H_sym, lam, VH[:, -1])

    if compute_gamma_loc:
        wL, VL = torch.linalg.eigh(L)
        gamma = wL[0]
        v0 = VL[:, 0]
        out["gamma_loc"] = float(gamma.item())
        out["eig_residual_L_star"] = _residual(L, gamma, v0)
        u_star, X_star = torch_unpack_tangent_fr(v0, N)
        out["_u_star"] = u_star.detach().cpu().numpy()
        out["_X_star"] = X_star.detach().cpu().numpy()
    out["eigh_runtime_seconds"] = float(time.time() - t_eig)

    # diagonal benchmark
    A_sym = 0.5 * (A + A.T)
    out["Lambda_hat_diag"] = float(torch.linalg.eigvalsh(A_sym)[-1].item())
    offdiag = A - torch.diag(torch.diagonal(A))
    out["diag_offdiag_norm"] = float(torch.linalg.norm(offdiag).item())

    # self-adjointness diagnostics (probes on dense matrices)
    out["self_adjoint_error_H_sym"] = _self_adjoint_error_matrix(
        H_sym, self_adjoint_probes, self_adjoint_seed)
    out["self_adjoint_error_H_raw"] = _self_adjoint_error_matrix(
        G, self_adjoint_probes, self_adjoint_seed)
    if compute_gamma_loc:
        out["self_adjoint_error_L_star"] = _self_adjoint_error_matrix(
            L, self_adjoint_probes, self_adjoint_seed)

    # raw-forward largest-real eigenvalue (diagnostic; small N only on GPU)
    do_raw = (N <= TORCH_RAW_EIG_MAX_N) if compute_raw_forward == "auto" else bool(compute_raw_forward)
    if do_raw:
        try:
            ev = torch.linalg.eigvals(G)
            out["Lambda_hat_raw_forward"] = float(ev.real.max().item())
        except Exception:
            out["Lambda_hat_raw_forward"] = float("nan")
    else:
        out["Lambda_hat_raw_forward"] = float("nan")

    del G, T_mat, A, H_sym, L
    return out


# ---------------------------------------------------------------------------
# Row assembly (mirrors estimator_suite.compute_row schema on the torch path)
# ---------------------------------------------------------------------------

def _safe_ratio(a, b):
    if b is None or not np.isfinite(b) or b == 0.0:
        return float("nan")
    return float(a) / float(b)


def _safe_diff(a, b):
    if a is None or b is None or not (np.isfinite(a) and np.isfinite(b)):
        return float("nan")
    return float(a) - float(b)


def compute_row_torch(cpu_potential, Z, point, opts, *, run_id="",
                      experiment_group="natural_gradient_local_rate"):
    """Torch/GPU analogue of :func:`estimator_suite.compute_row` (same schema + GPU columns).

    ``cpu_potential`` is the NumPy potential (parameters are copied to torch);
    ``Z`` is the NumPy sample bank (moved to the torch device). The separable
    exact benchmark reuses the cheap CPU quadrature on ``cpu_potential``.
    """
    torch = get_torch()
    # Resolve device/dtype *outside* the per-point try so a misconfigured CUDA
    # request fails fast and clearly rather than being swallowed per row.
    device = resolve_device(opts.get("device", "auto"))
    dtype = resolve_dtype(opts.get("dtype", "float64"))
    if not torch_supports_family(point["family"]):
        raise NotImplementedError(
            f"torch backend does not support family {point['family']!r}; "
            "use backend='numpy' (or backend='auto' to fall back per point).")

    md = cpu_potential.metadata()
    N = int(point["N_theta"])
    kappa = float(point["kappa_target"])
    beta = float(md.get("beta_target", 1.0))
    cs = opts.get("chunk_size")
    info = torch_device_info(opts.get("device", "auto"), opts.get("dtype", "float64"))

    row = {
        "run_id": run_id,
        "experiment_group": experiment_group,
        "potential_family": point["family"],
        "family": point["family"],
        "seed": int(point["seed"]),
        "N_theta": N,
        "kappa_target": kappa,
        "rho": float(md.get("rho", 0.0)),
        "M_mc": int(point["M_mc"]),
        "operator_estimator": opts.get("estimator", "symmetrized"),
        "backend": "torch",
        "chunk_size": cs if cs is not None else "",
        "quadrature_nodes": int(opts.get("quadrature_nodes", 80)),
        "device": str(device),
        "dtype": str(opts.get("dtype", "float64")),
        "eigensolver": str(opts.get("eigensolver", "torch_dense_eigh")),
        "explicit_dense_max_N_theta": int(opts.get("explicit_dense_max_N_theta", 64)),
        "basis_block_size": int(opts.get("basis_block_size", 32)),
        "torch_available": info["torch_available"],
        "torch_version": info["torch_version"],
        "cuda_available": info["cuda_available"],
        "cuda_version": info["cuda_version"],
        "device_name": info["device_name"],
        "gpu_available": info["gpu_available"],
        # numeric columns default to NaN until computed
        "Lambda_hat_raw_forward": float("nan"),
        "Lambda_hat_full_sym": float("nan"),
        "Lambda_hat_diag": float("nan"),
        "Lambda_hat_separable_exact": float("nan"),
        "diag_offdiag_norm": float("nan"),
        "full_over_diag": float("nan"),
        "full_minus_diag": float("nan"),
        "diag_minus_exact": float("nan"),
        "full_sym_minus_exact": float("nan"),
        "gamma_loc": float("nan"),
        "inverse_gamma_loc": float("nan"),
        "self_adjoint_error_H_raw": float("nan"),
        "self_adjoint_error_H_sym": float("nan"),
        "self_adjoint_error_L_star": float("nan"),
        "eig_residual_H": float("nan"),
        "eig_residual_L_star": float("nan"),
        "operator_matrix_dim": int(N + torch_sym_dim(N)),
        "dense_matrix_memory_mb": float("nan"),
        "matrix_construction_runtime_seconds": float("nan"),
        "eigh_runtime_seconds": float("nan"),
        "gpu_peak_memory_mb": float("nan"),
        "separable_exact_status": "",
        "status": "ok",
        "error_message": "",
    }
    row.update(diagnostics.true_benchmark_columns(point["family"]))

    dense_max = int(opts.get("explicit_dense_max_N_theta", 64))
    t0 = time.time()
    try:
        if N > dense_max:
            raise NotImplementedError(
                f"N_theta={N} exceeds explicit_dense_max_N_theta={dense_max}; the "
                "dense GPU eigensolver is the only supported torch path in this "
                "pass (raise explicit_dense_max_N_theta to allow, memory permitting).")
        if device.type == "cuda":
            torch.cuda.reset_peak_memory_stats(device)

        Zt = _as_device_tensor(Z, device, dtype)
        pot_t = torch_potential_from_cpu(cpu_potential, device, dtype)
        est = torch_operator_estimates(
            pot_t, Zt, chunk_size=cs,
            basis_block_size=opts.get("basis_block_size", 32),
            self_adjoint_probes=int(opts.get("self_adjoint_probes", 8)),
            self_adjoint_seed=int(opts.get("self_adjoint_seed", 0)),
            compute_gamma_loc=bool(opts.get("compute_gamma_loc", True)))
        row.update({k: v for k, v in est.items() if not k.startswith("_")})
        if "_u_star" in est:
            row["_u_star"] = est["_u_star"]
            row["_X_star"] = est["_X_star"]

        if np.isfinite(row["gamma_loc"]):
            row["inverse_gamma_loc"] = (1.0 / row["gamma_loc"]
                                        if row["gamma_loc"] != 0.0 else float("inf"))

        # separable exact via the cheap CPU quadrature on the CPU potential
        if bool(opts.get("separable_exact_benchmark", True)):
            lam_exact, exact_status = separable_exact_lambda(
                cpu_potential, n_nodes=int(opts.get("quadrature_nodes", 80)))
            row["Lambda_hat_separable_exact"] = lam_exact
            row["separable_exact_status"] = exact_status
            if row["baseline_type"] == "" and exact_status == "ok":
                row["baseline_type"] = "separable_exact"

        row["full_over_diag"] = _safe_ratio(row["Lambda_hat_full_sym"], row["Lambda_hat_diag"])
        row["full_minus_diag"] = _safe_diff(row["Lambda_hat_full_sym"], row["Lambda_hat_diag"])
        row["diag_minus_exact"] = _safe_diff(row["Lambda_hat_diag"], row["Lambda_hat_separable_exact"])
        row["full_sym_minus_exact"] = _safe_diff(row["Lambda_hat_full_sym"], row["Lambda_hat_separable_exact"])

        headline = (row["Lambda_hat_raw_forward"] if row["operator_estimator"] == "raw_forward"
                    else row["Lambda_hat_full_sym"])
        row["Lambda_hat"] = float(headline)
        ref = diagnostics.reference_columns(
            N, kappa, beta, Lambda_hat=headline,
            gamma_loc=(row["gamma_loc"] if np.isfinite(row["gamma_loc"]) else None))
        row.update(ref)
        for k in ("alpha_target", "beta_target", "L_A", "L_A_is_empirical",
                  "norm_mean_grad", "norm_mean_hess_minus_I",
                  "empirical_min_hess_eig", "empirical_max_hess_eig"):
            if k in md:
                row[k] = md[k]

        if device.type == "cuda":
            row["gpu_peak_memory_mb"] = float(torch.cuda.max_memory_allocated(device) / 1e6)
            del Zt, pot_t
            torch.cuda.empty_cache()
    except Exception as exc:
        row["status"] = "error"
        row["error_message"] = f"{type(exc).__name__}: {exc}"

    row["runtime_seconds"] = float(time.time() - t0)
    return row
