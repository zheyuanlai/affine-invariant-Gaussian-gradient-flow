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
    torch_unpack_tangent_fr, torch_outer_minus_I_symvec, _triu,
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
        return torch.diag_embed(self.batch_hess_diag(theta_batch))     # (B, N, N)

    def batch_hess_diag(self, theta_batch):
        """Diagonal entries of ``Hess V`` only, shape ``(B, N)`` (no dense matrix)."""
        d = self._fpp(theta_batch + self.c)            # (B, N) phi''(theta_i + c_i)
        return 1.0 + self.rho * (d - self.M_diag)      # (B, N)


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
# GPU-side centering statistics (the expensive bank reductions, on device)
# ---------------------------------------------------------------------------
#
# Profiling the low-dimensional high-M H200 path showed >90% of wall time is the
# *CPU* CenteredPotential construction: ``phi_grad``/``phi_hess`` are evaluated
# over the full multi-million-sample bank (twice -- once for E[Hess Phi], once
# for the eigenvalue extremes of Hess Phi - M), entirely on CPU, before the GPU
# operator pass even starts. These helpers move exactly those three reductions
#
#     b      = E[grad Phi(Z)]                              (N,)
#     M      = E[Hess Phi(Z)]            (symmetrized)     (N, N)
#     emin,  = min/max eigenvalue of Hess Phi(Z_j) - M over the bank
#     emax
#
# onto the device, reproducing the NumPy raw-feature math to float64. The cheap
# centering *algebra* (L_A, rho, alpha/beta_target, diagnostics) is unchanged and
# stays in CenteredPotential via its ``precomputed_stats`` path, so the resulting
# potential is byte-for-byte the CPU one. The numpy bank itself is reused as-is
# (its RNG stream is the experiment design; see notes in run_operator_*).

# cusolver's batched ``syevj`` (eigvalsh) overflows internal workspace above
# ~32k matrices in one call, so the per-chunk eigenvalue pass is sub-batched.
_EIG_SUBBATCH = 32768


def _batched_eig_extremes(Hb, emin, emax):
    """Update running (emin, emax) with eigenvalue extremes of a batch ``(B,N,N)``.

    ``eigvalsh`` is called in sub-batches of ``_EIG_SUBBATCH`` to stay within the
    cusolver batched-solver workspace limit.
    """
    torch = get_torch()
    Hs = 0.5 * (Hb + Hb.transpose(-1, -2))
    B = Hs.shape[0]
    for s in range(0, B, _EIG_SUBBATCH):
        w = torch.linalg.eigvalsh(Hs[s:s + _EIG_SUBBATCH])
        emin = min(emin, float(w[:, 0].min().item()))
        emax = max(emax, float(w[:, -1].max().item()))
    return emin, emax


class _TorchRawFeature:
    """On-device evaluators for a raw feature's ``mean grad Phi`` / ``Hess Phi``.

    Subclasses provide ``grad_phi_sum(Zc)`` (-> running-summable ``(N,)``) and
    ``hess_phi_batch(Zc)`` (-> ``(B, N, N)`` raw Hessians), matching the NumPy
    :class:`RawFeaturePotential` math exactly. The generic
    :func:`_centering_stats` driver does the chunked accumulation and eig pass.
    """

    N_theta: int

    def grad_phi_sum(self, Zc):  # pragma: no cover - abstract
        raise NotImplementedError

    def hess_phi_batch(self, Zc):  # pragma: no cover - abstract
        raise NotImplementedError


class _TorchRidgeRaw(_TorchRawFeature):
    """Raw ridge-sum feature: ``Hess Phi = W^T diag(coeff phi''(W theta + c)) W``."""

    def __init__(self, N_theta, W, c, coeff, phi_name, device, dtype):
        self.N_theta = int(N_theta)
        self.W = _to_t(W, device, dtype)                    # (r, N)
        self.c = _to_t(c, device, dtype).reshape(-1)        # (r,)
        self.coeff = _to_t(coeff, device, dtype).reshape(-1)  # (r,)
        self._fp = _phi_first_derivative(phi_name)
        self._fpp = _phi_second_derivative(phi_name)

    def grad_phi_sum(self, Zc):
        # grad Phi = (phi'(s) * coeff) @ W ; sum over the chunk's rows.
        s = Zc @ self.W.T + self.c                          # (B, r)
        g = (self._fp(s) * self.coeff) @ self.W             # (B, N)
        return g.sum(dim=0)

    def hess_phi_batch(self, Zc):
        torch = get_torch()
        s = Zc @ self.W.T + self.c                          # (B, r)
        d = self._fpp(s) * self.coeff                       # (B, r)
        return torch.einsum("br,ri,rj->bij", d, self.W, self.W)  # (B, N, N)


class _TorchSeparableRaw(_TorchRawFeature):
    """Raw separable feature: ``Hess Phi`` diagonal, ``grad Phi = phi'(theta + c)``."""

    def __init__(self, N_theta, c, phi_name, device, dtype):
        self.N_theta = int(N_theta)
        self.c = _to_t(c, device, dtype).reshape(-1)        # (N,)
        self._fp = _phi_first_derivative(phi_name)
        self._fpp = _phi_second_derivative(phi_name)

    def grad_phi_sum(self, Zc):
        return self._fp(Zc + self.c).sum(dim=0)             # (N,)

    def hess_phi_batch(self, Zc):
        torch = get_torch()
        d = self._fpp(Zc + self.c)                          # (B, N)
        return torch.diag_embed(d)                          # (B, N, N)


class _TorchRadialRaw(_TorchRawFeature):
    """Raw radial-tail feature: ``Hess Phi = 2 h'(R) I + 4 h''(R) theta theta^T``."""

    def __init__(self, N_theta, scale, shift, phi_name, device, dtype):
        self.N_theta = int(N_theta)
        self.scale = float(scale)
        self.shift = float(shift)
        self._fp = _phi_first_derivative(phi_name)
        self._fpp = _phi_second_derivative(phi_name)
        self._I = get_torch().eye(self.N_theta, device=device, dtype=dtype)

    def _hp_hpp(self, Zc):
        torch = get_torch()
        R = torch.sum(Zc * Zc, dim=1)                       # (B,)
        u = self.scale * (R / self.N_theta - self.shift)    # (B,)
        hp = (self.scale / self.N_theta) * self._fp(u)      # h'(R), (B,)
        hpp = (self.scale / self.N_theta) ** 2 * self._fpp(u)  # h''(R), (B,)
        return hp, hpp

    def grad_phi_sum(self, Zc):
        hp, _ = self._hp_hpp(Zc)
        return (2.0 * hp[:, None] * Zc).sum(dim=0)          # (N,)

    def hess_phi_batch(self, Zc):
        torch = get_torch()
        hp, hpp = self._hp_hpp(Zc)
        outer = torch.einsum("bi,bj->bij", Zc, Zc)          # (B, N, N)
        return (2.0 * hp[:, None, None] * self._I
                + 4.0 * hpp[:, None, None] * outer)         # (B, N, N)


def _torch_raw_feature_from_cpu(raw, device, dtype):
    """Build the on-device raw-feature evaluator from a NumPy raw feature."""
    from src.natural_gradient_local_rate.potentials import SeparablePotential
    from src.natural_gradient_local_rate.potentials.additive_index import RidgeSumFeature
    from src.natural_gradient_local_rate.potentials.radial_tail import RadialTailPotential
    if isinstance(raw, SeparablePotential):
        return _TorchSeparableRaw(raw.N_theta, raw.c, raw.phi_name, device, dtype)
    if isinstance(raw, RidgeSumFeature):
        return _TorchRidgeRaw(raw.N_theta, raw.W, raw.c, raw.coeff, raw.phi_name,
                              device, dtype)
    if isinstance(raw, RadialTailPotential):
        return _TorchRadialRaw(raw.N_theta, raw.scale, raw.shift, raw.phi_name,
                               device, dtype)
    raise NotImplementedError(
        f"GPU centering does not support raw feature {type(raw).__name__!r}.")


def centering_stats_on_device(raw, Z, device, dtype, chunk_size=None):
    """Compute ``(b, M, emin_dev, emax_dev, mean_Z)`` for a raw feature on device.

    ``Z`` is the NumPy sample bank (moved to ``device`` in chunks). Reproduces the
    NumPy :class:`CenteredPotential` reductions exactly:

    * ``b      = mean_j grad Phi(Z_j)``               (N,)
    * ``M      = sym(mean_j Hess Phi(Z_j))``          (N, N)
    * ``emin/emax`` = global min/max eigenvalue of ``Hess Phi(Z_j) - M``
    * ``mean_Z = mean_j Z_j``                         (N,)

    Two chunked passes (mean first, then eig-extremes against the finished ``M``)
    mirror the CPU code's two passes. Returns NumPy arrays / floats ready for
    :class:`CenteredPotential`'s ``precomputed_stats``.
    """
    torch = get_torch()
    rawt = _torch_raw_feature_from_cpu(raw, device, dtype)
    N, M_tot = Z.shape[1], Z.shape[0]
    b_acc = torch.zeros(N, device=device, dtype=dtype)
    M_acc = torch.zeros((N, N), device=device, dtype=dtype)
    z_acc = torch.zeros(N, device=device, dtype=dtype)
    for s, e in _chunk_bounds(M_tot, chunk_size):
        Zc = _as_device_tensor(Z[s:e], device, dtype)
        b_acc += rawt.grad_phi_sum(Zc)
        M_acc += rawt.hess_phi_batch(Zc).sum(dim=0)
        z_acc += Zc.sum(dim=0)
    b = b_acc / M_tot
    Mmat = M_acc / M_tot
    Mmat = 0.5 * (Mmat + Mmat.T)
    mean_Z = z_acc / M_tot

    emin, emax = float("inf"), float("-inf")
    for s, e in _chunk_bounds(M_tot, chunk_size):
        Zc = _as_device_tensor(Z[s:e], device, dtype)
        Hdev = rawt.hess_phi_batch(Zc) - Mmat
        emin, emax = _batched_eig_extremes(Hdev, emin, emax)

    return (b.detach().cpu().numpy(), Mmat.detach().cpu().numpy(),
            float(emin), float(emax), mean_Z.detach().cpu().numpy())


def build_centered_potential_gpu(family, point, cfg_potential, Z, device, dtype,
                                 *, chunk_size=None, safety_factor=2.0,
                                 phi="log_cosh", feature_multiplier=4):
    """Build a CPU :class:`CenteredPotential` whose bank reductions ran on device.

    Constructs the raw feature (single source of truth, shared with the CPU
    builder), computes the heavy ``(b, M, eig-extremes)`` reductions on ``device``
    via :func:`centering_stats_on_device`, then hands them to
    :class:`CenteredPotential` through ``precomputed_stats``. The returned object
    is the *NumPy* potential -- identical to ``build_potential(...)`` on the same
    bank -- so all downstream code (torch potential copy, separable-exact
    quadrature) is unchanged.
    """
    from src.natural_gradient_local_rate.potentials import (
        build_raw_feature, CenteredPotential,
    )
    raw = build_raw_feature(
        family, int(point["N_theta"]), int(point["seed"]),
        feature_multiplier=feature_multiplier, phi=phi)
    b, Mmat, emin_dev, emax_dev, mean_Z = centering_stats_on_device(
        raw, Z, device, dtype, chunk_size=chunk_size)
    stats = {
        "b": b, "M": Mmat, "emin_dev": emin_dev, "emax_dev": emax_dev,
        "mean_Z": mean_Z, "centering_samples": int(Z.shape[0]),
        "centering_seed": None,
    }
    return CenteredPotential(
        raw, float(point["kappa_target"]), safety_factor=safety_factor,
        precomputed_stats=stats)


# ---------------------------------------------------------------------------
# Chunking
# ---------------------------------------------------------------------------

def _to_t(x, device, dtype):
    """``torch.as_tensor`` on a writable contiguous float64 copy.

    ``np.array(..., copy=True)`` guarantees a writable, C-contiguous buffer even
    when ``x`` is a read-only view (e.g. a broadcast array or a slice of a
    memory-mapped bank), so ``torch.as_tensor`` never emits the "given NumPy
    array is not writable" UserWarning and never aliases a read-only buffer.
    """
    torch = get_torch()
    arr = np.array(x, dtype=np.float64, copy=True, order="C")
    return torch.as_tensor(arr, device=device, dtype=dtype)


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

def _accumulate_dense_generic(pot_t, Z, chunk_size=None):
    """Single chunked pass returning ``(G, T_mat, A)`` on device.

    ``G = (1/M) Vw^T Vq`` (p x p forward H_lin matrix), ``T_mat = (1/M) Z^T Vw``
    (N x p), ``A = (1/M) diag(W)^T (Z^2 - 1)`` (N x N diagonal benchmark).

    ``Vq[j] = sym_to_vec(Z_j Z_j^T - I)`` is built directly via
    :func:`torch_outer_minus_I_symvec` (no ``(B, N, N)`` outer-product tensor),
    and the upper-triangle indices are cached per ``(N, device)``.
    """
    torch = get_torch()
    N, M = Z.shape[1], Z.shape[0]
    p = torch_sym_dim(N)
    G = torch.zeros((p, p), device=Z.device, dtype=Z.dtype)
    T_mat = torch.zeros((N, p), device=Z.device, dtype=Z.dtype)
    A = torch.zeros((N, N), device=Z.device, dtype=Z.dtype)
    for s, e in _chunk_bounds(M, chunk_size):
        Zc = Z[s:e]
        W = pot_t.batch_hess(Zc)                          # (B, N, N)
        Vw = torch_sym_to_vec_batch(W)                    # (B, p)
        Vq = torch_outer_minus_I_symvec(Zc, N)            # (B, p)  == vec(ZcZc^T - I)
        G += Vw.T @ Vq
        T_mat += Zc.T @ Vw
        diagH = torch.diagonal(W, dim1=-2, dim2=-1)       # (B, N)
        A += diagH.T @ (Zc * Zc - 1.0)
    return G / M, T_mat / M, A / M


def _accumulate_dense_gaussian(pot_t, Z, chunk_size=None):
    """Analytic ``(G, T_mat, A)`` for the Gaussian potential (``W = I``).

    With ``W_j = I`` for every sample, ``Vw[j] = vec(I)`` is constant, so the
    Monte-Carlo accumulation collapses to two cheap reductions -- no per-sample
    Hessian batches and no ``(B, p)`` ``Vw``::

        S = (1/M) Z^T Z,    E = S - I,    m = mean_j Z_j
        G      = vec(I) vec(E)^T                       (p x p)
        H_sym  = 0.5 (vec(I) vec(E)^T + vec(E) vec(I)^T)
        T_mat  = m vec(I)^T   (column b is m * Tr(basis_b))
        A_{ik} = E_{kk}       (every row equals diag(E))

    These are exactly the limits of the generic accumulator for ``W = I``, so the
    Gaussian row matches the generic torch/CPU result to float64 round-off.
    """
    torch = get_torch()
    N, M = Z.shape[1], Z.shape[0]
    I = torch.eye(N, device=Z.device, dtype=Z.dtype)
    S = torch.zeros((N, N), device=Z.device, dtype=Z.dtype)
    zsum = torch.zeros(N, device=Z.device, dtype=Z.dtype)
    for s, e in _chunk_bounds(M, chunk_size):
        Zc = Z[s:e]
        S += Zc.T @ Zc
        zsum += Zc.sum(dim=0)
    S /= M
    mean_Z = zsum / M
    E = S - I
    vI = torch_sym_to_vec(I)                              # (p,)
    vE = torch_sym_to_vec(E)                              # (p,)
    G = torch.outer(vI, vE)                               # (p, p)
    T_mat = torch.outer(mean_Z, vI)                       # (N, p)
    diagE = torch.diagonal(E).contiguous()                # (N,)
    A = diagE.unsqueeze(0).repeat(N, 1)                   # (N, N), row i = diag(E)
    return G, T_mat, A


def _accumulate_dense_separable(pot_t, Z, chunk_size=None):
    """Fast ``(G, T_mat, A)`` for a separable potential (diagonal ``W``).

    ``W_j = diag(Wdiag_j)`` is diagonal, so ``Vw[j] = sym_to_vec(W_j)`` is nonzero
    only on the diagonal coordinates. This lets us skip materializing the full
    ``(B, N, N)`` Hessian and the ``(B, p)`` ``Vw``; we accumulate only the
    diagonal rows of ``G`` and diagonal columns of ``T_mat``::

        Wdiag      = diag(Hess V)                          (B, N)   (no dense W)
        G[diag, :] = (1/M) Wdiag^T Vq                      (N, p)
        T[:, diag] = (1/M) Z^T Wdiag                       (N, N)
        A          = (1/M) Wdiag^T (Z^2 - 1)               (N, N)

    where ``diag`` are the p-space indices of the diagonal matrix positions. The
    off-diagonal rows/columns of ``G`` / ``T_mat`` are exactly zero (a diagonal
    ``W`` has no off-diagonal vectorized component), so this reproduces the
    generic accumulator bit-for-bit.
    """
    torch = get_torch()
    N, M = Z.shape[1], Z.shape[0]
    p = torch_sym_dim(N)
    _, _, offdiag = _triu(N, Z.device)
    diag_idx = (~offdiag).nonzero(as_tuple=True)[0]       # (N,) p-space diag coords
    G_diagrows = torch.zeros((N, p), device=Z.device, dtype=Z.dtype)
    T_diagcols = torch.zeros((N, N), device=Z.device, dtype=Z.dtype)
    A = torch.zeros((N, N), device=Z.device, dtype=Z.dtype)
    for s, e in _chunk_bounds(M, chunk_size):
        Zc = Z[s:e]
        Wdiag = pot_t.batch_hess_diag(Zc)                 # (B, N)
        Vq = torch_outer_minus_I_symvec(Zc, N)            # (B, p)
        G_diagrows += Wdiag.T @ Vq                        # (N, p)
        T_diagcols += Zc.T @ Wdiag                        # (N, N)
        A += Wdiag.T @ (Zc * Zc - 1.0)                    # (N, N)
    G = torch.zeros((p, p), device=Z.device, dtype=Z.dtype)
    T_mat = torch.zeros((N, p), device=Z.device, dtype=Z.dtype)
    G[diag_idx, :] = G_diagrows / M
    T_mat[:, diag_idx] = T_diagcols / M
    return G, T_mat, A / M


def _resolve_accumulation_mode(pot_t):
    """Pick the accumulation fast path for a torch potential.

    ``gaussian_analytic`` for ``W = I``, ``separable_diagonal_fast`` for a
    diagonal-``W`` separable potential, else ``generic_dense``.
    """
    if isinstance(pot_t, TorchGaussianPotential):
        return "gaussian_analytic"
    if isinstance(pot_t, TorchSeparablePotential):
        return "separable_diagonal_fast"
    return "generic_dense"


_ACCUMULATORS = {
    "gaussian_analytic": _accumulate_dense_gaussian,
    "separable_diagonal_fast": _accumulate_dense_separable,
    "generic_dense": _accumulate_dense_generic,
}


def _accumulate_dense(pot_t, Z, chunk_size=None, mode=None):
    """Dispatch to the generic / Gaussian / separable accumulator.

    ``mode=None`` auto-selects via :func:`_resolve_accumulation_mode`; an explicit
    mode string (e.g. ``"generic_dense"``) forces a path -- used by tests that
    compare the fast paths against the generic one on the same bank.
    """
    if mode is None:
        mode = _resolve_accumulation_mode(pot_t)
    return _ACCUMULATORS[mode](pot_t, Z, chunk_size)


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
    mode = _resolve_accumulation_mode(pot_t)
    out["torch_accumulation_mode"] = mode
    G, T_mat, A = _accumulate_dense(pot_t, Z, chunk_size, mode=mode)
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
        "torch_accumulation_mode": "",
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
