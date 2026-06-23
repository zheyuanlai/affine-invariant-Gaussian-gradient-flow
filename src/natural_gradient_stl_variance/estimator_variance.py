"""Experiment A -- estimator-level variance comparison at fixed Gaussian states.

At a fixed state ``a = (m, C)`` we draw ``theta_j ~ N(m, C)`` and compare the
baseline and STL mean estimators

    b_base(theta) = score_post(theta),
    b_stl(theta)  = score_post(theta) + C^{-1}(theta - m),

reporting, for each estimator, the empirical mean, the bias against the
deterministic reference ``E_q[score_post]``, the Euclidean trace variance
``E||b - E b||^2`` and the Fisher--Rao mean-block variance
``E[(b - E b)^T C (b - E b)]``, together with the STL/baseline variance ratios.

For the well-specified Gaussian target the variances have closed forms used to
validate the Monte Carlo estimates: with ``D = C^{-1} - A``,

    Var_eucl(b_base) = Tr(A C A),     Var_FR(b_base) = Tr(C A C A),
    Var_eucl(b_stl)  = Tr(D C D),     Var_FR(b_stl)  = Tr(C D C D),

so at the optimum ``C = A^{-1}`` (``D = 0``) the STL variance is exactly zero.

The expensive sampling/reduction is batched over the ``M`` samples through the
:class:`~src.natural_gradient_stl_variance.linalg.ArrayBackend` (NumPy or a single
CUDA device) and chunked so that ``M`` may be made large without large memory.
"""
from __future__ import annotations

import numpy as np

from src.common.spd import symmetrize
from src.natural_gradient_stl_variance.linalg import ArrayBackend


def _score_backend(bk, kind, theta, a_bk, tau):
    """Batched ``score_post(theta)`` on the backend (``theta`` shape ``(M, d)``)."""
    base = -(theta * a_bk)
    if kind == "log_cosh" and tau != 0.0:
        base = base - tau * bk.tanh(theta)
    return base


def evaluate_state(target, state_name, m, C, distance, M, seed,
                   backend, chunk_size=200_000):
    """Estimator-level metrics for one ``(target, state, seed)`` over ``M`` samples.

    Returns a flat dict (one CSV row). Sampling and reductions use ``backend``;
    the reference expectation and the Gaussian closed forms use NumPy.
    """
    bk = backend
    m = np.asarray(m, dtype=np.float64)
    C = symmetrize(C)
    d = C.shape[0]

    C_bk = bk.asarray(C)
    m_bk = bk.asarray(m)
    a_bk = bk.asarray(target.a_diag)
    tau = float(target.tau)
    C_sqrt = bk.sqrtm_sym(C_bk)
    C_inv = bk.inv(C_bk)

    gen = bk.generator(seed)

    # Chunked accumulation of first and (Euclidean / FR) second moments.
    sum_base = bk.asarray(np.zeros(d))
    sum_stl = bk.asarray(np.zeros(d))
    sq_base = sq_stl = 0.0           # sum_j ||b_j||^2
    fr_base = fr_stl = 0.0           # sum_j b_j^T C b_j
    remaining = int(M)
    while remaining > 0:
        b = min(int(chunk_size), remaining)
        Z = bk.randn((b, d), gen)
        theta = m_bk + Z @ C_sqrt           # (b, d); C_sqrt symmetric
        score = _score_backend(bk, target.kind, theta, a_bk, tau)
        corr = (theta - m_bk) @ C_inv       # C^{-1}(theta - m), C_inv symmetric
        b_base = score
        b_stl = score + corr

        sum_base = sum_base + bk.sum(b_base, axis=0)
        sum_stl = sum_stl + bk.sum(b_stl, axis=0)
        sq_base += float(bk.to_numpy(bk.sum(b_base * b_base)))
        sq_stl += float(bk.to_numpy(bk.sum(b_stl * b_stl)))
        fr_base += float(bk.to_numpy(bk.sum(b_base * (b_base @ C_bk))))
        fr_stl += float(bk.to_numpy(bk.sum(b_stl * (b_stl @ C_bk))))
        remaining -= b

    mean_base = bk.to_numpy(sum_base) / M
    mean_stl = bk.to_numpy(sum_stl) / M
    # Var = E||b||^2 - ||E b||^2  (and the FR analogue with the C-inner product).
    var_eucl_base = sq_base / M - float(mean_base @ mean_base)
    var_eucl_stl = sq_stl / M - float(mean_stl @ mean_stl)
    var_fr_base = fr_base / M - float(mean_base @ C @ mean_base)
    var_fr_stl = fr_stl / M - float(mean_stl @ C @ mean_stl)
    # Numerical floor: tiny negative values from cancellation -> clamp at 0.
    var_eucl_base = max(var_eucl_base, 0.0)
    var_eucl_stl = max(var_eucl_stl, 0.0)
    var_fr_base = max(var_fr_base, 0.0)
    var_fr_stl = max(var_fr_stl, 0.0)

    ref = target.exp_score(m, C)
    bias_base = float(np.linalg.norm(mean_base - ref))
    bias_stl = float(np.linalg.norm(mean_stl - ref))

    def _ratio(num, den):
        return float(num / den) if den > 0 else float("nan")

    row = {
        "target_name": target.name, "kind": target.kind,
        "d": int(d), "kappa": float(target.kappa), "tau": tau,
        "state": state_name, "distance_to_optimum": float(distance),
        "M": int(M), "seed": int(seed),
        "backend": bk.backend, "device": bk.device_str(),
        "bias_base": bias_base, "bias_stl": bias_stl,
        "var_eucl_base": var_eucl_base, "var_eucl_stl": var_eucl_stl,
        "var_fr_base": var_fr_base, "var_fr_stl": var_fr_stl,
        "ratio_euclidean": _ratio(var_eucl_stl, var_eucl_base),
        "ratio_fr": _ratio(var_fr_stl, var_fr_base),
    }

    # Exact Gaussian closed forms (validation columns).
    if target.kind == "gaussian":
        A = np.diag(target.a_diag)
        Cinv_np = np.linalg.inv(C)
        D = Cinv_np - A
        AC = A @ C
        DC = D @ C
        row.update({
            "var_eucl_base_exact": float(np.trace(A @ C @ A)),
            "var_fr_base_exact": float(np.trace(AC @ AC)),
            "var_eucl_stl_exact": float(np.trace(D @ C @ D)),
            "var_fr_stl_exact": float(np.trace(DC @ DC)),
        })
    else:
        row.update({
            "var_eucl_base_exact": float("nan"),
            "var_fr_base_exact": float("nan"),
            "var_eucl_stl_exact": float("nan"),
            "var_fr_stl_exact": float("nan"),
        })
    return row


def make_backend(backend="numpy", device="cpu", dtype="float64"):
    """Convenience constructor mirroring the script CLI."""
    return ArrayBackend(backend=backend, device=device, dtype=dtype)
