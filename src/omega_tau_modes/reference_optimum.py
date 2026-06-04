"""
Reference Gaussian VI optimum for a log-concave target.

Minimises F(m, C) = E_{N(m,C)}[V(theta)] - 0.5 log det C  over all (m, C)
with C SPD, using a Cholesky parameterization C = L L^T.

Cholesky parameterization
-------------------------
L is lower-triangular with positive diagonal.  We store:
  - m            : shape (n,)          — mean
  - off-diagonal entries of L          — unconstrained reals
  - diagonal log-entries eta_i         — so L_ii = exp(eta_i) > 0

This gives a single unconstrained parameter vector of length n + n*(n+1)/2.

Objective (fixed-sample estimate)
----------------------------------
    F(m, L) = mean_j V(m + L z_j) - log det L
            = mean_j V(theta_j)   - sum_i eta_i

Gradients
---------
    grad_m F = mean_j grad V(theta_j)

    grad_L F = (1/K) sum_j grad V(theta_j) z_j^T - L^{-T}

For the packed vector:
  - off-diagonal L_ij   -> partial F / partial L_ij
  - diagonal eta_i      -> (partial F / partial L_ii) * L_ii
                        = (grad_L)_ii * exp(eta_i)

Saved file format (.npz)
------------------------
    m_star, C_star, L_star, F_star, metadata dict items as arrays.
"""
import json
import os

import numpy as np
import scipy.linalg
from scipy.optimize import minimize

from src.omega_tau_modes.qmc_samples import push_forward


# ---------------------------------------------------------------------------
# Cholesky packing / unpacking
# ---------------------------------------------------------------------------

def _pack(m, L):
    """Pack (m, L) into a 1-D unconstrained parameter vector.

    Layout: [ m | L_lower_offdiag | log(L_diag) ]
    """
    n = len(m)
    diag_log = np.log(np.diag(L))                   # (n,)
    # Lower-triangular off-diagonal entries, row-major
    rows, cols = np.tril_indices(n, k=-1)
    offdiag = L[rows, cols]                          # (n*(n-1)/2,)
    return np.concatenate([m, offdiag, diag_log])


def _unpack(params, n):
    """Unpack parameter vector back to (m, L).

    Returns m (n,) and L (n,n) lower-triangular with positive diagonal.
    """
    m = params[:n].copy()
    n_offdiag = n * (n - 1) // 2
    offdiag = params[n: n + n_offdiag]
    diag_log = params[n + n_offdiag:]

    L = np.zeros((n, n), dtype=np.float64)
    rows, cols = np.tril_indices(n, k=-1)
    L[rows, cols] = offdiag
    np.fill_diagonal(L, np.exp(diag_log))
    return m, L


# ---------------------------------------------------------------------------
# Objective and gradient
# ---------------------------------------------------------------------------

def _objective_and_grad(params, target, Z):
    """Return (F, grad_params) for the VI objective.

    F(m, L) = mean_j V(theta_j) - sum_i log L_ii
    """
    n = target.n
    K = Z.shape[0]
    m, L = _unpack(params, n)

    # Sample push-forward: theta_j = m + L z_j
    Theta = push_forward(m, L, Z)                    # (K, n)

    # Objective value
    V_vals = target.batch_value(Theta)               # (K,)
    F = float(np.mean(V_vals)) - float(np.sum(np.log(np.diag(L))))

    # Gradient w.r.t. m:  mean_j grad V(theta_j)
    GV = target.batch_grad(Theta)                    # (K, n)
    grad_m = np.mean(GV, axis=0)                     # (n,)

    # Gradient w.r.t. L:  (1/K) sum_j grad V(theta_j) z_j^T  -  L^{-T}
    # = grad_m z_j^T averaged, then subtract L^{-T}
    grad_L_data = GV[:, :, np.newaxis] * Z[:, np.newaxis, :]  # (K, n, n)
    grad_L = np.mean(grad_L_data, axis=0)            # (n, n)
    # Subtract L^{-T}: only needed on lower-triangular part
    # L^{-T} = (L^{-1})^T; scipy.linalg.solve_triangular is efficient
    L_inv = scipy.linalg.solve_triangular(L, np.eye(n), lower=True)
    grad_L -= L_inv.T                                # (n, n)

    # Pack gradient
    rows, cols = np.tril_indices(n, k=-1)
    grad_offdiag = grad_L[rows, cols]
    # Diagonal: dF/d eta_i = dF/d L_ii * L_ii
    grad_eta = np.diag(grad_L) * np.diag(L)

    grad_params = np.concatenate([grad_m, grad_offdiag, grad_eta])
    return F, grad_params


# ---------------------------------------------------------------------------
# Main optimiser
# ---------------------------------------------------------------------------

def compute_reference_optimum(target, Z_ref, maxiter: int = 2000, gtol: float = 1e-6):
    """Find the reference Gaussian VI optimum a_star = (m_star, C_star).

    Args:
        target  : LogCoshTarget instance
        Z_ref   : reference samples, shape (K_ref, n), standard normal
        maxiter : L-BFGS-B iteration budget
        gtol    : gradient tolerance for convergence

    Returns:
        dict with keys:
            m_star, C_star, L_star, F_star,
            grad_m_norm, cov_residual_norm,
            n, rho, m_features, seed, K_ref, converged
    """
    n = target.n
    K_ref = Z_ref.shape[0]

    # Initial point: m=0, L=I
    m0 = np.zeros(n, dtype=np.float64)
    L0 = np.eye(n, dtype=np.float64)
    params0 = _pack(m0, L0)

    result = minimize(
        fun=_objective_and_grad,
        x0=params0,
        args=(target, Z_ref),
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": maxiter, "gtol": gtol, "ftol": 1e-15},
    )

    m_star, L_star = _unpack(result.x, n)
    C_star = L_star @ L_star.T
    C_star = 0.5 * (C_star + C_star.T)              # symmetrise

    F_star = float(result.fun)

    # Stationarity check
    _, grad = _objective_and_grad(result.x, target, Z_ref)
    m_grad, _, _ = grad[:n], grad[n:n + n*(n-1)//2], grad[n + n*(n-1)//2:]
    grad_m_norm = float(np.linalg.norm(m_grad))

    # Covariance residual: how far is L_star L_star^T from a natural fixed point
    Theta_ref = push_forward(m_star, L_star, Z_ref)
    S_ref = np.mean(target.batch_hess(Theta_ref), axis=0)       # (n, n)
    cov_residual_norm = float(np.linalg.norm(C_star @ S_ref - np.eye(n), "fro"))

    return {
        "m_star":              m_star,
        "C_star":              C_star,
        "L_star":              L_star,
        "F_star":              F_star,
        "grad_m_norm":         grad_m_norm,
        "cov_residual_norm":   cov_residual_norm,
        "n":                   n,
        "rho":                 target.rho,
        "m_features":          target.m_features,
        "target_seed":         target.seed,
        "K_ref":               K_ref,
        "converged":           bool(result.success),
        "optim_message":       result.message,
    }


# ---------------------------------------------------------------------------
# Save / load
# ---------------------------------------------------------------------------

def save_reference_optimum(path: str, result: dict):
    """Save reference optimum to .npz (arrays) + .json (metadata)."""
    os.makedirs(os.path.dirname(path) if os.path.dirname(path) else ".", exist_ok=True)
    np.savez(
        path,
        m_star=result["m_star"],
        C_star=result["C_star"],
        L_star=result["L_star"],
    )
    meta = {k: v for k, v in result.items()
            if k not in ("m_star", "C_star", "L_star")}
    # Convert numpy scalars to Python for JSON serialisation
    meta = {k: (v.item() if hasattr(v, "item") else v) for k, v in meta.items()}
    meta_path = path.replace(".npz", "_meta.json")
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"  Saved reference optimum -> {path}")
    print(f"  Metadata               -> {meta_path}")


def load_reference_optimum(path: str) -> dict:
    """Load reference optimum previously saved by save_reference_optimum."""
    data = np.load(path)
    result = {
        "m_star": data["m_star"],
        "C_star": data["C_star"],
        "L_star": data["L_star"],
    }
    meta_path = path.replace(".npz", "_meta.json")
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            result.update(json.load(f))
    return result


def load_or_compute(path: str, target, Z_ref, force: bool = False,
                    maxiter: int = 2000, gtol: float = 1e-6) -> dict:
    """Load cached reference optimum or compute and save if missing/forced."""
    if not force and os.path.exists(path):
        print(f"  Loading cached reference optimum from {path}")
        return load_reference_optimum(path)
    print(f"  Computing reference optimum (K_ref={Z_ref.shape[0]}) ...")
    result = compute_reference_optimum(target, Z_ref, maxiter=maxiter, gtol=gtol)
    save_reference_optimum(path, result)
    return result
