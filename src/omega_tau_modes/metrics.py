"""
Metrics for evaluating the affine-invariant Gaussian gradient flow.

All metrics are relative to the target N(0, I_n).
Inputs (m, C) describe the current particle distribution N(m, C).
All computations use float64.
"""
import numpy as np
from src.omega_tau_modes.utils import spd_eigh, cosine_expectation


# ---------------------------------------------------------------------------
# KL divergence energy (primary convergence measure)
# ---------------------------------------------------------------------------

def kl_energy(m, C):
    """KL divergence KL(N(m, C) || N(0, I_n)).

    Exact formula:
        E = 0.5 * ( ||m||^2 + Tr(C) - logdet(C) - n )

    This equals zero iff m = 0 and C = I, and is always >= 0 (Gibbs inequality).
    """
    n = len(m)
    eigvals = np.linalg.eigvalsh(C)
    log_det_C = np.sum(np.log(np.maximum(eigvals, 1e-300)))
    return 0.5 * (float(np.dot(m, m)) + float(np.sum(eigvals)) - log_det_C - n)


# ---------------------------------------------------------------------------
# Full metric suite
# ---------------------------------------------------------------------------

def compute_all_metrics(m, C, E0, q, b):
    """Compute all diagnostic metrics for the current state (m, C).

    Args:
        m  : current mean, shape (n,)
        C  : current covariance, shape (n, n)
        E0 : initial KL energy (used for normalization; use 1.0 if E0 == 0)
        q  : test-function direction vector, shape (n,), unit norm
        b  : test-function phase offset (scalar)

    Returns:
        dict with keys:
            kl_energy, norm_energy, mean_error, cov_error, volume_error,
            shape_error, cosine_error, eig_min, eig_max, chi
    """
    n = len(m)
    eigvals, eigvecs = spd_eigh(C)  # eigvals >= 1e-300, sorted ascending

    log_eigvals = np.log(eigvals)
    trace_C     = float(np.sum(eigvals))
    log_det_C   = float(np.sum(log_eigvals))
    norm_m_sq   = float(np.dot(m, m))

    # 1. KL energy E = 0.5*(||m||^2 + Tr(C) - logdet(C) - n)
    E = 0.5 * (norm_m_sq + trace_C - log_det_C - n)
    E = max(E, 0.0)  # clamp tiny negative values from floating-point arithmetic

    # 2. Normalised energy E / E0
    norm_E = E / E0 if E0 > 1e-15 else float(E)

    # 3. Mean error ||m||_2
    mean_err = float(np.sqrt(norm_m_sq))

    # 4. Relative covariance error ||C - I||_F / sqrt(n)
    cov_err = float(np.linalg.norm(C - np.eye(n), 'fro') / np.sqrt(n))

    # 5. Volume error |logdet(C) / n|
    vol_err = abs(log_det_C / n)

    # 6. Shape error ||log C - (Tr(log C)/n) I||_F
    #    log C = Q diag(log lambda_i) Q^T
    #    The trace-free part of log C captures anisotropy independent of scale.
    mean_log_eig = log_det_C / n          # = Tr(log C) / n
    log_C = eigvecs @ np.diag(log_eigvals) @ eigvecs.T
    shape_mat = log_C - mean_log_eig * np.eye(n)
    shape_err = float(np.linalg.norm(shape_mat, 'fro'))

    # 7. Cosine test-function error
    #    E_{N(m,C)}[cos(q^T theta + b)] = exp(-0.5 q^T C q) cos(q^T m + b)
    #    True value under N(0,I):         exp(-0.5 ||q||^2) cos(b)
    cos_current = cosine_expectation(m, C, q, b)
    cos_true    = float(np.exp(-0.5 * float(np.dot(q, q))) * np.cos(b))
    cos_err     = abs(cos_current - cos_true)

    # 8. Eigenvalue extremes
    eig_min = float(eigvals[0])   # eigh returns sorted ascending
    eig_max = float(eigvals[-1])

    # 9. Trace dominance ratio chi for covariance residual
    #    Residuals r_i = 1 - lambda_i  (zero at target lambda_i = 1)
    #    chi = (sum r_i)^2 / (n * sum r_i^2)  in [0, 1]
    #    chi = 1  iff all residuals equal (maximally trace-dominated)
    #    chi = 1/n iff only one residual is non-zero (maximally shape-dominated)
    residuals = 1.0 - eigvals
    sum_r  = float(np.sum(residuals))
    sum_r2 = float(np.sum(residuals ** 2))
    chi = (sum_r ** 2 / (n * sum_r2)) if sum_r2 > 1e-300 else 1.0

    return {
        "kl_energy":    E,
        "norm_energy":  norm_E,
        "mean_error":   mean_err,
        "cov_error":    cov_err,
        "volume_error": vol_err,
        "shape_error":  shape_err,
        "cosine_error": cos_err,
        "eig_min":      eig_min,
        "eig_max":      eig_max,
        "chi":          chi,
    }
