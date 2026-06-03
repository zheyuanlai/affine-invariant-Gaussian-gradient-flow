"""
Metrics for evaluating the log-concave gradient flow experiments.

All metrics are expressed relative to the reference Gaussian VI optimum
a_star = (m_star, C_star).

Let  R = C_star^{-1/2} C C_star^{-1/2}  (covariance in whitened coordinates).

Metric catalogue
----------------
objective                 : F(m, C) ≈ mean_j V(theta_j) - 0.5 logdet C
objective_gap             : F - F_star   (raw, may be slightly negative due to MC)
normalized_objective_gap  : (F - F_star) / gap_0

whitened_mean_error       : ||C_star^{-1/2} (m - m_star)||_2

cov_error                 : ||R - I||_F / sqrt(n)
volume_error              : |logdet(R) / n|
shape_error               : ||logR - (Tr(logR)/n) I||_F

Stationarity residuals
  mean_residual     : ||g||_2      (g = E[grad V])
  cov_residual      : ||I - B||_F  where B = C^{1/2} S C^{1/2}
  trace_residual    : |Tr(I - B)| / sqrt(n)
  traceless_residual: ||I - B - (Tr(I-B)/n) I||_F
  chi               : (Tr(I-B))^2 / (n * ||I-B||_F^2)  — trace dominance

eig_min, eig_max    : eigenvalue extremes of C
cosine_error_to_star: |E cos(q^T theta + b)|_{m,C} - |_{m_star,C_star}|
"""
import numpy as np

from src.utils import spd_eigh, spd_sqrt, spd_invsqrt, cosine_expectation


def compute_lc_objective(m, L, Z, target):
    """Estimate the VI objective F(m, C) using fixed samples.

    F(m, C) = mean_j V(m + L z_j) - log det L
            = mean_j V(theta_j)   - sum_i log L_ii

    Args:
        m      : mean, shape (n,)
        L      : lower-triangular Cholesky factor of C, shape (n, n)
        Z      : standard-normal samples, shape (K, n)
        target : LogCoshTarget

    Returns:
        F : scalar float
    """
    from src.qmc_samples import push_forward
    Theta = push_forward(m, L, Z)
    V_vals = target.batch_value(Theta)
    return float(np.mean(V_vals)) - float(np.sum(np.log(np.diag(np.abs(L)))))


def compute_lc_metrics(m, C, g, S, m_star, C_star,
                        C_star_invsqrt, obj, F_star, gap0, q, b):
    """Compute all log-concave diagnostic metrics.

    Args:
        m              : current mean, (n,)
        C              : current covariance, (n, n) SPD
        g              : E[grad V(theta)] under N(m,C), (n,)
        S              : E[Hess V(theta)] under N(m,C), (n, n) SPD
        m_star         : reference mean, (n,)
        C_star         : reference covariance, (n, n) SPD
        C_star_invsqrt : C_star^{-1/2}, (n, n) — pre-computed once
        obj            : current VI objective F(m, C) (scalar)
        F_star         : reference objective F(m_star, C_star) (scalar)
        gap0           : initial objective gap F(m0,C0) - F_star (for normalisation)
        q              : test-function direction, (n,), unit norm
        b              : test-function phase offset (scalar)

    Returns:
        dict of metric name -> float
    """
    n = len(m)

    # ------------------------------------------------------------------
    # Eigendecomposition of C (used for multiple metrics)
    # ------------------------------------------------------------------
    eigvals_C, _ = spd_eigh(C)
    eig_min = float(eigvals_C[0])
    eig_max = float(eigvals_C[-1])

    # ------------------------------------------------------------------
    # Objective metrics
    # ------------------------------------------------------------------
    objective_gap = float(obj) - float(F_star)
    norm_gap = objective_gap / float(gap0) if abs(float(gap0)) > 1e-15 else float(objective_gap)

    # ------------------------------------------------------------------
    # Whitened mean error:  ||C_star^{-1/2} (m - m_star)||
    # ------------------------------------------------------------------
    dm = m - m_star
    wh_mean_err = float(np.linalg.norm(C_star_invsqrt @ dm))

    # ------------------------------------------------------------------
    # Whitened covariance:  R = C_star^{-1/2} C C_star^{-1/2}
    # ------------------------------------------------------------------
    R = C_star_invsqrt @ C @ C_star_invsqrt
    R = 0.5 * (R + R.T)

    eigvals_R, eigvecs_R = spd_eigh(R)
    log_eigvals_R = np.log(eigvals_R)
    log_det_R = float(np.sum(log_eigvals_R))

    # Relative covariance error ||R - I||_F / sqrt(n)
    cov_err = float(np.linalg.norm(R - np.eye(n), "fro") / np.sqrt(n))

    # Volume error |logdet(R) / n|
    vol_err = abs(log_det_R / n)

    # Shape error ||logR - (Tr(logR)/n) I||_F
    mean_log_R = log_det_R / n
    logR = eigvecs_R @ np.diag(log_eigvals_R) @ eigvecs_R.T
    shape_mat = logR - mean_log_R * np.eye(n)
    shape_err = float(np.linalg.norm(shape_mat, "fro"))

    # ------------------------------------------------------------------
    # Stationarity residuals
    # ------------------------------------------------------------------
    # Mean residual: ||g||_2  (zero at the fixed point)
    mean_residual = float(np.linalg.norm(g))

    # Whitened Hessian B = C^{1/2} S C^{1/2}
    C_sqrt = spd_sqrt(C)
    B = C_sqrt @ S @ C_sqrt
    B = 0.5 * (B + B.T)

    # Covariance residual matrix:  I - B  (zero at fixed point where B = I)
    res_mat = np.eye(n, dtype=np.float64) - B
    cov_residual = float(np.linalg.norm(res_mat, "fro"))

    # Trace residual: |Tr(res_mat)| / sqrt(n)
    tr_res = float(np.trace(res_mat))
    trace_residual = abs(tr_res) / np.sqrt(n)

    # Traceless residual: ||res_mat - (Tr/n) I||_F
    traceless_mat = res_mat - (tr_res / n) * np.eye(n)
    traceless_residual = float(np.linalg.norm(traceless_mat, "fro"))

    # Chi — trace dominance ratio
    cov_res_sq = float(np.sum(res_mat ** 2))           # = ||res_mat||_F^2
    chi = (tr_res ** 2 / (n * cov_res_sq)) if cov_res_sq > 1e-300 else 1.0

    # ------------------------------------------------------------------
    # Cosine test-function error relative to a_star
    # ------------------------------------------------------------------
    cos_now  = cosine_expectation(m,      C,      q, b)
    cos_star = cosine_expectation(m_star, C_star, q, b)
    cosine_error_to_star = abs(cos_now - cos_star)

    return {
        "objective":                  float(obj),
        "objective_gap":              objective_gap,
        "normalized_objective_gap":   norm_gap,
        "whitened_mean_error":        wh_mean_err,
        "cov_error":                  cov_err,
        "volume_error":               vol_err,
        "shape_error":                shape_err,
        "mean_residual":              mean_residual,
        "cov_residual":               cov_residual,
        "trace_residual":             trace_residual,
        "traceless_residual":         traceless_residual,
        "chi":                        float(chi),
        "eig_min":                    eig_min,
        "eig_max":                    eig_max,
        "cosine_error_to_star":       cosine_error_to_star,
    }
