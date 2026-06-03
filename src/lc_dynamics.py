"""
Affine-invariant gradient flow — one discrete update step for a general
strongly log-concave target.

Given pre-computed Monte Carlo estimates:
    g = E_{N(m,C)}[grad V(theta)]     shape (n,)
    S = E_{N(m,C)}[Hess V(theta)]     shape (n, n)

the Riemannian-distance discretization gives:

Mean update:
    m_next = m - dt * C @ g

Covariance update:
    Let  B   = C^{1/2} S C^{1/2}          (whitened Hessian)
         tau_tr = Tr(C S) = Tr(B)         (trace of whitened Hessian)
         alpha  = (omega + tau * tau_tr) / (omega + n * tau)

    Exponent matrix:  M = dt/(2*omega) * (-B + alpha * I)

    C_next = C^{1/2} expm(M) C^{1/2}

Unlike the Gaussian case, the exponent M is NOT generally a function of
C's eigenvectors, so scipy.linalg.expm is used for the matrix exponential.

Note on sign convention:
    This module expects g = E[grad V] and S = E[Hess V],
    i.e., the *negative* of E[grad log pi] and E[Hess log pi].
    The Gaussian target case writes E[grad log pi] = -m and E[Hess log pi] = -I,
    which corresponds to g = m and S = I here.

Validity constraints:
    omega > 0   and   omega + n * tau > 0
"""
import numpy as np
import scipy.linalg

from src.utils import validate_params, spd_sqrt


def logconcave_step(m, C, g, S, dt, omega, tau):
    """One discrete step of the affine-invariant gradient flow for a general target.

    Args:
        m     : current mean, shape (n,), float64
        C     : current covariance, shape (n, n), float64, SPD
        g     : MC estimate of E[grad V(theta)], shape (n,)
        S     : MC estimate of E[Hess V(theta)], shape (n, n), SPD
        dt    : positive time-step
        omega : positive float — covariance update rate parameter
        tau   : real float — trace-weighting; must satisfy omega + n*tau > 0

    Returns:
        (m_next, C_next) : updated mean and covariance, both float64

    Raises:
        ValueError : if omega <= 0 or omega + n*tau <= 0
    """
    n = len(m)

    if not validate_params(omega, tau, n):
        raise ValueError(
            f"Invalid parameters: omega={omega}, tau={tau}, n={n}. "
            f"Required: omega > 0 and omega + n*tau > 0 "
            f"(omega + n*tau = {omega + n * tau:.6g})."
        )

    # ------------------------------------------------------------------
    # Mean update:  m_next = m - dt * C * g
    # (compare Gaussian case: g = m, so m_next = m - dt * C @ m)
    # ------------------------------------------------------------------
    m_next = m - dt * (C @ g)

    # ------------------------------------------------------------------
    # Covariance update via matrix exponential
    # ------------------------------------------------------------------
    C_sqrt = spd_sqrt(C)                             # C^{1/2}

    # Whitened Hessian: B = C^{1/2} S C^{1/2}
    B = C_sqrt @ S @ C_sqrt
    B = 0.5 * (B + B.T)                             # symmetrise

    # Trace of whitened Hessian = Tr(C S)
    tau_tr = float(np.trace(B))

    # Equilibrium scalar: alpha = (omega + tau * tau_tr) / (omega + n * tau)
    alpha = (omega + tau * tau_tr) / (omega + n * tau)

    # Exponent matrix: M = dt/(2*omega) * (-B + alpha * I)
    scale = dt / (2.0 * omega)
    M = scale * (-B + alpha * np.eye(n, dtype=np.float64))

    # Matrix exponential (general, via Padé approximation)
    expM = scipy.linalg.expm(M)

    # C_next = C^{1/2} expm(M) C^{1/2}
    C_next = C_sqrt @ expM @ C_sqrt
    C_next = 0.5 * (C_next + C_next.T)              # symmetrise

    return m_next, C_next
