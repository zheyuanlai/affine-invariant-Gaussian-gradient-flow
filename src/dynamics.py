"""
Affine-invariant Gaussian gradient flow — one discrete update step.

Target distribution: pi = N(0, I_n).

Exact expectations under the current particle distribution q = N(m, C):
    E_q[ grad log pi(theta) ]  = E_q[ -theta ]  = -m
    E_q[ Hess log pi(theta) ]  =                  -I

Riemannian-distance discretization with parameters (omega, tau):

Mean update (affine-invariant gradient descent on the mean):
    m_{k+1} = m_k - dt * C_k @ m_k

Covariance update (matrix exponential form):
    C_{k+1} = C_k^{1/2}
                exp( dt/(2*omega) * [ -C_k + alpha * I ] )
              C_k^{1/2}

where
    alpha = (omega + tau * Tr(C_k)) / (omega + n * tau)

Because the exponent commutes with C_k (both are functions of C_k's
eigenvectors), the update reduces to a scalar rescaling of each eigenvalue:

    lambda_i_{k+1} = lambda_i * exp( dt/(2*omega) * (-lambda_i + alpha) )

with the same eigenvectors Q preserved.

Validity constraints:
    omega > 0   and   omega + n * tau > 0
"""
import numpy as np
from src.utils import validate_params, spd_eigh


def gaussian_step(m, C, dt, omega, tau):
    """One discrete step of the affine-invariant Gaussian gradient flow.

    Target: N(0, I_n).  All expectations are exact (no Monte Carlo).

    Args:
        m     : current mean, shape (n,), dtype float64
        C     : current covariance, shape (n, n), dtype float64, must be SPD
        dt    : positive time-step
        omega : positive float — overall covariance update rate
        tau   : real float — trace-weighting; must satisfy omega + n*tau > 0

    Returns:
        (m_next, C_next) : updated (mean, covariance), both float64

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
    # Mean update
    # dm/dt = C * E[grad log pi] = C * (-m)   =>   m_{k+1} = m_k - dt*C*m_k
    # ------------------------------------------------------------------
    m_next = m - dt * (C @ m)

    # ------------------------------------------------------------------
    # Covariance update via eigenvalues
    # ------------------------------------------------------------------
    eigvals, eigvecs = spd_eigh(C)

    # Trace-weighted equilibrium target eigenvalue
    trace_C = np.sum(eigvals)
    alpha = (omega + tau * trace_C) / (omega + n * tau)

    # Multiplicative eigenvalue update
    scale = dt / (2.0 * omega)
    eigvals_next = eigvals * np.exp(scale * (-eigvals + alpha))

    # Reconstruct SPD matrix; symmetrize to suppress floating-point skew
    C_next = eigvecs @ np.diag(eigvals_next) @ eigvecs.T
    C_next = 0.5 * (C_next + C_next.T)

    return m_next, C_next
