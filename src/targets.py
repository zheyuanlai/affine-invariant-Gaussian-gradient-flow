"""
Log-concave non-Gaussian target distribution.

V_rho(x) = 0.5 ||x||^2 + (rho / m) * sum_{ell=1}^m log cosh(a_ell^T x)

where a_ell are fixed random unit vectors in R^n and m = 4*n by default.

The posterior is  pi(x) ∝ exp(-V_rho(x)).

Sign convention:
    grad log pi  = -grad V
    Hess log pi  = -Hess V

Strong log-concavity:
    Hess V(x) = I + (rho/m) * sum_ell sech^2(a_ell^T x) a_ell a_ell^T  >= I
since sech^2 >= 0.  The minimum eigenvalue of Hess V is >= 1.

Smoothness:
    Hess V(x) = I + (rho/m) * sum_ell sech^2(a_ell^T x) a_ell a_ell^T  <= I + rho*I
since sech^2 <= 1 and A A^T <= I (rows are unit vectors).

Stable sech^2 implementation:
    sech^2(z) = 1/cosh(z)^2 = 1 - tanh^2(z)
avoids overflow for large |z| because tanh is bounded.
"""
import numpy as np


class LogCoshTarget:
    """Strongly log-concave target pi(x) ∝ exp(-V_rho(x)).

    Parameters
    ----------
    n          : int   — dimension
    rho        : float — coupling strength (>= 0)
    m_features : int   — number of measurement vectors (default 4*n)
    seed       : int   — RNG seed for reproducible A matrix

    Attributes
    ----------
    A : ndarray, shape (m_features, n)
        Rows are unit vectors drawn from N(0, I) and normalised.
    """

    def __init__(self, n: int, rho: float, m_features: int = None, seed: int = 0):
        self.n = int(n)
        self.rho = float(rho)
        self.m_features = int(m_features) if m_features is not None else 4 * self.n
        self.seed = int(seed)

        # Draw measurement directions and normalise rows to unit norm
        rng = np.random.default_rng(self.seed)
        A_raw = rng.standard_normal((self.m_features, self.n))
        norms = np.linalg.norm(A_raw, axis=1, keepdims=True)
        self.A = (A_raw / norms).astype(np.float64)   # shape (m, n)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _projections(self, x):
        """Return (m,) array of a_ell^T x."""
        return self.A @ x                              # (m,)

    @staticmethod
    def _sech2(z):
        """Numerically stable sech^2(z) = 1 - tanh^2(z)."""
        t = np.tanh(z)
        return 1.0 - t * t

    # ------------------------------------------------------------------
    # Single-point interface
    # ------------------------------------------------------------------

    def value(self, x):
        """Compute V_rho(x) = 0.5 ||x||^2 + (rho/m) sum_ell log cosh(a_ell^T x)."""
        proj = self._projections(x)                    # (m,)
        return 0.5 * float(np.dot(x, x)) + (self.rho / self.m_features) * float(
            np.sum(np.log(np.cosh(proj)))
        )

    def grad(self, x):
        """Compute grad V_rho(x) = x + (rho/m) sum_ell tanh(a_ell^T x) a_ell."""
        proj = self._projections(x)                    # (m,)
        weights = np.tanh(proj)                        # (m,)
        return x + (self.rho / self.m_features) * (self.A.T @ weights)

    def hess(self, x):
        """Compute Hess V_rho(x) = I + (rho/m) sum_ell sech^2(a_ell^T x) a_ell a_ell^T."""
        proj = self._projections(x)                    # (m,)
        w = self._sech2(proj)                          # (m,), values in [0,1]
        # Weighted outer-product sum: A^T diag(w) A
        wA = w[:, np.newaxis] * self.A                 # (m, n)
        return np.eye(self.n, dtype=np.float64) + (self.rho / self.m_features) * (self.A.T @ wA)

    # ------------------------------------------------------------------
    # Batch interface  (X has shape (K, n))
    # ------------------------------------------------------------------

    def batch_value(self, X):
        """Return (K,) array of V_rho values for each row of X."""
        # projections: (K, m)
        proj = X @ self.A.T                            # (K, m)
        quad = 0.5 * np.sum(X ** 2, axis=1)           # (K,)
        logcosh = (self.rho / self.m_features) * np.sum(np.log(np.cosh(proj)), axis=1)
        return quad + logcosh

    def batch_grad(self, X):
        """Return (K, n) array of grad V_rho for each row of X."""
        proj = X @ self.A.T                            # (K, m)
        # tanh: (K, m); each row weighted by tanh of its projections
        T = np.tanh(proj)                              # (K, m)
        return X + (self.rho / self.m_features) * (T @ self.A)

    def batch_hess(self, X):
        """Return (K, n, n) array of Hess V_rho for each row of X.

        Hess_k = I + (rho/m) * A^T diag(sech^2(A x_k)) A
        """
        K = X.shape[0]
        proj = X @ self.A.T                            # (K, m)
        W = self._sech2(proj)                          # (K, m)
        # For each k: (rho/m) * A^T diag(W[k]) A
        # Vectorised: einsum or loop — einsum is cleaner for memory
        # wA[k,ell,:] = W[k,ell] * A[ell,:]
        wA = W[:, :, np.newaxis] * self.A[np.newaxis, :, :]  # (K, m, n)
        cov_part = (self.rho / self.m_features) * np.einsum("kmi,kmj->kij", wA, self.A[np.newaxis])
        I_stack = np.broadcast_to(np.eye(self.n, dtype=np.float64), (K, self.n, self.n)).copy()
        return I_stack + cov_part
