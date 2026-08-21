"""The sharp bump-train counterexample of ``thm:sharp-disc`` (Appendix C).

This module builds the one-dimensional potential ``V_kappa`` of
:file:`manuscript/appendices/C-bump-train.tex` and exposes the two Gaussian
averages the scalar natural gradient schemes need,

    b(m, c) = E[V'(X)],   A(m, c) = E[V''(X)],   X ~ N(m, c),

together with the objective ``E(m, c) = E[V(X)] - 1/2 log c`` and its gap.

Construction (Appendix C, eq. bump-width-mass / bump-count-x0 / center-recursion)
--------------------------------------------------------------------------------
With a flat-top profile ``phi`` (``phi = 1`` on ``|u| <= 1/2``, ``0`` on
``|u| >= 1``), width ``w = (kappa-1) M_phi / LH`` and bump mass
``H = (kappa-1) w I_phi``, the centers solve

    x_{j+1} = x_j - s [ x_j + H (N - j + 1/2) ],     x_0 = Y kappa^2 / s,
    N = floor(T / s),

and ``V''(theta) = 1 + (kappa-1) sum_j [phi((theta-x_j)/w) + phi((theta+x_j)/w)]``.

**Generalized mean gain.** The appendix writes the recursion with the constant
``gamma/kappa^2``, which is exactly the per-step mean gain ``Delta t * c_kappa``
at the certified step ``Delta t = gamma/kappa`` and the matched covariance
``c_kappa = 1/kappa``. Here that constant is the free parameter ``s``, so the
same construction can be *retuned* to any stepsize: ``s = Delta t / kappa``. The
appendix family is ``s = gamma/kappa^2``. Every geometric statement of
``lem:bump-centers`` (``x_N >= 3/4 x_0``, spacing ``>= 3Y kappa^2 / 4``) is
``s``-free -- it only constrains ``T`` and ``Y`` -- and is asserted at
construction time, so the retuned family is the same counterexample aimed at a
different step. With ``N = T/s`` steps of blocked progress the family certifies
``Omega(T kappa / Delta t)`` iterations, which is ``Omega(kappa^2)`` at
``Delta t = Theta(1/kappa)`` and ``Omega(kappa)`` at ``Delta t = Theta(1)``.

All evaluations are exact up to Gauss--Hermite quadrature and the tabulated
profile antiderivatives; there is no Monte Carlo and no randomness.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import cumulative_trapezoid

from src.natural_gradient_discretization_stepsize.targets import gauss_hermite_nodes

# Resolution of the tabulated profile antiderivatives on [-1, 1].
_PROFILE_GRID_POINTS = 200_001


def _smoothstep(t):
    """``C^inf`` step: ``0`` for ``t <= 0``, ``1`` for ``t >= 1``, symmetric about ``1/2``."""
    t = np.asarray(t, dtype=np.float64)
    p = np.where(t > 0.0, np.exp(-1.0 / np.where(t > 0.0, t, 1.0)), 0.0)
    q = np.where(t < 1.0, np.exp(-1.0 / np.where(t < 1.0, 1.0 - t, 1.0)), 0.0)
    return p / (p + q)


def phi(u):
    """Flat-top profile of eq. (flat-bump): even, ``C^inf``, ``1`` on ``|u|<=1/2``, ``0`` on ``|u|>=1``."""
    a = np.abs(np.asarray(u, dtype=np.float64))
    out = np.zeros_like(a)
    out[a <= 0.5] = 1.0
    mid = (a > 0.5) & (a < 1.0)
    out[mid] = _smoothstep(2.0 * (1.0 - a[mid]))
    return out


class _Profile:
    """Tabulated ``phi``, ``Psi = int phi`` and ``Xi = int Psi``, with ``I_phi``, ``M_phi``.

    ``Psi(u) = 0`` for ``u <= -1`` and ``= I_phi`` for ``u >= 1``; ``Xi(u) = 0`` for
    ``u <= -1`` and ``= Xi(1) + I_phi (u - 1)`` for ``u >= 1``. Both are evaluated by
    linear interpolation on the tabulated grid, which is exact to ~1e-11 and enters
    only the bump contributions to ``V'`` and ``V`` (never the near-optimum regime,
    where ``V(theta) = theta^2/2`` identically).
    """

    def __init__(self, n_points=_PROFILE_GRID_POINTS):
        u = np.linspace(-1.0, 1.0, int(n_points))
        f = phi(u)
        psi = np.concatenate([[0.0], cumulative_trapezoid(f, u)])
        xi = np.concatenate([[0.0], cumulative_trapezoid(psi, u)])
        self.u = u
        self.psi = psi
        self.xi = xi
        self.I_phi = float(psi[-1])
        self.Xi1 = float(xi[-1])
        self.M_phi = float(np.max(np.abs(np.gradient(f, u))))

    def Psi(self, u):
        return np.interp(u, self.u, self.psi, left=0.0, right=self.I_phi)

    def Xi(self, u):
        u = np.asarray(u, dtype=np.float64)
        return np.where(u >= 1.0,
                        self.Xi1 + self.I_phi * (u - 1.0),
                        np.interp(u, self.u, self.xi, left=0.0, right=self.Xi1))


_PROFILE = _Profile()


class BumpTrain:
    """The potential ``V_kappa`` of Appendix C, tuned to per-step mean gain ``s``.

    Parameters
    ----------
    kappa : float
        Condition number; ``1 <= V'' <= kappa`` with both bounds attained.
    s : float
        Per-step mean gain ``Delta t * c_kappa`` the train is built against.
        The appendix family is ``s = gamma / kappa**2``.
    T, Y : float
        Geometry constants of eq. (bump-count-x0). ``T`` small and ``Y`` large
        enough that ``lem:bump-centers`` holds; both are asserted below.
    LH : float
        Hessian-Lipschitz budget, ``||V'''||_inf <= LH``, dimension-free.
    gh_nodes : int
        Gauss--Hermite order for the Gaussian averages.
    """

    def __init__(self, kappa, s, T=0.125, Y=4.0, LH=1.0, gh_nodes=32):
        if not (kappa > 1.0 and s > 0.0):
            raise ValueError("require kappa > 1 and s > 0")
        self.kappa = float(kappa)
        self.s = float(s)
        self.T = float(T)
        self.Y = float(Y)
        self.LH = float(LH)
        self.c_kappa = 1.0 / self.kappa

        p = _PROFILE
        self.I_phi, self.M_phi = p.I_phi, p.M_phi
        self.w = (self.kappa - 1.0) * self.M_phi / self.LH       # eq. bump-width-mass
        self.B = self.w * self.I_phi
        self.H = (self.kappa - 1.0) * self.B
        self.N = int(math.floor(self.T / self.s))                # eq. bump-count-x0
        self.x0 = self.Y * self.kappa ** 2 / self.s

        # Center recursion eq. (center-recursion); centers strictly decrease in j.
        x = np.empty(self.N + 1, dtype=np.float64)
        x[0] = self.x0
        for j in range(self.N):
            x[j + 1] = x[j] - self.s * (x[j] + self.H * (self.N - j + 0.5))
        self.centers = x

        # lem:bump-centers -- asserted, not assumed, at the tested (kappa, T, Y).
        gaps = -np.diff(x)
        self.min_spacing = float(gaps.min()) if gaps.size else math.inf
        self.x_N = float(x[-1])
        if self.x_N < 0.75 * self.x0:
            raise ValueError(
                f"lem:bump-centers fails: x_N/x_0 = {self.x_N / self.x0:.4f} < 0.75 "
                f"(kappa={kappa}, T={T}, Y={Y}); decrease T or increase Y")
        if self.min_spacing <= 2.0 * self.w:
            raise ValueError(
                f"bump supports overlap: min spacing {self.min_spacing:.4g} <= 2w "
                f"{2.0 * self.w:.4g} (kappa={kappa})")

        # Ascending centers + prefix sums: O(log N) evaluation of V, V', V''.
        self._asc = x[::-1].copy()
        self._csum = np.concatenate([[0.0], np.cumsum(self._asc)])

        self.gh_x, self.gh_w = gauss_hermite_nodes(int(gh_nodes))
        self.energy_star = self._energy_scalar(0.0, 1.0)

    # -- potential and derivatives (V is even; V' is odd) -------------------

    def _split(self, a):
        """``(n_below, csum_below, u_near)`` for ``a >= 0`` (array-valued).

        ``n_below`` counts centers at or below ``a - w`` (full bump mass passed);
        ``u_near`` is ``(a - x_near)/w`` for the single center within ``w`` of ``a``,
        set to ``-2`` when there is none -- left of the profile support, so
        ``phi = Psi = Xi = 0`` there and the term drops out. At most one such center
        exists because the spacing check above rules out overlap.
        """
        a = np.asarray(a, dtype=np.float64)
        n_below = np.searchsorted(self._asc, a - self.w, side="right")
        k = np.minimum(n_below, self._asc.size - 1)
        x_near = self._asc[k]
        u_near = np.where((n_below < self._asc.size) & (x_near < a + self.w),
                          (a - x_near) / self.w, -2.0)
        return n_below, self._csum[n_below], u_near

    def V2(self, theta):
        """``V''(theta) = 1 + (kappa-1) sum_j [phi(...) + phi(...)]``."""
        a = np.abs(np.asarray(theta, dtype=np.float64))
        _, _, u = self._split(a)
        return 1.0 + (self.kappa - 1.0) * phi(u)

    def V1(self, theta):
        """``V'(theta) = theta + (kappa-1) w sum_j Psi((theta - x_j)/w)``; odd in ``theta``."""
        t = np.asarray(theta, dtype=np.float64)
        a = np.abs(t)
        n_below, _, u = self._split(a)
        acc = self.I_phi * n_below + _PROFILE.Psi(u)
        return np.sign(t) * (a + (self.kappa - 1.0) * self.w * acc)

    def V0(self, theta):
        """``V(theta) = theta^2/2 + (kappa-1) w^2 sum_j Xi((theta - x_j)/w)``, ``V(0) = 0``."""
        a = np.abs(np.asarray(theta, dtype=np.float64))
        n_below, csum_below, u = self._split(a)
        # Centers fully below contribute Xi(v) = Xi1 + I_phi (v - 1) with v = (a - x_j)/w.
        acc = (n_below * (_PROFILE.Xi1 - self.I_phi)
               + self.I_phi * (n_below * a - csum_below) / self.w
               + _PROFILE.Xi(u))
        return 0.5 * a * a + (self.kappa - 1.0) * self.w ** 2 * acc

    # -- Gaussian averages and objective -----------------------------------

    def _nodes(self, m, c):
        return float(m) + math.sqrt(max(float(c), 0.0)) * self.gh_x

    def b_A(self, m, c):
        """``(b, A) = (E[V'(X)], E[V''(X)])`` for ``X ~ N(m, c)`` (Gauss--Hermite)."""
        y = self._nodes(m, c)
        return float(self.gh_w @ self.V1(y)), float(self.gh_w @ self.V2(y))

    def _energy_scalar(self, m, c):
        y = self._nodes(m, c)
        return float(self.gh_w @ self.V0(y)) - 0.5 * math.log(float(c))

    def energy(self, m, c):
        """``E(m, c) = E_{N(m,c)}[V] - 1/2 log c`` (eq. scalar-objective)."""
        return self._energy_scalar(m, c)

    def gap(self, m, c):
        """``DeltaE(m, c) = E(m, c) - E(a_star)``.

        ``a_star = (0, 1)``: ``V`` is even, and ``V(theta) = theta^2/2`` on the whole
        support of every ``N(0, c)`` with ``c = O(1)`` because the innermost bump
        sits at ``x_N ~ kappa^2/s``. This is asserted in the tests.
        """
        return self._energy_scalar(m, c) - self.energy_star

    def metadata(self):
        return {
            "kappa": self.kappa, "s": self.s, "T": self.T, "Y": self.Y, "LH": self.LH,
            "I_phi": self.I_phi, "M_phi": self.M_phi,
            "w": self.w, "H": self.H, "N": self.N, "x0": self.x0, "x_N": self.x_N,
            "x_N_over_x0": self.x_N / self.x0,
            "min_center_spacing": self.min_spacing,
            "c_kappa": self.c_kappa, "energy_star": self.energy_star,
            "gap_at_start": self.gap(self.x0, self.c_kappa),
        }
