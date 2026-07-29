"""Exact nonoptimal two-cycle: no fixed order-one step converges.

Companion counterexample to the bump train. Where the bump train blocks progress
for finitely many iterations, this one blocks it forever: for any fixed stepsize
``Delta t = gamma`` with ``2/kappa < gamma < 2`` it exhibits an admissible target
on which both scalar schemes have an exact period-two orbit, so the hitting time
is infinite rather than merely ``Omega(kappa^2)``.

Mechanism
---------
Both schemes share the mean update ``m+ = m - Delta t * c * b(m, c)``. The
covariance maps are stationary wherever ``c A(m, c) = 1``, so at such a state the
iteration reduces to the scalar recursion

    m+ = (1 - Delta t * c * b(m, c) / m) m .

The multiplier involves ``c`` times the *secant* slope ``b(m,c)/m``, while
stationarity of the covariance pins ``c`` to the reciprocal of the *tangent*
``A(m,c)``. Those two are independent: an even potential can have local curvature
``p`` at ``m = M`` while its averaged slope over ``[0, M]`` is any ``r`` with
``alpha <= r <= beta``. Choosing

    c = 1/p,     r/p = 2/gamma

gives multiplier ``1 - gamma * r/p = -1`` exactly, and evenness sends
``(M, c) -> (-M, c) -> (M, c)`` forever.

Realization with tight constants
--------------------------------
The construction is only meaningful if the potential's condition number really is
``kappa``. A potential with ``V'' in [p, P]`` and ``P/p = O(1)`` would reproduce
the cycle, but it would then just be the textbook explicit-step stability limit
``Delta t < 2/kappa_true`` restated at ``kappa_true = O(1)``. So ``V''`` here
*attains* both ``1`` and ``kappa``:

    V''(x) = kappa   on  [0, x1]           (max attained)
             1       on  [x1+d, x2]        (min attained)
             p       on  [x2+d, infinity)  (sets A(M, c) = p)

with ``C^inf`` transitions of width ``d``, and ``x1`` solved so that
``b(M, c) = r M`` exactly. The optimizer is ``(0, 1/kappa)``: ``V'' = kappa`` on a
neighbourhood of the origin of radius ``x1 = Theta(kappa)``.

The Hessian-Lipschitz constant is ``||V'''||_inf ~ (kappa-1) ||sigma'||/d``, so
``d = Theta(kappa/L_H)`` and ``M = Theta(kappa/L_H)``: the family is admissible for
any prescribed dimension-free ``L_H``, exactly as the bump train is.

Reading
-------
Combined with the bump train this closes the fixed-step question. For a fixed
``Delta t``:

* ``Delta t > 2/kappa``  -- this family gives an exact two-cycle, hitting time infinite;
* ``Delta t <= 2/kappa`` -- the bump train gives ``Omega(T kappa / Delta t) = Omega(kappa^2)``.

So the best fixed step is ``Theta(1/kappa)`` and the resulting ``Theta(kappa^2)`` is
intrinsic to fixed-step Fisher--Rao, not an artifact of the proof.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import brentq

from src.natural_gradient_discretization_stepsize.targets import gauss_hermite_nodes
from src.natural_gradient_fixed_step_barrier.bump_target import _smoothstep

_STEP_GRID_POINTS = 200_001


class _StepProfile:
    """Tabulated ``sigma`` (C^inf step), ``Sigma = int sigma`` and ``Xi = int Sigma``.

    ``sigma`` rises from ``0`` at ``t <= 0`` to ``1`` at ``t >= 1``. Outside ``[0,1]``
    the antiderivatives are continued analytically, so all three are exact there.
    """

    def __init__(self, n_points=_STEP_GRID_POINTS):
        t = np.linspace(0.0, 1.0, int(n_points))
        s = _smoothstep(t)
        Sg = np.concatenate([[0.0], cumulative_trapezoid(s, t)])
        Xi = np.concatenate([[0.0], cumulative_trapezoid(Sg, t)])
        self.t, self.s, self.Sg, self.Xi = t, s, Sg, Xi
        self.S1, self.X1 = float(Sg[-1]), float(Xi[-1])
        self.max_deriv = float(np.max(np.abs(np.gradient(s, t))))

    def sigma(self, u):
        return np.interp(u, self.t, self.s, left=0.0, right=1.0)

    def Sigma(self, u):
        u = np.asarray(u, dtype=np.float64)
        return np.where(u >= 1.0, self.S1 + (u - 1.0),
                        np.where(u <= 0.0, 0.0, np.interp(u, self.t, self.Sg)))

    def XiF(self, u):
        u = np.asarray(u, dtype=np.float64)
        return np.where(u >= 1.0, self.X1 + self.S1 * (u - 1.0) + 0.5 * (u - 1.0) ** 2,
                        np.where(u <= 0.0, 0.0, np.interp(u, self.t, self.Xi)))


_STEP = _StepProfile()


def cycle_constants(kappa, gamma):
    """``(q, p, r)`` of the two-cycle: ``q = 2/gamma``, ``p = (kappa+1)/(q+1)``, ``r = q p``.

    Requires ``2/kappa < gamma < 2`` so that ``1 < p, r < kappa``.
    """
    if not (2.0 / kappa < gamma < 2.0):
        raise ValueError(f"gamma={gamma} outside (2/kappa, 2) = ({2.0/kappa}, 2)")
    q = 2.0 / gamma
    p = (kappa + 1.0) / (q + 1.0)
    return float(q), float(p), float(q * p)


class TwoCycleTarget:
    """Even potential with ``min V'' = 1``, ``max V'' = kappa``, carrying the two-cycle.

    Parameters
    ----------
    kappa : float
        Condition number; both curvature bounds are attained (asserted below).
    gamma : float
        The fixed stepsize the cycle is built against, ``2/kappa < gamma < 2``.
    LH : float
        Hessian-Lipschitz budget, ``||V'''||_inf <= LH``.
    width_factor : float
        ``M0 / d``, the starting room for the plateaus. When the secant target ``r``
        sits close to ``kappa`` the mandatory ``1``-dip drags the achievable average
        below ``r``; the constructor then doubles ``width_factor`` (up to
        ``max_width_factor``) until the target is reachable, which only lengthens
        the plateaus and changes nothing else about the family.
    gh_nodes : int
        Gauss--Hermite order for the Gaussian averages.
    """

    name = "two_cycle"

    def __init__(self, kappa, gamma, LH=1.0, width_factor=24.0, gh_nodes=64,
                 max_width_factor=4096.0):
        self.kappa = float(kappa)
        self.gamma = float(gamma)
        self.LH = float(LH)
        self.q, self.p, self.r = cycle_constants(self.kappa, self.gamma)
        self.c = 1.0 / self.p

        # Transition width from the Hessian-Lipschitz budget: the steepest
        # transition is the kappa -> 1 one, of height kappa-1.
        self.d = (self.kappa - 1.0) * _STEP.max_deriv / self.LH
        # The Gaussian at M (std 1/sqrt(p) <= 1) must sit entirely in the p-region.
        self.Delta = 12.0 / math.sqrt(self.p)

        wf, solved = float(width_factor), False
        while wf <= float(max_width_factor):
            self.width_factor = wf
            self.M0 = wf * self.d
            self.M = self.M0 + self.Delta
            self.x2 = self.M0 - self.d
            lo, hi = 1e-9, self.x2 - 2.0 * self.d
            f = lambda x1: self._b_at(self.M, self.c, x1) - self.r * self.M
            if hi > lo and f(lo) * f(hi) < 0.0:
                self.x1 = float(brentq(f, lo, hi, xtol=1e-12, rtol=8.9e-16))
                solved = True
                break
            wf *= 2.0
        if not solved:
            raise ValueError(
                f"secant target r={self.r:.4g} unreachable at kappa={kappa}, gamma={gamma} "
                f"up to width_factor={max_width_factor}")

        self.gh_x, self.gh_w = gauss_hermite_nodes(int(gh_nodes))
        self.m_star = 0.0
        self.c_star = 1.0 / self.kappa      # V'' = kappa on a neighbourhood of 0
        self.energy_star = self._energy(self.m_star, self.c_star)

    # -- potential and derivatives (V even, V' odd) -------------------------

    def _u1(self, a):
        return (self.x1 + self.d - a) / self.d

    def _u2(self, a):
        return (self.x2 + self.d - a) / self.d

    def V2(self, x, x1=None):
        """``V''(a) = p + (1-p) sigma(u2) + (kappa-1) sigma(u1)``: kappa, then 1, then p."""
        x1 = self.x1 if x1 is None else x1
        a = np.abs(np.asarray(x, dtype=np.float64))
        u1 = (x1 + self.d - a) / self.d
        return (self.p + (1.0 - self.p) * _STEP.sigma(self._u2(a))
                + (self.kappa - 1.0) * _STEP.sigma(u1))

    def V1(self, x, x1=None):
        """``V'(a) = int_0^a V''``; odd in ``x``."""
        x1 = self.x1 if x1 is None else x1
        t = np.asarray(x, dtype=np.float64)
        a = np.abs(t)
        u1, u1_0 = (x1 + self.d - a) / self.d, (x1 + self.d) / self.d
        out = (self.p * a
               + (1.0 - self.p) * self.d * (_STEP.Sigma(self._u2(0.0)) - _STEP.Sigma(self._u2(a)))
               + (self.kappa - 1.0) * self.d * (_STEP.Sigma(u1_0) - _STEP.Sigma(u1)))
        return np.sign(t) * out

    def V0(self, x, x1=None):
        """``V(a) = int_0^a V'``, ``V(0) = 0``; even in ``x``."""
        x1 = self.x1 if x1 is None else x1
        a = np.abs(np.asarray(x, dtype=np.float64))
        u1, u1_0 = (x1 + self.d - a) / self.d, (x1 + self.d) / self.d
        u2, u2_0 = self._u2(a), self._u2(0.0)
        term2 = (1.0 - self.p) * self.d * (
            _STEP.Sigma(u2_0) * a - self.d * (_STEP.XiF(u2_0) - _STEP.XiF(u2)))
        term1 = (self.kappa - 1.0) * self.d * (
            _STEP.Sigma(u1_0) * a - self.d * (_STEP.XiF(u1_0) - _STEP.XiF(u1)))
        return 0.5 * self.p * a * a + term2 + term1

    # -- Gaussian averages and objective ------------------------------------

    def _nodes(self, m, c):
        return float(m) + math.sqrt(max(float(c), 0.0)) * self.gh_x

    def _b_at(self, m, c, x1):
        """``b(m,c)`` at a trial plateau length (used by the root solve at build time)."""
        x, w = gauss_hermite_nodes(64)
        y = float(m) + math.sqrt(max(float(c), 0.0)) * x
        return float(w @ self.V1(y, x1=x1))

    def b_A(self, m, c):
        """``(b, A) = (E[V'(X)], E[V''(X)])`` for ``X ~ N(m, c)``."""
        y = self._nodes(m, c)
        return float(self.gh_w @ self.V1(y)), float(self.gh_w @ self.V2(y))

    def _energy(self, m, c):
        return float(self.gh_w @ self.V0(self._nodes(m, c))) - 0.5 * math.log(float(c))

    def energy(self, m, c):
        return self._energy(m, c)

    def gap(self, m, c):
        return self._energy(m, c) - self.energy_star

    # -- diagnostics ---------------------------------------------------------

    def curvature_range(self, n=400_001):
        """``(min V'', max V'')`` on ``[0, 1.2 M]`` -- both bounds must be attained."""
        xs = np.linspace(0.0, 1.2 * self.M, int(n))
        v2 = self.V2(xs)
        return float(v2.min()), float(v2.max())

    def hessian_lipschitz(self, n=400_001):
        """``max |V'''|`` by finite differences on ``[0, 1.2 M]``."""
        xs = np.linspace(0.0, 1.2 * self.M, int(n))
        return float(np.max(np.abs(np.gradient(self.V2(xs), xs))))

    def metadata(self):
        v2min, v2max = self.curvature_range()
        b, A = self.b_A(self.M, self.c)
        return {
            "kappa": self.kappa, "gamma": self.gamma, "LH": self.LH,
            "q": self.q, "p": self.p, "r": self.r,
            "c_cycle": self.c, "M": self.M, "x1": self.x1, "x2": self.x2, "d": self.d,
            "V2_min": v2min, "V2_max": v2max,
            "kappa_realized": v2max / v2min,
            "hessian_lipschitz": self.hessian_lipschitz(),
            "A_at_M": A, "b_over_M_at_M": b / self.M, "cA_at_M": self.c * A,
            "mean_multiplier": 1.0 - self.gamma * self.c * b / self.M,
            "c_star": self.c_star, "energy_star": self.energy_star,
            "gap_at_cycle": self.gap(self.M, self.c),
        }
