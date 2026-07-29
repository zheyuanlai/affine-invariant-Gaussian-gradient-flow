"""Sharp bump-train counterexample: is the ``O(kappa^2)`` iteration count intrinsic?

The manuscript proves (``thm:sharp-disc``, Appendix C) that the globally certified
fixed step ``Delta t = gamma/kappa`` forces ``Omega(kappa^2)`` iterations for a
constant-factor reduction of the objective, on an explicit one-dimensional
bump-train potential with exact condition number ``kappa``, matched initial
covariance ``c_0 = 1/kappa`` and dimension-free ``||V'''||_inf``. The remark after
that theorem leaves open whether a ``kappa``-independent step -- fixed or adaptive
-- escapes the obstruction, which is where the gap to the Bures--Wasserstein
``O(kappa)`` complexity (arXiv:2304.05398) sits.

This group runs the counterexample itself. Two families of the same construction:

``manuscript``  the train of Appendix C, built against the certified step
                ``Delta t = gamma/kappa`` (per-step mean gain ``s = gamma/kappa^2``);
``retuned``     the same construction rebuilt against each arm's own step
                (``s = Delta t / kappa``), i.e. the adversary adapts to the method.

crossed with the two schemes (Riemannian, KL) and the three stepsize arms
(``theory``, ``const``, ``relcurv``; see :mod:`.runner`). The measured quantity is
the iteration count to a constant-factor gap reduction as a function of ``kappa``,
reported as a log-log slope: slope ``2`` reproduces ``thm:sharp-disc``, slope ``1``
is the Bures--Wasserstein scaling.

Everything is one-dimensional, deterministic and CPU-only: closed-form potential,
Gauss--Hermite Gaussian averages, no Monte Carlo.
"""
from src.natural_gradient_sharp_bump.bump_target import BumpTrain, phi
from src.natural_gradient_sharp_bump.runner import (
    ARMS, SCHEMES, nominal_dt, scheme_step, simulate, step_size,
)

__all__ = [
    "BumpTrain", "phi",
    "ARMS", "SCHEMES", "nominal_dt", "scheme_step", "simulate", "step_size",
]
