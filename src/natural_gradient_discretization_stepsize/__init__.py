"""Riemannian vs KL discretization of the Gaussian natural gradient flow.

This experiment group compares two time discretizations of the *same* Gaussian
natural gradient flow for the variational problem
``min_{m,C} KL(N(m, C) || rho_post)``:

* the **Riemannian-distance** discretization (geodesic / exponential-map step),
* the **KL/Bregman** discretization (rational inverse-covariance step).

Both schemes share the explicit mean update ``m_{n+1} = m_n + dt * C_n g_n`` and
differ only in the covariance update. The central question is whether the
stricter sufficient stepsize condition for the KL scheme in the current proof is
a genuine restriction or a proof artifact. All experiments are one- or
two-dimensional, deterministic, and CPU-only.

See ``reports/natural_gradient_discretization_stepsize_report.tex`` for the
write-up and the repository ``README.md`` for the reproduction commands.

Notation (matching the natural-gradient manuscript)::

    rho_a = N(m, C),  a = (m, C),  C in SPD(d)
    E(a)  = KL(rho_a || rho_post)
    g_n   = E_{rho_a}[grad log rho_post] = -E[grad V]
    H_n   = E_{rho_a}[grad^2 log rho_post] = -E[grad^2 V]
    dm/dt = C g,   dC/dt = C + C H C.
"""
from src.natural_gradient_discretization_stepsize.targets import (
    GaussianPosteriorTarget,
    LiteratureLogconcaveTarget,
    SmoothLogconcaveTarget,
    ScalarGaussianTarget,
    build_target,
    TARGET_NAMES,
)
from src.natural_gradient_discretization_stepsize.methods import (
    riemannian_cov_step,
    kl_cov_step,
    mean_step,
    discretization_step,
    METHOD_NAMES,
)
from src.natural_gradient_discretization_stepsize.metrics import (
    classify_run,
    theory_stepsize_bounds,
)

__all__ = [
    "GaussianPosteriorTarget",
    "LiteratureLogconcaveTarget",
    "SmoothLogconcaveTarget",
    "ScalarGaussianTarget",
    "build_target",
    "TARGET_NAMES",
    "riemannian_cov_step",
    "kl_cov_step",
    "mean_step",
    "discretization_step",
    "METHOD_NAMES",
    "classify_run",
    "theory_stepsize_bounds",
]
