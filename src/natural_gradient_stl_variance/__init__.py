"""Sticking-the-Landing (STL) variance reduction for stochastic Gaussian
natural-gradient variational inference.

Sixth experiment group. Stochastic Gaussian natural-gradient schemes estimate the
Gaussian expectations ``g = E_q[score_post]`` and ``H = E_q[Hess_log_post]`` from
samples ``theta ~ N(m, C)``. The *Sticking-the-Landing* (STL) trick replaces the
plain posterior-score mean estimator

    b_base(theta) = score_post(theta)

by the score-residual estimator

    b_stl(theta)  = score_post(theta) - grad_theta log q(theta)
                  = score_post(theta) + C^{-1}(theta - m),

which is unbiased for ``E_q[score_post]`` (because ``E_q[C^{-1}(theta-m)] = 0``)
and *vanishes pointwise* at a well-specified Gaussian optimum, where
``score_post(theta) = -C^{-1}(theta-m)`` exactly. For the *direct* Hessian
estimator used here, the analogous covariance-block residual rewrite
``K = S + C^{-1}`` is only an algebraic restatement of the covariance update (see
:mod:`src.natural_gradient_stl_variance.estimators`), so the experiment isolates
the mean-block STL estimator.

This package implements:

* :mod:`targets` -- a well-specified anisotropic Gaussian target and a
  misspecified smooth strongly log-concave separable ``log cosh`` target;
* :mod:`estimators` -- the baseline / STL mean estimators and the Riemannian and
  KL covariance updates, in both the direct and the Hessian-residual forms;
* :mod:`linalg` -- a small NumPy/torch batched-linear-algebra backend (float64,
  CPU or a single CUDA device);
* :mod:`states` -- the fixed Gaussian states used by the estimator-level study;
* :mod:`estimator_variance` -- Experiment A (estimator-level variance);
* :mod:`algorithm` -- Experiment B (algorithm-level stochastic trajectories);
* :mod:`metrics` -- tail / noise-floor metrics and W2;
* :mod:`plotting` -- display helpers and the report figure builders.

All expectations that define the reference optima are deterministic (closed form
for the Gaussian target, Gauss--Hermite quadrature / a fixed-point solve for the
``log cosh`` target). Every stochastic run is reproducible from explicit seeds and
uses float64 by default.
"""
from src.natural_gradient_stl_variance.targets import (
    GaussianTarget,
    LogCoshTarget,
    build_target,
    TARGET_NAMES,
)
from src.natural_gradient_stl_variance.estimators import (
    baseline_mean_estimator,
    stl_mean_estimator,
    riemannian_cov_step,
    kl_cov_step,
    riemannian_cov_step_residual,
    kl_cov_step_residual,
    METHODS,
    SCHEMES,
)

__all__ = [
    "GaussianTarget",
    "LogCoshTarget",
    "build_target",
    "TARGET_NAMES",
    "baseline_mean_estimator",
    "stl_mean_estimator",
    "riemannian_cov_step",
    "kl_cov_step",
    "riemannian_cov_step_residual",
    "kl_cov_step_residual",
    "METHODS",
    "SCHEMES",
]
