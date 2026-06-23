"""The fixed Gaussian states ``a = (m, C)`` used by the estimator-level study.

Every state is defined relative to the target's Gaussian VI optimum
``a_star = (m_star, C_star)`` (with ``m_star = 0`` for both targets):

* ``optimum``       -- ``(m_star, C_star)``;
* ``near``          -- mean offset ``rho_near * s`` along the marginal-std
  direction ``s = sqrt(diag(C_star))``, covariance ``C_star``;
* ``medium``        -- mean offset ``rho_medium * s``, covariance ``C_star``;
* ``far``           -- mean offset ``rho_far * s``, covariance ``C_star``;
* ``underdispersed``-- ``(m_star, under_scale * C_star)``;
* ``overdispersed`` -- ``(m_star, over_scale * C_star)``.

Measuring the mean offset in units of the marginal standard deviation makes the
"distance to optimum" comparable across dimensions and condition numbers. The
default multipliers are conservative and overridable from the config.
"""
from __future__ import annotations

import numpy as np

from src.common.spd import symmetrize

STATE_NAMES = ["optimum", "near", "medium", "far", "underdispersed", "overdispersed"]

DEFAULT_STATE_PARAMS = {
    "rho_near": 0.25,
    "rho_medium": 1.0,
    "rho_far": 4.0,
    "under_scale": 0.25,
    "over_scale": 4.0,
}


def build_states(m_star, C_star, params=None):
    """Return ``[(name, m, C, distance_to_optimum), ...]`` for the six states.

    ``distance_to_optimum`` is the mean-offset multiplier ``rho`` for the
    mean-perturbation states (0 for ``optimum``) and ``nan`` for the
    dispersion states (which perturb the covariance, not the mean).
    """
    p = dict(DEFAULT_STATE_PARAMS)
    if params:
        p.update({k: float(v) for k, v in params.items() if k in DEFAULT_STATE_PARAMS})

    m_star = np.asarray(m_star, dtype=np.float64)
    C_star = symmetrize(C_star)
    s = np.sqrt(np.clip(np.diag(C_star), 0.0, None))  # marginal std direction

    states = [
        ("optimum", m_star.copy(), C_star.copy(), 0.0),
        ("near", m_star + p["rho_near"] * s, C_star.copy(), p["rho_near"]),
        ("medium", m_star + p["rho_medium"] * s, C_star.copy(), p["rho_medium"]),
        ("far", m_star + p["rho_far"] * s, C_star.copy(), p["rho_far"]),
        ("underdispersed", m_star.copy(),
         symmetrize(p["under_scale"] * C_star), float("nan")),
        ("overdispersed", m_star.copy(),
         symmetrize(p["over_scale"] * C_star), float("nan")),
    ]
    return states, p
