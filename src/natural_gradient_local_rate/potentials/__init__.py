"""Potential families for the natural-gradient local-rate experiment.

A *potential* defines ``V`` with ``rho_post(theta) ~ exp(-V(theta))``. Use
:func:`build_potential` to construct one from a config-style description.
"""
from __future__ import annotations

from src.natural_gradient_local_rate.potentials.base import (
    BasePotential, GaussianPotential, CenteredPotential, RawFeaturePotential,
    get_nonlinearity,
)
from src.natural_gradient_local_rate.potentials.separable import SeparablePotential
from src.natural_gradient_local_rate.potentials.additive_index import (
    AdditiveIndexPotential, RidgeSumFeature,
)
from src.natural_gradient_local_rate.potentials.random_feature import RandomFeaturePotential
from src.natural_gradient_local_rate.potentials.radial_tail import RadialTailPotential

# Families that go through the centered-feature construction (exclude "gaussian").
FEATURE_FAMILIES = ["separable", "additive_index", "random_feature", "radial_tail"]
ALL_FAMILIES = ["gaussian"] + FEATURE_FAMILIES


def build_raw_feature(family, N_theta, seed, *, feature_multiplier=4,
                      phi="log_cosh", feature_scale=1.0, additive_r=None,
                      radial_scale=1.0, radial_shift=1.0):
    """Construct the *raw* (un-centered) feature for a coupled family.

    Returns the :class:`RawFeaturePotential` for ``separable`` /
    ``additive_index`` / ``random_feature`` / ``radial_tail``. ``gaussian`` has
    no raw feature and raises ``ValueError`` (callers handle it specially). This
    is the single source of truth for raw-feature parameters, shared by the CPU
    :func:`build_potential` and the GPU-centering builder so they stay in sync.
    """
    family = str(family)
    if family == "gaussian":
        raise ValueError("gaussian has no raw feature; use GaussianPotential")
    if family == "separable":
        return SeparablePotential(N_theta, seed=seed, feature_scale=feature_scale, phi=phi)
    if family == "additive_index":
        r = additive_r if additive_r is not None else feature_multiplier * N_theta
        return AdditiveIndexPotential(N_theta, r=r, seed=seed,
                                      feature_scale=feature_scale, phi=phi)
    if family == "random_feature":
        return RandomFeaturePotential(N_theta, r=additive_r, seed=seed,
                                      feature_multiplier=feature_multiplier,
                                      feature_scale=feature_scale, phi=phi)
    if family == "radial_tail":
        return RadialTailPotential(N_theta, seed=seed,
                                   scale=radial_scale, shift=radial_shift, phi=phi)
    raise ValueError(f"unknown potential family '{family}', have {ALL_FAMILIES}")


def build_potential(family, N_theta, kappa_target, seed, *,
                    feature_multiplier=4, phi="log_cosh", safety_factor=2.0,
                    centering_samples=8192, centering_chunk=2048,
                    feature_scale=1.0, additive_r=None,
                    radial_scale=1.0, radial_shift=1.0, la_mode="auto",
                    Z_ref=None):
    """Build a potential from a config-style description.

    Parameters mirror the YAML config: ``family`` is one of
    :data:`ALL_FAMILIES`; ``feature_multiplier`` sets ``r = mult * N_theta`` for
    the dense random-feature family; ``safety_factor`` multiplies the empirical
    ``L_A`` estimate for coupled families. Pass ``Z_ref`` (the run's operator /
    flow sample bank) so that centering shares it and the equilibrium
    ``a_star = (0, I)`` is the *exact* discrete fixed point.
    """
    family = str(family)
    if family == "gaussian":
        return GaussianPotential(N_theta, seed=seed)

    raw = build_raw_feature(
        family, N_theta, seed, feature_multiplier=feature_multiplier, phi=phi,
        feature_scale=feature_scale, additive_r=additive_r,
        radial_scale=radial_scale, radial_shift=radial_shift)

    return CenteredPotential(
        raw, kappa_target,
        centering_samples=centering_samples,
        centering_chunk=centering_chunk,
        safety_factor=safety_factor,
        la_mode=la_mode,
        Z_ref=Z_ref,
    )


__all__ = [
    "BasePotential", "GaussianPotential", "CenteredPotential", "RawFeaturePotential",
    "SeparablePotential", "AdditiveIndexPotential", "RandomFeaturePotential",
    "RadialTailPotential", "RidgeSumFeature", "get_nonlinearity",
    "build_potential", "build_raw_feature", "FEATURE_FAMILIES", "ALL_FAMILIES",
]
