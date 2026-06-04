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

    centering_kwargs = dict(
        centering_samples=centering_samples,
        centering_chunk=centering_chunk,
        safety_factor=safety_factor,
        la_mode=la_mode,
        Z_ref=Z_ref,
    )

    if family == "separable":
        raw = SeparablePotential(N_theta, seed=seed, feature_scale=feature_scale, phi=phi)
    elif family == "additive_index":
        r = additive_r if additive_r is not None else feature_multiplier * N_theta
        raw = AdditiveIndexPotential(N_theta, r=r, seed=seed,
                                     feature_scale=feature_scale, phi=phi)
    elif family == "random_feature":
        raw = RandomFeaturePotential(N_theta, r=additive_r, seed=seed,
                                     feature_multiplier=feature_multiplier,
                                     feature_scale=feature_scale, phi=phi)
    elif family == "radial_tail":
        raw = RadialTailPotential(N_theta, seed=seed,
                                  scale=radial_scale, shift=radial_shift, phi=phi)
    else:
        raise ValueError(f"unknown potential family '{family}', have {ALL_FAMILIES}")

    return CenteredPotential(raw, kappa_target, **centering_kwargs)


__all__ = [
    "BasePotential", "GaussianPotential", "CenteredPotential", "RawFeaturePotential",
    "SeparablePotential", "AdditiveIndexPotential", "RandomFeaturePotential",
    "RadialTailPotential", "RidgeSumFeature", "get_nonlinearity",
    "build_potential", "FEATURE_FAMILIES", "ALL_FAMILIES",
]
