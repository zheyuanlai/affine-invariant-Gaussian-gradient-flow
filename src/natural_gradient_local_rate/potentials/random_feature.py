"""Dense random-feature coupled potential (the main general ``V`` family).

Same ridge-sum form as :mod:`additive_index` but with a dense bank of random
directions and a larger default width ``r = feature_multiplier * N_theta``
(default multiplier 4). This is the primary stress test for general coupled
potentials.
"""
from __future__ import annotations

from src.natural_gradient_local_rate.potentials.additive_index import RidgeSumFeature


class RandomFeaturePotential(RidgeSumFeature):
    family = "random_feature"

    def __init__(self, N_theta, r=None, seed=0, feature_multiplier=4,
                 feature_scale=1.0, phi="log_cosh", coeff_dist="signs",
                 offset_scale=1.0):
        if r is None:
            r = int(feature_multiplier) * int(N_theta)
        self.feature_multiplier = int(feature_multiplier)
        super().__init__(N_theta, r, seed=seed, feature_scale=feature_scale,
                         phi=phi, coeff_dist=coeff_dist, offset_scale=offset_scale)

    def raw_metadata(self):
        md = super().raw_metadata()
        md["feature_multiplier"] = self.feature_multiplier
        return md
