"""Deterministic sharp-rate experiments for Gaussian Fisher--Rao flow."""

from .bump_train import BumpTrainTarget
from .local_targets import ChenLogBump, RidgeTarget, ShellTarget
from .spiral import SpiralValleyTarget

__all__ = [
    "BumpTrainTarget",
    "ChenLogBump",
    "RidgeTarget",
    "ShellTarget",
    "SpiralValleyTarget",
]
