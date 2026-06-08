"""Wasserstein--Fisher--Rao (WFR) gradient flow on the Gaussian manifold.

Section 5 experiment group. The WFR Gaussian flow augments the Fisher--Rao
(natural-gradient) flow with a Wasserstein transport contribution of strength
``lambda_t``; its forward--backward discretization interleaves one Wasserstein
step and one Fisher--Rao step per iteration. This package implements the targets,
the two half-steps, the transport schedules, and the run-level driver.
"""
