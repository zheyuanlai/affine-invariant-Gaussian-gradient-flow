"""Source package for the Gaussian gradient flow experiments.

Subpackages
-----------
common                       : shared numerical infrastructure
omega_tau_modes              : the (omega, tau) affine-invariant flow experiments
natural_gradient_local_rate  : the natural-gradient local convergence-rate study

This top-level package intentionally performs no eager imports so that importing
one experiment group never pulls in another (or matplotlib) as a side effect.
"""
