"""Natural-gradient local convergence-rate experiment.

Studies whether the *local* convergence rate of the Gaussian natural gradient
flow near equilibrium is essentially dimension-free. Notation follows the
natural-gradient manuscript:

* ``N_theta``                    : dimension
* ``rho_post(theta) ~ exp(-V)``  : target (posterior)
* ``rho_a = N(m, C)``, ``a=(m,C)``
* ``E(a) = KL(rho_a || rho_post)``
* equilibrium ``a_star = (0, I)`` in equilibrium-whitened coordinates.

See ``docs/specs/natural_gradient_local_rate_spec.md`` (source of truth) and
``reports/natural_gradient_local_rate_notes.md``.
"""
