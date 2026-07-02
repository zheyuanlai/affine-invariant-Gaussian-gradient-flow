"""Covariance-bootstrap experiment: enhanced global rates for the Gaussian
natural-gradient (Fisher--Rao) flow.

Modules
-------
targets     : exact diagonal Gaussian and smooth log-cosh targets (curvature
              convention A = E[grad^2 V]); exact energy gap, a_star, and the
              stochastic Hessian-fluctuation intensity Psi.
methods      : the Fisher--Rao Riemannian/KL covariance steps, the shared mean
              step, and the one-step Wasserstein/Bures bootstrap (all in the
              curvature convention; equal to the H_disc = -A discretization steps).
envelopes    : the KL / Riemannian covariance lower envelopes L_n (recurrence and
              closed form), the frozen lower bound, and dynamic/frozen contraction
              factors (reusing q_riem / q_kl from src.common.theory_constants).
metrics      : first-crossing (N_cov), hitting times, and tail-floor helpers.
runner       : deterministic trajectory simulation and the four experiment drivers.
stochastic   : batched stochastic STL trajectories (NumPy or one CUDA torch device).
plotting     : figure builders shared by plot_results.py and make_report_assets.py.
"""
