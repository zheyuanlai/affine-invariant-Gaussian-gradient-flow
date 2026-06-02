"""
AffineInvariantGaussianGradientFlow
====================================
Study of parameter effects (omega, tau) in affine-invariant Gaussian gradient
flows with Gaussian target N(0, I_n), where all expectations are exact.
"""
from src.dynamics import gaussian_step
from src.metrics import compute_all_metrics, kl_energy
from src.initializations import get_initialization, INIT_NAMES
from src.utils import validate_params, make_q_vector
