"""Self-adjointness of L_star in Fisher--Rao-whitened packed coordinates.

With ``H_sym`` in the covariance block and the metric-scaled packing
``y = (u, sym_to_vec(X)/sqrt2)``, the matrix representation of ``L_star`` must be
symmetric under the *Euclidean* dot product, which is the precondition for
``scipy.sparse.linalg.eigsh``.
"""
import numpy as np
import pytest

from src.common.monte_carlo import gaussian_samples
from src.natural_gradient_local_rate.potentials import build_potential
from src.natural_gradient_local_rate.operators import (
    make_L_star_operator, self_adjoint_error_L_star,
)
from src.natural_gradient_local_rate.linearized_rate import _to_matrix


@pytest.mark.parametrize("family", ["separable", "random_feature", "radial_tail"])
def test_L_star_self_adjoint_error_near_precision(family):
    N = 6
    pot = build_potential(family, N, kappa_target=12.0, seed=0, centering_samples=4096)
    Z = gaussian_samples(N, 8192, seed=1, antithetic=True)
    err = self_adjoint_error_L_star(pot, Z, n_probe=8)
    assert err < 1e-9, f"L_star self-adjointness error {err}"


@pytest.mark.parametrize("N", [4, 9])
def test_packed_matrix_is_symmetric(N):
    # dot(y1, A y2) == dot(A y1, y2): the dense packed matrix is symmetric.
    pot = build_potential("random_feature", N, kappa_target=8.0, seed=0,
                          centering_samples=4096)
    Z = gaussian_samples(N, 4096, seed=2, antithetic=True)
    op = make_L_star_operator(pot, Z, covariance_weight="half")
    Mat = _to_matrix(op)                       # no symmetrization applied
    asym = float(np.max(np.abs(Mat - Mat.T)))
    assert asym < 1e-9, f"asymmetry {asym}"
    rng = np.random.default_rng(0)
    for _ in range(5):
        y1 = rng.standard_normal(op.shape[0])
        y2 = rng.standard_normal(op.shape[0])
        assert float(y1 @ op.matvec(y2)) == pytest.approx(
            float(op.matvec(y1) @ y2), rel=1e-9, abs=1e-9)


def test_plain_coordinates_are_not_symmetric():
    # The unweighted 'plain' representation is generally asymmetric; this
    # documents why the Fisher--Rao 'half' packing is required for eigsh.
    N = 5
    pot = build_potential("random_feature", N, kappa_target=20.0, seed=0,
                          centering_samples=4096)
    Z = gaussian_samples(N, 4096, seed=1, antithetic=True)
    op = make_L_star_operator(pot, Z, covariance_weight="plain")
    Mat = _to_matrix(op)
    assert np.max(np.abs(Mat - Mat.T)) > 1e-6
