"""Fisher--Rao metric scaling of the symmetric-matrix vectorization.

These checks underpin the eigsh validity of the local-rate operator: the packed
tangent ``y = (u, sym_to_vec(X)/sqrt2)`` must turn the Euclidean dot product into
the Fisher--Rao inner product ``u.u + 0.5 Tr(X X)``.
"""
import numpy as np
import pytest

from src.common.symspace import (
    sym_to_vec, vec_to_sym, sym_inner, sym_norm, random_symmetric,
    pack_tangent_fr, unpack_tangent_fr, fisher_rao_inner, fisher_rao_norm,
    tangent_dim,
)


@pytest.mark.parametrize("N", [1, 2, 3, 5, 8])
def test_sym_vec_preserves_frobenius_inner_product(N):
    rng = np.random.default_rng(N)
    X = random_symmetric(N, rng)
    Y = random_symmetric(N, rng)
    assert float(np.dot(sym_to_vec(X), sym_to_vec(Y))) == pytest.approx(
        float(np.trace(X @ Y)), abs=1e-12)
    assert vec_to_sym(sym_to_vec(X), N) == pytest.approx(X, abs=1e-13)


@pytest.mark.parametrize("N", [1, 2, 4, 7])
def test_packed_tangent_preserves_fisher_rao_norm(N):
    rng = np.random.default_rng(10 + N)
    u = rng.standard_normal(N)
    X = random_symmetric(N, rng)
    y = pack_tangent_fr(u, X)
    assert y.shape == (tangent_dim(N),)
    # ||y||^2 (Euclidean) == ||u||^2 + 0.5 ||X||_F^2 (Fisher--Rao)
    fr_sq = float(u @ u) + 0.5 * sym_norm(X) ** 2
    assert float(np.dot(y, y)) == pytest.approx(fr_sq, abs=1e-12)
    assert fisher_rao_norm((u, X)) ** 2 == pytest.approx(fr_sq, abs=1e-12)


@pytest.mark.parametrize("N", [2, 4, 6])
def test_packed_dot_equals_fisher_rao_inner(N):
    rng = np.random.default_rng(99 + N)
    u1, u2 = rng.standard_normal(N), rng.standard_normal(N)
    X1, X2 = random_symmetric(N, rng), random_symmetric(N, rng)
    y1, y2 = pack_tangent_fr(u1, X1), pack_tangent_fr(u2, X2)
    lhs = float(np.dot(y1, y2))
    rhs = fisher_rao_inner((u1, X1), (u2, X2))
    assert lhs == pytest.approx(rhs, abs=1e-12)
    # and the explicit definition u.u + 0.5 Tr(X1 X2)
    assert rhs == pytest.approx(float(u1 @ u2) + 0.5 * sym_inner(X1, X2), abs=1e-12)


@pytest.mark.parametrize("N", [1, 3, 5])
def test_pack_unpack_roundtrip(N):
    rng = np.random.default_rng(7 * N)
    u = rng.standard_normal(N)
    X = random_symmetric(N, rng)
    u2, X2 = unpack_tangent_fr(pack_tangent_fr(u, X), N)
    assert u2 == pytest.approx(u, abs=1e-13)
    assert X2 == pytest.approx(X, abs=1e-13)
