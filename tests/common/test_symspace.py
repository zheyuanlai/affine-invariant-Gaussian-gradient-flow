"""Tests for the Frobenius-isometric symmetric-matrix vectorization."""
import numpy as np
import pytest

from src.common.symspace import (
    sym_dim, sym_to_vec, vec_to_sym, sym_inner, sym_norm,
    random_symmetric, pack_tangent, unpack_tangent, tangent_dim,
)


@pytest.mark.parametrize("N", [1, 2, 3, 5, 8])
def test_roundtrip(N):
    rng = np.random.default_rng(N)
    X = random_symmetric(N, rng)
    X2 = vec_to_sym(sym_to_vec(X), N)
    assert np.allclose(X2, X, atol=1e-14)
    assert sym_to_vec(X).shape == (sym_dim(N),)


@pytest.mark.parametrize("N", [1, 2, 3, 5, 8])
def test_frobenius_isometry(N):
    rng = np.random.default_rng(100 + N)
    X = random_symmetric(N, rng)
    Y = random_symmetric(N, rng)
    # dot(vec X, vec Y) == Tr(X Y)
    lhs = float(np.dot(sym_to_vec(X), sym_to_vec(Y)))
    rhs = float(np.trace(X @ Y))
    assert lhs == pytest.approx(rhs, abs=1e-12)
    # also equals sym_inner
    assert sym_inner(X, Y) == pytest.approx(rhs, abs=1e-12)


@pytest.mark.parametrize("N", [1, 2, 4, 7])
def test_norm_matches_vector_norm(N):
    rng = np.random.default_rng(7 * N)
    X = random_symmetric(N, rng)
    assert sym_norm(X) == pytest.approx(np.linalg.norm(sym_to_vec(X)), abs=1e-12)
    assert sym_norm(X) == pytest.approx(np.linalg.norm(X, "fro"), abs=1e-12)


def test_normalize_flag():
    rng = np.random.default_rng(0)
    X = random_symmetric(6, rng, normalize=True)
    assert sym_norm(X) == pytest.approx(1.0, abs=1e-12)
    assert np.allclose(X, X.T)


@pytest.mark.parametrize("N", [1, 3, 5])
def test_pack_unpack_tangent(N):
    rng = np.random.default_rng(N)
    u = rng.standard_normal(N)
    X = random_symmetric(N, rng)
    w = pack_tangent(u, X)
    assert w.shape == (tangent_dim(N),)
    u2, X2 = unpack_tangent(w, N)
    assert np.allclose(u2, u, atol=1e-14)
    assert np.allclose(X2, X, atol=1e-14)
    # packed Euclidean norm^2 == ||u||^2 + ||X||_F^2
    assert float(np.dot(w, w)) == pytest.approx(
        float(u @ u) + sym_norm(X) ** 2, abs=1e-12
    )


def test_diagonal_stored_directly_offdiag_scaled():
    X = np.array([[2.0, 3.0], [3.0, 5.0]])
    v = sym_to_vec(X)
    # upper triangle row-major: (0,0), (0,1), (1,1)
    assert v[0] == pytest.approx(2.0)
    assert v[1] == pytest.approx(np.sqrt(2.0) * 3.0)
    assert v[2] == pytest.approx(5.0)
