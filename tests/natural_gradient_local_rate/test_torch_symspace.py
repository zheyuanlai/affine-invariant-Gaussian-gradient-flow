"""Torch symmetric-space helpers must match the NumPy convention bit-for-bit."""
import numpy as np
import pytest

torch = pytest.importorskip("torch")

from src.common import symspace as ns
from src.common import torch_symspace as ts


@pytest.mark.parametrize("N", [1, 2, 3, 5, 8])
def test_sym_to_vec_matches_numpy(N):
    rng = np.random.default_rng(N)
    X = ns.random_symmetric(N, rng)
    v_np = ns.sym_to_vec(X)
    v_t = ts.torch_sym_to_vec(torch.as_tensor(X, dtype=torch.float64)).numpy()
    assert v_t == pytest.approx(v_np, abs=1e-12)
    assert ts.torch_sym_dim(N) == ns.sym_dim(N) == v_np.shape[0]


@pytest.mark.parametrize("N", [1, 2, 4, 7])
def test_vec_to_sym_matches_numpy(N):
    rng = np.random.default_rng(10 + N)
    v = rng.standard_normal(ns.sym_dim(N))
    X_np = ns.vec_to_sym(v, N)
    X_t = ts.torch_vec_to_sym(torch.as_tensor(v, dtype=torch.float64), N).numpy()
    assert X_t == pytest.approx(X_np, abs=1e-12)


@pytest.mark.parametrize("N", [2, 4, 6])
def test_roundtrip_and_isometry(N):
    rng = np.random.default_rng(99 + N)
    X = torch.as_tensor(ns.random_symmetric(N, rng), dtype=torch.float64)
    Y = torch.as_tensor(ns.random_symmetric(N, rng), dtype=torch.float64)
    # round-trip
    assert ts.torch_vec_to_sym(ts.torch_sym_to_vec(X), N).numpy() == pytest.approx(X.numpy(), abs=1e-12)
    # Frobenius isometry: dot(vec X, vec Y) == Tr(X Y)
    lhs = float(ts.torch_sym_to_vec(X) @ ts.torch_sym_to_vec(Y))
    rhs = float(torch.trace(X @ Y))
    assert lhs == pytest.approx(rhs, abs=1e-12)


@pytest.mark.parametrize("N", [2, 3, 5])
def test_batched_matches_loop(N):
    rng = np.random.default_rng(7 * N)
    Xb = np.stack([ns.random_symmetric(N, rng) for _ in range(6)])
    v_batch = ts.torch_sym_to_vec_batch(torch.as_tensor(Xb, dtype=torch.float64)).numpy()
    v_loop = np.stack([ns.sym_to_vec(Xb[i]) for i in range(6)])
    assert v_batch == pytest.approx(v_loop, abs=1e-12)


@pytest.mark.parametrize("N", [2, 4])
def test_fisher_rao_packing_matches_numpy(N):
    rng = np.random.default_rng(N)
    u = rng.standard_normal(N)
    X = ns.random_symmetric(N, rng)
    y_np = ns.pack_tangent_fr(u, X)
    y_t = ts.torch_pack_tangent_fr(torch.as_tensor(u), torch.as_tensor(X)).numpy()
    assert y_t == pytest.approx(y_np, abs=1e-12)
    u2, X2 = ts.torch_unpack_tangent_fr(torch.as_tensor(y_np), N)
    assert u2.numpy() == pytest.approx(u, abs=1e-12)
    assert X2.numpy() == pytest.approx(X, abs=1e-12)
    # packed Euclidean dot == Fisher--Rao inner product
    fr = float(ts.torch_fisher_rao_inner(torch.as_tensor(u), torch.as_tensor(X),
                                         torch.as_tensor(u), torch.as_tensor(X)))
    assert fr == pytest.approx(float(u @ u) + 0.5 * ns.sym_norm(X) ** 2, abs=1e-12)
