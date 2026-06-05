"""GPU-side centering must reproduce the CPU CenteredPotential byte-for-byte.

The torch/CUDA low-dimensional high-M path moves the three heavy bank reductions
(``E[grad Phi]``, ``E[Hess Phi]``, eigenvalue extremes of ``Hess Phi - M``) onto
the device via :func:`torch_backend.build_centered_potential_gpu`, then feeds them
to :class:`CenteredPotential` through its ``precomputed_stats`` path. The result
must be the *identical* NumPy potential that ``build_potential(...)`` produces on
the same bank, to float64 round-off.

These run on torch CPU when no GPU is present (the centering math is
device-agnostic) and additionally on CUDA when available.
"""
import numpy as np
import pytest

torch = pytest.importorskip("torch")

from src.common.monte_carlo import gaussian_samples
from src.natural_gradient_local_rate.potentials import build_potential
from src.natural_gradient_local_rate import torch_backend as tb

FAMILIES = ["separable", "additive_index", "random_feature", "radial_tail"]


def _devices():
    devs = [torch.device("cpu")]
    if torch.cuda.is_available():
        devs.append(torch.device("cuda"))
    return devs


@pytest.mark.parametrize("device", _devices())
@pytest.mark.parametrize("fam", FAMILIES)
@pytest.mark.parametrize("N", [4, 8])
def test_gpu_centering_matches_cpu_potential(fam, N, device):
    M = 4096
    Z = gaussian_samples(N, M, seed=7, antithetic=True)
    point = {"family": fam, "N_theta": N, "kappa_target": 8.0, "seed": 0}
    cpu = build_potential(fam, N, 8.0, 0, Z_ref=Z, feature_multiplier=4,
                          safety_factor=2.0)
    gpu = tb.build_centered_potential_gpu(
        fam, point, {}, Z, device, torch.float64,
        chunk_size=1024, safety_factor=2.0, feature_multiplier=4)

    # centering reductions
    assert gpu.M == pytest.approx(cpu.M, abs=1e-12)
    assert gpu.b == pytest.approx(cpu.b, abs=1e-12)
    # centering algebra derived from them
    assert gpu.L_A == pytest.approx(cpu.L_A, abs=1e-12)
    assert gpu.rho == pytest.approx(cpu.rho, abs=1e-12)
    assert gpu.alpha_target == pytest.approx(cpu.alpha_target, abs=1e-12)
    assert gpu.beta_target == pytest.approx(cpu.beta_target, abs=1e-12)
    assert gpu.L_A_is_empirical == cpu.L_A_is_empirical
    # Hessian-eigenvalue extremes and centering diagnostics
    assert gpu.empirical_min_hess_eig == pytest.approx(cpu.empirical_min_hess_eig, abs=1e-11)
    assert gpu.empirical_max_hess_eig == pytest.approx(cpu.empirical_max_hess_eig, abs=1e-11)
    assert gpu.norm_mean_grad == pytest.approx(cpu.norm_mean_grad, abs=1e-11)

    # Hess V on fresh samples must agree
    Zt = gaussian_samples(N, 256, seed=99)
    assert gpu.batch_hess(Zt) == pytest.approx(cpu.batch_hess(Zt), abs=1e-11)

    # metadata that downstream CSV columns read
    mc, mg = cpu.metadata(), gpu.metadata()
    for k in ("alpha_target", "beta_target", "rho", "L_A", "kappa_target",
              "N_theta", "seed"):
        assert mg[k] == pytest.approx(mc[k], abs=1e-12), k


@pytest.mark.parametrize("fam", FAMILIES)
def test_centering_stats_chunking_invariant(fam):
    """The on-device stats are identical regardless of the chunk size used."""
    N, M = 6, 2048
    Z = gaussian_samples(N, M, seed=3, antithetic=True)
    from src.natural_gradient_local_rate.potentials import build_raw_feature
    raw = build_raw_feature(fam, N, 0, feature_multiplier=4)
    dev = torch.device("cpu")
    s_full = tb.centering_stats_on_device(raw, Z, dev, torch.float64, chunk_size=M)
    s_chunk = tb.centering_stats_on_device(raw, Z, dev, torch.float64, chunk_size=257)
    for a, b in zip(s_full[:2], s_chunk[:2]):       # b, M arrays
        assert a == pytest.approx(b, abs=1e-12)
    for a, b in zip(s_full[2:4], s_chunk[2:4]):     # emin, emax floats
        assert a == pytest.approx(b, abs=1e-12)
