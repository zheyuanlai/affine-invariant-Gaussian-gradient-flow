"""Torch availability plumbing. Skips entirely if torch is not installed."""
import pytest

torch = pytest.importorskip("torch")

from src.common import torch_utils


def test_torch_available_true_here():
    assert torch_utils.torch_available() is True
    assert torch_utils.get_torch() is torch


def test_resolve_device_cpu():
    dev = torch_utils.resolve_device("cpu")
    assert dev.type == "cpu"


def test_resolve_device_auto_is_valid():
    dev = torch_utils.resolve_device("auto")
    assert dev.type in ("cpu", "cuda")


def test_resolve_device_cuda_errors_when_unavailable():
    if torch.cuda.is_available():
        pytest.skip("CUDA available; the error path is not exercised here")
    with pytest.raises(RuntimeError):
        torch_utils.resolve_device("cuda")


def test_resolve_dtype():
    assert torch_utils.resolve_dtype("float64") is torch.float64
    assert torch_utils.resolve_dtype("float32") is torch.float32
    with pytest.raises(ValueError):
        torch_utils.resolve_dtype("float16")


def test_resolve_backend():
    assert torch_utils.resolve_backend("numpy") == "numpy"
    assert torch_utils.resolve_backend("torch") == "torch"
    # auto picks torch only when a CUDA device resolves
    expected = "torch" if torch.cuda.is_available() else "numpy"
    assert torch_utils.resolve_backend("auto", "auto") == expected


def test_device_info_keys():
    info = torch_utils.torch_device_info("cpu", "float64")
    for k in ("torch_available", "torch_version", "cuda_available", "cuda_version",
              "device", "device_name", "dtype", "gpu_available"):
        assert k in info
    assert info["torch_available"] is True
