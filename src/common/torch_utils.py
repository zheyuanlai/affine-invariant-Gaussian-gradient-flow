"""Optional-PyTorch plumbing: lazy import, device/dtype resolution, info.

The repository runs fine without PyTorch installed; nothing here imports ``torch``
at module load. GPU code paths call :func:`get_torch` (which raises a clear error
if torch is missing) and :func:`resolve_device` (which raises if CUDA is
requested but unavailable -- it never silently falls back to CPU).
"""
from __future__ import annotations

import importlib


def torch_available():
    """True iff ``import torch`` succeeds (does not import it into this module)."""
    try:
        importlib.import_module("torch")
        return True
    except Exception:
        return False


def get_torch():
    """Return the imported ``torch`` module or raise a clear ImportError."""
    try:
        return importlib.import_module("torch")
    except Exception as exc:  # pragma: no cover - exercised only without torch
        raise ImportError(
            "PyTorch is required for the torch/GPU backend but is not importable. "
            "Install it (see requirements-gpu.txt and the README 'GPU backend' "
            "section), e.g. `pip install torch`, or use backend='numpy'."
        ) from exc


def resolve_device(device="auto"):
    """Resolve a device string to a ``torch.device``.

    * ``"auto"`` -> ``cuda`` if available, else ``cpu``;
    * ``"cuda"`` (or ``"cuda:N"``) -> error if CUDA is unavailable;
    * ``"cpu"``  -> cpu.

    Never silently downgrades a CUDA request to CPU.
    """
    torch = get_torch()
    dev = str(device).lower()
    if dev == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if dev == "cpu":
        return torch.device("cpu")
    if dev == "cuda" or dev.startswith("cuda:"):
        if not torch.cuda.is_available():
            raise RuntimeError(
                f"device={device!r} was requested but CUDA is not available "
                "(torch.cuda.is_available() is False). Use device='cpu' or "
                "device='auto', or run on a CUDA machine."
            )
        return torch.device(dev)
    raise ValueError(f"unknown device {device!r} (expected 'auto', 'cpu', or 'cuda')")


def resolve_backend(backend="numpy", device="auto"):
    """Resolve ``backend`` (``numpy`` / ``torch`` / ``auto``) to ``"numpy"`` or ``"torch"``.

    ``"auto"`` selects torch only when it is importable *and* the requested
    device resolves to CUDA (torch-on-CPU offers no advantage over NumPy);
    otherwise NumPy. ``"numpy"`` returns immediately without touching torch.
    """
    b = str(backend).lower()
    if b == "numpy":
        return "numpy"
    if b == "torch":
        return "torch"
    if b == "auto":
        if not torch_available():
            return "numpy"
        try:
            return "torch" if resolve_device(device).type == "cuda" else "numpy"
        except Exception:
            return "numpy"
    raise ValueError(f"unknown backend {backend!r} (expected numpy/torch/auto)")


def resolve_dtype(dtype="float64"):
    """Resolve a dtype string to a ``torch`` dtype (``float64`` or ``float32``)."""
    torch = get_torch()
    d = str(dtype).lower()
    if d in ("float64", "double", "f64"):
        return torch.float64
    if d in ("float32", "float", "f32"):
        return torch.float32
    raise ValueError(f"unknown dtype {dtype!r} (expected 'float64' or 'float32')")


def torch_device_info(device="auto", dtype="float64"):
    """Return a JSON-serializable dict describing the torch/CUDA environment.

    Safe to call without torch installed (returns ``torch_available=False`` and
    ``None``/``False`` fields). Resolution errors (e.g. cuda requested but
    unavailable) are reported in the dict rather than raised, so this is safe for
    CSV bookkeeping.
    """
    info = {
        "torch_available": False,
        "torch_version": None,
        "cuda_available": False,
        "cuda_version": None,
        "device": None,
        "device_name": None,
        "dtype": str(dtype),
        "gpu_available": False,
    }
    if not torch_available():
        return info
    torch = get_torch()
    info["torch_available"] = True
    info["torch_version"] = str(torch.__version__)
    info["cuda_available"] = bool(torch.cuda.is_available())
    info["gpu_available"] = bool(torch.cuda.is_available())
    info["cuda_version"] = getattr(torch.version, "cuda", None)
    try:
        dev = resolve_device(device)
        info["device"] = str(dev)
        if dev.type == "cuda":
            info["device_name"] = torch.cuda.get_device_name(dev)
        else:
            info["device_name"] = "cpu"
    except Exception as exc:
        info["device"] = f"unresolved({device})"
        info["device_name"] = f"error: {exc}"
    return info
