"""Root conftest: add project root to sys.path so 'src' is importable."""
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))


def pytest_configure(config):
    """Register custom markers (keeps ``-W error`` / strict-markers clean)."""
    config.addinivalue_line(
        "markers",
        "gpu: test requires a CUDA device; skipped unless torch.cuda is available.",
    )
