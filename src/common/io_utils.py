"""YAML / CSV / JSON I/O helpers and run-id utilities."""
from __future__ import annotations

import datetime
import json
import os

import numpy as np
import pandas as pd
import yaml


def load_yaml(path):
    """Load a YAML file into a Python object."""
    with open(path, "r") as f:
        return yaml.safe_load(f)


def ensure_dir(path):
    """Create ``path`` (a directory) if needed and return it."""
    if path:
        os.makedirs(path, exist_ok=True)
    return path


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def save_json(path, obj):
    """Write ``obj`` to ``path`` as pretty JSON (numpy-aware)."""
    ensure_dir(os.path.dirname(path))
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=_json_default)
    return path


def save_csv_append(path, row, columns=None):
    """Append ``row`` (a dict, list of dicts, or DataFrame) to a CSV.

    Writes the header only when the file is new/empty. When ``columns`` is given
    the frame is reindexed to that exact column order, which keeps appended rows
    aligned across calls.
    """
    ensure_dir(os.path.dirname(path))
    if isinstance(row, pd.DataFrame):
        df = row
    elif isinstance(row, dict):
        df = pd.DataFrame([row])
    else:
        df = pd.DataFrame(list(row))
    if columns is not None:
        df = df.reindex(columns=columns)
    write_header = (not os.path.exists(path)) or os.path.getsize(path) == 0
    df.to_csv(path, mode="a", header=write_header, index=False)
    return path


def save_dataframe(path, df, columns=None):
    """Write a DataFrame to CSV (overwrite), optionally reordering columns."""
    ensure_dir(os.path.dirname(path))
    if columns is not None:
        df = df.reindex(columns=columns)
    df.to_csv(path, index=False)
    return path


def timestamp():
    """Compact UTC-naive local timestamp ``YYYYmmdd_HHMMSS``."""
    return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")


def run_id(prefix="run"):
    """A human-readable run id, e.g. ``run_20260603_171530``."""
    return f"{prefix}_{timestamp()}"
