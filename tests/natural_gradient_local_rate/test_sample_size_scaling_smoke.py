"""Smoke test for the sample-size-scaling runner script.

Runs the real entry point on a tiny config and checks that the results CSV is
produced with the full required column schema and clean status.
"""
import os
import sys
import subprocess

import pandas as pd
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
_SCRIPTS = os.path.join(_ROOT, "scripts", "natural_gradient_local_rate")
sys.path.insert(0, _SCRIPTS)
import _common  # noqa: E402

SCRIPT = os.path.join(_SCRIPTS, "run_sample_size_scaling.py")

TINY_CONFIG = """\
experiment:
  name: tiny_sample_size_scaling
grid:
  N_theta: [4]
  kappa_target: [5]
  potential_family: ["separable"]
  seeds: [0]
monte_carlo:
  M_mc: [256, 512]
  antithetic: true
  chunk_size: 256
operator:
  estimator: "symmetrized"
  diagonal_benchmark: true
  separable_exact_benchmark: true
  compute_gamma_loc: true
  backend: "numpy"
  chunk_size: 256
  eigsh_tol: 1.0e-4
  eigsh_maxiter: 300
  quadrature_nodes: 40
outputs:
  base_dir: {base_dir}
"""


def test_sample_size_scaling_smoke(tmp_path):
    base_dir = tmp_path / "out"
    cfg_path = tmp_path / "tiny.yaml"
    cfg_path.write_text(TINY_CONFIG.format(base_dir=str(base_dir)))

    res = subprocess.run(
        [sys.executable, SCRIPT, "--config", str(cfg_path), "--overwrite"],
        cwd=_ROOT, capture_output=True, text=True, timeout=600,
    )
    assert res.returncode == 0, f"script failed:\nSTDOUT\n{res.stdout}\nSTDERR\n{res.stderr}"

    csv = base_dir / "sample_size_scaling" / "results_long.csv"
    assert csv.exists(), "results_long.csv was not written"
    df = pd.read_csv(csv)

    # two M_mc values x one base grid point -> two rows
    assert len(df) == 2
    assert sorted(df["M_mc"].tolist()) == [256, 512]

    # every required column is present, in schema order
    missing = [c for c in _common.REQUIRED_COLUMNS if c not in df.columns]
    assert not missing, f"missing required columns: {missing}"
    assert list(df.columns[: len(_common.REQUIRED_COLUMNS)]) == _common.REQUIRED_COLUMNS

    assert (df["status"] == "ok").all(), df[["M_mc", "status", "error_message"]].to_dict("records")
    # separable -> the exact benchmark is populated (finite)
    assert df["Lambda_hat_separable_exact"].notna().all()
