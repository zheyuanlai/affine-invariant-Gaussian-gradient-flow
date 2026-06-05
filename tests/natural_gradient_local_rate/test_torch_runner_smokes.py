"""Tiny torch-CPU smoke tests for the production runner entry points."""
import importlib.util
import os
import subprocess
import sys

import pandas as pd
import pytest

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("torch") is None,
    reason="torch not installed",
)

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
_SCRIPTS = os.path.join(_ROOT, "scripts", "natural_gradient_local_rate")

RUN_OPERATOR_GRID = os.path.join(_SCRIPTS, "run_operator_grid.py")
RUN_LINEARIZED_RATE_GRID = os.path.join(_SCRIPTS, "run_linearized_rate_grid.py")
RUN_OPERATOR_LINEARIZED_GRID = os.path.join(_SCRIPTS, "run_operator_linearized_grid.py")
RUN_SAMPLE_SIZE_SCALING = os.path.join(_SCRIPTS, "run_sample_size_scaling.py")


TINY_GRID_CONFIG = """\
experiment:
  name: tiny_torch_runner_grid
grid:
  N_theta: [3]
  kappa_target: [5]
  potential_family: ["radial_tail"]
  seeds: [0]
monte_carlo:
  M_mc: 128
  antithetic: true
  chunk_size: 32
operator:
  estimator: "symmetrized"
  diagonal_benchmark: true
  separable_exact_benchmark: true
  compute_gamma_loc: true
  backend: "numpy"
  device: "auto"
  dtype: "float32"
  eigensolver: "auto"
  explicit_dense_max_N_theta: 2
  basis_block_size: 99
  chunk_size: 32
  eigsh_tol: 1.0e-4
  eigsh_maxiter: 200
  quadrature_nodes: 20
outputs:
  base_dir: __BASE_DIR__
"""


TINY_SCALING_CONFIG = """\
experiment:
  name: tiny_torch_runner_scaling
grid:
  N_theta: [3]
  kappa_target: [5]
  potential_family: ["radial_tail"]
  seeds: [0]
monte_carlo:
  M_mc: [96, 128]
  antithetic: true
  chunk_size: 32
operator:
  estimator: "symmetrized"
  diagonal_benchmark: true
  separable_exact_benchmark: true
  compute_gamma_loc: true
  backend: "numpy"
  device: "auto"
  dtype: "float32"
  eigensolver: "auto"
  explicit_dense_max_N_theta: 2
  basis_block_size: 99
  chunk_size: 32
  eigsh_tol: 1.0e-4
  eigsh_maxiter: 200
  quadrature_nodes: 20
outputs:
  base_dir: __BASE_DIR__
"""


def _write_config(tmp_path, text, base_dir):
    cfg_path = tmp_path / "tiny_torch.yaml"
    cfg_path.write_text(text.replace("__BASE_DIR__", str(base_dir)))
    return cfg_path


def _run(script, cfg_path):
    cmd = [
        sys.executable, script,
        "--config", str(cfg_path),
        "--backend", "torch",
        "--device", "cpu",
        "--dtype", "float64",
        "--chunk-size", "64",
        "--basis-block-size", "2",
        "--explicit-dense-max-N-theta", "4",
        "--overwrite",
    ]
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("VECLIB_MAXIMUM_THREADS", "1")
    return subprocess.run(cmd, cwd=_ROOT, capture_output=True, text=True,
                          timeout=600, env=env)


def _assert_torch_cpu_csv(csv_path, expected_rows):
    assert csv_path.exists(), "results_long.csv was not written"
    df = pd.read_csv(csv_path)
    assert len(df) == expected_rows
    assert (df["status"] == "ok").all(), df[["status", "error_message"]].to_dict("records")
    assert (df["backend"] == "torch").all()
    assert (df["device"] == "cpu").all()
    assert (df["dtype"] == "float64").all()
    assert (df["chunk_size"] == 64).all()
    assert (df["basis_block_size"] == 2).all()
    assert (df["explicit_dense_max_N_theta"] == 4).all()
    assert (df["potential_family"] == "radial_tail").all()
    assert df["operator_matrix_dim"].notna().all()
    assert df["dense_matrix_memory_mb"].notna().all()
    assert df["Lambda_hat_full_sym"].notna().all()
    assert df["gamma_loc"].notna().all()
    return df


def test_run_operator_grid_torch_cpu_smoke(tmp_path):
    base_dir = tmp_path / "out"
    cfg_path = _write_config(tmp_path, TINY_GRID_CONFIG, base_dir)
    res = _run(RUN_OPERATOR_GRID, cfg_path)
    assert res.returncode == 0, f"script failed:\nSTDOUT\n{res.stdout}\nSTDERR\n{res.stderr}"
    _assert_torch_cpu_csv(base_dir / "operator_grid" / "results_long.csv", 1)


def test_run_linearized_rate_grid_torch_cpu_smoke(tmp_path):
    base_dir = tmp_path / "out"
    cfg_path = _write_config(tmp_path, TINY_GRID_CONFIG, base_dir)
    res = _run(RUN_LINEARIZED_RATE_GRID, cfg_path)
    assert res.returncode == 0, f"script failed:\nSTDOUT\n{res.stdout}\nSTDERR\n{res.stderr}"
    _assert_torch_cpu_csv(base_dir / "linearized_rate_grid" / "results_long.csv", 1)
    eig = base_dir / "linearized_rate_grid" / "eigenvectors"
    assert any(eig.glob("*.npz")), "linearized-rate smoke did not write an eigenvector"


def test_run_operator_linearized_grid_torch_cpu_smoke(tmp_path):
    base_dir = tmp_path / "out"
    cfg_path = _write_config(tmp_path, TINY_GRID_CONFIG, base_dir)
    res = _run(RUN_OPERATOR_LINEARIZED_GRID, cfg_path)
    assert res.returncode == 0, f"script failed:\nSTDOUT\n{res.stdout}\nSTDERR\n{res.stderr}"
    op = _assert_torch_cpu_csv(base_dir / "operator_grid" / "results_long.csv", 1)
    lr = _assert_torch_cpu_csv(base_dir / "linearized_rate_grid" / "results_long.csv", 1)
    assert op["gamma_loc"].iloc[0] == pytest.approx(lr["gamma_loc"].iloc[0])
    assert op["Lambda_hat_full_sym"].iloc[0] == pytest.approx(lr["Lambda_hat_full_sym"].iloc[0])
    assert (base_dir / "operator_grid" / "summary.csv").exists()
    assert (base_dir / "linearized_rate_grid" / "summary.csv").exists()
    eig = base_dir / "linearized_rate_grid" / "eigenvectors"
    assert any(eig.glob("*.npz")), "joint smoke did not write an eigenvector"


def test_run_sample_size_scaling_torch_cpu_smoke(tmp_path):
    base_dir = tmp_path / "out"
    cfg_path = _write_config(tmp_path, TINY_SCALING_CONFIG, base_dir)
    res = _run(RUN_SAMPLE_SIZE_SCALING, cfg_path)
    assert res.returncode == 0, f"script failed:\nSTDOUT\n{res.stdout}\nSTDERR\n{res.stderr}"
    df = _assert_torch_cpu_csv(base_dir / "sample_size_scaling" / "results_long.csv", 2)
    assert sorted(df["M_mc"].tolist()) == [96, 128]
