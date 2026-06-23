"""End-to-end smoke tests for the STL scripts and the figure builders.

Runs the estimator and algorithm smoke grids (CPU) once into a temporary output
directory, checks the expected CSVs and columns, then exercises the plotting
script and the report figure builders so they cannot silently break.
"""
import os
import subprocess
import sys

import pandas as pd
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_SMOKE_CFG = os.path.join(_ROOT, "configs", "natural_gradient_stl_variance",
                          "stl_variance_smoke.yaml")
_EST = os.path.join(_ROOT, "scripts", "natural_gradient_stl_variance",
                    "run_estimator_variance.py")
_ALG = os.path.join(_ROOT, "scripts", "natural_gradient_stl_variance",
                    "run_algorithm_grid.py")
_PLOT = os.path.join(_ROOT, "scripts", "natural_gradient_stl_variance",
                     "plot_results.py")


def _run(script, outdir, extra):
    cmd = [sys.executable, script, "--config", _SMOKE_CFG,
           "--outdir", outdir, "--overwrite"] + extra
    proc = subprocess.run(cmd, cwd=_ROOT, capture_output=True, text=True)
    assert proc.returncode == 0, f"{script} failed:\n{proc.stdout}\n{proc.stderr}"
    return proc


@pytest.fixture(scope="module")
def smoke_outputs(tmp_path_factory):
    outdir = str(tmp_path_factory.mktemp("stl_smoke"))
    _run(_EST, outdir, ["--device", "cpu"])
    _run(_ALG, outdir, ["--device", "cpu"])
    return outdir


def test_estimator_csvs_and_columns(smoke_outputs):
    long_path = os.path.join(smoke_outputs, "estimator_variance.csv")
    summ_path = os.path.join(smoke_outputs, "estimator_variance_summary.csv")
    assert os.path.exists(long_path) and os.path.exists(summ_path)
    df = pd.read_csv(long_path)
    required = {"target_name", "kind", "d", "kappa", "tau", "state", "seed",
                "bias_base", "bias_stl", "var_fr_base", "var_fr_stl",
                "ratio_fr", "ratio_euclidean"}
    assert required.issubset(df.columns)
    # Six states for both targets, two seeds each in the smoke grid.
    assert df["state"].nunique() == 6
    # STL is essentially zero-variance at the Gaussian optimum.
    opt = df[(df.kind == "gaussian") & (df.state == "optimum")]
    assert opt["var_fr_stl"].max() < 1e-12
    assert os.path.exists(os.path.join(smoke_outputs, "estimator_target_metadata.json"))
    assert os.path.exists(os.path.join(smoke_outputs, "estimator_run_metadata.json"))


def test_algorithm_csvs_and_columns(smoke_outputs):
    long_path = os.path.join(smoke_outputs, "algorithm_results_long.csv")
    summ_path = os.path.join(smoke_outputs, "algorithm_summary.csv")
    tail_path = os.path.join(smoke_outputs, "tail_noise_floor.csv")
    for p in (long_path, summ_path, tail_path):
        assert os.path.exists(p), p
    long_df = pd.read_csv(long_path)
    long_req = {"target_name", "method", "scheme", "stl", "seed", "d", "kappa",
                "tau", "dt", "batch_size", "step", "energy_gap", "sq_mean_error",
                "rel_cov_fro_error", "w2_sq", "min_eig_C", "max_eig_C", "spd_fail"}
    assert long_req.issubset(long_df.columns)
    assert set(long_df["method"].unique()) == {
        "riemannian", "riemannian_stl", "kl", "kl_stl"}
    summ_df = pd.read_csv(summ_path)
    assert {"gap_final", "tail_median_gap", "iter_to_1e_minus_3"}.issubset(summ_df.columns)
    tail_df = pd.read_csv(tail_path)
    assert {"tail_mean_gap", "tail_median_gap", "tail_std_gap",
            "final_gap"}.issubset(tail_df.columns)


def test_plot_results_script(smoke_outputs):
    proc = subprocess.run(
        [sys.executable, _PLOT, "--outdir", smoke_outputs],
        cwd=_ROOT, capture_output=True, text=True)
    assert proc.returncode == 0, f"plot_results failed:\n{proc.stdout}\n{proc.stderr}"
    figs = os.path.join(smoke_outputs, "figures")
    pngs = [f for f in os.listdir(figs) if f.endswith(".png")]
    assert pngs, "no figures produced"


def test_report_figure_builders(smoke_outputs):
    """The shared figure builders run on smoke outputs without crashing."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.natural_gradient_stl_variance import plotting as P

    est = pd.read_csv(os.path.join(smoke_outputs, "estimator_variance.csv"))
    long_df = pd.read_csv(os.path.join(smoke_outputs, "algorithm_results_long.csv"))
    tail = pd.read_csv(os.path.join(smoke_outputs, "tail_noise_floor.csv"))

    figs = [
        P.fig_estimator_variance_ratio(est, "gaussian"),
        P.fig_variance_by_distance(est, "gaussian"),
        P.fig_algorithm_gap(long_df, "gaussian"),
        P.fig_noise_floor(tail),
        P.fig_batchsize_effect(tail),
    ]
    for fig in figs:
        if fig is not None:
            plt.close(fig)
