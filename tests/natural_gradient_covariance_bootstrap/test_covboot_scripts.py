"""End-to-end smoke test for the covariance-bootstrap runner and figure builders.

Runs the smoke grid (CPU) once into a temporary output directory, checks the
required CSVs / metadata and their columns, then exercises the plotting script and
the shared figure builders so they cannot silently break.
"""
import json
import os
import subprocess
import sys

import pandas as pd
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_SMOKE_CFG = os.path.join(_ROOT, "configs", "natural_gradient_covariance_bootstrap",
                          "covariance_bootstrap_smoke.yaml")
_RUN = os.path.join(_ROOT, "scripts", "natural_gradient_covariance_bootstrap",
                    "run_experiments.py")
_PLOT = os.path.join(_ROOT, "scripts", "natural_gradient_covariance_bootstrap",
                     "plot_results.py")

REQUIRED_CSVS = [
    "results_long.csv", "covariance_bootstrap_summary.csv",
    "contraction_benchmark.csv", "wasserstein_bootstrap_summary.csv",
    "stl_floor_summary.csv",
]
REQUIRED_JSON = ["target_metadata.json", "run_metadata.json"]


@pytest.fixture(scope="module")
def smoke_outputs(tmp_path_factory):
    outdir = str(tmp_path_factory.mktemp("covboot_smoke"))
    cmd = [sys.executable, _RUN, "--config", _SMOKE_CFG, "--outdir", outdir,
           "--device", "cpu", "--overwrite"]
    proc = subprocess.run(cmd, cwd=_ROOT, capture_output=True, text=True)
    assert proc.returncode == 0, f"run_experiments failed:\n{proc.stdout}\n{proc.stderr}"
    return outdir


def test_required_outputs_exist(smoke_outputs):
    for name in REQUIRED_CSVS + REQUIRED_JSON:
        assert os.path.exists(os.path.join(smoke_outputs, name)), name


def test_results_long_columns(smoke_outputs):
    df = pd.read_csv(os.path.join(smoke_outputs, "results_long.csv"))
    required = {"experiment", "target", "scheme", "method", "kappa", "beta",
                "lambda0", "dt", "iteration", "energy_gap", "lambda_min_C", "L_n",
                "q_dynamic", "q_frozen", "status"}
    assert required.issubset(df.columns)
    assert set(df.experiment.unique()) == {
        "covariance_bootstrap_scaling", "dynamic_contraction_benchmark",
        "wasserstein_bootstrap_then_fr"}
    assert (df.status == "ok").all()


def test_scaling_summary(smoke_outputs):
    df = pd.read_csv(os.path.join(smoke_outputs, "covariance_bootstrap_summary.csv"))
    assert {"N_cov", "log_inv_beta_lambda0", "inv_lambda0",
            "N_cov_theory_slope"}.issubset(df.columns)
    # smaller lambda0 -> larger burn-in (monotone in log(1/(beta lambda0))).
    one = df[(df.scheme == "kl") & (df.kappa == df.kappa.max())].sort_values("lambda0")
    ncov = one.N_cov.values
    assert ncov[0] >= ncov[-1]         # smallest lambda0 first -> largest N_cov


def test_stl_floor_summary(smoke_outputs):
    df = pd.read_csv(os.path.join(smoke_outputs, "stl_floor_summary.csv"))
    assert {"target", "method", "dt", "tail_median_gap", "psi_star",
            "dt_psi_star"}.issubset(df.columns)
    # Gaussian target: Psi is exactly zero.
    g = df[df.target == "gaussian"]
    assert (g.psi_star == 0.0).all()
    # Gaussian STL floor is far below the Gaussian raw floor.
    for dt in sorted(g.dt.unique()):
        raw = g[(g.method == "kl_raw") & (g.dt == dt)].tail_median_gap.iloc[0]
        stl = g[(g.method == "kl_stl") & (g.dt == dt)].tail_median_gap.iloc[0]
        assert stl < raw
    # log-cosh: Psi is positive.
    lc = df[df.target == "log_cosh"]
    assert (lc.psi_star > 0).all()


def test_metadata(smoke_outputs):
    with open(os.path.join(smoke_outputs, "run_metadata.json")) as fh:
        meta = json.load(fh)
    assert meta["experiment_group"] == "natural_gradient_covariance_bootstrap"
    assert "cuda_visible_devices" in meta
    assert meta["deterministic_backend"] == "numpy"
    with open(os.path.join(smoke_outputs, "target_metadata.json")) as fh:
        tmeta = json.load(fh)
    assert tmeta["targets"]


def test_plot_results_script(smoke_outputs):
    proc = subprocess.run([sys.executable, _PLOT, "--outdir", smoke_outputs],
                          cwd=_ROOT, capture_output=True, text=True)
    assert proc.returncode == 0, f"plot_results failed:\n{proc.stdout}\n{proc.stderr}"
    figs = os.path.join(smoke_outputs, "figures")
    pngs = [f for f in os.listdir(figs) if f.endswith(".png")]
    assert pngs, "no figures produced"


def test_report_figure_builders(smoke_outputs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from src.natural_gradient_covariance_bootstrap import plotting as P

    long_df = pd.read_csv(os.path.join(smoke_outputs, "results_long.csv"))
    scaling = pd.read_csv(os.path.join(smoke_outputs, "covariance_bootstrap_summary.csv"))
    bench = pd.read_csv(os.path.join(smoke_outputs, "contraction_benchmark.csv"))
    wboot = pd.read_csv(os.path.join(smoke_outputs, "wasserstein_bootstrap_summary.csv"))
    floor = pd.read_csv(os.path.join(smoke_outputs, "stl_floor_summary.csv"))
    figs = [
        P.fig_covariance_envelope(long_df, scheme="kl"),
        P.fig_warmup_scaling(scaling),
        P.fig_dynamic_contraction(bench),
        P.fig_contraction_factors(bench),
        P.fig_wasserstein_bootstrap(wboot),
        P.fig_three_stage(long_df, scheme="kl"),
        P.fig_stl_noise_floor(floor),
    ]
    for fig in figs:
        if fig is not None:
            plt.close(fig)
