# Gaussian Gradient Flow Experiments

Numerical experiments for Gaussian gradient flows in variational inference. The
variational family is the non-degenerate Gaussians `N(m, C)` and the goal is to
minimize `KL(N(m, C) || target)`. The repository is organized into two
self-contained experiment groups, each with its own configs, source modules,
scripts, tests, and outputs.

## Experiment groups

### `omega_tau_modes`
The two-parameter affine-invariant Gaussian gradient flow family with parameters
`omega` and `tau`. Studies how `omega` and `tau` affect the mean,
covariance-volume, covariance-shape, and mixed convergence modes, for both an
exact Gaussian target and a strongly log-concave non-Gaussian target.

### `natural_gradient_local_rate`
The Gaussian natural gradient flow. Asks whether its **local** convergence rate
near equilibrium is essentially dimension-free — i.e. whether the rate depends on
the dimension `N_theta`, or only weakly on the conditioning `kappa` through
`log(kappa)`. Working in equilibrium-whitened coordinates with `a_star = (0, I)`,
it estimates the local rate `gamma_loc` (smallest eigenvalue of the linearized
positive generator `L_star`) and the operator norm `Lambda_hat` (largest
eigenvalue of the symmetric-matrix operator `H_lin`), and validates them against a
Riemannian-exponential-map flow simulation.

Notation used throughout this group: `N_theta` (dimension),
`rho_post(theta) ∝ exp(-V(theta))` (target), `rho_a = N(m, C)`, `a = (m, C)`,
`E(a) = KL(rho_a || rho_post)`, equilibrium `a_star = (0, I)`.

## Repository structure

```
configs/
  omega_tau_modes/                gaussian_target / logconcave_target configs
  natural_gradient_local_rate/    smoke + grid + production configs
src/
  common/                         spd, symspace, monte_carlo, io, plotting style
  omega_tau_modes/                omega/tau dynamics, targets, metrics, plotting
  natural_gradient_local_rate/    potentials, operators, linearized rate, flow
scripts/
  omega_tau_modes/                grid runners + plotting (+ back-compat shims in scripts/)
  natural_gradient_local_rate/    operator/rate/flow runners, run_all, plotting
tests/
  common/  omega_tau_modes/  natural_gradient_local_rate/
docs/
  references/                     local-only manuscript PDFs (git-ignored)
  specs/                          implementation specs (tracked)
reports/                          detailed mathematical write-ups
outputs/                          experiment outputs (git-ignored)
```

## Installation

```bash
pip install -r requirements.txt
```

CPU only, float64, no GPU/PyTorch. Python 3.9+.

## Running tests

```bash
pytest
```

## Running the omega/tau experiments

```bash
# fast smoke runs
python scripts/omega_tau_modes/run_gaussian_grid.py   --config configs/omega_tau_modes/gaussian_target.yaml   --smoke
python scripts/omega_tau_modes/run_logconcave_grid.py --config configs/omega_tau_modes/logconcave_target.yaml --smoke

# full grids (drop --smoke), then figures
python scripts/omega_tau_modes/run_gaussian_grid.py   --config configs/omega_tau_modes/gaussian_target.yaml
python scripts/omega_tau_modes/plot_gaussian_results.py   --indir outputs/omega_tau_modes/gaussian_grid
python scripts/omega_tau_modes/run_logconcave_grid.py --config configs/omega_tau_modes/logconcave_target.yaml
python scripts/omega_tau_modes/plot_logconcave_results.py --indir outputs/omega_tau_modes/logconcave_grid
```

The old entry points (`python scripts/run_gaussian_grid.py`, etc.) still work as
thin back-compat shims that delegate to the `scripts/omega_tau_modes/` versions.

## Running the natural-gradient local-rate experiments

```bash
# smoke (fast) — operator norm, linearized rate, flow validation
python scripts/natural_gradient_local_rate/run_operator_grid.py        --config configs/natural_gradient_local_rate/smoke.yaml
python scripts/natural_gradient_local_rate/run_linearized_rate_grid.py --config configs/natural_gradient_local_rate/smoke.yaml
python scripts/natural_gradient_local_rate/run_flow_validation.py      --config configs/natural_gradient_local_rate/smoke.yaml

# or all three in sequence
python scripts/natural_gradient_local_rate/run_all.py --smoke

# sample-size scaling (sweeps M_mc to separate real trends from MC noise)
python scripts/natural_gradient_local_rate/run_sample_size_scaling.py --config configs/natural_gradient_local_rate/sample_size_scaling.yaml

# figures
python scripts/natural_gradient_local_rate/plot_results.py              --input outputs/natural_gradient_local_rate   --outdir outputs/natural_gradient_local_rate/figures
python scripts/natural_gradient_local_rate/plot_estimator_diagnostics.py --input outputs/natural_gradient_local_rate   --outdir outputs/natural_gradient_local_rate/figures/estimator_diagnostics
```

The larger `operator_grid.yaml`, `linearized_rate_grid.yaml`,
`flow_validation.yaml`, and `production_all.yaml` configs are provided but are
**not** run automatically — they are expensive.

#### Operator estimator: modes and diagnostics

The operator norm has two estimator modes, selected by the `operator.estimator`
config key:

- **`symmetrized`** (default) — the self-adjoint `H_sym = ½(H_lin + H_lin*)`.
  **Use this for all eigenvalue computations**, including the `L_star`
  covariance block (`gamma_loc`), which is solved in Fisher–Rao-whitened
  coordinates so the matrix `eigsh` sees is genuinely symmetric.
- **`raw_forward`** — the uncorrected forward `H_lin`, kept only for diagnostics
  and backward comparison; it is *not* self-adjoint on a finite sample bank.

Before interpreting any dimension dependence, run the **separable sanity
checks**: the diagonal-restricted estimator and the Gauss–Hermite
`separable_exact` benchmark are dimension-free ground truths for separable
controls, and the **sample-size scaling** sweep shows whether a high-dimensional
trend shrinks as `M_mc` grows (finite-sample spectral noise) or persists (a real
effect). See [`reports/natural_gradient_local_rate_notes.md`](reports/natural_gradient_local_rate_notes.md).

#### GPU backend (optional PyTorch)

An optional PyTorch backend runs the same corrected estimators on a CUDA device
(for Colab Pro / A100); the NumPy/SciPy CPU path remains the default and the repo
works without torch installed.

- Select it with `operator.backend: torch` (or `--backend torch`) and
  `operator.device: cuda` (or `--device cuda`); `backend: auto` uses torch only
  when a CUDA device is available, else NumPy.
- The GPU path uses a dense `torch.linalg.eigh` eigensolver for
  `N_theta <= explicit_dense_max_N_theta` (default 64) and is numerically
  identical to the CPU path on the same bank.
- `device=cuda` raises a clear error if CUDA is unavailable (no silent CPU
  fallback); use `device=cpu` to exercise the torch path on a CPU-only machine.
- Install torch separately (it is not in base `requirements.txt`); see
  [`requirements-gpu.txt`](requirements-gpu.txt) and
  [`reports/gpu_colab_notes.md`](reports/gpu_colab_notes.md).

```bash
# GPU smoke (torch dense path) — runs on cuda, or cpu without a GPU
python scripts/natural_gradient_local_rate/run_sample_size_scaling.py \
    --config configs/natural_gradient_local_rate/gpu_smoke.yaml --backend torch --device cuda
```

For production diagnostics, prefer
`configs/natural_gradient_local_rate/gpu_baseline_scaling.yaml`, which runs
Gaussian, separable, and random-feature families side by side. Always include
Gaussian and separable baselines, run
`scripts/natural_gradient_local_rate/postprocess_baselines.py` before
interpreting plots, and do not interpret raw full-operator growth with
`N_theta` without baseline correction and sample-size convergence checks. The GPU
backend changes only the speed, not the meaning, of the estimates.

## Outputs

```
outputs/omega_tau_modes/gaussian_grid/        results_long.csv, summary.csv, figures/
outputs/omega_tau_modes/logconcave_grid/      results_long.csv, summary.csv, reference optimum, figures/
outputs/natural_gradient_local_rate/operator_grid/       results_long.csv, summary.csv
outputs/natural_gradient_local_rate/linearized_rate_grid/ results_long.csv, summary.csv, eigenvectors
outputs/natural_gradient_local_rate/flow_validation/     trajectories + summaries
outputs/natural_gradient_local_rate/figures/             PNG + PDF figures
```

The legacy directories `outputs/gaussian_grid/` and `outputs/logconcave_grid/`
from earlier runs are left in place for reference. `outputs/` is git-ignored.

## Reports and local references

Detailed mathematical write-ups live in [`reports/`](reports/):

- [`reports/omega_tau_modes_notes.md`](reports/omega_tau_modes_notes.md) — the
  full derivation, diagnostics, and result tables for the omega/tau experiments
  (the former root README).
- [`reports/natural_gradient_local_rate_notes.md`](reports/natural_gradient_local_rate_notes.md)
  — notes for the natural-gradient local-rate experiment.

Local manuscript PDFs in `docs/references/` are **not** committed (git-ignored);
the implementation spec in `docs/specs/` is the tracked source of truth.
