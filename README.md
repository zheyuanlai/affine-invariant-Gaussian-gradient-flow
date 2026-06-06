# Affine-Invariant Gaussian Gradient Flow Experiments

Numerical experiments for Gaussian variational inference, where the variational
family is the non-degenerate Gaussians `N(m, C)` and the objective is
`KL(N(m, C) || target)`. The repository contains two self-contained experiment
groups, each with its own configs, source modules, scripts, tests, and outputs.

The polished write-ups for both groups are the LaTeX reports in
[`reports/`](reports/); this file documents the code, the final outputs, and the
exact commands to reproduce them.

## Experiment groups

### 1. `omega_tau_modes` — affine-invariant `(omega, tau)` flows

A two-parameter family of affine-invariant Gaussian gradient flows. The scalar
`omega > 0` rescales the covariance dynamics uniformly; `tau > -omega/n` enters
only through a trace-weighting term and therefore acts solely on the
covariance-volume (trace) mode. The experiments isolate four convergence modes
(mean, covariance-volume, covariance-shape, mixed) and quantify the effect of
`omega` and `tau` on each, for an exact Gaussian target and a strongly
log-concave non-Gaussian target.

**Finding.** A negative `tau` accelerates volume-dominated transients (about a
2x speedup at `tau = -omega/2n`) and slows them by ~3/2 at `tau = +omega/2n`,
while the mean and shape modes are unaffected; the realized speedup tracks the
initial trace-dominance of the perturbation. See
[`reports/affine_invariant_omega_tau_report.tex`](reports/affine_invariant_omega_tau_report.tex).

### 2. `natural_gradient_local_rate` — local convergence rate

The Gaussian natural gradient flow near equilibrium. In equilibrium-whitened
coordinates with `a_star = (0, I)`, the local rate `gamma_loc` is the smallest
eigenvalue of the linearized positive generator `L_star` in the Fisher–Rao
metric. The question is whether `gamma_loc` genuinely depends on the dimension
`N_theta`, or only on the conditioning `kappa` through `log(kappa)`.

**Finding.** Over the final production grid the measured `gamma_loc` is
essentially flat in `N_theta` at every `kappa` and varies only with `kappa`,
indicating that the `N_theta` factor in the current proof bound is a proof
artifact rather than a property of the flow. See
[`reports/natural_gradient_local_rate_report.tex`](reports/natural_gradient_local_rate_report.tex).

## Which outputs are final

Only these directories are interpreted as evidence in the reports:

```
outputs/gaussian_grid/                              omega/tau, Gaussian target
outputs/logconcave_grid/                            omega/tau, log-concave target
outputs/natural_gradient_local_rate/operator_grid/        local rate: Lambda_hat + gamma_loc
outputs/natural_gradient_local_rate/linearized_rate_grid/ local rate: gamma_loc + eigenvectors
```

The local-rate final run is a single GPU production grid (`N_theta = 1..16`,
`kappa in {2,5,10,20,50,100}`, five potential families, three seeds,
`M_mc = 2^22 = 4,194,304`, torch/CUDA/float64, 1440 rows, all `status == ok`).
Exploratory smoke, baseline, and high-dimensional pilot runs are **not** evidence
and have been removed. `outputs/` is git-ignored except for the committed final
CSVs/summaries.

## Installation

```bash
pip install -r requirements.txt        # CPU, float64; NumPy/SciPy/matplotlib/pandas
```

The optional PyTorch GPU backend (used only for the local-rate production grid)
is not in the base requirements; see [`requirements-gpu.txt`](requirements-gpu.txt).
A CUDA build is needed for the production command below; a CPU torch build also
exercises the same code path (`--backend torch --device cpu`).

## Tests

```bash
pytest
```

## Reproducing the omega/tau experiments

The final outputs live at `outputs/gaussian_grid/` and
`outputs/logconcave_grid/`. Pass `--outdir` explicitly to write there (the config
default base directory is `outputs/omega_tau_modes/...`):

```bash
python scripts/omega_tau_modes/run_gaussian_grid.py \
    --config configs/omega_tau_modes/gaussian_target.yaml \
    --outdir outputs/gaussian_grid

python scripts/omega_tau_modes/run_logconcave_grid.py \
    --config configs/omega_tau_modes/logconcave_target.yaml \
    --outdir outputs/logconcave_grid
```

Add `--smoke` for fast reduced grids. Per-group figures (not required by the
reports) can be produced with
`scripts/omega_tau_modes/plot_gaussian_results.py` and
`plot_logconcave_results.py`.

## Reproducing the natural-gradient local-rate production run

The production grid runs on a CUDA GPU (developed on an NVIDIA H200). The joint
runner computes the shared dense accumulation once per grid point and writes both
the `operator_grid` and `linearized_rate_grid` stages plus the slow eigenvectors.
With `--outdir outputs/natural_gradient_local_rate` the two stages land directly
under that directory (the layout the reports read):

```bash
python scripts/natural_gradient_local_rate/run_operator_linearized_grid.py \
    --config configs/natural_gradient_local_rate/gpu_lowdim_operator_full.yaml \
    --backend torch --device cuda --dtype float64 \
    --chunk-size 1048576 \
    --outdir outputs/natural_gradient_local_rate \
    --overwrite
```

Notes:
- `--chunk-size 1048576` is the recommended value for H200-class GPUs (peak
  ~25 GB). Use `131072` (peak ~3.7 GB) or `65536` on smaller GPUs.
- `--device cuda` raises a clear error if CUDA is unavailable; it never silently
  falls back to CPU. A CPU torch build runs with `--device cpu`.
- If a few rows fail on a contended GPU (out-of-memory), re-run only the failed
  rows on a free device with
  `scripts/natural_gradient_local_rate/_patch_failed_rows.py` (re-runs every
  `status != ok` row in place and preserves the run id).

### GPU backend

Select the GPU path with `operator.backend: torch` / `--backend torch` and
`operator.device: cuda` / `--device cuda` (`backend: auto` uses torch only when a
CUDA device is available). The torch path uses a dense `torch.linalg.eigh`
eigensolver for `N_theta <= explicit_dense_max_N_theta` and is numerically
identical to the NumPy/SciPy CPU path on the same sample bank; it changes only
the speed, not the meaning, of the estimates. Potential centering and the
`H_sym` accumulation run on-device for the production grid.

## Reports

```bash
# 1. regenerate every report figure (PDF + PNG) and LaTeX table fragment
python reports/make_report_assets.py
#    -> reports/assets/figs/*.pdf, *.png  and  reports/assets/tab_*.tex

# 2. compile the two reports (tectonic resolves preamble.tex and assets/)
cd reports
tectonic affine_invariant_omega_tau_report.tex
tectonic natural_gradient_local_rate_report.tex
```

`make_report_assets.py` only reads the final CSVs and writes figures/tables; it
does not re-run any dynamics. The reports `\input` a shared
[`reports/preamble.tex`](reports/preamble.tex).

## Repository layout

```
configs/
  omega_tau_modes/                gaussian_target / logconcave_target configs
  natural_gradient_local_rate/    smoke + grid + production configs
src/
  common/                         spd, symspace, monte_carlo, io, plotting style,
                                  torch backend helpers
  omega_tau_modes/                (omega, tau) dynamics, targets, metrics, plotting
  natural_gradient_local_rate/    potentials, operators, linearized rate, torch backend
scripts/
  omega_tau_modes/                grid runners + plotting
  natural_gradient_local_rate/    operator/rate/flow runners, plotting, patch tool
tests/
  common/  omega_tau_modes/  natural_gradient_local_rate/
reports/                          LaTeX reports, shared preamble, asset generator, assets/
docs/specs/                       tracked implementation specs (source of truth)
outputs/                          experiment outputs (final CSVs committed)
```

## Specs

The tracked implementation source of truth lives in
[`docs/specs/`](docs/specs/):
[`affine_invariant_gradient_flow.md`](docs/specs/affine_invariant_gradient_flow.md)
(the `(omega, tau)` flow family) and
[`natural_gradient_local_rate_spec.md`](docs/specs/natural_gradient_local_rate_spec.md)
(the local-rate operators and bounds). The code is kept consistent with these.
