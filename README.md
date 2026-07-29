# Affine-Invariant Gaussian Gradient Flow Experiments

Numerical experiments for Gaussian variational inference, where the variational
family is the non-degenerate Gaussians `N(m, C)` and the objective is
`KL(N(m, C) || target)`. The repository contains seven self-contained experiment
groups, each with its own configs, source modules, scripts, tests, and outputs.

The polished write-ups for these groups are the LaTeX reports in
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

### 3. `natural_gradient_discretization_stepsize` — Riemannian vs KL discretization and stepsize stability

Two time discretizations of the same Gaussian natural gradient flow (`dm/dt =
C g`, `dC/dt = C + C H C`): the **Riemannian-distance** (geodesic) covariance
update `C' = e^{dt} C^{1/2} exp(dt C^{1/2} H C^{1/2}) C^{1/2}` and the **KL/Bregman**
update `C' = (1 + dt) (C^{-1} - dt H)^{-1}`. Both share the explicit mean step
`m' = m + dt C g`. Under the improved KL proof both schemes share the **same**
theorem-safe stepsize scale `dt <= 1/(beta lambda_max)`: the KL proof no longer
carries the obsolete cubic stepsize penalty `max{1, lambda_max^3 / (2 lambda_min^3)}`,
and the only remaining theoretical difference between the schemes is the per-step
contraction factor (`q_riem` vs `q_kl`), not the admissible stepsize. The
experiments now ask how conservative that shared bound is, how the observed
contraction compares against `q_riem(dt)` and `q_kl(dt)`, where the large-step
failures come from, and which scheme is more efficient at matched stepsizes.
Three deterministic 2-D targets (exact Gaussian posterior, non-smooth quartic
log-concave, smooth strongly log-concave) and a scalar `N(0, 1)` diagnostic are
swept over `dt in {0.001, ..., 10}` and each run is classified SPD-feasible /
stable / monotone / accurate.

**Finding.** With the shared theorem-safe scale `1/(beta lambda_max)`, both
schemes remain stable and monotone well beyond the sufficient bound on these
deterministic targets, with comparable empirical/theory ratios — there is no
KL-specific stepsize disadvantage from the theory. The scalar diagnostic isolates
the mechanism: the KL covariance update is unconditionally SPD for log-concave
targets, and the large-stepsize failures of both schemes are driven by the shared
explicit mean update, not by the covariance step. On the common theorem-safe rate
grid the observed contraction is faster than either theorem factor predicts, with
`q_kl` the more conservative of the two (now a genuine in-theorem rate, not an
outside-theorem benchmark). A matched-stepsize convergence-speed study adds the
practical counterpart: on both globally smooth targets the Riemannian (geodesic)
update reaches a smaller terminal energy gap at every convergent stepsize, with
the margin growing as the stepsize coarsens, while KL leads only on the non-smooth
quartic and only modestly. See
[`reports/natural_gradient_discretization_stepsize_report.tex`](reports/natural_gradient_discretization_stepsize_report.tex).

### 4. `wfr_gradient_flow` — Wasserstein–Fisher–Rao splitting and phase separation

The Wasserstein–Fisher–Rao (WFR) Gaussian flow combines the Bures–Wasserstein
transport flow with the Fisher–Rao (natural-gradient) flow at transport strength
`lambda_t`: `dm/dt = (C + lambda I) g`, `dC/dt = (C + C H C) + lambda (2I + C H +
H C)`. We discretize it by forward–backward operator splitting — one Wasserstein
half-step `C' = 1/2 (Ctilde + 2hI + [Ctilde(Ctilde + 4hI)]^{1/2})`,
`Ctilde = (I + hH) C (I + hH)`, then one Fisher–Rao (KL) half-step — and compare
five methods (`fr_only`, `w_only`, `wfr_fixed`, `wfr_theory`, `wfr_adaptive`) on
two ill-conditioned 2-D targets (an exact anisotropic Gaussian posterior `N(0,
diag(1, Lambda))` and a smooth nonseparable log-concave log-cosh posterior) from
a far, underdispersed start. Convergence is measured as the **energy gap versus
iteration**. All expectations are exact (closed form / Gauss–Hermite); CPU only.

**Finding.** The flows **phase-separate**. Wasserstein transport — whose additive
`2hI` term inflates an underdispersed covariance at a scale-independent rate —
wins the warmup (smallest gap in the first few iterations) but plateaus far from
the optimum; the Fisher–Rao step is slow during warmup but converges fastest in
the tail once the covariance is calibrated to the curvature. The WFR splitting
inherits both, reaching every tolerance in roughly half the iterations of pure
Fisher–Rao, and the curvature-adaptive schedule `h_n = h_max/(1 + (s_n/s0)^2)`
with `s_n = lambda_min(C^{1/2}(-H)C^{1/2})` — large transport while underdispersed
(`s_n << 1`), decaying as the covariance calibrates (`s_n -> 1`) — converges
fastest across every `(Lambda, epsilon)` regime. See
[`reports/wfr_gradient_flow_report.tex`](reports/wfr_gradient_flow_report.tex).

### 5. `natural_gradient_nonconvex_instability` — nonconvex FR instability and projected KL stationarity

A deterministic one-dimensional stress test for the smooth non-log-concave
target `V_R(x) = 0.5 x^2 - 2 R^2 log cosh(x/R)`, whose Hessian is globally
bounded by `-1 <= V_R'' <= 1`. The symmetric experiments use `m=0` and
Gauss-Hermite quadrature to compute `A_R(c)=E[V_R''(X)]`,
`F(c)=E[V_R(X)]-0.5 log c`. The Riemannian Fisher-Rao covariance step is tested
for a negative-curvature cascade, the KL/Bregman Fisher-Rao covariance step is
tested near its rational pole, and the Bures-Wasserstein forward-backward step is
checked against the nonconvex running-minimum BW-gradient stationarity envelope.
A fourth experiment adds the **projected (clipped) KL scheme** of the manuscript:
the KL covariance update is clipped back into a spectral interval
`[lambda_minus, lambda_plus]` and, at theorem-safe stepsizes
`dt = safety / L_clip` with `L_clip = 2*beta*lambda_plus` (so
`dt <= 1/(2*beta*lambda_plus)`), the running-minimum Bregman displacement
`D_min(N) = min_{0<=n<N} KL(rho_{a_n}||rho_{a_{n+1}})` is checked against the
energy-drop envelope `B_N = (dt/N)*(F(c_0) - F(c_N))`. The theorem-safe scale
depends on `lambda_plus` and `beta` only, not on `lambda_minus` (a sweep
confirms `dt_theory` is flat in `lambda_minus` and scales like `1/lambda_plus`).

**Finding.** Bounded Hessian smoothness alone does not numerically control the
trajectory-wise Fisher-Rao gradient norm independently of covariance: both
unprojected Fisher-Rao discretizations can produce very large FR gradients on
this target. With covariance clipping the pole disappears in the theorem-safe
regime and the projected KL scheme satisfies the constrained Bregman
stationarity check `D_min(N) <= B_N` (max violation `0`, `theorem_check_pass`
true for every `R`) across the safety fractions `0.25, 0.5, 0.9`. Once the upper
covariance constraint binds the covariance is pinned at `lambda_plus` and `D_n`
becomes zero — this is *constrained* stationarity on the feasible-set boundary,
not a small *unconstrained* Fisher-Rao gradient, and it does not imply
convergence to the unconstrained Gaussian VI optimum. A small comparison
additionally runs the outside-theorem stepsize `dt = 1/(beta*lambda_plus)`, which
is exactly twice the theorem edge `1/(2*beta*lambda_plus)` (so `dt*L_clip = 2`,
labelled `outside_theorem_2x`): the envelope still holds empirically because the
trajectory reaches the active constraint in a single step — a single-example
observation, not a relaxation of the theorem's hypotheses. The BW comparison is
deliberately different again: it uses the running minimum of the BW gradient
norm, matching the nonconvex stationarity theory, not a pointwise Fisher-Rao
norm. See
[`reports/natural_gradient_nonconvex_instability_report.tex`](reports/natural_gradient_nonconvex_instability_report.tex).

### 6. `natural_gradient_stl_variance` — Sticking-the-Landing variance reduction for stochastic Gaussian natural-gradient schemes

Stochastic Gaussian natural-gradient schemes estimate the Gaussian expectations
`g = E_q[score_post]` and `H = E_q[Hess_log_post]` from samples `theta ~ N(m, C)`.
The **Sticking-the-Landing** (STL) trick replaces the plain posterior-score mean
estimator `b_base = score_post(theta)` by the score-residual estimator
`b_stl = score_post(theta) + C^{-1}(theta - m)`, which subtracts the (zero-mean)
score of the variational Gaussian itself. It is unbiased and, at a well-specified
Gaussian optimum, vanishes pointwise. We isolate mean-block STL (the covariance
update is held fixed: for the direct Hessian sampler the covariance "residual"
rewrite `K = S + C^{-1}` is algebraically identical to the original update) with
an estimator-level variance study at six fixed states and an algorithm-level
trajectory study over `{riemannian, riemannian_stl, kl, kl_stl}`, on a
well-specified anisotropic Gaussian target and a misspecified smooth strongly
log-concave `log cosh` target. float64; single-GPU by default (the algorithm grid
may split across at most two GPUs).

**Finding.** For the Gaussian target the STL mean estimator has *exactly zero*
variance whenever `C = C_star` (the optimum and every mean-only perturbation), and
STL drops the algorithm's stochastic noise floor one-to-two orders of magnitude.
For the misspecified `log cosh` target STL still helps near the optimum but cannot
remove the floor (the nonlinear `tanh` part of the score is not captured by the
Gaussian residual). STL is not universally beneficial: under covariance
under-dispersion (`C` too concentrated, so `C^{-1}` is large) it *inflates* the
variance above the baseline. See
[`reports/natural_gradient_stl_variance_report.tex`](reports/natural_gradient_stl_variance_report.tex).

### 7. `natural_gradient_covariance_bootstrap` — covariance bootstrap and enhanced global rates

Numerical characterization of the manuscript note on enhanced global rates for the
Gaussian natural-gradient (Fisher–Rao) flow, in the curvature convention
`G = E[grad V]`, `A = E[grad^2 V]` with covariance updates
`C' = C^{1/2} exp(dt (I - C^{1/2} A C^{1/2})) C^{1/2}` (Riemannian) and
`C' = (1 + dt)(C^{-1} + dt A)^{-1}` (KL). Four experiments on an exact diagonal
Gaussian target and a smooth strongly log-concave `log cosh` target: (1) the
Fisher–Rao covariance burn-in `N_cov` to reach the curvature scale `1/(2 beta)`;
(2) dynamic (growing covariance lower envelope `L_n`) versus frozen energy-gap
contraction; (3) one Wasserstein/Bures bootstrap step
`C_b = 1/2 (Ctilde + 2 eta I + [Ctilde(Ctilde + 4 eta I)]^{1/2})`, `eta = c/beta`,
followed by a permanent Fisher–Rao tail; (4) the stochastic Sticking-the-Landing
(STL) noise floor and the intrinsic Hessian fluctuation
`Psi = E||C^{1/2}(hess V(X) - A) C^{1/2}||_F^2`. Experiments 1–3 are deterministic
(closed-form / Gauss–Hermite, CPU only); experiment 4 is stochastic (float64,
optional single CUDA device).

**Finding.** The old frozen covariance lower bound `min(lambda0, 1/beta)` is
pessimistic. The observed covariance burn-in is **logarithmic** in `1/(beta
lambda0)` with slope `1/log(1+dt)` (measured 2.39 vs. theory 2.47 for KL,
`R^2 > 0.99`) and flat against `1/lambda0` — not a `1/lambda0` law. The dynamic
`L_n`-based energy-gap contraction certifies progress that the frozen envelope
cannot (both remain valid, conservative upper bounds; the observed rate is
faster). One Wasserstein/Bures bootstrap step lifts `lambda_min(C)` to the
curvature scale `c/beta` in a single iteration and removes the
`lambda0`-dependent warm-up (40–50 pure-FR iterations vs. 13–17 with the
bootstrap at `lambda0 = 1e-8`, independent of `lambda0`). STL has **no**
asymptotic noise floor for the Gaussian target (machine zero) and an
`O(dt Psi)` floor for the non-Gaussian target (log–log slope ≈ 1.0). See
[`reports/natural_gradient_covariance_bootstrap_report.tex`](reports/natural_gradient_covariance_bootstrap_report.tex).

### 8. `natural_gradient_sharp_bump` — is the `O(kappa^2)` global localization intrinsic?

Runs the sharp counterexample of `thm:sharp-disc` (manuscript Appendix C) itself.
The 1-D bump-train potential `V_kappa` has exact condition number `kappa`,
`||V'''||_inf <= LH` dimension-free, and matched initial covariance `c_0 = 1/kappa`;
at the globally certified fixed step `dt = gamma/kappa` it blocks both schemes for
`N = T kappa^2 / gamma` iterations. The appendix writes the center recursion with
the constant `gamma/kappa^2`, which is exactly the per-step mean gain `dt * c_kappa`;
here that constant is a free parameter `s`, so the same construction can be
**retuned** to any stepsize (`s = dt / kappa`) — the adversary sees the stepsize rule
before choosing `V`. Two families (`manuscript`, `retuned`) are crossed with the two
schemes and two **fixed** stepsize arms: `theory` (`dt = gamma/kappa`, the
certified step) and `const` (`dt = gamma`, `kappa`-independent). Both are single
scalars held constant for the whole run — no line search, no state dependence. A
separate sweep varies the order-one constant `gamma` itself. The measured quantity
is iterations to `DeltaE <= DeltaE_0 / 2` versus `kappa`, reported as a log-log
slope.

**Finding.** The certified step reproduces the theorem: slope `1.93` over
`kappa in [128, 1024]` (`-> 2`), with the iterates shadowing the ideal centers `x_j`
to machine precision, `c_n A_n = 1` throughout, and `n_half / N_train -> 1.04`. An
order-one step is **not** destabilized by the construction — every run is monotone,
no failures — and on the manuscript train it defeats the obstruction outright
(slope `0.16`, essentially logarithmic). When the adversary retunes the train to
the order-one step the cost returns to slope `0.91` (Riemannian) / `0.90` (KL),
i.e. `Theta(kappa)`, matching the Bures–Wasserstein complexity of
arXiv:2304.05398. Both schemes agree to within two iterations everywhere — along
the shadowed trajectory `c_n A_n = 1`, at which the two covariance maps coincide
exactly. The order-one reading is not an artifact of `gamma = 0.5`: every
`gamma <= 1` is monotone and convergent at every `kappa` with exponent in
`[0.88, 0.92]`, and the fixed step breaks only at `gamma >= 2`, a `kappa`-free
threshold matching the explicit mean update's `dt * c_n A_n < 2` limit. So the `kappa^2` is
the price of committing to `dt = Theta(1/kappa)`, not an intrinsic property of the
Fisher–Rao scheme: the construction certifies `Omega(T kappa / dt)`, which is
`kappa^2` at the certified step and `kappa` at an order-one step. This is evidence
on the sharp example only — it does not supply a stability proof for `dt = Theta(1)`
outside the certified range `dt <= 1/(beta lambda_max)`. See
[`reports/natural_gradient_sharp_bump_report.tex`](reports/natural_gradient_sharp_bump_report.tex).

## Which outputs are final

Only these directories are interpreted as evidence in the reports:

```
outputs/gaussian_grid/                              omega/tau, Gaussian target
outputs/logconcave_grid/                            omega/tau, log-concave target
outputs/natural_gradient_local_rate/operator_grid/        local rate: Lambda_hat + gamma_loc
outputs/natural_gradient_local_rate/linearized_rate_grid/ local rate: gamma_loc + eigenvectors
outputs/natural_gradient_discretization_stepsize/         Riemannian vs KL stepsize stability
outputs/wfr_gradient_flow/                                WFR splitting: phase separation + schedules
outputs/natural_gradient_nonconvex_instability/           nonconvex FR instability + BW bound
outputs/natural_gradient_stl_variance/                   STL estimator + algorithm variance reduction
outputs/natural_gradient_covariance_bootstrap/           covariance bootstrap: burn-in, contraction, W-boot, STL floor
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

## Reproducing the discretization-stepsize experiments

Deterministic and CPU-only (closed-form / Gauss–Hermite expectations, no Monte
Carlo, no GPU). The final outputs live at
`outputs/natural_gradient_discretization_stepsize/`:

```bash
python scripts/natural_gradient_discretization_stepsize/run_stepsize_grid.py \
    --config configs/natural_gradient_discretization_stepsize/stepsize_grid.yaml \
    --outdir outputs/natural_gradient_discretization_stepsize --overwrite

python scripts/natural_gradient_discretization_stepsize/plot_results.py \
    --outdir outputs/natural_gradient_discretization_stepsize
```

The runner writes `results_long.csv`, `summary.csv`, `stepsize_summary.csv`,
`scalar_diagnostic.csv`, and `target_metadata.json`. Add `--smoke` for the fast
reduced grid (`outputs/natural_gradient_discretization_stepsize_smoke/`).

### Supplementary theoretical-rate benchmark

A supplementary experiment compares the theorem-predicted contraction factors of
both schemes against the observed contraction, on the common theorem-safe
reference stepsize `dt_ref = 1/(beta * lambda_max)`. Under the improved KL proof
both schemes are admitted on `dt <= dt_ref`, so `q_riem(dt)` and `q_kl(dt)` on
this grid are genuine in-theorem contraction rates. Only the two globally smooth
targets are run, since the non-smooth quartic has no global `beta`:

```bash
python scripts/natural_gradient_discretization_stepsize/run_rate_benchmark.py \
    --config configs/natural_gradient_discretization_stepsize/rate_benchmark.yaml \
    --outdir outputs/natural_gradient_discretization_stepsize --overwrite
```

This writes the prefixed files `rate_results_long.csv`, `rate_summary.csv`,
`rate_tolerance_summary.csv`, and `rate_metadata.json` into the same output
directory (it never overwrites the stepsize-study CSVs). Add `--smoke` for the
fast reduced grid (`outputs/natural_gradient_discretization_stepsize_rate_smoke/`).
The report assets (including the rate figures and
`tab_discretization_rate_summary.tex`) are regenerated by:

```bash
python reports/make_report_assets.py --only discretization
```

and the report compiled by:

```bash
cd reports
tectonic natural_gradient_discretization_stepsize_report.tex
```

## Reproducing the WFR gradient-flow experiments

Deterministic and CPU-only (closed-form / Gauss–Hermite expectations, no Monte
Carlo, no GPU). The final outputs live at `outputs/wfr_gradient_flow/`:

```bash
python scripts/wfr_gradient_flow/run_wfr_grid.py \
    --config configs/wfr_gradient_flow/wfr_grid.yaml \
    --outdir outputs/wfr_gradient_flow --overwrite
```

The runner writes `results_long.csv`, `summary.csv`, `hitting_times.csv`,
`schedule_sweep.csv`, `dt_sweep.csv`, `target_metadata.json`, and
`run_metadata.json`. Add `--smoke` for the fast reduced grid
(`outputs/wfr_gradient_flow_smoke/`). The report assets (five figures and the
`tab_wfr_*` tables) are regenerated by:

```bash
python reports/make_report_assets.py --only wfr
```

and the report compiled by:

```bash
cd reports
tectonic wfr_gradient_flow_report.tex
```

## Reproducing the nonconvex-instability experiments

Deterministic and CPU-only (Gauss-Hermite expectations, no Monte Carlo, no GPU).
The final outputs live at `outputs/natural_gradient_nonconvex_instability/`:

```bash
python scripts/natural_gradient_nonconvex_instability/run_experiments.py \
    --config configs/natural_gradient_nonconvex_instability/nonconvex_instability.yaml \
    --outdir outputs/natural_gradient_nonconvex_instability --overwrite

python scripts/natural_gradient_nonconvex_instability/plot_results.py \
    --outdir outputs/natural_gradient_nonconvex_instability
```

The runner writes `results_long.csv`, `summary.csv`, `kl_pole_summary.csv`,
`wasserstein_bound_summary.csv`, `clipped_kl_stationarity.csv`,
`clipped_kl_summary.csv`, `clipped_kl_sweep.csv`, `target_metadata.json`, and
`run_metadata.json`. Add `--smoke` for the fast reduced grid
(`outputs/natural_gradient_nonconvex_instability_smoke/`). The
`clipped_kl_stationarity.csv`/`clipped_kl_summary.csv` files carry a `dt_rule`
column distinguishing the theorem-safe stepsizes (`dt = safety / (2*beta*lambda_plus)`
at safety `0.25, 0.5, 0.9`) from the `outside_theorem_2x` comparison
`dt = 1/(beta*lambda_plus)` (`dt*L_clip = 2`); `clipped_kl_sweep.csv` records the
theorem-safe stepsize `1/(2*beta*lambda_plus)` as `lambda_minus`/`lambda_plus`
vary. `plot_results.py` adds the clipped KL figures
`fig_nonconvex_clipped_kl_stationarity.{pdf,png}`,
`fig_nonconvex_clipped_covariance.{pdf,png}`, the outside-theorem comparison
`fig_nonconvex_clipped_kl_largestep.{pdf,png}`, and the stepsize-scaling sweep
`fig_nonconvex_clipped_kl_theory_sweep.{pdf,png}` alongside the existing figures.
The report assets (including the clipped KL figures,
`tab_nonconvex_clipped_kl_summary.tex`,
`tab_nonconvex_clipped_kl_largestep.tex`, and
`tab_nonconvex_clipped_kl_sweep.tex`) are regenerated by:

```bash
python reports/make_report_assets.py --only nonconvex_instability
```

and the report compiled by:

```bash
cd reports
tectonic natural_gradient_nonconvex_instability_report.tex
```

## Reproducing the STL variance experiments

Two experiments share the output directory `outputs/natural_gradient_stl_variance/`:
an estimator-level variance comparison and an algorithm-level trajectory grid.
Each script's `--overwrite` clears only its own outputs, so the estimator and the
algorithm runs coexist. The smoke grids (CPU, fast) validate the scripts first:

```bash
# smoke (CPU): estimator then algorithm into the smoke directory
python scripts/natural_gradient_stl_variance/run_estimator_variance.py \
    --config configs/natural_gradient_stl_variance/stl_variance_smoke.yaml \
    --outdir outputs/natural_gradient_stl_variance_smoke --device cpu --overwrite

python scripts/natural_gradient_stl_variance/run_algorithm_grid.py \
    --config configs/natural_gradient_stl_variance/stl_variance_smoke.yaml \
    --outdir outputs/natural_gradient_stl_variance_smoke --device cpu --overwrite
```

The production grids run on a single CUDA GPU by default (developed on an NVIDIA
H200). The algorithm grid accepts `--max-gpus 2` to split its cells across at
most two GPUs (it never launches more than two GPU jobs):

```bash
# production (GPU)
python scripts/natural_gradient_stl_variance/run_estimator_variance.py \
    --config configs/natural_gradient_stl_variance/stl_variance.yaml \
    --outdir outputs/natural_gradient_stl_variance --device cuda --overwrite

python scripts/natural_gradient_stl_variance/run_algorithm_grid.py \
    --config configs/natural_gradient_stl_variance/stl_variance.yaml \
    --outdir outputs/natural_gradient_stl_variance --device cuda --max-gpus 2 --overwrite
```

The estimator runner writes `estimator_variance.csv`,
`estimator_variance_summary.csv`, `estimator_target_metadata.json`, and
`estimator_run_metadata.json`. The algorithm runner writes
`algorithm_results_long.csv`, `algorithm_summary.csv`, `tail_noise_floor.csv`,
`target_metadata.json`, and `run_metadata.json`. Per-group figures (written to
`outputs/natural_gradient_stl_variance/figures/`), the report figures, and the
report itself are produced by:

```bash
python scripts/natural_gradient_stl_variance/plot_results.py \
    --outdir outputs/natural_gradient_stl_variance

python reports/make_report_assets.py --only stl_variance

cd reports
tectonic natural_gradient_stl_variance_report.tex
```

## Reproducing the covariance-bootstrap experiments

Four experiments share the output directory
`outputs/natural_gradient_covariance_bootstrap/`. Experiments 1–3 are
deterministic (closed-form / Gauss–Hermite, CPU/NumPy only); experiment 4 is the
stochastic STL noise floor (float64). The single runner drives all four and is
resumable at block granularity (the deterministic block and the stochastic block
are each skipped if their outputs exist, unless `--overwrite` is passed):

```bash
# smoke (CPU, seconds): reduced grids into a separate directory
python scripts/natural_gradient_covariance_bootstrap/run_experiments.py \
    --config configs/natural_gradient_covariance_bootstrap/covariance_bootstrap_smoke.yaml \
    --outdir outputs/natural_gradient_covariance_bootstrap_smoke --overwrite

# production (CPU): ~1.5 min end to end on a workstation
python scripts/natural_gradient_covariance_bootstrap/run_experiments.py \
    --config configs/natural_gradient_covariance_bootstrap/covariance_bootstrap.yaml \
    --outdir outputs/natural_gradient_covariance_bootstrap --overwrite --backend auto
```

The runner writes `results_long.csv`, `covariance_bootstrap_summary.csv`,
`contraction_benchmark.csv`, `wasserstein_bootstrap_summary.csv`,
`stl_floor_summary.csv`, `target_metadata.json`, and `run_metadata.json`. The
`--smoke` flag also shrinks any config in place. Per-group figures are written by:

```bash
python scripts/natural_gradient_covariance_bootstrap/plot_results.py \
    --outdir outputs/natural_gradient_covariance_bootstrap
```

**Stochastic STL on a GPU (optional).** The deterministic experiments are always
CPU/NumPy. Experiment 4 accepts a single CUDA device, but for these small diagonal
problems the CPU path is faster (per-step kernel-launch/sync overhead dominates on
the GPU), so the committed evidence is the CPU run. To run the stochastic block on
one allowed GPU (physical index `>= 4`, never GPUs 0–3), first inspect
`nvidia-smi`, pick one idle GPU, and export its physical id:

```bash
nvidia-smi
export CUDA_VISIBLE_DEVICES=6      # one idle physical GPU with index >= 4
python scripts/natural_gradient_covariance_bootstrap/run_experiments.py \
    --config configs/natural_gradient_covariance_bootstrap/covariance_bootstrap.yaml \
    --outdir outputs/natural_gradient_covariance_bootstrap \
    --backend torch --device cuda --overwrite
```

Inside Python that GPU appears as `cuda:0`; the chosen physical id and
`CUDA_VISIBLE_DEVICES` are recorded in `run_metadata.json`. The report assets
(seven figures and four `tab_covboot_*` tables) are regenerated by:

```bash
python reports/make_report_assets.py --only covariance_bootstrap
```

and the report compiled by:

```bash
cd reports
tectonic natural_gradient_covariance_bootstrap_report.tex
tectonic natural_gradient_sharp_bump_report.tex
```

## Reproducing the sharp bump-train experiment

Deterministic, 1-D, CPU-only (~95 s for the full `kappa` grid up to 1024):

```bash
python scripts/natural_gradient_sharp_bump/run_sharp_bump.py \
    --config configs/natural_gradient_sharp_bump/sharp_bump.yaml \
    --outdir outputs/natural_gradient_sharp_bump --overwrite
```

```bash
python scripts/natural_gradient_sharp_bump/plot_sharp_bump.py \
    --indir outputs/natural_gradient_sharp_bump
```

This writes `sharp_bump_long.csv`, `sharp_bump_summary.csv`,
`sharp_bump_slopes.csv`, `sharp_bump_gamma_sweep.csv` and
`sharp_bump_metadata.json`, plus
`figures/sharp_bump_scaling.{png,pdf}`. Add `--smoke` for the reduced grid
(`outputs/natural_gradient_sharp_bump_smoke/`). The headline columns are
`slope_upper` in `sharp_bump_slopes.csv` (log-log exponent of `n_half` vs `kappa`
over the top half of the grid) and `shadow_mean_err_over_w` / `cA_max` in the
summary, which verify that the `theory` arm really is on the shadowed trajectory
of `lem:bump-shadowing` rather than merely slow.

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
#    (use --only {omega_tau,local_rate,discretization,wfr,nonconvex_instability,
#     stl_variance,covariance_bootstrap,sharp_bump} to build one group)

# 2. compile the reports (tectonic resolves preamble.tex and assets/)
cd reports
tectonic affine_invariant_omega_tau_report.tex
tectonic natural_gradient_local_rate_report.tex
tectonic natural_gradient_discretization_stepsize_report.tex
tectonic wfr_gradient_flow_report.tex
tectonic natural_gradient_nonconvex_instability_report.tex
tectonic natural_gradient_stl_variance_report.tex
tectonic natural_gradient_covariance_bootstrap_report.tex
tectonic natural_gradient_sharp_bump_report.tex
```

`make_report_assets.py` only reads the final CSVs and writes figures/tables; it
does not re-run any dynamics. The reports `\input` a shared
[`reports/preamble.tex`](reports/preamble.tex). A group whose final outputs are
absent on a given checkout is skipped (with a notice) so the others still build.

## Repository layout

```
configs/
  omega_tau_modes/                gaussian_target / logconcave_target configs
  natural_gradient_local_rate/    smoke + grid + production configs
  natural_gradient_discretization_stepsize/  stepsize_grid config
  wfr_gradient_flow/              wfr_smoke + wfr_grid configs
  natural_gradient_nonconvex_instability/  nonconvex_instability config
  natural_gradient_stl_variance/  stl_variance + stl_variance_smoke configs
  natural_gradient_covariance_bootstrap/  covariance_bootstrap + _smoke configs
src/
  common/                         spd, symspace, monte_carlo, io, plotting style,
                                  torch backend helpers
  omega_tau_modes/                (omega, tau) dynamics, targets, metrics, plotting
  natural_gradient_local_rate/    potentials, operators, linearized rate, torch backend
  natural_gradient_discretization_stepsize/  targets, methods, metrics, ode_reference,
                                  optimize_star, runner, plotting
  wfr_gradient_flow/              targets, methods (W/FR half-steps), schedules,
                                  runner, metrics, plotting
  natural_gradient_nonconvex_instability/  target, scalar methods, runner, plotting
  natural_gradient_stl_variance/  targets, estimators, linalg backend, states,
                                  estimator_variance, algorithm, grid, metrics, plotting
  natural_gradient_covariance_bootstrap/  targets, methods, envelopes, metrics,
                                  runner (det.), stochastic (STL), plotting
scripts/
  omega_tau_modes/                grid runners + plotting
  natural_gradient_local_rate/    operator/rate/flow runners, plotting, patch tool
  natural_gradient_discretization_stepsize/  stepsize grid runner + plotting
  wfr_gradient_flow/              WFR grid runner
  natural_gradient_nonconvex_instability/  nonconvex runner + plotting
  natural_gradient_stl_variance/  estimator + algorithm runners, plotting
  natural_gradient_covariance_bootstrap/  run_experiments + plot_results
tests/
  common/  omega_tau_modes/  natural_gradient_local_rate/
  natural_gradient_discretization_stepsize/  wfr_gradient_flow/
  natural_gradient_nonconvex_instability/  natural_gradient_stl_variance/
  natural_gradient_covariance_bootstrap/
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
(the local-rate operators and bounds), and
[`natural_gradient_covariance_bootstrap.md`](docs/specs/natural_gradient_covariance_bootstrap.md)
(the covariance envelopes, dynamic contraction, Bures bootstrap, and STL floor).
The code is kept consistent with these.
