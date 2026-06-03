# Affine-Invariant Gaussian Gradient Flow

A clean, reproducible Python experiment repo for studying the parameter effects
of (ω, τ) in the Riemannian-distance discretization of affine-invariant Gaussian
gradient flows.  The first version implements the **Gaussian target N(0, Iₙ)**,
where all expectations are exact — no Monte Carlo or quadrature needed.

---

## Scientific background

### Variational inference as gradient flow

We study variational inference as gradient flow of the KL divergence
KL(q ‖ π) over the manifold of Gaussians q = N(m, C), equipped with the
affine-invariant (Fisher–Rao-like) Riemannian metric parameterized by (ω, τ).

The resulting **Riemannian-distance discretization** gives the following update
at each step (derived in §2 of the associated paper):

**Mean update:**
```
m_{k+1} = m_k − Δt · C_k m_k
```

**Covariance update (matrix exponential form):**
```
C_{k+1} = C_k^{1/2}
             exp( Δt/(2ω) · [ −C_k + α I ] )
           C_k^{1/2}

where  α = (ω + τ Tr(C_k)) / (ω + n τ)
```

The matrix exponential is essential: without it the update is not the
Riemannian exponential-map step and does not preserve positive definiteness.

**Why the eigenvectors are preserved.**
Because both C_k and the exponent matrix are functions of the same eigenbasis,
they commute.  The update reduces to a scalar rescaling of each eigenvalue:

```
λᵢ_{k+1} = λᵢ · exp( Δt/(2ω) · (−λᵢ + α) )
```

with the eigenvectors Q of C_k unchanged.  This makes the implementation
exact and free of any matrix-exponential routine.

---

### Parameters ω and τ

| Symbol | Role | Constraint |
|--------|------|-----------|
| `ω > 0` | Overall covariance update rate — scales how fast all eigenvalues relax toward their equilibrium | `ω > 0` |
| `τ` | Trace-weighting — shifts the equilibrium target `α` up or down depending on whether Tr(C) is too large or too small | `ω + n τ > 0` |

**Intuition for ω.**
The mean update rate is fixed by Δt and C, independently of ω.
The covariance update rate scales as 1/(2ω).
Smaller ω → faster covariance convergence; larger ω → slower.
The choice ω = 1/2 with τ = 0 corresponds to the balanced Fisher–Rao flow.

**Intuition for τ.**
When τ = 0, the target eigenvalue is α = 1 for all i, independent of Tr(C).
Each eigenvalue independently drifts toward 1 — this is the standard Fisher–Rao flow.

When τ ≠ 0, the scalar α depends on Tr(C):
- If Tr(C) > n (covariance volume too large) and τ < 0,
  then α < 1, which *increases* the drift rate pushing eigenvalues down.
  This accelerates volume shrinkage.
- Conversely, τ > 0 reduces the drift rate and *slows* volume correction.

The specific choice τ = −ω/(2n) makes ω + nτ = ω/2, halving the denominator
and effectively doubling the volume correction speed in the local Gaussian theory.

**Key implication:** τ acts only on the *trace/volume* part of the covariance
error.  It does not accelerate traceless shape modes or mean convergence.

---

### Gaussian target: exact expectations

For target π = N(0, Iₙ):
```
∇_θ log π(θ) = −θ          ⟹  𝔼_{N(m,C)}[∇ log π] = −m
∇²_θ log π(θ) = −I          ⟹  𝔼_{N(m,C)}[∇² log π] = −I
```

Substituting into the general discrete scheme gives the closed-form updates
implemented in `src/dynamics.py`.

---

## Initializations

Five initial conditions are designed to **isolate distinct convergence modes**,
so that the roles of ω and τ can be disentangled.

### `mean_only` — pure mean error

```
m₀ = r · 1/√n,   C₀ = I,   r = 3
```

The covariance is already at the target; only the mean is displaced.

**What it tests:** whether ω and τ affect mean-dominated convergence.

**Expected finding:** all three τ values (τ₋, τ₀, τ₊) behave nearly identically,
because τ acts only on the covariance-volume dynamics and C₀ = I is already fixed.

---

### `volume_high` — pure volume expansion error

```
m₀ = 0,   C₀ = 4I
```

All eigenvalues are 4 (too large), so the entire covariance error is volume
(scalar scale).  There is no shape anisotropy and no mean error.

**What it tests:** the τ < 0 acceleration hypothesis.

**Expected finding:** τ₋ converges significantly faster than τ₀, and τ₊ is
slower.  This is the regime where τ < 0 provides the clearest benefit.

---

### `volume_low` — pure volume compression error

```
m₀ = 0,   C₀ = 0.25 · I
```

All eigenvalues are 0.25 (too small); same structure as `volume_high` but
the covariance must expand rather than contract.

**What it tests:** τ effect on volume expansion (not just contraction).

**Expected finding:** same qualitative pattern as `volume_high` — τ₋ is faster,
τ₊ is slower — confirming the τ acceleration is symmetric in direction.

---

### `shape_only` — pure shape (anisotropy) error

```
m₀ = 0,   C₀ = diag(e^r, e^{−r}, 1, …, 1),   r = 2
```

By construction det(C₀) = e^r · e^{−r} · 1 · … · 1 = 1, so the **volume is
exactly correct**.  All of the covariance error is anisotropy: the first axis
is too large, the second is too small, the rest are at target.

**What it tests:** whether τ helps when the error is purely in the traceless
(shape) part of the covariance.

**Expected finding:** τ₋ ≈ τ₀ ≈ τ₊.  Since there is no net volume error, the
trace shift introduced by τ buys nothing.  The relevant parameter here is ω,
which controls how fast each eigenvalue relaxes.

---

### `mixed` — simultaneous mean + volume + shape error

```
m₀ = 2 · 1/√n,   C₀ = s · diag(e^r, e^{−r}, 1, …, 1),   s = 2,   r = 1.5
```

This is the closest to a realistic scenario: nonzero mean, inflated volume
(det = s^n · 1 > 1), and anisotropic shape.

**What it tests:** whether τ < 0 gives a net benefit in the most common practical case.

**Expected finding:** τ₋ may help during the early volume-dominated transient,
but the final convergence rate is limited by mean and shape modes (which τ
does not accelerate).  Benefits are smaller and less guaranteed than in
`volume_high` / `volume_low`.

---

## Metrics

All metrics are computed against the target N(0, Iₙ) at each saved step.

### 1. KL energy  (`kl_energy`)

```
E = KL(N(m, C) ‖ N(0, I))
  = ½ ( ‖m‖² + Tr(C) − log det C − n )
```

The primary scalar convergence measure.  Equals zero iff m = 0 and C = I;
always ≥ 0 by the Gibbs inequality.

This corresponds to the energy gap E(aₙ) − E(a★) in the paper.

---

### 2. Normalised energy  (`norm_energy`)

```
Ê = E / E₀
```

Divides by the initial energy so all runs start at 1 and are comparable across
initializations with very different scales.  The time-to-tolerance thresholds
(1e-2, 1e-4, 1e-6) are defined on this quantity.

---

### 3. Mean error  (`mean_error`)

```
eₘ = ‖m‖₂
```

The Euclidean distance of the current mean from the target mean 0.
Corresponds to the first summary statistic in Figure 5 of the paper.

---

### 4. Relative covariance error  (`cov_error`)

```
e_C = ‖C − I‖_F / √n
```

Total covariance mismatch, normalized by √n so it is comparable across
dimensions.  Corresponds to the second summary statistic in Figure 5.

---

### 5. Volume error  (`volume_error`)

```
e_vol = |log det C / n|
      = |(1/n) Σᵢ log λᵢ|
```

Measures the per-dimension log-volume mismatch.  Zero iff det(C) = 1.

This diagnostic directly reveals whether τ is accelerating the trace/volume mode:
if τ < 0 helps, it should show up here first and most clearly.

---

### 6. Shape error  (`shape_error`)

```
log C = Q diag(log λᵢ) Qᵀ

e_shape = ‖log C − (Tr(log C)/n) I‖_F
```

This removes the scalar (volume) part of log C and retains only the
traceless anisotropy.  Zero iff C = s·I for any scalar s > 0.

If τ < 0 only accelerates volume modes, this metric should not benefit from τ.

Together, `volume_error` and `shape_error` decompose the full `cov_error` into
its two orthogonal components:  **volume** (scalar part of log C)
and **shape** (traceless part of log C).

---

### 7. Cosine test-function error  (`cosine_error`)

For θ ~ N(m, C), the exact identity is:
```
𝔼[cos(qᵀθ + b)] = exp(−½ qᵀCq) · cos(qᵀm + b)
```

The true value under N(0, I) is:
```
exp(−½ ‖q‖²) · cos(b)
```

The error is the absolute difference between these two quantities.

This is the third summary statistic from Figure 5.  The test vector is fixed as
q = (1, 2, …, n)ᵀ / ‖(1, 2, …, n)‖₂  and  b = 0.5.

---

### 8. Eigenvalue extremes  (`eig_min`, `eig_max`)

The smallest and largest eigenvalues of C over time.  Useful for spotting
near-singularity (eig_min → 0) or blow-up (eig_max → ∞), which can occur
with aggressive step sizes or extreme parameters.

---

### 9. Trace dominance ratio  (`chi`)

```
residuals rᵢ = 1 − λᵢ

χ = (Σᵢ rᵢ)² / (n · Σᵢ rᵢ²)
```

χ ∈ [1/n, 1]:
- χ = 1: all residuals are equal — the error is **pure volume** (maximally
  trace-dominated).  This is where τ < 0 helps most.
- χ = 1/n: only one residual is nonzero — the error is **pure shape**
  (maximally anisotropy-dominated).  τ gives no benefit here.

χ tracks whether the covariance error is "isotropic" (χ near 1, τ may help)
or "anisotropic" (χ near 1/n, only ω matters).

---

## Expected qualitative findings

### Effect of τ

| Initialization | τ < 0 vs τ = 0 | Explanation |
|---------------|----------------|-------------|
| `mean_only`   | No difference | τ only acts on covariance volume; C₀ = I is already the target |
| `volume_high` | τ < 0 faster (~2×) | Pure volume error; τ < 0 doubles the trace-mode convergence rate |
| `volume_low`  | τ < 0 faster (~2×) | Same as above; error is pure volume in the other direction |
| `shape_only`  | No difference | det(C₀) = 1; no volume error for τ to accelerate |
| `mixed`       | τ < 0 helps early | Speeds up volume phase; final rate still limited by mean / shape modes |

**τ > 0** is generally worse than τ = 0 for volume-dominated initializations,
and no better elsewhere.

**τ = 0** is the robust, parameter-free default choice.

### Effect of ω

- Smaller ω → faster covariance convergence (eigenvalues relax faster).
- Larger ω → slower covariance, mean convergence is unaffected.
- ω = 1/2 with τ = 0 is the balanced Fisher–Rao choice: it equates the mean
  and covariance natural gradient steps in a specific sense.
- For covariance-dominated initializations, ω < 1/2 can be faster.
- For mean-dominated initialization, varying ω makes no visible difference.

### Overall conclusion

> Smaller ω can accelerate covariance-dominated transients; τ < 0 can
> additionally accelerate *trace-dominated* covariance transients.  However,
> neither provides uniform improvement across mean, shape, and mixed modes.
> The choice (ω, τ) = (1/2, 0) remains the most robust parameter-free choice
> because it balances mean, covariance-shape, and covariance-volume dynamics.

---

## Installation

```bash
git clone <repo-url>
cd AffineInvariantGaussianGradientFlow
pip install -r requirements.txt
```

Python 3.9+ required.  All dynamics use NumPy/SciPy (CPU, float64).
No PyTorch or GPU dependencies.

---

## Running experiments

### Run tests first

```bash
pytest                 # all 93 tests
pytest -v tests/       # verbose output
```

### Smoke run (fast: n=2, T=2)

```bash
python scripts/run_gaussian_grid.py --n 2 --T 2 --dt 0.1 --outdir outputs/smoke
```

### Full default experiment (n ∈ {2, 5, 10}, T=20, dt=0.02)

```bash
python scripts/run_gaussian_grid.py
```

Writes:
- `outputs/gaussian_grid/results_long.csv` — ~45k rows, one per saved step
- `outputs/gaussian_grid/summary.csv` — 225 rows, one per run

Override any default:

```bash
python scripts/run_gaussian_grid.py --dt 0.01 --T 30 --n 5 10 --outdir outputs/fine
```

### Recompute summary from existing long CSV

```bash
python scripts/make_summary_tables.py
```

### Generate figures

```bash
python scripts/plot_gaussian_results.py
```

Figures for each dimension are written to `outputs/gaussian_grid/figures/n{N}/`:

```
outputs/gaussian_grid/figures/
├── n2/
│   ├── fig_tau_effect_omega_half_n2.{png,pdf}
│   ├── fig_omega_sweep_tau_zero_n2.{png,pdf}
│   ├── fig_time_to_tol_heatmap_n2.{png,pdf}
│   └── fig_tau_speedup_heatmap_n2.{png,pdf}
├── n5/
│   └── ...
└── n10/
    └── ...
```

To plot only specific dimensions:

```bash
python scripts/plot_gaussian_results.py --n 5 10
```

---

## Figure descriptions

| Figure | Description |
|--------|-------------|
| `fig_tau_effect_omega_half_n{N}` | 5 rows (inits) × 6 cols (metrics), comparing τ₋/τ₀/τ₊ for ω=0.5 |
| `fig_omega_sweep_tau_zero_n{N}`  | Normalised energy vs time for all ω values (τ=0), one panel per init |
| `fig_time_to_tol_heatmap_n{N}`  | Heatmap of time-to-1e-4 (rows=init, cols=ω, τ=0) |
| `fig_tau_speedup_heatmap_n{N}`  | Speedup ratio T(τ)/T(τ=0) for τ₋ and τ₊ (ω ∈ {1/4, 1/2, 1}) |

---

## Repository structure

```
├── configs/
│   └── gaussian_target.yaml    default experiment configuration
├── src/
│   ├── __init__.py
│   ├── dynamics.py             gaussian_step() — one closed-form update step
│   ├── metrics.py              compute_all_metrics(), kl_energy()
│   ├── initializations.py      get_initialization()
│   ├── plotting.py             figure generation (per-n subdirectories)
│   └── utils.py                SPD utilities, parameter validation
├── scripts/
│   ├── run_gaussian_grid.py    main grid runner → results_long.csv + summary.csv
│   ├── plot_gaussian_results.py  figure generation script
│   └── make_summary_tables.py  recompute summary.csv from results_long.csv
├── tests/
│   ├── test_gaussian_update.py
│   └── test_metrics.py
└── outputs/
    └── gaussian_grid/
        ├── results_long.csv
        ├── summary.csv
        └── figures/
            ├── n2/
            ├── n5/
            └── n10/
```

---

## Summary CSV columns (Gaussian target)

| Column | Description |
|--------|-------------|
| `n`, `omega`, `tau_type`, `tau_value`, `init_name` | Run identity |
| `dt`, `T` | Integration parameters |
| `final_energy` | KL energy at t = T |
| `final_normalized_energy` | E(T) / E(0) |
| `time_to_1e_minus_2/4/6` | First time normalised energy ≤ threshold (inf if not reached) |
| `monotone_energy_bool` | True if normalised energy is non-increasing throughout |
| `min_eig_min_over_time` | Minimum eigenvalue of C seen across all saved steps |
| `max_eig_max_over_time` | Maximum eigenvalue of C seen across all saved steps |

---

---

# Strongly log-concave non-Gaussian target

## Target definition

```
V_rho(x) = 0.5 ||x||² + (rho / m) * sum_{ell=1}^m  log cosh(a_ell^T x)
```

where:
- `a_ell` are **m = 4n** fixed random unit vectors in Rⁿ, drawn once per `target_seed`
- `rho >= 0` controls the coupling strength; rho = 0 recovers the Gaussian case

The posterior is `pi(x) ∝ exp(-V_rho(x))`.

**Gradient:**
```
grad V(x) = x + (rho/m) * sum_ell  tanh(a_ell^T x) * a_ell
```

**Hessian:**
```
Hess V(x) = I + (rho/m) * sum_ell  sech²(a_ell^T x) * a_ell a_ell^T
```

**Strong log-concavity:**
Since `sech²(z) >= 0` and `a_ell a_ell^T` is PSD, we have `Hess V(x) >= I` for all x.
The minimum eigenvalue of Hess V is ≥ 1 everywhere.

**Smoothness:**
Since `sech²(z) <= 1` and the rows of A have unit norm, `Hess V(x) <= (1 + rho) I`.

**Numerical note:**
`sech²(z) = 1 − tanh²(z)` is used throughout to avoid overflow for large |z|.

---

## Sign convention

The algorithm is written in terms of:
```
g(m, C) = E_{N(m,C)}[ grad V(theta) ]    (= -E[grad log pi])
S(m, C) = E_{N(m,C)}[ Hess V(theta) ]    (= -E[Hess log pi])
```

For the Gaussian target `N(0, I)`, `V(x) = 0.5||x||²`, so `g = m` and `S = I`,
recovering the closed-form Gaussian update exactly.

---

## Discrete update for general target

**Mean:**
```
m_{k+1} = m_k − dt * C_k * g(m_k, C_k)
```

**Covariance:**
```
B   = C_k^{1/2} S(m_k, C_k) C_k^{1/2}        (whitened Hessian)
alpha = (omega + tau * Tr(B)) / (omega + n * tau)
M   = dt/(2*omega) * (−B + alpha * I)
C_{k+1} = C_k^{1/2} expm(M) C_k^{1/2}
```

`scipy.linalg.expm` is used because M is not generally a function of C's eigenvectors —
the commutation shortcut that makes the Gaussian case exact does not apply here.

---

## Monte Carlo / QMC expectations

Expectations `g` and `S` are estimated using **K fixed samples** shared across all
(omega, tau, init) comparisons (common random numbers):

```
z_j ~ N(0, I_n),  j = 1,...,K
theta_j = m + L z_j,   where C = L L^T (Cholesky)

g ≈ (1/K) sum_j  grad V(theta_j)
S ≈ (1/K) sum_j  Hess V(theta_j)
```

Samples are generated via **Sobol quasi-Monte Carlo** (scipy.stats.qmc.Sobol),
transformed with the inverse normal CDF, clipped to [1e-12, 1-1e-12] before
the transform to avoid ±inf. For K not a power of two, or n > 21201, falls
back to a seeded NumPy generator.

Default `K = 4096` for dynamics, `K_ref = 8192` for the reference optimum.
Both seeds are fixed per experiment so all runs are reproducible.

---

## Reference Gaussian VI optimum

For the non-Gaussian target, there is no closed-form optimum. We compute the
**best Gaussian approximation** a★ = (m★, C★) by minimising:

```
F(m, C) = E_{N(m,C)}[V(theta)] − 0.5 log det C
```

This is the VI objective (negative ELBO up to an additive constant).

**Parameterisation:**
`C = L Lᵀ` with L lower-triangular, diagonal entries `L_ii = exp(η_i) > 0`.
The parameter vector is `[m | off-diagonal L entries | log-diagonal η]`,
entirely unconstrained.

**Fixed-sample objective:**
```
F(m, L) = mean_j V(m + L z_j) − sum_i log L_ii
```

**Gradients (analytic):**
```
grad_m F = mean_j grad V(theta_j)
grad_L F = mean_j grad V(theta_j) z_j^T  −  L^{-T}
```
Packed as: off-diagonal → direct; diagonal entry η_i → `(grad_L)_ii * L_ii`.

Optimised via `scipy.optimize.minimize(method="L-BFGS-B")` with up to 2000
iterations, starting from `m=0, L=I`.

The result is saved to `reference_optimum.npz` and reused on subsequent runs
unless `--force-optimize` is passed. Because V_rho is even, `||m★||` should be
very small (verified: ~ 1e-4 or less).

---

## Initializations relative to a★

All five initializations are defined in the **coordinate frame of the reference
optimum**, not relative to the identity. Let `C★_sqrt = C★^{1/2}`.

### `mean_only` — pure mean offset from a★
```
m0 = m★ + 3 * C★_sqrt @ (1/√n),    C0 = C★
```
Covariance is already optimal; only the mean is displaced.
**Expected:** τ has no effect (C0 = C★ is already fixed).

### `volume_high` — inflated covariance
```
m0 = m★,    C0 = 4 * C★
```
Volume is 4ⁿ times the optimal. In whitened coordinates R = 4I.
**Expected:** τ < 0 gives ~2× acceleration over τ = 0; τ > 0 is slower.

### `volume_low` — deflated covariance
```
m0 = m★,    C0 = 0.25 * C★
```
Volume is (1/4)ⁿ times the optimal. In whitened coordinates R = 0.25I.
**Expected:** same pattern as `volume_high`.

### `shape_only` — correct volume, wrong shape
```
m0 = m★,    C0 = C★_sqrt @ diag(e^r, e^{-r}, 1,...,1) @ C★_sqrt,    r = 2
```
In whitened coordinates `R = diag(e^2, e^{-2}, 1,...,1)`, so `det(R) = 1`.
Volume matches C★; error is purely in the eigenvector spread.
**Expected:** τ ≈ 0 independent of sign (no trace/volume error to accelerate).

### `mixed` — mean + volume + shape
```
m0 = m★ + 2 * C★_sqrt @ (1/√n)
C0 = C★_sqrt @ [2 * diag(e^1.5, e^{-1.5}, 1,...,1)] @ C★_sqrt
```
All three error modes present simultaneously.
**Expected:** τ < 0 may help the early volume phase; final convergence limited
by mean and shape modes. Benefit is smaller and less certain than `volume_high`.

---

## Metrics (log-concave target)

All metrics are expressed relative to the reference optimum a★ = (m★, C★).
Let `R = C★^{-1/2} C C★^{-1/2}` be the covariance in whitened coordinates.

| Metric | Formula | Interpretation |
|--------|---------|---------------|
| `objective` | `F(m,C) ≈ mean V(θⱼ) − 0.5 log det C` | VI objective |
| `objective_gap` | `F − F★` | Excess objective (raw; may be slightly negative from MC noise) |
| `normalized_objective_gap` | `(F − F★) / gap₀` | Primary convergence metric, starts at 1 |
| `whitened_mean_error` | `‖C★^{-1/2}(m − m★)‖₂` | Mean error in natural units |
| `cov_error` | `‖R − I‖_F / √n` | Total covariance mismatch in whitened coords |
| `volume_error` | `\|log det R / n\|` | Per-dim volume (scale) error |
| `shape_error` | `‖log R − (Tr log R / n) I‖_F` | Anisotropy error, independent of scale |
| `mean_residual` | `‖g‖₂` | Stationarity: mean equation residual |
| `cov_residual` | `‖I − B‖_F`,  B = C^{1/2} S C^{1/2} | Stationarity: covariance equation residual |
| `trace_residual` | `\|Tr(I − B)\| / √n` | Trace/volume part of covariance residual |
| `traceless_residual` | `‖(I−B) − (Tr(I−B)/n) I‖_F` | Shape part of covariance residual |
| `chi` | `(Tr(I−B))² / (n ‖I−B‖_F²)` | Trace dominance ratio ∈ [1/n, 1]; χ≈1 means τ<0 may help |
| `eig_min`, `eig_max` | eigenvalue extremes of C | Numerical health check |
| `cosine_error_to_star` | `\|E[cos(q^T θ+b)]_{m,C} − E[...]_{m★,C★}\|` | Test-function gap to reference |

---

## Expected qualitative findings (log-concave)

| Initialization | τ < 0 effect | Explanation |
|---------------|-------------|-------------|
| `mean_only` | None | C0 = C★; τ acts only on covariance volume |
| `volume_high` | ~2× speedup | Pure volume error; τ<0 doubles trace-mode rate |
| `volume_low` | ~2× speedup | Same, for volume expansion direction |
| `shape_only` | None | det(R0) = 1; no volume error for τ to exploit |
| `mixed` | Moderate, early only | Helps initial volume phase; shape and mean phases unaffected |

**χ as a predictor:**
The initial trace-dominance ratio χ = (Tr residual)² / (n ‖residual‖_F²) predicts
whether τ < 0 helps. When χ ≈ 1 (covariance residual is isotropic / volume-dominated),
τ < 0 accelerates convergence. When χ ≈ 1/n (shape-dominated), τ < 0 gives no benefit.
Figure 5 (`speedup_vs_chi`) tests this directly.

**Overall conclusion:**
> Smaller ω can accelerate covariance-dominated transients; τ < 0 can additionally
> accelerate *trace-dominated* covariance transients. Neither provides uniform
> improvement across all modes. (ω, τ) = (1/2, 0) remains the most robust choice.

---

## Running the log-concave experiment

### Smoke run (fast: n=3, rho=2, K=512, T=2)

```bash
python scripts/run_logconcave_grid.py \
  --n 3 --rho 2 --K 512 --K-ref 1024 --T 2 --dt 0.01 \
  --outdir outputs/logconcave_smoke
```

### Full default experiment (n=5, rho=5, K=4096, T=40)

```bash
python scripts/run_logconcave_grid.py
```

Writes to `outputs/logconcave_grid/`:
- `results_long.csv` — one row per saved step per run (45 runs × ~400 steps)
- `summary.csv` — one row per run
- `reference_optimum.npz` + `reference_optimum_meta.json` — cached Gaussian VI optimum
- `target_metadata.json` — target + sample parameters

Re-run with a fresh reference optimum:
```bash
python scripts/run_logconcave_grid.py --force-optimize
```

### Generate figures

```bash
python scripts/plot_logconcave_results.py
```

Figures are saved to `outputs/logconcave_grid/figures/`:

| File | Description |
|------|-------------|
| `fig_logconcave_tau_effect_omega_half_n{N}_rho{R}` | τ comparison (5 inits × 6 metrics) |
| `fig_logconcave_omega_sweep_tau_zero_n{N}_rho{R}` | ω sweep, normalised gap |
| `fig_logconcave_time_to_tol_heatmap_n{N}_rho{R}` | Time-to-1e-4 heatmap, 3 τ panels |
| `fig_logconcave_tau_speedup_heatmap_n{N}_rho{R}` | Speedup ratio T(τ)/T(τ=0) |
| `fig_logconcave_speedup_vs_chi_n{N}_rho{R}` | Scatter: initial χ vs τ speedup |

To plot for a specific (n, rho):
```bash
python scripts/plot_logconcave_results.py --n 5 --rho 5
```

### Validate reference optimum

```bash
python scripts/check_logconcave_reference.py
```

---

## Extended repository structure

```
├── configs/
│   ├── gaussian_target.yaml
│   └── logconcave_target.yaml        ← new
├── src/
│   ├── __init__.py
│   ├── dynamics.py                   (Gaussian, closed-form)
│   ├── lc_dynamics.py                ← new: general target, scipy.linalg.expm
│   ├── targets.py                    ← new: LogCoshTarget
│   ├── qmc_samples.py                ← new: Sobol QMC + push-forward
│   ├── reference_optimum.py          ← new: L-BFGS-B VI optimiser
│   ├── lc_initializations.py         ← new: initializations relative to a★
│   ├── lc_metrics.py                 ← new: metrics relative to a★
│   ├── lc_plotting.py                ← new: 5 log-concave figures
│   ├── initializations.py            (Gaussian, relative to N(0,I))
│   ├── metrics.py                    (Gaussian, KL divergence)
│   ├── plotting.py                   (Gaussian, 4 figures)
│   └── utils.py                      (shared SPD utilities)
├── scripts/
│   ├── run_gaussian_grid.py
│   ├── run_logconcave_grid.py        ← new
│   ├── plot_gaussian_results.py
│   ├── plot_logconcave_results.py    ← new
│   ├── make_summary_tables.py
│   └── check_logconcave_reference.py ← new
├── tests/
│   ├── test_gaussian_update.py
│   ├── test_metrics.py
│   └── test_logconcave.py            ← new (34 tests)
└── outputs/
    ├── gaussian_grid/
    └── logconcave_grid/
        ├── results_long.csv
        ├── summary.csv
        ├── reference_optimum.npz
        ├── reference_optimum_meta.json
        ├── target_metadata.json
        └── figures/
```

---

## Log-concave summary CSV columns

| Column | Description |
|--------|-------------|
| `n`, `rho`, `m_features`, `target_seed`, `sample_seed`, `K` | Target + sample identity |
| `omega`, `tau_type`, `tau_value`, `init_name`, `dt`, `T` | Run parameters |
| `final_objective_gap` | F(T) − F★ (raw) |
| `final_normalized_objective_gap` | (F(T)−F★) / gap₀ |
| `time_to_1e_minus_2/4/6` | First time normalised gap ≤ threshold (inf if not reached) |
| `monotone_objective_bool` | True if normalised gap is non-increasing |
| `min_eig_min_over_time` | Minimum eigenvalue of C across all saved steps |
| `max_eig_max_over_time` | Maximum eigenvalue of C across all saved steps |
| `initial_chi`, `final_chi` | Trace-dominance ratio at t=0 and t=T |
| `initial/final_volume_error` | Volume error at start and end |
| `initial/final_shape_error` | Shape error at start and end |
