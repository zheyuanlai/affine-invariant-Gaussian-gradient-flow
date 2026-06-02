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

## Summary CSV columns

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
