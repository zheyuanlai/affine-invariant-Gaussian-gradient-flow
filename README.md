# Affine-Invariant Gaussian Gradient Flows: Discretization and Numerical Mode Decomposition

This repository studies a two-parameter family of affine-invariant Gaussian
gradient flows for the variational problem

$$
\min_{m,\,C}\; F(m,C),\qquad
F(m,C)=\operatorname{KL}\!\big(\mathcal N(m,C)\,\Vert\,\pi\big),
$$

which, for a target $\pi(\theta)\propto e^{-V(\theta)}$ and up to additive
constants, equals

$$
F(m,C)=\mathbb E_{\theta\sim\mathcal N(m,C)}[V(\theta)]-\tfrac12\log\det C .
$$

The flows are parameterized by $(\omega,\tau)$. The experiments isolate four
convergence modes — mean, covariance-volume, covariance-shape, and mixed — and
quantify the effect of $\omega$ and $\tau$ on each. Two targets are studied: the
exact Gaussian target $\mathcal N(0,I)$, and a strongly log-concave non-Gaussian
target with Monte Carlo / quasi-Monte Carlo expectations and a numerically
computed reference Gaussian optimum. 

---

## 1. Variational problem

Let the variational family be the non-degenerate Gaussians

$$
q_{m,C}=\mathcal N(m,C),\qquad m\in\mathbb R^n,\quad C\in\mathbb S_{++}^n,
$$

and let the target be $\pi(\theta)\propto\exp(-V(\theta))$ with $V$ smooth. Define
the expected gradient and Hessian under the variational distribution,

$$
g(m,C)=\mathbb E_{q_{m,C}}[\nabla V(\theta)],\qquad
S(m,C)=\mathbb E_{q_{m,C}}[\nabla^2 V(\theta)].
$$

The first-order stationarity conditions of $F$ are

$$
g(m_\star,C_\star)=0,\qquad C_\star^{-1}=S(m_\star,C_\star).
$$

---

## 2. The $(\omega,\tau)$ family of flows

In the notation used by the code, the continuous-time affine-invariant Gaussian
flow is

$$
\dot m_t=-\,C_t\,g(m_t,C_t),
$$

$$
B_t=C_t^{1/2}\,S(m_t,C_t)\,C_t^{1/2},\qquad
\alpha_t=\frac{\omega+\tau\operatorname{Tr}(B_t)}{\omega+n\tau},
$$

$$
\dot C_t=\frac{1}{2\omega}\,C_t^{1/2}\big(-B_t+\alpha_t I\big)C_t^{1/2}.
$$

The matrix $B_t$ is the Hessian expectation in the natural (whitened) coordinates
of $C_t$, and $\alpha_t$ is a scalar trace-weighting term. The admissibility
condition is

$$
\omega>0,\qquad \omega+n\tau>0,
$$

which keeps $\alpha_t$ finite and the covariance equation well-posed.

---

## 3. Exponential-map discretization

With step $\Delta t$, the implemented update is

$$
m_{k+1}=m_k-\Delta t\,C_k\,g(m_k,C_k),
$$

$$
C_{k+1}=C_k^{1/2}\exp\!\Big(\frac{\Delta t}{2\omega}\big(-B_k+\alpha_k I\big)\Big)C_k^{1/2}.
$$

Because $-B_k+\alpha_k I$ is symmetric, its matrix exponential is symmetric
positive definite, and $C_{k+1}$ is a congruence transform
$C_k^{1/2}(\,\cdot\,)C_k^{1/2}$ of an SPD matrix; hence $C_{k+1}\in\mathbb S_{++}^n$
for every step. This is the discrete exponential map on the SPD manifold and
removes the positivity constraint without projection.

---

## 4. Trace-mode linearization: role of $\omega$ and $\tau$

The parameter $\omega$ rescales the covariance equation uniformly. The parameter
$\tau$ acts only through the scalar $\alpha_t$, i.e. only on the trace (volume)
component of the covariance residual. This can be made precise by linearizing
near stationarity.

Consider an isotropic perturbation of the whitened Hessian,
$B=(1+\varepsilon)I$. Then $\operatorname{Tr}(B)=n(1+\varepsilon)$ and

$$
\alpha=\frac{\omega+\tau n(1+\varepsilon)}{\omega+n\tau}
=1+\frac{n\tau}{\omega+n\tau}\,\varepsilon,
$$

so that

$$
-B+\alpha I=\Big(-1-\varepsilon+1+\tfrac{n\tau}{\omega+n\tau}\varepsilon\Big)I
=-\frac{\omega}{\omega+n\tau}\,\varepsilon\,I.
$$

Substituting into the covariance equation and evaluating near $C=I$ gives a
linear decay $\dot\varepsilon=-\kappa\,\varepsilon$ for the trace mode with rate

$$
\kappa=\frac{1}{2(\omega+n\tau)} .
$$

Consequences for the trace/volume mode:

- $\tau=0$: rate $\kappa=\dfrac{1}{2\omega}$.
- $\tau=-\dfrac{\omega}{2n}$: $\omega+n\tau=\dfrac{\omega}{2}$, so $\kappa=\dfrac{1}{\omega}$ — the rate doubles, i.e. the convergence time halves.
- $\tau=+\dfrac{\omega}{2n}$: $\omega+n\tau=\dfrac{3\omega}{2}$, so $\kappa=\dfrac{1}{3\omega}$ — the time scale increases by a factor $3/2$.

This linearization applies to trace/volume modes only. The mean equation does not
contain $\tau$, and traceless (shape) covariance modes have
$\operatorname{Tr}(R_{\mathrm{cov}})=0$, on which $\alpha$ has no leading-order
effect. The numerical results in §7 and §10 are consistent with these factors.

---

## 5. Covariance-error diagnostics

Let $R_{\mathrm{cov}}$ denote the covariance residual: for the Gaussian target it
is taken as $I-C$ (equivalently the eigenvalue residuals $1-\lambda_i$), and for
the general target as

$$
R_{\mathrm{cov}}=I-B,\qquad B=C^{1/2}S(m,C)\,C^{1/2}.
$$

Decompose it into trace and traceless parts:

$$
r_{\mathrm{tr}}=\frac{|\operatorname{Tr}(R_{\mathrm{cov}})|}{\sqrt n},
\qquad
r_{\mathrm{tf}}=\Big\|\,R_{\mathrm{cov}}-\tfrac{\operatorname{Tr}(R_{\mathrm{cov}})}{n}I\,\Big\|_F .
$$

The trace-dominance ratio is

$$
\chi=\frac{\big(\operatorname{Tr}R_{\mathrm{cov}}\big)^2}{n\,\|R_{\mathrm{cov}}\|_F^2},
\qquad 0\le\chi\le 1 .
$$

$\chi\approx 1$ indicates a residual that is almost purely trace/volume;
$\chi\approx 0$ indicates an almost purely traceless/shape residual. By the
linearization of §4, $\tau$ can only accelerate the trace component, so $\chi$ is
the diagnostic used to interpret which initializations admit a $\tau$-speedup.

Other recorded scalars (per saved step) include the objective gap, the whitened
mean error $\|C_\star^{-1/2}(m-m_\star)\|$, the relative covariance error
$\|C_\star^{-1/2}CC_\star^{-1/2}-I\|_F/\sqrt n$, the volume error
$|\log\det(C_\star^{-1/2}CC_\star^{-1/2})/n|$, the traceless shape error, the
covariance residual $\|I-B\|_F$, and the eigenvalue extremes of $C$. Column names
are documented in §13.

---

## 6. Gaussian target

For

$$
V(\theta)=\tfrac12\|\theta\|^2,\qquad \pi=\mathcal N(0,I),
$$

the expectations are exact:

$$
g(m,C)=m,\qquad S(m,C)=I.
$$

Hence $B=C$, and the covariance update is diagonal in the eigenbasis of $C$;
writing $C=Q\operatorname{diag}(\lambda_i)Q^\top$,

$$
\lambda_{i,k+1}=\lambda_{i,k}\exp\!\Big(\frac{\Delta t}{2\omega}\big(-\lambda_{i,k}+\alpha_k\big)\Big),
\qquad
\alpha_k=\frac{\omega+\tau\sum_j\lambda_{j,k}}{\omega+n\tau}.
$$

The stationary point is $m_\star=0$, $C_\star=I$, and the objective reduces to the
exact KL energy $F=\tfrac12(\|m\|^2+\operatorname{Tr}C-\log\det C-n)$.

**Initializations (mode decomposition).** With $\mathbf 1$ the all-ones vector,
$r$ a magnitude, and $\operatorname{diag}(\cdot)$ a diagonal covariance:

| Mode | $m_0$ | $C_0$ |
|---|---|---|
| mean-only | $\dfrac{r}{\sqrt n}\mathbf 1,\;r=3$ | $I$ |
| volume-high | $0$ | $4I$ |
| volume-low | $0$ | $0.25\,I$ |
| shape-only | $0$ | $\operatorname{diag}(e^{r},e^{-r},1,\dots,1),\;r=2$ |
| mixed | $\dfrac{2}{\sqrt n}\mathbf 1$ | $2\,\operatorname{diag}(e^{r},e^{-r},1,\dots,1),\;r=1.5$ |

The shape-only covariance has $\det C_0=1$, so its residual is purely traceless;
the volume initializations are purely isotropic ($\chi=1$); the mixed
initialization perturbs mean, volume, and shape simultaneously.

The grid is $n\in\{2,5,10\}$, $\omega\in\{0.125,0.25,0.5,1,2\}$, and for each
$\omega$ the three values $\tau\in\{-\omega/2n,\,0,\,+\omega/2n\}$, with
$\Delta t=0.02$, $T=20$.

---

## 7. Gaussian numerical results

Results below are for $n=5$, extracted from `outputs/gaussian_grid/summary.csv`.
Time-to-tolerance is the first time the normalized energy $F/F_0$ falls below the
threshold; `>20` denotes not reached within the horizon $T=20$.

**Time to $F/F_0\le 10^{-4}$, $\tau=0$.**

| Initialization | $\omega=0.125$ | $\omega=0.25$ | $\omega=0.5$ | $\omega=1$ | $\omega=2$ |
|---|---|---|---|---|---|
| mean-only | 4.6 | 4.6 | 4.6 | 4.6 | 4.6 |
| volume-high | 0.9 | 1.9 | 3.7 | 7.5 | 15.0 |
| volume-low | 1.4 | 2.8 | 5.6 | 11.2 | >20 |
| shape-only | 1.3 | 2.6 | 5.3 | 10.5 | >20 |
| mixed | 4.0 | 4.0 | 4.3 | 7.3 | 14.5 |

**$\tau$-speedup ratio $T(\tau)/T(0)$ at $\omega=0.5$** (tolerance $10^{-4}$):

| Initialization | $T(\tau{=}0)$ | $T(\tau_-)/T_0$ | $T(\tau_+)/T_0$ |
|---|---|---|---|
| mean-only | 4.6 | 1.00 | 1.00 |
| volume-high | 3.7 | 0.51 | 1.51 |
| volume-low | 5.6 | 0.50 | 1.50 |
| shape-only | 5.3 | 0.98 | 1.06 |
| mixed | 4.3 | 1.09 | 1.07 |

**$\tau$-speedup ratio across $\omega\in\{0.25,0.5,1\}$** (tolerance $10^{-4}$):

| Initialization | $\tau_-,\,\omega{=}0.25$ | $\tau_+,\,\omega{=}0.25$ | $\tau_-,\,\omega{=}0.5$ | $\tau_+,\,\omega{=}0.5$ | $\tau_-,\,\omega{=}1$ | $\tau_+,\,\omega{=}1$ |
|---|---|---|---|---|---|---|
| mean-only | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| volume-high | 0.47 | 1.47 | 0.51 | 1.51 | 0.49 | 1.49 |
| volume-low | 0.50 | 1.50 | 0.50 | 1.50 | 0.50 | 1.49 |
| shape-only | 1.00 | 1.08 | 0.98 | 1.06 | 0.99 | 1.07 |
| mixed | 1.05 | 0.95 | 1.09 | 1.07 | 1.07 | 1.22 |

**Observed numerical behavior.**

- The mean-only time-to-tolerance is $4.6$ for every $\omega$ and every $\tau$:
  the mean equation is independent of $\omega$ and $\tau$.
- The two volume initializations accelerate by a factor close to $2$ under
  $\tau_-$ (ratios $0.47$–$0.51$) and slow by a factor close to $3/2$ under
  $\tau_+$ (ratios $1.47$–$1.51$), uniformly in $\omega$. This matches the
  trace-mode rates $\kappa$ of §4.
- The shape-only ratios stay within $0.98$–$1.08$: $\tau$ does not substantially
  affect a purely traceless residual.
- The mixed ratios lie in $0.95$–$1.22$ with no consistent sign, since the
  limiting error combines mean and shape components that $\tau$ does not
  accelerate.
- Under $\tau=0$, covariance-dominated modes slow monotonically with $\omega$
  (e.g. volume-high from $0.9$ at $\omega=0.125$ to $15.0$ at $\omega=2$), while
  mean-only is flat, consistent with $\kappa\propto 1/\omega$.

---

## 8. Strongly log-concave non-Gaussian target

Define, with $m=4n$ fixed unit vectors $a_\ell\in\mathbb R^n$,

$$
V_\rho(x)=\tfrac12\|x\|^2+\frac{\rho}{m}\sum_{\ell=1}^{m}\log\cosh(a_\ell^\top x).
$$

Its gradient and Hessian are

$$
\nabla V_\rho(x)=x+\frac{\rho}{m}\sum_{\ell=1}^{m}\tanh(a_\ell^\top x)\,a_\ell,
$$

$$
\nabla^2 V_\rho(x)=I+\frac{\rho}{m}\sum_{\ell=1}^{m}\operatorname{sech}^2(a_\ell^\top x)\,a_\ell a_\ell^\top.
$$

Since $0\le\operatorname{sech}^2(\cdot)\le1$ and the $a_\ell$ are unit vectors,

$$
I\preceq\nabla^2 V_\rho(x)\preceq(1+\rho)I,
$$

so $V_\rho$ is $1$-strongly convex and $(1+\rho)$-smooth, giving a strongly
log-concave non-Gaussian target.

**Reference Gaussian optimum.** Expectations are not available in closed form. A
reference optimum $a_\star=(m_\star,C_\star)$ is computed by minimizing a
fixed-sample estimate of $F$ over a Cholesky factor $C=LL^\top$,

$$
F(m,L)=\frac1K\sum_{j=1}^{K}V_\rho(m+Lz_j)-\sum_{i=1}^{n}\log L_{ii},
\qquad z_j\sim\mathcal N(0,I),
$$

with analytic gradients

$$
\nabla_m F=\frac1K\sum_{j=1}^{K}\nabla V_\rho(m+Lz_j),\qquad
\nabla_L F=\frac1K\sum_{j=1}^{K}\nabla V_\rho(m+Lz_j)\,z_j^\top-L^{-\top},
$$

minimized by L-BFGS-B with common random numbers $\{z_j\}$. The same fixed sample
set is reused across all $(\omega,\tau,\text{init})$ runs so that comparisons share
identical Monte Carlo noise.

---

## 9. Log-concave numerical setup

The values below are read from `outputs/logconcave_grid/target_metadata.json` and
`outputs/logconcave_grid/reference_optimum_meta.json`.

| Quantity | Value |
|---|---|
| dimension $n$ | $5$ |
| coupling $\rho$ | $5.0$ |
| features $m=4n$ | $20$ |
| dynamics samples $K$ | $4096$ |
| reference samples $K_{\mathrm{ref}}$ | $8192$ |
| step $\Delta t$ | $0.005$ |
| horizon $T$ | $40$ |
| reference objective $F_\star$ | $3.980868$ |
| reference mean norm $\|m_\star\|$ | $2.92\times10^{-5}$ |
| reference gradient norm $\|\nabla_m F\|$ | $8.24\times10^{-10}$ |
| reference covariance residual $\|C_\star S_\star-I\|_F$ | $1.85\times10^{-3}$ |

For the production run reported here, the reference optimizer converged with
$F_\star=3.980868$, $\|\nabla_m F\|=8.24\times10^{-10}$, and covariance residual
$1.85\times10^{-3}$. Because $V_\rho$ is even, $m_\star$ is numerically zero
($\|m_\star\|=2.92\times10^{-5}$); it is optimized rather than assumed.

The grid is $\omega\in\{0.25,0.5,1\}$ with $\tau\in\{-\omega/2n,\,0,\,+\omega/2n\}$
over the same five initializations, now defined relative to $a_\star$ (e.g.
volume-high uses $C_0=4C_\star$, shape-only uses
$C_0=C_\star^{1/2}\operatorname{diag}(e^2,e^{-2},1,1,1)C_\star^{1/2}$).

---

## 10. Log-concave numerical results

Results are for $n=5,\ \rho=5$, extracted from
`outputs/logconcave_grid/summary.csv`. Time-to-tolerance uses the normalized
objective gap $(F-F_\star)/(F_0-F_\star)$.

**Objective-floor note.** The finite-sample objective gap has a numerical floor
near $10^{-4}$ for this run: the smallest observed final gap is
$\approx1.5\times10^{-4}$ (normalized $\approx1.0\times10^{-4}$). The $10^{-4}$
threshold therefore sits on the floor, and entries marked `>40` at that tolerance
(notably volume-low, whose normalized gap plateaus at $1.01\times10^{-4}$) should
not be read as divergence. The $10^{-2}$ tolerance is unaffected and is used for
the speedup ratios.

**Time to normalized gap $\le 10^{-2}$, $\tau=0$.**

| Initialization | $\omega=0.25$ | $\omega=0.5$ | $\omega=1$ |
|---|---|---|---|
| mean-only | 2.2 | 2.3 | 2.4 |
| volume-high | 0.9 | 1.9 | 3.6 |
| volume-low | 1.8 | 3.5 | 6.9 |
| shape-only | 1.4 | 2.9 | 5.7 |
| mixed | 1.8 | 2.0 | 3.1 |

**Time to normalized gap $\le 10^{-4}$, $\tau=0$** (at the objective floor):

| Initialization | $\omega=0.25$ | $\omega=0.5$ | $\omega=1$ |
|---|---|---|---|
| mean-only | 4.8 | 4.8 | 4.9 |
| volume-high | 2.4 | 4.7 | 9.3 |
| volume-low | >40 | >40 | >40 |
| shape-only | 3.0 | 5.9 | 11.7 |
| mixed | 4.2 | 4.6 | 8.1 |

**$\tau$-speedup ratio $T(\tau)/T(0)$ at $\omega=0.5$** (tolerance $10^{-2}$):

| Initialization | $T(\tau{=}0)$ | $T(\tau_-)/T_0$ | $T(\tau_+)/T_0$ |
|---|---|---|---|
| mean-only | 2.3 | 1.00 | 1.00 |
| volume-high | 1.9 | 0.51 | 1.49 |
| volume-low | 3.5 | 0.51 | 1.49 |
| shape-only | 2.9 | 1.00 | 0.98 |
| mixed | 2.0 | 1.10 | 1.02 |

**Initial trace dominance $\chi_0$ at $\omega=0.5$:**

| Initialization | $\chi_0$ |
|---|---|
| mean-only | 0.705 |
| volume-high | 0.998 |
| volume-low | 1.000 |
| shape-only | 0.133 |
| mixed | 0.295 |

**Observed numerical behavior.**

- Under $\tau=0$, the $\omega$-sweep slows covariance-dominated modes as $\omega$
  increases (volume-low $1.8\to6.9$, shape-only $1.4\to5.7$ as $\omega$ goes
  $0.25\to1$), while mean-only stays near $2.2$–$2.4$.
- The two volume initializations, which have $\chi_0\approx1$, accelerate by a
  factor close to $2$ under $\tau_-$ ($T(\tau_-)/T_0=0.51$) and slow by a factor
  close to $3/2$ under $\tau_+$ ($1.49$), matching the Gaussian case and the
  §4 rates.
- shape-only ($\chi_0=0.133$) is essentially unchanged under $\tau$
  (ratios $0.97$–$1.00$); mean-only ($\chi_0=0.705$, but mean-limited) is exactly
  $1.00$.
- mixed ($\chi_0=0.295$) shows ratios $1.06$–$1.10$ under $\tau_-$, i.e. no net
  speedup, because the limiting error is not trace dominated.
- Across the grid, the magnitude of the $\tau_-$ speedup tracks $\chi_0$: only
  the trace-dominated initializations realize the near-$2\times$ acceleration.

---

## 11. Figures

Figures are written by `scripts/plot_gaussian_results.py` and
`scripts/plot_logconcave_results.py` as paired PNG/PDF. Gaussian figures live in
`outputs/gaussian_grid/figures/n{n}/`; log-concave figures in
`outputs/logconcave_grid/figures/`.

**Main figures.**

- `fig_gaussian_tau_speedup_n{n}` / `fig_logconcave_tau_speedup_n{n}_rho{r}`:
  heatmaps of the observed ratio $T(\tau)/T(0)$, diverging colormap centered at
  $1$. Demonstrates the near-$2\times$ acceleration for volume modes under
  $\tau_-$ and the $\approx3/2$ slowdown under $\tau_+$.
- `fig_gaussian_omega_sweep_n{n}` / `fig_logconcave_omega_sweep_n{n}_rho{r}`:
  normalized objective gap versus time for $\tau=0$ across $\omega$, one panel per
  initialization. Demonstrates the covariance time-scale dependence on $\omega$.
- `fig_gaussian_tau_effect_n{n}` / `fig_logconcave_tau_effect_n{n}_rho{r}`:
  selected diagnostics for the informative modes, separating trace residual
  $r_{\mathrm{tr}}$ and traceless residual $r_{\mathrm{tf}}$ across the three
  $\tau$ values.
- `fig_logconcave_speedup_vs_chi_n{n}_rho{r}`: observed $\tau_-$ speedup against
  the initial trace dominance $\chi_0$, demonstrating that the speedup is realized
  only when $\chi_0\approx1$.

**Appendix / diagnostic figures.**

- `fig_gaussian_tau_effect_full_n{n}` / `fig_logconcave_tau_effect_full_n{n}_rho{r}`:
  the full diagnostic grid (all initializations $\times$ all recorded metrics).
- `fig_gaussian_time_to_tol_n{n}` / `fig_logconcave_time_to_tol_n{n}_rho{r}`:
  time-to-tolerance heatmaps; not-reached cells are masked and labeled `>T`.

---

## 12. Installation and reproducibility

Python 3.9+ with NumPy/SciPy (CPU, float64); no GPU or PyTorch. Determinism is
provided by fixed seeds for the target $a_\ell$ and for the QMC sample set.

```bash
pip install -r requirements.txt
pytest                                   # unit tests (127 tests)

python scripts/run_gaussian_grid.py      # Gaussian grid  -> outputs/gaussian_grid/
python scripts/plot_gaussian_results.py

python scripts/run_logconcave_grid.py    # log-concave grid -> outputs/logconcave_grid/
python scripts/plot_logconcave_results.py

python scripts/extract_readme_results.py # regenerate the tables in this file
```

The result tables in §7 and §10 are produced by
`scripts/extract_readme_results.py`, which only reads existing summary files and
prints Markdown; it does not run the dynamics. The grid runners accept
`--n`, `--dt`, `--T`, `--omega`, `--outdir`, and (log-concave) `--rho`, `--K`,
`--K-ref`, `--target-seed`, `--sample-seed`, `--force-optimize`.

---

## 13. Repository structure

```
configs/
  gaussian_target.yaml          Gaussian grid configuration
  logconcave_target.yaml        log-concave grid configuration
src/
  dynamics.py                   Gaussian eigenvalue update (exact)
  lc_dynamics.py                general update via scipy.linalg.expm
  targets.py                    LogCoshTarget (batched grad/Hess)
  qmc_samples.py                Sobol QMC sample generation
  reference_optimum.py          L-BFGS-B reference Gaussian optimum
  initializations.py            Gaussian initializations
  lc_initializations.py         log-concave initializations (relative to a_star)
  metrics.py                    Gaussian metrics (KL energy, residuals)
  lc_metrics.py                 log-concave metrics (objective gap, residuals, chi)
  plotting.py / lc_plotting.py  figure generation
  utils.py                      SPD utilities, parameter validation
scripts/
  run_gaussian_grid.py          run Gaussian grid -> CSV
  run_logconcave_grid.py        run log-concave grid -> CSV + reference optimum
  plot_gaussian_results.py      Gaussian figures
  plot_logconcave_results.py    log-concave figures
  make_summary_tables.py        recompute summary.csv from results_long.csv
  check_logconcave_reference.py inspect the cached reference optimum
  extract_readme_results.py     print the Markdown result tables in this README
tests/
  test_gaussian_update.py, test_metrics.py, test_logconcave.py
outputs/
  gaussian_grid/, logconcave_grid/   results_long.csv, summary.csv, figures/
```

### Recorded summary columns

Gaussian (`outputs/gaussian_grid/summary.csv`): `n, omega, tau_type, tau_value,
init_name, dt, T, final_energy, final_normalized_energy, time_to_1e_minus_2,
time_to_1e_minus_4, time_to_1e_minus_6, monotone_energy_bool,
min_eig_min_over_time, max_eig_max_over_time`.

Log-concave (`outputs/logconcave_grid/summary.csv`): the same run identifiers plus
`rho, m_features, target_seed, sample_seed, K`, the objective-gap fields
`final_objective_gap, final_normalized_objective_gap`, the tolerance times,
`monotone_objective_bool`, the eigenvalue extremes, and the diagnostic fields
`initial_chi, final_chi, initial_volume_error, initial_shape_error,
final_volume_error, final_shape_error`.
