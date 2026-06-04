# Natural Gradient Local Convergence Rate — Notes

> Placeholder for the detailed mathematical write-up of the
> `natural_gradient_local_rate` experiment group. The tracked implementation
> source of truth is [`docs/specs/natural_gradient_local_rate_spec.md`](../docs/specs/natural_gradient_local_rate_spec.md).

## Question

Does the **local** convergence rate of the Gaussian natural gradient flow near
equilibrium depend on the dimension `N_theta`, or only weakly on the conditioning
`kappa` through `log(kappa)`? The conjecture is that the rate is essentially
dimension-free.

## Setup and notation

- `N_theta` — dimension.
- `rho_post(theta) ∝ exp(-V(theta))` — target (posterior).
- `rho_a = N(m, C)`, `a = (m, C)` — Gaussian variational state.
- `E(a) = KL(rho_a || rho_post)` — objective.
- `a_star = (0, I)` — equilibrium in equilibrium-whitened coordinates.

The Gaussian natural gradient flow (with `g(a) = -E_{rho_a}[grad V]`,
`H(a) = -E_{rho_a}[Hess V]`):

```
dm/dt = C g(a)
dC/dt = C + C H(a) C
```

Local perturbation `m = u`, `C = I + X` with `X` symmetric. The linearized
positive generator is `L_star = -J_star`:

```
L_star(u, X) = ( u + 0.5 T[X],  X + T_star[u] + 0.5 H_lin[X] )
T[X]       = E[ Tr(X Hess V(Z)) Z ]
T_star[u]  = E[ Hess V(Z) (Z^T u) ]
H_lin[X]   = E[ Hess V(Z) (Z^T X Z - Tr(X)) ],     Z ~ N(0, I)
```

`L_star` is self-adjoint with respect to the Fisher–Rao inner product at
equilibrium, `<(u1,X1),(u2,X2)> = u1·u2 + 0.5 Tr(X1 X2)`.

## Key estimates

- `gamma_loc` = smallest eigenvalue of `L_star` — the local convergence rate
  (primary numerical estimate).
- `Lambda_hat` = largest eigenvalue of `H_lin` — the quantity behind the
  dimension-free conjecture.

Reference theoretical bounds (constants centralized in
`src/natural_gradient_local_rate/diagnostics.py`):

```
Gamma_current        = min( beta_target - 1,
                            N_theta * (3 + 4/sqrt(pi) * (1 + log(kappa_target))) )
current_bound_rate   = 1 / (4 + Gamma_current)
D_kappa              = (4 + 4/sqrt(pi)) * (1 + log(kappa_target))
conjecture_bound_rate = 1 / (4 + D_kappa)
```

## Validation

The Riemannian exponential-map discretization with small `Delta t`:

```
m_{n+1} = m_n + Delta_t * C_n * g_n
C_{n+1} = C_n^{1/2} exp( Delta_t * [ I + C_n^{1/2} H_n C_n^{1/2} ] ) C_n^{1/2}
```

initialized as a small perturbation along the slow eigenvector of `L_star`. The
fitted decay rate of `R^2 = ||m||^2 + 0.5 ||C - I||_F^2` should match `gamma_loc`.
This flow simulation validates the linearized rate; it is not the primary
evidence.

## Estimator corrections and diagnostics (CPU pass)

A pilot run exposed a problem: the matrix-free operator-norm estimate
`Lambda_hat` grew with dimension and the linearized rate `gamma_loc` went
**negative** at `N_theta = 64` even for *separable* potentials, which are
theoretically dimension-free. Before trusting any coupled-case conclusion the
estimator has to pass separable sanity checks. This pass (CPU only — no GPU)
adds the corrections and diagnostics below. The reported quantities and their
status are in the per-row CSV schema (`results_long.csv`); the figures are in
`outputs/natural_gradient_local_rate/figures/estimator_diagnostics/`.

### Why finite-sample self-adjointness matters

In exact expectation `H_lin[X] = E[Hess V(Z) (Zᵀ X Z − Tr X)]` is self-adjoint on
symmetric matrices under the Frobenius inner product. The plain Monte-Carlo
estimator `H_forward[X] = mean_j Hess V(Z_j) (Z_jᵀ X Z_j − Tr X)` is **not**: on a
finite bank its Frobenius adjoint is the *different* operator
`H_adjoint[Y] = mean_j Tr(Y Hess V(Z_j)) (Z_j Z_jᵀ − I)`. Feeding a non-self-adjoint
operator to `eigsh` (which assumes symmetry) silently corrupts the spectrum.

The fix is to estimate the explicitly symmetrized operator
`H_sym = ½(H_forward + H_adjoint)` on the *same* bank and use it everywhere
(`Lambda_hat`, the `L_star` covariance block, every eigensolve). The relative
self-adjointness error
`|⟨X, H[Y]⟩ − ⟨H[X], Y⟩| / max(1, |⟨X,H[Y]⟩|, |⟨H[X],Y⟩|)` drops from `O(1)` for
`raw_forward` to `~10⁻¹⁵` for `symmetrized` (figure `diag2`). `raw_forward` is
retained only as a diagnostic.

### Why Fisher–Rao metric scaling matters

`L_star` is self-adjoint with respect to the **local Fisher–Rao inner product**
`⟨(u₁,X₁),(u₂,X₂)⟩_⋆ = u₁·u₂ + ½ Tr(X₁X₂)`, not the naive Euclidean one. Handing
`eigsh` the operator in plain `(u, vec X)` coordinates presents a *non-symmetric*
matrix and yields a meaningless `gamma_loc`. We pack tangents in
metric-whitened coordinates `y = (u, sym_to_vec(X)/√2)` (with `sym_to_vec` a
Frobenius isometry that puts `√2` on the off-diagonals), so the Euclidean dot
product of packed vectors equals the Fisher–Rao inner product. In these
coordinates the `L_star` matrix is genuinely symmetric: the self-adjointness
error `self_adjoint_error_L_star` is `~10⁻¹⁵` (figure `diag6`). `gamma_loc` is
only reported when this check passes.

### Why a separable exact / diagonal benchmark is required

Even after symmetrization and correct metric scaling, the **full** `H_sym` acts
on the `N_theta(N_theta+1)/2`-dimensional space of symmetric matrices, so its
largest eigenvalue is biased upward by finite-sample spectral noise that grows
with that dimension. This is visible in the cleanest possible case — the
**Gaussian** target, where the true `Lambda = 0` and `gamma_loc = 1` — at fixed
`M_mc = 32768`:

| `N_theta` | `Lambda_hat_full_sym` | `Lambda_hat_diag` | `max_i A_ii` | `gamma_loc` |
|---:|---:|---:|---:|---:|
| 4  | 0.03 | 0.01 | 0.002 | 0.98 |
| 16 | 0.29 | 0.11 | 0.018 | 0.87 |
| 32 | 0.75 | 0.22 | 0.024 | 0.66 |
| 64 | 2.06 | 0.41 | 0.039 | 0.02 |

The full-operator estimate and `gamma_loc` are destroyed at `N_theta = 64`
purely by noise; the diagnostic that isolates this is the comparison against two
lower-noise references:

- **diagonal-restricted** `A = G − 11ᵀ`, `A_ij = E[Hess V_ii(Z)(Z_j² − 1)]`, an
  `N_theta × N_theta` operator (much smaller domain, much less noise); for
  separable potentials `A` is exactly diagonal in expectation, and `max_i A_ii`
  is the cleanest estimator;
- **separable exact**, `max_i E[V_i''(Z)(Z² − 1)]` by Gauss–Hermite quadrature,
  which has *no* Monte-Carlo error and is the dimension-free ground truth for
  separable controls (`NaN` for non-separable families).

If `Lambda_hat_full_sym` grows with dimension while the diagonal estimate tracks
the exact benchmark (small `diag_minus_exact`, small `diag_offdiag_norm`), the
full-operator growth is an estimator artifact, not a property of the flow.

### Why the pilot high-dimensional results must not yet be interpreted

The **sample-size scaling** sweep settles it. For a separable control the exact
benchmark is `M_mc`-independent, so any `M_mc`-dependence of `Lambda_hat_full_sym`
is pure estimator noise. In the smoke sweep (`separable`, `N_theta = 16`,
`kappa = 5`; exact `≈ 0.117` at every `M_mc`):

| `M_mc` | `Lambda_hat_full_sym` | `Lambda_hat_diag` | exact | `gamma_loc` |
|---:|---:|---:|---:|---:|
| 1024 | 1.41 | 0.51 | 0.117 | 0.30 |
| 2048 | 1.08 | 0.35 | 0.116 | 0.45 |
| 4096 | 0.66 | 0.15 | 0.116 | 0.63 |
| 8192 | 0.47 | 0.13 | 0.117 | 0.71 |

`Lambda_hat_full_sym` falls and `gamma_loc` climbs toward 1 as `M_mc` grows, with
the exact value flat throughout (figures `diag4`, `diag5`). The high-dimensional
"growth" and the negative `gamma_loc` in the pilot are therefore consistent with
**finite-sample spectral noise**, not a genuine dimension dependence. Drawing a
dimension-dependence conclusion — for separable controls or, a fortiori, for
general coupled potentials — requires `M_mc` large enough that this sweep has
plateaued and the diagonal/exact benchmarks agree. That production sweep has not
been run here; this pass only fixes and validates the CPU estimator.

## Results

_To be filled in after production runs._ Smoke runs only exercise the code paths;
they are **not** production evidence, and (per the diagnostics above) no
dimension-dependence conclusion should be drawn until the sample-size scaling has
plateaued and the separable diagonal/exact benchmarks agree.
