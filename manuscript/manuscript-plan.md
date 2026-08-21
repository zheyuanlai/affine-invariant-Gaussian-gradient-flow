# Manuscript plan: sharp Fisher--Rao Gaussian variational inference

Status: integrated full draft complete.  Remaining author-side inputs are the
title-page metadata and the decision whether to rerun selected experiments to
restore the missing run-level provenance files documented in Appendix E.

## Recommended working title

**Sharp Global-to-Local Theory for Fisher--Rao Gaussian Variational Inference**

Possible subtitle, if the stochastic contribution needs to be visible in the title:

**Deterministic and Price/Hessian Stochastic Algorithms**

The shorter title is preferable. The abstract and introduction can make the stochastic theory prominent without turning the title into a list.

## Central question and answer

The paper asks which geometric and statistical mechanisms control the convergence of full-covariance Gaussian variational inference under the Fisher--Rao natural gradient.

The answer is a three-scale rate landscape:

1. **Covariance initialization:** underdispersion creates a covariance burn-in whose Fisher--Rao cost depends on the initial covariance.
2. **Global localization:** strong log-concavity gives global convergence at the target-conditioning scale, and matching constructions show that this dependence cannot generally be improved.
3. **Local non-Gaussian stiffness:** after optimizer whitening, the local rate is governed by a score-operator spectral scale rather than the initial covariance; a growing-dimensional ridge family makes the universal bound sharp.

The stochastic Price/Hessian STL theory preserves the same geometric structure, with explicit variance floors, a decreasing-step floor-free result, and a finite-horizon local non-exit theory. The affine-invariant metric classification explains why Fisher--Rao is the balanced canonical choice. The non-log-concave section marks the boundary of the theory: the two unprojected schemes can be unstable, while a covariance-constrained KL proximal step admits a corrected stationarity guarantee.

## Scope boundary

Included:

- full-covariance Gaussian reverse-KL variational inference;
- the complete affine-invariant metric family on Gaussians, up to an overall scale;
- Fisher--Rao flow, Riemannian retraction, and KL/Bregman resolvent;
- deterministic global, local, and global-to-local theory;
- all three lower constructions: spiral, bump train, and convex ridge;
- Price/Hessian STL algorithms and their global and local stochastic theory;
- a concise Bures--Wasserstein warm start and comparison result;
- non-log-concave instability and covariance-constrained stationarity;
- focused numerical illustrations and reproducibility details.

Excluded:

- a general Wasserstein--Fisher--Rao theory or WFR algorithmic study;
- an extension of every global/local/stochastic rate theorem to every affine-invariant metric;
- an unrestricted dimension-free logarithmic local-rate conjecture;
- an unconditional infinite-horizon stochastic local-containment claim.

The general affine-invariant family is used to classify the geometry, derive its flow, identify the exact trace/traceless preconditioning, and prove that Fisher--Rao is the unique balanced choice on Gaussian targets. The sharp convergence theory remains centered on Fisher--Rao.

## Main-text architecture

### 1. Introduction and rate landscape

- Define the Gaussian reverse-KL problem in one paragraph.
- Explain affine invariance and the three bottlenecks.
- State the deterministic and stochastic contributions in two compact theorem tables.
- Distinguish throughout between the extrinsic condition number `kappa` and optimizer-whitened `kappa_star`.
- State the scope of every sharpness result: initialization class, dimension, target family, and whether the metric is energy gap or local parameter error.
- Relate the results to natural-gradient VI, Bures--Wasserstein Gaussian VI, STL estimators, and constrained nonconvex mirror descent.

### 2. Gaussian VI and affine-invariant geometry

#### 2.1 Objective and differential identities

- Gaussian family `N(m,C)`, reverse-KL energy, first-order conditions, averaged score `g`, averaged Hessian `A`, and residual `R=C^{-1}-A`.
- Price/Bonnet identities and the regularity convention used later.

#### 2.2 Classification of affine-invariant metrics

State and prove, up to an overall positive multiple,

```text
g_(m,C)((u,X),(v,Y))
 = eta u^T C^{-1}v
 + omega Tr(C^{-1}X C^{-1}Y)
 + tau Tr(C^{-1}X) Tr(C^{-1}Y),
```

with `eta>0`, `omega>0`, and `tau>-omega/d`. Normalize `eta=1` after the classification.

The proof is an isotropy argument at `(0,I)`: the tangent representation splits into mean, trace, and traceless covariance modes; the sign element `-I` eliminates mean--covariance cross terms.

#### 2.3 Why Fisher--Rao is balanced

- Give the general affine-invariant covariance flow and its exact trace/traceless residual decomposition.
- On a Gaussian target in optimizer-whitened coordinates, the linear rates are
  `1/eta` for the mean, `1/(2 omega)` for traceless covariance, and `1/[2(omega+d tau)]` for trace covariance.
- Conclude that, up to scale, Fisher--Rao (`eta=1`, `omega=1/2`, `tau=0`) is the unique metric that balances all three modes.
- Make clear that these are exact Gaussian/local modal rates, not general non-Gaussian convergence rates.

#### 2.4 Canonical Fisher--Rao dynamics and algorithms

- Fisher--Rao metric, gradient norm, and exact energy dissipation.
- Continuous flow.
- Riemannian exponential retraction scheme, explicitly called a retraction rather than the full Gaussian-manifold exponential map.
- KL/Bregman resolvent.
- Brief oracle taxonomy: exact expectations versus Price/Hessian minibatches.

### 3. Sharp global deterministic theory under strong log-concavity

#### 3.1 Continuous covariance bootstrap and intrinsic upper bound

- Precision Duhamel identity and covariance floor.
- Optimizer whitening and the affine-invariant rate

```text
Delta_t <= Delta_0 [1 + beta_star lambda_(0,star)(e^t-1)]^{-1/kappa_star}.
```

- Separate the covariance burn-in term from the global localization term.

#### 3.2 Continuous lower bound

- State the spiral family and the `Omega(kappa_star)` localization time in the main text.
- Explain its scope: dimension at least two, anisotropic admissible initialization, and a growing initial radius.
- Put the full Hessian certificate, rescaled four-dimensional shadow system, semigroup estimate, and Gaussian-smoothing bootstrap in an appendix.

#### 3.3 Deterministic discretizations

- Present the two covariance bootstraps side by side.
- Give the Riemannian one-step descent theorem and the KL forward-divergence descent theorem.
- State the common global complexity

```text
O(h^{-1} log_+(1/(beta lambda_0))
  + (kappa/h) log_+(Delta_0/delta)).
```

- Explain the different computational and stability behavior: the Riemannian update uses a matrix exponential; the KL update is a rational resolvent and has a wider covariance-stability range.

#### 3.4 Discrete lower bound

- State the one-dimensional bump-train construction and the `Omega(kappa/h)` constant-factor lower bound for both fixed-step maps.
- Emphasize that it starts after covariance burn-in and has a dimension-free Hessian-Lipschitz constant.

#### 3.5 Removing only the covariance burn-in

- Give the Bures--Wasserstein warm start and deterministic quadratic rescue as concise corollaries.
- Do not develop general WFR dynamics.
- State explicitly that these devices remove initialization cost but do not remove the sharp global-conditioning or local-stiffness obstructions.

### 4. Local spectral theory and global-to-local complexity

#### 4.1 Optimizer whitening and the linearized generator

- Define the first- and second-Hermite score operators.
- Use the newer quantity

```text
Gamma = min{beta_star-1,
            sqrt(d)(4 sqrt(2 t_0)+4 t_0+4)}
       = O(1+sqrt(d) log kappa_star).
```

- State the spectral sandwich

```text
max{alpha_star,1/(4+Gamma)}
 <= gamma_star <= Lambda_star
 <= min{beta_star,(4+Gamma)/2}.
```

#### 4.2 Gaussian-core nonlinear local region

- Calibrate first with Gaussian targets: exact continuous local threshold and exact KL invariant interval.
- Introduce energy coercivity, the covariance band, the non-Gaussian modulus, and the certified flow region.
- Be precise that the general discrete radius is constructive/existential through a supremum of second map derivatives; do not call it a closed-form practical radius.

#### 4.3 Deterministic local maps and online entry

- State local contractions for both maps.
- Give the computable residual gate and exact Gaussian spectral gates.

#### 4.4 Sharp local obstruction

- Integrate the repaired convex-ridge construction.
- State the combined conclusion

```text
gamma_star = Theta(kappa_star^{-1/2}),
tau_H^2    = Theta(kappa_star^{1/2})
```

on a family with `d=Theta(kappa_star^2)` and uniformly bounded Hessian-Lipschitz constant.
- Conclude that the universal `1/(4+Gamma)` branch is sharp up to constants on this family.
- State that this does not resolve fixed dimension or sharpness of the `sqrt(d) log kappa_star` branch.

#### 4.5 Deterministic global-to-local theorem

- Combine covariance burn-in, global localization, and the local contraction into a three-stage complexity statement.
- Keep the error metric transition explicit: global energy gap first, equilibrium parameter norm locally.

### 5. Price/Hessian stochastic Fisher--Rao algorithms

#### 5.1 Algorithm and covariance safety

- Define the joint gradient/Hessian oracle and the STL mean estimator.
- State quadratic rescue and the pathwise covariance bands.
- Explain that using the same sample for score and Hessian is allowed; the analysis does not assume independence between those two components.

#### 5.2 Intrinsic STL variance

- Introduce the intrinsic Hessian fluctuation `Psi`.
- Give both regimes:
  - self-bounding under strong log-concavity and smoothness alone;
  - Hessian-Lipschitz control, with anisotropic/effective-dimension form before the worst-case `d^2` bound.
- Highlight exact vanishing for Gaussian Hessians without overstating pathwise zero noise away from matched covariance.

#### 5.3 Global stochastic convergence

- Riemannian and KL constant-step results in one comparison theorem/table.
- Hessian-Lipschitz regime: `tilde O(kappa^2)` rounds to an explicit floor.
- assumption-minimal regime: `tilde O(kappa^3)` single-sample behavior, or `tilde O(kappa^2)` rounds with batch `Theta(kappa)` and `tilde O(kappa^3)` oracle pairs.
- Decreasing stepsizes: floor-free `O(1/N)` in the Hessian-Lipschitz regime.

#### 5.4 Local stochastic convergence

- State the local variance decomposition and perturbation lemmas.
- Present the local theorem honestly as a killed/non-exit mean-square contraction plus a separate finite-horizon exit-probability bound.
- With the universal spectral certificate, the local iteration bound is

```text
O((4+Gamma)^2 max{1,V_1/B} log(r_0^2/epsilon)),
```

not an unqualified linear-in-`log kappa` rate.
- State that high-probability containment over horizon `N` requires deep entry and a batch proportional to that horizon under the present argument.
- Give the stochastic global-to-local corollary with these qualifications visible in the theorem statement.

#### 5.5 Gaussian target calibration

- Quadratic rescue recovers the exact optimizer in one query.
- With matched covariance, later STL evolution is pathwise deterministic.
- With unmatched covariance, the additive floor vanishes but mean-score noise remains until the covariance matches.

### 6. Beyond log-concavity

#### 6.1 Analytic instability of the unprojected schemes

- Use the smooth double-well family

```text
V_R(x) = x^2/2 - 2 R^2 log cosh(x/R),
```

whose Hessian lies in `[-1,1]`.
- Prove a Riemannian covariance cascade and a KL pole, rather than presenting them only as numerical observations.

#### 6.2 Bures--Wasserstein comparison

- State the existing nonconvex running-minimum stationarity result for the Bures--Wasserstein forward--backward method and use it as a geometry comparison.
- Keep this concise; it is not a WFR section.

#### 6.3 Covariance-constrained KL stationarity

- Retain the current split mean/covariance algorithm and spectral covariance clipping.
- Use the corrected one-step stationarity measure

```text
S_n = (1/2)||m_(n+1)-m_n||^2_(C_n^{-1})
      + KL(N(0,C_n) || N(0,C_(n+1))).
```

- For `h <= 1/(2 beta lambda_+)`, prove

```text
min_{n<N} S_n <= (h/N)(E_0-E_N).
```

- Deduce the full reverse-Gaussian-KL version with the factor `lambda_+/lambda_-`.
- In the centered scalar experiment, the mean term vanishes, so the plotted factor-one envelope remains valid.
- Do not include the legacy stochastic nonconvex clipping theorem unless it is separately rederived.

### 7. Numerical illustrations

Use a small number of composite figures, each tied to a theorem or boundary claim:

1. **Three-stage deterministic behavior:** covariance burn-in, global localization, local acceleration, and Bures warm start.
2. **Riemannian versus KL discretization:** same-step speed, larger-step robustness, and computational distinction; avoid calling the retraction the exact Gaussian exponential map.
3. **STL behavior:** estimator variance versus distance, batch scaling, and constant-step floors.
4. **Affine metric modes:** Gaussian trace/traceless calibration plus one non-Gaussian illustration, interpreted as mode preconditioning rather than a general rate theorem.
5. **Non-log-concave boundary:** Riemannian cascade, KL pole, Bures running-minimum comparison, and clipped-KL constrained stationarity.

Move secondary heatmaps and parameter sweeps to a numerical appendix. Exclude the general WFR experiment suite from this paper.

### 8. Discussion

- What is genuinely sharp and what is only sufficient.
- Fixed-dimensional local-rate question.
- Closed-form discrete local region and sharper stochastic local stepsizes.
- General affine-invariant metrics beyond Gaussian modal analysis.
- General WFR as the natural second-paper direction.

## Appendix architecture

- **A.** Gaussian differential identities, metric classification, and affine covariance.
- **B.** Spiral Hessian and shadowing certificates.
- **C.** Deterministic one-step and global upper-bound details.
- **D.** Bump-train construction.
- **E.** Linear score operators and non-Gaussian local modulus.
- **F.** Exact Gaussian continuous and discrete regions.
- **G.** Convex-ridge local counterexample and smooth realization.
- **H.** Stochastic global proofs and decreasing stepsizes.
- **I.** Stochastic local perturbation, killed process, and containment.
- **J.** Non-log-concave instability and clipped-KL stationarity.
- **K.** Experimental protocols, metadata, and supplementary plots.

## Proof-readiness ledger

| Result package | Status | Manuscript treatment |
|---|---|---|
| Affine-invariant metric classification and Fisher--Rao balance | Ready after concise rewrite | Main theorem in Section 2 |
| Continuous global upper bound | Ready | Main theorem |
| Spiral continuous lower bound | Repaired and standalone-compiling | Main statement, full proof in Appendix B |
| Riemannian/KL deterministic global upper bounds | Ready | Main theorem pair |
| Bump-train fixed-step lower bound | Ready | Main statement, Appendix D |
| New `sqrt(d)` local spectral theorem | Ready | Main theorem |
| Gaussian-core nonlinear local region | Ready, with explicitness qualification | Main theorem |
| Convex-ridge counterexample | Repaired and standalone-compiling | Main sharpness corollary, Appendix G |
| Price/Hessian global stochastic theory | Audited, corrected, and standalone-compiling | Main Section 5 |
| Stochastic local theory | Audited; correct only in non-exit/finite-horizon form | Main theorem with prominent scope |
| Nonconvex instability | Analytic proof available | Main Section 6 |
| Clipped KL stationarity | Legacy full-KL statement needs corrected metric/factor | Corrected theorem in Section 6 |
| General WFR | Out of scope | Second paper |

## Source-assembly plan after approval

Do not continue growing a single monolithic source. Create a master file and modular inputs:

```text
manuscript/main.tex
manuscript/sections/01-introduction.tex
manuscript/sections/02-geometry.tex
manuscript/sections/03-global-deterministic.tex
manuscript/sections/04-local.tex
manuscript/sections/05-stochastic.tex
manuscript/sections/06-nonconvex.tex
manuscript/sections/07-numerics.tex
manuscript/sections/08-discussion.tex
manuscript/appendices/*.tex
manuscript/figures/*
manuscript/refs.bib
```

Use one namespaced label convention and one notation table before importing proofs. The old `natural-gradient.tex` is evidence and source material, not a trusted merge base. The improved deterministic and stochastic notes supply the primary proofs; the repaired counterexample and audited reports supply the extensions.

The journal currently accepts a submission PDF and asks for its EMS LaTeX template only after acceptance, so the drafting source can remain a clean portable `article` setup until the mathematics and exposition stabilize.
