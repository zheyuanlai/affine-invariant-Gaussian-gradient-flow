# Checkpoint 09: one-interface global completion

Date: 2026-07-16

## Status

No complete proof has yet passed audit. The reductions, lower-dimensional cases, approximation framework, and direct variational construction are no longer the bottleneck. The surviving bottleneck is a global theorem for one balanced Cheeger interface: either complete enough transverse mass into essentially disjoint cells, or exhibit a macroscopic competitor with a universal perimeter saving.

## 1. Direct-deficit theorem

For a Lipschitz witness F define

\[
 R(G)=\inf_c\int |G-c|\,d\mu,
 \qquad D_\psi(G)=|DG|_\mu-\psi R(G).
\]

Minimizing

\[
D_\psi(G)+\kappa\|M(G)-M(F)\|_*
\]

gives an exact minimizer with one median subgradient q and one compatible constant anisotropy H. Almost every level globally minimizes the same functional

\[
P_\Phi(B)-\psi\int_B q\,d\mu,
\qquad \Phi(v)=|v|+\langle H,v\rangle,
\]

and satisfies

\[
(1-\kappa)\psi m(v)\le P_\Phi(B)\le \psi m(v),
\qquad |\lambda|\le\psi.
\]

With the frozen parameters \(\kappa=10^{-6}\), \(\alpha=10^{-28}\), and \(\beta=10^{-14}\), all matrix losses before geometric completion are below \(3.664\cdot10^{-5}\). Conditional on regular completion, the remaining covariance coefficient is at least \(2.4524\cdot10^{-16}/\psi^3\).

The common objects are q and H. A common obstacle-free calibration does not follow: the obstacle multiplier is real and can carry the entire median jump.

## 2. Exact foliation limit

The compactness reduction is now stated with its indispensable uniform-integrability hypothesis. If the deficit tends to zero, the matrix moment converges, and median-centered approximate minimizers are uniformly integrable, then strict BV compactness produces a nonconstant exact minimizer whose almost-everywhere levels are exact Cheeger sets. Heat smoothing alone does not supply the missing uniform integrability, so this is an exact reduction rather than a closure theorem.

## 3. Diffuse levels versus the median jump

On a regular level contained in a region where q is locally constant, the forced minimizer is an unconstrained stable critical point of anisotropic perimeter plus or minus \(\psi\mu\). The Jacobi inequality with extrinsic cutoffs and the constant lapse forces vanishing second fundamental form and vanishing density curvature in the anisotropic normal direction. Consequently an exact diffuse regular band is flat and log-affine.

This mechanism does not touch a function with a single median jump. Every balanced Cheeger set can be encoded by such a two-valued function while preserving the normalized normal matrix exactly. The one-interface problem is therefore unavoidable.

## 4. Exact polyhedral closure

For an isotropic log-affine density on a convex polyhedron and a balanced polyhedral exact relative-perimeter minimizer, any nonflat interior ridge admits a chord/bevel replacement with first-order gain and only second-order volume repair. Hence no interior ridge survives. The remaining facets are full hyperplane slices.

Their associated smaller halfspaces form a central-cell covering, so their masses sum to at least one half. The one-dimensional isotropic log-concave halfline bound then yields

\[
P_\Phi(E)\ge \frac{\phi_-}{4\sqrt3},
\qquad
\psi\ge \frac{\phi_-}{2\sqrt3}.
\]

This is a genuine dimension-free theorem for the exact polyhedral/log-affine branch. Arbitrary triangulation does not transfer it: for a smooth mesh, aggregate bevel gain and approximation error are both quadratic in mesh size.

## 5. Audited local obstructions

The following candidate closures are now ruled out.

1. **Asymmetric quantile escape.** Off-center tent penalties lose \(4d/(1+2d)\) at displacement d from one half. Product-exponential max sets can be amplitude-scaled to match total variation and the full aggregate normal matrix while the cusp has smaller asymmetric deficit. Aggregate moment penalties cannot label levels.

2. **Jacobi or tube closure.** The matrix inequality \(Q\preceq I\), trace one, permits \(Q=I/n\). The isotropic radial measure with density proportional to \(\exp(-\sqrt{n+1}|x|)\) has a stable median sphere, \(Q=I/n\), bounded two-sided tube deficit, and dimension-free Cheeger constant. Thus local curvature, stability, tube, and high-rank data do not force the desired lower bound.

3. **Unweighted random-fiber incidence.** For the balanced coordinate cut of the isotropic cube, the random-direction fiber incidence is of order \(1/n\), despite a dimension-free Cheeger constant. The correct weighted fiber quantity involves inverse conditional scale; isotropy alone gives only a three-halves-power estimate and does not close the argument.

4. **Exact common calibration.** Diffuse max-of-affine switches create a singular divergence and are excluded, but a single jump can absorb the multiplier. The hoped-for Wulff-tube foliation therefore does not follow.

## 6. Live approach families

### A. Global calibrated incidence

Construct a measurable flow-cell or ray decomposition attached to the single interface. Each regular patch must be assigned transverse mass, with bounded reuse, so that either total assigned smaller-side mass is universal or cell collisions produce a competitor with universal saving. This must handle both flat cube cuts and spherical radial cuts.

### B. Spectral compatibility of localization or transport rays

Per-ray one-dimensional inequalities are exact, but their conditional scales can concentrate in many directions. The live task is to derive a global compatibility law from log-concavity and nonbranching, not merely from covariance trace bounds.

### C. Smooth-to-polyhedral completion

Seek a macroscopic chord replacement whose gain is linear in a geometrically meaningful collision or turning quantity while its volume repair is quadratic. Pure mesh refinement is blocked because both terms are quadratic.

## 7. Blocked list

- Any argument whose missing input is dimension-free control of the top covariance along stochastic localization.
- Thin shell alone, or a bootstrap retaining any logarithmic dimension loss.
- Pointwise or averaged Brascamp--Lieb inversion without a new degeneracy mechanism.
- Local Jacobi, tube, or normal-matrix estimates without a global incidence theorem.
- Aggregate matrix or total-variation penalties intended to prevent a median jump.
- Unweighted random-fiber incidence.
- Naive polyhedral approximation of a smooth interface.

## 8. Adversarial instantiations

- **Cube:** exact flat cut; defeats unweighted fiber incidence but satisfies completion.
- **Simplex:** flat cap behavior stresses small directional incidence.
- **Isotropic cross-polytope:** nonsmooth support and many competing facets stress approximation.
- **Product exponential:** max and cyclic-max interfaces have genuine bevel gains; aggregate moment penalties nevertheless fail to distinguish central and cusp encodings.
- **Radial exponential:** stable smooth sphere defeats every purely local curvature/tube completion premise.
- **Strongly asymmetric log-affine polyhedron:** covered by the exact polyhedral theorem.

## 9. Next audit gate

A candidate global theorem is admissible only if it states an explicit measurable assignment, proves bounded overlap without dimension dependence, and gives a dichotomy whose constants remain positive on the cube, simplex, product exponential, and radial exponential examples. Any proof that silently replaces this with a covariance trace estimate is rejected.
