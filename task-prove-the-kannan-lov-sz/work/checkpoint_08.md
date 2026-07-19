# Checkpoint 08: constant-anisotropy transfer, short Wulff amplification, and the direct Cheeger-deficit functional

## 1. Status

No complete KLS proof is asserted.  Since Checkpoint 07, the
near-minimizer transfer has advanced in three ways.

1. The spatial heat selector has been removed from the Euler equation.  An
   equimeasurable matrix-fidelity minimizer retains the **unweighted**
   coarea normal matrix and has one constant, uniformly elliptic anisotropy.
2. In the smooth positive-lapse branch, a near-Cheeger CMC leaf need not be
   globally isoperimetric.  Secants of its own concave nested perimeter
   control its multiplier, and a central-safe short Wulff excursion gives a
   fixed covariance amplification.
3. A still stronger candidate functional minimizes the direct Cheeger
   deficit itself, rather than imposing equimeasurability.  Its formal
   multiplier is prescribed in `[-psi,psi]`, including exact divided
   differences on jump interfaces.  This may remove the surrounding-
   foliation hypothesis altogether.

The first two items have complete smooth algebra and audited constant
chains.  The third is under a clean-room nonsmooth audit.  Two obstacles
remain load bearing even if it passes:

* singular/Cantor/free-boundary Wulff coverage for the retained matrix; and
* the balanced high-rank branch.  Local flatness, zero Jacobi charge, and
  small killed flux do not imply a global affine product.

## 2. Mass-preserving controlled localization

For a fixed witness set `S`, current covariance `A`, and

\[
 v=\operatorname {Cov}_{\mu_t}(1_S,X),\qquad
 w=A^{-1/2}v,\qquad q=v^TA^{-1}v,
\]

the exact rank-`n-1` control is

\[
 D=A^{-1/2}\left(I-{ww^T\over q}\right)A^{-1/2}.       \tag{2.1}
\]

It satisfies

\[
 Dv=0,\qquad
 dA\big|_{dt}=-A+{vv^T\over q},                       \tag{2.2}
\]

so the set mass is pathwise fixed and every whitened covariance mode except
one is damped.  The fixed-set perimeter is a nonnegative local martingale.

The discontinuous convention needed at `v=0` is not well posed.  In one
dimension, a Gaussian and a symmetric interval give an occupation-time
contradiction: “full control at zero and zero control off zero” need not
have a continuous weak solution.  Only mesh-relaxed controls are presently
valid.

Conditional on an actual terminal line disintegration, covariance closes
every fixed direction cap, and exact common terminal barycenters are
impossible in dimension at least two.  Thin shell gives the spectral
counting estimate

\[
 N_K(a)\le C a^{-9/2},
 \qquad K=\mathbb E[\sigma^2u\otimes u],               \tag{2.3}
\]

but only the dimension-dependent consequence
`E sigma^2<=C n^(7/9)`.  Product exponentials realize long dispersed
coordinate needles.  The multiscale small-eigenvalue terminal regime is
therefore still blocked.

## 3. Near-minimizers, fillings, and signed-distance reach

For the scalar proximal problem

\[
 \min_{\mu(B)=\mu(A)}
       P_\mu(B)+\Lambda\mu(A\mathbin\triangle B),      \tag{3.1}
\]

one has

\[
 P(B)\le P(A),\qquad
 \mu(A\triangle B)\le {\delta\over\Lambda},
 \qquad |H_\mu-\lambda|\le\Lambda                  \tag{3.2}
\]

on the regular boundary.  This does not preserve the quadratic normal
matrix.  Scalar divergence calibration controls only the first normal
moment; an explicit rotating divergence-free field destroys the proposed
projector transfer.  Adding a spatially selected matrix penalty preserves
the matrix but introduces the nonperturbative term `kappa grad omega`.

The exact phase-cut cost is

\[
 R_{cut}=2\operatorname {Fill}.                       \tag{3.3}
\]

Continuous max-flow/min-cut duality and corrected-flux Pythagoras hold, but
a thin tube around an arbitrarily turning embedded curve has an exact unit
flow and a minimum cross-section.  Thus flow alone cannot force parallel,
concurrent, or radial geometry.  Likewise, a spherical-cap minimizer in a
ball refutes an unconditional turning/Plateau estimate: a fabricated normal
phase trace can have positive Plateau excess while the isoperimetric and
bevel deficits are zero.

For a true T3 extremizer, a two-sided calibrated ray core has the sharp
reach inequality

\[
 |z-z'|^2-r^2|N-N'|^2
 \ge2r|\langle z-z',N+N'\rangle|,                    \tag{3.4}
\]

hence `|z-z'|>=r|N-N'|`.  Smooth charts satisfy

\[
 \|S\|_{op}\le r^{-1},\qquad
 \|S\|_{HS}^2\le C_\beta r^{-2}.                    \tag{3.5}
\]

Covariance excludes coherent packets, and translated thin shell excludes
exact concurrence.  The remaining packet has effective normal rank at
least `c r^2`; separated flat/product-like charts remain the obstruction.

## 4. Unweighted matrix retention

Let `F` be the centrally clipped heat function and define

\[
 M(F)=\int \sigma_F\sigma_F^T\,d|DF|_\mu,
 \qquad T=\operatorname {tr}M(F).                    \tag{4.1}
\]

The selected physical matrix is a soft thinning of this unweighted coarea
matrix.  Total variance gives

\[
 1-\operatorname {tr}Q^2>.004577971,
 \qquad Q={M(F)\over T}.                              \tag{4.2}
\]

For `0<kappa<1/3`, minimize among functions equimeasurable with `F`

\[
 \operatorname {TV}_\mu(G)
 +\kappa\|M(G)-M(F)\|_*.                             \tag{4.3}
\]

The direct method is valid on the relative interior of an arbitrary convex
support.  Nuclear/operator duality produces one constant tension

\[
 \Phi_H(\xi)=|\xi|+\kappa{\xi^TH\xi\over|\xi|},
 \qquad \|H\|_{op}\le1,                              \tag{4.4}
\]

with tangent ellipticity

\[
 (1-3\kappa)I\preceq D^2\Phi_H|_{\xi^\perp}
                   \preceq(1+3\kappa)I.              \tag{4.5}
\]

If

\[
 \Delta_F=\operatorname {TV}(F)
       -\int I_\mu(\mu(F>r))dr,
\]

then every minimizer obeys

\[
 \operatorname {TV}(G)\le\operatorname {TV}(F),
 \qquad \|M(G)-M(F)\|_*\le{\Delta_F\over\kappa}.    \tag{4.6}
\]

For normalized matrices `Q,Q'`, the sharp dimension-free estimate is

\[
 \|Q'-Q\|_*\le {\Delta_F\over T}(1+\kappa^{-1}),
 \qquad
 \operatorname {tr}(Q'^2)-\operatorname {tr}(Q^2)
 \le e+{e^2\over2}.                                  \tag{4.7}
\]

The audited hierarchy

\[
 \kappa=10^{-6},\qquad \alpha=10^{-28},qquad
 \beta=10^{-14}                                      \tag{4.8}
\]

retains central trace larger than `.0032827p`.  Matrix retention and the
Euclidean-to-Wulff direction conversion leave a fixed angular-variance
seed.  The required rank-500 threshold is about `1.81*10^47`; below this
fixed number Buser--Ledoux already supplies a universal conclusion.

## 5. What BV stationarity actually gives

There is one compatible constant nuclear subgradient for every
`mu`-preserving BV inner variation.  Equivalently, the entire coarea
varifold is stationary in aggregate.  This does **not** imply that almost
every threshold has one CMC multiplier.

For a finite-phase `SBV` function with values `a_i`, phase pressures `p_i`,
and an interface between phases `i,j`, the exact law is

\[
 H_{ij}={p_j-p_i\over a_j-a_i}.                       \tag{5.1}
\]

A threshold can contain interfaces with different divided differences.
An explicit two-disk configuration is stationary and locally stable while
one threshold has two different curvatures.  Moreover an indicator of an
isoperimetric set has zero integrated profile deficit and all of its normal
matrix on the jump part.  Therefore no estimate

\[
 \operatorname {tr}M^j(G)\le C\Delta_G              \tag{5.2}
\]

can hold without an additional hypothesis.

On a genuinely smooth noncritical band, exact quantile restoration does
give

\[
 -\operatorname {div}_\mu D\Phi_H(\nabla G)=\lambda(G). \tag{5.3}
\]

The remaining transfer cannot be called routine BV regularity; it must be a
no-contact theorem, a deficit charge for incompatible jumps, or a layered
tube theorem.

## 6. Jacobi charge and Wulff tubes

For a smooth positive-lapse CMC foliation, the Wulff-gauge lapse satisfies

\[
 -\lambda'(v)\int_{\Sigma_v}{1\over f}d\sigma_\Phi
 =\int\langle A_\Phi\nabla\log f,\nabla\log f\rangle d\sigma_\Phi
  +\int q_\Phi d\sigma_\Phi
  +\int b_\Phi d\tau_\Phi.                         \tag{6.1}
\]

Every term on the right is nonnegative.  Thus the chosen foliation's
perimeter is concave.  Under two-sided lapse control this gives a
dimension-free projector-drift estimate.  Integrated profile deficit alone
does not supply the lapse or slope hypotheses: nested moving-center caps in
a ball have zero profile deficit and still turn through curvature and
support contact.

For one regular hypersurface, put

\[
 z=D\Phi(N),\qquad B=D^2\Phi(N)|_{N^\perp}.
\]

Before first Wulff cut, focal time, or support contact,

\[
 j_x(t)=\det(I+tBS)
 e^{-V(x+tz)+V(x)}=e^{\lambda t-D_x(t)},             \tag{6.2}
\]

with

\[
 D_x(t)=\int_0^t(t-s)
 \left\{\operatorname {tr}[((I+sBS)^{-1}BS)^2]
       +\nabla^2V(x+sz)[z,z]\right\}ds\ge0.          \tag{6.3}
\]

The exact covariance output is for `z`, not `N`:

\[
 \operatorname {Cov}(\mu)\succeq
 {e^{-|\lambda|T-h}T^3\over12}
 \int_G z\otimes z\,d\sigma_{\Phi,\mu}.             \tag{6.4}
\]

There is no Loewner comparison between `z tensor z` and `N tensor N`; a
two-dimensional off-diagonal matrix gives an explicit counterexample.
Projective variance is stable, and for `kappa<=10^-3`

\[
 \|Q_z-Q_N\|_*\le21\kappa.                          \tag{6.5}
\]

## 7. Near-Cheeger CMC leaves need not be profile minimizers

Let `c m(v)` be a global anisotropic Cheeger lower bound.  If a smooth
nested CMC foliation has concave perimeter `P(v)` and one leaf satisfies

\[
 c m(v_0)\le P(v_0)\le(1+\delta)c m(v_0),            \tag{7.1}
\]

then secants over
`[v_0-rho m(v_0),v_0+rho m(v_0)]` give

\[
 |\lambda(v_0)|\le(1+\delta/\rho)c.                 \tag{7.2}
\]

Choose the leaf or its complement so that the outward multiplier is
nonpositive, and increase its volume by `gamma m(v_0)`.  The first tube time
`T` satisfies

\[
 {\gamma\over(1+\delta)c}\le T
 \le{\gamma\over(1-\gamma)c},                      \tag{7.3}
\]

while

\[
 {P_0-R(T-)\over P_0}
 \le{\delta+\gamma\over1+\delta}.                   \tag{7.4}
\]

This works on either arm and at the cusp `v=1/2`.  A tiny fixed `gamma`
still gives a universal multiple of `1/c`.

The heat construction supplies the correct direct Cheeger deficit.  For
the equimeasurable minimizer,

\[
 \int[P_E(G_r)-\psi m(v_r)]dr\le D_{co}(F_0).        \tag{7.5}
\]

Levels with relative gap above `epsilon` have total Euclidean surface trace
at most `D_co/epsilon`.  In addition,

\[
 \int_0^1|\mu(F_0>r)-1/2|dr
 \le\|F_0-1_S\|_1=2U(s),                            \tag{7.6}
\]

so all but negligible retained trace lies near half volume.

At

\[
 \kappa=10^{-6},\quad\alpha=10^{-28},\quad
 \epsilon=10^{-5},\quad h_0=10^{-2},
\]

the retained volumes lie in `[.49,.51]`, the anisotropic relative error is
`1.2000123*10^-5`, the short tube has length at least
`1.99997*10^-5/c`, and full bad-ray deletion is below
`8.3677*10^-5`.  The covariance coefficient is larger than
`2.4523*10^-16/c^3`.  The aggregate angular packet remains nonzero.

## 8. Direct Cheeger-deficit minimization

The following candidate removes the unknown value law and may remove the
foliation secants.  Minimize over `0<=G<=1`

\[
 \mathcal K(G)=D_\psi(G)+\kappa\|M(G)-M(F)\|_*,      \tag{8.1}
\]

where

\[
 D_\psi(G)=\operatorname {TV}(G)
 -\psi\int_0^1m(\mu(G>r))dr\ge0.                    \tag{8.2}
\]

The distribution term has the exact median representation

\[
 \int_0^1m(\mu(G>r))dr
 =\inf_{a\in\mathbb R}\int|G-a|d\mu.                \tag{8.3}
\]

Comparison with `F` formally gives the decisive pair

\[
 D_\psi(G)\le D_\psi(F),\qquad
 \|M(G)-M(F)\|_*\le{D_\psi(F)\over\kappa}.          \tag{8.4}
\]

The range constraint and (8.2) give BV coercivity.  On a smooth band, one
compatible constant matrix should give

\[
 -\operatorname {div}_\mu D\Phi_H(\nabla G)
 =\psi\,m'(v(G))\in[-\psi,\psi].                    \tag{8.5}
\]

Equivalently, away from a median the forcing is `psi sign(G-a)`.  On an
`SBV` jump, the curvature is the divided difference of `|s-a|`, and hence
also lies in `[-psi,psi]`.  A short Wulff tube only needs this pointwise
bound: even if different boundary components have different multipliers,
the factors `e^{lambda_x t}` stay between `e^{-psi t}` and `e^{psi t}`.

This is not yet certified.  The clean-room audit is checking lower
semicontinuity, simultaneous choice of the nuclear and median
subgradients, nonunique medians, atoms, the Cantor part, and whether the
first variation yields enough anisotropic quasiminimality for singular tube
coverage.  Until that audit passes, (8.5) is a candidate rather than a
lemma.

## 9. The high-rank obstruction is global

The proposed implication

\[
 \text{high flat rank + small local charge}
 \Longrightarrow\text{global affine product}         \tag{9.1}
\]

is false.  For

\[
 \Omega_{m,a}=\{x_i\ge0,\ x_i+x_{i+1}\ge a\},
 \qquad d\mu\propto e^{-\sum_i x_i}1_{\Omega_{m,a}}dx, \tag{9.2}
\]

the leaves `A_L={max_i x_i>=L}` have flat equal-CMC facets, zero Jacobi
charge, exact natural contact, and

\[
 Q(L)={I_m\over m}.                                   \tag{9.3}
\]

Their normalized killed loss to any fixed distance tends to zero as
`L-a-log m` tends to infinity.  Nevertheless the irredundant facet normals

\[
 e_i,\qquad e_i+e_{i+1}                              \tag{9.4}
\]

form a connected circuit and prevent every nontrivial affine product
decomposition.  The construction remains isotropic and uniformly elliptic
after a small cyclic conditioning and whitening.

The counterexample has vanishing tail flux, so it does not refute a
balanced global-saturation theorem.  At the exact median max level of the
product exponential measure,

\[
 p_m={m\over2}(2^{1/m}-1)
 \in[\tfrac12\log2,\tfrac12],                       \tag{9.5}
\]

and both inward and outward swept masses over a short tube are of the
expected fixed order.  This model has a universal perimeter, so it is a
calibration rather than a counterexample to KLS.

The correct remaining invariant must include global swept-flux coverage of
support normals and Hessian couplings.  Observed phase rank alone cannot
detect mixed support faces outside the tube.

## 10. Updated registry

### Active

1. Clean-room proof or refutation of the direct-deficit Euler and
   quasiminimality theorem, including jump/Cantor/free-boundary parts.
2. A layered/subgraph or finite-family construction if the direct
   functional does not supply pointwise curvature on its singular part.
3. Singular Wulff coverage for bounded-mean-curvature anisotropic BV
   boundaries, with exact convergence of the displacement matrix.
4. Balanced high-rank saturation: turn fixed swept mass plus small
   curvature/contact/collision loss into either a perimeter competitor or a
   globally controlled factor decomposition.
5. Mesh-relaxed mass-preserving localization, only if its terminal
   small-eigenvalue regime can be coupled to the same balanced-saturation
   invariant.

### Blocked

* arbitrary levelwise global minimization for an equimeasurable minimizer;
  the rectangle switch example disproves it;
* inference of levelwise CMC from aggregate BV stationarity;
* fixed `kappa=1/4` in a near-lossless Wulff profile comparison;
* Loewner comparison of Euclidean normals with Wulff displacements;
* rank-only product rigidity, local Jacobi rigidity, or flow-only filling
  rigidity;
* hard feedback at `v=0` in controlled stochastic localization;
* covariance or thin shell alone in the dispersed high-rank branch.

## 11. Audit ledger

* The sharper normalized-matrix constant and every retuned heat constant
  were independently recomputed.
* The Wulff Jacobian, hard-support edge convention, projective comparison,
  and covariance output were clean-room checked.
* The false Loewner assertion and the false levelwise-minimizer assertion
  were removed.
* The smooth near-Cheeger amplification was tested on Laplace, Gaussian,
  cube, ball-cap, one-sided exponential, and product-exponential models.
* The BV jump law and the cyclic support-coupling counterexample prevent
  the two main regularity/product shortcuts from re-entering implicitly.
* No dimension-free Poincare inequality, KLS-equivalent covariance-process
  bound, selected-boundary capacity, or log-Sobolev estimate is assumed.
* A full symbolic `n`-tracking audit and clean-room reproof of the eventual
  load-bearing singular and high-rank lemmas remain impossible until those
  lemmas are actually proved.
