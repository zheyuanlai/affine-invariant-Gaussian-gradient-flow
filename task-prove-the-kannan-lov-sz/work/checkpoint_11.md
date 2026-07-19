# Checkpoint 11: global uncertainty and transport-density rank

## 1. Candidate-proof status

There is still no complete proof.  No dimension-free conclusion is being
claimed.  This checkpoint records four exact completions, three decisive
no-go results, and two surviving global targets.  Every occurrence of a
universal constant below is independent of the dimension.

## 2. Median signed distance is an exact T3 formulation

For every nonpoint log-concave probability on its affine support,

\[
 D_1(\mu)=
 \sup_{\substack{E\text{ open}\\\mu(E)=1/2}}
       \int d(x,\partial E)d\mu(x).                       \tag{2.1}
\]

The proof in `signed_distance_extremizer.md` includes level atoms and zero
plateaus: a generic arbitrarily small linear tilt makes the pushforward law
atomless without changing the first moment.  Lower-dimensional support is
handled intrinsically.

For a smooth maximizer, write the normal-cell masses on the two sides of a
basepoint as \(a(y)\) and \(b(y)\).  The constrained first variation gives

\[
 b(y)-a(y)=c\rho(y).
\]

Integrating and using global half mass forces \(c=0\).  Hence every normal
Voronoi cell is individually bisected.  In a unique-nearest chart the exact
second variation is

\[
 J''(h)=2\int_\Sigma\rho h^2
 -\int_\Sigma\langle G_y\nabla h,\nabla h\rangle .        \tag{2.2}
\]

A codimension-one medial switch adds the nonpositive form

\[
 -\int_C{(\alpha_1-\alpha_2)^2\over
               |\nabla(d_1-d_2)|}\rho,d\mathcal H^{d-1}. \tag{2.3}
\]

It cannot be omitted.  Gaussian halfspaces, cube halfspaces, and the
isotropic radial exponential median sphere pass the formulas; the radial
model has degree-one equality.  An isotropic atomic regular simplex has a
global signed-distance maximum of order \(\sqrt d\), proving that isotropy
and global maximality without log-concavity are insufficient.  Equation
(2.1) is exact, but the desired upper bound on its right side remains T3
itself.

## 3. Smooth local chord aggregation is closed as a no-go

For disjoint trace-matched replacements, let

\[
 v_i=\Delta V_i,\qquad r_i=\Delta P_i-\lambda v_i.
\]

If the cells cancel volume, then exactly

\[
 \sum_i v_i=0\quad\Longrightarrow\quad
 \Delta P=\sum_i r_i.                                    \tag{3.1}
\]

Thus cancellation removes only the common CMC multiplier.  It cannot
reverse the nonnegative constrained residual of stable smooth cells.
Circular and spherical-cap residuals are strictly positive.  Cube and
quartic flat slices have zero residual.  Genuine angle ridges do yield a
negative first-order bevel residual with only a second-order volume error,
but this mechanism was already completed in the polyhedral branch.

The smooth survivor therefore requires a nonlocal change in adjacency or
reuse of material.  Summing more local chords is blocked.

## 4. Aggressive mass-preserving localization is stopped, not terminal

For

\[
 A=\operatorname {Cov}_{\mu_t}X,\quad
 v=\operatorname {Cov}_{\mu_t}(1_S,X),\quad
 q=A^{-1}v,\quad
 P=I-{q\otimes q\over|q|^2},
\]

the controller

\[
 D=A^{-1}PA^{-1}                                         \tag{4.1}
\]

satisfies

\[
 Dv=0,\qquad ADA=P.                                      \tag{4.2}
\]

The fixed set mass is pathwise constant, and the covariance has drift
\(-P\).  All stopped density, mean, label-centroid, covariance, inverse,
determinant, and perimeter Itô identities in
`aggressive_mass_preserving_localization.md` passed an independent
formula audit.  In particular

\[
 d\operatorname {tr}A_t=dM_t-(d-1)dt.                    \tag{4.3}
\]

This reaches the first covariance face in universal expected time, not a
rank-one posterior.  The protected direction rotates even for an
off-centre Gaussian ball; its exact initial angular quadratic-variation
rate is \((d-1)/r^2\).  In the Gaussian halfspace the natural clock
diverges logarithmically and the \(D\)-clock as the inverse terminal
eigenvalue.  A product-exponential maximum has angular rate of order
\((\log d)^2\).  There is no proved singular continuation, and no universal
survivor-variance estimate.  The controller is therefore blocked after
the stopped identity stage.

## 5. Random-line uncertainty is an exact global target

For a balanced Borel set, a random direction \(\theta\), projection
coordinate \(y\), conditional set mass \(p_{\theta,y}\), and conditional
line standard deviation \(\sigma_{\theta,y}\), define

\[
 \mathcal U_\mu(E)=\mathbb E_\theta\int
 {\min(p_{\theta,y},1-p_{\theta,y})\over\sigma_{\theta,y}}
 d(\pi_{\theta^\perp}\mu)(y).                            \tag{5.1}
\]

BV slicing, the sharp-scale one-dimensional log-concave Cheeger bound, and
Crofton give

\[
 {1\over2\sqrt3}\mathcal U_\mu(E)
 \le a_dP_\mu(E),\qquad
 a_d={\Gamma(d/2)\over\sqrt\pi\Gamma((d+1)/2)}
 \asymp d^{-1/2}.                                        \tag{5.2}
\]

Consequently the new assertion

\[
 \mathcal U_\mu(E)\ge {c\over\sqrt d}                   \tag{5.3}
\]

for every balanced set would close KLS.  The direction in (5.2) is the
opposite one, so (5.3) is not supplied by slicing.

The exact Blaschke--Petkantschin identity retains a constant cross-phase
pair mass but weights it by \(|s-t|^{d-1}\).  A direct conditional-moment
comparison introduces a fatal dimension-dependent factor.  All canonical
tests have the correct scale: Gaussian and cube halfspaces, a simplex cap,
the radial exponential median sphere, and the product-exponential maximum
give \(\Theta(d^{-1/2})\); a parity checkerboard is much larger.  No
counterexample is known, but (5.3) is a new variance-normalized hit-and-run
conductance theorem and is not yet proved.

## 6. Transport-density tensor and curvature pair

Let \(h=1_E-1_{E^c}\), let \(u\) be the oriented direction in an optimal
\(W_1\) ray decomposition between the two normalized halves, and let
\(\tau\) be half the Beckmann density:

\[
 -\operatorname {div}(\rho\tau u)=\rho h,
 \qquad D=\int\tau d\mu={1\over2}W_1(2\mu|_E,2\mu|_{E^c}). \tag{6.1}
\]

Raywise integration by parts gives the singularity-insensitive tensor
identity

\[
 M:=\int\tau u\otimes u,d\mu
   =\int hX\otimes u,d\mu.                              \tag{6.2}
\]

Isotropy implies

\[
 \|M\|_{HS}\le1,\qquad \operatorname {tr}M=D,
 \qquad \operatorname {rank}M\ge D^2.                   \tag{6.3}
\]

On smooth regular rays, Bochner and the Beckmann equation give the signed
identity

\[
 \int\tau\bigl(\|\nabla^2f\|_{HS}^2+
                    \nabla^2V(u,u)\bigr)d\mu
 =2\int_{\partial^*E}\langle u,n_{\rm in}\rangle
              dP_\mu
 \le2P_\mu(E),                                           \tag{6.4}
\]

provided the eikonal chart, endpoint traces, and integration by parts meet
the explicit hypotheses in `transport_density_tensor.md`.  The equality is
signed.  Replacing it by the positive boundary flux is false: an explicit
one-dimensional Gaussian set with three crossings has curvature budget
\(0.6372237\) while its positive boundary flux is \(1.3326089\).  Moreover,
the claim that arbitrary medial, potential-ridge, and support limits create
a finite nonnegative quadratic turning measure has not been proved.  A cube
spends the whole formal budget at its support, so the classical Hessian alone
is invalid.  The clean-room audit therefore validated (6.2)--(6.3), corrected
the sign in (6.4), and restricted its nonsmooth use.

Combining a tube lower bound \(D\ge c/P\) with (6.3) unconditionally forces
a bad cut to have effective direction rank at least \(c/P^2\).  The
additional \(CP^2\) normalized curvature conclusion is available only on
the audited smooth, monotone-crossing branch; it is not currently a theorem
for an arbitrary Cheeger cut.  Even on that branch, the missing global
theorem must exclude a spatially separated, high-rank collection of
almost-flat direction packets.  Applying an unproved Poincare inequality to
the transport-density law would be circular.

## 7. Active and blocked registry

### Active

1. Global compatibility from the audited transport-density tensor.  Any
   curvature use is restricted to the smooth signed-flux branch until a
   separate singular theorem is proved.
2. The variance-normalized random-line lower bound (5.3), with emphasis on
   a nonlocal uncertainty argument rather than another one-direction
   majority completion.
3. Normal-cell log-concavity: sharp one-dimensional cell scale identities
   combined with the global covariance and medial forms.
4. Boolean rigidity for the variance-normalized random-line form.

The normal-cell task completed immediately after this checkpoint was
written.  It proved the local log-concavity and the explicit comparisons

\[
 {1\over8r}\le\mathbb E|T|\le {1\over2r},\qquad
 {1\over48r^2}\le\operatorname {Var}T\le {1\over2r^2},
\]

including focal and support truncation.  It also proved that isotropy alone
aggregates these cells only at order \(\sqrt d\), and that the radial
exponential sphere saturates the remaining directional stability.  The
route is therefore retained as exact local infrastructure but removed from
the active completion list.

The simultaneous-polarization task also completed.  If \(\mu\) is already
reflection invariant, set-only two-point polarization preserves mass and
increases the signed-distance objective.  This does not extend to a
reduction for general log-concave measures: joint density/set sorting can
change half mass to \(3/4\), density polarization can destroy log-concavity
even after smooth strict-convex regularization, and Steiner symmetrization
followed by rebalancing and whitening has both monotonicity signs.  A
triangle cap loses a factor exactly \(1/\sqrt3\) after whitening.  The
symmetrization route is therefore closed outside already symmetric laws.

### Blocked

1. Disjoint smooth chord summation and scalar volume cancellation.
2. Quotient second variation without a global bridge; the radial exponential
   law is exactly sharp.
3. One-direction majority-fiber completion; max-flow reduces exactly to
   quasiminimality at the available error budget.
4. Aggressive mass-preserving localization beyond the first covariance
   face.
5. Pointwise support curvature to ruling, regular-chart Gauss Jacobians,
   and any omission of medial/support singular charges.
6. The direct high-moment use of the Radon identity.

## 8. Audit ledger

- Symbolic dimension tracking: every result above retains its explicit
  \(d\)-dependence; no dimension-free conclusion has been inferred.
- Degenerate support: (2.1) is intrinsic to the affine support; the new
  random-line and tensor statements are also formulated in that intrinsic
  dimension.
- Stochastic rigor: the aggressive controller is asserted only under the
  explicit bounded stops for which its Itô identities were audited.
- Singular geometry: medial, focal, potential-ridge, and support terms are
  separately flagged and are never replaced by almost-everywhere Hessians.
  No general nonsmooth curvature budget is currently asserted.
- Circularity: (5.3) and the missing global compatibility theorem after
  (6.4) are listed as unproved KLS-strength steps, not lemmas.
- Candidate proof: none yet; clean-room proof audit of a final chain has not
  begun.
