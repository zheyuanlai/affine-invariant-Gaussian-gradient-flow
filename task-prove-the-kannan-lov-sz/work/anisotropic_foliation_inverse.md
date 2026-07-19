# Constant-anisotropy CMC foliations: a Jacobi charge identity and a quantitative turning inverse

## 0. Outcome and exact limitation

This note analyzes the revised joint-replacement setup from
`equimeasurable_matrix_replacement.md`.  There is a useful exact theorem,
but the hypotheses currently available do **not** yet imply the requested
global parallel/concurrent/product classification.

The useful theorem is the following.  In a smooth nested foliation by
weighted constant-anisotropic-mean-curvature leaves, let `f>0` be the lapse
between neighboring leaves, expressed in the natural Wulff-normal gauge.
Then

\[
 -\lambda'(v)\int_{\Sigma_v}{1\over f}\,d\sigma_\Phi
 =\int_{\Sigma_v}\langle A_\Phi\nabla\log f,
              \nabla\log f\rangle\,d\sigma_\Phi
  +\int_{\Sigma_v}q_\Phi\,d\sigma_\Phi
  +\int_{\partial\Sigma_v}b_\Phi\,d\tau_\Phi .       \tag{0.1}
\]

Here `lambda(v)` is the constant weighted anisotropic mean curvature,
`A_Phi` is uniformly positive definite, and

\[
 q_\Phi=
 \operatorname {tr}\!\left[(D^2\Phi(N)S)^2\right]
       +\nabla^2V[D\Phi(N),D\Phi(N)]\ge0.              \tag{0.2}
\]

The boundary density `b_Phi` is the first-contact term and is nonnegative
for a convex support with the natural anisotropic free-boundary condition.
Formula (0.1) is the infinitesimal counterpart of the killed Wulff-tube
formula.  It charges, with one sign, all three smooth ways in which a CMC
foliation can turn:

* variation of the lapse along a leaf;
* anisotropic second fundamental form or convexity of the potential; and
* contact with the hard support.

At a critical value, a focal time, or a collision, the smooth identity stops
and the missing flux is exactly the killed-ray charge.  Thus representative
normal directions cannot rotate for free.

There is also a dimension-free quantitative consequence.  If on a volume
interval the lapse satisfies

\[
 {\ell\over P(v)}\le f\le {L\over P(v)},              \tag{0.3}
\]

the topology is fixed, `|lambda|/P` is bounded, and the CMC slope drop is
small, then the normalized normal-projector matrix `Q_v` satisfies

\[
 \int\|Q_v'\|_{HS}^2,dv
 \le C_\Phi\left\{
 {L^2\over\ell}{\lambda(a)-\lambda(b)\over P_{\min}}
 +L\sup {|\lambda|^2\over P^2}\int W(v)\,dv\right\}, \tag{0.4}
\]

where `W(v)=1-tr(Q_v^2)` is the within-leaf projector variance.  Hence a
family of coherent leaves with a nearly affine perimeter profile has one
essentially fixed representative direction.  Any change of direction is
paid by one of the terms in (0.1), failure of (0.3), or a
critical/contact/focal/collision event.

The limitation is equally precise.  The present assumption

\[
 \int [P(\{G>r\})-I(\mu\{G>r\})],dr\ll1             \tag{0.5}
\]

does not by itself bound the slope drop in (0.4), does not give (0.3), and
does not control topology changes.  Exact rotating spherical caps in a
Euclidean ball have zero profile deficit and nevertheless turn through
support contact; a monotone reparametrization of their common level
function can assign essentially arbitrary coarea weight to two separated
directions while preserving zero deficit.  This is the concurrent/contact
branch, not a counterexample to the intended trichotomy, but it proves that
contact cannot be estimated from (0.5).  A complete inverse still needs a
dimension-free argument which turns the heat-generated value law and joint
minimality into either the regular-lapse hypotheses of (0.4) or a
quantitatively large killed/contact/product charge.

## 1. Projector variance and its invariant dichotomy

Let `Sigma_r=partial{G>r}` on a regular central band and initially use
Euclidean weighted area `d sigma=e^{-V}dH^{n-1}`.  Put

\[
 P_r=\sigma(\Sigma_r),\qquad
 Q_r={1\over P_r}\int_{\Sigma_r}N\otimes N\,d\sigma,
 \qquad W_r=1-\operatorname {tr}(Q_r^2).              \tag{1.1}
\]

For the coarea probability

\[
 d\pi(r)={P_r,dr\over\int P_sds},
 \qquad \bar Q=\int Q_r,d\pi(r),                    \tag{1.2}
\]

the Hilbert-space law of total variance is

\[
 \boxed{
  1-\operatorname {tr}(\bar Q^2)
   =\int W_r,d\pi(r)
      +\int\|Q_r-\bar Q\|_{HS}^2d\pi(r).}            \tag{1.3}
\]

Thus the assumed lower bound `.003` is either within-leaf dispersion or
cross-leaf motion.  No regularity or choice of representative normal is
used in (1.3).

The constant anisotropy changes this only by a universal factor.  Let

\[
 d\sigma_\Phi=\Phi(N)d\sigma,
 \qquad 0<m\le\Phi(N)\le M.                           \tag{1.4}
\]

For any random Hilbert-space vector `Z`, the pair formula for variance gives

\[
 \operatorname {Var}_{\sigma_\Phi/P_\Phi}(Z)
 \ge {m^2\over M^2}
       \operatorname {Var}_{\sigma/P}(Z).             \tag{1.5}
\]

Indeed the density of each of two independent anisotropically reweighted
samples is at least `m/M` times the original density.  Consequently a fixed
Euclidean projector variance remains a fixed anisotropic projector
variance.  In what follows all constants may depend on the fixed numbers
`m,M` and the tangent ellipticity bounds of `Phi`, but never on dimension.

A useful elementary fact is

\[
 \lambda_{\max}(Q_r)\ge\operatorname {tr}(Q_r^2)=1-W_r. \tag{1.6}
\]

Thus `W_r` small means genuine projective coherence: for a top eigenvector
`u_r`,

\[
 {1\over P_r}\int_{\Sigma_r}
       \{1-(N\cdot u_r)^2\}\,d\sigma\le W_r.          \tag{1.7}
\]

The sign of `N` is immaterial, as it must be for disconnected planar
components.

## 2. The positive-lapse Jacobi identity

### 2.1 Euclidean formula, including the hard wall

Let `Omega` be a `C^2` convex domain, let `V in C^2(Omega)` be convex, and
let `dmu=Z^{-1}e^{-V}1_Omega dx`.  Suppose
`{Sigma_v}_{v in (a,b)}` is a `C^2` nested foliation, parametrized by
enclosed `mu`-volume, with no critical value or topology change.  Assume
every leaf meets `partial Omega` orthogonally and has constant weighted mean
curvature

\[
 H_\mu=\operatorname {tr}S-\nabla V\cdot N=\lambda(v). \tag{2.1}
\]

Choose the normal parametrization

\[
 \partial_vX=fN,
 \qquad f>0,
 \qquad \int_{\Sigma_v}f,d\sigma=1.                 \tag{2.2}
\]

The standard variation identities, with `S(Y)=D_YN`, are

\[
 \partial_vN=-\nabla_\Sigma f,
 \qquad
 \lambda'=-L_\Sigma f-qf,                            \tag{2.3}
\]

where

\[
 L_\Sigma=\Delta_\Sigma-\nabla_\Sigma V\cdot\nabla_\Sigma,
 \qquad q=|S|_{HS}^2+\nabla^2V(N,N)\ge0.             \tag{2.4}
\]

Preservation of the orthogonal contact condition gives the Robin condition

\[
 \partial_\eta f=\mathrm {II}_{\partial\Omega}(N,N)f
       \quad\hbox{on }\partial\Sigma_v,               \tag{2.5}
\]

where `eta=nu_Omega` is the outward conormal.  Convexity gives
`II_{partial Omega}>=0` with this convention.

**Lemma 2.1 (exact Jacobi charge).**  Under (2.1)--(2.5),

\[
 \boxed{
 -\lambda'(v)\int_{\Sigma_v}{1\over f}\,d\sigma
 =\int_{\Sigma_v}|\nabla_\Sigma\log f|^2d\sigma
  +\int_{\Sigma_v}q\,d\sigma
  +\int_{\partial\Sigma_v}
       \mathrm {II}_{\partial\Omega}(N,N)d\tau.}     \tag{2.6}
\]

**Proof.**  Divide the second identity in (2.3) by `f` and integrate.  The
weighted divergence theorem and (2.5) give

\[
 \int_{\Sigma_v}{L_\Sigma f\over f}d\sigma
 =\int_{\Sigma_v}|\nabla\log f|^2d\sigma
  +\int_{\partial\Sigma_v}{\partial_\eta f\over f}d\tau.
\]

Substitution proves (2.6).  Every term is finite under the stated smoothness
hypotheses and every term on the right is nonnegative.  QED.

The first variation of perimeter is

\[
 P'(v)=\lambda(v),                                    \tag{2.7}
\]

because (2.2) normalizes the volume speed.  Hence every such positive-lapse
CMC foliation has a concave perimeter function:

\[
 P''(v)=\lambda'(v)\le0.                              \tag{2.8}
\]

This conclusion uses the foliation, not merely the fact that each isolated
leaf is CMC.

### 2.2 Constant anisotropy

Let `Phi` be even, positive, convex, `C^3` off zero, and one-homogeneous.
Assume its spherical Hessian is uniformly elliptic.  Set

\[
 z(N)=D\Phi(N),\qquad G(N)=D^2\Phi(N)|_{N^\perp}.      \tag{2.9}
\]

Write a geometric deformation in Wulff gauge as

\[
 \partial_vX=fz(N).                                   \tag{2.10}
\]

Its normal volume speed is `f Phi(N)`, so volume parametrization is exactly

\[
 \int_{\Sigma_v}f,d\sigma_\Phi=1.                   \tag{2.11}
\]

The tangential part of (2.10) is not an error.  It gives the cancellation

\[
 D_vN=-\Phi(N)\nabla_\Sigma f.                        \tag{2.12}
\]

To verify (2.12), write `fz=fPhi N+fz_T`.  For a velocity `uN+T`,
`D_vN=-nabla u+ST`; Euler homogeneity gives
`nabla_Sigma Phi=S z_T`, and the two terms containing `fS z_T` cancel.

The linearization of weighted anisotropic mean curvature in this gauge is
self-adjoint for `d sigma_Phi`.  Denote it by

\[
 J_\Phi f=\operatorname {div}_{\sigma_\Phi}
                (A_\Phi\nabla_\Sigma f)+q_\Phi f.     \tag{2.13}
\]

Here `A_Phi` is symmetric and uniformly elliptic, with constants determined
only by `Phi`.  Evaluating the zero-order term on the constant Wulff offset
and using the exact Jacobian from Section 5 of
`equimeasurable_matrix_replacement.md` gives (0.2).  In particular it is
nonnegative.  The natural free-boundary linearization has Robin form

\[
 A_\Phi\nabla f\cdot\eta=b_\Phi f,\qquad b_\Phi\ge0, \tag{2.14}
\]

where `b_Phi` is the anisotropic support-contact curvature.  The inequality
is the differential form of loss at first contact for a convex support.

Since a CMC leaf has no tangential derivative of its mean curvature,
(2.10) gives

\[
 \lambda'=-J_\Phi f.                                  \tag{2.15}
\]

Dividing (2.15) by `f`, integrating by parts, and using (2.14) proves (0.1).
This argument also shows exactly which statement is needed in a nonsmooth
extension: the self-adjoint Jacobi measure, including its nonnegative
boundary part.  Merely writing an Euler equation on the regular set is not
enough to discard a singular contact measure.

## 3. A dimension-free quantitative turning estimate

We first prove the Euclidean estimate; the anisotropic statement follows
from (2.12)--(2.15) with a constant `C_Phi`.

Let

\[
 Q(v)={1\over P(v)}\int_{\Sigma_v}N\otimes N\,d\sigma,
 \qquad W(v)=1-\operatorname {tr}(Q(v)^2).             \tag{3.1}
\]

Differentiating the numerator, using (2.3), the local area variation
`D_vd sigma=lambda f d sigma`, and (2.7), gives the exact identity

\[
 Q'=-{1\over P}\int_{\Sigma_v}
       \operatorname {sym}(N\otimes\nabla f)d\sigma
 +{\lambda\over P}\int_{\Sigma_v}f(N\otimes N-Q)d\sigma. \tag{3.2}
\]

The second term is important: even if no normal turns along a trajectory,
different fixed planar components can be reweighted.  That is the local
signature of the product/multiphase branch.

Assume on `[a,b]`

\[
 {\ell\over P(v)}\le f\le {L\over P(v)},
 \qquad 0<\ell\le1\le L.                              \tag{3.3}
\]

By Cauchy--Schwarz and `int f d sigma=1`, (3.2) implies

\[
 \|Q'\|_{HS}^2
 \le {8L^2\over P^3}
         \int_{\Sigma_v}|\nabla\log f|^2d\sigma
       +{2L\lambda^2\over P^2}W(v).                  \tag{3.4}
\]

Indeed the squared Hilbert--Schmidt norm of the first integral is at most
`4P int|nabla f|^2`, while

\[
 \left\|\int f(N\otimes N-Q)d\sigma\right\|_{HS}^2
 \le \int f\,d\sigma
      \int f\|N\otimes N-Q\|_{HS}^2d\sigma
 \le LW(v).                                           \tag{3.5}
\]

The lower lapse bound in (3.3) and (2.6) give

\[
 \int|\nabla\log f|^2d\sigma
 \le-\lambda'(v)\int{1\over f}d\sigma
 \le-{P(v)^2\over\ell}\lambda'(v).                  \tag{3.6}
\]

Combining (3.4)--(3.6) proves:

**Proposition 3.1 (regular-lapse inverse).**  If `P(v)>=P_min>0` on
`[a,b]`, then

\[
 \boxed{
 \int_a^b\|Q'(v)\|_{HS}^2dv
 \le {8L^2\over\ell}
       {\lambda(a)-\lambda(b)\over P_{\min}}
 +2L\sup_{[a,b]}{\lambda^2\over P^2}
       \int_a^bW(v)dv.}                               \tag{3.7}
\]

For a fixed uniformly elliptic `Phi`, the same assertion holds with the
right side multiplied by `C_Phi`.  The proof uses (2.12) for the first term,
the anisotropic first-variation transport formula for the second, and
(0.1) in place of (2.6).  Every coefficient is a spherical derivative of
the one fixed `Phi`, so `C_Phi` is independent of `n`.

For any probability law `nu` on `[a,b]`, no regularity of its density is
needed:

\[
 \operatorname {Var}_\nu(Q)
 \le \operatorname {diam}\{Q(v):a\le v\le b\}^2
 \le(b-a)\int_a^b\|Q'(v)\|_{HS}^2dv.                 \tag{3.8}
\]

This observation avoids a hidden Poincare inequality for the value law of
`G`.  In particular it applies to the coarea weights inherited from the
heat function, however nonuniform they are.

Fix `nu_0=.003`.  If throughout an interval of length at most one

\[
 W(v)\le w,quad
 {\lambda(a)-\lambda(b)\over P_{\min}}\le\epsilon,
 \quad {|\lambda|\over P}\le\Lambda,                 \tag{3.9}
\]

then (3.7)--(3.8) give

\[
 \operatorname {Var}_\nu(Q)
 \le C_\Phi\left({L^2\over\ell}\epsilon
                         +L\Lambda^2w\right).          \tag{3.10}
\]

Taking `epsilon` and `w` below
`nu_0/[10 C_Phi max(L^2/ell,L Lambda^2)]` rules out a total coarea
projector variance `nu_0`.  Therefore, under a regular lapse and fixed
topology, the `.003` packet forces a definite curvature/contact slope drop
or a definite within-leaf variance.  If (3.3) fails, the exceptional set is
exactly a lapse/criticality alternative; it cannot be silently absorbed
into (3.10).

## 4. The coherent-leaf killed-tube branch

The regular-lapse estimate controls cross-leaf rotation.  A separate and
stronger observation handles a coherent individual leaf.

Suppose an anisotropic CMC leaf has perimeter `P asymp p`, its outward
Wulff tube survives to distance `T asymp1/p` after deleting at most
`eta P` of base flux, and its normalized projector variance is `W`.  From
(1.6), deleting `eta P` changes the operator norm of its normal matrix by at
most `eta P`.  For `W+eta` smaller than a fixed `Phi`-dependent constant,
the odd Cahn--Hoffman map `N mapsto z(N)` sends the coherent projective
packet to another coherent projective packet.  The killed-tube covariance
bound therefore gives

\[
 I=\operatorname {Cov}(\mu)
 \succeq c_\Phi T^3
     \int_{G}z(N)\otimes z(N)d\sigma_\Phi,
 \qquad
 \left\|\int_Gz\otimes z,d\sigma_\Phi\right\|_{op}
       \ge c_\Phi P.                                  \tag{4.1}
\]

Taking operator norms in (4.1) yields

\[
                         1\ge {c_\Phi\over p^2}.       \tag{4.2}
\]

Thus `p>=c_Phi`.  In the contradiction regime `p<c_Phi`, a coherent CMC
leaf cannot simultaneously have a near-linear long tube and small
contact/focal/collision loss.

This is the precise within-leaf inverse.  It is **coherence**, meaning
`W` small, which gives the immediate covariance contradiction.  Large
`W` is a high-rank escape and is not ruled out by (4.1): the operator norm
can be as small as the reciprocal of the rank.  The product of one-sided
exponentials realizes this obstruction.  Any argument claiming that an
arbitrary fixed positive lower bound on `W` alone contradicts isotropy has
reversed the operator-norm inequality.

For a CMC leaf which is only a near profile minimizer, the scalar killed-tube
proof acquires its relative perimeter deficit and the error in the slope
`lambda` as additional terms.  Formula (0.1) controls the latter over a
regular foliation, but (0.5) alone gives neither pointwise control.  This is
why the joint CMC foliation is useful but not yet a complete transfer.

## 5. What a small integrated deficit does and does not control

There is a clean one-dimensional statement once perimeter concavity has
been established.  Let `P` be concave on `[A,B]`, let `ell_0` be its endpoint
chord, and put `E=P-ell_0>=0`.  If `nu=-P''` is the curvature measure, the
Dirichlet Green kernel gives

\[
 \int_A^BE(v)dv
 ={1\over2}\int_A^B(v-A)(B-v)\,d\nu(v).              \tag{5.1}
\]

Consequently, for `0<h<(B-A)/2`,

\[
 \boxed{
 \nu([A+h,B-h])
 \le {2\over h^2}\int_A^B(P-\ell_0)dv.}              \tag{5.2}
\]

This is the desired slope-drop bound for (3.7), with no dimension factor.
If the isoperimetric profile is already known to be within `delta p` of the
same affine chord and

\[
 \int_A^B(P-I)dv\le\delta p,                          \tag{5.3}
\]

then (5.2) bounds the interior slope drop by `C delta p/h^2`.

The current coarea deficit is not (5.3).  If `r=r(v)` is the value
parameter of `G`, it is

\[
 \int(P-I)\,dr
 =\int(P(v)-I(v))r'(v)dv.                             \tag{5.4}
\]

No upper or lower bound for `r'(v)` follows from CMC.  Replacing `G` by
`h(G)` for a strictly increasing `h` leaves every leaf, its CMC equation,
and all contact geometry unchanged, but replaces `dr` by `h'(r)dr` in both
the coarea matrix and the integrated deficit.  Therefore (5.4) cannot imply
(5.3) without quantitative information from the prescribed heat value law
or the lapse.  This is a genuine invariance obstruction, not a regularity
technicality.

An exact profile foliation makes the point even more sharply.  Then
`P(v)=I(v)` and (5.4) is zero for every reparametrization, while the slope
drop in (5.2) is the curvature of `I` and may be nonzero.  Thus the right
quantity that charges smooth turning is profile curvature, not profile
deficit alone.  Near the putative Cheeger obstruction one expects a nearly
linear profile, but that separate input must be stated and quantitatively
linked to the central CMC leaves.

## 6. Exact rigidity and the possible escape geometries

Suppose on an open interval the right side of (0.1) vanishes and the
foliation has fixed topology.  Uniform ellipticity gives

\[
 \nabla_\Sigma f=0,qquad
 D^2\Phi(N)S=0,qquad
 \nabla^2V[z(N),z(N)]=0,qquad b_\Phi=0.              \tag{6.1}
\]

Since `D^2Phi(N)` is invertible on the tangent space, `S=0`.  Convexity of
`V` implies `nabla^2V z(N)=0`.  Equation (2.12) then shows that `N` is
constant along each material trajectory.  Every connected leaf component
is a planar piece and the components move by parallel Wulff translation.

This proves the following exact local classification.

**Proposition 6.1 (zero-charge classification).**  On every connected
spacetime component of a regular zero-charge CMC foliation, all leaves are
parallel planar sections with one fixed projective normal.  The potential
is affine in the corresponding Wulff direction throughout the swept slab,
and the support-contact curvature in that direction vanishes.  A change of
projective normal can occur only by:

1. creation or deletion of a component at a critical/support-contact value;
2. focal or cut-locus collision of components; or
3. coexistence and reweighting of several fixed planar phase families.

The third alternative is the local precursor of an affine product of
one-sided exponentials.  If the swept support and potential globally split
in the fixed phase directions, whitening makes those independent
exponential directions orthogonal.  The local conclusion (6.1), by itself,
does not prove that global splitting: support components may terminate on
different faces.  Calling all such configurations an orthogonal product
without proving the global support decomposition would be an unjustified
step.

Concurrent spherical leaves are not zero-charge.  Their second fundamental
form contributes to `q_Phi`; moving centers also move the contact set or
produce a cut event.  They belong to the concurrent/contact branch.  A
long-radius sphere can make the pointwise curvature small, but the
finite-distance killed-tube identity accumulates it over the corresponding
long scale.

## 7. Model audit

### 7.1 Parallel planes

For `Sigma_t={x dot u=t}`, `S=0` and `W=0`.  Weighted CMC says that
`partial_uV` is constant on each plane.  On a full slab this integrates to

\[
                         V(y+tu)=V_0(y)+v(t).           \tag{7.1}
\]

If the profile slope is exactly constant, (0.1) gives `v''=0` and the lapse
is constant on each plane.  This is the affine/exponential factor.  There is
no cross-level projector motion.

### 7.2 Moving-center spheres and ball caps

For a Euclidean sphere of radius `R`,

\[
 |S|^2={n-1\over R^2},                                \tag{7.2}
\]

so concurrence at its center is recorded in `q`.  For the uniform measure
on an isotropic Euclidean ball, the exact relative isoperimetric regions are
orthogonal spherical caps.  In a two-dimensional section, they are the
Euclidean circles orthogonal to the support circle.  Choose a continuously
nested family by choosing nested boundary arcs whose midpoints move while
their half-lengths increase.  The corresponding orthogonal circles are
disjoint, their cap regions are nested, and every leaf is an exact profile
minimizer with constant mean curvature and natural free-boundary contact.

Thus

\[
 P(v)-I(v)=0                                           \tag{7.3}
\]

on the whole family, although its centers and representative directions
move.  The motion is paid by (7.2) and the support-contact term.  By applying
a monotone reparametrization to the common cap function, one may put most
coarea weight near two caps whose projector matrices are separated; a lower
bound `.003` is immediate.  This example verifies all the static assumptions
of the revised setup and proves that the contact/concurrent alternative is
indispensable.  The isotropic ball itself has a universal Cheeger constant,
so it is not a KLS counterexample.

### 7.3 Cylinders

For a round cylinder with normal sphere of dimension `k-1`, symmetry gives

\[
 Q={1\over k}P_U,qquad W=1-{1\over k}.                \tag{7.4}
\]

Hence `k>=2` is a genuine within-leaf high-rank branch, not a coherent
cross-level branch.  Its curvature charge contains `(k-1)/R^2`.  If the
radius is so large that this is small, the killed tube carries a rank-`k`
packet and the covariance estimate only forces `k` to be correspondingly
large; it does not by itself give a contradiction.  A product or incidence
argument is still required in this high-rank case.

### 7.4 Product exponentials

For a product of one-sided exponential measures, the canonical phase
boundaries are flat coordinate sheets.  Each sheet has zero intrinsic
turning, while several sheets on one level give a high-rank `Q`.  Switching
which sheet is present occurs at a ridge, component birth, or collision.
After isotropic whitening, the independent factor directions are
orthogonal.  This exactly realizes item 3 of Proposition 6.1 and shows why
the reweighting term in (3.2) cannot be discarded.

### 7.5 Boxes

Coordinate slab cuts of a box are parallel and zero-charge.  A continuously
rotating family of planar cuts has constant interior mean curvature but does
not satisfy the natural contact condition except at isolated orientations.
Its endpoints slide across edges and the missing term is a contact/critical
measure.  In a smooth rounding of the box, the same cost appears as the
nonnegative boundary term in (2.6).  Thus a wavy sequence of planar chords
is not a counterexample to (0.1).

### 7.6 Wavy foliations

For a smooth wavy foliation with no contact, (2.6) gives

\[
 -\lambda'\int f^{-1}d\sigma
 \ge\int|\nabla\log f|^2d\sigma+\int|S|^2d\sigma.     \tag{7.5}
\]

Normal turning along physical normal distance is
`-nabla log f`.  Therefore a wavy foliation cannot have both a nearly
constant CMC slope and a regular lapse.  The only smooth escape is lapse
concentration; geometrically this is approach to a critical/focal/contact
set.  This is precisely the exceptional case excluded in Proposition 3.1,
not an untracked error.

## 8. What remains for the joint replacement route

The revised constant-anisotropy construction has removed the spatial
selector derivative.  If its `BV` minimizer can be approximated by a smooth
regular foliation so that (0.1) holds with all singular terms retained, the
remaining proof obligation can be stated narrowly:

1. derive, from joint equimeasurable minimality and the heat value law, a
   central subband on which the CMC perimeter is nearly affine and the
   Wulff lapse obeys (0.3), or charge the failure by a fixed amount of
   critical/contact/collision flux;
2. apply (3.7) to rule out cross-level rotation of coherent leaves;
3. apply (4.1) to any coherent long-surviving leaf;
4. prove that a high-rank within-leaf packet with small Jacobi and killed
   flux is a genuine global orthogonal-product factor, rather than merely
   a collection of locally flat components ending on unrelated support
   faces.

Items 1 and 4 are not consequences of CMC plus the integrated deficit as
currently stated.  Item 1 is a quantitative lapse/singular-measure theorem;
item 4 is the high-rank tube/product inverse.  Formula (0.1) and Proposition
3.1 provide the dimension-free analytic core once those two geometric
inputs are supplied.

## 9. Generality audit

1. **Dimension.**  All displayed constants depend only on the fixed
   ellipticity and spherical regularity of `Phi`; no trace estimate uses
   `n`.
2. **Hard support.**  The contact term is retained in (2.6), (0.1), and all
   model tests.  Dropping it makes the ball-cap and rotating-box examples
   false positives.
3. **Critical and singular levels.**  The proved identities apply only on
   a smooth fixed-topology band.  At excluded values the correct object is
   a nonnegative killed/contact/collision measure.  Its quantitative
   control for the `BV` minimizer is not asserted.
4. **Parameterization.**  Proposition 3.1 uses volume and (3.8) handles an
   arbitrary coarea law.  The deficit conversion failure (5.4) is stated
   explicitly.
5. **Circularity.**  No Poincare, Cheeger, ball-walk, or desired KLS bound
   is assumed.  The only covariance use is the elementary long-tube matrix
   inequality in (4.1).
6. **High rank.**  The note does not claim that covariance alone excludes a
   high-rank normal packet.  Cylinders and product exponentials explicitly
   audit that failure.

