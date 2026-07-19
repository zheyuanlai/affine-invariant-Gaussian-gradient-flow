# Facet completion from exact Cheeger minimality

## 0. Verdict

There is a complete, dimension-free closure of the **exact polyhedral
log-affine branch**.  In fact, neither the high-rank normal matrix nor the
killed-ray packet is needed in that branch.

Let \(\Omega\subset\mathbb R^n\) be a convex polyhedron, let

\[
 d\mu=Z^{-1}e^{-a\cdot x}{\bf1}_{\Omega}(x)\,dx,
\]

assume that \(\mu\) is isotropic, and let \(\Phi\) be an even, strictly
convex, one-homogeneous anisotropy
with \(\Phi(u)\ge\phi_-|u|\).  Suppose that a polyhedral set \(E\) of
mass \(1/2\) minimizes the relative anisotropic perimeter among all sets
of mass \(1/2\).  Then every component of its interior boundary is a
complete hyperplane slice of \(\Omega\).  If \(p=P_{\Phi,\mu}(E)\), then

\[
                 p\ge {\phi_-\over4\sqrt3},
 \qquad
                 {p\over\min(\mu(E),1-\mu(E))}
                    \ge {\phi_-\over2\sqrt3}.                 \tag{0.1}
\]

The mechanism is explicit.  At a non-flat interior ridge, replacing two
short facet segments by their chord saves

\[
 \varepsilon\bigl[
       \Phi(J\tau _1)+\Phi(J\tau _2)
       -\Phi(J(\tau _2-\tau _1))\bigr]
       \int_R Z^{-1}e^{-a\cdot x}\,d\mathcal H^{n-2}(x)
       +O(\varepsilon^2),                                    \tag{0.2}
\]

where the bracket is strictly positive.  The change of volume is only
\(O(\varepsilon^2)\), and can be repaired on a regular facet at perimeter
cost \(O(\varepsilon^2)\).  Thus exact minimality forbids every interior
ridge.  A polyhedral boundary facet can then end only on
\(\partial\Omega\), so it is a complete slice.  A central-cell incidence
argument shows that the smaller halfspaces cut off by all these slices
have total mass at least \(1/2\); the one-dimensional isotropic
log-concave bound then proves (0.1).

This is not yet a closure of the smooth high-rank branch.  The obstruction
to passing through arbitrary polyhedral approximation is quantitative and
real: for a mesh approximation of a smooth interface, dihedral angles and
facet clearances vanish together, and the total bevel gain is
\(O(\ell^2)\).  Section 8 gives a precise no-go calculation.  A successful
passage to the general branch must produce **macroscopic** facets/ridges
from the killed-ray packet, or replace local beveling by a global
union/intersection or calibration argument.

## 1. Setting

Assume throughout Sections 1--7 that \(n\ge2\), \(\Omega\) is a
full-dimensional open convex polyhedron, and

\[
                 \rho(x)=Z^{-1}e^{-a\cdot x}{\bf1}_\Omega(x) \tag{1.1}
\]

is an isotropic probability density.  For a finite-perimeter set \(F\),
write

\[
 P_{\Phi,\mu}(F)
   =\int_{\partial^*F\cap\Omega}\Phi(\nu_F)\rho\,
                                      d\mathcal H^{n-1}.       \tag{1.2}
\]

The support boundary is free and therefore is not counted.  The
anisotropy \(\Phi:\mathbb R^n\to[0,\infty)\) is assumed even,
one-homogeneous, strictly convex, and continuous, with

\[
                  \phi_-|q|\le\Phi(q)\le\phi_+|q|.            \tag{1.3}
\]

The small matrix perturbations used in the Wulff replacement satisfy
these assumptions; in the Euclidean case \(\Phi(q)=|q|\) and
\(\phi_-=1\).

A *polyhedral set relative to* \(\Omega\) means that, modulo null sets,
its interior reduced boundary is a finite union of relatively open convex
\((n-1)\)-polytopes.  We assume

\[
 \mu(E)=\frac12,
 \qquad
 P_{\Phi,\mu}(E)\le P_{\Phi,\mu}(F)
        \quad\text{whenever }\mu(F)=\frac12.                  \tag{1.4}
\]

If \(E\) is a balanced global Cheeger minimizer, (1.4) is automatic.

## 2. The local ridge formula

### Lemma 2.1 (strict bevel gain)

Let \(R\Subset\Omega\) be a relatively open piece of an
\((n-2)\)-dimensional ridge at which two non-coplanar facets of
\(\partial E\) meet and no other facet is incident.  Then there are
equal-volume competitors \(E_\varepsilon\), supported in an arbitrarily
small neighborhood of \(R\) together with one regular repair patch, such
that

\[
 P_{\Phi,\mu}(E)-P_{\Phi,\mu}(E_\varepsilon)
  =\varepsilon\beta_\Phi(\tau_1,\tau_2)
       \int_R\rho\,d\mathcal H^{n-2}+O_R(\varepsilon^2),      \tag{2.1}
\]

where

\[
 \beta_\Phi(\tau_1,\tau_2)
  :=\Phi(J\tau_1)+\Phi(J\tau_2)
       -\Phi(J(\tau_2-\tau_1))>0.                            \tag{2.2}
\]

Here \(N\) is the two-plane perpendicular to \(R\), \(J\) is a
quarter-turn in \(N\), and \(\tau_1,\tau_2\in N\) are the unit rays of
the two facets pointing away from the ridge.  The volume before repair
satisfies

\[
 |\mu(\widetilde E_\varepsilon)-\mu(E)|
 \le {\varepsilon^2\over2}
       |\det_N(\tau_1,\tau_2)|
       \int_R\rho\,d\mathcal H^{n-2}+O_R(\varepsilon^3).      \tag{2.3}
\]

#### Proof

First take a compact polyhedral subpatch \(R_0\Subset R\).  In the normal
two-plane, the old interface joins \(\varepsilon\tau_1\) to the origin
and the origin to \(\varepsilon\tau_2\).  Replace it by the chord joining
\(\varepsilon\tau_1\) to \(\varepsilon\tau_2\), choosing the side of the
chord consistent with the old phase.  Define

\[
                         g(v)=\Phi(Jv),\qquad v\in N.          \tag{2.4}
\]

This is an even strictly convex norm on \(N\).  Per unit ridge area the
old anisotropic length is

\[
                   \varepsilon[g(\tau_1)+g(\tau_2)],          \tag{2.5}
\]

whereas the chord has anisotropic length

\[
                         \varepsilon g(\tau_2-\tau_1).        \tag{2.6}
\]

The strict triangle inequality applied to
\(-\tau_1+\tau_2\) shows that (2.2) is positive unless
\(\tau_2=-\tau_1\).  The exceptional case is precisely a flat
hyperplane, not a genuine ridge.

Since

\[
 {\rho(x+y)\over\rho(x)}=e^{-a\cdot y}=1+O_R(\varepsilon)    \tag{2.7}
\]

on the bevel, (2.5)--(2.7), integrated over \(R_0\), give (2.1) before
volume repair.  Tapering the construction near \(\partial R_0\) adds
only \(O_R(\varepsilon^2)\): the taper has two-dimensional cross-sectional
area \(O(\varepsilon^2)\) over an \((n-3)\)-dimensional set.  For \(n=2\)
there is no taper.  The symmetric difference in each normal two-plane is
the triangle spanned by \(\varepsilon\tau_1,\varepsilon\tau_2\), proving
(2.3).

It remains to repair the mass.  Choose a compact patch \(U\) in the
relative interior of any boundary facet, disjoint from the bevel, and a
nonnegative piecewise smooth function \(\eta\) compactly supported in
\(U\), not identically zero.  Replacing this facet locally by its normal
graph of height \(s\eta\) changes volume by

\[
                    s\int_U\eta\rho\,d\mathcal H^{n-1}
                       +O_U(s^2),                              \tag{2.8}
\]

and changes perimeter by \(O_U(|s|)\).  Both signs of \(s\) are
admissible because \(U\Subset\Omega\).  Equations (2.3) and (2.8) give a
choice \(s=O_R(\varepsilon^2)\) restoring the mass exactly, at perimeter
cost \(O_R(\varepsilon^2)\).  This proves (2.1).  QED.

The same proof handles a ridge with more than two incident facets: in the
normal two-plane choose two consecutive boundary rays.  A single boundary
ray cannot terminate in the interior, because the indicator of a set has
an even number of jumps around a small circle.

### Corollary 2.2 (no interior ridges)

Under (1.4), every interior ridge of \(\partial E\) is flat.  After
merging coplanar facets, \(\partial E\cap\Omega\) has no interior
\((n-2)\)-faces.

Indeed, the positive \(O(\varepsilon)\) term in (2.1) contradicts (1.4)
for sufficiently small \(\varepsilon\).

### 2.1 A quantitative near-minimal version

If instead \(E\) has additive fixed-volume deficit at most \(\delta\),
the same calculation gives, for every isolated ridge patch and every
\(0<\varepsilon<\varepsilon_R\),

\[
 \delta\ge
   \varepsilon\beta_\Phi(\tau_1,\tau_2)
        \int_R\rho\,d\mathcal H^{n-2}-C_R\varepsilon^2.     \tag{2.9}
\]

This formula is useful only when the angle, ridge mass, and allowable
clearance \(\varepsilon_R\) have dimension-free lower control.  None of
those controls follows from a formal mesh approximation; Section 8 makes
the failure explicit.

## 3. Contact ridges and the anisotropic Young law

Although contact with \(\partial\Omega\) is compatible with slice
completion, it is useful to record the corresponding explicit competitor.
Let a flat interior facet with tangent ray \(\tau\) meet a flat support
face with exterior normal \(\nu_\Omega\).  In the normal two-plane let
\(s=J\nu_\Omega\) be the support tangent and put
\(z=D\Phi(J\tau)\).  Keep the interface point
\(\varepsilon\tau\) fixed and move its contact endpoint from the origin
to \(t\varepsilon s\).  The old and new relative perimeter costs per unit
contact-ridge area are

\[
                     \varepsilon g(\tau),
 \qquad
                     \varepsilon g(\tau-t s),                 \tag{3.1}
\]

and

\[
 {d\over dt}\bigg|_{t=0}g(\tau-t s)
                         =z\cdot\nu_\Omega.                   \tag{3.2}
\]

The swept volume is \(O(t\varepsilon^2)\), so it can again be repaired
at second-order cost.  Consequently an exact minimizer at a two-sided
flat support face must satisfy

\[
                         D\Phi(\nu_E)\cdot\nu_\Omega=0.        \tag{3.3}
\]

At a corner of \(\Omega\), the same argument gives the corresponding
normal-cone variational inequality.  Formula (3.3) is the anisotropic
orthogonality (Young) law.  A contact satisfying (3.3) is not a defect:
it is exactly how a complete hyperplane slice is allowed to end on the
free boundary.

## 4. No ridges implies full slices

### Lemma 4.1 (polyhedral completion)

After merging coplanar facets, every nonempty component of
\(\partial E\cap\Omega\) equals, modulo
\(\mathcal H^{n-1}\)-null sets,

\[
                              \Pi\cap\Omega                    \tag{4.1}
\]

for some affine hyperplane \(\Pi\).

#### Proof

Fix a facet plane \(\Pi\), and let \(S_\Pi\) be the union of the relative
interiors of all facets contained in \(\Pi\).  It is relatively open in
the connected convex set \(D=\Pi\cap\Omega\).  If its relative boundary
met the relative interior of \(D\), polyhedrality would give an
\((n-2)\)-dimensional face at a generic such point.  A boundary wall of a
set cannot terminate there by itself (intersect with a small normal circle
and count jumps), so another facet must be incident.  A coplanar incident
facet was already merged; a non-coplanar one is forbidden by Corollary
2.2.  Therefore \(S_\Pi\) has no relative boundary in the relative
interior of \(D\).  Its closure is both open and closed in the relative
interior of \(D\), and connectedness gives (4.1).  QED.

Write the resulting distinct slice hyperplanes as
\(\Pi_1,\ldots,\Pi_N\).  Their slices cannot cross in \(\Omega\), since a
crossing would be an interior non-flat ridge.  Intersections on
\(\partial\Omega\) or on sets of codimension at least two do not affect
perimeter.

## 5. The central-cell incidence bound

For each \(i\), choose the halfspace \(H_i^-\) bounded by \(\Pi_i\) with
smaller mass, and put

\[
 q_i=\mu(H_i^-\cap\Omega)\le\frac12,
 \qquad H_i^+=\mathbb R^n\setminus\operatorname{int}H_i^- .  \tag{5.1}
\]

### Lemma 5.1 (the smaller sides have total mass at least one half)

\[
                              \sum_{i=1}^Nq_i\ge\frac12.       \tag{5.2}
\]

#### Proof

Consider the central cell

\[
                         C=\Omega\cap\bigcap_{i=1}^NH_i^+.    \tag{5.3}
\]

The union bound gives

\[
                              \mu(C)\ge1-\sum_iq_i.            \tag{5.4}
\]

Its relative interior is convex and hence connected, and it meets none of
the slice hyperplanes.  Since the entire interior boundary of \(E\) is the
union of those slices, the polyhedral indicator \({\bf1}_E\) is constant
almost everywhere on \(C\).  If \(\sum_iq_i<1/2\), then (5.4) gives
\(\mu(C)>1/2\).  If the constant is one, then \(\mu(E)>1/2\); if it is
zero, then \(\mu(E^c)>1/2\).  Both contradict \(\mu(E)=1/2\).  This proves
(5.2).  QED.

This incidence lemma is the bounded-reuse statement that is missing from
a merely local killed-ray calculation: every complete slice is counted
once, while the central cell prevents all of its smaller sides from having
negligible total mass.

## 6. Dimension-free closure of the polyhedral branch

We use the following elementary one-dimensional fact.

### Lemma 6.1 (isotropic log-concave halfline bound)

If \(\nu\) is an isotropic log-concave probability on \(\mathbb R\), and
\(H\) is a halfline of mass \(q\le1/2\), then

\[
                              \nu^+(H)\ge {q\over2\sqrt3}.     \tag{6.1}
\]

#### Proof

Let \(f\) be the density and
\(I(t)=f(F^{-1}(t))\).  Log-concavity makes \(I\) concave on \([0,1]\),
with zero endpoint limits.  If \(M=\|f\|_\infty\), concavity gives
\(I(1/2)\ge M/2\).  A probability density bounded by \(M\) has variance
at least \(1/(12M^2)\), by symmetric decreasing rearrangement and the
bathtub principle.  Since the variance is one,
\(M\ge1/(2\sqrt3)\).  Finally, concavity and \(I(0)=0\) give

\[
 {I(q)\over q}\ge {I(1/2)\over1/2}\ge M\ge{1\over2\sqrt3},  \tag{6.2}
\]

which is (6.1).  QED.

### Theorem 6.2 (polyhedral facet completion inverse)

Under (1.1)--(1.4),

\[
                         P_{\Phi,\mu}(E)
                           \ge {\phi_-\over4\sqrt3}.           \tag{6.3}
\]

#### Proof

Let \(u_i\) be a unit normal to \(\Pi_i\).  Isotropy of \(\mu\) implies
that the marginal \(u_i\cdot X\) is an isotropic one-dimensional
log-concave probability.  Its density at the slice level is the Euclidean
weighted area of \(\Pi_i\cap\Omega\).  Hence Lemma 6.1 and (1.3) give

\[
 \int_{\Pi_i\cap\Omega}\Phi(u_i)\rho\,d\mathcal H^{n-1}
                         \ge {\phi_-q_i\over2\sqrt3}.          \tag{6.4}
\]

The full slices are precisely the components of the interior boundary;
they overlap only in null sets.  Summing (6.4) and using Lemma 5.1 yields
(6.3).  QED.

For Euclidean perimeter, a balanced Cheeger minimizer therefore has
\(\psi_\mu=2P_\mu(E)\ge1/(2\sqrt3)\).  For the small constant anisotropy
arising after whitening or matrix replacement, the same conclusion loses
only the lower ellipticity factor \(\phi_-\).

For a non-isotropic full-dimensional log-affine law with covariance
\(A\succ0\), apply the theorem after the whitening map
\(y=A^{-1/2}(x-\mathbb EX)\).  Euclidean weighted perimeter in the
original variables becomes the constant anisotropy

\[
                         \Phi_A(q)=|A^{-1/2}q|,                \tag{6.5}
\]

whose lower ellipticity is \(\|A\|_{\rm op}^{-1/2}\).  Thus the same
polyhedral branch gives

\[
                    \psi_\mu\ge {1\over
                        2\sqrt3\sqrt{\|A\|_{\rm op}}}.        \tag{6.6}
\]

If the support has affine dimension \(k\ge1\), the argument is applied
with its intrinsic Euclidean structure on that supporting subspace; its
intrinsic covariance is positive definite, and (6.6) uses the same
nonzero operator norm as the ambient covariance.  The case \(k=0\) is the
excluded point mass.

Notice that the hypotheses about high normalized normal variance,
bounded multiplier, and a fixed killed-ray packet were not used.  Exact
polyhedrality plus exact global minimality is already more rigid.

## 7. The two required stress tests

### 7.1 Product exponentials

Let \(Y_1,\ldots,Y_m\) be independent unit exponentials, let
\(d=1-e^{-L}\), and consider

\[
                        E_L=\{\max_iY_i\ge L\}.                \tag{7.1}
\]

At the median level \(d^m=1/2\), its perimeter is

\[
 p_m={m\over2}(2^{1/m}-1)
             \in\left[{\log2\over2},{1\over2}\right].       \tag{7.2}
\]

Thus this high-rank example already lies at universal Cheeger scale.  It
also illustrates exactly where the bevel branch fires.  For every pair
\(i<j\), the two facets meet on

\[
 R_{ij}=\{Y_i=Y_j=L,\ Y_k<L\ (k\ne i,j)\},                   \tag{7.3}
\]

whose weighted codimension-two area is

\[
                              q^2d^{m-2},\qquad q=e^{-L}.      \tag{7.4}
\]

The normal two-dimensional cross-section is a reentrant right corner.
Removing two segments of length \(\varepsilon\) and inserting the chord
saves

\[
             (2-\sqrt2)\varepsilon q^2d^{m-2}+O(\varepsilon^2)
                                                                    \tag{7.5}
\]

at that ridge, while changing mass by
\(\frac12\varepsilon^2q^2d^{m-2}+O(\varepsilon^3)\).  Summed over pairs,
the leading ridge mass is

\[
 R_m={m\choose2}(1-d)^2d^{m-2}
             \longrightarrow{(\log2)^2\over4}.              \tag{7.6}
\]

For a \((1\pm\kappa)\)-Euclidean anisotropy, the right-corner bracket is
at least

\[
                    2(1-\kappa)-\sqrt2(1+\kappa)>0            \tag{7.7}
\]

whenever \(\kappa<3-2\sqrt2\), in particular for the tiny perturbations
used in the matrix argument.

Therefore the maximum set is not an exact polyhedral fixed-volume
minimizer: its genuine internal ridges admit the explicit saving (7.5).
This is not a false rejection of the measure.  Bobkov--Houdre
tensorization gives a universal lower Cheeger bound for the product, and
(7.2) gives a universal upper bound.

### 7.2 Cyclically constrained exponentials

For

\[
 \Omega_{m,a}=\{y_i\ge0,\ y_i+y_{i+1}\ge a\},
 \qquad a=m^{-6},                                             \tag{7.8}
\]

with density proportional to \(e^{-\sum_i y_i}\), the median maximum
leaf has perimeter between `.33` and `.44`, as computed in
`balanced_tube_inverse.md`.  Since its median level satisfies \(L>a\),
every tie ridge \(Y_i=Y_j=L\) contains a positive-measure subpatch a
positive distance from the active support faces.  On that subpatch the
same right-angle calculation (7.5) applies.  A union bound on the remaining
path constraints changes the aggregate ridge weight from (7.6) by
\(O(m^{-11})\).  Thus this maximum leaf is likewise not an exact
polyhedral minimizer, and its perimeter is already universal.

The underlying conditioned measure also has a universal Cheeger constant,
so the stress test itself remains at scale one.  Indeed, the conditioning
event is convex and has product-exponential probability
\(1-O(m^{-11})\).  Dimension-free exponential concentration of the
product (from Bobkov--Houdre tensorization) passes to conditioning with a
universal change of constants: a set of conditional mass at least one
half has original mass at least `.49`, and the original concentration
bound applies after a universal initial enlargement.  The conditioned law
is log-concave, so E. Milman's concentration--Cheeger equivalence yields a
universal conditional Cheeger lower bound.  Whitening changes it only by
the ellipticity factor in (4.16) of `balanced_tube_inverse.md`.

These tests show the intended dichotomy cleanly.  Their incomplete facets
are accompanied by real, nondegenerate internal ridges; beveling rejects
the *candidate sets*, while independent dimension-free estimates confirm
that the *measures* have \(\psi\asymp1\).

## 8. Why arbitrary polyhedral approximation does not transfer the theorem

The exact result above is not stable under mesh refinement without a new
macroscopic input.  The obstruction can be quantified already in the
Euclidean plane.

Let two unit boundary rays make exterior turning angle \(\theta\), so that
the rays become opposite when \(\theta=0\).  The Euclidean triangle defect
in (2.2) is

\[
 \beta(\theta)=2-2\cos(\theta/2)
                   ={\theta^2\over4}+O(\theta^4).              \tag{8.1}
\]

Approximate a smooth unit circle by a regular \(N\)-gon.  Its mesh size
and exterior angle are both \(\ell\asymp N^{-1}\).  A bevel cannot be
taken deeper than \(C\ell\) without meeting the neighboring vertices.
Even beveling all \(N\) vertices therefore gives at most

\[
       N\,(C\ell)\,\beta(2\pi/N)=O(N^{-2})=O(\ell^2)          \tag{8.2}
\]

of perimeter gain.  The polygon-to-circle perimeter error is itself
\(O(\ell^2)\).  Meanwhile, if the polygon lies in a fixed larger convex
support, none of its edge patches is a complete support slice, and the sum
of their formal full-slice completion defects need not tend to zero (it
can grow with \(N\)).

The same scaling holds in higher dimensions.  For a shape-regular mesh of
a smooth hypersurface, total ridge measure is \(O(A/\ell)\), dihedral
turning is \(O(\ell|II|)\), and available bevel depth is \(O(\ell)\).
Consequently the summed local strict-triangle gain is

\[
 {A\over\ell}\cdot\ell\cdot O(\ell^2|II|^2)
                         =O(A\ell^2|II|^2).                   \tag{8.3}
\]

It vanishes at precisely the scale of the approximation error.  High
global normal rank does not repair (8.3): a smooth sphere has maximally
spread normals while every adjacent pair differs by only \(O(\ell)\).

Thus the following broad implication is false without additional
hypotheses:

\[
 \text{small fixed-volume deficit of a polyhedral approximation}
 \quad\Longrightarrow\quad
 \text{small full-slice completion defect}.                  \tag{8.4}
\]

Local bevel gains see squared adjacent turning times the local mesh scale;
full-slice completion is global incidence data.

## 9. Exact approximation requirement for the smooth high-rank branch

To turn Theorem 6.2 into the desired one-interface inverse, one needs one
of the following genuinely new statements, with constants independent of
dimension.

1. **Macroscopic polyhedral extraction.**  From the fixed killed-ray
   packet, extract boundary patches with tangential clearance
   \(r_0h/p\), and show that every noncompletion produces ridges with
   angle-weighted mass satisfying

   \[
    \sum_R r_R\,\beta_\Phi(R)
         \int_R\rho\,d\mathcal H^{n-2}\ge c p.                \tag{9.1}
   \]

   Formula (2.9) would then contradict the near-Cheeger deficit; surviving
   patches would be complete and Theorem 6.2 would apply.

2. **Global bounded-reuse competitor.**  Prove directly that missing
   portions of the full marginal slices can be assigned to
   union/intersection competitors with total physical perimeter charged at
   most \(C\) times.  This would establish
   \(R_{\rm comp}\le Cp\) on a packet sweeping fixed mass, and Proposition
   2.2 of `balanced_tube_inverse.md` would close the branch.

3. **Common-calibration incidence.**  Use the single bounded-divergence
   calibration shared by the direct Cheeger-deficit level sets to replace
   the polyhedral central cell.  One needs a measurable analogue of Lemma
   5.1 saying that persistent Wulff flow cells either carry more than half
   the mass or encounter a calibrated interface with bounded reuse.

Any approximation must also preserve, with the same universal constants,
the covariance normalization, the anisotropic ellipticity, the balanced
mass, the fixed swept mass, and the additive perimeter deficit.  Equations
(8.1)--(8.3) show that merely smoothing the density and triangulating the
interface cannot provide (9.1).

The positive conclusion is therefore exact but sharply delimited: **exact
polyhedral global minimality forces full-slice completion and yields a
dimension-free inverse; local beveling alone has no stable route from a
smooth high-rank minimizer to that rigid category.**
