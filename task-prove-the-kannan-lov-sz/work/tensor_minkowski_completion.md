# Tensor Minkowski completion for a balanced interface

## 0. Purpose and verdict

This note isolates an exact matrix identity that converts the smooth
one-interface problem into a support-contact problem.  It gives a complete
dimension-free inverse for a closed constant-mean-curvature interface once
its boundary-position second moment is controlled at the isotropic scale.
For a free-boundary interface, the only additional term is an explicit
contact tensor.  Thus the identity is a smooth analogue of the polyhedral
central-cell split: either the interior interface itself forces a universal
perimeter, or a dimension-sized tensor must be carried by its contact with
the convex support.

The contact-tensor inverse is not proved here.  The point of the calculation
is to identify exactly what a global flow-cell theorem must control, without
losing a factor of the dimension through a scalar trace estimate.

## 1. Exact tensor identity

Let \(K\subset\mathbb R^d\) be a bounded convex body with \(C^2\) boundary.
Let \(E\subset K\) be a relative isoperimetric region, and suppose for the
moment that the regular interior interface

\[
 \Sigma=\overline{\partial E\cap\operatorname{int}K}
\]

is a compact, oriented, \(C^2\) hypersurface with boundary
\(\Gamma=\partial\Sigma\subset\partial K\).  Write \(n\) for its unit
normal, \(H=\operatorname{div}_\Sigma n\) for its scalar mean curvature,
and \(P=\mathcal H^{d-1}(\Sigma)\).  First variation gives that \(H\) is
constant, and the free-boundary condition is

\[
                         n\cdot n_K=0\quad\hbox{on }\Gamma.       \tag{1.1}
\]

Let \(\nu\) be the outward unit conormal of \(\Gamma\) in \(\Sigma\).  Up to
the orientation of \(\Gamma\), (1.1) gives \(\nu=n_K\).  Fix
\(c\in\mathbb R^d\), and define

\[
 Q=\frac1P\int_\Sigma n\otimes n\,dA,
 \qquad
 X=\frac1P\int_\Sigma (x-c)\otimes n\,dA,
 \qquad
 B=\frac1P\int_\Gamma (x-c)\otimes\nu\,dS.                  \tag{1.2}
\]

Then

\[
                         \boxed{I-Q=H X+B.}                    \tag{1.3}
\]

To prove (1.3), take arbitrary constant vectors \(a,b\in\mathbb R^d\) and
the tangential vector field

\[
 Z=(a\cdot(x-c))\,[b-(b\cdot n)n].                           \tag{1.4}
\]

An orthonormal-frame computation gives

\[
 \operatorname{div}_\Sigma Z
 =a\cdot b-(a\cdot n)(b\cdot n)
   -H(a\cdot(x-c))(b\cdot n).                               \tag{1.5}
\]

The surface divergence theorem and \(n\cdot\nu=0\) give

\[
 \int_\Sigma\operatorname{div}_\Sigma Z\,dA
 =\int_\Gamma(a\cdot(x-c))(b\cdot\nu)\,dS.                 \tag{1.6}
\]

Equations (1.5)--(1.6), for every \(a,b\), are precisely (1.3), with the
sign of \(B\) changed if the opposite conormal convention is chosen.  Every
subsequent estimate uses only \(\|B\|\), so this convention is harmless.

The identity remains valid componentwise.  Singular strata of codimension
at least two require an exhaustion by regular patches and a proof that the
artificial boundary flux tends to zero; that limiting statement is part of
the regularity audit and is not silently assumed here.

## 2. The dimension is retained in Frobenius norm

The normal matrix is positive semidefinite and has trace one.  Hence

\[
 \|I-Q\|_F^2
 =d-2+\operatorname{tr}(Q^2)\ge d-2.                       \tag{2.1}
\]

This is the crucial cancellation missed by scalar normal estimates.  The
left side is of order \(d\), even when \(Q=I/d\).  Thus a high-dimensional
spread of normals strengthens rather than weakens the tensor identity.

For the position tensor, Cauchy--Schwarz gives

\[
 \|X\|_F^2
 \le \frac1P\int_\Sigma|x-c|^2\,dA.                        \tag{2.2}
\]

Combining (1.3)--(2.2) yields the exact dichotomy

\[
 \sqrt{d-2}
 \le |H|\left(\frac1P\int_\Sigma|x-c|^2dA\right)^{1/2}
       +\|B\|_F.                                             \tag{2.3}
\]

No operator-norm or trace relaxation has been used.

## 3. Closed-interface inverse

If \(\Gamma=\varnothing\), then \(B=0\).  Assume

\[
                    \frac1P\int_\Sigma|x-c|^2dA\le C_0d.    \tag{3.1}
\]

For a balanced relative Cheeger minimizer, the one-sided derivatives of the
concave isoperimetric profile give

\[
                             |H|\le\psi=2P                  \tag{3.2}
\]

after normalizing \(|K|=1\).  Equations (2.3), (3.1), and (3.2) imply, for
\(d\ge4\),

\[
                    P\ge\frac1{2\sqrt{2C_0}},
 \qquad
                    \psi\ge\frac1{\sqrt{2C_0}}.             \tag{3.3}
\]

The dimensions \(d\le3\) are absorbed into any fixed finite-dimensional
Cheeger bound.  Thus a closed interior component cannot be the source of a
vanishing KLS constant once (3.1) is available.

The estimate is sharp in its dimensional bookkeeping.  A sphere of radius
\(r\) has \(Q=I/d\), \(|H|=(d-1)/r\), and surface second moment \(r^2\);
both sides of (2.1)--(2.2) are of order \(d\).  There is no hidden
\(d^{-1/2}\) loss.

## 4. Free boundary: the exact remaining tensor

For \(\Gamma\ne\varnothing\), if both

\[
 \frac1P\int_\Sigma|x-c|^2dA\le C_0d,
 \qquad |H|\le2P,                                           \tag{4.1}
\]

hold and \(P<1/(4\sqrt{2C_0})\), then (2.3) forces

\[
                         \boxed{\|B\|_F\ge\tfrac14\sqrt d} \tag{4.2}
\]

for \(d\ge4\), after a harmless weakening of constants.  In unnormalised
form, the support contact carries a tensor of size

\[
 \left\|\int_\Gamma(x-c)\otimes n_K\,dS\right\|_F
                         \ge\tfrac14P\sqrt d.                \tag{4.3}
\]

For a flat complete slice with normal \(u\), \(H=0\), and (1.3) becomes

\[
                              B=I-u\otimes u.                \tag{4.4}
\]

Thus (4.2) correctly retains the cube and simplex slice models: a large
contact tensor is not an error term but the smooth encoding of slice
completion.

The desired global completion theorem can now be stated without geometric
ambiguity.

> **Contact-tensor inverse.**  For an isotropic convex body of volume one
> and a balanced globally perimeter-minimizing free-boundary interface,
> either \(P\ge c\), or the contact tensor in (1.2) can be decomposed into
> complete support-to-support cells \(C_j\) with bounded reuse so that the
> associated smaller-side masses \(q_j\) satisfy
> \(\sum_jq_j\ge c\) and each cell contributes at most \(C\) times its
> relative interface area.

Together with the one-dimensional isotropic halfline estimate, this
statement gives \(P\ge c\).  It is exactly the smooth counterpart of the
central-cell incidence lemma for full polyhedral slices.

## 5. Relation to normal-flow tubes

Suppose a regular subset \(G\subset\Sigma\) supports injective normal
segments of length \(R=c_1/P\) with Jacobian bounded above and below by
universal constants.  The tube occupies a fixed amount of volume.  For any
\(c\in\mathbb R^d\), integration along the segments gives

\[
 \frac1{\mathcal H^{d-1}(G)}\int_G|x-c|^2dA
 \le C\left(\frac1{|T_G|}\int_{T_G}|y-c|^2dy+R^2\right).    \tag{5.1}
\]

If \(T_G\) has fixed volume fraction in an isotropic convex body, its
second moment is at most \(C d\) by Paouris' tail estimate; hence

\[
                   \frac1{|G|}\int_G|x-c|^2dA\le C(d+R^2). \tag{5.2}
\]

When \(R^2\le d\), this is the isotropic-scale moment required in (3.1).
When \(R^2>d\), the tube itself already violates covariance in any direction
carrying a fixed fraction of its normal matrix.  What is not automatic is
passing (5.2) from the surviving set \(G\) to all of \(\Sigma\): a small
surface fraction can lie far out or be killed at the support and can carry
the contact flux.  Formula (1.3) shows that this exceptional part must be
handled together with \(B\), not discarded as a scalar error.

## 6. Stability supplies the correct contact curvature term

For completeness, the second variation of a smooth free-boundary CMC
interface in a convex body is

\[
 \mathcal Q(u)=\int_\Sigma(|\nabla_\Sigma u|^2-|A|^2u^2)dA
 -\int_\Gamma \mathrm{II}_{\partial K}(n,n)u^2dS\ge0        \tag{6.1}
\]

for every \(u\) with \(\int_\Sigma u=0\).  Put

\[
                         m=\frac1P\int_\Sigma n\,dA
\]

and use \(u_i=n_i-m_i\).  Since
\(\sum_i|\nabla n_i|^2=|A|^2\), summing (6.1) gives

\[
 \int_\Gamma\mathrm{II}_{\partial K}(n,n)|n-m|^2dS
 \le\int_\Sigma |A|^2[\,2n\cdot m-|m|^2\,]dA.              \tag{6.2}
\]

In particular,

\[
 \int_\Gamma\mathrm{II}_{\partial K}(n,n)|n-m|^2dS
 \le(2|m|+|m|^2)\int_\Sigma|A|^2dA.                        \tag{6.3}
\]

If the normal matrix has effective rank at least \(r\), then

\[
                         |m|^2\le\|Q\|_{\mathrm{op}}\le1/r. \tag{6.4}
\]

Thus the retained rank-\(235\) branch makes the support curvature in the
interface-normal direction quantitatively subordinate to the interior
curvature energy.  This is only a contact 2-jet statement.  It does **not**
imply that the support is ruled in the interface-normal direction.  For
example, the smooth convex graph

\[
 x_d=x_1^4+\sum_{j=2}^{d-1}x_j^2
\]

has zero second fundamental form in the \(e_1\) direction along the contact
stratum \(x_1=0\), while every nontrivial line in the \(e_1\) direction
immediately leaves the graph.  A valid cell theorem must therefore control
the integral of support curvature along a transverse trajectory, or obtain
an actual chord/reach estimate from global minimality.  Equation (6.3)
alone cannot be promoted to a ruled/product classification.

## 7. Audit and stress tests

1. **Cube flat cut.**  \(H=0\), \(Q=u\otimes u\), and the contact tensor is
   exactly \(I-u\otimes u\).  The theorem sends this case to completion,
   not to a false curvature contradiction.
2. **Simplex cap.**  A flat cap again has its full dimensional contribution
   in \(B\); the smaller-side marginal density is then controlled by the
   one-dimensional log-concave lemma.
3. **Radial exponential sphere.**  In the lifted uniform model a spherical
   component has a large mean-curvature tensor term.  It cannot satisfy
   simultaneously \(|H|\le2P\), isotropic-scale position moment, and
   \(P\to0\).
4. **Product exponential max interface.**  Its polyhedral ridges are
   excluded before the smooth identity is used.  Smoothing the ridges moves
   their first-order bevel defect into curvature/contact transition zones;
   it does not make \(B\) disappear.
5. **Slabs and cylinders.**  Stable cylindrical or unduloid interfaces can
   exist in high dimension.  Therefore stability alone is not used as a
   classification theorem.  Equations (1.3) and (6.2) are exact necessary
   identities, while global perimeter minimality is still needed for the
   contact-cell inverse.

## 8. Remaining load-bearing statements

The tensor mechanism closes the smooth branch if the following two claims
are established with universal constants:

1. a boundary-position moment/contact dichotomy that upgrades the good-tube
   estimate (5.2) to the part of the interface relevant in (1.3); and
2. the contact-tensor inverse of Section 4, supplemented by a genuinely
   global estimate that converts contact curvature along whole trajectories
   into complete cells or a perimeter-saving chord competitor.  The
   pointwise stability bound (6.3) is insufficient because of the quartic
   support counterexample above.

Neither statement is a covariance trace bound.  Both are global and both
remain valid on the flat-slice models where every local curvature argument
vanishes.
