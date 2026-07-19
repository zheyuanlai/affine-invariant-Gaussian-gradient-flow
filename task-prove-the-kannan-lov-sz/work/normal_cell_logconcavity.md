# Log-concavity of normal cells and the global compatibility barrier

## 0. Outcome

Let \(\mu\) have a log-concave density \(\rho=e^{-V}\) on its affine
support, and let \(\Sigma=\partial E\) be a smooth signed-distance
interface.  At \(y\in\Sigma\), orient the normal \(n(y)\) into \(E\),
write \(S_y=D_\Sigma n(y)\), and use the normal coordinate

\[
                 x=y+t n(y).
\]

On the maximal unique-nearest, pre-focal part of this normal line, set

\[
 w_y(t)=\rho(y+t n(y))\det(I+tS_y),                         \tag{0.1}
\]

and set it equal to zero after cut, focal, or support truncation.  The
first result of this report is the exact local statement

\[
                     \boxed{w_y\text{ is log-concave on }\mathbb R.} \tag{0.2}
\]

Abrupt cut-locus truncation and the jump of a uniform density at the
boundary of a convex body do not invalidate (0.2).  At a focal endpoint
the determinant tends to zero.  These three cases are treated separately
in Theorem 1.2.

Suppose now that the cell is bisected at \(t=0\).  Put

\[
 m(y)=\int w_y(t)\,dt,\qquad r(y)={\rho(y)\over m(y)},
 \qquad q_y(t)={w_y(t)\over m(y)}.                           \tag{0.3}
\]

Thus \(r(y)=q_y(0)\) is the conditional density at its median.  If
\(\delta_y=\mathbb E_y|T|\) and \(\sigma_y^2=\operatorname{Var}_yT\),
then the following completely explicit comparisons hold:

\[
 {1\over8r}\le\delta_y\le {1\over2r},\qquad
 {1\over48r^2}\le\sigma_y^2\le {1\over2r^2}.              \tag{0.4}
\]

They are sharp in scale: the symmetric Laplace law attains both upper
bounds, while a centered uniform interval has the same reciprocal-density
scaling.  In unnormalized form, if

\[
 L_y=\int |t|w_y(t)\,dt,qquad
 Q_y=\int(t-\mathbb E_yT)^2w_y(t)\,dt,
\]

then

\[
 {m^2\over8\rho}\le L_y\le {m^2\over2\rho},\qquad
 {m^3\over48\rho^2}\le Q_y\le {m^3\over2\rho^2}.         \tag{0.5}
\]

For a smooth half-mass stationary interface every regular normal cell is
bisected.  Since \(d\eta(y)=m(y)dA(y)\) is a probability, (0.5) gives

\[
 {1\over8}\mathcal I\le J_\mu(E)\le {1\over2}\mathcal I,
 \qquad
 \mathcal I:=\int_\Sigma {m(y)^2\over\rho(y)}\,dA(y)
             =\int r^{-1}\,d\eta.                          \tag{0.6}
\]

Consequently a dimension-free bound on the proposed aggregate
\(\mathcal I\) is not a weaker local lemma: for signed-distance
maximizers it is precisely the first-moment form of KLS, up to the
constants in (0.6).

Isotropy yields the exact matrix constraint

\[
 \int r^{-2}n\otimes n\,d\eta\preceq48I_n,                 \tag{0.7}
\]

but its trace gives only

\[
                         \mathcal I\le4\sqrt{3n}.          \tag{0.8}
\]

No operator-norm replacement can repair this loss.  For the isotropic
radial exponential law and its median sphere, \(\mathcal I\asymp1\) but
the operator norm of the matrix on the left of (0.7) is \(\asymp n^{-1}\).
This is an exact log-concave, balanced, second-variation-sharp
counterexample to directional-coherence aggregation.

The signed-distance second variation becomes a weighted Poincare
inequality on the space of normal cells, with a nonnegative medial-axis
Dirichlet form.  It does not turn (0.7) into a scalar estimate.  A new
local lemma below shows that every principal direction with
\(|\kappa_i|\le r/4\) has cell metric at least \((24r)^{-1}\).  The
remaining step is nevertheless global: one must control the quotient
geometry or its medial/contact graph by ambient isotropy.  The radial
exponential sphere saturates exactly this mechanism, so no strict local
improvement is available.

The cube, Gaussian halfspace, radial exponential sphere, regular-simplex
cap, and product-exponential maximum are audited explicitly.  The last
two show that free-boundary/contact strata can carry a fixed positive
amount of mass.  In the product-exponential maximum, the exact cut-locus
truncation changes every regular cell scale; replacing it by an uncut
coordinate fiber produces a spurious aggregate of order \(n\), whereas
the correctly truncated regular-cell aggregate stays of constant order.

## 1. The normal-cell density is log-concave

Throughout this section \(K=\{V<+\infty\}^{\rm cl}\) is the convex
support.  The density is represented by the upper-semicontinuous
log-concave version \(e^{-V}\), where \(V:\mathbb R^n\to
(-\infty,+\infty]\) is convex.

### 1.1 The nearest-point interval

Let \(\Sigma\) be a \(C^2\) embedded hypersurface without boundary and
fix \(y\in\Sigma\).  If \(y\) is the unique nearest point of
\(y+t n(y)\), then for every \(z\in\Sigma\setminus\{y\}\),

\[
 |y+t n-z|^2-t^2=|z-y|^2-2t\langle n,z-y\rangle>0.         \tag{1.1}
\]

For fixed \(z\), the right side is affine in \(t\) and is positive at
\(t=0\).  Hence uniqueness persists on the segment from \(0\) to \(t\).
The unique-nearest values therefore form an interval containing zero,
up to its endpoints.  Intersect it with the interval on which
\(I+tS_y\) is positive definite and with
\(\{t:y+tn(y)\in K\}\).  Denote the resulting interval by \(I_y\).
The last set is an interval because \(K\) is convex.

At an interior unique minimizer, the tangential Hessian of squared
distance is \(2(I+tS_y)\).  Thus loss of positive definiteness can occur
only at a focal endpoint; before that endpoint all the affine factors
below are positive.

### Theorem 1.2 (log-concavity through all three truncations)

Define

\[
 w_y(t)=1_{I_y}(t)e^{-V(y+tn(y))}\det(I+tS_y).              \tag{1.2}
\]

Then \(w_y\) is a one-dimensional log-concave density, allowing zero
values and allowing an open, closed, or half-open support interval.

#### Proof

Diagonalize the self-adjoint shape operator on \(T_y\Sigma\), with
principal curvatures \(\kappa_1,\ldots,\kappa_{n-1}\).  On \(I_y\),

\[
 \log\det(I+tS_y)=\sum_{i=1}^{n-1}\log(1+t\kappa_i),
\]

and

\[
 {d^2\over dt^2}\log\det(I+tS_y)
 =-\sum_{i=1}^{n-1}{\kappa_i^2\over(1+t\kappa_i)^2}\le0. \tag{1.3}
\]

The function \(-V(y+tn(y))\) is concave on its line support.  Their sum
is concave on \(I_y\), so the product in (1.2) is log-concave there.
Multiplication by the indicator of an interval preserves log-concavity
in the extended-value convention \(\log0=-\infty\).

There are three endpoint mechanisms.

1. At a cut endpoint another nearest branch ties the branch based at
   \(y\).  The determinant may remain positive, but setting the density
   abruptly to zero is exactly interval truncation and preserves
   log-concavity.
2. At a focal endpoint one factor \(1+t\kappa_i\) tends to zero.  The
   limiting value of the Jacobian is zero; again the zero extension is
   log-concave.
3. At a support endpoint, \(V\) becomes \(+\infty\).  This includes the
   jump of the uniform density on a convex body.  Since the intersection
   of a line with the convex support is an interval, the same extended
   log-concavity argument applies.

Endpoint values do not affect the induced measure.  This proves the
claim. \(\square\)

If a cell has positive mass on both sides of zero, then log-concavity
forces \(w_y(0)>0\).  In particular a bisected cell is based in the
relative interior of the support.  A basepoint outside the support whose
normal line meets the support only on one side cannot occur in a regular
stationary balanced decomposition.

## 2. Sharp-scale one-dimensional cell comparisons

### Lemma 2.1 (tail, first moment, and variance at a median)

Let \(q\) be a nondegenerate log-concave probability density on an
interval, suppose

\[
 \int_{-\infty}^0q=\int_0^\infty q={1\over2},
 \qquad r=q(0)>0,                                          \tag{2.1}
\]

and let \(T\sim q\).  Then for \(s\ge0\),

\[
 \mathbb P(|T|\ge s)\le e^{-2rs},                          \tag{2.2}
\]

and

\[
 \|q\|_\infty\le2r.                                      \tag{2.3}
\]

Consequently

\[
 {1\over8r}\le\mathbb E|T|\le{1\over2r},\qquad
 {1\over48r^2}\le\operatorname{Var}T\le{1\over2r^2}.    \tag{2.4}
\]

#### Proof

The distribution function \(F\) and survival function \(1-F\) of a
one-dimensional log-concave law are log-concave.  Their logarithmic
derivatives at the median are \(2r\) and \(-2r\), respectively.
Concavity below a supporting tangent gives

\[
 F(-s)\le {1\over2}e^{-2rs},\qquad
 1-F(s)\le {1\over2}e^{-2rs},                              \tag{2.5}
\]

which proves (2.2).

For completeness, normalize a mode to be at \(0\) and its height to be
one.  If the median is \(a>0\), put
\(b=q(a)\) and \(c=-\log b/a\).  Concavity of \(\log q\) gives

\[
 \int_0^a q(t)\,dt\ge{1-b\over c},\qquad
 \int_a^\infty q(t)\,dt\le{b\over c}.                    \tag{2.6}
\]

The two median half-masses are equal, so \(1-b\le b\), hence
\(b\ge1/2\).  Reflection handles a median to the other side of the mode;
plateau modes follow by a limit.  Undoing the normalization proves
(2.3).

Integrating (2.2) gives

\[
 \mathbb E|T|\le\int_0^\infty e^{-2rs}\,ds={1\over2r},
 \qquad
 \mathbb ET^2\le\int_0^\infty2s e^{-2rs}\,ds={1\over2r^2}.
                                                                    \tag{2.7}
\]

On the other hand, (2.3) implies

\[
 \mathbb P(|T|\le s)\le4rs.
\]

Integration on \([0,(4r)^{-1}]\) gives
\(\mathbb E|T|\ge(8r)^{-1}\).  More generally, a probability density
bounded by \(M\) satisfies, for every \(c\in\mathbb R\),

\[
 \mathbb E(T-c)^2
 \ge\int_0^{1/(2M)}2s(1-2Ms)\,ds={1\over12M^2}.           \tag{2.8}
\]

Take \(c=\mathbb ET\), use \(M\le2r\), and combine with
\(\operatorname{Var}T\le\mathbb ET^2\).  This proves (2.4).
\(\square\)

The mean need not be zero.  The centered exponential example in Section
6.6 has a bisected cell with \(\mathbb E T=1-\log2\).

### Corollary 2.2 (unnormalized form)

For a bisected normal cell, (0.5) holds.  Equivalently,

\[
 {1\over4\sqrt2}\sigma_y\le\delta_y\le2\sqrt3\,\sigma_y. \tag{2.9}
\]

#### Proof

Substitute \(r=\rho/m\) into Lemma 2.1 and multiply the first-moment and
variance inequalities by \(m\).  Combining the opposite ends of (2.4)
gives (2.9). \(\square\)

The constants are explicit rather than claimed optimal.  Their powers
of \(r\) are optimal.  The density
\(q(t)=r e^{-2r|t|}\) has

\[
 \mathbb E|T|={1\over2r},\qquad
 \operatorname{Var}T={1\over2r^2},                        \tag{2.10}
\]

while the uniform density \(r\) on
\([-(2r)^{-1},(2r)^{-1}]\) has first moment \((4r)^{-1}\) and variance
\((12r^2)^{-1}\).

## 3. The exact covariance constraint and the \(\sqrt n\) wall

Assume the smooth normal cells cover \(\mu\) up to a null set and are
bisected.  The normal area formula gives

\[
 d\mu(y+tn(y))=q_y(t)\,dt\,d\eta(y),
 \qquad d\eta(y)=m(y)dA(y),\qquad \eta(\Sigma)=1.          \tag{3.1}
\]

Let

\[
 \bar t_y=\mathbb E_yT,qquad z_y=y+\bar t_y n(y).
\]

Conditional covariance decomposition is an identity:

\[
 \operatorname{Cov}_\mu X
 =\int \sigma_y^2n(y)\otimes n(y)\,d\eta(y)
  +\operatorname{Cov}_\eta(z_y).                           \tag{3.2}
\]

If \(\mu\) is centered and isotropic, then \(\int z_y\,d\eta=0\) and

\[
 I_n=\int\left[\sigma_y^2n\otimes n+z_y\otimes z_y\right]d\eta.
                                                                    \tag{3.3}
\]

The lower variance estimate in (0.4) proves (0.7).  Taking the trace and
using Cauchy--Schwarz gives

\[
 \mathcal I^2=\left(\int r^{-1}d\eta\right)^2
 \le\int r^{-2}d\eta
 \le48n,                                                   \tag{3.4}
\]

which is (0.8).  Symbolically, the only dimension entry is
\(\operatorname{tr}I_n=n\); no smoothing or truncation parameter is
hidden in this calculation.

The matrix statement retains more information than its trace:

\[
 \int{\langle a,n(y)\rangle^2\over r(y)^2}\,d\eta(y)
 \le48|a|^2.                                               \tag{3.5}
\]

It controls long cells whose normals have a common direction.  It says
nothing dimension-free when long cells distribute their normals among
many directions.  Section 6.3 shows that replacing the trace by the
operator norm is impossible even for a perfectly symmetric stationary
interface.

## 4. What second variation adds

Define the normalized cell metric on \(T_y\Sigma\) by

\[
 \mathsf G_y={G_y\over m(y)}
 =\int |t|(I+tS_y)^{-1}q_y(t)\,dt.                         \tag{4.1}
\]

For a regular normal speed \(h\), the signed-distance second variation
has the form

\[
 J''(h)=2\int r h^2\,d\eta
 -\int\langle\mathsf G_y\nabla_\Sigma h,
                   \nabla_\Sigma h\rangle\,d\eta
 -\mathcal M(h),                                           \tag{4.2}
\]

where \(\mathcal M(h)\ge0\) is the exact contribution of moving
codimension-one nearest-branch switching sets.  Thus a local maximum
satisfies

\[
 2\int r h^2\,d\eta
 \le\int\langle\mathsf G_y\nabla h,\nabla h\rangle\,d\eta
   +\mathcal M(h),
 \qquad \int r h\,d\eta=0.                                \tag{4.3}
\]

The medial term is on the right of the stability inequality.  Dropping
it makes the assertion stronger and is not legitimate.

### Lemma 4.1 (a quantitative low-curvature cell metric)

Let \(e_i\) be a unit principal direction with curvature \(\kappa_i\).
If

\[
                         |\kappa_i|\le {r(y)\over4},       \tag{4.4}
\]

then

\[
                  \langle\mathsf G_y e_i,e_i\rangle
                  \ge {1\over24r(y)}.                     \tag{4.5}
\]

#### Proof

Lemma 2.1 gives \(\mathbb E|T|\ge(8r)^{-1}\).  Its tail bound gives,
with \(L=2/r\),

\[
 \mathbb E[|T|1_{\{|T|>L\}}]
 \le\left(L+{1\over2r}\right)e^{-2rL}
 ={5e^{-4}\over2r}<{1\over16r}.                           \tag{4.6}
\]

Therefore

\[
 \int_{|t|\le2/r}|t|q_y(t)\,dt\ge{1\over16r}.            \tag{4.7}
\]

On this interval, (4.4) gives
\(0<1+t\kappa_i\le3/2\).  The principal component of (4.1) is

\[
 \int {|t|\over1+t\kappa_i}q_y(t)\,dt,
\]

so restricting it to \(|t|\le2/r\) and using (4.7) proves (4.5).
\(\square\)

This produces a precise dichotomy: a tangent direction either contributes
at least \((24r)^{-1}\) to the quotient metric, or its principal
curvature is at least \(r/4\).  It still does not bound
\(\int r^{-1}d\eta\).  Inequality (4.3) is a Poincare inequality on the
actual quotient/medial cell complex.  Ambient covariance controls the
normal conditional variances in (3.3), but does not control the Poincare
constant of that quotient.  Proving such control with a universal
constant is the missing global compatibility theorem.

### 4.2 Exact equivalence of the proposed aggregate with the target

For every smooth bisected-cell interface, (0.6) is an exact two-sided
comparison.  For a global signed-distance maximizer,

\[
                         \mathcal I\le C                   \tag{GCC}
\]

would give \(D_1(\mu)\le C/2\), after the standard approximation of a
nonsmooth maximizing interface.  Conversely, a dimension-free first
moment bound gives \(J_\mu(E)\le C\) and hence
\(\mathcal I\le8C\) for every regular bisected-cell interface.  Thus
(GCC), including the approximation/contact-stratum passage, is
equivalent up to universal constants to T3.  It may not be inserted as a
``cell aggregation lemma'' without proof.

There is also a purely algebraic warning.  Take \(n\) abstract cells of
mass \(1/n\), normals \(e_i\), conditional symmetric Laplace laws with
\(r=n^{-1/2}\), and variance \(n/2\).  Then

\[
 \int\sigma^2 n\otimes n\,d\eta={1\over2}I_n,
 \qquad \int r^{-1}d\eta=\sqrt n.                         \tag{4.8}
\]

A medial graph form can be chosen to make (4.3) an equality.  This is not
a log-concave geometric example and is not a counterexample to KLS; it
shows exactly why (0.7), cell log-concavity, and the sign
\(\mathcal M\ge0\) do not algebraically imply (GCC).  The genuinely
geometric obstructions are next.

## 5. Exact failures of tempting aggregation steps

The following statements are all false with a universal constant.

1. **Directional coherence.**  One cannot assert
   \(\mathcal I^2\le C\|\int r^{-2}n\otimes n\,d\eta\|_{\rm op}\), nor
   the same estimate with normalized trace.  The radial exponential
   sphere in Section 6.3 has left side of constant order and operator
   norm of order \(1/n\).
2. **Cellwise centering.**  A median-bisected log-concave cell need not
   have \(\mathbb E_yT=0\).  The isotropic shifted exponential halfspace
   in Section 6.6 has mean \(1-\log2\).
3. **Global half mass implies regular-cell balance.**  It does not for a
   nonsmooth or free-boundary interface.  The simplex cap and the
   product-exponential maximum have half mass but a fixed amount of
   positive-volume normal cells based on lower-dimensional contact/ridge
   strata.  Their regular cells are not bisected.
4. **Singular strata are negligible because their base has zero surface
   area.**  Their normal cones can have full dimension.  In both examples
   just mentioned, their Voronoi region has mass tending to
   \((1-\log2)/2\).
5. **Small largest principal curvature forces coherent normals.**  On the
   radial exponential median sphere, every principal curvature is
   \(\asymp n^{-1/2}\), the normals are uniform on the sphere, and the
   degree-one second-variation inequality is an equality.  The summed
   quotient geometry, not a pointwise curvature bound, is what pays for
   stability.

These counterexamples do not disprove (GCC).  They eliminate proposed
shortcuts from local cell facts to (GCC).

## 6. Mandatory model audit

### 6.1 Isotropic cube and a coordinate halfspace

Let \(\mu\) be uniform on \([-a,a]^n\), \(a=\sqrt3\), and
\(E=\{x_1>0\}\).  For \(y=(0,x_2,\ldots,x_n)\),

\[
 w_y(t)=\rho 1_{[-a,a]}(t),\quad m=2a\rho,\quad
 r={1\over2a},\quad \delta={a\over2},\quad \sigma^2=1.   \tag{6.1}
\]

Thus

\[
 \mathcal I=2a=2\sqrt3,qquad J_\mu(E)={a\over2}={\sqrt3\over2},
 \qquad \int\sigma^2n\otimes n\,d\eta=e_1\otimes e_1.   \tag{6.2}
\]

The abrupt support cutoff in (6.1) is log-concave and verifies the support
part of Theorem 1.2.  Here \(S=0\) and
\(\mathsf G=(a/2)I_{T\Sigma}\); the cross-sectional cube Poincare
inequality gives strict second-variation stability.

### 6.2 Gaussian halfspace

For \(\gamma_n\) and \(E=\{x_1>0\}\),

\[
 q_y(t)=\varphi(t),\quad r=\varphi(0)={1\over\sqrt{2\pi}},
 \quad\delta=\sqrt{2\over\pi},\quad\sigma^2=1.            \tag{6.3}
\]

Consequently

\[
 \mathcal I=\sqrt{2\pi},\qquad
 J_{\gamma_n}(E)=\sqrt{2/\pi},\qquad
 \mathsf G=\sqrt{2/\pi}\,I_{T\Sigma}.                    \tag{6.4}
\]

The quotient is \(\gamma_{n-1}\), and (4.3) is exactly its Gaussian
Poincare inequality.  Linear quotient functions give equality.

### 6.3 Isotropic radial exponential and its median sphere

Let

\[
 d\mu(x)=c_ne^{-\lambda|x|}\,dx,qquad\lambda=\sqrt{n+1},
                                                                    \tag{6.5}
\]

so \(R=|X|\sim\operatorname{Gamma}(n,\lambda)\) and \(\mu\) is
isotropic.  Let \(r_0\) be the median of \(R\) and
\(E=\{|x|<r_0\}\).  With inward normal \(n(y)=-y/r_0\) and
\(t=r_0-R\),

\[
 q_y(t)=p_R(r_0-t)1_{\{t<r_0\}},qquad
 r(y)=p_R(r_0),qquad \sigma_y^2={n\over n+1}.             \tag{6.6}
\]

The focal endpoint is \(t=r_0\), where
\(\det(I+tS)=(1-t/r_0)^{n-1}\) vanishes.  Formula (1.3) proves
log-concavity all the way to that endpoint.  Lemma 2.1 and
\(n/(n+1)\asymp1\) give

\[
                         r(y)\asymp1,qquad\mathcal I\asymp1. \tag{6.7}
\]

Rotational symmetry yields the exact matrix

\[
 \int\sigma_y^2n(y)\otimes n(y)\,d\eta(y)
 ={1\over n+1}I_n.                                        \tag{6.8}
\]

The same is true, up to universal factors, with \(\sigma_y^2\) replaced
by \(r(y)^{-2}\).  Hence its operator norm is \(\asymp n^{-1}\) while
\(\mathcal I\asymp1\), proving the directional-coherence failure.

Every principal curvature has magnitude \(1/r_0\asymp n^{-1/2}\).
Nevertheless, writing \(G_y=g_yI_{T_y\Sigma}\), the degree-one spherical
mode satisfies the exact identity

\[
 {n-1\over r_0^2}\,g_y=2\rho(y),                           \tag{6.9}
\]

so (4.3) is an equality for infinitesimal translations of the sphere.
This model blocks any proposed strict gain from the local second
variation.

### 6.4 Regular-simplex half-volume cap and free-boundary contact mass

Use barycentric coordinates
\((U_1,\ldots,U_{n+1})\sim\operatorname{Dirichlet}(1,\ldots,1)\) on a
regular simplex.  Isotropic normalization is a scalar on its affine
span, so the following orthogonal-projection calculation is unchanged.
Put

\[
 s_n=1-2^{-1/n},\qquad E=\{U_1\ge s_n\}.                   \tag{6.10}
\]

Since \(\mathbb P(U_1\ge s)=(1-s)^n\), this cap has mass one half.  Its
relative interface is the simplex

\[
 \Sigma_s=\{U_1=s,\ U_i\ge0,\ \sum_{i=2}^{n+1}U_i=1-s\}.
\]

For a point with \(U_1=x<s\), orthogonal projection to the complete
cutting plane has coordinates

\[
 U'_1=s,qquad U'_i=U_i-{s-x\over n}\quad(i\ge2).          \tag{6.11}
\]

It lands in the regular face precisely when
\(U_i\ge(s-x)/n\) for every \(i\ge2\).  Conditional on \(U_1=x\), a
simplex translation gives

\[
 \mathbb P(\text{regular projection}\mid U_1=x)
 =\left({1-s\over1-x}\right)^{n-1}.                       \tag{6.12}
\]

Since the density of \(U_1\) is \(n(1-x)^{n-1}\), the complement-side
mass assigned to regular normal cells is

\[
                  n s(1-s)^{n-1}.                         \tag{6.13}
\]

All cap-side points have regular projection.  Therefore the mass of
Voronoi cells based on the lower-dimensional free boundary is exactly

\[
 c_n={1\over2}-n s_n(1-s_n)^{n-1}
 \longrightarrow {1-\log2\over2}>0.                       \tag{6.14}
\]

On a regular normal line based at barycentric point
\((s,z_2,\ldots,z_{n+1})\), the two uniform chord lengths, in the common
normal parameter, are

\[
                         s\quad\text{and}\quad n\min_{i\ge2}z_i.
                                                                    \tag{6.15}
\]

They are unequal for almost every basepoint.  Thus global half mass does
not imply regular-cell balance for this free-boundary cap.  The missing
mass is exactly carried by normal cones based on the boundary of
\(\Sigma_s\).  Treating that base as surface-null would delete the fixed
mass (6.14).

### 6.5 Product exponential maximum and ridge mass

Let \(Y_i\) be independent \(\operatorname{Exp}(1)\) variables and
\(X_i=Y_i-1\).  Then \(X=(X_i)\) is centered, isotropic, and log-concave.
Set

\[
 z_n=-\log(1-2^{-1/n}),\qquad
 E=\{\max_iY_i\le z_n\}.                                  \tag{6.16}
\]

This set has mass one half.  Write

\[
 p_n=e^{-z_n}=1-2^{-1/n}
\]

and let \(N\) count the coordinates exceeding \(z_n\).  Then
\(N\sim\operatorname{Binomial}(n,p_n)\),

\[
 \mathbb P(N=0)={1\over2},\qquad
 \mathbb P(N=1)=np_n(1-p_n)^{n-1}.                        \tag{6.17}
\]

Points with \(N\ge2\) have nearest boundary point on a ridge where at
least two facets meet.  The mass of these singular normal cones is

\[
 {1\over2}-np_n(1-p_n)^{n-1}
 \longrightarrow {1-\log2\over2}.                         \tag{6.18}
\]

On a regular facet based at \(Y_i=z_n\), write the other coordinates as
\(u_j\in[0,z_n]\), and put

\[
 M=\max_{j\ne i}u_j,\qquad d=z_n-M.                       \tag{6.19}
\]

Use \(t=Y_i-z_n\).  The normal branch ceases to be uniquely nearest when
\(Y_i\) ties the largest other coordinate, so the **exact** normal-cell
density, apart from a base-dependent constant \(A\), is

\[
                         w(t)=Ae^{-t}1_{[-d,\infty)}(t).   \tag{6.20}
\]

Thus \(m=Ae^d\), \(r=e^{-d}\), and the two masses are in the ratio

\[
 \int_{-d}^0e^{-t}dt:e^0=(e^d-1):1=(1-r):r.               \tag{6.21}
\]

They are equal only on the codimension-one base locus \(d=\log2\); they
are not equal almost everywhere for any \(n>1\).  This cell also gives a
clean instance of abrupt cut-locus truncation in Theorem 1.2.

The correctly truncated regular-cell aggregate can be evaluated exactly.
Since the base density is \(e^{-z_n-\sum u_j}\),

\[
 \begin{split}
 \mathcal I_{\rm reg}
 &=n e^{z_n}\int_{[0,z_n]^{n-1}}
       e^{-\sum_j u_j-2\max_j u_j}\,du\\
 &=n(n-1)e^{z_n}\int_0^{1-p_n}u^{n-2}(1-u)^2\,du.         \tag{6.22}
 \end{split}
\]

With \(v=n(1-u)\), dominated convergence gives

\[
 \mathcal I_{\rm reg}\longrightarrow
 {1\over\log2}\int_{\log2}^{\infty}v^2e^{-v}\,dv
 ={(\log2)^2+2\log2+2\over2\log2}.                        \tag{6.23}
\]

So the actual regular aggregate is of constant order.  If one retains
the mass of the assigned regular cells but incorrectly replaces their
cut endpoint \(-d\) by the support endpoint \(-z_n\), one assigns the
common false density \(r=p_n\) and obtains the spurious value

\[
 {\mathbb P(N\le1)\over p_n}
 ={\frac12+np_n(1-p_n)^{n-1}\over p_n}\asymp n.           \tag{6.24}
\]

That value is not a geometric cell integral: it is precisely the error
caused by ignoring nearest-branch switching.  The example is not a
counterexample to (GCC): the interface is nonsmooth, its regular cells
are not balanced, and it is not asserted to maximize \(J_\mu\).  It is
an exact warning that smoothing a ridge requires quantitative control of
the resulting high-curvature/contact packet.

### 6.6 A genuinely asymmetric balanced cell

For the same centered exponential law in one coordinate, take the
coordinate halfspace at its median.  In the normal coordinate
\(t=Y_1-\log2\),

\[
 q(t)={1\over2}e^{-t}1_{[-\log2,\infty)}(t),qquad
 r={1\over2},\qquad \sigma^2=1,qquad
 \mathbb ET=1-\log2.                                     \tag{6.25}
\]

The cell is exactly bisected but not centered.  Formula (3.2), rather
than the false substitution \(z_y=y\), is necessary in asymmetric
examples.

## 7. Formal separation of proved and unproved statements

The following are local or exact and have been proved above:

* normal-cell log-concavity through cut, focal, and support truncation;
* the universal median-density comparisons (0.4)--(0.5);
* the mixture covariance identity (3.2) and matrix bound (0.7);
* the dimension-dependent consequence (0.8);
* the low-curvature cell-metric estimate (4.5);
* all five model calculations, including the exact contact masses.

The following statement has **not** been proved:

\[
 \text{actual global maximality + isotropy + log-concavity}
 \quad\Longrightarrow\quad
 \int r^{-1}d\eta\le C.                                  \tag{7.1}
\]

By (0.6), (7.1), with the required smoothing and contact-stratum limit,
is KLS-equivalent.  The second variation rewrites the missing information
as control of a weighted quotient/medial Poincare problem, but does not
solve it.  The radial sphere proves that diffuse normals and sharp local
stability can coexist, while the simplex and product maximum prove that
singular normal cones cannot be discarded.  Any completion must therefore
establish a genuinely global compatibility theorem for the quotient and
contact geometry; none of the local log-concave-cell estimates above can
be promoted to that theorem by trace, operator norm, or omission of the
medial term.
