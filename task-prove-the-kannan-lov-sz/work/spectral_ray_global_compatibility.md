# Spectral blocks and global compatibility of balanced transport rays

## 0. Audit verdict

There are two rigorous improvements, but they do not close the inverse.

First, the global eikonal structure gives an exact pairwise packing law
which is stronger than the local focal estimate.  If two transport rays
cross the common zero level at \(y,z\), have directions \(u,v\), and both
contain the parameter interval \([-R,R]\), then

\[
                         \boxed{R|u-v|\le |y-z|.}     \tag{0.1}
\]

This follows from the two *cross-endpoint* Lipschitz inequalities and is
valid without smoothness.  It says that the orientation map on the long-ray
part of the zero level is \(R^{-1}\)-Lipschitz; equivalently, the basepoint
metric is \(R\)-co-Lipschitz with respect to orientation distance.
Transport cyclic monotonicity gives no additional pairwise constraint: each
cross edge already has length at least \(2R\) by the Kantorovich potential.

Second, the bounded spectral multiplicity theorem in
`terminal_needle_inverse.md` can be sharpened.  In an \(r\)-dimensional
spectral subspace the projected conditional standard deviation \(\tau\)
satisfies

\[
                            \mathbb E\tau^4\le C r,   \tag{0.2}
\]

not merely \(\mathbb E\tau^2\le r\).  Using this fourth-moment tail in the
inward-ball argument improves

\[
 N_K(\alpha)\le C\alpha^{-9/2}
 \quad\hbox{to}\quad
 \boxed{N_K(\alpha)\le C\alpha^{-17/5}.}             \tag{0.3}
\]

Consequently

\[
 \operatorname{tr}K\le Cn^{12/17},\qquad
 \mathbb E{1\over\sigma}\ge c n^{-6/17}.             \tag{0.4}
\]

This is only an audit checkpoint.  Dyadic summation still diverges, so
(0.3) is not a returnable KLS result.

The hoped-for operator packing step does not follow from the endpoint
Jacobian or cyclic monotonicity.  On a smooth long-ray chart, simultaneous
log-concavity at \(+R\) and \(-R\) gives exactly

\[
 R^2\|D_\Sigma u\|_{HS}^2\le C,                     \tag{0.5}
\]

which is the already known Frobenius focal bound.  It gives no lower
Jacobian for the Gauss map: flat orientation patches have
\(D_\Sigma u=0\).  Distinct flat patches can meet on singular transition
strata.

Section 6 gives a rigorous stress test with all of the following features:

* \(\mu\) is the isotropic uniform measure on a cube, hence globally
  log-concave;
* \(f\) is one globally defined 1-Lipschitz eikonal potential;
* its ray interiors are nonbranching and every conditional ray is exactly
  balanced and log-concave;
* \(K\) has \(m\) equal nonzero eigenvalues and normalized orientation rank
  \(m\); and
* the Gauss derivative is zero on every regular orientation patch.

The model does satisfy the desired inverse, very strongly:
\(\mathbb E(1/\sigma)=m\).  Its transition cells shorten the rays.  Thus it
does not contradict KLS; it disproves the narrower claim that diffuse
spectral orientation, global log-concavity, and a global eikonal congruence
alone force a Jacobian or log-concavity violation.  The fact not shared by
this stress test is that \(f\) is a global maximizer of the first-moment
functional.

The remaining load-bearing statement is therefore a **weighted transition
inequality**: substantial projective-orientation variance among rays on which
\(\sigma\) is large must force a positive amount of reciprocal-scale mass in
the transition strata.  Neither (0.1), (0.5), area formula, nor current
spectral multiplicity proves it.

## 1. Setting

Let \(\mu\) be isotropic and log-concave.  Let \(f\) be an extremal
1-Lipschitz potential for the mean-centered first-moment functional, shifted
so that \(\mu f=0\).  Write the transport-ray disintegration as

\[
 \mu=\int\nu_q\,\eta(dq),qquad
 X=y_q+T_q u_q,qquad f(X)=T_q,                       \tag{1.1}
\]

where \(y_q\in\{f=0\}\), \(|u_q|=1\), and the conditional law of \(T_q\)
is one-dimensional log-concave.  The optimal-cut argument gives one number
\(p\in(0,1)\) such that

\[
                         \nu_q(T_q>0)=p               \tag{1.2}
\]

for almost every active ray.  Put

\[
 m_q=\mathbb E_qT_q,qquad
 \sigma_q^2=\operatorname{Var}_qT_q,qquad
 b_q=y_q+m_qu_q.                                      \tag{1.3}
\]

Total covariance gives

\[
 I=\operatorname{Cov}_\eta(b_q)+K,qquad
 K=\mathbb E_\eta[\sigma_q^2u_qu_q^T]\preceq I.      \tag{1.4}
\]

The target inverse is

\[
                         \mathbb E_\eta{1\over\sigma_q}\ge c. \tag{1.5}
\]

When \(p\in[\delta,1-\delta]\), all constants in the geometric statements
below may depend on \(\delta\).  If one wants a universal proof through this
route, the dependence on the optimizing cut probability must ultimately be
removed or \(p\) must first be bounded away from the endpoints.

## 2. Exact cross-endpoint packing

### Theorem 2.1 (global eikonal nonintersection inequality)

Let \(f:\mathbb R^n\to\mathbb R\) be 1-Lipschitz.  Suppose

\[
 f(y+tu)=t,qquad f(z+tv)=t                         \tag{2.1}
\]

for all \(t\in[-R,R]\), where \(f(y)=f(z)=0\) and \(u,v\) are unit vectors.
Then

\[
                             R|u-v|\le|y-z|.          \tag{2.2}
\]

#### Proof

Apply the 1-Lipschitz inequality to the two pairs with opposite labels:

\[
\begin{aligned}
 |(y+Ru)-(z-Rv)|&\ge 2R,\\
 |(y-Ru)-(z+Rv)|&\ge 2R.                              \tag{2.3}
\end{aligned}
\]

Put \(d=z-y\).  Squaring (2.3) gives

\[
 |d|^2+R^2|u+v|^2\pm2R\langle d,u+v\rangle\ge4R^2. \tag{2.4}
\]

Take the smaller left side and use
\(|u+v|^2=4-|u-v|^2\).  Then

\[
 |d|^2-R^2|u-v|^2
 \ge2R|\langle d,u+v\rangle|\ge0,                   \tag{2.5}
\]

which is (2.2).  \(\square\)

The unequal-length version uses

\[
 R=\min\{R_y^+,R_y^-,R_z^+,R_z^-\},                 \tag{2.6}
\]

where \([-R_q^-,R_q^+]\) is contained in the parameter interval of the
\(q\)-th ray.

For a one-dimensional log-concave density of variance \(\sigma^2\), the
maximum density is at most \(C/\sigma\).  Under (1.2), if
\(p\in[\delta,1-\delta]\), each side of zero in the conditional support has
length at least \(c\delta\sigma\).  Therefore, on the set

\[
                         H_L=\{q:\sigma_q\ge L\},    \tag{2.7}
\]

Theorem 2.1 gives

\[
 |u_q-u_{q'}|\le {C_\delta\over L}|y_q-y_{q'}|
 \qquad(q,q'\in H_L).                                \tag{2.8}
\]

In particular, for two independent rays conditioned on \(H_L\),

\[
 L^2\operatorname{tr}\operatorname{Cov}(u_q\mid H_L)
 \le C_\delta\operatorname{tr}
                    \operatorname{Cov}(y_q\mid H_L). \tag{2.9}
\]

This is a genuine global compatibility relation.  It is nevertheless only
a trace estimate.  Since \(y_q=b_q-m_qu_q\) and
\(|m_q|\le C_\delta\sigma_q\), isotropy permits the right side of (2.9) to
be of order \(n/\eta(H_L)\).  Thus (2.9) still permits
\(L=n^{1/4}\) and a direction cloud spread through \(n\) coordinates.

### 2.1 Cyclic monotonicity is redundant here

Take endpoint pairs

\[
 x_i^+=y_i+Ru_i,qquad x_i^-=y_i-Ru_i.                \tag{2.10}
\]

The two-cycle inequality for the distance cost is

\[
 4R\le |x_1^+-x_2^-|+|x_2^+-x_1^-|.                  \tag{2.11}
\]

But every cross term on the right is already at least \(2R\), because
\(f(x_i^+)-f(x_j^-)=2R\) and \(f\) is 1-Lipschitz.  The same observation
applies to every edge of every higher transport cycle.  Hence cost-cyclic
monotonicity supplies no operator-valued strengthening of (2.2) for these
equal-label endpoint packets.

## 3. What the two endpoint Jacobians actually imply

Work on a \(C^2\) ray chart.  Let \(\Sigma=\{f=0\}\), \(u=N=\nabla f\), and

\[
                         S=D_\Sigma N.                \tag{3.1}
\]

The eikonal equation makes \(S\) self-adjoint on the tangent space and the
normal maps have differentials

\[
 D F_{\pm R}=I\pm RS.                                 \tag{3.2}
\]

If all rays in the chart survive to both \(\pm R\), every eigenvalue
\(\kappa_i\) of \(S\) obeys \(|R\kappa_i|<1\), and

\[
 J_+(y)J_-(y)
 =\det(I+RS)\det(I-RS)
 =\det(I-R^2S^2).                                     \tag{3.3}
\]

Let \(\rho=e^{-V}\).  The conditional density ratio is

\[
 {q_y(R)q_y(-R)\over q_y(0)^2}
 ={\rho(y+Ru)\rho(y-Ru)\over\rho(y)^2}
       \det(I-R^2S^2).                                \tag{3.4}
\]

Log-concavity gives

\[
                 \rho(y+Ru)\rho(y-Ru)\le\rho(y)^2.  \tag{3.5}
\]

For a balanced one-dimensional log-concave density, choose
\(R=c_\delta\sigma_y\) with \(c_\delta>0\) sufficiently small.  Standard
central-density and slope bounds give

\[
                         q_y(\pm R)\ge c'_\delta q_y(0). \tag{3.6}
\]

Combining (3.3)--(3.6),

\[
 (c'_\delta)^2\le\det(I-R^2S^2).                    \tag{3.7}
\]

Since \(-\log(1-x)\ge x\) for \(0\le x<1\),

\[
 \boxed{R^2\|S\|_{HS}^2
       \le-2\log c'_\delta.}                         \tag{3.8}
\]

Thus simultaneous control of the two endpoint Jacobians reproduces the
Frobenius estimate already obtained by integrating the one-dimensional
log-density curvature.  It does not yield

\[
                         |\det S|\ge cR^{-(n-1)}      \tag{3.9}
\]

or any lower bound on an exterior-power Jacobian of the Gauss map.  A flat
patch has \(S=0\), \(J_+=J_-=1\), and satisfies every inequality above.

The radial example tests the opposite extreme.  For a sphere of radius
\(r\simeq\sqrt n\),

\[
                         \|S\|_{HS}^2={n-1\over r^2}\simeq1. \tag{3.10}
\]

Equation (3.8) forces a balanced radial conditional scale to be \(O(1)\).
Thus the determinant mechanism correctly excludes a long spherical bundle.
It is powerless on flat polyhedral or cylindrical patches, exactly where
orientation changes can be concentrated on lower-dimensional transition
strata.

## 4. Fourth-moment improvement of spectral multiplicity

This section records the best unconditional spectral improvement found in
this audit.

### Theorem 4.1

Let an isotropic log-concave probability admit a line disintegration
satisfying (1.4).  Let \(K\) have \(r\) eigenvalues at least
\(\alpha\in(0,1]\).  Then

\[
                             r\le C\alpha^{-17/5}.    \tag{4.1}
\]

#### Proof

Let \(E\) be the corresponding spectral subspace and define

\[
 Y=P_EX,qquad y_q=P_Eb_q,qquad
 \tau_q=\sigma_q|P_Eu_q|.                            \tag{4.2}
\]

The marginal \(Y\) is isotropic and log-concave in \(E\simeq\mathbb R^r\).
As in Theorem 8.1 of `terminal_needle_inverse.md`, put

\[
 s=\sqrt{1-\alpha/2},qquad d=1-s,qquad
 q_0={\alpha\over2-\alpha}.                           \tag{4.3}
\]

The covariance identity gives

\[
 \mathbb P\{|y_q|\le s\sqrt r\}\ge q_0.             \tag{4.4}
\]

Conditionally on \(q\), \(Y\) lies on a line, has barycenter \(y_q\), and
has standard deviation \(\tau_q\).  The one-dimensional quadratic lemma
applied to \(|Y|^2\) gives

\[
 \operatorname{Var}(|Y|^2\mid q)\ge {1\over100}\tau_q^4. \tag{4.5}
\]

The law of total variance and the quadratic form of the thin-shell theorem
give

\[
 \boxed{\mathbb E\tau_q^4\le C_4r}                   \tag{4.6}
\]

with a universal \(C_4\).  Set

\[
 M=\left({2C_4r\over q_0}\right)^{1/4}.              \tag{4.7}
\]

Then \(\mathbb P\{\tau_q>M\}\le q_0/2\).  Hence the event

\[
 G=\{|y_q|\le s\sqrt r,\ \tau_q\le M\}              \tag{4.8}
\]

has probability at least \(q_0/2\).

There is a universal \(c_0>0\) such that a mean-zero variance-one
one-dimensional log-concave law gives mass at least
\(c_0\min(a,1)\) to \([-a,a]\).  Choose

\[
 a=\min\left\{1,{d\sqrt r\over2M}\right\}.           \tag{4.9}
\]

On \(G\), the central interval of radius \(a\tau_q\) about \(y_q\) is
contained in the ball of radius

\[
 (s+d/2)\sqrt r=(1-d/2)\sqrt r.                       \tag{4.10}
\]

Therefore

\[
 \mathbb P\{|Y|\le(1-d/2)\sqrt r\}
 \ge {c_0q_0\over2}\min(a,1).                        \tag{4.11}
\]

Thin shell and Chebyshev give the upper bound

\[
 \mathbb P\{|Y|\le(1-d/2)\sqrt r\}
 \le {4C_{TS}\over d^2r}.                            \tag{4.12}
\]

If \(a=1\), (4.11)--(4.12) give
\(r\le C d^{-2}q_0^{-1}\).  If \(a<1\), equations (4.7) and (4.9) give

\[
 a={d q_0^{1/4}r^{1/4}\over2(2C_4)^{1/4}}.           \tag{4.13}
\]

Substitution in (4.11)--(4.12) yields

\[
 r^{5/4}\le C d^{-3}q_0^{-5/4},qquad
 r\le C d^{-12/5}q_0^{-1}.                           \tag{4.14}
\]

For \(0<\alpha\le1\), \(d\ge\alpha/4\) and
\(q_0\ge\alpha/2\).  Both cases are bounded by
\(C\alpha^{-17/5}\), proving (4.1).  \(\square\)

### Corollary 4.2

In ambient dimension \(n\),

\[
\begin{aligned}
 \operatorname{tr}K
 &=\int_0^1N_K(\alpha)\,d\alpha\\
 &\le\int_0^1\min\{n,C\alpha^{-17/5}\}\,d\alpha
 \le C'n^{12/17}.                                    \tag{4.15}
\end{aligned}
\]

Since \(\mathbb E\sigma^2=\operatorname{tr}K\), Cauchy--Schwarz gives

\[
 \mathbb E{1\over\sigma}
 \ge {1\over\mathbb E\sigma}
 \ge {1\over\sqrt{\mathbb E\sigma^2}}
 \ge c n^{-6/17}.                                    \tag{4.16}
\]

This improves the exponent in the earlier terminal audit but remains
dimension dependent.

## 5. Why dyadic spectral blocks still do not sum

Let

\[
 E_j=\mathbf1_{(2^{-j-1},2^{-j}]}(K)\mathbb R^n,qquad
 r_j=\dim E_j.                                       \tag{5.1}
\]

Theorem 4.1 gives

\[
                         r_j\le C2^{17j/5}.           \tag{5.2}
\]

The trace allowed on this block is therefore

\[
 \operatorname{tr}(P_{E_j}K)\le2^{-j}r_j
                         \le C2^{12j/5},              \tag{5.3}
\]

which grows with the scale index.  Cutting off only when \(r_j\simeq n\)
reproduces (4.15).  Thus a threshold-by-threshold application of thin shell
cannot give a summable spectral series.

The global packing law (2.2) also does not by itself repair (5.3).  It says
that a block of rays of length \(L\) and orientation variance of order one
requires basepoint variance of order \(L^2\).  Isotropy has total basepoint
variance of order \(n\), so it permits \(L=n^{1/4}\) at \(n\) diffuse
orientations.  To improve (5.3), one needs an operator- or incidence-valued
statement which charges the *transition between orientation blocks*, not
only their pairwise separation.

## 6. A globally log-concave, globally eikonal checkerboard

This example is a countermodel to a purely Jacobian-based orientation
argument.  It also shows exactly how transition strata can pay the desired
reciprocal scale.

Let \(a=\sqrt3\), let \(\mu\) be uniform on the cube
\([-a,a]^n\), and fix \(2\le m\le n\).  This measure is isotropic and
log-concave.  Define, with the value zero on the coordinate hyperplanes,

\[
 f_m(x)=\left(\prod_{i=1}^m\operatorname{sgn}x_i\right)
             \min_{1\le i\le m}|x_i|.                \tag{6.1}
\]

On the region where \(j\) uniquely minimizes \(|x_i|\),

\[
 f_m(x)=\left(\prod_{i\ne j}\operatorname{sgn}x_i\right)x_j,
 \qquad
 \nabla f_m=\left(\prod_{i\ne j}\operatorname{sgn}x_i\right)e_j. \tag{6.2}
\]

The function is continuous and piecewise affine, and
\(|\nabla f_m|=1\) almost everywhere.  Hence it is globally 1-Lipschitz.

Fix \(j\) and the other coordinates.  Put

\[
 H=\min_{i\ne j,\ i\le m}|x_i|.                      \tag{6.3}
\]

The maximal ray through that point is obtained by varying \(x_j\) through
\([-H,H]\), with orientation chosen as in (6.2).  Along it,

\[
                         f_m=t,qquad T\sim\operatorname{Unif}[-H,H]. \tag{6.4}
\]

Ray interiors are disjoint; ties occur only at endpoints and on a null
transition set.  Every ray is exactly half-balanced, and reflection
\(t\mapsto-t\) is an optimal cross-cut coupling.  Thus \(f_m\) is an exact
Kantorovich potential for its own balanced cut.

The quotient law can be computed explicitly.  Conditional on the projective
orientation \(j\), the \(m-1\) other absolute coordinates are initially
uniform on \([0,a]\), and the ray quotient weights their minimum by its ray
length \(2H\).  Therefore

\[
                         {H\over a}\sim\operatorname{Beta}(2,m-1) \tag{6.5}
\]

under \(\eta\), and the projective direction is uniform on
\(\{e_1,\ldots,e_m\}\).  Since \(\sigma=H/\sqrt3=H/a\),

\[
 \mathbb E\sigma^2={6\over(m+1)(m+2)},qquad
 \mathbb E{1\over\sigma}=m.                          \tag{6.6}
\]

Moreover

\[
 \boxed{K={6\over m(m+1)(m+2)}P_{\operatorname{span}(e_1,\ldots,e_m)}.}
                                                               \tag{6.7}
\]

Thus \(K\) has \(m\) equal nonzero eigenvalues of order \(m^{-3}\), while
the orientation covariance has full rank \(m\).  On every regular patch
\(D_\Sigma u=0\), so every regular Gauss Jacobian vanishes.  All orientation
change is carried by the coordinate-hyperplane transition strata.

The reciprocal scale is large because a typical ray reaches a transition
after the minimum of \(m-1\) independent transverse distances.  In fact the
cut perimeter is

\[
 \mu^+(\{f_m>0\})={m\over2\sqrt3},qquad
 {\mu^+(\{f_m>0\})\over1/2}={m\over\sqrt3},           \tag{6.8}
\]

in agreement with the one-dimensional estimate and (6.6).

This is a fully global log-concave/eikonal/balanced model, so diffuse
spectral orientation alone does not create a log-concavity violation.  It
is not a counterexample to the desired theorem.  Indeed

\[
 \int|f_m|\,d\mu={a\over m+1},                        \tag{6.9}
\]

whereas a coordinate linear functional has mean absolute deviation
\(a/2\).  For \(m>1\), \(f_m\) is not the global maximizer of the
first-moment functional.  Any successful global transition theorem may
therefore have to use this extremality, rather than only balance and
Kantorovich optimality for the induced cut.

## 7. The exact remaining compatibility statement

Let

\[
                         A=\mathbb E_\eta{1\over\sigma_q}. \tag{7.1}
\]

If \(A\) were small, Markov's inequality would put most quotient mass on
rays with a fixed large lower bound on \(\sigma\).  Proposition 6.1 of the
terminal audit already closes the case in which a fixed positive fraction
of those directions lies in one fixed projective cap.  It would therefore
suffice to prove a statement of the following form.

For an event \(H\), define its projective dispersion by

\[
 Q_H=\mathbb E[u_qu_q^T\mid H],\qquad
 \mathcal V_{\rm proj}(H)
 =\mathbb E[\|u_qu_q^T-Q_H\|_{HS}^2\mid H]
 =1-\operatorname{tr}(Q_H^2).                            \tag{7.2}
\]

**Weighted transition statement needed.**  There are universal constants
\(c,C>0\) such that, for an extremal first-moment potential and every long-ray
event \(H\) of quotient probability at least \(c\),

\[
                    \mathcal V_{\rm proj}(H)\le C A.       \tag{7.3}
\]

It would also suffice to have any estimate which forces the left side to
tend to zero with \(A\).  The projector formulation is essential: a packet
with directions \(u=\pm e\) has large signed covariance but is already a
single projective direction.

If (7.3) holds and \(A\) is sufficiently small, then
\(\operatorname{tr}(Q_H^2)\) is close to one.  Since
\(\operatorname{tr}Q_H=1\), its top eigenvalue is close to one, and a
fixed fraction of \(H\) lies in a fixed projective cap.  The cap inverse
then gives \(A\ge c\), a contradiction.  The checkerboard model is
consistent with (7.3), since its right side equals \(m\), while
\(\mathcal V_{\rm proj}=1-1/m\).

The current geometric facts stop one step earlier:

\[
 R^2\mathcal V_{\rm proj}(H)
 \le C\operatorname{tr}\operatorname{Cov}(y\mid H),
 \qquad
 R^2\|D_\Sigma u\|_{HS}^2\le C.                    \tag{7.4}
\]

The first inequality follows from (2.2) and
\(\|uu^T-vv^T\|_{HS}\le\sqrt2\,|u-v|\).  Neither right side is bounded by
\(A\).  Area formula cannot supply the
missing implication because the Gauss map may have rank zero on every
regular high-mass patch, as in Section 6.  Cyclic monotonicity cannot supply
it because it is already saturated edge by edge by the Kantorovich
potential.  A radial bundle is controlled by (3.8), while a simplicial or
polyhedral bundle escapes locally through flat facets and must be charged
through its transition strata.

Thus the multiscale obstruction is now localized precisely: prove a
dimension-free lower bound on the reciprocal-scale mass of the transition
strata between long, nearly flat orientation blocks, using global
log-concavity **and the global extremality of \(f\)**.  Without that new
weighted transition theorem, the spectral series (5.3) remains
non-summable.

## 8. Exact support-contact compatibility

The contact-tensor identity in *tensor_minkowski_completion.md* makes the
support-contact branch load-bearing.  The eikonal endpoint inequalities do
give an exact local theorem in that branch.

### Theorem 8.1 (support plane versus a long segment)

Let \(K\subset\mathbb R^n\) be closed and convex, let \(x\in\partial K\),
and let \(\nu\) be an outward unit support normal at \(x\):

\[
                    \langle \nu,z-x\rangle\le0
                    \qquad(z\in K).                         \tag{8.1}
\]

Suppose \(y+[-R,R]u\subset K\), where \(R>0\) and \(|u|=1\), and put

\[
                         h_x(y)=\langle\nu,x-y\rangle.
\]

Then

\[
                       \boxed{R|\langle\nu,u\rangle|
                              \le h_x(y).}                  \tag{8.2}
\]

If \(u\perp\nu\) and \(x-su\in K\) for some \(s>0\), then

\[
 [x-su,x]\subset K\cap\{z:\langle\nu,z-x\rangle=0\}.       \tag{8.3}
\]

In particular, a tangent segment ending at \(x\) lies in the exposed face
of \(K\).  If that face is the singleton \(\{x\}\), such a segment has
zero length.

#### Proof

Apply (8.1) to \(y+Ru\) and \(y-Ru\):

\[
 -h_x(y)\pm R\langle\nu,u\rangle\le0.
\]

These two inequalities are exactly (8.2).  Under the hypotheses of
(8.3), every point of the segment lies in \(K\), while its scalar product
with \(\nu\) equals that of \(x\).  Hence every point lies in the exposed
face.  \(\square\)

Thus a free-boundary contact point, where the interface normal \(u\) is
tangent to \(\partial K\), cannot itself emit a nontrivial normal ray unless
the support has an actual face in direction \(u\).  More precisely, if
\(K\) is strictly convex, \(y_k\to x\), \(u_k\to u\perp\nu\), and
\(y_k+[-R_k,R_k]u_k\subset K\), then \(R_k\to0\); otherwise a subsequence
with \(R_k\ge R>0\) would converge to a forbidden tangent segment at
\(x\).  This is an exact face statement, not a
second-fundamental-form heuristic.

### Corollary 8.2 (local contact packet)

Let \(f\) be 1-Lipschitz.  Let

\[
 f(y_j+tu_j)=t\qquad(-R\le t\le R),\quad 1\le j\le N,      \tag{8.4}
\]

with all displayed segments contained in \(K\).  Let
\((x_i,\nu_i)\), \(1\le i\le M\), be support point-normal pairs satisfying
(8.1).  Assume

\[
                         |x_i-y_j|\le D
                         \qquad\text{for all }i,j.          \tag{8.5}
\]

Then, for all \(j,k,i\),

\[
 |u_j-u_k|\le {2D\over R},\qquad
 |\langle\nu_i,u_j\rangle|\le {D\over R}.                 \tag{8.6}
\]

Consequently, if \(u_0=u_1\) and
\(M_\Gamma=\sum_i a_i\nu_i\nu_i^T\), where \(a_i\ge0\) and
\(\sum_i a_i=1\), then

\[
 u_0^TM_\Gamma u_0\le {9D^2\over R^2},\qquad
 1-\operatorname{tr}\!\left[
       \left({1\over N}\sum_{j=1}^Nu_ju_j^T\right)^2
                         \right]\le {8D^2\over R^2}.        \tag{8.7}
\]

#### Proof

By (8.5), \(|y_j-y_k|\le2D\), so Theorem 2.1 gives the first
inequality in (8.6).  Theorem 8.1 and
\(h_{x_i}(y_j)\le|x_i-y_j|\) give the second.  Hence

\[
 |\langle\nu_i,u_0\rangle|
 \le|\langle\nu_i,u_j\rangle|+|u_j-u_0|
 \le {3D\over R},
\]

which proves the first part of (8.7).  For independent uniform indices
\(J,K\),

\[
\begin{aligned}
 2\left\{1-\operatorname{tr}\left[
       \left({1\over N}\sum_j u_ju_j^T\right)^2\right]\right\}
 &=\mathbb E\|u_Ju_J^T-u_Ku_K^T\|_{HS}^2\\
 &\le2\mathbb E|u_J-u_K|^2
 \le {8D^2\over R^2}.
\end{aligned}
\]

The displayed bound in (8.7) is a harmless weakening.  \(\square\)

For a contact patch with tensor

\[
 B_\Gamma={1\over P}\int_\Gamma(x-c)\otimes\nu_x\,dS(x), \tag{8.8}
\]

the same calculation gives, whenever one long packet lies within distance
\(D\) of every point in the patch,

\[
 |B_\Gamma u_0|
 \le {3D\over RP}\int_\Gamma|x-c|\,dS(x).                 \tag{8.9}
\]

Hence a spatially localized packet of \(R\)-long rays has one projective
direction and forces the nearby support-contact normal matrix, and the
nearby contact tensor on the right, to have that direction as an
approximate kernel.  In the flat-slice model the conclusion is exact:
\(B=I-u\otimes u\) and \(Bu=0\).

This identifies the remaining aggregation problem more sharply.  A
dimension-sized global contact tensor is compatible with one long-ray
packet because it can be large on \(u^\perp\).  To sustain many diffuse
long-ray directions, its contact mass must split into spatial patches with
different approximate kernels.  The pairwise eikonal inequality forces
patches whose kernels differ by a fixed angle to have zero-level
basepoints separated by order \(R\).  What is still missing is a
dimension-free, bounded-reuse estimate charging that spatial separation to
isotropic volume or to reciprocal ray scale.  This is precisely the
contact-cell completion theorem; neither the Frobenius lower bound on the
global contact tensor nor local support curvature alone supplies the
aggregation.
