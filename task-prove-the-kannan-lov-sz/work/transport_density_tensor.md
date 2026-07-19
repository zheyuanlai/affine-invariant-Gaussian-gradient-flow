# The transport-density tensor and its curvature budget

## 0. Scope and verdict

This note records one exact measure-theoretic tensor identity and one
exact identity under explicit smooth eikonal and endpoint hypotheses.
They do not prove KLS, but they sharpen the structure forced by a
hypothetical bad cut.

Let \(\mu\) be isotropic and log-concave on its affine support
\(F\simeq\mathbb R^d\), let \(E\) have \(\mu(E)=1/2\), and put

\[
 h=1_E-1_{E^c},\qquad \mu_+=2\,1_E\mu,\qquad
 \mu_-=2\,1_{E^c}\mu .                                     \tag{0.1}
\]

Choose an optimal \(W_1\) transport-ray decomposition from \(\mu_-\) to
\(\mu_+\).  If \(u\) is the oriented ray direction and \(\tau\) is half
the Beckmann transport density, then, with divergence understood
intrinsically on \(F\),

\[
 -\operatorname {div}(\rho\tau u)=\rho h,
 \qquad
 D:=\int\tau\,d\mu={1\over2}W_1(\mu_+,\mu_-).              \tag{0.2}
\]

The first identity is completely measure theoretic and remains valid in
the presence of focal and medial sets (subject only to the standard
transport-ray disintegration stated in Section 1):

\[
 \boxed{
 M:=\int\tau\,u\otimes u\,d\mu
   =\int h\,X\otimes u\,d\mu,
 \qquad \|M\|_{HS}\le1.}                                  \tag{0.3}
\]

In particular

\[
 M\succeq0,\qquad \operatorname {tr}M=D,
 \qquad \operatorname {rank}M\ge D^2.                      \tag{0.4}
\]

The second identity is the smooth eikonal--Bochner budget.  On regular
rays, with \(f\) a Kantorovich potential, \(u=\nabla f\),
\(|\nabla f|=1\), and \(L=\Delta-\nabla V\cdot\nabla\),

\[
 -\partial_u Lf
 =\|\nabla^2 f\|_{HS}^2+\nabla^2V(u,u).                     \tag{0.5}
\]

If \(f\) is globally \(C^3\) on the flux region (or admits an exhaustion
with no internal singular-stratum terms), all ray-end products needed in
the integration by parts vanish, the cut has finite weighted perimeter,
and \(n_E^{\rm in}\) denotes its measure-theoretic **inward** unit normal,
then the smooth identity is

\[
 \boxed{
 \int\tau\bigl(\|\nabla^2f\|_{HS}^2+
                       \nabla^2V(u,u)\bigr)d\mu
 =2\int_{\partial^*E}\langle u,n_E^{\rm in}\rangle\,dP_\mu(E)
 \le2\int_{\partial^*E}\langle u,n_E^{\rm in}\rangle_+\,dP_\mu(E)
 \le2\mu^+(E).}                                           \tag{0.6}
\]

Here \(dP_\mu(E)=e^{-V}d\mathcal H^{d-1}\) on the reduced boundary in
the smooth full-dimensional setting, and
\(P_\mu(E)\le\mu^+(E)\) by the relaxation characterization of perimeter
(with equality for the regular cuts used in the smooth calculation).
The positive part in (0.6) is an
upper bound, not an equality.  This distinction is necessary even in one
dimension: an optimal ray can cross \(E\) first from \(E^c\) to \(E\),
then from \(E\) to \(E^c\), and later back into \(E\).

For a convex-body support or a nonsmooth convex potential, smooth
regularizations suggest an additional nonnegative support/potential
curvature charge.  The present note does **not** prove convergence of
that charge for an arbitrary nonsmooth Kantorovich potential.  A cube cut
has zero classical interior curvature and is a diagnostic showing that
some endpoint/support term is indispensable; it is not by itself a proof
of the general regularization assertion.

Consequently, in the smooth regime of (0.6), if a balanced Cheeger cut
has perimeter \(P\) and its
transport cost obeys the standard tube lower bound \(D\ge c/P\), then its
transport-weighted direction law

\[
 d\nu={\tau\over D}\,d\mu
\]

satisfies

\[
 \left\|\int u\otimes u\,d\nu\right\|_{HS}\le {1\over D}
 \le C P,
 \qquad
 \int\bigl(\|\nabla^2f\|_{HS}^2+
                  \nabla^2V(u,u)\bigr)d\nu\le C P^2.
                                                               \tag{0.7}
\]

Thus a small Cheeger constant forces both high effective direction rank
and very small average turning in that regime.  The remaining
implication--even after one has separately justified the
nonsmooth/support limiting budget--that a
globally log-concave ray congruence cannot have both properties unless it
is controlled by a linear or translated-radial witness--is precisely the
unproved global compatibility step.  Thus (0.3) is an exact general
constraint, whereas (0.6)--(0.7) are conditional on the displayed
smoothness/end-flux hypotheses (or on a still-missing rigorous
regularization theorem).  Neither supplies a proof of KLS.

## 1. Raywise construction of the flux

The construction is first given in a form which does not mention a smooth
interface.  Fix the dual sign by requiring
\[
 f(x_+)-f(x_-)=|x_+-x_-|
\]
for almost every optimally coupled pair \(x_-\sim\mu_-\),
\(x_+\sim\mu_+\).  The standard nonbranching decomposition for this
optimal Kantorovich potential gives a quotient space \((Y,\eta)\), intervals
\(I_y=(a_y,b_y)\), basepoints \(z_y\), and unit vectors \(u_y\) such that

\[
 x=z_y+t u_y,\qquad
 d\mu(x)=q_y(t)\,dt\,d\eta(y).                            \tag{1.1}
\]

The conditional densities \(q_y\) are log-concave probability densities
on their intervals.  The orientation is \(u_y=\nabla f\), so \(f\)
increases at unit speed and optimal mass moves in the increasing-\(t\)
direction.  Every optimal pair lies on one ray; hence the quotient
marginals of \(\mu_-\) and \(\mu_+\) agree, and
\[
 \int_{I_y}h(z_y+t u_y)q_y(t)\,dt=0                     \tag{1.2}
\]
for \(\eta\)-almost every active ray.  The existence of a
one-dimensional coupling which moves only to the right also gives the
stochastic-order inequality

\[
 W_y(t):=-\int_{a_y}^t h(z_y+s u_y)q_y(s)\,ds\ge0.        \tag{1.3}
\]

Explicitly, the conditional source and target distribution functions
satisfy
\[
 F_{-,y}(t)-F_{+,y}(t)
 =2\int_{a_y}^t(1_{E^c}-1_E)q_y(s)\,ds
 =2W_y(t),
\]
so their factor \(2\) exactly matches the normalization in (0.1).

Thus \(W_y\) vanishes at both endpoints.  It is not true that an arbitrary
orientation can repair a cumulative flux which changes sign; the
nonnegativity follows from complementary slackness and the chosen dual
orientation.  Multiple sign intervals along one maximal ray are allowed,
but optimality constrains their cumulative function to remain
nonnegative.

The direction field is measurable: on the transported set take
\(u=\nabla f\), which exists \(\mu\)-almost everywhere, and choose any
Borel representative on the null complement.  Since \(\mu_-\) and
\(\mu_+\) are mutually singular, complementary
slackness places \(\mu\)-almost every point on a nontrivial active ray.

Define \(\tau\) on the ray by requiring that the flux per quotient area is
\(W_y\).  In coordinates with the geometric Jacobian absorbed into
\(q_y\), this means

\[
 \tau(z_y+t u_y)q_y(t)=W_y(t).                            \tag{1.4}
\]

For every compactly supported Lipschitz test function \(\phi\), integration
by parts on each interval gives

\[
\begin{aligned}
 \int\tau\,u\cdot\nabla\phi\,d\mu
 &=\int_Y\int_{I_y}W_y(t){d\over dt}
          \phi(z_y+t u_y)\,dt\,d\eta(y)\\
 &=\int_Y\int_{I_y}h(z_y+t u_y)q_y(t)
          \phi(z_y+t u_y)\,dt\,d\eta(y)\\
 &=\int h\phi\,d\mu .
\end{aligned}                                             \tag{1.5}
\]

This is (0.2) in distributions.  Notice the factor \(1/2\):
\(\mu_+-\mu_-=2h\mu\), whereas the displayed flux has divergence
\(h\mu\).  Thus it is half of the unit-demand Beckmann flux.  No
representative of a Minkowski sum and no regularity of the medial set
enters (1.5).

The one-dimensional layer-cake identity also gives

\[
 \int\tau\,d\mu
 =\int_Y\int_{I_y}W_y(t)\,dt\,d\eta(y)
 ={1\over2}W_1(\mu_+,\mu_-)=D.                           \tag{1.6}
\]

Indeed, on each ray the one-dimensional formula
\[
 W_1(\mu_{-,y},\mu_{+,y})
 =\int_{I_y}|F_{-,y}(t)-F_{+,y}(t)|\,dt
 =2\int_{I_y}W_y(t)\,dt
\]
uses (1.3); integration over the common quotient marginal proves the
last equality in (1.6).  This is the Beckmann factor \(1/2\) with the
normalization in (0.1).  In the special signed-distance threshold case,
after placing the
threshold at \(t=0\), the inner integral is exactly
\(\int |t|q_y(t)dt\).

## 2. Exact tensor identity

The tensor identity is most safely proved ray by ray.  Since \(u_y\) is
constant along its line, (1.2)--(1.3) and integration by parts imply

\[
 \int_{I_y}W_y(t)\,dt
 =\int_{I_y}t\,h(z_y+t u_y)q_y(t)\,dt.                  \tag{2.1}
\]

At an infinite endpoint this uses
\(tW_y(t)\to0\), not merely \(W_y(t)\to0\).  Indeed \(W_y(t)\) is bounded
by the corresponding tail of \(q_y\), and a one-dimensional
log-concave probability has a finite first moment (in fact exponential
tails), which gives the required product limit.

\[
 \int_{I_y}h(z_y+t u_y)q_y(t)\,dt=0.                   \tag{2.2}
\]

Therefore

\[
\begin{aligned}
 \int_{I_y}\tau\,u_y\otimes u_y\,q_y(t)dt
 &=\left(\int_{I_y}W_y(t)dt\right)u_y\otimes u_y\\
 &=\int_{I_y}h(z_y+t u_y)
       (z_y+t u_y)\otimes u_y\,q_y(t)dt.
\end{aligned}                                             \tag{2.3}
\]

Integrating in \(y\) proves the equality in (0.3).  Raywise balance also
proves, without differentiating \(u\) across a singular stratum,

\[
                         \int h u\,d\mu=0.              \tag{2.4}
\]

To prove the norm bound, let \(B\) be any real \(d\times d\) matrix.
Using isotropy and \(|u|=|h|=1\) almost everywhere on the transported
mass,

\[
\begin{aligned}
 |\langle B,M\rangle_{HS}|
 &=\left|\int h\,X^TBu\,d\mu\right|\\
 &\le\left(\int|B^TX|^2d\mu\right)^{1/2}
       \left(\int h^2|u|^2d\mu\right)^{1/2}\\
 &=\|B\|_{HS}.
\end{aligned}                                             \tag{2.5}
\]

Duality of the Hilbert--Schmidt norm gives \(\|M\|_{HS}\le1\).  Since
\(M\) is an integral of positive rank-one matrices,
\(M\succeq0\), and (1.6) gives \(\operatorname {tr}M=D\).  If
\(r=\operatorname {rank}M\), Cauchy--Schwarz yields

\[
 D^2=(\operatorname {tr}M)^2\le r\|M\|_{HS}^2\le r,     \tag{2.6}
\]

which is (0.4).

There is no direct affine-covariant version obtained by inserting
\(A^{-1/2}\) around this particular \(M\), because Euclidean unit
directions and \(W_1\) costs change under a nonsimilarity.  The safe
formulation is to first apply the whitening map on \(F\) and rebuild the
Euclidean optimal rays there.  Affine transport of an already chosen
\(W_1\) decomposition does not preserve its cost.

## 3. Smooth eikonal--Bochner identity

Assume in this section that \(d\mu=e^{-V}dx\), with
\(V\in C^2\) convex, and that on the active ray region the Kantorovich
potential \(f\) is \(C^3\), \(|\nabla f|=1\), and
\(u=\nabla f\).  We additionally assume that this region admits a compact
exhaustion for which all internal singular-stratum and ray-end terms in
the calculations below vanish.  Pointwise \(C^3\) regularity merely
almost everywhere does not imply this assumption.  Differentiating the
eikonal equation gives

\[
                         \nabla^2f\,u=0.                 \tag{3.1}
\]

The weighted Bochner identity is

\[
 {1\over2}L|\nabla f|^2
 =\|\nabla^2f\|_{HS}^2+
   \langle\nabla f,\nabla Lf\rangle+
   \nabla^2V(\nabla f,\nabla f).                         \tag{3.2}
\]

The left side vanishes, so (3.2) is (0.5).  Multiplying by \(\tau\) and
using (1.5) with \(\phi=Lf\), after truncating both the rays and \(Lf\),
gives

\[
 -\int\tau\bigl(\|\nabla^2f\|_{HS}^2+
                  \nabla^2V(u,u)\bigr)d\mu
 =\int hLf\,d\mu .                                      \tag{3.3}
\]

For a finite-perimeter set, use the standard BV convention
\[
 D1_E=-n_E^{\rm out}|D1_E|
      =n_E^{\rm in}|D1_E|.
\]
Since \(h=2\,1_E-1\), the weighted integration-by-parts formula gives

\[
 \int hLf\,d\mu
 =-2\int_{\partial^*E}
       \langle u,n_E^{\rm in}\rangle\,dP_\mu(E).         \tag{3.4}
\]

Equations (3.3)--(3.4) give the signed equality in (0.6).  At a crossing
from \(E^c\) to \(E\), the scalar
\(\langle u,n_E^{\rm in}\rangle\) is positive; at a later crossing from
\(E\) back to \(E^c\), it is negative.  Both crossings can occur on the
same active optimal ray.  Nonnegativity of the left side shows only that
the *total signed* boundary flux is nonnegative.  Dropping its negative
part gives the two upper bounds in (0.6), but not another equality.

There are three limit requirements in this calculation.

1. The truncated endpoint terms \(W_yLf\) must tend to zero.  This is
   automatic on a finite regular ray when \(W_y\to0\) and \(Lf\) has an
   integrable trace; it must be checked at a focal singularity.
2. The tensor \(\nabla^2f\) is interpreted only on regular ray charts.
   For a BV direction field a jump of the Gauss map has a signed matrix
   derivative, but there is no automatic finite measure which is the
   limit of the quadratic energy
   \(\tau\|\nabla^2f\|_{HS}^2d\mu\).  A smoothing can make this energy
   diverge or depend on how the eikonal defect is repaired.
3. For a convex support or a nonsmooth convex potential, a chosen convex
   regularization has a nonnegative Hessian charge.  To use it here one
   must also prove convergence of the optimal potentials, transport
   densities, endpoint products, and the weighted Hessian pairing.
   These facts do not follow just from convexity.

Thus (0.6) is proved here only under the displayed smooth no-end-flux
hypotheses.  The cube and ridge models below show what a successful
general extension must retain, but they do not establish such an
extension.  Keeping only the classical almost-everywhere Hessians is
incorrect.

## 4. Consequence for a bad balanced cut

Suppose \(E\) is a balanced near-Cheeger cut in the smooth regime of
Section 3, with weighted perimeter \(P\), and suppose the two-sided tube
estimate has supplied

\[
                         D\ge {c_0\over P}.             \tag{4.1}
\]

No claim about the approximation needed for nonattainment is proved in
this note.  The estimate is not being inferred from (0.3).

Normalize the transport density to a probability:

\[
                         d\nu={\tau\over D}d\mu .        \tag{4.2}
\]

Then (0.3), (0.6), and (4.1) give

\[
 \operatorname {tr}\int u\otimes u\,d\nu=1,
 \qquad
 \left\|\int u\otimes u\,d\nu\right\|_{HS}
     \le {1\over D}\le {P\over c_0}.                   \tag{4.3}
\]

\[
 \int\bigl(\|\nabla^2f\|_{HS}^2+
        \nabla^2V(u,u)\bigr)d\nu
     \le {2P\over D}\le {2P^2\over c_0}.               \tag{4.4}
\]

In particular the effective rank of the direction covariance is at least
\(c_0^2/P^2\).  Simultaneously, the total regular turning and potential
curvature per unit transported distance is \(O(P^2)\).

This pair is stronger than the trace-only covariance statement
\(\int\sigma_y^2u_y\otimes u_y\,d\eta\preceq I\): it controls the
Hilbert--Schmidt norm of the *first-scale transport weighting* directly.
It still does not contradict isotropy.  A matrix with
\(r\asymp P^{-2}\) eigenvalues all of order \(P^2\) satisfies every line
of (4.3).

## 5. Exact model checks

### 5.1 Gaussian halfspace

For \(\mu=\gamma_d\), \(E=\{x_1>0\}\), take \(f=x_1\) and
\(u=e_1\).  Writing \(\varphi\) for the one-dimensional Gaussian density,

\[
 \tau(t)={\Phi(-|t|)\over\varphi(t)},\qquad
 D=\int\tau d\gamma_d=\mathbb E|G|=\sqrt{2/\pi}.         \tag{5.1}
\]

Thus \(M=D e_1\otimes e_1\).  Here
\(\nabla^2f=0\), \(\nabla^2V(u,u)=1\), and

\[
 \int\tau\nabla^2V(u,u)d\mu=D
 =2\varphi(0)=2P.                                        \tag{5.2}
\]

Both identities are exact.

### 5.2 Gaussian multiple-crossing audit

The positive-part equality originally claimed in (0.6) is false even for
a one-dimensional Gaussian.  Let
\[
 a=\Phi^{-1}(1/5),\qquad b=\Phi^{-1}(3/10),\qquad
 c=\Phi^{-1}(3/5),
\]
and
\[
 E=(a,b)\cup(c,\infty).
\]
Then \(\gamma_1(E)=1/10+2/5=1/2\).  For \(f(t)=t\), \(u=1\), the cumulative
flux is
\[
 W(t)=
 \begin{cases}
   \Phi(t),&t<a,\\
   2/5-\Phi(t),&a<t<b,\\
   \Phi(t)-1/5,&b<t<c,\\
   1-\Phi(t),&t>c.
 \end{cases}                                             \tag{5.3}
\]
It is nonnegative, so \(f\) is an optimal dual potential with the
required orientation.  This single maximal ray crosses from \(E^c\) into
\(E\) at \(a\), back into \(E^c\) at \(b\), and into \(E\) again at
\(c\).  Directly,
\[
 D
 =\int h(t)t\varphi(t)\,dt
 =2\{\varphi(a)-\varphi(b)+\varphi(c)\}
 =2\int_{\partial^*E}
       \langle u,n_E^{\rm in}\rangle\,dP_{\gamma_1}(E).
                                                               \tag{5.4}
\]

\[
 2\int_{\partial^*E}
       \langle u,n_E^{\rm in}\rangle_+\,dP_{\gamma_1}(E)
 =2\{\varphi(a)+\varphi(c)\}>D.                          \tag{5.5}
\]
Since \(\nabla^2V=1\), the left side of (0.6) equals \(D\).  Thus the
signed identity and the positive-part upper bound are both sharp in their
logical roles.

### 5.3 Cube halfspace

For the isotropic cube \([ -a,a]^d\), \(a=\sqrt3\), and
\(E=\{x_1>0\}\), one has
\[
 \tau(t)=a-|t|,\qquad D={a\over2},\qquad
 M={a\over2}e_1\otimes e_1,\qquad P={1\over2a}.          \tag{5.6}
\]
The tensor identity is therefore exact.  The classical interior terms
in (0.6) vanish, while the signed interface term is \(2P=1/a>0\).
Consequently (0.6) cannot be extended to hard supports without an
endpoint/support contribution.  For example, symmetric soft-wall
densities proportional to
\(\exp[-k(|t|-a)_+^2]\), followed by an arbitrarily small smoothing,
converge to the uniform marginal; their smooth identity puts
\(2q_k(0)\to1/a\) into the wall Hessian term.  This verifies the support
charge for this explicit approximation, but does not establish a
canonical charge for arbitrary supports and ray fields.

### 5.4 Isotropic radial exponential median sphere

For \(d\ge2\), let \(d\mu\propto e^{-\lambda|x|}dx\), with
\(\lambda=\sqrt{d+1}\), take \(E=\{|x|>r_0\}\) for a median radius
\(r_0\), and orient \(u\) radially outward.
Rotational symmetry gives

\[
                         M={D\over d}I,qquad
                         \|M\|_{HS}={D\over\sqrt d}.     \tag{5.7}
\]

The convex potential satisfies \(\nabla^2V(u,u)=0\), whereas for the
Kantorovich potential \(f(x)=|x|\),

\[
                         \|\nabla^2f\|_{HS}^2={d-1\over r^2}. \tag{5.8}
\]

The raywise integration-by-parts identity is

\[
 \int\tau{d-1\over r^2}d\mu=2P.                         \tag{5.9}
\]

This is the same degree-one equality which appears in the quotient second
variation.  At the singular endpoint \(r=0\), the radial density is
\(O(r^{d-1})\), \(W(r)=O(r^d)\), and
\(Lf=(d-1)/r-\lambda\), so \(W(r)Lf(r)\to0\).  Hence the displayed
integration by parts is valid despite that removable endpoint.  The
curvature budget therefore does not falsely exclude the radial model.

### 5.5 Product exponential coordinate cut

For a product whose first marginal is
\[
 q(t)={\alpha\over2}e^{-\alpha|t|}
\]
and \(E=\{x_1>0\}\), the regular Hessians vanish away from the ridge
\(x_1=0\).  Here \(W(0)=1/2\), \(\tau(0)=1/\alpha\), and the
distributional Hessian is \(D^2V=2\alpha\delta_0\).  Its natural
regularized pairing is
\[
 \tau(0)q(0)\,2\alpha=\alpha=2q(0)=2P.                  \tag{5.10}
\]
Thus this explicit regularization supplies exactly the right side of
(0.6).  Treating a piecewise-linear convex potential as having zero
curvature almost everywhere would invalidate the identity.

## 6. The remaining global lemma and its circular forms

### 6.1 The monotone signed-distance branch

The multiple-crossing correction does not destroy the intended
signed-distance application.  Assume that \(E\) is a smooth stationary
maximizer of

\[
 J_\mu(E)=\int d(x,\partial E)\,d\mu
\]

under the constraint \(\mu(E)=1/2\), and assume that its regular normal
Voronoi cells exhaust the measure.  The first variation in
signed_distance_extremizer.md proves that every cell is bisected.  If
\(t\) is signed normal distance, positive on \(E\), then on each cell

\[
 h(t)=\operatorname {sign}t,\qquad
 W_y(t)=
 \begin{cases}
  \displaystyle\int_{a_y}^tq_y(s)\,ds,&t<0,\\[4pt]
  \displaystyle\int_t^{b_y}q_y(s)\,ds,&t>0.
 \end{cases}                                               \tag{6.1}
\]

Thus \(W_y\ge0\), and

\[
 \int_{I_y}W_y(t)\,dt=\int_{I_y}|t|q_y(t)\,dt.             \tag{6.2}
\]

The resulting normal flux is feasible for the two normalized halves and
is calibrated by the signed-distance potential.  Hence it is an optimal
\(W_1\) flux and

\[
 D=J_\mu(E).                                                \tag{6.3}
\]

At the unique sign crossing, \(u=n_E^{\rm in}\).  Therefore, whenever the
smooth endpoint hypotheses of (0.6) hold, its signed boundary integral is
exactly \(2P_\mu(E)\); the Gaussian multiple-crossing cancellation from
Section 5.2 cannot occur on this branch.

There is also a dimension-free operator bound for the cell basepoints.
Let \(\eta\) be the quotient probability and write

\[
 X=Y+T N,\qquad
 \bar t(Y)=\mathbb E[T\mid Y],\qquad
 \sigma^2(Y)=\operatorname {Var}(T\mid Y).                \tag{6.4}
\]

Since zero is a conditional median, Cantelli's inequality applied on the
appropriate side gives

\[
 |\bar t(Y)|^2\le\sigma^2(Y).                              \tag{6.5}
\]

The conditional barycenter is \(B=Y+\bar tN\).  The law of total covariance
and isotropy give

\[
 \mathbb E[B\otimes B]\preceq I,\qquad
 \mathbb E[\sigma^2N\otimes N]\preceq I.                  \tag{6.6}
\]

Consequently, for every \(a\in\mathbb R^d\),

\[
 \mathbb E\langle a,Y\rangle^2
 \le2\mathbb E\langle a,B\rangle^2
   +2\mathbb E[\bar t^2\langle a,N\rangle^2]
 \le4|a|^2.                                                \tag{6.7}
\]

Thus the basepoint law has a universal covariance-operator bound even if
its trace is of order \(d\).

Finally, suppose the maximizer is within a fixed factor of \(D_1(\mu)\),
and write \(s(Y)=1/q_Y(0)\) for its conditional median-density scale.
Milman's universal implication \(C_P(\mu)\le C D_1(\mu)^2\), applied to the
one-Lipschitz signed distance, gives

\[
 \mathbb E T^2\le C D_1(\mu)^2.                            \tag{6.8}
\]

The sharp cell comparisons in normal_cell_logconcavity.md give
\(\mathbb E_Y|T|\asymp s(Y)\) and
\(\operatorname {Var}_Y T\asymp s(Y)^2\).  Hence

\[
 \mathbb E_\eta s\asymp D_1(\mu),\qquad
 \mathbb E_\eta s^2\le C D_1(\mu)^2.                      \tag{6.9}
\]

Paley--Zygmund and Markov therefore leave a fixed quotient mass of cells
whose scale is comparable to \(D_1(\mu)\).  On this mass the tensor identity
forces effective direction rank of order \(D_1(\mu)^2\), while (6.7) keeps
each fixed projection of the basepoints at universal scale.  This removes
the multiscale escape.  It does not remove the remaining spatial escape:
different direction packets may use different, nearly orthogonal basepoint
coordinates, whose total trace can still be of order \(d\).
A completion from (4.3)--(4.4), after first proving the required
nonsmooth/support passage, would require a theorem of the following kind.

> **Global eikonal compatibility target.**  There is a universal
> \(c>0\) such that no transport-ray field of an isotropic log-concave
> probability can simultaneously satisfy
> \(\|\int u\otimes u\,d\nu\|_{HS}\le\varepsilon\) and total regular plus
> singular turning at most \(c\varepsilon^2\), unless the transported
> cost is controlled by a finite combination of affine and
> translated-radial functions.

Here “singular turning” is only a placeholder until a lower-semicontinuous
quantity compatible with the smooth quadratic energy has actually been
constructed.  The affine and translated-radial conclusions would be harmless: isotropy
controls affine functions, and the translated thin-shell theorem controls
\(|X-z|\) uniformly in \(z\).  But the target is not proved here.

Several tempting substitutes are circular or false.

* A Poincare inequality for \(\nu\) applied to the coordinate functions of
  \(u\) would give the conclusion immediately, but no dimension-free
  Poincare constant for the transport-density law is available.
* Small classical \(\nabla u\) does not control jumps of \(u\) across a
  medial set.  Cube checkerboards and Gaussian fans place all turning on
  transition strata.
* High effective rank alone is compatible with isotropy and with radial
  congruences.
* The support/ridge charge in (0.6) cannot be replaced by an almost-everywhere
  Hessian.

Accordingly, (0.3)--(0.4) are retained as exact general structural
constraints.  Equations (0.6)--(0.7) are retained only as smooth
conditional constraints, with the signed-flux correction above; they are
not a completed proof of KLS.
