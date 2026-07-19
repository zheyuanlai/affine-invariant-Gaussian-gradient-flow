# Transport rays, normal-flow geometry, and the exact missing compatibility

This note treats the mean-centered first-moment constant

\[
  \mathcal D(\mu):=\sup_{\operatorname{Lip}(f)\leq 1}
       \int |f-\mu f|\,d\mu .
\]

The median-centered version differs from this one by at most a factor two.

## 1. Exact cut identity and balanced rays

For every probability measure with a finite first moment,

\[
 \boxed{\quad \mathcal D(\mu)
 =2\sup_{0<\mu(E)<1}\mu(E)\mu(E^c)
     W_1(\mu_E,\mu_{E^c}).\quad}                 \tag{1}
\]

Indeed, for `g=f-mu f`,

\[
 \int |g|\,d\mu=2\sup_E\int_Eg\,d\mu,
\]

and, writing `p=mu(E)`,

\[
 \int_E(f-\mu f)\,d\mu
 =p(1-p)\left(\int f\,d\mu_E-\int f\,d\mu_{E^c}\right).
\]

Now use Kantorovich--Rubinstein duality and interchange the two suprema.

If `mu` has finite first moment, an optimizer `f` for `D(mu)` exists after
fixing `f(0)=0`: a maximizing sequence has a locally uniformly convergent
subsequence, and `|f(x)|<=|x|` gives passage to `L1(mu)`.  Put

\[
 E=\{f>\mu f\},\qquad p=\mu(E).
\]

Then `f` is also an optimal Kantorovich potential between `mu_E` and
`mu_Ec` (with the appropriate orientation).  Let `pi` be an optimal
coupling.  Equality in the dual problem gives

\[
 f(x)-f(y)=|x-y|\qquad\text{for }\pi\text{-a.e. }(x,y).       \tag{2}
\]

The nonbranching transport-ray decomposition of a Euclidean
1-Lipschitz function partitions the interiors of the strained set into
maximal line segments on which `f` is affine with slope one.  Let `Q` be
the quotient map to rays.  Since (2) implies `Q(x)=Q(y)`,

\[
 Q_\#(\mu|_E/p)=Q_\#(\mu|_{E^c}/(1-p)).                     \tag{3}
\]

Disintegrate `mu=int nu_q eta(dq)` along rays and let
`r_q=nu_q(E)`.  Equality (3) reads

\[
 \frac{r_q}{p}=\frac{1-r_q}{1-p}\quad\text{for }\eta\text{-a.e. }q,
\]

so

\[
 \boxed{\nu_q(E)=p\quad\text{on almost every active ray}.} \tag{4}
\]

Thus balance is exact, not merely an averaged property.  This argument
also handles nonsmooth potentials and singular ray endpoints; only the
normal-Jacobian formulas below require a regular chart.

## 2. Smooth eikonal chart and exact conditional density

Assume in this section that `rho=exp(-V)` is positive and `C2`, `V` is
convex, and that on an open transport region `f` is `C2` and
`|grad f|=1`.  Shift so that `mu f=0`, put

\[
 \Sigma=\{f=0\},\qquad N=\nabla f|_\Sigma,
\]

and let `S_y=D_\Sigma N(y)` be the shape operator.  The eikonal equation
implies `nabla^2 f N=0`; consequently the normal trajectories are straight:

\[
 F(y,t)=y+tN(y),\qquad f(F(y,t))=t.                          \tag{5}
\]

Before the first focal time,

\[
 J(y,t)=\det(I+tS_y),                                       \tag{6}
\]

and the ray conditional has density

\[
 q_y(t)=Z_y^{-1}e^{-V(y+tN(y))}\det(I+tS_y),
 \quad
 Z_y=\int_{I_y}e^{-V(y+tN(y))}\det(I+tS_y)\,dt.             \tag{7}
\]

The quotient probability is `eta(dy)=Z_y dH^{n-1}(y)`.  Formula (7)
is one-dimensional log-concave, because

\[
 -\frac{d^2}{dt^2}\log q_y(t)
 =N^T\nabla^2V(y+tN)N+
   \operatorname{tr}\!\left[S_y^2(I+tS_y)^{-2}\right]\geq0. \tag{8}
\]

For the cut `E={f>0}`, (4) becomes

\[
 \int_0^{\sup I_y}q_y(t)\,dt=p\quad\text{for }\eta\text{-a.e. }y. \tag{9}
\]

The weighted exterior perimeter is exactly

\[
 \mu^+(E)=\int_\Sigma e^{-V(y)}\,dH^{n-1}(y)
          =\int q_y(0)\,\eta(dy).                           \tag{10}
\]

These formulas extend chart by chart; the branching endpoints have
zero conditional one-dimensional measure.

## 3. Three rigorous consequences

Write

\[
 m_y=\int t q_y(t)dt,\qquad
 \sigma_y^2=\int(t-m_y)^2q_y(t)dt.
\]

### 3.1 Per-ray perimeter

The isoperimetric profile `I(s)=q(F^{-1}(s))` of a one-dimensional
log-concave density is concave.  Together with the standard
one-dimensional estimate `q(median)>=c/sigma`, this gives

\[
 q_y(0)\geq c\,\frac{\min(p,1-p)}{\sigma_y}.                 \tag{11}
\]

Hence

\[
 \frac{\mu^+(E)}{\min(p,1-p)}
 \geq c\int\frac{1}{\sigma_y}\,\eta(dy).                   \tag{12}
\]

### 3.2 Covariance matrix constraint

Conditionally on `y`, `X=y+TN(y)`, and therefore

\[
 \operatorname{Cov}(X\mid y)=\sigma_y^2N(y)N(y)^T.
\]

The total-covariance identity gives, for isotropic `mu`,

\[
 \boxed{\int\sigma_y^2N(y)N(y)^T\,\eta(dy)\preceq I.}       \tag{13}
\]

Its trace only gives `int sigma_y^2<=n`; it is dimension-free only when
the long rays have a common direction.

### 3.3 Thin shell gives exactly a fourth-moment constraint

We first record an explicit one-dimensional algebraic lemma.

**Lemma.**  If `T` has a log-concave density and variance `sigma^2`, then
for every real `a`,

\[
 \operatorname{Var}(T^2+2aT)\geq \frac{1}{100}\sigma^4.     \tag{14}
\]

**Proof.**  Standardize `T=m+sigma Z`, so `EZ=0`, `EZ^2=1`.
The elementary sharp density estimate for a one-dimensional isotropic
log-concave law is `||density(Z)||_infty<=1`.  For every `b`, set
`P_b(z)=z^2+bz-1`.  With `epsilon=1/8`, completing the square shows

\[
 \operatorname{Leb}\{|P_b|\leq\epsilon\}
 =\frac{4\epsilon}{
 \sqrt{1+b^2/4+\epsilon}+\sqrt{1+b^2/4-\epsilon}}
 <0.251.
\]

Thus `E P_b(Z)^2 > epsilon^2(1-0.251)>1/100`.  Since
`Var(Z^2+bZ)=E P_b(Z)^2`, rescaling proves (14).  QED.

Along a ray,

\[
 |y+tN|^2=|y|^2+2\langle y,N\rangle t+t^2.
\]

Apply (14) conditionally and then the law of total variance:

\[
 \operatorname{Var}_\mu|X|^2
 \geq\frac1{100}\int\sigma_y^4\,\eta(dy).                  \tag{15}
\]

The dimension-free thin-shell theorem in its quadratic form,
`Var|X|^2<=C_TS n`, therefore yields

\[
 \boxed{\int\sigma_y^4\,\eta(dy)\leq100C_{TS}n.}           \tag{16}
\]

Combining (12), Jensen, and (16) only yields

\[
 \frac{\mu^+(E)}{\min(p,1-p)}\geq c n^{-1/4}.              \tag{17}
\]

There is an equally direct first-moment statement.  Balance (9) makes
the quotient marginals on the two sides identical.  On each ray the
optimal transport distance between the normalized halves is the
difference of their means, `d_y`.  Variance decomposition gives

\[
 \sigma_y^2\geq p(1-p)d_y^2.
\]

Consequently the exact identity (1) gives

\[
 \mathcal D(\mu)=2p(1-p)\int d_y\,\eta(dy)
 \leq 2\sqrt{p(1-p)}\int\sigma_y\,\eta(dy)
 \leq (100C_{TS}n)^{1/4}.                                  \tag{18}
\]

Thus isotropy plus thin shell, used only through conditional radial
variance, has a precise `n^(1/4)` ceiling.

## 4. Focal-curvature dichotomy

Suppose additionally `p in [delta,1-delta]`.  A one-dimensional
log-concave density of standard deviation `sigma` satisfies
`||q||_infty<=1/sigma`.  Hence each side of zero in its support has
length at least `delta sigma`.  Positivity of every factor
`1+t kappa_i` in (6) on the full ray interval then gives

\[
 \boxed{\|S_y\|_{op}\leq\frac{1}{\delta\sigma_y}.}          \tag{19}
\]

Therefore long balanced needles are quantitatively flat.  Parallel
long needles are ruled out by (13), but (13) permits long flat bundles
whose directions are dispersed among many orthogonal directions.

There is a stronger Frobenius-curvature version, which is useful for
testing the proposed focal/Jacobian mechanism.

**Lemma (long balanced rays have small total curvature).**  For every
`delta in (0,1/2]` there is an explicit `C_delta<infinity` such that, if
`p in [delta,1-delta]`, then

\[
 \boxed{\|S_y\|_{HS}^2\leq \frac{C_\delta}{\sigma_y^2}.}     \tag{19a}
\]

**Proof.**  Write `q=exp(-W)` after absorbing the normalizing constant.
Let `t_-` and `t_+` be the `delta/2` and `1-delta/2` quantiles.  Since
`0` is between the `delta` and `1-delta` quantiles, Chebyshev's
inequality gives

\[
 |t_-|+|t_+|\leq C_\delta\sigma.
\]

The interval `[t_-,t_+]` has `q`-mass `1-delta`; the bound
`||q||_infty<=1/sigma` therefore gives

\[
 t_+-t_-\geq(1-\delta)\sigma.                               \tag{19b}
\]

Convexity of `W` and the two tail masses `delta/2` imply

\[
 W'_+(t_+)\leq \frac{2}{\delta\sigma},\qquad
 W'_-(t_-)\geq-\frac{2}{\delta\sigma}.                      \tag{19c}
\]

For example, if `W'_+(t_+)>0`, convexity bounds the right tail by
`q(t_+)/W'_+(t_+)`, and `q(t_+)<=1/sigma`; the other inequality is
identical.

The proof of (19) bounds `|kappa_i|<=1/(delta sigma)`.  Together with
the quantile-location bound, this makes every factor
`1+t kappa_i` on `[t_-,t_+]` at most a constant depending only on
`delta`.  Formula (8), interpreted in the sense of convex second
derivative measures if necessary, now gives

\[
 C_\delta^{-1}(t_+-t_-)\sum_i\kappa_i^2
 \leq \int_{t_-}^{t_+}W''(t)dt
 \leq \frac4{\delta\sigma}.
\]

Use (19b) to obtain (19a).  QED.

Thus a spherical family cannot carry long conditionals: curvature in
`n-1` directions enters through the Hilbert--Schmidt norm in (8).
The only remaining escape is large, almost-flat orientation patches
whose changes of direction occur in transition regions of very small
quotient weight.

## 5. Explicit countermodel to all local/averaged constraints

The following labelled mixture saturates (16) and shows exactly what
the preceding information fails to encode.  Let

\[
 s^2=\sqrt n,\quad L=\sqrt3\,n^{1/4},\quad
 a=\frac{n-\sqrt n}{n-1}.
\]

Choose `J` uniformly from `{1,...,n}`.  Conditional on `J=j`, let

\[
 T\sim\operatorname{Unif}[-L,L],\qquad
 Y\sim N(0,a(I-e_je_j^T)),\qquad X=Y+Te_j,
\]

independently.  Then

\[
 \operatorname{Cov}(X)=\frac1n\sum_j
 \left[a(I-e_je_j^T)+\sqrt n\,e_je_j^T\right]=I,            \tag{20}
\]

and

\[
 \operatorname{Var}|X|^2
 =2a^2(n-1)+\operatorname{Var}(T^2)
 =2a^2(n-1)+\frac45n<3n.                                   \tag{21}
\]

Every labelled component is log-concave.  Its exact ray chart is the
flat normal flow

\[
 F_j(y,t)=y+te_j,\qquad S=0,\qquad J=1,
\]

and every conditional needle is balanced, log-concave, and has

\[
 \sigma_y=n^{1/4},\qquad
 q_y(0)=\frac1{2\sqrt3}n^{-1/4}.                             \tag{22}
\]

Moreover

\[
 \int\sigma_y^2N_yN_y^T= n^{-1/2}I\preceq I,qquad
 \int\sigma_y^4=n.                                         \tag{23}
\]

Thus (11), (13), (16), and the strongest possible focal estimate all
hold, yet the averaged reciprocal ray scale is `n^(-1/4)`.

The push-forward density on `R^n` is the equal mixture

\[
 \rho(x)=\frac1n\sum_{j=1}^n
 \frac{e^{-\sum_{k\ne j}x_k^2/(2a)}}
 {2L(2\pi a)^{(n-1)/2}}\,1_{\{|x_j|\leq L\}}.              \tag{24}
\]

It is not log-concave: at `Le_1` and `Le_2` a single summand is
unsuppressed, whereas at their midpoint every summand is suppressed by
at least `exp(-L^2/(8a))`; for all sufficiently large `n` this violates
the midpoint inequality.  Also, after forgetting the label, the flat
ray charts overlap and are not the fibers of one globally defined
eikonal potential.

This identifies the missing compatibility precisely:

> one must use, in a quantitative global way, that all orientation
> bundles are fibers of a single injective normal congruence and that
> their transverse weights are restrictions of one globally
> log-concave density.

Neither covariance, thin shell, one-dimensional log-concavity,
balancedness, the exact normal Jacobian on each chart, nor focal
curvature separately records this compatibility.  A proposed
"widely-varying rays imply Jacobian concentration" lemma must therefore
contain a global transition estimate between orientation charts.  Such
an estimate cannot be replaced by the local bound (19): the
countermodel has `S=0` on every high-mass chart and puts all changes of
orientation into the absent global gluing data.

### 5.1 A stronger Euclidean countermodel with one global potential

The label can be made geometrically visible, so that overlapping charts
are not the issue.  This sharper version isolates global log-concavity
as the missing input.

Assume `n>=2`, use cyclic indices, and set

\[
 R=10\sqrt n,\qquad r=\sqrt{n+1},\qquad
 L=\sqrt3\,n^{1/4}.
\]

For `j in {1,...,n}` and `epsilon in {+1,-1}`, let

\[
 b_{j,\epsilon}=\epsilon R e_{j+1}
\]

and let `K_{j,epsilon}` be the cylinder

\[
 K_{j,\epsilon}=b_{j,\epsilon}
  +\{z\in e_j^\perp:|z|\leq r\}+[-L,L]e_j.                 \tag{26}
\]

The `2n` cylinders are pairwise disjoint.  Indeed, distinct centers are
at distance at least `R sqrt(2)`, while every cylinder lies in the ball
of radius `sqrt(r^2+L^2)` about its center, and

\[
 R\sqrt2-2\sqrt{r^2+L^2}>2L.
\]

Let `mu_0` be the equal mixture of the uniform probabilities on these
cylinders.  On `K_{j,epsilon}`, define

\[
 f_0(b_{j,\epsilon}+z+te_j)=t.
\]

This is 1-Lipschitz on each cylinder, and the displayed separation
shows it is 1-Lipschitz on their union.  Therefore it has a global
1-Lipschitz McShane extension to `R^n`.  On the support of `mu_0`, all
its transport rays are the disjoint flat fibers
`b_{j,epsilon}+z+[-L,L]e_j`.

The mixture is centered.  If `Z` is uniform on the radius-`r` ball in
`e_j^perp`, then `Cov(Z)=I-e_je_j^T`; if `T` is uniform on `[-L,L]`,
then `Var(T)=sqrt(n)`.  Consequently

\[
 \operatorname{Cov}(\mu_0)=\lambda I,
 \qquad
 \lambda=100+\frac{n-1}{n}+\frac1{\sqrt n}.                 \tag{27}
\]

After the scalar rescaling `X'=X/sqrt(lambda)`, the resulting measure
`mu` is exactly isotropic and

\[
 f'(x')=\lambda^{-1/2}f_0(\sqrt\lambda x')
\]

is globally 1-Lipschitz.  Its cut `{f'>0}` has mass one half on every
ray.  Reflection `t -> -t` gives a cross-cut optimal coupling, because
`f'(x)-f'(y)=|x-y|` for every paired pair.  Thus `f'` is an exact
Kantorovich potential for this cut, not merely a chartwise surrogate.

The thin-shell-sized radial variance also holds.  With
`X=b+Z+Te_j`, symmetry and independence give

\[
 \operatorname{Var}|X|^2
 =\operatorname{Var}|Z|^2+\operatorname{Var}(T^2)
   +4\operatorname{Var}\langle b,Z\rangle.
\]

For the `(n-1)`-dimensional ball of radius `sqrt(n+1)`,

\[
 \operatorname{Var}|Z|^2=\frac{4(n-1)}{n+3}<4,
\]

while `Var(T^2)=4n/5` and
`4 Var<b,Z>=4R^2=400n`.  Hence

\[
 \operatorname{Var}_\mu|X'|^2
 <\frac{401n+4}{\lambda^2}<0.05n.                           \tag{28}
\]

Nevertheless every ray has

\[
 \sigma_y=\frac{n^{1/4}}{\sqrt\lambda},\qquad
 q_y(0)=\frac{\sqrt\lambda}{2\sqrt3}\,n^{-1/4},            \tag{29}
\]

and shape operator zero.  Thus Euclidean injectivity and a single
global 1-Lipschitz Kantorovich potential still do not improve the
`n^{-1/4}` reciprocal scale.

The sole failed hypothesis is global log-concavity: `mu_0` is a mixture
of uniform measures on disjoint convex cylinders, and its support is
not convex.  Therefore any successful transition estimate must use
log-concavity *across* orientation bundles; neither the existence of a
global potential nor absence of ray overlap is enough.

## 6. What a closing lemma would have to say

A sufficient genuinely new statement would be, for the *single*
normal congruence arising from an optimizer and a globally log-concave
isotropic density,

\[
 \int \sigma_y\,\eta(dy)\leq C                              \tag{25}

\]

(or the stronger `int sigma_y^2<=C`).  By (18), (25) immediately gives
the desired dimension-free first-moment concentration.  The labelled
mixture proves that (25) cannot be derived from the currently available
raywise and moment constraints; the load-bearing new input must rule
out orientation mixing without paying through a low-weight transition
region.  Establishing that transition estimate is the unresolved
geometric core of this route.
