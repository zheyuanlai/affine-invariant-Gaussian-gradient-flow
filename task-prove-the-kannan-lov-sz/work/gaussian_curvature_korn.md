# Failure of a post-Gaussian curvature--Korn estimate

## 0. Verdict

The dimension-free slice estimate

\[
 \mathbb E_q\langle D^2U\,F,F\rangle
 \le C\left(\mathbb E_q\|DF\|_{\rm HS}^2
                  +|\mathbb E_qF|^2\right),
 \qquad U=-\log q,
 \tag{0.1}
\]

is false even in one dimension when \(q=\nu*\gamma\), with \(\nu\)
log-concave, and \(F\) is a smooth gradient field.  More strongly, it is
false for the exact centered Poisson field arising from the conditional
slices of one fixed isotropic log-concave law after the required unit
Gaussian convolution.

The counterexample is a long interval convolved with a standard Gaussian.
The Gaussian curvature is of order one in a boundary layer having mass of
order \(a^{-1}\), whereas the slice deformation field changes on scale
\(a\).  Thus the curvature energy is of order \(a^{-1}\) and the Korn
energy is of order \(a^{-2}\).

This does **not** refute the integrated weighted-Fisher target.  In the
exact joint example, the bad slices lie in an exponentially small tail of
the scalar marginal, and the entire weighted curvature term is in fact
finite with a universal bound.  It refutes only a pointwise-in-slice
closure of that target.

## 1. A smooth one-dimensional counterexample

For \(a\ge2\), let

\[
 \nu_a={\rm Unif}[-a,a],\qquad q_a=\nu_a*\gamma,
 \qquad U_a=-\log q_a,
\]

where \(\gamma\) is the standard Gaussian law.  Write \(\phi\) and
\(\Phi\) for the standard Gaussian density and distribution function.
Then

\[
 q_a(x)=\frac{\Phi(x+a)-\Phi(x-a)}{2a}.                 \tag{1.1}
\]

The density \(q_a\) is positive, analytic, even, and log-concave.  Take

\[
 F_a(x)=\frac{x}{a}
       =\frac{d}{dx}\left(\frac{x^2}{2a}\right).       \tag{1.2}
\]

Evenness gives

\[
 \mathbb E_{q_a}F_a=0,
 \qquad \mathbb E_{q_a}|F_a'|^2=\frac1{a^2}.           \tag{1.3}
\]

We next lower-bound the curvature energy.  Put \(x=a+t\), where
\(0\le t\le1/2\), and set

\[
 A=\Phi(2a+t)-\Phi(t),\qquad
 P=\phi(t)-\phi(2a+t).
\]

Differentiating (1.1) gives

\[
 q_a=\frac A{2a},\qquad
 q_a'=-\frac P{2a},\qquad
 q_a''=\frac{t\phi(t)-(2a+t)\phi(2a+t)}{2a}.
\]

Consequently

\[
 U_a''(a+t)q_a(a+t)
 =\frac1{2a}\left{
       \frac{P^2}{A}-t\phi(t)+(2a+t)\phi(2a+t)
   \right}.                                           \tag{1.4}
\]

For \(a\ge2\) and \(0\le t\le1/2\),

\[
 A\le\frac12,\qquad
 P\ge\phi(1/2)-\phi(4),\qquad
 t\phi(t)\le\frac12\phi(1/2).
\]

Hence the expression in braces in (1.4) is at least

\[
 c_*:=2\big(\phi(1/2)-\phi(4)\big)^2
             -\frac12\phi(1/2)>\frac1{16}.             \tag{1.5}
\]

For completeness, the elementary estimates
\(0.352<\phi(1/2)<0.353\) and \(\phi(4)<0.00014\) make the
right side of (1.5) larger than \(0.071>1/16\).  Integrating (1.4) over
this half-unit boundary layer yields

\[
 \int_a^{a+1/2}U_a''q_a\,dx\ge\frac1{64a}.             \tag{1.6}
\]

Since \(F_a^2\ge1\) on the same interval and \(U_a''\ge0\),

\[
 \mathbb E_{q_a}[U_a''F_a^2]\ge\frac1{64a}.            \tag{1.7}
\]

Combining (1.3) and (1.7), the ratio of the left side of (0.1) to its
right-side quadratic form is at least \(a/64\to\infty\).  Thus no
universal constant in (0.1) exists, even for analytic post-Gaussian
densities and smooth gradient fields.

## 2. Realization as an exact Poisson slice field

The preceding field was arbitrary.  The same obstruction occurs for the
field selected by the conditional Poisson equation in the weighted-Fisher
construction.

### 2.1 One fixed isotropic log-concave law

On \(\mathbb R^2\), consider

\[
 p(r,x)=\frac12e^{-r}\mathbf 1_{\{r>0,\ |x|<r\}}.       \tag{2.1}
\]

Its support is a convex cone and its potential is affine there, so this is
a log-concave probability density.  Its first marginal is
\(R\sim\Gamma(2,1)\).  Conditional on \(R=r\), \(X\) is uniform on
\([-r,r]\).  Therefore

\[
 \mathbb ER=2,\quad {\rm Var}(R)=2,\quad
 \mathbb EX=0,\quad
 \mathbb EX^2=\frac{\mathbb ER^2}{3}=2,\quad
 {\rm Cov}(R,X)=0.                                     \tag{2.2}
\]

It follows that

\[
 S=\frac{R-2}{\sqrt2},\qquad Y=\frac X{\sqrt2}          \tag{2.3}
\]

has an isotropic log-concave joint law.  Given \(S=s\),

\[
 Y\mid S=s\sim{\rm Unif}[-a,a],qquad
 a=s+\sqrt2,\qquad s>-\sqrt2.                           \tag{2.4}
\]

After adding an independent \(G\sim N(0,1)\) in the transverse variable,
the conditional law of \(Z=Y+G\) is exactly \(q_a\) from (1.1).

### 2.2 Continuity equation and Poisson field

Let \(u_a(y)=(2a)^{-1}\mathbf 1_{[-a,a]}(y)\).  In distributions,

\[
 \partial_a u_a
 =\partial_y\left(u_a(y)\frac{-y}{a}\right).           \tag{2.5}
\]

Because \(da/ds=1\), convolution of (2.5) with the Gaussian gives

\[
 \partial_sq_s=\partial_z(q_sF_s),\qquad
 F_s(z)=-\frac{M_a(z)}a,qquad
 M_a(z)=\mathbb E[Y\mid Y+G=z].                         \tag{2.6}
\]

In one dimension \(F_s\) is a gradient.  If \(g_s'=F_s\), then

\[
 q_sL_sg_s=(q_sg_s')'=\partial_sq_s,
 \qquad L_s=\partial_{zz}-U_a'\partial_z,
\]

and hence

\[
 L_sg_s=\partial_s\log q_s.                             \tag{2.7}
\]

Thus (2.6) is the exact field of the centered Poisson solution, up to the
irrelevant additive constant in \(g_s\).  It is centered because

\[
 \mathbb E_{q_a}F_s=-\frac{\mathbb EY}{a}=0.            \tag{2.8}
\]

Tweedie's identity and its derivative give

\[
 M_a(z)=z+(\log q_a)'(z),\qquad
 M_a'(z)={\rm Var}(Y\mid Y+G=z).                         \tag{2.9}
\]

The posterior law is a standard Gaussian likelihood restricted to the
interval \([-a,a]\).  The one-dimensional Brascamp--Lieb inequality,
valid by approximation also for the hard endpoints, gives
\(0\le M_a'\le1\).  Consequently the deformation energy satisfies

\[
 C_s:=\mathbb E_{q_a}|F_s'|^2\le\frac1{a^2}.            \tag{2.10}
\]

### 2.3 Boundary-layer curvature energy

Fix \(z=a+t\) with \(0\le t\le1/2\), and write \(W=a-Y\).  Conditional on
\(Z=z\), \(W\) has density proportional to

\[
 e^{-tw-w^2/2}\mathbf 1_{[0,2a]}(w).                   \tag{2.11}
\]

Deleting the part \(w>2a\) from the corresponding law on
\([0,\infty)\) can only lower its mean.  The untruncated mean is decreasing
in \(t\), since its derivative is minus its variance, and at \(t=0\) it
equals \(\sqrt{2/\pi}<1\).  Therefore

\[
 \mathbb E[W\mid Z=a+t]<1,qquad
 M_a(a+t)>a-1\ge\frac a2.                               \tag{2.12}
\]

It follows that \(|F_s|\ge1/2\) throughout the boundary layer.  Using
(1.6),

\[
 B_s:=\mathbb E_{q_a}[U_a''F_s^2]
 \ge\frac14\int_a^{a+1/2}U_a''q_a\,dz
 \ge\frac1{256a}.                                      \tag{2.13}
\]

Equations (2.8), (2.10), and (2.13) imply

\[
 \frac{B_s}{C_s+|\mathbb E_{q_a}F_s|^2}
 \ge\frac a{256}\longrightarrow\infty.                \tag{2.14}
\]

This is a failure of (0.1) for the exact Poisson field on conditional
slices of a fixed isotropic log-concave distribution, not merely for a
freely chosen test field.

## 3. Mean/residual decomposition

For a general conditional family, put

\[
 \bar F_s=\mathbb E_sF_s=-m'(s),\qquad
 F_s^0=F_s-\bar F_s,\qquad \mathbb E_sF_s^0=0,
\]

where \(m(s)=\mathbb E_sZ\), and let \(H_s=D^2_{zz}U\).  Then exactly

\[
 \begin{aligned}
 B_s={}&\bar F_s^T(\mathbb E_sH_s)\bar F_s
       +2\bar F_s^T\mathbb E_s(H_sF_s^0)
       +\mathbb E_s\langle H_sF_s^0,F_s^0\rangle.
 \end{aligned}                                         \tag{3.1}
\]

For a unit Gaussian convolution, \(0\preceq H_s\preceq I\).  Hence, with

\[
 B_s^0=\mathbb E_s\langle H_sF_s^0,F_s^0\rangle,
\]

the elementary inequality \(|u+v|^2\le2|u|^2+2|v|^2\) in the
\(H_s\)-seminorm gives

\[
 B_s\le2|m'(s)|^2+2B_s^0.                              \tag{3.2}
\]

Pure translations account only for the mean sector.  The cone example is
more decisive: every conditional slice in (2.4) is centered, so

\[
 m(s)=0,\qquad \bar F_s=0,
 \qquad F_s^0=F_s.                                      \tag{3.3}
\]

Thus its divergence is entirely in the residual sector:

\[
 B_s^0\ge\frac1{256a},\qquad
 \mathbb E_s|D_zF_s^0|^2\le\frac1{a^2}.                \tag{3.4}
\]

Subtracting the conditional mean therefore cannot repair a local
curvature--Korn estimate.  Any successful weighted-Fisher proof must use
relations between different \(s\)-slices.  In the present identity ledger,
the additional nonnegative datum is the material-curvature charge

\[
 \mathcal Q_s=D^2U[(1,-F_s),(1,-F_s)],
 \qquad
 \Phi''(s)=\mathbb E_s\mathcal Q_s
              +\mathbb E_s\|D_zF_s\|_{\rm HS}^2,       \tag{3.5}
\]

together with its Stein-weighted integral and the centroid identities.
A bound using only \(q_s,F_s\) at one fixed slice cannot close the route.

## 4. Why the integrated WFI target survives this example

The scalar marginal in the cone example is, in terms of
\(a=s+\sqrt2\),

\[
 \rho(s)=2a e^{-\sqrt2a}\mathbf 1_{\{a>0\}}.           \tag{4.1}
\]

The canonical Stein kernel of the standardized gamma marginal is

\[
 \tau(s)=\frac R2=\frac a{\sqrt2}.                      \tag{4.2}
\]

Indeed, the gamma integration-by-parts identity
\(\mathbb E[(R-2)h(R)]=\mathbb E[R h'(R)]\), followed by the scaling
\(S=(R-2)/\sqrt2\), proves (4.2).  Moreover, \(|M_a|\le a\), so
\(|F_s|\le1\); and the Gaussian channel gives \(0\le U_a''\le1\).
Therefore \(B_s\le1\), and

\[
 \int\rho(s)\tau(s)^2B_s\,ds
 \le\mathbb E\tau(S)^2
 =\frac{\mathbb ER^2}{4}=\frac32.                      \tag{4.3}
\]

Thus the same example that makes the pointwise ratio (2.14) arbitrarily
large obeys a universal integrated bound.  Its role is to force any proof
to retain the joint \(s\)-geometry and the weighted material-curvature
budget, rather than replacing them by a slice-by-slice Korn inequality.
