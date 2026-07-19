# Continuous Gaussian posterior reduction

## 0. The dimension-free comparison

Let \(\mu\) be a centered isotropic log-concave probability measure on
\(\mathbb R^n\), let \(X\sim\mu\), let \(G\sim N(0,I_n)\) be independent,
and, for a fixed \(t\in(0,1)\), put

\[
 Z_t=\sqrt t\,X+\sqrt{1-t}\,G,
 \qquad \nu_t=\mathcal L(Z_t).
\]

Both laws are isotropic.  Write

\[
 a=\lambda _1(\mu),\qquad b_t=\lambda _1(\nu_t).
\]

Then

\[
 \boxed{
 a\ge \frac{t b_t}{1+(1-t)b_t}.}                    \tag{0.1}
\]

If \(a<t/(1-t)\), the equivalent forward test-function comparison is

\[
 \boxed{
 b_t\le \frac{a}{t-(1-t)a}.}                        \tag{0.2}
\]

At \(t=1/2\), (0.1) is the fixed-unit-noise estimate
\(a\ge b_{1/2}/(2+b_{1/2})\).  The continuous version has two useful
consequences.

1. For every fixed \(t>0\), a universal gap for all \(\nu_t\) transfers
   to every isotropic log-concave input with no dimension loss.
2. The output potential is positive analytic and satisfies
   \[
   0\preceq D^2V_t\preceq\frac1{1-t}I.               \tag{0.3}
   \]
   Hence, for every fixed \(\varepsilon>0\), KLS is quantitatively
   equivalent to proving a universal gap for the analytic isotropic class
   \(0\preceq D^2V\preceq(1+\varepsilon)I\): take
   \(t=\varepsilon/(1+\varepsilon)\).

The last statement is an upper-curvature target, not a Bakry--Emery lower
curvature theorem.  It remains KLS-strength even when \(\varepsilon\) is
an arbitrarily small fixed number.

## 1. Posterior variance and energy

Take first a bounded smooth \(f\in W^{1,2}(\mu)\) such that

\[
 E f(X)=0,\qquad E f(X)^2=1,
 \qquad q=E|\nabla f(X)|^2.
\]

Set

\[
 F_t(z)=E[f(X)\mid Z_t=z].                            \tag{1.1}
\]

Conditionally on \(Z_t=z\), the potential of \(X\) is

\[
 V(x)+\frac{|z-\sqrt t\,x|^2}{2(1-t)}+
 \text{constant}.
\]

It is \(t/(1-t)\)-strongly convex on the intrinsic convex support.
Posterior Brascamp--Lieb therefore gives, with
\(A_z=\operatorname{Cov}(X\mid Z_t=z)\),

\[
 A_z\preceq\frac{1-t}{t}I,                           \tag{1.2}
\]

\[
 \operatorname{Var}(f(X)\mid Z_t=z)
 \le\frac{1-t}{t}
 E(|\nabla f(X)|^2\mid Z_t=z).                       \tag{1.3}
\]

Averaging and using total variance gives

\[
 d_t:=E\operatorname{Var}(f(X)\mid Z_t)
 \le\frac{1-t}{t}q,                                 \tag{1.4}
\]

\[
 \operatorname{Var}_{\nu_t}(F_t)=1-d_t
 \ge1-\frac{1-t}{t}q.                               \tag{1.5}
\]

Differentiation of the posterior exponential family gives

\[
 \nabla F_t(z)=\frac{\sqrt t}{1-t}
 \operatorname{Cov}(X,f(X)\mid Z_t=z).              \tag{1.6}
\]

For every unit vector \(u\), covariance Cauchy--Schwarz and (1.2) imply

\[
 \begin{aligned}
 |u\cdot\nabla F_t(z)|^2
 &\le\frac{t}{(1-t)^2}
 \operatorname{Var}(u\cdot X\mid Z_t=z)
 \operatorname{Var}(f(X)\mid Z_t=z)\\
 &\le\frac1{1-t}
 \operatorname{Var}(f(X)\mid Z_t=z).
 \end{aligned}                                      \tag{1.7}
\]

Taking the supremum over \(u\) introduces no trace or dimension factor.
After averaging, (1.4) and (1.7) give

\[
 \boxed{
 \int|\nabla F_t|^2d\nu_t
 \le\frac{d_t}{1-t}\le\frac qt.}                  \tag{1.8}
\]

If \(q<t/(1-t)\), (1.5), (1.8), and the variational definition of
\(b_t\) yield

\[
 b_t\le
 \frac{q/t}{1-(1-t)q/t}
 =\frac{q}{t-(1-t)q}.                                \tag{1.9}
\]

## 2. Passage to the spectral edge

If \(a<t/(1-t)\), take centered unit form-domain functions with energies
\(q_k\downarrow a\) in (1.9).  This proves (0.2), and rearrangement gives

\[
 b_t[t-(1-t)a]\le a
 \quad\Longrightarrow\quad
 a\ge\frac{t b_t}{1+(1-t)b_t}.                       \tag{2.1}
\]

If \(a\ge t/(1-t)\), (0.1) is automatic.  Indeed isotropy gives
\(b_t\le1\), and therefore

\[
 \frac{t b_t}{1+(1-t)b_t}
 \le\frac{t}{2-t}<\frac{t}{1-t}\le a.
\]

Thus (0.1) holds without spectral attainment and for every \(t\in(0,1)\).

For general form-domain tests, truncate and approximate intrinsically.
Disintegration places the test in the posterior Sobolev space almost
surely; (1.3) applies there, while Gaussian convolution makes \(F_t\)
smooth.  Equations (1.5) and (1.8) pass by \(L^2\) convergence, weak lower
semicontinuity, and closure of the weighted Sobolev form.  Extended convex
potentials and hard convex supports follow from the standard strongly
log-concave Poincare theorem, or by convex smoothing preserving the
constant \(t/(1-t)\).

## 3. Output curvature and affine supports

Let \(V_t=-\log(d\nu_t/dz)\).  Posterior differentiation gives

\[
 D^2V_t(z)
 =\frac1{1-t}I-\frac{t}{(1-t)^2}A_z.                 \tag{3.1}
\]

The posterior covariance bound (1.2) gives the lower Hessian bound in
(0.3), while the elementary covariance bound \(A_z\succeq0\) gives the
upper Hessian bound.  The density is positive analytic.

If \(\mu\) is supported on a proper affine subspace, translate to the
barycenter, identify the supporting space isometrically with
\(\mathbb R^k\), and take \(G\) intrinsically in that space.  All gradients,
covariances, and Hessians above are then intrinsic.  A point mass is the
excluded case \(k=0\).

## 4. Exact scope

The reduction permits the upper-curvature constant to be chosen as close
to one as desired, provided it is fixed before the dimension and the input
are chosen.  It does not show that a law satisfying (0.3) has a universal
gap.  In particular, any perturbative argument whose constants deteriorate
with dimension, or which uses the unproved coercivity of the Witten
one-form operator, does not close the target.
