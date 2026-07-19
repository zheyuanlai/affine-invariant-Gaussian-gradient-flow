# Posterior reverse smoothing and a valid fixed-Gaussian reduction

## 0. The theorem

Let \(\mu\) be a centered isotropic log-concave probability measure on
\(\mathbb R^n\), let \(X\sim\mu\), let \(G\sim N(0,I_n)\) be independent,
and put

\[
 S=\frac{X+G}{\sqrt2},\qquad \nu=\mathcal L(S).
\]

Write

\[
 a=\lambda_1(\mu),\qquad b=\lambda_1(\nu).
\]

Then

\[
 \boxed{
 b\le \frac{2a}{1-a}\quad (a<1),
 \qquad
 a\ge\frac{b}{2+b}.}                                  \tag{0.1}
\]

The second inequality also holds when \(a=1\).  In particular, if every
fixed unit-Gaussian regularization has gap at least \(c_0>0\), then every
isotropic log-concave input has gap at least

\[
 \boxed{\lambda_1(\mu)\ge\frac{c_0}{2+c_0}.}          \tag{0.2}
\]

This is the valid fixed-Gaussian reverse reduction.  It is logically
independent of the forward two-law inequality

\[
 b(1-b)\le4(b-a),
\]

whose correct rearrangement is \(a\le b(3+b)/4\) and which by itself has
the wrong direction for a reverse reduction.

## 1. Posterior inequalities

First take \(f\) bounded, smooth, centered, and normalized in
\(L^2(\mu)\).  Write

\[
 q=\int |\nabla f|^2\,d\mu,
 \qquad Y=X+G,
 \qquad F(y)=E[f(X)\mid Y=y].                           \tag{1.1}
\]

Let \(\mu_y=\mathcal L(X\mid Y=y)\), and let

\[
 A_y=\operatorname{Cov}_{\mu_y}(X).
\]

The posterior potential is

\[
 V_y(x)=V(x)+\frac12|x-y|^2+\text{constant}
\]

on the convex support of \(\mu\).  Thus the posterior is
\(1\)-strongly log-concave.  The Brascamp--Lieb inequality, first for a
smooth finite convex potential and then by monotone convex approximation,
gives

\[
 A_y\preceq I,                                         \tag{1.2}
\]

and

\[
 \operatorname{Var}_{\mu_y}(f)
 \le E_{\mu_y}|\nabla f|^2.                            \tag{1.3}
\]

Averaging (1.3) over \(Y\) yields

\[
 E\operatorname{Var}(f(X)\mid Y)\le q.                \tag{1.4}
\]

The law of total variance therefore gives

\[
 \boxed{\operatorname{Var}(F(Y))
 =1-E\operatorname{Var}(f(X)\mid Y)
 \ge1-q.}                                             \tag{1.5}
\]

Differentiating the posterior exponential family gives the exact identity

\[
 \nabla F(y)=\operatorname{Cov}_{\mu_y}(X,f).          \tag{1.6}
\]

For any unit vector \(u\), covariance Cauchy--Schwarz and (1.2) imply

\[
 |u\cdot\nabla F(y)|^2
 \le \operatorname{Var}_{\mu_y}(u\cdot X)
       \operatorname{Var}_{\mu_y}(f)
 \le \operatorname{Var}_{\mu_y}(f).
\]

Taking the supremum over unit \(u\) is legitimate because the right side
does not depend on \(u\), and gives the dimension-free pointwise bound

\[
 |\nabla F(y)|^2\le\operatorname{Var}_{\mu_y}(f).      \tag{1.7}
\]

Combining (1.3), (1.7), and disintegration gives

\[
 \boxed{
 \int|\nabla F|^2d\mathcal L(Y)
 \le E\operatorname{Var}(f(X)\mid Y)
 \le q.}                                              \tag{1.8}
\]

Notice that (1.8) is a genuine posterior \(H^1\) contraction.  It is not
the reverse conditional projection under self-convolution which fails on
the interval model: posterior strong log-concavity from the Gaussian
likelihood is the additional ingredient.

## 2. Spectral-gap comparison

Define the test function on the isotropic output by

\[
 \widetilde F(s)=F(\sqrt2s).
\]

It is centered under \(\nu\), has the same variance as \(F(Y)\), and by
(1.8) satisfies

\[
 \int|\nabla\widetilde F|^2d\nu
 =2\int|\nabla F|^2d\mathcal L(Y)
 \le2q.                                                \tag{2.1}
\]

When \(q<1\), (1.5) and the variational definition of \(b\) give

\[
 b\le\frac{2q}{1-q}.                                  \tag{2.2}
\]

The spectral edge need not be attained.  Since isotropy gives \(a\le1\),
if \(a<1\) choose centered normalized form-domain functions \(f_k\) with
energies \(q_k\downarrow a\), apply (2.2), and let \(k\to\infty\).  This
proves the first inequality in (0.1).  Rearrangement gives

\[
 b(1-a)\le2a
 \quad\Longleftrightarrow\quad
 a\ge\frac{b}{2+b}.
\]

If \(a=1\), the latter conclusion is automatic because \(b\le1\) by
isotropy.  This proves (0.1) in all cases.

## 3. Form-domain and support audit

For a positive smooth log-concave input, (1.6) follows by differentiating
the ratio defining the posterior expectation.  The estimates above then
extend from bounded smooth tests to \(W^{1,2}(\mu)\) by truncation and
closure: (1.5) is stable under \(L^2\) convergence, while (1.8) supplies a
uniform output-form bound and weak lower semicontinuity.

For an extended-valued convex potential or hard convex support, approximate
the posterior potential by smooth finite convex potentials with the same
quadratic lower-curvature constant, apply Brascamp--Lieb, and pass to the
limit.  Equivalently, (1.2)--(1.3) are the standard Poincare and covariance
consequences of \(1\)-strong log-concavity, valid on a convex affine
support with its intrinsic gradient.  Log-concave weighted Sobolev forms
have bounded smooth functions as a core after this intrinsic
approximation.

If the original measure is supported on a proper affine subspace \(E\),
translate to its barycenter, identify \(E\) isometrically with
\(\mathbb R^k\), and add the Gaussian intrinsically in \(E\).  Every
identity above then holds in dimension \(k\).  A point mass has \(k=0\)
and is excluded.  Thus no ambient null directions enter the posterior
covariance argument.

## 4. Curvature of the regularized output

Let \(U=-\log(d\mathcal L(Y)/dy)\).  The Gaussian posterior formula gives

\[
 D^2U(y)=I-A_y.
\]

Together with (1.2) this yields

\[
 0\preceq D^2U\preceq I.
\]

For the potential \(\widetilde U(s)=U(\sqrt2s)+\text{constant}\) of
\(S\),

\[
 \boxed{0\preceq D^2\widetilde U\preceq2I.}           \tag{4.1}
\]

Gaussian convolution makes the output density positive and analytic.
Consequently KLS is quantitatively equivalent to proving a universal gap
for the isotropic analytic log-concave subclass (4.1): one direction is
restriction to a subclass, and the reverse direction is (0.2).

## 5. Canonical checks

For \(\mu=N(0,I)\), the input and normalized output gaps are both one;
(0.1) reads \(1\ge1/3\), with slack because the posterior estimate was
designed for arbitrary slow modes rather than Gaussian linear equality.

For products, the posterior factorizes but the proof never sums posterior
coordinate variances: (1.7) takes an operator norm through a unit-vector
test.  Hence no factor of \(n\) appears.

For a hard-support body such as the cube, simplex, or crosspolytope, the
posterior quadratic likelihood gives the same unit strong-convexity
constant on the convex support.  For the one-sided exponential and its
products, truncation justifies (1.3), and the averaged input energy in
(1.4) remains finite.  Non-symmetry therefore causes no change in the
constant chain.
