# Near-linear first eigenfunctions: exact Bochner defects and the missing rigidity lemma

## 0. Verdict

Let \(\mu\) be isotropic and log-concave, let \(A=-L\) be its weighted
Neumann generator, and let

\[
 Af=\lambda f,\qquad \mathbb E f=0,\qquad \mathbb E f^2=1.
\]

Write the orthogonal projection of \(f\) onto the linear functions as

\[
 P_Uf=a\cdot x,\qquad r=f-a\cdot x,\qquad
 \delta=\|r\|_2.
\]

There is a useful exact second-order identity which seems not to have been
recorded in the preceding Target B notes:

\[
 \boxed{
 \left(\mathbb E\|D^2f\|_{\mathrm{HS}}^2
       -\lambda\operatorname {Var}(\nabla f)\right)
 +\mathcal R_V(f)
 =\lambda^3(1-\delta^2).}
 \tag{0.1}
\]

Here \(\mathcal R_V(f)\ge0\) is the curvature-plus-Reilly-boundary term.
The two terms on the left are separately nonnegative.  Thus, if \(\lambda\)
is small and \(\delta<1\), every centered derivative of \(f\) is collectively
an almost extremizer of the same bottom Poincare inequality, and the
curvature/boundary part of Bochner is also small.

This is genuine extra structure, but it does not by itself close Target B.
Without a second spectral gap, the derivative almost-extremizers can occupy
a large bottom spectral cluster.  Curl-freeness does not presently turn that
cluster into a contradiction.  The natural inequality which would do so for
all functions,

\[
 \operatorname {Var}(\nabla g)
 \le C\mathbb E\|D^2g\|_{\mathrm{HS}}^2,
 \tag{0.2}
\]

is quantitatively equivalent to KLS.  A seemingly weaker Landau-type
interpolation,

\[
 \operatorname {Var}(\nabla g)
 \le C\operatorname {dist}_{L^2}(g,\mathrm{Aff})
       \|D^2g\|_2,
 \tag{0.3}
\]

is false for general log-concave functions.  Section 3 gives a smooth,
full-support, radial convex counterexample for which the left side is \(4\)
and the right side is \(O(n^{-1/2})\).

It is tempting to retain the factor \(\delta\) only for first
eigenfunctions:

\[
 \boxed{
 \operatorname {Var}(\nabla f)
 \le C\,\delta\,\|D^2f\|_2
 \quad\text{for a first eigenfunction }f,}
 \tag{EIG-INT}
\]

but this eigenfunction-specific version is also false.  The first
degree-one Neumann mode of the isotropic Euclidean ball has

\[
 \delta\sim n^{-1},\qquad
 \operatorname {Var}(\nabla f)\sim n^{-1},\qquad
 \|D^2f\|_2\sim\sqrt{3/n},
\]

so the ratio in (EIG-INT) grows like \(\sqrt{n/3}\).  This is an explicit
convex-body counterexample consisting of genuine first eigenfunctions, not
arbitrary test functions.

The precise live lemma left by this audit is instead the residual
mean-gradient bound

\[
 \boxed{
 |\mathbb E\nabla r|
 \le C\big(\|r\|_2+\|D^2r\|_2\big).}
 \tag{MG}
\]

It closes Target B without using the normalized KKT equation.  It is proved
below in dimension one, and it tensorizes rigorously to arbitrary products
of isotropic one-dimensional log-concave laws with no dimension loss.  It
also survives the Gaussian, interval/cube, regular-simplex quadratic, radial
exponential, and Euclidean-ball tests.  It remains open for irreducible
log-concave laws.  Even its quadratic restriction contains a nontrivial
Hilbert--Schmidt third-moment tensor estimate, so (MG) should not be treated
as an elementary interpolation theorem.

The net conclusion is therefore:

* Bochner gives the exact defect decomposition (0.1), a stronger starting
  point than the raw estimate \(\|D^2f\|_2\le\lambda\).
* The full curl-free/Korn closure is KLS in disguise.
* The natural multiplicative interpolation fails even after restricting to
  true first eigenfunctions; isotropic balls are the countermodel.
* The narrower residual mean-gradient bound (MG) would close Target B.  It
  holds in one dimension and for product laws, but no irreducible
  high-dimensional proof was found.
* The missing step is a boundary/score or spectral-cluster rigidity statement,
  not another Bochner calculation.

## 1. Exact identities for the affine residual

### 1.1 Setting

The formulas apply to a smooth convex domain \(\Omega\), with

\[
 d\mu=Z^{-1}e^{-V(x)}\mathbf 1_\Omega(x)\,dx,
 \qquad V\text{ convex},
\]

and Neumann boundary condition.  Full space is obtained by taking
\(\Omega=\mathbb R^n\).  The identities involving only the closed form

\[
 \mathcal E(g,h)=\mathbb E\langle\nabla g,\nabla h\rangle
\]

continue to hold for nonsmooth convex bodies.  Smoothness is used only to
write the pointwise Hessian and Reilly terms.

Isotropy and the normalization of \(f\) give

\[
 a=\mathbb E[Xf(X)],\qquad
 \mathbb E r=0,\qquad \mathbb E[Xr(X)]=0,
 \tag{1.1}
\]

and Pythagoras gives

\[
 \boxed{|a|^2=1-\delta^2.}
 \tag{1.2}
\]

Testing the Rayleigh quotient with a linear function shows \(0<\lambda\le1\).

### 1.2 Mean gradients and residual energy

The weak eigenfunction equation tested against \(b\cdot x\) gives, for every
\(b\in\mathbb R^n\),

\[
 b\cdot\mathbb E\nabla f
 =\mathcal E(f,b\cdot x)
 =\lambda\mathbb E[f(b\cdot X)]
 =\lambda b\cdot a.
\]

Consequently

\[
 \boxed{\mathbb E\nabla f=\lambda a,\qquad
        \mathbb E\nabla r=(\lambda-1)a.}
 \tag{1.3}
\]

Since \(\mathbb E|\nabla f|^2=\lambda\), direct expansion gives

\[
\begin{aligned}
 \mathbb E|\nabla r|^2
 &=\lambda-2\lambda|a|^2+|a|^2\\
 &=\boxed{1-\lambda+(2\lambda-1)\delta^2}.
\end{aligned}
\tag{1.4}
\]

The centered gradients of \(f\) and \(r\) agree.  Their exact common
variance is

\[
\begin{aligned}
 \operatorname {Var}(\nabla f)
 &=\mathbb E|\nabla f-\lambda a|^2\\
 &=\boxed{\lambda-\lambda^2|a|^2
 =\lambda\big(1-\lambda(1-\delta^2)\big).}
\end{aligned}
\tag{1.5}
\]

Finally,

\[
 \boxed{D^2r=D^2f.}
 \tag{1.6}
\]

Equations (1.3)--(1.6) display the cancellation a hypothetical small gap
would require.  If \(\lambda\ll1\) and \(\delta\ll1\), then \(r\) is small
in \(L^2\), \(D^2r\) is small in \(L^2\), but

\[
 \mathbb E\nabla r\simeq-a,
 \qquad
 \nabla r\simeq-a\quad\text{in }L^2,
 \tag{1.7}
\]

because \(\nabla r+a=\nabla f\) has squared norm \(\lambda\).  Ruling out
exactly this cancellation is Target B in second-order form.

### 1.3 Bochner--Reilly and the exact derivative-defect split

Put

\[
 B_f=\mathbb E\|D^2f\|_{\mathrm{HS}}^2.
\]

For a smooth convex weighted domain, the integrated Bochner--Reilly identity
is

\[
 \lambda^2
 =B_f+\mathcal R_V(f),
 \tag{1.8}
\]

where

\[
\begin{aligned}
 \mathcal R_V(f)
 ={}&\mathbb E\big[D^2V(\nabla f,\nabla f)\big]\\
 &+Z^{-1}\int_{\partial\Omega}
   \mathrm {II}(\nabla_Tf,\nabla_Tf)e^{-V}\,d\mathcal H^{n-1}
 \ge0.
\end{aligned}
\tag{1.9}
\]

Thus

\[
 \boxed{B_f\le\lambda^2.}
 \tag{1.10}
\]

Apply the same bottom Poincare inequality componentwise to
\(\partial_i f\).  Since Poincare subtracts the component mean, (1.5) gives

\[
 B_f
 =\sum_i\mathbb E|\nabla\partial_i f|^2
 \ge\lambda\sum_i\operatorname {Var}(\partial_i f)
 =\lambda\operatorname {Var}(\nabla f).
 \tag{1.11}
\]

Define the total derivative Poincare deficit

\[
 \mathcal D_\nabla(f)
 :=B_f-\lambda\operatorname {Var}(\nabla f)\ge0.
 \tag{1.12}
\]

Substituting (1.5) and (1.8) gives the exact identity

\[
\begin{aligned}
 \mathcal D_\nabla(f)+\mathcal R_V(f)
 &=\lambda^2-\lambda
   \big(\lambda-\lambda^2|a|^2\big)\\
 &=\boxed{\lambda^3|a|^2
 =\lambda^3(1-\delta^2).}
\end{aligned}
\tag{1.13}
\]

In particular,

\[
 \lambda^2\big(1-\lambda(1-\delta^2)\big)
 \le B_f\le\lambda^2,
 \tag{1.14}
\]

and both the derivative deficit and the curvature/boundary energy are at
most \(\lambda^3(1-\delta^2)\).

This is stronger than the bare Hessian bound.  If \(\lambda\to0\) while
\(\delta\le\delta_0<1\), then

\[
 B_f=(1+O(\lambda))\lambda^2,
 \qquad
 \operatorname {Var}(\nabla f)=(1+O(\lambda))\lambda.
 \tag{1.15}
\]

Thus the centered gradient has size \(\sqrt\lambda\), its derivative has
size \(\lambda\), and it is an almost-bottom spectral field.

### 1.4 Spectral-cluster formulation

Let

\[
 h_i=\partial_i f-\lambda a_i,
 \]

so that every \(h_i\) is centered.  If \(E_A(dt)\) is the spectral
resolution of \(A\) on \(1^\perp\), then

\[
 \mathcal D_\nabla(f)
 =\sum_i\int_{[\lambda,\infty)}(t-\lambda)
       \,d\langle E_A(t)h_i,h_i\rangle.
 \tag{1.16}
\]

Consequently, for every \(\eta>0\),

\[
 \sum_i\left\|
 \mathbf 1_{[(1+\eta)\lambda,\infty)}(A)h_i
 \right\|_2^2
 \le {\mathcal D_\nabla(f)\over\eta\lambda}
 \le {\lambda^2(1-\delta^2)\over\eta}.
 \tag{1.17}
\]

Relative to the total mass (1.5), only \(O(\lambda/\eta)\) of the centered
gradient lies above \((1+\eta)\lambda\).  This is a dimension-free low-mode
statement which does not assume that the first eigenvalue is isolated.

What is missing is a rigidity theorem for a *curl-free vector of many
clustered low modes*.  A bound on the second scalar eigenvalue would turn
(1.17) into proximity to the first eigenspace, but such a bound is not
available and the first eigenspace may have large multiplicity.

### 1.5 Interior PDE identities and the uncontrolled forcing

On full space, with enough regularity,

\[
 A(a\cdot x)=a\cdot\nabla V,
\]

and hence

\[
 \boxed{Ar=\lambda r+\lambda a\cdot x-a\cdot\nabla V.}
 \tag{1.18}
\]

On a bounded Neumann domain this holds in the interior, while

\[
 \partial_\nu r=-a\cdot\nu
 \tag{1.19}
\]

is the nonhomogeneous boundary condition.  Differentiating the eigenvalue
equation on full space gives the one-form equation

\[
 A\nabla f+D^2V\nabla f=\lambda\nabla f.
 \tag{1.20}
\]

After integration,

\[
 \mathbb E[D^2V\nabla f]=\lambda\mathbb E\nabla f=\lambda^2a.
 \tag{1.21}
\]

Equations (1.18)--(1.21) do not bound the score forcing
\(a\cdot\nabla V\), nor the boundary trace in (1.19).  These are precisely
the terms that invalidate a formal unweighted integration-by-parts proof of
a Landau inequality.  For uniform convex bodies, the entire obstruction is
on the boundary; for full-space measures it is carried by the score.

## 2. Which second-order statements are circular?

### 2.1 Affine Korn is KLS-equivalent

Consider the dimension-free estimate

\[
 \operatorname {Var}_\mu(\nabla g)
 \le C\mathbb E\|D^2g\|_{\mathrm{HS}}^2
 \qquad\text{for every smooth }g.
 \tag{2.1}
\]

KLS implies (2.1) immediately by applying scalar Poincare to every
\(\partial_i g\) and summing.  Conversely, apply (2.1) to a normalized first
eigenfunction.  From (1.5), \(|a|\le1\), and (1.10),

\[
 \lambda(1-\lambda)
 \le\operatorname {Var}(\nabla f)
 \le C B_f\le C\lambda^2.
\]

Therefore \(\lambda\ge1/(C+1)\).  The usual spectral-band argument gives
the same implication when the bottom of the spectrum is not attained.
Thus a universal bound in (2.1) is quantitatively equivalent to KLS.

The same warning applies to the Bochner-energy version

\[
 \operatorname {Var}(\nabla g)
 \le C\left(
 \mathbb E\|D^2g\|_{\mathrm{HS}}^2+\mathcal R_V(g)
 \right).
 \tag{2.2}
\]

The optimal constant in (2.2) is within an additive one of
\(C_P(\mu)\).  Calling (2.1) or (2.2) curl-free rigidity does not make it
weaker than the original problem.

### 2.2 A targeted eigenfunction interpolation which also fails

The following statement uses the extra information that \(f\) is both a
first eigenfunction and close to the affine kernel of the Hessian:

\[
 \operatorname {Var}(\nabla f)
 \le C\delta B_f^{1/2}.
 \tag{2.3}
\]

Combining (1.5), (1.10), and (2.3) gives

\[
 \lambda\big(1-\lambda(1-\delta^2)\big)
 \le C\delta\lambda.
\]

After cancelling \(\lambda\),

\[
 \boxed{
 \lambda\ge {1-C\delta\over1-\delta^2}}
 \qquad(C\delta<1).
 \tag{2.4}
\]

Thus (2.3), if true, would close Target B, in fact forcing
\(\lambda\to1\) as \(\delta\to0\).  It is not formally the
KLS-equivalent estimate (2.1): (2.3) is asserted only for a bottom
eigenfunction and gains the small angle \(\delta\).

The defect identity (1.13) initially makes (2.3) look natural.  It says that
the centered gradient already sits almost entirely at the bottom spectral
scale; (2.3) would say that a bottom curl-free field produced by a nearly
affine potential cannot have order-one normalized mass.  Section 3.3 shows
that this inference is false: the first degree-one Neumann mode of an
isotropic Euclidean ball violates (2.3) by a factor asymptotic to
\(\sqrt{n/3}\).  Thin radial concentration makes the eigenfunction
\(O(n^{-1})\)-close to linear while its angular gradient continues to
fluctuate at scale \(n^{-1/2}\).

### 2.3 Mean-gradient coercivity for the affine-orthogonal residual

A different targeted route only tries to control the nonzero mean in (1.3).
For every \(g\perp\mathrm{Aff}\), consider

\[
 |\mathbb E\nabla g|
 \le C\big(\|g\|_2+\|D^2g\|_2\big).
 \tag{2.5}
\]

Applied to \(g=r\), equations (1.3), (1.6), and (1.10) give

\[
 (1-\lambda)\sqrt{1-\delta^2}
 \le C(\delta+\lambda).
 \tag{2.6}
\]

For any \(\delta_0\) satisfying
\(C\delta_0<\sqrt{1-\delta_0^2}\), this yields the numerical lower bound

\[
 \lambda\ge
 {\sqrt{1-\delta_0^2}-C\delta_0
  \over \sqrt{1-\delta_0^2}+C}
 \qquad(\delta\le\delta_0).
 \tag{2.7}
\]

Two stronger variants are

\[
 |\mathbb E\nabla g|^2
 \le C\|g\|_2\|D^2g\|_2,
 \tag{2.8}
\]

and

\[
 |\mathbb E\nabla g|^2
 \le C\|D^2g\|_2^2.
 \tag{2.9}
\]

Either also closes Target B.  For example, (2.8) gives

\[
 (1-\lambda)^2(1-\delta^2)\le C\delta\lambda.
 \tag{2.10}
\]

These mean-gradient estimates are not the same as affine Korn.  KLS does
imply (2.9): if \(m=\mathbb E\nabla g\) and \(g\perp\mathrm{Aff}\), apply
Poincare first to each derivative and then to \(g-m\cdot x\):

\[
\begin{aligned}
 \operatorname {Var}(\nabla g)&\le C_P B_g,\\
 \|g-m\cdot x\|_2^2
 &\le C_P\operatorname {Var}(\nabla g)\le C_P^2B_g.
\end{aligned}
\]

Since \(g\perp m\cdot x\), the last left side equals
\(\|g\|_2^2+|m|^2\).  Hence \(|m|^2\le C_P^2B_g\).

The converse is not immediate: (2.9) only sees functions with a nonzero
mean gradient, and for an eigenfunction it is coercive only through its
linear projection.  It is therefore a plausible intermediate statement
for Target B, not an established reformulation of all of KLS.

### 2.4 A complete proof in one dimension and for product laws

The mean-gradient statement is true, in the stronger form (2.9), for every
one-dimensional isotropic log-concave law.  Let \(\nu\) be centered with
variance one.  The standard one-dimensional log-concave Poincare estimate
may be taken as

\[
 \operatorname {Var}_\nu u\le C_0\int (u')^2d\nu,
 \qquad C_0=12.
 \tag{2.11}
\]

Let \(h\perp\{1,x\}\), and put

\[
 m=\int h'\,d\nu.
\]

Apply (2.11) first to \(h'\):

\[
 \|h'-m\|_2^2\le C_0\|h''\|_2^2.
 \tag{2.12}
\]

Now set \(k=h-mx\).  Since \(\nu\) and \(h\) are centered,
\(\int k\,d\nu=0\).  A second application of (2.11) gives

\[
 \|k\|_2^2
 \le C_0\|h'-m\|_2^2
 \le C_0^2\|h''\|_2^2.
 \tag{2.13}
\]

The affine orthogonality and isotropy give

\[
 \langle k,x\rangle
 =\langle h,x\rangle-m\|x\|_2^2=-m.
 \tag{2.14}
\]

Therefore Cauchy--Schwarz and (2.13) yield the explicit estimate

\[
 \boxed{
 \left|\int h'\,d\nu\right|
 \le12\left(\int(h'')^2d\nu\right)^{1/2}.}
 \tag{2.15}
\]

No boundary condition is imposed on \(h\) or \(h'\).  Hence the proof
covers compact intervals and nonsmooth log-concave densities through their
closed Poincare forms.  General \(W^{2,2}\) functions follow by the usual
one-dimensional approximation.

The estimate tensorizes exactly.  Let

\[
 \mu=\nu_1\otimes\cdots\otimes\nu_n,
\]

where every factor is centered, variance one, and log-concave.  For
\(g\perp\mathrm{Aff}\), define

\[
 g_i(x_i)=\mathbb E[g(X)\mid X_i=x_i].
\]

Then \(g_i\perp\{1,x_i\}\), and Fubini gives

\[
 \mathbb E_\mu\partial_i g=\mathbb E_{\nu_i}g_i',
 \qquad
 g_i''=\mathbb E[\partial_{ii}g\mid X_i].
 \tag{2.16}
\]

Apply (2.15), square, and sum over the coordinates.  Conditional Jensen
then gives

\[
\begin{aligned}
 |\mathbb E_\mu\nabla g|^2
 &=\sum_i|\mathbb E_{\nu_i}g_i'|^2\\
 &\le144\sum_i\|g_i''\|_{L^2(\nu_i)}^2\\
 &\le144\sum_i\mathbb E_\mu(\partial_{ii}g)^2\\
 &\le144\mathbb E_\mu\|D^2g\|_{\mathrm{HS}}^2.
\end{aligned}
\tag{2.17}
\]

Thus product laws satisfy (2.9) with \(C=144\), and (MG) with constant
\(12\), without the \(\|g\|_2\) term.  The argument is a genuine tensor
reduction: it only uses the one-coordinate conditional projections and the
diagonal Hessian entries, so there is no hidden factor of \(n\).

The missing issue is therefore irreducible dependence between coordinates,
not one-dimensional endpoints or tensor bookkeeping.

### 2.5 The quadratic third-moment diagnostic

Even (2.5)--(2.9) contain nontrivial high-dimensional moment information.
For a symmetric matrix \(B\), put

\[
 q_B=X^TBX-\mathbb E[X^TBX],
 \qquad d_B=\mathbb E[Xq_B],
 \qquad g_B=q_B-d_B\cdot X.
 \tag{2.18}
\]

Then \(g_B\perp\mathrm{Aff}\) and

\[
 \mathbb E\nabla g_B=-d_B,
 \qquad D^2g_B=2B,
 \qquad
 \|g_B\|_2^2=\operatorname {Var}(q_B)-|d_B|^2.
 \tag{2.19}
\]

This has an exact operator formulation.  On the Hilbert space
\(\mathrm {Sym}_n\) with the Hilbert--Schmidt inner product, define

\[
 \mathsf T B=\mathbb E[X(X^TBX-\operatorname {Tr}B)]
 \tag{2.20}
\]

and define the covariance operator \(\mathsf C\) by

\[
 \langle B,\mathsf C B\rangle_{\mathrm{HS}}
 =\operatorname {Var}(X^TBX).
 \tag{2.21}
\]

Since \(\mathsf T B\cdot X\) is exactly the linear projection of \(q_B\),

\[
 \mathsf S:=\mathsf C-\mathsf T^*\mathsf T\succeq0,
 \qquad
 \|g_B\|_2^2=\langle B,\mathsf SB\rangle_{\mathrm{HS}}.
 \tag{2.22}
\]

Consequently the restriction of (MG) to quadratic functions is **exactly**
the operator estimate

\[
 \boxed{
 \|\mathsf TB\|_2
 \le C\left(
 \langle B,\mathsf SB\rangle_{\mathrm{HS}}^{1/2}
 +2\|B\|_{\mathrm{HS}}
 \right)
 \quad\text{for every }B\in\mathrm {Sym}_n.}
 \tag{QMG}
\]

Thus the graph norm contains both the part of the quadratic orthogonal to
affine functions and its Euclidean Hessian.  The only dangerous situation
is a quadratic which is very close in \(L^2(\mu)\) to a large linear
function while its coefficient matrix has small Hilbert--Schmidt norm.

Thus (2.9), even only on quadratics, would imply

\[
 |\mathbb E[X(X^TBX)]|\le C\|B\|_{\mathrm{HS}}.
 \tag{2.23}
\]

By Hilbert--Schmidt duality, this is equivalent to

\[
 \sup_{|u|=1}
 \left\|\mathbb E[(u\cdot X)(XX^T-I)]\right\|_{\mathrm{HS}}
 \le C.
 \tag{2.24}
\]

Indeed, Poincare applied to \(q_B\) gives

\[
 \|\mathsf TB\|_2^2
 \le\operatorname {Var}(q_B)
 \le4C_P(\mu)\|B\|_{\mathrm{HS}}^2.
\]

Thus KLS implies (2.24), while any explicit family with
\(\sup_u\|M_u\|_{\mathrm{HS}}\to\infty\) would force
\(C_P(\mu)\to\infty\) and would itself disprove KLS.  A dimension-free proof
of (2.24) from bare log-concavity is therefore a substantive third-moment
tensor theorem.  The regular simplex and product exponential tests give a
constant, not a growing counterexample, but they also show that the quantity
need not be small.

### 2.6 Thin shell gives only a \(\sqrt{\log n}\) tensor bound

Fix a unit vector \(u\) and write

\[
 M_u=\mathbb E[(u\cdot X)(XX^T-I)].
\]

Dimension-free thin shell for every marginal gives, for every rank-\(r\)
orthogonal projection \(P\),

\[
 \operatorname {Var}(|PX|^2)\le C r.
 \tag{2.25}
\]

Therefore covariance Cauchy yields

\[
 |\operatorname {Tr}(PM_u)|
 =|\mathbb E[(u\cdot X)(|PX|^2-r)]|
 \le C\sqrt r.
 \tag{2.26}
\]

Take \(P\) successively onto the positive and negative top eigenspaces of
\(M_u\).  The Ky Fan sums in (2.26) imply that the decreasing absolute
eigenvalues obey

\[
 |\lambda_j(M_u)|\le {C\over\sqrt j}.
 \tag{2.27}
\]

Consequently thin shell proves only

\[
 \|M_u\|_{\mathrm{HS}}^2
 \le C\sum_{j=1}^n{1\over j}
 \le C\log(en).
 \tag{2.28}
\]

This logarithm cannot be removed from the projection inequalities alone:
the abstract matrix

\[
 \operatorname {diag}(1,2^{-1/2},\ldots,n^{-1/2})
\]

satisfies all bounds of the form (2.26) and has squared Hilbert--Schmidt
norm comparable to \(\log n\).  A slicing bound controls normalization and
density height, but supplies no missing compatibility among these
eigendirections.  Hence neither slicing nor thin shell, by itself, proves
(2.24).  One needs a genuinely tensorial property of third moments.

### 2.7 An asymmetric radial-cap test remains bounded

A simple asymmetric body of revolution does not realize the logarithmic
matrix obstruction.  Let \(Y\) be uniform on the half-ball

\[
 K_n=\{y\in\mathbb R^n:|y|\le1,\ y_1\ge0\}.
 \tag{2.29}
\]

Put

\[
 m_n=\mathbb EY_1
 ={\Gamma((n+2)/2)\over
   \sqrt\pi\,\Gamma((n+3)/2)},
 \qquad
 \sigma_{1,n}^2={1\over n+2}-m_n^2,
 \qquad
 \sigma_{\perp,n}^2={1\over n+2}.
 \tag{2.30}
\]

The affine normalization

\[
 X_1={Y_1-m_n\over\sigma_{1,n}},
 \qquad
 X_i={Y_i\over\sigma_{\perp,n}}\quad(i\ge2)
 \tag{2.31}
\]

is exactly centered and isotropic.  Rotational symmetry in the last
\(n-1\) coordinates shows that the only nonzero entries of its symmetric
third-moment tensor are

\[
 \alpha_n=\mathbb EX_1^3,
 \qquad
 \beta_n=\mathbb E[X_1X_i^2]\quad(i\ge2),
 \tag{2.32}
\]

and permutations of the latter.  Conditional on \(Y_1=t\), a transverse
coordinate has second moment \((1-t^2)/(n+1)\).  Hence, exactly,

\[
 \beta_n
 =-{n+2\over n+1}
 {\mathbb E[(Y_1-m_n)Y_1^2]\over\sigma_{1,n}}.
 \tag{2.33}
\]

Let

\[
 \mu_0=\sqrt{2/\pi},
 \qquad \sigma_0=\sqrt{1-2/\pi}.
\]

The beta moments in (2.30)--(2.33) give

\[
 \alpha_n\longrightarrow
 {\mu_0(4/\pi-1)\over\sigma_0^3},
 \qquad
 \beta_n=-{\mu_0\over\sigma_0n}+O(n^{-2}).
 \tag{2.34}
\]

For a unit vector \(u=(u_1,u_\perp)\), the directional third-moment matrix
is therefore

\[
 M_u=
 \begin{pmatrix}
  \alpha_nu_1&\beta_nu_\perp^T\\
  \beta_nu_\perp&\beta_nu_1I_{n-1}
 \end{pmatrix}.
 \tag{2.35}
\]

Consequently

\[
 \|M_u\|_{\mathrm{HS}}^2
 =\alpha_n^2u_1^2+2\beta_n^2|u_\perp|^2
  +(n-1)\beta_n^2u_1^2
 \le C.
 \tag{2.36}
\]

Thus a maximally elementary asymmetric radial cap, with affine
normalization tracked exactly, satisfies the desired directional tensor
bound.  The same is true for products of asymmetric one-dimensional laws:
independence makes \(M_u\) diagonal with entries
\(u_i\mathbb EX_i^3\), and one-dimensional log-concave third moments are
uniformly bounded.  Neither family supplies a counterexample to quadratic
MGI.

### 2.8 A dependent multiscale wedge also remains bounded

One can try to manufacture many third-moment directions by letting a common
coordinate tilt many transverse widths.  Pairing the slopes keeps the
one-dimensional marginal centered.  Namely, for numbers \(a_i\ge0\), consider
the convex body with coordinates \((t,x_1,y_1,\ldots,x_m,y_m)\)

\[
 K_a=\left\{
 |x_i|\le1+a_it,\quad |y_i|\le1-a_it\quad(1\le i\le m)
 \right\},
 \tag{2.36a}
\]

on the interval where all displayed widths are nonnegative.  The \(t\)-slice
volume is proportional to

\[
 p(t)=\prod_{i=1}^m(1-a_i^2t^2),
\]

so the marginal of the original coordinate \(\tau=t\) is even.  Write
\(v=\mathbb E\tau^2\), and set \(T=\tau/\sqrt v\).  After also rescaling
each transverse coordinate to variance one, the two coordinates in the
\(i\)-th pair have opposite directional third moments

\[
 \lambda_i={2a_i\sqrt v\over1+a_i^2v},
 \qquad -\lambda_i.
 \tag{2.36b}
\]

Moreover,

\[
 (-\log p)''(t)
 =\sum_i {2a_i^2(1+a_i^2t^2)\over(1-a_i^2t^2)^2}
 \ge2\sum_i a_i^2.
\]

The one-dimensional Brascamp--Lieb inequality therefore gives

\[
 v\sum_i a_i^2\le {1\over2},
 \qquad
 2\sum_i\lambda_i^2\le 8v\sum_i a_i^2\le4.
 \tag{2.36c}
\]

Sign symmetry leaves no other third-tensor coefficients except permutations
of \(\mathbb E[T X_i^2]=\lambda_i\) and
\(\mathbb E[T Y_i^2]=-\lambda_i\).  Since also
\(|\lambda_i|\le1\), the block matrix for an arbitrary unit direction \(u\)
satisfies

\[
 \|M_u\|_{\mathrm{HS}}^2
 =2u_t^2\sum_i\lambda_i^2
  +2\sum_i\lambda_i^2(u_{x_i}^2+u_{y_i}^2)
 \le4.
 \tag{2.36d}
\]

Thus even this genuinely dependent, many-slope construction obeys the
quadratic MGI tensor bound.  The mechanism is informative: attempting to
accumulate many directional third moments makes the common slice marginal
strongly log-concave, and its variance shrinks by exactly the amount needed
to keep their Hilbert--Schmidt sum bounded.

### 2.9 Why conditioning on the distinguished direction stalls

Fix \(u\), write \(T=u\cdot X\), \(Y=P_{u^\perp}X\), and let

\[
 m(t)=\mathbb E[Y\mid T=t],
 \qquad K(t)=\mathbb E[YY^T\mid T=t].
\]

In the splitting \(\mathbb Ru\oplus u^\perp\), the three blocks of the
directional third-moment matrix are exactly

\[
 \mathbb ET^3,
 \qquad \mathbb E[T^2m(T)],
 \qquad \mathbb E[TK(T)].
 \tag{2.37}
\]

One-dimensional log-concavity bounds the first block.  Isotropy and total
covariance control only

\[
 \mathbb E K(T)+\mathbb E[m(T)m(T)^T]=I_{u^\perp};
 \tag{2.38}
\]

they do not control the Hilbert--Schmidt covariance
\(\mathbb E[TK(T)]\).  Applying the one-dimensional result to the
conditional projection of a quadratic does not repair this.  Differentiating
a conditional expectation creates

\[
 {d\over dt}\mathbb E[g(T,Y)\mid T=t]
 =\mathbb E[\partial_ug\mid T=t]
  -\operatorname {Cov}
   (g,\partial_uV\mid T=t),
 \tag{2.39}
\]

with the analogous moving-boundary term for a convex body.  No
dimension-free estimate for that conditional score covariance follows from
convexity.

Here is the second derivative explicitly.  Put

\[
 s_t(y)=\partial_uV(t,y)-\mathbb E_t\partial_uV.
\]

For any slice function \(F\),

\[
 {d\over dt}\mathbb E_tF
 =\mathbb E_t(\partial_uF-Fs_t).
 \tag{2.40}
\]

Since

\[
 \partial_ts_t
 =\partial_{uu}V-\mathbb E_t\partial_{uu}V
   +\operatorname {Var}_t(\partial_uV),
 \tag{2.41}
\]

differentiating (2.39) once more gives

\[
\boxed{
\begin{aligned}
 h''(t)={}&\mathbb E_t\partial_{uu}g
 -2\mathbb E_t[(\partial_ug)s_t]\\
 &+\mathbb E_t\!\left[g
   \left(s_t^2-\partial_ts_t\right)\right].
\end{aligned}}
 \tag{2.42}
\]

The last multiplier is centered.  Equivalently it is the centered version
of \(s_t^2-\partial_{uu}V\).  If
\(\bar V(t)=-\log\int e^{-V(t,y)}dy\) is the marginal potential, then

\[
 \bar V''(t)
 =\mathbb E_t\partial_{uu}V
  -\operatorname {Var}_t(\partial_uV)\ge0
 \tag{2.43}
\]

by Prekopa.  This only yields

\[
 \operatorname {Var}_t(\partial_uV)
 \le\mathbb E_t\partial_{uu}V.
 \tag{2.44}
\]

The right side of (2.44) has no dimension-free bound from isotropy; it is a
directional Fisher/curvature quantity and can diverge along smooth
approximations of hard convex supports.  On a hard support it is replaced
by the velocity and curvature of the moving slice boundary.  Thus separately
estimating either covariance in (2.39) or the last two terms of (2.42) must
lose an uncontrolled score/boundary norm.  Any successful one-dimensional
reduction has to exploit cancellation among all three terms in (2.42),
together with the global affine orthogonality of \(g\).  Joint convexity
alone supplies no sign for that cancellation.

The same obstruction is visible directly in (2.37): conditional covariance
matrices may rotate and change scale even though their average is fixed by
(2.38).  Thus the most direct one-dimensional localization/conditional
covariance proof stops at the same tensor quantity \(M_u\); it does not
establish (2.24).

## 3. The unrestricted interpolation is false

### 3.1 Uniform isotropic ball

Let \(X\) be uniform on the isotropic ball

\[
 B^n_{\sqrt{n+2}}.
\]

For

\[
 g(x)={|x|^2-n\over\sqrt n},
\]

radial symmetry makes \(g\perp\mathrm{Aff}\).  Direct calculation gives

\[
 \operatorname {Var}(|X|^2)={4n\over n+4},
 \qquad
 \|g\|_2^2={4\over n+4},
 \tag{3.1}
\]

while

\[
 \operatorname {Var}(\nabla g)=4,
 \qquad
 \|D^2g\|_2=2.
 \tag{3.2}
\]

Therefore the right side of (0.3) is \(O(n^{-1/2})\), while the left side
is \(4\).  This already disproves a general dimension-free interpolation.
Here the convex potential is the indicator potential of a ball.

### 3.2 A smooth full-support convex-potential counterexample

The same failure is not a boundary artifact.  In dimension \(n\), let

\[
 V_n(x)=\left({|x|\over s_n}\right)^{2n},
 \qquad
 s_n^2={n\Gamma(1/2)\over\Gamma(1/2+1/n)},
 \tag{3.3}
\]

and let \(d\mu_n\propto e^{-V_n}dx\).  The potential is a smooth convex
polynomial and the choice of \(s_n\) makes \(\mu_n\) isotropic.  Indeed, if
\(R=|X|\), then \((R/s_n)^{2n}\) has the Gamma law with shape \(1/2\), and
\(\mathbb E R^2=n\).

Again take

\[
 g_n(x)={|x|^2-n\over\sqrt n}.
\]

Then

\[
 g_n\perp\mathrm{Aff},\qquad
 \operatorname {Var}(\nabla g_n)=4,
 \qquad \|D^2g_n\|_2=2.
 \tag{3.4}
\]

The radial moment formula gives

\[
 \operatorname {Var}(R^2)
 =n^2\left[
 {\Gamma(1/2)\Gamma(1/2+2/n)
  \over\Gamma(1/2+1/n)^2}-1
 \right].
 \tag{3.5}
\]

Using \(\psi_1(1/2)=\pi^2/2\),

\[
 \|g_n\|_2^2
 ={\operatorname {Var}(R^2)\over n}
 ={\pi^2\over2n}+O(n^{-2}).
 \tag{3.6}
\]

Thus

\[
 {\operatorname {Var}(\nabla g_n)
  \over \operatorname {dist}(g_n,\mathrm{Aff})\|D^2g_n\|_2}
 \sim {2\sqrt{2n}\over\pi}\longrightarrow\infty.
 \tag{3.7}
\]

Adding an arbitrarily small \(\varepsilon_n|x|^2\) and re-isotropizing
produces a smooth strictly convex version with the same divergence.

This example does **not** disprove (2.5)--(2.9), because radial symmetry
gives \(\mathbb E\nabla g_n=0\).  It sharply isolates why a targeted
mean-gradient statement might survive even though full gradient
interpolation fails.

### 3.3 The eigenfunction restriction is false on isotropic balls

The factor \(\delta\) cannot be rescued by restricting to first
eigenfunctions.  Let \(\mu_n\) be uniform on

\[
 B^n_R,\qquad R=\sqrt{n+2}.
 \tag{3.8}
\]

This ball is isotropic.  Its first nonzero Neumann eigenspace is the
degree-one spherical-harmonic sector.  A normalized member has the form

\[
 f_n(x)=c_nh_n(s)\theta_1,
 \qquad s={|x|\over R},\quad \theta={x\over|x|},
 \tag{3.9}
\]

where, up to a scalar multiple,

\[
 h_n(s)=s^{-(n-2)/2}J_{n/2}(z_ns).
 \tag{3.10}
\]

The Neumann condition is

\[
 J_{n/2}(z_n)-z_nJ_{n/2+1}(z_n)=0,
 \tag{3.11}
\]

with \(z_n\) the first positive root, and

\[
 \lambda_n={z_n^2\over n+2}.
 \tag{3.12}
\]

For a radial function \(q\), write

\[
 \langle q\rangle_n=n\int_0^1q(s)s^{n-1}ds,
 \qquad I_0=\langle h_n^2\rangle_n,
 \qquad I_1=\langle sh_n\rangle_n.
\]

Spherical integration gives the exact linear projection and Hessian energy

\[
 |a_n|^2={n+2\over n}{I_1^2\over I_0},
 \tag{3.13}
\]

and

\[
 B_{f_n}
 ={1\over(n+2)^2I_0}
 \left\langle
 (h_n'')^2+3(n-1)
 \left({sh_n'-h_n\over s^2}\right)^2
 \right\rangle_n.
 \tag{3.14}
\]

For completeness, write \(m=(n+2)/2\), \(x=z_n^2\), and

\[
 S_m(x)=\sum_{k\ge0}{(-x/4)^k\over k!(m)_k}.
\]

The Bessel power series turns (3.11) into

\[
 1={x\over2m}{S_{m+1}(x)\over S_m(x)}.
\]

A standard large-order expansion, obtained here by expanding \((m)_k\) in
the absolutely convergent series, gives uniformly for \(x/m\) in a fixed
compact subset of \([0,\infty)\)

\[
 {S_{m+1}(x)\over S_m(x)}
 =1+{x\over4m^2}+O(m^{-2}).
\]

The same expansion first localizes the first positive solution to
\(x=2m+O(1)\), and substitution back into the displayed ratio then gives
\(x=2m-1+O(m^{-1})\).  Thus

\[
 z_n^2=n+1+O(n^{-1}),
 \qquad
 \lambda_n=1-{1\over n}+O(n^{-2}).
 \tag{3.15}
\]

Normalize \(h_n(1)=1\).  The ODE and Neumann condition give

\[
 h_n'(1)=0,\qquad h_n''(1)=n-1-z_n^2=-2+O(n^{-1}),
 \tag{3.16}
\]

so boundary-layer integration in (3.13)--(3.14) yields

\[
 I_0=1-{4\over n^2}+O(n^{-3}),
 \qquad
 I_1=1-{1\over n}-{1\over n^2}+O(n^{-3}).
\]

Indeed the radial law has
\(\langle(1-s)^k\rangle_n=k!/[(n+1)\cdots(n+k)]\), and the radial ODE
together with (3.16) gives
\(h_n(s)=1-(1-s)^2+O((1-s)^3+n(1-s)^4)\) on the contributing boundary
layer.  Substitution gives

\[
 |a_n|^2=1-{1\over n^2}+O(n^{-3}),
 \qquad
 \delta_n={1\over n}+O(n^{-2}),
 \tag{3.17}
\]

and

\[
 B_{f_n}={3\over n}+O(n^{-2}).
 \tag{3.18}
\]

Using the exact identity (1.5),

\[
 \operatorname {Var}(\nabla f_n)
 ={1\over n}+O(n^{-2}).
 \tag{3.19}
\]

Consequently

\[
 \boxed{
 {\operatorname {Var}(\nabla f_n)
  \over\delta_n\sqrt{B_{f_n}}}
 \sim\sqrt{n\over3}\longrightarrow\infty.}
 \tag{3.20}
\]

This disproves (EIG-INT) for genuine normalized first eigenfunctions of
isotropic convex bodies.  Numerically, the ratio in (3.20) is
\(2.945\) at \(n=20\), \(4.314\) at \(n=50\), and \(5.937\) at
\(n=100\).

The mean-gradient gate survives.  For the affine residual,

\[
 |\mathbb E\nabla r_n|
 =(1-\lambda_n)|a_n|\sim n^{-1},
 \]

whereas \(\|r_n\|_2+\|D^2r_n\|_2\sim\sqrt{3/n}\).  The ball therefore
separates the false full-gradient interpolation from the still-live
mean-gradient statement.

## 4. Canonical eigenfunction and moment stress tests

### 4.1 Gaussian

For \(\gamma_n\), a normalized first eigenfunction is

\[
 f(x)=u\cdot x,\qquad \lambda=1.
\]

Hence \(\delta=0\), \(r=0\), \(B_f=0\), and every identity above is exact
with zero defect.  More generally, Gaussian integration by parts gives
\(\mathbb E\nabla g=\mathbb E[Xg]=0\) whenever
\(g\perp\mathrm{Aff}\).  Thus the Gaussian is the equality model for
Target B and places no positive lower bound on the constant in (MG).

### 4.2 Isotropic interval and cube

On \([ -\sqrt3,\sqrt3]\) with normalized Lebesgue measure, put

\[
 f(x)=\sqrt2\sin\left({\pi x\over2\sqrt3}\right),
 \qquad \lambda={\pi^2\over12}.
 \tag{4.1}
\]

Its linear projection coefficient is

\[
 a=\mathbb E[Xf(X)]={4\sqrt6\over\pi^2},
 \tag{4.2}
\]

so

\[
 \delta^2=1-{96\over\pi^4},
 \qquad
 \delta=0.120273\ldots.
 \tag{4.3}
\]

Also

\[
 \|f''\|_2=\lambda,
 \qquad
 \mathbb E f'=\lambda a,
 \qquad
 \mathbb E r'=(\lambda-1)a.
 \tag{4.4}
\]

Numerically,

\[
 \lambda=0.822467\ldots,\qquad
 a=0.992741\ldots,
 \qquad |\mathbb E r'|=0.176244\ldots.
\]

For the eigenfunction interpolation (2.3), the required constant is

\[
 {\operatorname {Var}(f')\over\delta\|f''\|_2}
 =1.574999\ldots.
 \tag{4.5}
\]

Thus the interval is a nontrivial constant-order test, not a degenerate
one.

On the isotropic cube \([ -\sqrt3,\sqrt3]^n\), every normalized combination

\[
 f_u(x)=\sum_{i=1}^n u_i\sqrt2
 \sin\left({\pi x_i\over2\sqrt3}\right),
 \qquad |u|=1,
 \tag{4.6}
\]

is a first Neumann eigenfunction with the same \(\lambda\), \(\delta\), and
ratios as the interval.  Moreover \(B_{f_u}=\lambda^2\); all classical
Reilly curvature vanishes on the flat facets.  Hence tensor multiplicity
does not falsify the targeted inequalities, but it prevents any argument
which assumes a simple first eigenspace.

### 4.3 Regular simplex: a sharp quadratic test of mean-gradient rigidity

Let \(X\) be uniform on an isotropic regular \(n\)-simplex and let \(u\)
point from its center toward a vertex.  The standardized coordinate
\(Z=u\cdot X\) is the centered, variance-one version of a
\(\mathrm {Beta}(1,n)\) variable.  Its third and fourth moments are

\[
 m_3={2(n-1)\sqrt{n+2}\over(n+3)\sqrt n},
 \tag{4.7}
\]

and

\[
 m_4={3(n+2)(3n^2-n+2)\over n(n+3)(n+4)}.
 \tag{4.8}
\]

Put

\[
 g=Z^2-1-m_3Z.
 \tag{4.9}
\]

Simplex symmetry gives \(g\perp\mathrm{Aff}\), and

\[
 \mathbb E\nabla g=-m_3u,
 \qquad \|D^2g\|_2=2,
 \tag{4.10}
\]

while

\[
 \|g\|_2^2=m_4-1-m_3^2
 ={4(n+1)^4\over n(n+3)^2(n+4)}.
 \tag{4.11}
\]

As \(n\to\infty\),

\[
 m_3\to2,\qquad \|g\|_2\to2,
 \tag{4.12}
\]

so both (2.8) and (2.9) require a constant at least one asymptotically.
The simplex therefore survives, and nearly saturates at constant order, the
quadratic mean-gradient tests.

The full third-moment matrix is also harmless but nonzero.  The stabilizer
of \(u\) forces

\[
 M_u:=\mathbb E[(u\cdot X)(XX^T-I)]
 =m_3uu^T-{m_3\over n-1}P_{u^\perp}.
 \tag{4.13}
\]

Thus

\[
 \|M_u\|_{\mathrm{HS}}^2
 =m_3^2\left(1+{1\over n-1}\right)=O(1).
 \tag{4.14}
\]

No growing simplex counterexample to (2.24) results.

These are moment tests, not formulas for the first Neumann eigenfunction of
the simplex.  The latter has no comparably simple closed form.  In
particular, affine geometry alone cannot be used to declare a linear
function a Neumann eigenfunction: its normal derivative is nonzero on the
facets.

### 4.4 Isotropic radial exponential: an exact near-linear first mode

For \(n\ge2\), consider

\[
 d\mu_n(x)=Z_n^{-1}e^{-\alpha|x|}\,dx,
 \qquad \alpha=\sqrt{n+1}.
 \tag{4.15}
\]

The radius has the \(\mathrm {Gamma}(n,\alpha)\) law and
\(\mathbb E|X|^2=n\), so \(\mu_n\) is isotropic.

The ground-state transform sends \(A\) to the Coulomb Hamiltonian

\[
 -\Delta+{\alpha^2\over4}-{\alpha(n-1)\over2|x|}.
 \tag{4.16}
\]

Its first excited shell has eigenvalue

\[
 \boxed{\lambda_n={n\over n+1}.}
 \tag{4.17}
\]

One of its degree-one modes is explicitly

\[
 F_n(x)=c_n x_1e^{t_n|x|},
 \qquad t_n={\alpha\over n+1}={1\over\sqrt{n+1}},
 \tag{4.18}
\]

where \(c_n\) normalizes the \(L^2(\mu_n)\) norm.  Direct substitution gives
\(AF_n=\lambda_nF_n\).

The only linear projection is onto \(x_1\).  Gamma moments give the exact
coefficient

\[
 |a_n|^2
 =\left(1-{1\over n^2}\right)^{n+2},
 \qquad
 \delta_n^2
 =1-\left(1-{1\over n^2}\right)^{n+2}.
 \tag{4.19}
\]

Hence

\[
 \delta_n^2={1\over n}+O(n^{-2}),
 \tag{4.20}
\]

so this is an explicit sequence of genuinely near-linear first
eigenfunctions.  Its gap tends to one, as Target B predicts.

The Hessian can also be integrated exactly:

\[
 B_{F_n}
 ={3n^2-3n+1\over n(n+1)^2}
 ={3\over n}+O(n^{-2}).
 \tag{4.21}
\]

Together with

\[
 \operatorname {Var}(\nabla F_n)
 =\lambda_n-\lambda_n^2|a_n|^2
 ={2\over n}+O(n^{-2}),
 \tag{4.22}
\]

this yields

\[
 {\operatorname {Var}(\nabla F_n)
  \over\delta_n\sqrt{B_{F_n}}}
 \longrightarrow {2\over\sqrt3}.
 \tag{4.23}
\]

Thus the radial exponential gives a constant-order value for the false
eigenfunction interpolation (2.3); unlike the much thinner uniform ball,
it does not make that ratio grow.  Its residual mean-gradient is
\((1-\lambda_n)|a_n|=O(n^{-1})\), so it also survives (MG).

The potential in (4.15) is nonsmooth only at the origin.  Replacing
\(|x|\) by \(\sqrt{|x|^2+\varepsilon^2}\), then re-isotropizing, gives smooth
convex approximants; the displayed eigenpair and ratios are recovered in
the \(\varepsilon\downarrow0\) limit.

## 5. Why the exact defect identity does not yet close

### 5.1 Near equality only gives a spectral cluster

Equation (1.13) says that every row of \(D^2f\) almost saturates Poincare
after the derivative mean is removed.  If there were a dimension-free
separation between the first and second positive eigenvalues, (1.17) would
place every \(h_i\) near the first eigenspace.  No such separation is
available.  Cubes already have first-eigenvalue multiplicity \(n\), and
radial potentials have representation-theoretic multiplicities without
product structure.

Even exact membership of all \(h_i\) in a large first eigenspace would leave
a compatibility problem: one must use

\[
 \partial_jh_i=\partial_ih_j
 \tag{5.1}
\]

to classify curl-free tuples of first modes.  No dimension-free
classification of this kind is known.  A statement that every large bottom
cluster comes from a Gaussian factor is false without additional
hypotheses; irreducible radial measures already have large spherical-harmonic
multiplicities.

There is also a concrete regularity loss in trying to use curl-freeness.
Let

\[
 p_i=\mathbf 1_{[\lambda,(1+\eta)\lambda]}(A)h_i.
\]

Equation (1.17) controls \(h_i-p_i\) only in \(L^2\).  The exact relation
\(\partial_jh_i=\partial_ih_j\) does not pass to \((p_i)\): scalar spectral
projection need not commute with differentiation, and controlling
\(\partial_j(h_i-p_i)\) would require third-derivative information absent
from Bochner.  Applying Bochner again to \(h_i\) is invalid because
\(h_i\) is not a scalar eigenfunction; differentiating the PDE creates the
\(D^2V\nabla f\) commutator and new boundary data.

The standard countermodel to the naive inference “large bottom
multiplicity implies a Gaussian/product factor” is

\[
 V_\varepsilon(x)={|x|^2\over2}
 +{\varepsilon\over4n}|x|^4,
 \qquad \varepsilon>0.
 \tag{5.2}
\]

This smooth strictly convex radial law is irreducible and has no nontrivial
orthogonal product factor, yet its degree-one spherical-harmonic sector has
exact multiplicity \(n\).  As \(\varepsilon\downarrow0\) it can be made
arbitrarily close to the Gaussian first cluster without acquiring an exact
factor.  The uniform-ball calculation in Section 3.3 gives a complementary
countermodel: a genuinely near-linear first mode can retain angular
centered-gradient variance too large for the proposed \(\delta\)-gain.

Neither example has a small spectral gap, so neither is a counterexample to
Target B.  Rather, they formally rule out the two gap-free classification
steps which the derivative-tower argument would need.  Producing an
isotropic log-concave countermodel with \(\lambda\to0\) would itself be
KLS-relevant.

### 5.2 The residual is not an eigenfunction

The small object is \(r=f-a\cdot x\), but (1.18) shows that it solves a
forced equation.  Its forcing contains the directional score
\(a\cdot\nabla V\), and on a convex body its boundary data are
\(-a\cdot\nu\).  Neither is controlled by \(\delta\) and \(B_f\).

Formally applying the unweighted identity

\[
 \|\nabla r\|_2^2=-\langle r,\Delta r\rangle
\]

therefore drops exactly the term which carries the cancellation in (1.3).
On a polytope it is a facet trace; on full space it is the score.  The smooth
radial example in Section 3 shows that log-concavity alone does not restore
a dimension-free general Landau inequality.

### 5.3 Curvature smallness is state-specific

For full-space measures, (1.13) also gives

\[
 \mathbb E[D^2V(\nabla f,\nabla f)]
 \le\lambda^3(1-\delta^2).
 \tag{5.3}
\]

This says that the eigenfunction gradient avoids curvature, not that
\(D^2V\) is small or large in the affine direction \(a\).  In a hypothetical
small-gap regime, \(\nabla r\) cancels \(a\) in both Euclidean and curvature
metrics.  Uniform convex bodies show the limiting obstruction sharply:
\(D^2V=0\) throughout the interior, so all coercion must come from the
boundary/normal cone.  Bochner cannot control that term by itself.

### 5.4 Relation to the normalized KKT equations

Nothing in (1.1)--(1.21) uses KKT stationarity.  Consequently, a proof of
(MG) for all isotropic log-concave laws would remove Target B entirely after
tensor near-linearity has supplied \(\delta\le\delta_0\).  The ball
counterexample shows that this conclusion cannot instead be obtained from
the full centered-gradient interpolation (EIG-INT).

The failure to prove these statements also explains why the KKT equations
may still be essential.  At flat potentials and faceted boundaries, the
missing score/trace coercion is represented by the complementary PSD normal
cone.  A successful KKT argument would have to control that particular
state-selected term; mere positivity and the scalar Bochner budget do not.

## 6. Recommended precise target

The ball calculation rules out the natural full-gradient interpolation,
even on first eigenfunctions.  The sole surviving second-order replacement
for the unresolved directional \(H^{-1}\) estimate is the residual
graph-norm estimate

\[
 |\mathbb E\nabla g|
 \le C(\|g\|_2+\|D^2g\|_2),
 \qquad g\perp\mathrm{Aff},
 \tag{6.1}
\]

(6.1) closes Target B through (2.7).  Section 2.4 proves the stronger
Hessian-only form with \(C=12\) in one dimension and for every product of
isotropic one-dimensional log-concave laws.  Its quadratic restriction
already contains the third-moment map (2.18)--(2.24).  A general proof would
therefore have to use more than generic second-order interpolation; a
counterexample should be sought among
asymmetric thin-shell bodies where a centered quadratic is almost linear.

At present, (1.13), (1.17), the 1D/product theorem, and the explicit model
audits are the positive output.  The irreducible high-dimensional version of
(6.1) remains open.
