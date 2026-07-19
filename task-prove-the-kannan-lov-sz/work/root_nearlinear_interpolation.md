# Near-linear first eigenfunctions: an exact mean-gradient gate

## 1. Purpose

This note isolates a strictly narrower analytic statement that would close
the quantitative part of the tensor-extremizer route.  It does not prove
that statement for every log-concave measure.

Let \(\mu=e^{-V}dx\) be centered, isotropic, and log-concave, with the usual
Neumann reading on a convex support.  Let \(f\) be a normalized first
eigenfunction:

\[
 \mathcal A_\mu f=\lambda f,\qquad
 \int f\,d\mu=0,\qquad \int f^2d\mu=1,\qquad
 \int|\nabla f|^2d\mu=\lambda.                         \tag{1.1}
\]

Put

\[
 a=\int Xf(X)d\mu,\qquad \beta=|a|^2,\qquad
 \ell(x)=a\cdot x,\qquad r=f-\ell.                     \tag{1.2}
\]

The tensor factorization theorem makes \(1-\beta\) small once the mixed
tensor stress is close to the quadratic moment space.  The question is
whether that information, together with the eigenfunction equation, already
forces \(\lambda\) to be numerical.

## 2. Exact algebra

Isotropy and the definition of \(a\) give

\[
 \int r\,d\mu=0,\qquad \int Xr\,d\mu=0,\qquad
 \|r\|_2^2=1-\beta.                                    \tag{2.1}
\]

Testing (1.1) against the coordinate functions gives

\[
                  \int\nabla f\,d\mu=\lambda a.         \tag{2.2}
\]

Consequently

\[
 \boxed{\int\nabla r\,d\mu=-(1-\lambda)a.}              \tag{2.3}
\]

Testing against \(\ell\) also yields

\[
 \int\nabla f\cdot a\,d\mu=\lambda\beta,
\]

and therefore

\[
 \boxed{\int|\nabla r|^2d\mu
       =\lambda+\beta(1-2\lambda).}                    \tag{2.4}
\]

Thus, when \(\beta\) is close to one and \(\lambda\) is small, the residual
has small \(L^2\) norm but its gradient is close in mean and in energy to
\(-a\).  This is the precise cancellation that ordinary \(L^2\)
factorization does not see.

The integrated Bochner--Reilly identity gives

\[
\begin{aligned}
 \lambda^2={}&\int\|D^2f\|_{\rm HS}^2d\mu
 +\int D^2V(\nabla f,\nabla f)d\mu\\
 &+\int_{\partial\operatorname {supp}\mu}
       {\rm II}(\nabla_\tau f,\nabla_\tau f)e^{-V}dS/Z.
                                                               \tag{2.5}
\end{aligned}
\]

All terms on the right are nonnegative.  Since \(D^2r=D^2f\),

\[
                       \|D^2r\|_2\le\lambda.            \tag{2.6}
\]

No Poincare inequality has been used in (2.1)--(2.6).

## 3. A sufficient dimension-free interpolation statement

The following statement would close Target B.

> **Mean-gradient interpolation (MGI).** There is a numerical \(C_M\) such
> that for every isotropic log-concave \(\mu\) and every \(C^2\) function
> \(g\) with
> \[
>       \int g\,d\mu=0,\qquad \int Xg\,d\mu=0,
> \]
> one has
> \[
> \left|\int\nabla g\,d\mu\right|
> \le C_M\left(\|g\|_{L^2(\mu)}
>              +\|D^2g\|_{L^2(\mu;{\rm HS})}\right).    \tag{3.1}
> \]

Indeed, applying (3.1) to \(r\) and using (2.3), (2.6) gives

\[
 (1-\lambda)\sqrt\beta
       \le C_M\left(\sqrt{1-\beta}+\lambda\right).       \tag{3.2}
\]

For example, if

\[
 \sqrt{1-\beta}\le {1\over8C_M},\qquad
 \lambda\le\min\left({1\over2},{1\over8C_M}\right),
\]

then the left side of (3.2) is at least
\(\frac12\sqrt{1-(8C_M)^{-2}}>3/8\), whereas the right side is at most
\(1/4\), a contradiction.  Thus tensor factorization with a sufficiently
small numerical error and (3.1) imply a numerical lower bound for
\(\lambda\).

Statement (3.1) is weaker than the full Poincare inequality in one visible
respect: it controls only the finite-dimensional functional
\(g\mapsto\int\nabla g\), and only after \(g\) is orthogonal to constants
and all linear functions.  A slow eigenfunction orthogonal to the linear
space makes both sides of the eigenfunction version of (2.2) zero, so
(3.1) alone does not immediately exclude a general KLS witness.

## 4. Product measures: exact tensor reduction

Suppose

\[
                 \mu=\nu_1\otimes\cdots\otimes\nu_n,
\]

where every \(\nu_i\) is centered, variance one, and log-concave.  For
\(g\in C^2\), define the one-coordinate projection

\[
                 g_i(x_i)=\mathbb E[g(X)\mid X_i=x_i].
\]

If \(g\) is orthogonal to constants and coordinates, then every \(g_i\) is
orthogonal in \(L^2(\nu_i)\) to \(1\) and \(x_i\).  Moreover,

\[
 \mathbb E_\mu\partial_i g=\mathbb E_{\nu_i}g_i',\qquad
 g_i''=\mathbb E[\partial_{ii}g\mid X_i].               \tag{4.1}
\]

The centered one-coordinate subspaces are mutually orthogonal.  Hence

\[
 \sum_i\|g_i\|_{L^2(\nu_i)}^2\le\|g\|_{L^2(\mu)}^2,\qquad
 \sum_i\|g_i''\|_{L^2(\nu_i)}^2
       \le\int\|D^2g\|_{\rm HS}^2d\mu.                  \tag{4.2}
\]

Consequently a universal one-dimensional estimate

\[
 \left|\int h'd\nu\right|^2
 \le C_1\left(\int h^2d\nu+\int(h'')^2d\nu\right),
 \quad h\perp\{1,x\},                                  \tag{4.3}
\]

tensorizes without a dimension loss and proves (3.1) for product
log-concave measures:

\[
\begin{aligned}
 \left|\int\nabla g\,d\mu\right|^2
 &=\sum_i\left|\int g_i'd\nu_i\right|^2\\
 &\le C_1\sum_i(\|g_i\|_2^2+\|g_i''\|_2^2)\\
 &\le C_1(\|g\|_2^2+\|D^2g\|_2^2).                     \tag{4.4}
\end{aligned}
\]

This explains why the cube and products of centered exponentials do not
stress (3.1): exact Hoeffding orthogonality removes the apparent sum over
coordinates.

## 5. Why the known Stein trace estimate does not prove MGI

Let \(\tau\) be a Stein kernel:

\[
             \mathbb E[Xg]=\mathbb E[\tau\nabla g].
\]

For a function in the constraint space of (3.1),

\[
             0=\mathbb E[\tau\nabla g].
\]

Writing \(c=\mathbb E\nabla g\) gives

\[
             c=-\mathbb E[(\tau-I)(\nabla g-c)].         \tag{5.1}
\]

The known Hilbert--Schmidt estimate
\(\mathbb E\|\tau\|_{\rm HS}^2\le Cn\) loses \(\sqrt n\) in (5.1), and
controlling \(\nabla g-c\) by \(D^2g\) through Poincare would be circular.
An operator bound on the Stein kernel would still leave the same
gradient-to-Hessian issue.  Therefore (3.1) is not obtained by merely
restating the trace theorem.

## 6. Remaining audit

The unresolved tasks are:

1. prove the one-dimensional estimate (4.3) uniformly for arbitrary
   log-concave laws, including interval endpoints and nonsmooth potentials;
2. determine whether the product orthogonality in (4.2) has a replacement
   for irreducible log-concave measures;
3. test (3.1) on isotropic simplices, crosspolytopes, and radial exponential
   laws with boundary-layer and degree-one spherical-harmonic functions;
4. construct a smooth isotropic log-concave counterexample if (3.1) is
   false; and
5. if (3.1) is true, prove approximation to affine supports with its
   numerical constant unchanged.

Until these tasks are completed, (3.1) is a precise candidate gate, not a
proved lemma and not a KLS proof.
