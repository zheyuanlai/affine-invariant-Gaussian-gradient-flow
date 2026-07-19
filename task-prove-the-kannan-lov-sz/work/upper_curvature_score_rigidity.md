# Upper-curvature score rigidity and the near-affine spectral branch

## 0. Scope and conclusion

Let \(\mu(dx)=e^{-V(x)}dx\) be a centered isotropic probability measure
on \(\mathbb R^n\).  Assume first that \(V\in C^2\), that the integrations
by parts below have no boundary term, and that

\[
 0\preceq D^2V\preceq \beta I                                      \tag{0.1}
\]

for some \(\beta\ge1\).  These hypotheses hold for every fixed
unit-Gaussian regularization from
\(\mathcal G\mu=\mathcal L((X+G)/\sqrt2)\), with \(\beta=2\).

The score identity below proves a complete dimension-free estimate for
the near-affine branch of the bottom spectrum.  It does not control a
bottom spectral sequence whose affine residual stays bounded away from
zero, and therefore does not by itself prove KLS.

## 1. Score covariance ledger

Put \(Y=\nabla V(X)\), where \(X\sim\mu\), and set

\[
 J=E[YY^T].
\]

Integration by parts gives, for all \(i,j\),

\[
 E Y_i=0,\qquad E[X_iY_j]=\delta_{ij},\qquad
 E[Y_iY_j]=E[\partial_{ij}V(X)].                       \tag{1.1}
\]

Consequently

\[
 J=E D^2V(X),qquad I\preceq J\preceq\beta I,          \tag{1.2}
\]

and, exactly,

\[
 E[(Y-X)(Y-X)^T]=J-I\preceq(\beta-1)I.                \tag{1.3}
\]

The lower bound in (1.2) also follows directly from (1.3), since its
left-hand side is positive semidefinite.  Thus no lower curvature bound
is being assumed.

Let

\[
 \mathrm{Aff}=\operatorname{span}\{1,x_1,\ldots,x_n\}
 \subset L^2(\mu).
\]

For every \(g\in H^1(\mu)\) orthogonal to \(\mathrm{Aff}\), integration
by parts and (1.3) give

\[
 E\nabla g=E[gY]=E[g(Y-X)]
\]

and hence

\[
 \boxed{
 |E\nabla g|\le \sqrt{\beta-1}\,\|g\|_{L^2(\mu)}.}   \tag{1.4}
\]

Indeed, test the vector on an arbitrary unit vector \(u\), apply
Cauchy--Schwarz, and use
\(E\langle u,Y-X\rangle^2\le\beta-1\).

For nonsmooth convex \(V\) satisfying the upper Hessian bound in the
distributional sense, (1.1)--(1.4) follow by convolution of \(V\) and a
vanishing quadratic confinement, followed by local dominated convergence
and truncation.  The fixed-Gaussian class needs no such passage: its
density is positive and analytic, its score is square-integrable by
(1.2), and cutoff integration by parts is justified by first proving the
identities for compactly supported tests and then taking the cutoff limit
in \(L^2\).

## 2. Attained eigenfunctions

Let \(A=-L\) be the nonnegative self-adjoint operator associated with the
Dirichlet form \(\mathcal E(f,h)=\int\langle\nabla f,\nabla h\rangle d\mu\).
Suppose that a centered normalized \(f\) satisfies

\[
 Af=\lambda f,qquad \|f\|_2=1.
\]

Put

\[
 \ell=E[Xf(X)],\qquad g=f-\ell\cdot x,qquad
 \delta=\|g\|_2.
\]

Isotropy and the definition of \(\ell\) imply

\[
 g\perp\mathrm{Aff},qquad |ell|^2=1-\delta^2.        \tag{2.1}
\]

Testing the weak eigenfunction equation against coordinate functions
gives \(E\nabla f=\lambda\ell\).  Therefore

\[
 E\nabla g=-(1-\lambda)\ell.
\]

Combining this identity with (1.4) yields the rigidity estimate

\[
 \boxed{
 (1-\lambda)\sqrt{1-\delta^2}
 \le\sqrt{\beta-1}\,\delta.}                          \tag{2.2}
\]

In the fixed-Gaussian class \(\beta=2\).  If \(\delta\le1/2\), then

\[
 \boxed{\lambda\ge1-\frac1{\sqrt3}.}                 \tag{2.3}
\]

## 3. Continuous bottom spectral edge

The preceding conclusion does not require the spectral edge to be an
eigenvalue.  Let

\[
 b=\inf\operatorname{spec}(A|_{L^2_0(\mu)}).
\]

Isotropy gives \(b\le1\).  For every \(\varepsilon>0\), choose a centered
unit vector

\[
 f_\varepsilon\in
 \mathbf 1_{[b,b+\varepsilon]}(A)L^2_0(\mu).
\]

Set

\[
 z_\varepsilon=(A-b)f_\varepsilon,qquad
 \ell_\varepsilon=E[Xf_\varepsilon],qquad
 g_\varepsilon=f_\varepsilon-\ell_\varepsilon\cdot x,
\]

and \(\delta_\varepsilon=\|g_\varepsilon\|_2\).  Spectral calculus and
isotropy give

\[
 \|z_\varepsilon\|_2\le\varepsilon,qquad
 |E[Xz_\varepsilon]|\le\varepsilon,qquad
 |\ell_\varepsilon|^2=1-\delta_\varepsilon^2.         \tag{3.1}
\]

Testing \(Af_\varepsilon=bf_\varepsilon+z_\varepsilon\) weakly against
the coordinates gives

\[
 E\nabla g_\varepsilon
 =-(1-b)\ell_\varepsilon+E[Xz_\varepsilon].           \tag{3.2}
\]

Equations (1.4), (3.1), and (3.2) imply

\[
 \boxed{
 (1-b)\sqrt{1-\delta_\varepsilon^2}
 \le\sqrt{\beta-1}\,\delta_\varepsilon+\varepsilon.} \tag{3.3}
\]

In particular, if a sequence \(\varepsilon_k\downarrow0\) can be chosen
with \(\delta_{\varepsilon_k}\le1/2\), then for \(\beta=2\)

\[
 b\ge1-\frac1{\sqrt3}.                                \tag{3.4}
\]

Thus any possible small-gap sequence in the fixed-Gaussian class must
have every sufficiently accurate bottom spectral vector genuinely
non-affine:

\[
 \liminf_{\varepsilon\downarrow0}
 \operatorname{dist}_{L^2(\mu)}
 \left(f_\varepsilon,\mathrm{Aff}\right)>\frac12
\]

for every choice of bottom spectral vectors.  Closing this nonlinear
branch is the remaining proof task.

## 4. Sanity checks

For a standard Gaussian, \(\beta=1\), \(Y=X\), and (1.4) gives
\(E\nabla g=0\) exactly for every affine-orthogonal \(g\), as expected
from Hermite orthogonality.

For tensor products, (1.3) is block diagonal, so the estimate does not
acquire a dimension factor.  For hard-support inputs after fixed Gaussian
regularization, the posterior covariance formula proves
\(0\preceq D^2V\preceq2I\), so the argument remains valid without any
boundary integration on the original support.
