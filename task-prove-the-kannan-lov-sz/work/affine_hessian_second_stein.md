# Affine--Hessian rigidity and second-order Stein kernels

## 0. Verdict

Let \(\mu\) be a centered isotropic log-concave probability on
\(\mathbb R^n\).  The proposed affine--Hessian inequality is

\[
 \boxed{\quad
 \left\|\mathbb E\left[
 X\otimes\nabla u+\nabla u\otimes X\right]\right\|_{\rm HS}^2
 \le C\,\mathbb E\|D^2u\|_{\rm HS}^2.
 \quad}                                                \tag{AH}
\]

This report does not prove (AH) for every log-concave law, but it gives the
following exact results.

1. KLS with Poincare constant \(C_P(\mu)\) implies (AH) with
   \(C=4C_P(\mu)\).  Conversely, (AH) sees only the projection of the
   centered gradient onto linear coordinates and has a large kernel.  For a
   near-linear first eigenfunction it controls a quadratic moment of the
   eigenfunction, but it is completely blind to the load-bearing mean
   gradient of the affine-orthogonal residual.  Thus (AH) alone does not
   close the near-linear eigenfunction route.

2. (AH) is exactly equivalent, by Hilbert-space duality, to the existence
   for every symmetric matrix \(B\) of a \(B\)-specific second-order Stein
   field \(K_B\) satisfying

   \[
    \mathbb E[K_B:D^2u]
      =2\mathbb E[(BX)\cdot\nabla u],
    \qquad
    \mathbb E\|K_B\|_{\rm HS}^2\le C\|B\|_{\rm HS}^2.  \tag{0.1}
   \]

   A common first-order Stein kernel \(\tau\) gives
   \(K_B=2\operatorname {sym}(B\tau)\).  A uniform \(L^2\) operator-norm
   bound on \(\tau\) is sufficient but unnecessary.

3. There is a complete dimension-free proof for every unconditional
   isotropic log-concave measure, including dependent unconditional convex
   bodies.  Conditional one-dimensional Stein kernels give a diagonal
   \(\tau\) with

   \[
      \sup_i\mathbb E\tau_i^2\le C,
   \]

   and hence (AH).  For the isotropic cube the resulting constant is
   \(24/5\); for the isotropic crosspolytope it is at most \(16\).

4. Explicit common kernels also prove (AH) with the following constants:

   \[
   \begin{array}{c|c}
   \text{law}&\text{valid constant in (AH)}\\ \hline
   \text{standard Gaussian}&4\\
   \text{isotropic Euclidean ball}&<8\\
   \text{isotropic regular simplex}&<8\\
   \text{isotropic circular cone, and every affine shear of it}&\le12.
   \end{array}                                         \tag{0.2}
   \]

   The cone calculation gives a genuinely non-centrally-symmetric,
   nonproduct test.

5. Three tempting universal sublemmas are false.

   * Central symmetry does not make
     \(\mathbb E[X\otimes\nabla u]\) symmetric.  An explicit even quartic
     on the cube gives unequal transposed entries.
   * The canonical simplex Stein kernel has
     \(\mathbb E\|\tau\|_{\rm op}^2\gtrsim(\log n)^2\), even though its
     \(B\)-specific averaged action proves (AH) with constant below eight.
   * The universal iid-chord representation produces a valid second-order
     kernel, but on the cube its \(L^2\) norm grows exponentially in the
     dimension.  Thus density domination or an \(L^2\) estimate for that
     particular chord kernel cannot prove (AH).

The positive results show that \(B\)-specific or symmetry-averaged kernels
are materially stronger than a common pointwise operator bound.  The
remaining general obstruction is controlling those kernels when the
conditional fiber centers move and no irreducible symmetry averages their
large eigenvectors.

## 1. Exact operator relations

Write

\[
 \mathcal T_\mu u
 =\mathbb E[X\otimes\nabla u+\nabla u\otimes X]
 \in\mathbb S^n.                                      \tag{1.1}
\]

### 1.1 KLS implies (AH)

Put \(M=\mathbb E[X\otimes\nabla u]\).  For each \(j\), the coordinates
\(X_1,\ldots,X_n\) are an orthonormal family in \(L^2(\mu)\).  Bessel's
inequality gives

\[
\begin{aligned}
 \sum_i\bigl(\mathbb E[X_i\partial_j u]\bigr)^2
 &\le\operatorname {Var}_\mu(\partial_j u).
\end{aligned}                                         \tag{1.2}
\]

Consequently

\[
 \|M\|_{\rm HS}^2
 \le\sum_j\operatorname {Var}_\mu(\partial_j u)
 \le C_P(\mu)\,\mathbb E\|D^2u\|_{\rm HS}^2.           \tag{1.3}
\]

Since \(\|M+M^T\|_{\rm HS}^2\le4\|M\|_{\rm HS}^2\),

\[
 \|\mathcal T_\mu u\|_{\rm HS}^2
 \le4C_P(\mu)\mathbb E\|D^2u\|_{\rm HS}^2.             \tag{1.4}
\]

No converse Poincare estimate follows from Bessel: (AH) retains only the
linear-coordinate projection of each centered derivative.

The numerical constant in any universal (AH) must be at least four.  For

\[
 u(x)=\frac12x^TBx,\qquad B=B^T,
\]

isotropy gives

\[
 \mathcal T_\mu u=2B,\qquad D^2u=B.                    \tag{1.5}
\]

Thus equality holds in (AH) with \(C=4\), for every isotropic law and every
quadratic \(u\).

### 1.2 Exact first-eigenfunction specialization

Let \(f\) be a normalized first eigenfunction:

\[
 -L_\mu f=\lambda f,\qquad
 \mathbb Ef=0,\qquad\mathbb Ef^2=1,\qquad
 \mathbb E|\nabla f|^2=\lambda.                        \tag{1.6}
\]

Testing the weak eigenfunction equation against \(x_ix_j\) gives

\[
 \mathbb E[X_i\partial_jf+X_j\partial_if]
 =\lambda\mathbb E[fX_iX_j].                           \tag{1.7}
\]

Hence

\[
 \boxed{\quad
 \mathcal T_\mu f
 =\lambda\,\mathbb E[f(XX^T-I)].
 \quad}                                                \tag{1.8}
\]

The Bochner--Reilly identity gives

\[
 \mathbb E\|D^2f\|_{\rm HS}^2\le\lambda^2.             \tag{1.9}
\]

Therefore (AH) would imply only

\[
 \left\|\mathbb E[f(XX^T-I)]\right\|_{\rm HS}^2\le C.  \tag{1.10}
\]

To compare this with the near-linear gate, put

\[
 a=\mathbb E[Xf],\qquad
 \ell=a\cdot x,\qquad r=f-\ell.
\]

Then

\[
 \mathbb Er=0,\qquad\mathbb E[Xr]=0,\qquad
 \mathbb E\nabla r=-(1-\lambda)a,\qquad D^2r=D^2f.    \tag{1.11}
\]

Because \(\mathbb EX=0\), adding an affine function does not change
\(\mathcal T_\mu\):

\[
                         \mathcal T_\mu r=\mathcal T_\mu f. \tag{1.12}
\]

Thus (AH) is insensitive to the mean gradient
\(-(1-\lambda)a\), which is the exact obstruction in (1.11).  Even if
\(\|r\|_2\to0\), (AH) does not yield the mean-gradient estimate needed to
force \(\lambda\) away from zero.  Its use in the non-Gaussian slice route
is different: it converts a Brascamp--Lieb deficit into a covariance
velocity after solving a Poisson equation.

For a centrally symmetric law, every odd \(u\) has even gradient, so
\(\mathcal T_\mu u=0\).  This parity kernel further demonstrates why (AH)
cannot by itself be read as a full Poincare inequality.

## 2. \(B\)-specific second-order Stein kernels

For \(B\in\mathbb S^n\), define

\[
 \Lambda_B(u)
 =\langle B,\mathcal T_\mu u\rangle_{\rm HS}
 =2\mathbb E[(BX)\cdot\nabla u].                       \tag{2.1}
\]

The functional annihilates every affine \(u\).  On the quotient by affine
functions, equip smooth functions with the Hessian seminorm

\[
                         \|u\|_{\dot H^2(\mu)}
 =\bigl(\mathbb E\|D^2u\|_{\rm HS}^2\bigr)^{1/2}.       \tag{2.2}
\]

By Riesz representation on the closure of the Hessian range, the following
are equivalent with the same constant \(C\).

* For every smooth \(u\),
  \[
     |\Lambda_B(u)|^2\le C\|B\|_{\rm HS}^2
                   \mathbb E\|D^2u\|_{\rm HS}^2.
  \]
* There is a symmetric matrix field \(K_B\in L^2(\mu;\mathbb S^n)\) such
  that
  \[
   \mathbb E[K_B:D^2u]=\Lambda_B(u),\qquad
   \mathbb E\|K_B\|_{\rm HS}^2\le C\|B\|_{\rm HS}^2.   \tag{2.3}
  \]

For a smooth density \(p\) on a domain \(\Omega\), (2.3) has the weak PDE
form

\[
 -\operatorname {div}(pK_B)=2pBx                         \tag{2.4}
\]

up to addition of a divergence-free matrix field, with zero normal flux
\(K_Bn=0\) on \(\partial\Omega\) when there is a boundary.  Equation (2.4)
is the \(B\)-specific second-order Stein problem.

Suppose now that \(\tau\) is a first-order Stein kernel, not necessarily
symmetric:

\[
 \mathbb E[X_i\phi]
 =\mathbb E\sum_k\tau_{ik}\partial_k\phi.              \tag{2.5}
\]

Applying (2.5) to \(\phi=\partial_j u\) gives

\[
 \mathbb E[X\otimes\nabla u]=\mathbb E[\tau D^2u].
\]

For symmetric \(B\), a valid choice in (2.3) is

\[
 \boxed{\quad K_B=2\operatorname {sym}(B\tau).\quad}   \tag{2.6}
\]

In particular, if

\[
                    \mathbb E[\tau\tau^T]\preceq K I, \tag{2.7}
\]

then

\[
\begin{aligned}
 \mathbb E\|K_B\|_{\rm HS}^2
 &\le4\mathbb E\|B\tau\|_{\rm HS}^2\\
 &=4\operatorname {Tr}\bigl(B^2\mathbb E[\tau\tau^T]\bigr)
 \le4K\|B\|_{\rm HS}^2.                               \tag{2.8}
\end{aligned}
\]

This averaged matrix condition is much weaker than
\(\mathbb E\|\tau\|_{\rm op}^2\le K\).  The simplex in Section 5 separates
the two conditions by a factor of order \((\log n)^2\).

## 3. Uniform polytopes and boundary bookkeeping

If \(\mu\) is uniform on a centered convex body \(K\), integration by parts
gives, entrywise,

\[
\begin{aligned}
 |K|\,\mathcal T_\mu u{}_{ij}
 &=
 \int_{\partial K}u(x)(x_in_j+x_jn_i)\,dS(x)
 -2\delta_{ij}\int_Ku(x)\,dx.                          \tag{3.1}
\end{aligned}
\]

Equivalently,

\[
 |K|\,\Lambda_B(u)
 =2\int_{\partial K}u(x)\,Bx\cdot n\,dS
 -2\operatorname {Tr}B\int_Ku.                        \tag{3.2}
\]

The separate terms in (3.2) may grow with the number of facets or with
surface-to-volume ratio.  Their affine cancellations are essential:
\(\Lambda_B\) annihilates affine functions exactly.  Taking absolute values
facet by facet therefore loses the structure needed for (AH).

For a polytope, a field \(K_B\) satisfying (2.4) incorporates all facet
cancellations at once.  The zero-flux condition is the precise replacement
for informal cancellation of boundary terms.  Sections 4 and 5 construct
such fields for unconditional polytopes and simplices.

## 4. A complete proof for unconditional dependent measures

Assume that \(\mu\) is unconditional:

\[
 (X_1,\ldots,X_n)\stackrel d=
 (\varepsilon_1X_1,\ldots,\varepsilon_nX_n)
 \quad\text{for every }\varepsilon_i\in\{-1,1\}.       \tag{4.1}
\]

No independence is assumed.

### 4.1 Conditional diagonal Stein kernel

Fix \(i\).  Conditional on \(X_{-i}=x_{-i}\), the law of \(X_i\) is an
even one-dimensional log-concave probability on an interval.  Let
\(\tau_i(x)\) be its canonical one-dimensional Stein kernel.  Then

\[
 \mathbb E[X_i\phi(X)\mid X_{-i}]
 =\mathbb E[\tau_i(X)\partial_i\phi(X)\mid X_{-i}].    \tag{4.2}
\]

Consequently

\[
 \tau(X)=\operatorname {diag}(\tau_1(X),\ldots,\tau_n(X)) \tag{4.3}
\]

is a first-order Stein kernel for the full dependent law.

We use the following elementary one-dimensional estimate.

**Lemma 4.1.**  If \(Z\) is centered and log-concave with variance \(v\),
and \(\tau_Z\) is its canonical Stein kernel, then

\[
                         \mathbb E\tau_Z^2\le C_0v^2.  \tag{4.4}
\]

**Proof.**  After scaling, take \(v=1\).  The canonical formula

\[
 \tau_Z(z)
 ={1\over p(z)}\int_z^\infty x p(x)\,dx
\]

and the increasing-hazard and reverse-hazard bounds for a one-dimensional
log-concave density give

\[
                         0\le\tau_Z(z)\le C(1+|z|).
\]

A centered variance-one log-concave variable has a numerical exponential
tail.  Squaring the last display and integrating proves (4.4) for \(v=1\).
The general statement follows from the scaling identity
\(\tau_{\sqrt v Z}(\sqrt vz)=v\tau_Z(z)\). \(\square\)

Let

\[
 v_i(X_{-i})=\mathbb E[X_i^2\mid X_{-i}].
\]

Lemma 4.1, Jensen's inequality, and the universal fourth-moment bound for a
one-dimensional isotropic log-concave marginal give

\[
\begin{aligned}
 \mathbb E\tau_i^2
 &\le C_0\mathbb E v_i^2\\
 &\le C_0\mathbb E X_i^4
 \le C_1.                                               \tag{4.5}
\end{aligned}
\]

Here the middle inequality is pointwise:
\((\mathbb E[X_i^2\mid X_{-i}])^2
\le\mathbb E[X_i^4\mid X_{-i}]\).

It follows from (4.3)--(4.5) that

\[
 \mathbb E[\tau\tau^T]
 =\operatorname {diag}(\mathbb E\tau_1^2,\ldots,
                        \mathbb E\tau_n^2)
 \preceq C_1I.                                         \tag{4.6}
\]

Equations (2.6)--(2.8) prove (AH) for every unconditional isotropic
log-concave probability, with the numerical constant \(4C_1\).  This proof
includes nonsmooth densities and convex-body boundaries by conditional
one-dimensional approximation.

### 4.2 Explicit uniform-body fibers

Let \(K\) be an unconditional convex body.  For fixed \(x_{-i}\), its
\(i\)-th coordinate fiber is

\[
                         [-\rho_i(x_{-i}),\rho_i(x_{-i})].
\]

The canonical kernel of the uniform law on that interval is

\[
 \tau_i(x)={\rho_i(x_{-i})^2-x_i^2\over2}.             \tag{4.7}
\]

Thus (4.3) is completely explicit for every unconditional polytope.  At
almost every boundary point, if the outer normal has \(n_i\ne0\), that
point is an endpoint of the \(i\)-th coordinate fiber and \(\tau_i=0\).
Therefore \(\tau n=0\), verifying the polyhedral zero-flux condition in
(2.4) directly.

### 4.3 Cube

For the isotropic cube

\[
                         K=[-\sqrt3,\sqrt3]^n,
\]

\[
 \tau_i(x)={3-x_i^2\over2}.
\]

Since \(\mathbb EX_i^2=1\) and \(\mathbb EX_i^4=9/5\),

\[
 \mathbb E\tau_i^2
 ={1\over4}\left(9-6+{9\over5}\right)
 ={6\over5}.                                          \tag{4.8}
\]

Hence

\[
 \boxed{\quad
 \|\mathcal T_\mu u\|_{\rm HS}^2
 \le {24\over5}\,\mathbb E\|D^2u\|_{\rm HS}^2.
 \quad}                                                \tag{4.9}
\]

### 4.4 Crosspolytope

Let

\[
 K=\left\{x:\sum_{i=1}^n|x_i|\le R\right\},
\qquad
 R^2={(n+1)(n+2)\over2},                               \tag{4.10}
\]

so that the uniform law is isotropic.  Put

\[
 L_i=R-\sum_{j\ne i}|X_j|.
\]

Conditional on \(X_{-i}\), \(X_i\) is uniform on \([-L_i,L_i]\), and

\[
 \tau_i={L_i^2-X_i^2\over2}.                           \tag{4.11}
\]

The normalized absolute coordinates together with the slack have the
Dirichlet\((1,\ldots,1)\) law.  Therefore \(L_i/R\) has the
\(\operatorname {Beta}(2,n-1)\) law.  Conditional integration gives

\[
 \mathbb E[\tau_i^2\mid L_i]
 ={2\over15}L_i^4,                                    \tag{4.12}
\]

while

\[
\begin{aligned}
 \mathbb EL_i^4
 &=R^4\,{(2)_4\over(n+1)_4}\\
 &=30\,{(n+1)(n+2)\over(n+3)(n+4)}.                   \tag{4.13}
\end{aligned}
\]

Consequently

\[
 \mathbb E\tau_i^2
 =4\,{(n+1)(n+2)\over(n+3)(n+4)}
 <4.                                                   \tag{4.14}
\]

The crosspolytope therefore satisfies

\[
 \boxed{\quad
 \|\mathcal T_\mu u\|_{\rm HS}^2
 \le16\,\mathbb E\|D^2u\|_{\rm HS}^2.
 \quad}                                                \tag{4.15}
\]

This is a dependent, nonproduct example in every dimension.

## 5. Ball and simplex: averaged kernels beat operator norm

### 5.1 Isotropic Euclidean ball

Let \(X\) be uniform on the ball of radius

\[
                         R=\sqrt{n+2}.
\]

The scalar matrix field

\[
 \tau(x)={R^2-|x|^2\over2}\,I                         \tag{5.1}
\]

is a Stein kernel.  Indeed, its scalar coefficient vanishes on the boundary
and has gradient \(-x\), so integration by parts gives

\[
                         \mathbb E[X\phi]
 =\mathbb E[\tau\nabla\phi].
\]

The radial moments

\[
 \mathbb E|X|^2=n,\qquad
 \mathbb E|X|^4={n(n+2)^2\over n+4}
\]

give

\[
 \mathbb E\left({R^2-|X|^2\over2}\right)^2
 ={2(n+2)\over n+4}<2.                                 \tag{5.2}
\]

Thus

\[
 \boxed{\quad
 \|\mathcal T_\mu u\|_{\rm HS}^2
 <8\,\mathbb E\|D^2u\|_{\rm HS}^2.
 \quad}                                                \tag{5.3}
\]

### 5.2 Regular simplex

Put \(N=n+1\).  Let

\[
 P=(P_1,\ldots,P_N)\sim\operatorname {Dirichlet}(1,\ldots,1)
\]

on the standard \(N\)-vertex simplex, and let \(U:\mathbb R^n\to\mathbb
R^N\) be an isometric embedding onto \({\bf1}^\perp\).  The isotropic regular
simplex may be represented as

\[
 X=\sqrt{N(N+1)}\,U^T\left(P-{1\over N}{\bf1}\right).  \tag{5.4}
\]

The Dirichlet integration-by-parts identity says that

\[
 \tau_P(P)={1\over N}\left(\operatorname {diag}P-PP^T\right) \tag{5.5}
\]

is a Stein kernel on \({\bf1}^\perp\).  After the affine transformation
(5.4), a Stein kernel for \(X\) is

\[
 \tau_X=(N+1)\,U^T
       \left(\operatorname {diag}P-PP^T\right)U.       \tag{5.6}
\]

Permutation symmetry acts irreducibly on \({\bf1}^\perp\).  Hence

\[
                         \mathbb E\tau_X^2=c_NI_n.     \tag{5.7}
\]

The constant is explicit.  If
\(C(P)=\operatorname {diag}P-PP^T\), then

\[
 \operatorname {Tr}C(P)^2
 =\sum_iP_i^2-2\sum_iP_i^3+\left(\sum_iP_i^2\right)^2. \tag{5.8}
\]

The Dirichlet moments give

\[
\begin{aligned}
 \mathbb E\sum_iP_i^2&={2\over N+1},\\
 \mathbb E\sum_iP_i^3&={6\over(N+1)(N+2)},\\
 \mathbb E\left(\sum_iP_i^2\right)^2
 &={4(N+5)\over(N+1)(N+2)(N+3)}.
\end{aligned}                                         \tag{5.9}
\]

Therefore

\[
 \mathbb E\operatorname {Tr}C(P)^2
 ={2(N-1)\over(N+1)(N+3)},                            \tag{5.10}
\]

and

\[
 c_N={1\over N-1}\mathbb E\operatorname {Tr}\tau_X^2
 ={2(N+1)\over N+3}<2.                                \tag{5.11}
\]

Equations (2.6)--(2.8) prove

\[
 \boxed{\quad
 \|\mathcal T_\mu u\|_{\rm HS}^2
 \le {8(N+1)\over N+3}\,
       \mathbb E\|D^2u\|_{\rm HS}^2
 <8\,\mathbb E\|D^2u\|_{\rm HS}^2.
 \quad}                                                \tag{5.12}
\]

This gives a direct second-order Stein proof for the simplex, without using
its Poincare constant.

### 5.3 A false common-kernel shortcut

The averaged estimate (5.7) does not come from a dimension-free operator
bound.  Let

\[
                         P_{\max}=\max_iP_i.
\]

On the event \(P_{\max}\le1/2\), testing the categorical covariance
\(C(P)\) on \((e_i-e_j)/\sqrt2\), where \(i\) maximizes \(P_i\) and
\(P_j\le P_i\), gives

\[
                         \|C(P)\|_{\rm op}\ge {P_{\max}\over4}. \tag{5.13}
\]

Using the exponential representation

\[
                         P_i={E_i\over\sum_jE_j},
\qquad E_i\stackrel{\rm iid}{\sim}\operatorname {Exp}(1),
\]

one has, with probability tending to one,

\[
                         P_{\max}\ge {c\log N\over N},
\qquad P_{\max}\le{1\over2}.                           \tag{5.14}
\]

For example,
\[
 \mathbb P\{\max_iE_i<\tfrac12\log N\}
 =(1-N^{-1/2})^N\to0,
\]
and a Chernoff bound gives
\(\mathbb P\{\sum_iE_i>2N\}\le e^{-cN}\).
Combining (5.6), (5.13), and (5.14) yields

\[
 \boxed{\quad
 \mathbb E\|\tau_X\|_{\rm op}^2\ge c(\log N)^2.
 \quad}                                                \tag{5.15}
\]

Thus the implication

\[
 \text{dimension-free (AH)}
 \Longrightarrow
 \text{a common Stein kernel with bounded }L^2\text{ operator norm}
\]

is false even for regular simplices.  The \(B\)-specific averaged action in
(5.7) is the correct scale.

## 6. A non-symmetric cone and an inherited Stein kernel

Let the total dimension be \(n\ge2\), put \(d=n-1\), and let \(W\) be
uniform on the isotropic \(d\)-dimensional Euclidean ball.  Let \(S\) be
independent of \(W\), with density

\[
                         ns^{n-1}{\bf1}_{[0,1]}(s).
\]

The uniform law on the cone with apex at the origin and base
\(\{1\}\times\operatorname {supp}W\) is represented by \((S,SW)\).
Set

\[
\begin{aligned}
 \mu_S&={n\over n+1},&
 \sigma_S^2&={n\over(n+1)^2(n+2)},\\
 Y&={S-\mu_S\over\sigma_S},&
 Z&=cSW,\qquad c=\sqrt{n+2\over n}.
\end{aligned}                                         \tag{6.1}
\]

Then \(X=(Y,Z)\) is isotropic.

### 6.1 Explicit nonsymmetric kernel

The canonical Stein kernel of \(Y\) is

\[
 a(S)={(n+1)(n+2)\over n}S(1-S).                      \tag{6.2}
\]

Indeed, the centered Stein kernel of \(S\) is
\(S(1-S)/(n+1)\), and standardization divides it by
\(\sigma_S^2\).

Let the ball Stein kernel be

\[
 \tau_W(W)=\theta(W)I_d,\qquad
 \theta(W)={d+2-|W|^2\over2}.                          \tag{6.3}
\]

Define the upper-triangular matrix field

\[
 \tau_X=
 \begin{pmatrix}
  a(S)&r(S,W)^T\\
  0&Q(S,W)
 \end{pmatrix},                                       \tag{6.4}
\]

where

\[
 r(S,W)=c^2S(1-S)W,\qquad
 Q(S,W)=c^2S^2\theta(W)I_d.                            \tag{6.5}
\]

This is a Stein kernel for \(X\).  For the transverse rows, condition on
\(S\) and apply the Stein identity of \(W\).  For the axial row, if
\(\phi=\phi(Y,cSW)\), differentiation with respect to \(Y\) gives

\[
 {d\over dY}\phi(Y,cSW)
 =\partial_Y\phi+{\sigma_S\over S}Z\cdot\nabla_Z\phi.
\]

The off-diagonal row in (6.4) is exactly
\[
 a(S){\sigma_S\over S}Z
 =c^2S(1-S)W.
\]

### 6.2 Uniform \(L^2\) matrix bound

The beta moments give

\[
 \mathbb E[S^2(1-S)^2]
 ={2n\over(n+2)(n+3)(n+4)},\qquad
 \mathbb ES^4={n\over n+4}.                            \tag{6.6}
\]

Consequently

\[
\begin{aligned}
 \mathbb Ea(S)^2
 &={2(n+1)^2(n+2)\over n(n+3)(n+4)}<2,\\
 \mathbb E|r(S,W)|^2
 &={2(n-1)(n+2)\over n(n+3)(n+4)}<1.                  \tag{6.7}
\end{aligned}
\]

For the isotropic \(d\)-ball,

\[
 \mathbb E\theta(W)^2={2(d+2)\over d+4}
 ={2(n+1)\over n+3}<2.                                \tag{6.8}
\]

Thus

\[
 \mathbb E[Q Q^T]
 ={(n+2)^2\over n(n+4)}\,
   \mathbb E\theta(W)^2\,I_d
 \prec {8\over3}I_d.                                  \tag{6.9}
\]

Rotational symmetry gives \(\mathbb E[rQ^T]=0\).  Equations
(6.7)--(6.9) imply

\[
                         \mathbb E[\tau_X\tau_X^T]\preceq3I_n. \tag{6.10}
\]

Therefore

\[
 \boxed{\quad
 \|\mathcal T_\mu u\|_{\rm HS}^2
 \le12\,\mathbb E\|D^2u\|_{\rm HS}^2.
 \quad}                                                \tag{6.11}
\]

The isotropic image of any affine shear of this cone differs only by an
orthogonal transformation, so (6.11) covers skew affine cones as well.
This example has neither central symmetry nor product structure.

The construction also gives an inheritance statement.  If an isotropic
base \(W\) admits a Stein kernel with
\(\mathbb E[\tau_W\tau_W^T]\preceq K I\), the same upper-triangular formula
gives a cone kernel with an averaged matrix bound depending only on \(K\),
not on \(n\).
