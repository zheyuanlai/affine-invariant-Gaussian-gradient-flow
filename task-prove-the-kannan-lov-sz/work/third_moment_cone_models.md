# Directional third moments in cones and one-parameter slice models

## 1. Question and outcome

Let \(X\) be centered and isotropic in \(\mathbb R^n\), and define

\[
 M_u=\mathbb E\big[(u\cdot X)(XX^T-I)\big]
     =\mathbb E\big[(u\cdot X)XX^T\big],\qquad |u|=1.       \tag{1.1}
\]

The second equality uses \(\mathbb E X=0\).  The issue is whether

\[
                \sup_{|u|=1}\|M_u\|_{\rm HS}\le C          \tag{1.2}
\]

holds for every isotropic log-concave law.  This note does **not** prove
(1.2) in complete generality.  It does prove it, without invoking KLS or a
high-dimensional Poincare inequality, for four fairly broad and strongly
asymmetric families:

1. all homothetic bodies over a centrally symmetric base, including cones,
   pyramids, frusta, and bodies of revolution;
2. all uniform affine box-slice bodies, with arbitrarily many unequal
   slopes;
3. all log-concave Gaussian slice models whose conditional covariance is
   an affine positive-definite matrix pencil; and
4. every Dirichlet density on a simplex whose parameters are at least one.

The Dirichlet result is exact: the bound is \(2\), and it is asymptotically
sharp.  The cone and a balanced anisotropic box family also have exact
formulas tending to \(2\).  Thus none of these natural cone, pyramid, or
common-latent-scale constructions gives a growing counterexample.

If a full-dimensional law is first centered and then whitened, the choice
of whitening is immaterial here: any two whitening maps differ by an
orthogonal map, and (1.2) is orthogonally invariant.  The same observation
applies on the affine hull of a lower-dimensional law.

The mechanism in the slice models is explicit.  After whitening, the
Frobenius norm in (1.2) is the Frobenius norm of the linear regression of
the conditional slice covariance on the axial coordinate.  Log-concavity
forces the marginal axial potential to pay a curvature cost equal to the
sum of squares of all simultaneous covariance changes.  Isotropy then
prevents that curvature cost from growing.

## 2. Two elementary one-dimensional facts

We use only the following standard one-dimensional consequences of
log-concavity.  They are recorded with enough detail to make clear that no
form of high-dimensional KLS is hidden here.

### Lemma 2.1 (isotropic core and third moment)

There are numerical constants \(b,K_0,K_3>0\) with the following property.
If \(Z\) has a log-concave density \(e^{-V}\) on an interval, with
\(\mathbb EZ=0\) and \(\mathbb EZ^2=1\), then

\[
 [-b,b]\subset\operatorname {int}(\operatorname {supp}Z),\qquad
 \operatorname {osc}_{[-b,b]}V\le K_0,qquad
 |\mathbb EZ^3|\le K_3.                                \tag{2.1}
\]

One may take, for example, \(b=1/20\), with nonoptimized numerical
\(K_0,K_3\).

**Justification.**  The elementary one-dimensional height estimate for a
variance-one log-concave density is \(\|e^{-V}\|_\infty\le1\), and
Grunbaum's inequality gives at least \(e^{-1}\) mass on either side of its
mean.  Chebyshev's inequality therefore gives

\[
 \mathbb P\{-3\le Z\le-b\},\ 
 \mathbb P\{b\le Z\le3\}
 \ge q:=e^{-1}-b-\frac19>0.                             \tag{2.2}
\]

There are points \(x_-\in[-3,-b]\) and \(x_+\in[b,3]\) at which the
density is at least \(q/3\).  Log-concavity makes it at least \(q/3\)
throughout \([x_-,x_+]\).  Together with the height upper bound this proves
the first two assertions, for instance with
\(K_0=\log(3/q)\).  The usual one-dimensional log-concave tail estimate
\(\mathbb P(|Z|>t)\le C e^{-ct}\), obtained by applying concavity of
\(\log f\) past fixed quantiles, gives the last assertion by integrating
the tail.  All constants here are one-dimensional.

The following direct consequence will be used repeatedly.

### Lemma 2.2 (core curvature budget)

Under the hypotheses of Lemma 2.1, if \(V\) is twice differentiable in the
interior and

\[
                         V''(z)\ge H
       \quad\text{for }|z|\le b/2,                     \tag{2.3}
\]

then

\[
                         H\le \frac{4K_0}{b^2}.         \tag{2.4}
\]

The same statement holds for distributional second derivatives.

**Proof.**  Convexity and the oscillation bound imply

\[
 V'_+(b/2)\le \frac{2K_0}{b},\qquad
 V'_-(-b/2)\ge-\frac{2K_0}{b}.
\]

The difference of these two derivatives is at least \(bH\), proving
(2.4).  Notice that this is a local, one-dimensional argument rather than
a spectral-gap estimate.

## 3. Exact block reduction for centered slices

The following calculation is the common algebra behind the first three
families.

### Lemma 3.1 (slice block formula)

Let \(X=(Z,U)\in\mathbb R\times\mathbb R^d\) be centered and isotropic.
Assume that, conditional on \(Z\), the law of \(U\) is centrally symmetric.
Put

\[
 \kappa=\mathbb EZ^3,\qquad
 C(z)=\mathbb E[UU^T\mid Z=z],\qquad
 D=\mathbb E[ZC(Z)].                                  \tag{3.1}
\]

Then, for \(u=(a,v)\), \(a^2+|v|^2=1\),

\[
 M_{(a,v)}=
 \begin{pmatrix}
       a\kappa &(Dv)^T\\
       Dv      &aD
 \end{pmatrix},                                      \tag{3.2}
\]

and consequently

\[
 \|M_{(a,v)}\|_{\rm HS}^2
 =a^2\big(\kappa^2+\|D\|_{\rm HS}^2\big)+2|Dv|^2,  \tag{3.3}
\]

\[
 \sup_{a^2+|v|^2=1}\|M_{(a,v)}\|_{\rm HS}^2
 =\max\left\{\kappa^2+\|D\|_{\rm HS}^2,
                    2\|D\|_{\rm op}^2\right\}.      \tag{3.4}
\]

**Proof.**  Conditional central symmetry kills every conditional odd
moment of \(U\).  Thus the upper-left, off-diagonal, and lower-right blocks
of \(\mathbb E[(aZ+v\cdot U)XX^T]\) are respectively
\(a\kappa\), \(Dv\), and \(aD\).  Squaring the blocks gives (3.3), and
optimizing the resulting quadratic form gives (3.4).

This lemma shows exactly what a counterexample inside a centered-slice
family would require: a one-dimensional log-concave marginal for which
\(\|\mathbb E[ZC(Z)]\|_{\rm HS}\) grows.

## 4. Affine positive-definite Gaussian covariance pencils

This family allows all eigenvalues of the conditional covariance to move
at different rates, so it is a useful test of possible Frobenius
accumulation.

### Theorem 4.1 (anisotropic Gaussian slices)

Let \(I\subset\mathbb R\) be an interval, let \(W:I\to(-\infty,+\infty]\)
be convex, and let

\[
                       R(t)=R_0+tR_1\succ0\quad(t\in I) \tag{4.1}
\]

be an affine symmetric positive-definite matrix pencil on \(\mathbb R^d\).
Assume that

\[
 d\mu(t,y)=Z_\mu^{-1}
 \exp\left[-W(t)-\frac12y^TR(t)^{-1}y\right]dt\,dy   \tag{4.2}
\]

is a probability law with nonzero axial variance.  After centering and
whitening \(\mu\), its directional third moments satisfy

\[
                   \sup_{|u|=1}\|M_u\|_{\rm HS}\le C \tag{4.3}
\]

for a numerical constant independent of \(d,W,R_0,R_1\).

**Log-concavity.**  The matrix-fractional function
\((R,y)\mapsto y^TR^{-1}y\) is jointly convex on positive-definite \(R\).
Since \(R(t)\) is affine and \(W\) is convex, (4.2) is log-concave.

**Whitening and exact third-moment matrix.**  Let

\[
 m=\mathbb ET,\quad \sigma^2=\operatorname {Var}T,\quad
 \bar R=\mathbb ER(T)=R(m),
\]

and set

\[
 Z={T-m\over\sigma},\qquad U=\bar R^{-1/2}Y,\qquad
 S=\sigma\bar R^{-1/2}R_1\bar R^{-1/2}.               \tag{4.4}
\]

Then \((Z,U)\) is isotropic and

\[
                  C(z)=\mathbb E[UU^T\mid Z=z]=I+zS.  \tag{4.5}
\]

In particular, \(D=S\) in Lemma 3.1.  Thus this model tests the full
Hilbert--Schmidt size of an arbitrary symmetric covariance velocity, not
only a scalar radial change.

**Curvature estimate.**  The marginal density of \(Z\) is proportional to

\[
 e^{-\widehat W(z)}\det(I+zS)^{1/2},                  \tag{4.6}
\]

where \(\widehat W\) is convex.  If \(V\) is its negative logarithm, then

\[
 V''(z)\ge\frac12\operatorname {tr}
       \left[S^2(I+zS)^{-2}\right].                   \tag{4.7}
\]

Lemma 2.1 implies that \(I+zS\succ0\) for every \(|z|\le b\).  Hence every
eigenvalue \(s\) of \(S\) satisfies \(|s|<1/b\).  For \(|z|\le b/2\),
\(|1+zs|\le3/2\), and therefore

\[
                  V''(z)\ge {2\over9}\|S\|_{\rm HS}^2. \tag{4.8}
\]

Lemma 2.2 yields

\[
                  \|S\|_{\rm HS}^2\le {18K_0\over b^2}. \tag{4.9}
\]

Finally \(|\kappa|\le K_3\), and (3.4) proves (4.3), for example with
\(C^2=K_3^2+36K_0/b^2\).  No Poincare inequality was used.

The determinant factor in (4.6) is the precise obstruction to a proposed
counterexample in which many conditional variances change by order one
while the axial variable retains order-one variance: its second derivative
charges the sum of the squared changes.

## 5. Uniform affine box slices

The analogous uniform model is a convex body with independently and
unequally moving pairs of facets.

Let \(I\) be an interval and let \(r_i(t)=\alpha_i+\beta_i t>0\) on its
interior.  Consider the convex body

\[
 K=\left\{(t,y)\in I\times\mathbb R^d:
                    |y_i|\le r_i(t),\ 1\le i\le d\right\}. \tag{5.1}
\]

Assume that \(K\) has finite positive volume and that its axial variance is
nonzero.  Convexity follows because all defining inequalities are affine.  Let
\((T,Y)\) be uniform on \(K\).  The density of \(T\) is proportional to
\(\prod_i r_i(t)\), and, conditional on \(T=t\), the coordinates \(Y_i\)
are independent uniforms on \([-r_i(t),r_i(t)]\).

### Theorem 5.1 (all affine box slices)

After centering and whitening the uniform law on (5.1),

\[
                   \sup_{|u|=1}\|M_u\|_{\rm HS}\le C \tag{5.2}
\]

with a numerical constant independent of the number and slopes of the
facets.

**Exact normalization.**  Put \(m=\mathbb ET\),
\(\sigma^2=\operatorname {Var}T\), \(Z=(T-m)/\sigma\), and

\[
 U_i={\sqrt3Y_i\over\sqrt{\mathbb Er_i(T)^2}},\qquad
 s_i={\sigma\beta_i\over r_i(m)}.                     \tag{5.3}
\]

The vector \((Z,U)\) is isotropic.  Since \(r_i\) is affine,

\[
 {r_i(T)\over r_i(m)}=1+s_iZ,\qquad
 {\mathbb Er_i(T)^2\over r_i(m)^2}=1+s_i^2.           \tag{5.4}
\]

Writing \(\kappa=\mathbb EZ^3\), the conditional covariance is diagonal,
and (3.1) gives \(D=\operatorname {diag}(c_1,\ldots,c_d)\), where

\[
              c_i={2s_i+\kappa s_i^2\over1+s_i^2}.    \tag{5.5}
\]

Thus (3.4) is an exact formula for this family.

**Why the slopes cannot accumulate.**  The standardized marginal density
of \(Z\) has the form

\[
             f(z)=C\prod_{i=1}^d(1+s_i z)             \tag{5.6}
\]

on an interval where every factor is positive.  Its potential satisfies

\[
             V''(z)=\sum_i{s_i^2\over(1+s_i z)^2}.    \tag{5.7}
\]

For \(|z|\le b/2\), Cauchy--Schwarz gives

\[
 (1+s_i z)^2\le(1+s_i^2)(1+z^2)\le2(1+s_i^2).
\]

Consequently, with

\[
                         B=\sum_i{s_i^2\over1+s_i^2}, \tag{5.8}
\]

one has \(V''\ge B/2\) on the core.  Lemma 2.2 gives

\[
                         B\le {8K_0\over b^2}.         \tag{5.9}
\]

Moreover, by (5.5) and the triangle inequality in \(\ell_2\),

\[
 \left(\sum_i c_i^2\right)^{1/2}
 \le
 \left(\sum_i{4s_i^2\over(1+s_i^2)^2}\right)^{1/2}
 +|\kappa|\left(\sum_i{s_i^4\over(1+s_i^2)^2}\right)^{1/2}
 \le(2+K_3)\sqrt B.                                  \tag{5.10}
\]

Equations (3.4), (5.9), and (5.10) prove the theorem.

This argument is useful beyond the result itself: even perfectly balanced
positive and negative facet slopes, for which the first derivative of the
slice volume cancels, produce the positive curvature sum in (5.7).

### Exact balanced family

Let \(d=2k\), \(I=[-1,1]\), and take

\[
 r_i(t)=1+t\quad(1\le i\le k),\qquad
 r_i(t)=1-t\quad(k<i\le2k).                            \tag{5.11}
\]

Then the axial density is proportional to \((1-t^2)^k\), and

\[
 \mathbb ET=0,\qquad \mathbb ET^2={1\over2k+3},\qquad \kappa=0. \tag{5.12}
\]

Formula (5.5) gives

\[
 c_i= {2\sqrt{2k+3}\over2k+4}\quad(i\le k),\qquad
 c_i=-{2\sqrt{2k+3}\over2k+4}\quad(i>k).             \tag{5.13}
\]

Therefore

\[
 \sup_{|u|=1}\|M_u\|_{\rm HS}^2
 =\sum_{i=1}^{2k}c_i^2
 ={2k(2k+3)\over(k+2)^2}\longrightarrow4.            \tag{5.14}
\]

There are \(2k\) nonzero diagonal third-moment coefficients, each of order
\(k^{-1/2}\), but their squared sum tends to four rather than growing.
This is the cleanest explicit stress test of the hoped-for universal
bound.

## 6. Homothetic cones, pyramids, and bodies of revolution

Let \(B\subset\mathbb R^d\) be a centrally symmetric convex body whose
uniform law is centered and has covariance \(I_d\).  Let \(r:I\to[0,\infty)\)
be concave, and consider

\[
 K_{r,B}=\{(t,r(t)y):t\in I,\ y\in B\}.               \tag{6.1}
\]

Equivalently, \(K_{r,B}=\{(t,x):\|x\|_B\le r(t)\}\), so it is convex.
Assume throughout that it has finite positive volume and nonzero axial
variance.
This includes ordinary bodies of revolution when \(B\) is a Euclidean
ball and generalized cones or pyramids when \(r\) is affine and vanishes
at an endpoint.

Let \(T\) have density proportional to \(r(t)^d\), and let \(Y\) be uniform
on \(B\), independently of \(T\).  A uniform point of (6.1) is
\((T,r(T)Y)\).  With

\[
 Z={T-\mathbb ET\over\sqrt{\operatorname {Var}T}},\qquad
 U={r(T)Y\over\sqrt{\mathbb Er(T)^2}},                \tag{6.2}
\]

the vector \((Z,U)\) is isotropic and Lemma 3.1 applies with

\[
 D=cI_d,\qquad
 c={\mathbb E[Zr(T)^2]\over\mathbb Er(T)^2}.          \tag{6.3}
\]

The one-dimensional Berwald--Borell moment inequality for a nonnegative
concave function on an interval gives

\[
 {\left(\int_I r^{d+4}\right)\left(\int_I r^d\right)
  \over\left(\int_I r^{d+2}\right)^2}
 \le{(d+3)^2\over(d+1)(d+5)}
 =1+{4\over(d+1)(d+5)}.                               \tag{6.4}
\]

For affine \(r\), (6.4) follows directly from the concavity in \(p\) of
\(\log((p+1)\int r^p)\); the general concave case is the standard
one-dimensional Berwald extension, with the same sharp constant.
Hence

\[
 \operatorname {Var}\left({r(T)^2\over\mathbb Er(T)^2}\right)
 \le {4\over(d+1)(d+5)}.                              \tag{6.5}
\]

Cauchy--Schwarz and \(\mathbb EZ^2=1\) imply

\[
              d c^2\le {4d\over(d+1)(d+5)}\le4.       \tag{6.6}
\]

Together with \(|\kappa|\le K_3\), the exact formula (3.4) proves a
universal bound for all (6.1).

### Exact cone over an arbitrary symmetric base

Take ambient dimension \(n=d+1\), \(I=[0,1]\), and \(r(t)=t\).  Then
\(T\sim\operatorname {Beta}(n,1)\),

\[
 \mathbb ET={n\over n+1},\quad
 \operatorname {Var}T={n\over(n+1)^2(n+2)},\quad
 \mathbb ET^2={n\over n+2}.                           \tag{6.7}
\]

The standardized skewness and covariance coefficient are

\[
 \kappa=-{2(n-1)\sqrt{n+2}\over(n+3)\sqrt n},\qquad
 c={2\sqrt{n+2}\over(n+3)\sqrt n},\qquad
 \kappa=-(n-1)c.                                      \tag{6.8}
\]

Thus, for \(u=(a,v)\),

\[
 M_u=\begin{pmatrix}a\kappa&cv^T\\cv&acI_{n-1}\end{pmatrix}, \tag{6.9}
\]

and for \(n\ge2\),

\[
 \sup_{|u|=1}\|M_u\|_{\rm HS}^2
 ={4(n-1)(n+2)\over(n+3)^2}\longrightarrow4.         \tag{6.10}
\]

The answer is independent of the dimension and geometry of the symmetric
base.  In particular, making the base itself a badly conditioned-looking
cube, crosspolytope, or smooth body does not change the result after affine
normalization.

## 7. Exact Dirichlet computation on arbitrary log-concave simplices

The preceding slice families have a distinguished axis.  Dirichlet laws
give a different, irreducibly coupled family and allow a complete tensor
calculation.

Let \(P=(P_1,\ldots,P_N)\) have the Dirichlet law with parameters
\(\alpha_1,\ldots,\alpha_N\), and write

\[
 A=\sum_i\alpha_i,\qquad q_i={\alpha_i\over A}.        \tag{7.1}
\]

The density on the affine simplex is proportional to
\(\prod_i p_i^{\alpha_i-1}\).  It is log-concave exactly in the range
\(\alpha_i\ge1\) for every \(i\), which is the range assumed below.

### Theorem 7.1 (Dirichlet third-moment theorem)

Let \(X\in\mathbb R^{N-1}\) be any affine whitening of \(P-q\).  Then, for
every unit vector \(u\), there is a unique vector
\(a=(a_1,\ldots,a_N)\) in the \(q\)-mean-zero subspace satisfying

\[
       \sum_iq_i a_i=0,\qquad \sum_iq_i a_i^2=1,       \tag{7.2}
\]

such that

\[
 \boxed{\displaystyle
 \|M_u\|_{\rm HS}^2
 ={4(A+1)\over(A+2)^2}\left(\sum_{i=1}^N a_i^2-2\right).} \tag{7.3}
\]

Consequently, if every \(\alpha_i\ge1\),

\[
                         \sup_{|u|=1}\|M_u\|_{\rm HS}<2. \tag{7.4}
\]

**Proof.**  Work in the finite Hilbert space \(L^2(q)\), and let
\(H=\{a:\mathbb E_q a=0\}\).  The covariance formula for a Dirichlet
vector shows that every standardized linear functional of \(P-q\) is

\[
 Z_a=\sqrt{A+1}\sum_i a_i(P_i-q_i),\qquad
 a\in H,\quad\|a\|_{L^2(q)}=1.                        \tag{7.5}
\]

For \(a,b,c\in H\), the third Dirichlet moment formula gives

\[
 \mathbb E[Z_aZ_bZ_c]
 ={2\sqrt{A+1}\over A+2}\,\mathbb E_q[abc].          \tag{7.6}
\]

Indeed, in the expansion of \(\mathbb E[(a\cdot P)(b\cdot P)(c\cdot P)]\),
all terms containing a factor \(\mathbb E_q a\), \(\mathbb E_q b\), or
\(\mathbb E_q c\) vanish, leaving

\[
 {2\sum_i\alpha_i a_ib_ic_i\over A(A+1)(A+2)}
 ={2\mathbb E_q[abc]\over(A+1)(A+2)}                 \tag{7.7}
\]

before the three standardizing factors are inserted.

Let \(P_H\) be orthogonal projection onto \(H\), and let \(D_a\) denote
multiplication by \(a\) in \(L^2(q)\).  Formula (7.6) says that the matrix
\(M_u\), as a bilinear form on the whitened tangent space, is

\[
 {2\sqrt{A+1}\over A+2}\,P_HD_aP_H.                 \tag{7.8}
\]

If \(Q\) is projection onto the constants, \(P_H=I-Q\).  Since
\(\mathbb E_qa=0\) and \(\mathbb E_qa^2=1\), a direct trace expansion gives

\[
\begin{aligned}
 \|P_HD_aP_H\|_{\rm HS}^2
 &=\operatorname {tr}(D_a(I-Q)D_a(I-Q))\\
 &=\operatorname {tr}(D_a^2)-2\langle1,D_a^21\rangle
   +\langle1,D_a1\rangle^2\\
 &=\sum_i a_i^2-2.                                    \tag{7.9}
\end{aligned}
\]

This proves (7.3).  Finally \(\alpha_i\ge1\) implies \(q_i\ge1/A\), so

\[
               \sum_i a_i^2\le A\sum_iq_i a_i^2=A.   \tag{7.10}
\]

Substitution into (7.3) gives a number strictly below four, proving (7.4).

For the uniform simplex, \(\alpha_i=1\), \(A=N\), and
\(\sum_i a_i^2=N\) for every unit tangent direction.  Therefore

\[
 \|M_u\|_{\rm HS}^2={4(N+1)(N-2)\over(N+2)^2}        \tag{7.11}
\]

for every \(u\), recovering the regular-simplex calculation and tending
to four.  The constant two is also approached by highly asymmetric beta
laws, for example with one Dirichlet parameter tending to infinity while
another remains equal to one.  Thus (7.4) has the correct scale and is
asymptotically sharp within this class.

## 8. Consequences for the general gate

These calculations rule out several plausible ways of falsifying (1.2).

* A single common scale variable cannot make all transverse diagonal
  coefficients order one: the slice-volume factor concentrates that
  variable, as (6.4)--(6.6) quantify.
* Balancing expanding and contracting directions cancels the first
  derivative of slice volume, but not its curvature.  Equations
  (4.7) and (5.7) charge the Hilbert--Schmidt square of the covariance
  velocity.
* Arbitrarily asymmetric simplex weights do not help.  In the exact
  Dirichlet formula, the potentially large unweighted norm
  \(\sum a_i^2\) is cancelled by the factor
  \(4(A+1)/(A+2)^2\); the log-concavity threshold \(\alpha_i\ge1\) is
  precisely what makes the cancellation dimension-free.

### Harmonic-spectrum audit

The projection thin-shell estimate by itself permits, abstractly, a
matrix whose positive or negative eigenvalues have the borderline profile

\[
                         |\lambda_j|\asymp j^{-1/2}.   \tag{8.1}
\]

Such a profile has bounded rank-\(r\) partial traces of order \(\sqrt r\)
but \(\sum_{j\le d}\lambda_j^2\asymp\log d\).  None of the families above
realizes it:

* in the Gaussian pencil family, the transverse eigenvalues are exactly
  those of \(S\), and (4.9) gives \(\sum_j\lambda_j^2\le C\);
* in the affine box family, they are the \(c_i\), and (5.10) gives
  \(\sum_i c_i^2\le C\), even for arbitrary multiscale positive and
  negative slopes;
* in a cone, all \(n-1\) transverse eigenvalues equal
  \(c\asymp n^{-1}\), while the single axial eigenvalue is of constant
  order; and
* for every log-concave Dirichlet law, (7.3) bounds the full Frobenius
  square by four before any eigenvalue ordering is used.

In particular, assigning raw box or Gaussian slice slopes
\(s_j\asymp j^{-1/2}\) does not create (8.1) after isotropic normalization.
The log-volume curvature is already of order
\(\sum_{j\le d}s_j^2\asymp\log d\), so the axial variance shrinks by the
reciprocal scale and the normalized spectrum has bounded square sum.  A
genuine harmonic-spectrum example, if one exists, must therefore use
non-affine, noncommuting slice changes or nonzero conditional third
moments in a way not charged by these log-determinant curvature identities.

The remaining general difficulty is now sharper.  For a general
log-concave law and a direction \(u\), disintegrate over
\(Z=u\cdot X\).  Even after removing conditional means by an affine shear,
the slice covariance need not be an affine matrix pencil, the conditional
third moments need not vanish, and rotating/nonhomothetic slices introduce
terms absent from Lemma 3.1.  A complete proof of (1.2) would need a
matrix-valued curvature or mixed-volume inequality controlling all of
those terms.  The calculations above establish the desired estimate for
large nonproduct classes, but they do not supply that missing general
inequality and therefore do not by themselves resolve the KLS task.
