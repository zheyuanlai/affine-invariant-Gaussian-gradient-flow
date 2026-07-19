# Noncommuting Gaussian slices: exact convexity and curvature compensation

## 1. Main result

Consider a density on \(\mathbb R\times\mathbb R^d\) of the form

\[
 p(t,y)=Z^{-1}\exp\left[-W(t)-\frac12y^TQ(t)y\right],
 \qquad Q(t)\succ0.                                    \tag{1.1}
\]

Write \(R(t)=Q(t)^{-1}\), so that, conditional on \(T=t\), \(Y\) is a
centered Gaussian with covariance \(R(t)\).  The principal result of this
note is:

> **Theorem.**  Suppose (1.1) is log-concave, has finite second moments,
> and \(\operatorname {Var}T>0\).  Center and whiten \(T\), and whiten
> \(Y\) by its unconditional covariance.  For the resulting isotropic
> vector \(X\),
> \[
>             \sup_{|u|=1}
>             \left\|\mathbb E[(u\cdot X)(XX^T-I)]\right\|_{\rm HS}
>             \le C,                                  \tag{1.2}
> \]
> where \(C\) is numerical.  It is independent of \(d,W\), of the
> condition numbers of the matrices, and of the rate at which their
> eigenspaces rotate.

This strictly extends the affine covariance-pencil calculation.  The
covariance path may be nonlinear, and the matrices \(R(t)\), \(R'(t)\),
and \(R''(t)\) need not commute.

The proof does not use KLS or a high-dimensional Poincare inequality.  It
uses:

1. the exact equivalence between joint convexity and Loewner concavity of
   \(R\);
2. the curvature identity for the one-dimensional marginal;
3. the canonical one-dimensional Stein kernel; and
4. a simple growth lemma for a positive matrix-concave function with mean
   \(I\).

## 2. Exact joint-convexity condition

Let

\[
                  \mathcal V(t,y)=W(t)+\frac12y^TQ(t)y. \tag{2.1}
\]

Assume first that \(W,Q\) are \(C^2\).  Its block Hessian is

\[
 D^2\mathcal V=
 \begin{pmatrix}
 W''+\frac12y^TQ''y &(Q'y)^T\\
 Q'y &Q
 \end{pmatrix}.                                      \tag{2.2}
\]

Since \(Q\succ0\), the Schur complement shows that (2.2) is positive
semidefinite for every \(y\) if and only if

\[
 W''\ge0,\qquad
 \frac12Q''-Q'Q^{-1}Q'\succeq0.                       \tag{2.3}
\]

Differentiating \(Q=R^{-1}\) gives

\[
\begin{aligned}
 Q'&=-QR'Q,\\
 Q''&=2QR'QR'Q-QR''Q,
\end{aligned}
\]

and hence

\[
 \frac12Q''-Q'Q^{-1}Q'=-\frac12QR''Q.                 \tag{2.4}
\]

We have proved the exact equivalence

\[
 \boxed{\quad \mathcal V\text{ is jointly convex}
 \quad\Longleftrightarrow\quad
 W\text{ is convex and }R''\preceq0.\quad}            \tag{2.5}
\]

For nonsmooth paths, (2.5) has the same distributional meaning:
\(R\) is matrix concave when \(v^TR(t)v\) is concave for every \(v\).
Approximation by convolution on compact subintervals gives the formulas
below; hard endpoints are obtained by monotone convex approximation.

### Genuine rotation is allowed

Matrix concavity does not force commutation.  For example, in dimension
two let

\[
 A=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad
 B=\begin{pmatrix}1&\rho\\ \rho&1\end{pmatrix},
 \quad 0<|\rho|<1,
\]

and

\[
                         R(t)=I+tA-t^2B.              \tag{2.6}
\]

On a sufficiently small interval around zero, \(R(t)\succ0\), while
\(R''=-2B\preceq0\).  Since \([A,B]\ne0\), the eigenspaces rotate with
\(t\), and covariance matrices at distinct times generally do not
commute.  Thus the theorem genuinely covers rotating nonlinear paths.

There is, however, no pure volume-preserving rotation.  If the eigenvalues
of \(R(t)\) are all constant, then \(\operatorname {tr}R\) is constant.
Matrix concavity gives \(R''\preceq0\), and
\(\operatorname {tr}R''=0\), hence \(R''=0\).  The path is affine.
Constancy of \(\operatorname {tr}R^2\) then forces its affine velocity to
vanish.  Therefore a nonconstant rotation must pay by changing the
eigenvalues as well.

## 3. Whitening and the exact third-moment block

Let

\[
 m=\mathbb ET,\qquad \sigma^2=\operatorname {Var}T,\qquad
 \bar R=\mathbb ER(T).
\]

Since \(R\) is not assumed affine, \(\bar R\) need not equal \(R(m)\).
Define

\[
 Z=\frac{T-m}{\sigma},\qquad
 U=\bar R^{-1/2}Y,\qquad
 C(z)=\bar R^{-1/2}R(m+\sigma z)\bar R^{-1/2}.         \tag{3.1}
\]

Then

\[
 \mathbb EZ=0,\quad\mathbb EZ^2=1,\quad
 \mathbb EU=0,\quad\mathbb E UU^T=I,\quad
 \mathbb E[ZU]=0,                                     \tag{3.2}
\]

so \(X=(Z,U)\) is isotropic.  Moreover,

\[
 C(z)=\mathbb E[UU^T\mid Z=z],\qquad
 \mathbb EC(Z)=I,                                     \tag{3.3}
\]

and \(C\) is a positive matrix-concave function.

Put

\[
 \kappa=\mathbb EZ^3,\qquad D=\mathbb E[ZC(Z)].        \tag{3.4}
\]

Conditional Gaussian symmetry gives, for \(u=(a,v)\),
\(a^2+|v|^2=1\),

\[
 M_{(a,v)}=
 \begin{pmatrix}
 a\kappa &(Dv)^T\\
 Dv&aD
 \end{pmatrix},                                      \tag{3.5}
\]

and consequently

\[
 \sup_{|u|=1}\|M_u\|_{\rm HS}^2
 =\max\left\{\kappa^2+\|D\|_{\rm HS}^2,
                    2\|D\|_{\rm op}^2\right\}.        \tag{3.6}
\]

Thus it remains only to bound \(\|D\|_{\rm HS}\).

## 4. Three one-dimensional lemmas

Let \(\nu(dz)=f(z)dz\) be centered, variance one, and log-concave on an
interval \(J\).

### 4.1 Isotropic core

There are numerical \(b,c_0,C_0>0\) such that

\[
 [-b,b]\subset\operatorname {int}J,\qquad
                c_0\le f(z)\le C_0\quad(|z|\le b).    \tag{4.1}
\]

This follows from the elementary one-dimensional height bound, Grunbaum's
inequality on both sides of the mean, Chebyshev's inequality, and
log-concavity.  One may take \(b=1/20\) after changing the other numerical
constants.

### 4.2 Canonical Stein kernel

Define

\[
 \tau(z)=\frac1{f(z)}\int_z^{\sup J}x f(x)\,dx
        =-\frac1{f(z)}\int_{\inf J}^z x f(x)\,dx.      \tag{4.2}
\]

Then \(\tau\ge0\), \(\mathbb E\tau=1\), and

\[
                   \mathbb E[Zg(Z)]=\mathbb E[\tau(Z)g'(Z)] \tag{4.3}
\]

for every locally absolutely continuous \(g\) for which the two sides are
integrable.  Standard one-dimensional log-concave tail and hazard bounds
also give

\[
                  \tau(z)\le C(1+|z|),\qquad
 \mathbb E\big[\tau(Z)(1+|Z|)^2\big]\le C.             \tag{4.4}
\]

For completeness, on \(z\ge0\), write

\[
 \int_z^\infty x f(x)dx
 =z\,\mathbb P(Z\ge z)+\int_z^\infty\mathbb P(Z\ge s)ds.
\]

The increasing-hazard property of a log-concave density bounds both tail
terms divided by \(f(z)\) by a numerical multiple of \(1+z\).  The
left-hand estimate is identical using the reverse hazard.  The last part
of (4.4) follows from the exponential tails of an isotropic one-dimensional
log-concave law.

### 4.3 Growth of a normalized matrix-concave function

Suppose \(C:J\to\mathbb S_+^d\) is matrix concave and
\(\mathbb EC(Z)=I\).  Then

\[
                         C(z)\preceq L(1+|z|)I         \tag{4.5}
\]

for a numerical \(L\).

**Proof.**  Matrix Jensen gives \(C(0)\succeq\mathbb EC=I\).  For
\(|z|\le b/2\), concavity between \(0\) and either \(b\) or \(-b\), and
positivity at the latter point, gives

\[
                         C(z)\succeq\frac12C(0).
\]

Integrating over this core and using (4.1) and \(\mathbb EC=I\) yields
\(C(0)\preceq L_0I\).

For a unit vector \(v\), put \(h(z)=v^TC(z)v\).  This is nonnegative and
concave.  For \(z>0\), if \(h(z)>h(0)\), concavity gives

\[
 {h(z)-h(0)\over z}\le {h(0)-h(-b)\over b}
                    \le {h(0)\over b}\le {L_0\over b}.
\]

If \(h(z)\le h(0)\), the desired estimate is immediate.  The argument for
\(z<0\) is symmetric.  Taking the supremum over \(v\) proves (4.5).

## 5. Marginal curvature pays for noncommuting velocity

The marginal density of \(Z\) is

\[
 f(z)=\text{const}\cdot
       e^{-\widehat W(z)}\det C(z)^{1/2},             \tag{5.1}
\]

where \(\widehat W(z)=W(m+\sigma z)\) up to an irrelevant affine
normalization.  Its convex potential is

\[
                  \Phi(z)=\widehat W(z)-\frac12\log\det C(z)+\text{const}.
                                                                    \tag{5.2}
\]

Define the affine-invariant covariance velocity

\[
                  A(z)=C(z)^{-1/2}C'(z)C(z)^{-1/2}.   \tag{5.3}
\]

Direct differentiation, with no commutativity assumption, gives

\[
\begin{aligned}
 \Phi''(z)
 &=\widehat W''(z)
   +\frac12\operatorname {tr}\big(A(z)^2\big)
   +\frac12\operatorname {tr}\big(C(z)^{-1}[-C''(z)]\big)\\
 &\ge\frac12\|A(z)\|_{\rm HS}^2.                     \tag{5.4}
\end{aligned}
\]

The last term is nonnegative because \(C''\preceq0\).  Formula (5.4) is
the exact noncommuting log-volume curvature identity.

Apply the Stein identity to \(g=\Phi'\).  If the density vanishes at the
ends of its support,

\[
 \mathbb E[\tau\Phi'']=\mathbb E[Z\Phi']=1.           \tag{5.5}
\]

For hard finite endpoints, integration by parts leaves the nonnegative
boundary term

\[
 [zf(z)]_{\inf J}^{\sup J},
\]

so the left side of (5.5) is at most one.  Convex approximation gives in
all cases

\[
              \mathbb E[\tau(Z)\|A(Z)\|_{\rm HS}^2]\le2. \tag{5.6}
\]

On the other hand, applying (4.3) entrywise to \(C\) yields

\[
                         D=\mathbb E[\tau(Z)C'(Z)].    \tag{5.7}
\]

Since

\[
 C'=C^{1/2}AC^{1/2},\qquad
 \|C'\|_{\rm HS}\le\|C\|_{\rm op}\|A\|_{\rm HS},
\]

Cauchy--Schwarz, (4.5), (4.4), and (5.6) give

\[
\begin{aligned}
 \|D\|_{\rm HS}
 &\le \mathbb E[\tau\|C'\|_{\rm HS}]\\
 &\le
 \left(\mathbb E[\tau\|A\|_{\rm HS}^2]\right)^{1/2}
 \left(\mathbb E[\tau\|C\|_{\rm op}^2]\right)^{1/2}\\
 &\le C.                                               \tag{5.8}
\end{aligned}
\]

Finally, \(|\kappa|\le C\) by the one-dimensional log-concave moment
bound.  Substitution in (3.6) proves (1.2).

### What prevents a harmonic spectrum

At each time, (5.4) charges

\[
 \sum_j s_j(z)^2,
\]

where \(s_j(z)\) are the eigenvalues of the relative velocity \(A(z)\).
Equation (5.6) says that this charge has total Stein-weighted budget at
most two.  Thus noncommuting rotations cannot hide a borderline
\(j^{-1/2}\) velocity spectrum: Frobenius energy, rather than only its
trace, appears in the exact determinant identity.  Matrix concavity then
prevents large covariance from amplifying a small relative velocity in
the tails, through (4.5).

## 6. Conditional means: a necessary quotient

The centered-slice assumption is essential in its literal form.  If

\[
             \mathcal V(t,y)=V_0(y-tb)
\]

on a bounded \(t\)-interval, then \(\mathcal V\) is jointly convex and all
slice partition functions are equal.  Marginal curvature is zero, while
the raw conditional second moment changes because the mean is \(tb\).
After the affine shear \((t,y)\mapsto(t,y-tb)\), the covariance is constant
and the apparent change disappears.

Thus any extension to noncentered slices must quotient out affine
translations and work with conditional covariance, not raw second moment.
Nonlinear conditional means would require an additional estimate and are
not covered by Theorem 1.

## 7. Test of the Brascamp--Lieb-deficit extension

The Gaussian proof suggests a possible extension to general log-concave
slices.  This section identifies exactly what is needed and why the usual
Brascamp--Lieb inequality alone does not provide it.

Let

\[
                         p_t(y)=Z_t^{-1}e^{-V(t,y)}
\]

be a smooth strictly log-concave conditional law.  Put
\(H=\nabla_y^2V\), \(g=V_t-\mathbb E_tV_t\), and let
\(\Phi(t)=-\log Z_t\).  Differentiation gives

\[
 \Phi''(t)=\mathbb E_tV_{tt}-\operatorname {Var}_t(V_t). \tag{7.1}
\]

Joint convexity and Brascamp--Lieb split this into two nonnegative terms:

\[
\begin{aligned}
 \Phi''
={}&\mathbb E_t\left[
 V_{tt}-\langle H^{-1}\nabla_yV_t,\nabla_yV_t\rangle\right]\\
 &+\underbrace{\left\{
 \mathbb E_t\langle H^{-1}\nabla g,\nabla g\rangle
 -\operatorname {Var}_t(g)\right\}}_{\mathcal D_{BL}(g)}.
                                                               \tag{7.2}
\end{aligned}
\]

The first line is the Schur-complement slack; the second is the
Brascamp--Lieb deficit.

If the conditional law is centered with covariance \(C(t)\), then

\[
 C'(t)=-\mathbb E_t\left[
 ((Y-\mathbb E_tY)(Y-\mathbb E_tY)^T-C(t))\,g(Y)\right]. \tag{7.3}
\]

Therefore the direct non-Gaussian analogue of (5.4) would be the stability
estimate

\[
\boxed{\quad
 \left\|C^{-1/2}
   \operatorname {Cov}_t\!\left(g,(Y-m)(Y-m)^T\right)
   C^{-1/2}\right\|_{\rm HS}^2
 \le C\,\mathcal D_{BL}(g).
 \quad}                                                \tag{7.4}
\]

For a standard Gaussian, (7.4) holds sharply with constant \(2\).  Indeed,
the covariance on the left extracts the second Hermite chaos, while
\(\mathcal D_{BL}\) multiplies chaos of degree \(k\) by \(k-1\).

### Hörmander reduction of the missing estimate

Normalize a fixed conditional law to be centered and isotropic.  Let

\[
                  L=-\Delta+\nabla V\cdot\nabla
\]

and solve \(Lu=g\), with \(\mathbb Eu=0\).  The integrated Bochner identity
is

\[
 \mathbb Eg^2
 =\mathbb E\|D^2u\|_{\rm HS}^2+
  \mathbb E\langle H\nabla u,\nabla u\rangle.          \tag{7.5}
\]

Also,

\[
 \mathbb Eg^2=\mathbb E\langle\nabla g,\nabla u\rangle.
\]

Cauchy--Schwarz in the \(H^{-1},H\) metrics, followed by (7.5), gives the
quantitative deficit bound

\[
                         \mathcal D_{BL}(g)
 \ge\mathbb E\|D^2u\|_{\rm HS}^2.                     \tag{7.6}
\]

On the other hand, integration by parts gives

\[
\begin{aligned}
 \operatorname {Cov}(g,YY^T-I)
 &=\mathbb E[g(YY^T-I)]\\
 &=\mathbb E\left[Y\otimes\nabla u+
                  \nabla u\otimes Y\right].           \tag{7.7}
\end{aligned}
\]

Consequently (7.4) would follow from the dimension-free affine-Hessian
rigidity estimate

\[
\boxed{\quad
 \left\|\mathbb E\left[
 Y\otimes\nabla u+\nabla u\otimes Y\right]\right\|_{\rm HS}^2
 \le C\,\mathbb E\|D^2u\|_{\rm HS}^2.
 \quad}                                                \tag{7.8}
\]

Estimate (7.8) is true with the sharp order for Gaussians, for quadratic
\(u\), and by one-dimensional tensorization for product measures.  It is
not a consequence of ordinary Bessel alone: Bessel gives

\[
 \|\mathbb E[Y\otimes\nabla u]\|_{\rm HS}^2
 \le\sum_j\operatorname {Var}(\partial_j u),
\]

and bounding the right side by the Hessian energy for a general
log-concave law is a Poincare step.  No dimension-free proof of (7.8) for
arbitrary log-concave measures is supplied here.

There is a second general-slice issue.  Matrix concavity of \(C(t)\), used
in the growth estimate (4.5), follows exactly in the Gaussian quadratic
case but is not automatic for arbitrary log-concave slices.  Even a proof
of (7.4) would therefore need either a replacement for (4.5) or a direct
Stein-weighted estimate of the conditional covariance growth.

Thus the Brascamp--Lieb decomposition is informative but does not, by
itself, extend the theorem to arbitrary non-Gaussian slices.  It isolates
two formal load-bearing statements:

1. the affine-Hessian rigidity (7.8), which converts BL deficit into
   covariance velocity; and
2. a dimension-free growth bound for the conditional covariance path.

Assuming either without proof would conceal the same sort of
high-dimensional rigidity sought in the main KLS problem.

## 8. Audit and scope

* All matrix identities above are noncommutative; no simultaneous
  diagonalization is used.
* The only spectral inequality is the one-dimensional Stein identity.
* Hard axial endpoints and nonsmooth concave covariance paths follow by
  convex approximation; the boundary term only decreases the curvature
  budget in (5.6).
* Degenerate transverse support is handled on its affine hull before
  whitening.
* The theorem proves the directional third-moment gate for the broadest
  centered conditional-Gaussian class compatible with exact joint
  log-concavity.
* The non-Gaussian extension remains conditional on the explicit new
  estimates (7.8) and the covariance-growth replacement; neither is
  treated as established.

## 9. Universal iid symmetrization

There is a canonical way to turn the general directional third moment
into a conditional-covariance regression.  Let \(X,Y\) be independent
copies of a centered isotropic log-concave vector, and put

\[
                  S={X+Y\over\sqrt2},\qquad
                  D={X-Y\over\sqrt2}.                 \tag{9.1}
\]

The orthogonal image \((S,D)\) is log-concave, both marginals are
isotropic, and \(\mathbb E[S D^T]=0\).  Exchangeability under
\(X\leftrightarrow Y\) makes the conditional law of \(D\) given \(S=s\)
centrally symmetric.  Hence

\[
 \mathbb E[D\mid S]=0,\qquad
 C(s):=\operatorname {Cov}(D\mid S=s),\qquad
 \mathbb EC(S)=I.                                    \tag{9.2}
\]

The last identity is also the law of total covariance.  A direct expansion
using independence gives the exact identity

\[
\boxed{\quad
 M_u(X)=\sqrt2\,\mathbb E\left[
             (u\cdot S)(C(S)-I)\right].
\quad}                                                \tag{9.3}
\]

Indeed,

\[
\begin{aligned}
 \mathbb E[(u\cdot S)C(S)]
 &=\mathbb E[(u\cdot S)DD^T]\\
 &=\frac1{2\sqrt2}\mathbb E\left[
 (u\cdot(X+Y))(X-Y)(X-Y)^T\right]\\
 &=\frac1{\sqrt2}M_u(X).
\end{aligned}
\]

All mixed terms vanish by centering and independence.  Formula (9.3) is
useful, but by itself it is exactly equivalent to the directional
third-moment gate rather than an estimate of it.

### 9.1 Conditional Prékopa formula

If \(X\) has a smooth density \(e^{-V}\), the joint potential in
sum--difference coordinates is

\[
 F(s,d)=V\left({s+d\over\sqrt2}\right)
       +V\left({s-d\over\sqrt2}\right).               \tag{9.4}
\]

Let \(\Psi(s)=-\log\int e^{-F(s,d)}\,dd\) be the potential of \(S\).
For \(H_\pm=\nabla^2V((s\pm d)/\sqrt2)\),

\[
 F_{ss}=F_{dd}={H_++H_-\over2},\qquad
 F_{sd}={H_+-H_-\over2}.                              \tag{9.5}
\]

For a direction \(u\), put
\[
 g_u(d)=\partial_uF(s,d)-\mathbb E_s\partial_uF(s,D).
\]
Then

\[
 \partial_{uu}\Psi(s)
 =\mathbb E_s[u^TF_{ss}u]-\operatorname {Var}_s(g_u). \tag{9.6}
\]

Brascamp--Lieb decomposes the right side into the Schur slack

\[
 \mathbb E_s\left[
 u^T(F_{ss}-F_{sd}F_{dd}^{-1}F_{ds})u\right]          \tag{9.7}
\]

and the BL deficit of \(g_u\).  The matrix in (9.7) is the
noncommutative parallel sum

\[
 F_{ss}-F_{sd}F_{dd}^{-1}F_{ds}
       =2(H_+^{-1}+H_-^{-1})^{-1}.                    \tag{9.8}
\]

The covariance derivative is

\[
 \partial_u C(s)
 =-\operatorname {Cov}_s(DD^T,g_u(D)).                \tag{9.9}
\]

Thus the same missing stability estimate (7.4) appears, now for a
centrally symmetric conditional law and an even score \(g_u\).  Neither
the nonnegative scalar curvature in (9.6) nor the law-of-total-covariance
identity (9.2) directly bounds the Frobenius norm in (9.9).

There are two possible ways to reduce (9.3) to one dimension, and each
loses the special structure:

1. applying a multivariate Stein identity to \(S\) requires an operator
   bound for a Stein kernel of \(S\), whereas the general known estimate
   is only a trace/Hilbert--Schmidt bound of order \(n\);
2. conditioning first on \(T=u\cdot S\) gives the one-dimensional
   regression
   \[
   \mathbb E[T(\bar C(T)-I)],\qquad
   \bar C(t)=\mathbb E[DD^T\mid u\cdot S=t],
   \]
   but \(\bar C(t)\) is not a matrix-concave Gaussian covariance path.

Consequently, Prékopa and Brascamp--Lieb positivity alone do not complete
the Frobenius estimate.

### 9.2 Smallest structural obstruction: one-dimensional uniforms

Take \(X\) uniform on \([-\sqrt3,\sqrt3]\).  Then \(S\) is supported on
\([-\sqrt6,\sqrt6]\), and, conditional on \(S=s\),

\[
               D\ \text{ is uniform on }\
 [-(\sqrt6-|s|),\,\sqrt6-|s|].
\]

Therefore

\[
                       C(s)={(\sqrt6-|s|)^2\over3}.   \tag{9.10}
\]

On either side of zero,

\[
                              C''(s)={2\over3}>0.      \tag{9.11}
\]

Thus the conditional covariance in the universal iid symmetrization is
not matrix concave even in dimension one.  This is the smallest explicit
obstruction to importing Theorem 1 directly.  The example has zero third
moment by symmetry, but the next example shows that the same phenomenon
can carry a nonzero regression.

### 9.3 Smallest nonzero example: centered exponentials

Let \(X=E-1\), where \(E\) is exponential of mean one.  This law is
centered, variance one, and log-concave.  If

\[
                         r=\sqrt2s+2,
\]

then, conditional on \(S=s\), the pair \((E_1,E_2)\) is uniform on the
segment \(E_1+E_2=r\).  Hence

\[
 D\mid S=s\ \text{ is uniform on }[-r/\sqrt2,r/\sqrt2],
\qquad
 C(s)={r^2\over6}.                                    \tag{9.12}
\]

Again \(C''=2/3>0\).  Since \(\mathbb EX^3=2\), (9.3) gives

\[
                         \mathbb E[S C(S)]=\sqrt2.    \tag{9.13}
\]

The marginal \(S\) is a shifted and rescaled Gamma\((2,1)\) law.  Its
potential satisfies

\[
 \Psi''(s)={2\over r^2}
 ={1\over4}\left({C'(s)\over C(s)}\right)^2.          \tag{9.14}
\]

Thus Prékopa curvature does pay for the scalar covariance velocity in this
example, but through moving-support geometry rather than covariance
concavity.

This also exposes a boundary failure of a naive smooth BL computation.
Inside its support, the exponential potential is affine and the
conditional potential in \(d\) is constant.  The interior Hessian and
interior BL deficit vanish, while \(C'(s)\ne0\).  All of the curvature is
carried by the moving endpoints.  Smooth approximations concentrate the
missing deficit in a boundary layer.  Therefore any non-Gaussian
Brascamp--Lieb proof must retain Reilly/boundary terms uniformly; simply
inverting the interior Hessian cannot be correct for cubes, exponentials,
or their limits.

### 9.4 Conditional-variance identities alone allow a harmonic spectrum

To see why (9.2), conditional symmetry, and one-dimensional concentration
are algebraically insufficient, let \(S\sim N(0,1)\), let

\[
 B_n=\operatorname {diag}(1,2^{-1/2},\ldots,n^{-1/2}),
\qquad
 C_n(s)=I+\frac12\tanh(s)B_n,                         \tag{9.15}
\]

and conditionally let \(D\mid S=s\) be centered Gaussian with covariance
\(C_n(s)\).  Then \(C_n(s)\succeq\frac12I\),
\(\mathbb EC_n(S)=I\), \(D\mid S\) is symmetric, and
\(\operatorname {Cov}D=I\).  Nevertheless,

\[
\begin{aligned}
 \left\|\mathbb E[S(C_n(S)-I)]\right\|_{\rm HS}^2
 &=\frac14\big(\mathbb E[S\tanh S]\big)^2
   \sum_{j=1}^n{1\over j}\\
 &\asymp\log n.                                       \tag{9.16}
\end{aligned}
\]

This model is not asserted to be an iid sum--difference law or jointly
log-concave; in fact the curvature theorem above explains why a
conditional-Gaussian realization of this harmonic profile cannot be
jointly log-concave.  Its purpose is precise: no manipulation using only
conditional symmetry, \(\mathbb EC=I\), and total covariance can prove the
desired bound.  A successful iid-symmetrization proof must exploit the
full Prékopa/BL geometry, including its boundary part, and must establish
a genuinely matrix-valued stability estimate beyond the identities
themselves.
