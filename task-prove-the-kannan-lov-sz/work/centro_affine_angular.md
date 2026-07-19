# The centro-affine angular route: exact radial form and a conditioned counterexample

## 1. Result of the audit

Let \(K\subset\mathbb R^n\), \(n\geq 2\), contain the origin in its interior
and have \(C^2\) boundary with everywhere positive Gauss curvature.  Write

\[
 G(x)=\|x\|_K,\qquad \rho(\theta)=G(\theta)^{-1},\qquad
 d\nu_K(\theta)=\frac{G(\theta)^{-n}}{n|K|}\,d\sigma(\theta).
\]

Thus \(\nu_K\) is radial cone measure.  Define the radial centro-affine
metric

\[
 \mathsf C_G
 :=\frac{\nabla_S^2G+Gg_S}{G}
 =g_S+\nabla_S^2\log G
       +d\log G\otimes d\log G.                         \tag{1.1}
\]

The tensor is positive definite under the stated curvature hypothesis.  The
classical Brunn--Minkowski inequality, in its Hilbert--Colesanti form, gives
the exact inequality

\[
 \boxed{
 \operatorname{Var}_{\nu_K} f
 \leq \frac1{n-1}\int_{S^{n-1}}
     \big\langle \mathsf C_G^{-1}\nabla_Sf,\nabla_Sf\big\rangle
     \,d\nu_K .}                                           \tag{1.2}
\]

There is no suppressed dimensional constant in (1.2).  For the Euclidean
ball, \(G\) is constant on the sphere, \(\mathsf C_G=g_S\), and (1.2) is the
sharp spherical Poincare inequality.

A comparison

\[
 \int\langle\mathsf C_G^{-1}\nabla f,\nabla f\rangle d\nu_K
 \leq A\int |\nabla f|^2d\nu_K                              \tag{1.3}
\]

with universal \(A\) would immediately imply

\[
 \operatorname{Var}_{\nu_K}f\leq \frac{A}{n-1}
       \int|\nabla_Sf|^2d\nu_K
\]

and hence \(D_{\rm ang}(\nu_K)\leq C/\sqrt n\).  The principal conclusion of
this note is that (1.3) is false even under the exact condition

\[
 \int\theta\otimes\theta\,d\nu_K(\theta)=\frac1n I,        \tag{1.4}
\]

and even when \(K\) is the Ball body of an isotropic log-concave probability.
In fact, for every even integer \(p\geq4\), set

\[
 G_{p,\varepsilon}(x)=\|x\|_p+\varepsilon\|x\|_2,
 \qquad K_{p,\varepsilon}=\{G_{p,\varepsilon}\leq1\}.
                                                               \tag{1.5}
\]

Then \(K_{p,\varepsilon}\) is smooth and strongly convex, (1.4) holds by
signed-permutation symmetry, but

\[
 \mathsf C_{G_{p,\varepsilon}}(e_1)|_{e_1^\perp}
       =\frac{\varepsilon}{1+\varepsilon}I.                 \tag{1.6}
\]

Consequently the best constant in (1.3) is at least
\((1+\varepsilon)/\varepsilon\), and is unbounded as
\(\varepsilon\downarrow0\).  A homothety makes the uniform law on
\(K_{p,\varepsilon}\) isotropic without changing either its angular law or
the tensor (1.1).  For a normalized uniform density, the Ball ray-integral
body is the body itself.  Thus (1.5) is a counterexample in precisely the
isotropic/Ball-body class relevant here.

The missing input is therefore not angular second-moment conditioning.  It is
a lower bound, or a genuinely nonlocal substitute for a lower bound, on the
centro-affine curvature tensor

\[
 \mathsf C_G=\frac{D^2G|_{TS^{n-1}}}{G}.                    \tag{1.7}
\]

Angular covariance is a zeroth-order integral constraint and does not control
(1.7).  This counterexample rules out the direct Hilbert-to-round Dirichlet
comparison; it does **not** disprove the desired round angular Poincare
inequality itself.

## 2. The three exact Brunn--Minkowski forms

This section fixes the constants and shows that the Colesanti, Hilbert, and
centro-affine formulations used below are the same statement.

### 2.1 Gauss-normal notation

Let \(u\in S^{n-1}\) denote an outer unit normal and let

\[
 h(u)=h_K(u),\qquad
 A(u)=\nabla_S^2h(u)+h(u)g_S.                                \tag{2.1}
\]

The inverse Gauss parametrization is

\[
 x(u)=\nabla_Sh(u)+h(u)u.
\]

The positive-definite endomorphism \(A\) is the reverse Weingarten map.  At
\(x(u)\), the Weingarten map is \(W=A^{-1}\), and

\[
 d\mathcal H^{n-1}_{\partial K}=\det A\,d\sigma(u),\qquad
 d\bar V_K(u)=\frac{h\det A}{n|K|}\,d\sigma(u).              \tag{2.2}
\]

Here \(\bar V_K\) is the cone-volume probability pushed to the normal sphere.

### 2.2 Colesanti's boundary form

For every \(\psi\in C^1(\partial K)\) satisfying
\(\int_{\partial K}\psi\,d\mathcal H^{n-1}=0\), the local
Brunn--Minkowski inequality is

\[
 \int_{\partial K}\operatorname{tr}(W)\psi^2\,d\mathcal H^{n-1}
 \leq
 \int_{\partial K}
 \langle W^{-1}\nabla_{\partial K}\psi,
                 \nabla_{\partial K}\psi\rangle
 \,d\mathcal H^{n-1}.                                      \tag{2.3}
\]

This is Colesanti's sharp inequality.  Equality is attained by infinitesimal
translations, \(\psi(x)=\langle\nu_K(x),a\rangle\).

Put \(\phi(u)=\psi(x(u))\).  Since
\(\nabla_S\phi=A\nabla_{\partial K}\psi\), (2.3) becomes

\[
 \int \operatorname{tr}(A^{-1})\phi^2\det A\,d\sigma
 \leq
 \int\langle A^{-1}\nabla_S\phi,\nabla_S\phi\rangle
            \det A\,d\sigma,                                \tag{2.4}
\]

under the condition \(\int\phi\det A\,d\sigma=0\).

### 2.3 From Colesanti to Hilbert

Set \(\phi=hz\).  Then the mean-zero condition in (2.4) is exactly

\[
 \int z\,d\bar V_K=0.                                       \tag{2.5}
\]

The following identity contains the entire conversion of (2.4):

\[
\begin{aligned}
 &\int\!\left[
  \langle A^{-1}\nabla(hz),\nabla(hz)\rangle
  -\operatorname{tr}(A^{-1})h^2z^2\right]\det A\,d\sigma \\
 &\quad=\int\!\left[
 h^2\langle A^{-1}\nabla z,\nabla z\rangle
 -(n-1)hz^2\right]\det A\,d\sigma .                         \tag{2.6}
\end{aligned}
\]

For completeness, let \(C^{ij}=\det(A)(A^{-1})^{ij}\) be the cofactor
matrix.  The Cheng--Yau divergence identity gives \(\nabla_jC^{ij}=0\).
After expanding the left side of (2.6), integrate

\[
 2hzC^{ij}h_i z_j=hC^{ij}h_i(z^2)_j
\]

by parts.  The resulting \(-C^{ij}h_ih_jz^2\) cancels the matching term,
and

\[
 C^{ij}h_{ij}=C^{ij}(A_{ij}-hg_{ij})
 =(n-1)\det A-h\operatorname{tr}(C).
\]

This proves (2.6).  Dividing the consequence of (2.4) by \(n|K|\) gives

\[
 \operatorname{Var}_{\bar V_K}z
 \leq\frac1{n-1}\int
 h\langle A^{-1}\nabla_Sz,\nabla_Sz\rangle\,d\bar V_K.    \tag{2.7}
\]

This is Hilbert's sharp spectral-gap form.  Equivalently, the
centro-affine metric in the Gauss-normal parametrization is

\[
 g_K^{\rm normal}=\frac{A}{h},                                \tag{2.8}
\]

and the quadratic form in (2.7) is its inverse metric.  The first nonzero
eigenvalue is exactly \(n-1\).  Its eigenspace is

\[
 z_a(u)=\frac{\langle u,a\rangle}{h(u)},\qquad a\in\mathbb R^n.
                                                               \tag{2.9}
\]

One may define the nonnegative Hilbert--Brunn--Minkowski operator without any
sign convention ambiguity by

\[
 \mathsf H_Kz
 =-\frac1{h\det A}\operatorname{div}_S
       \big(h^2\operatorname{cof}(A)\nabla_Sz\big).           \tag{2.10}
\]

Then

\[
 \int z\mathsf H_Kz\,d\bar V_K
 =\int h\langle A^{-1}\nabla z,\nabla z\rangle d\bar V_K,
 \qquad \lambda_1(\mathsf H_K)=n-1.                            \tag{2.11}
\]

## 3. Exact pullback from normal to radial coordinates

Let

\[
 x(\theta)=\rho(\theta)\theta,\qquad
 \varphi=\log\rho=-\log G,
 \qquad q=(1+|\nabla_S\varphi|^2)^{1/2}.                       \tag{3.1}
\]

For tangent vectors \(v,w\in T_\theta S^{n-1}\), direct differentiation of
the radial graph gives

\[
 u=\frac{\theta-\nabla_S\varphi}{q},\qquad
 h=\langle x,u\rangle=\frac{\rho}{q},                          \tag{3.2}
\]

and

\[
 \mathrm{II}_{\partial K}(dx(v),dx(w))
 =\frac{\rho}{q}\left[
 g_S+d\varphi\otimes d\varphi-\nabla_S^2\varphi
 \right](v,w).                                                 \tag{3.3}
\]

On \(\partial K\), the centro-affine metric is

\[
 g_K^{\partial K}=\frac{\mathrm{II}_{\partial K}}h.            \tag{3.4}
\]

Substituting \(\varphi=-\log G\) in (3.3)--(3.4) yields

\[
 (x(\cdot))^*g_K^{\partial K}
 =g_S+\nabla_S^2\log G+d\log G\otimes d\log G
 =\frac{\nabla_S^2G+Gg_S}{G}=\mathsf C_G.                      \tag{3.5}
\]

This also follows from the dual identity
\(g_K^{\rm radial}=D^2h_{K^\circ}/h_{K^\circ}\), since
\(h_{K^\circ}=G\).  Formula (3.5) shows explicitly the Hessian and the
rank-one gradient term which are both lost in a naive pullback of (2.8).

The radial surface Jacobian is

\[
 d\mathcal H^{n-1}_{\partial K}=\rho^{n-1}q\,d\sigma,
 \qquad h\,d\mathcal H^{n-1}_{\partial K}=\rho^n\,d\sigma.
                                                               \tag{3.6}
\]

Consequently cone measure becomes

\[
 d\nu_K(\theta)=\frac{\rho(\theta)^n}{n|K|}\,d\sigma
 =\frac{G(\theta)^{-n}}{n|K|}\,d\sigma.                        \tag{3.7}
\]

Pulling (2.7) through the radial map and using (3.5)--(3.7) proves (1.2).
The equality functions (2.9) become

\[
 z_a(\theta)=\langle\nabla_{\mathbb R^n}G(\theta),a\rangle,
                                                               \tag{3.8}
\]

because \(u/h=\nabla_{\mathbb R^n}G\).  In particular, except for the
Euclidean ball (or after the projective identification which turns an
ellipsoid into a ball), the first Hilbert modes are not the ordinary round
linear functions \(\langle\theta,a\rangle\).
Angular covariance controls the latter and says nothing directly about
(3.8).

## 4. What comparison would be sufficient, and what it really requires

Let

\[
 \mathcal E_{\rm ca}(f)
 =\int\langle\mathsf C_G^{-1}\nabla f,\nabla f\rangle d\nu_K,
 \qquad
 \mathcal E_{\rm rd}(f)=\int|\nabla f|^2d\nu_K.                \tag{4.1}
\]

The Hilbert inequality says

\[
 (n-1)\operatorname{Var}_{\nu_K}f\leq\mathcal E_{\rm ca}(f).
                                                               \tag{4.2}
\]

Therefore an upper comparison
\(\mathcal E_{\rm ca}\leq A\mathcal E_{\rm rd}\) is the useful direction.
For all smooth tests, such a comparison is equivalent to

\[
 \mathsf C_G^{-1}\leq A g_S^{-1}
 \quad\nu_K\text{-a.e.},
 \qquad\text{or equivalently}\qquad
 \mathsf C_G\geq A^{-1}g_S.                                   \tag{4.3}
\]

The implication from the form inequality to the pointwise inequality follows
by localization (or by localized high-frequency tests for the principal
symbol).  Thus a direct comparison is, exactly, a uniform lower
centro-affine-curvature estimate.

In contrast, angular conditioning is only

\[
 \frac cn I\leq\int\theta\otimes\theta\,d\nu_K
                    \leq\frac CnI.                              \tag{4.4}
\]

It contains no derivatives of \(G\).  The family in the next section shows
that even equality in (4.4) does not imply any positive constant in (4.3).

## 5. Smooth angularly isotropic counterexample

### 5.1 Geometry and angular covariance

Fix an even \(p\geq4\) and \(\varepsilon>0\), and use (1.5).  Both summands
of \(G_{p,\varepsilon}\) are convex and one-homogeneous.  On
\(\mathbb R^n\setminus\{0\}\) the function is smooth.  Moreover, for
\(v\perp\theta\),

\[
 D^2\|\cdot\|_2(\theta)[v,v]=|v|^2,
\]

so

\[
 D^2G_{p,\varepsilon}(\theta)[v,v]\geq\varepsilon|v|^2.
                                                               \tag{5.1}
\]

Hence the boundary of \(K_{p,\varepsilon}\) is smooth with positive Gauss
curvature.

The density of radial cone measure is proportional to
\(G_{p,\varepsilon}(\theta)^{-n}\).  It is invariant under every signed
permutation.  Sign changes kill the off-diagonal entries of its angular
second-moment matrix, permutations make all diagonal entries equal, and the
trace is one.  Therefore

\[
 \int\theta\otimes\theta\,d\nu_{p,\varepsilon}=\frac1nI.      \tag{5.2}
\]

### 5.2 Degeneration of the centro-affine tensor

For a one-homogeneous \(C^2\) function,

\[
 (\nabla_S^2G+Gg_S)(v,v)=D^2_{\mathbb R^n}G(v,v),
 \qquad v\perp\theta.                                         \tag{5.3}
\]

At \(\theta=e_1\), all transverse second derivatives of the \(\ell_p\) norm
vanish when \(p>2\), while those of the Euclidean norm equal one.  Since
\(G_{p,\varepsilon}(e_1)=1+\varepsilon\), (5.3) gives the exact formula
(1.6).  If (1.3) held with constant \(A\), localization at \(e_1\), with
gradient in the \(e_2\) direction, would give

\[
 A\geq\frac{1+\varepsilon}{\varepsilon}.                       \tag{5.4}
\]

Thus no constant depending only on the angular covariance constants in
(4.4) exists.

There is also a failure of the weaker estimate one might try only for round
one-Lipschitz functions.  Take \(n=2,p=4\), write
\(\theta=(\cos t,\sin t)\), and let \(f(\theta)=\theta_2=\sin t\).  On the
circle,

\[
 G(t)=(\cos^4t+\sin^4t)^{1/4}+\varepsilon,
 \qquad \mathsf C_G(t)=\frac{G''(t)+G(t)}{G(t)}.                \tag{5.5}
\]

Taylor expansion at \(t=0\) gives

\[
 (\cos^4t+\sin^4t)^{1/4}
 =1-\frac12t^2+\frac7{24}t^4+O(t^6),
\]

and hence, on a fixed small interval,

\[
 \mathsf C_G(t)\leq C_0(\varepsilon+t^2).                      \tag{5.6}
\]

The cone density is bounded below there uniformly in \(\varepsilon\), and
\(|f'(t)|=|\cos t|\geq1/2\).  It follows that

\[
 \mathcal E_{\rm ca}(f)
 \geq c_0\int_{-t_0}^{t_0}\frac{dt}{\varepsilon+t^2}
 \geq\frac{c_1}{\sqrt\varepsilon},                             \tag{5.7}
\]

whereas \(f\) is round one-Lipschitz and
\(\mathcal E_{\rm rd}(f)\leq1\).  Thus even a universal bound on the
Hilbert energy of round one-Lipschitz tests is false under exact angular
isotropy.

### 5.3 Realization by an isotropic log-concave law

Let \(Y\) be uniform on \(K_{p,\varepsilon}\).  Signed-permutation symmetry
gives

\[
 \mathbb EY=0,\qquad \operatorname{Cov}(Y)=a_{p,\varepsilon}I
\]

for some \(a_{p,\varepsilon}>0\).  Then
\(X=Y/\sqrt{a_{p,\varepsilon}}\) is isotropic and log-concave.  A homothety
does not change \(\Theta=X/|X|\), and it multiplies \(G\) by a constant, which
does not change \(\mathsf C_G\).

Finally, if \(f=|K|^{-1}1_K\), then \(f(0)=|K|^{-1}\) and, on every ray,

\[
 \frac n{f(0)}\int_0^\infty r^{n-1}f(r\theta)\,dr
 =n|K|\frac{\rho_K(\theta)^n}{n|K|}
 =\rho_K(\theta)^n.                                            \tag{5.8}
\]

So Ball's body \(K_n(f)\) is exactly \(K\).  Equations (5.2), (5.4), and
(5.8) prove that the counterexample lies in the required paired class.

## 6. The requested model tests

### 6.1 Cube

For the cube, up to scale,

\[
 G(\theta)=\|\theta\|_\infty.
\]

On the interior of every signed max-coordinate chamber, \(G\) is the
restriction of a linear functional.  Therefore

\[
 D^2_{\mathbb R^n}G|_{T_\theta S^{n-1}}=0,
 \qquad \mathsf C_G=0                                         \tag{6.1}
\]

there.  All curvature of the gauge is singular and sits on the chamber
walls.  The smooth gauges \(\|x\|_p+\varepsilon\|x\|_2\), with
\(p\to\infty\) and \(\varepsilon\downarrow0\), preserve exact angular
isotropy and make (6.1) a smooth degeneration.  Thus the cube is the basic
obstruction to a lower bound for \(\mathsf C_G\).

This does not contradict the good round angular behavior of the cube.  It
says that the Hilbert form becomes prohibitively large on gradients lying in
nearly flat centro-affine regions, so Hilbert's inequality cannot be converted
to the round one by form domination.

### 6.2 Regular simplex

Let \(u_0,\ldots,u_n\in S^{n-1}\) be the vertices of a centered regular
simplex:

\[
 \sum_{i=0}^n u_i=0,\qquad
 \langle u_i,u_j\rangle=-\frac1n\quad(i\ne j).
\]

For \(p>2\), define

\[
 H_p(x)=\left(\sum_{i=0}^n\langle u_i,x\rangle_+^p\right)^{1/p},
 \qquad G_{p,\varepsilon}^{\Delta}(x)=H_p(x)+\varepsilon|x|.
                                                               \tag{6.2}
\]

The function \(H_p\) is positive away from the origin, convex, and
one-homogeneous.  Taking, for example, an integer \(p\geq4\), (6.2) is
\(C^2\) away from the origin; the Euclidean summand makes its transverse
Hessian positive definite.  The simplex symmetry group acts irreducibly, so
both cone angular covariance and the covariance of the uniform law are scalar
matrices.  The former is exactly \(I/n\) because its trace is one; a homothety
makes the latter exactly \(I\).

At \(x=u_0\), the quantities \(\langle u_i,x\rangle\), \(i\ne0\), are
negative.  Hence on a whole neighborhood of \(u_0\),
\(H_p(x)=\langle u_0,x\rangle\) is linear.  Consequently

\[
 \mathsf C_{G_{p,\varepsilon}^{\Delta}}(u_0)|_{u_0^\perp}
 =\frac{\varepsilon}{1+\varepsilon}I.                           \tag{6.3}
\]

As \(p\to\infty\), \(H_p\) tends to
\(\max_i\langle u_i,x\rangle\), the gauge of a regular simplex (up to
polarity and scale).  Thus the simplex gives the same conditioned obstruction
as the cube: flat chambers, with curvature concentrated on their interfaces.

### 6.3 Ellipsoids

Let \(K=A B_2^n\) and put \(B=A^{-T}A^{-1}\).  Then

\[
 G(\theta)=(\theta^TB\theta)^{1/2}.
\]

For \(v\perp\theta\), direct differentiation gives

\[
 \mathsf C_G(\theta)[v,v]
 =\frac{v^TBv}{\theta^TB\theta}
  -\frac{(v^TB\theta)^2}{(\theta^TB\theta)^2}.                 \tag{6.4}
\]

Writing \(\lambda_-\) and \(\lambda_+\) for the extreme eigenvalues of
\(B\), one has

\[
 \frac{\lambda_-}{\lambda_+}|v|^2
 \leq\mathsf C_G(\theta)[v,v]
 \leq\frac{\lambda_+}{\lambda_-}|v|^2.                         \tag{6.5}
\]

Indeed, the bracket before division by \(\theta^TB\theta\) in (6.4) is

\[
 \min_{a\in\mathbb R}(v-a\theta)^TB(v-a\theta),
\]

which lies between \(\lambda_-|v|^2\) and \(\lambda_+|v|^2\) because
\(v\perp\theta\).  Thus an ellipsoid of bounded eccentricity has a bounded
comparison, with constant at most
\(\lambda_+/\lambda_-=\operatorname{cond}(A)^2\).  A centered uniform
ellipsoid which is isotropic has \(A\) scalar, so \(\mathsf C_G=g_S\) after
scale.  Ellipsoids are therefore a positive calibration, not a counterexample
after isotropic normalization.

### 6.4 A smoothed high-curvature spike

The covariance condition also gives no upper control on centro-affine
curvature.  Fix \(1<q<2\) and \(\eta>0\), and define

\[
 r_i(x)=\big(x_i^2+\eta^2|x|^2\big)^{1/2},\qquad
 G_{q,\eta}^{\rm sp}(x)=\left(\sum_{i=1}^n r_i(x)^q\right)^{1/q}.
                                                               \tag{6.6}
\]

This is a smooth, one-homogeneous, signed-permutation invariant gauge.  Each
\(r_i\) is a smooth strictly convex quadratic norm.  Its Hessian is positive
on every direction not parallel to \(x\); since all first derivatives of the
outer \(\ell_q\) norm are positive, the transverse Hessian of (6.6) is
positive definite.  Thus its unit body is smooth and strongly convex.  As
\(\eta\downarrow0\) and \(q\downarrow1\), it approaches the pointy
cross-polytope model.

To calculate the curvature at a tip, put \(x(t)=e_1+te_2\) and
\(S(t)=G_{q,\eta}^{\rm sp}(x(t))^q\).  Then \(S'(0)=0\),

\[
\begin{aligned}
 S(0)&=(1+\eta^2)^{q/2}+(n-1)\eta^q,\\
 S''(0)&=q\eta^2(1+\eta^2)^{q/2-1}
       +q(1+\eta^2)\eta^{q-2}
       +q(n-2)\eta^q.
\end{aligned}                                                   \tag{6.7}
\]

Consequently

\[
 \mathsf C_{G_{q,\eta}^{\rm sp}}(e_1)[e_2,e_2]
 =\frac{S''(0)}{qS(0)}\asymp_{n,q}\eta^{q-2}
 \longrightarrow\infty.                                      \tag{6.8}
\]

Its angular covariance is nevertheless exactly \(I/n\), again by
signed-permutation symmetry, and its uniform law becomes isotropic by a
homothety.  Thus angular covariance controls neither the lower nor the upper
side of the centro-affine tensor.  For the one-sided comparison useful in
(1.3), large curvature is harmless because it makes
\(\mathsf C_G^{-1}\) small; the fatal missing input is the lower-curvature
side exhibited by (1.5) and (6.2).

## 7. Consequences for the angular KLS route

The exact centro-affine theorem supplies a universal spectral gap, but for a
body-dependent metric.  To obtain the round angular estimate one must compare
that metric with the round metric in the direction (4.3).  The examples above
prove:

1. The stable angular condition
   \(cI/n\leq\mathbb E\Theta\Theta^T\leq CI/n\) does not give such a
   comparison; exact equality is insufficient.
2. This remains true for smooth, strongly convex Ball bodies paired with
   isotropic log-concave laws.
3. In polyhedral limits, the centro-affine tensor vanishes in open chambers
   and the positive curvature is a singular measure on walls.  An inverse of
   the absolutely continuous tensor cannot encode the wall gluing needed for
   round mixing.
4. The rank-one term
   \(d\log G\otimes d\log G\) in (1.1) is an essential part of the exact
   radial pullback, but retaining it does not repair the comparison: the full
   tensor still degenerates in the conditioned examples.
5. Hilbert's first eigenspace consists of the adapted conormal functions
   \(\langle\nabla G,a\rangle\), whereas angular covariance controls the
   round functions \(\langle\theta,a\rangle\).  Bridging these two families
   also requires derivative/curvature information absent from covariance.

Therefore the Hilbert--Colesanti inequality does not independently prove
\(D_{\rm ang}(\nu)\leq C/\sqrt n\) by a dimension-free metric comparison.
A surviving centro-affine route would need a new nonlocal coercivity theorem
which uses the full curvature measure, including singular wall curvature, and
which remains stable as smooth curvature collapses onto walls.  On the paired
isotropic class, a theorem strong enough to yield the round gap at scale
\(n\) is itself a KLS-strength input; none has been proved in this note.

## 8. Primary references used for normalization

* A. Colesanti, *From the Brunn--Minkowski inequality to a class of Poincare
  type inequalities*, Comm. Contemp. Math. 10 (2008), 765--772,
  [arXiv:math/0703584](https://arxiv.org/abs/math/0703584).
* E. Milman, *Centro-affine differential geometry and the log-Minkowski
  problem*, J. Eur. Math. Soc. 27 (2025), 709--772,
  [doi:10.4171/JEMS/1386](https://doi.org/10.4171/JEMS/1386).  In particular,
  its normal metric \(D^2h_K/h_K\), radial metric
  \(D^2h_{K^\circ}/h_{K^\circ}\), cone-volume measure, and Hilbert spectral
  gap \(n-1\) agree with (2.7)--(3.7).
