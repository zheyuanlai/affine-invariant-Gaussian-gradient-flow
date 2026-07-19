# Mean-gradient coercivity: directional third moments, wedges, and cones

## 0. Verdict

Let \(X\) have an isotropic log-concave law \(\mu\) on
\(\mathbb R^n\).  For a smooth function \(g\), write

\[
 g\perp\mathrm {Aff}
 \quad\Longleftrightarrow\quad
 \mathbb Eg=0,\qquad \mathbb E[Xg(X)]=0.               \tag{0.1}
\]

The residual mean-gradient targets are

\[
 |\mathbb E\nabla g|
 \le C\big(\|g\|_2+\|D^2g\|_2\big)                    \tag{MG}
\]

and the stronger Hessian-only estimate

\[
 |\mathbb E\nabla g|
 \le C\|D^2g\|_2.                                     \tag{HMG}
\]

This report stress-tests both statements under exact affine normalization.
It does not prove either inequality for every isotropic log-concave law,
and it finds no growing counterexample.  The concrete conclusions are:

1. On quadratic functions, (HMG) is exactly the directional third-moment
   estimate
   \[
    \sup_{|u|=1}
    \left\|\mathbb E[(u\cdot X)(XX^T-I)]\right\|_{\rm HS}
    \le C.                                             \tag{0.2}
   \]
   This is the stochastic-localization tensor parameter
   \(\kappa_n\).  Projection thin-shell estimates alone give only
   \(O(\sqrt{\log n})\), so a universal proof is already a substantive
   tensor theorem.

2. A multiscale family of dependent polytopes, with half of its transverse
   widths increasing and half decreasing along an axial coordinate, is
   completely benign.  For arbitrary slopes \(a_i\), exact whitening
   gives
   \[
    \sup_{|u|=1}\|M_u\|_{\rm HS}\le2.                  \tag{0.3}
   \]
   The slice-volume potential has curvature
   \(2\sum_i a_i^2\); this forces the axial variance \(v\) to satisfy
   \(v\sum_i a_i^2\le1/2\).  It exactly cancels every attempted
   multiscale accumulation.

3. The same wedge polytope is a \(4\)-Lipschitz image of a product of
   isotropic one-dimensional log-concave laws.  Therefore its isotropic
   uniform law has
   \[
     C_P\le192,
   \]
   and every smooth \(g\perp\mathrm {Aff}\) satisfies the full
   Hessian-only estimate
   \[
     |\mathbb E\nabla g|\le192\|D^2g\|_2.              \tag{0.4}
   \]
   This is a genuinely dependent, non-centrally-symmetric polyhedral
   class, not an affine product in the displayed coordinates.

   A one-sided version with widths \(1+a_it\), \(0\le a_i\le1/2\), has no
   pair-exchange symmetry.  It still has a universal directional
   third-moment bound, is a \((\sqrt{28}+3)\)-Lipschitz product image, and
   satisfies
   \(C_P\le12(\sqrt{28}+3)^2<972\), with the same constant in (HMG).

4. A pyramid over an arbitrary isotropic base has an exact third-tensor
   block decomposition.  The axial slice has Hilbert--Schmidt norm
   \(O(1)\), while transverse slices inherit the base tensor with a factor
   at most \(8/(5\sqrt2)<1.14\), plus universal low-rank terms.  Thus a
   single coning operation cannot create a growing tensor from a benign
   base.

5. The isotropic half-ball, a strongly non-symmetric curved example, also
   has \(\sup_u\|M_u\|_{\rm HS}=O(1)\).  Its one axial eigenvalue stays
   constant, while its \(n-1\) tangential eigenvalues are \(O(1/n)\).

6. Gaussian convolution of the wedge examples gives smooth, full-support,
   genuinely dependent log-concave laws.  Whitening only decreases their
   directional third tensor, and their Poincare constants remain
   universal, so (0.4) survives smoothing.

The attempted counterexample mechanism is now precise.  To make
\(\|M_u\|_{\rm HS}\) grow, one wants many conditional transverse variances
to change with one standardized coordinate.  In every explicit convex
slice family below, Brunn--Minkowski volume growth concentrates that
coordinate by exactly the amount needed to keep the squared changes
summable.  A genuine counterexample would have to vary many covariance
directions while paying substantially less determinant or slice-volume
curvature than these models.

## 1. Quadratic reduction and boundary bookkeeping

### 1.1 The exact third-moment diagnostic

Let \(B\in\mathrm {Sym}_n\) and put

\[
\begin{aligned}
 q_B(x)&=x^TBx-\operatorname {Tr}B,\\
 d_B&=\mathbb E[Xq_B(X)],\\
 g_B(x)&=q_B(x)-d_B\cdot x.
\end{aligned}                                         \tag{1.1}
\]

Isotropy gives \(g_B\perp\mathrm {Aff}\), and direct differentiation gives

\[
 \mathbb E\nabla g_B=-d_B,\qquad
 D^2g_B=2B.                                            \tag{1.2}
\]

For a unit vector \(u\), define

\[
 M_u=\mathbb E[(u\cdot X)(XX^T-I)].                    \tag{1.3}
\]

Then

\[
 u\cdot d_B
 =\mathbb E[(u\cdot X)(X^TBX-\operatorname {Tr}B)]
 =\langle M_u,B\rangle_{\rm HS}.                       \tag{1.4}
\]

Hilbert--Schmidt duality therefore gives

\[
 \sup_{B\ne0}{|d_B|\over\|B\|_{\rm HS}}
 =\sup_{|u|=1}\|M_u\|_{\rm HS}.                        \tag{1.5}
\]

Thus (HMG), even only on quadratics, implies (0.2), and (0.2) is exactly
the quadratic instance of (HMG), up to the factor two in (1.2).

The \(\|g\|_2\) term in (MG) changes the diagnostic.  Indeed,

\[
 \|g_B\|_2^2
 =\operatorname {Var}(X^TBX)-|d_B|^2.                 \tag{1.6}
\]

Consequently a growing \(M_u\) would disprove Hessian-only coercivity, but
it would not by itself disprove (MG): the quadratic residual in (1.6)
could pay for the same growth.

### 1.2 What thin shell proves, and no more

Let \(P\) be a rank-\(r\) orthogonal projection.  The dimension-free
thin-shell theorem applied to the marginal \(PX\) gives

\[
 \operatorname {Var}(|PX|^2)\le C r.                  \tag{1.7}
\]

Therefore

\[
 |\operatorname {Tr}(PM_u)|
 =|\mathbb E[(u\cdot X)(|PX|^2-r)]|
 \le C\sqrt r.                                        \tag{1.8}
\]

Applying (1.8) to the positive and negative spectral projections of
\(M_u\) gives, for the decreasing absolute eigenvalues,

\[
 |\lambda_j(M_u)|\le {C\over\sqrt j}.                 \tag{1.9}
\]

Hence

\[
 \|M_u\|_{\rm HS}\le C\sqrt{\log(en)}.                 \tag{1.10}
\]

The abstract matrix
\(\operatorname {diag}(1,2^{-1/2},\ldots,n^{-1/2})\)
satisfies all Ky Fan estimates (1.8) and has squared Hilbert--Schmidt norm
comparable to \(\log n\).  Removing the logarithm requires compatibility
among the eigendirections which projection thin shell does not encode.

### 1.3 Boundary and score forms

For a smooth full-support density \(d\mu=e^{-V}dx/Z\), integration by
parts gives

\[
 \mathbb E[\partial_u g]
 =\mathbb E[g\,\partial_uV].                           \tag{1.11}
\]

For the uniform law on a convex body \(K\),

\[
 \mathbb E_K[\partial_u g]
 ={1\over|K|}\int_{\partial K}g(x)\,u\cdot n(x)\,
 d\mathcal H^{n-1}(x).                                \tag{1.12}
\]

Thus the mean-gradient functional is carried entirely by the score or by
the boundary.  Thin wedges can make individual raw facet terms large.
Those terms are not affine invariant.  If a raw body has covariance
\(\Sigma\), the relevant body is
\(\Sigma^{-1/2}(K-\mathbb EX)\), and the normal slopes transform with the
same matrix.  Every example below is whitened explicitly before its
boundary or third-moment size is assessed.

## 2. A multiscale paired box wedge

### 2.1 Definition and exact whitening

Fix \(m\ge1\) and parameters

\[
 0\le a_i<1,\qquad
 S=\sum_{i=1}^m a_i^2.
\]

In dimension \(n=2m+1\), define the polytope

\[
\begin{aligned}
 K_a=\{(t,x,y):\;&-1\le t\le1,\\
 &|x_i|\le1+a_it,\quad
 |y_i|\le1-a_it\quad(1\le i\le m)\}.
\end{aligned}                                         \tag{2.1}
\]

It is an intersection of affine halfspaces and hence convex.  Conditional
on \(T=t\), all transverse coordinates are independent centered uniforms.
The slice volume is

\[
 |(K_a)_t|=4^m\prod_{i=1}^m(1-a_i^2t^2),              \tag{2.2}
\]

so the axial density is

\[
 \rho(t)=Z^{-1}\mathbf1_{[-1,1]}(t)
 \prod_{i=1}^m(1-a_i^2t^2).                           \tag{2.3}
\]

It is even.  Put

\[
 v=\mathbb ET^2.
\]

The interior potential

\[
 W(t)=-\sum_i\log(1-a_i^2t^2)
\]

satisfies

\[
 W''(t)
 =2\sum_i {a_i^2(1+a_i^2t^2)\over(1-a_i^2t^2)^2}
 \ge2S.                                                \tag{2.4}
\]

The one-dimensional Brascamp--Lieb inequality on the interval, obtained
also by smooth convex approximation of the hard walls, gives

\[
 \boxed{\quad vS\le{1\over2}.\quad}                   \tag{2.5}
\]

When \(S=0\), all formulas below are read by continuity and the law is a
product.

The raw covariance is diagonal.  The two coordinates in the \(i\)-th pair
have the same variance

\[
 \sigma_i^2
 ={1\over3}\mathbb E(1+a_iT)^2
 ={1+a_i^2v\over3}.                                   \tag{2.6}
\]

The exact isotropic coordinates are

\[
 U={T\over\sqrt v},\qquad
 X_i={x_i\over\sigma_i},\qquad
 Y_i={y_i\over\sigma_i}.                              \tag{2.7}
\]

All mixed covariances vanish by the independent sign symmetries and the
involution

\[
 (t,x,y)\longmapsto(-t,y,x).
\]

Thus (2.7), not a heuristic diagonal rescaling, is the full covariance
whitening.

### 2.2 Complete directional third-moment calculation

Define

\[
 \lambda_i
 ={2a_i\sqrt v\over1+a_i^2v}.                          \tag{2.8}
\]

Conditional second moments and the evenness of \(T\) give

\[
\begin{aligned}
 \mathbb E[UX_i^2]
 &={\mathbb E[T(1+a_iT)^2]\over
       \sqrt v(1+a_i^2v)}
 =\lambda_i,\\
 \mathbb E[UY_i^2]&=-\lambda_i.                       \tag{2.9}
\end{aligned}
\]

Every other third moment vanishes: \(U^3\) vanishes by the axial
involution, and every transverse coordinate has its own sign symmetry.

Let

\[
 u=\alpha e_0+\sum_i(p_i e_i+q_i f_i),
\qquad
 \alpha^2+\sum_i(p_i^2+q_i^2)=1,                      \tag{2.10}
\]

where \(e_0\) is axial and \(e_i,f_i\) are the two coordinates of pair
\(i\).  The only nonzero entries of \(M_u\) are

\[
\begin{array}{lll}
 (M_u)_{e_i e_i}=\alpha\lambda_i,&
 (M_u)_{f_i f_i}=-\alpha\lambda_i,\\
 (M_u)_{0e_i}=(M_u)_{e_i0}=p_i\lambda_i,&
 (M_u)_{0f_i}=(M_u)_{f_i0}=-q_i\lambda_i.
\end{array}                                           \tag{2.11}
\]

Consequently

\[
 \|M_u\|_{\rm HS}^2
 =2\alpha^2\sum_i\lambda_i^2
  +2\sum_i\lambda_i^2(p_i^2+q_i^2).                  \tag{2.12}
\]

Equations (2.5) and (2.8) yield

\[
 \sum_i\lambda_i^2
 \le4v\sum_i a_i^2
 \le2.                                                 \tag{2.13}
\]

Using (2.10) in (2.12),

\[
 \boxed{\quad
 \sup_{|u|=1}\|M_u\|_{\rm HS}^2\le4.
 \quad}                                                \tag{2.14}
\]

This permits arbitrary multiscale slopes.  For example, the attempted
harmonic choice \(a_i=i^{-1/2}\) has
\(\sum_i a_i^2\asymp\log m\), but (2.5) forces
\(v=O(1/\log m)\).  The normalized eigenvalues become
\(O((i\log m)^{-1/2})\), whose squares sum to \(O(1)\).

### 2.3 A universal Poincare bound for the whole dependent class

There is an exact product parametrization.  Let

* \(U=T/\sqrt v\), with the isotropic one-dimensional log-concave law
  induced by (2.3);
* \(R_i,Q_i\) be independent uniforms on
  \([-\sqrt3,\sqrt3]\), independent of \(U\).

Define

\[
\begin{aligned}
 F_0&=U,\\
 F_{e_i}
 &= {1+a_i\sqrt v\,U\over\sqrt{1+a_i^2v}}R_i,\\
 F_{f_i}
 &= {1-a_i\sqrt v\,U\over\sqrt{1+a_i^2v}}Q_i.
\end{aligned}                                         \tag{2.15}
\]

The law of \(F\) is exactly the uniform law on the whitened polytope
(2.7).  This follows from the slice density (2.3) and the conditional
uniform product structure.

Put

\[
 c_i={a_i\sqrt v\over\sqrt{1+a_i^2v}}.
\]

The first column of \(DF\) has squared norm

\[
 1+\sum_i c_i^2(R_i^2+Q_i^2)
 \le1+6vS\le4.                                        \tag{2.16}
\]

Every other column has norm at most two, because
\[
 0\le1\pm a_i\sqrt v\,U=1\pm a_iT\le2.
\]

Writing \(DF\) as its first column plus its transverse diagonal part
therefore gives the pointwise operator bound

\[
 \boxed{\quad\|DF\|_{\rm op}\le4.\quad}                \tag{2.17}
\]

Every isotropic one-dimensional log-concave probability has Poincare
constant at most \(12\), and the same holds for the isotropic interval.
Tensorization gives Poincare constant at most \(12\) for the product
input in (2.15).  Pulling a test function back through the
\(4\)-Lipschitz map gives

\[
 \boxed{\quad C_P(\lambda_{K_a}^{\rm iso})\le
 12\cdot4^2=192.\quad}                                 \tag{2.18}
\]

This proof uses no KLS input in growing dimension; it uses only the
one-dimensional log-concave Poincare bound.

Let \(g\perp\mathrm {Aff}\) be smooth on the whitened body and put

\[
 m=\mathbb E\nabla g.
\]

Apply (2.18) componentwise to the derivatives:

\[
 \mathbb E|\nabla g-m|^2
 \le192\,\mathbb E\|D^2g\|_{\rm HS}^2.                \tag{2.19}
\]

Apply (2.18) once more to \(g-m\cdot x\).  Its mean is zero, and

\[
\begin{aligned}
 \|g-m\cdot x\|_2^2
 &\le192\,\mathbb E|\nabla g-m|^2\\
 &\le192^2\,\mathbb E\|D^2g\|_{\rm HS}^2.
\end{aligned}                                         \tag{2.20}
\]

Affine orthogonality and isotropy give

\[
 \|g-m\cdot x\|_2^2=\|g\|_2^2+|m|^2.
\]

Therefore

\[
 \boxed{\quad
 |\mathbb E\nabla g|
 \le192\|D^2g\|_2.
 \quad}                                                \tag{2.21}
\]

Thus the full Hessian-only mean-gradient estimate holds on this dependent
polyhedral class, not merely on its quadratic subspace.

### 2.4 Boundary interpretation

Before whitening, the slanted facets have axial slopes \(a_i\).  A naive
sum of their absolute boundary contributions in (1.12) sees
\(\sum_i|a_i|\), which can grow like \(\sqrt m\) or faster.  In the
raw coordinates, the two \(x_i\)-facets have equations

\[
 \pm x_i-a_it=1
\]

and outward normals proportional to
\((-a_i,\pm e_i)\); the two \(y_i\)-facets have normals proportional to
\((a_i,\pm f_i)\).  In the divergence formula, the normalization of the
normal cancels the same factor in the graph surface element.  Thus the
apparently dangerous axial boundary coefficients really are the raw
numbers \(\pm a_i\).

After substituting \(t=\sqrt v\,U\) and \(x_i=\sigma_iX_i\), an
\(x_i\)-facet becomes

\[
 \pm X_i-{\sqrt3\,a_i\sqrt v\over
                 \sqrt{1+a_i^2v}}\,U
 ={1\over\sigma_i}.
\]

Hence the normalized width slopes are

\[
 \sqrt3\,c_i
 ={\sqrt3\,a_i\sqrt v\over\sqrt{1+a_i^2v}},
\]

and

\[
 \sum_i 2(\sqrt3\,c_i)^2
 \le6vS\le3.                                          \tag{2.22}
\]

Equation (2.22) is the boundary version of the third-moment cancellation.
The growing raw facet sum is an artifact of measuring normals before
affine normalization.

### 2.5 A genuinely one-sided wedge

The paired construction makes the algebra exact but has an involution
which exchanges its two transverse groups.  The same compensation survives
without that symmetry.  Let

\[
 0\le a_i\le{1\over2},\qquad S=\sum_i a_i^2,
\]

and define

\[
 K_a^+=\{(t,x):-1\le t\le1,\ |x_i|\le1+a_it\}.        \tag{2.23}
\]

The axial density is proportional to

\[
 \prod_i(1+a_it).
\]

Let \(\mu=\mathbb ET\), \(v=\operatorname {Var}T\), and
\(\xi=T-\mu\).  Since

\[
 {d^2\over dt^2}\left[-\sum_i\log(1+a_it)\right]
 =\sum_i{a_i^2\over(1+a_it)^2}
 \ge {4\over9}S,
\]

one-dimensional Brascamp--Lieb gives

\[
 vS\le {9\over4}.                                     \tag{2.24}
\]

Put

\[
 L_i=1+a_i\mu,\qquad
 D_i=L_i^2+a_i^2v.
\]

Here \(1/2\le L_i\le3/2\).  In isotropic coordinates

\[
 U={\xi\over\sqrt v},\qquad
 X_i={\sqrt3\,x_i\over\sqrt{D_i}},
\]

the only nonzero third moments, apart from \(\mathbb EU^3\), are

\[
\begin{aligned}
 \lambda_i
 :=\mathbb E[UX_i^2]
 ={2L_ia_iv+a_i^2\mathbb E\xi^3\over
        \sqrt v\,D_i}.                                \tag{2.25}
\end{aligned}
\]

Every centered variance-one one-dimensional log-concave variable has a
universal third absolute moment.  We may take the explicit (deliberately
nonoptimized) bound

\[
 \mathbb E|U|^3\le C_3
 :=8+{96\over[\log(4/e)]^3}.                           \tag{2.26}
\]

Indeed, Grünbaum's one-dimensional lemma gives
\(\mathbb P(U\ge0),\mathbb P(U\le0)\ge e^{-1}\), while Chebyshev gives
\(\mathbb P(U\ge2),\mathbb P(U\le-2)\le1/4\).  Each one-sided survival
function is log-concave by Prékopa--Leindler.  After dividing by its value
at zero, concavity of the logarithm therefore gives, for \(t\ge2\),

\[
 \mathbb P(|U|\ge t)
 \le2\exp\left[-{t\over2}\log(4/e)\right].
\]

Integrating \(3t^2\mathbb P(|U|\ge t)\) over \([0,2]\) and \([2,\infty)\),
and enlarging the latter integral to \([0,\infty)\), gives (2.26).

Using \(D_i\ge1/4\), \(L_i\le3/2\), and (2.26),

\[
 |\lambda_i|
 \le12a_i\sqrt v+4C_3a_i^2v.                          \tag{2.27}
\]

Since \(v\le1\), \(a_i^2\le1/4\), and (2.24),

\[
\begin{aligned}
 \sum_i\lambda_i^2
 &\le288vS+32C_3^2v^2\sum_i a_i^4\\
 &\le648+18C_3^2.                                     \tag{2.28}
\end{aligned}
\]

Also \(|\mathbb EU^3|\le C_3\).  Individual transverse sign symmetries
now give, for
\(u=\alpha e_0+\sum_i p_ie_i\),

\[
 \|M_u\|_{\rm HS}^2
 =\alpha^2\left((\mathbb EU^3)^2+\sum_i\lambda_i^2\right)
  +2\sum_i p_i^2\lambda_i^2.                          \tag{2.29}
\]

Equations (2.28)--(2.29) prove a universal directional third-moment bound
for this one-sided family.

There is again a product parametrization.  Let \(R_i\) be independent
uniforms on \([-\sqrt3,\sqrt3]\), independent of \(U\), and put

\[
 F_0=U,\qquad
 F_i={L_i+a_i\sqrt v\,U\over\sqrt{D_i}}R_i.           \tag{2.30}
\]

The first column of \(DF\) has squared norm at most

\[
 1+3\sum_i{a_i^2v\over D_i}
 \le1+12vS\le28,
\]

and every transverse multiplier is at most three.  Hence

\[
 \|DF\|_{\rm op}\le\sqrt{28}+3<9.                     \tag{2.31}
\]

The product input has Poincare constant at most \(12\), so

\[
 C_P((\lambda_{K_a^+})^{\rm iso})
 \le12(\sqrt{28}+3)^2
 =444+72\sqrt{28}<972.                                \tag{2.32}
\]

Repeating (2.19)--(2.21) proves

\[
 |\mathbb E\nabla g|
 \le(444+72\sqrt{28})\|D^2g\|_2
\]

for every smooth affine-orthogonal \(g\) on the isotropic one-sided
wedge.  Thus the positive result is not caused by the pair-exchange
symmetry of (2.1).

## 3. Pyramids inherit rather than create third-moment growth

Let \(W\in\mathbb R^{n-1}\) be centered and isotropic, and let \(S\) be
independent with density

\[
 n s^{n-1}\mathbf1_{[0,1]}(s).
\]

The uniform law on the pyramid with apex at the origin and base
\(\{1\}\times L\), where \(W\) is uniform on \(L\), is represented by

\[
 (S,SW).
\]

Put

\[
\begin{aligned}
 \mu_S&={n\over n+1},&
 \sigma_S^2&={n\over(n+1)^2(n+2)},\\
 Y&={S-\mu_S\over\sigma_S},&
 Z&=\sqrt{{n+2\over n}}\,SW.
\end{aligned}                                         \tag{3.1}
\]

Then \((Y,Z)\) is isotropic.  Direct beta-moment calculation gives

\[
\begin{aligned}
 A_n&=\mathbb EY^3
 =-{2(n-1)\sqrt{n+2}\over\sqrt n(n+3)},\\
 \beta_n&=\mathbb E[YZ_iZ_j]\delta_{ij}
 ={2\sqrt{n+2}\over\sqrt n(n+3)},\\
 \gamma_n&=
 \left({n+2\over n}\right)^{3/2}\mathbb ES^3
 ={(n+2)^{3/2}\over\sqrt n(n+3)}.
\end{aligned}                                         \tag{3.2}
\]

More explicitly,

\[
 \mathbb E[YZ_iZ_j]=\beta_n\delta_{ij},
\qquad
 \mathbb E[Z_iZ_jZ_k]
 =\gamma_n\mathbb E[W_iW_jW_k],                       \tag{3.3}
\]

and all moments with two \(Y\)'s and one \(Z\) vanish.

For \(u=(a,v)\in\mathbb R\oplus\mathbb R^{n-1}\), the directional
third-moment matrix is therefore

\[
 M_u=
 \begin{pmatrix}
  aA_n&\beta_n v^T\\
  \beta_n v&a\beta_nI+\gamma_nM_v^W
 \end{pmatrix}.                                       \tag{3.4}
\]

The coefficients satisfy, for \(n\ge2\),

\[
 |A_n|\le2,\qquad
 \beta_n\sqrt{n-1}\le1,\qquad
 \gamma_n\le {8\over5\sqrt2}<1.14.                    \tag{3.5}
\]

It follows from the triangle inequality in Hilbert--Schmidt space that

\[
 \sup_{|u|=1}\|M_u\|_{\rm HS}
 \le4+{8\over5\sqrt2}
 \sup_{|v|=1}\|M_v^W\|_{\rm HS}.                      \tag{3.6}
\]

In the axial direction the sharper exact formula is

\[
 \|M_{(1,0)}\|_{\rm HS}^2
 =A_n^2+(n-1)\beta_n^2\le5.                            \tag{3.7}
\]

Thus the cone direction and the common transverse dilation create only a
constant tensor.  Any growth in a single pyramid must already be present
in the base.  Formula (3.4) also explains the regular simplex
cancellation: its large vertex skew is one axial eigenvalue, while the
\(n-1\) tangential eigenvalues are of order \(1/n\).

## 4. Half-ball: an asymmetric curved stress test

Let \(Y\) be uniform on

\[
 H_n=\{y\in\mathbb R^n:|y|\le1,\ y_1\ge0\}.
\]

Put

\[
 m_n=\mathbb EY_1
 ={\Gamma((n+2)/2)\over
   \sqrt\pi\,\Gamma((n+3)/2)},                         \tag{4.1}
\]

\[
 \sigma_1^2={1\over n+2}-m_n^2,\qquad
 \sigma_\perp^2={1\over n+2},
\]

and whiten by

\[
 X_1={Y_1-m_n\over\sigma_1},\qquad
 X_i={Y_i\over\sigma_\perp}\quad(i\ge2).              \tag{4.2}
\]

Rotational invariance in the last \(n-1\) coordinates implies that the
only nonzero third moments are

\[
 \alpha_n=\mathbb EX_1^3,\qquad
 \beta_n=\mathbb E[X_1X_i^2]\quad(i\ge2).             \tag{4.3}
\]

The exact half-ball moments

\[
 \mathbb EY_1^2={1\over n+2},\qquad
 \mathbb EY_1^3={2m_n\over n+3}
\]

give

\[
\begin{aligned}
 \alpha_n
 &={2m_n/(n+3)-3m_n/(n+2)+2m_n^3\over\sigma_1^3},\\
 \beta_n
 &=-{m_n\over\sigma_1(n+3)}.                          \tag{4.4}
\end{aligned}
\]

For \(u=(u_1,u_\perp)\),

\[
 M_u=
 \begin{pmatrix}
 \alpha_nu_1&\beta_nu_\perp^T\\
 \beta_nu_\perp&\beta_nu_1I_{n-1}
 \end{pmatrix},                                       \tag{4.5}
\]

and hence

\[
 \|M_u\|_{\rm HS}^2
 =\alpha_n^2u_1^2
  +2\beta_n^2|u_\perp|^2
  +(n-1)\beta_n^2u_1^2.                               \tag{4.6}
\]

As \(n\to\infty\),

\[
 \sqrt n\,m_n\to\sqrt{2/\pi},\qquad
 \sqrt n\,\sigma_1\to\sqrt{1-2/\pi},
\]

so

\[
\begin{aligned}
 \alpha_n&\longrightarrow
 {\sqrt{2/\pi}\,(4/\pi-1)\over(1-2/\pi)^{3/2}},\\
 \beta_n&=
 -{\sqrt{2/\pi}\over\sqrt{1-2/\pi}}\,
 {1\over n}+O(n^{-2}).                                \tag{4.7}
\end{aligned}
\]

Equations (4.6)--(4.7) give

\[
 \sup_{|u|=1}\|M_u\|_{\rm HS}=O(1).                   \tag{4.8}
\]

Thus a sharp one-sided boundary and a nonzero axial skew do not produce
many constant-size tangential third moments.

## 5. Smooth full-support versions

Let \(X\) be the isotropic whitened wedge vector from Section 2, let
\(G\sim N(0,I_n)\) be independent, and put

\[
 X_\varepsilon={X+\varepsilon G\over\sqrt{1+\varepsilon^2}}.
\]

This law is isotropic, log-concave, smooth, and has full support.  All
mixed third moments involving \(G\) vanish because \(X\) and \(G\) are
centered and the Gaussian third moments vanish.  Therefore

\[
 M_u(X_\varepsilon)
 ={1\over(1+\varepsilon^2)^{3/2}}M_u(X),              \tag{5.1}
\]

and (2.14) remains valid.

Poincare constants are subadditive under convolution and scale
quadratically.  Hence

\[
 C_P(X_\varepsilon)
 \le {C_P(X)+\varepsilon^2\over1+\varepsilon^2}
 \le192.                                               \tag{5.2}
\]

The argument (2.19)--(2.21) gives, for every smooth
\(g\perp\mathrm {Aff}\) under the convolved law,

\[
 |\mathbb E\nabla g|
 \le192\|D^2g\|_2.                                    \tag{5.3}
\]

Thus neither polyhedral corners nor compact support are responsible for
the positive wedge result.

## 6. What remains open

The tested families instantiate three different compensation mechanisms.

* Product measures place their third moments on coordinate-diagonal tensor
  entries, giving \(\|M_u\|_{\rm HS}=O(1)\).
* Pyramids concentrate their axial variable at scale \(1/n\), making the
  common tangential eigenvalue \(O(1/n)\).
* Paired multiscale wedges keep slice volume nearly balanced to first
  order, but the second-order log-volume curvature is
  \(2\sum_i a_i^2\); axial whitening then divides every slope by
  \((\sum_i a_i^2)^{1/2}\).

The regular simplex, product exponential, half-ball, arbitrary paired
wedge, and a pyramid over any base all fail to realize the abstract
harmonic spectrum allowed by (1.9).  No explicit isotropic log-concave law
with

\[
 \sup_u\|M_u\|_{\rm HS}\longrightarrow\infty
\]

is obtained.

A universal proof of (0.2) would settle the quadratic restriction of
(HMG) and improve the generic thin-shell tensor bound.  It must use more
than rank-by-rank projection variance.  Even such a proof would not
automatically establish (MG) for nonquadratic functions; the boundary
functional (1.12) can couple to higher-order modes.

The wedge theorem supplies one genuinely dependent class where this last
problem is completely closed.  Beyond product-like Lipschitz
parametrizations, the active obstruction is a dimension-free second-order
trace estimate for the score or boundary functional after isotropic
normalization.  Proving that estimate for every log-concave law remains
conjecture-strength; it is not inserted here.
