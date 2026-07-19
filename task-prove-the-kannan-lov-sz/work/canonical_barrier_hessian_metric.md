# Canonical barriers and Hessian metrics: exact identities, stress tests, and the transfer obstruction

## 0. Verdict

There are two different objects that are often both called the
Cheng--Yau or canonical barrier, and they must not be conflated.

1. On every proper open convex domain \(\Omega\subset\mathbb R^d\), the
   direct canonical potential is the unique convex solution
   \[
      \det D^2F=e^{2F},\qquad F(x)\longrightarrow+\infty
      \quad(x\to\partial\Omega).
      \tag{0.1}
   \]
2. For a bounded convex body \(K\subset\mathbb R^n\), one may instead
   apply (0.1) to the homogenization cone
   \[
      \mathcal C_K=\{(tx,t):t>0,\ x\in\operatorname{int}K\}
      \subset\mathbb R^{n+1}
      \tag{0.2}
   \]
   and restrict its canonical potential to the slice \(t=1\).  This is
   the cone-induced canonical barrier of \(K\).  It is generally not the
   direct solution of (0.1) on \(K\).  Already on an interval the direct
   potential is a logarithm of a secant, whereas the cone-induced
   barrier is a logarithmic facet barrier.

Both constructions give natural Hessian metrics.  The main conclusions
of this audit are as follows.

* The canonical potential is strongly self-concordant.  This is a
  dimension-free Frobenius estimate on the *relative variation* of its
  Hessian.  It gives no absolute comparison between the Hessian metric
  and the Euclidean covariance metric.

* The Monge--Ampere equation does not give positive Bakry--Emery
  curvature for the measure one first encounters.  In fact, for
  \(g=D^2F\), the metric-measure space
  \((\Omega,g,e^{2F}dx)\) has non-positive Bakry--Emery tensor.  The
  strong self-concordance estimate is exactly the tensor inequality
  behind this non-positivity.

* For uniform measure on a bounded domain, the Riemannian volume is
  \(d\operatorname{vol}_g=e^Fdx\), so uniform Lebesgue measure is
  \(e^{-F}d\operatorname{vol}_g\).  Its Bakry--Emery tensor has an
  additional drift term of no fixed sign.  Even on an interval the
  curvature-to-metric ratio tends to zero at the boundary.  Therefore
  the Bakry--Emery/Lichnerowicz criterion supplies no positive uniform
  spectral gap.

* Most decisively, **strong self-concordance plus isotropy is
  insufficient even for a dimension-free weighted Poincare inequality**.
  On the isotropic cube there is a smooth barrier whose Hessian is
  strongly self-concordant, whose log determinant is convex, and whose
  standard barrier parameter is at most \(2n-1\), but whose weighted
  Poincare constant is at least \(\pi^2n/6\).  Thus even the optimal
  order \(O(n)\) barrier parameter does not repair the implication.

* The analogous barrier-weighted residual mean-gradient estimate fails
  on the same example.  There are affine-orthogonal boundary-layer
  functions \(r_n\) for which
  \[
     |\mathbb E\nabla r_n|\ge c,
     \qquad
     \|r_n\|_2+
     \left(\mathbb E\operatorname{Tr}
       [D^2r_n\,H^{-1}D^2r_n]\right)^{1/2}
     \le Cn^{-1/4}.
     \tag{0.3}
  \]

* Exact canonical normalization avoids the artificial anisotropic
  rescaling in the preceding counterexample.  A dimension-free
  Poincare inequality in the *direct canonical metric* remains a
  mathematically meaningful possible intermediate theorem.  It does
  hold in the elementary product models below.  However, no local
  curvature proof gives it, and even if it were proved it would not by
  itself imply Euclidean KLS.  The required metric-to-Euclidean energy
  comparison is false on the isotropic simplex, on a skew cone, and for
  product exponential measure.

* The entropic barrier has the missing probabilistic meaning: its third
  derivative is the standardized third-moment tensor of an exponential
  tilt.  Dimension-free strong self-concordance of the entropic barrier
  is precisely the directional third-moment estimate appearing in the
  residual mean-gradient program.  The canonical barrier has strong
  self-concordance unconditionally, but its third derivative is not a
  moment tensor.  Comparing the two third derivatives with a universal
  residual bound is already the missing third-moment theorem, not a
  free consequence of canonicality.

Consequently, the canonical-barrier route is not closed, but its precise
missing mechanism is now isolated: it needs a new global transfer from
canonical intrinsic energy to Euclidean low modes (or from the canonical
cubic tensor to the entropic cubic tensor).  Pointwise metric comparison,
Bakry--Emery curvature, and strong self-concordance alone cannot provide
that transfer.

## 1. Definitions and exact invariances

### 1.1 Proper domains and the direct canonical potential

An open convex domain \(\Omega\subset\mathbb R^d\) is **proper** if it
contains no complete affine line.  The Cheng--Yau existence and
uniqueness theorem, in its real Monge--Ampere formulation, gives a unique
smooth strictly convex function
\(F_\Omega:\Omega\to\mathbb R\) satisfying (0.1).  The boundary condition
is understood as divergence to \(+\infty\) along every sequence tending
to the Euclidean boundary.  This formulation applies without a smoothness
assumption on \(\partial\Omega\); all differential statements are in the
interior.

Write
\[
   H=D^2F,\qquad g=H_{ij}\,dx_i\,dx_j.
   \tag{1.1}
\]
If \(A\in GL_d\), \(b\in\mathbb R^d\), and
\(\widetilde\Omega=A\Omega+b\), uniqueness gives
\[
 F_{\widetilde\Omega}(Ax+b)
   =F_\Omega(x)-\log|\det A|.
 \tag{1.2}
\]
Indeed the Hessian determinant on the left changes by
\(|\det A|^{-2}\).  Thus the Hessian metric is affinely covariant.

If \(\Omega=\mathcal C\) is a proper cone of dimension \(d\), uniqueness
also gives logarithmic homogeneity:
\[
 F_{\mathcal C}(tx)=F_{\mathcal C}(x)-d\log t,
 \qquad t>0.
 \tag{1.3}
\]
Differentiating (1.3) yields the useful identities
\[
 H(x)x=-\nabla F(x),
 \qquad
 \nabla F(x)^TH(x)^{-1}\nabla F(x)=d.
 \tag{1.4}
\]

For a bounded full-dimensional convex body \(K\subset\mathbb R^n\),
the cone (0.2) is proper.  Define the cone-induced slice barrier
\[
   \phi_K(x)=F_{\mathcal C_K}(x,1).
   \tag{1.5}
\]
It is affine-covariant, convex, smooth in \(\operatorname{int}K\), and
diverges at \(\partial K\).  It need not solve
\(\det D^2\phi=e^{2\phi}\) in dimension \(n\).

For a measure supported on a lower-dimensional affine subspace, all of
these definitions are made after identifying its affine hull isometrically
with \(\mathbb R^k\).  This does not create an additional barrier issue;
the covariance and Euclidean gradients are then taken in the affine hull.

### 1.2 Self-concordance and strong self-concordance

Let \(\phi\in C^3(\operatorname{int}K)\) be convex with
\(H_\phi=D^2\phi\succ0\).  It is self-concordant in the standard
normalization if
\[
  |D^3\phi(x)[h,h,h]|
  \le 2\big(h^TH_\phi(x)h\big)^{3/2}
  \tag{1.6}
\]
for all \(x,h\).  It is a \(\nu\)-self-concordant barrier if, in addition,
\[
  \nabla\phi(x)^TH_\phi(x)^{-1}\nabla\phi(x)\le\nu
  \tag{1.7}
\]
and it has the barrier boundary behavior.

Following Laddha--Lee--Vempala, a positive-definite matrix field \(H\)
is **strongly self-concordant** if
\[
 \left\|H(x)^{-1/2}DH(x)[h]H(x)^{-1/2}\right\|_F
 \le 2\sqrt{h^TH(x)h}
 \tag{1.8}
\]
for all \(x,h\).  This is stronger than (1.6), because it controls all
relative eigenvalue velocities in Frobenius norm rather than only one
cubic contraction.

Hildebrand proved that the canonical barrier on a regular convex cone is
a self-concordant barrier with parameter no larger than the dimension.
His canonical-barrier construction, together with the curvature
calculation below, gives (1.8).  Laddha--Lee--Vempala explicitly record
both convexity of \(\log\det H=2F\) and strong self-concordance for this
barrier.  Klartag--Kolesnikov's hyperbolic Monge--Ampere curvature theorem
also gives (1.8) for the direct solution (0.1) on a proper convex domain.

The exact primary references used here are:

* R. Hildebrand, *Canonical Barriers on Convex Cones*, Mathematics of
  Operations Research 39 (2014), 841--850: equation (0.1), affine
  invariance, logarithmic homogeneity, and barrier parameter at most the
  cone dimension.
* A. Laddha, Y. T. Lee, S. Vempala, *Strong Self-Concordance and
  Sampling*, STOC 2020, arXiv:1911.05656: Definition 3 is (1.8), and the
  canonical barrier is listed as strongly self-concordant with convex log
  determinant.
* B. Klartag, A. Kolesnikov, *Remarks on Curvature in the Transportation
  Metric*, arXiv:1604.04165: for every proper open convex domain and the
  boundary-blowup solution of the hyperbolic equation, the relevant
  Bakry--Emery tensor is non-positive.
* D. J. F. Fox, *A Schwarz Lemma for Kahler Affine Metrics and the
  Canonical Potential of a Proper Convex Cone*, Annali di Matematica Pura
  ed Applicata 194 (2015), 1--42: canonical-potential and
  smooth-metric-measure formulations.
* A. Kolesnikov, E. Milman, *Riemannian Metrics on Convex Sets with
  Applications to Poincare and Log-Sobolev Inequalities*, Calculus of
  Variations and Partial Differential Equations 55 (2016), Article 77:
  their generalized Brascamp--Lieb and Bakry--Emery implications require
  the stated geometric-convexity hypotheses and an appropriate positive
  generalized-Ricci lower bound.  The paper explicitly warns that these
  geometric hypotheses are not automatic for a general Hessian metric.

These results concern the intrinsic barrier metric.  None asserts a
dimension-free comparison with covariance.

## 2. Exact Hessian-metric calculus

### 2.1 Riemannian volume and the uniform weighted Laplacian

For the direct potential (0.1),
\[
  d\operatorname{vol}_g
  =\sqrt{\det H}\,dx=e^Fdx.
  \tag{2.1}
\]
If \(K\) is bounded and \(\mu_K=dx/|K|\), then
\[
  d\mu_K=|K|^{-1}e^{-F}\,d\operatorname{vol}_g.
  \tag{2.2}
\]
The intrinsic Dirichlet form is therefore
\[
 \mathcal E_H(f)
 =\int_K\langle H^{-1}\nabla f,\nabla f\rangle\,d\mu_K.
 \tag{2.3}
\]
Its formal generator on \(C_c^\infty(K)\) is
\[
 L_Hf=\partial_i(H^{ij}\partial_jf).
 \tag{2.4}
\]

The cofactor matrix of a Hessian is divergence-free (the Piola identity):
\[
  \partial_i\big((\det H)H^{ij}\big)=0.
  \tag{2.5}
\]
Using \(\partial_i\log\det H=2F_i\) gives
\[
  \boxed{\quad \partial_iH^{ij}=-2H^{ij}F_i.\quad}
  \tag{2.6}
\]
Consequently,
\[
 L_Hf=H^{ij}f_{ij}-2H^{ij}F_if_j.
 \tag{2.7}
\]
For compactly supported smooth \(f\), (2.6) also gives the canonical
integration-by-parts identity
\[
 \int H^{-1}\nabla f\,dx
 =2\int fH^{-1}\nabla F\,dx.
 \tag{2.8}
\]
This is a Stein identity for the vector field
\(2H^{-1}\nabla F\), not for the Euclidean coordinate \(x\).  On a cone,
(1.4) makes this vector field equal to \(-2x\), but Lebesgue measure on
the whole cone is infinite; taking a bounded slice introduces additional
Schur-complement and moving-slice terms.

### 2.2 Curvature: the sign is not the desired one

Let
\[
  T_{ijk}=F_{ijk},
  \qquad
  \mathcal C_{ij}
   =T_{iab}H^{ac}H^{bd}T_{jcd}.
  \tag{2.9}
\]
Then
\[
 h^T\mathcal C h
 =\left\|H^{-1/2}DH[h]H^{-1/2}\right\|_F^2.
 \tag{2.10}
\]
Klartag--Kolesnikov's formula, specialized to
\(\det D^2F=e^{2F}\), is
\[
 \operatorname{Ric}_{g,e^{2F}dx}
 =\frac14\mathcal C-H.
 \tag{2.11}
\]
Their tensor maximum-principle theorem says that (2.11) is non-positive.
Equivalently,
\[
  \mathcal C\preceq4H,
  \tag{2.12}
\]
which is exactly (1.8).  Thus strong self-concordance here encodes a
non-positive Bakry--Emery statement; it must not be cited as a positive
curvature bound.

Relative to Riemannian volume, \(e^{2F}dx=e^F d\operatorname{vol}_g\)
has potential \(-F\), whereas uniform measure has potential \(+F\).
Since
\[
 (\operatorname{Hess}_gF)_{ij}
 =H_{ij}-\frac12T_{ij\ell}H^{\ell k}F_k,
 \tag{2.13}
\]
the uniform-measure tensor is
\[
 \boxed{\quad
 \operatorname{Ric}_{g,dx}
 =\frac14\mathcal C+H
   -T_{ij\ell}H^{\ell k}F_k.
 \quad}
 \tag{2.14}
\]
The last term has no sign.  Strong self-concordance and the usual barrier
gradient estimate give only
\[
 \left|T[h,h,H^{-1}\nabla F]\right|
 \le2\|\nabla F\|_{H^{-1}}\,\|h\|_H^2
 \le2\sqrt\nu\,\|h\|_H^2,
 \tag{2.15}
\]
which loses \(\sqrt\nu\), hence generally \(\sqrt n\).

There is a sharper obstruction: even when (2.14) is positive, it need not
be bounded below by a positive multiple of \(g\).  On the interval
\((-a,a)\), put
\[
 k=\frac\pi{2a},
 \qquad
 F(x)=\log k-\log\cos(kx).
 \tag{2.16}
\]
Then \(F''=e^{2F}=k^2\sec^2(kx)\).  A direct one-dimensional calculation
gives
\[
 \operatorname{Ric}_{g,dx}=k^2\,dx^2,
 \qquad
 \frac{\operatorname{Ric}_{g,dx}}g=\cos^2(kx)\longrightarrow0
 \quad(x\to\partial K).
 \tag{2.17}
\]
Hence the Bakry--Emery criterion has lower-bound constant zero even in
dimension one.

Accordingly, the generalized Hessian-metric Brascamp--Lieb machinery of
Kolesnikov--Milman cannot be invoked here merely from the fact that
\(g\) is Hessian: its positive generalized-Ricci hypothesis fails in the
uniform form needed for Lichnerowicz, already in (2.17), and its separate
geometric-convexity assumptions would also have to be checked.

In intrinsic coordinate
\[
 r=\operatorname{arsinh}(\tan(kx)),
 \tag{2.18}
\]
the metric is \(dr^2\) and uniform measure has density proportional to
\(\operatorname{sech}r\).  This measure does have a positive spectral
gap, but it comes from global exponential confinement, not from a
positive local curvature lower bound.

## 3. What weighted Poincare would say

For a probability measure \(\mu\) on \(K\) and a positive-definite field
\(H\), define
\[
 C_P^H(\mu)
 =\sup_{f\not\equiv\mathrm{const}}
 \frac{\operatorname{Var}_\mu f}
 {\int\langle H^{-1}\nabla f,\nabla f\rangle\,d\mu}.
 \tag{3.1}
\]
A dimension-free bound for the exactly normalized direct canonical metric
would be an interesting intrinsic theorem.  Strong self-concordance by
itself does not imply it.

### 3.1 A sharp counterexample to the implication from barrier axioms

Let
\[
 K_n=(-\sqrt3,\sqrt3)^n,
 \qquad
 \mu_n=\operatorname{Unif}(K_n).
 \tag{3.2}
\]
This measure is isotropic.  Put
\[
 a=\sqrt3,
 \quad k=\frac\pi{2a},
 \quad F_0(s)=\log k-\log\cos(ks),
 \tag{3.3}
\]
and define the anisotropically weighted barrier
\[
 \phi_n(x)=nF_0(x_1)+\sum_{i=2}^nF_0(x_i).
 \tag{3.4}
\]
Writing \(\alpha_1=n\) and \(\alpha_i=1\) for \(i\ge2\), its Hessian is
\[
 H_{ii}(x)=\alpha_i k^2\sec^2(kx_i),
 \qquad H_{ij}=0\quad(i\ne j).
 \tag{3.5}
\]

This example satisfies all of the following.

1. It is a smooth barrier on an isotropic body.
2. It is strongly self-concordant.  Indeed
   \[
    \left\|H^{-1/2}DH[h]H^{-1/2}\right\|_F^2
    =4k^2\sum_i h_i^2\tan^2(kx_i)
    \le4h^THh.
    \tag{3.6}
   \]
3. Its log determinant is convex:
   \[
    \log\det H=\mathrm{const}+2\sum_i\log\sec(kx_i).
    \tag{3.7}
   \]
4. It has optimal-order barrier parameter:
   \[
    \nabla\phi_n^TH^{-1}\nabla\phi_n
    =\sum_i\alpha_i\sin^2(kx_i)\le2n-1.
    \tag{3.8}
   \]

Nevertheless, for \(f(x)=x_1\),
\[
 \operatorname{Var}_{\mu_n}f=1,
 \qquad
 \int\langle H^{-1}\nabla f,\nabla f\rangle\,d\mu_n
 =\frac1{nk^2}\mathbb E\cos^2(kX_1)
 =\frac6{\pi^2n}.
 \tag{3.9}
\]
Therefore
\[
 \boxed{\quad C_P^H(\mu_n)\ge\frac{\pi^2}{6}n.\quad}
 \tag{3.10}
\]

This proves that isotropy, convex \(\log\det H\), strong
self-concordance with the universal constant in (1.8), and even
\(\nu=O(n)\) do not imply a dimension-free intrinsic spectral gap.
Exact canonical PDE normalization is genuinely extra information.

### 3.2 The unweighted-to-weighted energy comparison is pointwise

For a continuous positive-definite field on an open set of positive
Lebesgue density,
\[
 \sup_{0\ne f\in C_c^\infty(K)}
 \frac{\int\langle H^{-1}\nabla f,\nabla f\rangle\,d\mu}
 {\int|\nabla f|^2\,d\mu}
 =\operatorname*{ess\,sup}_{x\in K}
   \lambda_{\max}(H(x)^{-1}).
 \tag{3.11}
\]
The upper bound is immediate.  For the lower bound, localize near a
Lebesgue point of \(H\) and use a high-frequency cutoff plane wave in a
maximizing eigendirection; cutoff gradients are lower order in frequency.

Thus the black-box transfer
\[
 \mathcal E_H(f)\le C\int|\nabla f|^2d\mu
 \quad\hbox{for every }f
 \tag{3.12}
\]
is exactly the pointwise comparison \(H\succeq C^{-1}I\).  Section 5
shows that (3.12) is false for the cone-induced canonical barrier of an
isotropic simplex and for the natural canonical metric of product
exponential measure.  Any successful transfer must therefore be
low-mode- or capacity-sensitive; it cannot be an all-functions energy
domination.

## 4. The residual mean-gradient gate

### 4.1 The desired Euclidean statement

For isotropic \(\mu\), write
\[
 g\perp\mathrm{Aff}
 \quad\Longleftrightarrow\quad
 \mathbb Eg=0,
 \qquad
 \mathbb E[Xg(X)]=0.
 \tag{4.1}
\]
The surviving mean-gradient target in the main proof program is
\[
 |\mathbb E\nabla g|
 \le C\big(\|g\|_2+\|D^2g\|_2\big),
 \qquad g\perp\mathrm{Aff}.
 \tag{4.2}
\]
Its Hessian-only strengthening is
\[
 |\mathbb E\nabla g|\le C\|D^2g\|_2.
 \tag{4.3}
\]
The eigenfunction-residual calculation in the main registry shows that
(4.2) with a universal constant forces a universal spectral gap.  It is
therefore a conjecture-strength gate, not a general consequence one may
assume about a canonical metric.

Already on quadratic functions, (4.3) is the directional third-moment
problem.  For \(B\in\operatorname{Sym}_n\), set
\[
 q_B(x)=x^TBx-\operatorname{Tr}B,
 \qquad
 d_B=\mathbb E[Xq_B(X)],
 \qquad
 g_B=q_B-d_B\cdot x.
 \tag{4.4}
\]
Then \(g_B\perp\mathrm{Aff}\),
\[
 \mathbb E\nabla g_B=-d_B,
 \qquad D^2g_B=2B,
 \tag{4.5}
\]
and
\[
 u\cdot d_B
 =\left\langle
   \mathbb E[(u\cdot X)(XX^T-I)],B
  \right\rangle_{HS}.
 \tag{4.6}
\]
Hence the quadratic restriction of (4.3) is exactly
\[
 \sup_{|u|=1}
 \left\|\mathbb E[(u\cdot X)(XX^T-I)]\right\|_{HS}\le C.
 \tag{4.7}
\]

### 4.2 The natural barrier-weighted replacement is false

A canonical-metric attempt would naturally replace the Hessian norm in
(4.2) by
\[
 \|D^2g\|_{H^{-1},2}^2
 :=\int\operatorname{Tr}
      (D^2g\,H^{-1}D^2g)\,d\mu
 =\sum_i\int
   \langle H^{-1}\nabla\partial_i g,
           \nabla\partial_i g\rangle\,d\mu.
 \tag{4.8}
\]
Even the version with the \(L^2\) residual term fails under the barrier
hypotheses of Section 3.

Use the cube and barrier (3.2)--(3.5).  Fix a smooth
\(\eta:[0,\infty)\to[0,1]\) with
\[
 \eta(0)=1,
 \qquad \eta(s)=0\quad(s\ge1),
 \qquad \|\eta''\|_\infty<\infty.
 \tag{4.9}
\]
For \(0<\varepsilon<1\), put
\[
 b_\varepsilon(x)=
 \eta\left(\frac{a-x_1}{\varepsilon}\right),
 \quad
 c_\varepsilon=\mathbb Eb_\varepsilon,
 \quad
 d_\varepsilon=\mathbb E[X_1b_\varepsilon],
 \tag{4.10}
\]
and
\[
 r_\varepsilon(x)=b_\varepsilon(x)-c_\varepsilon-d_\varepsilon x_1.
 \tag{4.11}
\]
Since \(\mathbb EX_1=0\), \(\mathbb EX_1^2=1\), and the other
coordinates are independent and centered,
\[
 r_\varepsilon\perp\mathrm{Aff}.
 \tag{4.12}
\]
Moreover,
\[
 \mathbb E\partial_1b_\varepsilon
 =\frac{b_\varepsilon(a)-b_\varepsilon(-a)}{2a}
 =\frac1{2a},
 \qquad
 |d_\varepsilon|\le\frac\varepsilon2.
 \tag{4.13}
\]
Thus, for all sufficiently small \(\varepsilon\),
\[
 |\mathbb E\nabla r_\varepsilon|\ge\frac1{4a}.
 \tag{4.14}
\]
On the other hand,
\[
 \|r_\varepsilon\|_2\le C\sqrt\varepsilon.
 \tag{4.15}
\]

On the support of \(b_\varepsilon''\), write
\(\delta=a-x_1\in[0,\varepsilon]\).  Since
\[
 (H^{-1})_{11}
 =\frac{\cos^2(kx_1)}{nk^2}
 =\frac{\sin^2(k\delta)}{nk^2}
 \le\frac{\delta^2}{n},
 \tag{4.16}
\]
and \(|b_\varepsilon''|\le C\varepsilon^{-2}\),
\[
 \|D^2r_\varepsilon\|_{H^{-1},2}^2
 \le\frac{C}{n\varepsilon^4}
       \int_0^\varepsilon\delta^2d\delta
 \le\frac{C}{n\varepsilon}.
 \tag{4.17}
\]
Taking \(\varepsilon=n^{-1/2}\) gives
\[
 \|r_\varepsilon\|_2+
 \|D^2r_\varepsilon\|_{H^{-1},2}
 \le Cn^{-1/4},
 \tag{4.18}
\]
while (4.14) stays bounded below.  This proves (0.3).

The exactly normalized direct canonical cube does not have the extra
factor \(n\) in its first coordinate, so this counterexample does not
disprove a special theorem for the exact canonical solution.  It does
prove that strong self-concordance, convex log determinant, isotropy, and
\(\nu=O(n)\) cannot be the proof of such a theorem.

## 5. Required model tests

### 5.1 Isotropic cube: the exact direct metric is benign

For
\(K=\prod_{i=1}^n(-a_i,a_i)\), the direct solution of (0.1) is
\[
 F(x)=\sum_{i=1}^n
 \left(\log k_i-\log\cos(k_ix_i)\right),
 \qquad k_i=\frac\pi{2a_i}.
 \tag{5.1}
\]
Indeed its Hessian is diagonal and
\(\det D^2F=e^{2F}\).  On the isotropic cube \(a_i=\sqrt3\),
\[
 H^{-1}\preceq \frac{12}{\pi^2}I.
 \tag{5.2}
\]

The intrinsic weighted gap is also dimension-free.  In one coordinate,
with \(u=kx\), the Dirichlet form is
\[
 \frac1\pi\int_{-\pi/2}^{\pi/2}
       \cos^2u\,|f'(u)|^2du.
 \tag{5.3}
\]
The one-dimensional weighted Hardy--Muckenhoupt criterion gives
\[
 C_P\le4\max(B_-,B_+),
 \tag{5.4}
\]
where, for example,
\[
 B_+=\sup_{0<u<\pi/2}
 \frac{\pi/2-u}{\pi}
 \int_0^u\frac{\pi\,dt}{\cos^2t}
 =\sup_{0<u<\pi/2}(\pi/2-u)\tan u\le1.
 \tag{5.5}
\]
Tensorization of this product Dirichlet form gives
\[
 C_P^H(\operatorname{Unif}K)\le4
 \tag{5.6}
\]
for every dimension.  Thus the exact direct metric passes the cube test;
the anisotropic rescaling (3.4), not the cube itself, exposes the logical
insufficiency of the barrier axioms.

### 5.2 Isotropic Euclidean ball: cone-induced metric dominates Euclidean

Let \(K=B_2^n(R)\) with \(R=\sqrt{n+2}\), so uniform measure is
isotropic.  Its homogenization is a Lorentz cone of dimension
\(m=n+1\).  Up to the additive constant fixed by (0.1), its canonical
potential is
\[
 F_{\mathcal C}(x,t)
 =-\frac m2\log\left(t^2-\frac{|x|^2}{R^2}\right)+c_m.
 \tag{5.7}
\]
The slice barrier is
\[
 \phi_K(x)=-\frac m2\log\left(1-\frac{|x|^2}{R^2}\right)+c_m,
 \tag{5.8}
\]
and
\[
 D^2\phi_K(x)
 =\frac m{R^2-|x|^2}I
  +\frac{2m}{(R^2-|x|^2)^2}xx^T
 \succeq\frac{n+1}{n+2}I.
 \tag{5.9}
\]
Hence
\[
 (D^2\phi_K)^{-1}\preceq\frac{n+2}{n+1}I\preceq2I.
 \tag{5.10}
\]
There is no metric-transfer loss in this model.  Notice that (5.8) is
the cone-induced barrier; the direct solution of (0.1) on the ball is a
different radial function.

### 5.3 Isotropic simplex: pointwise transfer loses \(n\)

Let
\[
 \Delta_n=\{\lambda_i>0:\ \lambda_1+\cdots+\lambda_{n+1}=1\}
 \tag{5.11}
\]
in its \(n\)-dimensional affine hull.  Its homogenization cone is linearly
isomorphic to \(\mathbb R_+^{n+1}\), whose canonical potential is
\(-\sum_i\log z_i\).  Thus the cone-induced barrier is exactly
\[
 \phi(\lambda)=-\sum_{i=1}^{n+1}\log\lambda_i.
 \tag{5.12}
\]

For uniform measure, \(S=\lambda_1\sim\operatorname{Beta}(1,n)\), so
\[
 \mathbb ES=\frac1{n+1},
 \qquad
 \sigma^2:=\operatorname{Var}S
 =\frac{n}{(n+1)^2(n+2)}.
 \tag{5.13}
\]
Along the permutation-symmetric axial line
\[
 \lambda_1=s,
 \qquad
 \lambda_2=\cdots=\lambda_{n+1}=\frac{1-s}{n},
 \tag{5.14}
\]
put \(z=(s-\mathbb ES)/\sigma\).  Permutation symmetry makes this axis
an eigendirection of the barrier Hessian, and
\[
 H_{zz}(s)
 =\sigma^2\left(\frac1{s^2}+\frac n{(1-s)^2}\right).
 \tag{5.15}
\]
At the barycenter,
\[
 H_{zz}\left(\frac1{n+1}\right)
 =\frac{n+1}{n+2}\asymp1.
 \tag{5.16}
\]
At \(s=1/2\), however,
\[
 H_{zz}(1/2)
 =\frac{4n}{(n+1)(n+2)}\asymp\frac1n,
 \tag{5.17}
\]
so
\[
 \lambda_{\max}(H^{-1})\ge
 \frac{(n+1)(n+2)}{4n}\asymp n.
 \tag{5.18}
\]
This is an interior point and the estimate persists on an open
neighborhood.  By (3.11), an all-functions transfer from canonical energy
to Euclidean energy necessarily loses a factor \(n\).

The bad region is rare:
\[
 \mathbb P\{S\ge1/2\}=2^{-n}.
 \tag{5.19}
\]
Thus (5.18) does not rule out a low-mode or capacity-sensitive transfer;
it rules out only the tempting pointwise or averaged-all-gradients
comparison.

### 5.4 Product exponential measure: intrinsic gap but unbounded transfer

Let \(Y_i\) be independent \(\operatorname{Exp}(1)\) variables and
\(X_i=Y_i-1\).  Then \(X\) is isotropic and log-concave.  On the support
\((-1,\infty)^n\), the positive-orthant canonical barrier is
\[
 \phi(x)=-\sum_{i=1}^n\log(x_i+1),
 \qquad
 H^{-1}(x)=\operatorname{diag}((x_i+1)^2).
 \tag{5.20}
\]
It satisfies the strong self-concordance inequality with equality in the
coordinate directions and has parameter \(n\).  At the mean, \(H=I\),
but \(\|H^{-1}\|_{op}\) is unbounded in the tails.  Therefore (3.12)
fails with every finite constant.

The intrinsic weighted spectral gap nevertheless tensorizes with a
universal constant.  In one dimension the measure is
\(e^{-y}dy\) on \((0,\infty)\) and the energy is
\[
 \int_0^\infty y^2|f'(y)|^2e^{-y}dy.
 \tag{5.21}
\]
The weighted Hardy--Muckenhoupt quantities are finite uniformly:
near zero,
\[
 (1-e^{-x})\int_x^{\log2}\frac{e^t}{t^2}dt=O(1),
 \tag{5.22}
\]
and at infinity,
\[
 e^{-x}\int_{\log2}^x\frac{e^t}{t^2}dt=O(1).
 \tag{5.23}
\]
Equivalently, in intrinsic coordinate \(r=\log y\), the measure has
density proportional to \(e^{r-e^r}\) and the energy is Euclidean.

Thus the weighted geometry sees a perfectly good product spectral gap,
while black-box transfer to Euclidean energy is impossible.  This is the
cleanest unbounded-support version of the obstruction.

### 5.5 A skew cone: canonical scaling still has a rare \(n\)-loss

Let the total dimension be \(n\) and consider the truncated circular cone
\[
 K_n=\{(t,y)\in\mathbb R\times\mathbb R^{n-1}:
       0<t<1,\ |y|<t\}.
 \tag{5.24}
\]
For its uniform law,
\[
 T\sim\operatorname{Beta}(n,1),
 \qquad
 \mathbb ET=\frac n{n+1},
 \qquad
 \sigma_T^2=\frac{n}{(n+1)^2(n+2)}.
 \tag{5.25}
\]
Moreover
\[
 \operatorname{Var}(Y_i)
 =\frac{n}{(n+1)(n+2)}.
 \tag{5.26}
\]
Therefore scaling \(t\) by \(\sigma_T^{-1}\) and every transverse
coordinate by the reciprocal square root of (5.26) puts the body in
isotropic position.

Consider the cone-plus-cap barrier
\[
 \phi_n(t,y)
 =2\left[-\log(1-t)-\frac n2\log(t^2-|y|^2)\right].
 \tag{5.27}
\]
The Lorentz term has the canonical \(n/2\) scaling; the cap is the
one-dimensional logarithmic barrier.  Each summand is strongly
self-concordant, and a factor two makes their sum strongly
self-concordant by the triangle inequality and
\(\sqrt a+\sqrt b\le\sqrt{2(a+b)}\).  Its standard parameter is
\(O(n)\).  This explicit composite is used only as a stress-test metric;
the exact canonical barrier of the homogenized capped cone is not known
in closed form.

On the axis \(y=0\),
\[
 \partial_{tt}\phi_n
 =2\left(\frac1{(1-t)^2}+\frac n{t^2}\right).
 \tag{5.28}
\]
In the standardized axial coordinate
\(z=(t-\mathbb ET)/\sigma_T\), at \(t=1/2\),
\[
 H_{zz}
 =\sigma_T^2\partial_{tt}\phi_n
 =\frac{8n}{(n+1)(n+2)}\asymp\frac1n.
 \tag{5.29}
\]
Hence the inverse metric is of order \(n\) in a deep interior region,
although
\[
 \mathbb P\{T\le1/2\}=2^{-n}.
 \tag{5.30}
\]
At the typical location \(t=1-O(1/n)\), the cap curvature in (5.28)
is of order \(n^2\), and after multiplication by \(\sigma_T^2\) the
standardized Hessian is of order one.  As for the simplex, the metric
automatically normalizes the high-mass region and becomes weak only in
an exponentially rare deep-cone region.

For axial functions the intrinsic weighted Poincare constant is still
universal.  This assertion uses the full inverse Hessian, not merely its
value on the axis.  A block inversion of the Lorentz Hessian gives
\[
 (H^{-1})_{tt}(t,y)
 =\left(
    \frac{2n}{t^2+|y|^2}+\frac2{(1-t)^2}
   \right)^{-1}.
 \tag{5.31}
\]
Since \(|y|<t\), this is between the reciprocal of (5.28) and twice that
reciprocal.  The axial marginal density is
\(\rho_n(t)=nt^{n-1}\).  If \(a_n(t)\) denotes the reciprocal of (5.28),
then its median is \(m_n=2^{-1/n}\), and the one-dimensional criterion
uses
\[
 \frac1{a_n(t)\rho_n(t)}
 =2t^{-n-1}+\frac2n\frac{t^{1-n}}{(1-t)^2}.
 \tag{5.32}
\]
For \(x<m_n\), multiplication by \(\mu(0,x)=x^n\) makes the first term
at most \(2/n\); splitting the second integral at \(1/2\) makes it
universally bounded.  For \(x>m_n\), use
\(1-x^n\le n(1-x)\) and \(t^{1-n}\le m_n^{1-n}<4\) to get a universal
bound.  Thus the rare inverse-metric growth does not create an axial
intrinsic bottleneck, but it again defeats pointwise Euclidean transfer.

## 6. Canonical versus entropic barriers

### 6.1 The entropic cubic tensor is exactly the third moment

For a bounded convex body \(K\), define the log-Laplace transform
\[
 A(\theta)=\log\int_K e^{\langle\theta,x\rangle}dx
 \tag{6.1}
\]
and the entropic barrier \(E=A^*\).  Let \(\mu_\theta\) be the exponential
tilt, with mean \(m_\theta\) and covariance \(\Sigma_\theta\).  Then
\[
 \nabla A(\theta)=m_\theta,
 \qquad
 D^2A(\theta)=\Sigma_\theta,
 \qquad
 D^2E(m_\theta)=\Sigma_\theta^{-1}.
 \tag{6.2}
\]
Let
\[
 Z=\Sigma_\theta^{-1/2}(X-m_\theta),
 \qquad
 a=\Sigma_\theta^{-1/2}h.
 \tag{6.3}
\]
Differentiating the inverse covariance gives the exact identity
\[
 \left\|(D^2E)^{-1/2}
       D(D^2E)[h]
       (D^2E)^{-1/2}\right\|_F
 =\left\|
    \mathbb E_{\mu_\theta}[(a\cdot Z)(ZZ^T-I)]
  \right\|_F.
 \tag{6.4}
\]
Also
\[
 h^TD^2E(m_\theta)h=|a|^2.
 \tag{6.5}
\]
Consequently, dimension-free strong self-concordance of the entropic
barrier is precisely a universal directional third-moment bound for all
isotropic exponential tilts.  Laddha--Lee--Vempala prove their entropic
and universal-barrier strong-self-concordance estimates by invoking the
known KLS-dependent third-moment bound, and state that their argument is,
up to logarithmic factors, equivalent to KLS.  That result cannot be
reversed here as a proof of KLS.

### 6.2 What a canonical-to-entropic comparison would have to prove

Suppose \(K\) is isotropic with barycenter \(b=0\).  At \(b\),
\[
 D^2E(b)=I,
 \qquad
 D^3E(b)[u,\cdot,\cdot]
 =-\mathbb E[(u\cdot X)(XX^T-I)].
 \tag{6.6}
\]
Let \(\phi\) be either canonical barrier.  Its strong
self-concordance controls
\[
 (D^2\phi)^{-1/2}
 D^3\phi[u]
 (D^2\phi)^{-1/2},
 \tag{6.7}
\]
but this tensor has no known moment interpretation.

A proposed comparison needs at least the following two ingredients:
\[
 cI\preceq D^2\phi(b)\preceq CI,
 \tag{6.8}
\]
and, for the residual \(R=E-\phi\),
\[
 \sup_{|u|=1}\|D^3R(b)[u,\cdot,\cdot]\|_{HS}\le C.
 \tag{6.9}
\]
If (6.8) holds, canonical strong self-concordance bounds the canonical
term in Euclidean Hilbert--Schmidt norm.  Equation (6.6) then shows that
(6.9) is, up to that already bounded term, exactly the missing
directional third-moment estimate (4.7).  Thus (6.9) is not a harmless
barrier-comparison lemma; it is the quadratic residual mean-gradient
gate itself.

Even a uniform Hessian comparison at one point would not control the
derivative difference.  Two uniformly comparable positive Hessian
fields can have arbitrarily different derivatives.  Any valid
canonical-to-entropic argument must exploit an additional PDE,
variational, or probabilistic identity controlling the residual cubic
tensor; Loewner comparability alone is insufficient.

## 7. Exact blocked and live statements

### 7.1 Blocked statements

The following statements cannot be used as automatic consequences of
canonical barriers.

1. **Strong self-concordance \(+\) isotropy \(\Rightarrow\) weighted
   Poincare.**  This is disproved by (3.2)--(3.10), even with convex log
   determinant and parameter \(O(n)\).

2. **Bakry--Emery \(\Rightarrow\) canonical weighted gap.**  The
   hyperbolic tensor has the wrong sign, and the uniform tensor has no
   positive lower bound even in dimension one.

3. **Pointwise canonical-to-Euclidean energy domination.**  By (3.11),
   this is equivalent to a uniform lower bound on the barrier Hessian.
   It fails by factors \(n\) on the simplex and skew cone and fails
   without any finite constant for product exponentials.

4. **Canonical strong self-concordance \(\Rightarrow\) entropic strong
   self-concordance.**  The missing derivative residual is the third
   moment tensor itself.

5. **Barrier-weighted residual mean-gradient coercivity.**  The natural
   version (4.8) is disproved by the boundary layer (4.9)--(4.18).

### 7.2 Live statements, with their exact burden

Two genuinely new statements remain logically possible.

**(CB-P)** For the *exact direct* canonical potential on every bounded
convex body,
\[
 \operatorname{Var}_{\mu_K}f
 \le C\int_K\langle(D^2F_K)^{-1}\nabla f,\nabla f\rangle d\mu_K
 \tag{7.1}
\]
with universal \(C\).

This holds in the product interval/cube model and is compatible with the
one-dimensional and cone tests above.  It is not implied by local
curvature or abstract barrier axioms.  Even if proved, (7.1) is only an
intrinsic gap; it does not imply KLS without a new low-mode transfer.

**(CB-LM)** A low-mode/capacity comparison which avoids the false
all-functions estimate (3.12).  One exact sufficient form is: for some
universal \(\alpha\) and \(\beta<C^{-1}\), where \(C\) is the constant in
(7.1),
\[
 \mathcal E_H(f)
 \le\alpha\int|\nabla f|^2d\mu_K
    +\beta\operatorname{Var}_{\mu_K}f
 \quad\text{for all mean-zero }f.
 \tag{7.2}
\]
Combining (7.1) and (7.2) yields
\[
 \operatorname{Var}f
 \le\frac{C\alpha}{1-C\beta}
      \int|\nabla f|^2d\mu_K.
 \tag{7.3}
\]
Thus a universal proof of (7.2) is already KLS-strength.  Unlike (3.12),
it could in principle exploit that the large inverse metric in the
simplex and skew-cone tests lives in exponentially rare regions.  No
such estimate follows from strong self-concordance.

The alternative live burden is the canonical-to-entropic residual
estimate (6.9), which is the third-moment part of the residual
mean-gradient route.  A full functional extension from (6.9) to (4.2)
would still be needed to complete KLS.

## 8. Final registry entry

**Family:** Cheng--Yau/canonical barrier, Hessian metric, and
canonical-to-entropic comparison.

**Exact unconditional assets:** affine covariance; logarithmic
homogeneity on cones; \(\det D^2F=e^{2F}\); Piola identity (2.6);
strong self-concordance (2.12); explicit product, Lorentz, orthant, and
cone models.

**New no-go theorem:** isotropy + strong self-concordance + convex
\(\log\det H\) + \(\nu=O(n)\) does not give a dimension-free weighted
Poincare constant or barrier-weighted residual mean-gradient estimate.
The isotropic cube examples (3.4) and (4.11) give losses \(n\) and
\(n^{1/4}\), respectively.

**Geometric obstruction:** canonical inverse metrics can be of order
\(n\) in exponentially rare interior regions of the isotropic simplex
and skew cone, and unbounded in exponential tails.  Intrinsic spectral
gaps may remain universal, so only a capacity-sensitive transfer could
work.

**Conjecture-strength comparisons:** a universal low-mode transfer such
as (7.2), the Euclidean residual mean-gradient estimate (4.2), or the
canonical-to-entropic cubic residual estimate (6.9).  None may be
assumed.

**Status:** blocked as a direct self-concordance/curvature proof.  Reopen
only with a new global canonical intrinsic isoperimetry theorem together
with a rare-region capacity transfer, or with a new PDE identity that
controls the canonical-to-entropic residual cubic tensor.
