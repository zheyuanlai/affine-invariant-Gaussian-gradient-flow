# Tensor mixed modes and difference-of-convex buffers

## Executive conclusion

Let \(\mu\) be isotropic, let \(f\) be a normalized first eigenfunction,
and write
\[
 \int f\,d\mu=0,\qquad \int f^2\,d\mu=1,\qquad
 \int |\nabla f|^2\,d\mu=\lambda .
\]
On the tensor square put
\[
 h(x,y)=f(x)f(y).
\]
There is a clean negative result for the proposed difference-of-convex
(DC) buffer.

The Bochner estimate does imply
\[
 \|D^2h\|_{L^2(\mu\otimes\mu;\mathrm{HS})}\le 2\lambda.
\]
However, this is the correct scale, not an overestimate.  If
\[
 h-P=k_+-k_-
\]
for a polynomial \(P\) of total degree at most two and convex functions
\(k_\pm\), then, with
\(S=D^2k_++D^2k_-\),
\[
 \boxed{
 \frac1\lambda\|S\|_{L^2(\mu\otimes\mu;\mathrm{HS})}
 \ge
 \sqrt{2\bigl(1-\lambda^2|a|^4\bigr)}
 \ge \sqrt{2(1-\lambda^2)},
 }
 \tag{0.1}
\]
where
\[
 a=\int x f(x)\,d\mu(x),\qquad |a|\le1.
\]
Consequently the best normalized \(L^2\) DC curvature tends to at least
\(\sqrt2\), rather than to zero, if \(\lambda\to0\).  The obstruction is
the mixed Hessian block
\(\nabla f(x)\otimes\nabla f(y)\); a quadratic can remove only its mean.

There is an even sharper obstruction at the level of the actual KKT
action.  Let
\[
 q=h-\Pi_{\mathcal P_2}h
\]
be the part of \(h\) orthogonal to degree-at-most-two polynomials.  If a
normalized PSD normal-cone multiplier \(N\) carries the mixed source
modulo a polynomial, in the sense that for some \(p\in\mathcal P_2\),
\[
 \int W(p-h)\,d(\mu\otimes\mu)=-\int D^2W:N             \tag{0.2}
\]
for all tests \(W\), then every exact DC representation
\(q=k_+-k_-\) satisfies
\[
 \boxed{
 \int (D^2k_++D^2k_-):N
 \ge \|q\|_2^2.
 }
 \tag{0.3}
\]
Thus the buffer action is order one whenever the mixed residual is order
one.  Proving that this action is \(o(1)\) would already prove that
\(q=o(1)\), which is precisely the desired mixed full-matrix conclusion.
It does not follow from Bochner control.

On unbounded support, compact truncation produces honest global DC
decompositions, but their pointwise convexifying constants are not
controlled by the Bochner \(L^2\) estimate, and (0.1) survives the
truncation.  Hence truncation does not create a vanishing buffer.

The conclusion is not that a state-specific KKT argument is impossible.
It is that the proposed generic implication
\[
 \text{Bochner }L^2\text{ control}
 \Longrightarrow
 \text{DC buffer with }o(1)\text{ normalized normal-cone action}
\]
is false.  Any successful version must use additional structure that
controls the actual multiplier on zero-curvature sets, or must obtain
cancellation special to the spectral KKT equation.  Positivity and
complementarity alone do not do this.

## 1. Setting and two elementary eigenfunction identities

Let \(\Omega\subseteq\mathbb R^n\) be a convex support, possibly all of
\(\mathbb R^n\), and let
\[
 d\mu=Z^{-1}e^{-V}\,dx
\]
with \(V\) convex.  On a bounded support use the natural Neumann
realization.  Assume enough regularity for the displayed Hessians; the
Sobolev formulation is given below.

The weak eigenfunction equation is
\[
 \int \nabla f\cdot\nabla\phi\,d\mu
 =\lambda\int f\phi\,d\mu.                               \tag{1.1}
\]
Testing with \(\phi=x_i\), using cutoffs first on unbounded support, gives
\[
 m:=\int\nabla f\,d\mu
 =\lambda a,
 \qquad
 a:=\int xf\,d\mu.                                      \tag{1.2}
\]
Since the coordinate functions form an orthonormal family in
\(L^2(\mu)\), Bessel's inequality gives
\[
 |a|^2=\sum_i\left(\int x_i f\,d\mu\right)^2\le1.        \tag{1.3}
\]
The coordinate Rayleigh quotients also give
\[
 0<\lambda\le1.                                         \tag{1.4}
\]

The weighted Bochner--Reilly identity gives
\[
 B:=\int\|D^2f\|_{\mathrm{HS}}^2\,d\mu\le\lambda^2.    \tag{1.5}
\]
Boundary and curvature terms have the favorable sign.  Only (1.1),
(1.3), and (1.5) are used below.

## 2. Exact Hessian calculation for the tensor mode

Write \(\rho=\mu\otimes\mu\), and let
\(g=\nabla f\).  The block Hessian of
\(h(x,y)=f(x)f(y)\) is
\[
 D^2h(x,y)=
 \begin{pmatrix}
  f(y)D^2f(x) & g(x)g(y)^T\\
  g(y)g(x)^T & f(x)D^2f(y)
 \end{pmatrix}.                                         \tag{2.1}
\]
Independence and \(\|f\|_2=1\) yield the exact identity
\[
 \begin{aligned}
 \int\|D^2h\|_{\mathrm{HS}}^2\,d\rho
 &=2\int\|D^2f\|_{\mathrm{HS}}^2\,d\mu
   +2\left(\int|\nabla f|^2\,d\mu\right)^2\\
 &=2B+2\lambda^2
 \le4\lambda^2.                                        \tag{2.2}
 \end{aligned}
\]
Thus \(D^2h=O(\lambda)\) in \(L^2\).  The mixed block alone already has
exact squared norm \(\lambda^2\), so this scale cannot in general be
improved.

### Quadratic corrections cannot remove the fluctuation

For \(P\in\mathcal P_2(\mathbb R^{2n})\), its \(xy\) Hessian block is a
constant matrix, say \(C\).  Therefore
\[
 D^2_{xy}(h-P)=g(x)g(y)^T-C.                              \tag{2.3}
\]
The best constant approximation in the Hilbert space
\(L^2(\rho;\mathbb R^{n\times n})\) is its mean
\[
 C_*=\mathbb E[g(X)g(Y)^T]=mm^T.
\]
Consequently
\[
 \begin{aligned}
 \inf_C\int\|g(x)g(y)^T-C\|_{\mathrm{HS}}^2\,d\rho
 &=\left(\int|g|^2\,d\mu\right)^2-|m|^4\\
 &=\lambda^2-\lambda^4|a|^4.                             \tag{2.4}
 \end{aligned}
\]
Both off-diagonal Hessian blocks occur in (2.1), so for every
\(P\in\mathcal P_2\),
\[
 \|D^2(h-P)\|_{L^2(\rho;\mathrm{HS})}^2
 \ge2\lambda^2\bigl(1-\lambda^2|a|^4\bigr).             \tag{2.5}
\]
This lower bound is independent of dimension.  In particular,
\[
 \inf_{P\in\mathcal P_2}
 \frac{\|D^2(h-P)\|_2}{\lambda}
 \ge\sqrt{2(1-\lambda^2)}
 \longrightarrow\sqrt2                                  \tag{2.6}
\]
along any hypothetical sequence with \(\lambda\to0\).

The appearance of \(\lambda a\), rather than \(a\), in (2.4) matters.
Even if \(f\) is very close in value norm to a linear function, its mean
gradient is only \(\lambda a\).  At small \(\lambda\), the gradient of
the nonlinear remainder must cancel almost all of the gradient of that
linear part.  The tensor mixed Hessian detects this cancellation.

## 3. A quantitative DC lower bound

For a twice weakly differentiable convex function \(k\), its Hessian is
positive semidefinite almost everywhere.  Suppose
\[
 h-P=k_+-k_-                                              \tag{3.1}
\]
and put
\[
 A_\pm=D^2k_\pm\succeq0,
 \qquad S=A_++A_-.
\]
Then
\[
 D^2(h-P)=A_+-A_-.
\]
Pointwise,
\[
 \begin{aligned}
 \|A_++A_-\|_{\mathrm{HS}}^2
 -\|A_+-A_-\|_{\mathrm{HS}}^2
 =4\operatorname{tr}(A_+A_-)\ge0.                        \tag{3.2}
 \end{aligned}
\]
Combining (3.2) with (2.5) proves the following statement.

### Theorem 3.1 (no vanishing normalized \(L^2\) DC buffer)

Under the assumptions of Section 1, every representation (3.1) with
\(D^2k_\pm\in L^2(\rho)\) obeys
\[
 \frac1\lambda
 \|D^2k_++D^2k_-\|_{L^2(\rho;\mathrm{HS})}
 \ge\sqrt{2\bigl(1-\lambda^2|a|^4\bigr)}
 \ge\sqrt{2(1-\lambda^2)}.                               \tag{3.3}
\]
If no such global DC representation exists, the left side is understood
as \(+\infty\).

Equivalently, define the polynomial-quotiented DC curvature
\[
 \mathfrak B_2(h):=
 \inf_{\substack{P\in\mathcal P_2,\ k_\pm\text{ convex}\\
                   h-P=k_+-k_-}}
 \frac1\lambda
 \|D^2k_++D^2k_-\|_{L^2(\rho;\mathrm{HS})}.              \tag{3.4}
\]
Then
\[
 \mathfrak B_2(h)\ge\sqrt{2(1-\lambda^2)}.               \tag{3.5}
\]

This is already true for the relaxed algebraic problem in which
\(A_\pm\) are merely PSD matrix fields and are not required to be
compatible Hessians.  For a symmetric matrix \(H\), the optimal
pointwise algebraic split is
\[
 H=H_+-H_-,\qquad H_++H_-=|H|,
\]
and \(\||H|\|_{\mathrm{HS}}=\|H\|_{\mathrm{HS}}\).  Hessian
compatibility can only increase the cost.

For nonsmooth convex \(k_\pm\), their distributional Hessians are PSD
matrix-valued measures.  Since the difference equals the smooth Hessian
of \(h-P\), their singular parts must agree.  That common positive
singular part adds to the buffer and cannot help.  The same lower bound
applies to the absolutely continuous parts whenever the displayed
\(L^2\) norm is finite; otherwise the DC cost is infinite.

### Approximate decompositions

The estimate is stable in the norm relevant to a Hessian multiplier.  If
\[
 h-P=k_+-k_-+r,
 \qquad D^2r\in L^2(\rho),                                \tag{3.6}
\]
then
\[
 \frac{\|D^2k_++D^2k_-\|_2}{\lambda}
 \ge
 \sqrt{2(1-\lambda^2|a|^4)}
 -\frac{\|D^2r\|_2}{\lambda}.                            \tag{3.7}
\]
Thus an \(o(\lambda)\) Hessian error does not change the obstruction.

Approximation only in \(L^2(\rho)\) is different.  It does not control
the Hessian distribution and hence does not control pairing with a KKT
matrix measure.  On compact sets, smooth functions can be approximated
in value norm by functions having arbitrarily large narrow curvature
spikes.  An \(L^2\)-only DC approximation therefore cannot supply the
normal-cone estimate needed in Target A.

## 4. The normalized KKT action is not small for free

The preceding curvature bound is a geometric obstruction.  There is
also an exact action obstruction that explains why the proposed route is
circular unless one obtains new information on the multiplier.

### 4.1 The polynomial residual

Let \(U=\operatorname{span}\{x_1,\ldots,x_n\}\subset L^2(\mu)\).  Since
\(f\) is centered, all terms in a degree-at-most-two polynomial depending
only on \(x\), or only on \(y\), are orthogonal to \(h=f\otimes f\).
Only bilinear terms survive.  Hence
\[
 \Pi_{\mathcal P_2}h=(a\cdot x)(a\cdot y),                 \tag{4.1}
\]
and, with \(q=h-\Pi_{\mathcal P_2}h\),
\[
 \|q\|_2^2=1-|a|^4.                                      \tag{4.2}
\]

### 4.2 Action barrier

Suppose a normalized PSD potential multiplier \(N\) realizes the mixed
source modulo the polynomial gauge: for some \(p\in\mathcal P_2\),
\[
 \int W(p-h)\,d\rho=-\int D^2W:N.                          \tag{4.3}
\]
for every admissible smooth test \(W\).  This is the normalized
\(\log\lambda\) version of the potential KKT identity one would need for
the mixed tensor entry.

If \(q=k_+-k_-\), set
\[
 \mathcal A_N(k):=\int D^2k:N\ge0
\]
for convex \(k\).  Then
\[
 \begin{aligned}
 \mathcal A_N(k_+)+\mathcal A_N(k_-)
 &\ge
 \left|\mathcal A_N(k_+)-\mathcal A_N(k_-)\right|\\
 &=\left|\int D^2q:N\right|\\
 &=\left|\int q(p-h)\,d\rho\right|
 =\int q^2\,d\rho.
                                                               \tag{4.4}
 \end{aligned}
\]
Therefore
\[
 \boxed{
 \mathcal A_N(k_++k_-)
 \ge1-|a|^4.
 }                                                           \tag{4.5}
\]

For an approximate value decomposition
\[
 q=k_+-k_-+r,
\]
the same argument gives the general estimate
\[
 \mathcal A_N(k_++k_-)
 \ge \|q\|_2^2-\|r\|_2\|p-h\|_2.                         \tag{4.6}
\]
Indeed, apply (4.3) to \(q-r\) and use Cauchy--Schwarz.  In the
orthogonal quotient gauge \(p=\Pi_{\mathcal P_2}h\), this becomes
\[
 \mathcal A_N(k_++k_-)
 \ge \|q\|_2^2-\|r\|_2\|q\|_2.                           \tag{4.7}
\]

Equations (4.5)--(4.6) have two implications.

1. If the mixed residual is bounded away from zero, the normalized DC
   buffer action is bounded away from zero for every decomposition.
2. A theorem proving \(\mathcal A_N(k_++k_-)=o(1)\) would immediately
   prove \(\operatorname{dist}(h,\mathcal P_2)=o(1)\).  The desired
   action bound is therefore essentially the mixed conclusion itself;
   Bochner does not provide it.

This does not negate the logical usefulness of an action estimate.  It
identifies exactly what new coercivity such an estimate would have to
contain.

### 4.3 Why flat complementarity permits a large action

For the original \(\lambda\)-objective, write its multiplier as \(M\).
After dividing the stationarity equation by \(\lambda\), the normalized
multiplier is \(N=M/\lambda\).  Complementarity says
\[
 \int D^2V:N=0.                                           \tag{4.8}
\]
If \(D^2V=0\) on an open set, (4.8) imposes no size bound on a PSD
multiplier supported there.

This can be seen directly from any DC buffer.  Let
\(S=D^2k_++D^2k_-\succeq0\) be nonzero and square-integrable on a flat
region.  Then
\[
 N_S:=\frac{S}{\|S\|_2^2}\,\rho                           \tag{4.9}
\]
is a finite PSD matrix measure, is complementary to a flat potential,
and satisfies
\[
 \int S:N_S=1.                                            \tag{4.10}
\]
If a proposed construction had \(\|S\|_2=O(\lambda)\), the dual norm of
\(N_S\) would be of order at least \(1/\lambda\).  Nothing in positivity
or flat complementarity excludes this scaling.

The measure in (4.9) is not asserted to be the multiplier of an actual
spectral extremizer.  It proves the precise functional-analytic point:
the normal cone contains directions that turn an \(O(\lambda)\) Hessian
buffer into an order-one normalized action.  To rule them out one must
use more than the conic KKT conditions.

The same issue is worse for singular multipliers.  If normalized PSD
measures are controlled only in total variation, their dual Hessian norm
is pointwise \(L^\infty\), not \(L^2\).  Bochner gives no pointwise
Hessian bound at all.

## 5. Unbounded support and honest pointwise convexity

An \(L^2\) lower bound on a Hessian is not a pointwise convexity
statement.  Conversely, an \(L^2\) upper bound does not give a global DC
decomposition.

### 5.1 What is always available locally

If \(u\in C^2\) on a compact convex set \(K\), then
\[
 L_K:=\sup_K\bigl(-\lambda_{\min}(D^2u)\bigr)_+<\infty.
\]
On \(K\), one has the semiconvex decomposition
\[
 u=\left(u+\frac{L_K}{2}|z|^2\right)
   -\frac{L_K}{2}|z|^2.                                  \tag{5.1}
\]
Both terms are convex.  On all of \(\mathbb R^{2n}\), the same formula
works only if the negative Hessian has a global pointwise lower bound.
Bochner's \(L^2\) estimate does not give such a bound.

Moreover, the buffer Hessian in (5.1) contains \(2L_KI\).  Even if
\(\|D^2u\|_2\) is small, \(L_K\) can be arbitrarily large because of a
narrow curvature spike.  Thus a semiconvexity constant cannot be
estimated from (1.5) alone, uniformly in dimension or support size.

### 5.2 Compact truncation

Assume now that \(f\in W^{2,2}(\mu)\) and is smooth locally.  Choose
standard cutoffs \(\chi_R\) with
\[
 \chi_R=1\text{ on }B_R,\qquad
 \chi_R=0\text{ outside }B_{2R},\qquad
 |D^j\chi_R|\le C_jR^{-j},\quad j=1,2.
\]
Put
\[
 f_R=\chi_Rf,
 \qquad h_R(x,y)=f_R(x)f_R(y).
\]
Tail convergence of \(f,\nabla f,D^2f\), together with the derivative
bounds on \(\chi_R\), gives
\[
 f_R\longrightarrow f\quad\text{in }W^{2,2}(\mu),
 \qquad
 h_R\longrightarrow h\quad\text{in }W^{2,2}(\rho).       \tag{5.2}
\]
Each \(h_R-P\) is smooth and compactly supported up to the polynomial
part.  Its Hessian is globally bounded, so (5.1) supplies an honest
global DC decomposition.

This does not make the buffer quantitative.  If
\[
 e_R=\int|\nabla f_R|^2\,d\mu,
 \qquad m_R=\int\nabla f_R\,d\mu,
\]
then the same mixed-block projection calculation gives, for every
quadratic \(P_R\),
\[
 \|D^2(h_R-P_R)\|_2^2
 \ge2(e_R^2-|m_R|^4).                                    \tag{5.3}
\]
By (5.2), \(e_R\to\lambda\) and \(m_R\to m=\lambda a\).
Hence
\[
 \liminf_{R\to\infty}
 \inf_{P_R\in\mathcal P_2}
 \frac{\|D^2(h_R-P_R)\|_2}{\lambda}
 \ge\sqrt{2(1-\lambda^2|a|^4)}.                          \tag{5.4}
\]
Every genuine DC buffer for the truncation has at least this much
normalized \(L^2\) curvature.  Meanwhile its pointwise convexifying
constant may diverge with \(R\).  Thus truncation resolves existence but
not smallness.

## 6. A smooth isotropic stress test

An actual sequence of isotropic log-concave measures with
\(\lambda_1\to0\) would itself disprove KLS, so constructing such a
sequence cannot be used as an auxiliary counterexample without solving
the main problem negatively.  The lower bound (3.3) is instead
conditional and universal: if such a sequence existed, its normalized
DC curvature would converge to at least \(\sqrt2\).

There is nevertheless a concrete smooth family showing that the
obstruction is not an artifact of rough convex bodies.  For integers
\(r\ge2\), let
\[
 d\mu_r(x)=Z_r^{-1}
 \exp\left[-\left(\frac{x}{s_r}\right)^{2r}\right]dx,    \tag{6.1}
\]
where \(s_r\) is chosen so that \(\operatorname{Var}(\mu_r)=1\).
These are smooth, strictly positive, isotropic log-concave densities on
\(\mathbb R\), with confining smooth convex potentials.  As
\(r\to\infty\), they converge to the uniform probability on
\([-\sqrt3,\sqrt3]\).  The one-dimensional Dirichlet forms converge in
the Mosco sense to the Neumann form on this interval, so their first
eigenpairs converge to the interval first eigenpair.

For the limiting interval, write \(A=\sqrt3\).  The normalized first odd
eigenfunction and eigenvalue are
\[
 f_0(x)=\sqrt2\sin\left(\frac{\pi x}{2A}\right),
 \qquad
 \lambda_0=\frac{\pi^2}{12}.                              \tag{6.2}
\]
Its linear coefficient is
\[
 a_0=\int xf_0(x)\,d\mu_0(x)
 =\frac{4\sqrt6}{\pi^2}.                                 \tag{6.3}
\]
Therefore the sharp mixed-block contribution to (3.3) tends to
\[
 \sqrt{2(1-\lambda_0^2a_0^4)}
 =\sqrt{2\left(1-\frac{64}{\pi^4}\right)}
 \approx0.8282.                                          \tag{6.4}
\]
This is a numerical order-one normalized DC cost in a limit of smooth
full-support log-concave states.

The example also separates value closeness from curvature cost.  Here
\[
 \operatorname{dist}_{L^2(\mu_0\otimes\mu_0)}
 (f_0\otimes f_0,\mathcal P_2)^2
 =1-a_0^4
 =1-\frac{9216}{\pi^8}
 \approx0.02872,                                         \tag{6.5}
\]
so the tensor is already rather close to a bilinear polynomial in value
norm, while its normalized DC curvature remains about \(0.83\).  This is
exactly why an \(H^2\) buffer estimate and the desired \(L^2\) mixed
rigidity estimate are not interchangeable.

## 7. Consequences for Target A

The DC-buffer proposal yields three audited conclusions.

### 7.1 What Bochner really gives

Bochner gives the correct absolute Hessian scale
\[
 \|D^2h\|_2=O(\lambda).
\]
It does not give a buffer of size \(o(\lambda)\).  Modulo every quadratic,
the total positive DC curvature divided by \(\lambda\) stays bounded
below by a numerical constant as \(\lambda\to0\).

### 7.2 What would still suffice

An \(O(\lambda)\) buffer could still have \(o(1)\) action if the actual
normalized KKT multiplier were uniformly bounded in a dual norm, or if
it had a special cancellation against the mixed Hessian.  Neither fact
follows from positivity or complementarity.  Flat potentials admit
normalized PSD multipliers with reciprocal-\(\lambda\) size, and
singular multipliers require pointwise rather than \(L^2\) Hessian
control.

Thus (3.3) does not by itself prove that the actual KKT action is large.
It proves that small action cannot be obtained from a generic
small-curvature argument.

This analysis also grants the DC proposal its most favorable spectral
premise: that a normalized PSD normal functional carrying the mixed
source has already been identified.  The independent nonsmooth issue
that the scalar spectral multiplier \(\Theta\) may discard the mixed
eigenspace entry remains.  A DC curvature estimate does not select that
entry.

### 7.3 The exact remaining statement

If the normalized mixed KKT identity (4.3) is present, positivity gives
the lower bound (4.5).  Therefore an \(o(1)\) buffer-action theorem is
equivalent in strength to forcing
\[
 \|h-\Pi_{\mathcal P_2}h\|_2=o(1).
\]
This is the desired mixed full-matrix rigidity, not a preliminary
consequence of the Bochner estimate.

Accordingly, the quantitative DC route does not close Target A on its
own.  A successful proof must add at least one genuinely new ingredient:

* a dimension-free bound on the normalized KKT multiplier on flat
  potential directions;
* a structural theorem forcing that multiplier to annihilate most of
  the mixed Hessian;
* a regularization whose complementary multiplier remains controlled in
  a dual norm strong enough for (3.7); or
* a different two-sided variation that bypasses the convex-potential
  normal cone.

Without such input, pointwise convexity, unbounded support, and flat
complementarity all preserve an order-one normalized obstruction.
