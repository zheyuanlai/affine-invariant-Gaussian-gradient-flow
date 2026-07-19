# Cyclic binary dilation: measurable laminar closure and a gap-charging potential

## 0. Scope and verdict

Let \(K\subset\mathbb R^d\) be a convex body whose uniform probability is
centered and isotropic.  Let \(R\subset K\) have volume \(|K|/2\), and let
\(b:R\to\mathbb R^d\) be locally constant on finitely many measurable
plateaux.  Assume that

\[
 K=(I-b)(R)\mathbin{\dot\cup}(I+b)(R)                  \tag{0.1}
\]

modulo null sets, with the two branch maps injective and volume preserving.
Assume, in addition, actual cyclic firm realizability: \(b\) is the
restriction to \(R\) of the gradient of a convex \(C^1\) function
\(\psi:\mathbb R^d\to\mathbb R\) whose gradient is \(1\)-Lipschitz.  Thus

\[
 0\preceq D^2\psi\preceq I
 \quad\hbox{in the Alexandrov sense},                  \tag{0.2}
\]

and \(b=\nabla\psi\) on \(R\).  These are the convex-displacement subbranch
of the finite zero-strain data produced by cyclic Brenier
piecewise-translation maps.  Actual Brenier cyclicity alone gives a
nonexpansive midpoint gradient; it need not give the firm inequality or a
convex midpoint potential.  Section 5.1 audits this distinction for a
general balanced Brenier map.

Write \(Z\) for a uniform point of \(R\), \(B=b(Z)\), and let
\(\varepsilon\) be a fair sign independent of \(Z\).  Then

\[
 U=Z+\varepsilon B
\]

is uniform on \(K\).  In particular,

\[
 \mathbb EZ=0,\qquad
 I=\mathbb E(ZZ^T)+\mathbb E(BB^T),\qquad
 \mathbb E(BB^T)\preceq I,\qquad |\mathbb EB|^2\le1.  \tag{0.3}
\]

The desired cyclic binary-dilation theorem is

\[
 \mathbb E|B|^2\le C.                                  \tag{0.4}
\]

This note does not prove (0.4) in complete generality.  It gives the
following rigorous advances and obstructions.

1. The global-corridor theorem remains true for **arbitrary measurable**
   cores, not only convex cores.  If a half-volume set can be enlarged
   through a sequential family of genuine Minkowski segment extrusions,
   then the total corridor length is at most \(\sqrt{24}\), and every
   label in the resulting zonotope satisfies
   \(\mathbb E|B|^2\le13\).  A quantitative version for a core of volume
   \(\alpha|K|\) is proved, as is a laminar-forest version.

2. A finite polyhedral complex whose midpoint set has \(q\) positive-gap
   connected components obeys
   \(\mathbb E|B|^2\le q\).  In particular, a connected midpoint set forces
   a single label and has energy at most one.  Distinct labels can occur
   only across genuine physical gaps.

3. A homothetic radial class has a dimension-free bound.  If
   \(R\subset(1-\delta)K\) and
   \(b(R)\subset \frac{\delta}{2}(K-K)\), then
   \[
      \mathbb E|B|^2\le3(\log2)^2.
   \]
   The proof tracks \(d\) explicitly and uses the sharp
   centroid--covariance containment of a convex body.

4. There is an exact global functional which sees all physical gaps.  If
   \(K_\pm=(I\pm b)(R)\) and
   \(\sigma=\mathbf1_{K_+}-\mathbf1_{K_-}\), then
   \[
    \boxed{\quad
      {3\over4}\mathbb E_R|b|^2
      \le \mathbb E_K[\sigma\psi]
      \le {5\over4}\mathbb E_R|b|^2 .
    \quad}                                             \tag{0.5}
   \]
   This is map-adaptive, uses the actual endpoint partition, and never
   discards physical gaps.  A dimension-free variance estimate for this
   special convex transport potential would close the firm zero-strain
   theorem.  For a completely general balanced Brenier map there is an
   exact analogue with a \(1\)-smooth, generally nonconvex Cayley
   potential.  Poincare restricted to all of those general Cayley
   potentials is quantitatively equivalent to KLS.  Restricting further
   to the convex members is not known to be equivalent and cannot be
   justified by an entropy-strain perturbation.

5. The most natural symmetric mixed-volume repair is false.  The symmetric
   Jacobian excess
   \[
    \det(I+D^2\psi)+\det(I-D^2\psi)-2
   \]
   contains only even elementary symmetric functions and vanishes for
   rank-one Hessians.  An explicit actual cyclic Brenier binary dilation
   on every isotropic cube has rank-one gap Hessian, zero symmetric
   Jacobian excess, zero source-wall charge, but
   \[
      \mathbb E|B|^2={1\over4},
      \qquad \operatorname {tr}\operatorname {Cov}B={1\over54}.
   \]
   Thus one must use an oriented, map-adaptive charge; averaging the two
   branches destroys the first-order Hessian term which records rank-one
   physical gaps.

The remaining obstruction is global and precise.  Local facet or gap
budgets on small components scale like
\(|v|\,|A|^{3/2}\).  That scaling permits \(d\) rare orthogonal labels of
energy \(d\) while spending only constant formal local volume.  Actual
cyclic complementarity must therefore manufacture a common large core,
an overlap charge, or an equivalent global coupling.  None of those three
conclusions is inferred from local wall data here.

## 1. Binary identities and the hypothesis hierarchy

For a Borel set \(A\subset K\), branchwise change of variables in (0.1)
gives

\[
\begin{aligned}
 \mathbb P(U\in A)
 &= {1\over2}{2\over|K|}
    \left(|(I-b)^{-1}(A)|+|(I+b)^{-1}(A)|\right)\\
 &= {|A|\over|K|}.
\end{aligned}                                         \tag{1.1}
\]

Also \(\mathbb E[U\mid Z]=Z\), so this is a binary martingale dilation.
Taking its first and second moments proves (0.3).  Put

\[
 M=\mathbb E(BB^T),\qquad m=\mathbb EB.
\]

Then

\[
 M\preceq I,\qquad |m|^2\le1,\qquad
 \operatorname {Cov}B\preceq I.                       \tag{1.2}
\]

Consequently, if the centered labels span an \(r\)-dimensional affine
space,

\[
 \mathbb E|B|^2
 =|m|^2+\operatorname {tr}\operatorname {Cov}B
 \le 1+r.                                              \tag{1.3}
\]

This elementary rank bound will be used below, but it cannot decide a
growing-rank complex.

The cyclic hypothesis is stronger than the pairwise inequality

\[
 (b(z)-b(z'))\cdot(z-z')\ge |b(z)-b(z')|^2.            \tag{1.4}
\]

In particular, (1.4) implies

\[
 |b(z)-b(z')|\le |z-z'|.                               \tag{1.5}
\]

The arguments in Sections 2--4 use only consequences explicitly stated
there.  Section 5 uses the full cyclic extension \(b=\nabla\psi\).  The
counterexample in Section 6 is not merely pairwise firm: it is produced by
an increasing one-dimensional Brenier map and then tensored with identity
maps.

## 2. A measurable laminar extrusion theorem

The earlier corridor argument was stated for convex cores.  Convexity is
unnecessary.  What matters is a genuine global extrusion inclusion.
All sets in this section are taken Borel; their projections and sums with
compact segments are analytic and hence Lebesgue measurable.  The same
proof applies to any completed-measurable data for which the displayed
Minkowski sums are measurable.

### 2.1 Projection thickness

For a unit vector \(u\), write \(P_u=P_{u^\perp}\).  Fibers of \(K\)
parallel to \(u\) are intervals.  Let \(L(y)\) be their lengths.  Since
the integral of \(t^2\) on any interval of length \(L\) is at least
\(L^3/12\), isotropy in direction \(u\) gives

\[
 {1\over12}\int_{P_uK}L(y)^3\,dy
 \le \int_K\langle x,u\rangle^2\,dx
 =|K|.                                                 \tag{2.1}
\]

If \(A\subset K\) is measurable and \(Q=P_uA\), then

\[
\begin{aligned}
 |A|
 &\le\int_Q L(y)\,dy\\
 &\le\left(\int L(y)^3dy\right)^{1/3}|Q|^{2/3}.
\end{aligned}
\]

Thus

\[
 \boxed{\quad
 |P_uA|\ge {|A|^{3/2}\over\sqrt{12|K|}}.
 \quad}                                                \tag{2.2}
\]

This proof uses neither convexity of \(A\) nor an isoperimetric inequality.

### 2.2 Measurable one-segment growth

Let \(v=\ell u\), \(\ell\ge0\).  On every nonempty \(u\)-fiber,
the one-dimensional Brunn--Minkowski inequality says

\[
 |A_y+[0,\ell]|\ge |A_y|+\ell .
\]

After integrating over \(P_uA\),

\[
 \boxed{\quad
 |A+[0,v]|\ge |A|+|v|\,|P_uA|.
 \quad}                                                \tag{2.3}
\]

For convex \(A\), (2.3) is the exact one-segment Steiner formula.  For
arbitrary measurable \(A\), the inequality is all that is needed.

### 2.3 Sequential measurable corridors

**Theorem 2.1 (measurable global-corridor bound).**  Let
\(0<\alpha<1\).  Suppose measurable sets

\[
 A_0\subset A_1\subset\cdots\subset A_N\subset K       \tag{2.4}
\]

satisfy

\[
 |A_0|\ge\alpha|K|,
 \qquad
 A_{j-1}+[0,v_j]\subset A_j
 \quad(1\le j\le N).                                  \tag{2.5}
\]

Then

\[
 \boxed{\quad
 \sum_{j=1}^N|v_j|
 \le L(\alpha):={\sqrt{12}(1-\alpha)\over\alpha^{3/2}}.
 \quad}                                                \tag{2.6}
\]

In particular, \(L(1/2)=\sqrt{24}\).

**Proof.**  Equations (2.2)--(2.5), and
\(|A_{j-1}|\ge\alpha|K|\), imply

\[
\begin{aligned}
 |A_j|-|A_{j-1}|
 &\ge |v_j|\,|P_{v_j^\perp}A_{j-1}|\\
 &\ge {|v_j|\alpha^{3/2}\over\sqrt{12}}\,|K|.
\end{aligned}                                         \tag{2.7}
\]

Sum in \(j\).  The left side telescopes and is at most
\((1-\alpha)|K|\), proving (2.6). \(\square\)

No angle condition, commutation assumption, boundary regularity, or
polyhedral structure occurs in this theorem.

**Corollary 2.2 (measurable zonotopal laminar class).**  In addition to
the hypotheses of Theorem 2.1, suppose the label range is contained in

\[
 b_0+\sum_{j=1}^N[0,v_j].                              \tag{2.8}
\]

Then

\[
 \operatorname {tr}\operatorname {Cov}B
 \le {L(\alpha)^2\over2}
 ={6(1-\alpha)^2\over\alpha^3},                        \tag{2.9}
\]

and

\[
 \boxed{\quad
 \mathbb E|B|^2
 \le1+{6(1-\alpha)^2\over\alpha^3}.
 \quad}                                                \tag{2.10}
\]

For \(\alpha=1/2\), the right side is \(13\).

**Proof.**  The zonotope in (2.8) has diameter at most
\(\sum_j|v_j|\).  If \(B'\) is an independent copy of \(B\), then

\[
 \operatorname {tr}\operatorname {Cov}B
 ={1\over2}\mathbb E|B-B'|^2
 \le {1\over2}\left(\sum_j|v_j|\right)^2.              \tag{2.11}
\]

Use Theorem 2.1 and \(|\mathbb EB|^2\le1\). \(\square\)

The gain over the convex-core statement is substantive for midpoint
complexes: a core may contain holes and physical target gaps.  The theorem
charges them correctly as long as the *whole measurable core* can actually
be extruded.  Calling a collection of unrelated local facet prisms a global
extrusion would still be invalid.

### 2.4 A laminar-forest version

Let \(S\) be a finite cluster variable measurable with respect to \(B\),
with probabilities \(p_s\).  Suppose that, conditional on \(S=s\), the
label range is contained in a zonotope supplied by a measurable corridor
chain with starting fraction \(\alpha_s\).  Let

\[
 r=\dim\operatorname {span}\{
       \mathbb E[B\mid S=s]-\mathbb EB:p_s>0\}.
\]

The law of total covariance and Corollary 2.2 give

\[
\begin{aligned}
 \operatorname {tr}\operatorname {Cov}B
 &=\operatorname {tr}\operatorname {Cov}(\mathbb E[B\mid S])
   +\sum_s p_s\operatorname {tr}\operatorname {Cov}(B\mid S=s)\\
 &\le r+\sum_s p_s\,{6(1-\alpha_s)^2\over\alpha_s^3}.
\end{aligned}                                         \tag{2.12}
\]

Hence

\[
 \boxed{\quad
 \mathbb E|B|^2
 \le1+r+\sum_s p_s\,{6(1-\alpha_s)^2\over\alpha_s^3}.
 \quad}                                                \tag{2.13}
\]

This is a rigorous laminar-forest theorem.  It is dimension free whenever
the between-cluster rank and the displayed weighted core penalty are
universal.  It also exposes the exact defect of a forest made of many tiny
cores: the power \(\alpha_s^{-3}\) cannot be discarded.

## 3. Polyhedral connectivity forces physical gaps

Assume \(b\) takes finitely many values.  By (1.5), \(b\) is continuous on
the metric subspace \(R\).  A continuous map from a connected set into a
finite set is constant.  Therefore:

**Proposition 3.1 (component-rank bound).**  If \(R\) has \(q<\infty\)
connected components, then \(b\) takes at most \(q\) values and

\[
 \boxed{\quad
 \mathbb E|B|^2\le 1+\min(d,q-1).
 \quad}                                                \tag{3.1}
\]

In particular, if \(R\) is connected, then \(b\) is constant and

\[
 \mathbb E|B|^2=|\mathbb EB|^2\le1.                   \tag{3.2}
\]

The same conclusion holds with \(q\) equal to the number of components of
the graph obtained by joining two plateau cells whenever their set
distance is zero.  Indeed, if
\(\operatorname {dist}(R_i,R_j)=0\), then (1.5), applied to sequences
approaching a common limit, forces \(b_i=b_j\).

Thus a nonconstant finite label field is impossible without actual
positive-width gaps in the midpoint set.  A source-wall calculation which
only records adjacent occupied cells necessarily misses the only places
where labels can change.

Proposition 3.1 is sharp at the level of the covariance information
\(\operatorname {Cov}B\preceq I\): \(q\) points can have centered affine
rank \(q-1\).  Improving (3.1) for large \(q\) requires the placement and
cyclic ordering of the physical components, not merely their count.

## 4. A radial and homothetic class

The following is the sharp Kannan--Lovasz--Simonovits
inertia-ellipsoid inclusion; its one-dimensional localization proof is
recalled to fix the normalization.

**Lemma 4.1 (centroid--covariance containment).**  If \(K\) is centered
and isotropic, then

\[
 K\subset \sqrt{d(d+2)}\,B_2^d.                       \tag{4.1}
\]

One proof fixes \(u\in S^{d-1}\).  The density of
\(\langle X,u\rangle\), \(X\) uniform on \(K\), has an
\((d-1)\)-concave root.  Among such densities with a prescribed upper
endpoint \(h\) and mean zero, the smallest second moment is attained by
the cone density

\[
 c\,(h-t)^{d-1}\mathbf1_{[-h/d,h]}(t),
\]

whose variance is \(h^2/[d(d+2)]\).  Since the variance is one,
\(h\le\sqrt{d(d+2)}\).  Apply the same statement to every \(u\).
Equality is attained in vertex directions of the isotropic regular
simplex.

**Theorem 4.2 (homothetic radial binary dilations).**  In the setup of
Section 0, suppose that for some \(0\le\delta<1\),

\[
 R\subset(1-\delta)K,
 \qquad
 b(R)\subset {\delta\over2}(K-K).                     \tag{4.2}
\]

Then

\[
 \boxed{\quad
 \mathbb E|B|^2\le3(\log2)^2.
 \quad}                                                \tag{4.3}
\]

**Proof.**  The half-volume condition and the first inclusion give

\[
 {1\over2}|K|=|R|
 \le(1-\delta)^d|K|.
\]

Therefore

\[
 \delta\le1-2^{-1/d}\le{\log2\over d}.                 \tag{4.4}
\]

By Lemma 4.1, any \(x,y\in K\) satisfy

\[
 \left|{x-y\over2}\right|\le\sqrt{d(d+2)}.
\]

The second inclusion in (4.2) consequently gives

\[
 |B|^2\le\delta^2d(d+2)
 \le(\log2)^2\left(1+{2\over d}\right)
 \le3(\log2)^2.                                       \tag{4.5}
\]

Average (4.5). \(\square\)

The theorem covers a broad homothetic inward-core/outward-label regime.
Neither condition in (4.2) follows merely from \(Z\pm B\in K\): a central
point can lie on a long chord.  The conditions encode the radial mechanism
which makes the crosspolytope and related fans benign.

Two additional radial polyhedral models have exact dimension-free
calculations.

* In the isotropic crosspolytope
  \(K=\{x:\|x\|_1\le R_d\}\), a symmetric outward coordinate fan with
  displacement length \(t\) draws all source mass from
  \((1-t/R_d)K\).  Half volume forces
  \[
     (1-t/R_d)^d\ge{1\over2},
     \qquad
     t\le {R_d\log2\over d}.
  \]
  With \(R_d=\sqrt{(d+1)(d+2)/2}\), its midpoint-label energy is at most
  \((\log2)^2\).

* In the isotropic regular simplex, labels
  \(b_i=\lambda v_i\) on \(d+1\) equally weighted vertex pieces can be
  feasible only where at most one barycentric coordinate is smaller than
  \(c=2\lambda/(d+1)\).  For a uniform Dirichlet
  \((1,\ldots,1)\) vector, the exact probability of this event is
  \[
   (d+1)\left(1-{2\lambda d\over d+1}\right)_+^d
   -d(1-2\lambda)_+^d.                                \tag{4.6}
  \]
  For \(d\ge2\), feasibility first forces \(\lambda<1/2\), so the
  positive parts may then be removed.  Requiring (4.6) to be at least
  \(1/2\) gives
  \(\lambda d<3\), and hence
  \[
     \mathbb E|B|^2
     =\lambda^2d(d+2)\le18.                            \tag{4.7}
  \]

The exact overlap term in (4.6) is essential; a union bound loses
\(\log d\).  These fan calculations and Theorem 4.2 are radial mechanisms,
whereas Theorem 2.1 is a laminar mixed-volume mechanism.

## 5. An exact potential functional which charges every physical gap

Let

\[
 K_+=(I+b)(R),\qquad K_-=(I-b)(R),
\qquad
 \sigma=\mathbf1_{K_+}-\mathbf1_{K_-}.                \tag{5.1}
\]

Both sets have volume \(|K|/2\), and they partition \(K\).

**Theorem 5.1 (signed transport-potential certificate).**  Let
\(\psi\) be the convex \(1\)-smooth extension with
\(\nabla\psi=b\) on \(R\).  Then

\[
 \boxed{\quad
 {3\over4}\mathbb E_R|b|^2
 \le \mathbb E_K[\sigma\psi]
 \le {5\over4}\mathbb E_R|b|^2.
 \quad}                                                \tag{5.2}
\]

**Proof.**  Convexity and the \(1\)-smooth descent lemma imply, for
\(z\in R\) and \(h=b(z)=\nabla\psi(z)\),

\[
\begin{array}{rcl}
 \psi(z+h)&\ge&\psi(z)+|h|^2,\\
 \psi(z+h)&\le&\psi(z)+\frac32|h|^2,\\
 \psi(z-h)&\ge&\psi(z)-|h|^2,\\
 \psi(z-h)&\le&\psi(z)-\frac12|h|^2.
\end{array}                                           \tag{5.3}
\]

Consequently

\[
 {3\over4}|h|^2
 \le {\,\psi(z+h)-\psi(z-h)\over2}
 \le {5\over4}|h|^2.                                  \tag{5.4}
\]

Each branch pushes the uniform law on \(R\) to the uniform conditional law
on \(K_\pm\).  Hence

\[
\begin{aligned}
 \mathbb E_K[\sigma\psi]
 &= {1\over |K|}
   \left(\int_{K_+}\psi-\int_{K_-}\psi\right)\\
 &= {1\over2}\mathbb E_R[
      \psi(Z+B)-\psi(Z-B)].
\end{aligned}                                         \tag{5.5}
\]

Average (5.4). \(\square\)

This functional is insensitive to the choice of additive constant in
\(\psi\).  It uses values on the actual endpoint partition and therefore
cannot overlook a physical gap separating midpoint components.

There is also an unsigned convex-order identity.  For every convex
\(\varphi\),

\[
\begin{aligned}
 \mathcal J_\varphi
 &:=\mathbb E_K\varphi-\mathbb E_R\varphi\\
 &=\mathbb E_R\left[
   {\varphi(Z+B)+\varphi(Z-B)\over2}-\varphi(Z)\right]
 \ge0.                                                 \tag{5.6}
\end{aligned}
\]

If \(\varphi\) is \(C^2\), then

\[
 \mathcal J_\varphi
 ={1\over2}\mathbb E_R\int_{-1}^1
 (1-|s|)\,
 B^T D^2\varphi(Z+sB)B\,ds.                            \tag{5.7}
\]

For \(\varphi=\psi\), \(1\)-smoothness gives only

\[
 0\le\mathcal J_\psi\le{1\over2}\mathbb E|B|^2.       \tag{5.8}
\]

Unlike the signed certificate (5.2), the unsigned certificate may vanish
for a constant nonzero label.

The signed certificate yields a concrete restricted analytic target.  By
the \(1\)-Lipschitz property of \(\nabla\psi\),

\[
 |\nabla\psi(z\pm b(z))|
 \le|\nabla\psi(z)|+|b(z)|=2|b(z)|.
\]

Changing variables on the two branches gives

\[
 \int_K|\nabla\psi|^2\,d\lambda_K
 \le4\mathbb E_R|b|^2.                                \tag{5.9}
\]

Since \(\mathbb E_K\sigma=0\), Cauchy--Schwarz and (5.2) show that a bound

\[
 \operatorname {Var}_{\lambda_K}\psi
 \le C_0\int_K|\nabla\psi|^2\,d\lambda_K               \tag{5.10}
\]

for this special class of cyclic binary-dilation potentials would imply

\[
 {3\over4}E
 \le\sqrt{\operatorname {Var}_K\psi}
 \le2\sqrt{C_0E},
 \qquad E=\mathbb E_R|b|^2,
\]

and therefore

\[
 E\le {64\over9}C_0.                                  \tag{5.11}
\]

Equation (5.10) is not assumed.  Full Poincare would prove it and would be
circular in the KLS task.  The point of (5.2) is narrower: it replaces a
vague request to "charge the gaps" by one exact map-adaptive scalar
functional and a precisely delimited restricted variance problem.

### 5.1 General balanced Brenier maps: exact certificate, but nonconvex

The preceding convex potential is special to the cyclically **firm**
subbranch.  There is nevertheless an exact signed certificate for every
balanced Brenier map.

Let \(\mu\) be any centered log-concave probability with finite second
moment, let \(E\) have \(\mu(E)=1/2\), and put

\[
 \mu_+=2\mathbf1_E\mu,\qquad
 \mu_-=2\mathbf1_{E^c}\mu.                             \tag{5.12}
\]

Let \(T=\nabla\phi\) be the Brenier map from \(\mu_+\) to \(\mu_-\).
For \(X\sim\mu_+\), write

\[
 Y=T(X),\qquad Z={X+Y\over2},\qquad B={Y-X\over2}.     \tag{5.13}
\]

Let

\[
 e_\phi(w)=\inf_x\left\{\phi(x)+{1\over2}|w-x|^2\right\}
\]

be the Moreau envelope, and define the Cayley midpoint potential

\[
 q(z)={1\over2}e_\phi(2z)-{1\over2}|z|^2.              \tag{5.14}
\]

The proximal optimality condition at \(2z=x+y\) gives

\[
 \operatorname {prox}_\phi(2z)=x,\qquad
 \nabla q(z)=y-z={y-x\over2}=B.                        \tag{5.15}
\]

The map \(\nabla q\) is \(1\)-Lipschitz.  Indeed, if
\(x=\operatorname {prox}_\phi(2z)\) and
\(x'=\operatorname {prox}_\phi(2z')\), monotonicity of
\(\partial\phi\) gives

\[
 |x-x'|^2\le2\langle x-x',z-z'\rangle,
\]

and therefore

\[
 |(z-x)-(z'-x')|^2\le|z-z'|^2.                        \tag{5.16}
\]

At a regular point where \(A=DT(X)\),

\[
 D^2q(Z)=(A-I)(A+I)^{-1}.                              \tag{5.17}
\]

Its eigenvalues lie in \([-1,1]\), but they need not be nonnegative.
Thus \(q\) is \(1\)-smooth and is generally **not convex**.

Put

\[
 \sigma=\mathbf1_{E^c}-\mathbf1_E.
\]

The two-sided Taylor remainder for a function with \(1\)-Lipschitz
gradient, applied at \(z\) with \(h=\nabla q(z)\), gives

\[
 {1\over2}|h|^2
 \le {q(z+h)-q(z-h)\over2}
 \le {3\over2}|h|^2.                                  \tag{5.18}
\]

Since \(z-h=x\), \(z+h=y\), and \(T_\#\mu_+=\mu_-\),

\[
 \mathbb E_\mu[\sigma q]
 ={1\over2}\mathbb E_{\mu_+}[q(Y)-q(X)].
\]

Consequently every balanced Brenier map satisfies the exact
dimension-free certificate

\[
 \boxed{\quad
 {1\over8}W_2^2(\mu_+,\mu_-)
 \le\mathbb E_\mu[\sigma q]
 \le {3\over8}W_2^2(\mu_+,\mu_-).
 \quad}                                                \tag{5.19}
\]

There is no positive-strain or entropy error in (5.19).  Cayley curvature
and physical-gap geometry change convexity, not the signed energy
comparison.

The endpoint gradient is controlled by the same displacement.  From
(5.15)--(5.16),

\[
 |\nabla q(X)|\le2|B|,
 \qquad
 |\nabla q(Y)|\le2|B|,
\]

and hence

\[
 \int|\nabla q|^2d\mu
 \le4\mathbb E|B|^2
 =W_2^2(\mu_+,\mu_-).                                  \tag{5.20}
\]

This identifies the exact logical strength of a restricted Poincare
estimate.  Define

\[
 C_{\rm Cay}(\mu)
 =\sup_{\mu(E)=1/2}
 { \operatorname {Var}_\mu q_E
  \over \int|\nabla q_E|^2d\mu},                       \tag{5.21}
\]

where \(q_E\) is (5.14), and zero denominators are omitted.  Ordinary
Poincare gives

\[
 C_{\rm Cay}(\mu)\le C_P(\mu).                         \tag{5.22}
\]

Conversely, (5.19), Cauchy--Schwarz, and (5.20) give, for every balanced
set,

\[
 {W_2^2\over8}
 \le\sqrt{\operatorname {Var}_\mu q_E}
 \le\sqrt{C_{\rm Cay}}\;W_2,
\]

so

\[
 \sup_{\mu(E)=1/2}W_2^2(\mu_+,\mu_-)
 \le64C_{\rm Cay}(\mu).                                \tag{5.23}
\]

The balanced-half transport implication

\[
 D_1(\mu)\le{1\over2}
 \sup_{\mu(E)=1/2}W_2(\mu_+,\mu_-)
\]

and E. Milman's equivalence for log-concave probabilities then imply

\[
 C_P(\mu)\le C\,C_{\rm Cay}(\mu)                       \tag{5.24}
\]

with a universal numerical \(C\).  Thus a dimension-free Poincare
inequality for **all Cayley midpoint potentials \(q_E\)** is
quantitatively equivalent to full KLS.  It is not a weaker target.

By contrast, an inequality only for the convex members of this class
controls only the firm subbranch.  There is no entropy-deficit
convexification which repairs this distinction.  The normalized
derivative strain does satisfy

\[
 \mathbb E_{\mu_+}
 \left\|(A-I)(A+I)^{-1}\right\|_{HS}^2
 \le2\{\log2-\operatorname {Ent}_\mu(\nu_{1/2})\},     \tag{5.25}
\]

but the negative curvature of \(q\) can lie entirely in a physical
midpoint gap and be invisible to the expectation in (5.25).

An exact zero-deficit example already occurs on \(K=[0,4]\).  Translation
to its barycenter and the one-dimensional isotropic rescaling preserve
all sign and zero-deficit assertions, so the unnormalized coordinates are
used for transparency.  Take

\[
 E=[0,1]\cup[3,4],\qquad E^c=[1,3].
\]

The increasing Brenier map is \(T(x)=x+1\) on \([0,1]\) and
\(T(x)=x-1\) on \([3,4]\).  It has \(T'=1\) at every source density
point, so both the matrix strain and the midpoint entropy deficit vanish.
The midpoint labels, however, are

\[
 b=+{1\over2}\quad\hbox{on }[1/2,3/2],
 \qquad
 b=-{1\over2}\quad\hbox{on }[5/2,7/2].
\]

On the physical midpoint gap \([3/2,5/2]\), the canonical Cayley extension
has

\[
 \nabla q(z)=2-z,\qquad D^2q=-1.                       \tag{5.26}
\]

Thus \(q\) is nonconvex even when the right side of (5.25) is zero.
This is not merely a defect of the canonical choice of \(q\).  For every
convex function \(r:\mathbb R\to\mathbb R\),

\[
 \mathbb E_{\mu_+}r
 ={1\over2}\int_0^1[r(u)+r(4-u)]\,du
 \ge {1\over2}\int_0^1[r(1+u)+r(3-u)]\,du
 =\mathbb E_{\mu_-}r.                                  \tag{5.26a}
\]

Indeed, both pairs in the integrand have midpoint \(2\), and the first
pair has the larger separation.  Hence

\[
 \mathbb E_\mu[\sigma r]\le0
\]

for every convex \(r\), whereas \(W_2^2(\mu_+,\mu_-)=1\) and the entropy
deficit is zero.  Therefore no inequality of the form

\[
 cW_2^2\le \mathbb E_\mu[\sigma r]+C\,
 \{\log2-\operatorname {Ent}_\mu(\nu_{1/2})\}
\]

can hold with \(c>0\) and a convex \(r\), even if \(r\) is allowed to
depend on the transport, **for this prescribed source-to-target
orientation**.  Reversing the transport replaces the Cayley potential by
its negative and happens to repair this one-dimensional concave example;
it does not provide a convexification of the original \(q\), and a
mixed-sign Cayley Hessian is repaired by neither orientation.  Thus a
canonical signed convex-potential extension of (5.19), with the same
orientation and an entropy-controlled error, is impossible.

The most direct convexification also displays an uncontrolled term.  The
function

\[
 r(z)={1\over2}\left(q(z)+{1\over2}|z|^2\right)
 ={1\over4}e_\phi(2z)                                  \tag{5.27}
\]

is convex with \(1\)-Lipschitz gradient, but

\[
 \mathbb E_\mu[\sigma r]
 ={1\over2}\mathbb E_\mu[\sigma q]
 +{1\over8}\left(
   \mathbb E_{\mu_-}|Y|^2-\mathbb E_{\mu_+}|X|^2
  \right).                                             \tag{5.28}
\]

The conditional second-moment imbalance in (5.28) is not controlled by
the entropy deficit; in general it can carry dimension dependence.
Therefore (5.27) does not transfer (5.19) to the convex class with a
universal positive-strain error.

## 6. Why symmetric mixed volume loses rank-one gaps

For \(0\preceq H\preceq I\),

\[
\begin{aligned}
 &\det(I+H)+\det(I-H)-2\\
 &\qquad
 =2\sum_{\substack{k\ge2\\k\ {\rm even}}}e_k(H),
\end{aligned}                                         \tag{6.1}
\]

where \(e_k\) is the \(k\)-th elementary symmetric polynomial of the
eigenvalues.  In particular,

\[
 \operatorname {rank}H\le1
 \quad\Longrightarrow\quad
 \det(I+H)+\det(I-H)-2=0.                              \tag{6.2}
\]

Thus adding the two branch Jacobian excesses cancels the trace term.  The
following actual family shows that this is fatal, not a formal concern.

### 6.1 An actual cyclic rank-one-gap family in every dimension

Put \(s=\sqrt3\) and

\[
 K_d=[-s,s]^d.                                         \tag{6.3}
\]

This cube is centered and isotropic.  In its first coordinate define

\[
\begin{array}{lll}
 P_1=[-s,-2/s],&
 Q_1=[-2/s,-1/s],&
 a_1=1/s,\\[2mm]
 P_2=[-1/s,1/s],&
 Q_2=[1/s,s],&
 a_2=2/s.
\end{array}                                            \tag{6.4}
\]

Translation by \(a_i\) sends \(P_i\) to \(Q_i\).  The source
\(P_1\cup P_2\) and target \(Q_1\cup Q_2\) are complementary half-length
subsets of \([-s,s]\).  Because \(a_1<a_2\) and \(P_1\) lies to the left
of \(P_2\), this is the increasing one-dimensional Brenier map.

The midpoint intervals and labels are

\[
\begin{array}{lll}
 R_1=[-5/(2s),-3/(2s)],& b_1=1/(2s),\\[1mm]
 R_2=[0,2/s],& b_2=1/s.
\end{array}                                            \tag{6.5}
\]

Let

\[
\begin{aligned}
 R&=(R_1\cup R_2)\times[-s,s]^{d-1},\\
 b(x)&=
 \begin{cases}
 b_1e_1,&x_1\in R_1,\\
 b_2e_1,&x_1\in R_2.
 \end{cases}
\end{aligned}                                         \tag{6.6}
\]

Tensoring the one-dimensional Brenier map with the identity on the last
\(d-1\) coordinates proves actual cyclic Brenier realizability.  Directly,

\[
 K_d=(I-b)(R)\mathbin{\dot\cup}(I+b)(R).               \tag{6.7}
\]

The two midpoint components have weights \(1/3\) and \(2/3\).  Therefore

\[
\begin{aligned}
 \mathbb E|B|^2
 &= {1\over3}{1\over4s^2}
   +{2\over3}{1\over s^2}
 ={1\over4},\\
 |\mathbb EB|^2
 &=\left({5\over6s}\right)^2={25\over108},\\
 \operatorname {tr}\operatorname {Cov}B
 &= {1\over4}-{25\over108}
 ={1\over54}.                                         \tag{6.8}
\end{aligned}
\]

A canonical cyclic firm extension is

\[
 \nabla\psi(x)=\beta(x_1)e_1,
\qquad
 \beta(t)=
 \begin{cases}
 1/(2s),&t\le-1/s,\\
 t+s/2,&-1/s\le t\le-1/(2s),\\
 1/s,&t\ge-1/(2s).
 \end{cases}                                          \tag{6.9}
\]

Thus

\[
 D^2\psi(x)
 =\mathbf1_{[-1/s,-1/(2s)]}(x_1)\,e_1e_1^T
 \quad\hbox{a.e.}                                     \tag{6.10}
\]

and has rank at most one everywhere.  Equations (6.1)--(6.2) give

\[
 \int_{K_d}
 \left[\det(I+D^2\psi)+\det(I-D^2\psi)-2\right]dx=0,   \tag{6.11}
\]

while (6.8) is strictly positive.

There is no source--source wall: the two occupied midpoint slabs are
separated by a physical gap.  Hence every local functional supported only
on interfaces of occupied source cells is also zero.  This proves, in every
dimension and within the exact cyclic class, that neither source-wall
volume nor symmetric branch mixed volume can control even the centered
label energy without an additive term.

The oriented determinant retains the missing trace:

\[
 \det(I+H)-1=\operatorname {tr}H
 \quad\hbox{when }\operatorname {rank}H\le1.           \tag{6.12}
\]

But an ambient inclusion for the oriented extension is false.  In (6.9),
points with \(2/s<x_1\le s\) lie in the physical gap \(K_d\setminus R\)
and satisfy

\[
 (I+\nabla\psi)_1(x)=x_1+1/s>s
\]

on a positive-volume terminal slab.  Hence

\[
 (I+\nabla\psi)(K_d\setminus R)
 \not\subset K_d\setminus(I+b)(R).                    \tag{6.13}
\]

This is the exact dichotomy:

* symmetrizing the two Jacobians gives an ambiently natural expression but
  cancels the rank-one trace charge;
* keeping one orientation sees the gap, but its canonical extension need
  not remain in the ambient body.

The signed potential in Theorem 5.1 retains orientation without extending
an image of all physical gaps.

## 7. A no-go theorem for map-independent convex-order tests

Identity (5.6) may tempt one to choose a single convex function
\(\varphi_K\), depending only on \(K\), whose Jensen gap controls every
possible binary displacement.  A globally strongly convex choice does
charge the energy:

\[
 D^2\varphi_K\succeq cI
 \quad\Longrightarrow\quad
 \mathcal J_{\varphi_K}\ge {c\over2}\mathbb E|B|^2.   \tag{7.1}
\]

It cannot have a dimension-free ambient budget.  On the isotropic cube
\([-s,s]^d\), for each coordinate \(j\) there is an actual constant-label
binary dilation

\[
 R^{(j)}
 =[-s,s]^{j-1}\times[-s/2,s/2]\times[-s,s]^{d-j},
 \qquad
 b^{(j)}={s\over2}e_j.                                \tag{7.2}
\]

The two translates by \(\pm b^{(j)}\) partition the cube, and the label is
cyclic firm.  Any map-independent function which is \(c\)-strongly convex
on every coordinate line satisfies, by conditioning successively on the
coordinates,

\[
 \mathbb E_K\varphi_K-\varphi_K(0)
 \ge {c\over2}\sum_{j=1}^d\mathbb EX_j^2
 ={cd\over2}.                                         \tag{7.3}
\]

Thus no fixed globally strongly convex scalar test has both a universal
energy charge in all realizable directions and an \(O(1)\) expectation
budget.  A successful convex-order functional must be map-adaptive or
anisotropic.  The transport potential \(\psi\) in Theorem 5.1 is precisely
map-adaptive.

This no-go statement does not exclude a subtler \(\varphi=\varphi_{K,b}\)
whose curvature follows the realized label covariance.  It only rules out
the most direct universal strong-convexity completion.

## 8. Why a local laminar budget cannot close by itself

The exponent in the measurable projection estimate (2.2) is sharp for the
information being used.  The following numerical model isolates the
failure.  It is deliberately not claimed to be Brenier realizable.

For \(1\le i\le d\), take formal component masses and corridor vectors

\[
 p_i=\alpha_i={1\over d},
 \qquad
 v_i=\sqrt d\,e_i.                                    \tag{8.1}
\]

Then the first-order local projection cost predicted by (2.2) is

\[
 |v_i|\alpha_i^{3/2}={1\over d},
 \qquad
 \sum_{i=1}^d |v_i|\alpha_i^{3/2}=1.                  \tag{8.2}
\]

Yet if \(B=v_I\) for a uniform random index \(I\), then

\[
 \mathbb E(BB^T)=I,\qquad
 |\mathbb EB|^2=1,\qquad
 \operatorname {tr}\operatorname {Cov}B=d-1,\qquad
 \mathbb E|B|^2=d.                                    \tag{8.3}
\]

Thus all covariance operator constraints in (0.3) and a constant total
formal local mixed-volume budget coexist with trace energy \(d\).

This calculation identifies what a valid local-to-global theorem must add.
It must prove, from actual cyclic complementarity, at least one of:

1. many small components share a common core of fixed volume;
2. the cross-corridors overlap in a way whose exact union charge is
   coercive, as in the regular-simplex calculation;
3. the rare orthogonal pattern (8.1) is incompatible with cyclic power
   cells and complementary endpoint tilings.

Replacing \(\alpha_i^{3/2}\) by \(\alpha_i\), or simply summing local facet
areas as if every component were a half-volume core, would assume the
missing conclusion.

## 9. Audit and remaining theorem

The dimension dependence of every positive result is explicit:

\[
\begin{array}{c|c}
\text{class}&\text{bound on }\mathbb E|B|^2\\ \hline
\text{one measurable half-core extrusion chain}&13\\
\text{connected midpoint set}&1\\
\text{\(q\)-component finite complex}&q\\
\text{homothetic radial class}&3(\log2)^2\\
\text{symmetric crosspolytope fan}&(\log2)^2\\
\text{symmetric regular-simplex fan}&18.
\end{array}                                           \tag{9.1}
\]

The measurable-chain proof tracks the only auxiliary parameter explicitly:

\[
 1+{6(1-\alpha)^2\over\alpha^3}.                      \tag{9.2}
\]

It degenerates when the starting core becomes small; no hidden limiting
argument suppresses that dependence.  The rank-one family of Section 6 is
valid for \(d=1\) and for every higher \(d\), is non-symmetric in its active
coordinate, has unboundedly many passive dimensions, and retains the full
physical endpoint gaps.

The exact unresolved statement can now be posed in either of two
nonlocal forms.

**Geometric form.**  Prove that every actual cyclic complementary
polyhedral complex admits a distribution of measurable extrusion chains
whose starting fractions and between-chain ranks make the right side of
(2.13) universal.

**Analytic form.**  Prove the restricted potential estimate (5.10), with a
universal constant, for the convex \(1\)-smooth potential generated by
every actual cyclic binary dilation.

This analytic form closes the firm zero-strain branch only.  If it is
enlarged to the nonconvex Cayley potentials of every balanced Brenier map,
then (5.21)--(5.24) show that it is quantitatively equivalent to full KLS.
The zero-deficit gap example (5.26) prevents an entropy-controlled passage
from the latter family to the convex family.

The geometric form cannot be deduced from local projection costs because
of (8.1)--(8.3).  The analytic form cannot be replaced by full Poincare,
which would be circular.  Symmetric determinant volume cannot bridge the
two forms because the actual family in Section 6 makes it identically zero.

Accordingly, the report proves the cyclic binary-dilation theorem on the
largest presently justified measurable laminar, bounded-component
polyhedral, and homothetic radial classes, and identifies an exact
gap-sensitive potential certificate.  It does not claim the arbitrary
complex theorem.
