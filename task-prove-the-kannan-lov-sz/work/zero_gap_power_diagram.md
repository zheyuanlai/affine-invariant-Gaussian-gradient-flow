# The exact zero-gap branch: translations, walls, and rank

## 0. Verdict

Let \(K\subset\mathbb R^d\) be a convex body, let \(\mu\) be its
isotropic uniform probability, let \(|E|=|K|/2\), and let \(T\) be the
Brenier map from \(2\,1_E\mu\) to \(2\,1_{K\backslash E}\mu\).  Suppose

\[
 D_aT=I\qquad\text{Lebesgue-a.e. on }E.                  \tag{0.1}
\]

This note records the strongest verified consequences of (0.1).

First, a necessary correction: (0.1) does **not** imply that
\(T=I+\nabla h\) with \(h\) convex and purely singular Hessian.  The
Brenier potential \(\phi\) is convex, but
\(h=\phi-|x|^2/2\) is only \(1\)-semiconvex in general and can have
negative absolutely continuous Hessian on target gaps.  Section 1 gives
an exact one-dimensional counterexample to the implication.

There is a useful stronger subcase:

\[
 T(x)=x+\nabla h(x)\quad(x\in E),\qquad h\text{ convex},  \tag{0.2}
\]

with \(h\) polyhedral at first.  In this subcase:

* \(T\) is expansive on \(E\);
* the displacement label on the midpoint cells is firmly
  nonexpansive;
* if the translations span an \(r\)-dimensional subspace, then

  \[
  \boxed{\int_E |T(x)-x|^2\,{2\,dx\over|K|}\le4r;}       \tag{0.3}
  \]

* every source--source power wall opens a disjoint midpoint prism, and
  the exact total wall budget is at most \(|K|/2\);
* for arbitrary rank, expansiveness and the thin-shell theorem give only
  the dimension-dependent fallback

  \[
  W_2^2(2\,1_E\mu,2\,1_{K\backslash E}\mu)\le C\sqrt d.  \tag{0.4}
  \]

The one-cell and collinear cases are therefore dimension-free.  The
unresolved finite-cell case is a high-rank collection of low-strain
translations separated either by disjoint wall prisms or by physical
target gaps.  No growing-dimensional example satisfying convex support,
complementarity, isotropy, and Brenier optimality was found.  Conversely,
the wall budget by itself is an \(L^1\) control of the jump Hessian and
does not yet control the \(L^2\) translation amplitudes without a global
connectivity inequality.  Thus this note does not claim the desired
dimension-free zero-gap lemma.

## 1. Zero strain does not imply a convex displacement potential

Take \(K=[0,4]\) with uniform measure and

\[
 E=[0,1]\cup[3,4],\qquad K\backslash E=[1,3].            \tag{1.1}
\]

The increasing rearrangement is

\[
 T(x)=
 \begin{cases}
 x+1,&0\le x\le1,\\
 x-1,&3\le x\le4.
 \end{cases}                                             \tag{1.2}
\]

It pushes the uniform restriction of \(E\) to the uniform restriction of
its complement and is the unique one-dimensional Brenier map.  Moreover
\(T'=1\) source-almost everywhere, so both the potential and determinant
midpoint deficits vanish.

The monotone extension of \(T=\phi'\) across the target gap is
\(T(x)=2\) on \([1,3]\).  Consequently

\[
 h=\phi-{x^2\over2}
 \quad\Longrightarrow\quad
 h''=
 \begin{cases}
 0,&x\in E,\\
 -1,&x\in(1,3).
 \end{cases}                                             \tag{1.3}
\]

The source translations decrease from \(+1\) to \(-1\), so \(h\) is not
convex.  After centering and scaling \(K\) to variance one, the model
remains isotropic and has \(W_2^2=3/4\).

More generally, for **every** Borel half-partition of a uniform interval,
the increasing rearrangement has derivative one at almost every source
density point.  All geometric information is then carried by jumps,
singular continuous increments, and target gaps.  Any proof of the whole
zero-gap branch must retain those parts.

Thus (0.2) is an additional expansiveness hypothesis, not a consequence
of (0.1).

## 2. Finite convex-displacement model

Assume now the stronger hypothesis (0.2) with

\[
 h(x)=\max_{1\le i\le N}\{a_i\cdot x+b_i\}.              \tag{2.1}
\]

Let

\[
 C_i=\{x:a_i\cdot x+b_i\ge a_j\cdot x+b_j
             \text{ for all }j\},\qquad P_i=E\cap C_i,   \tag{2.2}
\]

up to null tie sets.  Then

\[
 T(x)=x+a_i\quad(x\in P_i),\qquad Q_i=P_i+a_i,           \tag{2.3}
\]

and complementarity says

\[
 K=\mathop{\dot\bigcup}_{i=1}^N P_i
   \mathop{\dot\bigcup}_{i=1}^N Q_i
 \quad\text{modulo null sets}.                           \tag{2.4}
\]

For \(x\in P_i\), \(x'\in P_j\), the power inequalities give

\[
 (a_i-a_j)\cdot(x-x')\ge0.                               \tag{2.5}
\]

Consequently

\[
 |T(x)-T(x')|^2
 =|x-x'|^2+2(a_i-a_j)\cdot(x-x')+|a_i-a_j|^2
 \ge|x-x'|^2.                                            \tag{2.6}
\]

Thus \(T\) is expansive.  Its inverse on the target is a
volume-preserving \(1\)-Lipschitz map, piecewise given by
\(y\mapsto y-a_i\).

Put

\[
 R_i=P_i+{a_i\over2},\qquad R=\mathop{\dot\bigcup}_iR_i. \tag{2.7}
\]

The midpoint law is uniform on \(R\), and \(|R|=|K|/2\).  If
\(z=x+a_i/2\in R_i\) and \(z'=x'+a_j/2\in R_j\), (2.5) becomes

\[
 (a_i-a_j)\cdot(z-z')\ge{1\over2}|a_i-a_j|^2.            \tag{2.8}
\]

Equivalently, the label

\[
 b(z)={a_i\over2}\quad(z\in R_i)                         \tag{2.9}
\]

is firmly nonexpansive:

\[
 (b(z)-b(z'))\cdot(z-z')\ge|b(z)-b(z')|^2.              \tag{2.10}
\]

In particular \(b\) is \(1\)-Lipschitz on \(R\).  It has zero
Lebesgue derivative on each midpoint cell; all variation occurs across
gaps and singular walls.

The same statement holds without polyhedrality.  If \(h\) is any proper
closed convex function and \(z=x+\nabla h(x)/2\), then
\(b(z)=\nabla h(x)/2\) is the gradient of the Moreau envelope of
\(h/2\); it extends to a firmly nonexpansive map on the ambient space.
This extension fact alone gives only
\(\mathbb E|b-\mathbb Eb|^2\le O(d)\), not a universal trace bound.

## 3. Covariance and expansiveness identities

Let \(X\) be uniform on \(E\), \(Y=T(X)\), \(A=Y-X\), and

\[
 m=\mathbb EX=-\mathbb EY,\qquad
 A_c=A-\mathbb EA=A+2m.                                 \tag{3.1}
\]

Write

\[
 S=\mathbb E|A|^2,\qquad S_c=\mathbb E|A_c|^2
   =S-4|m|^2.                                            \tag{3.2}
\]

Isotropy of the full mixture gives

\[
 \mathbb E[A\otimes A]\preceq4I,\qquad |m|\le1.          \tag{3.3}
\]

Taking two independent source points in (2.5) yields

\[
 \operatorname {tr}\operatorname {Cov}(X,A)\ge0.         \tag{3.4}
\]

Since \(Y=X+A\),

\[
\begin{split}
 \operatorname {tr}\operatorname {Cov}(Y)
 -\operatorname {tr}\operatorname {Cov}(X)
 &=2\operatorname {tr}\operatorname {Cov}(X,A)+S_c\\
 &\ge S_c.
\end{split}                                               \tag{3.5}
\]

Thus the target half is at least as spread as the source half, and the
centered translation energy is bounded by their trace-variance gap.
This is the precise clustered-to-spread orientation imposed by convexity
of \(h\).

The displacement effective-rank bound is

\[
 {\operatorname {tr}\operatorname {Cov}(A)
  \over\|\operatorname {Cov}(A)\|_{\mathrm{op}}}
 \ge {S-4\over4}                                         \tag{3.6}
\]

whenever \(S>4\).  Hence a hypothetical large-cost example must involve
many displacement directions; one long translation is impossible.

## 4. Exact dimension-free bound for low translation rank

The following lemma does not require the convexity of \(h\), only the
fiber-preserving property.

**Lemma 4.1 (translation-span bound).**  Suppose all displacements
\(T(x)-x\) belong to a fixed \(r\)-dimensional linear subspace
\(U\subset\mathbb R^d\).  Then

\[
 W_2^2(2\,1_E\mu,2\,1_{K\backslash E}\mu)\le4r.          \tag{4.1}
\]

**Proof.**  Let \(y=P_{U^\perp}x\) and
\(K_y=K\cap(y+U)\).  Since \(T\) preserves \(y\), disintegration of the
pushforward identity gives

\[
 |E\cap K_y|=|(K\backslash E)\cap K_y|={1\over2}|K_y|    \tag{4.2}
\]

for almost every nonempty fiber.  The restriction of \(T\) to that
fiber is an optimal quadratic coupling; otherwise replacing it on a
positive set of fibers would lower the global cost.

Let \(c_y\) and \(\Sigma_y\) be the centroid and covariance of the
uniform probability on \(K_y\).  Couple the two conditional half-laws
independently.  If their means relative to \(c_y\) are \(q_y\) and
\(-q_y\), its cost is

\[
 2\operatorname {tr}\Sigma_y+2|q_y|^2
 \le4\operatorname {tr}\Sigma_y.                         \tag{4.3}
\]

Optimality of the conditional Brenier coupling gives the same upper
bound for its cost.  The source mixing weight of the fiber is
\(|K_y|/|K|\).  Therefore

\[
\begin{split}
 W_2^2
 &\le {4\over|K|}\int_{U^\perp}
       |K_y|\operatorname {tr}\Sigma_y\,dy\\
 &\le4\operatorname {tr}(P_U\operatorname {Cov}(\mu)P_U)
 =4r.
\end{split}                                               \tag{4.4}
\]

The second inequality is the conditional-variance decomposition. \(\square\)

For a single translation, complementarity gives
\(a=\mathbb EY-\mathbb EX=-2m\), hence the sharper estimate
\(|a|^2\le4\).  For collinear translations Lemma 4.1 gives \(S\le4\).
Thus every fixed-rank finite power diagram is dimension-free.  The lemma
does not control a high-rank family of individually tiny translations.

## 5. The exact wall-prism budget

For two adjacent power cells, let

\[
 W_{ij}=\{x:a_i\cdot x+b_i=a_j\cdot x+b_j
             \ge a_k\cdot x+b_k\ \forall k\}.            \tag{5.1}
\]

Its unit normal is

\[
 n_{ij}={a_i-a_j\over|a_i-a_j|}.                         \tag{5.2}
\]

On a relatively open source--source facet
\(\Gamma_{ij}\subset W_{ij}\) for which both adjacent traces belong to
\(E\), the distributional Hessian is

\[
 D^2h\!\restriction_{\Gamma_{ij}}
 =|a_i-a_j|\,n_{ij}\otimes n_{ij}\,
   \mathcal H^{d-1}\!\restriction_{\Gamma_{ij}}.          \tag{5.3}
\]

At \(x\in\Gamma_{ij}\), the midpoint subgradient segment is

\[
 \Pi_{ij}(x)=x+{1\over2}[a_i,a_j].                       \tag{5.4}
\]

Both endpoints are limits of regular midpoint points and lie in \(K\);
convexity of \(K\) puts the whole segment in \(K\).  Its relative normal
length is \(|a_i-a_j|/2\), so

\[
 \left|\bigcup_{x\in\Gamma_{ij}}\Pi_{ij}(x)\right|
 ={1\over2}|a_i-a_j|\,\mathcal H^{d-1}(\Gamma_{ij}).      \tag{5.5}
\]

These open wall prisms are pairwise disjoint, including across different
walls.  Indeed they are subgradient images under the strongly convex
function

\[
 g(x)={|x|^2\over2}+{h(x)\over2};                        \tag{5.6}
\]

strong monotonicity makes the inverse subgradient single-valued.
They are also disjoint from the regular midpoint image \(R\).  Since
\(|R|=|K|/2\), summing (5.5) yields

\[
 \boxed{\sum_{i<j}|a_i-a_j|\,
       \mathcal H^{d-1}(\Gamma_{ij})\le|K|.}              \tag{5.7}
\]

At a codimension-\(k\) junction, the corresponding \(k\)-dimensional
polytope of active slopes combines with the junction face to give an
additional nonnegative \(d\)-volume.  Strong monotonicity again makes all
such strata disjoint.  Their total, together with the wall-prism volumes,
is at most \(|K|/2\).

Equation (5.7) is an exact total-variation budget for the singular Hessian
on source adjacency walls.  It does not see two source components
separated by a positive-volume target corridor.  Nor does an upper bound
on jump total variation imply an \(L^2\) bound on the cell labels without
a Poincare/connectivity inequality for the cell complex.  The one-
dimensional examples show both defects sharply.

## 6. Convex scalar witness and the thin-shell fallback

For \(x\in P_i\), convexity of \(h\) gives

\[
 h(Tx)-h(x)\ge a_i\cdot(Tx-x)=|a_i|^2.                  \tag{6.1}
\]

After integration,

\[
 \boxed{
 \mathbb E_{\mu_-}h-\mathbb E_{\mu_+}h\ge S.}            \tag{6.2}
\]

Thus a large-cost convex-displacement example automatically supplies a
convex scalar function whose two conditional means are separated by at
least the squared transport cost.  A dimension-free convex Poincare or
first-moment bound strong enough to compare the left side with
\((\mathbb E_E|\nabla h|^2)^{1/2}=\sqrt S\) would close the branch.
Such an inequality is not presently derived from isotropy here and may
not be assumed.

There is a rigorous dimension-dependent fallback.  From (3.5),

\[
 S_c\le
 \mathbb E_{\mu_-}|X|^2-\mathbb E_{\mu_+}|X|^2
 \le2\,\mathbb E_\mu\big||X|^2-d\big|.                  \tag{6.3}
\]

The thin-shell theorem gives
\(\mathbb E(|X|-\sqrt d)^2\le C\) for isotropic log-concave laws.
Cauchy--Schwarz and
\(\mathbb E(|X|+\sqrt d)^2\le4d\) imply

\[
 \mathbb E\big||X|^2-d\big|\le C\sqrt d.                 \tag{6.4}
\]

Since \(S=S_c+4|m|^2\) and \(|m|\le1\), (0.4) follows.  This estimate is
not returnable for KLS; it only verifies that the expansive branch cannot
have cost larger than order \(\sqrt d\).

## 7. Approximation audit

Suppose \(h_k\) are polyhedral convex functions converging locally
uniformly to a convex \(h\).  Their resolvents and Moreau-envelope
gradients converge locally uniformly, so the firm midpoint inequality
(2.10) passes to the limit.  The Hessian measures \(D^2h_k\) converge
weakly-* to \(D^2h\), and the full subgradient-image volume in Section 5
is lower semicontinuous.

This does **not** by itself approximate a complementary half-transport.
One must simultaneously choose source sets \(E_k\) such that

\[
 (I+\nabla h_k)_\#(1_{E_k}dx)=1_{K\backslash E_k}dx,     \tag{7.1}
\]

with no loss of mass through cell walls or the boundary of \(K\).
Local uniform convergence of \(h_k\) does not preserve (7.1).  Therefore
finite-cell estimates may pass to a general singular \(h\) only after a
separate complementarity-stable approximation theorem.  No such theorem
is assumed here.

## 8. Remaining finite-cell question

The verified facts reduce the stronger convex-displacement problem to the
following concrete high-rank statement.

> Let a convex body \(K\) be partitioned into translated pairs
> \(P_i,Q_i=P_i+a_i\), with the source pieces contained in the power cells
> of the convex polyhedral function
> \(\max_i(a_i\cdot x+b_i)\).  Assume \(K\) is isotropic and the two
> unions have equal volume.  Prove
> \[
> {2\over|K|}\sum_i|P_i|\,|a_i|^2\le C.
> \]

Low-dimensional translation span is handled by Lemma 4.1.  Large cost
would force effective displacement rank at least \((S-4)/4\), while all
source wall prisms and higher junction prisms have total volume at most
\(|K|/2\).  The missing implication is a dimension-free conversion of
this wall/gap volume into \(L^2\) control of a high-rank label field.

The premise \(D_aT=I\) alone is broader still because of Section 1.  Any
eventual proof must first handle negative displacement curvature on
target gaps, then pass from finite translations to singular continuous
Brenier increments.  No growing counterexample satisfying all of the
listed conditions is certified by the present analysis.
