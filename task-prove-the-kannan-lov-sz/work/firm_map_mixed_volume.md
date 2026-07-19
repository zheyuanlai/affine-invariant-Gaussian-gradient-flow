# Firm midpoint maps, martingale dilation, and the mixed-volume gap

## 0. Verdict

This note audits the proposed dimension-free mixed-volume/Hessian theorem
for the exact zero-entropy-gap branch of balanced transport in a convex
body.

Let \(K\subset\mathbb R^d\) be a convex body, let \(R\subset K\) have
volume \(|K|/2\), and let \(b:R\to\mathbb R^d\) be locally constant on
measurable cells.  Assume that

\[
 K=(I-b)(R)\mathbin{\dot\cup}(I+b)(R)                 \tag{0.1}
\]

modulo null sets, with both branch maps injective and volume preserving.
These are the midpoint data of a zero-gap piecewise-translation coupling.
In the convex-displacement Brenier subbranch, \(b\) is the restriction of
the gradient of a convex \(1\)-smooth function and is firmly nonexpansive:

\[
 (b(z)-b(z'))\cdot(z-z')\ge |b(z)-b(z')|^2.           \tag{0.2}
\]

The audit produces one new exact global reduction and one decisive
counterexample to the most natural local-to-global proof.

1.  Equation (0.1) makes the uniform probability on \(K\) a **binary
    martingale dilation** of the uniform probability on \(R\): if
    \(Z\) is uniform on \(R\), \(B=b(Z)\), and \(\varepsilon\) is a fair
    sign, then

    \[
    U=Z+\varepsilon B
    \]

    is uniform on \(K\) and \(\mathbb E[U\mid Z]=Z\).

2.  Iterating this binary step until it first exits \(R\) gives a
    martingale from the uniform law on \(R\) to the uniform law on
    \(K\setminus R\).  Its lifetime has the exact geometric tail
    \(\mathbb P(\tau>n)=2^{-n}\), and

    \[
    \boxed{\mathbb E|Z_\tau-Z_0|^2=2\mathbb E_R|b|^2.} \tag{0.3}
    \]

    Equivalently, if \(K\) is isotropic,

    \[
    \operatorname {Cov}(R)=I-\mathbb E(BB^T),\qquad
    \operatorname {Cov}(K\setminus R)=I+\mathbb E(BB^T). \tag{0.4}
    \]

    This retains physical target gaps automatically and is stronger than
    the local wall-prism budget.  It does not by itself bound the quadratic
    variation: doing so is exactly the new global step.

3.  Local singular-Hessian mass cannot supply that step.  The realizable
    one-dimensional convex-body example

    \[
    K=[0,6],\qquad E=[0,1]\cup[2,4]                    \tag{0.5}
    \]

    has translations \(+1\) and \(+2\).  It satisfies the full cyclic
    Brenier, firm-midpoint, complementarity, and zero-gap conditions, but
    its source--source wall functional is zero.  Nevertheless, after
    isotropic normalization,

    \[
    \mathbb E|B|^2=\frac14,qquad
    \operatorname {Var}(R)=\frac34,qquad
    -\log\det\operatorname {Cov}(R)=\log\frac43>0.     \tag{0.6}
    \]

    Thus entropy deficit, Cayley strain, and the exact local source-wall
    prism volume can all vanish while the covariance determinant drops.
    Any valid functional must also charge physical gaps between source
    components.

4.  Extending \(b\) firmly across those gaps does not finish the mixed-
    volume argument.  In the same example the canonical cyclic firm
    extension sends part of \(K\setminus R\) outside \(K\) under
    \(I+b\).  Therefore the tempting bound

    \[
    |(I+b)(K\setminus R)|\le |K|-|(I+b)(R)|            \tag{0.7}
    \]

    is false.  Without (0.7), the determinant expansion of
    \(I+D b\) has no ambient-volume budget.

No proof of \(\mathbb E|B|^2\le C\) for every actual zero-gap Brenier
complex is obtained, and no growing-dimensional counterexample to that
absolute statement is constructed.  The new martingale formulation shows
precisely what remains: prove a dimension-free quadratic-variation bound
for this special cyclic, firmly nonexpansive, binary martingale
decomposition of an isotropic convex body.  Pairwise firm
nonexpansiveness alone is not the full Brenier condition; cyclic
monotonicity must be retained.

Finally, the **general** absolute midpoint determinant bound, over all
balanced cuts and not merely the exact zero-gap subclass, is quantitatively
equivalent to the balanced-half \(W_2\) bound and hence to KLS.  The
zero-gap-only statement is a distinguished KLS branch; no reduction from
the general positive-strain case to that branch is currently proved.

## 1. Actual midpoint data and the hierarchy of hypotheses

In the finite convex-displacement model, start with source cells \(P_i\),
translations \(a_i\), target cells \(Q_i=P_i+a_i\), and midpoint cells

\[
 R_i=P_i+\frac{a_i}{2},\qquad b_i=\frac{a_i}{2}.       \tag{1.1}
\]

Then

\[
 P_i=R_i-b_i,qquad Q_i=R_i+b_i,                       \tag{1.2}
\]

and complementarity is exactly (0.1).  If the displacement is the
gradient of a convex polyhedral function, the midpoint label is the
gradient of its Moreau envelope and extends to a cyclically firmly
nonexpansive map on \(\mathbb R^d\).  In particular, (0.2) holds.

It is important to distinguish three levels.

* Pairwise firm nonexpansiveness is the metric inequality (0.2).
* Cyclic firm nonexpansiveness means that \(b=\nabla\psi\) for a convex
  function \(\psi\) with \(0\preceq D^2\psi\preceq I\), in the smooth
  reading.
* Actual zero-gap Brenier realizability additionally requires the
  volume-preserving complementary partition (0.1), with the correct
  cellwise translations.

The first condition does not imply the second.  For example, let

\[
 J=\begin{pmatrix}0&-1\\1&0\end{pmatrix},
 \qquad M=(I+J)^{-1}=\frac12(I-J).                     \tag{1.3}
\]

Then

\[
 M^TM=\frac12I=\frac{M+M^T}{2},                       \tag{1.4}
\]

so \(z\mapsto Mz\) is firmly nonexpansive, with equality in (0.2).
But \(M\) is not symmetric and hence is not the gradient of a convex
function.  The failure is already finite.  On the five clockwise vertices
\(z_k=(\cos(2\pi k/5),-\sin(2\pi k/5))\), the cyclic sum is

\[
 \sum_{k=0}^4 (Mz_k)\cdot(z_k-z_{k+1})
 =\frac52\left(1-\cos\frac{2\pi}{5}
                  -\sin\frac{2\pi}{5}\right)<0,       \tag{1.5}
\]

with indices modulo five.  Thus these five points give pairwise firm data
which are not cyclic Brenier data.  A counterexample built only from
(0.2) would not decide the actual zero-gap branch.

## 2. Binary martingale dilation

Let \(\lambda_R\) and \(\lambda_K\) denote the uniform probabilities on
\(R\) and \(K\).  The branch maps

\[
 F_-(z)=z-b(z),\qquad F_+(z)=z+b(z)                    \tag{2.1}
\]

are translations on each plateau cell, hence have Jacobian one there.
Assumption (0.1) says their images are disjoint and cover \(K\).

The argument below is not intrinsically polyhedral.  In an arbitrary
uniform zero-gap transport, the midpoint density equals \(2/|K|\) on its
regular image and zero elsewhere, so it is uniform on a measurable
half-volume set \(R\).  Conditional on a midpoint, its two deterministic
endpoints have midpoint \(z\) and can be written \(z\pm b(z)\).  Their
equal mixture is uniform on \(K\).  These measure-theoretic facts are all
that Lemmas 2.1 and 3.1 use.

**Lemma 2.1 (binary dilation).**  If \(Z\sim\lambda_R\) and
\(\varepsilon\) is an independent fair sign, then

\[
 U=Z+\varepsilon b(Z)                                  \tag{2.2}
\]

has law \(\lambda_K\), and \((Z,U)\) is a martingale coupling:

\[
 \mathbb E[U\mid Z]=Z.                                 \tag{2.3}
\]

**Proof.**  For a Borel set \(A\subset K\), cellwise change of variables
and (0.1) give

\[
\begin{aligned}
 \mathbb P(U\in A)
 &=\frac12\frac{2}{|K|}
   \left(|F_-^{-1}(A)|+|F_+^{-1}(A)|\right)\\
 &=\frac{|A|}{|K|}.
\end{aligned}                                           \tag{2.4}
\]

The conditional mean identity follows immediately from the symmetric
sign. \(\square\)

Suppose now that \(K\) is centered and isotropic.  Write
\(B=b(Z)\).  Taking the first two moments in (2.2) gives

\[
 \mathbb EZ=0,qquad
 I=\mathbb E[ZZ^T]+\mathbb E[BB^T].                    \tag{2.5}
\]

Thus

\[
 \mathbb E[BB^T]\preceq I,qquad |\mathbb EB|\le1.   \tag{2.6}
\]

The desired theorem is the trace upgrade

\[
 \operatorname {tr}\mathbb E[BB^T]\le C.             \tag{2.7}
\]

Equation (2.5) shows that it is exactly a dimension-free covariance loss
for a half-volume midpoint set which admits the special binary dilation
(2.2).

The martingale property also gives convex order:

\[
 \lambda_R\preceq_{cx}\lambda_K.                      \tag{2.8}
\]

This is stronger than mere inclusion \(R\subset K\), but convex order
alone does not visibly bound the trace in (2.7); the firm/cyclic geometry
must enter.

## 3. The killed midpoint chain

The dilation can be iterated because \(R\subset K\).  Starting with
\(Z_0\sim\lambda_R\), choose independent fair signs \(\varepsilon_n\)
and, while \(Z_n\in R\), set

\[
 Z_{n+1}=Z_n+\varepsilon_{n+1}b(Z_n).                  \tag{3.1}
\]

Let

\[
 \tau=\inf\{n\ge1:Z_n\notin R\}.                     \tag{3.2}
\]

After time \(\tau\), keep the chain fixed.  Put \(G=K\setminus R\).

**Lemma 3.1 (geometric survival and exact quadratic variation).**  The
chain above satisfies

\[
 \mathbb P(\tau>n)=2^{-n},                             \tag{3.3}
\]

and, conditional on \(\tau>n\), \(Z_n\) is uniform on \(R\).  The
terminal point \(Z_\tau\) is uniform on \(G\), and

\[
 \boxed{
 \mathbb E|Z_\tau-Z_0|^2
 =\mathbb E\sum_{j=0}^{\tau-1}|b(Z_j)|^2
 =2\mathbb E_R|b|^2.}                                  \tag{3.4}
\]

**Proof.**  If the current state is uniform on \(R\), Lemma 2.1 says the
next, unkilled state is uniform on \(K\).  It lands in \(R\) with
probability \(1/2\), and conditional on that event it is again uniform on
\(R\).  Induction proves (3.3) and the conditional-law assertion.

For a Borel set \(A\subset G\), the probability of surviving \(n-1\)
steps and then exiting into \(A\) is

\[
 2^{-(n-1)}\frac{|A|}{|K|}.
\]

Summing over \(n\ge1\) gives \(2|A|/|K|\), the uniform probability on
\(G\).

The stopped process is a bounded martingale.  Its increments are
orthogonal in \(L^2\), so

\[
 \mathbb E|Z_\tau-Z_0|^2
 =\sum_{j\ge0}\mathbb E\left[
   \mathbf1_{\{\tau>j\}}|b(Z_j)|^2\right].            \tag{3.5}
\]

Conditional on \(\tau>j\), the state is uniform on \(R\), and the event
has probability \(2^{-j}\).  The geometric sum is two, proving (3.4).
\(\square\)

Since \(R\) and \(G\) have equal volume and the full body is centered,
both conditional means are zero.  The total covariance identity and
(2.5) then give

\[
 \operatorname {Cov}(\lambda_R)=I-M,qquad
 \operatorname {Cov}(\lambda_G)=I+M,qquad
 M=\mathbb E[BB^T].                                    \tag{3.6}
\]

Taking traces in (3.6) recovers (3.4): the terminal martingale cost is
the covariance increase \(2\operatorname {tr}M\).

This formulation retains every physical gap, including gaps on which the
displacement potential has negative curvature in the broader nonconvex-
displacement zero-gap branch.  What it does not supply is a universal
bound on martingale quadratic variation.

## 4. What firm and cyclic geometry say along the chain

Assume (0.2).  If one step survives, write

\[
 z'=z+\varepsilon b(z)\in R,qquad
 d=b(z')-b(z).
\]

Firm nonexpansiveness gives

\[
 d\cdot(\varepsilon b(z))\ge|d|^2.                    \tag{4.1}
\]

For a forward, \(\varepsilon=+1\), transition,

\[
 |b(z')|^2\ge |b(z)|^2+3|d|^2,                        \tag{4.2}
\]

while for a backward transition,

\[
 |b(z')|^2\le |b(z)|^2-|d|^2.                         \tag{4.3}
\]

Thus labels cannot wander arbitrarily along a surviving path: the sign of
the branch imposes a strict Lyapunov change unless the label stays fixed.
This recovers the one-dimensional corridor ordering and excludes naive
parity tensorizations, whose inverse branches reverse label order.

Equations (4.2)--(4.3) do not yet control \(|b(z)|\) itself.  A path has
mean lifetime two, but it can start with a high-rank label of large norm;
the transition inequalities only charge changes of that label.  This is
the pathwise version of the distinction between wall jumps and cellwise
translation constants.

Cyclic monotonicity supplies all higher cycle inequalities in addition to
(4.1).  Those inequalities are indispensable: the linear map (1.3) obeys
every pairwise firm inequality but is not a Brenier gradient.  No argument
below replaces cyclic data by pairwise data.

## 5. Why local wall prisms do not control the determinant

In a polyhedral cyclic model, a source--source wall \(\Gamma_{ij}\)
between labels \(b_i,b_j\) opens the disjoint midpoint prism

\[
 \{x+[b_i,b_j]:x\in\Gamma_{ij}\}.
\]

Consequently the exact normalized source-wall functional is

\[
 \mathfrak W_{\rm src}
 =\frac1{|K|}\sum_{i<j}|b_i-b_j|
       \mathcal H^{d-1}(\Gamma_{ij})
 \le\frac12.                                           \tag{5.1}
\]

Higher junction strata give additional nonnegative mixed volumes, also
with total at most \(1/2\).  This is a universal singular-Hessian budget,
but only on walls where both traces are source cells.

The following example proves that no coercive determinant estimate can
use only (5.1).

**Proposition 5.1 (realizable physical-gap obstruction).**  Let

\[
 K=[0,6],\qquad E=[0,1]\cup[2,4],qquad
 K\setminus E=[1,2]\cup[4,6].                         \tag{5.2}
\]

The increasing Brenier map is

\[
 T(x)=
 \begin{cases}
  x+1,&x\in[0,1],\\
  x+2,&x\in[2,4].
 \end{cases}                                           \tag{5.3}
\]

Its midpoint data are

\[
\begin{array}{c|c|c|c}
R_i&b_i&R_i-b_i&R_i+b_i\\ \hline
[1/2,3/2]&1/2&[0,1]&[1,2]\\
[3,5]&1&[2,4]&[4,6].
\end{array}                                             \tag{5.4}
\]

They satisfy cyclic firm nonexpansiveness and (0.1), but
\(\mathfrak W_{\rm src}=0\).  Moreover the entropy deficit and Cayley
strain vanish for every displacement time.  Nevertheless, after centering
at three and scaling by \(1/\sqrt3\) to isotropy,

\[
 \mathbb E_R B^2=\frac14,
 \qquad \operatorname {Var}(R)=\frac34,
 \qquad -\log\det\operatorname {Cov}(R)=\log\frac43.
\tag{5.5}
\]

**Proof.**  The map (5.3) is increasing and pushes Lebesgue measure on
\(E\) to Lebesgue measure on its complement, hence is the Brenier map.
It has derivative one at every source density point.  The two midpoint
labels differ by \(1/2\), while the distance between their midpoint cells
is at least \(3-3/2=3/2\); hence (0.2) holds.  The power wall separating
the two slopes lies in the physical target gap \([1,2]\), not on a
source--source interface, so (5.1) is zero.

Before isotropic scaling, the label weights are \(1/3,2/3\), whence

\[
 \mathbb EB^2=\frac13\frac14+\frac23=\frac34.          \tag{5.6}
\]

The midpoint set has mean three and variance \(9/4\); the ambient interval
has variance three.  Division by three under isotropic scaling proves
(5.5). \(\square\)

Thus an inequality of the form

\[
 -\log\det\operatorname {Cov}(R)
 \le C\left(\text{entropy deficit}+\text{Cayley strain}
             +\mathfrak W_{\rm src}\right)            \tag{5.7}
\]

is false, even in one dimension and even for actual cyclic Brenier data.
A universal additive baseline would make Proposition 5.1 harmless, but
then the unproved content is exactly the desired absolute bound.

## 6. The failed global Hessian-volume completion

Let \(b=\nabla\psi\) be the canonical cyclic firm extension, so in the
smooth reading

\[
 0\preceq D^2\psi\preceq I.                            \tag{6.1}
\]

On a physical gap \(G=K\setminus R\), the natural full Hessian charge is

\[
 \mathfrak H_1
 =\frac1{|K|}\int_G\operatorname {tr}D^2\psi.         \tag{6.2}
\]

The associated mixed-volume charge is

\[
 \mathfrak H_{\det}
 =\frac1{|K|}\int_G
   \left[\det(I+D^2\psi)-1\right].                    \tag{6.3}
\]

Because all elementary symmetric polynomials of a positive semidefinite
matrix are nonnegative,

\[
 \mathfrak H_{\det}\ge\mathfrak H_1.                 \tag{6.4}
\]

If one had the global inclusion

\[
 (I+b)(G)\subset K\setminus(I+b)(R),                  \tag{6.5}
\]

the area formula would bound \(\mathfrak H_{\det}\) by the unused half
of \(K\).  Since \(|G|=|K|/2\), it would actually force
\(\mathfrak H_{\det}=0\); a fixed-dilation variant would still give a
dimension-free trace-Hessian budget.  This is the most direct global
mixed-volume completion of the local wall-prism proof.

Inclusion (6.5) is false for actual data.  In Proposition 5.1, the convex
displacement can be chosen with its power wall at \(x=3/2\).  Its canonical
midpoint firm extension is

\[
 b(z)=
 \begin{cases}
  1/2,&z\le2,\\
  z-3/2,&2\le z\le5/2,\\
  1,&z\ge5/2.
 \end{cases}                                           \tag{6.6}
\]

It agrees with the plateau labels on \(R\).  But on the physical gap
\((5,6]\subset K\setminus R\),

\[
 (I+b)(z)=z+1\in(6,7],                                 \tag{6.7}
\]

which lies outside \(K\).  Similarly, the opposite branch exits at the
left edge.  Thus the Jacobian expansion on all physical gaps cannot be
charged to unused volume inside the original convex body.

One can restore a valid mixed-volume estimate by imposing a common convex
core or a sequential family of global extrusions.  Under those additional
hypotheses the projection-thickness argument gives

\[
 \sum_j|v_j|\le\sqrt{24},\qquad
 \mathbb E|B|^2\le13.                                  \tag{6.8}
\]

The issue is not positivity of mixed volumes; (6.4) has the correct sign.
The missing statement is a range/localization theorem which replaces the
false global inclusion (6.5) by a laminar collection of large convex cores
while retaining every physical gap.

## 7. The exact candidate theorem after the audit

The martingale formulation isolates a clean global statement.

> **Cyclic binary-dilation theorem.**  Let \(K\) be an isotropic convex
> body.  Let \(R\subset K\) have half its volume, and suppose that a
> locally constant cyclically firmly nonexpansive map \(b:R\to\mathbb R^d\)
> satisfies the complementary, injective, volume-preserving identity
> \(K=(I-b)(R)\dot\cup(I+b)(R)\).  Then
> \[
> \mathbb E_R|b|^2\le C.
> \tag{7.1}
> \]

By Lemma 3.1, this is equivalently a dimension-free quadratic-variation
bound for the stopped martingale whose stationary survival law is uniform
on \(R\) and whose terminal law is uniform on \(K\setminus R\).

All currently verified cases fit one of the following mechanisms.

* If the labels span an \(r\)-dimensional space, fiber disintegration gives
  \(\mathbb E|B|^2\le r\).
* Orthogonal product corridors lose volume multiplicatively and have
  universal energy.
* Sequential noncommuting corridors sharing a half-volume core satisfy
  (6.8).
* Symmetric simplex and crosspolytope fans have universal energy by exact
  overlap or homothetic-core calculations.
* In one dimension the physical gaps order all labels, and isotropy gives
  the required bound.

A hypothetical counterexample must therefore have all of the following:
growing centered label rank, no common large convex core, small local wall
patches, substantial physical-gap variation, and genuinely cyclic (not
merely pairwise firm) label data.  No realizable complex with these
properties is constructed here.

The exact unproved local-to-global estimate can be stated as a replacement
for (6.5): decompose the canonical Hessian measure on physical gaps into
random laminar extrusion chains \((C_j,v_j)\), each starting from a convex
core of volume at least \(c|K|\), so that

\[
 \mathbb E_R|B-B'|^2
 \le C\,\mathbb E_{\rm chains}\left(\sum_j|v_j|\right)^2. \tag{7.2}
\]

Projection thickness and mixed-volume positivity would then prove (7.1).
Neither the local wall budget (5.1) nor the martingale identity (3.4)
constructs these chains.

## 8. Relation to the midpoint determinant and KLS

For an actual balanced-half displacement, the midpoint covariance is

\[
 C_{1/2}=I-\frac14K=I-\mathbb E(BB^T).                 \tag{8.1}
\]

The domination \(\nu_{1/2}\le2\mu\), together with one-dimensional
log-concavity of projections, gives a universal spectral floor

\[
 c_0I\preceq C_{1/2}\preceq I                         \tag{8.2}
\]

with, for example, \(c_0=1/4800\).  Hence

\[
 \mathbb E|B|^2
 \le-\log\det C_{1/2}
 \le c_0^{-1}\mathbb E|B|^2.                          \tag{8.3}
\]

Thus the cyclic binary-dilation theorem is exactly the absolute determinant
bound on the zero-gap subclass.

For **all** balanced cuts, a universal bound on
\(-\log\det C_{1/2}\) is quantitatively equivalent by (8.3) to a universal
bound on

\[
 W_2^2(2\mathbf1_E\mu,2\mathbf1_{E^c}\mu).
\]

Poincare implies this balanced-half transport bound, and the reverse bound
controls the median \(L^1\) deviation of every 1-Lipschitz function; for
log-concave measures that scale is equivalent to the Cheeger/Poincare
scale.  Therefore the general absolute determinant bound is exactly
KLS-equivalent up to numerical constants.

The qualifier “general” matters.  Proving (7.1) only for exact zero-gap
piecewise translations would close the singular-translation branch exposed
by entropy, but it would not by itself prove KLS.  One would still need a
stability theorem converting the bounded Cayley/potential strain of a
general transport into a controlled perturbation of cyclic binary-dilation
data.  No such reduction is assumed here.
