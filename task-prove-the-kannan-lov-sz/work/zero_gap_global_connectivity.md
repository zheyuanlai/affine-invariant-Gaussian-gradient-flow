# Zero-gap power diagrams: global corridor inequalities and their limit

## 0. Scope and verdict

Let \(K\subset\mathbb R^d\) be a convex body whose uniform probability is
isotropic.  Let

\[
 R=\mathop{\dot\bigcup}_{i=1}^N R_i,
 \qquad |R|={|K|\over2},
 \qquad b(z)=b_i\quad(z\in R_i),                         \tag{0.1}
\]

and assume that

\[
 (b_i-b_j)\cdot(z-z')\ge |b_i-b_j|^2
 \quad(z\in R_i,\ z'\in R_j).                          \tag{0.2}
\]

Thus \(b\) is firmly nonexpansive on its plateau set.  Assume also that
the endpoint pieces form a partition

\[
 K=\mathop{\dot\bigcup}_i(R_i-b_i)
   \mathop{\dot\bigcup}_i(R_i+b_i)                     \tag{0.3}
\]

modulo null sets.  These are the finite midpoint data of the convex
piecewise-translation branch of a zero-entropy-gap Brenier map.

This note proves three dimension-free subcases.

1.  In an orthogonal product-corridor complex, the sum of all coordinate
    corridor widths is at most \(\log 2\), and consequently
    \(\mathbb E|b|^2\le 1+3(\log2)^2\).
2.  More generally, if the local corridors can be organized as sequential
    Minkowski extrusions of one convex half-volume core, then their total
    translation length is at most \(\sqrt{24}\), with no assumption that
    the extrusion directions commute.  Any label contained in the resulting
    zonotope has \(\mathbb E|b|^2\le13\).
3.  The symmetric regular-simplex fan and the symmetric crosspolytope fan
    both have universal label energy.  The simplex proof uses an exact
    Dirichlet union calculation; a union bound alone would incorrectly lose
    \(\log d\).

The missing step is now precise.  A general power complex supplies only
**local** wall prisms.  It need not contain a convex half-volume core to
which every wall jump can be applied as a global Minkowski extrusion.
Alexandrov--Fenchel positivity therefore cannot simply sum the local wall
data.  An arbitrary-complex theorem would have to manufacture a laminar
family of large convex cores, or prove a substitute mixed-volume inequality
for local facet patches.  Neither statement is proved below.  Accordingly,
this note does not claim the desired bound for every finite power diagram.

## 1. Identities that remove the mean-displacement issue

Let \(Z\) be uniform on \(R\), let \(B=b(Z)\), and let \(\varepsilon\) be an
independent symmetric sign.  By (0.3),

\[
 U=Z+\varepsilon B                                      \tag{1.1}
\]

is uniform on \(K\).  Since \(K\) is isotropic,

\[
 \mathbb EZ=0,
 \qquad
 I=\mathbb E(Z\otimes Z)+\mathbb E(B\otimes B).         \tag{1.2}
\]

In particular,

\[
 \mathbb E(B\otimes B)\preceq I.                       \tag{1.3}
\]

If \(m=\mathbb EB\ne0\), apply (1.3) to
\(u=m/|m|\).  Then

\[
 |m|^2=(\mathbb E\langle B,u\rangle)^2
 \le\mathbb E\langle B,u\rangle^2\le1.                \tag{1.4}
\]

Thus a proof only has to control \(\operatorname {tr}\operatorname {Cov}B\);
the mean label always costs at most one.  Notice also that

\[
 \|\operatorname {Cov}B\|_{\rm op}\le1.                \tag{1.5}
\]

Consequently a label field with centered rank \(r\) obeys

\[
 \mathbb E|B|^2\le r+1.                                 \tag{1.6}
\]

This disposes of every fixed-valence fan, but not of a complex whose label
rank tends to infinity.

## 2. A projection-thickness lemma from isotropy

The following elementary estimate is the input that turns mixed volume into
a dimension-free length bound.

**Lemma 2.1 (half-volume projection thickness).**  Let \(K\) be isotropic,
let \(A\subset K\) be measurable with \(|A|=\alpha|K|\), and let \(u\) be a
unit vector.  Then

\[
 { |A|\over |P_{u^\perp}A|}
 \le \sqrt{12\over\alpha}.                              \tag{2.1}
\]

In particular, if \(|A|\ge|K|/2\), its mean projection thickness in every
direction is at most \(\sqrt{24}\).

**Proof.**  Write the fibers of \(K\) parallel to \(u\) as intervals of
length \(L(y)\), \(y\in P_{u^\perp}K\).  On any interval of length \(L\),
the integral of \(t^2\) is at least \(L^3/12\), independently of the
location of the interval.  Isotropy in direction \(u\) therefore gives

\[
 \int L(y)^3\,dy\le12|K|.                               \tag{2.2}
\]

Put \(Q=P_{u^\perp}A\).  Since every \(A\)-fiber has length at most the
corresponding \(K\)-fiber, Holder's inequality yields

\[
 |A|\le\int_Q L(y)\,dy
 \le(12|K|)^{1/3}|Q|^{2/3}.                             \tag{2.3}
\]

Hence

\[
 |Q|\ge {|A|^{3/2}\over\sqrt{12|K|}},
\]

which is (2.1). \(\square\)

No isoperimetric or Poincare inequality is used in this lemma.

## 3. Sequential Brunn--Minkowski closes every global corridor system

For a convex body \(C\) and a vector \(v\), Cavalieri's principle gives the
exact one-segment Steiner formula

\[
 |C+[0,v]|=|C|+|v|\,|P_{v^\perp}C|.                    \tag{3.1}
\]

It is also the first mixed-volume coefficient
\(dV(C[d-1],[0,v])\).  All higher mixed-volume coefficients in a sum of
segments are nonnegative.

**Theorem 3.1 (sequential global-corridor bound).**  Suppose there are convex
bodies

\[
 C_0\subset C_1\subset\cdots\subset C_m\subset K,
 \qquad |C_0|={|K|\over2},                              \tag{3.2}
\]

and vectors \(v_1,\ldots,v_m\) such that

\[
 C_{j-1}+[0,v_j]\subset C_j\qquad(1\le j\le m).        \tag{3.3}
\]

Then

\[
 \boxed{\sum_{j=1}^m|v_j|\le\sqrt{24}.}                \tag{3.4}
\]

The directions \(v_j\) may be arbitrary and need not be orthogonal or
commuting.

**Proof.**  By (3.1), (3.3), and Lemma 2.1,

\[
 |C_j|-|C_{j-1}|
 \ge |v_j|\,|P_{v_j^\perp}C_{j-1}|
 \ge {|v_j|\over\sqrt{24}}|C_{j-1}|
 \ge {|v_j|\over2\sqrt{24}}|K|.                        \tag{3.5}
\]

Sum in \(j\).  The left side telescopes and is at most
\(|K|-|C_0|=|K|/2\), proving (3.4). \(\square\)

**Corollary 3.2 (zonotopal labels).**  Under the hypotheses of Theorem 3.1,
suppose that, for some \(b_0\), every label belongs to

\[
 b_0+\sum_{j=1}^m[0,v_j].                               \tag{3.6}
\]

Then

\[
 \operatorname {tr}\operatorname {Cov}B\le12,
 \qquad \mathbb E|B|^2\le13.                           \tag{3.7}
\]

Indeed, (3.4) makes the diameter of the zonotope in (3.6) at most
\(\sqrt{24}\), and

\[
 \operatorname {tr}\operatorname {Cov}B
 ={1\over2}\mathbb E|B-B'|^2\le12.                    \tag{3.8}
\]

The last unit in (3.7) is (1.4).

### 3.1 Simultaneous Minkowski addition

There is an equivalent one-step version.  If

\[
 C_0+\sum_{j=1}^m[0,v_j]\subset K,
 \qquad |C_0|=|K|/2,                                   \tag{3.9}
\]

then the Steiner polynomial and positivity of all mixed volumes imply

\[
 \sum_j |v_j|\,|P_{v_j^\perp}C_0|
 \le |K|-|C_0|.                                        \tag{3.10}
\]

Lemma 2.1 again gives \(\sum_j|v_j|\le\sqrt{24}\).
This is the exact Alexandrov--Fenchel mechanism suggested by the product
examples: first-order corridor costs add, while all noncommuting cross terms
have the favorable sign.

The limitation is the common-core premise (3.9), not the angles between the
corridors.

## 4. Orthogonal product corridors

The product calculation makes the multiplicative volume mechanism fully
explicit.

Let

\[
 K=[-\sqrt3,\sqrt3]^d,
 \qquad R=R_1\times\cdots\times R_d,                    \tag{4.1}
\]

where each \(R_j\subset[-\sqrt3,\sqrt3]\) is a finite union of intervals.
Put

\[
 \delta_j=1-{|R_j|\over2\sqrt3}.                       \tag{4.2}
\]

Assume \(|R|=|K|/2\).  Let

\[
 b(z)=(\beta_1(z_1),\ldots,\beta_d(z_d)),               \tag{4.3}
\]

where \(\beta_j\) is constant on every component of \(R_j\) and satisfies
the scalar firm inequality

\[
 (\beta_j(s)-\beta_j(t))(s-t)
 \ge(\beta_j(s)-\beta_j(t))^2.                         \tag{4.4}
\]

Assume finally that the endpoint mixture (1.1) is uniform on \(K\).

**Theorem 4.1 (orthogonal corridor bound).**  Under these assumptions,

\[
 \operatorname {tr}\operatorname {Cov}B
 \le3(\log2)^2,
 \qquad
 \boxed{\mathbb E|B|^2\le1+3(\log2)^2.}                \tag{4.5}
\]

**Proof.**  The scalar firm inequality says that each \(\beta_j\) is
nondecreasing and 1-Lipschitz on \(R_j\).  Since it is constant on the
occupied interval components, its total increase is at most the total length
of the omitted gaps:

\[
 \operatorname {range}(\beta_j)
 \le2\sqrt3\,\delta_j.                                 \tag{4.6}
\]

Consequently Popoviciu's elementary variance bound gives

\[
 \operatorname {Var}(\beta_j(Z_j))
 \le {\operatorname {range}(\beta_j)^2\over4}
 \le3\delta_j^2.                                       \tag{4.7}
\]

The half-volume identity is

\[
 \prod_{j=1}^d(1-\delta_j)={1\over2}.                  \tag{4.8}
\]

Since \(\delta\le-\log(1-\delta)\),

\[
 \sum_j\delta_j\le\log2.                              \tag{4.9}
\]

Therefore

\[
 \operatorname {tr}\operatorname {Cov}B
 =\sum_j\operatorname {Var}(\beta_j(Z_j))
 \le3\sum_j\delta_j^2
 \le3(\log2)^2.                                       \tag{4.10}
\]

Add (1.4). \(\square\)

Thus the canonical rank-many cube attempt cannot work: if coordinate
corridors have widths of order \(1/d\), their total missing volume is of
constant order but their squared label energy is only of order \(1/d\).
More uneven widths only reduce \(\ell_2\) relative to their fixed \(\ell_1\)
budget.

### 4.1 Why tensoring a one-dimensional example does not evade (4.8)

Take a one-dimensional central-source translation model and tensor it in
\(d\) coordinates.  The set on which all \(d\) forward branches are
available has mass \(2^{-d}\), not \(1/2\).  Replacing it by an even-parity
half uses inverse branches in some coordinates.  In one dimension the inverse
branch reverses the order of the displacement labels, so the resulting field
violates monotonicity and hence cannot be a Brenier gradient.

There is also a purely set-theoretic obstruction.  If
\(E=\prod E_j\) and \(F=\prod F_j\) are complementary product sets, then at
most one coordinate can have both \(E_j\) and \(F_j\) of positive measure.
Otherwise a mixed point, using an \(E\)-coordinate in one position and an
\(F\)-coordinate in another, belongs to neither product.  Thus an exact
product complement has translation rank one.  Rank-many constructions must
use genuine nonproduct corridors, which is exactly where (4.8) is lost.

## 5. A symmetric regular-simplex fan

This section tests the most acute high-rank fan.  Let
\(v_0,\ldots,v_d\) be the vertices of the regular isotropic simplex, with

\[
 |v_i|^2=d(d+2),
 \qquad v_i\cdot v_j=-(d+2)\quad(i\ne j).               \tag{5.1}
\]

Write every point as \(x=\sum_{i=0}^d p_i v_i\), where \(p\) is uniform on
the probability simplex.  Suppose there are \(d+1\) equally weighted source
pieces \(P_i\), and the label on \(P_i\) is

\[
 b_i=\lambda v_i,
 \qquad \lambda\ge0.                                   \tag{5.2}
\]

Assume \(P_i\) and \(P_i+2b_i\) lie in \(K\), the \(P_i\) are disjoint,
and their union has volume \(|K|/2\).  These conditions are weaker than full
complementarity, so the resulting bound applies to the symmetric fan a
fortiori.

Adding \(2\lambda v_i\) changes barycentric coordinates by

\[
 p_i\mapsto p_i+{2\lambda d\over d+1},
 \qquad
 p_j\mapsto p_j-{2\lambda\over d+1}\quad(j\ne i).      \tag{5.3}
\]

Thus

\[
 P_i\subset A_i:=\left\{p:p_j\ge c\ \hbox{for all }j\ne i\right\},
 \qquad c={2\lambda\over d+1}.                         \tag{5.4}
\]

The union \(\cup_iA_i\) is the event that at most one barycentric coordinate
is smaller than \(c\).  Put \(m=d+1\), and let
\(N=\#\{i:p_i<c\}\).  A direct translation of the simplex gives

\[
 \mathbb P(N=0)=(1-mc)_+^d,                             \tag{5.5}
\]

and integration in the exceptional coordinate gives

\[
 \mathbb P(N=1)
 =m\left[(1-(m-1)c)_+^d-(1-mc)_+^d\right].             \tag{5.6}
\]

Therefore

\[
 \left|\bigcup_iA_i\right|/|K|
 =m(1-(m-1)c)_+^d-(m-1)(1-mc)_+^d.                    \tag{5.7}
\]

Since the disjoint source pieces of total volume \(|K|/2\) lie in this
union, the right side is at least \(1/2\).

Set

\[
 a=1-{2\lambda d\over d+1},
 \qquad q=1-2\lambda.
\]

If \(\lambda\ge1/2\), then \(q\le0\), and (5.7) is at most
\[
 (d+1)\left(1-{d\over d+1}\right)^d=(d+1)^{1-d}<1/2
 \tag{5.7a}
\]
for \(d\ge2\).  Hence the half-volume premise forces \(\lambda<1/2\), so
the positive parts in (5.7) may now be removed.  By the mean-value theorem,

\[
\begin{split}
 (d+1)a^d-dq^d
 &=a^d+d(a^d-q^d)\\
 &\le a^{d-1}\left(a+{2\lambda d^2\over d+1}\right)\\
 &\le (1+x)e^{-x},
 \qquad x={2\lambda d(d-1)\over d+1}.
\end{split}
\tag{5.8}
\]

Since \(3e^{-2}<1/2\), the last expression being at least \(1/2\) forces
\(x<2\), and therefore \(\lambda d<3\) for every \(d\ge2\).  The case
\(d=1\) is one-dimensional.  Hence, for example,

\[
 \mathbb E|B|^2=\lambda^2d(d+2)\le18.                  \tag{5.9}
\]

This calculation identifies a real audit trap.  Bounding each
\(|A_i|\) separately and then taking a union bound gives only
\(\lambda d\lesssim\log d\).  The exact overlap term in (5.7) removes that
logarithm.  In geometric language, all vertex corridors compete for the same
central simplex volume; they are not independent low-area escapes.

## 6. A symmetric crosspolytope fan

Let

\[
 K=\{x:\|x\|_1\le R\},
 \qquad R=\sqrt{(d+1)(d+2)/2},                           \tag{6.1}
\]

so that the uniform law on \(K\) is isotropic.  Consider the symmetric
outward fan with labels

\[
 b_{\sigma i}={t\over2}\sigma e_i,
 \qquad \sigma\in\{-1,1\},\quad1\le i\le d,            \tag{6.2}
\]

and source pieces satisfying \(\sigma x_i\ge0\).  If both \(x\) and
\(x+t\sigma e_i\) lie in \(K\), then

\[
 \|x\|_1\le R-t.                                       \tag{6.3}
\]

Thus the union of all source pieces is contained in
\((1-t/R)K\).  If it has half the volume of \(K\),

\[
 (1-t/R)^d\ge{1\over2},
 \qquad
 t\le R(1-2^{-1/d})\le {R\log2\over d}.                \tag{6.4}
\]

Since \(R/d\le2\) for all \(d\ge1\),

\[
 \mathbb E|B|^2={t^2\over4}\le(\log2)^2.              \tag{6.5}
\]

Again the high fan valence is harmless: convexity makes all outward pieces
draw their mass from one homothetic inner core.

## 7. A noncommuting three-wall fan

Three labels \(b_1,b_2,b_3\) need not be collinear.  Nevertheless their
centered covariance has rank at most two.  Equations (1.4)--(1.5) give the
sharp dimension-free audit

\[
 \mathbb E|B|^2
 =\operatorname {tr}\operatorname {Cov}B+|\mathbb EB|^2
 \le2+1=3.                                              \tag{7.1}
\]

Thus noncommutativity at a single triple junction is not an obstruction.
Growing cost requires growing label rank and hence a globally coupled family
of junctions.

It is important, however, not to infer a false determinant estimate from the
triple junction.  Take

\[
 b_1=0,
 \qquad b_2=Le_1,
 \qquad b_3=L(e_1+\varepsilon e_2).                     \tag{7.2}
\]

The two-dimensional junction prism has label area
\(\varepsilon L^2/2\), while the label diameter is comparable to \(L\).
Choosing \(\varepsilon=L^{-3}\) makes the junction volume tend to zero as
the diameter tends to infinity.  The codimension-one wall prisms must then
carry the missing charge, unless their physical facets are correspondingly
short.  Hence the sum of higher-stratum determinant volumes alone cannot
control label energy.  A valid global inequality must combine all wall
strata with the masses and connectivity of their adjacent cells.

The sequential theorem does handle this fan if a common half-volume core
\(C\) satisfies

\[
 C+[0,b_2-b_1]+[0,b_3-b_1]\subset K.                   \tag{7.3}
\]

In that case (3.10) controls both edge lengths, irrespective of
\(\varepsilon\).  General power data provide only cellwise prisms and do not
imply (7.3).

## 8. What the wall-prism budget does and does not give

For a source--source wall \(\Gamma_{ij}\), the exact midpoint prism has
volume

\[
 |b_i-b_j|\,\mathcal H^{d-1}(\Gamma_{ij}),              \tag{8.1}
\]

and all open wall prisms are disjoint.  Thus

\[
 \sum_{i<j}|b_i-b_j|\,
 \mathcal H^{d-1}(\Gamma_{ij})\le {|K|\over2}.          \tag{8.2}
\]

This is a local first-variation budget.  Compare it with the global
first-variation term in (3.10):

\[
 |v|\,|P_{v^\perp}C|.                                  \tag{8.3}
\]

Lemma 2.1 makes (8.3) coercive because \(C\) has half the volume of \(K\).
In (8.2), the facet \(\Gamma_{ij}\) may be an arbitrarily small patch.  No
known consequence of isotropy lower-bounds its area in terms of either
adjacent cell mass.  Replacing it by such a lower bound is a relative
isoperimetric assertion and cannot be assumed.

The exact arbitrary-complex statement needed by this route can be formulated
as follows.

> **Local-to-global corridor conjecture.**  From every finite data set
> satisfying (0.1)--(0.3), construct a probability distribution over laminar
> chains of convex sets \(C_0\subset\cdots\subset C_m\subset K\), with
> \(|C_0|\ge c|K|\), and corridor vectors \(v_j\) satisfying
> \(C_{j-1}+[0,v_j]\subset C_j\), such that for the random label \(B\),
> \[
> \mathbb E|B-B'|^2
> \le C\,\mathbb E_{\rm chains}\left(\sum_j|v_j|\right)^2. \tag{8.4}
> \]

Theorem 3.1 would immediately make the right side universal.  The content is
the construction of the chains: (8.2) alone does not supply them.  Equivalently,
one may seek a localized Alexandrov--Fenchel inequality that replaces the
full projection area in (8.3) by the collection of facet patches in (8.2)
without losing a dimension or a graph Cheeger constant.

This is not a harmless reformulation.  A long chain of nearly parallel cells
is controlled by (1.5), an orthogonal product complex is controlled by
Theorem 4.1, and a symmetric fan is controlled by Sections 5--6.  A putative
bad complex must mix these geometries: it needs growing centered label rank,
small local wall patches, no common convex core, and repeated changes of
corridor direction.  Those are precisely the hypotheses under which a
laminar extraction theorem is presently missing.

## 9. Model registry

* **Cube down-sets and coordinate staircases.**  Product corridors satisfy
  Theorem 4.1.  Tensoring forward one-dimensional branches gives the wrong
  mass, while parity repair violates monotonicity.  Exact enumeration of
  congruent-cell matchings on the \(2\times2\), \(2\times3\), and
  \(4\times2\) grids yields only one-direction slab translations; this is a
  finite audit, not a proof for every grid.
* **Regular simplex.**  For the symmetric vertex fan, the exact event
  \(N\le1\) in (5.7), rather than a union bound over vertices, forces
  \(\lambda d=O(1)\) and hence universal energy.
* **Crosspolytope.**  Every symmetric outward fan draws from a common inner
  homothetic crosspolytope, and the half-volume condition gives (6.4).
* **Zonotopal tilings.**  Whenever the tiling corridors arise by Minkowski
  addition to one half-volume convex core, all mixed terms are nonnegative
  and Theorem 3.1 gives the bound.  A generic zonotopal *local* tiling need
  not exhibit such a core; calling every wall segment a Minkowski summand
  would be the invalid step.
* **Three-wall noncommuting junction.**  Fixed rank already gives (7.1), but
  the degenerating triangle (7.2) disproves any attempt to use only the
  top-dimensional label-polytope volume.

## 10. Conclusion

The global inequality is proved whenever rank-many corridors incur
multiplicative volume loss on a common large core.  The constants are
explicit:

\[
 \sum |v_j|\le\sqrt{24},
 \qquad
 \operatorname {tr}\operatorname {Cov}B\le12,
 \qquad
 |\mathbb EB|^2\le1.                                   \tag{10.1}
\]

The simplex calculation shows that substantial overlap among apparent
escape corridors can itself substitute for a common core, and it removes a
spurious \(\log d\) loss exactly.  What remains is to prove that every
Brenier-valid complementary power complex has one of these two global
features, or to exhibit a genuinely growing-rank counterexample avoiding
both.  Local wall-prism volume, fixed-rank covariance, and pairwise overlap
bounds do not by themselves decide that alternative.
