# Signed-distance polarization: the exact monotonicity and its stopping point

## 0. Result of the audit

For a probability measure \(\mu\) on a Euclidean affine space and an open
set \(E\), write

\[
 J_\mu(E)=\int d(x,\partial E)\,d\mu(x).
\]

There is one useful exact polarization statement. If reflection in a
hyperplane preserves \(\mu\), then polarizing **the set only** toward either
side cannot decrease \(J_\mu\). More precisely, polarization commutes with
inner parallel sets in the favorable direction, and the same assertion
applied to the complement controls the two signs of the signed distance:

\[
 \boxed{\quad J_\mu(P_HE)\ge J_\mu(E),\qquad
        \mu(P_HE)=\mu(E).\quad}                         \tag{0.1}
\]

This is proved in Sections 1--2 without a perimeter or smoothness
assumption. The covariance does not change, because the measure does not
change. The result is sharp: cube halfspaces and radial spheres are fixed
points.

It does not give a symmetrization proof of KLS. The reasons are exact.

1. A general log-concave measure has no reflection group on which (0.1)
   can operate. Two-point polarization of the density need not preserve
   log-concavity, even for a smooth strictly log-concave density on the
   line; Proposition 4.1 gives an explicit convex potential for which the
   polarized potential has a downward jump in its derivative.
2. Polarizing a density and a balanced set in the same direction does not
   preserve balance. For the one-sided exponential, a median halfline of
   mass \(1/2\) becomes a set of mass exactly \(3/4\) (or \(1/4\) in the
   opposite orientation).
3. Steiner rearrangement of a log-concave density does preserve
   log-concavity, but joint rearrangement and subsequent rebalancing has no
   monotonicity after covariance normalization. A one-parameter family of
   truncated log-affine laws gives both strict signs. Smooth strictly
   log-concave approximations retain both signs.
4. Even for a uniform triangle, a natural simultaneous Steiner
   symmetrization of the body and a half-volume cap reduces the whitened
   signed-distance integral by the exact factor \(1/\sqrt3\).
5. Iteration can only use genuine symmetries. Coordinate polarizations of
   an unconditional law preserve the occupancy count on every sign orbit;
   origin-hyperplane polarizations of a radial law preserve the angular
   measure of every radial section. They therefore cannot reduce an
   arbitrary balanced set to a halfspace, a ball, or its complement.

Thus polarization is a valid monotone improvement inside a fixed symmetry
class, but there is no dimension-free reduction from general log-concave
measures to unconditional, radial, or halfspace structure. Sections 7--9
track covariance, whitening, compactness, and all requested model tests.

## 1. Definitions and the erosion identity for \(J\)

Let \(L\) be an affine hyperplane, let \(\sigma\) be reflection in \(L\),
and let \(H\) be one of the two closed halfspaces bounded by \(L\). Values
on \(L\) are immaterial in all measure statements below. The two-point
polarization of a Borel set is

\[
 P_HE=\bigl((E\cup\sigma E)\cap H\bigr)
       \cup\bigl((E\cap\sigma E)\cap H^c\bigr).          \tag{1.1}
\]

Thus, on each two-point orbit \(\{x,\sigma x\}\), a lone occupied point is
moved into \(H\). Write

\[
 E^{\ominus t}=\{x:B(x,t)\subset E\},\qquad t>0,         \tag{1.2}
\]

where open balls may be replaced by closed balls after changing strict to
nonstrict inequalities. For an open set \(E\),

\[
 \{d(x,\partial E)>t\}=E^{\ominus t}\,\uplus\,(E^c)^{\ominus t}.  \tag{1.3}
\]

Consequently Tonelli's theorem gives

\[
 J_\mu(E)=\int_0^\infty
   \left[\mu(E^{\ominus t})+\mu((E^c)^{\ominus t})\right]dt.       \tag{1.4}
\]

Formula (1.4) remains valid for any set for which the boundary and the two
metric interiors are interpreted using a fixed topological representative.

### Lemma 1.1 (polarization improves every inner parallel set)

For every set \(E\subset\mathbb R^n\) and every \(t>0\),

\[
              P_H(E^{\ominus t})\subset (P_HE)^{\ominus t}.       \tag{1.5}
\]

#### Proof

Put coordinates so that \(H=\{x_1\ge0\}\). We use the elementary
reflection comparison

\[
\begin{aligned}
 |u-v|&\le|\sigma u-v| &&\text{if \(u,v\) lie on the same side of \(L\)},\\
 |\sigma u-v|&\le|u-v| &&\text{if \(u,v\) lie on opposite sides of \(L\)}
\end{aligned}                                             \tag{1.6}
\]

when \(\sigma u\) is on the side of \(v\). It follows immediately by
comparing the squared normal coordinates.

First let \(x\in H\) and suppose \(x\in P_H(E^{\ominus t})\). At least one
of \(x,\sigma x\) is the center of a radius-\(t\) ball contained in \(E\).
Assume first that \(B(x,t)\subset E\). If \(z\in B(x,t)\cap H\), then
\(z\in E\), hence \(z\in P_HE\). If \(z\in B(x,t)\cap H^c\), then both
\(z\) and \(\sigma z\) lie in \(B(x,t)\), the latter by (1.6); hence both
belong to \(E\), and again \(z\in P_HE\). If instead
\(B(\sigma x,t)\subset E\), reflect this argument. For a point
\(z\in H^c\cap B(x,t)\), both \(z\) and \(\sigma z\) are within distance
\(t\) of \(\sigma x\): one equality follows from reflection and the other
is (1.6). Thus \(B(x,t)\subset P_HE\).

Now let \(x\in H^c\). Membership in \(P_H(E^{\ominus t})\) says that both
\(B(x,t)\) and \(B(\sigma x,t)\) lie in \(E\). If
\(z\in B(x,t)\cap H^c\), then \(z\in E\) and
\(\sigma z\in B(\sigma x,t)\subset E\). If \(z\in B(x,t)\cap H\), one of
\(z,\sigma z\) belongs to \(E\) by the same reflection contraction. These
are exactly the intersection and union conditions in (1.1). Hence
\(B(x,t)\subset P_HE\). This proves (1.5). \(\square\)

The complement of a polarization is the oppositely oriented polarization
of the complement:

\[
                    (P_HE)^c=P_{H^c}(E^c).               \tag{1.7}
\]

This is why both terms in (1.4), rather than only the interior term, have
the same favorable monotonicity.

## 2. The exact fixed-measure theorem

### Theorem 2.1 (reflection-invariant signed-distance polarization)

Let \(\mu\) be a probability measure invariant under \(\sigma\). Assume
that \(J_\mu(E)<\infty\); in particular this holds when \(\mu\) has a finite
first moment and \(\partial E\ne\varnothing\). Then

\[
 \mu(P_HE)=\mu(E),\qquad J_\mu(P_HE)\ge J_\mu(E).         \tag{2.1}
\]

The assertion holds for open sets by direct use of topological boundary,
and for measurable sets after fixing the metric-interior representative.

#### Proof

On each reflection pair, polarization preserves the number of occupied
points. Reflection invariance therefore gives

\[
                         \mu(P_HA)=\mu(A)                 \tag{2.2}
\]

for every measurable \(A\). Lemma 1.1 and (2.2) imply

\[
 \mu((P_HE)^{\ominus t})
 \ge\mu(P_H(E^{\ominus t}))=\mu(E^{\ominus t}).          \tag{2.3}
\]

Apply the same argument to \(E^c\), oriented toward \(H^c\), and use
(1.7):

\[
 \mu(((P_HE)^c)^{\ominus t})\ge\mu((E^c)^{\ominus t}).  \tag{2.4}
\]

Integration in \(t\) using (1.4) proves (2.1). \(\square\)

### Corollary 2.2 (no covariance or isotropy loss)

If \(\mu\) is centered or isotropic, the same is true before and after the
operation in Theorem 2.1, because only \(E\) changes. In particular, a
finite chain of polarizations through genuine symmetry hyperplanes of
\(\mu\) preserves the half-mass constraint and increases \(J_\mu\) at
every step with **no affine factor**.

This is the maximal exact positive statement available from two-point
polarization. Every extension in which the measure is also changed must
separately reprove balance, log-concavity, covariance control, and the
metric comparison.

## 3. What the fixed-measure theorem can and cannot symmetrize

### 3.1 Unconditional measures

Suppose \(\mu\) is invariant under all coordinate sign changes. For
\(r=(|x_1|,\ldots,|x_n|)\), let

\[
 N_E(r)=\#\{\varepsilon\in\{-1,1\}^n:
                   (\varepsilon_1r_1,\ldots,\varepsilon_nr_n)\in E\}. \tag{3.1}
\]

Every coordinate polarization merely compresses the occupied vertices of
the Boolean sign cube and preserves \(N_E(r)\) for almost every \(r\).
Repeated coordinate compressions can produce an increasing Boolean family
on each sign orbit, but they cannot change the integer-valued function
\(N_E(r)\). A halfspace, an orthant, and a radial set generally have
different functions (3.1). Thus even in the unconditional class,
polarization does not reduce all balanced interfaces to one canonical
shape.

### 3.2 Radial measures

Suppose \(\mu\) is radial. Its available reflection symmetries are the
hyperplanes through the origin. If

\[
 p_E(r)=\frac{\mathcal H^{n-1}(E\cap rS^{n-1})}
                    {\mathcal H^{n-1}(rS^{n-1})},         \tag{3.2}
\]

then every allowed polarization preserves \(p_E(r)\) for almost every
radius. A dense sequence of angular polarizations may turn each angular
section into a spherical cap about a common pole, but the arbitrary
function \(r\mapsto p_E(r)\) survives. In particular,

\[
 p_{\{x_1>0\}}(r)=\tfrac12,\qquad
 p_{\{|x|<r_0\}}(r)=\mathbf1_{(0,r_0)}(r),                \tag{3.3}
\]

so a median sphere cannot be polarized into a halfspace, or conversely.
Theorem 2.1 therefore does not reduce the radial problem to either one.

### 3.3 General measures

For a generic log-concave density the symmetry group is trivial. Theorem
2.1 then has no nontrivial direction. The tempting response is to
polarize the density along with the set. Sections 4--6 show that this
fails before any KLS estimate is reached.

## 4. Two-point polarization of the density is inadmissible

For a nonnegative density \(\rho\), its same-side two-point polarization is

\[
 \rho^H(x)=
 \begin{cases}
  \max\{\rho(x),\rho(\sigma x)\},&x\in H,\\
  \min\{\rho(x),\rho(\sigma x)\},&x\in H^c.
 \end{cases}                                             \tag{4.1}
\]

It preserves total mass and all integrals of reflection-invariant test
functions. It does **not** preserve log-concavity.

### Proposition 4.1 (explicit log-concavity failure)

On \(\mathbb R\), define the convex coercive potential

\[
 V(x)=
 \begin{cases}
  -3(x+1),&x\le-1,\\
  0,&-1\le x\le0,\\
  x,&x\ge0.
 \end{cases}                                             \tag{4.2}
\]

Let \(\rho=Z^{-1}e^{-V}\), reflect in the origin, and polarize the larger
density toward \(H=[0,\infty)\). Then \(\rho^H\) is not log-concave.

#### Proof

Writing \(W=-\log(Z\rho^H)\), for \(r>0\) one has

\[
 W(r)=\min\{V(r),V(-r)\}
 =\begin{cases}
    0,&0<r<1,\\
    3(r-1),&1<r<3/2,\\
    r,&r>3/2.
   \end{cases}                                           \tag{4.3}
\]

The left derivative of \(W\) at \(3/2\) is \(3\), whereas the right
derivative is \(1\). Hence \(W\) is not convex, so \(\rho^H\) is not
log-concave. \(\square\)

This is not an artifact of corners. Let \(\eta_\varepsilon\) be an even
standard \(C^\infty\) mollifier supported in
\((-\varepsilon,\varepsilon)\), and put

\[
 V_\varepsilon=V*\eta_\varepsilon+\varepsilon x^2.       \tag{4.4}
\]

Then \(V_\varepsilon\) is smooth, strictly convex, and coercive. If
\(0<\varepsilon<1/4\), both points in a neighborhood of \(3/2\), and their
negatives, remain farther than the mollifier radius from the corners of
\(V\). Evenness of the mollifier makes convolution exact on those affine
pieces. The two branches cross at \(r=3/2\), with derivatives
\(3+2\varepsilon r\) on the left and \(1+2\varepsilon r\) on the right.
Their lower envelope has a downward derivative jump of exactly \(2\).
Thus every choice \(0<\varepsilon<1/4\) in (4.4) gives an explicit smooth
strictly log-concave counterexample.

The issue is stronger than mere loss of smoothness. Polarization keeps
convexity of all individual superlevel sets in this one-dimensional
example, but log-concavity requires the quantitative interpolation between
different levels; (4.3) violates exactly that interpolation.

## 5. Joint two-point polarization destroys the half-mass constraint

Even in cases where a particular polarized density happens to remain
log-concave, sorting the density and the set in the same direction changes
their association. On a reflection pair with density values \(a\ge b\)
and exactly one occupied point, same-side sorting changes the occupied
mass from either \(a\) or \(b\) to \(a\); opposite-side sorting changes it
to \(b\). There is no pairwise conservation law for \(\mu(E)\).

### Proposition 5.1 (the exact \(3/4\)--\(1/4\) exponential test)

Let \(X\) have density \(e^{-x}\mathbf1_{[0,\infty)}(x)\). Its median is
\(m=\log2\). In the coordinate \(y=x-m\), let

\[
 E=\{y>0\},\qquad H=\{y>0\}.                             \tag{5.1}
\]

Then \(\mu(E)=1/2\). Polarize the density so that its larger value on each
pair \(\{y,-y\}\) lies in \(H\). Although \(P_HE=E\),

\[
                          \mu^H(E)=\frac34.               \tag{5.2}
\]

Putting the larger value in \(H^c\) instead gives mass \(1/4\).

#### Proof

The translated density is

\[
 \rho(y)=\tfrac12e^{-y}\mathbf1_{[-m,\infty)}(y).        \tag{5.3}
\]

For \(0<y<m\), \(\rho(-y)=\tfrac12e^y\) is larger than
\(\rho(y)=\tfrac12e^{-y}\); for \(y>m\), the reflected value is zero.
Therefore

\[
 \mu^H(E)=\frac12\int_0^m e^y\,dy
             +\frac12\int_m^\infty e^{-y}\,dy
          =\frac12(2-1)+\frac14=\frac34.                 \tag{5.4}
\]

The opposite orientation has the complementary mass. \(\square\)

Notice that \(d(y,\partial E)=|y|\) is reflection invariant. Hence its
integral is unchanged by density polarization and equals \(\log2\), even
though the admissibility constraint has been lost. A subsequent shift of
the boundary to restore half mass is a new global operation and has no
polarization monotonicity.

## 6. Steiner symmetrization: admissible density, indefinite normalized gain

Two-point density polarization should not be confused with Steiner
rearrangement. Write points as \((z,t)\in e^\perp\times\mathbb R\). The
Steiner rearrangement \(\rho^*\) replaces, for each fixed \(z\), the
function \(t\mapsto\rho(z,t)\) by its even nonincreasing equimeasurable
rearrangement.

For a log-concave \(\rho\), \(\rho^*\) is log-concave. One proof applies
one-dimensional symmetric rearrangement to the fibers of every superlevel
set and uses

\[
 S_e(\lambda A+(1-\lambda)B)
 \supset \lambda S_eA+(1-\lambda)S_eB                 \tag{6.1}
\]

for convex sets \(A,B\); the log-concave interpolation of the original
superlevels then passes to the rearranged superlevels. Total mass and the
\(z\)-marginal are preserved.

The set constraint and Euclidean metric remain problematic. The following
one-dimensional calculation gives both signs after the standard
"symmetrize, rebalance at the median, and whiten" repair.

### Proposition 6.1 (both monotonicity directions in a log-affine family)

For \(a>0\), let

\[
 \rho_a(x)=\frac{a e^{ax}}{e^a-1}\mathbf1_{[0,1]}(x).    \tag{6.2}
\]

Its Steiner decreasing rearrangement is

\[
 \rho_a^*(y)=\frac{a e^{a(1-2|y|)}}{e^a-1}
                     \mathbf1_{[-1/2,1/2]}(y).           \tag{6.3}
\]

Let \(m_a\) be the median of \(\rho_a\), and compare the balanced
halflines \(E_a=(m_a,\infty)\) and \(E_a^*=(0,\infty)\). Define the
whitened signed-distance values

\[
 R(a)=\frac{J_{\rho_a}(E_a)}{\sqrt{\operatorname{Var}_{\rho_a}X}},
 \qquad
 R^*(a)=\frac{J_{\rho_a^*}(E_a^*)}
                  {\sqrt{\operatorname{Var}_{\rho_a^*}Y}}.       \tag{6.4}
\]

Then \(R^*(a)<R(a)\) for all sufficiently small positive \(a\), whereas
\(R^*(a)>R(a)\) for all sufficiently large \(a\).

#### Proof

Direct integration gives

\[
 m_a=\frac1a\log\frac{e^a+1}{2},                         \tag{6.5}
\]

\[
 M_a:=J_{\rho_a}(E_a)
 =\frac{a e^a-(e^a+1)\log((e^a+1)/2)}{a(e^a-1)},          \tag{6.6}
\]

\[
 v_a:=\operatorname{Var}_{\rho_a}X
 =\frac1{a^2}-\frac{e^a}{(e^a-1)^2},                     \tag{6.7}
\]

and

\[
 M_a^*=\frac{e^a-1-a}{2a(e^a-1)},\qquad
 v_a^*=\frac{2e^a-2-2a-a^2}{4a^2(e^a-1)}.                \tag{6.8}
\]

Expanding the difference of the squared normalized values at \(a=0\)
gives

\[
 \frac{M_a^2}{v_a}-\frac{(M_a^*)^2}{v_a^*}
       =\frac a{16}-\frac{23a^2}{960}+O(a^3).             \tag{6.9}
\]

It is positive for all sufficiently small \(a>0\). At the other end,
after the scaling \(a(1-X)\), the original law converges with moments to a
one-sided exponential, while \(2aY\) converges with moments to the
symmetric Laplace law of density \(e^{-|u|}/2\). Hence

\[
 \lim_{a\to\infty}R(a)=\log2,\qquad
 \lim_{a\to\infty}R^*(a)=\frac1{\sqrt2}.                 \tag{6.10}
\]

Since \(\log2<1/\sqrt2\), the reverse strict inequality holds for large
\(a\). \(\square\)

The counterexamples can be pinned to concrete parameters. Substitution in
(6.6)--(6.8) gives

\[
\begin{array}{c|cc}
 a&R(a)&R^*(a)\\ \hline
 1&0.85239936\ldots&0.82932334\ldots\\
 10&0.69428866\ldots&0.70775039\ldots .
\end{array}                                               \tag{6.11}
\]

For example, the displayed strict intervals are certified directly by
bounding the exponential and logarithm series with their first omitted
terms. Thus \(a=1\) and \(a=10\), rather than only limiting choices, give
the two signs.

Both signs persist for smooth strictly log-concave full-support laws.
Indeed, with

\[
 \phi_\varepsilon(s)=\varepsilon\log(1+e^{s/\varepsilon}),
\quad
 V_{a,M,\varepsilon}(x)
 =-ax+M\phi_\varepsilon(-x)
       +M\phi_\varepsilon(x-1)+\varepsilon x^2,           \tag{6.12}
\]

the normalized density \(e^{-V_{a,M,\varepsilon}}\) is smooth and strictly
log-concave. Taking first \(\varepsilon\downarrow0\) with
\(M\varepsilon\downarrow0\), and then \(M\to\infty\), gives (6.2) in
total variation with convergence of every fixed moment. Its symmetric
decreasing rearrangement converges in the same way to (6.3), because the
distribution functions of the density values converge at every nonflat
level. Medians, the integrals in (6.4), and variances are therefore
continuous. The strict inequalities supplied by (6.9) and (6.10) hold
for finite choices of \(M,\varepsilon\).

Equivalently, one may convolve (6.2) with a centered Gaussian of variance
\(\delta^2\). Convolution preserves log-concavity and produces a positive
real-analytic density. Symmetric rearrangement and all quantities in
(6.4) converge to their values above as \(\delta\downarrow0\); hence the
convolutions of the concrete \(a=1\) and \(a=10\) laws retain their
respective strict signs for every sufficiently small positive \(\delta\).

Proposition 6.1 also isolates the affine loss. Before whitening, Steiner
rearrangement decreases the displayed median-halfspace value in both the
small- and large-\(a\) regimes. In the large-\(a\) limit it changes
variance from \(a^{-2}\) to \((2a^2)^{-1}\); expansion by \(\sqrt2\) during
whitening reverses the normalized inequality.

## 7. Covariance and affine bookkeeping

### 7.1 Two-point density polarization

Assume the reflecting hyperplane is \(e^\perp\) through the origin. Pair
the points \(z+te\) and \(z-te\), \(t>0\). Two-point density
polarization preserves every reflection-even raw moment. In particular,

\[
 \int |x|^2\rho^H(x)dx=\int |x|^2\rho(x)dx,              \tag{7.1}
\]

as well as all purely tangential second moments and the normal raw second
moment. It need not preserve the barycenter or normal--tangential cross
moments. The new normal mean is

\[
 \langle b_{\rho^H},e\rangle
 =\int_{e^\perp}\int_0^\infty
       t\,|\rho(z+te)-\rho(z-te)|\,dt\,dz,                \tag{7.2}
\]

when the larger value is placed on the positive side. Thus a centered
density usually becomes noncentered, and its covariance trace drops by
the square of the new barycenter length (with the unchanged tangential
mean understood).

For the centered exponential

\[
 \rho(x)=e^{-(x+1)}\mathbf1_{[-1,\infty)}(x),             \tag{7.3}
\]

which has variance one, polarization toward the positive halfline gives

\[
 b_{\rho^H}=4e^{-2},\qquad
 \operatorname{Var}(\rho^H)=1-16e^{-4}.                  \tag{7.4}
\]

Even this elementary operation therefore needs a nontrivial whitening
factor. More importantly, Proposition 4.1 says the resulting law is not
in general in the log-concave class at all.

### 7.2 Steiner covariance

Steiner rearrangement preserves the tangential marginal, makes the normal
conditional distributions even, and minimizes their normal second moment
among equimeasurable fibers. If the symmetrizing hyperplane passes through
the barycenter, then

\[
 A_{\rho^*}=
 \begin{pmatrix}
  A_{\rm tan}&0\\
  0&v_e^*
 \end{pmatrix},\qquad
 v_e^*\le \int\langle x,e\rangle^2\,d\mu(x).             \tag{7.5}
\]

In particular, if \(\mu\) is isotropic before one Steiner step,

\[
                         A_{\rho^*}=I_{e^\perp}\oplus[\lambda],
                         \qquad0<\lambda\le1.             \tag{7.6}
\]

The operation is not isotropy preserving.

The covariance operator itself is monotone in the favorable scalar sense:
the largest eigenvalue of the original block matrix is at least both
\(\|A_{\rm tan}\|_{\rm op}\) and its normal diagonal entry, while
\(v_e^*\) is no larger than that entry. Hence

\[
                    \|A_{\rho^*}\|_{\rm op}
                    \le \|A_\rho\|_{\rm op}.              \tag{7.7}
\]

This does not control the smallest eigenvalue needed for whitening.

For an invertible affine map \(T\), Euclidean distance obeys

\[
 s_{\min}(T)d(x,\partial E)
 \le d(Tx,\partial(TE))
 \le\|T\|_{\rm op}d(x,\partial E).                       \tag{7.8}
\]

Consequently re-isotropization by \(A^{-1/2}\) can cost any factor between
\(\lambda_{\max}(A)^{-1/2}\) and
\(\lambda_{\min}(A)^{-1/2}\). Applying (7.7) separately at \(N\)
successive steps produces the unaudited factor

\[
                         \prod_{j=1}^N\lambda_j^{-1/2},   \tag{7.9}
\]

where \(\lambda_j\) is the newly contracted eigenvalue after the \(j\)-th
step. If \(N\) grows with dimension, even a universal one-step bound
\(\lambda_j\ge c_0<1\) gives \(c_0^{-N/2}\), not a KLS constant.

The product example of centered one-sided exponentials makes the danger
concrete. Symmetrizing one coordinate changes that coordinate's variance
from \(1\) to \(1/2\); whitening it costs \(\sqrt2\). Multiplying the
one-step estimate over \(n\) coordinate symmetrizations gives
\(2^{n/2}\). If all affine maps commute and whitening is postponed, the
actual final map is only \(\sqrt2 I\), so this particular exponential loss
is bookkeeping rather than geometry. It nevertheless proves that a
stepwise polarization proof is invalid unless it controls the **single
final product of affine maps**. With changing, nonorthogonal directions
the maps do not commute, and no such dimension-free control follows from
the one-step inequalities.

## 8. The two-dimensional simplex test

Let

\[
 K=\{(x,y):x\ge0,\ y\ge0,\ x+y\le1\}                    \tag{8.1}
\]

and let \(\mu\) be uniform probability on \(K\), with density \(2\).
Set

\[
 m=1-\frac1{\sqrt2},\qquad E=\{x<m\}.                    \tag{8.2}
\]

Since the \(x\)-marginal has density \(2(1-x)\), \(\mu(E)=1/2\), and

\[
 J_\mu(E)=2\int_0^1(1-x)|x-m|\,dx
          =\frac{2-\sqrt2}{3}.                           \tag{8.3}
\]

Steiner symmetrize vertical fibers in the \(x\)-direction. The body
becomes

\[
 K^*=\{(x,y):0\le y\le1,\ |x|\le(1-y)/2\}.              \tag{8.4}
\]

The sections of \(E\cap K\) have length \(\min\{m,1-y\}\). Centering
these sections gives \(K^*\cap E^*\), where the natural ambient extension
is the strip

\[
                         E^*=\{|x|<m/2\}.                \tag{8.5}
\]

This set still has probability \(1/2\). A direct fiber integration yields

\[
 J_{\mu^*}(E^*)
 =2\int_{K^*}\bigl||x|-m/2\bigr|\,dxdy
 =\frac{2-\sqrt2}{6}=\frac12J_\mu(E).                    \tag{8.6}
\]

The covariance matrices are

\[
 A_\mu=
 \begin{pmatrix}1/18&-1/36\\-1/36&1/18\end{pmatrix},
 \qquad
 A_{\mu^*}=\begin{pmatrix}1/24&0\\0&1/18\end{pmatrix}. \tag{8.7}
\]

After exact whitening, distance to the original vertical boundary is
multiplied by \(\sqrt{18}\), while distance to the two rearranged vertical
boundaries is multiplied by \(\sqrt{24}\). Therefore

\[
 \frac{J_{A_{\mu^*}^{-1/2}\#\mu^*}
              (A_{\mu^*}^{-1/2}E^*)}
      {J_{A_\mu^{-1/2}\#\mu}(A_\mu^{-1/2}E)}
 =\frac{\sqrt{24}}{2\sqrt{18}}=\frac1{\sqrt3}.           \tag{8.8}
\]

Thus even simultaneous volume-preserving Steiner symmetrization can
strictly decrease the whitened functional.

There is a representative issue hidden in every body-only formulation.
The set is only specified \(\mu\)-almost everywhere on \(K\), whereas
\(J\) uses its ambient boundary. Formula (8.5) is the canonical
fiber-centered ambient extension; choosing \(E^*\cap K^*\) instead adds
support-boundary pieces and decreases \(J\) further. A valid reduction
must prescribe this extension and prove that it is compatible with limits.

The strict example can be made smooth. Let

\[
 V_{M,\varepsilon}(x,y)=M\{\phi_\varepsilon(-x)
 +\phi_\varepsilon(-y)+\phi_\varepsilon(x+y-1)\}
 +\varepsilon(x^2+y^2),                                  \tag{8.9}
\]

where \(\phi_\varepsilon\) is as in (6.12). The normalized density
\(e^{-V_{M,\varepsilon}}\) is smooth and strictly log-concave. Along, for
example, \(M\to\infty\) and \(\varepsilon=M^{-2}\), it converges to
uniform measure on \(K\), with second moments. The associated fiber
rearranged densities converge to uniform measure on \(K^*\), and smooth
sets whose sections converge to those in (8.2) and (8.5) have the
signed-distance limits (8.3) and (8.6). Hence the numerical obstruction is
stable under smooth log-concave regularization. This last approximation is
not claimed to define a unique ambient-set Steiner operation: the
representative ambiguity in the preceding paragraph is precisely why no
such canonical operation is available. Any proposed relative-section
operation that is continuous under this regularization inherits the fixed
gap \(1/\sqrt3<1\).

## 9. Mandatory invariant-model tests

### 9.1 Tilted Gaussian nonlinear cut: strict increase

Let \((U,V)\) be any orthogonal rotation of a standard Gaussian pair, and
let

\[
                         E=\{UV>0\}.                     \tag{9.1}
\]

This is a balanced nonlinear two-cone cut. Reflect in the tilted line
\(\{U=0\}\) and polarize toward \(\{U>0\}\). Exactly one point in every
nondegenerate reflection pair belongs to \(E\), hence

\[
                         P_HE=\{U>0\}.                   \tag{9.2}
\]

The Gaussian is reflection invariant, so Theorem 2.1 applies. Here the
gain is explicit:

\[
 J_\gamma(E)=\mathbb E\min\{|U|,|V|\}
 =\frac{2(\sqrt2-1)}{\sqrt\pi},\qquad
 J_\gamma(P_HE)=\mathbb E|U|=\sqrt{\frac2\pi}.           \tag{9.3}
\]

Thus the exact fixed-measure monotonicity can be strict and can turn a
nonlinear tilted cut into a halfspace when a matching reflection symmetry
is present.

### 9.2 Isotropic cube halfspace: equality

For \(\mu\) uniform on \([-\sqrt3,\sqrt3]^n\) and
\(E=\{x_1>0\}\),

\[
                         J_\mu(E)=\mathbb E|X_1|
                                  =\frac{\sqrt3}{2}.      \tag{9.4}
\]

Coordinate polarization either fixes \(E\) or reflects it, so equality
holds in (0.1). Coordinate polarizations of a tilted halfspace only sort
the signs of its normal coefficients; cube symmetry again leaves the
value unchanged. This model supplies no strict improvement to iterate.

### 9.3 Isotropic radial exponential median sphere: equality

Let

\[
 d\mu_n(x)=c_n e^{-\sqrt{n+1}|x|}\,dx.                   \tag{9.5}
\]

Then \(R=|X|\) has the gamma law with shape \(n\) and rate
\(\sqrt{n+1}\), and \(\operatorname{Cov}(\mu_n)=I\). If \(r_n\) is the
median of \(R\) and \(E=\{|x|<r_n\}\), then

\[
                         J_{\mu_n}(E)=\mathbb E|R-r_n|.   \tag{9.6}
\]

Every reflection through an origin hyperplane fixes both \(\mu_n\) and
\(E\), so all allowed polarizations give equality. This is the radial
sphere obstruction in its exact isotropic normalization.

### 9.4 One-dimensional asymmetric log-affine test

Propositions 5.1 and 6.1 supply the two distinct one-dimensional audits:
two-point joint sorting loses exact balance by \(1/4\), while admissible
Steiner symmetrization followed by rebalancing and whitening has both
possible strict signs. The smooth potentials (4.4) and (6.12) show that
neither phenomenon depends on nonsmooth support.

### 9.5 Two-dimensional simplex

Section 8 gives exact mass, signed-distance, covariance, and whitening
constants. The normalized loss \(1/\sqrt3\) is dimension-free but in the
wrong direction for a monotone reduction.

## 10. Iterative compactness and why it does not restore the route

Theorem 2.1 is valid at every finite stage. Passing to infinitely many
polarizations requires more than weak or \(L^1\) convergence of indicator
functions, because \(J\) depends on the topological boundary. For example,
under uniform measure on \([0,1]\),

\[
 E=(-\infty,1/2),\qquad
 E'=(-\infty,1/2)\cup(\mathbb Q\cap(1/2,\infty))          \tag{10.1}
\]

are equal almost everywhere, but

\[
 J(E)=\frac14,\qquad J(E')=\frac18,                       \tag{10.2}
\]

because every point of \([1/2,1]\) belongs to \(\partial E'\). Open
approximations obtained by adding a fine family of intervals of total
length tending to zero exhibit the same discontinuity in the limit.

There is a more suitable compactness device for a **fixed** tight measure.
For balanced open sets \(E_j\), choose signed-distance functions

\[
 s_j=d(\,\cdot\,,E_j^c)-d(\,\cdot\,,E_j).                \tag{10.3}
\]

They are 1-Lipschitz. Balance and tightness anchor a zero of \(s_j\) in a
common bounded set: choose a ball carrying more than half the mass; it
meets both signs, and a segment between two such points meets the zero set.
Arzela--Ascoli gives a locally uniform subsequential limit \(s\). A finite
first moment gives uniform integrability and hence

\[
                         \int|s_j|\,d\mu\longrightarrow
                         \int|s|\,d\mu.                  \tag{10.4}
\]

The limit may have a positive-mass zero plateau and need not be the signed
distance of the open set \(\{s>0\}\). A generic linear tilt recovers the
set formulation approximately, but that tilt need not preserve any of the
symmetries produced by the iteration. Thus compactness can retain the
functional value, but it does not manufacture the missing canonical
halfspace/radial structure.

For density Steiner symmetrization there is an additional incompatibility.
After a step, whitening changes Euclidean reflections to affine
involutions; the next Euclidean Steiner step is not conjugate to the
previous one. Iterating without whitening loses isotropy, while iterating
with whitening incurs the noncommuting affine products discussed in
(7.9). Neither weak compactness nor the slicing bound controls those
products direction by direction.

## 11. Formal conclusion

The exact reusable lemma from this route is Theorem 2.1:

\[
 \text{reflection symmetry of \(\mu\)}
 \quad\Longrightarrow\quad
 \text{set polarization preserves mass and increases \(J\)}.     \tag{11.1}
\]

It gives a legitimate improvement operation for an already symmetric
measure, and the Gaussian cone calculation shows that it can be powerful.
It cannot be promoted to a reduction for all log-concave measures:

\[
 \begin{array}{c|c}
 \text{proposed extension}&\text{exact obstruction}\\ \hline
 \text{polarize density and set}&\text{balance changes: }1/2\to3/4\\
 \text{two-point polarize density}&\text{log-concavity can fail smoothly}\\
 \text{Steiner symmetrize, rebalance, whiten}&
       \text{both normalized monotonicity signs occur}\\
 \text{iterate symmetry polarizations}&
       \text{radial-section/sign-orbit invariants survive}\\
 \text{re-isotropize after every step}&
       \text{uncontrolled noncommuting affine products.}
 \end{array}                                               \tag{11.2}
\]

Accordingly, no dimension-free reduction to unconditional, radial, or
halfspace structure follows from polarization alone. Any successful use
of (11.1) in the KLS problem must obtain new reflection symmetries without
changing the balanced signed-distance functional, or must couple the two
interior-core gains in (1.4) to a separate global mechanism that controls
the covariance and the rebalancing error.
