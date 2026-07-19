# Multiway midpoint incidence at the natural-tilt scale

## Executive verdict

There is a clean quantitative midpoint dichotomy after one adds a genuine
packet hypothesis.  Suppose that the natural-tilt law is \(\nu_t\), that

\[
 t={\tau\over \mathsf K},
\]

and that disjoint phase cells \(A_i\) satisfy

\[
 \nu_t(A_i)\ge {a\over M},\qquad
 A_i\subset B(c_i,R\sqrt t),\qquad
 |c_i-c_j|>4R\sqrt t\quad(i\ne j).                 \tag{0.1}
\]

Define the posterior-overlap weights

\[
 w_{ij}=\exp\left[-{(|c_i-c_j|+2R\sqrt t)^2\over8t}\right]
 =\exp\left[-{1\over8}
 \left(\sqrt{\mathsf K\over\tau}|c_i-c_j|+2R\right)^2\right]
                                                               \tag{0.2}
\]

and their normalized density

\[
 \rho={1\over M^2}\sum_{i,j}w_{ij}.                  \tag{0.3}
\]

For every measurable proposed high-defect separator \(H\), one of the
following alternatives holds.

1. \(\nu_t(H)\ge a\rho/4\).
2. On a set of \(\nu_t\)-mass at least \(a\rho/4\), a single midpoint
   \(x\) approximately reflects at least \(a\rho M/4\) packet centers in
   disjoint pairs:

   \[
    |c_i+c_j-2x|\le2R\sqrt t,
   \]

   and every retained pair obeys

   \[
    |c_i-c_j|+2R\sqrt t
    \le\sqrt{8t\log {4\over a\rho}}.                 \tag{0.4}
   \]

Thus a bounded-overlap branch really does force a dimension-free
high-defect charge when \(\rho\) is numerical.  The parameter \(\rho\),
not the unweighted number \(M^2\), is the correct effective pair density.
It retains exactly the exponential cost required to compare two posterior
states.

The conclusion cannot be extended to arbitrary cells, and even (0.4) does
not organize the phase normals.  Two decisive obstructions are proved
below.

* Minkowski midpoint multiplicity is unstable under null changes even when
  the cells are kept exactly disjoint.  Distinct null spheres can make all
  \(M^2\) midpoint sets contain one common positive-mass annulus.
* For the Gaussian natural-tilt law, the canonical separated winner cones

  \[
   A_i(h)=\{c:c_i\ge c_k+h\text{ for every }k\ne i\}
  \]

  have masses comparable to \(1/M\), orthogonal projective labels, and

  \[
   {A_i(h)+A_j(h)\over2}=\mathbb R^M\qquad(i\ne j).   \tag{0.5}
  \]

  Hence their midpoint overlap is maximal, while their phase lines are
  neither common affine lines nor concurrent radial lines.  The cells are
  convex and canonical, so this is not merely a null-representative
  pathology.

Accordingly, multiway Brunn--Minkowski incidence proves a useful
**bounded-packet reflection theorem**, but membership in Minkowski midpoint
sets cannot prove a common affine or radial organization of posterior phase
normals.  That last step needs a null-invariant measure on actual midpoint
representations and a theorem tying the representation to the same
posterior certificate.

## 1. Natural tilts and the unavoidable posterior-overlap cost

Let \(\mu\) be log-concave on its affine support and, on that support, write

\[
 d\pi_c(x)={1\over Z(c)}
 \exp\{c\cdot x-t|x|^2/2\}\,d\mu(x).                 \tag{1.1}
\]

The natural parameter

\[
 C_t=tX+\sqrt t\,G,
 \qquad X\sim\mu,\quad G\sim N(0,I)\text{ independent},       \tag{1.2}
\]

has a log-concave law \(\nu_t\).  Its density on the relevant affine
space is proportional to

\[
                         e^{-|c|^2/(2t)}Z(c).          \tag{1.3}
\]

No isoperimetric estimate for \(\nu_t\) will be used.

Put \(b(c)=\log Z(c)\).  Brascamp--Lieb for the
\(t\)-strongly log-concave posterior gives

\[
                         \nabla^2b(c)=\operatorname{Cov}_{\pi_c}(X)
                         \preceq t^{-1}I.             \tag{1.4}
\]

For \(c,d\) and \(m=(c+d)/2\), direct calculation gives

\[
 \begin{aligned}
 \mathbb E_{\pi_c}\left({d\pi_m\over d\pi_c}\right)^2
 &= {Z(c)Z(d)\over Z(m)^2}\\
 &\le \exp\left\{{|c-d|^2\over4t}\right\}.          \tag{1.5}
 \end{aligned}
\]

Indeed, the identity follows by integrating the square of the exponential
likelihood ratio, and the inequality is the centered second-difference
bound furnished by (1.4).  Consequently, for every measurable
\(F:\mathbb R^n\to[0,1]\),

\[
 \boxed{\quad
 \mathbb E_{\pi_m}F
 \le \exp\left\{{|c-d|^2\over8t}\right\}
       \sqrt{\mathbb E_{\pi_c}F}.
 \quad}                                                  \tag{1.6}
\]

Thus a certificate transported from an endpoint to its midpoint loses the
factor

\[
                         \exp\{-|c-d|^2/(8t)\}.        \tag{1.7}
\]

At \(t=\tau/\mathsf K\), the exponent is

\[
                         {|c-d|^2\over8t}
 ={\mathsf K|c-d|^2\over8\tau}.                       \tag{1.8}
\]

For packets as in (0.1), every \(c\in A_i,d\in A_j\) obeys

\[
 |c-d|\le |c_i-c_j|+2R\sqrt t.                       \tag{1.9}
\]

This is the origin of the safe, uniform weight (0.2).  Using the distance
between the two cells instead of the right side of (1.9) is not safe: a
point in a Minkowski midpoint set may be represented by endpoint pairs much
farther apart than the closest pair of cells.

## 2. The exact weighted multiplicity identity

For Borel cells \(A_1,\ldots,A_M\), set

\[
 S_{ij}={A_i+A_j\over2},\qquad
 N_w(x)=\sum_{i,j=1}^M w_{ij}\mathbf1_{S_{ij}}(x),     \tag{2.1}
\]

where \(w_{ij}\ge0\) are arbitrary deterministic weights.

**Lemma 2.1 (weighted midpoint incidence).**  If \(\nu\) is log-concave
and

\[
                         \nu(A_i)\ge {a\over M},       \tag{2.2}
\]

then

\[
 \boxed{\quad
 \int N_w\,d\nu
 =\sum_{i,j}w_{ij}\nu(S_{ij})
 \ge {a\over M}\sum_{i,j}w_{ij}.
 \quad}                                                  \tag{2.3}
\]

More generally, with \(a_i=\nu(A_i)\), the right side is
\(\sum_{i,j}w_{ij}\sqrt{a_i a_j}\).

**Proof.**  Log-concavity gives

\[
                         \nu(S_{ij})\ge\sqrt{a_i a_j}.
\]

Multiply by \(w_{ij}\), sum, and use Tonelli.  The statement for Borel
sets follows from inner approximation by compact sets, with an arbitrarily
small loss followed by a limit.  \(\square\)

For the unweighted multiplicity \(N\), (2.3) says

\[
                         \int N\,d\nu\ge aM.          \tag{2.4}
\]

This is only a first-moment statement.  Let \(H\) be a high-defect set and
\(G=H^c\).  If \(N_w\le L\) on \(G\), then

\[
 \boxed{\quad
 \sum_{i,j}w_{ij}\nu(S_{ij}\cap H)
 =\int_HN_w\,d\nu
 \ge {a\over M}\sum_{i,j}w_{ij}-L.
 \quad}                                                  \tag{2.5}
\]

This is a large **incidence-weighted** separator charge.  It need not be a
large ordinary mass.  Since the unweighted multiplicity can be \(M^2\),
(2.5) alone gives at best

\[
                         \nu(H)\gtrsim {1\over M}.     \tag{2.6}
\]

An upper bound of order \(M\) for the pointwise multiplicity is precisely
what converts (2.5) into a numerical separator mass.  The next section
identifies a transparent geometric hypothesis giving that bound.

## 3. A bounded-packet reflection dichotomy

The following result is the positive part of the midpoint program.

**Lemma 3.1 (midpoint graphs of separated packets are matchings).**
Suppose

\[
 A_i\subset B(c_i,r),\qquad |c_i-c_j|>4r\quad(i\ne j). \tag{3.1}
\]

For fixed \(x\), form a bipartite graph \(G_x\) on two copies of
\(\{1,\ldots,M\}\) by putting an edge \((i,j)\) when
\(x\in(A_i+A_j)/2\).  Then \(G_x\) is a matching.  In particular,

\[
                         N(x)\le M.                   \tag{3.2}
\]

For every edge,

\[
                         |c_i+c_j-2x|\le2r.           \tag{3.3}
\]

**Proof.**  If \((i,j)\) is an edge, take \(p\in A_i,q\in A_j\) with
\(p+q=2x\).  Then

\[
 |c_j-(2x-c_i)|
 \le |c_j-q|+|p-c_i|\le2r,                            \tag{3.4}
\]

which proves (3.3).  If both \((i,j)\) and \((i,k)\) were edges, (3.4)
would give \(|c_j-c_k|\le4r\), contrary to (3.1).  The same argument on
the other bipartite side excludes two edges ending at the same vertex.
Thus \(G_x\) is a matching.  \(\square\)

**Theorem 3.2 (overlap-weighted separator or common reflection).**
Assume (0.1), use the weights (0.2), and let \(\rho\) be (0.3).  For every
measurable \(H\), at least one of the following holds:

\[
                         \nu_t(H)\ge {a\rho\over4};    \tag{3.5}
\]

or there is a measurable set \(B\subset H^c\) satisfying

\[
                         \nu_t(B)\ge {a\rho\over4},    \tag{3.6}
\]

such that for every \(x\in B\) there are subsets \(I_x,J_x\) and a
bijection \(\sigma_x:I_x\to J_x\) with

\[
 |I_x|=|J_x|\ge {a\rho M\over4},                      \tag{3.7}
\]

\[
 |c_i+c_{\sigma_x(i)}-2x|\le2R\sqrt t,                \tag{3.8}
\]

and

\[
 |c_i-c_{\sigma_x(i)}|+2R\sqrt t
 \le\sqrt{8t\log {4\over a\rho}}.                   \tag{3.9}
\]

**Proof.**  Lemma 2.1 gives

\[
                         \int N_w\,d\nu_t\ge a\rho M.\tag{3.10}
\]

Lemma 3.1 and \(w_{ij}\le1\) give \(N_w\le M\).  Suppose (3.5) fails and
put

\[
 B=\{x\in H^c:N_w(x)\ge a\rho M/2\}.                 \tag{3.11}
\]

Then

\[
 \begin{aligned}
 a\rho M
 &\le\int N_w\,d\nu_t\\
 &\le M\nu_t(H)+{a\rho M\over2}+M\nu_t(B),
 \end{aligned}
\]

so \(\nu_t(B)>a\rho/4\), proving (3.6).

At a fixed \(x\in B\), the incident edges form a matching and have total
weight at least \(a\rho M/2\).  Edges of weight less than
\(a\rho/4\) contribute less than \(a\rho M/4\).  Therefore at least
\(a\rho M/4\) matching edges have weight at least \(a\rho/4\).  Restrict
to them.  Lemma 3.1 gives (3.8), while the definition (0.2) gives (3.9).
\(\square\)

The theorem has an immediate defect-charge form.  If \(D(c)\ge e\) on
\(H\), then the first branch gives

\[
                         \int D\,d\nu_t\ge {ea\rho\over4}.     \tag{3.12}
\]

If the phase profile uses a factor \(I(g)\) and the states in question are
central, \(I(g)\ge i_\delta\), then the same conclusion holds with the
right side multiplied by \(i_\delta\).  This uses no expansion theorem.

There is also a direct bounded-overlap corollary.  If
\(N_w\le L\) on \(H^c\), then (3.10) and \(N_w\le M\) give

\[
                         \nu_t(H)\ge a\rho-{L\over M}. \tag{3.13}
\]

Thus a genuinely bounded low-defect overlap, \(L=O(1)\), forces a
numerical separator mass whenever \(\rho\) is numerical and \(M\) is large.

The theorem also identifies its own limitation.  It reflects **packet
centers**, not phase normals.  Brunn--Minkowski and the multiplicity
identity are invariant under an arbitrary reassignment of projective labels
to the cells.  Hence no statement organizing those labels can follow from
(2.3) without an additional analytic relation between the actual midpoint
representation and the posterior phase certificate.

## 4. Two decisive obstructions outside the packet regime

### 4.1 Exact disjointness does not repair null instability

The usual null-set objection can be made while preserving literal
disjointness.

**Lemma 4.1 (disjoint null-sphere inflation).**  Let \(n\ge2\) and let
\(\nu\) be an absolutely continuous probability.  Suppose there are two
disjoint open balls \(U,V\) of positive \(\nu\)-mass and that \(\nu\) has
positive density almost everywhere on \(V\).  Given \(M\) disjoint compact
positive-mass sets \(B_i\subset U\), there are exactly disjoint compact
sets \(A_i\) such that

\[
                         \nu(A_i)=\nu(B_i)             \tag{4.1}
\]

and the \(M^2\) midpoint sets \((A_i+A_j)/2\) all contain one common
positive-\(\nu\)-mass annulus.

**Proof.**  Choose \(z,R,\varepsilon\) so that the closed ball containing
all the spheres below lies in \(V\).  Choose distinct radii

\[
                         R_i\in[R,R+\varepsilon]
\]

and put \(H_i=z+R_iS^{n-1}\).  Arrange that the \(B_i\) avoid all
\(H_j\), and set \(A_i=B_i\cup H_i\).  The sets are exactly disjoint and,
by absolute continuity, (4.1) holds.  For every \(i,j\),

\[
 {H_i+H_j\over2}
 =\left\{z+y:{|R_i-R_j|\over2}\le|y|
                    \le {R_i+R_j\over2}\right\}.     \tag{4.2}
\]

The equality follows by varying the angle between two vectors of lengths
\(R_i,R_j\).  Every set in (4.2) contains

\[
                         \{z+y:\varepsilon/2\le|y|\le R\},    \tag{4.3}
\]

which has positive \(\nu\)-mass.  \(\square\)

Thus neither exact disjointness nor fixed cell masses makes the indicator
multiplicity a measure-theoretic functional.  Any load-bearing inverse
theorem must use canonical representatives or an actual representation
measure.

### 4.2 Canonical convex winner cells have maximal overlap

Null modifications are not the only obstruction.  Let

\[
 C_i=\{x\in\mathbb R^M:x_i\ge x_k\text{ for all }k\},
 \qquad
 A_i(h)=h e_i+C_i.                                    \tag{4.4}
\]

Equivalently,

\[
                         A_i(h)=\{x:x_i-x_k\ge h\ (k\ne i)\}.
                                                               \tag{4.5}
\]

For \(h>0\), these are pairwise disjoint closed convex cells.

**Lemma 4.2 (winner-cone sum).**  If \(i\ne j\), then

\[
                         C_i+C_j=\mathbb R^M,
 \qquad
                         A_i(h)+A_j(h)=\mathbb R^M.    \tag{4.6}
\]

Moreover,

\[
                         \operatorname{dist}(A_i(h),A_j(h))
                         =\sqrt2\,h.                  \tag{4.7}
\]

**Proof.**  Fix \(z\in\mathbb R^M\).  Set \(x_j=0\),

\[
 x_k=z_k-z_j\quad(k\notin\{i,j\}),
\]

and choose

\[
 x_i\ge\max\{0,z_i-z_j,\max_{k\notin\{i,j\}}(z_k-z_j)\}.
\]

Then \(x\in C_i\), and \(y=z-x\) satisfies
\(y_j\ge y_k\) for all \(k\), so \(y\in C_j\).  This proves the first
identity in (4.6); the shifted identity follows because
\(h(e_i+e_j)+\mathbb R^M=\mathbb R^M\).

If \(x\in A_i(h),y\in A_j(h)\), then

\[
 (x-y)_i-(x-y)_j=(x_i-x_j)+(y_j-y_i)\ge2h,
\]

so \(|x-y|\ge\sqrt2h\).  Equality is attained by taking \(x_i=h/2\),
all other coordinates of \(x\) equal to \(-h/2\), and the analogous point
with \(j\) as winner for \(y\).  \(\square\)

Now take a genuinely natural Gaussian tilt law.  Let

\[
                         \mu=N(0,\mathsf K I_M),
 \qquad t={\tau\over\mathsf K}.
\]

Then \(\mathsf K\) is the Poincare constant of \(\mu\), and (1.2) gives

\[
                         \nu_t=N(0,P I_M),
 \qquad P=t(1+\tau).                                  \tag{4.8}
\]

If \(\Delta_M\) is the gap between the largest and second-largest of
\(M\) independent \(N(0,P)\) variables, permutation symmetry gives

\[
                         \nu_t(A_i(h))
 ={1\over M}\mathbb P\{\Delta_M\ge h\}.              \tag{4.9}
\]

The elementary Gaussian gap estimate

\[
 \mathbb P\{\Delta_M\le s\sqrt P\}
 \le C s\sqrt{\log M},
 \qquad 0\le s\le {c\over\sqrt{\log M}},             \tag{4.10}
\]

follows from

\[
 \mathbb P\{\Delta_M>r\sqrt P\}
 =M\int_{\mathbb R}\varphi(x)\Phi(x-r)^{M-1}\,dx
                                                               \tag{4.11}
\]

after splitting where \(M\bar\Phi(x)\) is respectively below, within,
and above fixed numerical bounds and applying the two-sided Mills ratios.
In particular, for a sufficiently small universal \(\varepsilon>0\),

\[
 h={\varepsilon\sqrt P\over\sqrt{\log M}}
 \quad\Longrightarrow\quad
 {1\over2M}\le\nu_t(A_i(h))\le{1\over M}.            \tag{4.12}
\]

Assign the projective label \(P_i=e_ie_i^T\) to \(A_i(h)\).  Then

\[
                         \|P_i-P_j\|_{HS}=\sqrt2       \tag{4.13}
\]

for \(i\ne j\), while Lemma 4.2 gives

\[
 \sum_{i,j}\mathbf1_{(A_i(h)+A_j(h))/2}(x)
 \ge M(M-1)\qquad\text{for every }x.                 \tag{4.14}
\]

There is no point \(z_0\) for which the lines
\(c+\mathbb Re_i\), for all \(c\) in all the full-dimensional cells
\(A_i(h)\), pass through \(z_0\).  Nor do the labels have one common affine
direction.  Thus maximal additive overlap does not organize the phase
normals.

Even the nearest-cell posterior penalty remains visible:

\[
 \exp\left[-{\operatorname{dist}(A_i(h),A_j(h))^2\over8t}\right]
 =\exp\left[-{h^2\over4t}\right]
 =\exp\left[-{\varepsilon^2(1+\tau)\over4\log M}\right].       \tag{4.15}
\]

If one incorrectly used nearest-cell distance as though it were a uniform
representation bound, the resulting pseudo-density would be

\[
 \widetilde\rho_{\min}
 ={1\over M}+\left(1-{1\over M}\right)
   \exp\left[-{\varepsilon^2(1+\tau)\over4\log M}\right].
                                                               \tag{4.16}
\]

For \(\tau=\kappa\log M\), this tends to
\(e^{-\varepsilon^2\kappa/4}\).  By contrast, the actual packet density
\(\rho\) in Theorem 3.2 is not available: the cones have infinite diameter,
so the safe uniform endpoint-distance bound is infinite and the
corresponding off-diagonal retention is zero.  In this precise sense,
\(\widetilde\rho_{\min}\) is numerical while
\(\rho_{\mathrm{safe}}=0\).

Thus, for \(\tau=O(\log M)\), (4.15) is numerical but cannot be attached
to all the incidences in (4.14): the cells are unbounded, and a given
midpoint may use endpoints arbitrarily farther apart than the closest pair.
The safe uniform weight (0.2) is zero in the unbounded-packet limit.  This
is an explicit demonstration that set membership loses the posterior
overlap cost.

Even the elementary angular-Lipschitz constraint does not remove this
model.  On the union of the cells, the assignment
\(c\mapsto P_i=e_ie_i^T\) has Lipschitz constant at most \(h^{-1}\), by
(4.7) and (4.13).  Kirszbraun's theorem extends it to a
Hilbert--Schmidt-valued map on all of \(\mathbb R^M\) with the same
constant.  When \(\tau\asymp\log M\), (4.12) has \(h\asymp\sqrt t\), so
this is exactly the \(O(t^{-1/2})\) angular scale.  The extension need not
be rank one in the tie region, which is precisely where a high-defect
phase transition is allowed.  This observation is not a posterior
realization; it shows that cell masses, separation, midpoint incidence, and
the known Lipschitz scale still do not organize the labels.

The winner cones are the exact phase-cell geometry of maximum/facet models.
Lemma 4.2 is a countermodel to any inverse theorem whose assumptions use
only the natural log-concave tilt law, equal cell masses, separated
projective labels, and Minkowski midpoint multiplicity.  If “low defect” is
required to mean realization by the posterior centroid of one fixed cut,
that analytic realization is an additional hypothesis and must be used in
the proof; it is invisible to the midpoint count.

## 5. Mandatory model tests

### 5.1 Orthogonal bounded cells

Let

\[
                         c_i=L\sqrt t\,e_i,
 \qquad A_i\subset B(c_i,r\sqrt t).                   \tag{5.1}
\]

The enclosing midpoint balls have centers
\((c_i+c_j)/2\) and radius \(r\sqrt t\).  Distinct unordered pairs of
orthogonal centers are separated by at least \(L\sqrt t/\sqrt2\).  Hence,
if

\[
                         r<{L\over2\sqrt2},            \tag{5.2}
\]

the midpoint sets are disjoint apart from the symmetry
\(S_{ij}=S_{ji}\), and \(N\le2\).  Lemma 2.1 would then give

\[
                         aM\le\int N\,d\nu_t\le2,      \tag{5.3}
\]

which is impossible for \(M>2/a\).  Thus a log-concave law cannot assign
mass \(a/M\) to many genuinely small orthogonal packets.  It must create
midpoint overlap or destroy one of the packet hypotheses.  The posterior
cost between two centers is

\[
                         {|c_i-c_j|^2\over8t}={L^2\over4}.      \tag{5.4}
\]

This test is exactly captured by \(\rho\): it is useful when \(L=O(1)\)
and exponentially weak when \(L\gg1\).

### 5.2 Radial Gaussian shells

For \(\nu_t=N(0,P I_d)\), a radial phase has

\[
                         u(c)={c\over|c|}.
\]

Angular packets on a common shell have a genuine common center at the
origin.  Thus the intended high-overlap inverse branch is correct in this
model.  However, on the typical shell \(R^2\simeq dP=d\,t(1+\tau)\), two
directions separated by angle \(\theta\) pay

\[
 {|c-d|^2\over8t}
 \simeq {d(1+\tau)\theta^2\over8}.                   \tag{5.5}
\]

Opposite packets therefore pay \(\exp[-\Theta(d(1+\tau))]\).  Only angles
\(\theta=O((d(1+\tau))^{-1/2})\) have numerical posterior overlap.  This
is the same scale as spherical concentration.  Large set-level overlap of
opposite radial sectors is consequently not a usable posterior bridge.

### 5.3 Gaussian parity

For the Gaussian posterior

\[
                         \pi_c=N(c/t,t^{-1}I_2)
\]

and the parity cut \(S=\{x_1x_2\ge0\}\), the tilts
\(c=(a\sqrt t,0)\) and \(d=(0,a\sqrt t)\) have orthogonal active phases.
The normalized centroid defect is a constant multiple of

\[
                         \Phi(-a)\asymp {e^{-a^2/2}\over1+a}.  \tag{5.6}
\]

Their squared separation is \(|c-d|^2/t=2a^2\), so (1.6) loses

\[
                         e^{a^2/4}.                    \tag{5.7}
\]

The square root of (5.6) cancels (5.7), up to a polynomial factor.  Thus
the overlap estimate is sharp and supplies no small midpoint defect.  The
four parity arms have only \(M=4\), their midpoint multiplicity is bounded,
and the central state has zero centroid.  This is a high-defect separator,
but its Fisher/profile weight can be exponentially smaller than the arm
defect.  Multiway counting gives no amplification in this fixed-\(M\)
model.

### 5.4 Product-exponential maximum cells

The coordinate-winner geometry is exactly (4.4)--(4.6), independently of
whether the underlying one-dimensional factor is Gaussian or exponential.
Moving a distance \(h\) into a winner cone replaces \(C_i\) by
\(h e_i+C_i\), but every cross midpoint set remains all of \(\mathbb R^M\).
The coordinate phase normals stay orthogonal and the tie regions carry the
phase-switching defect.

For a gap \(h=L\sqrt t\), the nearest-pair posterior retention is

\[
                         e^{-h^2/(4t)}=e^{-L^2/4}.      \tag{5.8}
\]

This number must be retained.  If saturation against all \(M-1\)
competitors requires \(L^2\gtrsim\log M\), it consumes a power of \(M\);
whether any multiway gain survives depends on the exact constant.  The
unweighted identity (4.14) cannot answer that question, because its
representations have no uniform endpoint-distance bound.  Truncating the
cells restores a valid weight and places the problem under Theorem 3.2, but
destroys the all-space sum identity.

This model is therefore decisive: high Minkowski overlap may be maximal for
the very phase geometry one wanted to exclude, without yielding affine or
radial concurrence of its normals.

## 6. Extremality-dependent test: canonical coarea tubes

The preceding countermodels use only the tilt-space cells.  A near-Cheeger
extremizer supplies additional canonical geometry which must be audited
separately.

Let \(S\) have mass \(1/2\) and perimeter

\[
                         p\le{\psi_\mu\over2}+\epsilon,
\]

and let \(F:\mathbb R^n\to[0,1]\) be a Lipschitz smoothing with
\(\mathbb EF=1/2\).  Put

\[
 U={1\over2}\mathbb E|F-\mathbf1_S|,
 \qquad E_r=\{F>r\}.
\]

For the relaxed weighted relative \(BV\) perimeter, the exact coarea
cascade is

\[
 \mathcal D_{\mathrm{co}}(F)
 :=\int_0^1\left[
 P_\mu(E_r)-\psi_\mu\min\{\mu(E_r),1-\mu(E_r)\}\right]\,dr
 \le\epsilon+2\psi_\mu U.                            \tag{6.1}
\]

For posterior-resampling smoothing at \(s=\alpha\mathsf K\), the bad-scale
estimates give

\[
                         U\le C\sqrt\alpha.
\]

Consequently, if \(\epsilon\le Cp\sqrt\alpha\) and
\(p\asymp\psi_\mu\), then

\[
                         \mathcal D_{\mathrm{co}}(F)
                         \le Cp\sqrt\alpha.           \tag{6.2}
\]

Also,

\[
                         \int_0^1\mu(E_r\triangle S)\,dr=2U.   \tag{6.3}
\]

For almost every regular value \(r\), the reduced boundary of \(E_r\)
has the canonical normal

\[
                         \theta={\nabla F\over|\nabla F|}.
\]

One may therefore define canonical boundary packets by imposing Borel
conditions on \(F\), \(\theta\), centrality, and defect, and then taking
their perimeter-essential representatives.  This removes Lemma 4.1:
adjoining arbitrary null spheres no longer changes the chosen packet.
It does **not** remove Lemma 4.2.  Canonical winner cells are already
essentially closed convex cells, and regular level surfaces can have many
winner facets.

There is a further mismatch of spaces.  The packets in Theorem 3.2 live in
the natural-parameter variable \(c\), whereas the regular tubes above live
in the physical boundary variable \(x\).  The exact Jensen identity aligns
the posterior direction and \(\theta(x)\) only in the joint
boundary-flux law.  It does not give an injective map from tilt packets to
disjoint physical boundary patches.  Such a disintegration would be an
additional lemma.

### 6.1 Gaussian winner cones fail through the coarea deficit

For the Gaussian maximum cut

\[
 S_M=(-\infty,q_M]^M,\qquad \Phi(q_M)^M={1\over2},
\]

the relative Gaussian perimeter is

\[
 p_M=M\varphi(q_M)\Phi(q_M)^{M-1}
 ={M\varphi(q_M)\over2\Phi(q_M)}
 \asymp\sqrt{\log M}.                                 \tag{6.4}
\]

Gaussian halfspaces show

\[
                         {\psi_{\gamma_M}\over2}=I(1/2)
                         ={1\over\sqrt{2\pi}}.
\]

Therefore

\[
 p_M-{\psi_{\gamma_M}\over2}\ge c\sqrt{\log M}
 \ge c'p_M.                                           \tag{6.5}
\]

For the \(BV\) function \(F=\mathbf1_{S_M}\), every nontrivial superlevel
is \(S_M\), so

\[
                         \mathcal D_{\mathrm{co}}(F)
 =p_M-{\psi_{\gamma_M}\over2}\asymp p_M.              \tag{6.6}
\]

The same limit holds along any strict \(BV\) smoothing of
\(\mathbf1_{S_M}\).  Hence these canonical Gaussian winner facets cannot
satisfy (6.2) as \(\alpha\downarrow0\).  This is the exact
extremality-dependent reason that the Gaussian winner-cone countermodel is
inadmissible for a near-Cheeger extremizer; midpoint incidence itself does
not exclude it.

This conclusion is not uniform in a simultaneous \(M\to\infty\),
\(\alpha\to0\) limit without a quantitative strict-\(BV\) convergence
rate.  The initial perimeter gap (6.5), however, already shows that the cut
does not meet the assumed near-minimizer hypothesis
\(\epsilon=O(p\sqrt\alpha)\) for small \(\alpha\).

### 6.2 The product-exponential winner family survives this test

For \(M\) independent rate-one one-sided exponentials, the nested inner
boxes

\[
                         B_q=[0,q]^M
\]

have mass \(m=(1-e^{-q})^M\) and relative perimeter

\[
 P_\mu(B_q)
 =M e^{-q}(1-e^{-q})^{M-1}
 =M(1-m^{1/M})m^{(M-1)/M}.                            \tag{6.7}
\]

At fixed \(m\in(0,1)\),

\[
                         P_\mu(B_q)\longrightarrow -m\log m.  \tag{6.8}
\]

In particular, the balanced maximum cut has perimeter tending to
\((\log2)/2\).  No known comparison in this argument gives a strictly
smaller balanced perimeter, and the exact Cheeger profile of this product
is not being assumed.  Therefore there is no proved lower bound

\[
 \mathcal D_{\mathrm{co}}(F)\ge c\,p
\]

for its canonical nested winner tubes.  Unlike the Gaussian case, the
coarea test does not exclude the product-exponential maximum family.

### 6.3 The coarea budget does not yet give a multiway charge

Write

\[
 d(r)=P_\mu(E_r)-\psi_\mu
                  \min\{\mu(E_r),1-\mu(E_r)\}\ge0.
\]

Equation (6.2) gives only

\[
 \left|\{r:d(r)>\lambda p\}\right|
 \le {C\sqrt\alpha\over\lambda}.                      \tag{6.9}
\]

It has no pair index and no cross term between normals on distinct boundary
patches.  In particular, \(M^2\) midpoint pairs may reuse one tie skeleton;
the perimeter and \(d(r)\) count that skeleton once, whereas
\(N_w\) counts it with multiplicity.

Combining Theorem 3.2 with coarea would require a new inequality of the
form

\[
 \mathcal D_{\mathrm{co}}(F)
 \ge c p\int_H d\nu_t
 \quad\hbox{or}\quad
 \mathcal D_{\mathrm{co}}(F)
 \ge {cp\over M}\int_HN_w\,d\nu_t,                   \tag{6.10}
\]

under a canonical coupling of tilt packets to regular boundary tubes.
Neither inequality follows from (6.1), flux-normal alignment, or
log-concavity.  If the first inequality in (6.10) were proved, Theorem 3.2
would give

\[
                         \mathcal D_{\mathrm{co}}(F)
 \ge {cpa\rho\over4}
\]

on the separator branch, contradicting (6.2) when
\(\sqrt\alpha\ll a\rho\).  This is only a conditional implication, not an
available estimate.

Thus regular nested level-set tubes cure the representative pathology and
exclude Gaussian winner cones by extremality, but they yield no proved
finite multiway separator charge.  The product-exponential winner family
remains the mandatory survivor.

## 7. What remains after the incidence theorem

The rigorous output is the following.

1. For arbitrary cells, Brunn--Minkowski gives only the weighted first
   moment (2.3).  A bounded low-defect multiplicity forces the
   incidence-weighted separator charge (2.5), not an ordinary mass.
2. For separated bounded packets, midpoint graphs are matchings.  Theorem
   3.2 then gives a dimension-free separator mass or a common approximate
   reflection matching, with the exact effective density \(\rho\) and with
   \(t=\tau/\mathsf K\) tracked.
3. Reflection of packet centers does not imply concurrence of the phase
   normals.  The incidence count is blind to the labels and to the
   same-posterior or same-ray matching.
4. Without bounded canonical packets, both null-sphere inflation and the
   canonical winner cones refute the desired inverse theorem.
5. Extremal coarea tubes canonically fix representatives and exclude the
   Gaussian winner family through a deficit of order \(p\), but the coarea
   budget has no incidence multiplicity.  It neither yields (6.10) nor
   excludes the product-exponential winner family.

A closing lemma would have to replace the indicator
\(\mathbf1_{x\in(A_i+A_j)/2}\) by a null-invariant density on actual
representations \((c,d,x)\), retain the factor
\(e^{-|c-d|^2/(8t)}\), and prove that a positive density of reflected
representations uses the posterior matching which defines the phase
normals.  It would then still need a multiway rigidity theorem showing that
the resulting pairwise matchings synchronize into one affine direction or
one radial/concurrent center.

Neither synchronization is a consequence of log-concavity or of the
midpoint multiplicity identity.  The bounded-packet theorem is therefore a
valid new incidence lemma, while the proposed unrestricted
midpoint-to-normal dichotomy is decisively false.
