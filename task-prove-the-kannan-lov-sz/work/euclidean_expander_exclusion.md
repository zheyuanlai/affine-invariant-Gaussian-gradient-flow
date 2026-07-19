# Euclidean expander gluing: midpoint bridges, calibration slack, and the remaining obstruction

## Executive conclusion

No complete Euclidean exclusion of the long-ray expander branch is proved
here, and no counterexample satisfying the full hypotheses is constructed.
The reason for this distinction is substantive.  A positive quotient-mass
family of balanced rays with conditional scale tending to infinity for a
global `T3` extremizer is already a counterexample to the desired
dimension-free first-moment bound.  Conversely, excluding precisely that
configuration is the load-bearing KLS step.

What can be proved, without KLS, concentration, or a Poincare inequality, is
the following.

1. Cross-packet chords which are still in the interior of their calibrated
   rays pay an explicit additive slack of order
   `s |N-N'|^2`.  Exact cross calibration in the long core is necessarily
   the original same-ray matching.
2. The transport proof of Brunn--Minkowski gives an actual, null-invariant
   midpoint law, not merely membership in a Minkowski sum.  For `M^2`
   packet pairs, a midpoint loses at least one of the two `log M` packet
   labels.  This is sharp at the entropy scale `log M ~ s` and still does
   not select the calibrated matching.
3. Simultaneously barycentering many packet laws produces an exact
   log-determinant gain.  It is linear in `s` for ideal dispersed tubes, but
   long conditional variance does not force the required synchronized
   Brenier differentials.  An explicit common-source example has dispersed
   scale-`s` target variances and only `Theta(log s)` determinant gain.
4. A regular simplex under a Gaussian law realizes a complete,
   equal-conductance medial graph with high-rank normals for one global
   signed-distance function.  Thus line incidence, global log-concavity,
   and complete medial expansion alone do not give the desired obstruction.
   The model fails exactly at long conditional scale and raywise balance.
5. Any correct inverse theorem needs a third harmless branch besides
   parallel and concurrent families: constant cross distances give spheres
   in orthogonal affine subspaces (the Clifford branch).  Interior
   calibration rules this branch out in the core, so it can arise only in
   the endpoint/focal gluing which remains uncontrolled.

The new results therefore identify a precise missing lemma: endpoint
conductance must force either a strong-convexity-like bridge gain, synchronized
multi-map determinant diversity, or proximity to a parallel, concurrent, or
orthogonal-spherical model.  None of those implications follows from the
present scalar facts.

## 1. Setup and why a full countermodel would be a KLS countermodel

Let `f` be one-Lipschitz and let its active rays be written

\[
        x=z_y+tN_y,\qquad f(x)=t,\qquad |N_y|=1.
\]

Write `q_y(t)dt` for the conditional law and suppose the exact balance is

\[
 q_y(( -\infty,0])=q\in[\delta,1-\delta],\qquad p=1-q.
                                                               \tag{1.1}
\]

The following elementary observation fixes the logical strength of the
requested construction.  A one-dimensional log-concave density of standard
deviation `sigma` obeys

\[
                    \|q_y\|_\infty\le {1\over\sigma}.           \tag{1.2}
\]

Let `Q_y(r)` be its quantile.  Set

\[
 t_-=Q_y(q/2),\qquad t_+=Q_y((1+q)/2).
\]

The interval between them has probability `1/2`, so (1.2) gives

\[
                         t_+-t_-\ge\sigma_y/2.                  \tag{1.3}
\]

The positive quantile band between probabilities `(1+q)/2` and
`(3+q)/4` has mass `p/4`, and the negative band between `q/4` and
`q/2` has mass `q/4`.  Since all values in the former are at least
`t_+` and all values in the latter are at most `t_-`, the difference
`d_y` of the two sign-conditional means satisfies

\[
                              d_y\ge\sigma_y/8.                  \tag{1.4}
\]

For a global extremizer the exact cut--transport formula is

\[
             {\cal D}(\mu)=2pq\int d_y\,d\eta(y).                \tag{1.5}
\]

Consequently, if a quotient set `Omega` has mass `alpha` and
`sigma_y>=s` on it, then

\[
             {\cal D}(\mu)\ge {\delta^2\alpha s\over4}.          \tag{1.6}
\]

Thus an actual model with fixed `alpha,delta>0` and `s->infinity` is not a
local countertest; it is an unbounded first-moment example.  This report does
not claim to construct one.

## 2. Interior cross-calibration forces the same-ray matching

The next lemma is purely Euclidean and is the sharp useful statement about
same-ray matching.

**Lemma 2.1 (two-ended interior calibration).**  Let `f` be one-Lipschitz,
let `x,y in R^n`, `u,v in S^{n-1}`, and put

\[
                         L=f(x)-f(y)>0.
\]

Assume that, for some `h>0`, both endpoints continue away from the chord
inside calibrated rays:

\[
 f(x+tu)=f(x)+t,\qquad f(y-tv)=f(y)-t,qquad 0\le t\le h.          \tag{2.1}
\]

Then

\[
 \boxed{\quad
 |x-y|-L\ge {hL\over8(L+h)}|u-v|^2.
 \quad}                                                          \tag{2.2}
\]

In particular, equality `|x-y|=L` forces `u=v=(x-y)/|x-y|`; the two
pieces belong to one maximal calibrated ray.

**Proof.**  Put `d=|x-y|>=L` and `w=(x-y)/d`.  Comparing `x+hu` with
`y` gives

\[
 (L+h)^2\le |dw+hu|^2,
\]

and hence

\[
 1-u\cdot w\le{(d-L)(d+L+2h)\over2hd}.                         \tag{2.3}
\]

The comparison of `x` with `y-hv` gives the same bound with `v`.  Since

\[
 |u-v|^2\le2|u-w|^2+2|v-w|^2,
\]

we obtain

\[
 |u-v|^2\le {4(d-L)(d+L+2h)\over hd}.                           \tag{2.4}
\]

Finally,

\[
 {d\over d+L+2h}\ge {L\over2(L+h)},
\]

because the difference after cross multiplication is
`(d-L)(L+2h)`.  Rearranging (2.4) proves (2.2).  In the equality case,
(2.3) forces `u=v=w`; adjoining the chord to the two continuations proves
the last assertion.  QED.

There is a uniform positive-mass version on long balanced rays.

**Corollary 2.2 (quantile bands have linear cross slack).**  Suppose
`sigma_y>=s` and (1.1) hold.  Define the negative and positive bands

\[
 C_y=[Q_y(q/4),Q_y(q/2)],\qquad
 B_y=[Q_y((1+q)/2),Q_y((3+q)/4)].                       \tag{2.5}
\]

Their conditional masses are respectively `q/4` and `p/4`, hence at least
`delta/4`.  On at least half of any given quotient family, one can choose
one of two types so that every cross choice `b in B_y`, `c in C_{y'}` of
the same type obeys

\[
 L=f(b)-f(c)\ge s/4,                                    \tag{2.6}
\]

and both endpoints have a calibrated outward continuation of length

\[
 h=\delta s/4.                                          \tag{2.7}
\]

Consequently

\[
 \boxed{\quad
 |b-c|-\bigl(f(b)-f(c)\bigr)
 \ge {\delta s\over64}|N_y-N_{y'}|^2.
 \quad}                                                  \tag{2.8}
\]

**Proof.**  From (1.3), either `t_+>=sigma_y/4` or
`-t_- >=sigma_y/4`; one of the two alternatives holds on at least half
the quotient mass.  In the first type every `b` in (2.5) is at least
`s/4`; in the second every `c` is at most `-s/4`.  This proves (2.6).
There is tail mass at least `delta/4` beyond each selected band.  Equation
(1.2) therefore forces the one-dimensional support to continue by at least
`delta sigma_y/4>=delta s/4`.  Take exactly that much continuation.  Since
`h<=L`, Lemma 2.1 gives a coefficient at least `h/16`, proving (2.8).
QED.

Thus separated packets cannot be joined by alternative calibrated chords in
the long core.  All expander-like switching must occur after one loses the
outward continuation, namely at medial or focal endpoints.

The midpoint consequence is also explicit.  With `m=(b+c)/2`, the
parallelogram identity and (2.8) give

\[
\boxed{\quad
 |m|^2\le {|b|^2+|c|^2\over2}-{L^2\over4}
 -{\delta s^2\over512}|N_y-N_{y'}|^2.
\quad}                                                   \tag{2.9}
\]

Indeed, `|b-c|^2-L^2>=2L(|b-c|-L)` and `L>=s/4`.
The extra squared-radius deficit is only order `s^2`; after conversion to
radius at the isotropic scale `sqrt(n)`, it is order `s^2/sqrt(n)`.  This is
exactly too small to give a dimension-free contradiction from radial
second-moment information.

## 3. Actual midpoints from Brenier transport

Indicator membership in `(A+B)/2` is null-set unstable.  The next
proposition supplies actual endpoint representations and is the appropriate
replacement.

**Proposition 3.1 (transported midpoint domination).**  Let

\[
                   d\mu=Z^{-1}e^{-V(x)}dx=\rho(x)dx
\]

be log-concave, let `A,B` have masses `a,b>0`, and let `T` be the Brenier
map from `mu_A` to `mu_B`.  If

\[
             \nu=\left({\operatorname{Id}+T\over2}\right)_\#\mu_A,
\]

then

\[
 \boxed{\qquad {d\nu\over d\mu}\le {1\over\sqrt{ab}}.\qquad}   \tag{3.1}
\]

At `z=(x+T(x))/2` the sharper smooth pointwise identity is

\[
 {d\nu\over d\mu}(z)
 = {e^{-\Delta_V(x,T(x))}\over\sqrt{ab}\,Q(DT(x))},             \tag{3.2}
\]

where

\[
 \Delta_V={V(x)+V(Tx)\over2}-V(z)\ge0,
 \qquad
 Q(A)={\det((I+A)/2)\over\sqrt{\det A}}\ge1.                     \tag{3.3}
\]

**Proof.**  Monge--Ampere and midpoint change of variables give

\[
 \det DT(x)={b\rho(x)\over a\rho(Tx)},\qquad
 {d\nu\over dx}(z)
 ={\rho(x)\over a\det((I+DT(x))/2)}.                            \tag{3.4}
\]

Substitution of the first identity into the second, followed by division by
`rho(z)`, proves (3.2).  Convexity of `V` gives `Delta_V>=0`, and the
scalar arithmetic--geometric mean inequality applied to the eigenvalues of
the positive matrix `DT` gives `Q>=1`.  Standard regularization and
Alexandrov differentiation give (3.1) without smoothness.  QED.

**Corollary 3.2 (midpoint label entropy).**  Let `A_i,C_j`,
`1<=i,j<=M`, all have mass at least `beta/M`.  Choose the midpoint law
`nu_ij` from Proposition 3.1, choose `(I,J)` uniformly, and conditionally
sample `Z` from `nu_IJ`.  Then

\[
 D({\cal L}(Z)\Vert\mu)+I((I,J);Z)\le\log(M/\beta),              \tag{3.5}
\]

and hence

\[
                         H(I,J\mid Z)\ge\log(\beta M).           \tag{3.6}
\]

**Proof.**  Proposition 3.1 gives
`D(nu_ij||mu)<=log(M/beta)`.  Average this identity and use the exact
chain rule

\[
 {1\over M^2}\sum_{i,j}D(\nu_{ij}\Vert\mu)
 =D({\cal L}(Z)\Vert\mu)+I((I,J);Z).                             \tag{3.7}
\]

Subtract the mutual information from `H(I,J)=2log M`.  QED.

This says that a genuine transported midpoint is compatible, on average,
with at least `beta M` packet pairs.  It is a measure-theoretic theorem and
survives all null modifications.  Nevertheless it does not select the
same-ray matching:

* the map in Proposition 3.1 is `W_2`-optimal, whereas the calibrated ray
  plan is `W_1`-optimal;
* (2.8) makes separated cross chords longer, but mere convexity of `V`
  attaches no penalty to their length;
* translations through an affine part of `V` have `Delta_V=0`, `DT=I`,
  and `Q=1` at arbitrary separation; and
* for `M=exp(Theta(s))`, the allowance `log M` in (3.5) is already of
  order `s`.

For comparison, if `V` is `lambda`-strongly convex and
`dist(A,B)>=r`, then `Delta_V>=lambda r^2/8`.  Integrating (3.2) yields

\[
                 r^2\le {4\over\lambda}\log{1\over ab}.          \tag{3.8}
\]

Thus the cross-midpoint strategy closes immediately in the strongly
log-concave case.  Its exact missing input for a general convex potential is
a substitute for the vanished curvature term.

## 4. Simultaneous barycenters and the determinant audit

Pairwise determinants give only `log s` for orthogonal scale-`s` tubes.  A
simultaneous construction is stronger.

**Proposition 4.1 (multi-packet density and entropy inequality).**  Let
`A_i` have masses `a_i`, let `theta_i>0` sum to one, and take an absolutely
continuous common source `gamma` with Brenier maps

\[
                         (T_i)_\#\gamma=\mu_{A_i}.
\]

Put

\[
 S=\sum_i\theta_iT_i,\qquad \nu=S_\#\gamma,
\]

and, at points of Alexandrov differentiability,

\[
 \begin{split}
 \Delta_V&=\sum_i\theta_iV(T_i)-V(S),\\
 Q_\theta&={\det(\sum_i\theta_iDT_i)
       \over\prod_i\det(DT_i)^{\theta_i}}.
 \end{split}                                                    \tag{4.1}
\]

Then `Delta_V>=0`, `Q_theta>=1`, and, in the smooth injective setting,

\[
 {d\nu\over d\mu}(Sx)
 = {e^{-\Delta_V(x)}\over
 Q_\theta(x)\prod_i a_i^{\theta_i}}.                            \tag{4.2}
\]

Consequently

\[
 \boxed{\quad
 D(\nu\Vert\mu)+E_\gamma[\Delta_V+\log Q_\theta]
 =\sum_i\theta_i\log{1\over a_i}.
 \quad}                                                         \tag{4.3}
\]

**Proof.**  If `g` is the density of `gamma`, Monge--Ampere gives

\[
 g={\rho(T_i)\over a_i}\det DT_i,
 \qquad {d\nu\over dx}(Sx)={g\over\det DS}.                     \tag{4.4}
\]

Take the weighted geometric mean of the first identities.  Convexity gives

\[
 \rho(S)\ge e^{\Delta_V}\prod_i\rho(T_i)^{\theta_i},
\]

and concavity of `log det` gives `Q_theta>=1`.  Cancelling `g` proves
(4.2); integration proves (4.3).  For nonsmooth data, regularization and
lower semicontinuity retain the corresponding one-sided inequality, which is
all later bounds require.  QED.

The matrix scaling is exactly right in the ideal tube model.  Let

\[
 A_i=I+(s-1)u_i u_i^T,qquad
 C=\sum_i\theta_i u_i u_i^T.                                  \tag{4.5}
\]

For two orthogonal directions with equal weights, the gap is only

\[
 Q(A_i,A_j)={(s+1)^2\over4s}
 ={s+2+s^{-1}\over4},\qquad \log Q=\log s+O(1).                 \tag{4.5a}
\]

More generally, for two positive matrices `A,B`, if the eigenvalues of
`A^{-1/2}BA^{-1/2}` are `e^{t_k}`, then

\[
 {\det((A+B)/2)\over\sqrt{\det A\det B}}
       =\prod_k\cosh(t_k/2).                                   \tag{4.5b}
\]

Thus a pair sees only logarithmic distortion of the relative stretches; the
linear gain comes from averaging many dispersed rank-one stretches at once.

Then `det A_i=s` and

\[
                  \log Q_\theta
 =\log\det(I+(s-1)C)-\log s.                                  \tag{4.6}
\]

Since `Tr C=1`, if

\[
                         \|C\|_{op}\le {\kappa\over s^2},       \tag{4.7}
\]

then `Tr C^2<=kappa/s^2` and `log(1+x)>=x-x^2/2` give

\[
 \boxed{\quad
 \log Q_\theta\ge s-1-\log s-\kappa/2.
 \quad}                                                         \tag{4.8}
\]

For orthonormal `u_i=e_i`, equal weights, and `M>=s^2`, this specializes
to

\[
 \log Q=M\log(1+(s-1)/M)-\log s\ge s-2-\log s.                 \tag{4.9}
\]

The ray covariance constraint supplies precisely a bound of the form
(4.7) for the *normal directions*.  The unresolved step is to transfer it
to the Brenier derivatives in (4.5).  Conditional variance does not do so.

**Example 4.2 (desynchronized stretches).**  Let `gamma` be uniform on
`[0,1]^M`, let `I=[1/2-1/(2M),1/2+1/(2M)]`, and define the increasing map

\[
 \tau(t)=t+sM\,|[0,t]\cap I|.                                  \tag{4.10}
\]

For `1<=i<=M`, let `T_i` act as `tau` in coordinate `i` and as the identity
in every other coordinate.  Each `T_i` is a Brenier map.  Its target has
variance at least `(1-1/M)^2s^2/4` in direction `e_i`: the lower and upper
source blocks have masses `(1-1/M)/2` and their images are separated by at
least `s`.

At a source point set

\[
 K(x)=\#\{i:x_i\in I\}.
\]

For equal weights,

\[
 \det\left({1\over M}\sum_iDT_i\right)=(1+s)^{K(x)},
 \qquad
 \left(\prod_i\det DT_i\right)^{1/M}
 =(1+sM)^{K(x)/M}.                                           \tag{4.11}
\]

Since `E K=1`,

\[
 \boxed{\quad
 E_\gamma\log Q
 =\log(1+s)-{1\over M}\log(1+sM)
 =\Theta(\log s),
 \quad}                                                       \tag{4.12}
\]

even when `M>>s^2`.  Thus dispersed scale-`s` variances and positive
Brenier differentials do not imply the linear gain (4.8): the long stretches
may occur on disjoint source events.  The targets in this example have a
low-density bridge and are not log-concave; this is relevant rather than
incidental, because ray-packet restrictions and unions of their two endpoint
bands need not be log-concave either.  A successful use of (4.3) must prove
synchronization from the fact that all packets come from one global
log-concave density and one extremal ray congruence.

The determinant loss can in fact be made arbitrarily smaller than
`log s`, even with exact balance.

**Example 4.3 (balanced rare-stretch maps with vanishing determinant
gain).**  Let `G` be standard Gaussian in `R^M`, let `u_i=e_i`, and, for
`epsilon>0`, define the odd increasing function

\[
 \theta_\varepsilon(t)=
 \begin{cases}
 -1,&t\le-\varepsilon,\\
 t/\varepsilon,&|t|<\varepsilon,\\
 1,&t\ge\varepsilon.
 \end{cases}                                                   \tag{4.13}
\]

For `L>0`, put

\[
 T_i(x)=x+L\theta_\varepsilon(x_i)e_i
 =\nabla\left({|x|^2\over2}
       +L\int_0^{x_i}\theta_\varepsilon(r)dr\right).            \tag{4.14}
\]

These are Brenier maps.  Their target laws are exactly balanced in the long
coordinate, and

\[
 \operatorname{Var}\langle T_i(G),e_i\rangle
 \ge L^2\mathbb P(|G_i|\ge\varepsilon)
 \ge L^2(1-\sqrt{2/\pi}\,\varepsilon),                          \tag{4.15}
\]

while all transverse coordinates remain standard Gaussian.  On the other
hand,

\[
 DT_i=I+{L\over\varepsilon}{\bf1}_{\{|x_i|<\varepsilon\}}
             e_ie_i^T.                                         \tag{4.16}
\]

If `K(x)=#\{i:|x_i|<epsilon\}`, then the equal-weight gap is

\[
 \log Q(x)=K(x)\left[
 \log\left(1+{L\over M\varepsilon}\right)
 -{1\over M}\log\left(1+{L\over\varepsilon}\right)
 \right].                                                       \tag{4.17}
\]

Writing `p_epsilon=P(|G_1|<epsilon)<=sqrt(2/pi)epsilon` gives

\[
 0\le E\log Q
 \le \sqrt{2/\pi}\,M\varepsilon
       \log\left(1+{L\over M\varepsilon}\right)
 \longrightarrow0.                                             \tag{4.18}
\]

At the singular limit `epsilon=0`,
`T_i(x)=x+L sign(x_i)e_i` has `DT_i=I` almost everywhere and `Q=1`,
although its target has two balanced lobes separated at scale `L` in a
different orthogonal direction for every `i`.  Thus even balance, long
variance, and maximal direction dispersion do not force any determinant
gain for non-log-concave packet laws.  What must rule out this example is
not matrix algebra but the global geometry of packet restrictions and their
same-ray central bridges.

There is also no free gain from iterating pairwise midpoint operations along
an expander path.  Log-determinant Jensen gaps telescope under a binary
averaging tree:

\[
 \log\det\Big(\sum_i\theta_iA_i\Big)
 -\sum_i\theta_i\log\det A_i
\]

equals the root gap plus the appropriately weighted gaps inside its child
groups.  It is the single global gap (4.1), not an unweighted sum of every
pair encountered.  Reusing laws also reintroduces the corresponding entropy
costs.  Any path-iteration proof therefore needs an additional
synchronization or rank-growth invariant; graph expansion by itself does not
supply one.

## 5. Exact midpoint incidence and the third inverse branch

At the atomic level Euclidean space does rule out the abstract cyclic
many-center pattern.

**Lemma 5.1 (one full matching center).**  Let
`P={p_1,...,p_M}` and `C={c_1,...,c_M}` be sets of distinct points.  For a
fixed `z`, the bipartite graph

\[
                  G_z=\{(i,j):p_i+c_j=2z\}                       \tag{5.1}
\]

is a matching.  If it has `M` edges, then `C=2z-P`.  There cannot be two
distinct centers with `M` edges.

**Proof.**  A fixed `p_i` determines `c_j=2z-p_i`, and conversely, so no
vertex has degree two.  A matching with `M` edges is full and gives the
reflection identity.  If both `z,z'` were full centers, the finite set `P`
would be invariant under translation by `2(z-z')`.  Maximizing a linear
functional in that translation direction shows that a finite set cannot be
invariant under a nonzero translation.  QED.

This exact fact is not stable enough for positive-volume bands, and even its
strongest conclusion is reflection, not concurrence of the calibrated rays.
The deranged cyclic endpoint configuration gives an explicit failure of that
last implication.

Moreover, reflection/concurrence is not the only equality geometry which a
robust inverse theorem must allow.

**Lemma 5.2 (constant cross-distance classification).**  If nonempty
sets `B,C` satisfy

\[
                         |b-c|=d\quad(b\in B,c\in C),             \tag{5.2}
\]

then there are orthogonal subspaces `U,V`, affine centers `b_0,c_0`, and
`w=b_0-c_0` orthogonal to `U+V`, such that

\[
 B\subset b_0+r_BS(U),\qquad C\subset c_0+r_CS(V),
 \qquad r_B^2+r_C^2+|w|^2=d^2.                                 \tag{5.3}
\]

**Proof.**  Subtracting the four squared instances of (5.2) gives

\[
       \langle b-b',c-c'\rangle=0.                               \tag{5.4}
\]

Hence `U=span(B-B)` and `V=span(C-C)` are orthogonal.  Choose the closest
affine centers so `w` is orthogonal to both.  Writing `b=b_0+u` and
`c=c_0+v` turns (5.2) into

\[
                 |w|^2+|u|^2+|v|^2=d^2.
\]

Fixing first `u` and then `v` proves that both radii are constant.  QED.

For independent endpoint laws `B,C`, the robust algebraic diagnostic is the
exact ANOVA identity

\[
\begin{split}
 \operatorname{Var}|B-C|^2={}&
 \operatorname{Var}(|U|^2+2\langle a,U\rangle)
 +\operatorname{Var}(|V|^2-2\langle a,V\rangle)\\
 &+4\operatorname{Tr}(\Sigma_B\Sigma_C),                       \tag{5.5}
\end{split}
\]

where `a=EB-EC`, `U=B-EB`, and `V=C-EC`.  Thus nearly constant product
cross distances force nearly radial endpoint laws and nearly orthogonal
covariance ranges.  What is missing is a bound on the left side of (5.5):
neither (3.1) nor midpoint entropy controls the product law of the endpoints.

Lemma 2.1 clarifies where the Clifford model can occur.  If the common
cross distance also equals the `f`-value gap and both endpoints retain a
positive calibrated continuation, then all the normals agree and the model
collapses to the same-ray/parallel case.  A nondegenerate orthogonal-spherical
branch is therefore necessarily an endpoint or focal phenomenon, exactly
where the expander energy is currently hidden.

## 6. A globally log-concave complete medial graph

The following construction disproves any exclusion based only on Euclidean
facet incidence and log-concavity.

**Theorem 6.1 (Gaussian regular-simplex countertest).**  Let `M>=3` and

\[
 H=\{x\in\mathbb R^M:\sum_i x_i=0\},\qquad d=M-1,
\]

and choose regular-simplex unit vectors `u_1,...,u_M in H` satisfying

\[
 \langle u_i,u_j\rangle=-{1\over M-1}\quad(i\ne j),
 \qquad {1\over M}\sum_i u_iu_i^T={1\over M-1}I_H.              \tag{6.1}
\]

For `a>0`, put

\[
 K_a=\{x\in H:\langle u_i,x\rangle\le a\ \hbox{for all }i\},  \tag{6.2}
\]

and let `F` be signed Euclidean distance to `partial K_a`, positive inside.
Under the standard Gaussian on `H`:

1. In `K_a`,
   \[
             F(x)=\min_i\{a-\langle u_i,x\rangle\}.              \tag{6.3}
   \]
   Every pair of branches has a relatively open medial interface.  Hence
   the medial adjacency graph is `K_M`.
2. For facet-height variables `h_i`, the exact envelope junction form is
   \[
       {\cal E}_{med}(h)=w_{M,a}\sum_{i<j}(h_i-h_j)^2,            \tag{6.4}
   \]
   with the same `w_{M,a}>0` on every edge.  Its normalized graph gap is
   `M/(M-1)>=1`, and the normal law has rank `M-1` by (6.1).
3. Every open facet-normal ray has conditional variance at most one.
   If `z=a u_i+w`, `w perpendicular u_i`, and `N=-u_i`, then
   \[
      F(z+tN)=t\quad(-\infty<t<r_i(z)),\qquad
      q_z(t)\propto e^{-(a-t)^2/2}{\bf1}_{t<r_i(z)},              \tag{6.5}
   \]
   where
   \[
      r_i(z)=\min_{j\ne i}{M-1\over M}
                   (a-\langle u_j,z\rangle)\in(0,a].             \tag{6.6}
   \]
   Direct one-dimensional integration by parts for the truncated Gaussian
   gives `Var_qz(T)<=1`.
4. Even if `a` is chosen so the simplex has Gaussian mass `1/2`, exact
   raywise balance fails.  In fact
   \[
    q_z(T>0)=
    {\Phi(r_i(z)-a)-\Phi(-a)\over\Phi(r_i(z)-a)},                 \tag{6.7}
   \]
   which varies from zero as `r_i(z)` tends to zero to
   `1-2\Phi(-a)` at `r_i(z)=a`.

**Proof.**  Formula (6.3) is the distance from an interior point to the
union of the complementary facet halfspaces.  For every `i!=j`, the set

\[
 \Gamma_{ij}=\{x:\langle u_i,x\rangle=\langle u_j,x\rangle
       >\langle u_k,x\rangle\ (k\ne i,j)\}                       \tag{6.8}
\]

has a nonempty `(d-1)`-dimensional part in `K_a`; take
`x=t(u_i+u_j)` with
`0<t<a(M-1)/(M-2)`.  The two-branch envelope formula gives edge weight

\[
 w_{M,a}=left({2M\over M-1}\right)^{-1/2}
      \int_{\Gamma_{ij}}\varphi_H\,d\mathcal H^{d-1}>0,          \tag{6.9}
\]

because `|u_i-u_j|=sqrt(2M/(M-1))`.  Permutation symmetry makes all weights
equal, and the complete-graph identity gives the claimed gap.

Along the ray in (6.5),

\[
 |z-tu_i|^2=|w|^2+(a-t)^2,
\]

and the flat normal Jacobian is one.  Solving the first equality with a
competing branch gives (6.6).  The truncated Gaussian variance bound and
the elementary integral in (6.7) finish the proof.  QED.

This is a connected, globally log-concave, full-dimensional Euclidean model
with a single global signed-distance function, zero smooth facet curvature,
high-rank normals, and a complete equal-conductance medial graph.  It is not
a counterexample to the full proposed exclusion: its ray conditionals have
scale at most one and are not balanced.  It proves that any successful
obstruction must use the long--balanced--extremal conjunction, not merely
Euclidean realizability of an expander graph.

## 7. Exact remaining statement

The interior of the long core is now understood: by (2.8), alternative
cross-packet calibration cannot occur there.  The entire unresolved problem
is at the completed endpoint graph.

A sufficient new theorem would take a positive-mass family of long balanced
rays from one global extremizer and prove one of the following alternatives:

1. the focal/medial conductance has too little spectral gap, contradicting
   the extremal second variation;
2. transported packet laws have synchronized Brenier stretches, so (4.3)
   and (4.8) produce a linear-in-`s` determinant gain exceeding the packet
   entropy;
3. the endpoint product cross-distance variance in (5.5) is small, yielding
   an approximately orthogonal-spherical model; or
4. the rays are approximately parallel or concurrent.

The regular-simplex example rules out item 1 from incidence alone.
Example 4.2 rules out item 2 from conditional variances and direction
dispersion alone.  Proposition 3.1 gives actual midpoint representations but
does not control item 3 or identify the `W_1` same-ray matching.  Establishing
one of these implications requires genuinely new use of the fact that the
packet laws are restrictions of one globally log-concave density and arise
from one jointly extremal signed-distance congruence.

That is the Euclidean/log-concave realization obstruction still missing from
the global `T3` route.
