# Global compatibility of calibrated rays

This note isolates three pieces of information which are genuinely global and
which are absent from the labelled-cylinder countermodel:

1. a sharp positive-reach inequality between any two two-sided calibrated
   rays of one Lipschitz potential;
2. the singular turning cost carried by a connected polyhedral fan, even
   though its classical shape operator vanishes almost everywhere; and
3. a midpoint-overlap constraint forced by global log-concavity across many
   ray bundles.

It also records the exact deterministic dichotomy behind the set-mass
preserving stochastic-localization control.  None of these facts alone proves
the desired dimension-free estimate, but together they specify more sharply
what a closing compatibility theorem has to control.

## 1. The two-sided core has positive reach

Let `f:R^n -> R` be 1-Lipschitz.  For `r>0`, define the oriented two-sided
`r`-core

\[
 \mathcal R_r(f)=\{(z,u): f(z)=0,\ |u|=1,\
                   f(z+tu)=t\text{ for every }|t|\le r\}.
\]

No differentiability of `f` or regularity of its zero set is assumed.

**Theorem 1 (sharp global reach inequality).**  If `(z,u)` and `(z',u')`
belong to `R_r(f)`, then, with `d=z-z'`,

\[
 \boxed{|d|^2-r^2|u-u'|^2
        \ge 2r\,|\langle d,u+u'\rangle|.}                 \tag{1}
\]

In particular,

\[
 |z-z'|\ge r|u-u'|.                                      \tag{2}
\]

Consequently the oriented direction is unique at each base point of the
two-sided core and the Gauss map `N(z)=u` is `1/r`-Lipschitz there.  Moreover,
for `|t|<r`, the normal map

\[
 F_t(z)=z+tN(z)
\]

is injective and satisfies

\[
 (1-|t|/r)|z-z'|
 \le |F_t(z)-F_t(z')|
 \le (1+|t|/r)|z-z'|.                                    \tag{3}
\]

**Proof.**  Calibration and the global Lipschitz inequality give the two
cross-ray inequalities

\[
 |(z+ru)-(z'-ru')|\ge 2r,
 \qquad |(z'+ru')-(z-ru)|\ge 2r.                         \tag{4}
\]

Writing `s=u+u'`, squaring (4), and using
`|s|^2=4-|u-u'|^2` gives

\[
 |d|^2-r^2|u-u'|^2+2r\langle d,s\rangle\ge0,
\]

\[
 |d|^2-r^2|u-u'|^2-2r\langle d,s\rangle\ge0.
\]

These are exactly (1), and (2) follows by discarding the nonnegative
right-hand side.  If `z=z'`, (2) gives `u=u'`.  Finally,

\[
 |d+t(u-u')|\ge |d|-|t||u-u'|
                \ge(1-|t|/r)|d|,
\]

and the reverse triangle inequality in the other direction gives the upper
bound in (3).  The lower bound proves injectivity for `|t|<r`.  QED.

At every smooth point, (2) implies `||D N||_op<=1/r`, but Theorem 1 is
strictly stronger than that pointwise statement: it remains valid between
different smooth charts and across a nonsmooth zero set.  Thus a family of
long rays with substantially different directions must either have base
points separated on the same scale or place the direction change outside the
two-sided core, at a focal, medial, or singular endpoint set.

There is also a precise, albeit unweighted, normal-variation consequence.

**Lemma 1.1 (direction rank forces full turning complexity).**  Let `nu` be a
probability on the unit sphere satisfying

\[
 \int uu^T\,d\nu(u)\preceq M I.                           \tag{4a}
\]

If a connected rectifiable set `K` in the unit sphere contains
`supp(nu)`, then

\[
 \mathcal H^1(K)\ge \frac1{8M}-\frac12.                  \tag{4b}
\]

Accordingly, any connected completed ray graph whose Gauss map contains this
direction support has total absolute normal variation at least the right-hand
side of (4b).

**Proof.**  If `|u-v|<=1`, then `u dot v>=1/2`, so (4a) gives

\[
 \nu\{u:|u-v|\le1\}\le4M.                                \tag{4c}
\]

Choose a maximal `1`-separated subset `v_1,...,v_k` of `supp(nu)`.  Its unit
balls cover the support, hence `k>=1/(4M)`.  The open half-unit balls about
the `v_i` are disjoint.  If `k>=2`, connectedness forces `K` to join each
`v_i` to the complement of its half-unit ball, contributing at least `1/2`
of one-dimensional Hausdorff measure inside that ball.  Summing over all but,
conservatively, one ball gives

\[
 \mathcal H^1(K)\ge(k-1)/2\ge1/(8M)-1/2.
\]

For a rectifiable graph map, total absolute variation dominates the Hausdorff
measure of its image, proving the final assertion.  QED.

For long rays of scale `s` occupying quotient mass `alpha`, covariance gives
`M<=(alpha s^2)^(-1)`; hence a connected normal completion must carry
`Omega(alpha s^2)` total turning.  This is the clean geometric uncertainty
principle suggested by the local focal estimate.  Its unresolved weighted
version is crucial: one must show that global log-concavity and joint
extremality prevent all of this connecting graph from being hidden in regions
of negligible density.  If the completion has many connected components,
Lemma 2 below supplies a different, additive bridge constraint.

## 2. A connected strictly log-concave fan with singular turning cost

The next example disproves any assertion that connected support, strict
log-concavity, exact balance, and one global Kantorovich potential force the
orientation changes to occur through a positive-measure region of smooth
curvature.

Fix an integer `m>=2`, put `alpha=pi/m`, and in `R^2` let

\[
 Z=\bigcup_{k=0}^{m-1}L_k,
 \qquad L_k=\mathbb R(\cos(k\alpha),\sin(k\alpha)).        \tag{5}
\]

The complement consists of `2m` sectors of angle `alpha`.  Give consecutive
sectors alternating signs and define `f_m` to be the signed Euclidean distance
to `Z` with those signs.

### 2.1 Global Lipschitzness and eikonal structure

The function `f_m` is globally 1-Lipschitz.  For two points in sectors of the
same sign this follows from the 1-Lipschitz property of distance to `Z`.  For
points of opposite sign, their joining segment crosses `Z`; if `z` is a
crossing point, then

\[
 |f_m(x)-f_m(y)|=d(x,Z)+d(y,Z)
 \le |x-z|+|z-y|=|x-y|.                                  \tag{6}
\]

It has `|grad f_m|=1` away from the zero rays and the angular bisectors, which
form a Lebesgue-null set.

Rotation by `alpha` preserves the standard Gaussian `gamma_2` and changes the
sign of `f_m`.  Hence

\[
 \gamma_2(f_m>0)=\gamma_2(f_m<0)=1/2,
 \qquad \int f_m\,d\gamma_2=0.                            \tag{7}
\]

### 2.2 Exact balanced-ray disintegration

Consider a zero ray and a base point `z=s tau` on it, where `s>0` and `tau`
is the radial unit vector.  Let `N` be the unit normal pointing into the
adjacent positive sector.  The maximal calibrated normal segment is

\[
 z+tN,\qquad |t|<T_s,
 \qquad T_s=s\tan(\alpha/2).                              \tag{8}
\]

Indeed, its endpoints are the intersections with the two neighboring angular
bisectors.  On the segment `f_m(z+tN)=t`.  The normal chart is flat and has
Jacobian one.  Since `z` is orthogonal to `N`, the Gaussian conditional density
is

\[
 q_s(t)=\frac{e^{-t^2/2}\mathbf 1_{(-T_s,T_s)}(t)}
              {\int_{-T_s}^{T_s}e^{-v^2/2}dv}.            \tag{9}
\]

It is exactly symmetric, and therefore every active ray assigns mass `1/2`
to each sign.  Reflection `t -> -t` preserves (9), interchanges the two signs,
and satisfies

\[
 f_m(z+tN)-f_m(z-tN)=2t=|(z+tN)-(z-tN)|.                 \tag{10}
\]

Integrating these conditional reflections gives a coupling between the two
normalized sign restrictions.  Equality (10) proves that `f_m` is an optimal
Kantorovich potential for that cut.

All classical shape operators on the open zero rays vanish.  The conditional
standard deviations nevertheless satisfy

\[
 \sigma_s\le \min\{1,T_s\}.                              \tag{11}
\]

Thus increasing the number of orientations shortens the available rays rather
than creating a large-scale witness.

### 2.3 The missing curvature is a singular turning charge

Let `N_j` be the oriented normal on the `j`-th zero ray, pointing into its
neighboring positive sector.  The sign alternation gives

\[
 \angle(N_j,N_{j+1})=\pi-\alpha,
 \qquad |N_j-N_{j+1}|=2\cos(\alpha/2).                   \tag{12}
\]

Each of the `m` positive sector components has two boundary rays meeting at
the origin.  Complete each such boundary through its vertex and sum absolute
turning, rather than allowing the vector first variations of different
components to cancel.  The resulting full focal/turning measure has chordal
normal variation at least

\[
 2m\cos(\pi/(2m)),                                       \tag{13}
\]

and the corresponding angular turning measure is

\[
 m(\pi-\pi/m)=\pi(m-1).                                  \tag{14}
\]

The Gaussian weight at the vertex is `(2pi)^(-1)`, so this weighted absolute
singular charge is still of order `m`.  This is a measure on the completed
ray-boundary graph; it is deliberately stronger than the classical shape
operator on the reduced boundary, which does not see the vertex.  Equations
(8) and (13) exhibit the uncertainty
relation in this model:

\[
 \text{typical two-sided reach}\asymp m^{-1},
 \qquad \text{singular normal variation}\asymp m.        \tag{15}
\]

This cost is invisible to `|II|` almost everywhere.  By contrast, the
disconnected-cylinder construction has no junction inside the support; a
McShane extension may connect the charts only through a region where the
measure has zero density, so there is no weighted singular turning charge.

Finally, the fan potential is not a first-moment extremizer when `m` is large.
If `G` is standard Gaussian in `R^2`, then

\[
 |f_m(G)|\le |G|\frac{\pi}{2m},
 \qquad
 \mathbb E|f_m(G)|\le\frac{\pi}{2m}\sqrt{\frac\pi2}.     \tag{16}
\]

A linear Gaussian coordinate has first absolute moment `sqrt(2/pi)`.  Thus
for `m>=3` the fan is quantitatively far from the global `D_1` maximizer.  The
example proves that being a Kantorovich potential for one's own cut is not a
substitute for joint extremality over all cuts and potentials.

## 3. Global log-concavity forces multibundle bridge overlap

The following elementary lemma is a useful exact way to encode the gluing
which the disconnected mixture lacks.

**Lemma 2 (midpoint-overlap inequality).**  Let `mu` be log-concave and let
`B_1,...,B_m,C_1,...,C_m` be compact sets such that

\[
 \mu(B_i)\ge\beta/m,\qquad \mu(C_j)\ge\beta/m             \tag{17}
\]

for all `i,j`.  Put

\[
 M_{ij}=(B_i+C_j)/2,
 \qquad K(x)=\sum_{i,j=1}^m\mathbf1_{M_{ij}}(x).           \tag{18}
\]

Then

\[
 \boxed{\int K\,d\mu\ge\beta m.}                        \tag{19}
\]

In particular, if `K<=K_0` almost everywhere, then `m<=K_0/beta`.

**Proof.**  Log-concavity gives

\[
 \mu(M_{ij})\ge\sqrt{\mu(B_i)\mu(C_j)}\ge\beta/m.
\]

Sum this inequality over all `m^2` pairs and apply Tonelli.  QED.

The compactness assumption is only to match the set formulation of
log-concavity.  Inner approximation gives the same conclusion, with an
arbitrarily small loss, for Borel sets.

Here is its consequence for hypothetical long-ray data.  Suppose a quotient
set `Omega` of mass at least `alpha` consists of rays with conditional scale
at least `s`, and suppose the covariance decomposition gives

\[
 \int_\Omega \sigma_y^2 N_yN_y^T\,d\eta(y)\preceq I.      \tag{20}
\]

For the normalized direction law `nu=N_#(eta|Omega)/eta(Omega)`,

\[
 \int uu^T\,d\nu(u)\preceq\frac1{\alpha s^2}I.            \tag{21}
\]

Thus every atom of `nu` has mass at most `(alpha s^2)^(-1)`, and every
Euclidean unit cap `{u:|u-v|<=1}` has mass at most `4/(alpha s^2)`.  Therefore
the long rays necessarily occupy order `alpha s^2` direction packets (atoms
may be grouped and the nonatomic part may be split so that packet masses are
comparable).

If, on every such ray, positive and negative bands at distance comparable to
`s` have conditional mass at least `kappa`, divide the quotient into
`m=floor(c alpha s^2)` packets of comparable mass and let `B_i,C_i` be the
unions of those two bands.  Then

\[
 \mu(B_i),\mu(C_i)\ge c\kappa\alpha/m.                    \tag{22}
\]

Lemma 2 says that the cross-bundle midpoint sets have average overlap at least
`c kappa alpha m`, hence of order `s^2`.  This is a genuine global constraint:
the disconnected-cylinder model violates exactly the midpoint lower bounds,
while a globally log-concave counterexample must realize very high additive
overlap among different orientation bundles.

What is not yet proved is the needed inverse statement: high bridge overlap,
together with Theorem 1, should force either a common cylindrical factor or a
large full-BV/focal charge.  A universal bounded-overlap assertion would close
the argument immediately, but it is false without additional geometric input;
high overlap can arise from approximate reflection or additive structure.

## 4. Set-mass-preserving stochastic localization

Let `p_t` be an Eldan localization with

\[
 dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle,
 \qquad Q_T=\int_0^T C_t^2dt.                             \tag{23}
\]

For a fixed Borel set `S`, put

\[
 g_t=\mu_t(S),\qquad
 v_t=\operatorname{Cov}_{\mu_t}(\mathbf1_S,X).
\]

When `v_t!=0`, set `e_t=v_t/|v_t|` and

\[
 C_t=P_{e_t^\perp}=I-e_te_t^T.                            \tag{24}
\]

When `v_t=0`, take `C_t=I` until a direction emerges.  Put
`R_t=e_te_t^T` when `v_t!=0` and `R_t=0` when `v_t=0`.  All identities below
are first applied to a bounded stopped process.

Since

\[
 dg_t=v_t^TC_t\,dW_t=0,                                  \tag{25}
\]

the set mass is pathwise constant.  Also

\[
 Q_T=T I-\int_0^T R_tdt.                                  \tag{26}
\]

The following spectral alternative is deterministic.

**Lemma 3 (full curvature or one pathwise line).**  For `0<epsilon<1`, either

\[
 \lambda_{\min}(Q_T)\ge\epsilon T,                        \tag{27}
\]

or there is a unit vector `u`, depending on the localization path, such that

\[
 \int_0^T\|R_t-P_u\|_F^2dt<2\epsilon T.                   \tag{28}
\]

**Proof.**  Put `R=int_0^T R_tdt`.  If (27) fails, then
`lambda_max(R)>(1-epsilon)T`.  For a maximizing unit eigenvector `u`,

\[
 \int_0^T u^TR_tu\,dt>(1-\epsilon)T.
\]

Since `R_t` is either zero or a rank-one orthogonal projector,

\[
 \|R_t-P_u\|_F^2
 =\operatorname{tr}(R_t^2)+1-2u^TR_tu
 \le2(1-u^TR_tu).
\]

Integration proves (28).  QED.

On paths satisfying (27), `p_T` is `(epsilon T)`-strongly log-concave in the
ordinary Euclidean metric.  The strongly-log-concave Cheeger inequality and
(25) give

\[
 P_{\mu_T}(S)\ge c\sqrt{\epsilon T}\min(g_0,1-g_0).       \tag{29}
\]

Unlike covariance-normalized localization, there is no random metric in
(29).  The fixed-boundary transfer is exact:

\[
 \mathbb E P_{\mu_T}(S)=P_\mu(S),                         \tag{30}
\]

first for smooth finite-perimeter boundaries by Tonelli and
`E p_T(x)=p(x)`, and then by BV approximation.

Thus only the near-one-line paths in (28) remain.  If their vector `u` were a
single deterministic direction, transverse localization would disintegrate
`mu` into parallel one-dimensional conditionals, every one retaining mass
`g_0`; the covariance identity would then give

\[
 \int\sigma_y^2d\eta(y)\le\operatorname{Var}\langle X,u\rangle=1,
\]

and the one-dimensional log-concave Cheeger bound would close the argument.
The difficulty is precisely that `u` in Lemma 3 is path-dependent.

This path dependence cannot be removed by symmetry or connectedness alone.
For the Gaussian fan of Section 2, the initial set covariance is zero.  Indeed,
rotation by `alpha` changes the sign set to its complement, while its action on
vectors has no `-1` eigenvector when `m>=2`; hence

\[
 \operatorname{Cov}_{\gamma_2}(\mathbf1_{\{f_m>0\}},X)=0. \tag{31}
\]

The control initially has no preferred line, and its law is equivariant under
the fan rotations.  Any pathwise residual direction selected later can
therefore have a nontrivial rotation-invariant mixture law.  The fan is not a
bad-scale example because all its one-dimensional variances are bounded by
one, but it rules out the inference "near rank-one accumulated precision
implies one deterministic direction" without a scale-sensitive global
compatibility estimate.

## 5. Exact remaining statement

Both the ray route and the mass-preserving localization route reduce to the
same missing phenomenon.  A sufficient theorem would say that a globally
log-concave isotropic density cannot decompose, through one jointly extremal
cut, into a positive-mass collection of balanced one-dimensional conditionals
of scale `s>>1` whose directions vary by path or bundle.  Theorem 1 forces the
variation either into base separation or into the focal complement of the
`s`-core.  The fan shows that the latter must be measured by full BV/turning
curvature, including singular junctions.  Lemma 2 shows that log-concavity
forces extensive midpoint overlap across the separated bundles.

A closing result must therefore prove one of the following quantitative
alternatives with universal constants:

1. the full weighted BV/focal charge created outside the long core yields a
   competing cut or potential with strictly larger first-moment objective;
2. the high midpoint overlap has an inverse theorem giving a common
   cylindrical factor, after which covariance and tensorization close the
   estimate; or
3. in the localization formulation, every path-dependent residual direction
   of large conditional variance incurs a universal amount of accumulated
   curvature before it is selected.

The Gaussian fan satisfies all local balance and calibration hypotheses but
pays the singular charge and has scale at most one.  The disconnected-cylinder
model has large scale but evades both the charge and the bridge overlap only by
failing global log-concavity.  This is the sharp distinction established in
this note.
