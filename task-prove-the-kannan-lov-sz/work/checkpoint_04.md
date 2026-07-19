# Checkpoint 04: endpoint radial rank, Gaussian-profile rigidity, and adaptive no-go tests

## 1. Exact results added since Checkpoint 03

### 1.1 Projection-wise endpoint rank and the far-endpoint branch

For a boundary-weighted family `F` of balanced rays with conditional scale
at least `s`, put

\[
 P_F=\int_Fq_y(0)d\eta(y),\qquad
 d\nu={q_y(0)\over P_F}d\eta.
\]

For every orthogonal projection `P` of rank `d`, covariance and the
dimension-free thin-shell theorem for the marginal give

\[
 P_Fs^3E_\nu|PN|^2\le Cd,
 \qquad
 P_Fs^5E_\nu|PN|^4\le Cd.                         \tag{1.1}
\]

Consequently the span of the active normals has rank

\[
 r_N\ge cP_Fs^5.                                      \tag{1.2}
\]

On a single variance band of quotient mass `alpha`, this becomes
`r_N>=c alpha s^4`.  Ambient padding is therefore irrelevant.

Under the additional, explicitly marked hypothesis that the ideal
two-endpoint heat-bath gap is at least `1-epsilon`, the product endpoint
slack

\[
 G(b,c)=|b-c|^2-(r(b)+ell(c))^2
\]

obeys

\[
 E_{\mu_B\otimes\mu_C}G
 \le {2\epsilon\over1-\epsilon}E_\nu L^2,             \tag{1.3}
\]

and

\[
 E_\nu L^2\ge c(1-\epsilon)P_Fs^5.                   \tag{1.4}
\]

If both endpoint densities are a fixed fraction of the central ray density,
one-dimensional log-concavity gives `L=O(s)` and forces `s=O(1)` on a
fixed-mass band.  Thus any remaining ideal-gap survivor must have high
normal rank, RMS chord length at least order `s^2`, and exponentially weak
endpoint density or large half-length variation.

The ideal heat-bath hypothesis cannot be inferred from the translation
charge.  The geometric medial form has conductances
`rho/|N_i-N_j|`, whereas a conditional-variance heat bath requires the
rank-one factorization `w_ij=c pi_i pi_j`.  The translation trace fixes only
the normal-jump-weighted total charge and has no canonical Markov time
normalization.

### 1.2 A genuine smooth Clifford-cone survivor

There is an exact globally log-concave isotropic model in
`R^m direct-sum R^m` with density

\[
 \exp[-\kappa_m(|u|+|v|)^2]
\]

and signed distance

\[
 f(u,v)={|u|-|v|\over\sqrt2}.                           \tag{1.5}
\]

On a ray with radial parameter `R`,

\[
 q_R(t)={ (1-t^2/R^2)^{m-1}\over R B(1/2,m)},
 \qquad
 \sigma_R^2={R^2\over2m+1}.                            \tag{1.6}
\]

Both focal endpoint densities vanish, while

\[
 \sigma_R^2\|S_R\|_{HS}^2={2(m-1)\over2m+1}<1.        \tag{1.7}
\]

The genuine smooth normal-height form, not an idealized endpoint form, has
the exact product coefficients

\[
 {\mathcal S(h)\over P}
 =E_\nu\left[{R^2\over2m}|\partial_Rh|^2
 +{1\over m-1}(|\nabla_\xi h|^2+|\nabla_\zeta h|^2)\right], \tag{1.8}
\]

and satisfies the full centered stability inequality.  Therefore the claim
“positive endpoint-family weight forces active endpoints at distance
`O(s)`” is false.

The escape is quantitatively sharp: on the family `sigma_R>=s`,

\[
 P_F\le C\exp(-c m s^2).                               \tag{1.9}
\]

It does not refute a fixed-bulk-mass theorem.  The exact surviving geometric
target is coherence of the curvature centers on such a fixed-mass band.

A pointwise spectral lemma now shows that if

\[
 c\le\sigma^2\operatorname {tr}S^2\le C,
 \quad \operatorname {rank}S\le Ak,
 \quad \sigma^2\|S\|_{op}^2\le A/k,
\]

then at least `ck` principal focal radii are comparable to
`sigma sqrt(k)`.  Codazzi gives exact relative-nullity and splitting-tensor
identities, and a repeated-curvature block has a leafwise constant focal
center.  What is not proved is that the nullity splitting tensor vanishes or
that the leafwise center maps and core projections are coherent across the
fixed-mass family.  Under that additional coherence, projected Paouris and
the covariance decomposition immediately force `s=O(1)`.

### 1.3 Gaussian-profile localization and a new equality decomposition

For ordinary stochastic localization of a set, let

\[
 g_t=\mu_t(S),\qquad v_t=\operatorname {Cov}_t(1_S,X),
 \qquad \mathcal I(a)=\varphi(\Phi^{-1}(a)).
\]

Strong log-concavity gives the sharp posterior centroid bound

\[
 \sqrt t|v_t|\le\mathcal I(g_t),                       \tag{1.10}
\]

and Ito calculus gives

\[
 d(\sqrt t\,\mathcal I(g_t))
 =\sqrt t\,\mathcal I'(g_t)v_t\cdot dW_t
 +{\mathcal I(g_t)^2-t|v_t|^2
   \over2\sqrt t\,\mathcal I(g_t)}dt.                 \tag{1.11}
\]

Thus `F(t)=sqrt(t) E I(g_t)` is nondecreasing, bounded above by the
original perimeter, and converges to that perimeter.  The exponent `1/2`
is sharp.  The identity starts from `F(0)=0`, so it supplies no seed.

The sharp centroid deficit has now been decomposed exactly.  In the active
one-dimensional marginal, let `q(y)=P(S|Y=y)`, let `H` be the upper
`g`-quantile threshold, let `T` transport a standard Gaussian to `Y`, and
write `T(z)=z/sqrt(t)+R(z)` with `R` nonincreasing.  Then

\[
 {\mathcal I(g)\over\sqrt t}-|v|
 =\int(y-c)(H-q)d\rho
 -E[R(Z)1_{Z\ge z_g}],                                 \tag{1.12}
\]

and both terms on the right are nonnegative.  On central posterior states,
small deficit gives an explicit threshold-error bound away from the
interface and an explicit near-maximal-slope bound for `T` across central
quantile blocks.  Equality forces a Gaussian active marginal and a halfspace
cut.  The profile drift (1.11) is precisely the time integral of these two
defects, up to a factor in `[1,2]`.

The live Gaussian route is therefore a quantitative gluing problem: control
the transverse part of

\[
 D_t=E_t[(1_S-g_t)(X-a_t)(X-a_t)^T]                   \tag{1.13}
\]

and the rotation of `v_t/|v_t|` from the scalar defects in (1.12), using
overlap of nearby posterior tilts.  A trace estimate is insufficient because
near equality is an adaptive top-variance event.

### 1.4 Exact limits of low-SNR and two-phase localization

For the Gaussian channel, the first two expansions of binary chi-square
information and mutual information are now exact.  With
`a=E[hX]` and `M=E[hXX^T]`,

\[
 \chi(t)=t|a|^2+t^2\left({1\over2}\|M\|_{HS}^2-|a|^2\right)+o(t^2), \tag{1.14}
\]

and

\[
 I(t)={t\over2}|a|^2+t^2\left({1\over4}\|M\|_{HS}^2
 -{1\over2}|a|^2+{1\over4}|a|^4\right)+o(t^2).        \tag{1.15}
\]

Indicator structure gives `|a|<=1`, `||M||op<=1`, while projection thin
shell gives only `||M||HS^2<=C(1+log n)`.  This is infinitesimal and has no
uniform finite-time remainder.  For centered exponential input the
likelihood ratio is outside `L^2(gamma)` for every positive channel time and
the origin density has zero Taylor radius.  Hermite summation from
`Psi_1` moment estimates is therefore invalid.

For the informative-then-preserving control, Phase A has exact scalar SDEs

\[
 dg=|v|d\beta,
 \quad dv=De\,d\beta-A v\,dt,
 \quad dA=\mathcal T(e)d\beta-Aee^TA\,dt.              \tag{1.16}
\]

If `|v_0|>=eta`, dimension-free Hilbert-martingale estimates give, with
positive probability, the directional seed

\[
 e_0^TQ_Ae_0\ge c\eta^2.                               \tag{1.17}
\]

This is only a Rayleigh bound, not a Loewner rank-one bound.  At the switch,
a small eigenvalue of the total precision is exactly a line which Phase A
missed and Phase B subsequently protected for almost all of its duration.
Gaussian sign parity has `v_0=D_0=0` and all label-correlated polynomial
moments below degree `n` equal to zero.  Thus no dimension-independent
bounded-degree informative seed can handle arbitrary half-sets.

### 1.5 Canonical random-block descent fails when its factors are separated

For a split `R^n=U direct-sum V`, the two-block heat-bath form has exact gap

\[
 {1-\rho(X_U,X_V)\over2}.                              \tag{1.18}
\]

For the uniform ball with block dimensions `k,l`, maximal correlation is

\[
 \rho=\sqrt{kl\over(k+2)(l+2)}.                        \tag{1.19}
\]

At a half split the gap is `2/(n+4)`, with the slow mode supplied by the
radial block constraint.  More decisively, for iid shifted exponentials in
dimension two, points `(t,t)` and directions within `O(1/t)` of the
anti-diagonal have conditional fiber variance `Omega(t^2)`.  Haar averaging
leaves an `Omega(t)` conditional-covariance weight; smooth localized tests
make its ratio to Euclidean Dirichlet energy unbounded.  Differentiating a
conditional expectation also produces a conditional score and boundary flux,
not merely the conditional mean of the projected gradient.

Therefore a recurrence obtained by separately lower-bounding the heat-bath
gap and upper-bounding conditional covariance weights is false.  Any block
route must compare its slow modes and energy jointly and must remove the
radial modes using more than a scalar gap estimate.

## 2. Newly blocked families

1. **Literal active-endpoint coercivity.**  The Clifford cone has zero focal
   endpoint densities and exact smooth stability.  Reopen only with a fixed
   bulk/boundary fraction and a mechanism globalizing curvature centers.
2. **Translation charge implies endpoint heat bath.**  Conductance
   normalization is not fixed by the normal-jump trace.
3. **Finite-order Gaussian-channel summation.**  Exponential tails destroy
   `L^2(gamma)` analyticity at every positive time.
4. **Universal bounded-degree informative seed.**  Gaussian sign parity
   annihilates every bounded-degree correlated moment.
5. **Separated random-block recurrence.**  Both the half-block gap and the
   averaged conditional-covariance weight fail in explicit isotropic
   log-concave models.
6. **Convolution residual stability.**  A proposed all-test-function reverse
   bound forces a Gaussian Stein identity and fails for an explicit smooth
   uniformly convex one-dimensional perturbation of a Gaussian.  Only the
   corrected affine-Stein defect identity survives.

## 3. Active incompatible mechanisms

1. **Fixed-bulk focal coherence.**  Use Codazzi/relative-nullity and the
   exact `sigma sqrt(k)` focal scale to prove that a fixed-mass stable ray
   band has one coherent core projection or a finite improving competitor.
   The leafwise-to-global center step is the load-bearing lemma.
2. **Gaussian-centroid gluing.**  Convert the two scalar near-equality
   defects into transverse splitting, small angular quadratic variation, and
   overlap-compatible posterior halfspaces.  A successful occupation bound
   would seed the exact Gaussian profile formula.
3. **Non-Gaussian convex-cell localization.**  Audit random slab or
   simultaneous-bisection posteriors, which preserve log-concavity without
   Gaussian extreme tilts.  The tests must distinguish a real covariance
   gain from classical needle localization with one long exceptional
   direction.

## 4. Updated extremal constraint set

A hypothetical sequence with T3 scale `s->infinity` must now satisfy all of
the following.

1. A fixed quotient-mass variance band has `sigma_y comparable s`.
2. On that band the active normal span has rank at least `c alpha s^4`, and
   coherent direction caps have exponentially small mass.
3. The completed endpoint geometry either uses exponentially weak endpoint
   densities/large half-length variation or avoids the ideal heat-bath
   normalization entirely through a genuine smooth curvature form.
4. Every rank-balanced smooth piece has order `sigma sqrt(k)` focal radii,
   but the leafwise focal centers and core projections must fail to cohere;
   otherwise projected covariance/Paouris bounds force `s=O(1)`.
5. Under ordinary localization, any small-profile-drift central state must be
   close to a posterior halfspace with an almost maximally Gaussian active
   marginal.  A counterexample must rotate or relabel these halfspaces across
   posterior states without paying the transverse `D_t` energy.
6. It cannot obtain this rotation from any fixed collection of bounded-degree
   label moments, and it cannot be certified by a separated block gap plus a
   conditional covariance estimate.

## 5. Audit status

There is still no candidate complete proof.  Every displayed implication in
this checkpoint has been kept on the proved side of its stated hypotheses;
in particular, ideal endpoint gap, coherent focal centers, transverse
posterior splitting, and covariance shrinkage under convex-cell localization
are not conclusions.  Symbolic dimension tracking therefore stops at these
three explicit load-bearing statements rather than concealing them in an
`O(1)` term.
