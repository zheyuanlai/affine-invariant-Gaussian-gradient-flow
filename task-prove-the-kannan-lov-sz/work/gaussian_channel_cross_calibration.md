# Gaussian-channel overlap: exact normalization and the calibration gap

## 1. Setup

Let `X~mu`, let `B=1_S`, `p=mu(S)`, and observe

\[
 C=tX+\sqrt t\,G,
 \qquad
 q_t(c|x)=(2\pi t)^{-n/2}
       \exp\{-|c-tx|^2/(2t)\}.                       \tag{1.1}
\]

Write

\[
 \nu(dc)=p_t(c)dc,\qquad
 p_t(c)=\int q_t(c|x)d\mu(x),\qquad
 \ell_c(x)={q_t(c|x)\over p_t(c)}.                  \tag{1.2}
\]

Thus `dmu_c=ell_c dmu` is the localization posterior.  Put

\[
 g(c)=\mu_c(S),\quad v(c)=\operatorname {Cov}_{\mu_c}(B,X),
 \quad
 \mathsf B(c)={t\,v(c)v(c)^T\over g(c)(1-g(c))}.     \tag{1.3}
\]

At states where `v!=0`, write `u=v/|v|` and

\[
 r(c)=\operatorname {tr}\mathsf B(c)
      ={t|v(c)|^2\over g(c)(1-g(c))};                \tag{1.4}
\]

put `r=0` otherwise.  The binary Fisher matrix is

\[
 R=\int\mathsf B(c)d\nu(c)=\int r(c)u(c)u(c)^Td\nu(c),
 \qquad Z=\operatorname {tr}R.                       \tag{1.5}
\]

For an isotropic input, conditional covariance and total covariance give

\[
                         R\preceq tI.                \tag{1.6}
\]

The issue addressed here is whether the high effective rank forced by
`Z asymp 1` and `t=K^{-1}` can be converted into the positive-density
*cross-calibration* required in Sections 8--9 of
`finite_medial_competitor.md`.

## 2. Independent Fisher phases really are almost orthogonal

Define the normalized active-phase law

\[
                         d\rho(c)={r(c)\over Z}d\nu(c). \tag{2.1}
\]

If `C,C'` are independent with law `rho`, then

\[
 \boxed{
 \mathbb E_\rho\langle u(C),u(C')\rangle^2
 =\operatorname {tr}\left({R\over Z}\right)^2
 \le {t\over Z}.}                                    \tag{2.2}
\]

More generally, on a measurable good set `E`, put
`Z_E=int_E r dnu` and normalize there.  Since
`int_E r uu^T dnu preceq R preceq tI`, the right side of (2.2) is at most
`t/Z_E`.  This is the exact product-phase orthogonality statement.

## 3. A common Brenier source gives positive cross incidence

Assume the posteriors are `t`-strongly log-concave.  Let `T_c` be the
Brenier contraction from `N(0,t^{-1}I)` to `mu_c`, and, for `G~N(0,I)`,
write

\[
 X_c(G)=T_c(G/\sqrt t),\qquad
 A_c=\{G:X_c(G)\in S\}.                              \tag{3.1}
\]

Suppose `g(c) in [delta,1-delta]` and the relative centroid defect is

\[
 \epsilon_c=1-{sqrt t|v(c)|\over\mathcal I(g(c))}.  \tag{3.2}
\]

The estimates in Sections 2--4.1 of `spatial_phase_coherence.md` imply the
following source-space formulation.  With

\[
 \zeta_\delta(s)=C_\delta s\sqrt{\log(e/s)},\qquad
 \alpha_\delta(s)=C_\delta\{\sqrt s+\zeta_\delta(s)^{1/3}\}, \tag{3.3}
\]

there is `b_c=Phi^{-1}(1-g(c))` such that

\[
 \boxed{
 \gamma_n\bigl(A_c\mathbin\triangle
       \{\langle G,u(c)\rangle\ge b_c\}\bigr)
 \le\alpha_\delta(\epsilon_c).}                    \tag{3.4}
\]

Here is the short deduction.  The target cut differs from its active
halfspace by `O_delta(sqrt(epsilon_c))`.  If
`Y=<u,sqrt(t)(X_c-EX_c)>`, the Brenier derivative estimate and Gaussian
Poincare give

\[
       \mathbb E|Y-\langle G,u\rangle|^2\le\zeta_\delta(\epsilon_c).
\]

The elementary smoothing inequality

\[
 \sup_a|P(Y\le a)-P(\langle G,u\rangle\le a)|
 \le \inf_{h>0}\{\zeta/h^2+Ch\}\le C\zeta^{1/3}
\]

and the lower Gaussian density on central quantiles give (3.4).

If two good phases obey

\[
 |\langle u(c),u(c')\rangle|\le\delta^2/4,
 \qquad
 \alpha_\delta(\epsilon_c)+\alpha_\delta(\epsilon_{c'})
       \le\delta^2/4,                                \tag{3.5}
\]

then

\[
 \boxed{\gamma_n(A_c\cap A_{c'}^c)\ge\delta^2/2.}   \tag{3.6}
\]

Indeed, at correlation zero the two Gaussian halfspace events have cross
mass at least `delta^2`.  For correlation `a` with `|a|<=1/2`, Pinsker's
inequality and
`D(N(0,[[1,a],[a,1]])||N(0,I_2))=-log(1-a^2)/2`
change the probability of any event by at most `|a|`; then use (3.4).

Combining (2.2), Markov, and (3.6), a product-positive proportion of good
Fisher phase pairs has a positive common-source cross incidence whenever
`t/Z_E` and the defects are sufficiently small.  The exact normalized
overlap weight is

\[
 d\Omega(c,c',z)=d\rho_E(c)d\rho_E(c')d\gamma_n(z),
 \qquad
 \omega_{+-}(c,c')=\gamma_n(A_c\cap A_{c'}^c).       \tag{3.7}
\]

Thus there is no missing normalization at the level of *cross-labelled
incidence*.

## 4. The Gaussian channel supplies two other exact normalizations

### 4.1 Common latent input

Generate two observations independently conditional on the same `X`.  The
joint output law `J` has Radon--Nikodym derivative

\[
 \boxed{
 {dJ\over d(\nu\otimes\nu)}(c,c')
 =\Lambda(c,c'):=\int\ell_c(x)\ell_{c'}(x)d\mu(x).}  \tag{4.1}
\]

This overlap is not merely globally normalized; it is doubly stochastic:

\[
 \boxed{
 \int\Lambda(c,c')d\nu(c')=1
 \quad\hbox{for `nu`-a.e. `c`, and symmetrically}.}   \tag{4.2}
\]

This follows from `int ell_c'(x)dnu(c')=int q_t(c'|x)dc'=1`.

However, this coupling need not preserve (2.2).  Put

\[
 M(x)=\int \mathsf B(c)q_t(c|x)dc.                  \tag{4.3}
\]

Then `int M dmu=R` and conditional independence gives the exact identity

\[
 \boxed{
 \mathbb E_J[r(C)r(C')\langle u(C),u(C')\rangle^2]
 =\mathbb E_\mu\|M(X)\|_{HS}^2
 =\|R\|_{HS}^2+\mathbb E_\mu\|M(X)-R\|_{HS}^2.}    \tag{4.4}
\]

The independent-phase value is only `||R||_HS^2`.  Thus the normalized
channel overlap adds a nonnegative, presently uncontrolled phase-clustering
term.  The matrix bound `R preceq tI` gives no upper bound on that term.

### 4.2 One observation and two posterior resamples

The base-space posterior-resampling kernel is

\[
 k_t(x,y)=\int {q_t(c|x)q_t(c|y)\over p_t(c)}dc,
 \qquad K_t(x,dy)=k_t(x,y)d\mu(y).                  \tag{4.5}
\]

It is symmetric and exactly Markov:

\[
 \int k_t(x,y)d\mu(y)=1,
 \qquad
 d\mu(x)K_t(x,dy)=\int d\nu(c)d\mu_c(x)d\mu_c(y).  \tag{4.6}
\]

Its cross-label mass is

\[
 \boxed{
 \int_S K_t(x,S^c)d\mu(x)=\int g(c)(1-g(c))d\nu(c).} \tag{4.7}
\]

This gives normalized positive/negative pairs, but both samples use the
same phase `c`; it contains no pairwise orthogonality from (2.2).

### 4.3 Opposite-label independent channels

There is also an exact product coupling which retains independent phases.
If `X_+~mu(.|S)` and `X_-~mu(.|S^c)` are independent and each is passed
through the observation channel, the output laws are

\[
 d\nu_+(c)={g(c)\over p}d\nu(c),\qquad
 d\nu_-(c)={1-g(c)\over1-p}d\nu(c).                 \tag{4.8}
\]

Conditional on their outputs, the endpoint laws are exactly `mu_c(.|S)`
and `mu_c'(.|S^c)`.  On central near-equality states, (1.4) gives

\[
 0<c_\delta\le r(c)le C_\delta<\infty,             \tag{4.9}
\]

so the Fisher phase law, `nu_+`, `nu_-`, and `nu` are mutually comparable
there (assuming the original label is balanced).  Hence this coupling also
produces a positive density of independently tagged, near-orthogonal
positive/negative endpoint pairs.  Its normalization is exact, not an
overlap estimate.

## 5. Why none of these identities is cross-calibration

For the finite endpoint-defect/Clifford argument one needs, on positive
mass, not only `x in S,y notin S`, but small metric calibration slack such
as

\[
 |x-y|-\{f(x)-f(y)\}=o(t^{-1/2})                    \tag{5.1}
\]

at the long scale `t^{-1/2}` (or the corresponding endpoint feature
concentration in (8.4) of `finite_medial_competitor.md`).  Equations
(2.2)--(4.9) do not bound (5.1).

The incompatibility is precise:

* product Fisher phases preserve angular dispersion and the common Brenier
  source gives (3.6), but two different Brenier maps have no cross-map
  distance or calibration inequality;
* common-latent and posterior-resampling couplings have exact overlap/Markov
  normalization, but (4.4) shows that the former can cluster directions,
  while the latter uses only one direction;
* opposite-label independent channels give actual endpoint pairs and
  independent phase tags, but no small-slack estimate.

Even a proper common-potential Gaussian case shows the missing implication.
Take `W(x)=lambda|x|^2/2`, so the input is a centered Gaussian and
`tau=t+lambda`.  Then

\[
 \mu_c=N(c/\tau,\tau^{-1}I),\qquad
 X_c(G)=c/\tau+G/\sqrt\tau,                          \tag{5.2}
\]

and the common-source coupling satisfies

\[
                         |X_c(G)-X_{c'}(G)|={|c-c'|\over\tau} \tag{5.3}
\]

identically.  For the fixed parity set
`S={x_1x_2>=0}`, tilts with means `(a,0)` and `(0,a)` have asymptotically
sharp orthogonal centroid directions as `a sqrt(tau)->infinity`, while
(5.3) is `sqrt(2)a`, an arbitrarily large multiple of the posterior scale
`tau^{-1/2}`.  Taking `lambda/t->0` also makes the centroid ratio normalized
with the localization curvature `t` arbitrarily sharp.  Thus near-halfspace
structure plus orthogonality and common
Brenier-source overlap do not imply cross calibration, even for one common
smooth convex potential and one fixed set.

**Verdict.**  The Gaussian channel does add the exact normalization which
simple Hellinger or midpoint tilt overlap lacks: (4.2), (4.6), and (4.8)
are dimension-free identities.  But normalization and Fisher angular
dispersion live in incompatible couplings, and the product coupling that
retains both labels and dispersion still lacks (5.1).  The remaining input
is exactly a positive-density endpoint rematching/calibration theorem; the
channel identities do not prove the finite endpoint-defect lemma.

## 6. Audit status of the spatial note

The algebraic channel identities (2.2), (4.1)--(4.9) above have complete
one-line conditional-expectation proofs.  The following pieces of
`spatial_phase_coherence.md` have **not** received an independent clean-room
audit and should not be treated as load-bearing without one:

1. Lemma 2.1's quantitative strip argument, especially the pairing step in
   (2.7).
2. Lemma 3.1's tail modulus and its use in the Brenier product estimate.
3. The smooth-potential midpoint transfer (4.6)--(4.9), which is only a
   bulk certificate and in any event misses hard support facets.
4. The fixed-dimensional halfspace line estimate (5.5); the asserted
   uniform lower density/central-quantile geometry was sketched, not fully
   constant-tracked.
5. The equality splitting assertion in Section 6 was cited at the level of
   the equality case of Lichnerowicz, not reproved under all approximation
   and lower-dimensional-support cases.

Accordingly, the spatial note's safe current output is a candidate local
phase/product mechanism plus an explicit global obstruction, not an audited
KLS lemma.
