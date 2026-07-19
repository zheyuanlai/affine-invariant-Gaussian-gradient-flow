# Checkpoint 05: general angular stability, Fisher--ray constraints, and exact global obstructions

## 1. Exact results added since Checkpoint 04

### 1.1 Clean-room posterior stability

Let `pi` be `t`-strongly log-concave on its affine hull, let
`g=pi(S) in [delta,1-delta]`, and put

\[
 v=\operatorname {Cov}_\pi(1_S,X),\qquad
 \varepsilon=1-\frac{\sqrt t|v|}{I(g)}.
\]

An independent proof now establishes all of the following, including hard
support and lower-dimensional affine hulls.

1. If `H` is the active halfspace with the same mass as `S`, then
   \[
   \pi(S\mathbin\triangle H)
       \le \sqrt{\frac{2}{\pi\delta}}\sqrt\varepsilon .       \tag{1.1}
   \]
2. If `u=v/|v|`, then
   \[
   1-t\operatorname {Var}\langle u,X\rangle
       \le C_\delta\varepsilon\sqrt{\log(e/\varepsilon)}.    \tag{1.2}
   \]
3. If `T` is the standardized Brenier contraction from a standard Gaussian
   and `H_T` is its Hessian, then
   \[
   E|H_Tu-u|^2
       \le C_\delta\varepsilon\sqrt{\log(e/\varepsilon)}.    \tag{1.3}
   \]
   The corresponding product-map approximation loses at most twice this
   amount.
4. In fixed dimension at most two, two central halfspaces under an isotropic
   log-concave law have oriented normal/offset distance at most
   `C_delta` times their symmetric-difference mass.

The first estimate has an exact one-dimensional quantile proof.  The third
must use the standardized map; its unstandardized analogue is false.

### 1.2 Dimension-free general angular stability

Write, with `m=E_pi X`,

\[
 D=E[(1_S-g)(X-m)(X-m)^T],\qquad P=I-uu^T.
\]

Combining the preceding clean-room theorem with the independently proved
Gaussian-halfspace pullback theorem gives

\[
 \boxed{\|PD\|_{HS}^2\le C_\delta t^{-2}
       \Omega_\delta(\varepsilon)},                 \tag{1.4}
\]

where

\[
 \rho_\delta(r)=C_\delta\{\sqrt r+
       (r\sqrt{\log(e/r)})^{1/3}\},\qquad
 \Omega_\delta(r)=r+\sqrt{\rho_\delta(r)}.          \tag{1.5}
\]

In particular
`Omega_delta(r)=O_delta(r^(1/6)log(e/r)^(1/12))`.
This lemma has passed an independent audit.  The two possible hidden
dimension losses are absent:

* the pullback of the target threshold differs from a source Gaussian
  halfspace by at most `2(zeta/h^2+2 phi(0)h)`;
* for every symmetric `M` with `||M||HS=1`, target Poincare gives
  `Var(X^T M X)<=4`, so the correlated second moments are stable in
  Hilbert--Schmidt norm with error `2 sqrt(rho)`.

The change of active direction costs no `sqrt(n)`, because the difference of
the two rank-one projections has rank at most two and the comparison matrix
has operator norm at most one.  At a central state, (1.4) gives

\[
 \frac{\|PD\|_{HS}^2}{|v|^2}
       \le \frac{C_\delta}{t}\Omega_\delta(\varepsilon).      \tag{1.6}
\]

Thus the formerly missing *pointwise* transverse estimate is now proved for
an arbitrary set and arbitrary strongly log-concave posterior.  Its weak
modulus does not by itself resolve the long-time scale mismatch.

### 1.3 Exact stochastic and heat-flow defect identities

For ordinary Gaussian localization, let
`r=|v|`, `u=v/r`, `P=I-uu^T`, and
`eta=sqrt(t)r/I(g)`.  With `M=sqrt(t)I(g)` and
`Delta=M-r`, the exact semimartingale identity is

\[
 d\Delta=d\mathcal M-left\{
 \frac{\|PD\|_{HS}^2}{2r}
 +r\,u^T(t^{-1}I-A)u
 +\frac{I(g)(1-\eta)^2}{2t^{3/2}}
 \right\}dt.                                        \tag{1.7}
\]

All three drift terms are nonnegative.  This simultaneously records
transverse rotation, active marginal non-Gaussianity, and scalar amplitude
defect.  The heat-flow picture has an independently derived pointwise
Bernstein identity whose integral is an exact sum of squares, together with
temporal gluing and the exact comparison to (1.7).

The global-minimizer argument yields a genuine posterior-resampling gap and
a global Cheeger contraction.  It also exposes a sharp obstruction: the
natural profile/Dirichlet comparison pays an amplitude of order
`sqrt(alpha)`, whereas the available deficit is order `alpha`.  This is not
an omitted estimate; Gaussian halfspaces attain the square-root scale.

### 1.4 Binary Fisher information and original calibrated rays

For the Gaussian observation channel, the binary Fisher matrix satisfies

\[
 R_t=E[\nabla g_t\nabla g_t^T/(g_t(1-g_t))]\preceq tI,        \tag{1.8}
\]

and the binary I--MMSE identity is exact.  On a plateau with fixed Fisher
trace, (1.8) forces effective rank of order the observation scale.

The posterior observation can be disintegrated over the *original*
calibrated rays `X=z_q+T N_q`.  For `B=1_{T>0}`, posterior ray bias `b_q`,
`k_q=Cov(B,T|c,q)`, synergy variance
`zeta=Var_q b_q`, and posterior active direction `u`, one has

\[
\begin{split}
 &[I(g)-E I(b_q)]
 +E[I(b_q)-\sqrt t\,k_q]\\
 &\qquad+E[\sqrt t\,k_q(1-|u\cdot N_q|)]
 \le \delta+\sqrt\zeta .                            \tag{1.9}
\end{split}
\]

Hence small centroid defect plus small conditional ray-identity information
forces both one-dimensional equality and projective alignment on the
original quotient.  Constant channel information at `t=1/K` also forces a
fixed quotient mass of original rays with conditional variance at least
`cK`; a Fisher plateau supplies the required information.  Pairwise phase
orthogonality is quantified by

\[
 E_{\lambda\otimes\lambda}\langle u,u'\rangle^2
   ={\operatorname {tr}R_t^2\over(\operatorname {tr}R_t)^2}
   \le {t\over\operatorname {tr}R_t}.                \tag{1.10}
\]

The conversion from posterior phase mass back to an extremal surface can
lose a likelihood-overlap factor which may be exponential.  That loss has
not been removed.

### 1.5 Exact synergy height and its no-go model

On calibrated rays of length scale `L`, the centered height

\[
 h_c(q)=Lw_c(q)(b_q-g(c))                            \tag{1.11}
\]

has a rigorously controlled finite switching cost.  If the relevant median
ray densities are at least `beta/L`, true T3 extremality implies

\[
 D_f(h_c)\ge {\kappa\beta L\over8}\zeta(c)^2.        \tag{1.12}
\]

Its average is exactly

\[
 E_Ch_C(q)=L\{1/2-E[g(C)|Q=q]\}
          =L E[(I-A^*A)B|Q=q].                       \tag{1.13}
\]

An isotropic log-concave diamond model has positive synergy but zero
average height.  Uniform dilation preserves the dimensionless synergy while
scaling the least defect linearly.  Therefore channel identities plus a
single quotient-height competitor cannot alone produce an `o(L)` cost.
The model is not a T3 extremizer, so an extremality-dependent phase-capacity
theorem remains logically possible.

### 1.6 A high-rank gradient wedge inequality

Let a probability `nu` have Poincare constant `C`, let `F` be a Sobolev
vector field, set `m=EF`, `R=E[FF^T]`, and `H=grad F`.  Exterior-product
Poincare gives the exact inequality

\[
\begin{split}
 &(\operatorname {tr}R)^2-\operatorname {tr}(R^2)
 -\{|m|^2\operatorname {tr}R-m^TRm\}\\
 &\qquad\le C E\operatorname {tr}
       [H^T(\operatorname {tr}R\,I-R)H].             \tag{1.14}
\end{split}
\]

If `m=0` and `||R||op<=kappa tr R` with `kappa<1/2`, then

\[
 E\|H\|_{HS}^2\ge{(1-2\kappa)\operatorname {tr}R\over C}.    \tag{1.15}
\]

For the smoothed observation law `Y=X+sqrt(s)G`, its Poincare constant is at
most `s+C_P(mu)`, and the field `F=sqrt(s) grad h` has second moment equal to
the binary Fisher matrix.  At `s` comparable to a hypothetical bad scale
`K=C_P(mu)`, high Fisher rank therefore forces Hessian energy of order
`K^{-2}`.  The longitudinal term in the eikonal change of variables has
exactly the same order, so (1.15) does not yet yield positive excess cost.

### 1.7 Exact convex-nullity rigidity and failure of its quantitative inverse

For a convex potential/Brenier contraction, genuine zero variance defect in
a direction forces a global affine direction.  A weighted family of exact
zero directions forces a global cylinder and Gaussian factor containing the
span of those directions, including nonsmooth hard support.

Two plausible approximate upgrades are false without an overlap hypothesis.
First, a bounded bulk resolvent misses singular curvature: smooth
approximations of `|x|` give an explicit cusp counterexample.  Second,
polyhedral maxima have many exact affine cells whose local fixed directions
do not globalize.  Thus a quantitative cylinder-or-focus theorem must pay a
posterior overlap/incidence cost; local nullity alone cannot supply it.

## 2. Newly blocked families

1. **Scalar additive-window baseline.**  Scalar admissibility does not
   select a Gaussian profile.  Compact perturbations preserve every proposed
   scalar constraint while changing the defect decay, so a unique-baseline
   argument cannot close the scale integral.
2. **Ordinary resampling/amplitude closure.**  The sharp spatial competitor
   has `sqrt(alpha)` size when the profile deficit is `alpha`.
3. **Projection-height synergy.**  Positive conditional phase information
   can average to zero on the quotient, and scaling defeats any channel-only
   `o(L)` estimate.
4. **Channel normalization without metric calibration.**  Fisher effective
   rank controls phase diversity but not the physical overlap/conductance
   needed to compare posterior cells.
5. **Approximate convex nullity from local cells.**  Polyhedral affine cells
   refute globalization without incidence or common posterior mass.
6. **Bulk Hessian resolvent through singular limits.**  A convex cusp is
   invisible to the proposed bounded resolvent.
7. **Moment-map or bounded-degree phase seeds.**  Gaussian parity and
   exponential-tail nonanalyticity remain decisive countertests.
8. **Random slabs/ordinary block descent.**  The conditional score, boundary
   flux, radial slow mode, and shifted-exponential fiber tail prevent the
   required separated recurrence.

## 3. Active incompatible mechanisms

1. **Phase partition/capacity.**  Combine (1.4), (1.7), and high Fisher
   effective rank.  The precise live target is a dimension-free theorem that
   either charges the interfaces between many posterior directions to the
   profile defect or forces a coherent Gaussian/product factor.  Every
   scale parameter `t= tau/K` must remain explicit.
2. **Gradient multiplicity rigidity.**  Determine whether the wedge lower
   bound can be localized to the transverse Bernstein square, rather than
   being cancelled by the longitudinal eikonal term.  Equality and
   convolution-Poincare saturation are being classified.
3. **Fisher-weighted focal coherence.**  Use the fixed original long-ray
   mass, projective phase alignment, and exact-nullity cylinder theorem to
   obtain the missing overlap/incidence.  The polyhedral-cell no-go must be
   excluded using true extremality, not local convexity.
4. **Temporal integration of angular stability.**  Seek a weighted
   occupation identity that integrates the weak modulus in (1.4) without
   selecting an increasingly long time window.  A bare pigeonhole is known
   to lose the required scale.

## 4. Updated extremal constraint set

A hypothetical isotropic T3 sequence with scale `L` tending to infinity must
now satisfy all of the following.

1. It has central near-minimizing cuts with fixed upper and lower tail mass,
   and calibrated original rays of variance comparable to `L^2` on a fixed
   quotient fraction.
2. At every low-defect strongly log-concave posterior, the cut is close to a
   halfspace, its active marginal is nearly Gaussian, its Brenier map nearly
   splits, and its transverse correlated Hessian obeys the dimension-free
   modulus (1.4).
3. The localization direction can rotate only by paying the first square in
   (1.7).  Consequently a bad sequence must distribute its defect over time
   so that no fixed window simultaneously has enough profile mass and small
   angular cost.
4. Whenever binary Fisher information stays bounded below at scale
   `t=1/L^2`, its phases have effective rank of order `L^2`, while a fixed
   fraction of the original calibrated rays are genuinely length `L`.
5. Those many phases cannot be reconciled by a single quotient height, a
   bounded-degree moment seed, a scalar Gaussian baseline, or merely local
   affine cells.
6. Exact zero transverse defects would force a genuine Gaussian cylinder.
   Hence every nonproduct survivor must retain positive defects on an
   incidence set, but present estimates do not lower-bound its mass.
7. Linear, translated-radial, product, and coherent Gaussian-factor branches
   remain dimension-free by isotropy, translated thin shell,
   Bobkov--Houdre tensorization, and one-dimensional Gaussian isoperimetry.
   A survivor must be nonlinear, nonradial, irreducible, and phase-incoherent.

## 5. Audit status

There is still no candidate complete proof.  The clean-room posterior
stability theorem and the general angular-stability composition have been
independently rederived, and the wedge inequality has received an independent
line audit.  Hard support and lower-dimensional affine hulls are included in
those statements.  The stochastic and heat-flow identities have exact sign
and scaling checks on Gaussian halfspaces, nonhalfspaces, product
exponentials, simplices, and parity/radial tests.

Symbolic dimension tracking presently stops at the following load-bearing
global statements, none of which is being assumed:

1. a phase-partition/capacity estimate converting high Fisher rank to
   interface cost;
2. a temporal occupation estimate integrating the weak angular modulus with
   no growing-window loss;
3. a Fisher-weighted focal-incidence theorem excluding the polyhedral-cell
   escape;
4. a transverse coercivity theorem separating the wedge Hessian energy from
   the longitudinal eikonal term.

Until one of these supplies a dimension-free seed and the other audit stages
are passed, the KLS conclusion has not been obtained.
