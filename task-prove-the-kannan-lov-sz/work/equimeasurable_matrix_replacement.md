# Removing the spatial selector: unweighted coarea variance and joint matrix fidelity

## 0. Outcome

The analytic selector is needed to *discover* a ranked physical packet, but
it is not needed to retain a fixed amount of angular variance.  The selected
normal law is a thinning of the full coarea normal law.  Conditional-variance
monotonicity therefore gives a completely unweighted angular-variance lower
bound.  This removes the `nabla omega` term which blocked the spatially
weighted matrix-fidelity problem.

A joint equimeasurable `BV` minimization then preserves the unweighted normal
matrix of the entire heat function with error `Delta/kappa`, where `Delta`
is the integrated profile deficit.  The integrand has a constant anisotropy,
not a spatially varying one.  These two statements are proved below.

This is not yet a KLS proof.  The joint minimizer is a nested family, and
global equimeasurable minimality has not been shown to make almost every
level an unconstrained anisotropic isoperimetric minimizer.  Independently
minimizing the levels destroys the heat-phase incidence.  This exact nesting
versus stationarity issue is the remaining transfer problem for this route.

## 1. Thinning cannot create projector variance

Let `(Omega,Lambda)` be a finite measure space, let `Z:Omega->H` take values
in a Hilbert space, and assume `|Z|=1` almost everywhere.  Let
`0<=omega<=1`, put

\[
 T=\Lambda(\Omega),\qquad t=\int\omega\,d\Lambda,
\]

and suppose `t>0`.  Denote by `P` the normalized law `Lambda/T` and by
`P_s` the selected law `omega Lambda/t`.

**Lemma 1.1 (soft thinning).**

\[
 \boxed{
 \operatorname {Var}_{P}(Z)
 \ge {t\over T}\operatorname {Var}_{P_s}(Z).}       \tag{1.1}
\]

**Proof.**  Enlarge the probability space by a Bernoulli variable `I` with
`P(I=1|x)=omega(x)`.  Then `P(I=1)=t/T` and the law of `Z` conditional on
`I=1` is `P_s`.  The Hilbert-space law of total variance gives

\[
 \operatorname {Var}(Z)
 =E\operatorname {Var}(Z\mid I)
  +\operatorname {Var}(E[Z\mid I])
 \ge P(I=1)\operatorname {Var}(Z\mid I=1).
\]

QED.

Apply this with the central truncated heat function.  Let
`C=[r_-,r_+]` be the audited central level interval and let `h_C` be the
nondecreasing clipping map with `h_C'=1_C`.  Put `F=h_C(F_0)`.  Then

\[
 d\Lambda=|\nabla F|d\mu,
 \qquad Z=\theta\theta^T,
 \qquad \theta={\nabla F\over|\nabla F|},           \tag{1.2}
\]

viewing symmetric matrices with the Hilbert--Schmidt norm.  Since `Z` is a
rank-one projector, if

\[
 Q={1\over T}\int Z\,d\Lambda,
 \qquad Q_s={1\over t}\int\omega Z\,d\Lambda,
\]

then

\[
 \operatorname {Var}_{P}(Z)=1-\operatorname {tr}Q^2,
 \qquad
 \operatorname {Var}_{P_s}(Z)=1-\operatorname {tr}Q_s^2. \tag{1.3}
\]

The audited floored analytic selector on the central interval satisfies

\[
 t>.0048495p,\qquad {t\over\|tQ_s\|_{op}}>17.86.    \tag{1.4}
\]

Because `Q_s` is positive semidefinite with trace one,
`tr Q_s^2<=||Q_s||op<1/17.86`.  Heat contraction gives
`T=int|nabla F|dmu<=int R dmu<=p`.  Hence (1.1)--(1.4) yield the explicit
unweighted statement

\[
 \boxed{
        1-\operatorname {tr}Q^2
        >.0048495(1-1/17.86)>.0045779.}             \tag{1.5}
\]

No derivative, lower bound, or regularity of `omega` occurs in (1.5).

## 2. The unweighted `BV` normal matrix

Work first on a bounded convex support with a positive continuous density.
For `G in BV(mu)` let

\[
 DG=\sigma_G|DG|,qquad |\sigma_G|=1\quad |DG|\text{-a.e.},
\]

and define

\[
 \operatorname {TV}_\mu(G)=|DG|_\mu(E),
 \qquad
 M(G)=\int\sigma_G\sigma_G^T\,d|DG|_\mu.           \tag{2.1}
\]

Thus `tr M(G)=TV_mu(G)`.  Coarea gives

\[
 M(G)=\int_{\mathbb R}M(\{G>r\})\,dr,              \tag{2.2}
\]

with the same orientation convention on both sides; the projector removes
the harmless sign.

Fix `0<kappa<1/3`.  Among all `[0,1]`-valued `G` equimeasurable with `F`,
minimize

\[
 \mathcal J_\kappa(G)=\operatorname {TV}_\mu(G)
       +\kappa\|M(G)-M(F)\|_*.                     \tag{2.3}
\]

Equimeasurability means
`mu(G>r)=mu(F>r)` for every continuity point of the common distribution.

**Lemma 2.1 (existence and lower semicontinuity).**  Problem (2.3) has a
minimizer `G_kappa`.

**Proof.**  By nuclear/operator duality,

\[
 \mathcal J_\kappa(G)
 =\sup_{H=H^T,\ \|H\|_{op}\le1}
 \left\{\int\Phi_H(dDG)
       -\kappa\operatorname {tr}(HM(F))\right\},   \tag{2.4}
\]

where the one-homogeneous integrand is

\[
 \Phi_H(\xi)=|\xi|+
       \kappa{\xi^TH\xi\over|\xi|}.                \tag{2.5}
\]

For a tangent vector `h perp n`, the spherical Hessian of (2.5) is

\[
 |h|^2+\kappa\{2h^THh-(n^THn)|h|^2\}
 \ge(1-3\kappa)|h|^2.                              \tag{2.6}
\]

Thus every `Phi_H` is convex, one-homogeneous, and uniformly elliptic.
Each integral in (2.4) is `BV` lower semicontinuous, and so is their
supremum.  A minimizing sequence has uniformly bounded total variation.
`BV` compactness gives an `L^1(mu)` limit after passage to a subsequence.
Convergence in measure preserves the pushforward distribution, so the limit
is equimeasurable with `F`.  QED.

On an unbounded log-concave support, local `BV` compactness plus tightness of
the probability gives the same conclusion.  For an extended-valued convex
potential one first works in the relative interior and then uses monotone
convex truncation.  A final proof would need to spell out the boundary trace
in this last passage.

## 3. Exact matrix retention from the profile deficit

Let

\[
 \Delta_F=\operatorname {TV}_\mu(F)
 -\int_{\mathbb R}I_\mu(\mu(F>r))\,dr.              \tag{3.1}
\]

This is nonnegative by coarea.  It is precisely the integrated
isoperimetric deficit.  With `F=h_C(F_0)`, it integrates only the central
levels, and the fixed-scale physical-splicing audit bounds it by
`6.02*10^-5p`.  If the heat fraction is retuned, the trace/rank and deficit
estimates must be rederived together; no independent limiting assertion is
being used here.

**Lemma 3.1 (joint rank preservation).**  Every minimizer in Lemma 2.1
satisfies

\[
 \boxed{
 \operatorname {TV}_\mu(G_\kappa)\le
       \operatorname {TV}_\mu(F),
 \qquad
 \|M(G_\kappa)-M(F)\|_*\le{\Delta_F\over\kappa}.}  \tag{3.2}
\]

**Proof.**  Comparison with `F` gives

\[
 TV(G_\kappa)+\kappa\|M(G_\kappa)-M(F)\|_*
 \le TV(F).                                        \tag{3.3}
\]

Equimeasurability and coarea give

\[
 TV(G_\kappa)
 \ge\int I_\mu(\mu(G_\kappa>r))dr
 =\int I_\mu(\mu(F>r))dr=TV(F)-\Delta_F.          \tag{3.4}
\]

Substitute (3.4) in (3.3).  QED.

There is a direct variance-stability estimate.  Put `T=tr M(F)`,
`T'=tr M(G_kappa)`, `Q=M(F)/T`, and `Q'=M(G_kappa)/T'`.  If
`Delta_F<T/2`, then (3.2) implies

\[
 \|Q'-Q\|_*
 \le {\Delta_F\over T}\left(1+{1\over\kappa}\right),           \tag{3.5}
\]

and hence

\[
 |\operatorname {tr}(Q'^2)-\operatorname {tr}(Q^2)|
 \le {\Delta_F\over T}\left(1+{1\over\kappa}\right).          \tag{3.6}
\]

Indeed `0<=T-T'<=Delta_F` and

`Q'-Q=(M'-M)/T+(1/T'-1/T)M'`,

which proves (3.5).  Both normalized matrices are positive contractions of
trace one, so `E=Q'-Q` has trace zero and

`tr(Q'^2-Q^2)=tr[E(Q'+Q-I)]`.

The operator norm of `Q'+Q-I` is at most one, proving (3.6).  Thus the fixed lower bound (1.5)
survives whenever the right side of (3.6) is below `.0045779`.

At the currently frozen constants this test does not pass.  Since
`T>=t>.0048495p` and `Delta_F<=6.02*10^-5p`, the audited ratio can be as
large as `.012414`.  Even the endpoint `kappa=1/3` makes the right side of
(3.6) about `.049656`, far above `.0045779`.  Thus the joint functional is a
structural removal of `nabla omega`, not a numerical closure with the
present parameters.  It requires a fresh fixed-scale run in which
`Delta_F/T` is below approximately `1.14*10^-3`, while the selected trace
and variance lower bounds remain fixed.  Sections 7 and
`equimeasurable_audit.md` give audited retunings which pass this matrix
test; the stricter Wulff-profile test requires the smaller anisotropy used
in Section 7.

## 4. Constant anisotropy and the nesting obstruction

At a minimizer, choose a nuclear-norm subgradient `H`, `||H||op<=1`.
For first variations, the objective is the anisotropic total variation with
the **constant** surface tension (2.5).  In particular, there is no term
containing the spatial derivative of the heat selector.  This is the main
gain over the weighted fidelity in `rank_preserving_eikeland.md`.

It is tempting to conclude that for almost every `r`, the set
`{G_kappa>r}` minimizes the anisotropic perimeter `P_H` at its prescribed
volume.  That conclusion is not automatic.  The admissible family in (2.3)
is nested because it is generated by one function.  Replacing one level by
an arbitrary same-volume anisotropic minimizer can violate nesting with all
neighboring levels.  Global equimeasurable stationarity initially gives only
variations realized by a common deformation of the function.

For a smooth `G` with nonvanishing gradient, a level-dependent normal flow
localized to a small interval of values suggests the Euler equation

\[
 H_{\Phi_H}(\partial\{G>r\})=\lambda(r)             \tag{4.1}
\]

for almost every regular `r`.  A rigorous derivation must restore the
distribution after the flow, control crossings of neighboring levels, and
pass through critical values and `BV` singularities.  No such theorem is
proved here.  Even (4.1) would be anisotropic CMC; one must still extend the
killed-tube identity to the uniformly elliptic constant tension (2.5), or
compare it quantitatively to a fixed affine metric.

Independently minimizing each level by the unweighted matrix-fidelity
functional does preserve its normal matrix without a selector derivative.
It also destroys the cross-level phase labels, collision structure, and
common-selector softmax law.  That operation therefore cannot be substituted
silently for the joint problem.

The exact next lemma required by this route is a dimension-free statement
that a minimizer of (2.3) has, outside a profile-deficit charge, levelwise
constant-anisotropic CMC geometry with a killed-tube formula and preserves
the physical incidence of neighboring levels.  This is narrower than the
old spatial-selector transfer but remains load bearing.

## 5. Exact Wulff-normal tube algebra

The constant anisotropy admits an exact tube calculus; it need not be
treated as an error in the Euclidean CMC equation.  Let `Phi` be a positive,
even, `C^2` away from zero, convex one-homogeneous function whose Hessian is
positive definite on every tangent space.  For a smooth hypersurface with
unit normal `n`, put

\[
 z=D\Phi(n),\qquad G=D^2\Phi(n)|_{n^\perp}.          \tag{5.1}
\]

Euler homogeneity gives `z cdot n=Phi(n)` and `G n=0`.  If `S` is the
Euclidean shape operator, the anisotropic parallel map is

\[
                         T_t(x)=x+t z(x).             \tag{5.2}
\]

Its tangential differential is `I+tGS`.  The normal remains `n` along the
ray, because `n cdot G=0`.  Up to the first focal, support-contact, or cut
time, the weighted flux Jacobian relative to
`Phi(n)e^{-V(x)}dH^{k-1}(x)` is

\[
 j_x(t)=\det(I+tGS)
       \exp\{-V(x+tz)+V(x)\}.                       \tag{5.3}
\]

Since `GS` is similar to the symmetric matrix `G^(1/2)SG^(1/2)`, its
eigenvalues are real and

\[
 {d^2\over dt^2}\log j_x(t)
 =-\operatorname {tr}
   [((I+tGS)^{-1}GS)^2]
   -\langle\nabla^2V(x+tz)z,z\rangle\le0.           \tag{5.4}
\]

The initial logarithmic slope is

\[
 {d\over dt}\bigg|_{t=0}\log j_x(t)
 =\operatorname {tr}(GS)-\langle\nabla V,z\rangle
 =:H_{\Phi,\mu}(x).                                 \tag{5.5}
\]

Thus, if a level has constant weighted anisotropic mean curvature
`H_{Phi,mu}=lambda`, then

\[
                         j_x(t)=e^{\lambda t-D_x(t)},
 \qquad D_x(t)\ge0,                                  \tag{5.6}
\]

with exactly the same killed-ray monotonicity as in the Euclidean tube
theorem.  Collisions, focal loss, and hard-support contact are stopped in
the same way.  The volume derivative of the anisotropic offset is its
anisotropic perimeter because the normal flux is
`z cdot n=Phi(n)`.

For (2.5),

\[
 (1-\kappa)|\xi|\le\Phi_H(\xi)\le(1+\kappa)|\xi|,
 \qquad
 |D\Phi_H(n)-n|\le C\kappa                         \tag{5.7}
\]

with a numerical `C`, and (2.6) gives the required ellipticity.  Hence
anisotropic perimeter and Wulff displacement length are comparable to their
Euclidean counterparts.  Matrix comparison is only in nuclear or
Hilbert--Schmidt norm: in general `z tensor z` and `n tensor n` do not
Loewner-dominate one another.  The bound `|z-n|<=C kappa` implies that their
normalized projector laws and angular variances differ by `O(kappa)`, which
is useful only when `kappa` is chosen much smaller than the fixed variance
seed.

Equations (5.3)--(5.6) are pointwise smooth identities.  A complete use
still needs a rectifiable killed-tube theorem for singular anisotropic
minimizers and the levelwise stationarity assertion discussed in Section 4.

## 6. Formal equimeasurable Euler equation

There is a precise smooth candidate for the missing stationarity theorem.
Suppose `G` is smooth, `|nabla G|>0` on a central value band, and its value
law has a positive smooth density there.  The first-order tangent space to
the equimeasurability constraint consists of perturbations `h` satisfying

\[
                         E[h\mid G]=0.                \tag{6.1}
\]

Indeed differentiation of `int q(G+t h)dmu` for every smooth `q` gives
necessity, and monotone distributional rearrangement of `G+t h` back to the
law of `G` has derivative `h-E[h|G]`, giving sufficiency at the formal
level.  Stationarity of the constant-anisotropic total variation then says

\[
 \int h\,[-\operatorname {div}_\mu D\Phi_H(\nabla G)]d\mu=0
 \quad\hbox{whenever }E[h\mid G]=0.                 \tag{6.2}
\]

The orthogonal-complement identity in `L^2(mu)` yields a measurable scalar
function `lambda` with

\[
 -\operatorname {div}_\mu D\Phi_H(\nabla G)
                         =\lambda(G).                \tag{6.3}
\]

Because `D Phi` is zero-homogeneous and
`D^2Phi(n)n=0`, the left side restricted to the level `G=r` is precisely
its weighted anisotropic mean curvature.  Thus (6.3) implies
`H_{Phi,mu}=lambda(r)` on every regular level in the band.

This calculation explains why nesting need not destroy CMC at smooth
regular levels.  It is not yet a theorem for the `BV` minimizer: exact
distribution-restoring paths, integrability of the multiplier, critical
levels, atoms created by clipping, and anisotropic partial regularity all
require proof.  Those are now the well-defined analytic tasks, rather than
a spatial-selector derivative estimate.

## 7. Candidate fixed-scale retuning

The numerical failure after (3.6) comes from the frozen value
`sqrt(alpha)=10^-5`, not from a power mismatch.  The explicit estimates in
`fixed_scale_wedge_extremality.md` retain their parameters:

\[
 {\mathcal D_{co}(F_s)\over p}
 \le\beta+{4\sqrt\alpha\over c_G(1-\beta)},
 \qquad c_G=\sqrt{2/\pi},                           \tag{7.1}
\]

whereas the good physical flux trace is
`a_-(alpha)/(8 pi)` times `p` up to the displayed alignment losses.  Here

\[
 M_\alpha=\sqrt{2\log {8I_0\over c_0\sqrt\alpha}},
 \qquad
 a_-(\alpha)\ge {3\over4}{M_\alpha\over1+M_\alpha^2}.             \tag{7.2}
\]

Thus the trace deteriorates only like `1/sqrt(log(1/alpha))`, while the
deficit is `O(sqrt(alpha))`.

The earlier provisional choice `kappa=1/4` is not compatible with the tube
step.  Euclidean and anisotropic profiles are known here only within the
multiplicative factors `1+-O(kappa)`, which would insert an order-`p`
profile loss at that value.  A viable hierarchy must make `kappa` tiny and
then make the coarea deficit much tinier still.  For example, take
provisionally

\[
                         \alpha=10^{-28},
 \qquad \beta=10^{-14},\qquad \kappa=10^{-6}.        \tag{7.3}
\]

Then (7.1) is less than `6.02*10^-14`.  Formula (7.2) gives
`M_alpha<8.46`.  The audited Mills estimate gives
`a_-(alpha)>.0875`; after the analytic floor, the endpoint deletion, and
the alignment deficits this leaves

\[
                         T\ge t>.0032827p.            \tag{7.4}
\]

Consequently

\[
 {\Delta_F\over T}\left(1+{1\over\kappa}\right)
 <{6.02\,10^{-14}\over.0032827}(1+10^6)
 <1.84\,10^{-5}.                                    \tag{7.5}
\]

The first inequality in (7.5) uses the sharper normalization estimate
`||Q'-Q||_*<=(Delta_F/T)(1+1/kappa)`; the purity change is of the same
order because both normalized matrices have trace one.  The audited
thinning variance at the same parameters is above `.0032688` once the
posterior effective rank is at least `500`.  Even using a four-times looser
matrix-retention loss and then the explicit Wulff direction comparison
leaves variance above `.0031535`; the subsequent `h=1` killed-tube deletion
leaves more than `.003036`.  The anisotropic-profile loss is `O(10^-6)p`
provided the relevant Euclidean leaves have relative profile error at most
`10^-6`.  The rank condition holds for
`K>=1.81*10^47`, an enormous but universal threshold.  Below that
threshold, Buser--Ledoux already gives a universal lower bound for `p`; its
size is irrelevant to dimension-freeness.

This constant chain has been independently checked in
`wulff_tube_audit.md`, but it is conditional rather than a KLS proof.  The
audit does not derive the required leafwise profile error from the
coarea-integrated deficit, nor does it establish singular Wulff coverage.
The point of (7.1)--(7.5) is only that the old fixed-power and fixed-
anisotropy numerical obstructions disappear after the spatial-selector
derivative is removed.
