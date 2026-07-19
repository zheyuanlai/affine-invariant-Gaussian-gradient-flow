# Locked transverse localization and balanced one-dimensional fibers

## Executive conclusion

This note proves the fixed-line branch of hard mass-preserving stochastic
localization, with explicit constants.

Let `S` have mass `1/2` under an isotropic log-concave probability `mu`, and
run the hard control

\[
 C_t=P_{e_t^\perp},\qquad
 e_t=\frac{\operatorname {Cov}_{\mu_t}(1_S,X)}
          {|\operatorname {Cov}_{\mu_t}(1_S,X)|}.
\]

Suppose that, for a deterministic line `Ru`, almost surely throughout a time
interval of length `T`,

\[
 \sqrt{1-\langle e_t,u\rangle^2}\leq\delta<1/\sqrt{2} .
 \tag{0.1}
\]

Then

\[
 \mu^+(S)\ge \frac{1}{\sqrt{96}}
 \left[
 \frac{1}{T(1-\delta^2)}+
 \frac{1-\delta^2}{1-2\delta^2}
 \left(1+
 \frac{\delta}{(1-\delta^2)\sqrt{T}}\right)^2
 \right]^{-1/2}.                                      \tag{0.2}
\]

In particular,

\[
 \mu^+(S)\ge \frac{1}{\sqrt{96}}-o_{\delta\downarrow0,
 T\uparrow\infty}(1).                                  \tag{0.3}
\]

The same result holds conditionally after a stopping time.  If `U` is the
line selected by the past, `K=Var_{mu_tau}<X,U>`, and the future locked-line
event has conditional probability `p`, the exact lower bound is

\[
 \frac{1}{\sqrt{96}}
 \frac{p}{\sqrt{
       \alpha^{-1}+\kappa^2
       (\sqrt{K/p}+r/\sqrt{\alpha})^2}},               \tag{0.4}
\]

where

\[
 \alpha=T(1-\delta^2),\qquad
 r=\frac{\delta}{\sqrt{1-\delta^2}},\qquad
 \kappa=\sqrt{\frac{1-\delta^2}{1-2\delta^2}}.        \tag{0.5}
\]

Thus the exact historical quantity needed to deal with early direction
selection is, in the ideal limit,

\[
 \boxed{\quad
 \mathbb E\,\frac{p_\tau^{3/2}}{\sqrt{K_\tau}}\geq c.
 \quad}                                                \tag{0.6}
\]

For exact locking one may replace `K_tau` by the smaller mean conditional
fiber variance.  Estimate (0.6), or its finite-`T`, finite-`delta` version
(0.4), is not supplied by isotropy: the line is selected from the same
history as the posterior.  A generic adaptive selected-variance estimate is
false for product exponentials.  What is still needed is a theorem coupling
large selected variance to the probability of *subsequent persistent
label-locking*.

There is also a clean geometric statement.  If the protected line is
exactly fixed for all future observation outcomes, then the original
conditional measure of `S` on almost every parallel line is exactly
one-half.  This follows from injectivity of Gaussian convolution, not from
an informal limiting argument.  Integrating the sharp one-dimensional
log-concave bound then gives

\[
 P_\nu(S)\geq \frac{1}{\sqrt{48V_u(\nu)}},\qquad
 V_u(\nu)=\mathbb E_\nu
 \operatorname {Var}(\langle X,u\rangle\mid P_{u^\perp}X).
 \tag{0.7}
\]

For the initial isotropic law, `V_u(mu)<=1`, including for every fixed `u`.
This proves the desired branch.  The surviving obstruction is only the
early, history-dependent selection of `u`.

All perimeters below are relative weighted `BV` perimeters on the affine
hull.  They are bounded above by exterior Minkowski content, so every lower
bound for `P_mu(S)` is also a lower bound for `mu^+(S)`.

If the affine hull has dimension one, the one-dimensional estimate in
Section 2 directly gives `mu^+(S)>=1/sqrt(48)` in isotropic scale.  The
spectral discussion below is therefore only needed in dimension at least
two.

---

## 1. Posterior and perimeter conventions

Work first with a full-dimensional, compactly supported log-concave density
and a bounded stopped process.  For a bounded predictable symmetric control
`C_t`, stochastic localization may be coupled to a sample `X~mu` and an
observation process by

\[
 dY_t=C_tX\,dt+dB_t.                                   \tag{1.1}
\]

The control is predictable from the observation history.  With the
innovation Brownian motion

\[
 dW_t=dY_t-C_ta_tdt,
 \qquad a_t=\mathbb E[X\mid\mathcal F_t],
\]

the conditional law `mu_t=L(X|F_t)` satisfies

\[
 dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle.          \tag{1.2}
\]

Its likelihood has quadratic precision

\[
 Q_t=\int_0^tC_s^2ds.                                  \tag{1.3}
\]

For

\[
 g_t=\mu_t(S),\qquad
 v_t=\operatorname {Cov}_{\mu_t}(1_S,X),
\]

one has

\[
 dg_t=v_t^TC_t\,dW_t.                                  \tag{1.4}
\]

Consequently `C_t=P_{e_t^perp}`, `e_t` parallel to `v_t`, preserves
`g_t=1/2` pathwise.  At a zero of `v_t`, one must either specify a
predictable line or use the relaxed-control construction; no direction SDE
is asserted at such a zero.

For every deterministic `x`, `p_t(x)/p_0(x)` is a nonnegative likelihood
martingale.  Tonelli on the reduced boundary gives, for bounded stopping
times,

\[
 \mathbb E[P_{\mu_t}(S)]=P_\mu(S).                     \tag{1.5}
\]

More precisely, if `sigma<=tau` are bounded stopping times, then

\[
 \mathbb E[P_{\mu_\tau}(S)\mid\mathcal F_\sigma]
 =P_{\mu_\sigma}(S).                                   \tag{1.5a}
\]

Indeed, conditionally on `F_sigma`, the future observation law is a Markov
kernel of total mass one for every fixed `x`.  Hence the future likelihood
ratio has conditional mean one pointwise in `x`, including at a bounded
stopping time.  For a continuous likelihood ratio `L(x)`, the weighted `BV`
representation reads

\[
 P_{L\mu}(S)=\int_{\partial^*S}L(x)p(x)
 \,d\mathcal H^{n-1}(x).
\]

Conditional Tonelli proves (1.5a).  For a nonsmooth log-concave density the
same identity follows from its precise representative on the relative
reduced boundary, or by monotone interior truncation and `BV` lower
semicontinuity.  The identity is relative to the affine hull and does not
count an artificial boundary of the support.

The corresponding exterior-content statement has the useful direction

\[
 \mu^+(S)\geq\mathbb E[\mu_t^+(S)].                    \tag{1.6}
\]

Equation (1.5) is first literal for a continuous density and a finite-
perimeter representative.  General log-concave laws follow by truncation on
the affine hull, convolution there, stopping, and lower semicontinuity of
weighted `BV` perimeter.  No equality for arbitrary setwise representatives
of exterior content is used.

The posterior-mixture identity is equally important.  For a fixed vector
`u`, or for an `F_tau`-measurable vector after conditioning on `F_tau`,

\[
 \mathbb E\left[
  \operatorname {Var}_{\mu_t}\langle X,u\rangle
  \mid\mathcal F_\tau\right]
 \leq\operatorname {Var}_{\mu_\tau}\langle X,u\rangle.
 \tag{1.7}
\]

This is just conditional total variance.  It is false with the right side
replaced by one when `u` is selected from the history after time zero.

---

## 2. A quantitative fiber lemma

Let `nu` be a log-concave probability on an affine space `E`, let `u` be a
unit vector in its translation space, and write

\[
 x=y+zu,\qquad y=P_{u^\perp}x.
\]

Disintegrate

\[
 \nu(dx)=\eta(dy)\nu_y(dz).
\]

The conditional laws `nu_y` are one-dimensional log-concave probabilities
for `eta`-almost every `y`.  Put

\[
 b(y)=\nu_y(S_y),\quad
 w(y)=\min(b(y),1-b(y)),\quad
 W=\int w\,d\eta,
\]

and

\[
 \sigma_y^2=\operatorname {Var}_{\nu_y}(z),\qquad
 V=\int\sigma_y^2d\eta.                               \tag{2.1}
\]

### Lemma 2.1 (fiber perimeter)

For every finite-perimeter `S`,

\[
 \boxed{\quad
 P_\nu(S)\geq \frac{W^{3/2}}{\sqrt{6V}}.
 \quad}                                                \tag{2.2}
\]

If `b(y)=1/2` almost everywhere, this becomes

\[
 P_\nu(S)\geq \frac{1}{\sqrt{48V}}.                   \tag{2.3}
\]

#### Proof

For a one-dimensional log-concave probability of variance `sigma^2`, its
Cheeger constant is at least

\[
 h\geq \frac{1}{\sqrt{12}\,\sigma}.                   \tag{2.4}
\]

Here is a constant-only verification.  If `rho` is the density and `m` a
median, log-concavity gives

\[
 \|\rho\|_\infty\leq2\rho(m).
\]

The bathtub principle gives

\[
 \sigma^2\geq\frac{1}{12\|\rho\|_\infty^2},
\]

and the monotonicity of the one-dimensional log-concave hazard and reverse
hazard gives `h>=2rho(m)`.  This proves (2.4).  Equivalently, one may use the
exact one-dimensional half-line isoperimetric profile.

Slicing a `BV` set and applying the coarea formula in the `u` direction
gives

\[
 P_\nu(S)\geq
 \int P_{\nu_y}(S_y)d\eta(y)
 \geq \frac{1}{\sqrt{12}}
       \int\frac{w(y)}{\sigma_y}d\eta(y).              \tag{2.5}
\]

Normalize the measure `w d eta` by its mass `W`.  Since
`r -> r^{-1/2}` is convex,

\[
 \int\frac{w}{\sigma}d\eta
 \geq \frac{W^{3/2}}{
          \sqrt{\int w\sigma^2d\eta}}.
\]

As `w<=1/2`, the denominator is at most `sqrt(V/2)`.  This proves
(2.2).  Formula (2.3) follows by setting `W=1/2`.  QED.

Conditional total variance gives

\[
 V\leq\operatorname {Var}_\nu\langle X,u\rangle.       \tag{2.6}
\]

Thus (2.3) is at least `1/sqrt(48)` for every fixed direction of an
isotropic law.

---

## 3. Exact locking forces exact conditional half-mass

The word "locking" must refer to the whole conditional experiment, not to
one realized path of probability zero.

### Theorem 3.1 (Gaussian completeness of a locked line)

Let `nu` be log-concave, `nu(S)=1/2`, and fix a line `Ru`.  Starting from
`nu`, suppose the hard mass-preserving process uses

\[
 C_s=P_{u^\perp}                                      \tag{3.1}
\]

for all future observation outcomes on some nontrivial time interval.  Then

\[
 \nu(S\mid P_{u^\perp}X)=1/2\qquad\nu\text{-a.s.}      \tag{3.2}
\]

Consequently,

\[
 P_\nu(S)\geq\frac{1}{\sqrt{48V_u(\nu)}}.             \tag{3.3}
\]

#### Proof

Let `h=1_S-1/2`, put `y=P_{u^perp}x`, and run (3.1) for a fixed time
`s>0`.  Given the natural parameter `c in u^perp`, the unnormalized signed
posterior mass is

\[
 F_s(c)=\int h(x)
   \exp\left(\langle c,y\rangle-\frac{s}{2}|y|^2\right)
   d\nu(x).                                            \tag{3.4}
\]

The hard control annihilates `v_t`, hence the posterior mass is one-half for
every realized observation.  Under the joint observation representation,

\[
 c=sP_{u^\perp}X+\sqrt{s}\,G,                         \tag{3.5}
\]

where `G` is a standard Gaussian on the span of the projected support.
Therefore `c` has a strictly positive density there, and `F_s(c)=0` for
Lebesgue-almost every `c`.

The Gaussian factor in (3.4) makes `F_s` an entire Laplace transform.  It
vanishes identically.  Uniqueness of the Fourier transform then implies

\[
 (P_{u^\perp})_\#(h\nu)=0.                            \tag{3.6}
\]

This is exactly (3.2).  Apply Lemma 2.1 and (2.6) to obtain (3.3).  QED.

It is enough that `v_s` remain parallel to `u` for every outcome, because
then the hard controller is (3.1).  Notice that no limit `s to infinity` is
needed.  The frequently used limiting explanation is nevertheless correct:
the infinite-time observation reveals `P_{u^perp}X`, while the label remains
independent of the observation.

The theorem does not say that a locked event selected by looking at the
whole future balances every original fiber.  Such a future event biases the
observation outcomes.  That distinction is the source of the factor `p` in
the next theorem.

---

## 4. Quantitative near-locking: the endpoint theorem

We first record the anisotropic endpoint estimate used below.

### Lemma 4.1 (one flat direction)

Let `rho=exp(-V)` be log-concave and suppose, distributionally,

\[
 \nabla^2V\succeq\alpha P_{w^\perp},\qquad\alpha>0.
\]

For every half-mass finite-perimeter set `A`,

\[
 P_\rho(A)\geq \frac{1}{\sqrt{96}}
 \frac{1}{\sqrt{\alpha^{-1}+
       \operatorname {Var}_\rho\langle X,w\rangle}}.  \tag{4.1}
\]

#### Proof

Apply Euclidean needle localization to `1_A-1/2`.  Almost every needle has
half of its conditional mass in `A`.  If its unit direction is `theta` and
its arclength variance is `tau^2`, curvature along the needle is at least
`alpha(1-<theta,w>^2)`.  If
`<theta,w>^2<1/2`, one-dimensional Brascamp--Lieb gives
`tau^2<=2/alpha`.  Otherwise

\[
 \tau^2\leq2\operatorname {Var}_{\mathrm{needle}}
                \langle X,w\rangle.
\]

Thus always

\[
 \tau^2\leq2(\alpha^{-1}+s^2),
\]

where `s^2` is the within-needle variance in direction `w`.  A half-mass set
under a one-dimensional log-concave probability has perimeter at least
`1/(sqrt(48)tau)`.  Integrating, using Jensen, and using
`E s^2<=Var_rho<X,w>` proves (4.1).  All inequalities pass through smooth
approximation in the affine hull.  QED.

### Theorem 4.2 (locked-cap event)

Assume the affine dimension is at least two and fix a bounded stopping time
`tau`.  Let `U` be an `F_tau`-measurable unit
vector, let `A` be an event measurable at time `tau+T`, and put

\[
 p=\mathbb P(A\mid\mathcal F_\tau),\qquad
 K=\operatorname {Var}_{\mu_\tau}\langle X,U\rangle.  \tag{4.2}
\]

Assume that on `A`, for almost every `s in [tau,tau+T]`,

\[
 \sqrt{1-\langle e_s,U\rangle^2}\leq\delta<1/\sqrt{2}. \tag{4.3}
\]

Define `alpha,r,kappa` by (0.5).  Then, conditionally on `F_tau`,

\[
 \boxed{\quad
 P_{\mu_\tau}(S)\geq \frac{1}{\sqrt{96}}
 \frac{p}{\sqrt{
       \alpha^{-1}+\kappa^2
       (\sqrt{K/p}+r/\sqrt{\alpha})^2}}.
 \quad}                                                \tag{4.4}
\]

The right side is zero by convention when `p=0`.

#### Proof

On `A`, put

\[
 M=\int_\tau^{\tau+T}e_se_s^Tds,\qquad Q=TI-M.
\]

Let `w` be a top unit eigenvector of `M`, and write
`q=T-lambda_1(M)`.  Since `Tr(M)=T` and `M` is positive semidefinite,

\[
 q=\sum_{j\geq2}\lambda_j(M).
\]

Together with (4.3), this gives

\[
 q\leq T\delta^2,
 \qquad \lambda_2(M)\leq q.                            \tag{4.5}
\]

Hence

\[
 Q\succeq\alpha P_{w^\perp},
 \qquad \alpha=T(1-\delta^2).                         \tag{4.6}
\]

Moreover,

\[
 U^TQU\leq T\delta^2.
\]

Since the eigenvalues of `Q` on `w^perp` are at least `alpha`,

\[
 \sin\angle(U,w)\leq
 \frac{\delta}{\sqrt{1-\delta^2}}=r.                  \tag{4.7}
\]

The top eigenline is unique under `delta<1/sqrt2`, and can be chosen
measurably.

At the endpoint let

\[
 z_U=\operatorname {Var}_{\mu_{\tau+T}}\langle X,U\rangle,
 \qquad z_w=\operatorname {Var}_{\mu_{\tau+T}}\langle X,w\rangle.
\]

Brascamp--Lieb and (4.6) give variance at most `1/alpha` in every direction
orthogonal to `w`.  Expressing `U` in the `w` decomposition and using the
triangle inequality in `L^2(mu_{tau+T})` gives

\[
 \sqrt{z_w}\leq\kappa
       \left(\sqrt{z_U}+\frac{r}{\sqrt{\alpha}}\right),
 \qquad
 \kappa=(1-r^2)^{-1/2}.                                \tag{4.8}
\]

Lemma 4.1 therefore yields, on `A`,

\[
 P_{\mu_{\tau+T}}(S)\geq\frac{1}{\sqrt{96}}D^{-1/2},
\]

where

\[
 D=\alpha^{-1}+\kappa^2
       (\sqrt{z_U}+r/\sqrt{\alpha})^2.                \tag{4.9}
\]

Conditional total variance gives

\[
 \mathbb E[1_Az_U\mid\mathcal F_\tau]\leq K,
 \qquad
 \mathbb E[1_A\sqrt{z_U}\mid\mathcal F_\tau]
 \leq\sqrt{pK}.                                       \tag{4.10}
\]

It follows that

\[
 \mathbb E[D\mid A,\mathcal F_\tau]
 \leq\alpha^{-1}+\kappa^2
       (\sqrt{K/p}+r/\sqrt{\alpha})^2.                \tag{4.11}
\]

The function `x -> x^{-1/2}` is convex.  Jensen, restriction to `A`, and
the exact conditional perimeter transfer (1.5) prove (4.4).  QED.

If (4.3) holds almost surely after time zero with a fixed deterministic
`u`, then `p=1`, `K=1`, and (4.4) is exactly (0.2).  Thus the claimed
dimension-free locked-line result does not require a quantitative theorem
about convergence of every transverse fiber.

---

## 5. The exact historical estimate

Suppose a line `U_tau` is selected at a stopping time and the process then
enters a cap around it.  Averaging (4.4) back to time zero shows that this
route closes if and only if it supplies a universal lower bound on

\[
 \mathbb E\left[
 \frac{p_\tau}{\sqrt{
       \alpha^{-1}+\kappa^2
       (\sqrt{K_\tau/p_\tau}+r/\sqrt{\alpha})^2}}
 \right].                                             \tag{5.1}
\]

In the exact, infinite-transverse-precision limit this is (0.6).  A simple
sufficient condition is the existence of universal `p_0,K_0,rho_0>0` such
that

\[
 \mathbb P\{p_\tau\geq p_0, K_\tau\leq K_0\}
 \geq\rho_0.                                          \tag{5.2}
\]

For an exact all-outcome lock, Theorem 3.1 improves the relevant variance
to

\[
 V_\tau=\mathbb E_{\mu_\tau}
 \operatorname {Var}(\langle X,U_\tau\rangle
       \mid P_{U_\tau^\perp}X),                       \tag{5.3}
\]

and the historical functional is `E[V_tau^{-1/2}]`.  At `tau=0`, isotropy
gives `V_0<=1` for every line, even if an external randomizer chooses the
line.  After the history has selected the line, isotropy no longer gives
`V_tau<=1` pathwise.

There is one easy favorable case.  Define the effective precision of `Q_tau`
in the selected direction by

\[
 \rho_\tau(U_\tau)=
 \sup\{r\geq0:Q_\tau\succeq rU_\tau U_\tau^T\}.
 \tag{5.4}
\]

When `Q_tau` is positive definite, this is

\[
 \rho_\tau(U_\tau)=
 \big(U_\tau^TQ_\tau^{-1}U_\tau\big)^{-1}.
\]

If `rho_tau(U_tau)>=r_0>0`, then semidefinite Brascamp--Lieb gives
`K_tau<=1/r_0`.  Thus only a line selected before it receives universal
*effective* longitudinal precision can be bad.  The Rayleigh quantity
`U_tau^TQ_tau U_tau` alone is not sufficient unless `U_tau` is an
eigenvector: off-diagonal precision can make the corresponding Schur
complement arbitrarily smaller.

The needed assertion cannot be replaced by a generic theorem about an
adaptively chosen line.  For the isotropic product of one-sided
exponentials, expose all coordinates for time

\[
 t=\frac{1}{2\log n}.
\]

The largest natural parameter is asymptotically at the critical exponential
slope.  Its scalar posterior variance is `Theta(log n)`, although its
accumulated quadratic precision is `Theta(1/log n)`.  Thus a history can
select a high-variance line very early with order-one probability.

This does not defeat (5.1).  For the maximum-tail label, the selected
coordinate does not yet carry order-one posterior label mass and does not
remain persistently protected.  Therefore the genuinely new estimate must
couple all three quantities

\[
 \boxed{\quad
 \text{selected variance }K_\tau,\quad
 \text{effective past precision }\rho_\tau(U_\tau),\quad
 \text{future lock probability }p_\tau.
 \quad}                                                \tag{5.5}
\]

An estimate omitting `p_tau` is false; one omitting `K_tau` does not control
the one remaining fiber; and one assuming a universal positive lower bound
for the effective precision in (5.4) simply assumes away early winner
selection.

---

## 6. Why finite precision alone does not yield a fiber-balance modulus

Exact locking gives exact balance by Theorem 3.1.  There is no uniform
statement of the form

\[
 \text{small label information at finite transverse precision}
 \Longrightarrow
 \mathbb E|\nu(S\mid P_{u^\perp}X)-1/2|\text{ small}   \tag{6.1}
\]

without using boundary variation or the full hard-control geometry.

Here is a two-dimensional Gaussian test.  Let `(Z,Y)` be standard Gaussian,
fix `a>0`, and put

\[
 S_k=\{Z\geq a\sin(kY)\}.                              \tag{6.2}
\]

Symmetry gives `gamma_2(S_k)=1/2`, while the masses on the `Z`-fibers are

\[
 b_k(y)=\Phi(-a\sin(ky)).                              \tag{6.3}
\]

They stay a fixed positive `L^1` distance from one-half as `k to infinity`.
Observe `Y` through the finite-SNR channel

\[
 H=\sqrt{T}\,Y+G.                                    \tag{6.4}
\]

Conditionally on `H`, the variance of `Y` is `1/(1+T)`.  The nonconstant
Fourier modes of the periodic function in (6.3) are therefore multiplied by

\[
 \exp\left(-\frac{j^2k^2}{2(1+T)}\right).
\]

For every fixed `T`,

\[
 \mathbb E\left|
  \mathbb P(S_k\mid H)-\frac{1}{2}\right|\longrightarrow0,
 \qquad k\longrightarrow\infty,                      \tag{6.5}
\]

although the true fiber imbalance does not tend to zero.  The Gaussian
perimeter of (6.2) is

\[
 \int\varphi(y)\varphi(a\sin ky)
       \sqrt{1+a^2k^2\cos^2ky}\,dy=\Theta_a(k),        \tag{6.6}
\]

so the KLS conclusion is safe: unresolved oscillations pay boundary.

This example is not a counterexample to Theorem 4.2 and is not asserted to
follow the exact hard feedback.  It rules out only the proposed generic
Bayesian inference from finite precision to fiber balance.  A proof that
insists on approximate fibers must quantitatively charge the unresolved
oscillations to `P_mu(S)`; in low-frequency configurations that charge is a
lower-dimensional Cheeger problem.  The endpoint proof avoids this
circularity.

---

## 7. Canonical tests

### 7.1 Product and coordinate halfspaces

Let `mu=nu tensor eta`, let `u` be the `nu` coordinate, and let
`S={z>=m_nu}`.  Transverse tilts change only `eta`, so `v_t` is always
parallel to `u`.  Every `u`-fiber has mass one-half.  Theorem 3.1 applies
with

\[
 V_u=\operatorname {Var}_\nu Z=1.
\]

For a centered rate-one exponential the actual perimeter is `1/2`.  For
the isotropic interval `[-sqrt3,sqrt3]` it is `1/(2sqrt3)`.  Gaussian
halfspaces and coordinate half-cubes are equally exact locked models.

### 7.2 A random protected line need not cause a cap-counting loss

For a Gaussian product take the parity cut

\[
 S=\{\operatorname {sign}X_1\operatorname {sign}X_2=1\}.
\]

At `v_0=0`, an external convention may protect either coordinate.  Once
one is protected and the other observed, `v_t` lies on the protected line.
Conditionally on the observed coordinate, the cut is a halfline of Gaussian
mass one-half on the protected fiber.  Thus `V=1` whichever line was chosen.
Conditioning at the line-selection time, rather than covering the sphere by
caps and summing `p_i^{3/2}`, gives the correct constant with no factor
depending on the number of possible lines.

### 7.3 The cube away from a product direction

For a diagonal halfspace of the cube, a line parallel to its normal meets a
typical transverse fiber in an interval whose two endpoints are not
symmetric about the cutting plane.  The conditional halfspace mass varies
with the base point.  Theorem 3.1 therefore shows that the hard process
cannot remain exactly locked to that diagonal for every transverse tilt.
The coordinate halfspace is the product exception and is covered by
Section 7.1.

### 7.4 The simplex

In the standard simplex

\[
 \{x_i\geq0,\ \sum_{i=1}^n x_i\leq1\},
\]

consider `S={x_1>=m}`.  Given
`y=(x_2,...,x_n)`, the conditional law of `x_1` is uniform on
`[0,1-sum_{i>=2}y_i]`, and

\[
 \mathbb P(S\mid y)=
 \frac{(1-\sum_{i\geq2}y_i-m)_+}
      {1-\sum_{i\geq2}y_i}.                            \tag{7.1}
\]

This is not constant.  To check the Euclidean fibers after isotropic
whitening, write the covariance of the standard simplex as a positive
multiple of `(n+1)I-11^T`.  The normal of the whitened cut pulls back to
parallel lines in direction

\[
 ((n+1)I-\mathbf1\mathbf1^T)e_1
 =ne_1-\sum_{i\geq2}e_i.
\]

Their intersection lengths with the simplex, and the fractions lying above
`x_1=m`, vary with the parallel-line label.  In dimension two these are
lines of direction `(2,-1)` in the triangle; fibers near different vertices
already give different fractions.  Thus the whitened conditional masses
are not constant either.  Although symmetry aligns `v_0` with the whitened
vertex axis, generic transverse tilts must rotate it.  Hence the simplex
does not satisfy the exact locked-line hypothesis except in a genuinely
fiber-balanced cut.

### 7.5 Maximum cuts of one-sided exponential products

Let `Z_i` be independent rate-one exponentials and choose `L` by

\[
 (1-e^{-L})^n=1/2,
 \qquad S=\{\max_iZ_i\geq L\}.
\]

Then

\[
 P_\mu(S)=ne^{-L}(1-e^{-L})^{n-1}\longrightarrow
 \frac{\log 2}{2}.                                    \tag{7.2}
\]

Initially the protected line is the diagonal.  However, if `e_t` is its
unit direction, the exact initial angular quadratic-variation rate is

\[
 \frac{d}{dt}\operatorname {Tr}[e]_t\bigg|_{t=0}
 =(1+o(1))(\log n)^2.                                  \tag{7.3}
\]

Moreover, conditional masses on diagonal fibers depend on the transverse
base point.  Thus exact diagonal locking is impossible, and even a
near-lock argument must pay for a very fast initial selection layer.  This
is exactly the situation measured by `p_tau` and `K_tau` in (5.1).  The
perimeter in (7.2) is already universal, so it is a stress test rather than
a counterexample.

---

## 8. Formal status

The following statements are proved without KLS:

1. exact all-outcome locking to a fixed line forces exact conditional
   half-mass on almost every parallel fiber;
2. balanced log-concave fibers give the explicit bound (0.7);
3. projective cap-locking for time `T` gives the quantitative endpoint
   bounds (0.2) and (0.4), including all event-probability factors;
4. perimeter transfer has the correct direction, and the selected-line
   covariance is used only after conditioning at the stopping time;
5. the unique remaining historical functional is (5.1), asymptotically
   (0.6).

The route is therefore complete once one proves a persistence-weighted
early-selection theorem of the form (5.1).  Neither isotropy, a deterministic
direction variance average, nor Brascamp--Lieb at small longitudinal
precision proves that theorem.  Assuming it would be the same adaptive
exceptional-variance step isolated by the hard-localization reduction.
