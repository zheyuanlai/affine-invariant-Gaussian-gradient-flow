# Inertial protected-line localization

## Executive conclusion

This note analyzes a bounded-speed protected direction.  Let

\[
 g_t=\mu_t(S),\qquad G_t=g_t(1-g_t),\qquad
 v_t=\operatorname {Cov}_{\mu_t}(1_S,X),
\]

let `u_t` be a predictable unit vector of finite variation, and use

\[
 C_t=P_{u_t^\perp}=I-u_tu_t^T.                         \tag{0.1}
\]

The exact information cost of lagging behind the instantaneous label
direction is

\[
 \mathcal E_T=\int_0^T
 \frac{|P_{u_t^\perp}v_t|^2}{G_t}\,dt.                \tag{0.2}
\]

Indeed,

\[
 d\langle g\rangle_t=|P_{u_t^\perp}v_t|^2dt,
 \qquad
 \mathbb E\mathcal E_T
 =2\{H(g_0)-\mathbb E H(g_T)\}.                       \tag{0.3}
\]

A concrete bounded-speed tracking law is

\[
 Z_t={P_{u_t^\perp}v_t\over\sqrt{G_t}},\qquad
 \dot u_t={Z_t\over\sqrt{1+|Z_t|^2/\Lambda^2}},       \tag{0.4}
\]

used while `g_t` stays in a fixed compact subinterval of `(0,1)`.  It is
globally regular at `v=0`, has speed below `Lambda`, and obeys the exact kinetic
identity

\[
 { |P_{u_t^\perp}v_t|^2\over G_t}
 ={|\dot u_t|^2\over1-|\dot u_t|^2/\Lambda^2}
 \geq|\dot u_t|^2.                                    \tag{0.5}
\]

Let

\[
 Q_T=\int_0^TP_{u_t^\perp}dt,
 \qquad q_T=\lambda_{\min}(Q_T),
\]

and, when `q_T<T/2`, let `w_T` be its exceptional eigenline.  A deterministic
trace inequality gives the promised inertia tradeoff:

\[
 \boxed{\quad
 d_{\mathbb{RP}}(u_0,w_T)^2
 \leq \frac{\pi^2q_T}{2T}+T\mathcal E_T.
 \quad}                                                \tag{0.6}
\]

Thus an exceptional endpoint line can move far from the deterministic
initial line only by paying either Schur precision `q_T` in that line or
binary-entropy clock.  Since `w_T` is an eigenvector, its effective
Schur-complement precision is exactly `q_T`; no invalid replacement by the
Rayleigh quotient is made.

On an event where `g_T` remains central, `q_T<=epsilon T`, and the entropy
clock is at most `B`, with

\[
 \theta_*^2=\frac{\pi^2\epsilon}{2}+TB<\frac{\pi^2}{4}, \tag{0.7}
\]

the endpoint selected variance is controlled by the posterior variance in
the deterministic direction `u_0` and the transverse curvature.  Conditional
total variance then gives a completely explicit, dimension-free endpoint
perimeter lower bound.  The complementary full-curvature event
`q_T>=epsilon T` is also good.

What is not proved is that either good event has universal probability.
The information identity alone does not do this.  A burst with normalized
binary correlation `K`, transverse to `u_t`, can spend entropy clock `B` in
time `B/K`, before a speed-`Lambda` line rotates appreciably.  The exact clock,
kinetic, precision, and endpoint inequalities permit such bursts.  Excluding
them requires a new assertion:

> a large label-correlated posterior variance cannot arise, in a direction
> not already effectively exposed, before it loses persistent correlation
> with the fixed set `S`.

This is the same early calibration/persistence obstruction found for the
hard protected-line process, now in an especially sharp form.  No bound on
this obstruction is proved here, so the inertial driver does not currently
close KLS.

The mandatory product-exponential maximum test does **not** realize the
naive bad burst in the exactly computable initial state or in the
full-exposure large-deviation proxy.  Its instantaneous angular rate is
`(log n)^2`, but its normalized binary correlation is only
`Theta((log n)^2/n)`.  At the subconstant proxy time at which a high-variance
coordinate is statistically selected, that coordinate still carries
vanishing posterior mass of the fixed maximum-tail label.  Transferring the
proxy uniformly to the inertial feedback is not proved; the calculation is
a mandatory cancellation test, not a good-event theorem.

---

## 1. Stopped construction and well-posedness

Work first with a smooth log-concave density on a compact convex support.
Stop before the natural parameters or the first four moments leave a fixed
compact set and, for `0<delta<1/2`, before

\[
 \tau_\delta=\inf\{t:g_t\notin(\delta,1-\delta)\}.      \tag{1.1}
\]

For a predictable symmetric control `C_t`, stochastic localization is

\[
 dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle,
 \qquad a_t=\mathbb E_tX.                              \tag{1.2}
\]

Equivalently,

\[
 p_t(x)={1\over Z_t}
 \exp\left(c_t\cdot x-{1\over2}x^TQ_tx\right)p_0(x),
 \qquad dQ_t=C_t^2dt.                                  \tag{1.3}
\]

Put

\[
 Y=X-a_t,\quad h=1_S-g_t,\quad
 A_t=\mathbb E_tYY^T,\quad D_t=\mathbb E_t[hYY^T],
\]

and

\[
 \mathcal T_t(z)=\mathbb E_t[YY^T\langle Y,z\rangle].
\]

For `C_t=P_{u_t^perp}`, direct Itô differentiation gives

\[
 \boxed{\begin{aligned}
 da_t&=A_tC_t\,dW_t,\\
 dg_t&=v_t^TC_t\,dW_t,\\
 dv_t&=D_tC_t\,dW_t-A_tC_tv_t\,dt,\\
 dA_t&=\mathcal T_t(C_t\,dW_t)-A_tC_tA_t\,dt.
 \end{aligned}}                                       \tag{1.4}
\]

The drift in `dv_t` is essential.  It is `-A_tC_t^2v_t`; here `C_t^2=C_t`.
It disappears only for the hard instantaneous choice `C_tv_t=0`.

On the stopped band, define the coupled tracking law

\[
 Z(c,Q,u)={P_{u^\perp}v(c,Q)\over\sqrt{g(c,Q)(1-g(c,Q))}},
 \qquad F_\Lambda(Z)={Z\over\sqrt{1+|Z|^2/\Lambda^2}}, \tag{1.5}
\]

and set

\[
 du_t=F_\Lambda(Z_t)dt.                                \tag{1.6}
\]

Take `u_0` to be the deterministic line of `v_0` when `v_0` is nonzero,
and use any fixed deterministic convention when `v_0=0`.  Isotropy then
always gives `Var_mu<X,u_0>=1`.

Because `Z_t` is orthogonal to `u_t`, (1.6) preserves `|u_t|=1`.
On compact support, the moment maps `(c,Q) -> (g,v,A,D)` are smooth.  The
denominator in (1.5) is bounded away from zero, `F_Lambda` is globally smooth,
and `u -> P_{u^perp}` is smooth on the sphere.  Hence the coupled finite-
dimensional SDE/ODE is locally Lipschitz and has a unique strong solution up
to the stated stops.  No division by `|v|` occurs, so `v=0` causes no
singularity.

For a general log-concave law, truncate on increasing compact convex sets in
its affine hull, smooth there, apply the stopped identities, and take a weak
subsequential limit.  All controls are bounded, the entropy identities are
localized on `g in [delta,1-delta]`, and the likelihood representation gives
uniform integrability at bounded parameter stops.  The report makes no
claim of a global unstopped strong solution after `g` approaches zero or
one.  One may set `C=0` after `tau_delta`; then all identities remain valid
with active time `t wedge tau_delta`.

The law (1.6) tracks an oriented lift of `v`.  Since the actual control
depends only on `u_tu_t^T`, changing the lift has no effect when `u_t` is
already parallel or antiparallel to `v_t`.  A fully projective smooth
variant is obtained by replacing `Z` with

\[
 Z_\eta={\langle u,v\rangle P_{u^\perp}v
 \over \sqrt G\sqrt{\langle u,v\rangle^2+\eta G}},
 \qquad\eta>0.                                        \tag{1.7}
\]

Then `|Z_eta|^2<=|P_{u^perp}v|^2/G`, so every tradeoff below remains valid
with inequality in place of equality in (0.5).  This removes any choice of
orientation at the cost of slowing at an exactly orthogonal line.

---

## 2. Exact mass and entropy clocks

From (1.4),

\[
 d\langle g\rangle_t
 =|P_{u_t^\perp}v_t|^2dt=G_t|Z_t|^2dt.                \tag{2.1}
\]

For binary entropy

\[
 H(s)=-s\log s-(1-s)\log(1-s),
\]

Itô's formula gives

\[
 dH(g_t)=H'(g_t)dg_t-{1\over2}|Z_t|^2dt.              \tag{2.2}
\]

Thus, up to bounded stops,

\[
 \boxed{\quad
 \mathbb E\mathcal E_t
 =2\{H(g_0)-\mathbb E H(g_t)\},
 \qquad
 \mathcal E_t=\int_0^t|Z_s|^2ds.
 \quad}                                                \tag{2.3}
\]

Starting from `g_0=1/2`, (2.3) says that the expected entropy clock is at
most `2 log 2`.  If the process is stopped
when it first exits `[delta,1-delta]`, the sharper bound is

\[
 \mathbb E\mathcal E_{t\wedge\tau_\delta}
 \leq2\{\log 2-H(\delta)\}.                           \tag{2.4}
\]

For `delta=1/4`, the right side is approximately `0.2616`.

The tracking speed satisfies

\[
 |\dot u_t|^2={|Z_t|^2\over1+|Z_t|^2/\Lambda^2},
\]

which is equivalent to (0.5).  In particular,

\[
 \int_0^T|\dot u_t|^2dt\leq\mathcal E_T.              \tag{2.5}
\]

This is the exact sense in which changing the protected direction spends
binary information.

There is a useful intrinsic-time identity.  In entropy-clock time `s`,

\[
 dg_s=\sqrt{g_s(1-g_s)}\,dB_s.                         \tag{2.6}
\]

Writing `G_s=g_s(1-g_s)`, Itô gives

\[
 d(e^sG_s)=e^s(1-2g_s)\sqrt{G_s}\,dB_s.               \tag{2.7}
\]

Consequently, for every bounded entropy-time stopping time `sigma<=B`,

\[
 \mathbb E G_\sigma
 \geq e^{-B}\mathbb E[e^\sigma G_\sigma]
 =\frac14e^{-B}.                                      \tag{2.8}
\]

This leaves universal expected central mass if one can cap the entropy
clock pathwise.  A pure projected controller does not provide such a cap
over a fixed physical time; stopping the whole process when the cap is hit
may leave arbitrarily little accumulated curvature.

---

## 3. Precision spectrum and the inertial trace inequality

On a path active for a deterministic time `T`, put

\[
 M_T=\int_0^Tu_tu_t^Tdt,
 \qquad Q_T=TI-M_T.                                    \tag{3.1}
\]

Both matrices are positive semidefinite and `Tr M_T=T`.  Let `w` be a top
unit eigenvector of `M_T` and put

\[
 q=T-\lambda_1(M_T)=\lambda_{\min}(Q_T).              \tag{3.2}
\]

Since

\[
 q=\sum_{j\geq2}\lambda_j(M_T),
\]

one has

\[
 Q_T\succeq(T-q)P_{w^\perp}                           \tag{3.3}
\]

whenever `q<T`.  Also

\[
 q=w^TQ_Tw=
 \int_0^T\{1-\langle u_t,w\rangle^2\}dt.             \tag{3.4}
\]

Thus `q` is both the smallest eigenvalue and the effective Schur precision
of the exceptional eigenline.

Let

\[
 \phi_t=d_{\mathbb{RP}}(u_t,w)
 =\arccos|\langle u_t,w\rangle|\in[0,\pi/2].          \tag{3.5}
\]

Its metric derivative is at most `|dot u_t|`.  For every `t`,

\[
 \phi_0\leq\phi_t+\int_0^t|\dot u_s|ds.
\]

Square, average in `t`, use Cauchy--Schwarz, and then use
`phi<=(pi/2)sin phi`.  This yields

\[
\begin{aligned}
 \phi_0^2
 &\leq {2\over T}\int_0^T\phi_t^2dt
 +{2\over T}\int_0^Tt\int_0^t|\dot u_s|^2ds\,dt\\
 &\leq {\pi^2\over2T}q
 +T\int_0^T|\dot u_s|^2ds\\
 &\leq {\pi^2\over2T}q+T\mathcal E_T.                \tag{3.6}
\end{aligned}
\]

This proves (0.6).  The constants are independent of dimension and of the
speed cap `Lambda`.  The cap matters separately: a rotation through
projective angle `theta` takes at least `theta/Lambda` physical time.

For comparison, if `u_t` follows the shortest planar geodesic from a line at
angle `theta` to `w` at maximal speed `Lambda` and then stays at `w`, the deposited
Rayleigh precision is exactly

\[
 A(\theta,\Lambda)={1\over \Lambda}
 \int_0^\theta\sin^2s\,ds
 ={1\over \Lambda}\left({\theta\over2}-{\sin2\theta\over4}\right).
 \tag{3.7}
\]

The off-diagonal block has norm `sin^2(theta)/(2Lambda)`.  If the total active
time in this two-dimensional rotation plane is `R` (so the other diagonal
entry is `R-A(theta,Lambda)`), the effective Schur precision is

\[
 A(\theta,\Lambda)-
 {\sin^4\theta\over4\Lambda^2(R-A(\theta,\Lambda))}. \tag{3.8}
\]

For small `theta` and fixed positive `R`, this is
`theta^3/(3Lambda)+O(theta^4/(Lambda^2 R))`.  Thus bounded-speed rotation deposits a
genuine cubic seed, but a seed that vanishes with the angle and can be
formed after the label has already spent its information clock.

---

## 4. Endpoint perimeter and selected variance

For a finite-perimeter fixed set, posterior perimeter transfers exactly:

\[
 P_\mu(S)=\mathbb E P_{\mu_T}(S).                     \tag{4.0}
\]

Indeed the endpoint likelihood is continuous in `x` and has pointwise mean
one; Tonelli on the relative reduced boundary proves (4.0).  Exterior
Minkowski content dominates the relative weighted `BV` perimeter, so a
lower bound obtained from (4.0) is also a lower bound for `mu^+(S)`.

The following arbitrary-mass version of the one-flat-direction lemma is
useful.

### Lemma 4.1

Let `rho=exp(-V)` be log-concave and suppose

\[
 \nabla^2V\succeq\alpha P_{w^\perp}
\]

distributionally.  For every finite-perimeter `A`,

\[
 P_\rho(A)\geq {1\over\sqrt{24}}
 {\min(\rho(A),1-\rho(A))
  \over\sqrt{\alpha^{-1}+
       \operatorname {Var}_\rho\langle X,w\rangle}}.  \tag{4.1}
\]

#### Proof

Localize `1_A-rho(A)` into one-dimensional needles.  Every needle gives
`A` the same mass `rho(A)`.  A one-dimensional log-concave law of variance
`tau^2` has Cheeger constant at least `1/(sqrt(12)tau)`.  Curvature and
Brascamp--Lieb give

\[
 \tau^2\leq2\{\alpha^{-1}+
 \operatorname {Var}_{\rm needle}\langle X,w\rangle\}.
\]

Integrate and apply Jensen and conditional total variance.  QED.

At time `T`, on the exceptional branch `q<T/2`, (3.3) and Lemma 4.1 give

\[
 P_{\mu_T}(S)\geq {m_T\over\sqrt{24}}
 {1\over\sqrt{(T-q)^{-1}+z_w}},                       \tag{4.2}
\]

where

\[
 m_T=\min(g_T,1-g_T),\qquad
 z_w=\operatorname {Var}_{\mu_T}\langle X,w\rangle.
\]

Suppose on an event `A` that

\[
 q\leq\epsilon T,\qquad
 \mathcal E_T\leq B,\qquad
 m_T\geq\delta,
\]

and put

\[
 \theta_*=\sqrt{\frac{\pi^2\epsilon}{2}+TB}<\frac\pi2,
 \quad r=\sin\theta_*,\quad\kappa=(\cos\theta_*)^{-1},
 \quad\alpha=(1-\epsilon)T.                           \tag{4.3}
\]

Let

\[
 z_0=\operatorname {Var}_{\mu_T}\langle X,u_0\rangle.
\]

Brascamp--Lieb in `w^perp` and (3.6) give

\[
 \sqrt{z_w}\leq\kappa
 \left(\sqrt{z_0}+{r\over\sqrt\alpha}\right).        \tag{4.4}
\]

If `u_0` is deterministic and the starting law is isotropic, posterior
mixture gives

\[
 \mathbb Ez_0\leq1.                                   \tag{4.5}
\]

Writing `p=P(A)`, conditional Jensen and perimeter transfer therefore give

\[
\boxed{\quad
 P_\mu(S)\geq {\delta\over\sqrt{24}}
 {p\over\sqrt{
 \alpha^{-1}+\kappa^2
 (p^{-1/2}+r/\sqrt\alpha)^2}}.
\quad}                                                 \tag{4.6}
\]

This is dimension free whenever `epsilon,T,B,delta,p` are universal and
`theta_*<pi/2`.

On the complementary spectral branch `q>=epsilon T`, the endpoint is
`epsilon T`-strongly log-concave.  The standard strongly-log-concave
Cheeger bound gives

\[
 P_{\mu_T}(S)\geq c\sqrt{\epsilon T}\,m_T.             \tag{4.7}
\]

Equations (4.6)--(4.7) show exactly what a successful probability estimate
would buy.  They do not prove that `m_T` is central on either branch.

---

## 5. Why the entropy identity does not close the probability estimate

The mean identity

\[
 \mathbb E\mathcal E_T
 =2\{\log 2-\mathbb E H(g_T)\}\leq2\log 2           \tag{5.1}
\]

does not provide a good event satisfying (4.3).  Large entropy clock is
exactly the mechanism that makes `m_T` small.  Markov's inequality can bound
the probability of a large clock only at a threshold comparable to
`2 log 2`; the trace estimate then leaves no useful joint control when the
same paths are selected by small endpoint mass.

Stopping at the first exit from `[delta,1-delta]` keeps `m=delta`, but it
does not solve the issue.  The active time, and hence every nonzero
eigenvalue of `Q`, can then be arbitrarily small.  Continuing with `C=0`
adds no curvature.  Switching instantaneously to the hard line `e_t`
preserves mass but reintroduces a path-selected unexposed direction and
abandons the inertial premise.

There is a sharp statewise model showing that the clocks and precision
geometry alone cannot rule out the failure.  Fix large `K` and a clock
budget `B`.  During an interval of length

\[
 \Delta t={B\over K},                                  \tag{5.2}
\]

take a posterior state with

\[
 \frac{|v|^2}{G}=K,\qquad e\perp u.                   \tag{5.3}
\]

Then the entropy clock spent is exactly `B`.  Any speed-`Lambda` controller
rotates by at most

\[
 {\Lambda B\over K},                                  \tag{5.4}
\]

and the precision deposited in the high-correlation direction is only
`B/K`.  This is the sharp scaling allowed by Brascamp--Lieb: a log-concave
one-dimensional posterior can have variance `Theta(K)` at quadratic
precision `Theta(1/K)` and a median halfspace has normalized binary
correlation `Theta(K)`.

In entropy time the label follows (2.6).  At deterministic clock `B`, (2.7)
gives

\[
 \mathbb EG_B={1\over4}e^{-B},
 \qquad
 \mathbb E\min(g_B,1-g_B)\leq{1\over2}e^{-B}.          \tag{5.5}
\]

After this burst one may keep `u` fixed and arrange zero further label
clock while accumulating transverse precision.  The exceptional line then
remains deterministic, but its endpoint isoperimetric contribution is
multiplied by the arbitrarily small quantity in (5.5).  Taking, for example,
`B=log K`, the burst duration and line rotation both tend to zero while the
endpoint mass contribution tends to zero.

This is a countermodel to any proof using only (0.3), (0.5), the spectrum of
`Q`, and Brascamp--Lieb.  It is not asserted to be a complete trajectory from
an initially isotropic log-concave law.  Realizing or excluding the burst
from such an initial law is precisely the missing historical theorem.  At
time zero isotropy gives `|v_0|^2/G_0<=1`; the issue is whether adaptive
tilting can create (5.3) before effective exposure while retaining
correlation with the *fixed* label.

The obstruction is independent of the particular speed formula (0.4).
Freeze, as a scalar comparison block, a direction `e` and a normalized
binary correlation `K`.  For any protected-line trajectory define the
precision deposited in `e` during the block by

\[
 R_e=\int |P_{u_t^\perp}e|^2dt.
\]

On this block the two exact clocks satisfy

\[
 \boxed{\quad \mathcal E=KR_e.\quad}                  \tag{5.6}
\]

Thus changing the tracking speed merely chooses where to sit on the same
information--precision curve.  A one-dimensional log-concave posterior of
variance `Theta(K)` can saturate the scale

\[
 \operatorname {Var}_{\rm end}(e)
 \asymp {1\over K^{-1}+R_e}.                           \tag{5.7}
\]

In the intrinsic label diffusion, a deterministic clock `KR_e` leaves
expected central mass at most a constant times `exp(-KR_e)`.  The strongest
endpoint scale obtainable from these two scalar quantities is consequently

\[
 \sup_{R_e\geq0}
 e^{-KR_e}\sqrt{K^{-1}+R_e}
 ={1\over\sqrt K}
 \sup_{x\geq0}e^{-x}\sqrt{1+x}
 ={1\over\sqrt K}.                                    \tag{5.8}
\]

For a very fast tracker, `R_e=O(1/K)` preserves label mass but leaves
variance `Theta(K)`.  For a slow tracker, `R_e` is larger but the label
polarizes before the curvature is useful.  Therefore no choice of local
tracking speed, bounded or unbounded, closes the endpoint by a scalar
statewise tradeoff.  A successful inertial proof must show that the frozen
high-`K` calibrated block cannot be created from the isotropic history.

---

## 6. Deterministic rotation deposits Schur precision, but not enough by itself

It is tempting to use the speed bound alone.  If the eventual exceptional
line is at angle `theta` from `u_0`, (3.7) shows that a direct rotation pays
`Theta(theta^3/Lambda)` effective precision.  Brascamp--Lieb then bounds the
eventual variance by `O(Lambda/theta^3)`.  This is useful for a fixed order-one
angle.

It does not control the information spent before the rotation.  In the
burst model, the label clock is spent in time `B/K`, while the line moves
only `Lambda B/K`; the later cubic rotation seed arrives after `g` has nearly
polarized.  Perimeter transfer weights the endpoint Cheeger inequality by
`min(g_T,1-g_T)`, so late curvature cannot restore the lost factor.

Conversely, making `Lambda` arbitrarily large removes inertia and tends back to
the hard controller `u=e`.  The endpoint then again contains a direction
selected from the same history as its posterior variance.  Thus no fixed
choice of `Lambda` eliminates both obstructions using the identities proved here.

A closing inertial theorem would have to prove a persistence-weighted burst
exclusion, for example universal constants `K_0,B_0,p_0` such that with
probability at least `p_0`, before the label spends clock `B_0`, either

1. every direction with normalized binary correlation above `K_0` has
   effective precision at least `c/K_0`, or
2. the eventual exceptional line remains in a fixed projective cap around
   `u_0` and the endpoint label mass stays central.

Neither alternative follows from one-dimensional log-concavity or from an
operator-covariance bound, and no such bound is assumed here.

---

## 7. Mandatory product-exponential maximum test

Let `Z_1,...,Z_n` be independent rate-one exponentials, centered to make the
product isotropic.  Choose `L=L_n` by

\[
 (1-e^{-L})^n={1\over2},
 \qquad S=\{\max_iZ_i\geq L\}.                         \tag{7.1}
\]

Put `q=e^{-L}=1-2^{-1/n}`.  Initially

\[
 (v_0)_i={Lq\over2(1-q)},
 \qquad
 |v_0|=\frac{\sqrt{n}Lq}{2(1-q)}.                     \tag{7.2}
\]

Since `G_0=1/4` and `nq` tends to `log 2`,

\[
 k_0=\frac{|v_0|^2}{G_0}
 =(1+o(1))\frac{(\log 2)^2L^2}{n}\longrightarrow0.   \tag{7.3}
\]

The initial label direction is the diagonal.  In contrast, its exact
angular quadratic-variation rate under the hard process is

\[
 {d\over dt}\operatorname {Tr}[e]_t\bigg|_{t=0}
 =(1+o(1))L^2.                                         \tag{7.4}
\]

Thus the direction can begin moving extremely fast while an inertial line
lags at negligible entropy cost.  This is favorable, not fatal: the
protected line remains close to the deterministic diagonal, and the fixed-
direction posterior mixture bound remains available.

The more delicate test is a possible early winning coordinate.  In the
full-exposure proxy at time

\[
 t={\kappa\over L},
\]

the largest scalar natural parameter is

\[
 c_*(\kappa)=
 \begin{cases}
 \sqrt{2\kappa}+o_{\mathbb P}(1),&\kappa\leq1/2,\\
 \kappa+1/2+o_{\mathbb P}(1),&\kappa\geq1/2.
 \end{cases}                                          \tag{7.5}
\]

The posterior tail mass above the fixed threshold `L` has exponential rate

\[
 J_\kappa(c)=
 \begin{cases}
 1+\kappa/2-c,&c\leq1,\\
 (1+\kappa-c)^2/(2\kappa),&1<c<1+\kappa,\\
 0,&c\geq1+\kappa.
 \end{cases}                                          \tag{7.6}
\]

At the winning value,

\[
 J_\kappa(c_*)=
 \begin{cases}
 1+\kappa/2-\sqrt{2\kappa},&\kappa\leq1/2,\\
 1/(8\kappa),&\kappa\geq1/2.
 \end{cases}                                          \tag{7.7}
\]

This is positive for every fixed `kappa`.  Hence even the winning coordinate
has posterior label-tail probability `n^{-J_kappa+o_P(1)}`.  At
`kappa=1/2` its scalar posterior variance is `Theta(L)`, but it does not have
normalized binary correlation `Theta(L)` with the fixed set `S`.  It is not
the burst (5.3).

In the proxy, making a coordinate carry order-one posterior mass of `S`
requires `kappa` of order `L`, so the physical time is order one.  An
inertial line which remained near the diagonal would then have accumulated
order-one precision in every coordinate except for its `1/sqrt n` diagonal
component; that component is controlled in expectation by the deterministic
direction mixture identity.

The full-exposure computation does not by itself control the evolving
off-diagonal quadratic tilt or the feedback line of the actual inertial
process.  It therefore rules out the naive early-winner story for this
model, but it is not a proof of the good event in (4.3).

Its actual perimeter is already universal:

\[
 P_\mu(S)=nq(1-q)^{n-1}\longrightarrow\frac{\log 2}{2}. \tag{7.8}
\]

The lesson is exact: fast angular selection is harmless when binary
correlation is small; high posterior variance is harmless when the fixed
label is not calibrated to it; and order-one calibrated label correlation
appears only after order-one exposure in this model.

---

## 8. Other canonical tests

### Product and cube coordinate halfspaces

For a product halfspace, `v_t` stays on its deterministic coordinate line
under transverse tilts.  Taking `u_0` on that line gives `Z_t=0`, so the
inertial line does not move, `g_t=1/2`, and

\[
 Q_T=TP_{u_0^\perp}.
\]

The exact balanced-fiber theorem gives a universal perimeter.  This includes
the centered one-sided exponential coordinate cut and a coordinate
half-cube.

### Gaussian halfspaces and radial origin cuts

The same exact locking occurs for a Gaussian halfspace and, by reflection
symmetry on every parallel fiber, for an origin halfspace under a
rotation-invariant log-concave law.  The inertia mechanism is inactive and
loses no constant.

### Simplex and nonproduct cube cuts

For a vertex-axis simplex cut or a diagonal cube cut, conditional masses on
parallel fibers are not constant.  Transverse tilts therefore rotate `v_t`.
The trace inequality (3.6) applies, but no model-specific estimate proved
here bounds the joint event of entropy loss, rotation, and selected
variance.  These examples neither refute nor close the driver.

### Dimension one

There is no transverse direction.  The one-dimensional log-concave Cheeger
bound directly gives the conjectured constant, so the inertial construction
is unnecessary.

---

## 9. Formal status

The following are proved with universal constants and without covariance-op
assumptions:

1. strong stopped well-posedness of the regularized SDE/ODE;
2. the exact SDEs (1.4), entropy clock (2.3), and kinetic identity (0.5);
3. the exact spectrum (3.3)--(3.4) and the dimension-free inertial trace
   inequality (0.6);
4. the endpoint selected-variance and perimeter estimate (4.6);
5. a sharp clock/precision burst countermodel showing that these identities
   alone cannot close the endpoint probability estimate;
6. the product-exponential maximum cancellation, including its distinct
   angular and binary-correlation scales.

No universal probability bound for the good event in (4.3) is proved.  In
particular, this note does not establish KLS and does not assume the missing
burst-exclusion theorem.  Any continuation must use the fixed-label
calibration or boundary geometry to rule out the statewise burst; a scalar
tracking law, entropy identity, and accumulated-precision spectrum are not
enough.
