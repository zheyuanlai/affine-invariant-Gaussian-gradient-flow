# Soft mass localization: the exact information--curvature tradeoff

## Executive conclusion

Let

\[
 g_t=\mu_t(S),\qquad
 v_t=\operatorname {Cov}_{\mu_t}(\mathbf 1_S,X),\qquad
 e_t={v_t\over |v_t|},
\]

and, until `g_t` exits `[delta,1-delta]`, use the full-rank control

\[
 C_t=P_{e_t^\perp}+\alpha_t e_te_t^T,
 \qquad
 C_t^2=P_{e_t^\perp}+\alpha_t^2e_te_t^T.              \tag{0.1}
\]

There is an exact tradeoff, not merely an inequality.  Put

\[
 G_t=g_t(1-g_t),\qquad
 r_t={|v_t|^2\over G_t},\qquad
 \kappa_t={|v_t|\over G_t}
 =\left|\mathbb E_t[X\mid S]-\mathbb E_t[X\mid S^c]\right|,
\]

and let

\[
 R_t=\int_0^t\alpha_s^2ds,
 \qquad
 U_t=\left\langle\log{g\over1-g}\right\rangle_t .
\]

Then, pathwise before the stop,

\[
 \boxed{\qquad dU_t=\kappa_t^2\,dR_t,
 \qquad
 dE_t=r_t\,dR_t=G_t\,dU_t,\qquad}                    \tag{0.2}
\]

where `E_t=int d<g>/(g(1-g))` is the binary-entropy clock.  Thus `alpha`
changes the speed of the experiment but cannot improve the information cost
per unit of curvature in the currently exceptional direction.

The accumulated quadratic potential is

\[
 Q_t=\int_0^t C_s^2ds
     =tI-\int_0^t(1-\alpha_s^2)e_se_s^Tds,             \tag{0.3}
\]

and for every unit vector `u`,

\[
 u^TQ_tu=R_t+\int_0^t(1-\alpha_s^2)
                    |P_{u^\perp}e_s|^2ds.             \tag{0.4}
\]

In particular `lambda_min(Q_t)>=R_t`, and all but possibly one eigenvalue
of `Q_t` are at least `(t+R_t)/2`.

Two pointwise-optimal feedbacks follow immediately.  If the allowed
entropy-clock rate is `b`, the maximal possible longitudinal curvature rate
is

\[
 \alpha_t^2=\min\left\{1,{b\over r_t}\right\}.        \tag{0.5}
\]

If the sharper log-odds clock is capped at rate `b`, the maximal rate is

\[
 \alpha_t^2=\min\left\{1,{b\over\kappa_t^2}\right\}
 =\min\left\{1,{bG_t\over r_t}\right\}.               \tag{0.6}
\]

No other scalar feedback satisfying the same pointwise clock cap gives more
curvature at the same posterior state.  Smooth versions, for example

\[
 \alpha_t^2={b\over b+\kappa_t^2},                    \tag{0.7}
\]

lose at most a factor two and obey the exact affine identity

\[
                       R_t+{U_t\over b}=t             \tag{0.8}
\]

on every unstopped interval.

Stopping the log-odds diffusion at `[delta,1-delta]` gives the sharp expected
budget

\[
 \mathbb E U_{T\wedge\tau_\delta}
 \le 2\Lambda(1-2\delta),
 \qquad \Lambda=\log{1-\delta\over\delta},             \tag{0.9}
\]

and any deterministic bound `U<=u` leaves a strictly positive, explicit
survival probability.  These facts solve the mass-survival part.

They do **not** by themselves create a universal lower eigenvalue.  Under a
clock cap, failure of `R_t` to grow is exactly persistent large conditional-
barycenter separation, equivalently persistent large protected variance.
Brascamp--Lieb shows that after a positive seed has formed this obstruction
decays at the optimal exponential rate.  The estimate degenerates at seed
zero, however, and a fixed survival budget can amplify a seed by only a
fixed factor.

Consequently the entire soft-control proposal reduces to the following
initial-seed assertion:

> There are numerical `rho,p>0` such that, for every isotropic log-concave
> law and every half-mass set, the process reaches `R=rho` before
> `g` exits `[1/4,3/4]` with probability at least `p`.

This assertion would itself imply the dimension-free Cheeger inequality by
the endpoint transfer in Section 7.  It is not proved by the entropy clock,
the quantization gap, or the product-exponential calculations below.  No
smooth isotropic log-concave counterexample is produced either.  What is
proved here is the exact reduction and the strongest scalar feedback
available without that new seed theorem.

## 1. Stopped construction and exact SDEs

Let `p_0=exp(-V_0)` be a full-dimensional log-concave density with finite
second moment.  Initially work with a smooth positive density on a compact
convex support.  For a predictable symmetric control `C_t`, define

\[
 p_t(x)={1\over Z_t}
 \exp\left(c_t\cdot x-{1\over2}x^TQ_tx\right)p_0(x),   \tag{1.1}
\]

where

\[
 dc_t=C_t\,dW_t+C_t^2a_t\,dt,
 \qquad dQ_t=C_t^2dt,
 \qquad a_t=\mathbb E_tX.                              \tag{1.2}
\]

Stop before the natural parameters or the first four moments leave a fixed
compact set, and also at

\[
 \tau_\delta=\inf\{t:g_t\notin(\delta,1-\delta)\}.     \tag{1.3}
\]

After this stop one may set `C=0`.  Ito's formula gives

\[
 dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle .        \tag{1.4}
\]

Write

\[
 Y=X-a_t,\qquad h=\mathbf1_S-g_t,\qquad
 A_t=\mathbb E_tYY^T,
 \qquad D_t=\mathbb E_t[hYY^T].                        \tag{1.5}
\]

For `v_t!=0`, put `ell_t=|v_t|`, `e_t=v_t/ell_t`, and
`P_t=I-e_te_t^T`.  Choose

\[
 C_t=P_t+\alpha_te_te_t^T,
 \qquad 0<\alpha_t\le1.                                \tag{1.6}
\]

The exact stopped equations are

\[
 \boxed{\begin{aligned}
 da_t&=A_tC_t\,dW_t,\\
 dg_t&=\alpha_t\ell_t\,d\beta_t,\\
 dv_t&=D_tC_t\,dW_t-\alpha_t^2A_tv_t\,dt,\\
 dA_t&={\cal T}_t(C_t\,dW_t)-A_tC_t^2A_t\,dt,
 \end{aligned}}                                       \tag{1.7}
\]

where

\[
 {\cal T}_t(z)=\mathbb E_t[YY^T\langle Y,z\rangle]
\]

and `beta` is the Brownian motion in the current `e_t` direction.  The drift
in `dv_t` is absent for the hard mass-preserving projection, but it is
essential for every genuinely soft control.

Here is a verification of that drift.  If
`m_{S,t}=E_t[1_SX]`, then

\[
 dm_{S,t}=(g_tA_t+D_t+a_tv_t^T)C_t\,dW_t.
\]

On the other hand, the martingale part of `d(g_ta_t)` is
`(g_tA_t+a_tv_t^T)C_tdW_t`, while

\[
 d[g,a]_t=A_tC_t^2v_t\,dt.
\]

Subtracting proves the third line of (1.7).

For completeness, let

\[
 K_t=D_tC_t^2D_t.
\]

While `ell_t>0`, normalization of the vector semimartingale in (1.7) gives

\[
\begin{aligned}
 d\ell_t={}&e_t^TD_tC_t\,dW_t
 -\alpha_t^2\ell_t e_t^TA_te_t\,dt
 +{\operatorname {tr}(P_tK_t)\over2\ell_t}\,dt,\\
 de_t={}&{1\over\ell_t}P_tD_tC_t\,dW_t
 -\alpha_t^2P_tA_te_t\,dt
 -{1\over\ell_t^2}P_tK_te_t\,dt
 -{\operatorname {tr}(P_tK_t)\over2\ell_t^2}e_t\,dt,\\
 d[e]_t={}&{1\over\ell_t^2}P_tK_tP_t\,dt .            \tag{1.8}
\end{aligned}
\]

These formulas display the two effects of softening: longitudinal
information moves `g`, and the term `-alpha^2 P A e` changes the selected
direction itself.

## 2. The two intrinsic information clocks

Put

\[
 G_t=g_t(1-g_t),\qquad
 r_t={\ell_t^2\over G_t},\qquad
 \kappa_t={\ell_t\over G_t}.                           \tag{2.1}
\]

Since

\[
 v_t=G_t\left(\mathbb E_t[X\mid S]
                   -\mathbb E_t[X\mid S^c]\right),    \tag{2.2}
\]

`kappa_t` is exactly the distance between the two posterior conditional
barycenters and `r_t=G_t kappa_t^2`.

### 2.1 Binary entropy

Let

\[
 H(x)=-x\log x-(1-x)\log(1-x).
\]

Because `H''(x)=-1/[x(1-x)]`, (1.7) gives

\[
 dH(g_t)=H'(g_t)\,dg_t-{1\over2}\alpha_t^2r_t\,dt.    \tag{2.3}
\]

Define the entropy clock

\[
 E_t=\int_0^{t\wedge\tau_\delta}\alpha_s^2r_s\,ds
     =\int_0^{t\wedge\tau_\delta}
             {d\langle g\rangle_s\over G_s}.           \tag{2.4}
\]

Optional stopping in the bounded stopped construction yields the exact
identity

\[
 \boxed{\qquad
 \mathbb E E_t=2\{H(g_0)-\mathbb EH(g_{t\wedge\tau_\delta})\}.
 \qquad}                                               \tag{2.5}
\]

In particular, when `g_0=1/2`,

\[
 \mathbb E E_t\le
 2\{\log2-H(\delta)\}.                                 \tag{2.6}
\]

For `delta=1/4`, the right side is approximately `0.2616`.

The unweighted quadratic variation has the independent exact identity

\[
 \mathbb E\langle g\rangle_{t\wedge\tau_\delta}
 =\mathbb E(g_{t\wedge\tau_\delta}-g_0)^2.             \tag{2.6a}
\]

For `g_0=1/2` this is at most `(1/2-delta)^2`, and at
`delta=1/4` it is at most `1/16`.

### 2.2 Log odds

Let

\[
 L_t=\log{g_t\over1-g_t},\qquad
 U_t=\langle L\rangle_{t\wedge\tau_\delta}.
\]

Then

\[
 dU_t=\alpha_t^2\kappa_t^2dt,
 \qquad dE_t=G_t\,dU_t.                                \tag{2.7}
\]

In its own quadratic-variation time, `L` satisfies the universal diffusion

\[
 dL_u=dB_u+{1\over2}\tanh(L_u/2)\,du,                 \tag{2.8}
\]

killed at `+-Lambda`, where

\[
 \Lambda=\log{1-\delta\over\delta}.                   \tag{2.9}
\]

Indeed `2g-1=tanh(L/2)`, and Ito's formula applied to the log odds gives
(2.8).

The mean intrinsic exit time is explicit.  The function

\[
 F_\delta(x)=2\left\{
 \Lambda\tanh(\Lambda/2)-x\tanh(x/2)\right\}          \tag{2.10}
\]

vanishes at `+-Lambda` and satisfies

\[
 \left({1\over2}{d^2\over dx^2}
 +{1\over2}\tanh(x/2){d\over dx}\right)F_\delta=-1. \tag{2.11}
\]

Consequently

\[
 \boxed{\qquad
 \mathbb E U_{t\wedge\tau_\delta}
 \le F_\delta(L_0).
 \qquad}                                               \tag{2.12}
\]

For a half-mass initial set this becomes

\[
 \mathbb E U_{t\wedge\tau_\delta}
 \le2\Lambda\tanh(\Lambda/2)
 =2\Lambda(1-2\delta).                                 \tag{2.13}
\]

Equality holds if the intrinsic diffusion is run all the way to its exit.
For `delta=1/4`, (2.13) is `log 3`.

## 3. Explicit survival from a clock cap

Suppose first that `U_{T\wedge\tau_delta}<=u` pathwise.  Since the drift in
(2.8) has absolute value at most `1/2`, the reflection principle gives, for
`u<2Lambda`,

\[
 \mathbb P\{\tau_\delta\le T\}
 \le4\Phi\left(-{\Lambda-u/2\over\sqrt u}\right).     \tag{3.1}
\]

There is also a bound directly from `g`.  Since

\[
 d\langle g\rangle=G^2dU\le {1\over16}dU,
\]

Doob's inequality yields

\[
 \mathbb P\{\tau_\delta\le T\}
 \le {u\over16(1/2-\delta)^2}.                         \tag{3.2}
\]

Both bounds are dimension free.

For example, when `delta=1/4`, (3.2) reduces to

\[
                    \mathbb P\{\tau_{1/4}\le T\}\le u.
                                                               \tag{3.2a}
\]

Thus the pathwise log-odds budget `U_T<=1/2` leaves survival probability at
least `1/2`.

If only the entropy clock is capped, say `E_T<=w`, then on the stopped band
`G>=G_delta:=delta(1-delta)`, so

\[
                         U_T\le {w\over G_\delta}.      \tag{3.3}
\]

Thus a pathwise finite entropy-rate budget over a fixed time leaves a
strictly positive universal survival probability.  An even simpler, though
weaker, proof uses Dambis--Dubins--Schwarz directly: if
`dot E<=b`, then

\[
 \langle g\rangle_T=\int_0^TG_t\,dE_t\le {bT\over4},
\]

and the event that a Brownian motion stays in
`(-(1/2-delta),1/2-delta)` until time `bT/4` forces survival.

The sharp budgets (2.5) and (2.12) should not be confused with a curvature
estimate.  They say how much information can be spent before absorption;
they do not say that a fixed portion of physical or curvature time has been
spent.

## 4. Accumulated curvature: exact spectrum

On an active interval of length `t`, put

\[
 R_t=\int_0^t\alpha_s^2ds,
 \qquad
 M_t=\int_0^t(1-\alpha_s^2)e_se_s^Tds.                 \tag{4.1}
\]

If the whole process is frozen at `tau_delta`, every occurrence of `t` in
this section means the active time `t wedge tau_delta`.

Then

\[
 Q_t=tI-M_t.                                           \tag{4.2}
\]

For every unit `u`,

\[
\begin{aligned}
 u^TQ_tu
 &=\int_0^t\{1-(1-\alpha_s^2)\langle u,e_s\rangle^2\}ds\\
 &=R_t+\int_0^t(1-\alpha_s^2)|P_{u^\perp}e_s|^2ds.
                                                               \tag{4.3}
\end{aligned}
\]

This proves

\[
                         Q_t\succeq R_tI.              \tag{4.4}
\]

It also isolates the only other source of minimum curvature: rotation of
the protected direction.

Let `q_1<=q_2<=...` be the eigenvalues of `Q_t`.  Since

\[
 \operatorname {tr}M_t=t-R_t,
\]

if `rho_1>=rho_2>=...` are the eigenvalues of `M_t`, then

\[
 \rho_2\le(t-R_t)-\rho_1=q_1-R_t.
\]

Therefore

\[
 \boxed{\qquad
 q_2\ge\max\{q_1,t-q_1+R_t\}
       \ge {t+R_t\over2}.
 \qquad}                                               \tag{4.5}
\]

Thus soft localization retains the deterministic all-but-one-direction
curvature of hard mass preservation and adds the scalar floor `R_t`.

If `q_1` is small and `u` is the corresponding eigenvector, (4.3) gives

\[
 \int_0^t(1-\alpha_s^2)|P_{u^\perp}e_s|^2ds
 =q_1-R_t.                                             \tag{4.6}
\]

Hence low full curvature requires both small longitudinal expenditure and
persistence of one pathwise direction.  Angularly diffuse product models do
not automatically obstruct this control.

## 5. Pointwise-optimal scalar feedbacks

### 5.1 Entropy-clock cap

Fix `b>0`.  Requiring

\[
 {dE_t\over dt}=\alpha_t^2r_t\le b                  \tag{5.1}
\]

at every possible posterior state forces

\[
 \alpha_t^2\le a_b(r_t):=min\left\{1,{b\over r_t}\right\},
                                                               \tag{5.2}
\]

with `a_b(0)=1`.  Taking equality is therefore the unique pointwise maximal
longitudinal curvature rate under the cap.  Since only the `e_t` eigenvalue
of `C_t^2` depends on `alpha`, this is also Loewner-optimal at the current
state.

For this hard-cap feedback,

\[
 \dot E_t=\min\{r_t,b\},\qquad
 \dot R_t=\min\{1,b/r_t\},                             \tag{5.3}
\]

and

\[
 1\le \dot R_t+{\dot E_t\over b}\le2.                 \tag{5.4}
\]

The smooth feedback

\[
 \widetilde a_b(r)={b\over b+r}                       \tag{5.5}
\]

is within a factor two of (5.2) and satisfies the exact identity

\[
             \dot R_t+{\dot E_t\over b}=1.             \tag{5.6}
\]

It also has a canonical nonsingular matrix form:

\[
 C_t^2=I-{v_tv_t^T\over bG_t+|v_t|^2}.                 \tag{5.7}
\]

The right side tends smoothly to `I` at `v=0`.

### 5.2 Log-odds-clock cap

For the sharper constraint

\[
                       \dot U_t\le b,                  \tag{5.8}
\]

the exact optimal rate is

\[
 \alpha_t^2=\min\left\{1,{b\over\kappa_t^2}\right\}
 =\min\left\{1,{bG_t\over r_t}\right\}.               \tag{5.9}
\]

The smooth version

\[
 \alpha_t^2={b\over b+\kappa_t^2}
 ={bG_t\over bG_t+r_t}                                 \tag{5.10}
\]

satisfies

\[
                       \dot R_t+{\dot U_t\over b}=1.    \tag{5.11}
\]

It has the nonsingular matrix representation

\[
 C_t^2=I-{v_tv_t^T\over bG_t^2+|v_t|^2}.              \tag{5.11a}
\]

If the feedback is required to depend on `r_t` but not on the current mass,
the conservative choice

\[
 \alpha_t^2=\min\left\{1,{bG_\delta\over r_t}\right\} \tag{5.12}
\]

guarantees (5.8) throughout the stopped band.

### 5.3 Protected-variance feedback is not pointwise stronger

Put

\[
                         s_t=e_t^TA_te_t.              \tag{5.13}
\]

The law of total covariance and (2.2) give

\[
                         r_t\le s_t.                   \tag{5.14}
\]

There is in fact a universal strict gap.  The two-point quantization lemma
for log-concave laws says that, for a numerical `eta>0`,

\[
 v_t^TA_t^{-1}v_t\le(1-\eta)G_t.                       \tag{5.15}
\]

Covariance Cauchy applied to `e_t` and `v_t` then yields

\[
                         r_t\le(1-\eta)s_t.            \tag{5.16}
\]

Thus

\[
 \alpha_t^2=\min\{1,b/s_t\}                           \tag{5.17}
\]

does cap the entropy clock, but it is never larger than the `r`-feedback
(5.2) with the same `b`.  A protected-variance rule can be useful only if
one proves additional dynamics for `s_t`; as a statewise solution of the two
objectives it is dominated by the normalized-correlation rule.

## 6. What Brascamp--Lieb gives after a seed

The endpoint potential in (1.1) has distributional Hessian at least `Q_t`.
By (4.4), once `R_t>0`, Brascamp--Lieb gives

\[
 A_t\preceq Q_t^{-1}\preceq R_t^{-1}I.                \tag{6.1}
\]

Combining (5.16) and (6.1),

\[
 r_t\le {1-\eta\over R_t},
 \qquad
 \kappa_t^2\le {1-\eta\over G_tR_t}
 \le {1-\eta\over G_\delta R_t}.                     \tag{6.2}
\]

For the entropy-optimal feedback (5.2), this implies

\[
 \dot R_t\ge
 \min\left\{1,{bR_t\over1-\eta}\right\}.             \tag{6.3}
\]

For the log-odds-optimal feedback (5.9),

\[
 \dot R_t\ge
 \min\left\{1,{bG_\delta R_t\over1-\eta}\right\}.    \tag{6.4}
\]

Consequently, if `R_{t_0}=epsilon>0`, then before the right side reaches one,

\[
 R_t\ge\varepsilon
 \exp\left({bG_\delta\over1-\eta}(t-t_0)\right)       \tag{6.5}
\]

for the log-odds rule; after that point it grows at unit speed.  This is the
best conclusion available from the quantization gap and
Brascamp--Lieb.

The degeneracy at `R=0` is genuine.  More importantly, keeping the total
log-odds budget `bT` universal makes the exponent in (6.5) at most

\[
                       {G_\delta\over1-\eta}\,bT,      \tag{6.6}
\]

also universal.  Hence this mechanism amplifies an initial seed by only a
fixed factor; it cannot turn a dimension-dependent seed tending to zero
into a universal one.

The same point can be written as an occupation estimate.  Under any control,
for every `K>0`, (0.2) gives

\[
 R_T
 \le |\{t\le T:\kappa_t\le K\}|+{U_T\over K^2}.      \tag{6.7}
\]

For the pointwise-optimal capped control, a lower curvature theorem is
therefore exactly a theorem that the process spends a fixed amount of time
with moderate posterior conditional-barycenter separation.  Equation (6.2)
controls this occupation only after the desired seed already exists.

### 6.1 The seed statement is load-bearing

Fix `delta=1/4`.  The remaining assertion is

\[
 \boxed{\quad
 \exists\rho_0,p_0>0:\quad
 \mathbb P\{\tau_{\rho_0}^R<\tau_{1/4}\}\ge p_0
 \quad}                                                \tag{6.8}
\]

for every isotropic log-concave initial law and every half-mass Borel set,
where

\[
 \tau_{\rho_0}^R=\inf\{t:R_t\ge\rho_0\}.
\]

Neither (2.12) nor (6.2) proves (6.8): the former bounds information spent,
whereas the latter has the nonunique zero solution at the comparison-ODE
level.  A uniform proof of (6.8) must control creation of a large protected
variance in an arbitrarily short initial layer.  This is the same adaptive
winner-selection obstruction isolated for the hard mass-preserving flow.

Section 7 shows that (6.8) immediately gives a universal Cheeger bound.
Thus it should not be treated as a routine stochastic-continuity lemma.

### 6.2 The entropy-plus-quantization stopping argument is insufficient

There is a precise logical obstruction to proving (6.8) from only the
intrinsic diffusion and the bound in (6.2).  Write curvature itself as the
time variable `q=R`, and put

\[
 c={1-\eta\over G_\delta},
 \qquad K_0=4(1-\eta).                                  \tag{6.9}
\]

The initial isotropic quantization bound permits `kappa_0^2<=K_0`, while
(6.2) permits `kappa_q^2<=c/q` for `q>0`.  For every small `epsilon>0`, the
continuous deterministic profile

\[
 \kappa_q^2=
 \begin{cases}
 K_0+(c/\varepsilon-K_0)q/\varepsilon,
      &0\le q\le\varepsilon,\\
 c/q,&\varepsilon\le q\le\rho
 \end{cases}                                            \tag{6.10}
\]

obeys both restrictions, provided `epsilon` is small enough.  Indeed
`q kappa_q^2` is increasing to `c` on the first interval and equals `c`
afterward.  But its log-odds expenditure before curvature `rho` is

\[
 \int_0^\rho\kappa_q^2dq
 ={c\over2}+O(\varepsilon)
   +c\log{\rho\over\varepsilon}\longrightarrow\infty.
                                                               \tag{6.11}
\]

The diffusion (2.8) exits every bounded interval with probability tending
to one as its intrinsic running time tends to infinity.  Hence the proposed
two inequalities allow the central mass to be lost with probability tending
to one before any fixed `rho` is reached.

Profile (6.10) is an abstract coefficient path, not a localization
trajectory.  Its role is to prove that entropy, isotropic initial
quantization, continuity, and the Brascamp--Lieb estimate do not logically
imply the seed theorem.  A successful proof must add a dynamical restriction
which prevents the rapid rise to the critical profile `kappa_q^2 about 1/q`.
Producing a smooth isotropic log-concave trajectory which actually realizes
that profile would disprove (6.8); no such realization is established here.

## 7. Endpoint transfer and why a seed proves KLS

For every fixed `x`, (1.4) makes `p_t(x)/p_0(x)` a nonnegative stopped
martingale.  In the compact stopped construction it is a true martingale,
so

\[
                         \mathbb E p_T(x)=p_0(x).       \tag{7.1}
\]

For a finite-perimeter set, Tonelli on the reduced boundary gives equality
of expected weighted perimeters.  For an arbitrary Borel set, exterior
Minkowski content and Fatou give the direction needed here:

\[
 \mathbb E\,\mu_T^+(S)\le\mu_0^+(S).                  \tag{7.2}
\]

Indeed, apply Fatou to

\[
 {\mu_T(S_\varepsilon)-\mu_T(S)\over\varepsilon}
\]

and use (7.1) before taking `epsilon downarrow0`.

On a path for which `Q_T\succeq rho I`, the posterior is
`rho`-strongly log-concave.  The strongly-log-concave Cheeger inequality
therefore gives

\[
 \mu_T^+(S)\ge c\sqrt\rho\min(g_T,1-g_T).              \tag{7.3}
\]

If the path reaches `R=rho` before leaving `[delta,1-delta]`, stop the whole
process at that time.  Equations (4.4), (7.2), and (7.3) imply

\[
 \boxed{\qquad
 \mu_0^+(S)\ge
 c\,\delta\sqrt\rho\,
 \mathbb P\{\tau_\rho^R<\tau_\delta\}.
 \qquad}                                               \tag{7.4}
\]

Thus (6.8) is already a dimension-free isoperimetric theorem.  Conversely,
the present calculation provides no independent proof of its probability
lower bound.

For an unbounded hitting time, apply the argument first at
`(tau_rho^R wedge tau_delta wedge m)`.  Retain in (7.3) only the paths on
which `tau_rho^R<=tau_delta wedge m`, and then let `m` increase to infinity.
The events increase to `{tau_rho^R<tau_delta}`, so monotone convergence gives
(7.4) without assuming uniform integrability at the unbounded stopping time.

## 8. Sharp tests and countermodels

### 8.1 A Gaussian halfspace rules out deterministic success

Take `n=1`, `mu_0=N(0,1)`, and `S=[0,infinity)`.  In one dimension the
curvature clock is exactly

\[
                         q=R=Q.                        \tag{8.1}
\]

Changing `alpha` only reparametrizes this clock.  At curvature time `q`, the
posterior is

\[
 N\left({c_q\over1+q},{1\over1+q}\right),
 \qquad
 g_q=\Phi\left({c_q\over\sqrt{1+q}}\right).            \tag{8.2}
\]

Under the unstopped path law,

\[
                         \Phi^{-1}(g_q)\sim N(0,q).     \tag{8.3}
\]

Therefore, with `z_delta=Phi^{-1}(1-delta)`,

\[
 \mathbb P\{\tau_\delta\le q_0\}
 \ge2\Phi\left(-{z_\delta\over\sqrt{q_0}}\right)>0.  \tag{8.4}
\]

For `delta=1/4` and `q_0=1`, the right side is exactly `1/2`.  Hence no
policy can guarantee a fixed curvature pathwise before mass exit.  This
example is compatible with a small-curvature, positive-probability seed and
does not refute (6.8).

### 8.2 An exponential product is a sharp statewise obstruction

Let `Z` be rate-one exponential and let

\[
 X_1=\sigma(Z_1-1),
\]

with the remaining coordinates independent centered rate-one exponentials.
Take

\[
                         S=\{Z_1\ge\log2\}.            \tag{8.5}
\]

Then `g=1/2`, `e=e_1`, and direct integration gives

\[
 s=e^TAe=\sigma^2,
 \qquad
 v_1={\sigma\log2\over2},
 \qquad
 r=\sigma^2(\log2)^2.                                  \tag{8.6}
\]

Thus any feedback which caps the entropy-clock rate by `b` must have

\[
                 \alpha^2\le {b\over\sigma^2(\log2)^2}.
                                                               \tag{8.7}
\]

This is a genuine log-concave product and shows that balance plus posterior
log-concavity alone cannot give an instantaneous curvature floor.  It is
anisotropic.  After whitening, (8.6) has `s=1`, so it is not a counterexample
to the isotropic initial-seed assertion (6.8).

The support edge is inessential for this statewise example.  Convolving
each exponential factor with a centered Gaussian of variance `epsilon^2`
produces a positive smooth log-concave product.  Choose the first-coordinate
median as the new threshold and then let `epsilon` decrease to zero.  Its
values of `g`, `s`, and `r/s` converge respectively to
`1/2`, `sigma^2`, and `(log2)^2`, so the bound (8.7) persists with an
arbitrarily small relative error.  The obstruction can therefore be made
smooth, but still not isotropic.

It can also be viewed dynamically.  Because the set depends only on the
first factor, the product structure keeps `e=e_1`.  After rescaling the first
coordinate back by `sigma`, the longitudinal natural parameter and
quadratic parameter are respectively multiplied by `sigma` and `sigma^2`.
On any fixed stopped neighborhood on which the rescaled one-dimensional
correlation stays bounded above and below, the capped feedback therefore
has an order-one rescaled clock while its curvature in the original
coordinate is of order `sigma^{-2}`.  Hence no theorem uniform over arbitrary
anisotropic starting posteriors is possible.  The missing content in (6.8)
is precisely that every such posterior must be reached from one isotropic
initial law along the same localization path.

### 8.3 The isotropic maximum-tail product is not a counterexample

Let `Z_1,...,Z_n` be independent rate-one exponentials, set
`X_i=Z_i-1`, and choose `L` from

\[
 (1-e^{-L})^n={1\over2}.
\]

For

\[
                         S=\{\max_iZ_i\ge L\},          \tag{8.8}
\]

write `q=e^{-L}=1-2^{-1/n}`.  Exchangeability gives

\[
 |v_0|={\sqrt nLq\over2(1-q)},
 \qquad
 \kappa_0={|v_0|\over1/4}
 ={2\sqrt nLq\over1-q}
 \sim {2(\log2)\log n\over\sqrt n}.                  \tag{8.9}
\]

Consequently full longitudinal exposure initially costs

\[
                         \dot U_0=\kappa_0^2=o(1),      \tag{8.10}
\]

even though the exact angular quadratic-variation rate of `e_t` under the
hard projection is asymptotic to `(log n)^2`.  Indeed the transverse
eigenvalue of `D_0` is

\[
 d_n={L^2q\over2(1-q)^2},
\]

and hence

\[
 {d\over dt}\operatorname {tr}[e]_t\bigg|_{t=0}
 ={(n-1)d_n^2\over|v_0|^2}
 ={n-1\over n}{L^2\over(1-q)^2}.                       \tag{8.11}
\]

Fast angular motion is helpful in (4.3), and by itself does not obstruct a
soft seed.

The full-exposure large-deviation calculation for this example further
shows that a single coordinate does not acquire order-one posterior tail
mass on the `1/log n` curvature scale.  At time `t=zeta/L`, the largest
one-coordinate natural parameter satisfies

\[
 c_*(\zeta)=
 \begin{cases}
  \sqrt{2\zeta}+o_{\mathbb P}(1),&\zeta\le1/2,\\
  \zeta+1/2+o_{\mathbb P}(1),&\zeta\ge1/2.
 \end{cases}                                            \tag{8.12}
\]

For a coordinate with natural parameter `c`, its posterior probability of
exceeding `L` has exponential rate

\[
 J_\zeta(c)=
 \begin{cases}
  1+\zeta/2-c,&c\le1,\\
  (1+\zeta-c)^2/(2\zeta),&1<c<1+\zeta,\\
  0,&c\ge1+\zeta.
 \end{cases}                                            \tag{8.13}
\]

At (8.12),

\[
 J_\zeta(c_*)=
 \begin{cases}
  1+\zeta/2-\sqrt{2\zeta},&\zeta\le1/2,\\
  1/(8\zeta),&\zeta\ge1/2,
 \end{cases}                                            \tag{8.14}
\]

which is positive for every fixed `zeta`.  These formulas follow by a
one-dimensional Laplace principle applied first to
`c=tZ+sqrt(t)G` and then to the posterior density
`exp((c-1)z-tz^2/2)`.

What remains unproved is a uniform occupation estimate preventing a diffuse
collection of candidates from creating and then maintaining `kappa_t^2` of
order `1/R_t`.  The large-deviation calculation is for full exposure and
does not control the evolving transverse projection in (0.1).  Thus this
model is a sharp stress test for (6.8), not a counterexample to it.

## 9. Existence at `v=0`, removal of stops, and approximation

The hard mass-preserving projection is discontinuous at `v=0`.  The soft
`r`-feedbacks above can avoid that defect completely.

For the hard entropy cap (5.2), `alpha=1` on the open neighborhood
`|v|^2<=bG` of `v=0`; hence `C=I` there.  Across the threshold the radial
matrix coefficient is locally Lipschitz.  The smooth feedback (5.5) has the
explicit form (5.7), which is smooth in `(v,g)` while
`g in [delta,1-delta]`.  Its positive square root is also locally Lipschitz:

\[
 C=I-(1-\sqrt{bG/(bG+|v|^2)}){v v^T\over|v|^2},       \tag{9.1}
\]

where the final term has a removable zero at `v=0`.

The hard log-odds cap (5.9) is identically `I` on
`|v|^2<=bG^2`, and its smooth version has the same locally Lipschitz
extension by (5.11a).  Thus the sharper clock requires no relaxed convention
at `v=0` either.

For a compactly supported density the natural-parameter moment maps
`a(c,Q)`, `A(c,Q)`, `g(c,Q)`, and `v(c,Q)` are smooth.  The finite-dimensional
SDE (1.2) with either of the preceding feedbacks therefore has a unique
strong solution up to the parameter and mass stops.  All stochastic
integrals in Sections 1--8 are then bounded stopped martingales, so optional
stopping and Tonelli are legitimate.

A feedback based on `s=e^TAe` alone need not approach `I` as `v->0`, because
the direction `e` is then undefined.  Such a rule requires either an
explicit convention plus a weak-existence argument, or a relaxed-control
limit.  This is another reason to prefer the `r`-feedback, which is also
pointwise stronger for the stated objectives.

For an unbounded log-concave density, first truncate to increasing convex
compact sets, smooth by convolution, recenter and whiten, and impose the
parameter stop.  For an initially isotropic law the recentering and
whitening maps converge to the identity.  The moment maps converge locally
uniformly on every stopped parameter set, and the corresponding solution
laws are tight.  Equations (2.5), (2.12), (4.3), and (7.2) are stable under
weak convergence and Fatou.  One then takes an increasing sequence of
parameter stops and finally removes the truncation.  If a quantitative seed
estimate such as (6.8) were proved uniformly in these approximations, (7.4)
would pass to the original law.  Without uniformity, pointwise likelihood
local martingales must not be treated as automatically uniformly integrable.

## 10. Exact remaining theorem

The soft-control route is now reduced without hidden operator norms or
unstated endpoint transfers.  A closing result must prove one of the
following comparable alternatives for a half-mass cut of an isotropic
log-concave law:

1. a positive-probability initial seed as in (6.8);
2. a dimension-free lower bound on the occupation time
   `|{t:kappa_t<=K}|` before the central-mass exit, for some numerical `K`;
3. a direct perimeter gain on paths where `kappa_t^2` grows like `1/R_t`
   during the initial layer.

The exact identities show why merely tuning `alpha` cannot supply this
input.  Entropy and log-odds clocks solve survival, Brascamp--Lieb solves
post-seed growth, and rotation solves every direction except a persistent
pathwise line.  Uniformly creating curvature in that final line is the new
content.  Establishing it would prove the desired dimension-free
isoperimetric bound; disproving it requires a smooth isotropic log-concave
trajectory with persistent early winner-locking, not just a posterior state
with large variance or a large instantaneous angular rate.
