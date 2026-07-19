# The protected direction: exact SDEs and a sharp product stress test

This note studies the one uncontrolled variance in the mass-preserving
localization

\[
 C_t=P_t:=I-e_te_t^T,\qquad
 e_t={v_t\over |v_t|},\qquad
 v_t=\operatorname {Cov}_{\mu_t}(1_S,X).
\]

It gives exact stopped equations, including the equation for the exceptional
eigendirection of the accumulated curvature.  It also gives a product
example in which the angular noise and the selection drift of the protected
variance diverge although every one-dimensional moment in the protected
direction is bounded.  Consequently, a scalar Lyapunov argument which uses
only the protected variance and the one-dimensional third/fourth moments
cannot close with a universal constant.  The calculation does **not** refute
the endpoint exceptional-variance estimate: a final subsection explains the
additional persistence issue that an actual counterexample must solve.

## 1. Stopped identities

Work first with a smooth positive log-concave density on a compact convex
set.  Fix a Borel set `S`, put `g=mu_t(S)`, and stop before

\[
 |v_t|\leq\varepsilon,
 \quad \|A_t\|+|a_t|+\|c_t\|+\|Q_t\|\geq R.
\]

All coefficients below are then bounded and the identities follow from
ordinary Ito calculus.  They therefore also hold up to the corresponding
localizing times in any weak limit of compactly supported approximations.
Put

\[
 Y=X-a_t,\qquad h=1_S-g,
\]

and define the symmetric matrix and third-moment map

\[
 D_t=\mathbb E_t[hYY^T],\qquad
 {\cal T}_t(z)=\mathbb E_t[YY^T\langle Y,z\rangle].       \tag{1}
\]

The density equation is

\[
 d p_t(x)=p_t(x)\langle x-a_t,P_t,dW_t\rangle .          \tag{2}
\]

Since `P_t v_t=0`, the set mass is constant.  Direct differentiation gives

\[
 \boxed{\begin{aligned}
 dg_t&=0,\\
 da_t&=A_tP_t,dW_t,\\
 dv_t&=D_tP_t,dW_t,\\
 dA_t&={\cal T}_t(P_t,dW_t)-A_tP_tA_t,dt.
 \end{aligned}}                                          \tag{3}
\]

For the third line, differentiating
`E_t[1_SX]-g_ta_t` initially produces
`(D_t+a_tv_t^T)P_t dW_t`; the second term vanishes because
`v_t^TP_t=0`.  Thus there is no omitted finite-variation term in `dv_t`.

### 1.1 Length and direction of the binary correlation

Write `r=|v|`, `e=v/r`, and `P=I-ee^T`.  Set

\[
 H=DPD,
 \qquad
 \Gamma=\operatorname {Tr}((PDP)^2).
\]

Ito's formula gives

\[
 dr=e^TDP,dW+
 {\operatorname {Tr}H-e^THe\over2r},dt,                 \tag{4}
\]

and

\[
 \boxed{
 de={1\over r}PDP,dW
 -{1\over r^2}PDPDe,dt
 -{\Gamma\over2r^2}e,dt .}                             \tag{5}
\]

In particular,

\[
 d[e]_t={1\over r^2}(PDP)^2dt,
 \qquad
 {d\over dt}\operatorname {Tr}[e]_t={\Gamma\over r^2}. \tag{6}
\]

The last two drift terms in (5) are both necessary: the first is the
tangential Ito correction and the second is the radial correction which
keeps `|e|=1`.

### 1.2 Instantaneous protected variance

Let

\[
 s_t=e_t^TA_te_t,
 \qquad G_t={1\over r_t}P_tD_tP_t,
 \qquad \tau_t=\mathbb E_t[(e_t^TY)^2Y].
\]

For an orthonormal coordinate system `(f_k)`, put
`N_k={\cal T}(Pf_k)`.  Expanding all the quadratic covariations in
`e^TAe` gives

\[
 ds=\beta^TdW+{\cal B}\,dt,                              \tag{7}
\]

where

\[
 \boxed{
 \beta=P\tau+{2\over r}PDP,Ae}                         \tag{8}
\]

and

\[
\boxed{\begin{aligned}
 {\cal B}={}&-|PAe|^2
 -{2\over r^2}e^TAPDPDe
 -{s\over r^2}\Gamma\\
 &+{1\over r^2}\operatorname {Tr}(A(PDP)^2)
 +2\sum_k (Gf_k)^TN_ke .
\end{aligned}}                                           \tag{9}
\]

The final summand is the covariation between the covariance martingale and
the angular martingale.  Omitting it gives a false evolution equation; the
product calculation in Section 3 makes it diverge like `log n`.

A useful moment fact does hold without dimension loss.  Covariance Cauchy
and the one-dimensional fourth-moment bound for log-concave laws give

\[
 \tau^TA^{-1}\tau
 \leq\operatorname {Var}((e^TY)^2)
 \leq C(e^TAe)^2.                                       \tag{10}
\]

However, (8) is measured in the Euclidean quadratic variation.  Passing
from (10) to `|P tau|^2` costs a transverse covariance eigenvalue.  Formula
(9) shows that this is not merely a poor estimate: direction selection also
appears through a genuine third-moment covariation.

## 2. Accumulated curvature and its exceptional eigendirection

Let

\[
 M_t=\int_0^t e_se_s^Tds,\qquad Q_t=tI-M_t.              \tag{11}
\]

Suppose the largest eigenvalue `rho_t` of `M_t` is simple.  Choose a unit
eigenvector `u_t`, let `(u_j,rho_j)`, `j>=2`, be the other eigenpairs, and
put

\[
 \alpha_t=\langle u_t,e_t\rangle,
 \qquad w_t=P_{u_t^\perp}e_t,
 \qquad
 R_t=\sum_{j\geq2}{u_ju_j^T\over\rho_t-\rho_j(t)}.
\]

Ordinary perturbation theory (there is no stochastic differential in
`M_t`) gives

\[
 \boxed{\begin{aligned}
 d\rho_t&=\alpha_t^2dt,\\
 du_t&=\alpha_tR_tw_tdt.
 \end{aligned}}                                         \tag{12}
\]

The sign-free spectral projector `U_t=u_tu_t^T` satisfies

\[
 dU_t=(R_te_te_t^TU_t+U_te_te_t^TR_t)dt.                \tag{13}
\]

Put `q_t=t-rho_t`, the curvature accumulated in this running exceptional
direction.  Since `Tr M_t=t`,

\[
 \rho_2(t)\leq q_t,
 \qquad
 \rho_t-\rho_2(t)\geq t-2q_t.                           \tag{14}
\]

Consequently, throughout the genuinely exceptional regime `q_t<t/2`, the
top eigendirection is automatically simple and

\[
 \boxed{\dot q_t=|w_t|^2,
 \qquad \|R_t\|\leq(t-2q_t)^{-1}.}                     \tag{15}
\]

Thus `q_T` is exactly the time-integrated squared misalignment of the
exceptional direction and the protected direction.

For

\[
 z_t=u_t^TA_tu_t,
 \qquad b_t=P_{u_t^\perp}A_tu_t,
 \qquad \tau_{u,t}=\mathbb E_t[(u_t^TY)^2Y],
\]

the finite variation of `u_t` eliminates all `dA du` covariations, and
(3), (12) give the particularly clean equation

\[
\boxed{\begin{aligned}
 dz_t={}&(P_t\tau_{u,t})^TdW_t
 -|P_tA_tu_t|^2dt\\
 &+2\alpha_t\,b_t^TR_tw_tdt.
\end{aligned}}                                          \tag{16}
\]

The last term is the exact adaptive-selection term.  When `q_t<t/4`, its
denominator is at least `t/2`, but this only controls late rotation.  A
direction can still be selected during an initial interval comparable to
the final small curvature `q_T`; this is precisely the possible
winner-locking scenario.

Equations (5), (9), and (16) are valid after adding a stopping time at which
the relevant eigenvalue gap falls below a fixed `gamma>0`.  Formula (14)
removes that auxiliary stop whenever `q_t<t/2`.  At a zero of `v_t` one must
use the relaxed-control formulation from `mass_preserving_rankone.md`; no
SDE for `e_t` itself exists there.

### 2.1 Dimension-free control after the exceptional direction has formed

There is a useful bound which avoids `||A_t||_op`.  In the regime
`0<q_t<t/2`, the spectral facts above imply

\[
 Q_t\succeq q_tu_tu_t^T+(t-q_t)P_{u_t^\perp}.            \tag{16a}
\]

Brascamp--Lieb, first with a harmless positive regularization if needed,
therefore gives

\[
 z_t\leq q_t^{-1},\qquad
 \|P_{u_t^\perp}A_tP_{u_t^\perp}\|_{\rm op}
 \leq(t-q_t)^{-1},                                      \tag{16b}
\]

and positivity of `A_t` gives

\[
 |b_t|^2\leq {z_t\over t-q_t}.                          \tag{16c}
\]

Covariance Cauchy and the one-dimensional log-concave fourth-moment bound
give

\[
 \tau_{u,t}^TA_t^{-1}\tau_{u,t}\leq Cz_t^2.
\]

Writing `beta_t=|P_{u_t^perp}e_t|=sqrt(dot q_t)` and decomposing
`tau_{u,t}` into its `u_t` and `u_t^perp` components yields

\[
 \boxed{
 |P_t\tau_{u,t}|^2
 \leq C\left({z_t^2\over t-q_t}+\beta_t^2z_t^3\right).} \tag{16d}
\]

Indeed, dual covariance Cauchy bounds the squared transverse component by
`Cz_t^2/(t-q_t)`, the squared longitudinal component by `Cz_t^3`, and
projection orthogonal to `e_t` retains only a `beta_t` fraction of the
longitudinal component.  Likewise the selection term in (16) obeys

\[
 \boxed{
 |2\alpha_tb_t^TR_tw_t|
 \leq {2\sqrt{z_t}\,\beta_t
       \over\sqrt{t-q_t}\,(t-2q_t)}.}                  \tag{16e}
\]

Thus the late-time stochastic terms are dimension free.  These estimates
do not close the argument at `t=0`: if the endpoint exceptional variance is
`R`, Brascamp--Lieb only forces `q_T<=1/R`, and a direction may in principle
be selected during the initial layer of that same length.  The product
calculation below shows that this initial layer can have arbitrarily large
angular quadratic variation.

One can make the last assertion sharper.  Fix `R>0` and consider any time
interval contained in

\[
 \{t\geq4/R,\ q_t\leq1/R,\ z_t\geq R/2\}.
\]

The positive part of the selection drift divided by `z_t` has universally
bounded total integral.  Indeed, (16e) and `dot q=beta^2` give

\[
\begin{aligned}
 \int { (2\alpha b^TRw)_+\over z}\,dt
 &\leq {C\over\sqrt R}
       \int_{4/R}^T{\beta_t\over t^{3/2}}dt\\
 &\leq {C\over\sqrt R}
       \left(\int_0^T\beta_t^2dt\right)^{1/2}
       \left(\int_{4/R}^\infty t^{-3}dt\right)^{1/2}
 \leq C.                                                \tag{16f}
\end{aligned}
\]

Consequently, on every high-variance excursion after time `4/R`, multiplying
`z_t` by the exponential of minus this positive drift produces a local
supermartingale (up to the exit time).  Thus adaptive rotation can amplify a
high variance by at most a universal multiplicative cost during the late
part of the flow.  What remains uncontrolled is exactly the creation of
variance `R` by time `O(1/R)`; (21)--(23) show why that initial creation
cannot be discarded as a routine estimate.

### 2.2 Exponential smallness of every coherent high-variance packet

Although the exceptional direction is adaptive, a high-variance event
cannot put much probability in any fixed direction cap.  Fix a numerical
time `T>0`, a unit vector `theta`, and `delta<=1/2`.  Let `E` be any
path event on which the top eigenvector `u_T` can be oriented so that

\[
 |u_T-\theta|\leq\delta,
 \qquad u_T^TA_Tu_T\geq R.                              \tag{16g}
\]

Then, for `R>=C_T`,

\[
 \boxed{\quad \mathbb P(E)\leq C_Te^{-c\sqrt R}.\quad} \tag{16h}
\]

Here is a proof.  Couple the path with a posterior sample:

\[
 d\Pi(\omega,x)=d\mathbb P(\omega)d\mu_T^\omega(x).
\]

The `x`-marginal of `Pi` is the original isotropic law `mu_0`.  On each
path in `E`, take the upper and lower posterior quartile tails of
`<u_T,x-a_T>`, denoted `U_omega,L_omega`.  A one-dimensional log-concave
law of variance `z` has density at most `C/sqrt(z)`.  Hence its first and
third quartiles, and therefore the conditional means of its lower and upper
quartile tails, are separated by at least `c sqrt(z)`.

The endpoint curvature satisfies

\[
 Q_T\succeq (T/2)P_{u_T^\perp}.
\]

Brascamp--Lieb therefore bounds every posterior variance in `u_T^perp` by
`2/T`.  Writing `theta=alpha u_T+w`, `w perpendicular u_T`, and using
Cauchy--Schwarz on each quartile tail shows

\[
 \left\langle\theta,
 E_T[X\mid U_\omega]-E_T[X\mid L_\omega]\right\rangle
 \geq c\sqrt R-C_T\delta\geq c'\sqrt R.                \tag{16i}
\]

Define joint events

\[
 B=\{\omega\in E,x\in U_\omega\},\qquad
 C=\{\omega\in E,x\in L_\omega\}.
\]

Both have `Pi`-mass `epsilon/4`, where
`epsilon=P(E)`, and the path distribution conditioned on either is exactly
`P(.|E)`.  Averaging (16i) therefore gives

\[
 \langle\theta,E[X\mid B]-E[X\mid C]\rangle
 \geq c'\sqrt R.                                       \tag{16j}
\]

For an isotropic log-concave random vector and any (possibly randomized)
event of mass `eta`, the one-dimensional exponential-tail estimate and
layer cake give

\[
 |E[X\mid\text{event}]|\leq C(1+\log(1/\eta)).          \tag{16k}
\]

The randomized form follows because the event selects a subprobability
density bounded by one relative to `mu_0`.  Apply (16k) to `B` and `C` in
(16j) and rearrange to obtain (16h).

Thus, if exceptional variance `R` occurs with fixed positive probability,
its direction law must occupy at least `exp(c sqrt(R))` separated
unit-scale caps.  This is a genuine restriction on a counterexample, but it
does not by itself contradict high ambient dimension.

## 3. A half-mass product cut with angular rate `(log n)^2`

Let `Z_1,...,Z_n` be independent rate-one exponentials and put
`X_i=Z_i-1`.  The product law is isotropic.  Choose `L=L_n` by

\[
 (1-e^{-L})^n={1\over2},
 \qquad q=e^{-L}=1-2^{-1/n},                             \tag{17}
\]

and take

\[
 S=\{\max_i Z_i\geq L\}.
\]

Then `mu(S)=1/2` and `L=log n-log(log 2)+o(1)`.  The exterior boundary
measure is

\[
 \mu^+(S)=nq(1-q)^{n-1}longrightarrow {\log2\over2}.   \tag{18}
\]

Thus this is a stress test for the stochastic reduction, not a putative KLS
counterexample.

Exchangeability gives `e_0=n^{-1/2}(1,...,1)`.  Direct integration yields

\[
 (v_0)_i=Lq(1-q)^{n-1}={Lq\over2(1-q)},
 \qquad
 r_0={\sqrt nLq\over2(1-q)}.                            \tag{19}
\]

The matrix `D_0` is exchangeable.  Its eigenvalue on `e_0^perp` is

\[
 d_n=(D_0)_{ii}-(D_0)_{ij}
 ={L^2q\over2(1-q)^2}.                                  \tag{20}
\]

For completeness, conditional on `S^c` the coordinates are independent
exponentials truncated at `L`, and

\[
 \operatorname {Var}(Z_i\mid Z_i<L)
 =1-{L^2q\over(1-q)^2}.
\]

Using
`d_n=(1/2)E[(1_S-1/2)(X_i-X_j)^2]` gives (20).
Consequently `P_0D_0P_0=d_nP_0`, and (6) gives the exact identity

\[
 \boxed{
 {d\over dt}\operatorname {Tr}[e]_t\bigg|_{t=0}
 ={(n-1)d_n^2\over r_0^2}
 ={n-1\over n}{L^2\over(1-q)^2}
 =(1+o(1))(\log n)^2.}                                  \tag{21}
\]

This divergence is invisible in the protected one-dimensional marginal:

\[
 \operatorname {Var}(e_0^TX)=1,
 \quad E(e_0^TX)^3={2\over\sqrt n},
 \quad E(e_0^TX)^4=3+{6\over n}.                        \tag{22}
\]

There is an equally sharp failure for the naive protected-variance
generator.  For the centered exponential product,

\[
 {\cal T}_0(z)=2\operatorname {diag}(z_1,...,z_n).
\]

At time zero `A=I`, the martingale coefficient (8) vanishes, the third and
fourth terms in (9) cancel, and the only surviving drift is the last
covariation term.  Substitution gives

\[
 \boxed{
 {d\over dt}\mathbb E[e_t^TA_te_t]\bigg|_{t=0}
 ={4d_n(n-1)\over r_0\sqrt n}
 ={4L(n-1)\over n(1-q)}
 =(4+o(1))\log n.}                                      \tag{23}
\]

Therefore no inequality of the form

\[
 {\cal L}(e^TAe)\leq C\,F(e^TAe,
 \text{one-dimensional standardized moments})
\]

with dimension-free `C` and locally bounded `F` can be true.  More
generally, any proposed Lyapunov functional must explicitly pay for the
small binary correlation `r`, the transverse conditional-covariance
difference `PDP`, or the accumulated misalignment `q`; log-concave
third/fourth moment bounds alone do not control the adaptive direction.

## 4. Why Section 3 is not yet an endpoint counterexample

The example above suggests a tempting story: briefly expose all
coordinates, select the coordinate with the largest exponential tilt, and
then protect it.  The following exact large-deviation calculation shows why
that story is incomplete.

Replace the rank-one-controlled observation during the proposed exploration
phase by full scalar Gaussian observations.  For one unshifted exponential
coordinate the natural parameter after time

\[
 t={\kappa\over L}
\]

has the signal representation

\[
 c=tZ+\sqrt t\,G,
 \qquad G\sim N(0,1).                                   \tag{24}
\]

For fixed `kappa,s>0`, Laplace's method gives

\[
 \mathbb P\{c\approx s\}
 =\exp[-L I_\kappa(s)+o(L)],
\]

where

\[
 I_\kappa(s)=
 \begin{cases}
 s^2/(2\kappa),&s\leq1,\\
 (s-1/2)/\kappa,&s\geq1.
 \end{cases}                                            \tag{25}
\]

Since `log n=L+O(1)`, the maximum natural parameter among the `n`
coordinates is

\[
 c_*(\kappa)=
 \begin{cases}
 \sqrt{2\kappa}+o_P(1),&\kappa\leq1/2,\\
 \kappa+1/2+o_P(1),&\kappa\geq1/2.
 \end{cases}                                            \tag{26}
\]

The scalar posterior is proportional on the half-line to

\[
 \exp((c-1)z-\kappa z^2/(2L)).
\]

Another one-dimensional Laplace calculation shows that its posterior tail
mass above `L` has exponential rate

\[
 -{1\over L}\log\mu_{c,t}\{Z\geq L\}\longrightarrow
 J_\kappa(c)=
 \begin{cases}
 1+\kappa/2-c,&c\leq1,\\
 (1+\kappa-c)^2/(2\kappa),&1<c<1+\kappa,\\
 0,&c\geq1+\kappa.
 \end{cases}                                            \tag{27}
\]

At the winning value (26),

\[
 J_\kappa(c_*)=
 \begin{cases}
 1+\kappa/2-\sqrt{2\kappa},&\kappa\leq1/2,\\
 1/(8\kappa),&\kappa\geq1/2.
 \end{cases}                                            \tag{28}
\]

This is strictly positive for every fixed `kappa`.  Hence at every time
`t=kappa/L` with fixed `kappa`, even the largest coordinate tail probability
is `n^{-J_kappa+o_P(1)}` and tends to zero.  A single coordinate does not
dominate on the `1/log n` time scale.  To obtain a coordinate with
order-one posterior tail mass one needs `kappa` of order `L`, namely a
constant amount of accumulated curvature, which already bounds its scalar
posterior variance by a universal constant.

Equations (25)--(28) are for the full-exposure proxy; transferring them to
the rank-one mass-preserving process would require controlling the
off-diagonal quadratic tilt and the evolving projection.  They nevertheless
identify the missing requirement for a genuine counterexample: not merely
fast angular motion, but **persistent locking before any candidate direction
accumulates constant curvature**.  The max-tail product has the former and
the full-exposure calculation rules out the naive version of the latter.

## 5. Conclusion

The endpoint estimate

\[
 \mathbb P\{\operatorname {Var}_{\mu_T}\langle
 X,u_T\rangle\leq C\}\geq c
\]

is neither proved nor refuted by the identities above.  What is proved is:

1. the exact SDE/ODE system (3)--(16), including all selection terms;
2. dimension-free control of the late-time noise and rotation terms, plus
   exponential smallness of every coherent packet of high-variance paths;
3. a fixed isotropic log-concave product and a half-mass set for which the
   angular rate is `(log n)^2` and the instantaneous selected-variance drift
   is `4 log n`, despite bounded protected one-dimensional moments;
4. a sharp large-deviation cancellation showing that this example does not
   lock onto one coordinate at subconstant accumulated curvature in the
   full-exposure proxy.

Thus any successful continuation must couple the accumulated-misalignment
identity `dot q=|P_u e|^2` to the third-moment selection term in (16).  A
bound on `||A_t||_op`, or a bound using only the marginal moments in the
currently protected direction, simply restates or misses the obstruction.
