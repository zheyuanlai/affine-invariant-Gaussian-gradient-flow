# Centroid-deficit occupation and gluing in ordinary localization

## 1. Conclusion

For ordinary stochastic localization there is an exact answer to the
occupation part of the angular-rigidity problem.  Let

\[
 g_t=\mu_t(S),\qquad
 v_t=\operatorname {Cov}_t(\mathbf 1_S,X),\qquad
 r_t=|v_t|,
\]

and, when `r_t>0`, put `u_t=v_t/r_t`, `P_t=I-u_tu_t^T`.  With

\[
 Y=X-a_t,\quad A_t=\mathbb E_tYY^T,\quad
 D_t=\mathbb E_t[(\mathbf1_S-g_t)YY^T],                 \tag{1.1}
\]

define the absolute and relative sharp-centroid defects

\[
 \Delta_t={\mathcal I(g_t)\over\sqrt t}-r_t,
 \qquad
 \eta_t={\sqrt t\,r_t\over\mathcal I(g_t)}\in[0,1].   \tag{1.2}
\]

Then, up to the usual bounded stops,

\[
\boxed{
 d\Delta_t=dM_t-\left\{
 {\|P_tD_t\|_{HS}^2\over2r_t}
 +r_tu_t^T(t^{-1}I-A_t)u_t
 +{\mathcal I(g_t)(1-\eta_t)^2\over2t^{3/2}}
 \right\}dt.}                                           \tag{1.3}
\]

All three displayed finite-variation terms are nonnegative.  In
particular, the very quantity which rotates the active direction is paid
for by decrease of the nonnegative centroid defect.  This is an integrated
statement for every strongly log-concave posterior, and requires neither a
pointwise transverse-splitting theorem nor a covariance-survival estimate.
At every fixed state, the exact two-deficit identity identifies
`Delta=D_cut+D_map`; thus (1.3) is precisely a dissipation identity for
the sum of the threshold defect and the active-marginal transport defect.

On a fixed multiplicative time window, (1.3) implies the dimension-free
rigidity statement

\[
 \mathbb E_s\int_s^\tau
 {\|P_tD_t\|_{HS}^2\over r_t^2}\,dt
 \le C_{\delta,\alpha,\Lambda}(1-\eta_s),              \tag{1.4}
\]

where `tau<=Lambda s` stops when `g_t` leaves
`[delta,1-delta]` or `eta_t` falls below `alpha>0`.  Moreover

\[
 \mathbb E_s\sup_{s\le t\le\tau}|u_t-u_s|
 \le C_{\delta,\alpha,\Lambda}\sqrt{1-\eta_s}.         \tag{1.5}
\]

Thus near equality is stable in its active direction along the actual
localization path.  This is the needed time-gluing statement.  A static
pointwise inequality is also true for Gaussian posteriors and for posteriors
which split in the active direction, but is not needed to obtain (1.4)--(1.5).
No pointwise inequality is asserted here for an arbitrary nonsplit
strongly log-concave posterior.

## 2. Exact SDEs for `v` and `D`

Let

\[
 p_t(x)={1\over Z_t}
 \exp\left(c_t\cdot x-{t\over2}|x|^2\right)p_0(x),
 \qquad dc_t=dW_t+a_tdt.                                \tag{2.1}
\]

Then

\[
 dp_t(x)=p_t(x)\langle x-a_t,dW_t\rangle,
 \qquad da_t=A_t\,dW_t.                               \tag{2.2}
\]

Write `h=1_S-g_t` and define the two third-moment contractions

\[
 \begin{aligned}
 \mathcal T_t(z)&=\mathbb E_t[YY^T\langle Y,z\rangle],\\
 \mathcal K_t(z)&=\mathbb E_t[hYY^T\langle Y,z\rangle].
 \end{aligned}                                         \tag{2.3}
\]

Direct Ito differentiation gives

\[
\boxed{
 dg_t=v_t^TdW_t,\qquad
 dv_t=D_t\,dW_t-A_tv_t\,dt.}                           \tag{2.4}
\]

\[
\boxed{\begin{aligned}
 dD_t={}&\mathcal K_t(dW_t)-A_t\langle v_t,dW_t\rangle
          -v_t(A_tdW_t)^T-(A_tdW_t)v_t^T\\
 &-\{\mathcal T_t(v_t)+A_tD_t+D_tA_t\}\,dt.
\end{aligned}}                                         \tag{2.5}
\]

There is no omitted drift in (2.4).  For (2.5), one may differentiate
`hY_iY_j` before multiplying by the likelihood.  Its own Ito drift has
zero expectation; its covariation with the likelihood is exactly
`-\mathcal T(v)-AD-DA`.  The four displayed martingale terms are respectively
the binary third central moment, the fluctuation of `g`, and the two
fluctuations of the barycenter.

The equations imply the following exact length and direction dynamics.  Set
`r=|v|`, `u=v/r`, and `P=I-uu^T`.  Then

\[
\boxed{
 dr=u^TD\,dW-r\,u^TAu\,dt
       +{\|PD\|_{HS}^2\over2r}\,dt.}                   \tag{2.6}
\]

\[
\boxed{
 du={PD\over r}\,dW-PAu\,dt
       -{PD^2u\over r^2}\,dt
       -{\|PD\|_{HS}^2\over2r^2}u\,dt.}                \tag{2.7}
\]

In particular,

\[
 d[u]_t={P_tD_t^2P_t\over r_t^2}\,dt,
 \qquad {d\over dt}\operatorname {Tr}[u]_t
 ={\|P_tD_t\|_{HS}^2\over r_t^2}.                     \tag{2.8}
\]

These formulas hold first while `r` and all relevant moments stay in a
compact set.  Convex truncation, smoothing, and removal of the stops give
the usual local statements.

## 3. The exact defect-dissipation identity

Let `z_t=Phi^{-1}(g_t)`.  The Gaussian profile satisfies

\[
 \mathcal I'(g)=-z,qquad \mathcal I''(g)=-1/\mathcal I(g).
\]

Consequently (2.4) gives

\[
 d\left({\mathcal I(g_t)\over\sqrt t}\right)
 =-{z_tv_t^T\over\sqrt t}dW_t
 -\left\{{r_t^2\over2\sqrt t\,\mathcal I(g_t)}
          +{\mathcal I(g_t)\over2t^{3/2}}\right\}dt.   \tag{3.1}
\]

Subtract (2.6), and put

\[
 C_t=t^{-1}I-A_t\succeq0.                              \tag{3.2}
\]

The positivity is Brascamp--Lieb for the `t`-strongly log-concave
posterior.  Since `r=\eta\mathcal I/\sqrt t`, the scalar terms combine
exactly as

\[
 {r\over t}-{r^2\over2\sqrt t\,\mathcal I}
 -{\mathcal I\over2t^{3/2}}
 =-{\mathcal I(1-\eta)^2\over2t^{3/2}}.                \tag{3.3}
\]

This proves (1.3), with

\[
 dM_t=\left(-{z_tv_t\over\sqrt t}-D_tu_t\right)^TdW_t. \tag{3.4}
\]

The sharp centroid inequality says `Delta_t>=0`.  Hence `Delta` is a
nonnegative local supermartingale, and for every bounded stopped interval
`[s,tau]`,

\[
\boxed{\begin{aligned}
 &\mathbb E_s\Delta_\tau
 +\mathbb E_s\int_s^\tau {\|P_tD_t\|_{HS}^2\over2r_t}dt\\
 &\quad+\mathbb E_s\int_s^\tau r_tu_t^TC_tu_tdt
 +\mathbb E_s\int_s^\tau
 {\mathcal I(g_t)(1-\eta_t)^2\over2t^{3/2}}dt
 =\Delta_s.                                             \tag{3.5}
\end{aligned}}
\]

For an initially `kappa`-strongly log-concave law, replace `t` everywhere
in (1.2)--(3.5) by `kappa+t`.  Its derivative is still one, so the proof is
unchanged.

## 4. A dimension-free occupation and direction-rigidity lemma

We first record what near equality already says at the initial state of a
window.  This is useful because it prevents the direction in the following
lemmas from being supported by a degenerate marginal.

**Lemma 4.0 (central marginal separation).**  For every
`delta in (0,1/2)` there are `epsilon_delta,c_delta>0` with the following
property.  Suppose at time `s`

\[
 g_s\in[\delta,1-\delta],\qquad
 \epsilon_s=1-{\sqrt s\,r_s\over\mathcal I(g_s)}
 \le\epsilon_\delta.                                  \tag{4.0}
\]

Let `Y=<X-a_s,u_s>`, let `rho` be its posterior law, and let `q(y)` be
the conditional probability of `S`.  Then, for every `b>0`,

\[
 \int_{|y-c|\ge b/\sqrt s}|q(y)-\mathbf1_{\{y\ge c\}}|d\rho(y)
 \le {\mathcal I(g_s)\over b}\epsilon_s,              \tag{4.0a}
\]

and

\[
 \boxed{s\operatorname {Var}_s(Y)\ge c_\delta.}       \tag{4.0b}
\]

**Proof.**  Formula (4.0a) is the threshold-error consequence of
`D_cut<=Delta_s`, with strip width `b/sqrt(s)`.  For (4.0b), choose
`a_-<z_g<a_+` so that each Gaussian interval
`[a_-,z_g]` and `[z_g,a_+]` has mass `delta/2`.  The two outer Gaussian
tails then also have mass at least `delta/2`.  If `T` is the monotone map
to `Y`, the map-defect estimate gives

\[
 T(a_+)-T(a_-)
 \ge {a_+-a_-\over\sqrt s}-{\Delta_s\over(\delta/2)^2}. \tag{4.0c}
\]

For `epsilon_delta` small enough, the right side is at least
`c_delta/sqrt(s)`.  Independent copies of `Y`, one in each outer tail,
then give `Var(Y)>=c_delta/s`.  QED.

Fix `delta in (0,1/2)`, `alpha in (0,1)`, `Lambda>1`, and a deterministic
`s>0`.  Let

\[
 \tau=\inf\{t\ge s:g_t\notin[\delta,1-\delta]
                 \text{ or }\eta_t<\alpha\}\wedge\Lambda s. \tag{4.1}
\]

We also stop at `r=0`, which is redundant before the `eta<alpha` stop.
Put `I_delta=min_{[delta,1-delta]} I>0` and

\[
 \epsilon_s=1-\eta_s={\sqrt s\,\Delta_s\over\mathcal I(g_s)}. \tag{4.2}
\]

**Lemma 4.1 (central near-equality occupation).**  Conditionally on the
state at time `s`,

\[
\boxed{
 \mathbb E_s\int_s^\tau
 {\|P_tD_t\|_{HS}^2\over r_t^2}dt
 \le C_{\delta,\alpha}\sqrt\Lambda\,\epsilon_s.}       \tag{4.3}
\]

\[
\boxed{
 \mathbb E_s\int_s^\tau(1-t\,u_t^TA_tu_t){dt\over t}
 \le C_{\delta,\alpha}\sqrt\Lambda\,\epsilon_s.}       \tag{4.4}
\]

\[
\boxed{
 \mathbb E_s\int_s^\tau(1-\eta_t)^2{dt\over t}
 \le C_{\delta}\sqrt\Lambda\,\epsilon_s.}              \tag{4.5}
\]

**Proof.**  On `[s,tau]`,

\[
 r_t\ge {\alpha I_\delta\over\sqrt t},
 \qquad {1\over r_t}\le {\sqrt{\Lambda s}\over\alpha I_\delta}. \tag{4.6}
\]

The first dissipation in (3.5), followed by (4.6), gives

\[
 \mathbb E_s\int_s^\tau {\|PD\|_{HS}^2\over r^2}dt
 \le {2\sqrt{\Lambda s}\over\alpha I_\delta}\Delta_s.
\]

Use (4.2) and the universal upper bound on `I` to get (4.3).  Similarly,

\[
 1-tu^TAu=t\,u^TCu,
\]

so (4.4) follows from the second dissipation and the same bound on `1/r`.
Finally

\[
 {(1-\eta)^2\over t}
 \le {2\sqrt{\Lambda s}\over I_\delta}
 {\mathcal I(g)(1-\eta)^2\over2t^{3/2}},
\]

which proves (4.5).  QED.

The third estimate is a genuine occupation statement.  For example, for
any `beta>0`,

\[
 \mathbb E_s\int_s^\tau
 \mathbf1_{\{1-\eta_t\ge\beta\}}{dt\over t}
 \le {C_\delta\sqrt\Lambda\over\beta^2}\epsilon_s.    \tag{4.7}
\]

Thus a central state which nearly saturates the centroid bound cannot spend
much subsequent logarithmic time far from saturation before the central
stop.

We next include the finite-variation terms in (2.7).

**Lemma 4.2 (stability of the active direction).**  Under the same stops,

\[
\boxed{
 \mathbb E_s\sup_{s\le t\le\tau}|u_t-u_s|
 \le C_{\delta,\alpha,\Lambda}
       (\sqrt{\epsilon_s}+\epsilon_s).}                \tag{4.8}
\]

In particular, for `epsilon_s<=1` the right side is
`C sqrt(epsilon_s)`.

**Proof.**  First note the dimension-free static bound

\[
 \boxed{\|D_t\|_{HS}\le {2\sqrt{g_t(1-g_t)}\over t}\le {1\over t}.} \tag{4.9}
\]

Indeed, for every symmetric `B` with `||B||HS=1`, Poincare and
`A_t\preceq t^{-1}I` give

\[
 \operatorname {Var}_t(Y^TBY)
 \le {4\over t}\mathbb E_t|BY|^2
 ={4\over t}\operatorname {Tr}(B^2A_t)\le {4\over t^2}.
\]

Covariance Cauchy followed by Hilbert--Schmidt duality proves (4.9).

The martingale in (2.7) has trace quadratic variation equal to the left
side of (4.3).  The vector-valued BDG inequality therefore bounds its
expected maximal displacement by `C sqrt(epsilon_s)`.

For the first tangential drift, `C=t^{-1}I-A` satisfies
`0\preceq C\preceq t^{-1}I`, and

\[
 |PAu|^2=|PCu|^2\le u^TC^2u\le t^{-1}u^TCu.            \tag{4.10}
\]

Consequently

\[
 \int_s^\tau|PAu|dt
 \le\left(\int_s^\tau r u^TCu\,dt\right)^{1/2}
     \left(\int_s^\tau{dt\over tr}\right)^{1/2},      \tag{4.11}
\]

and the second factor is bounded deterministically using (4.6).  Taking
conditional expectations and applying (3.5) gives
`C_{delta,alpha,Lambda} sqrt(epsilon_s)`.

For the other tangential drift, (4.9) gives

\[
 \int_s^\tau{|PD^2u|\over r^2}dt
 \le\left(\int_s^\tau{\|PD\|_{HS}^2\over r^2}dt\right)^{1/2}
     \left(\int_s^{\Lambda s}{dt\over t^2r^2}\right)^{1/2}. \tag{4.12}
\]

The last factor is at most
`sqrt(log Lambda)/(alpha I_delta)`.  Equations (4.3) and Cauchy--Schwarz
again give `C sqrt(epsilon_s)`.  The radial Ito correction in (2.7) has
total norm one half of the integral in (4.3), hence costs
`C epsilon_s`.  Combining the four contributions proves (4.8).  QED.

This proves stability along one realized path, not merely closeness of
threshold descriptions at unrelated posterior states.

As an exact rigidity corollary, if `Delta_s=0`, then (3.5) forces, almost
surely up to every bounded central stop,

\[
 P_tD_t=0,\qquad A_tu_t=t^{-1}u_t,qquad \eta_t=1.
\]

Equation (2.7) then gives `du_t=0`.  This dynamic conclusion agrees with
the static equality classification: the active marginal is Gaussian, the
cut is its threshold, and the active factor splits from the transverse
posterior.

## 5. The static pointwise question

At a fixed natural tilt `c`, differentiation gives

\[
 \nabla_cg=v,\qquad \nabla_cv=D,
 \qquad \nabla_cu={PD\over|v|}.                         \tag{5.1}
\]

Thus the proposed static quantity is exactly the squared Hilbert--Schmidt
angular derivative in tilt space.

### 5.1 Exact Gaussian posterior

For `X~N(0,t^{-1}I)`, every central set satisfies

\[
\boxed{
 \|PD\|_{HS}^2
 \le {C_\delta\over t^{3/2}}
 \left\{{\mathcal I(g)\over\sqrt t}-|v|\right\}.}      \tag{5.2}
\]

Here is the short audit of the proof.  Scale to `t=1`, rotate
`u=v/|v|` to the first coordinate, and let `H` be the upper Gaussian
`g`-halfspace.  For each transverse point `z`, let `m_0(z)` be the net
mass of flips from `H` to `S`, and `m_1(z)` their `|y-c|`-weighted mass.
The sign of a flip gives

\[
 \mathbb E m_1=\mathcal I(g)-|v|,
 \quad \mathbb E m_0=0,
 \quad \mathbb E[Zm_0]=0.                              \tag{5.3}
\]

If `M=(2pi)^{-1/2}`, layer cake gives, pointwise in `z`,

\[
 m_1(z)\ge {m_0(z)^2\over4M}.                           \tag{5.4}
\]

Indeed, if `A` is the total flip mass, then

\[
 \int|y-c|a(y)d\gamma_1(y)
 \ge\int_0^{A/(2M)}(A-2Ms)ds={A^2\over4M}.
\]

The transverse-transverse block is the second Gaussian-chaos coefficient
of `m_0`, and the transverse-active block is the first-chaos coefficient
of `m_1`.  Bessel gives

\[
 \|B\|_{HS}^2\le2\mathbb E m_0^2,
 \qquad |b|^2\le\mathbb E m_1^2
 \le C_\delta\mathbb E m_1.                            \tag{5.5}
\]

This proves (5.2).  Notice that the Gaussian second-chaos step really is a
Bessel inequality; bounding each quadratic form separately would not by
itself justify summing its coefficients.

### 5.2 A posterior split in the active direction

The same estimate holds if

\[
 \pi=\rho\otimes\nu
\]

in active/transverse coordinates, with both factors `t`-strongly
log-concave.  The one-dimensional transport identity

\[
 t^{-1/2}|z-z_g|=|T(z)-c|+|R(z)-R(z_g)|               \tag{5.6}
\]

and the two-deficit decomposition imply

\[
 \mathbb E_\nu m_0^2
 \le C_\delta\sqrt t(D_{cut}+D_{map}),
 \qquad
 \mathbb E_\nu m_1^2\le {C_\delta\over\sqrt t}D_{cut}. \tag{5.7}
\]

Poincare for the transverse factor gives, uniformly for
`||B||HS=1`,

\[
 \operatorname {Var}_\nu(W^TBW)\le4/t^2.              \tag{5.8}
\]

Hilbert--Schmidt duality then proves (5.2).  Thus non-Gaussianity of the
one-dimensional active marginal is not the obstruction; nonsplitting is.

For a general nonsplit posterior the exact two scalar defects do not by
themselves identify the variation of the conditional transverse law.  No
static extension of (5.2) is used in Sections 3--4.  Identity (1.3) is the
appropriate tilt-averaged substitute: it proves exactly the integrated
bound needed for localization.

## 6. Overlap of nearby posterior tilts

Fix the quadratic time `t` and write

\[
 \pi_c(dx)=\exp\{c\cdot x-b_t(c)\}\,
            e^{-t|x|^2/2}p_0(x)dx.                     \tag{6.1}
\]

Since `\nabla^2b_t(c)=A_{t,c}\preceq t^{-1}I`, two same-time tilts have the
dimension-free overlap bounds

\[
\boxed{\begin{aligned}
 \chi^2(\pi_{c+h}\|\pi_c)
 &=\exp\{b_t(c+2h)-2b_t(c+h)+b_t(c)\}-1\\
 &\le e^{|h|^2/t}-1.
\end{aligned}}                                         \tag{6.2}
\]

\[
\boxed{\begin{aligned}
 \int\sqrt{d\pi_c\,d\pi_{c+h}}
 &=\exp\left\{b_t(c+h/2)-{b_t(c)+b_t(c+h)\over2}\right\}\\
 &\ge e^{-|h|^2/(8t)}.
\end{aligned}}                                         \tag{6.3}
\]

These estimates safely transfer the threshold error of a fixed set between
two tilts separated by `O(sqrt(t))`.  They do not, however, glue an entire
multiplicative time interval.  Even in the Gaussian model,

\[
 \operatorname {Aff}\big(N(0,s^{-1}I),N(0,t^{-1}I)\big)
 =\left({2\sqrt{st}\over s+t}\right)^{n/2}.             \tag{6.4}
\]

For `t=(1+theta)s`, this is `exp(-Theta(n theta^2))`.  Also a typical
Brownian tilt increment in `n` dimensions has squared norm of order
`n dt`, so applying (6.2) step by step incurs the same trace loss.  Thus
full posterior overlap across times is not dimension free.  The
defect-dissipation identity avoids that loss because it tracks only the
fixed set's active centroid and exactly cancels all irrelevant trace
directions.

Equations (4.3) and (4.8) are therefore the correct gluing statements:
nearby states along the actual path share an active direction in quadratic
mean even though their full densities can have exponentially small
high-dimensional overlap.

## 7. Mandatory model tests

### 7.1 Gaussian halfspace

Let the posterior be `N(m,lambda^{-1}I)` and let `S` be a halfspace with
unit normal `u`.  If `z=Phi^{-1}(g)`, then

\[
 v={\mathcal I(g)\over\sqrt\lambda}u,
 \qquad
 D=-{z\mathcal I(g)\over\lambda}uu^T.                  \tag{7.1}
\]

Thus `Delta=0`, `PD=0`, `A=lambda^{-1}I`, and all three dissipations in
(1.3) vanish.  The active direction is exactly constant.  For an initially
`kappa`-Gaussian law one must use `lambda=kappa+t`; using only the weaker
curvature `t` correctly leaves a positive scalar defect but still no
angular noise.

### 7.2 Gaussian parity and radial cuts

For `S={x_1x_2>=0}` under `N((a,b),I)`, put `p=Phi(a)` and `q=Phi(b)`.
Then

\[
 g=1-p-q+2pq,qquad
 v=(\varphi(a)(2q-1),\varphi(b)(2p-1)).                 \tag{7.2}
\]

Along the central line `b=0`,

\[
 g={1\over2},\quad
 v=\varphi(0)(2\Phi(a)-1)e_2,
 \quad D_{12}=D_{21}=2\varphi(a)\varphi(0),             \tag{7.3}
\]

and the other entries of `D` vanish.  Therefore

\[
 \mathcal I(1/2)-|v|=2\varphi(0)\Phi(-|a|),
 \qquad \|PD\|_{HS}^2=4\varphi(a)^2\varphi(0)^2.       \tag{7.4}
\]

The two far arms `(|a|,0)` and `(0,|a|)` have nearly saturated, orthogonal
active directions, but they are separated by a region containing the
state `v=0`.  Their Gaussian tilt overlap is exponentially small in
`a^2`.  Hence this example forbids a single global direction for all
near-equality tilts, while fully agreeing with local occupation gluing.

For the radial cut `S={|x|<=R}`, write
`g(m)=G(|m|)` under `N(m,I)`.  At `m ne0`,

\[
 v=G'(|m|){m\over|m|},qquad
 D=G''ee^T+{G'\over|m|}(I-ee^T).                        \tag{7.5}
\]

Consequently

\[
 {\|PD\|_{HS}^2\over|v|^2}={n-1\over|m|^2}.            \tag{7.6}
\]

The Gaussian static estimate forces a central near-equality radial state
to have `|m|^2` much larger than `n`; the collective rotation of its
`n-1` tangential directions is paid exactly by its centroid defect.  At
`m=0`, symmetry gives `v=0`, so the state is outside the critical regime.

### 7.3 Product exponentials and the maximum cut

Let `Z_i` be iid rate-one exponentials, `X_i=Z_i-1`, and choose

\[
 (1-e^{-L})^n={1\over2},\qquad q=e^{-L}.
\]

For `S={max_i Z_i>=L}`, exchangeability gives

\[
 r_0={\sqrt nLq\over2(1-q)},
 \qquad
 P_0D_0P_0={L^2q\over2(1-q)^2}P_0.                     \tag{7.7}
\]

Hence

\[
 {\|P_0D_0\|_{HS}^2\over r_0^2}
 ={n-1\over n}{L^2\over(1-q)^2}
 =(1+o(1))(\log n)^2.                                  \tag{7.8}
\]

This large angular rate occurs at curvature zero.  At a small positive
ordinary-localization time `s`, its normalized centroid coefficient begins
as

\[
 {\sqrt s\,r_0\over\mathcal I(1/2)}
 =O\left({\sqrt s\log n\over\sqrt n}\right),           \tag{7.9}
\]

so the symmetric maximum state is very far from centroid saturation.  The
large-deviation calculation for the full-observation proxy is sharper: at
`s=kappa/L`, the winning coordinate still has posterior tail mass
`n^{-J_kappa+o(1)}` with `J_kappa>0` for every fixed `kappa`.  Thus the
maximum cut has fast angular exploration but no subconstant-curvature
locking.  It is a stress test passed by (1.3): the factor `Delta_s` is large
precisely where (7.8) is large.

More generally, for an iid centered factor of variance `sigma^2`, with
`F(L)^n=1/2`, let `sigma_L^2=Var(Z|Z<L)` and
`m_L=E[(Z-EZ)1_{Z>=L}]`.  The exact exchangeable formulas are

\[
 r={\sqrt n\,m_L\over2F(L)},
 \qquad
 P D P={\sigma^2-\sigma_L^2\over2}P.                   \tag{7.10}
\]

They apply in particular to Gaussian-damped exponential factors, which
are strongly log-concave.

### 7.4 Uniform simplex

Let `U` be uniform on the `N-1` dimensional probability simplex and

\[
 X=\sqrt{N(N+1)}\left(U-{\mathbf1\over N}\right).
\]

Split the indices into two equal groups and let `S` record which group
contains `argmax_i U_i`.  This is a balanced cut.  Since
`E max_i U_i=H_N/N`, permutation symmetry gives

\[
 \left|\mathbb E[(2\mathbf1_S-1)X]\right|
 ={\sqrt{N+1}(H_N-1)\over N-1}
 =\Theta\left({\log N\over\sqrt N}\right).             \tag{7.11}
\]

Thus `|v|` is one half of (7.11) and tends to zero.  The uniform simplex
has curvature zero in its interior, so its symmetric max-cell partition is
not a near-equality strongly-log-concave posterior state.  Adding the
ordinary quadratic tilt makes it `t`-strongly log-concave, but (7.11) by
itself gives no control at a later dimension-dependent time.  What it does
show exactly is that the rapid-decoding argument for an atomic simplex
cannot be transferred using vertex covariance alone.  An atomic regular
simplex is not log-concave and is therefore not an admissible
counterexample.

The parallel-facet halfspace cut is the complementary calibration.  By
simplex symmetry its centroid direction is a fixed vertex direction and
the active-transverse column of `D` vanishes at the symmetric state.
The transverse-transverse block need not vanish: the conditional simplex
slice scales with the active barycentric coordinate.  This is precisely a
nonsplit effect, carried by the contact facets.  At positive quadratic
time it is covered by the integrated identity (1.3), not by an assumed
product decomposition.

## 8. What the lemma does and does not prove

The exact result is the following local-to-pathwise implication:

\[
 \text{central near centroid equality at time }s
 \quad\Longrightarrow\quad
 \text{dimension-free directional stability on }[s,\Lambda s]
\]

up to the first loss of centrality or near equality.  The Gaussian-profile
drift is

\[
 {\mathcal I(g_t)\over2\sqrt t}(1-\eta_t^2)
 ={1+\eta_t\over2}\Delta_t,                            \tag{8.1}
\]

so a small total profile drift puts most weighted central occupation time
in states with small `Delta`; Lemmas 4.1--4.2 then glue each such state
forward without a dimension factor.

This is not a seed.  It does not show that the process reaches a central
near-equality state at a universal deterministic time, and it does not
bound the probability of remaining central.  Either assertion would be a
global Cheeger/KLS input.  What has been removed is the separate angular
compatibility gap: once a central near-equality state occurs, its active
direction cannot split or rotate cheaply along ordinary localization.
