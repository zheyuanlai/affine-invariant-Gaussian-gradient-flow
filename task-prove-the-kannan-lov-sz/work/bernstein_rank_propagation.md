# Bernstein rank propagation: an exact temporal lemma and a scale obstruction

## 0. Verdict

Let `mu` be isotropic and log-concave, let `K=C_P(mu)`, let `S` be a
balanced near-Cheeger set, and use the heat notation of
`heatflow_bernstein.md`.  Put

\[
 B(s)=\sqrt s\,A_*(s),\qquad
 R(s)=s\int q_s\nabla h_s\nabla h_s^T,
 \qquad h_s=2\arcsin\sqrt{g_s}.                    \tag{0.1}
\]

There are two exact evolution formulas.  First,

\[
 B'(s)=D_{\rm Hess}(s)+D_{\rm curv}(s)+D_{\rm eik}(s),       \tag{0.2}
\]

where all three terms are explicit nonnegative squares.  Second, if
`M(s)=R(s)/s` is the binary Fisher matrix, then `tr M(s)` has its own exact
three-square dissipation identity.  The matrix `M(s)` itself is not
monotone in Loewner order: a curvature anticommutator can rotate its
eigenspaces.

The main positive result of this report is a nonlocal rank-propagation
lemma.  Suppose that at time `s` a set `G` of boundary states is central,
has `e_s>=1-tau`, and carries a high-rank Fisher matrix.  For `r=Ls`, let

\[
 d_{s,r}={B(r)-B(s)\over \sqrt s\,J(s)}.                       \tag{0.3}
\]

There is a reverse-boundary coupling for which

\[
 E|U_{\lambda,s}-U_{\lambda,r}|
 \le 20\left\{\sqrt{L d_{s,r}/\lambda}
                     +d_{s,r}/\lambda\right\},
 \qquad
 U_{\lambda,t}={\sqrt t\nabla z_t
                   \over\sqrt{\lambda+t|\nabla z_t|^2}}.     \tag{0.4}
\]

Consequently, if the right side is smaller than a fixed fraction of the
initial marked boundary mass, then

\[
 \operatorname {tr}R(r)\ge c\sqrt\alpha,\qquad
 \operatorname {rank}_{\rm eff}R(r)\ge c r\sqrt\alpha,        \tag{0.5}
\]

starting from the fixed-scale theorem at `s=alpha K`.  In particular,
rank cannot disappear without paying

\[
 B(r)-B(s)\ge {c\over L}\sqrt s\,J(s)
              \ge {c\alpha\over L}\sqrt K.                   \tag{0.6}
\]

This is a genuine propagation theorem with universal constants once
`alpha`, `delta`, and `tau` are frozen.

It does not close KLS.  Near-Cheeger global minimality gives, for
`r=cK`,

\[
 B(r)\le Cc^{3/2}\sqrt K.                                    \tag{0.7}
\]

Combining (0.6) and (0.7) prevents collapse only while

\[
                         c^{5/2}\lesssim\alpha^2.              \tag{0.8}
\]

Even if one grants perfect rank persistence and also grants a global
conversion of all wedge-Poincare energy into Bernstein energy, the forced
integrated charge is only

\[
              c\sqrt\alpha\,c^{3/2}\sqrt K,                  \tag{0.9}
\]

which is smaller than (0.7) by `sqrt(alpha)`.  If the natural Fisher trace
grows like `sqrt(r/K)`, the forced charge is `c c^2 sqrt(K)`, still smaller
than (0.7) by `sqrt c`.  Balanced Gaussian balls realize exactly the latter
power while their Fisher rank remains `n` at every small scale.  Orthogonal
polyhedral phase patterns, including one-sided exponential maximum boxes,
have an even larger corner charge of order `c^{3/2}`, exactly matching the
available near-minimality budget.  Thus persistence is not the missing
ingredient: a new extremality/phase charge is still required.

No covariance-process estimate, KLS-strength Poincare statement, or
pointwise angular-stability modulus is used below.

## 1. Heat notation and the complete Bernstein derivative

For `s>0`, let

\[
 q_s=P_s\mu,\qquad q_sg_s=P_s(1_S\mu),\qquad
 z_s=\Phi^{-1}(g_s),\qquad \rho_s=q_s\varphi(z_s).              \tag{1.1}
\]

Write

\[
 a=\nabla z,\quad K_z=\nabla^2z,\quad \ell=\log q,\quad
 e=s|a|^2,\quad
 J(s)=\int\rho_s,\quad A_*(s)=\int\rho_s(1-e).                 \tag{1.2}
\]

Gaussian convolution makes all quantities smooth.  The posterior is
`s^{-1}`-strongly log-concave on the affine support of `mu`, and the sharp
Gaussian centroid inequality gives

\[
 0\le e\le1,\qquad -s^{-1}I\preceq\nabla^2\ell\preceq0.        \tag{1.3}
\]

The exact heat equations are

\[
 \partial_s z={1\over2}\Delta z+\langle\nabla\ell,\nabla z\rangle
                    -{z\over2}|\nabla z|^2,
 \qquad
 \partial_s\rho={1\over2}\Delta\rho+{1\over2}\rho|a|^2.    \tag{1.4}
\]

Integration by parts gives

\[
 {d\over ds}{J(s)\over\sqrt s}
        =-{A_*(s)\over2s^{3/2}}=-{B(s)\over2s^2}.              \tag{1.5}
\]

The pointwise Bochner calculation, with no discarded term, yields

\[
\boxed{
\begin{aligned}
 B'(s)={}&s^{3/2}\int\rho_s\|K_z\|_{HS}^2\\
 &+2s^{3/2}\int\rho_s
          \big|(-\nabla^2\ell)^{1/2}a\big|^2\\
 &+{1\over2\sqrt s}\int\rho_s(1-e)^2.
\end{aligned}}                                                   \tag{1.6}
\]

Thus the complete list of Bernstein squares is

\[
\begin{aligned}
 D_{\rm Hess}(s)&=s^{3/2}\int\rho_s\|\nabla^2z_s\|_{HS}^2,\\
 D_{\rm curv}(s)&=2s^{3/2}\int\rho_s
       |(-\nabla^2\log q_s)^{1/2}\nabla z_s|^2,\\
 D_{\rm eik}(s)&={1\over2\sqrt s}\int\rho_s
       (1-s|\nabla z_s|^2)^2.
\end{aligned}                                                    \tag{1.7}
\]

In particular, `B` is nondecreasing and, for `0<s<r`,

\[
 B(r)-B(s)=\int_s^r
       (D_{\rm Hess}+D_{\rm curv}+D_{\rm eik})(t)\,dt.         \tag{1.8}
\]

For nonsmooth `mu` and Borel `S`, first replace the label by
`epsilon+(1-2epsilon)1_S`, use spatial cutoffs, and send `epsilon` to zero.
Every term in (1.7) is nonnegative, so Fatou gives (1.8) as an equality of
nonnegative Radon measures and all inequalities below.  Ambient Gaussian
convolution also makes the formulas valid for lower-dimensional support.

## 2. A second exact identity: binary Fisher dissipation

Put

\[
 h=2\arcsin\sqrt g,\qquad p=\nabla h,\qquad H=\nabla^2h,\qquad
 c(h)=\cot h.                                                   \tag{2.1}
\]

Since `g=(1-cos h)/2`, the quotient heat equation for `g` is equivalent to

\[
 \partial_sh={1\over2}\Delta h+\langle\nabla\ell,p\rangle
                         +{c(h)\over2}|p|^2.                    \tag{2.2}
\]

Define the matrix Fisher information

\[
 \mathsf M(s)=\int q_s pp^T,\qquad R(s)=s\mathsf M(s).          \tag{2.3}
\]

Differentiating (2.3), using `partial_s q=(1/2)Delta q`, and integrating
the two Laplacians by parts gives the exact matrix equation

\[
\boxed{
\begin{aligned}
 \mathsf M'(s)=\int q_s\{&-(H-cpp^T)^2-|p|^2pp^T\\
 &+(\nabla^2\ell\,p)p^T+p(\nabla^2\ell\,p)^T\}.
\end{aligned}}                                                   \tag{2.4}
\]

Here `(H-cpp^T)^2` means the product of the symmetric matrix with itself.
For completeness, differentiating (2.2) gives

\[
 p_s={1\over2}\Delta p+H\nabla\ell+\nabla^2\ell\,p
       +cHp-{1\over2}\csc^2(h)|p|^2p.                           \tag{2.5}
\]

The heat derivative of `pp^T` contributes `HH^T`.  Integration by parts
turns the two `Delta p` terms into `-2HH^T` and cancels both
`H nabla ell` terms.  Finally

\[
 -H^2+c\{(Hp)p^T+p(Hp)^T\}-\csc^2(h)|p|^2pp^T
 =-(H-cpp^T)^2-|p|^2pp^T,                                    \tag{2.6}
\]

because `csc^2 h-cot^2 h=1`.  This proves (2.4).

Taking the trace produces another exact sum of nonnegative squares:

\[
\boxed{
 -{d\over ds}\operatorname {tr}\mathsf M(s)
 =\int q_s\left\{
 \|H-cpp^T\|_{HS}^2+|p|^4
 +2|(-\nabla^2\ell)^{1/2}p|^2\right\}.}                        \tag{2.7}
\]

Thus total binary Fisher information decreases.  The matrix does not in
general decrease in Loewner order.  Even when `nabla^2 ell preceq0`, the
anticommutator

\[
 (\nabla^2\ell\,p)p^T+p(\nabla^2\ell\,p)^T                    \tag{2.8}
\]

can have one positive eigenvalue when `p` is not an eigenvector of
`nabla^2 ell`.  Consequently (2.7) controls loss of trace but does not by
itself control rotation or effective rank.

At a posterior state let

\[
 v=\operatorname {Cov}(1_S,X\mid Y),\qquad G=g(1-g).
\]

Then

\[
 R(s)={1\over s}E\,{vv^T\over G}.                              \tag{2.9}
\]

Conditional covariance Cauchy--Schwarz and total covariance give the
crucial operator cap

\[
                         \boxed{R(s)\preceq s^{-1}I.}           \tag{2.10}
\]

This cap, rather than matrix monotonicity, is what makes temporal rank
propagation possible.

## 3. Reverse boundary transport

Fix `0<s<r` and put `L=r/s`.  Set

\[
 x_t=\sqrt t\nabla z_t,\qquad
 U_{\lambda,t}={x_t\over\sqrt{\lambda+|x_t|^2}},\qquad
 0<\lambda\le1.                                                \tag{3.1}
\]

The weight `rho_t` is the alive density of an exact reverse killed
diffusion.  Starting at time `r` with unnormalised density `rho_r`, run

\[
 dY_\tau=\nabla\log\rho_{r-\tau}(Y_\tau)d\tau+dW_\tau          \tag{3.2}
\]

and kill at rate `|nabla z_{r-tau}|^2/2`.  Its alive density at heat time
`t=r-tau` is `rho_t`.  The law of the pair of endpoints on surviving paths
therefore defines a measure `Gamma_{s,r}` of mass `J(s)` such that

\[
 (\operatorname {pr}_s)_\#\Gamma_{s,r}=\rho_s\,dy,\qquad
 (\operatorname {pr}_r)_\#\Gamma_{s,r}\le\rho_rdy.             \tag{3.3}
\]

The inequality in the second marginal is just survival probability at
most one.  No Poincare inequality is used.

Let

\[
 D=B(r)-B(s),\qquad d={D\over\sqrt sJ(s)}.                     \tag{3.4}
\]

The equation

\[
 \left(\partial_t-{1\over2}\Delta-
       \langle\nabla\log\rho_t,\nabla\rangle\right)x_t
 =\nabla^2\ell_t x_t+{1-|x_t|^2\over2t}x_t                     \tag{3.5}
\]

and the chain rule for (3.1) split the motion of `U_lambda` into a
martingale, a curvature drift, a scalar drift, and a Hessian drift.  The
derivative bounds

\[
 \|DU_\lambda\|\le\lambda^{-1/2},\qquad
 \|D^2U_\lambda\|\le6\lambda^{-1}                              \tag{3.6}
\]

and (1.7) give, after integrating over the surviving paths,

\[
\begin{aligned}
 E_\Gamma\langle M\rangle&\le {D\over\lambda\sqrt s},\\
 E_\Gamma\left(\int|b_{\rm curv}|\right)^2
   +E_\Gamma\left(\int|b_{\rm scal}|\right)^2
   &\le {2(L-1)D\over\lambda\sqrt s},\\
 E_\Gamma\int|b_{\rm Hess}|&\le {3D\over\lambda\sqrt s}.
\end{aligned}                                                    \tag{3.7}
\]

For example, `|nabla U_lambda|^2<=t||K_z||^2/lambda`, whose ratio to
the Hessian square in (1.7) is at most `1/(lambda sqrt(s))`.  For the two
finite-variation squares, Cauchy--Schwarz contributes the interval length
`r-s=(L-1)s`; the curvature bound

\[
 \|\nabla^2\ell\,u\|^2
 \le t^{-1}\{-\langle\nabla^2\ell\,u,u\rangle\}               \tag{3.8}
\]

follows from `-t^{-1}I preceq nabla^2 ell preceq0`.  These observations
prove every line of (3.7) directly from (1.7).

Divide by the path mass `J(s)`, use Ito isometry and Cauchy--Schwarz, and
absorb the four numerical terms.  One obtains the promised estimate

\[
\boxed{
 {1\over J(s)}\int |U_{\lambda,s}(y)-U_{\lambda,r}(y')|
                  \,d\Gamma_{s,r}(y,y')
 \le20\left\{\sqrt{Ld/\lambda}+d/\lambda\right\}.}            \tag{3.9}
\]

The constant 20 is deliberately loose.  The dyadic calculation gives 12;
the displayed value covers every `L>=1` in (3.7).

## 4. Explicit marked-rank propagation

Let `G_s` be a Borel set of states at time `s` such that

\[
 \delta\le g_s\le1-\delta,\qquad e_s\ge1-\tau,\qquad 0<\tau<1.
                                                                    \tag{4.1}
\]

Define

\[
 R_{G_s}(s)=\int_{G_s}q_s\,s\nabla h_s\nabla h_s^T,\qquad
 M_\delta=\max_{\delta\le u\le1-\delta}
             {I(u)\over u(1-u)}.                                \tag{4.2}
\]

Since

\[
 q_s s|\nabla h_s|^2
 =\rho_s e_s{I(g_s)\over g_s(1-g_s)},                            \tag{4.3}
\]

we have

\[
 \rho_s(G_s)\ge {\operatorname {tr}R_{G_s}(s)\over M_\delta}.
                                                                    \tag{4.4}
\]

Put

\[
 m={\rho_s(G_s)\over J(s)},\qquad
 a_0^2={1-\tau\over\lambda+1}.                                  \tag{4.5}
\]

On `G_s`, `|U_{lambda,s}|^2>=a_0^2`.  Let

\[
 \epsilon=20\{\sqrt{Ld/\lambda}+d/\lambda\}.                    \tag{4.6}
\]

If

\[
                         \epsilon\le {ma_0^2\over4},             \tag{4.7}
\]

then `||U||<=1` and `||aa^T-bb^T||<=2|a-b|` give

\[
 {1\over J(s)}\int_{\{y\in G_s\}}|U_{\lambda,r}(y')|^2
             d\Gamma_{s,r}(y,y')
 \ge {ma_0^2\over2}.                                             \tag{4.8}
\]

Because `|U_lambda|^2=e/(lambda+e)`, one has
`e>=lambda|U_lambda|^2`.  Also the elementary Gaussian-profile inequality
`I(u)>=c_Gu(1-u)`, with `c_G=sqrt(2/pi)`, gives

\[
\begin{aligned}
 \operatorname {tr}R(r)
 &=\int\rho_r e_r{I(g_r)\over g_r(1-g_r)}\\
 &\ge c_G\int\rho_r e_r\\
 &\ge {c_G\lambda ma_0^2\over2}J(s)
  ={c_G\lambda a_0^2\over2}\rho_s(G_s).
\end{aligned}                                                    \tag{4.9}
\]

The third line uses the second-marginal domination in (3.3), restricted to
paths marked by `y in G_s`.  Combining (2.10) and (4.9) proves:

**Proposition 1 (marked rank propagation).**  Under (4.1)--(4.7),

\[
\boxed{
\begin{aligned}
 \operatorname {tr}R(r)&\ge
 {c_G\lambda a_0^2\over2M_\delta}
       \operatorname {tr}R_{G_s}(s),\\
 \operatorname {rank}_{\rm eff}R(r)&\ge
 {c_G\lambda a_0^2\over2M_\delta}
       r\operatorname {tr}R_{G_s}(s).
\end{aligned}}                                                   \tag{4.10}
\]

This conclusion concerns the full endpoint Fisher matrix, not a
non-measurable population of marked paths.  It remains valid if all marked
paths arrive in one small part of observation space, because the global
operator cap (2.10) prevents that part from carrying a large eigenvalue.

### 4.1 Instantiation with the fixed-scale theorem

At `s=alpha K`, the fixed-scale result supplies

\[
 \operatorname {tr}R_{G_s}(s)\ge{\sqrt\alpha\over8\pi},
 \qquad J(s)\ge c_0\sqrt\alpha,                                \tag{4.11}
\]

where `G_s` satisfies (4.1).  It also gives `H(s)=J(s)/sqrt(s)<=p` and

\[
 p\le{1\over(1-\beta)\sqrt K},                                  \tag{4.12}
\]

so

\[
 J(s)\le{\sqrt\alpha\over1-\beta},\qquad
 m\ge {1-\beta\over8\pi M_\delta}.                              \tag{4.13}
\]

Freeze `lambda=1/2`, `beta<=beta_0<1`, and `tau<=tau_0<1`, replace the
right side of (4.13) by the uniform lower bound

\[
 m_*={1-\beta_0\over8\pi M_\delta},
\]

and put

\[
 \bar a^2={1-\tau_0\over3/2},\qquad
 d_*:=\min\left\{
 {\lambda m_*^2\bar a^4\over25600},
 {\lambda m_*\bar a^2\over160}\right\}.                      \tag{4.14}
\]

Then `d<=d_*/L` implies (4.7).  Proposition 1 and (4.11) give

\[
\boxed{
 \operatorname {tr}R(r)\ge c_*\sqrt\alpha,\qquad
 \operatorname {rank}_{\rm eff}R(r)
       \ge c_*r\sqrt\alpha,}                                   \tag{4.15}
\]

where

\[
 c_*={c_G\lambda\bar a^2\over16\pi M_\delta}>0.               \tag{4.16}
\]

Equivalently, failure of the trace conclusion in (4.15) forces

\[
\boxed{
 B(r)-B(s)>{d_*\over L}\sqrt sJ(s)
             \ge {d_*c_0\alpha\over L}\sqrt K.}               \tag{4.17}
\]

All constants are independent of `n`, `mu`, `K`, and the near-minimizer.
The dependence of `d_*` on `alpha` occurs only through the once-for-all
central cutoff `delta_alpha`; after `alpha` is frozen it is universal.

## 5. What near-minimality can pay

Let `p=P_mu(S)`, `p<=psi/2+epsilon_0`, and `beta=epsilon_0/p`.  The exact
heat-profile bounds are

\[
 H(t)={J(t)\over\sqrt t}\le p,\qquad
 H(t)\ge\psi\{1/2-2U(t)\},\qquad
 U(t)\le{\sqrt t\,p\over c_G}.                                  \tag{5.1}
\]

Since `p(1-beta)<=psi/2`,

\[
 p-H(2t)\le p\{\beta+4U(2t)\}.                                 \tag{5.2}
\]

Monotonicity of `B` and (1.5) imply

\[
 H(t)-H(2t)=\int_t^{2t}{B(u)\over2u^2}du
              \ge {B(t)\over4t}.                               \tag{5.3}
\]

Therefore

\[
\boxed{
 B(t)\le4tp\left\{\beta+{4\sqrt{2t}\,p\over c_G}\right\}.}  \tag{5.4}
\]

At `t=cK`, use (4.12) to obtain

\[
\boxed{
 B(cK)\le {4c\over1-\beta}
 \left\{\beta+{4\sqrt{2c}\over c_G(1-\beta)}\right\}\sqrt K.}
                                                                    \tag{5.5}
\]

The near-minimizer error may be sent to zero after all constants are fixed;
the structural scale is consequently

\[
                         B(cK)\lesssim c^{3/2}\sqrt K.          \tag{5.6}
\]

Take `r=cK` in (4.17), so `L=c/alpha`.  Rank collapse costs at least

\[
 {d_*c_0\alpha^2\over c}\sqrt K.                               \tag{5.7}
\]

Comparison with (5.6) excludes collapse only under

\[
 c^{5/2}\le c_{\rm prop}\alpha^2,
 \qquad
 c_{\rm prop}>0\text{ universal after }\alpha\text{ is fixed}. \tag{5.8}
\]

Thus the exact temporal argument genuinely propagates the fixed-scale rank
to a larger heat scale, for example to

\[
 c=c_1\alpha^{4/5}                                                   \tag{5.9}
\]

with a sufficiently small fixed `c_1`.  It cannot reach a scale independent
of `alpha` using the near-minimality budget alone.

### 5.1 The optimistic integrated-energy audit

It is useful to grant more than has been proved and check the powers.  Assume
throughout `[alpha K,cK]` that

1. `tr R(t)>=a sqrt(alpha)` and `||R(t)||op<=1/t`;
2. all wedge-Poincare derivative energy occurs on central states;
3. it converts to Bernstein energy with the best scale factor
   `B'(t)>=b sqrt(t) int q_t||nabla(sqrt(t)nabla h_t)||^2`.

The second and third assumptions are deliberately **stronger than the
proved inputs** and are false without a cutoff theorem.  They are used only
for the no-go power audit in this subsection, not in Proposition 1 or in
(5.4)--(5.8).  Thus any failure under these assumptions is decisive for
this scalar chain.  The
wedge-Poincare inequality and `C_P(q_t)<=K+t` give

\[
 \int q_t\|\nabla(\sqrt t\nabla h_t)\|_{HS}^2
       \ge {a\sqrt\alpha\over2(K+t)}.                            \tag{5.10}
\]

Consequently

\[
 B(cK)-B(\alpha K)
 \ge c_2{\sqrt\alpha\over K}
       \int_{\alpha K}^{cK}\sqrt t\,dt
 \ge c_3\sqrt\alpha\,c^{3/2}\sqrt K                            \tag{5.11}
\]

when `c>=2alpha`.  This is smaller than (5.6) by the factor `sqrt(alpha)`.
If instead the natural boundary mass grows and one grants

\[
                         \operatorname {tr}R(t)\gtrsim\sqrt{t/K}, \tag{5.12}
\]

the same calculation gives only

\[
                         B(cK)-B(\alpha K)\gtrsim c^2\sqrt K,    \tag{5.13}
\]

which is smaller than (5.6) by `sqrt c`.  Effective rank removes the mean
term in wedge Poincare, but it cannot multiply (5.10) by the rank: the exact
wedge inequality is proportional to `tr R`, not to
`rank_eff(R) tr R`.

Equations (5.11)--(5.13) show that even ideal temporal persistence cannot
produce the missing contradiction.  A new charge must use global
near-extremality in a way not expressible as the scalar Dirichlet energy of
the Fisher feature.

## 6. Model audits

### 6.1 Gaussian halfspaces

For `mu=N(0,I)` and `S={x_1>=0}`,

\[
 z_s(y)={y_1\over\sqrt{s(1+s)}},\qquad
 e_s={1\over1+s},\qquad \nabla^2z_s=0,\qquad
 \nabla^2\log q_s=-{I\over1+s}.                                  \tag{6.1}
\]

Moreover

\[
 J(s)=I_0\sqrt{s\over1+s},\quad
 A_*(s)=I_0{s^{3/2}\over(1+s)^{3/2}},\quad
 B(s)=I_0{s^2\over(1+s)^{3/2}},                                  \tag{6.2}
\]

and

\[
 B'(s)=I_0{2s+s^2/2\over(1+s)^{5/2}}.                            \tag{6.3}
\]

The curvature and eikonal squares in (1.7) add exactly to (6.3), while the
Hessian square vanishes.  The Fisher matrix has rank one at every time.
This is the affine equality branch and verifies that no rank charge should
be present.

### 6.2 Balanced Gaussian balls: persistent rank with only `s^2` charge

Let `S={|x|<=R_n}` have Gaussian mass one half.  Rotational symmetry gives

\[
                         R(s)=r_s I_n,\qquad
 \operatorname {rank}_{\rm eff}R(s)=n                              \tag{6.4}
\]

for every `s>0` for which `r_s>0`.  Thus there is no rank collapse at all.

The dimensional normalization of the weighted perimeter is

\[
 p_n=(2\pi)^{-n/2}e^{-R_n^2/2}
       |\mathbb S^{n-1}|R_n^{n-1}
 ={R_n^{n-1}e^{-R_n^2/2}\over
       2^{n/2-1}\Gamma(n/2)}.                                    \tag{6.5}
\]

This is exactly the density of the chi distribution with `n` degrees of
freedom at its median.  The chi central limit theorem (or Stirling applied
at the median) gives

\[
 R_n^2/n\longrightarrow1,\qquad p_n\longrightarrow1/\sqrt\pi.   \tag{6.6}
\]

A standard tubular-coordinate calculation gives the exact leading powers
as `s downarrow0`.  If

\[
 C_F=\int_{-\infty}^{\infty}
 {\varphi(z)^2\over\Phi(z)(1-\Phi(z))}\,dz\in(0,\infty),          \tag{6.7}
\]

then

\[
 \operatorname {tr}R(s)=C_Fp_n\sqrt s+o(\sqrt s),                \tag{6.8}
\]

and (1.6) gives

\[
 B'(s)=p_n\left(2+{n-1\over R_n^2}\right)s+o(s).                 \tag{6.9}
\]

Indeed, put `y=(R_n+sqrt(s)t)\theta`.  Gaussian conditioning gives
`X|Y=y ~ N(y/(1+s),sI/(1+s))`, and expansion of the radial norm gives

\[
 g_s(y)=\Phi(-t)+\sqrt s\,\varphi(t)
 \left(R_n-{n-1\over2R_n}\right)+O(s)
\]

locally in `C^2` in the tubular variables (differentiate the Gaussian
integral under the sign).  Consequently
`z_s(y)=-t+sqrt(s)(R_n-(n-1)/(2R_n))+O(s)`.  Thus
`nabla z_s=-s^{-1/2}\theta+O(s^{1/2})`, so
`1-e_s=O(s)`; the eikonal square is `O(s^2)=o(s)` after boundary
integration.  The tangential Hessian has `n-1` eigenvalues
`1/(R_n sqrt(s))+O(s^{1/2})`, and `rho_sdy/sqrt(s)` converges to
`varphi(t)dt` times Gaussian surface measure.  Hence the curvature
square contributes `2p_ns` and the Hessian square contributes
`p_n(n-1)s/R_n^2`.  Gaussian tails dominate the local expansion,
proving (6.8)--(6.9).  Since `B(0)=0`, integration gives the more
explicit audit

\[
 B(c)={p_n\over2}\left(2+{n-1\over R_n^2}\right)c^2+o(c^2)
       \asymp c^2.                                                \tag{6.10}
\]

for each fixed `n` as `c downarrow0`, with leading coefficients bounded
above and below uniformly as `n to infinity`.  Equivalently, a diagonal
choice `c_n downarrow0` realizes (6.10) with uniform constants.  No
uniform-in-`n` remainder for a fixed positive `c` is being asserted.
This realizes the optimistic power (5.13).  Perfect rank persistence
therefore does not improve `c^2` to the `c^{3/2}` scale needed to
contradict (5.6).

### 6.3 Symmetric Laplace halfline

For `dmu=(1/2)e^{-|x|}dx` and `S=[0,infinity)`, every superlevel of every
smoothed label is a halfline.  The exact one-dimensional isoperimetric
profile is `min(u,1-u)`, so the coarea deficit is zero at every heat time.
Nevertheless `h_s'=d(2arcsin sqrt(g_s))/dy` is nonconstant, hence

\[
                         \int q_s|h_s''|^2>0.                       \tag{6.11}
\]

The Fisher rank is one.  This exact near-minimizer shows that longitudinal
and curvature/eikonal Bernstein dissipation can be positive while global
minimality has zero scalar deficit.  Such energy must be removed as an
affine branch, not charged wholesale to near-minimality.

### 6.4 One-sided exponential maximum boxes

Let, and translate the product by its mean so that it is isotropic,

\[
 d\mu_n=e^{-\sum x_i}1_{\{x_i\ge0\}}dx,\qquad
 S=[0,q_n]^n,\qquad \mu_n(S)=1/2.                                  \tag{6.12}
\]

The Poincare constant is exactly four.  At small heat time the `n` top
faces are locally flat and carry equal flux, so

\[
 {R(s)\over\operatorname {tr}R(s)}\longrightarrow{I_n\over n}
 \qquad(s\downarrow0).                                           \tag{6.13}
\]

The codimension-two weighted edge content is

\[
 E_{2,n}=\binom n2e^{-2q_n}(1-e^{-q_n})^{n-2}
             \longrightarrow{(\log2)^2\over4}.                   \tag{6.14}
\]

Blow up an orthogonal edge by `y=y_0+sqrt(s)xi`.  The local posterior label
is

\[
 g_\angle(\xi)=\Phi(-\xi_1)\Phi(-\xi_2),\qquad
 z_\angle=\Phi^{-1}(g_\angle).                                  \tag{6.15}
\]

The normal Jacobian is `s`, while
`||nabla_y^2z||^2=s^{-2}||nabla_xi^2z_angle||^2` and
`e=|nabla_xi z_angle|^2`.  Thus both the Hessian square and the
eikonal square contribute at order `sqrt(s)`.  Each unit of weighted
edge content contributes

\[
 C_\angle\sqrt s+o(\sqrt s),\qquad
 C_\angle=\int_{\mathbb R^2}\varphi(z_\angle)
 \left\{\|\nabla^2z_\angle\|_{HS}^2+
 {1\over2}(1-|\nabla z_\angle|^2)^2\right\}d\xi
 \in(0,\infty),                                                  \tag{6.16}
\]

to `B'(s)`.  Gaussian tail bounds at infinity make the integral finite,
and nonaffinity makes it positive.  The local exponential density has
affine logarithm, so there is no leading curvature square.  Triple and
higher intersections are `O(s)`, while isolated flat faces have zero
leading Hessian and eikonal squares.  Hence

\[
 B'(s)=C_\angle E_{2,n}\sqrt s+O(s),\qquad
 B(s)-B(0)={2C_\angle E_{2,n}\over3}s^{3/2}+O(s^2).              \tag{6.17}
\]

This is a realizable high-rank phase pattern whose corner Bernstein charge
has exactly the `s^{3/2}` power allowed by (5.6).  Its bounded Poincare
constant means it is not a large-`K` counterexample.  It is nevertheless a
sharp local certificate: no propagation-plus-Bernstein lemma can demand a
larger power solely from the presence of many orthogonal faces.  Whether
these boxes are asymptotically Euclidean-Cheeger extremal is a separate,
unproved question and is not used here.

### 6.5 Cube

For the isotropic cube `[-sqrt(3),sqrt(3)]^n`, the label of a coordinate
half-cube and every heat posterior factor through one coordinate; `R(s)`
has rank one and the affine audit applies.
For a concentric half-volume inner cube, the small-time flux is equally
distributed among the coordinate faces, so `R/tr R` tends to `I_n/n`.
The codimension-two edge blow-up is again (6.15), and (6.17) holds with the
appropriate total Euclidean edge content.  In dimension two this is the
inner-square example: the normalized Bernstein dissipation tends to zero
although the limiting projector is `I_2/2`.  The set is far from a Cheeger
minimizer, so it certifies the necessity of global extremality but does not
contradict Proposition 1.

### 6.6 Simplex

Put the simplex in isotropic position.  For a facet-parallel cut, the
small-time label has one active normal and is asymptotically affine.  For a
concentric homothetic
half-volume simplex, permutation symmetry makes the `n+1` facet normals a
tight frame on the `n`-dimensional affine hull, hence

\[
                         R(s)/\operatorname {tr}R(s)\to I_n/n.   \tag{6.18}
\]

Pairwise facet intersections again have codimension two, and their
Gaussian blow-up gives a finite positive `sqrt(s)` term in `B'(s)`.  Thus
the simplex exhibits the same flat-phase/corner mechanism as the cube,
without a product decomposition.  No claim of Cheeger optimality for the
homothetic simplex is made.

## 7. Formal status

The following statements are proved in this report:

1. the complete Bernstein derivative (1.6)--(1.8);
2. the exact matrix Fisher equation (2.4) and scalar three-square identity
   (2.7);
3. the reverse-boundary transport estimate (3.9);
4. the marked rank-propagation Proposition 1, including the explicit
   fixed-scale consequence (4.15)--(4.17);
5. the exact near-minimality upper budget (5.4)--(5.6);
6. the propagation range (5.8) and the optimistic power obstruction
   (5.11)--(5.13).

The route does not prove KLS.  Its failure is stronger than a lack of
temporal control: rank does propagate for a nontrivial scale range, but
even perfect persistence supplies a Bernstein charge with the wrong power.
The Gaussian-ball and exponential-box calculations show that both relevant
powers are realized by honest log-concave models.  Reopening this family
requires an additional functional which distinguishes globally near-Cheeger
phase patterns from radial curvature and from codimension-two polyhedral
corner layers; another scalar use of `B` cannot do so.
