# Heat-flow Bernstein identity for Gaussian-profile slack

## 1. Setup and normalization

Let `mu` be a log-concave probability on `R^n`, possibly supported on an
affine subspace, and let `S` be Borel with

\[
        0<\mu(S)<1.
\]

For `s>0`, let

\[
 \gamma_s(x)=(2\pi s)^{-n/2}\exp(-|x|^2/(2s)),
 \qquad P_s\nu(y)=\int\gamma_s(y-x)\,d\nu(x),
\]

and define

\[
 q=q_s=P_s\mu,\qquad h=P_s(\mathbf 1_S\mu),\qquad
 g={h\over q},\qquad z=\Phi^{-1}(g).                 \tag{1.1}
\]

Thus `q` is a probability density with respect to ambient Lebesgue
measure `dy`, and every unlabelled spatial integral below is with respect
to `dy`.  Write

\[
 I(g)=\varphi(z),\qquad \rho=qI(g),\qquad
 \ell=\log q,qquad a=\nabla z,qquad K=\nabla^2z,
 \qquad w=|a|^2,qquad e=sw.                         \tag{1.2}
\]

Because the heat kernel is strictly positive and both `S` and its
complement have positive `mu`-mass, `q,h,q-h` are positive and smooth for
every `s>0`; hence `g` takes values in `(0,1)` and `z` is smooth.
Prékopa's theorem gives

\[
                 \nabla^2\ell\preceq0.              \tag{1.3}
\]

The posterior at `(s,y)` is

\[
 d\pi_{s,y}(x)={\gamma_s(y-x)\over q_s(y)}\,d\mu(x).
                                                               \tag{1.4}
\]

On the affine support of `mu`, it is `s^{-1}`-strongly log-concave.
The sharp centroid inequality for a strongly log-concave probability,
applied to `1_S`, gives

\[
                  \boxed{0\le e=s|\nabla z|^2\le1.} \tag{1.5}
\]

For completeness, differentiation of (1.4) gives, with
`m=E_pi X`, `A=Cov_pi(X)`, and
`v=E_pi[(1_S-g)(X-m)]`,

\[
 \nabla_y g={v\over s},\qquad
 \nabla_y\ell={m-y\over s},\qquad
 \nabla_y^2\ell={A-sI\over s^2}.                    \tag{1.6}
\]

Brascamp--Lieb on the affine support gives `A preceq sI`; in ambient
normal directions `A=0`.  Thus, in addition to (1.3),

\[
             -s^{-1}I\preceq\nabla^2\ell\preceq0.   \tag{1.7}
\]

## 2. Exact scalar PDEs

Both `q` and `h` solve the heat equation with generator `(1/2)Delta`.
Taking their quotient yields

\[
 \boxed{
 \partial_sg={1\over2}\Delta g+\langle\nabla\ell,\nabla g\rangle.}
                                                               \tag{2.1}
\]

Since `nabla g=I(g)nabla z` and

\[
 \Delta g=I(g)\{\Delta z-z|\nabla z|^2\},
\]

equation (2.1) is equivalent to

\[
 \boxed{
 \partial_sz={1\over2}\Delta z+\langle\nabla\ell,\nabla z\rangle
              -{z\over2}|\nabla z|^2.}              \tag{2.2}
\]

The boundary weight `rho=q varphi(z)` obeys a particularly useful
linear equation with a source:

\[
 \boxed{
       \partial_s\rho={1\over2}\Delta\rho+{1\over2}\rho w.}
                                                               \tag{2.3}
\]

Indeed, `nabla log varphi(z)=-z nabla z` and

\[
 \Delta\varphi(z)=\varphi(z)\{(z^2-1)w-z\Delta z\};
\]

substitution of (2.2) in `partial_s(q varphi(z))` proves (2.3), with
all mixed first-order terms cancelling.

Define

\[
 J(s)=\int\rho_s(y)\,dy,
 \qquad A_*(s)=\int\rho_s(y)(1-e_s(y))\,dy.           \tag{2.4}
\]

Here `A_*` denotes profile slack and is unrelated to the posterior
covariance `A` in (1.6).  Integrating (2.3) gives

\[
 J'(s)={1\over2}\int\rho w={1\over2s}\int\rho e,
 \qquad
 \boxed{
 {d\over ds}{J(s)\over\sqrt s}
       =-{A_*(s)\over2s^{3/2}}.}                    \tag{2.5}
\]

Thus `J(s)/sqrt(s)` is decreasing, and its loss is exactly the integrated
Gaussian-profile slack.

## 3. Pointwise Bochner identity

Set

\[
 \widetilde b=\nabla\ell-z\nabla z=\nabla\log\rho,
 \qquad
 \mathcal M=\partial_s-{1\over2}\Delta-
              \langle\widetilde b,\nabla\rangle.     \tag{3.1}
\]

Then

\[
 \boxed{
 \mathcal Mw=-w^2-\|K\|_{HS}^2
       +2\langle\nabla^2\ell\,a,a\rangle.}          \tag{3.2}
\]

To verify every sign, put
`L=(1/2)Delta+<nabla ell,nabla>`.  The ordinary Bochner formula is

\[
 Lw=\|K\|_{HS}^2+2\langle a,\nabla Lz\rangle
       -2\langle\nabla^2\ell\,a,a\rangle.           \tag{3.3}
\]

On the other hand, (2.2) gives

\[
 \partial_sw=2\langle a,\nabla Lz\rangle-w^2
                    -z\langle a,\nabla w\rangle.    \tag{3.4}
\]

Subtracting (3.3) from (3.4), then moving the final transport term into
the drift, gives (3.2).  Since `mathcal M s=1`, multiplication by `s`
gives

\[
 \boxed{
 \mathcal Me={e(1-e)\over s}-s\|K\|_{HS}^2
       +2s\langle\nabla^2\ell\,a,a\rangle.}         \tag{3.5}
\]

The same operator makes the vector gradient equation unusually simple:

\[
 \boxed{
 \mathcal M a=\nabla^2\ell\,a-{w\over2}a.}          \tag{3.6}
\]

Indeed, `mathcal Mz=zw/2`; commuting `nabla` with
`(1/2)Delta+<widetilde b,nabla>` and using
`nabla widetilde b=nabla^2ell-a otimes a-zK` cancels both `zKa`
terms.  Formula (3.6) separates rotation caused by the heat-density
curvature from rotation caused by spatial noise through `K`.

## 4. Exact integrated sum-of-squares identity

For a smooth integrable test function `F`, (2.3) and integration by
parts give

\[
 {d\over ds}\int\rho F
   =\int\rho\left(\partial_sF+{1\over2}\Delta F+{w\over2}F\right)
   =\int\rho\mathcal MF+{1\over2}\int\rho wF.       \tag{4.1}
\]

The second equality uses

\[
 \int\rho\langle\nabla\log\rho,\nabla F\rangle
       =-\int\rho\Delta F.                          \tag{4.2}
\]

Apply (4.1) to `F=e` and use (3.5).  If
`B(s)=int rho e`, then

\[
 B'(s)={1\over s}\int\rho(e-\tfrac12e^2)
       -s\int\rho\|K\|_{HS}^2
       +2s\int\rho\langle\nabla^2\ell a,a\rangle. \tag{4.3}
\]

Since `A_*=J-B` and `J'=B/(2s)`, subtraction gives

\[
 A_*'(s)=s\int\rho\|K\|_{HS}^2
       -2s\int\rho\langle\nabla^2\ell a,a\rangle
       -{1\over2s}\int\rho e(1-e).                 \tag{4.4}
\]

Adding `A_*/(2s)` to both sides produces the main identity.

**Theorem 1 (heat-flow Bernstein identity).**  Under the setup above,

\[
 \boxed{
 \begin{aligned}
 {d\over ds}\{\sqrt s\,A_*(s)\}
  ={}&s^{3/2}\int\rho\left(
       \|\nabla^2z\|_{HS}^2
       -2\langle\nabla^2\log q\,\nabla z,\nabla z\rangle\right)\\
    &+{1\over2\sqrt s}\int\rho(1-s|\nabla z|^2)^2.
 \end{aligned}}                                      \tag{4.5}
\]

Every term on the right is nonnegative by (1.3).  In particular,

\[
 \boxed{
 \int_{s_0}^{s_1}r^{3/2}\int\rho_r
            \|\nabla^2z_r\|_{HS}^2\,dy\,dr
 \le \sqrt{s_1}A_*(s_1)-\sqrt{s_0}A_*(s_0)}         \tag{4.6}
\]

for `0<s_0<s_1`.  There is no dimension factor, and the full Hessian,
including all mixed active--transverse and transverse--transverse
entries, occurs with coefficient one.

Combining (2.5) and the monotonicity in (4.5) gives a useful delayed
estimate.  Put `H(s)=J(s)/sqrt(s)`.  For every `s>0`,

\[
 \boxed{
 \begin{aligned}
 \int_s^{2s}r^{3/2}\int\rho_r\|\nabla^2z_r\|_{HS}^2\,dy\,dr
 &\le \sqrt{2s}A_*(2s)\\
 &\le 8s\{H(2s)-H(4s)\}.
 \end{aligned}}                                      \tag{4.7}
\]

Indeed `sqrt(r)A_*(r)` is nondecreasing, and hence

\[
 H(2s)-H(4s)
 ={1\over2}\int_{2s}^{4s}{A_*(r)\over r^{3/2}}dr
 ={1\over2}\int_{2s}^{4s}{\sqrt rA_*(r)\over r^2}dr
 \ge {\sqrt{2s}A_*(2s)\over8s}.
\]

Thus a profile loss on the *following* dyadic heat interval controls all
Hessian energy on the current interval.  This avoids any unjustified
pointwise conversion from an average profile deficit.

## 5. Posterior second moments and the nonsplitting term

At a fixed `(s,y)`, let `pi=pi_{s,y}` and use the posterior notation from
(1.6).  Define the correlated centered second moment

\[
 D=\mathbb E_\pi[(\mathbf1_S-g)(X-m)(X-m)^T].        \tag{5.1}
\]

Differentiating in the natural parameter `c=y/s` gives

\[
 \nabla_cg=v,\qquad \nabla_cv=D.                    \tag{5.2}
\]

Since `nabla_y=(1/s)nabla_c`, `nabla_y z=v/(sI(g))`, and
`I'(g)=-z`, one obtains the exact matrix relation

\[
 \boxed{
 \nabla_y^2z={D\over s^2I(g)}+z\,\nabla_yz\otimes\nabla_yz.}
                                                               \tag{5.3}
\]

If `u=v/|v|=nabla z/|nabla z|` and `P=I-u otimes u`, then

\[
 \boxed{
 P\nabla_y^2z={PD\over s^2I(g)},\qquad
 \nabla_yu={P\nabla_y^2z\over|\nabla_yz|}.}         \tag{5.4}
\]

Consequently the mixed conditional-law term that obstructs a
pointwise extension of split-posterior angular stability is not discarded:
it is exactly the projected Hessian in (4.5).  Heat-time integration
charges its Hilbert--Schmidt square to the derivative of scalar profile
slack.

On the set `{e>=kappa}`, (5.4) and `w=e/s` give

\[
 \|\nabla u\|_{HS}^2={\|P\nabla^2z\|_{HS}^2\over w}
       \le {s\over\kappa}\|\nabla^2z\|_{HS}^2.      \tag{5.5}
\]

Thus (4.6) implies the dimension-free angular-energy estimate

\[
 \boxed{
 \int_{s_0}^{s_1}\sqrt r\int_{\{e_r\ge\kappa\}}
       \rho_r\|\nabla u_r\|_{HS}^2\,dy\,dr
 \le {\sqrt{s_1}A_*(s_1)-\sqrt{s_0}A_*(s_0)\over\kappa}.}
                                                               \tag{5.6}
\]

There is also a curvature-rotation estimate.  From (1.7), for every unit
vector `u`,

\[
 \|P\nabla^2\ell\,u\|^2
 \le\|\nabla^2\ell\,u\|^2
 \le {1\over s}\{-\langle\nabla^2\ell u,u\rangle\}. \tag{5.7}
\]

On `{e>=kappa}`, the curvature term in (4.5) therefore yields

\[
 \boxed{
 \int_{s_0}^{s_1}r^{3/2}\int_{\{e_r\ge\kappa\}}
  \rho_r\|P_r\nabla^2\log q_r\,u_r\|^2\,dy\,dr
 \le {\sqrt{s_1}A_*(s_1)-\sqrt{s_0}A_*(s_0)\over2\kappa}.}
                                                               \tag{5.8}
\]

Together, (5.6) and (5.8) control the two mechanisms that rotate the
active direction in the gradient equation (3.6).

### 5.1 Reverse boundary diffusion

The weight `rho` has an exact stochastic interpretation.  Fix
`0<s_0<s_1`, set `s(tau)=s_1-tau`, and consider the time-inhomogeneous
diffusion

\[
 dY_\tau=\nabla\log\rho_{s(\tau)}(Y_\tau)d\tau+dB_\tau,          \tag{5.9}
\]

killed at instantaneous rate `w_{s(tau)}(Y_tau)/2`.  Its generator
before killing is

\[
 \widetilde L_s={1\over2}\Delta+\langle\nabla\log\rho_s,\nabla\rangle.
\]

The adjoint identity

\[
 \widetilde L_s^*\rho_s
 ={1\over2}\Delta\rho_s-
        \operatorname {div}(\rho_s\nabla\log\rho_s)
 =-{1\over2}\Delta\rho_s
\]

and (2.3) show that

\[
 -\partial_s\rho_s=\widetilde L_s^*\rho_s-{w_s\over2}\rho_s.  \tag{5.10}
\]

Hence, if the initial unnormalised density at `tau=0` is `rho_{s_1}`,
the alive density at reverse time `tau` is exactly `rho_{s(tau)}`.  This
is not merely an analogy: (5.10) is the forward equation of the killed
process.  Since `w<=1/s`, its survival mass on a dyadic interval obeys

\[
 {J(s_0)\over J(s_1)}\ge\sqrt{s_0\over s_1}.                    \tag{5.11}
\]

Thus for `s_1<=2s_0` it gives a common subcoupling of mass at least
`2^{-1/2}` between the normalized boundary laws at the two times.  The
earlier marginal of the survivor measure is exactly
`J(s_0)rho_{s_0}/J(s_0)`; the later marginal is survival-biased, so it
would be incorrect to call this a full coupling of both normalized laws.

### 5.2 A quantitative temporal gluing estimate

Direction is singular where `nabla z=0`.  A regularization avoids any
stopping argument.  Put

\[
 x_s=\sqrt s\,\nabla z_s,\qquad
 U_{\lambda,s}={x_s\over\sqrt{\lambda+|x_s|^2}},
 \qquad 0<\lambda\le1.                              \tag{5.12}
\]

Then `|x_s|^2=e_s<=1`, `|U_lambda|<=1`, and on `{e>=kappa}` the vector
`U_lambda` has the same direction as `u` and length at least
`sqrt(kappa/(lambda+kappa))`.  Equation (3.6) gives

\[
 \boxed{
 \mathcal Mx=\nabla^2\ell\,x+{1-e\over2s}x.}        \tag{5.13}
\]

Let `F_lambda(x)=x/(lambda+|x|^2)^{1/2}`.  Elementary differentiation
gives

\[
 \|DF_\lambda\|_{op}\le\lambda^{-1/2},\qquad
 \|D^2F_\lambda\|_{bil}\le {6\over\lambda},
 \qquad \|DF_\lambda(x)x\|\le1.                   \tag{5.14}
\]

The chain rule for `mathcal M` and `nabla x=sqrt(s)K` therefore gives

\[
 \mathcal MU_\lambda
 =DF_\lambda(x)\left(\nabla^2\ell x+{1-e\over2s}x\right)
  -{1\over2}D^2F_\lambda(x):\{sKK^T\}.             \tag{5.15}
\]

Let

\[
 \mathfrak D_{01}=\sqrt{s_1}A_*(s_1)-\sqrt{s_0}A_*(s_0)\ge0,
 \qquad s_1\le2s_0.                                \tag{5.16}
\]

Along the killed reverse diffusion, Itô's formula reads

\[
 dU_{\lambda,s(\tau)}(Y_\tau)
 =-\mathcal MU_{\lambda,s(\tau)}(Y_\tau)d\tau
   +\nabla U_{\lambda,s(\tau)}(Y_\tau)dB_\tau.     \tag{5.17}
\]

The exact occupation law (5.10), (4.5), and (5.14)--(5.15) imply the
following estimates, where expectations are *unnormalized* by integrating
the starting point against `rho_{s_1}dy` and paths are stopped at killing:

\[
 \mathbb E\langle M\rangle_{s_1-s_0}
 \le {\mathfrak D_{01}\over\lambda\sqrt{s_0}},      \tag{5.18}
\]

for the martingale in (5.17),

\[
 \mathbb E\left(\int |R_{\rm curv}|d\tau\right)^2
 \le {\mathfrak D_{01}\over2\lambda\sqrt{s_0}},
 \qquad
 \mathbb E\left(\int |R_{\rm scal}|d\tau\right)^2
 \le {\mathfrak D_{01}\over2\sqrt{s_0}},          \tag{5.19}
\]

and

\[
 \mathbb E\int |R_{\rm Hess}|d\tau
 \le {3\mathfrak D_{01}\over\lambda\sqrt{s_0}}.   \tag{5.20}
\]

Here the three drifts are respectively the three terms on the right of
(5.15).  To check the weights explicitly:

\[
 \begin{aligned}
 \|\nabla U_\lambda\|_{HS}^2
   &\le {s\over\lambda}\|K\|_{HS}^2,\\
 |R_{\rm curv}|^2
   &\le {e\over\lambda s}
       \{-\langle\nabla^2\ell u,u\rangle\},\\
 |R_{\rm scal}|^2&\le{(1-e)^2\over4s^2},\\
 |R_{\rm Hess}|&\le {3s\over\lambda}\|K\|_{HS}^2.
 \end{aligned}                                      \tag{5.21}
\]

For the curvature line, (1.7) was used.  Equations (5.18)--(5.20) then
follow by comparing (5.21) term by term with (4.5), using
`s>=s_0` and `s_1-s_0<=s_0`.

Conditioning the common subcoupling on survival and using (5.11) yields
the concrete temporal coherence bound

\[
 \boxed{
 \mathbb E_{\rm surv}|U_{\lambda,s_0}(Y_{s_0})
              -U_{\lambda,s_1}(Y_{s_1})|
 \le C\left\{
 \sqrt{{\mathfrak D_{01}\over
              \lambda\sqrt{s_0}J(s_0)}}
 +{\mathfrak D_{01}\over
              \lambda\sqrt{s_0}J(s_0)}\right\}.}  \tag{5.22}
\]

The constant `C` is numerical (the displayed estimates give, for
example, `C=12`).  Thus if the dimensionless relative endpoint defect

\[
       {\mathfrak D_{01}\over\sqrt{s_0}J(s_0)}
\]

is small, a constant fraction of the boundary laws at the two times have
regularized active directions coupled with small angular displacement.
Together with (4.7), the relative defect is controlled by the Gaussian
profile loss on the following dyadic interval.

More explicitly, for `s_0=s` and `s_1=2s`, (4.7) and
`J(s)=sqrt(s)H(s)` give

\[
 {\mathfrak D_{01}\over\sqrt sJ(s)}
 \le 8\,{H(2s)-H(4s)\over H(s)}.                   \tag{5.23}
\]

Thus (5.22) is small whenever the next-dyadic profile loss is a small
fraction of the current profile, with no scale or dimension loss.

This is genuine temporal gluing and uses no Poincaré inequality.  It does
not yet show that two *different* surviving trajectories at the same time
have the same direction; a spatial mixing or overlap argument is still
needed for that stronger conclusion.

### 5.3 Exact comparison with ordinary stochastic localization

The preceding calculation is the heat-coordinate version of an exact
ordinary-localization dissipation.  Put `t=1/s` and let

\[
 C_t=tX+\sqrt tG,qquad Y_s={C_t\over t}=X+\sqrt sG.              \tag{5.24}
\]

The posterior of `X` given either observation is (1.4), with natural tilt
`c=C_t`.  The law of `Y_s` is `q_sdy`.  At a posterior state define

\[
 r=|v|,quad u=v/r,quad P=I-u\otimes u,quad
 \eta={\sqrt t,r\over I(g)},quad
 \Delta={I(g)\over\sqrt t}-r,quad
 C=t^{-1}I-A.                                        \tag{5.25}
\]

Then the exact dictionary is

\[
 e=\eta^2,qquad
 P\nabla_y^2z={t^2PD\over I(g)},qquad
 \nabla_y^2\log q_s=-t^2C,                           \tag{5.26}
\]

and

\[
 \sqrt sA_*(s)=\mathbb E\{(1+\eta)\Delta\}.         \tag{5.27}
\]

For comparison, ordinary localization has the exact SDEs

\[
 dg=v^TdW,qquad dv=DdW-Avdt.                         \tag{5.28}
\]

Itô's formula therefore gives

\[
 dr=u^TDdW-r\,u^TAu\,dt+{\|PD\|_{HS}^2\over2r}dt
\]

and

\[
 d\left({I(g)\over\sqrt t}\right)
 =-{zv^T\over\sqrt t}dW
 -\left\{{r^2\over2\sqrt tI(g)}+{I(g)\over2t^{3/2}}\right\}dt.
\]

Subtracting and completing the scalar square yields

\[
 \boxed{
 d\Delta=dM-\left\{
 {\|PD\|_{HS}^2\over2r}
 +r\,u^TCu
 +{I(g)(1-\eta)^2\over2t^{3/2}}
 \right\}dt,}                                       \tag{5.29}
\]

where `dM=(-zv/sqrt(t)-Du)^TdW`.  Thus `Delta` is a nonnegative local
supermartingale.  Equations (5.26), the change `ds=-dt/t^2`, and
(5.27) show term by term that (5.29) and (4.5) have the same angular,
curvature, and scalar dissipations after averaging over the observation
law.  This independently audits all signs and powers of `s` and `t`.

## 6. Justification of differentiation and integration by parts

No smoothing of `mu` or of the set is actually needed before applying the
heat kernel.  To keep the inverse Gaussian CDF uniformly regular, first use

\[
 f_\epsilon=\epsilon+(1-2\epsilon)\mathbf1_S,
 \qquad 0<\epsilon<1/2.                              \tag{6.1}
\]

The quotient `g_epsilon=epsilon+(1-2epsilon)g` still solves (2.1), so
`z_epsilon` stays in the fixed compact interval
`[Phi^{-1}(epsilon),Phi^{-1}(1-epsilon)]`.  For an arbitrary finite
measure `mu`, Gaussian convolution makes `q` and `h` `C^infinity` and
strictly positive.  The formulas are therefore pointwise classical.

Here are global bounds sufficient to remove spatial cutoffs.  If `f` is
any `[0,1]`-valued posterior observable and

\[
 D_f=\mathbb E_\pi[(f-\mathbb E_\pi f)(X-m)(X-m)^T],
\]

then

\[
                         \|D_f\|_{HS}\le s.          \tag{6.2}
\]

Indeed, for a symmetric matrix `C` with `||C||HS=1`, Lichnerowicz for
the `s^{-1}`-strongly log-concave posterior gives

\[
 \operatorname {Var}_\pi((X-m)^TC(X-m))
 \le s\,\mathbb E_\pi|2C(X-m)|^2
 \le4s^2.
\]

Covariance Cauchy--Schwarz and `Var_pi(f)<=1/4`, followed by
Hilbert--Schmidt duality, prove (6.2).  Equations (5.3), (1.5), and the
compact range of `z_epsilon` consequently bound `nabla z_epsilon` and
`nabla^2z_epsilon` uniformly in space on every compact heat-time interval.
Higher derivatives are bounded there by the analogous Gaussian-kernel
formulas and centered posterior moments; strong log-concavity makes every
such moment finite with a bound depending only on its order, `n`, and the
time interval.  Dimension dependence in these temporary domination bounds
does not enter any inequality above.

Moreover, the score of a Gaussian convolution satisfies

\[
 \int q_s|\nabla\log q_s|^2dy\le {n\over s},          \tag{6.3}
\]

because it is the conditional expectation of the Gaussian noise divided
by `s`.  Hence `nabla rho_epsilon` is integrable: use
`nabla log rho=nabla log q-z nabla z`, the boundedness of `z_epsilon`,
and (6.3).  Multiplying all identities by standard cutoffs
`chi_R` with `|nabla chi_R|<=C/R`, `|Delta chi_R|<=C/R^2`, the boundary
terms are dominated by integrable functions times these vanishing
coefficients.  Equivalently, (4.2) is the distributional integration by
parts formula for `rho in W^{1,1}` and a bounded smooth `nabla F`.  This
proves (2.1)--(4.5) for fixed `epsilon` without any assumption on the
regularity or dimension of the support of `mu`.

Finally send `epsilon downarrow0`.  On every compact subset of
`(0,infinity) times R^n`, `z_epsilon` and all its derivatives converge to
`z`.  The terms on the right of (4.5) are nonnegative, so Fatou's lemma
gives (4.6), (4.7), (5.6), and (5.8).  The scalar endpoint quantities
converge by dominated convergence because
`0<=I(g_epsilon)<=1/sqrt(2pi)` and `0<=e_epsilon<=1`.  Whenever the
nonnegative spacetime terms are finite, the same cutoff argument gives
equality (4.5); in full generality it holds as an equality of nonnegative
Radon measures in `s`, while its integrated inequalities are sufficient
for all subsequent uses.

Lower-dimensional support causes no extra term: ambient Gaussian
convolution is full-dimensional, and (1.6)--(1.7) hold with zero posterior
covariance in normal directions.

## 7. Sign and scaling audits

### 7.1 Flat prior and a halfspace

Formally take Lebesgue measure as an improper flat prior and
`S={x_1>=0}`.  Then

\[
 q\equiv1,\qquad g(y)=\Phi(y_1/\sqrt s),\qquad
 z(y)=y_1/\sqrt s.
\]

Hence `e=1`, `nabla^2z=0`, `nabla^2log q=0`, and `A_*=0` pointwise.
Both sides of (4.5) vanish.  Although the transverse integrals are
infinite for this improper prior, this is an exact local audit of every
sign and factor.

### 7.2 Gaussian prior and a halfspace

Let `mu=N(0,sigma^2 I)` and `S={x_1>=0}`.  Direct Gaussian conditioning
gives

\[
 q_s=N(0,(\sigma^2+s)I),\quad
 z={\sigma y_1\over\sqrt{s(\sigma^2+s)}},\quad
 e={\sigma^2\over\sigma^2+s},\quad \nabla^2z=0,
 \quad\nabla^2\ell=-{I\over\sigma^2+s}.              \tag{7.1}
\]

Writing `c_0=(2pi)^(-1/2)`, one finds

\[
 J(s)=c_0\sqrt{s\over\sigma^2+s},\qquad
 A_*(s)=c_0{s^{3/2}\over(\sigma^2+s)^{3/2}}.         \tag{7.2}
\]

The derivative of the left side of (4.5) is

\[
 c_0{2\sigma^2s+s^2/2\over(\sigma^2+s)^{5/2}}.      \tag{7.3}
\]

The curvature term on the right is

\[
 c_0{2\sigma^2s\over(\sigma^2+s)^{5/2}},
\]

and the final scalar-square term is

\[
 c_0{s^2/2\over(\sigma^2+s)^{5/2}}.
\]

They add exactly to (7.3).  This audit detects both the factor `2` in the
curvature term and the factor `1/2` in the scalar-square term.

### 7.3 One-dimensional non-halfspace

For the same flat prior but the bounded set `S=[-a,a]`,

\[
 g_s(y)=\Phi((a-y)/\sqrt s)-\Phi((-a-y)/\sqrt s)     \tag{7.4}
\]

is nonlinear, `nabla^2z=z''` is nonzero, and the curvature term vanishes.
The identity becomes the nontrivial one-dimensional formula

\[
 {d\over ds}\left[\sqrt s\int\varphi(z)(1-sz'^2)dy\right]
 =s^{3/2}\int\varphi(z)z''^2dy
  +{1\over2\sqrt s}\int\varphi(z)(1-sz'^2)^2dy.    \tag{7.5}
\]

This provides a direct numerical audit with no posterior-curvature
contribution.  With `a=1`, adaptive quadrature on `[-8,8]` and a centered
relative step `10^{-4}` in `s` gives

\[
\begin{array}{c|c|c|c}
s&A_*(s)&\hbox{left side of (7.5)}&\hbox{right side of (7.5)}\\ \hline
0.5&0.2335725982&0.6924462790&0.6924462804\\
1&0.5618991786&0.8401834670&0.8401834677\\
2&0.9850855940&0.8049407664&0.8049407684
\end{array}
\]

The largest absolute discrepancy is `2.1 times 10^{-9}`.  The
derivatives used in the quadrature were analytic:

\[
 g'={\varphi(B)-\varphi(A)\over\sqrt s},\qquad
 g''={B\varphi(B)-A\varphi(A)\over s},\qquad
 z'={g'\over\varphi(z)},\qquad z''={g''\over\varphi(z)}+z(z')^2,
\]

where `A=(a-y)/sqrt(s)` and `B=(-a-y)/sqrt(s)`.

## 8. Model-class audit and present boundary

No symmetry or product assumption entered the proof.

* For the cube and a coordinate cut, `q_s` factorizes and `z` is
  one-dimensional; every transverse Hessian entry vanishes, consistently
  with (4.5).
* For a product of one-sided exponentials and a coordinate cut, the same
  one-dimensional reduction holds even though the convolved coordinate
  law is non-Gaussian.  The curvature and scalar-square terms pay for the
  strict centroid deficit.
* For a simplex, `q_s` is nonproduct and mixed Hessian entries generally
  occur.  Equation (4.5) retains their full Hilbert--Schmidt square; no
  conditional-independence or trace bound is inserted.
* For a radial exponential and a radial set, `z=z(r)` and
  `||nabla^2z||_HS^2=z''(r)^2+(n-1)(z'(r)/r)^2`.  Thus the angular
  multiplicity `n-1` is explicitly charged on the right of (4.5), rather
  than lost in a dimension-dependent estimate.

The identity closes the *conditional-law variation* gap in an integrated
heat-time sense.  It is not by itself a proof of KLS.  In particular,
small total profile loss controls angular energy on a preceding dyadic
interval via (4.7), but turning weighted Dirichlet energy of `u` into a
single coherent direction still requires a legitimate gluing mechanism.
The measure `rho=qI(g)` is not known to have a dimension-free Poincare
constant, and assuming one would reintroduce the conjecture.  Also, an
endpoint value `A_*(s_1)` cannot simply be declared small; (4.7) explains
the forward dyadic interval needed to control it.

## 9. A rigorous scale-`C_P` seed

One useful global input can nevertheless be proved without KLS.  Assume
`mu` has Poincaré constant `K<infinity`, and let `S` have mass `1/2`.
Then

\[
                         \boxed{J(K)\ge c_0}          \tag{9.1}
\]

for a numerical `c_0>0`; one may take

\[
 c_0=\min\left\{{\sqrt{2/\pi}\over8},
                    {1\over16I(1/2)}\right\}.
\]

Indeed the convolution `q_s=mu*gamma_s` satisfies

\[
 C_P(q_s)\le K+s.                                    \tag{9.2}
\]

To prove (9.2), decompose the variance of `F(X+G_s)` first in the
Gaussian variable and then in `X`; Jensen applied to the averaged
gradient gives the sum `s+K`.

Let `Y=X+G_K` and `g(Y)=P(S|Y)`.  Total variance gives

\[
 {1\over4}=\mathbb E_q[g(1-g)]+\operatorname {Var}_q(g).          \tag{9.3}
\]

If the first term is at least `1/8`, the one-dimensional Gaussian
Cheeger inequality

\[
 I(r)\ge\sqrt{2/\pi}\min(r,1-r)
      \ge\sqrt{2/\pi}\,r(1-r)
\]

gives the first lower bound in (9.1).  Otherwise `Var_q(g)>=1/8`,
and (9.2), (1.5), and `I(g)^2<=I(1/2)I(g)` give

\[
 {1\over8}\le\operatorname {Var}_q(g)
 \le2K\int q|\nabla g|^2
 =2K\int qI(g)^2w
 \le2I(1/2)J(K),                                     \tag{9.4}
\]

which proves the second lower bound.

Consequently

\[
                         H(K)={J(K)\over\sqrt K}ge{c_0\over\sqrt K}.
                                                               \tag{9.5}
\]

If a balanced set with exterior boundary at most `C/sqrt(K)` is chosen,
the small-time perimeter identity and monotonicity of `H` show that
`H(s)` remains within a universal multiplicative range throughout
`0<s<=K`.  Therefore, among any prescribed finite collection of dyadic
scales, at least one has a relative next-dyadic profile loss bounded by
the logarithmic total loss divided by the number of scales.  Equations
(4.7), (5.22), and (5.23) then give a genuinely dimension-free
near-equality scale.

There are two cautions.

1. Producing the balanced small-boundary set uses concavity of the
   isoperimetric profile (or an equivalent exact balanced reduction) and
   the small-time heat/perimeter limit; these must be supplied in any
   complete argument.
2. The tempting inequality `I(g)^2/[g(1-g)]>=c` is false: its left side
   is asymptotic to `g log(1/g)` as `g downarrow0`.  Thus covariance
   lower bounds cannot ignore tail states.  The seed (9.1) supplies a
   fixed amount of central mass at `s=K`, but a good near-equality scale
   found far below `K` need not inherit that centrality without an
   additional argument.

## 10. Quantified failure of the direct spatial-Poincaré closure

This section records a useful no-go calculation.  It prevents the exact
Bernstein identity from being turned into a claimed proof by hiding a bad
cutoff constant.

Suppose a dyadic pigeonhole argument finds a scale `r` with

\[
 r\ge\alpha K,qquad J(r)\ge c\sqrt\alpha,qquad
 {A_*(r)\over J(r)}\le C\varepsilon,qquad
 \int\rho_r\|\nabla^2z_r\|_{HS}^2
       \le {C\varepsilon J(r)\over r^2}.             \tag{10.1}
\]

These are exactly the bounds supplied by (4.7), (9.5), and a relative
profile drop `epsilon` on the next dyadic interval.  Let

\[
 U_\lambda={\sqrt r\nabla z\over
                  \sqrt{\lambda+r|\nabla z|^2}}.
\]

Choose a cutoff `chi(z)` which equals one on `|z|<=M`, vanishes outside
`|z|<=M+R`, and has `|chi'|<=C/R`.  On its support
`I(g)>=i_{M+R}:=varphi(M+R)`.  Hence

\[
 \int q_r|\nabla(\chi U_\lambda)|^2
 \le {C\varepsilon J(r)\over
             \lambda r i_{M+R}}+{C\over R^2r}.       \tag{10.2}
\]

Since `C_P(q_r)<=K+r`, ordinary Poincaré gives only

\[
 \operatorname {Var}_{q_r}(\chi U_\lambda)
 \le C(1+\alpha^{-1})
       \left\{{\varepsilon\over\lambda i_{M+R}}
                         +{1\over R^2}\right\}.      \tag{10.3}
\]

To make the cutoff-gradient term a small numerical constant requires
`R^2>=c/alpha`.  But then

\[
 i_{M+R}\le C\exp(-c/\alpha),                         \tag{10.4}
\]

so the Hessian term requires

\[
                   \varepsilon\le C\exp(-c/\alpha). \tag{10.5}
\]

The profile pigeonhole over the range `[alpha K,K]` supplies at best

\[
                   \varepsilon\le {C\over\log(1/\alpha)}.       \tag{10.6}
\]

The bounds (10.5) and (10.6) are incompatible in the regime in which
one hopes to improve a large `K`.  In particular, taking
`alpha=K^{a-1}` gives `R^2>=K^{1-a}` and an exponentially small
`i_{M+R}`, whereas the available profile loss is only `O(1/log K)`.
Thus this cutoff-plus-`C_P(q)` scheme does not prove an inequality
`K<=CK^a` for any `a<1`.

The obstruction is structural, not just algebraic.  A set with several
well-separated flat faces can have boundary weight concentrated on phases
with distinct normals; each phase has almost zero local angular energy and
the transition occurs where `I(g)` is tiny.  The weight `rho=qI(g)` sees
the faces but not a Poincaré bridge between them.  Isoperimetric
near-minimality may rule out such a multiphase configuration, but that
requires a separate stability or splitting theorem.  Neither (4.5) nor
`C_P(q)<=K+s` supplies it automatically.

Here is a concrete asymptotic counterexample to an unqualified coherence
lemma.  Let `mu` be uniform on `[-1,1]^2` and let

\[
                         S=[-2^{-1/2},2^{-1/2}]^2.
\]

Then `mu(S)=1/2`.  As `s downarrow0`, the normalized boundary law
`nu_s=rho_sdy/J(s)` gives asymptotically equal weight to the four faces.
Consequently

\[
 \mathbb E_{\nu_s}[u_su_s^T]\longrightarrow {1\over2}I_2,       \tag{10.7}
\]

so even the unoriented normal line has no coherent top eigendirection.
Nevertheless the normalized dissipation in (4.5) tends to zero.  On every
face portion a fixed positive distance from a corner and from the outer
square, the rescaled problem converges exponentially fast to the flat
halfspace, for which `e=1` and `nabla^2z=0`.  A corner neighborhood has
area `O(s)` after the change `y=y_0+sqrt(s)xi`, while
`nabla_y^2z=s^{-1}nabla_xi^2z_corner`.  Hence each corner contributes

\[
 s^2\int_{\rm corner}\rho_s\|\nabla^2z_s\|^2dy=O(s).
\]

Since `J(s)~sqrt(s)mu^+(S)`, the normalized Hessian contribution is
`O(sqrt(s))`; the scalar-square term has the same order, and the curvature
term from `q_s` is exponentially small because the inner square is a fixed
distance from the support boundary.  This proves that log-concavity and
small Bernstein dissipation alone cannot yield projector coherence.
The inner square is not an isoperimetric minimizer, pinpointing the extra
hypothesis a successful clustering lemma must exploit.

## 11. Formal exact-attainment calculation (superseded by Section 12)

This preliminary calculation assumes that a balanced Cheeger minimizer is
attained.  It is retained only as a short guide to the identities.  The
audited statement in Section 12 does **not** assume attainment: it constructs
a half-mass `varepsilon`-minimizer and retains the error in every inequality.

Let `f=1_S`, let `A_sf(y)=E[f(X)|Y=y]=g_s(y)` for
`Y=X+G_s`, and let

\[
 T_s=A_s^*A_s,qquad
 T_sf(x)=\mathbb E[g_s(Y)\mid X=x].                  \tag{11.1}
\]

The operator `T_s` is a self-adjoint Markov contraction on `L^2(mu)`.  Put

\[
 U(s)=\mathbb E_{q_s}[g_s(1-g_s)].                   \tag{11.2}
\]

For balanced binary `f`, self-adjointness gives the exact identity

\[
 \boxed{\mathbb E_\mu|f-T_sf|=2U(s).}                \tag{11.3}
\]

Indeed `<f,T_sf>=||A_sf||_2^2=E_qg^2`, while both `f` and `T_sf` have
mean `1/2`.

Assume now that `mu^+(S)=p` and that the Cheeger constant is `psi=2p`.
The `L^1` Cheeger inequality and (11.3) imply

\[
 \begin{aligned}
 \int|\nabla T_sf|d\mu
 &\ge\psi\inf_c\int|T_sf-c|d\mu\\
 &\ge2p\left\{{1\over2}-\mathbb E|T_sf-f|\right\}
 =p\{1-4U(s)\}.                                      \tag{11.4}
 \end{aligned}
\]

On the other hand differentiation of the Gaussian convolution gives

\[
 \nabla T_sf(x)=\mathbb E[\nabla g_s(Y)\mid X=x],
\]

and hence

\[
 \boxed{
 p\{1-4U(s)\}
 \le\int|\nabla T_sf|d\mu
 \le\mathbb E_q|\nabla g_s|
 \le H(s)\le p.}                                    \tag{11.5}
\]

The last inequality is the monotonicity (2.5) together with the
small-time perimeter limit.  Since

\[
 I(r)\ge\sqrt{2/\pi}\,r(1-r),
\]

one also has

\[
 U(s)\le {J(s)\over\sqrt{2/\pi}}
       ={\sqrt sH(s)\over\sqrt{2/\pi}}
       \le {\sqrt s\,p\over\sqrt{2/\pi}}.           \tag{11.6}
\]

Thus if `p<=C/sqrt(K)` and `s=alpha K`, the entire chain (11.5) loses
only `O(sqrt(alpha))` in relative terms.

The triangle inequality in (11.5) also has an exact angular stability
interpretation.  Let `W=grad g_s(Y)`, and, when the conditional mean is
nonzero, set

\[
 \theta(X)={\mathbb E[W\mid X]\over|\mathbb E[W\mid X]|}.
\]

Then

\[
 \begin{aligned}
 &\mathbb E|W|-\mathbb E|\mathbb E[W\mid X]|\\
 &\qquad=\mathbb E\{|W|(1-\langle W/|W|,\theta(X)\rangle)\}
 \le4pU(s),
 \end{aligned}                                       \tag{11.7}
\]

and therefore

\[
 \boxed{
 \mathbb E\{|W|\,|W/|W|-\theta(X)|^2\}\le8pU(s).} \tag{11.8}
\]

This aligns the posterior active direction with a direction field on the
original sample space.  It does not say that `theta(X)` is constant.

For reference, the Gibbs kernel `T_s` has a quantitative spectral gap
which follows only from `C_P(mu)=K`:

\[
 \boxed{
 \langle h,(I-T_s)h\rangle
 \ge {s\over K+2s}\operatorname {Var}_\mu(h).}       \tag{11.9}
\]

To prove it, set `D=E Var(h(X)|Y)` and `g_h=E[h|Y]`.  Posterior
`s^{-1}`-strong log-concavity gives

\[
 |\nabla g_h|^2={|\operatorname {Cov}(h,X\mid Y)|^2\over s^2}
 \le {\operatorname {Var}(h\mid Y)\over s}.
\]

Together with `C_P(q_s)<=K+s`, this yields
`Var_q(g_h)<=(K+s)D/s`; total variance then proves (11.9).
For `h=f`, the left side is exactly `U(s)`.  The gap (11.9) explains why
(11.5) alone recovers the usual Cheeger--Poincaré equivalence but does not
yet improve its exponent: iterating long enough to mix also accumulates
the one-step uncertainty.

## 12. Audited global-minimality and posterior-resampling lemma

This section gives the version used in the checkpoint.  In contrast with
Section 11, it neither assumes that the Cheeger infimum is attained nor drops
the near-minimizer error.

### 12.1 Obtaining a balanced near-minimizer

All perimeter statements are relative to the affine hull `E` of `mu`.  Let
`P_mu` denote relaxed weighted `BV` perimeter and set

\[
 \mathcal I_\mu(v)=\inf\{P_\mu(B):\mu(B)=v\},\qquad 0<v<1.   \tag{12.1}
\]

We use the concavity theorem for the isoperimetric profile under
`CD(0,infinity)`: on a complete smooth weighted Riemannian manifold with
`Ric+Hess V>=0`, the lower-semicontinuous profile is concave.  For a
log-concave probability `e^{-V}dx` on a convex subset of Euclidean space,
the same statement follows by smooth positive convex approximation,
exhaustion of the convex support, and lower semicontinuity of weighted
perimeter.  This is the nonsmooth Euclidean specialization of the
Bavard--Pansu/Bayle--Morgan--Johnson concavity theorem; see E. Milman,
*Invent. Math.* 177 (2009), Theorem 1.8 and its stated approximation
convention.  It applies on `E`, including when `E` is lower-dimensional.

The profile is symmetric because complementation preserves perimeter.
Concavity and `mathcal I_mu(0)=0` therefore give

\[
 {\mathcal I_\mu(v)\over\min(v,1-v)}
 \ge 2\mathcal I_\mu(1/2),\qquad 0<v<1.              \tag{12.2}
\]

The exterior-Minkowski and relaxed-perimeter definitions of the Cheeger
constant agree by Lipschitz relaxation, coarea, and distance cutoffs.  Hence

\[
                  \psi_\mu=2\mathcal I_\mu(1/2).    \tag{12.3}
\]

It follows that, for every `varepsilon>0`, there is a finite-perimeter set
`S=S_varepsilon` satisfying

\[
 \mu(S)=1/2,\qquad p:=P_\mu(S)le {\psi_\mu\over2}+\varepsilon. 
                                                               \tag{12.4}
\]

No minimizer is asserted to exist.  Weighted `BV` blow-up at the reduced
boundary gives the heat-profile perimeter formula

\[
 \lim_{r\downarrow0}{1\over\sqrt r}
 \int q_r I\left({P_r(\mathbf1_S\mu)\over q_r}\right)=P_\mu(S)=p. 
                                                               \tag{12.5}
\]

The normalization in (12.5) is exact: the tangent label is a Euclidean
halfspace and its normal integral is `int_R varphi(t)dt=1`.  On a
lower-dimensional affine support the ambient-normal Gaussian factor
integrates to one, reducing the formula to `E`.  Equivalently, prove (12.5)
first for smooth cuts and pass through a strict weighted-`BV` approximation.
This is also the approximation that justifies using the set in all formulas
below.  If profile concavity is not invoked, the rest of Section 12 remains
valid for any half-mass finite-perimeter `S`; only (12.4) must then be
supplied separately.

### 12.2 Conditional-expectation operators and the exact `L^1` identity

Let `f=mathbf1_S` and consider the joint law

\[
 X\sim\mu,\qquad Y=X+G_s,\qquad G_s\sim N(0,sI).
\]

Define

\[
 A_sh(y)=\mathbb E[h(X)\mid Y=y],\qquad
 A_s^*k(x)=\mathbb E[k(Y)\mid X=x].                 \tag{12.6}
\]

Then `A_sf=g_s` and

\[
 T_s=A_s^*A_s,\qquad T_sf(x)=\mathbb E[g_s(Y)\mid X=x].      \tag{12.7}
\]

Conditional Jensen shows that `A_s` and `A_s^*` are `L^2` contractions;
they preserve positivity and constants.  Thus `T_s` is a self-adjoint
positive Markov contraction.  Set

\[
                 U(s)=\mathbb E_{q_s}[g_s(1-g_s)].   \tag{12.8}
\]

The following identity has no inequality hidden in it.  Since
`0<=T_sf<=1`, both `f` and `T_sf` have mean `1/2`, and

\[
 \langle f,T_sf\rangle_{L^2(\mu)}
 =\langle A_sf,A_sf\rangle_{L^2(q_s)}=\mathbb E_{q_s}g_s^2,
\]

splitting the absolute value on `S` and `S^c` gives

\[
 \begin{aligned}
 \mathbb E_\mu|f-T_sf|
 &=\mathbb E\{f(1-T_sf)+(1-f)T_sf\}\\
 &=1-2\langle f,T_sf\rangle
 =2\mathbb E_{q_s}(g_s-g_s^2)
 =\boxed{2U(s)}.                                    \tag{12.9}
 \end{aligned}
\]

The precise normalization of the `L^1` Cheeger inequality is

\[
 \int_E|\nabla_EF|d\mu
 \ge\psi_\mu\inf_{c\in\mathbb R}\int_E|F-c|d\mu.  \tag{12.10}
\]

Layer-cake and coarea prove it first for locally Lipschitz `F`, and `BV`
relaxation gives the displayed generality.  Since
`inf_c int|f-c|dmu=1/2`, the triangle inequality and (12.9) yield

\[
 \begin{aligned}
 \int|\nabla T_sf|d\mu
 &\ge\psi_\mu\inf_c\int|T_sf-c|d\mu\\
 &\ge\psi_\mu\left\{{1\over2}-\int|T_sf-f|d\mu\right\}\\
 &=\psi_\mu\{1/2-2U(s)\}.                          \tag{12.11}
 \end{aligned}
\]

There is also no regularity gap in differentiating `T_sf`.  Indeed

\[
 T_sf(x)=\int g_s(x+z)\gamma_s(z)dz,
\]

and the posterior centroid bound gives the global estimate

\[
 |\nabla g_s|=I(g_s)|\nabla z_s|
       \le {I(1/2)\over\sqrt s}.                    \tag{12.12}
\]

Dominated differentiation therefore gives

\[
 \nabla_ET_sf(x)=\mathbb E[\nabla_Eg_s(Y)\mid X=x].           \tag{12.13}
\]

Write `W=grad_E g_s(Y)`.  Conditional Jensen, (1.5), (2.5), and the
small-time limit (12.5) give

\[
 \begin{aligned}
 \int|\nabla T_sf|d\mu
 &\le\mathbb E|W|
 =\int q_sI(g_s)|\nabla z_s|\\
 &\le {J(s)\over\sqrt s}=H(s)\le p.                \tag{12.14}
 \end{aligned}
\]

Combining the two sides gives the audited global contraction chain

\[
 \boxed{
 \psi_\mu\{1/2-2U(s)\}
 \le\int|\nabla T_sf|d\mu
 \le\mathbb E|W|\le H(s)\le p.}                    \tag{12.15}
\]

For an attained minimizer, `varepsilon=0` and `psi_mu=2p`, but neither fact
is used in (12.15).  Finally, concavity of the Gaussian profile, symmetry,
and `I(0)=0` imply `I(r)>=2I(1/2)min(r,1-r)`.  Since
`min(r,1-r)>=r(1-r)`, this gives the elementary one-variable bound

\[
 I(r)\ge\sqrt{2/\pi}\,r(1-r),\qquad 0\le r\le1,    \tag{12.16}
\]

implies

\[
 \boxed{
 U(s)\le {J(s)\over\sqrt{2/\pi}}
 ={\sqrt sH(s)\over\sqrt{2/\pi}}
 \le {\sqrt s\,p\over\sqrt{2/\pi}}.}              \tag{12.17}
\]

### 12.3 Jensen deficit and a bounded regularized direction feature

Define the norm-Jensen deficit

\[
 \delta_J(s)=\mathbb E|W|-\mathbb E|\mathbb E[W\mid X]|.    \tag{12.18}
\]

Equations (12.15) and (12.4), with no attainment assumption, give

\[
 \begin{aligned}
 0\le\delta_J(s)
 &\le p-\psi_\mu\{1/2-2U(s)\}\\
 &=p-\psi_\mu/2+2\psi_\mu U(s)\\
 &\le\boxed{\varepsilon+4pU(s)}.                   \tag{12.19}
 \end{aligned}
\]

For `m(X)=E[W|X]`, choose a measurable unit vector
`theta(X)=m(X)/|m(X)|` on `{m!=0}` and any fixed unit vector on `{m=0}`.
Put `u(Y)=W/|W|` on `{W!=0}`.  On the null set `{W=0}` choose `u`
arbitrarily.  Conditional expectation gives the exact identity

\[
 \begin{aligned}
 \delta_J(s)
 &=\mathbb E\{|W|(1-\langle u(Y),\theta(X)\rangle)\}\\
 &={1\over2}\mathbb E\{|W||u(Y)-\theta(X)|^2\}.     \tag{12.20}
 \end{aligned}
\]

This remains true on `{m=0}`, because the conditional expectation of the
inner-product term is zero there.  Thus

\[
 \boxed{
 \mathbb E\{|W||u(Y)-\theta(X)|^2\}
 \le2\varepsilon+8pU(s).}                          \tag{12.21}
\]

The estimate aligns the posterior active direction with a direction field
on the original sample space; it does not make that field constant.

For a bounded-feature formulation, (12.12) supplies

\[
                    |W|\le M_s:={I(1/2)\over\sqrt s}.         \tag{12.22}
\]

For `lambda>0`, define

\[
 F_\lambda(Y)={W\over\sqrt{\lambda+|W|^2}}
 =\kappa_\lambda(|W|)u(Y),\qquad
 \kappa_\lambda(r)={r\over\sqrt{\lambda+r^2}}.     \tag{12.23}
\]

This regularized direction feature has norm at most one.  Given `X`, let
`Y,Y'` be conditionally independent channel outputs and use primes for the
corresponding quantities.  The polarization identity gives the exact split

\[
 \begin{aligned}
 &{1\over2}\mathbb E[|F_\lambda(Y)-F_\lambda(Y')|^2\mid X]\\
 &\quad={1\over2}\mathbb E[(\kappa_\lambda(|W|)
                  -\kappa_\lambda(|W'|))^2\mid X]\\
 &\qquad\quad+{1\over2}\mathbb E[\kappa_\lambda(|W|)
             \kappa_\lambda(|W'|)|u-u'|^2\mid X].  \tag{12.24}
 \end{aligned}
\]

Call the expectation of the last line `D_phase(lambda)`.  If
`a(X)=E[|W||X]` and `m(X)=E[W|X]`, conditional independence gives

\[
 \begin{aligned}
 {1\over2}\mathbb E[|W||W'|\,|u-u'|^2\mid X]
 &=a(X)^2-|m(X)|^2\\
 &\le2M_s\{a(X)-|m(X)|\}.                           \tag{12.25}
 \end{aligned}
\]

Since `kappa_lambda(r)<=r/sqrt(lambda)`, averaging (12.25) and using
(12.19) yields

\[
 \boxed{
 D_{\rm phase}(\lambda)
 \le {2M_s\over\lambda}\delta_J(s)
 \le {2M_s\over\lambda}\{\varepsilon+4pU(s)\}.}   \tag{12.26}
\]

At the natural regularization `lambda=M_s^2`, this becomes

\[
 D_{\rm phase}(M_s^2)
 \le {2\sqrt s\over I(1/2)}\{\varepsilon+4pU(s)\}. \tag{12.27}
\]

The first term in (12.24) is an amplitude Dirichlet energy and is not
controlled by norm-Jensen stability: arbitrary positive multiples of one
fixed vector are equality cases for the triangle inequality.  This is an
exact obstruction to replacing `D_phase` by the full conditional variance
of `F_lambda` in (12.26).

### 12.4 Full proof of the posterior-resampling spectral gap

Assume `C_P(mu)=K<infinity`.  Then, for every `h in L^2(mu)`,

\[
 \boxed{
 \langle h,(I-T_s)h\rangle_{L^2(\mu)}
 \ge {s\over K+2s}\operatorname {Var}_\mu(h).}      \tag{12.28}
\]

It is enough first to take bounded `h`, and then use `L^2` truncation.  Set

\[
 g_h(y)=\mathbb E[h(X)\mid Y=y],\qquad
 D=\mathbb E\operatorname {Var}(h(X)\mid Y).
\]

Adjointness and conditional variance give

\[
 D=\mathbb E h^2-\mathbb E g_h^2
  =\langle h,(I-T_s)h\rangle.                       \tag{12.29}
\]

Differentiating the posterior formula yields

\[
 \nabla_Eg_h(y)={1\over s}\operatorname {Cov}(h(X),X\mid Y=y).
                                                               \tag{12.30}
\]

The posterior covariance satisfies `A_{s,y}preceq sI_E`.  Conditional
Cauchy--Schwarz, tested against every unit vector, therefore gives

\[
 |\operatorname {Cov}(h,X\mid Y)|^2
 \le s\operatorname {Var}(h\mid Y),\qquad
 \int|\nabla g_h|^2dq_s\le {D\over s}.              \tag{12.31}
\]

The convolution `q_s=mu*gamma_s` has Poincare constant at most `K+s`.
Indeed, apply the product Poincare inequality to `k(x+z)`, first in `x`
with constant `K` and then in `z` with Gaussian constant `s`; its two
gradients are both `grad k`.  Consequently

\[
 \operatorname {Var}_{q_s}(g_h)
 \le(K+s)\int|\nabla g_h|^2dq_s
 \le {K+s\over s}D.                                 \tag{12.32}
\]

The total-variance identity

\[
 \operatorname {Var}_\mu(h)=D+\operatorname {Var}_{q_s}(g_h)
\]

and (12.32) prove (12.28).  Equations (12.30)--(12.32) pass from bounded
truncations to every `L^2` function by closedness of the weak gradient and
lower semicontinuity.  For `h=f`, (12.29) is exactly `D=U(s)`.

The dual Gibbs kernel `widetilde T_s=A_sA_s^*` on `L^2(q_s)` has the same
gap.  More explicitly, (12.28) says that the squared operator norm of
`A_s` from the mean-zero part of `L^2(mu)` to the mean-zero part of
`L^2(q_s)` is at most `1-delta_s`, where

\[
                     \delta_s={s\over K+2s}.         \tag{12.33}
\]

The adjoint has the same norm.  Hence, componentwise for scalar or
Euclidean-vector-valued `F in L^2(q_s)`,

\[
 \boxed{
 \mathbb E\operatorname {Var}(F(Y)\mid X)
 =\langle F,(I-\widetilde T_s)F\rangle
 \ge\delta_s\operatorname {Var}_{q_s}(F).}          \tag{12.34}
\]

For `F=F_lambda`, the left side is exactly the sum of the amplitude and
phase terms in (12.24).

### 12.5 Sharp scaling obstruction

Suppose `p<=C_0/sqrt(K)`, put `s=alpha K` with `0<alpha<=1`, and let the
near-minimizer error tend to zero.  Equations (12.17) and (12.19) give

\[
 U(s)\le C_1\sqrt\alpha,\qquad
 {\delta_J(s)\over p}\le4U(s)\le C_2\sqrt\alpha.    \tag{12.35}
\]

The posterior-resampling gap, however, is only

\[
              \delta_s={\alpha\over1+2\alpha}\asymp\alpha.   \tag{12.36}
\]

The bounded feature makes the power mismatch exact.  Since
`I(1/2)=1/sqrt(2pi)`, (12.17) and (12.27) imply

\[
 D_{\rm phase}(M_s^2)
 \le {8\sqrt s\,pU(s)\over I(1/2)}
 \le C_3sp^2
 \le C_4{s\over K}.                                 \tag{12.37}
\]

Thus the available phase-Dirichlet upper bound is of exactly the same order
as the Gibbs gap, not smaller.  Before regularization this is the equivalent
mismatch `U(s)~sqrt(s/K)` versus `delta_s~s/K`.  Ordinary Poincare
contraction therefore cannot force phase variance to vanish.  Moreover,
(12.24) contains an uncontrolled amplitude energy.  A successful use of
(12.21) needs a phase-specific coercivity or a multiscale cancellation that
improves one of these powers; merely iterating `T_s` cannot do so.
