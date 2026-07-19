# Fixed-scale Fisher rank, angular wedge energy, and the missing extremality charge

## 0. Verdict

There is a completely quantitative fixed-scale statement which does not
require a dyadic pigeonhole.  Let `mu` be isotropic and log-concave, put
`K=C_P(mu)`, choose a half-mass near-Cheeger set `S`, and smooth its label at

\[
                         s=\alpha K.                 \tag{0.1}
\]

For every sufficiently small fixed numerical `alpha`, one obtains all of the
following at this *prescribed* scale:

1. the posterior-resampling uncertainty is bounded above and below with the
   correct powers of `alpha`;
2. the Gaussian-profile centroid defect, normalized by the profile flux, is
   `O(sqrt(alpha))`;
3. a set of central posteriors of `q_s`-mass `Omega(sqrt(alpha))` has
   pointwise centroid defect `O(sqrt(alpha))`;
4. the binary Fisher matrix carried by those posteriors has trace at least

   \[
                         {\sqrt\alpha\over8\pi},      \tag{0.2}
   \]

   and effective rank at least

   \[
                         {\alpha^{3/2}K\over8\pi}.    \tag{0.3}
   \]

All constants below are explicit.  In particular, the fixed-scale seed is
not the missing step.

The attempted wedge closure fails for a precise reason.  General angular
stability controls the good-state angular derivative by
`C_delta Omega_delta(tau)/s`, whereas the rank-sensitive Poincare inequality
forces total derivative energy only of order `sqrt(alpha)/K`.  Since
`s=alpha K` and the available pointwise defect is
`tau=O(sqrt(alpha))`, the ratio of the available upper scale to the required
lower scale is

\[
                 {C_\delta\Omega_\delta(O(\sqrt\alpha))
                    \over \alpha^{3/2}}.             \tag{0.4}
\]

It becomes worse, not better, as `alpha` decreases.  Moreover the wedge
energy may be paid entirely by longitudinal amplitude changes or by the
small bad-state set.  The coarea deficit is a scalar integral over nested
levels and contains no cross-normal term capable of excluding this.

Section 6 states a formal spatial phase-charge lemma which would be sufficient
for dimension-free KLS and proves that implication.  The lemma is not proved
by the four audited inputs.  The product of one-sided exponentials with its
balanced maximum box is the exact unresolved stress test: all perimeter,
normal-projector, and limiting profile calculations are explicit, but its
Euclidean Cheeger optimality is not available.  Thus it cannot honestly be
declared either an admissible counterexample or excluded.  An unqualified
weighted-projector lemma (one which omits global near-minimality) is already
refuted by the inner-square example in `heatflow_bernstein.md`; the genuinely
extremality-dependent version remains the load-bearing new theorem.

## 1. Setup and constants

Let `mu` be a full-dimensional isotropic log-concave probability on
`R^n`.  The same argument works on its affine hull, but isotropy already
forces full-dimensional support.  Put

\[
                         K=C_P(\mu).                  \tag{1.1}
\]

Since linear functions have variance one, `K>=1`.  Let `psi=psi_mu`.  Fix a
half-mass finite-perimeter set `S` such that

\[
 P_\mu(S)=p\le {\psi\over2}+\varepsilon_0,
 \qquad \beta={\varepsilon_0\over p}.                \tag{1.2}
\]

Such sets exist by concavity and symmetry of the log-concave isoperimetric
profile and strict `BV` approximation.  We will assume

\[
             0<\alpha\le {1\over4},\qquad
             0\le\beta\le\sqrt\alpha.               \tag{1.3}
\]

There is no existence issue in imposing the second condition: after `alpha`
is fixed, choose the near-minimizer error smaller than
`sqrt(alpha) psi/2`.

Let

\[
 X\sim\mu,\qquad Y=X+\sqrt sG,\qquad s=\alpha K,
 \qquad q_s={\cal L}(Y),                              \tag{1.4}
\]

and write

\[
 g(y)=P(X\in S\mid Y=y),\quad z=\Phi^{-1}(g),
 \quad I(g)=\varphi(z),\quad \rho=q_s I(g).          \tag{1.5}
\]

Set

\[
 e=s|\nabla z|^2,\quad
 J(s)=\int\rho,\quad H(s)={J(s)\over\sqrt s},
 \quad A_*(s)=\int\rho(1-e).                        \tag{1.6}
\]

The sharp posterior centroid inequality gives `0<=e<=1`.  We shall use the
two numerical Gaussian constants

\[
 I_0=I(1/2)={1\over\sqrt{2\pi}},\qquad
 c_G=\sqrt{2/\pi}=2I_0,                              \tag{1.7}
\]

and the explicit scale-`K` seed from `heatflow_bernstein.md`,

\[
 J(K)\ge c_0,\qquad
 c_0=\min\left\{{c_G\over8},{1\over16I_0}\right\}
     ={I_0\over4}.                                   \tag{1.8}
\]

The equality on the right follows directly from `c_G=2I_0` and
`I_0^2=1/(2pi)`.

For self-containment, here is the seed proof.  At `s=K`, total variance of
the balanced binary label gives

\[
 {1\over4}=E[g(1-g)]+\operatorname {Var}_{q_K}(g).   \tag{1.9}
\]

If the first term is at least `1/8`, (3.12) below gives
`J(K)>=c_G/8`.  Otherwise `Var(g)>=1/8`.  Since
`C_P(q_K)<=2K`, `e=K|\nabla z|^2<=1`, and `I(g)^2<=I_0I(g)`,

\[
 {1\over8}\le\operatorname {Var}(g)
 \le2K\int q_K I(g)^2|\nabla z|^2
 =2\int q_KI(g)^2e
 \le2I_0J(K).                                       \tag{1.10}
\]

This proves (1.8), including its exact constants.

## 2. Fixed-scale uncertainty and profile defect

### 2.1 Two-sided control of `U`

Let `f=1_S`, let `A_s f=E[f(X)|Y]`, and put `T_s=A_s^*A_s`.  Define

\[
 U(r)=E_{q_r}[g_r(1-g_r)]
     =\langle f,(I-T_r)f\rangle.                    \tag{2.1}
\]

The posterior-resampling gap proved in `heatflow_bernstein.md` gives

\[
 U(r)\ge {r\over K+2r}\operatorname {Var}(f)
        ={r\over4(K+2r)}.                            \tag{2.2}
\]

The heat-profile upper bound gives

\[
                         U(r)\le {\sqrt r\,p\over c_G}.       \tag{2.3}
\]

Cheeger's inequality in the normalization of the task is
`K<=4/psi^2`; hence `psi<=2/sqrt(K)`.  From (1.2),

\[
 p(1-\beta)\le {\psi\over2}\le {1\over\sqrt K}.    \tag{2.4}
\]

Consequently

\[
 \boxed{
 {\alpha\over4(1+2\alpha)}\le U(s)
 \le {\sqrt\alpha\over c_G(1-\beta)},
 \qquad
 U(2s)\le u_2:={\sqrt{2\alpha}\over c_G(1-\beta)}.} \tag{2.5}
\]

No Buser--Ledoux estimate is used in (2.5).

### 2.2 An endpoint estimate for the centroid defect

The exact profile identity is

\[
 H'(r)=-{A_*(r)\over2r^{3/2}},                       \tag{2.6}
\]

and the Bernstein sum-of-squares identity says that

\[
                         r\longmapsto\sqrt r A_*(r) \quad
                         \hbox{is nondecreasing}.    \tag{2.7}
\]

Therefore, for `r in [s,2s]`,
`A_*(r)>=sqrt(s/r)A_*(s)`, and hence

\[
\begin{aligned}
 H(s)-H(2s)
 &=\int_s^{2s}{A_*(r)\over2r^{3/2}}dr\\
 &\ge {A_*(s)\sqrt s\over2}\int_s^{2s}{dr\over r^2}
 ={A_*(s)\over4\sqrt s}.                            \tag{2.8}
\end{aligned}
\]

The audited global-minimality chain gives, at every `r>0`,

\[
                         H(r)\ge\psi(1/2-2U(r)),
 \qquad H(r)\le p.                                  \tag{2.9}
\]

It follows that

\[
 p-H(2s)
 \le p-\psi/2+2\psi U(2s)
 \le\varepsilon_0+4pU(2s).                          \tag{2.10}
\]

For `s<=K`, monotonicity of `H` and (1.8) give

\[
 H(s)\ge H(K)={J(K)\over\sqrt K}
              \ge {c_0\over\sqrt K}.               \tag{2.11}
\]

Combining (2.4), (2.8)--(2.11) proves the explicit normalized defect bound

\[
 \boxed{
 \eta_s:={A_*(s)\over J(s)}
 \le {4\{\beta+4u_2\}\over c_0(1-\beta)}.}         \tag{2.12}
\]

Notice that this is a bound at the left endpoint `s`, not merely at an
unspecified scale.  The forward interval `[s,2s]` and monotonicity (2.7) are
both essential.

For later reference define

\[
 C_*={8\over c_0}\left(1+{8\sqrt2\over c_G}\right),
 \qquad
 \alpha_*=min\left\{{1\over4},{1\over(8C_*)^2}\right\}.      \tag{2.13}
\]

When `alpha<=alpha_*` and `beta<=sqrt(alpha)`, equations (2.5) and (2.12)
give

\[
                         \boxed{\eta_s\le C_*\sqrt\alpha
                                      \le {1\over8}.}          \tag{2.14}
\]

Indeed, `1/(1-beta)<=2`, so `u_2<=2sqrt(2alpha)/c_G`, and substitution in
(2.12) gives exactly the first inequality in (2.14).  Numerically,
`c_0=I_0/4`, `C_*<1220`, and one may safely take `alpha=10^{-10}`; the
symbolic definition (2.13) is the audited constant.

Finally, (2.11) supplies the profile mass

\[
                         \boxed{J(s)\ge c_0\sqrt\alpha.}       \tag{2.15}
\]

## 3. Central good posteriors and an explicit Fisher-rank seed

At a posterior state put

\[
 m=E[X|Y],\qquad
 v=\operatorname {Cov}(f(X),X|Y),\qquad
 u={v\over|v|}\quad(v\ne0).                         \tag{3.1}
\]

Then

\[
 \nabla g={v\over s},\qquad
 e={|v|^2\over sI(g)^2}.                             \tag{3.2}
\]

Choose

\[
 M_\alpha=\sqrt{2\log{8I_0\over c_0\sqrt\alpha}},
 \qquad \delta_\alpha=\Phi(-M_\alpha).             \tag{3.3}
\]

Thus

\[
 I(\delta_\alpha)=\varphi(M_\alpha)
                  ={c_0\sqrt\alpha\over8}.         \tag{3.4}
\]

Put `tau=8eta_s` and define the central good set

\[
 {\cal G}=\{y:\delta_\alpha\le g(y)\le1-\delta_\alpha,
                         \ 1-e(y)\le\tau\}.         \tag{3.5}
\]

When `alpha<=alpha_*`, (2.14) gives `tau<=1`.  Since `I` is symmetric and
increasing on `[0,1/2]`,

\[
 \int_{\{g\notin[\delta_\alpha,1-\delta_\alpha]\}}
                   q_s I(g)e
 \le I(\delta_\alpha)\le {J(s)\over8}.              \tag{3.6}
\]

If `eta_s>0`, Markov's inequality in the *profile-flux measure* gives

\[
 \int_{\{1-e>\tau\}}q_sI(g)e
 \le\int_{\{1-e>\tau\}}q_sI(g)
 \le {A_*(s)\over\tau}={J(s)\over8}.                \tag{3.7}
\]

If `eta_s=0`, then `A_*(s)=0`; since its integrand is nonnegative,
`1-e=0` for `rho`-almost every state, and (3.7) means the same zero bound
without division by `tau=0`.

As `int q_sI(g)e=J(s)-A_*(s)`, equations (2.14), (3.6), and (3.7) give

\[
 \boxed{
                    \int_{\cal G}q_sI(g)e\ge {J(s)\over2}.}   \tag{3.8}
\]

In particular, because `eI(g)<=I_0`,

\[
 \boxed{q_s({\cal G})\ge {J(s)\over2I_0}
                       \ge {\sqrt\alpha\over8}.}    \tag{3.9}
\]

The last equality uses `c_0=I_0/4`.

Define the binary Fisher feature and its good-state second-moment matrix by

\[
 h=2\arcsin\sqrt g,\qquad
 F=\sqrt s\,\nabla h={v\over\sqrt{s\,g(1-g)}},     \tag{3.10}
\]

\[
 R_{\cal G}=E_{q_s}[FF^T1_{\cal G}],\qquad
 R=E_{q_s}[FF^T].                                   \tag{3.11}
\]

The elementary Gaussian profile inequality

\[
                         I(a)\ge c_Ga(1-a)           \tag{3.12}
\]

implies

\[
\begin{aligned}
 \operatorname {tr}R_{\cal G}
 &=\int_{\cal G}q_s\,{eI(g)^2\over g(1-g)}\\
 &\ge c_G\int_{\cal G}q_s eI(g)
 \ge {c_GJ(s)\over2}
 \ge {\sqrt\alpha\over8\pi}.                     \tag{3.13}
\end{aligned}
\]

The last constant is exact:
`c_Gc_0/2=(2I_0)(I_0/4)/2=I_0^2/4=1/(8pi)`.

Conditional covariance Cauchy--Schwarz gives

\[
 {vv^T\over g(1-g)}\preceq\operatorname {Cov}(X|Y).
\]

Averaging and using total covariance and isotropy yields

\[
                         R_{\cal G}\preceq R\preceq {1\over s}I.       \tag{3.14}
\]

Consequently

\[
 \boxed{
 {\operatorname {tr}R_{\cal G}\over\|R_{\cal G}\|_{op}}
 \ge s\operatorname {tr}R_{\cal G}
 \ge {\alpha^{3/2}K\over8\pi}.}                    \tag{3.15}
\]

Equations (3.8)--(3.15) are the promised central Fisher seed.  They use no
conditioning on `cal G`; all matrices are unnormalized integrals under the
original log-concave observation law.  Thus no non-log-concave conditional
measure is introduced.

Finally, on `cal G` the pointwise deficit used by
`general_angular_stability.md` is

\[
 \epsilon_y=1-{|v|\over\sqrt s\,I(g)}=1-\sqrt e
 \le1-e\le\tau.                                     \tag{3.16}
\]

At the non-strict endpoint where `tau=1`, states with `v=0` are omitted
before invoking angular stability.  The fixed choice (3.17) has
`tau_0<0.098`, so every state in `cal G_0` has `e>0` and hence `v ne0`.

### 3.1 One fixed universal instantiation

To remove any ambiguity about auxiliary choices, take once and for all

\[
 \boxed{\alpha_0=10^{-10},\qquad \beta_0=10^{-5},
 \qquad \tau_0=8C_*10^{-5}<0.098,}                  \tag{3.17}
\]

and

\[
 \boxed{
 M_0=\sqrt{2\log(3.2\cdot10^6)},\qquad
 \delta_0=\Phi(-M_0).}                              \tag{3.18}
\]

Here `C_*<1220`, and (3.18) is (3.3) with
`sqrt(alpha_0)=10^{-5}` because `8I_0/c_0=32`.  Given a balanced
near-minimizer, choose its error so that `beta<=beta_0`, put
`s=alpha_0K`, and define

\[
 {\cal G}_0=\{\delta_0\le g\le1-\delta_0,
                         \ 1-e\le\tau_0\}.          \tag{3.19}
\]

Since `eta_s<=C_*10^{-5}` and `tau_0=8C_*10^{-5}`, the proof of
(3.7)--(3.15) applies verbatim (with a possibly larger good set) and gives

\[
 \boxed{
 q_s({\cal G}_0)\ge1.25\cdot10^{-6},\qquad
 \operatorname {tr}R_{{\cal G}_0}\ge {10^{-5}\over8\pi},
 \qquad
 \operatorname {rank}_{eff}(R_{{\cal G}_0})
       \ge {10^{-15}K\over8\pi}.}                  \tag{3.20}
\]

Thus `alpha`, `tau`, and `delta` can all be frozen before the dimension,
measure, or near-minimizing sequence is chosen.

## 4. Exact angular derivative and the wedge power audit

Define the posterior second label moment

\[
 D=E[(f-g)(X-m)(X-m)^T\mid Y].                      \tag{4.1}
\]

Differentiation with respect to the natural parameter `c=y/s` gives

\[
                         \nabla_c v=D.               \tag{4.2}
\]

Indeed, differentiate `v=E[fX]-gm`; the two product-rule terms combine to
the centered expression (4.1).  Therefore, in the heat coordinate,

\[
                         \nabla_yu={P_uD\over s|v|},
 \qquad P_u=I-uu^T.                                  \tag{4.3}
\]

Since `F=|v|u/sqrt(sg(1-g))`, projecting its output derivative removes the
amplitude derivative and gives the exact identity

\[
 \boxed{
 P_u\nabla F={P_uD\over s^{3/2}\sqrt{g(1-g)}}.}      \tag{4.4}
\]

The posterior is `s^{-1}`-strongly log-concave.  Applying the general
angular stability theorem with `t=1/s`, (3.16), and
`g in [delta_alpha,1-delta_alpha]` yields on `cal G`

\[
 \|P_uD\|_{HS}^2
 \le C_{\delta_\alpha}s^2
                    \Omega_{\delta_\alpha}(\tau).  \tag{4.5}
\]

Equations (4.4)--(4.5) imply the pointwise and integrated bounds

\[
 \boxed{
 \int_{\cal G}\|P_u\nabla F\|_{HS}^2dq_s
 \le {C_{\delta_\alpha}\over s}
                    \Omega_{\delta_\alpha}(\tau).} \tag{4.6}
\]

Here and below the factor `1/[delta_alpha(1-delta_alpha)]` is absorbed into
`C_delta`; no dimension enters.

On the other hand, the rank-sensitive wedge Poincare inequality applies to
the full feature `F`.  Since `C_P(q_s)<=K+s`, (3.13)--(3.14) imply, whenever

\[
                         s\operatorname {tr}R>4,     \tag{4.7}
\]

that

\[
 \boxed{
 E_{q_s}\|\nabla F\|_{HS}^2
 \ge {\operatorname {tr}R\over2(K+s)}
 \ge {\sqrt\alpha\over16\pi(1+\alpha)K}.}          \tag{4.8}
\]

Condition (4.7) holds once
`K>32pi/alpha^(3/2)`, by (3.13).  If it fails, `K` is already bounded by a
universal constant because `alpha` has been fixed.

Comparing (4.6) and (4.8), with `s=alpha K`, gives the exact obstruction:

\[
 {\hbox{available good angular upper scale}
  \over\hbox{forced total lower scale}}
 \le C_{\delta_\alpha}(1+\alpha)
       {\Omega_{\delta_\alpha}(\tau)\over\alpha^{3/2}}.        \tag{4.9}
\]

The direction of this estimate is not itself a contradiction: (4.6) is an
upper bound on only part of the energy, whereas (4.8) is a lower bound on
the whole energy.  Formula (4.9) merely shows that even if one could remove
that logical mismatch, the powers do not close.  The audited modulus is,
up to a logarithm, `Omega_delta(r)=O_delta(r^(1/6))`, while
`tau=O(sqrt(alpha))`.  Thus the right side of (4.9) behaves no better than
`alpha^(-17/12)` as `alpha` decreases.  Even a hypothetical linear angular
modulus would leave a factor `alpha^{-1}`.

The orthogonal decomposition

\[
 \|\nabla F\|_{HS}^2
 =\|P_u\nabla F\|_{HS}^2+\|u^T\nabla F\|^2         \tag{4.10}
\]

identifies the first missing term: longitudinal amplitude energy.  The
second missing term is the unrestricted derivative energy on
`cal G^c`.  Neither is controlled by the norm-Jensen direction identity.
That identity weights the joint `(X,Y)` law by `|\nabla g|` and permits the
direction field attached to `X` to vary from one phase cell to another.

## 5. What the coarea deficit does and does not control

Let

\[
 F_s=T_s1_S,
 \qquad E_r=\{F_s>r\}.                               \tag{5.1}
\]

The audited coarea cascade gives

\[
 \mathcal D_{co}(F_s)
 :=\int_0^1\left[P_\mu(E_r)-
       \psi\min\{\mu(E_r),1-\mu(E_r)\}\right]dr
 \le\varepsilon_0+2\psi U(s).                      \tag{5.2}
\]

Since `psi<=2p`, (2.5) gives

\[
 \boxed{
 {\mathcal D_{co}(F_s)\over p}
 \le\beta+4U(s)
 \le\beta+{4\sqrt\alpha\over c_G(1-\beta)}.}      \tag{5.3}
\]

Thus the relative coarea deficit is `O(sqrt(alpha))`.  It is important that
(5.2) is an exact integral of nonnegative scalar deficits.  Boolean
combinations of the nested levels produce only disjoint level bands, whose
perimeters add.  They produce no term containing

\[
                         1-\langle u(y),u(y')\rangle \tag{5.4}
\]

for two spatially separated phase packets.  Consequently neither (5.2) nor
Markov extraction from it controls the two missing terms in (4.10).

This is not a removable regularity issue.  The canonical reduced boundaries
of regular levels fix the representative of each level set, but their
scalar perimeter is still counted only once even when many phase pairs reuse
the same transition skeleton.

### 5.1 An exact near-minimizer refutes longitudinal absorption

There is a particularly clean reason that one cannot repair (4.10) by
charging all longitudinal energy to `D_co`.  Let

\[
                         d\lambda(x)={1\over2}e^{-|x|}dx,
 \qquad S=[0,\infty).                                \tag{5.5}
\]

The symmetric exponential law is log-concave, `lambda(S)=1/2`, and its
exact one-dimensional profile is

\[
                         \mathcal I_\lambda(a)=\min(a,1-a).   \tag{5.6}
\]

Thus `S` is an exact balanced Cheeger minimizer.  For every heat time, the
posterior probability `g_s(y)=P(X>=0|X+sqrt(s)G=y)` is strictly increasing.
This follows directly from the monotone likelihood ratio of the Gaussian
kernel: for `y_2>y_1`, the ratio
`gamma_s(y_2-x)/gamma_s(y_1-x)` is strictly increasing in `x`.  Its further
Gaussian average `T_s1_S(x)` is also strictly increasing.  Hence every
nontrivial superlevel of `T_s1_S` is a half-line, and (5.6) gives

\[
                         \boxed{\mathcal D_{co}(T_s1_S)=0.}   \tag{5.7}
\]

If isotropic normalization is desired, replace `X` by `X/sqrt(2)` and
rescale `s`; zero coarea deficit and positivity of the energy below are
unchanged.

On the other hand, the one-dimensional Fisher feature

\[
                         F(y)=\sqrt s\,{d\over dy}
                          \left(2\arcsin\sqrt{g_s(y)}\right)  \tag{5.8}
\]

has strictly positive derivative energy.  If `F'` vanished almost
everywhere under `q_s`, positivity and smoothness of `q_s` would make `F`
constant on `R`.  Then `2arcsin sqrt(g_s)` would be affine.  It is bounded
between zero and `pi` and nonconstant, an impossibility.  Therefore

\[
                         \int |F'(y)|^2q_s(y)dy>0.             \tag{5.9}
\]

Equations (5.7)--(5.9) give an explicit log-concave *exact near-minimizer*
counterexample to every proposed inequality of the form

\[
 \int\|u^T\nabla F\|^2dq_s
 \le {C\over s}\,{\mathcal D_{co}(T_s1_S)\over P_\lambda(S)}. \tag{5.10}
\]

This is precisely the affine branch.  It does not refute a theorem which
correctly removes affine energy before charging incompatible phases, but it
proves that such a branch is logically mandatory and that longitudinal
energy cannot simply be absorbed into the scalar coarea deficit.

## 6. A formal phase-charge lemma sufficient for KLS

**Status boundary.**  Sections 1--5, including the fixed-scale theorem,
angular derivative identity, power audit, and symmetric-exponential
counterexample, are proved.  The lemma in this section is a proposed new
theorem and is **not proved** here.  Only the implication
"proposed lemma implies KLS" is proved below.

The preceding calculations isolate a self-contained theorem which would
close the route.  It is stated here so that its strength and its exceptional
branches cannot be hidden.

**Proposed spatial phase-charge lemma.**  There are numerical constants
`c_sp>0`, `N_sp<infinity`, `K_sp<infinity`, `C_rad<infinity`, and a
nondecreasing modulus `Xi(r)->0` as `r->0` with the
following property.  Let `nu` be an isotropic log-concave probability with
`C_P(nu)=K>=K_sp`, let `S` be a half-mass finite-perimeter set with
`P_nu(S)<=psi_nu/2+epsilon_0`, and form the heat objects at `s=alpha K`.
Here `0<alpha<=alpha_*`.  Suppose that for some `delta,tau>0`, some
`N>=1`, and a Borel state set `G`,

\[
\begin{gathered}
 G\subset\{\delta\le g\le1-\delta,\ 1-e\le\tau\},\\
 \operatorname {tr}R_G\ge {\sqrt\alpha\over8\pi},\qquad
 {\operatorname {tr}R_G\over\|R_G\|_{op}}\ge N.
                                                               \tag{6.1}
\end{gathered}
\]

Then at least one of the following holds:

1. **affine branch:** there is a line `L` such that
   `tr(P_{L^perp}R_G)<=Xi(tau)trR_G`;
2. **radial branch:** the flux-normal field is, outside relative flux
   `Xi(tau)`, within squared angle `Xi(tau)` of the rays from one center,
   and the resulting radial separator satisfies `K<=C_rad`;
3. **phase charge:**

   \[
   {\mathcal D_{co}(F_s)\over P_\nu(S)}
   \ge c_{sp}-\Xi(\tau)-{N_{sp}\over N}.             \tag{6.2}
   \]

One may and will replace `c_sp` by `min(c_sp,1)`, so assume
`0<c_sp<=1`.

The radial conclusion deliberately includes its own dimension-free
Poincare bound; merely naming a configuration radial would not close KLS.
The affine conclusion is stated in the same Fisher matrix which appears in
(6.1), so it is quantitatively incompatible with high effective rank.

Here is the complete implication from this lemma to KLS.  Fix `alpha`
smaller than `alpha_*` and then smaller, if necessary, so that

\[
 C_*\sqrt\alpha\le1/8,\qquad
 \Xi(8C_*\sqrt\alpha)\le c_{sp}/4,\qquad
 \sqrt\alpha+{4\sqrt\alpha\over c_G(1-\sqrt\alpha)}
                       \le c_{sp}/4.                 \tag{6.3}
\]

Choose the near-minimizer so that `beta<=sqrt(alpha)`.  Sections 2--3
produce (6.1) with

\[
 \delta=\delta_\alpha,\qquad \tau\le8C_*\sqrt\alpha,
 \qquad N={\alpha^{3/2}K\over8\pi}.                \tag{6.4}
\]

If `K` is larger than

\[
 \max\left\{K_{sp},{16\pi\over\alpha^{3/2}},
 {32\pi N_{sp}\over c_{sp}\alpha^{3/2}}\right\},  \tag{6.5}
\]

then `N>2` and the last term in (6.2) is at most `c_sp/4`.  The affine branch is impossible:
if `tr(P_{L^perp}R_G)<=Xi trR_G` with `Xi<c_sp/4<1/2`, then
`rank_eff(R_G)<=1/(1-Xi)<2`, contradicting (6.4) for the same large `K`.
The radial branch already gives `K<=C_rad`.  The phase-charge branch gives
`D_co/p>=c_sp/2`, contradicting (5.3) and (6.3).  Thus `K` is bounded by the
maximum of the constants in (6.5) and `C_rad`.  This is dimension-free KLS.

The implication is rigorous, but the proposed lemma is not a consequence of
the audited inputs.  In particular, replacing (6.2) by a bound on a
tilt-space separator or an incidence-weighted midpoint multiplicity is not
enough; one needs an actual charge to the physical coarea deficit.

## 7. Countermodels and the product-exponential survivor

### 7.1 Why global near-minimality is indispensable

For the uniform measure on `[-1,1]^2`, let

\[
                         S=[-2^{-1/2},2^{-1/2}]^2.   \tag{7.1}
\]

It has mass one half.  As heat time tends to zero, the boundary-flux normal
projector tends to `I_2/2`, while the normalized Bernstein dissipation tends
to zero.  The four flat faces are distinct phase cells and all direction
changes occur in corner layers of vanishing normalized mass.  Thus any
weighted-projector inequality inferred only from small local centroid or
Bernstein defect is false.  This is a fully explicit log-concave example.
It is not a near-Cheeger set: a coordinate half-square has much smaller
perimeter.  Therefore it refutes the unqualified lemma, not the genuinely
extremality-dependent statement (6.2).

Gaussian maximum cuts give the analogous high-rank winner geometry, but
their perimeter is `Theta(sqrt(log n))` while a halfspace has constant
perimeter.  Their coarea deficit is therefore a fixed fraction of their
perimeter and (6.2) is compatible with them.

### 7.2 Exact one-sided exponential calculations

Let

\[
 d\mu_n(x)=e^{-\sum_{i=1}^nx_i}1_{\{x_i\ge0\}}dx
\]

and let

\[
 B_q=[0,q]^n,\qquad m=\mu_n(B_q)=(1-e^{-q})^n.      \tag{7.2}
\]

Translation by the mean makes `mu_n` isotropic.  Its exact Poincare constant
is four in every dimension, by the one-dimensional exponential constant and
tensorization.  Thus the large-`K` restriction in the proposed lemma
deliberately prevents this bounded-`K` family from being a literal
counterexample to (6.2).  It remains a mandatory test for every *local*
phase-charge or splicing estimate used to prove that lemma.

The relative Euclidean boundary consists of `n` disjoint coordinate faces,
so its perimeter is exactly

\[
\begin{aligned}
 P_{\mu_n}(B_q)
 &=ne^{-q}(1-e^{-q})^{n-1}\\
 &=n(1-m^{1/n})m^{(n-1)/n}.                         \tag{7.3}
\end{aligned}
\]

For every fixed `m in (0,1)`,

\[
                         P_{\mu_n}(B_q)\longrightarrow-m\log m.       \tag{7.4}
\]

At `m=1/2`,

\[
                         p_n\longrightarrow{\log2\over2}.            \tag{7.5}
\]

The normalized boundary-flux projector is exactly

\[
                         {1\over n}\sum_{i=1}^ne_ie_i^T={1\over n}I.
                                                               \tag{7.6}
\]

Thus this is the canonical high-rank, nonconcurrent flat-phase pattern.
The coordinate median competitor has perimeter `1/2`, so it does not beat
(7.5).  Unlike the Gaussian and uniform-cube maximum cuts, no elementary
competitor in the present argument produces a fixed relative gap below
(7.5).

For a nested box with mass `m=1/2+x`, expansion of (7.4) gives

\[
 -m\log m={\log2\over2}+(\log2-1)x+O(x^2).          \tag{7.7}
\]

If the Euclidean Cheeger constant were asymptotic to `log 2`, the scalar
Cheeger deficit of these levels would be only linear in `|x|`.  Therefore a
heat family whose level masses deviate from one half by total amount
`O(U)` would spend only `O(pU)`, exactly the budget in (5.2), while retaining
the projector (7.6).  This explains why the example survives every scalar
coarea test in this route.

There is an exact theorem that boxes are asymptotically extremal for
one-sided exponential products under *uniform/supremum enlargement*, at
least in the monotone class (S. G. Bobkov, *Studia Math.* **123** (1997),
81--95, for product measures with the supremum metric and monotone Borel
sets).  It does not imply the Euclidean isoperimetric
statement required here: because `d_infinity<=d_2`, uniform neighborhoods
are larger and their boundary measure is an upper, not a lower, control on
Euclidean boundary.  Importing that theorem would therefore be a metric
error.

Accordingly, (7.2)--(7.7) are an explicit stress test but not a proved
near-minimizer counterexample to (6.2).  Proving either

\[
 \psi_{\mu_n}\le(1-c)\log2
 \quad\hbox{for a universal }c>0,                   \tag{7.8}
\]

or

\[
                         \psi_{\mu_n}\longrightarrow\log2           \tag{7.9}
\]

would decide whether this family is excluded by or refutes the proposed
phase charge.  Neither conclusion follows from the standard
Bobkov--Houdre tensorization estimate, which is only up to universal
constants.

## 8. Final status

The fixed-scale part is complete and formalizable:

\[
\begin{gathered}
 U(s)\in\left[{\alpha\over4(1+2\alpha)},
       {\sqrt\alpha\over c_G(1-\beta)}\right],\\
 {A_*(s)\over J(s)}\le C_*\sqrt\alpha,\qquad
 q_s({\cal G})\ge {\sqrt\alpha\over8},\\
 \operatorname {tr}R_{\cal G}\ge {\sqrt\alpha\over8\pi},\qquad
 \operatorname {rank}_{eff}(R_{\cal G})
       \ge {\alpha^{3/2}K\over8\pi}.               \tag{8.1}
\end{gathered}
\]

Every estimate is on the original observation law and keeps `alpha`, the
near-minimizer error, and the central cutoff explicit.
The fully frozen choice is (3.17)--(3.20):
`alpha=10^{-10}`, `tau<0.098`, and
`delta=Phi(-sqrt(2log(3.2*10^6)))`.

What is not proved is the spatial phase-charge lemma (6.2).  The exact
angular derivative (4.4), the wedge lower bound (4.8), and the coarea budget
(5.3) show that no rearrangement of the existing scalar estimates can prove
it: the powers mismatch and the longitudinal/bad-state energies remain
uncontrolled.  The next viable step must construct a physical competitor
which splices distinct canonical boundary patches and proves that the
resulting perimeter saving cannot be reused by many phase pairs, while
classifying the equality cases as affine or genuinely radial.  Any proposed
version must be run first on (7.2)--(7.7).
