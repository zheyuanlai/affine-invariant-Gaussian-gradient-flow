# Projection-wise radial exclusion and the endpoint trichotomy

## Executive conclusion

There are two dimension-free statements which survive variable endpoint
lengths and do not use geometric product conductances.

First, let `F` be a boundary-weighted family of balanced log-concave rays,
let

\[
 P_F=\int_F q_y(0)\,d\eta(y),\qquad
 d\nu(y)={q_y(0)\over P_F}\,d\eta(y),                    \tag{0.1}
\]

and suppose `sigma_y>=s` on `F`.  For every orthogonal projection `P` of
rank `d`,

\[
 \boxed{\quad
 P_Fs^3 E_\nu|PN|^2\le C d,\qquad
 P_Fs^5 E_\nu|PN|^4\le C d.
 \quad}                                                   \tag{0.2}
\]

The first inequality is total covariance; the second is thin shell for the
isotropic log-concave marginal `PX`.  In particular, if `r_N` is the
dimension of the linear span of the normals, then

\[
 r_N\ge cP_Fs^5.                                        \tag{0.3}
\]

If the family has actual ray mass `alpha` and lies in one scale band
`s<=sigma_y<=A_\sigma s`, then

\[
 P_F\asymp_{\delta,A_\sigma}{\alpha\over s},
 \qquad r_N\ge {c_\delta\alpha\over A_\sigma}s^4.        \tag{0.4}
\]

Thus padding ambient directions which do not occur in the normals does
nothing: the projection onto the normal span removes them.  A genuine
high-scale escape needs at least order `s^4` active normal dimensions, not
merely ambient dimension.

Second, suppose the ideal two-endpoint heat-bath form on `(F,nu)` has gap
at least `1-epsilon`.  Write the finite endpoint half-lengths as

\[
 B=Z+r(B)N,\qquad C=Z-\ell(C)N,\qquad
 L=r(B)+\ell(C)=|B-C|.                                  \tag{0.5}
\]

The fact that `r` is a function of `B` and `ell` a function of `C` follows
from the global signed-distance function.  For independent endpoint
marginals define the nonnegative Lipschitz slack

\[
 G(b,c)=|b-c|^2-(r(b)+\ell(c))^2\ge0.                   \tag{0.6}
\]

Then

\[
 \boxed{\quad
 E_{\mu_B\otimes\mu_C}G
 \le {2\epsilon\over1-\epsilon}E_\nu L^2.
 \quad}                                                   \tag{0.7}
\]

This estimate has no correlation-rank loss.  It works because `G` is
nonnegative on the full product and is exactly zero on the calibrated
joint law.  In the exact-gap case, the lifted endpoint differences

\[
 (B,r(B))-(B',r(B')),\qquad
 (C,\ell(C))-(C',\ell(C'))                              \tag{0.8}
\]

span orthogonal subspaces of `R^{n+1}`.  They lie on complementary
Lorentz quadrics.  If the two half-lengths are constant, this reduces to
the ordinary orthogonal-spherical/Clifford classification.  Large
half-length variation is therefore a genuine additional branch, not a
technical nuisance.

High gap also turns a chord-length bound into a common-center bound.  One
obtains

\[
 \boxed{\quad
 E_\nu L^2\ge c(1-\epsilon)P_Fs^5.
 \quad}                                                   \tag{0.9}
\]

Consequently a positive-mass, single-scale family must have root-mean-square
endpoint distance at least order `s^2`.  If finite endpoint densities remain
a fixed fraction of `q_y(0)`, one-dimensional log-concavity instead gives
`L=O_\delta(s)`, and (0.9) forces `s=O(1)`.  Hence every putative large-scale
family falls into the following quantitative trichotomy:

1. its normals are aligned, in which case covariance already forces
   `P_Fs^3=O(1)`;
2. its endpoints are (lifted) orthogonal-radial, and if half-length variation
   is small this is the usual dimension-free Clifford branch, again
   incompatible with positive mass at large `s`; or
3. it uses at least `cP_Fs^5` normal dimensions and endpoints at distance at
   least `c sqrt(P_F)s^(5/2)`, with very small endpoint density or substantial
   scale/length variation.

The third branch is real at the level of all projection and radial
inequalities.  A padded projective-plane construction below has normal rank
`R much larger than s^4`, ideal gap tending to one, and unbounded absolute
Clifford defect, while its conditional density at the finite endpoints is
`exp[-Theta(R/s^2)]` times its boundary density.  It is not a global
log-concave realization.  It shows exactly why radial projections alone
cannot finish the inverse theorem.

The exact translation-charge identity cannot replace the ideal-gap
hypothesis in (0.7): it fixes a coefficient times a Dirichlet energy, not
the heat-bath normalization.  Under translation charge alone, (0.2) is the
strongest general conclusion available here.

## 1. Setup and the two relevant ray measures

Let `mu` be isotropic and log-concave on `R^n` and let a signed-distance
localization give

\[
 X=Z_y+T_yN_y,\qquad |N_y|=1,\qquad y\sim\eta.          \tag{1.1}
\]

The conditional density `q_y` of `T_y` is log-concave, has variance
`sigma_y^2`, and is balanced at zero:

\[
 \delta\le P(T_y\ge0\mid y)\le1-\delta.                 \tag{1.2}
\]

Assume the ray has two finite endpoints

\[
 B_y=Z_y+r_yN_y,\qquad C_y=Z_y-\ell_yN_y,
 \qquad r_y,\ell_y>0.                                   \tag{1.3}
\]

Along a calibrated ray the signed-distance function satisfies

\[
 f(B_y)=r_y,\qquad f(C_y)=-\ell_y,\qquad
 |B_y-C_y|=r_y+\ell_y.                                  \tag{1.4}
\]

If several rays meet at the same endpoint `B`, the value `f(B)` is unique.
Thus `r_y=r(B_y)` is endpoint-measurable; similarly
`ell_y=ell(C_y)`.

Put `b_y=q_y(0)`.  On a measurable finite-endpoint family `F` define the
absolute boundary weight and its normalized boundary law by (0.1).
If `P=\int b_y\,d\eta` is the total perimeter weight and `F` has boundary
fraction `beta`, then `P_F=beta P`; every estimate below retains this
absolute factor.
The distinction between `eta` and `nu` matters.  The standard one-dimensional
log-concave density estimates give

\[
 {c_\delta\over\sigma_y}\le b_y\le {C\over\sigma_y}.
                                                               \tag{1.5}
\]

The upper bound alone implies

\[
 d\eta={P_F\over b_y}\,d\nu
 \ge cP_F\sigma_y\,d\nu.                               \tag{1.6}
\]

If `alpha=eta(F)` and `s<=sigma_y<=A_\sigma s` on `F`,
then both sides of (1.5) give

\[
 {c_\delta\alpha\over A_\sigma s}
 \le P_F\le {C\alpha\over s}.                           \tag{1.7}
\]

Thus a positive *fraction* of the perimeter is not by itself positive bulk
mass.  If the total perimeter is small, or if the selected rays occur over
many scales, `P_F` can be too small for a dimension-free conclusion.  This
is the first unavoidable weight escape.

## 2. Thin shell for every orthogonal marginal

The key point is that every orthogonal marginal of an isotropic log-concave
law is again isotropic and log-concave on its range.

**Theorem 2.1 (projection budgets).**  Let `P` be an orthogonal projection
of rank `d`.  Under the assumptions above and `sigma_y>=s` on `F`,

\[
 P_Fs^3 E_\nu|PN|^2\le C d,                             \tag{2.1}
\]

and

\[
 P_Fs^5 E_\nu|PN|^4\le C d.                             \tag{2.2}
\]

More locally, fix `z in ran P` and `R>0`, and let

\[
 F(P,z,R)=\{y\in F:|PB_y-z|\le R,\ |PC_y-z|\le R\}.
                                                               \tag{2.3}
\]

Then

\[
 P_Fs^5 E_\nu\!\left[
 |PN|^4{\bf1}_{F(P,z,R)}\right]\le C R^2.               \tag{2.4}
\]

**Proof.**  Total covariance gives

\[
 \int_F\sigma_y^2PN_yN_y^TP\,d\eta(y)\preceq P.
\]

Taking traces, using (1.6), and using `sigma_y>=s` proves
(2.1).

Conditionally on `y`,

\[
 |PX|^2=|PZ_y|^2+2\langle PZ_y,PN_y\rangle T_y
             +|PN_y|^2T_y^2.                           \tag{2.5}
\]

The one-dimensional quadratic anti-concentration lemma says that for a
log-concave variable of variance `sigma^2`,

\[
 \operatorname {Var}(T^2+2aT)\ge{\sigma^4\over100}
 \quad\hbox{for every }a\in R.                          \tag{2.6}
\]

Factoring `|PN_y|^2` in (2.5) therefore gives

\[
 \operatorname {Var}(|PX|^2\mid y)
 \ge{\sigma_y^4|PN_y|^4\over100}.                       \tag{2.7}
\]

The quadratic thin-shell bound for the rank-`d` marginal is

\[
 \operatorname {Var}|PX|^2\le C_{TS}d.                 \tag{2.8}
\]

The law of total variance, (1.6), and `sigma_y>=s` prove
(2.2).

For `y` in (2.3), the whole projected segment lies in the radius-`R` ball
about `z`.  If `Q=|PX-z|^2`, an independent-copy identity gives

\[
 \operatorname {Var}(Q\mid y)
 \le4R^2\operatorname {Var}(|PX-z|\mid y).              \tag{2.9}
\]

Combining (2.7), now centered at `z`, with (2.9) yields

\[
 \operatorname {Var}(|PX-z|\mid y)
 \ge{\sigma_y^4|PN_y|^4\over400R^2}.                   \tag{2.10}
\]

Translated thin shell for the marginal says

\[
 \sup_z\operatorname {Var}|PX-z|\le C_{TS}.             \tag{2.11}
\]

Integrating (2.10) and using (1.6) proves (2.4).  QED.

### 2.1 Rank, Ky Fan, and alignment consequences

Let

\[
 A_N=E_\nu NN^T
\]

and let `lambda_1>=lambda_2>=...` be its eigenvalues.  Applying (2.1) and
(2.2) to the projection onto the first `d` eigenvectors, with Jensen for
the fourth moment, gives

\[
 \boxed{\quad
 \sum_{j\le d}\lambda_j
 \le C\min\left\{
 {d\over P_Fs^3},
 \sqrt{d\over P_Fs^5}
 \right\}.
 \quad}                                                   \tag{2.12}
\]

If `W=span{N_y:y in F}` has dimension `r_N`, take `P=P_W` in
(2.2).  Since `|PN|=1`,

\[
 \boxed{\quad r_N\ge cP_Fs^5.\quad}                     \tag{2.13}
\]

For the covariance effective rank

\[
 r_{\rm eff}(A_N)={\operatorname {Tr}A_N\over\|A_N\|_{op}}
 ={1\over\lambda_1},
\]

the case `d=1` in (2.12) gives the distinct bound

\[
 r_{\rm eff}(A_N)
 \ge c\max\{P_Fs^3,\sqrt{P_Fs^5}\}.                    \tag{2.13a}
\]

Thus (2.13) is a support-rank statement and can be stronger by two powers
of the scale.  In the single-scale situation (1.7), the two conclusions are

\[
 r_N\ge {c_\delta\alpha\over A_\sigma}s^4,
 \qquad
 r_{\rm eff}(A_N)\ge
 c_\delta\sqrt{\alpha\over A_\sigma}\,s^2.             \tag{2.14}
\]

At the opposite extreme, suppose for some unit `v`

\[
 E_\nu|N-v|^2\le\tau.                                   \tag{2.15}
\]

Then `E<N,v>=1-\tau/2`, so Jensen and the rank-one case of (2.1)
give

\[
 \boxed{\quad
 P_Fs^3(1-\tau/2)^2\le C.
 \quad}                                                   \tag{2.16}
\]

Under (1.7), any fixed amount of alignment therefore gives

\[
 \alpha s^2\le {C A_\sigma\over\delta(1-\tau/2)^2}.
                                                               \tag{2.17}
\]

Thus aligned positive-mass rays are a harmless dimension-free branch; they
cannot carry `s` to infinity.

## 3. High ideal gap with variable endpoint half-lengths

Let `gamma=(B,C)_#nu` and let `mu_B,mu_C` be its marginals.  The ideal
two-endpoint form is

\[
 {\cal D}_{HB}(h)
 =E_\nu\operatorname {Var}(h\mid B)
  +E_\nu\operatorname {Var}(h\mid C).                  \tag{3.1}
\]

Assume

\[
 {\cal D}_{HB}(h)\ge(1-\epsilon)\operatorname {Var}_\nu h
 \qquad(E_\nu h=0).                                    \tag{3.2}
\]

Equivalently, the maximal correlation between `B` and `C` is at most
`epsilon`.  No assertion is made here that the geometric medial form equals
(3.1).

Define the lifted endpoints

\[
 \widetilde B=(B,r(B)),\qquad
 \widetilde C=(C,\ell(C))\in R^{n+1}.                  \tag{3.3}
\]

Assume their second moments are finite; otherwise infinite endpoint moment
is already the far-endpoint escape.

**Theorem 3.1 (dimension-free lifted product rigidity).**  Let

\[
 G(b,c)=|b-c|^2-(r(b)+\ell(c))^2.                       \tag{3.4}
\]

Then `G>=0` on `supp(mu_B) times supp(mu_C)` and

\[
 \Delta:=E_{\mu_B\otimes\mu_C}G
 \le2\epsilon
 \sqrt{\operatorname {Tr}\Sigma_{\widetilde B}\,
       \operatorname {Tr}\Sigma_{\widetilde C}}.        \tag{3.5}
\]

Moreover,

\[
 \operatorname {Tr}\Sigma_{\widetilde B}
 +\operatorname {Tr}\Sigma_{\widetilde C}
 \le {2-\epsilon\over1-\epsilon}E_\nu L^2,             \tag{3.6}
\]

and hence

\[
 \boxed{\quad
 \Delta\le {2\epsilon\over1-\epsilon}E_\nu L^2.
 \quad}                                                   \tag{3.7}
\]

For independent copies `B,B'` and `C,C'`,

\[
 E\left|
 \langle\widetilde B-\widetilde B',
         \widetilde C-\widetilde C'\rangle
 \right|
 \le2\Delta
 \le {4\epsilon\over1-\epsilon}E_\nu L^2.              \tag{3.8}
\]

**Proof.**  Global one-Lipschitzness gives

\[
 r(b)+\ell(c)=f(b)-f(c)\le|b-c|,
\]

so `G>=0` on the full product.  Calibration gives `G(B,C)=0` under
`gamma`.

Write

\[
 Q(b,r)=|b|^2-r^2.
\]

Then

\[
 G(b,c)=Q(b,r(b))+Q(c,\ell(c))
          -2\langle\widetilde b,\widetilde c\rangle.    \tag{3.9}
\]

The first two terms have the same expectations under `gamma` and under the
product of its marginals.  Therefore

\[
 \Delta
 =2E_\gamma\left\langle
 \widetilde B-E\widetilde B,\,
 \widetilde C-E\widetilde C
 \right\rangle.                                        \tag{3.10}
\]

Maximal correlation at most `epsilon`, applied coordinate by coordinate
and summed with Cauchy--Schwarz, proves (3.5).

For the spatial endpoint covariances, maximal correlation gives

\[
 \left|E\langle B-EB,C-EC\rangle\right|
 \le{\epsilon\over2}
 \left(E|B-EB|^2+E|C-EC|^2\right).                     \tag{3.11}
\]

Expanding `E|B-C|^2` yields

\[
 E|B-EB|^2+E|C-EC|^2
 \le {E L^2\over1-\epsilon}.                            \tag{3.12}
\]

Also

\[
 \operatorname {Var}r+\operatorname {Var}\ell
 \le E(r^2+\ell^2)\le E L^2.                            \tag{3.13}
\]

Equations (3.12)--(3.13) prove (3.6); the arithmetic-geometric mean in
(3.5) proves (3.7).

Finally, taking the double difference of (3.9) gives

\[
 G(b,c)+G(b',c')-G(b,c')-G(b',c)
 =-2\langle\widetilde b-\widetilde b',
            \widetilde c-\widetilde c'\rangle.          \tag{3.14}
\]

Since all four slacks are nonnegative and have mean `Delta`, (3.8)
follows.  QED.

### 3.1 Exact and approximately constant-length geometry

If `epsilon=0`, (3.7) and nonnegativity imply `G=0` for almost every
product pair.  Equation (3.14) then says that

\[
 U=\operatorname {span}(\widetilde B-\widetilde B),
 \qquad
 V=\operatorname {span}(\widetilde C-\widetilde C)
\]

are orthogonal subspaces of `R^{n+1}`.  Fixing
`widetilde c_0 in supp(mu_C)` in (3.9) shows that the positive lifted
endpoint set lies on one level set of

\[
 Q(\widetilde b)-2\langle\widetilde b,\widetilde c_0\rangle;
                                                               \tag{3.15}
\]

the negative set satisfies the analogous equation.  These are complementary
Lorentz quadrics on the orthogonal affine spans.  This is the exact
classification with variable half-lengths.

If `r` and `ell` are constant, the last coordinates disappear from the
difference spans.  The quadrics become Euclidean spheres in orthogonal
spatial subspaces, recovering the standard orthogonal-radial/Clifford
branch.

There is a useful quantitative version.  Suppose

\[
 \operatorname {Var}r+\operatorname {Var}\ell
 \le\tau^2 E_\nu L^2,                                  \tag{3.16}
\]

and put `d_0=Er+Eell`.  Then (3.7)--(3.8) imply

\[
 E\left|
 \langle B-B',C-C'\rangle
 \right|
 \le\left({4\epsilon\over1-\epsilon}+\tau^2\right)
 E_\nu L^2,                                             \tag{3.17}
\]

and

\[
 E_{\mu_B\otimes\mu_C}
 \left||B-C|^2-d_0^2\right|
 \le C\left({\epsilon\over1-\epsilon}+\tau\right)
 E_\nu L^2.                                             \tag{3.18}
\]

Indeed, subtract the product
`(r-r')(ell-ell')` from the lifted inner product in (3.8), and use
Cauchy--Schwarz.  For (3.18), write the defect as `G` plus
`(r+ell)^2-d_0^2` and use (3.16).

Equations (3.17)--(3.18) are a dimension-free `L^1` version of
orthogonal-radial rigidity.  An `L^2` cross-distance theorem still needs
bounded defect or effective correlation rank; high gap alone cannot provide
it.

## 4. High gap plus short active endpoints forces bounded scale

The projection inequality becomes stronger after high gap places the
endpoints in one deterministic ball.

**Theorem 4.1 (far-endpoint necessity).**  Under (3.2),

\[
 \boxed{\quad
 E_\nu L^2\ge c(1-\epsilon)P_Fs^5.
 \quad}                                                   \tag{4.1}
\]

Equivalently, if

\[
 E_\nu L^2\le A_L^2s^2,                                 \tag{4.2}
\]

then

\[
 \boxed{\quad
 P_Fs^3\le {C A_L^2\over1-\epsilon}.
 \quad}                                                   \tag{4.3}
\]

**Proof.**  Let `b_0=E B`, `c_0=E C` and
`z=(b_0+c_0)/2`.  The same maximal-correlation calculation as (3.12)
gives

\[
 E\bigl(|B-z|^2+|C-z|^2\bigr)
 \le {E L^2\over1-\epsilon}.                            \tag{4.4}
\]

Under (4.2), Markov's inequality shows that a set of `nu`-mass at least
one half has both endpoints in the ball about `z` of squared radius

\[
 R^2={2A_L^2s^2\over1-\epsilon}.                        \tag{4.5}
\]

Apply (2.4) with `P=I`.  On this set `|PN|=1`, so

\[
 {1\over2}P_Fs^5\le CR^2.
\]

This is (4.3), and rearranging gives (4.1).  QED.

In the single-scale situation (1.7), (4.3) gives

\[
 \alpha s^2
 \le {C A_\sigma A_L^2\over
          \delta(1-\epsilon)}.                          \tag{4.6}
\]

Thus fixed actual mass and `L=O(s)` force `s=O(1)`.

### 4.1 Endpoint density records the only way to make the chord much longer

Let

\[
 e_y^+=q_y(r_y-),\qquad e_y^-=q_y(-\ell_y+).            \tag{4.7}
\]

The pointwise tail bound for a one-dimensional log-concave density, together
with balance at zero, gives

\[
 {e_y^+\over b_y}\le C_\delta
       \exp\left(-c_\delta {r_y\over\sigma_y}\right),
 \qquad
 {e_y^-\over b_y}\le C_\delta
       \exp\left(-c_\delta {\ell_y\over\sigma_y}\right).
                                                               \tag{4.8}
\]

Consequently, if both finite endpoints are `kappa`-active,

\[
 e_y^+\ge\kappa b_y,\qquad e_y^-\ge\kappa b_y,          \tag{4.9}
\]

then

\[
 L_y\le C_\delta\sigma_y
             \log{C_\delta\over\kappa}.                 \tag{4.10}
\]

If also `sigma_y<=A_\sigma s` throughout the family, take

\[
 A_L=C_\delta A_\sigma\log(C_\delta/\kappa)
\]

in (4.6).  For fixed `alpha,delta,A_\sigma,kappa` this bounds `s`
universally.

Conversely, (4.1) says that a large-scale family must have RMS chord length
at least

\[
 (E_\nu L^2)^{1/2}
 \ge c\sqrt{(1-\epsilon)P_F}\,s^{5/2}.                 \tag{4.11}
\]

When `P_F asymp alpha/s`, this is `c sqrt(alpha)s^2`.  On a
single-scale family, (4.8) then forces at least one of the following:

* a substantial long-chord tail, rather than uniform endpoint control;
* an endpoint density exponentially smaller than `b_y`; or
* a failure of the single-scale assumption through rays with
  `sigma_y much larger than s`.

The endpoint densities in (4.7) are exactly the finite-endpoint terms in
the one-dimensional smooth/focal charge.  Thus the far-endpoint escape also
makes those charges negligible.  Treating every finite endpoint as an
equally weighted heat-bath block loses this information.

## 5. The quantitative trichotomy

The preceding results can be packaged without hiding any missing
hypothesis.

**Theorem 5.1 (radial endpoint trichotomy).**  Fix a family `F` as in
Section 1 and assume the ideal gap (3.2).  Then:

1. **Aligned branch.**  If (2.15) holds, then (2.16) holds.  For positive
   actual mass in one scale band this forces `s=O(1)`.

2. **Lifted orthogonal-radial branch.**  The full product slack and lifted
   cross-angle defects obey (3.7)--(3.8), independently of dimension and
   correlation rank.  If the half-length variation obeys (3.16), then the
   ordinary spatial endpoint laws satisfy the dimension-free approximate
   Clifford estimates (3.17)--(3.18).  If, in addition, endpoints are active
   and the family is single-scale, (4.6)--(4.10) force `s=O(1)`.

3. **High-rank/far-endpoint branch.**  Any remaining large-`s` family must
   satisfy simultaneously

   \[
    r_N\ge cP_Fs^5,\qquad
    E_\nu L^2\ge c(1-\epsilon)P_Fs^5.                  \tag{5.1}
   \]

   In a positive-mass single-scale family these become

   \[
    r_N\ge c_{\delta,A_\sigma}\alpha s^4,\qquad
    E_\nu L^2\ge
       c_{\delta,A_\sigma}(1-\epsilon)\alpha s^4.       \tag{5.2}
   \]

   It must also exhibit low endpoint density, large half-length variation,
   or scale spread as described after (4.11).

This is the strongest conclusion obtainable from orthogonal-marginal thin
shell and ideal endpoint mixing alone.  Notice that ambient padding outside
`span N` is absent from (5.1); only active normal rank counts.

## 6. A padded incidence survivor with rank much larger than `s^4`

The high-rank/far-endpoint branch is not empty at the level of all the
preceding necessary inequalities.

Let `Q` be a prime power, let

\[
 M=Q^2+Q+1,\qquad m=M-1=Q(Q+1),                         \tag{6.1}
\]

and use the projective-plane endpoint vectors `b_i,c_ell` in an
`m`-dimensional space `H`:

\[
 |b_i|^2=|c_\ell|^2=m,
\]

\[
 |b_i-c_\ell|^2=
 \begin{cases}
 D_0^2=2m-2Q^{3/2},&i\in\ell,\\
 D_0^2+2M/\sqrt Q,&i\notin\ell.
 \end{cases}                                            \tag{6.2}
\]

Let `U,V` be two further orthogonal `m`-dimensional spaces.  For independent
uniform `u in S(U)` and `v in S(V)`, define

\[
 B=(b_i,\sqrt m\,u,0),\qquad
 C=(c_\ell,0,\sqrt m\,v).                               \tag{6.3}
\]

The joint endpoint law chooses a uniform incidence `i in ell` and independent
`u,v`.  Its marginals choose the point, line, and sphere variables
independently.

Every joint edge has the same squared length

\[
 D^2=D_0^2+2m=4m-2Q^{3/2}\asymp m,                     \tag{6.4}
\]

whereas every nonedge is longer by `2M/sqrt Q`.  Hence assigning endpoint
values `D/2` and `-D/2` is one-Lipschitz and calibrates every joint edge.
The independent sphere padding adds no endpoint correlation, so the ideal
gap remains

\[
 \lambda_{HB}=1-{\sqrt Q\over Q+1}.                    \tag{6.5}
\]

The padding also does not repair the shared projective-plane covariance
block.  Under the product endpoint law,

\[
 \operatorname {Var}|B-C|^2=4m,\qquad
 4\operatorname {Tr}(\Sigma_B\Sigma_C)=4m.             \tag{6.6}
\]

Thus the absolute Clifford defect diverges.  The nonnegative mean product
slack from Theorem 3.1 is exactly

\[
 E_{\mu_B\otimes\mu_C}G
 =\left(1-{Q+1\over M}\right){2M\over\sqrt Q}
 =2Q^{3/2}\asymp \epsilon D^2.                          \tag{6.6a}
\]

Hence the dimension-free `L^1` estimate (3.7) is sharp in scale, while its
absolute error can be arbitrarily large when the endpoint chord is far
beyond the conditional standard deviation.

For an edge put `N=(B-C)/D`.  The old incidence covariance estimate and
the spherical covariance identities give

\[
 E NN^T\preceq {C\over m}I_{H\oplus U\oplus V}.         \tag{6.7}
\]

In particular its covariance effective rank, as well as its support rank,
is at least `cm`.

The corresponding midpoint-plus-conditional mixture also satisfies
`Cov(X)\preceq I` for `Q` large and `m\gg s^2`: its midpoint covariance is
at most `(1+O(Q^{-1/2}))/2` on `H` and `1/4` on each of `U,V`, while the
conditional covariance is bounded by `Cs^2/m`.

Consequently, for every rank-`d` projection,

\[
 E|PN|^4\le E|PN|^2\le {Cd\over m}.                    \tag{6.8}
\]

Now put on every chord the symmetric truncated Gaussian

\[
 q_{s,D}(t)=Z^{-1}\exp\left(-{t^2\over2a^2}\right)
 {\bf1}_{[-D/2,D/2]}(t),                               \tag{6.9}
\]

where `a` is chosen so that the variance is exactly `s^2`.  When
`D/s` tends to infinity, `a asymp s`,

\[
 q_{s,D}(0)\asymp {1\over s},\qquad
 {q_{s,D}(D/2-)\over q_{s,D}(0)}
 =\exp\left[-\Theta\left({m\over s^2}\right)\right].
                                                               \tag{6.10}
\]

If the whole ray family has quotient mass one, then `P_F asymp1/s`.
Choose

\[
 m\gg s^4.                                               \tag{6.11}
\]

Equations (6.8) and (6.11) give, for every `P`,

\[
 P_Fs^5E|PN|^4\lesssim {s^4d\over m}\le Cd,             \tag{6.12}
\]

so all global projection budgets pass.  The local radial tests pass as
well.  If both projected endpoints lie in a radius-`R` ball, then

\[
 D|PN|=|PB-PC|\le2R.                                    \tag{6.13}
\]

For `R<=D` this gives

\[
 P_Fs^5E[|PN|^4{\bf1}_{F(P,z,R)}]
 \lesssim{s^4R^4\over D^4}\le C R^2,                   \tag{6.14}
\]

because `D^2 asymp m>=s^4`.  For `R>=D` the same conclusion follows
from `P_Fs^5 asymp s^4<=D^2<=R^2`.  Thus (2.4) also sees no
contradiction.

All endpoints have radius `sqrt(2m)`.  The equal-radius radial lower bound
is only

\[
 {\sigma^4\over C m}\asymp{s^4\over m}=o(1),            \tag{6.15}
\]

which is compatible with translated thin shell.  Taking, for example,
`m=s^{10}` gives active normal rank much larger than `s^4`.  Taking
`m>=exp(Cs)` also makes the covariance bound on every fixed coherent cap
smaller than `exp(-cs)`; tensorizing the incidence core gives the same
qualitative escape.

The mixture in (6.3)--(6.9) is not asserted to be one full-dimensional
log-concave measure or a signed-distance extremizer.  It has the correct
calibration, balance, ideal gap, covariance capacity, endpoint-density
degeneracy, and every projection/radial scaling derived above.  Therefore
no argument using only those inequalities can eliminate the
high-rank/far-endpoint branch.

## 7. Why exact translation charge does not yield this endpoint theorem

The smooth-plus-junction translation identity controls

\[
 {\cal S}_N+
 {1\over2}\int|N-N'|^2\,dW(N,N')
 =P-{K\over2}.                                         \tag{7.1}
\]

It does not control the scalar mass of `W` or identify the fiber
conductances with conditional endpoint probabilities.  On a two-sheet
junction with angular gap `theta`, the geometric heat-bath coefficient is
of order `1/theta` while its contribution to (7.1) is of order `theta`.
More abstractly, multiplying any endpoint Dirichlet form by an arbitrary
constant changes the Markov time normalization without changing the form
after the inverse rescaling.

Concretely, if an ideal heat bath has any gap `lambda in (0,1)`, then

\[
 {\cal J}={P\over\lambda}{\cal D}_{HB}                  \tag{7.2}
\]

satisfies `J(h) >= P Var(h)` for every centered
height, and it saturates that bound on a bottom eigenfunction, while the
unscaled endpoint gap remains the arbitrary number `lambda`.  A single
translation trace cannot distinguish these normalizations.

Consequently (7.1), even together with stability, cannot imply
`lambda_HB>=1-epsilon`.  The projection inequalities of Section 2 remain
valid because they use only isotropy, log-concavity, and the actual ray
conditionals.  The product-slack conclusions of Section 3 genuinely require
the ideal endpoint-gap hypothesis or a new comparison between geometric
conductance and the boundary-weighted conditional laws.

## 8. Final verdict

Thin shell for every orthogonal marginal does substantially sharpen the
radial exclusion:

* it replaces ambient dimension by active normal rank;
* boundary weighting produces the exact scale `P_Fs^5`;
* high ideal endpoint gap forces both product near-calibration in the lifted
  geometry and RMS endpoint distance at least `c sqrt(P_F)s^(5/2)`;
* active finite endpoints and a single variance scale collapse this to
  `s=O(1)`.

What remains is an explicit high-rank/far-endpoint escape.  It has at least
order `s^4` normal dimensions at positive bulk mass, chord length of order
`s^2` or larger, and exponentially weak endpoint density.  Padded incidence
models realize every numerical feature of that escape while retaining a
near-maximal ideal gap.  Eliminating it requires a global log-concave
realization or a geometric endpoint-conductance theorem; neither follows
from radial projections or exact translation charge alone.
