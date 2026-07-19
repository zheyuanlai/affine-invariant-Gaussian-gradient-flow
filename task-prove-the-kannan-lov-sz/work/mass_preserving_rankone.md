# Mass-preserving rank-one localization

This note isolates what is rigorous in the control

```
C_t=P_{v_t^\perp},\qquad
v_t=\operatorname{Cov}_{\mu_t}(1_S,X),\qquad \mu(S)=1/2,
```

and what remains conjecture-strength.  All assertions below are on the
minimal affine hull of the initial law.

## 1. An anisotropic Cheeger lemma with one flat direction

### Proposition 1

Let `mu=exp(-V)dx` be a full-dimensional log-concave probability on
`R^n`.  Suppose, in the sense of distributions,

```
Hess V >= alpha P_{u^perp}
```

for a unit vector `u` and `alpha>0`.  Then every Borel set `S` with
`mu(S)=1/2` satisfies

```
mu^+(S) >= c_0 / sqrt(alpha^{-1}+Var_mu <u,X>),       (1)
```

where `c_0>0` is numerical.  The same conclusion holds for a log-concave
law on an affine subspace, with all objects interpreted intrinsically.

### Proof

Apply the Euclidean localization/needle-decomposition theorem to

```
h=1_S-1/2.
```

In the integrable-function version of that theorem there are a probability
space `(Omega,pi)`, a measurable partition, modulo null sets, into line
segments

```
I_omega={z_omega+t w_omega:t in J_omega}, |w_omega|=1,
```

and conditional probabilities `mu_omega` such that

```
mu=int mu_omega pi(domega),       mu_omega(S)=1/2       (2)
```

for `pi`-almost every nontrivial needle.  In Euclidean space the conditional
density in arclength has the form

```
rho_omega(t)=exp(-V(z_omega+t w_omega))J_omega(t)/Z_omega,
```

where `log J_omega` is concave.  This is the curvature-preserving form of
the needle theorem (equivalently, the conditional measures satisfy the
one-dimensional `CD(0,infinity)` inequality).  Consequently, writing
`W_omega=-log rho_omega`,

```
W_omega'' >= alpha(1-<u,w_omega>^2)                  (3)
```

distributionally.  The theorem applies directly to integrable `h`; one may
alternatively approximate `h` in `L^1(mu)` by bounded continuous functions
and use tightness of the disintegrations.

Put

```
tau_omega^2=Var_{mu_omega}(t),
s_omega^2=Var_{mu_omega}<u,X>
             =<u,w_omega>^2 tau_omega^2.             (4)
```

If `<u,w_omega>^2<1/2`, the one-dimensional Brascamp--Lieb inequality
applied to (3) gives `tau_omega^2<=2/alpha`.  If
`<u,w_omega>^2>=1/2`, (4) gives
`tau_omega^2<=2s_omega^2`.  Thus, in all cases,

```
tau_omega^2 <= 2(alpha^{-1}+s_omega^2).              (5)
```

For a one-dimensional log-concave probability of variance `tau^2`, every
half-mass Borel set has exterior perimeter at least `c_1/tau`.  Here is a
constant-only verification.  One-dimensional log-concave isoperimetry says
that a half-line minimizes perimeter.  If `m` is a median and `rho` the
density, log-concavity gives `||rho||_infinity<=2rho(m)`.  Among probability
densities bounded by `K`, the centered interval density `K 1_I` minimizes
variance, so `tau^2>=1/(12K^2)`.  Hence
`rho(m)>=1/(sqrt(48)tau)` and one can take `c_1=1/sqrt(48)`.

Intrinsic enlargement on a needle is contained in Euclidean enlargement:

```
(S cap I_omega)^{I_omega}_epsilon
       subset S_epsilon cap I_omega.
```

Disintegrating and applying Fatou to the nonnegative difference quotients
therefore gives

```
mu^+(S) >= int mu_omega^+(S cap I_omega) pi(domega)
         >= (c_1/sqrt(2)) int
              (alpha^{-1}+s_omega^2)^{-1/2} pi(domega).   (6)
```

The map `r -> (alpha^{-1}+r)^{-1/2}` is convex and decreasing.  Jensen and
the conditional-variance identity give

```
int s_omega^2 pi(domega)
 <= Var_mu<u,X>.
```

Substitution in (6) proves (1), with `c_0=1/sqrt(96)`.  Approximation by
smooth positive densities and bounded convex supports preserves the
distributional curvature inequality and yields the stated generality; the
perimeter step already uses exterior Minkowski content, so no boundary
regularity of `S` is needed.

## 2. Exact stochastic identities

For a bounded predictable symmetric control `C_t`, put

```
p_t(x)=Z_t^{-1}exp(c_t*x-x^TQ_tx/2)p(x),
dc_t=C_t dW_t+C_t^2a_tdt,       dQ_t=C_t^2dt.
```

After the usual parameter stopping, Ito's formula gives

```
dp_t(x)=p_t(x)<x-a_t,C_tdW_t>.                       (7)
```

Consequently, for `g_t=mu_t(S)` and
`v_t=Cov_t(1_S,X)`,

```
dg_t=v_t^TC_t dW_t.                                 (8)
```

Choose a predictable unit vector `e_t` parallel to `v_t` when `v_t!=0`
and set

```
C_t=I-e_te_t^T.                                     (9)
```

At a zero of `v_t`, any predictable unit vector is algebraically admissible.
Then (8) has zero quadratic variation, so `g_t=1/2` pathwise.  Moreover

```
Q_T=T I-M_T,      M_T=int_0^T e_te_t^Tdt,            (10)
Tr M_T=T.
```

If `u_T` is a measurable top eigenvector of `M_T`, the second eigenvalue of
`M_T` is at most `T/2`; hence

```
Q_T >= (T/2)P_{u_T^perp}.                            (11)
```

The endpoint potential has distributional Hessian at least `Q_T`.
Proposition 1 therefore gives, pathwise,

```
mu_T^+(S) >= c_0/
 sqrt(2/T+Var_{mu_T}<u_T,X>).                        (12)
```

For every fixed `x`, (7) makes `p_t(x)/p(x)` a nonnegative likelihood
martingale.  Thus `E mu_T(B)=mu(B)` for every Borel `B`.  For exterior
Minkowski content, Fatou gives the correctly oriented transfer

```
E mu_T^+(S)
 <= liminf_{epsilon downarrow0}
       E[mu_T(S_epsilon)-mu_T(S)]/epsilon
 =mu^+(S).                                          (13)
```

For a finite-perimeter set and a continuous density, equality holds by
integrating the pointwise likelihood identity over the reduced boundary;
only (13) is needed here.

Combining (12)--(13) reduces the route exactly to

```
E (2/T+Var_{mu_T}<u_T,X>)^{-1/2} >= c_T>0            (14)
```

for one fixed numerical `T`.  It is not enough to bound the posterior
variance in each deterministic direction: the mixture identity gives
`E A_T<=A_0=I` in Loewner order, while `u_T` is selected from the same path
as `A_T`.

### Well-posedness at `v=0`

On the open set `v(c,Q)!=0`, (9) is a smooth bounded coefficient after the
usual moment/parameter stopping, so strong existence and uniqueness hold up
to the first zero of `v`.  There is no canonical projection at `v=0`, and a
claim of global strong well-posedness without a convention is false as a
statement of the problem.

There is, however, a rigorous relaxed weak formulation.  First take `mu`
with compact convex support; then the natural-parameter moment maps and all
their derivatives needed on a finite time interval are bounded.  On a mesh,
hold fixed during `[t_k,t_{k+1})` the projection annihilating `v_{t_k}` (and
choose any rank-`n-1` projection when `v_{t_k}=0`).  The laws are tight.
The semimartingale modulus of continuity of `v`, uniformly on the compact
support, gives

```
E int_0^T |P^{(mesh)}_t v_t|^2dt ->0.                (22)
```

After passage to a subsequence, the projection-valued controls converge in
the relaxed sense to a predictable matrix `R_t` satisfying

```
0<=R_t<=I,  Tr R_t=n-1,  R_tv_t=0                   (23)
```

almost everywhere.  The limiting density martingale is driven by
`R_t^{1/2}` and its quadratic potential obeys `dQ_t=R_tdt`.  Equation (23)
makes (8) identically zero.  Moreover

```
Q_T=T I-M_T,\qquad M_T=int_0^T(I-R_t)dt,
```

where `M_T>=0` and `Tr M_T=T`.  Thus the second eigenvalue of `M_T` is still
at most `T/2`, and (11)--(14) remain unchanged.  A literal projection-valued
weak solution, if it exists, is a special case; uniqueness is not needed for
a bound uniform over all relaxed limits.

For an unbounded log-concave law one must first truncate on increasing convex
compact sets and later pass the final Cheeger estimate to the limit.  This is
the same quantitative approximation issue required by any stochastic-
localization proof; it cannot be bypassed by treating the pointwise
likelihood local martingales as automatically uniformly integrable.
Alternatively one can stop forever at `v=0`; then mass is preserved but the
rank-`n-1` curvature assertion (11) is lost.

## 3. Product-exponential median-sum stress test

Let `Z_1,...,Z_n` be independent rate-one exponentials, put
`G=sum_i Z_i`, let `m_n` be a median of `Gamma(n,1)`, and take

```
S={G>=m_n}.
```

Translation by the all-ones vector makes the product law isotropic.  At time
zero exchangeability gives

```
v_i=E[(Z_i-1)1_S]
    =(1/n)E[(G-n)1_{G>=m_n}]=r_n/n,                 (15)
```

where `r_n/sqrt(n)->1/sqrt(2pi)`.  Hence `e_0` is the normalized all-ones
vector and `|v_0|->1/sqrt(2pi)`; in particular the initial protected
direction is not small.

The martingale part of `v_t` under a mass-preserving control is

```
dv_t=D_tC_tdW_t,
D_t=E_t[(1_S-1/2)(X-a_t)(X-a_t)^T].                 (16)
```

At time zero `D_0=d_n I+b_n 11^T`.  Its eigenvalue on the hyperplane
`1^perp` can be computed exactly from Dirichlet conditional moments:

```
d_n=(D_0)_{ii}-(D_0)_{ij}
   =E[(1_S-1/2)G^2]/[n(n+1)].                       (17)
```

The central limit theorem with uniform integrability gives

```
sqrt(n)d_n -> sqrt(2/pi).                            (18)
```

Therefore the total transverse angular quadratic variation is order one per
unit time, not order `n`:

```
Tr(P_{1^perp}D_0^2P_{1^perp})
       =(n-1)d_n^2 ->2/pi.                          (19)
```

At the time `t_n=asymp1/log n` when the largest of the transverse Brownian
tilts first approaches the critical exponential slope, (19) predicts total
angular displacement only `O(1/sqrt(log n))`.  Coordinatewise Gaussian
maximal estimates give

```
max_i |e_{t_n,i}|=O_P(1/sqrt(n)),                    (20)
```

up to a numerical factor (the largest coordinate perturbation is comparable
to, not larger in order than, the initial `1/sqrt(n)`).  Hence every
coordinate has already accumulated

```
int_0^{t_n}(1-e_{s,i}^2)ds=(1-o_P(1))t_n.            (21)

```

A coordinate posterior can indeed have variance of order `log n` at this
instant, but it is not protected.  By the order-one time on which `e_t` can
rotate through an order-one angle, that coordinate has accumulated
order-one quadratic curvature, which bounds its variance by a numerical
constant.  Thus the proposed `1/log n` winner-locking mechanism for this
particular set has a cancellation: the halfspace covariance is spread over
all `n` coordinates with norm of order one, and feedback cannot rotate it to
one coordinate on the extreme-value time scale.

Equations (17)--(21) are a local-in-time stress test, not a proof of (14).
To make (20)--(21) fully uniform beyond a stopped neighborhood of the
initial parameters would require moment bounds for `D_t`; no such bounds are
being assumed here.  They do, however, rule out the naive claim that the
initial extreme coordinate is immediately locked merely because its scalar
posterior variance reaches `log n`.

## 4. Status of the route

Proposition 1 and the reduction (12)--(14) are dimension free and do not use
KLS.  The missing estimate (14) is an adaptive conditional-variance theorem.
Neither deterministic-direction total variance nor the spectral fact (11)
controls it.  A proof must use the special identity that the localization
filtration is independent of the label `1_S` (because
`P(S|F_t)=1/2` pathwise), together with global log-concavity.  Without such a
new input, replacing `u_T` by a deterministic direction or applying Jensen
to `E A_T<=I` is invalid.
