# Covariance-normalized stochastic localization: exact identities and a transfer obstruction

## 1. Stopped construction

Let `p` be a full-dimensional log-concave density with finite second moment.
For an adapted symmetric positive-definite matrix process `C_t`, put

```
p_t(x) = Z_t^{-1} exp(c_t dot x - x^T Q_t x/2) p(x),
d c_t = C_t dW_t + C_t^2 a_t dt,
d Q_t = C_t^2 dt,
a_t = int x p_t(x) dx,
A_t = int (x-a_t)(x-a_t)^T p_t(x) dx.
```

For rigor, first stop at

```
tau_R = inf{t: ||A_t|| + ||A_t^{-1}|| + |a_t| + |c_t| + ||Q_t|| >= R}
```

and set `C_t=0` after `tau_R`. On the stopped region all coefficients are
locally bounded. The usual finite-dimensional parameter SDE for `(c,Q)` is
well posed (one can first take a smooth compactly supported positive density;
all identities below are stopped identities and survive monotone removal of
the stop whenever the process is nonexplosive).

Ito's formula gives, pointwise in `x`,

```
d p_t(x) = p_t(x) (x-a_t)^T C_t dW_t.                 (1)
```

Consequently `E p_t(x)=p(x)` and `mu_t(S)` is a martingale for every Borel
set `S`.

In the covariance-normalized choice

```
C_t=A_t^{-1/2},                                      (2)
```

the covariance equation is

```
dA_t = T_t(A_t^{-1/2}dW_t) - A_t dt,                (3)
```

where `T_t(u)=int (x-a_t)(x-a_t)^T <x-a_t,u> p_t(x)dx`.
Also

```
Q_T = int_0^T A_s^{-1}ds,
Q_T^{-1} <= T^{-2} int_0^T A_s ds.                  (4)
```

The last assertion is the matrix arithmetic-harmonic mean inequality,
obtained from operator convexity of inversion applied to the average of
`A_s^{-1}`.

## 2. The set-mass estimate does work

Let `g_t=mu_t(S)` and `v_t=Cov_t(1_S,X)`. Equations (1)--(2) give

```
dg_t = v_t^T A_t^{-1/2}dW_t,
d<g>_t/dt = v_t^T A_t^{-1}v_t <= g_t(1-g_t).         (5)
```

The inequality is exactly Bessel's inequality in `L^2(mu_t)`: the squared
norm of the orthogonal projection of `1_S-g_t` on the span of the centered
coordinate functions is `v_t^T A_t^{-1}v_t`.

Applying Ito to `g_t(1-g_t)` and using (5),

```
E[g_T(1-g_T)] >= exp(-T) g_0(1-g_0),
E min(g_T,1-g_T) >= (exp(-T)/2) min(g_0,1-g_0).      (6)
```

This estimate is dimension free and remains true for the stopped process.

## 3. Endpoint anisotropic isoperimetry

Conditioned on the localization path, `p_T=exp(-V_T)` has generalized
Hessian at least `Q_T`. After the affine change `y=Q_T^{1/2}x`, it is
1-strongly log-concave. The standard strongly-log-concave Cheeger theorem
therefore yields

```
P_{Q_T^{-1},mu_T}(S)
 := int_{partial* S} sqrt(N_S(x)^T Q_T^{-1}N_S(x)) p_T(x)dH^{n-1}(x)
 >= c min(g_T,1-g_T).                                (7)
```

The displayed surface weight follows exactly from the area formula under
`y=Q_T^{1/2}x`.

Combining (6)--(7) would prove KLS if one had the averaged transfer estimate

```
E P_{Q_T^{-1},mu_T}(S) <= C_T P_mu(S)                (8)
```

for all finite-perimeter sets, with `C_T` universal for one fixed `T>0`.
The next sections show that (8) is false already in one dimension.

## 4. Exact posterior/change-of-measure identities

Write `L_t(x)=p_t(x)/p(x)`. For every fixed `x`, (1) says that `L_t(x)` is
the likelihood martingale

```
dL_t(x)=L_t(x)(x-a_t)^T C_t dW_t.                   (9)
```

Let `P^x` be the path law tilted by `L_t(x)`. Girsanov gives a `P^x`-Brownian
motion

```
dB_t=dW_t-C_t(x-a_t)dt.
```

Under `P^x`, the localization parameters satisfy the exact filtering
identities

```
dc_t=C_t dB_t+C_t^2x dt,
dQ_t=C_t^2dt,
d(c_t-Q_tx)=C_t dB_t.                               (10)
```

The posterior mean has an especially simple conditional equation. Since
`da_t=A_tC_t dW_t` under the original path law, (2) and Girsanov give

```
da_t=A_t^{1/2}dB_t+(x-a_t)dt,                        (10a)
a_t=(1-exp(-t))x+exp(-t)a_0
    +int_0^t exp(-(t-s))A_s^{1/2}dB_s.
```

Thus a tail signal pulls the posterior mean a fixed fraction of the way into
the tail in every positive physical time. The exponential example below
shows that the accompanying conditional posterior variance can be of order
the signal distance, even though the unconditional covariance starts at one.

Moreover, as functions of the signal location,

```
grad_x log L_t(x)=c_t-Q_tx,
Hess_x log L_t(x)=-Q_t.                             (11)
```

Thus every density/metric correlation on a boundary is a posterior
expectation. In particular, for every nonnegative path functional `F`,

```
E[p_T(x)F]=p(x) E^x[F].                              (12)
```

Equations (10)--(12), rather than `E p_T=p` alone, are the correlations that
must be retained in any endpoint transfer.

## 5. One-sided exponential counterexample to (8)

Let the initial law be the rate-one exponential density

```
p(z)=exp(-z) 1_{z>=0}.
```

It has variance one (translation by `-1` makes it isotropic in dimension
one). For scalar parameters `(c,q)`, the posterior density is

```
pi_{c,q}(z) proportional to exp((c-1)z-qz^2/2) 1_{z>=0},
a(c,q)=Var_{pi_{c,q}}(Z).
```

With `C_t=A_t^{-1/2}`, under `P^y` equation (10) becomes

```
dq_t=dt/a(c_t,q_t),
dc_t=dB_t/sqrt(a(c_t,q_t))+y dq_t.                  (13)
```

After the quadratic-variation time change there is a standard Brownian
motion `beta_q` for which

```
c(q)=yq+beta_q,
t(q)=int_0^q a(yu+beta_u,u)du.                      (14)
```

The physical-time value `Q_T` is the inverse clock determined by
`t(Q_T)=T`.

This representation also verifies nonexplosion for this example, rather
than assuming it. For `P^y`-almost every Brownian path,
`(yq+beta_q)/q -> y`. For every `y>0`, completing the square shows
`a(yq+beta_q,q)~1/q` as `q->infinity`, and hence `t(q)->infinity`. The clock
therefore has a finite inverse at every finite physical time. Since an
exponential signal has `y>0` almost surely, disintegrating over `y`
constructs the unstopped filter globally; all preceding stopped identities
agree with it before the stops and the stops tend to infinity.

Fix `T>0` and let `y` tend to infinity. Put `q=s/y`. On every bounded
`s`-interval,

```
sup_{0<=s<=M}|beta_{s/y}| -> 0
```

in probability. For `s<1`, bounded away from one, the posterior has a
positive exponential decay rate bounded away from zero, so
`a(s+beta_{s/y},s/y)/y -> 0`. For `s>1`, bounded away from one, complete the
square:

```
pi(z) proportional to exp(-(s/y)(z-y(s-1)/s+o(y))^2/2) 1_{z>=0}.
```

The truncation point is a diverging number of standard deviations below the
mean. Hence, uniformly on compact subintervals of `(1,infinity)`,

```
a(s+beta_{s/y},s/y)/y -> 1/s.                       (15)
```

For the small transition interval around `s=1`, Brascamp--Lieb on the
half-line gives `a(c,s/y)<=y/s`. This dominates the rescaled variance near
one. Splitting the integral into intervals below, around, and above one and
then shrinking the middle interval proves, for every fixed `s>1`,

```
t(s/y)=int_0^s a(r+beta_{r/y},r/y)dr/y -> log s      (16)
```

in probability. (Below `1-delta`, monotone-likelihood-ratio comparison with
an exponential law gives a uniform second-moment bound; above `1+delta`,
(15) is uniform.) Since the clock is increasing, (16) at
`s=exp(T-epsilon)` and `s=exp(T+epsilon)` implies

```
y Q_T -> exp(T)                                     (17)
```

in `P^y`-probability. In particular, for all sufficiently large `y`,

```
E^y[Q_T^{-1/2}] >= c_T sqrt(y).                     (18)
```

No uniform-integrability assertion is needed: on the event
`yQ_T<=2exp(T)`, whose probability tends to one, the integrand is at least
`sqrt(y/(2exp(T)))`.

Now take the half-line `S_y=[y,infinity)`. Its initial exterior perimeter is
`P_mu(S_y)=p(y)`. The one-dimensional anisotropic endpoint perimeter is
`p_T(y)Q_T^{-1/2}`. Using (12) and (18),

```
E P_{Q_T^{-1},mu_T}(S_y) / P_mu(S_y)
 = E[p_T(y)Q_T^{-1/2}]/p(y)
 = E^y[Q_T^{-1/2}]
 >= c_T sqrt(y) -> infinity.                        (19)
```

Thus (8), and also any pointwise boundary-density estimate intended to imply
it, is false even for a variance-one log-concave law on the line. Taking
products and the coordinate halfspace `{x_1>=y}` gives the same obstruction
for the isotropic product of shifted one-sided exponentials in every
dimension.

For comparison, with a standard Gaussian initial law, `A_t=exp(-t)` and
`Q_T=exp(T)-1` deterministically, so the ratio in (19) is
`(exp(T)-1)^{-1/2}`. The obstruction is therefore the posterior covariance
inflation along rare exponential-tail signals, not a defect in the Gaussian
normalization.

## 6. Consequence for the proof search

The normalized driver genuinely solves the balanced-mass survival problem,
but strong convexity at the endpoint is expressed in a random anisotropic
metric whose boundary weight is positively correlated with the posterior
density. This correlation is unbounded on exponential-tail boundary points.
Therefore the direct chain (6) + (7) + (8) cannot prove KLS. A viable repair
would have to separate boundary points already possessing a direct
one-dimensional tail expansion estimate, or use a different endpoint
functional in which the likelihood/metric correlation cancels. Merely
applying (4), Jensen, or an operator norm only worsens the obstruction.

## 7. Restricting to half-mass cuts: exact product recursion

The one-sided example in Section 5 has very small mass. Restricting (8) to
half-mass sets removes that particular one-dimensional counterexample, but it
does not create a new dimension-free mechanism. This can be seen exactly on
the recursive product set from the boundary-stability calculation.

For a measure `theta`, define

```
R_T^theta(D)=E P_{Q_T^{-1},theta_T}(D).
```

Let `nu` be isotropic log-concave on `R^d`, let `nu(B)=1/2`, and put
`p=P_nu(B)`. For standard Gaussian `gamma` on the first coordinate, set

```
A_L = ((-infinity,L] x B) union ((-infinity,-L] x B^c),
alpha_L=Phi(L)-Phi(-L).
```

The covariance-normalized localization of a product remains a product:
the posterior covariance and `C=A^{-1/2}` are block diagonal, and the two
blocks are driven by independent Brownian motions. The initial and endpoint
perimeters therefore factor exactly:

```
P_{gamma tensor nu}(A_L)=phi(L)+alpha_L p,            (20)
R_T^{gamma tensor nu}(A_L)
 = phi(L)/sqrt(exp(T)-1) + alpha_L R_T^nu(B).         (21)
```

For (21), the two horizontal faces contribute the first term because the
Gaussian accumulated precision is deterministically `exp(T)-1` and
`E nu_T(B)=1/2`. The vertical face is `(-L,L) x partial*B`; independence
factors its expectation into `E gamma_T((-L,L))=alpha_L` and
`R_T^nu(B)`.

Letting `L->infinity` in (20)--(21) gives

```
P_{gamma tensor nu}(A_L) -> p,
R_T^{gamma tensor nu}(A_L) -> R_T^nu(B).             (22)
```

Hence a universal half-mass version of (8), even assumed only for these
product sets, already implies the identical endpoint-transfer inequality in
one lower dimension. By (6)--(7), `R_T^nu(B)>=c_T`; consequently (22) would
give `p>=c_T/C`. The term left after removing the Gaussian tail faces is
exactly the lower-dimensional KLS assertion, not a tail estimate.

## 8. Sharp tail weight for an isotropic Laplace coordinate

The Gaussian recursion hides the density/metric correlation because its
precision is deterministic. Replace the first factor by the variance-one
symmetric Laplace law

```
rho_lambda(z)=(lambda/2)exp(-lambda|z|),  lambda=sqrt(2).
```

For a fixed positive signal `L`, the posterior under `P^L` has density

```
pi_{c,q}(z) proportional to
 exp(-lambda|z|+cz-qz^2/2),
c(q)=Lq+beta_q.
```

Putting `q=s/L`, the same clock argument as in Section 5 gives

```
Var_{pi_{s+beta_{s/L},s/L}}(Z)/L -> 0       for s<lambda,
Var_{pi_{s+beta_{s/L},s/L}}(Z)/L -> 1/s     for s>lambda,
t(s/L) -> log(s/lambda)                     for s>lambda,
LQ_T -> lambda exp(T)                       in probability.       (23)
```

Here the first limit uses the two exponential decay rates
`lambda-s` and `lambda+s`; the second follows by completing the square on
the positive half-line, while the negative half-line has negligible mass.
Brascamp--Lieb, `Var<=1/q`, controls the transition around `s=lambda`.

The inverse square roots in (23) are uniformly integrable. Indeed, fix
`r<lambda/4`. If `Q_T<r/L`, then the clock has accumulated time `T` before
`r/L`. On the event

```
sup_{q<=r/L}|c(q)|<=lambda/2,
```

the posterior has two-sided exponential decay rate at least `lambda/2`, so
its variance is at most a universal multiple of `lambda^{-2}` and the clock
time before `r/L` is `O(r/L)<T` for large `L`. Thus `Q_T<r/L` forces
`sup_{q<=r/L}|beta_q|>=lambda/4`, and the reflection principle gives

```
P^L(LQ_T<r) <= 4 exp(-c lambda^2 L/r).                (24)
```

This controls the lower tail strongly enough for uniform integrability.
Therefore the boundary-point weight has the sharp asymptotic

```
w_T(L):=E^L[Q_T^{-1/2}],
w_T(L)/sqrt(L) -> 1/sqrt(lambda exp(T)).              (25)
```

The same holds at `-L` by symmetry. In terms of the one-dimensional upper
tail `delta_L=(1/2)exp(-lambda L)` and its hazard
`rho_lambda(L)/delta_L=lambda`, (25) reads

```
w_T(L) asymptotic to
 (1/lambda) exp(-T/2) sqrt(log(1/(2delta_L))).        (26)
```

Thus a constant hazard does not control the endpoint weight: the missing
factor is the square root of the log-tail depth (equivalently the square root
of the Mahalanobis distance in this example). For a Gaussian the hazard is
of order `L`, and the corresponding weight is constant, consistent with the
heuristic scale `sqrt(L/hazard(L))`; no general multidimensional upper bound
of this form is asserted here.

For the Laplace analogue of `A_L`, with
`alpha_L=eta_lambda((-L,L))=1-exp(-lambda L)` for the measure
`deta_lambda=rho_lambda dx`, product factorization gives

```
P_{rho tensor nu}(A_L)=rho_lambda(L)+alpha_L p,       (27)
R_T^{rho tensor nu}(A_L)
 =rho_lambda(L)w_T(L)+alpha_L R_T^nu(B).              (28)
```

The high-weight horizontal term is, by (25),

```
rho_lambda(L)w_T(L)
 asymptotic to C_T rho_lambda(L)
                 sqrt(log(lambda/(2rho_lambda(L)))). (29)
```

The ordinary one-dimensional tail expansion is only
`rho_lambda(L)=lambda delta_L`; it contains no factor capable of absorbing
the square root in (29). If `rho_lambda(L)` is chosen comparable to a small
transverse perimeter `p`, charging (29) to the initial perimeter loses
`sqrt(log(1/p))`. Removing the tail slabs changes the mass by only
`delta_L`, and the remaining central vertical boundary is precisely the
second term in (28), namely the normalized-localization endpoint functional
for the lower-dimensional cut `B`.

Near-optimality does not by itself remove this tail. Along a hypothetical
sequence with `p->0`, choose `L` so that

```
rho_lambda(L)/p=(log(1/p))^{-1/4}.
```

Then (27) is `(1+o(1))p`, so if `B` is a near-minimizing half-mass cut the
lifted `A_L` is asymptotically just as good (up to the product tensorization
constant), but its horizontal endpoint term divided by `p` is of order
`(log(1/p))^{1/4}`. Thus even an `o(1)`-near-Cheeger restriction does not
justify dropping the likelihood/metric correlation. This statement is
conditional only in the proper proof-by-contradiction sense: if a bad KLS
sequence existed, the proposed variational restriction would still fail to
control its harmless tail decorations.

This is a quantified no-go for the proposed repair: splitting off boundary
points with large posterior weight either incurs a `sqrt(log(1/p))` loss on
Laplace tails, or discards those tails and leaves exactly the
lower-dimensional KLS/endpoint-transfer problem. The calculation does not
rule out a new global cancellation involving both terms of (28), but neither
a pointwise hazard bound nor a direct tail-perimeter estimate supplies such a
cancellation.

## 9. The clipped metric removes the tail correlation but restates the core

There is one exact way to make the endpoint functional transferable. Define

```
bar R_T^mu(S)=E int_{partial* S}
 p_T(x) min{1,sqrt(N^TQ_T^{-1}N)} dH^{n-1}(x).
```

Since the multiplier is at most one and `E p_T(x)=p(x)`, one has without any
uncorrelatedness assumption

```
bar R_T^mu(S)<=P_mu(S).                               (30)
```

For the Laplace tail in Section 8, the conditional multiplier is eventually
one, so clipping changes the amplified term `rho(L)w_T(L)` back to order
`rho(L)`. The price is that strong log-concavity only proves (7) for the
unclipped metric. A dimension-free lower bound

```
bar R_T^mu(S)>=c min(mu(S),1-mu(S))                   (31)
```

would close KLS immediately by (30), and is not supplied by endpoint
curvature.

The product recursion shows exactly what remains. For the Gaussian `A_L`,
put `b_T=min{1,(exp(T)-1)^{-1/2}}`. Then

```
bar R_T^{gamma tensor nu}(A_L)
 =b_T phi(L)+alpha_L bar R_T^nu(B).                  (32)
```

Letting `L->infinity` reduces (31) in dimension `d+1` to (31) for `B` in
dimension `d`. Thus clipping is a clean tail-safe reformulation, but proving
its averaged endpoint lower bound is the same unresolved nonlinear,
lower-dimensional bottleneck. Neither the covariance mixture identity
`E A_T<=A_0` nor the set-mass quadratic-variation upper bound controls the
boundary-conditioned clipped functional; an additional global cancellation
would be required.
