# Posterior calibration decomposition

This note records an exact operator identity suggested by the failure of
scalar-weighted endpoint perimeter.  It does not yet prove the required
uniform flux estimate.

Fix a time `T` in a stochastic localization scheme and write `mu_omega` for
the random posterior.  Define the posterior-resampling kernel

```
K_T h(x) = E^x[ mu_omega h ],
```

where `E^x` is likelihood tilt by `d mu_omega/d mu` at `x`.  Equivalently,
sample the localization observation conditional on `X=x` and then sample an
independent `Y` from its posterior; `K_T h(x)=E[h(Y)|X=x]`.

## Lemma 1 (reversibility and exact Dirichlet form)

`K_T` is a positive self-adjoint Markov operator on `L^2(mu)`, and

```
<h,(I-K_T)h>_mu = E Var_{mu_omega}(h).               (1)
```

Indeed,

```
mu(dx)K_T(x,dy)=E[mu_omega(dx)mu_omega(dy)],
```

which is symmetric, and expanding the right side of (1) proves the identity.

For covariance-normalized localization, Bessel's inequality gives, for
`V_t=Var_{mu_t}(h)`,

```
d E V_t/dt
 = -E[Cov_t(h,X)^T A_t^{-1}Cov_t(h,X)]
 >= -E V_t.
```

Consequently

```
I-K_T >= exp(-T) I                                  (2)
```

on mean-zero `L^2(mu)`, and the spectrum of `K_T` there lies in
`[0,1-exp(-T)]`.

## Lemma 2 (averaging posterior fluxes)

Let `h` have mean zero.  Suppose that for almost every posterior there is a
vector field `F_omega[h]` satisfying

```
-div_{mu_omega} F_omega[h]
   = h-mu_omega h.                                  (3)
```

Define the likelihood average

```
B_T h(x)=E^x[F_omega[h](x)].                         (4)
```

Then, distributionally,

```
-div_mu(B_T h)=(I-K_T)h.                            (5)
```

To prove this, multiply (3) by `mu_omega`, average, and use
`E mu_omega=mu` together with
`E[mu_omega(x)mu_omega h]=mu(x)K_T h(x)`.

It follows from (2), (5), and the Neumann series that

```
F[h] = sum_{j>=0} B_T K_T^j h                       (6)
```

formally satisfies `-div_mu F[h]=h`.

## Exact remaining estimate

If one could choose the posterior calibrations in (3) so that

```
||B_T h||_{L^2(mu;R^n)} <= C_T ||h||_{L^2(mu)}       (7)
```

with `C_T` universal, then (2) and (6) would give

```
||F[h]||_2 <= C_T exp(T)||h||_2.
```

Weighted divergence duality would then yield the dimension-free Poincare
inequality.  Thus (7) is a precise face-specific replacement for the false
scalar perimeter transfer.  Strong convexity of each posterior only gives an
anisotropic bound `int F_omega^T Q_T F_omega`, so converting that bound to
(7), while retaining the likelihood correlation, is the load-bearing open
step.  An `L^infinity` version would imply Cheeger directly, but the `L^2`
version already suffices through T2 and avoids any log-Sobolev claim.
