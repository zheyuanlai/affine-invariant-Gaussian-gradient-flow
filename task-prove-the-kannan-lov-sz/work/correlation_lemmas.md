# Linear predictability of a log-concave cut

This note records a dimension-free fact useful for covariance-normalized
localization.  It is not the central KLS estimate.

## Lemma 1 (a Bernoulli cut cannot be almost linear)

There is a numerical `eta>0` with the following property.  Let `mu` be a
full-dimensional log-concave probability with covariance `A>0`, let `S` be
Borel with `p=mu(S)` in `(0,1)`, and put

```
v = Cov_mu(1_S,X).
```

Then

```
v^T A^{-1}v <= (1-eta) p(1-p).                    (1)
```

The same statement on the affine hull covers singular ambient covariance.

### Proof

Whiten and translate, so `A=I`.  If `v=0` there is nothing to prove.  Put
`u=v/|v|` and `Y=<u,X>`.  The law of `Y` is a one-dimensional log-concave
probability with mean zero and variance one.  A standard one-dimensional
normalization estimate gives

```
||rho_Y||_infinity <= M_0                            (2)
```

for a universal numerical `M_0`.

Let `a=E[Y|S]` and `b=E[Y|S^c]`.  The law of total variance gives

```
1 = p Var(Y|S) + (1-p) Var(Y|S^c)
      + p(1-p)(a-b)^2,

|v|^2/[p(1-p)] = p(1-p)(a-b)^2.                    (3)
```

For every real `y`, the squared error made by the two conditional means is at
least `min{(y-a)^2,(y-b)^2}`.  Consequently

```
p Var(Y|S)+(1-p)Var(Y|S^c)
 >= E min{(Y-a)^2,(Y-b)^2}.                         (4)
```

The union of the two intervals of radius `r` about `a,b` has Lebesgue length
at most `4r`; by (2) it has probability at most `4M_0 r`.  With
`r=(8M_0)^{-1}`, at least half the probability lies outside that union, and
therefore the right side of (4) is at least

```
eta := 1/(128 M_0^2).
```

Substitution in (3) proves (1).  Notice that the argument allows `S` to
depend on all coordinates; no conditional log-concavity of `S` is used.

## Sharpness and limitation

For the isotropic uniform law on an interval and its median half-interval,
the ratio on the left of (1), divided by `p(1-p)`, equals `3/4`.  Thus no
version of (1) can have `eta>1/4`.

In normalized stochastic localization, (1) improves

```
d<g>_t/dt <= g_t(1-g_t)
```

to `(1-eta)g_t(1-g_t)` at every nondegenerate stopped time.  This is a real
dimension-free gain, but by itself it does not close KLS: anisotropic
strong-convexity improves Euclidean perimeter at a square-root curvature
rate, and `1-eta` can exceed `1/2` (the interval example gives `3/4`).  A
successful use must add a separate treatment of the highly linearly
predictable regime or use a sharper likelihood-aware endpoint profile.
