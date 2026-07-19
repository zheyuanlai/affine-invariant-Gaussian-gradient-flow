# Audit of the averaged-calibration `L^2` estimate

The averaged posterior calibration is a useful face-specific formulation,
but its proposed `L^2` estimate is quantitatively equivalent to the full
Poincare/KLS bound.  This note gives the exact Hilbert-space statement and
tests the construction on the Gaussian model.  It also identifies precisely
where the posterior strong-convexity estimate loses the likelihood
correlation.

## 1. Minimal flux norm

Let `mu` be a full-dimensional log-concave probability.  On mean-zero
`L^2(mu)`, let

```
D f = grad f,
D^* F = -div_mu F,
L=D^*D.
```

All statements can first be made on the closure of the range of `L` and
then read through the spectral calculus.  For a mean-zero source `s`, define

```
beta_mu(s)=inf{||F||_{L^2(mu)}: D^*F=s}.             (1)
```

The orthogonal projection theorem gives

```
beta_mu(s)^2=<s,L^{-1}s>,
F_min=D L^{-1}s.                                    (2)
```

Indeed, `D^*D L^{-1}s=s`; every other feasible field differs from
`D L^{-1}s` by an element of `ker D^*`, which is orthogonal to the closure
of `ran D`.  Consequently

```
sup_{s != 0} beta_mu(s)^2/||s||_2^2=C_P(mu).         (3)
```

This is also valid with value `+infinity` if no Poincare inequality holds.

## 2. Exact equivalence for the posterior-resampling source

Let `K=K_T` be the posterior-resampling operator and `R=I-K`.  The identities
in `work/calibration_decomposition.md` give, on mean-zero functions,

```
exp(-T) I <= R <= I.                                (4)
```

Let

```
M_T(mu)=sup_{h != 0} beta_mu(Rh)/||h||_2
       =||L^{-1/2}R||_{2->2}.                       (5)
```

Then

```
exp(-T) sqrt(C_P(mu)) <= M_T(mu) <= sqrt(C_P(mu)).  (6)
```

The upper bound follows from `||R||<=1`.  For the lower bound, (4) gives
`||R^{-1}||<=exp(T)`, and

```
L^{-1/2}=(L^{-1/2}R)R^{-1}.
```

Taking operator norms and using (3) proves the other half of (6).

Now let `B_T` be any linear choice of averaged posterior fields satisfying

```
D^*B_T=R.                                           (7)
```

Pointwise minimality in (2) gives

```
||B_T||_{2->L^2(mu;R^n)} >= M_T(mu).                (8)
```

Conversely, if `C_P(mu)<infinity`, the choice

```
B_min=D L^{-1}R                                     (9)
```

satisfies (7) and has norm exactly `M_T(mu)`.  Therefore the assertion

```
||B_T h||_2 <= C_T||h||_2                           (10)
```

uniformly over isotropic log-concave `mu` is quantitatively equivalent, up
to the fixed factor `exp(T)`, to KLS in Poincare form.  This does not make
the calibration route circular: a new posterior argument could prove
(10).  It does mean that posterior strong convexity alone has to supply the
entire missing KLS content; the decomposition and Neumann series do not
weaken the load-bearing estimate.

There is an equivalent pseudo-Poincare form.  Since `R` is self-adjoint,

```
M_T=||R L^{-1/2}||,
||Rf||_2 <= M_T ||grad f||_2.                       (11)
```

Conversely, (4) and (11) give

```
||f-Ef||_2 <= exp(T)M_T||grad f||_2.                (12)
```

Thus an attempted coupling proof of (11) is another exact formulation of
the same endpoint.

## 3. What strong convexity proves, and the failed conversion

For a posterior `mu_omega` with curvature matrix `Q_omega`, the gradient
solution of

```
-div_{mu_omega}F_omega=h-mu_omega h
```

obeys the Bochner estimate

```
int F_omega^T Q_omega F_omega dmu_omega
 <= int (h-mu_omega h)^2 dmu_omega.                 (13)
```

Averaging (13) gives a dimension-free anisotropic energy bound.  Applying
conditional Jensen and Cauchy--Schwarz to convert it to Euclidean energy
introduces the matrix

```
M_T(x)=E^x[Q_T^{-1}].                               (14)
```

This matrix is not uniformly bounded.  For the variance-one Laplace signal
`x=L`, the calculation in `work/normalized_localization.md` gives

```
E^L[Q_T^{-1}] asymptotic to L/[lambda exp(T)].       (15)
```

Hence any proof of (10) which first discards the correlation between
`F_omega[h](x)` and `Q_omega` fails already in one dimension.  A successful
calibration proof must use that the posterior flux is small in directions
and at points where the residual `h-mu_omega h` is irrelevant; this is the
same face-specific cancellation exhibited by the product no-go in
`work/likelihood_scalar_no_go.md`.

## 4. Gaussian model

Let `mu=gamma_n`.  Covariance-normalized localization at time `T` has
posterior covariance `exp(-T)I`, and posterior resampling is the
Ornstein--Uhlenbeck operator with correlation

```
rho=1-exp(-T).
```

On the Hermite chaos of degree `k>=1`,

```
K=rho^k,
R=1-rho^k,
L=k.
```

Thus

```
M_T(gamma_n)=sup_{k>=1}(1-rho^k)/sqrt(k).            (16)
```

Writing `delta=1-rho=exp(-T)` and using
`1-(1-delta)^k<=min{1,k delta}` gives

```
c sqrt(delta) <= M_T(gamma_n) <= sqrt(delta)         (17)
```

with a numerical `c>0` (choose `k` comparable to `1/delta` for the lower
bound).  The bound is independent of dimension.  On each chaos the Neumann
sum is exact:

```
[(1-rho^k)/sqrt(k)] sum_{j>=0}rho^{kj}=1/sqrt(k),    (18)
```

which is the standard Gaussian minimal flux.  In particular nonlinear
tests such as a smooth approximation to `max_i x_i` do not create a hidden
dimension loss in the Gaussian model; their Hermite components are all
controlled by (16).

For one-dimensional exponential or Laplace factors, the divergence equation
has a unique finite-energy flux and the usual one-dimensional Hardy
inequality bounds it.  The standard tensorized Poincare inequality therefore
bounds the *minimal* flux (9) for product exponentials with a universal
constant, including for max-coordinate tests.  This does not by itself show
that an arbitrary posterior-by-posterior choice in (3) averages to the
minimal field.  These solvable models confirm that the face-specific source
has the right scaling, but they do not address the irreducible case: by (6),
doing so uniformly is exactly the remaining KLS estimate.
