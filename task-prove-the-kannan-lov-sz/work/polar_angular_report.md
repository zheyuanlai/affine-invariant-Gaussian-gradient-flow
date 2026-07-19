# Polar/angular route: exact structure, a radial lemma, and curvature no-go results

Throughout this note the ambient dimension is `n >= 2`, `mu` has full-dimensional
log-concave density `f = exp(-V)`, `int f = 1`, and `X ~ mu` is isotropic.  Write

```
X = R Theta,       R = |X|,       Theta = X/|X| in S^{n-1}.
```

The origin is in the interior of the support of `f`, so the usual upper-semicontinuous
version of `f` satisfies `f(0) > 0`.  All formulas below therefore apply at the
barycentric origin; no translation to a mode is being made.

## 1. Exact disintegration

Let `sigma` denote unnormalized spherical area measure and set

```
w(theta) = int_0^infty r^{n-1} f(r theta) dr.
```

Then `int_S w dsigma = 1`, the law `nu` of `Theta` is

```
dnu(theta) = w(theta) dsigma(theta),
```

and, conditionally on `Theta = theta`, the radial law is

```
q_theta(dr) = w(theta)^{-1} r^{n-1} exp(-V(r theta)) dr.            (1.1)
```

In particular, `q_theta` is a one-dimensional log-concave law on `(0,infty)`.

## 2. The angular marginal is exactly a cone measure

For `p > 0`, define Ball's ray-integral body by its radial function

```
rho_p(theta)^p = (p/f(0)) int_0^infty r^{p-1} f(r theta) dr.         (2.1)
```

Ball's ray-integral theorem says: if `f` is integrable and log-concave and
`f(0)>0`, then the positively homogeneous function whose restriction to the
sphere is `1/rho_p` is convex; equivalently, (2.1) is the radial function of a
convex body `K_p(f)` containing the origin.  These are exactly the hypotheses
verified above.

For `p=n`,

```
rho_n(theta)^n = n w(theta)/f(0),       |K_n(f)| = 1/f(0).          (2.2)
```

The radial projection to `S^{n-1}` of cone measure on `partial K_n(f)` has
density

```
rho_n(theta)^n/(n |K_n(f)|) = w(theta).
```

Thus **the angular marginal of every full-dimensional log-concave law is
exactly the cone measure of the convex Ball body `K_n(f)`**.

Let `G` be the gauge of `K_n(f)`.  Up to an additive constant, the angular
potential is

```
H(theta) := -log w(theta) = n log G(theta) + constant.             (2.3)
```

When `G` is twice differentiable at `theta`, convexity and one-homogeneity give

```
D^2 G(theta)[u,u] = Hess_S G(theta)[u,u] + G(theta)|u|^2 >= 0
```

for `u perpendicular theta`.  Consequently

```
Hess_S H
 = n { G^{-1} Hess_S G - d(log G) tensor d(log G) }
 >= -n { g_S + d(log G) tensor d(log G) }.                         (2.4)
```

This is the convexity that actually survives.  It is one-sided in the wrong
direction for Bakry--Emery, and it is sharp on every linearity chamber of a
polyhedral gauge.

### 2.1 Exact gnomonic-chart formulation

There is a second exact form of the same structure which may be more useful
than the spherical Hessian.  On the hemisphere `theta_1>0`, write

```
x(y)=(1,y) in R x R^{n-1},
s(y)=sqrt(1+|y|^2),
theta(y)=x(y)/s(y).
```

Spherical area has Jacobian `s(y)^{-n}dy`.  The homogeneous extension of `w`
has degree `-n`, so (2.2) gives the exact cancellation

```
dnu(theta(y)) = (f(0)/n) G(1,y)^{-n} dy.                           (2.5)
```

Thus, in every projective chart, angular measure is a critical convex measure

```
phi(y)^{-n}dy,       phi(y)=G(1,y) convex,       dim(y)=n-1.        (2.6)
```

It is generally not log-concave; its `(-1/n)`-power is convex.  The pullback of
the round metric and its inverse are

```
g_y = s^{-2} I - s^{-4} y tensor y,
g_y^{-1} = s^2 (I+y tensor y),                                    (2.7)
```

so, for `h(theta(y))=h_tilde(y)`,

```
|grad_S h|^2
 = (1+|y|^2) { |grad h_tilde|^2 + (y dot grad h_tilde)^2 }.        (2.8)
```

Equations (2.5)--(2.8) turn the desired angular Poincare inequality into a
weighted Poincare inequality for a critical convex measure.  Both the exponent
`n=(n-1)+1` and the metric weight are essential.  For example, in a
max-coordinate chart of the cube, `G(1,y)=1` on `[-1,1]^{n-1}`: the chart law
is ordinary uniform measure on a cube, while the factor in (2.8), typically of
order `n`, supplies the target spherical scaling.  A Euclidean Poincare bound
for `phi^{-n}` alone cannot supply that factor.

Partitioning the sphere (up to null walls) into the `2n` cells

```
C_{j,sign}={theta: sign*theta_j=||theta||_infty}
```

puts every cell into the bounded chart `[-1,1]^{n-1}` with a density of the
form (2.6).  A possible polar continuation would need two estimates with
universal constants: the weighted within-cell inequality dictated by (2.8),
and a gluing inequality across the tie walls.  Neither follows from local
curvature: on the cube the absolutely continuous curvature is negative and
all compensating convexity is carried by those walls.

## 3. Smooth geodesic Hessian formula

Assume temporarily that `V` is `C^2` and that differentiation and radial
integration by parts are justified.  Fix `theta in S^{n-1}` and a unit tangent
vector `u`, and put `theta_s = theta cos(s) + u sin(s)`.  Under `q_theta`, set

```
a(r) = r <grad V(r theta),u>.
```

Twice differentiating `H(theta_s)=-log w(theta_s)` gives

```
Hess_S H(theta)[u,u]
 = E_theta[ r^2 D^2 V(r theta)[u,u] - r <grad V(r theta),theta> ]
   - Var_theta(a(R)).                                               (3.1)
```

Radial integration by parts gives

```
E_theta[ R <grad V(R theta),theta> ] = n,                           (3.2)
```

and hence

```
Hess_S H(theta)[u,u]
 = E_theta[R^2 D^2V(R theta)[u,u]] - n
   - Var_theta(R <grad V(R theta),u>).                              (3.3)
```

The nonnegative transverse Hessian term does not dominate either negative
term.  Formula (3.3) is consistent with (2.4), while (2.4) remains meaningful
for nonsmooth/polyhedral examples.

## 4. A dimension-free radial reduction for arbitrary Lipschitz tests

This part needs neither thin shell nor KLS.

### Lemma 4.1 (weighted radial Brascamp--Lieb)

For every direction for which (1.1) is defined, every locally Lipschitz `h`,
and `n >= 2`,

```
Var_{q_theta}(h)
 <= (1/(n-1)) E_{q_theta}[ R^2 h'(R)^2 ].                           (4.1)
```

Indeed, the potential of `q_theta` is

```
U_theta(r) = V(r theta) - (n-1) log r + constant,
```

so, in the distributional sense,

```
U_theta''(r) >= (n-1)/r^2.
```

The one-dimensional Brascamp--Lieb inequality proves (4.1) in the smooth,
strict case.  Approximation of the convex function `r -> V(r theta)` by smooth
convex functions on compact subintervals, followed by monotone truncation at
`0` and `infty`, proves the displayed inequality in the stated generality.

If `F:R^n -> R` is 1-Lipschitz and

```
g(theta) = E[ F(X) | Theta=theta ],
```

then `|d F(r theta)/dr| <= 1` almost everywhere.  Integrating (4.1) in `theta`
and using isotropy gives the exact estimate

```
E_nu Var(F(X) | Theta)
 <= (1/(n-1)) E|X|^2
 = n/(n-1).                                                         (4.2)
```

Therefore

```
E |F(X)-g(Theta)| <= sqrt(n/(n-1)).                                 (4.3)
```

Thus all radial fluctuation of every 1-Lipschitz test is already controlled by
a universal constant.  The unresolved term is the oscillation of the
conditional mean `g` across directions.  Notice that `g` is not generally
Lipschitz for the round metric: differentiating it introduces the conditional
score covariance from (3.1).

A weaker pointwise consequence, useful for checking examples, is

```
Var_{q_theta}(R) <= E_theta R^2/(n-1).
```

For `n >= 3`, rearranging this also gives

```
Var_{q_theta}(R) <= (E_theta R)^2/(n-2).                            (4.4)
```

The coefficient of variation is small, but the conditional mean and variance
need not be uniformly controlled in `theta`; see Section 6.

## 5. Exact equivalence of round-angular first-moment concentration and T3

Let

```
D_ang(nu) = sup { E_nu |h-E_nu h| : h is 1-Lipschitz for chord distance on S }.
```

Set

```
delta(mu) = ( E (R-sqrt(n))^2 )^{1/2}.
```

### Proposition 5.1

For every isotropic law (log-concavity is not needed for the comparison),

```
| D_1(mu) - sqrt(n) D_ang(nu) | <= 2 delta(mu),                    (5.1)
```

where `D_1(mu)` uses centering by the mean.

Proof.  If `F` is 1-Lipschitz, put `h(theta)=F(sqrt(n) theta)`.  Then `h` is
`sqrt(n)`-Lipschitz and

```
E |F(X)-h(Theta)| <= delta(mu).
```

For any integrable random variables `Y,Z`,

```
E|Y-EY| <= E|Z-EZ| + 2 E|Y-Z|.
```

Taking the supremum over `F` proves

```
D_1(mu) <= sqrt(n) D_ang(nu) + 2 delta(mu).                         (5.2)
```

Conversely, let `h` be chord-1-Lipschitz.  On the sphere of radius `sqrt(n)`,
the function `sqrt(n) theta -> h(theta)` is `1/sqrt(n)`-Lipschitz.  Its McShane
extension `H` to `R^n` has the same Lipschitz constant.  Comparing `H(X)` with
`H(sqrt(n)Theta)=h(Theta)` gives

```
D_ang(nu) <= D_1(mu)/sqrt(n) + 2 delta(mu)/sqrt(n),                 (5.3)
```

which completes the proof.

The dimension-free thin-shell theorem gives `delta(mu) <= C_TS`.  Therefore

Indeed, since `E R^2=n`, one has the exact identity

```
delta(mu)^2
= 2 sqrt(n)(sqrt(n)-E R)
= 2 sqrt(n) Var(R)/(sqrt(n)+E R)
<= 2 Var(R),                                                        (5.4a)
```

so the usual dimension-free thin-shell variance bound has precisely the
required form.  Therefore

```
D_1(mu) <= C   for every isotropic log-concave mu
```

is equivalent, up to an additive universal constant, to

```
D_ang(nu) <= C/sqrt(n) for every corresponding angular marginal.  (5.4)
```

In particular, a round-spherical Poincare inequality

```
Var_nu(h) <= (C/n) int |grad_S h|^2 dnu                            (5.5)
```

would prove KLS immediately, but (5.1) shows that the first-moment content of
such an assertion is precisely conjecture-strength.  Polar decomposition by
itself does not turn the missing estimate into a weaker theorem.  The known
`sqrt(log n)` KLS bound translates into the same `sqrt(log n)` loss in (5.4).

## 6. Explicit isotropic stress tests

### 6.1 Isotropic cube: negative curvature of order `n^2`

Let `mu` be uniform on

```
K = [-sqrt(3),sqrt(3)]^n.
```

It is isotropic.  Its radial endpoint and angular potential are

```
rho(theta) = sqrt(3)/||theta||_infty,
H(theta) = n log ||theta||_infty + constant.                        (6.1)
```

On the smooth chamber where `theta_j` is the unique coordinate of largest
absolute value (fix its sign), every unit tangent `u` satisfies

```
Hess_S H(theta)[u,u] = -n (1 + (u_j/theta_j)^2).                   (6.2)
```

Since `Ric_{S^{n-1}}=(n-2)g`, the Bakry--Emery tensor is at most `-2g` even in
directions with `u_j=0`.  There is no nonnegative angular curvature after
isotropization.

The failure is much larger than order `n`.  Take

```
theta_1 = sqrt(2/n),
theta_2 = ... = theta_n = sqrt((n-2)/(n(n-1))),
u = (e_1-theta_1 theta)/sqrt(1-theta_1^2).
```

Then `theta_1` is the unique maximum coordinate, `u` is a unit tangent, and

```
(u_1/theta_1)^2 = n/2-1,
Hess_S H(theta)[u,u] = -n^2/2.                                    (6.3)
```

Also `|grad_S H(theta)|` is of order `n^{3/2}`.  Thus neither a lower curvature
bound `-C n` nor an `O(n)` score bound is available.

Conditionally on `theta`,

```
R/rho(theta) ~ Beta(n,1),
E_theta R = n rho/(n+1),
Var_theta R = n rho^2/((n+2)(n+1)^2).                              (6.4)
```

The conditional mean ranges from asymptotically `sqrt(3)` on coordinate axes
to `sqrt(3n)` on diagonal directions.  Hence pointwise replacement of the
conditional radius by `sqrt(n)` loses order `sqrt(n)` on open angular sets,
even though its integrated loss is controlled by thin shell.

### 6.2 Isotropic regular simplex: radial variance of order one on vertex rays

Work in the `n`-dimensional subspace

```
H_n = { x in R^{n+1} : sum_i x_i = 0 }
```

and let

```
a = sqrt((n+2)/(n+1)),
K = { x in H_n : x_i >= -a for every i }.
```

If `P` is uniform on the probability simplex and

```
X = sqrt((n+1)(n+2)) (P - (1/(n+1)) 1),
```

then `X` is uniform on `K` and `Cov(X)=I_{H_n}`.  Thus `K` is exactly isotropic.

For `theta in S(H_n)`,

```
rho(theta) = a/(-min_i theta_i),
H(theta) = n log(-min_i theta_i) + constant.                       (6.5)
```

On a chamber with unique minimizing coordinate `j`, (6.2) holds verbatim.
There are explicit chamber points with `theta_j^2=2/(n+1)` and a tangent vector
with `u_j^2=(n-2)/(n+1)`, giving again

```
Hess_S H(theta)[u,u] = -n^2/2.                                    (6.6)
```

For completeness, one such point has coordinates

```
theta_j = -alpha,    alpha=sqrt(2/(n+1)),
theta_k = b          for one k != j,
theta_i = c          for the remaining n-1 coordinates,
R0 = sqrt(1-2/n),
b = alpha/n + sqrt((n-1)/n) R0,
c = alpha/n - R0/sqrt(n(n-1)).
```

It has sum zero, norm one, and `theta_j` is its unique minimum.  Projecting
`e_j` onto `H_n intersect theta^perp` gives the stated tangent.

Along the direction of a vertex,

```
rho_vertex = sqrt(n(n+2)).
```

The conditional law is again (6.4), now with this value of `rho`, and hence

```
Var(R | vertex direction) = n^2/(n+1)^2 -> 1.                     (6.7)
```

Thus no pointwise estimate `Var(R|Theta=theta) <= C/n` is possible.  The
conditional mean ranges from order one to order `n`.

### 6.3 Product Laplace: conditional radial variance can be order `n`

Let

```
f(x) = product_i [ (1/sqrt(2)) exp(-sqrt(2)|x_i|) ].
```

This law is isotropic and log-concave.  Direct integration gives

```
w(theta) = constant * ||theta||_1^{-n},
H(theta) = n log ||theta||_1 + constant,                            (6.8)
```

and

```
R | Theta=theta ~ Gamma(shape=n, rate=sqrt(2)||theta||_1).
```

Consequently

```
E_theta R = n/(sqrt(2)||theta||_1),
Var_theta R = n/(2||theta||_1^2).                                  (6.9)
```

On a coordinate axis the conditional variance is `n/2`; on a diagonal it is
`1/2`.  Hence neither the conditional variance nor the conditional radial
Poincare constant is pointwise dimension-free.  Estimate (4.2) succeeds only
after averaging with the angular law.

Inside a fixed orthant, (6.2) again holds with the active linear functional
`sum_i sign(theta_i) theta_i`.  At a diagonal point the Bakry--Emery curvature
equals `-2`; approaching a coordinate axis through the interior of an orthant
produces curvature approaching `-n^2+O(n)`.

### 6.4 Ellipsoids are only a positive baseline

A centered ellipsoid with semiaxes `a_i` has covariance
`diag(a_i^2/(n+2))`.  Exact isotropy forces every `a_i=sqrt(n+2)`, so the body is
the Euclidean ball.  Its angular law is uniform, `H` is constant, and the
spherical gap is `n-1`.  Thus ellipsoids do not support a negative example
after whitening; they merely show that the desired angular mechanism is exact
in the radial case.

The cube, regular simplex, and product Laplace angular laws are themselves
exactly angular-isotropic by their irreducible symmetry:

```
int theta tensor theta dnu = I/n.
```

Thus adding angular second-moment isotropy to Ball-body convexity still does
not restore any pointwise curvature lower bound.

## 7. Status of the polar route

What is proved without KLS:

1. The angular marginal is the cone measure of a genuine convex Ball body.
2. Its gauge has the distributional convexity (2.4).
3. The radial part of every 1-Lipschitz test costs at most a universal constant,
   by (4.2)--(4.3).
4. Ordinary spherical Bakry--Emery, log-concavity along geodesics, uniformly
   bounded angular score, uniform conditional-radius bounds, and naive
   pointwise replacement by `sqrt(n)` all fail on explicit isotropic examples.
5. A round-angular first-moment estimate of scale `1/sqrt(n)` is quantitatively
   equivalent to T3, by (5.1), once the dimension-free thin-shell theorem is
   used.

Therefore a viable continuation must use the **global positive curvature
measure of the convex gauge** (the singular chamber walls are essential in the
cube and simplex), together with angular isotropy, to prove a genuinely new
global mixing inequality.  Any argument using only the absolutely continuous
Hessian inside smooth chambers is ruled out by (6.2)--(6.6).  Proving (5.5), or
even only (5.4), without such a new global mechanism is simply a restatement of
the missing KLS step.
