# Critical convex charts and the max-cell gluing problem

This note continues `polar_angular_report.md`.  Put `d=n-1`.  If `G` is the
gauge of the Ball body and

```
  dnu(theta)=Z^{-1}G(theta)^{-n} dsigma(theta),
```

then on the signed max-coordinate cell `C_a`, `a in {+/-1,...,+/-n}`,
the gnomonic coordinate `theta=(e_a+y)/sqrt(1+|y|^2)` identifies the cell
with `Q=[-1,1]^d` and

```
  dnu = c_G phi_a(y)^(-n)dy,
  phi_a(y)=G(e_a+y).
```

The function `phi_a` is convex.  For a locally Lipschitz `h` on the sphere,
write `h_a(y)=h(theta(y))`.  The exact Dirichlet form is

```
  E_nu(h,h)
   = sum_a c_G int_Q (1+|y|^2)
       ( |grad h_a|^2+(y dot grad h_a)^2 ) phi_a(y)^(-n)dy.       (1)
```

No constants are suppressed in (1); walls have zero `nu`-mass and hence are
counted only once after choosing any measurable convention.

## 1. Angular covariance is automatically well conditioned

The following observation removes a possible escape from the gluing problem.

**Lemma 1 (angular second moments).**  There are universal `0<c<C<infty`
such that, for every dimension `n>=2`, every isotropic log-concave random
vector `X`, `R=|X|`, `Theta=X/R`, and every unit vector `u`,

```
               c/n <= E <Theta,u>^2 <= C/n.                       (2)
```

**Proof.**  Write `Z=<X,u>` and `q=E[Z^2/R^2]`.  The standard one-dimensional
moment estimate for isotropic log-concave laws gives `E Z^4<=C_4`.  Paouris'
negative-moment estimate, used only with exponent four, gives
`E R^(-4)<=C_- n^(-2)` for `n>=n_0`; enlarging the constants handles the
finitely many smaller dimensions.  Cauchy--Schwarz therefore gives

```
  q <= (E Z^4)^(1/2)(E R^(-4))^(1/2) <= C/n.                       (3)
```

For the reverse inequality, apply Cauchy--Schwarz to
`Z^2=sqrt((Z^2/R^2)(Z^2R^2))`:

```
  1=(E Z^2)^2 <= q E[Z^2R^2]
    <= q (E Z^4)^(1/2)(E R^4)^(1/2).                              (4)
```

Paouris' positive fourth-moment bound gives `E R^4<=C_+ n^2`.
Equations (4) and the one-dimensional fourth-moment bound give `q>=c/n`.
This proves (2).  Notice that no thin-shell estimate is used.  QED.

Thus the angular covariance matrix `Q_nu=int theta tensor theta dnu` obeys

```
                    (c/n)I <= Q_nu <= (C/n)I.                     (5)
```

Exact equality holds in the cube, regular-simplex and product-Laplace stress
tests by irreducible symmetry.  Consequently a gluing theorem formulated
under the stable angular-isotropy condition (5) is aimed at the full class
rather than a special subclass.  Condition (5), unlike exact equality, is
inherited directly and requires no projective change of coordinates.

## 2. Exact wall measure

Let `F_{ab}` be the common wall of two compatible signed cells (opposite
signed copies of the same coordinate are not adjacent).  In the chart of
`C_a`, choose coordinates so that `F_{ab}={y_k=epsilon}`.  Write
`y=(z,epsilon)`.  Since

```
  g_y^(-1)=(1+|y|^2)(I+y tensor y),
```

the coarea factor of the coordinate `y_k` is

```
  |grad_S y_k|=sqrt((g_y^(-1))_kk)
              =sqrt((1+|y|^2)(1+y_k^2)).                          (6)
```

It follows that the weighted spherical `(n-2)`-measure of the interface is

```
  b_ab = c_G int_[-1,1]^(d-1)
           sqrt(2(2+|z|^2)) phi_a(z,epsilon)^(-n) dz.              (7)
```

This is invariant under using the chart from the other side.  Formula (7) is
the correct conductance for every BV gluing argument.  Omitting the coarea
factor loses `sqrt(n)` on the basic examples.

### Cube calibration

For the isotropic cube, `phi_a=1`, every cell has mass `p_a=1/(2n)`, and

```
  b_ab = (1/(4n)) E sqrt(2(2+|Z|^2)),                              (8)
```

where `Z` is uniform on `[-1,1]^(n-2)`.  Since
`E|Z|^2=(n-2)/3` and `E|Z|^4<=C n^2`, Paley--Zygmund, followed by the trivial
upper bound, gives universal `c,C` with

```
                 c/sqrt(n) <= b_ab <= C/sqrt(n).                  (9)
```

In particular `b_ab/p_a` is of order `sqrt(n)`.  The often quoted
`1/sqrt(n)` interface-to-cell ratio has the coarea factor backwards.

### Product-Laplace calibration

For the isotropic product Laplace law, write

```
  Theta = (epsilon_i E_i)_i / (sum_i E_i^2)^(1/2),                (10)
```

where the signs are independent uniform signs and the `E_i` are independent
unit exponentials.  This is exactly its angular law.  Across the wall where
the `i`th sign changes,

```
 lim_(eps downarrow 0) eps^(-1) P(0<=Theta_i<=eps | fixed signs)
       = E (sum_(j != i) E_j^2)^(1/2).                            (11)
```

The last expectation is between `c sqrt(n)` and `C sqrt(n)` by the second and
fourth moments of an exponential variable.  Hence each orthant wall has
conductance comparable to `2^(-n)sqrt(n)`, again a `sqrt(n)`
interface-to-cell ratio.  This calculation includes the singular Hessian
charge of the `ell_1` gauge at the coordinate wall.

### Regular-simplex calibration

Use the exactly isotropic realization from `polar_angular_report.md`,

```
 K={x in H_n:x_i>=-a},   a^2=(n+2)/(n+1),
 H_n={sum_i x_i=0}.
```

Let `F_i={x_i=-a}` and `R_ij=F_i intersect F_j`.  The cell `C_i` is the cone
over `F_i`, and symmetry gives `p_i=1/(n+1)`.  Put
`d_F=dist(0,aff F_i)` and `d_R=dist(0,aff R_ij)`.  Orthogonal projection in
`H_n` gives

```
            d_F^2=(n+2)/n,       d_R^2=2(n+2)/(n-1).              (12)
```

If `A_F` and `A_R` are the Euclidean facet and ridge areas and `Y` is uniform
on `R_ij`, radial area change on the affine ridge gives

```
  b_ij/p_i = (d_R A_R)/(d_F A_F) E|Y|.                            (13)
```

Indeed `d_R A_R=int_Gamma rho^(n-1)dsigma`, while the extra round coarea
factor changes this to `int_Gamma rho^n dsigma=d_R A_R E|Y|`; division by
the corresponding cone-volume formula for `F_i` proves (13).

The simplex side length is `sqrt(2(n+1)(n+2))`.  The regular-simplex area
formula therefore yields

```
 (d_R A_R)/(d_F A_F)
     = sqrt(2)(n-1)/sqrt((n+1)(n+2)) asymp 1.                      (14)
```

Writing `Y=sum_(k != i,j) P_k v_k`, with
`P~Dirichlet(1,...,1)`, gives

```
       E|Y|^2=(n+2)^2/n,             E|Y|^4<=C n^2.               (15)
```

The first identity follows from
`|v_k|^2=n(n+2)`, `<v_k,v_l>=-(n+2)` for `k!=l`, and the second Dirichlet
moments.  The fourth bound follows from the fourth Dirichlet moments and the
same Gram matrix.  Cauchy--Schwarz gives `E|Y|<=C sqrt(n)`;
Paley--Zygmund applied to `|Y|^2` gives the reverse bound.  Combining
(13)--(15) proves

```
                         b_ij/p_i asymp sqrt(n).                   (16)
```

## 3. A precise counterexample to chartwise coercivity

Convexity at the critical exponent does not give the required gap on one
max-coordinate chart.  This remains false even when the chart function is the
restriction of a genuine global gauge.

Fix `n>=3`, put `m=n-2`, and on `Q=[-1,1]^(n-1)` write `y=(t,z)` with
`z in R^m`.  For `0<eps<1` and `a>1`, set

```
       phi_(a,eps)(t,z)=1+eps|t|+a||z||_1.                         (17)
```

This is the restriction to a signed max cell of the norm

```
 G_(a,eps)(x)=|x_1|+eps|x_2|+a sum_(j=3)^n |x_j|.                 (18)
```

Let `lambda_(a,eps)` be the probability law on `Q` proportional to
`phi_(a,eps)^(-n)`.  For `c>0`, scaling `u=az` and monotone convergence give

```
 a^m int_[-1,1]^m (c+a||z||_1)^(-(m+2))dz
   -> int_R^m (c+||u||_1)^(-(m+2))du
    = [2^m/(m+1)!]c^(-2).                                        (19)
```

Consequently, as `a->infty`, the `t`-marginal converges in total variation
to the law on `[-1,1]` with density proportional to
`(1+eps|t|)^(-2)`.  Also `z->0` in probability and

```
 E_lambda |z|^2 = O((log a)/a^2).                                 (20)
```

The last estimate follows by the same scaling: the radial `ell_1` integral
for the second moment grows only logarithmically at its upper cutoff.

For the test function `h(t,z)=t`, the exact round form (1) is

```
 E_lambda[(1+t^2+|z|^2)(1+t^2)].                                 (21)
```

First let `a->infty` and then `eps->0`.  Equations (19)--(21) give

```
 Var_lambda(h)->1/3,
 E_lambda[(1+t^2+|z|^2)(1+t^2)]->E_Unif[-1,1](1+t^2)^2<=4.        (22)
```

Hence the local weighted spectral gap is at most `12+o(1)`, rather than
`c n`.  The counterexample is compatible with every local hypothesis
`phi>0`, `phi` convex, and exponent `n=d+1`.

It is excluded only by the global angular condition (5).  Indeed the gauge
(18) is the angular gauge of a product of one-dimensional exponentials with
coordinate scales proportional to `1`, `eps^(-1)`, and `a^(-1)`; as
`eps->0` its angular covariance acquires an eigenvalue tending to one in the
second coordinate.  Thus the counterexample pinpoints the necessary missing
input: angular conditioning must couple different cells.  It cannot be used
after a separate chartwise Poincare estimate, because that estimate is false.

There is also a universality reason that no covariance-normalized local
theorem is an easier substitute.  Let `L` be any convex body of positive
volume contained in the interior of `Q`, and define

```
                    phi_A(y)=1+A dist(y,L).                        (23)
```

Then `phi_A` is positive and convex.  Dominated convergence gives

```
   phi_A^(-(d+1))dy / int_Q phi_A^(-(d+1))dy
                   -> 1_L dy/|L|                                  (24)
```

in total variation as `A->infty`.  For every fixed Lipschitz test function,
both its variance and its round-chart energy converge to the corresponding
quantities under the uniform law on `L`.

Moreover (23) is not an artificial chart function.  Any positive convex
function `phi` on `Q` which is the restriction of a finite globally Lipschitz
convex function has a gauge extension.  Write

```
       phi(y)=sup_(alpha,beta) {alpha+beta dot y}
```

over its affine supporting minorants and put

```
 P(t,x)=sup_(alpha,beta){alpha t+beta dot x},
 G(t,x)=max{P(t,x),||x||_infty,-t}.                                (25)
```

The bounded subgradients make `P` finite.  The function `G` is finite,
convex, positively homogeneous and positive away from the origin.  On the
max cell `t>=||x||_infty`, one has

```
             G(t,x)=P(t,x)=t phi(x/t),                            (26)
```

because `phi>=1`.  Thus it is the gauge of a genuine convex body and restricts
to `phi` in that chart.  The distance function in (23) is globally
Lipschitz, so the construction applies.

It follows from (24)--(26) that the class of individual critical convex
charts contains uniform measures on arbitrary convex bodies as limits.  A
local Poincare theorem normalized by the chart covariance would therefore
contain convex-body KLS itself: translate and dilate an arbitrary convex body
into `B(0,1/2) subset Q`, where the inverse round metric in (1) lies between
`I` and `(5/4)(I+yy^T)<= (25/16)I`; pass to the limit in (24), and undo the
dilation.  The stronger `c n` local gap is already disproved by (17)--(22).

## 4. Why the remaining gluing estimate is load-bearing

Define the global critical-chart assertion

```
  (CG2)  Var_nu(h) <= (C/n) E_nu(h,h)                              (27)
```

for every Ball gauge paired with an isotropic log-concave law.  By (1), this
is exactly the requested weighted Poincare inequality after all max cells are
glued; it is not a comparison with a different quadratic form.

Let `F:R^n->R` be one-Lipschitz and set
`h(theta)=F(sqrt(n)theta)`.  Then `|grad_S h|<=sqrt(n)`.  Assuming (27),

```
  E_nu |h-E_nu h| <= sqrt(Var_nu h) <= sqrt(C).                    (28)
```

The dimension-free thin-shell theorem and
`|F(X)-F(sqrt(n)Theta)|<=|R-sqrt(n)|` now imply

```
                  E|F(X)-EF(X)| <= C'.                            (29)
```

Milman's first-moment equivalence for log-concave measures turns (29) into
KLS.  Lemma 1 verifies that the apparently additional angular-isotropy input
is already available universally.  Therefore (CG2) is a complete KLS-strength
statement, not an intermediate estimate.

The exact first-moment version

```
  (CG1)  sup_(Lip_S h<=1) E_nu|h-E_nu h| <= C/sqrt(n)              (30)
```

is quantitatively equivalent to KLS for the paired law: Proposition 5.1 of
`polar_angular_report.md` proves both implications, with an additive
`2(E(R-sqrt(n))^2)^(1/2)` which is universal by thin shell.  Thus any gluing
argument whose final step is (30) has merely renamed target (T3).

The wall calibrations (9)--(11) do rule out a *spurious* chartwise loss: on
the three mandatory polyhedral tests, the missing factor `sqrt(n)` is present
exactly in the round coarea factor.  They do not prove a global graph
inequality.  A cut need not be a union of max cells; it may cross every cell,
and its two one-sided BV traces on a tie wall agree unless the cut itself has
a jump there.  Consequently the fixed adjacency graph and the numbers
`(p_a,b_ab)` do not determine its perimeter.  A valid gluing theorem must
control simultaneously

```
  (i)  within-cell BV variation,
  (ii) trace deviation from cell representatives, and
  (iii) jumps of those representatives across walls.                              (31)
```

The required coercivity constant in (31) is precisely the Cheeger constant of
`nu`: summing the distributional BV derivative over the partition gives the
identity

```
 Per_nu(A)
  = sum_a Per_nu(A; interior C_a)
    + sum_(a<b) int_F_ab |Tr_a 1_A-Tr_b 1_A| dnu_F.                (32)
```

Here (32) is the standard BV partition formula; `dnu_F` is exactly (7).
Taking the infimum of the right-hand side divided by
`min(nu(A),1-nu(A))` is therefore not a lower bound for the angular Cheeger
constant but its definition written in cell coordinates.  Any proposed
"formal gluing lemma" must supply an estimate beyond (32); if that estimate
has constant `c sqrt(n)`, then (28)--(29) prove KLS.

## 5. Curvature-measure separation is not stable

Two complementary examples prevent assigning the gluing gain to only one
part of the distributional Hessian of `G`.

* For `G=||.||_infty`, `D^2G` vanishes in every open max cell and all of its
  positive measure is carried by tie walls.
* For `G=||.||_2`, `D^2G` is absolutely continuous away from the origin and
  gives zero mass to every artificial max-coordinate tie wall.

Both gauges are angularly isotropic and both have the target spherical scale.
Moreover the hyperoctahedrally symmetric interpolation

```
             G_eps=(1-eps)||.||_2+eps||.||_infty                  (33)
```

is angularly isotropic for every `eps`; its singular wall charge is exactly
`eps` times that of the cube and can tend to zero without changing the
dimension or the symmetry condition.  Thus a uniform estimate based on a
positive lower bound for singular wall curvature is false, while an estimate
using only the absolutely continuous Hessian fails at `eps=1`.  Any surviving
curvature functional must use the full distributional convexity measure and
must prove a coercivity theorem stable under (33).  By Section 4, such a
coercivity theorem at scale `sqrt(n)` (Cheeger) or `n` (Poincare) is already
KLS-strength.

## 6. Status

The projective calculation corrects the apparent chartwise `n`-loss on the
cube, simplex and product-Laplace models: their interfaces have the right
`sqrt(n)` conductance once the inverse round metric is used.  What remains is
not a finite-state gluing problem.  It is the continuum trace coercivity in
(31), and its sharp constant is the global angular Cheeger/Poincare constant.
The first-moment form is exactly equivalent to KLS, while the Poincare form
implies KLS and may be strictly stronger because the angular measure is not
log-concave.  No dimension-free gluing lemma follows from the mass vector,
adjacency graph, or the singular wall charge alone.
