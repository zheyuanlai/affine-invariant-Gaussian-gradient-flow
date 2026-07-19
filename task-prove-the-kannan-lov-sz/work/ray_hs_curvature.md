# Long balanced rays force small full shape curvature

This note proves a one-dimensional lemma for the Jacobian factors in a smooth
normal-ray chart.  Unlike the elementary reach estimate, it controls the full
Hilbert--Schmidt norm of the shape operator.

## Proposition 1 (polynomial-Jacobian curvature bound)

Let `kappa_1,...,kappa_d` be real numbers and let

\[
 I=\{t\in\mathbb R:1+\kappa_i t>0\text{ for every }i\}.
\]

Let `phi:I->R` be convex and suppose that

\[
 q(t)={1\over Z}e^{-\phi(t)}\prod_{i=1}^d(1+\kappa_i t)
                                                                  \tag{1}
\]

is a probability density.  If, for some `delta in (0,1/2]`,

\[
 \delta\le\int_{I\cap(-\infty,0]}q(t)dt\le1-\delta,               \tag{2}
\]

and `sigma^2=Var_q(T)`, then

\[
 \boxed{\quad
 \sigma^2\sum_{i=1}^d\kappa_i^2\le C_\delta .\quad}              \tag{3}
\]

The constant is independent of `d`, the interval, the curvatures, and
`phi`.

### Proof

Put

\[
 K=\sum_i\kappa_i^2.
\]

There is nothing to prove if `K=0`.  Choose a subgradient `a` of `phi` at
zero and write

\[
 \psi(t)=\phi(t)-\phi(0)-at\ge0,
 \qquad
 h(t)=e^{-at}\prod_i(1+\kappa_i t).
\]

Then `h(0)=1` and

\[
 {q(t)\over q(0)}=e^{-\psi(t)}h(t)\le h(t).                       \tag{4}
\]

On `I`, set `ell=log h`.  Its first two derivatives are

\[
 \ell'(t)=-a+\sum_i{\kappa_i\over1+\kappa_i t},
 \qquad
 \ell''(t)=-\sum_i{\kappa_i^2\over(1+\kappa_i t)^2}.              \tag{5}
\]

Let `r=(2sqrt(K))^{-1}`.  Since
`max_i|kappa_i|<=sqrt(K)`, the whole interval `[-r,r]` is contained
in `I`, and on it

\[
 -\ell''(t)\ge {4\over9}K.                                      \tag{6}
\]

Write `b=ell'(0)`.  If `b<=0`, (6) gives, for `0<=t<=r`,

\[
 \ell(t)\le-{2\over9}Kt^2,
 \qquad
 \ell'(r)\le-{2\over9}\sqrt K.                                 \tag{7}
\]

Since `ell` is concave, after `r` it lies below its tangent there.
Consequently

\[
 \int_{I\cap[0,\infty)}h(t)dt\le {C\over\sqrt K}.                \tag{8}
\]

If `b>=0`, the identical argument on the negative half-line gives

\[
 \int_{I\cap(-\infty,0]}h(t)dt\le {C\over\sqrt K}.               \tag{9}
\]

In the first case, (2), (4), and (8) imply

\[
 1-\delta\ge\int_{t>0}q(t)dt
 \quad\hbox{and, more importantly,}\quad
 \delta\le\int_{t>0}q(t)dt
       \le q(0){C\over\sqrt K},                                 \tag{10}
\]

where the lower bound by `delta` follows from (2).  Thus
`q(0)>=c delta sqrt(K)`.  In the second case the same conclusion follows
from the left mass and (9).  Hence in all cases

\[
 q(0)\ge c\delta\sqrt K.                                        \tag{11}
\]

For every one-dimensional log-concave probability density of variance
`sigma^2`,

\[
 \|q\|_\infty\le {C_0\over\sigma}.                              \tag{12}
\]

For completeness, (12) follows by translating a mode to zero and using the
two supporting exponential tails of `log q`: normalization fixes their
decay lengths to be at most universal multiples of
`1/||q||_infinity`, and integration of the two tail second moments gives
`Var(T)<=C_0^2/||q||_infinity^2`.  Since `q(0)<=||q||_infinity`,
(11)--(12) give

\[
 \sigma\sqrt K\le C_0/(c\delta),
\]

which is (3).

## 2. Geometric interpretation

Let `Sigma` be a smooth oriented hypersurface, let `z in Sigma`, let
`N(z)` be its unit normal, and let `kappa_i(z)` be its principal
curvatures with the convention

\[
 \det D(z,t)\mapsto z+tN(z)=\prod_i(1+t\kappa_i(z)).               \tag{13}
\]

Suppose a log-concave density `e^{-V}` is disintegrated in this normal
chart before its first focal endpoints.  The conditional density on the ray
has exactly the form

\[
 q_z(t)\propto e^{-V(z+tN(z))}
                 \prod_i(1+t\kappa_i(z)),                         \tag{14}
\]

and `t->V(z+tN(z))` is convex.  If the separator `t=0` has conditional
mass on each side at least `delta`, Proposition 1 gives

\[
 \boxed{\quad
 \operatorname {Var}_{q_z}(T)\,|\mathrm {II}(z)|_{HS}^2
 \le C_\delta .\quad}                                           \tag{15}
\]

For a sphere of radius `R` in `d+1` dimensions, (15) has the correct scale:
`|II|_{HS}^2=d/R^2`, and the Jacobian alone confines the normal coordinate
to scale `R/sqrt(d)`.  Thus a family of rays with conditional scale `s`
forces full curvature `|II|_{HS}=O(1/s)`, not merely
`||II||_op=O(1/s)`.

## 3. Remaining audit issues

Equation (15) is pointwise on a smooth normal chart.  A global argument still
has to address:

1. ray endpoints at which one or more Jacobian factors vanish;
2. polyhedral or medial junctions, where curvature is a singular full-BV
   turning measure rather than an almost-everywhere shape operator; and
3. conversion of small weighted full curvature on the long-ray quotient into
   a parallel or concurrent structure without invoking a hidden Poincare
   inequality on that quotient.

The Gaussian fan shows that item 2 is essential.  Proposition 1 nevertheless
removes the earlier `sqrt(d)` loss on every smooth chart and supplies the
correct load-bearing local estimate for a global focal/turning argument.

