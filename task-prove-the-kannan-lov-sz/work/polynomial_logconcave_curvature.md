# A polynomial-times-log-concave curvature lemma

This note proves the one-dimensional estimate needed on a smooth normal
chart.  It allows mixed principal-curvature signs, nonsmooth convex
potentials, and focal factors which vanish at endpoints of the support.

## The lemma

**Lemma.**  Fix `delta in (0,1/2]`.  Let `J` be an interval containing
zero, let `phi:J->R union {+infinity}` be proper lower-semicontinuous and
convex, and let `kappa_1,...,kappa_d` be real numbers such that

\[
 1+t\kappa_i>0\quad\text{in the interior of }J             \tag{1}
\]

for every `i`.  Endpoint equality in (1) is allowed.  Suppose

\[
 q(t)=Z^{-1}e^{-\phi(t)}\prod_{i=1}^d(1+t\kappa_i)
        {\bf1}_J(t)                                       \tag{2}
\]

is a probability density with finite, nonzero variance `sigma^2` and

\[
 \mathbb P(T\le0)\ge\delta,
 \qquad \mathbb P(T\ge0)\ge\delta.                       \tag{3}
\]

Then

\[
 \boxed{\qquad
   \sigma^2\sum_{i=1}^d\kappa_i^2\le C_\delta,
 \qquad}                                                  \tag{4}
\]

where, for example, one may take

\[
 C_\delta={400\over\delta(1-\delta)}
       \left(1+{40\over\delta^{3/2}}\right)^2.            \tag{5}
\]

No dependence on `d` occurs.

### Preliminary density bound

We use the following elementary estimate: if `h` is a one-dimensional
log-concave probability density with variance `s^2`, then

\[
 \|h\|_\infty\le {10\over s}.                             \tag{6}
\]

Here is a proof with deliberately loose constants.  Translate a mode to
zero and write `M=h(0)`.  On each side of zero, let `a` be the first point
at which `h` is at most `M/e`, using the support endpoint if the density
does not reach that value before the endpoint.  Log-concavity gives
`a<=e/M`, since `h>=M/e` before `a`.  Beyond a genuine crossing,

\[
 h(x)\le M e^{-x/a};                                      \tag{7}
\]

if `a` is a support endpoint there is no tail.  Consequently the second
moment about the mode on either side is at most

\[
 M a^3\left({1\over3}+\int_1^\infty u^2e^{-u}\,du\right)
 < 47M^{-2}.                                              \tag{8}
\]

Both sides contribute less than `94M^{-2}`.  Variance is no larger than
the second moment about the mode, so `s^2<94M^{-2}`, proving (6).

### Proof of the lemma

Put

\[
 W(t)=\phi(t)-\sum_{i=1}^d\log(1+t\kappa_i)+\log Z.       \tag{9}
\]

It is convex on the interior of the support and `q=e^{-W}` there.  In the
sense of distributional second derivatives,

\[
 W''\ge \sum_{i=1}^d {\kappa_i^2\over(1+t\kappa_i)^2}\,dt. \tag{10}
\]

Formula (10) remains valid when `phi` is nonsmooth and when a factor
vanishes at a support endpoint.

Let `m=ET`.  Cantelli's inequality and (3) imply

\[
 |m|\le {\sigma\over\sqrt\delta}.                         \tag{11}
\]

Indeed, if `m>0`, then
`delta<=P(T-m<=-m)<=sigma^2/(sigma^2+m^2)`; the other sign is
identical.

Let `t_-` and `t_+` be respectively a `delta/2` and a
`1-delta/2` quantile.  Chebyshev's inequality together with (11) gives

\[
 |t_-|,|t_+|\le {3\sigma\over\sqrt\delta}.                \tag{12}
\]

The interval `[t_-,t_+]` has probability at least `1-delta`.  Applying
(6) to `q` therefore gives

\[
 t_+-t_-\ge {(1-\delta)\sigma\over10}.                    \tag{13}
\]

We next keep every focal denominator uniformly bounded above on this
quantile interval.  If `kappa_i>0`, its zero lies on the negative side.
The negative support length is at least `delta/||q||_infty`, by (3), and
the zero cannot lie in the interior of the support.  Thus

\[
 {1\over\kappa_i}\ge {\delta\sigma\over10}.               \tag{14}
\]

For `kappa_i<0`, the same argument on the positive side gives (14) with
`|kappa_i|`.  Hence

\[
 |\kappa_i|\le {10\over\delta\sigma}.                    \tag{15}
\]

Combining (12) and (15), throughout `[t_-,t_+]`,

\[
 0<1+t\kappa_i\le
 B_\delta:=1+{30\over\delta^{3/2}}.                       \tag{16}
\]

It remains to bound the total increase of the convex slope of `W`.  If
`W'_+(t_+)>0`, convexity gives

\[
 {\delta\over2}\le\int_{t_+}^\infty q(t)dt
 \le {q(t_+)\over W'_+(t_+)},                             \tag{17}
\]

and (6) yields `W'_+(t_+)<=20/(delta sigma)`.  If that
derivative is nonpositive, the same upper bound is automatic.  Similarly,

\[
 W'_-(t_-)\ge-{20\over\delta\sigma}.                      \tag{18}
\]

With the usual one-sided derivative convention for a convex function,
(17)--(18) imply

\[
 W''((t_-,t_+))\le {40\over\delta\sigma}.                 \tag{19}
\]

Integrating (10), and using (13) and (16), now gives

\[
 { (1-\delta)\sigma\over10B_\delta^2}
       \sum_i\kappa_i^2
 \le {40\over\delta\sigma}.                             \tag{20}
\]

This proves (4), with a constant slightly smaller than (5); (5) absorbs
all endpoint and quantile conventions.  QED.

## Checks and consequences

1. **Repeated spherical factor.**  If all `kappa_i=1/R`, (4) gives

   \[
    \sigma\le \sqrt{C_\delta}\,{R\over\sqrt d}.          \tag{21}
   \]

   Thus the model density proportional to
   `exp(-phi(t))(1+t/R)^d` cannot have scale comparable to `R` when
   `d` is large.

2. **Mixed signs.**  No cancellation is possible: every factor contributes
   the positive term `kappa_i^2/(1+t kappa_i)^2` to (10).  The two-sided
   quantile hypothesis separately keeps positive and negative focal roots
   out to distance `c_delta sigma`.

3. **Singular focal endpoints.**  A factor may vanish at an endpoint of
   `J`; the proof only integrates over the interior quantile interval, where
   it is positive.  The blow-up of `-log(1+t kappa_i)` at the endpoint can
   only increase the convex second-derivative measure and hence helps (19).

4. **Normal charts.**  For
   `q_z(t) proportional to exp(-V(z+tN)) det(I+tS_z)`, take
   `phi(t)=V(z+tN)`, whose convexity follows from that of `V`, and take the
   eigenvalues of `S_z` as the `kappa_i`.  If zero is between fixed
   `delta` and `1-delta` quantiles, (4) gives

   \[
    \sigma_z^2\|S_z\|_{HS}^2\le C_\delta.                \tag{22}
   \]

The estimate is entirely local along a ray.  It proves that a long balanced
conditional can change orientation only across base distances of order its
scale, and that spherical curvature in many directions shortens the
conditional by the square root of the curvature rank.  It does not by itself
control orientation jumps across focal or medial endpoint sets; that remains
the global compatibility issue.
