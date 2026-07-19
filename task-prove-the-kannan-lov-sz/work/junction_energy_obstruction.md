# What signed-distance stability does and does not force

The smooth and singular second variations naturally produce a Dirichlet form
on normal heights.  This note proves the exact low-energy-to-parallel lemma
and gives a sharp expander obstruction to the converse.  The obstruction
satisfies the covariance, long-ray, coherent-packet, and stability scalings;
only Euclidean log-concave realizability is absent.

## 1. An abstract completed-ray form

Let `(Omega,eta)` be a probability space of rays.  On ray `y`, let
`N_y in S^{n-1}`, let the conditional standard deviation be `sigma_y`,
and let `b_y=q_y(0)` be its boundary density.  Suppose, for fixed constants
`0<c_0<=C_0`,

\[
 \sigma_y\ge s,
 \qquad {c_0\over s}\le b_y\le {C_0\over s}              \tag{1}
\]

almost everywhere, and suppose isotropy gives

\[
 \int\sigma_y^2N_yN_y^T d\eta(y)\preceq I.               \tag{2}
\]

Let `E(h)` denote the complete second-variation cost: the smooth ray-metric
term from signed-distance deformation plus every medial/focal junction term.
For an extremal potential, the stability calculation requires

\[
 \int b_y(h-\eta h)^2d\eta\le E(h)                       \tag{3}
\]

for all admissible heights.

Define the total normal energy

\[
 E_N:=\sum_{k=1}^n E(\langle e_k,N\rangle).               \tag{4}
\]

**Lemma 1 (low energy forces the bounded-scale branch).**  Under
(1)--(3),

\[
 E_N\ge {c_0\over s}\left(1-{1\over s^2}\right).         \tag{5}
\]

In particular, if `s>=sqrt(2)`, then `E_N>=c_0/(2s)`.

**Proof.**  Put `m=int N d eta`.  From (2) and Jensen,

\[
 |m|^2\le\left\|\int N N^T d\eta\right\|_{op}
 \le {1\over s^2}.                                       \tag{6}
\]

Sum (3) for the `n` heights `<e_k,N>`.  Equations (1) and (6) give

\[
 E_N\ge\int b_y|N_y-m|^2d\eta
 \ge {c_0\over s}\int|N_y-m|^2d\eta
 ={c_0\over s}(1-|m|^2),                                 \tag{7}
\]

which is (5).  QED.

Thus a completed normal congruence whose smooth-plus-junction normal energy
is `o(1/s)` cannot support long rays.  More quantitatively, without using
(2), if

\[
 E_N\le {c_0\varepsilon\over s},                          \tag{8}
\]

then (7) gives `int|N-m|^2<=epsilon`; hence
`|m|>=sqrt(1-epsilon)` and, with `v=m/|m|`,

\[
 \int|N-v|^2d\eta=2(1-|m|)\le2\varepsilon.               \tag{9}
\]

This is the rigorous approximate-parallel branch.  Combining (9) with (2)
immediately bounds `s` by a universal constant when `epsilon` is a fixed
small constant.

## 2. A saturated expander obstruction

Lemma 1 is sharp at exactly the scale supplied by the one-dimensional trace
identity.  Fix a large `s` and an integer

\[
 M\ge e^{c s}.                                            \tag{10}
\]

Let `Omega={1,...,M}`, `eta_i=1/M`, and take orthonormal normal
directions `N_i=e_i in R^M`.  Set

\[
 b_i={c_0\over s}.                                       \tag{11}
\]

On the complete graph define

\[
 E(h)={c_0\over sM^2}\sum_{1\le i<j\le M}(h_i-h_j)^2.   \tag{12}
\]

The elementary identity

\[
 \sum_{i<j}(h_i-h_j)^2=M\sum_i(h_i-\bar h)^2             \tag{13}
\]

shows that (3) is an equality for every height.  Also

\[
 E_N={c_0(M-1)\over sM}\asymp {1\over s},                \tag{14}
\]

so (5) is saturated.  If every conditional variance is `s^2`, then

\[
 \int\sigma_i^2N_iN_i^Td\eta(i)={s^2\over M}I_M\preceq I_M \tag{15}
\]

for `M>=s^2`.  The coherent-cap mass is `1/M<=e^{-cs}`, exactly
consistent with the logarithmic barycenter constraint.

The geometry can be chosen nonparallel and nonconcurrent at the level of
calibrated components.  For example, use the cyclic endpoint pairing

\[
 p_i=L e_i,\qquad c_i=-L e_{i+1},                         \tag{16}
\]

whose calibrated directions are `(e_i+e_{i+1})/sqrt(2)` and whose
lines have no common point.  Alternatively, use separated flat cylinders
with axes `e_i` and generic bases.  The smooth shape energy is zero; (12)
is an abstract medial-junction expander energy.

This construction is not a globally log-concave Euclidean measure.  That is
the point: every inequality currently obtained from raywise log-concavity,
isotropy, the exponential coherent-packet bound, Hilbert--Schmidt curvature,
and extremal second variation is compatible with (10)--(16).  A constant-gap
junction graph absorbs the boundary gain without forcing either alignment or
concurrence.

## 3. Sharp remaining geometric statement

The signed-distance mechanism therefore reduces the inverse problem to the
following genuinely geometric assertion.

> **Euclidean gluing exclusion.**  A single globally defined signed-distance
> potential for a globally log-concave density cannot realize, on a
> positive-mass family of balanced rays of scale `s>>1`, a completed medial
> graph whose node masses are at most `e^{-cs}` and whose junction Dirichlet
> form has a universal spectral gap at total weight `Theta(1/s)`, unless the
> normal lines are approximately concurrent.

This assertion is not a consequence of the abstract graph stability
inequality: the complete-graph model proves that.  Nor does indicator
midpoint multiplicity prove it, because that multiplicity is null-set
unstable and does not weight actual reflected representations.  A proof must
relate Euclidean junction conductances to positive-density log-concave bridges
and to the *same-ray* endpoint matching.  That relation is the precise new
content still absent from the parallel/concurrent inverse theorem.
