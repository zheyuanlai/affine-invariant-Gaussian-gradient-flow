# Coarea cascade for a smoothed Cheeger near-minimizer

## 1. Exact integrated deficit

Let `mu` be a log-concave probability on its affine hull `E`, let `psi`
be its Cheeger constant, and let `P_mu` denote relaxed weighted `BV`
perimeter on `E`.  The exterior-Minkowski and relaxed-perimeter definitions
of `psi` agree by the distance-cutoff relaxation stated and audited in
Section 12.1 of `heatflow_bernstein.md`.  Let `S` have mass `1/2`, and let
`F` be a `[0,1]`-valued locally Lipschitz function with finite weighted
total variation and

\[
                         E F={1\over2}.
\]

Put

\[
 U={1\over2}E|F-\mathbf1_S|,
 \qquad A_r=\{F>r\},\quad 0<r<1.                    \tag{1.1}
\]

Because `EF=mu(S)`, the two one-sided errors are exactly equal:

\[
 \int_S(1-F)d\mu=\int_{S^c}F d\mu=U.               \tag{1.2}
\]

The weighted `BV` coarea theorem and layer cake give the two exact identities

\[
 \int|\nabla F|d\mu=\int_0^1P_\mu(A_r)dr,           \tag{1.3}
\]

\[
 \inf_c\int|F-c|d\mu
 =\int_0^1\min\{\mu(A_r),1-\mu(A_r)\}dr.             \tag{1.4}
\]

To verify (1.4), choose a median `m` of `F`.  For `r<m` the smaller
side is `{F<=r}`, while for `r>m` it is `{F>r}`; integrating the two
tails gives `E|F-m|`.  Atom conventions do not affect the integral.

Consequently the `L^1` Cheeger deficit has the exact level-set
decomposition

\[
\boxed{
 \int|\nabla F|d\mu-
 \psi\inf_c\int|F-c|d\mu
 =\int_0^1\left[P_\mu(A_r)-
 \psi\min\{\mu(A_r),1-\mu(A_r)\}\right]dr.}          \tag{1.5}
\]

There is a second exact cascade:

\[
\boxed{
 \int_0^1\mu(A_r\mathbin{\triangle} S)dr=2U.}        \tag{1.6}
\]

Indeed the integrand identity holds pointwise after integrating in `r`.

Suppose now that `P_mu(S)=p<=psi/2+varepsilon` and

\[
                         \int|\nabla F|d\mu\le p.    \tag{1.7}
\]

The triangle inequality gives

\[
 \inf_c\int|F-c|d\mu
 \ge {1\over2}-E|F-\mathbf1_S|
 ={1\over2}-2U.                                      \tag{1.8}
\]

Equations (1.5), (1.7), and (1.8) therefore imply

\[
\boxed{
 \int_0^1\left[P_\mu(A_r)-
 \psi\min\{\mu(A_r),1-\mu(A_r)\}\right]dr
 \le\varepsilon+2\psi U.}                           \tag{1.9}
\]

Thus most thresholds are simultaneously close to `S` and close to the
global Cheeger ratio whenever `U` and `varepsilon/psi` are small.  More
explicitly, for positive `lambda_1,lambda_2`, the exceptional set of
thresholds on which either

\[
 \mu(A_r\mathbin{\triangle} S)>\lambda_1
 \quad\hbox{or}\quad
 P_\mu(A_r)-\psi\min(\mu(A_r),1-\mu(A_r))>lambda_2
\]

has Lebesgue measure at most

\[
                         {2U\over\lambda_1}
 +{\varepsilon+2\psi U\over\lambda_2}.             \tag{1.10}
\]

No attainment of the Cheeger infimum is used.

## 2. A sharp transition-layer lower bound

The preceding integral statement has a geometric complement which keeps
the smoothing scale explicit.

**Lemma 2.1.**  Suppose `||grad F||_infinity<=L_F`.  For every
`a in (0,1/2)`,

\[
\boxed{
 U\ge {a\over2}\left[1-
 \exp\left\{-{\psi(1-2a)\over2L_F}\right\}\right].} \tag{2.1}
\]

**Proof.**  Put

\[
 L=\{F\le a\},\qquad H=\{F\ge1-a\}.
\]

Equation (1.2) and Markov's inequality give

\[
 \mu(L),\mu(H)\ge m:={1\over2}-{U\over a}.           \tag{2.2}
\]

If `m<=0`, (2.1) is automatic.  Otherwise Lipschitzness gives

\[
                         \operatorname {dist}(L,H)
 \ge d:={1-2a\over L_F}.                              \tag{2.3}
\]

For a set of mass at least `m<=1/2`, the Cheeger differential inequality
for its open neighborhoods shows that its `r`-neighborhood has mass at
least `m exp(psi r)` until it first reaches mass `1/2`, and strictly more
than `1/2` afterwards.  (This is obtained by integrating the lower Dini
derivative inequality for the neighborhood-volume function; equivalently
it is the standard Cheeger-to-exponential-concentration implication.)  For
every `r<d/2`, the open `r`-neighborhoods of `L` and `H` are disjoint.  They
therefore cannot both have mass greater than `1/2`.  Letting `r` increase
to `d/2` gives

\[
                         m e^{\psi d/2}\le {1\over2}. \tag{2.4}
\]

Substitute (2.2)--(2.3) and rearrange. QED.

The proof also gives the crude middle-layer estimate

\[
 \mu\{a<F<1-a\}\le {2U\over a}.                     \tag{2.5}
\]

For the Gaussian posterior-resampling regularization in
`heatflow_bernstein.md`, one has

\[
                         L_F\le {I(1/2)\over\sqrt s}. \tag{2.6}
\]

Hence

\[
\boxed{
 U(s)\ge {a\over2}\left[1-
 \exp\left\{-{\psi(1-2a)\sqrt s\over2I(1/2)}\right\}\right].} \tag{2.7}
\]

When `psi sqrt(s)<=1`, this has the sharp scale

\[
                         U(s)\ge c_a\psi\sqrt s.     \tag{2.8}
\]

The already proved heat-content upper bound for a balanced near-minimizer
is `U(s)<=C sqrt(s) p`.  Since `p=psi/2+o(psi)`, the two estimates have
the same order.  Thus the square-root amplitude is a real transition-layer
cost, not a loss caused by posterior resampling.

## 3. Normal field and the exact Jensen gap

For the Gaussian channel notation, let

\[
 W=\nabla g_s(Y),\qquad
 m(X)=E[W\mid X]=\nabla F(X),\qquad
 \theta(X)={m(X)\over|m(X)|}                         \tag{3.1}
\]

where the value on `{m=0}` is arbitrary.  If
`u(Y)=W/|W|`, the norm-Jensen deficit is exactly

\[
 \delta_J=E|W|-E|m(X)|
 ={1\over2}E\{|W|\,|u(Y)-\theta(X)|^2\}.             \tag{3.2}
\]

For the balanced Cheeger near-minimizer,

\[
                         \delta_J\le
 \varepsilon+4pU(s).                                  \tag{3.3}
\]

Thus the posterior phase direction aligns, in the boundary-flux weight,
with the ordinary normal field of every regular level set of `F`.  Equations
(1.5)--(1.10) show that the same nested level sets are near-global
isoperimetric sets for most thresholds.

What is not supplied by these identities is a finite competitor which
turns variation of `theta` between different boundary patches into a
strict improvement of (1.9).  The coarea deficit is additive over patches;
it has no cross term comparing their normals.  Such a cross term must come
from the common nested geometry, the full second variation, or a multiway
incidence argument.

## 4. Model audit and exact limitation

* For a Gaussian halfspace, the transition layer is one-dimensional and
  (2.7) has the correct `psi sqrt(s)` order.  The Jensen direction is
  constant.
* For a centered Gaussian ball, the normal field is radial and the nested
  levels are concentric.  Any valid rigidity conclusion must retain this
  branch rather than demand a single direction.
* For the inner square in a cube and for Gaussian or exponential maximum
  cuts, the small-heat level sets inherit many flat phase patches and the
  covariant angular energy tends to zero.  Their integrated quantity in
  (1.9), however, is not small relative to the global Cheeger optimum; they
  are therefore not counterexamples to an extremality-dependent theorem.
* Conversely, (1.5) alone cannot distinguish a disjoint union of exact
  one-dimensional Cheeger pieces when the isoperimetric profile is linear:
  every component pays exactly `psi` times its smaller mass and the deficits
  add to zero.  Log-concavity of the ambient measure or incidence of the
  pieces must be used to rule out that formal model.

The rigorous output is the exact cascade (1.5)--(1.10), the transition
lower bound (2.1), and the flux-normal alignment (3.2)--(3.3).  They show
precisely where near-global minimality enters the phase problem.  They do
not by themselves prove affine/radial coherence, and no Cheeger expansion
for a conditioned phase law has been inserted.
