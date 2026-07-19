# Stability decomposition for the sharp posterior centroid bound

## 1. Statement

Let `t>0`.  Let `pi` be a `t`-strongly log-concave probability on
`R^n`, let `S` be Borel, and write

\[
 g=\pi(S),\qquad
 v=\int_S(x-\mathbb E_\pi X)\,d\pi(x).
\]

Assume `v ne 0`, put `u=v/|v|`, and let `Y=<X-E X,u>`.  Denote by
`rho` the law of `Y` and put

\[
 q(y)=\pi(S\mid Y=y).
\]

Let `c` be an upper `g`-quantile of `rho`, chosen so that
`H(y)=1_{y>=c}` has mass `g`; the usual randomized convention at an atom
may be used, although strong log-concavity gives a continuous density.
Let `Z` be standard Gaussian, let `T` be the increasing transport from
`Z` to `Y`, and set

\[
 L=t^{-1/2},\qquad R(z)=T(z)-Lz,qquad
 z_g=\Phi^{-1}(1-g).
\]

The one-dimensional Caffarelli theorem says that `R` is nonincreasing.
Define

\[
 \Delta={\mathcal I(g)\over\sqrt t}-|v|,
 \qquad \mathcal I(g)=\varphi(\Phi^{-1}(g)).
\]

**Proposition 1 (exact two-deficit identity).**  One has

\[
 \boxed{\Delta=D_{\rm cut}+D_{\rm map},}                 \tag{1.1}
\]

where

\[
 \begin{aligned}
 D_{\rm cut}
   &=\int (y-c)(H(y)-q(y))\,d\rho(y)\ge0,\\
 D_{\rm map}
   &=-\mathbb E[R(Z)1_{\{Z\ge z_g\}}]\\
   &=g(1-g)\left(
       \mathbb E[R(Z)\mid Z<z_g]
       -\mathbb E[R(Z)\mid Z\ge z_g]
     \right)\ge0.
 \end{aligned}                                          \tag{1.2}
\]

Thus near equality in the sharp centroid bound simultaneously forces the
cut to be close to a threshold in its active marginal and forces the
one-dimensional Caffarelli map to be close to its maximal slope across the
two relevant quantile blocks.  No ambient dimension occurs.

## 2. Proof of the identity

Translate `pi` so that its barycenter is zero.  The definition of `u`
gives

\[
 |v|=\int yq(y)\,d\rho(y).                              \tag{2.1}
\]

Because `H` and `q` have the same mean,

\[
 \begin{aligned}
 \int y(H-q)d\rho
 &=\int (y-c)(H-q)d\rho\\
 &=\int_{y\ge c}(y-c)(1-q(y))d\rho(y)
   +\int_{y<c}(c-y)q(y)d\rho(y)\ge0.                  \tag{2.2}
 \end{aligned}
\]

This is exactly `D_cut`.  It is the quantitative bathtub-principle
defect.

Under the monotone transport, `H(T(Z))=1_{Z\ge z_g}`.  Since `EY=0`
and `EZ=0`, one has `ER(Z)=0`.  Therefore

\[
 \begin{aligned}
 {\mathcal I(g)\over\sqrt t}-\int yH(y)d\rho(y)
 &=L\,\mathbb E[Z1_{\{Z\ge z_g\}}]
   -\mathbb E[T(Z)1_{\{Z\ge z_g\}}]\\
 &=-\mathbb E[R(Z)1_{\{Z\ge z_g\}}].                 \tag{2.3}
 \end{aligned}
\]

Writing the last covariance by conditioning on the two sides gives the
second formula for `D_map`.  It is nonnegative because `R` is
nonincreasing.  Adding (2.2) and (2.3), and using (2.1), proves (1.1).

## 3. Explicit consequences on central states

Fix `delta in (0,1/2)` and suppose

\[
 \delta\le g\le1-\delta.                               \tag{3.1}
\]

### 3.1 Threshold error away from the interface

For every `a>0`, (1.2) gives

\[
 \boxed{
 \int_{|y-c|\ge a}|H(y)-q(y)|\,d\rho(y)
 \le {\Delta\over a}.}                                \tag{3.2}
\]

Indeed, the two integrands in (2.2) equal
`|y-c||H-q|` on their respective half-lines.  Thus all disagreement not
confined to an `a`-strip about the threshold is paid linearly by the
centroid deficit.

### 3.2 Quantile-slope error of the transport

Let `a_-<z_g<a_+` be such that

\[
 \mathbb P\{a_-\le Z<z_g\}\ge\eta,
 \qquad
 \mathbb P\{z_g\le Z\le a_+\}\ge\eta                 \tag{3.3}
\]

for some `eta>0`.  Monotonicity of `R` and (1.2) imply

\[
 \boxed{
 0\le L(a_+-a_-)-\{T(a_+)-T(a_-)\}
 =R(a_-)-R(a_+)
 \le {\Delta\over\eta^2}.}                            \tag{3.4}
\]

To verify the last inequality, let `Z_-` and `Z_+` have the Gaussian law
conditioned below and above `z_g`, respectively.  On the event
`Z_- in [a_-,z_g)` and `Z_+ in [z_g,a_+]`, whose conditional probability
is at least `eta^2/[g(1-g)]`, monotonicity gives

\[
 R(Z_-)-R(Z_+)\ge R(a_-)-R(a_+).
\]

Since

\[
 D_{\rm map}=g(1-g)\mathbb E[R(Z_-)-R(Z_+)],
\]

equation (3.4) follows.  Under (3.1), fixed central Gaussian quantiles may
be chosen with `eta=c_delta`, so the loss is a constant depending only on
`delta`.

### 3.3 Integrated form for stochastic localization

For ordinary stochastic localization of a fixed set, put

\[
 g_t=\mu_t(S),\quad v_t=\operatorname {Cov}_{\mu_t}(1_S,X),
 \quad
 r_t={\sqrt t|v_t|\over\mathcal I(g_t)}\in[0,1].       \tag{3.5}
\]

The Gaussian-profile drift identity gives, for `0<s<T`,

\[
 F(T)-F(s)
 ={1\over2}\mathbb E\int_s^T
 {\mathcal I(g_t)\over\sqrt t}(1-r_t^2)dt,             \tag{3.6}
\]

where `F(t)=sqrt(t) E I(g_t)`.  On the central event
`g_t in [delta,1-delta]`,

\[
 {\mathcal I(g_t)\over\sqrt t}(1-r_t^2)
 =(1+r_t)\left({\mathcal I(g_t)\over\sqrt t}-|v_t|\right).
                                                               \tag{3.7}
\]

Consequently (3.6) is exactly the time integral of the two nonnegative
defects (1.2), up to the factor `1+r_t in [1,2]`.  A small total profile
drift therefore forces posterior cuts to be nearly one-dimensional
thresholds and their active marginals to be nearly Gaussian dilations for
most *weighted central occupation time*.

This does not by itself prove a seed.  The unproved global step is a
dimension-free compatibility theorem gluing the active directions and
near-Gaussian marginal coordinates across different posterior states and
times.  Treating (3.2)--(3.4) as if they supplied that compatibility would
be circular.

## 4. Equality audit

If `Delta=0`, then `D_cut=0` and (2.2) gives `q=H` away from the threshold.
Moreover `D_map=0`.  Since `R` is nonincreasing and both half-lines have
positive Gaussian mass, `R` is constant almost everywhere across their
product; hence `T(z)=Lz+b` almost everywhere.  The active marginal is a
Gaussian of variance `1/t`, and `S` is, modulo `pi`-null sets, a halfspace
orthogonal to `u` after conditioning on that marginal.  This recovers the
sharp equality case without invoking a multidimensional stability theorem.

For a posterior of an initially `kappa`-strongly log-concave Gaussian
halfspace, the correct curvature in the sharp bound is `kappa+t`, and the
same calculation has zero deficit.  With only the posterior curvature `t`,
the ratio in (3.5) is strictly below one, as it must be.  Thus the identity
does not falsely classify a finite-time Gaussian posterior as an equality
case for the weaker curvature parameter.
