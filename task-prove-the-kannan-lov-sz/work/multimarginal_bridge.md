# A null-invariant multi-packet barycenter inequality

The indicator midpoint count is unstable under null modifications.  This
note replaces it by an exact density and entropy estimate for actual
transported representations.  The multi-packet form has a determinant gain
which is invisible in pairwise midpoint estimates.

## 1. Common-source displacement barycenters

Let

\[
 d\mu(x)=Z^{-1}e^{-V(x)}dx=:\rho(x)dx
\]

be a full-dimensional log-concave probability.  Let `A_1,...,A_M` be Borel
sets of positive masses `a_i=mu(A_i)`, and write

\[
 d\mu_i={1_{A_i}\rho\over a_i}dx.
\]

Let `lambda_i>0`, `sum_i lambda_i=1`.  Take any absolutely continuous
common source `gamma` and Brenier maps

\[
 T_i=\nabla\phi_i,\qquad (T_i)_\#\gamma=\mu_i.
\]

Put

\[
 S=\sum_i\lambda_iT_i,qquad \nu=S_\#\gamma.             \tag{1}
\]

The map `S` is again a gradient of a convex function.  First assume all
densities and maps are smooth and `DT_i` are positive definite.  At
`z=S(x)`, define

\[
 \Delta_V(x)=\sum_i\lambda_iV(T_i(x))-V(S(x))\ge0        \tag{2}
\]

and

\[
 Q(x)={\det(\sum_i\lambda_iDT_i(x))
       \over\prod_i\det(DT_i(x))^{\lambda_i}}\ge1.       \tag{3}
\]

The last inequality is concavity of `log det`.

**Proposition 1 (multi-packet density bound).**  At almost every such point,

\[
 \boxed{\quad
 {d\nu\over d\mu}(S(x))
 \le {e^{-\Delta_V(x)}\over
           Q(x)\prod_i a_i^{\lambda_i}}.
 \quad}                                                  \tag{4}
\]

Consequently

\[
 \boxed{\quad
 D(\nu\Vert\mu)
 +E_\gamma[\Delta_V+\log Q]
 \le\sum_i\lambda_i\log {1\over a_i}.
 \quad}                                                  \tag{5}
\]

In particular,

\[
 E_\gamma[\Delta_V+\log Q]
 \le\sum_i\lambda_i\log {1\over a_i}.                  \tag{6}
\]

**Proof.**  Write the density of `gamma` as `g`.  Monge--Ampere gives

\[
 g(x)={\rho(T_i(x))\over a_i}\det DT_i(x),
 \qquad
 {d\nu\over dx}(S(x))={g(x)\over\det DS(x)}.            \tag{7}
\]

Convexity of `V` says

\[
 \rho(S(x))\ge e^{\Delta_V(x)}
                 \prod_i\rho(T_i(x))^{\lambda_i}.       \tag{8}
\]

Substitute (7) into (8), cancel the common source density, and use
`DS=sum lambda_i DT_i`; this is exactly (4).  Taking logarithms and
integrating against `gamma`, then using nonnegativity of relative entropy,
proves (5)--(6).  QED.

The statement extends to nonsmooth log-concave densities and Alexandrov
derivatives of Brenier maps by Gaussian regularization, truncation to convex
compact supports, and lower semicontinuity of relative entropy.  Any use in a
final proof must spell out that approximation; the smooth formula is the
only part needed for the present mechanism audit.

For equal packet masses `a_i=beta/M` and equal weights, the right side of
(6) is

\[
 \log(M/\beta),                                         \tag{9}
\]

not `M log(M/beta)`.  Thus simultaneous averaging retains exactly the packet
entropy scale.

## 2. Why the multi-map determinant can be stronger than pairs

The following constant-matrix model calibrates (3).  Work in `R^M`, take
orthonormal vectors `e_i`, and put

\[
 A_i=I+(s-1)e_ie_i^T\qquad(s\ge1).                       \tag{10}
\]

Then `det A_i=s`, while

\[
 {1\over M}\sum_iA_i=\left(1+{s-1\over M}\right)I.
\]

Therefore

\[
 \log Q_M
 =M\log\left(1+{s-1\over M}\right)-\log s.            \tag{11}
\]

If `M>=s^2`, then

\[
 \log Q_M\ge s-2-\log s.                               \tag{12}
\]

By contrast, the two-map Jensen gap between `A_i` and `A_j` is only of
order `log s`.  Hence a simultaneous barycenter is capable, at the purely
matrix level, of producing the linear-in-`s` gain needed to compete with
`log M` when `M=exp(Theta(s))`.

This is the correct calibration for packets whose conditional laws have one
scale-`s` direction and whose long directions are dispersed.  The determinant
gain vanishes for translated parallel packets (`DT_i` identical), exactly as
it should.

## 3. The unresolved derivative-diversity lemma

To use (6), one would need a lower bound of the following kind.  For
positive- or negative-tail packet laws produced by balanced rays of scale
`s`, and for a common source chosen canonically,

\[
 E\log Q\ge c s-C                                      \tag{13}
\]

unless the packet laws are approximately translates, concurrent radial
images, or the two-block orthogonal-radial equality model.

Long one-dimensional conditional variance alone does not yet prove (13).
Using a Gaussian common source and Brenier maps gives the valid covariance
estimate

\[
 \operatorname {Cov}(T_i(G))\preceq E[DT_i(G)^2]         \tag{14}
\]

by Gaussian Poincare, but (14) is only an `L2` statement.  The pointwise
log-determinant gap in (3) can in principle be evaded if the large derivatives
of different maps occur on disjoint source events.  No dimension-free
synchronization estimate has been proved for the nonconvex ray packets.

Thus Proposition 1 is a genuinely stronger, null-invariant bridge than the
Minkowski midpoint count, and (11) shows the desired scaling in the canonical
tube model.  The load-bearing next lemma is synchronization of the long
Brenier stretches.  It cannot be replaced by covariance comparison or by a
pairwise determinant estimate.
