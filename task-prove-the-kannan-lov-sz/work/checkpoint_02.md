# Checkpoint 02: full junction stability and the Clifford equality branch

## 1. Exact results added since Checkpoint 01

### 1.1 Full signed-distance second variation

For a globally maximizing T3 witness, smooth normal deformation gives the
exact boundary-gain inequality in `signed_distance_stability.md`.  At a
generic medial interface between two distance charts, the missing second
epi-derivative is

\[
 \int_{M_{ij}}\frac{(h_i-h_j)^2}{|N_i-N_j|}\rho\,dH^{n-1}.
\]

For the normal-coordinate heights `h=a dot N`, summation over an ambient
orthonormal basis charges precisely

\[
 \int_{M_{ij}}\rho|N_i-N_j|\,dH^{n-1}.
\]

The global translation family `f(x-epsilon a)` bypasses arbitrary focal-set
regularity and yields the full trace identity

\[
 \int_{f=0}|N-m_N|^2\rho\,dH^{n-1}
 +\frac12E[(sign f-(p-q))<grad V,grad f>]
 \le P_mu(f>0).
\]

The second term is nonnegative by ray balance and one-dimensional convexity.
This identity includes every smooth, medial, multiway and focal contribution.

### 1.2 Sharp limitation of infinitesimal stability

If all relevant ray scales are at least `s` and their boundary densities are
comparable to `1/s`, stability forces total normal energy at least `c/s`.
Energy `o(1/s)` gives `L2` normal alignment and isotropy then gives `s=O(1)`.

This branch is sharp.  An abstract complete graph on
`M>=exp(c s)` orthogonal normal packets, with node boundary density `c/s` and
edge conductance `c/(sM^2)`, saturates the full height stability inequality,
the covariance constraint, the Hilbert--Schmidt curvature constraint and the
small-packet entropy constraint.  Therefore no further abstract manipulation
of these inequalities can force alignment.

### 1.3 Gaussian fan audit

For the standard-Gaussian `2m`-sector fan, the smooth charge is zero and the
singular charge is

\[
 P\cos(\pi/(2m))<P.
\]

Hence the fan is correctly excluded as a global T3 extremizer.  Its deficit
is only of relative order `m^{-2}` (absolute order `m^{-1}`), so this example
does not supply a uniform gap.

### 1.4 Orthogonal-spherical (Clifford) equality branch

If two endpoint sets have constant cross-distance `d`, their difference
spans are orthogonal and each set lies on a sphere in its affine span.  For
independent endpoint laws, the exact identity

\[
\begin{aligned}
 Var|B-C|^2={}&Var(|U|^2+2<a,U>)
 +Var(|V|^2-2<a,V>)\\
 &+4Tr(Cov(B)Cov(C))
\end{aligned}
\]

is a robust measure-level version of the same classification.

This third branch is harmless.  Two radial coordinates in orthogonal
marginals, together with the one possible separating linear coordinate,
form a one-Lipschitz feature map.  Translated thin shell controls the radial
coordinates and isotropy controls the linear coordinate.  Thus two
constant-mass tail bands concentrated near the classified feature values
cannot be separated by more than a universal constant.  The correct target
classification is therefore parallel, radial-on-a-subspace, or
orthogonal-radial, rather than merely parallel or concurrent.

## 2. Routes closed or merged

1. **Known stochastic bootstrap.**  The exact recurrence is
   `K_n<=C sqrt(K_n log n)` and has the fixed point `K_n=O(log n)`.
   Product exponentials have universal Poincare constant but genuinely
   develop top posterior variance of order `log n` at time `1/log n`, so the
   covariance horizon cannot be replaced by `log K_n`.
2. **Mass-preserving localization.**  Late-time rotation of the exceptional
   direction is dimension-free; only creation of large variance in the
   initial `O(1/R)` layer remains.  Every coherent high-variance cap has
   probability at most `C exp(-c sqrt R)`, but diffuse adaptive directions
   remain.  This is the same global packet-gluing obstruction as the ray
   route.
3. **Polar charts.**  Radial fluctuations are completely controlled and the
   angular law is an exactly conditioned critical convex measure.  Local
   chart coercivity is false, while the required global gluing inequality is
   quantitatively KLS itself.  No independent bridge emerged.
4. **Smooth or singular trace energy alone.**  Smooth curvature has only the
   perimeter scale `1/s`; medial charge has no known upper bound and can pay
   the entire stability inequality.  A finite competitor or Euclidean
   realization theorem is indispensable.

## 3. Current load-bearing statement

Let a true global T3 extremizer have a positive-mass family of balanced rays
of scale `s`.  Prove that the full Euclidean, density-weighted medial gluing
obeys one of the following alternatives with universal quantitative errors:

1. low gluing energy, which gives normal alignment;
2. many alternative near-calibrated same-ray rematchings, whose endpoint
   geometry is close to the orthogonal-spherical classification; or
3. a finite admissible signed-distance competitor with strictly larger
   centered first moment.

The missing implication must be null-invariant and must use actual reflected
representations, not membership in Minkowski midpoint sets.  Indicator
midpoint multiplicity is unstable under null modifications and does not
respect the calibrated matching.

## 4. Active independent tests

- A Euclidean/log-concave exclusion of expander-like medial realizations,
  based on positive-density bridges and same-ray matching.
- A finite normal-offset/rematching competitor which remains effective when
  the infinitesimal junction graph has a constant spectral gap.
- An adversarial construction search for a genuine log-concave realization;
  finding one would invalidate the proposed inverse statement and force a
  different mechanism.

No complete central lemma has yet passed audit.

