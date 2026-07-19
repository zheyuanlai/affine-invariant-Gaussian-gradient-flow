# Root audit of the local spatial phase lemmas

## Scope and verdict

This is a line-by-line algebraic audit of the candidate lemmas in
`spatial_phase_coherence.md`. It is not the independent clean-room audit,
which is being run separately. Lemmas 2.1, 3.1, and the Brenier product
estimate (4.0f) survive the checks below. The midpoint resolvent transfer is
correct in the smooth bulk setting but is not a hard-support statement. The
two-dimensional halfspace lemma and the equality splitting assertion still
need fully self-contained approximation proofs before they can be marked
load-bearing.

## 1. Threshold approximation

After scaling to `t=1`, write the active marginal as `T(Z)=Z+R(Z)`, where
`T` is the increasing Gaussian contraction, `R` is nonincreasing, and the
threshold is `z_0=Phi^{-1}(1-g)`. If `H=1_{Z>=z_0}`, then direct centering
gives

\[
 D_{map}=I(g)-E[T(Z)H]
 =g(1-g)\{E(R\mid Z<z_0)-E(R\mid Z>z_0)\}.          \tag{1.1}
\]

Let `p` be the marginal mass mapped into the target strip of width `b`
around the threshold. If its upper part has mass at least `p/2`, divide
that upper Gaussian preimage interval into equal-mass near and far pieces.
The near piece has Gaussian-coordinate length at least `c p`, since the
standard Gaussian density is at most `(2pi)^(-1/2)`, while its full image
has length at most `b`. Hence between `z_0` and every point in the far
piece, `R` drops by at least `c p-b`. In the double-integral form of (1.1),
pair the far piece, of unconditional mass at least `p/4`, with the entire
lower half, of mass `1-g>=delta`. This gives

\[
 D_{map}\ge c_delta p(cp-b)_+.                     \tag{1.2}
\]

The lower-strip case is symmetric. Therefore
`p<=C_delta(b+sqrt(D_map))`. Outside the strip, positivity of
`(y-a)(H-q)` gives threshold error at most `D_cut/b`. Since
`Delta=I(g)epsilon`, choosing `b=sqrt(epsilon)` proves
`pi(S triangle H)<=C_delta sqrt(epsilon)`. Restoring scale replaces both
`b` and the two defects by their dimensionless `sqrt(t)` versions and does
not change the conclusion.

The pairing uses unconditional masses; this is why no missing factor
`g(1-g)` occurs in (1.2). Generalized quantiles handle atoms, although a
strongly log-concave one-dimensional marginal is in fact absolutely
continuous on the relative interior of its support.

## 2. Variance rigidity

With `q=-R'` in the a.e. sense, `0<=q<=1`. Gaussian integration by parts
and the threshold-flux representation give

\[
 D_{map}=E[k_g(Z)q(Z)],\qquad
 k_g(z)\ge {c_delta\over1+|z|}.                     \tag{2.1}
\]

For `L>=1`, split at `|Z|=L`:

\[
 Eq\le C_delta(1+L)D_{map}+P(|Z|>L).               \tag{2.2}
\]

Taking `L=sqrt(2log(e/D_map))` yields
`Eq<=C_delta D_map sqrt(log(e/D_map))`. Since `ER=0`,

\[
 1-Var(T(Z))=2Eq-ER^2\le2Eq.                       \tag{2.3}
\]

This proves the claimed logarithmic modulus. The formula remains valid for
nonsmooth contractions by approximating the monotone 1-Lipschitz map in
Gaussian Sobolev space.

## 3. Brenier product estimate

For the standardized Brenier contraction, its a.e. Jacobian `H` is
symmetric with `0<=H<=I`. The variance bound and Gaussian Poincare imply
`E|Hu|^2>=1-zeta`. The spectral inequality

\[
 |Hu-u|^2\le1-|Hu|^2                               \tag{3.1}
\]

follows from `H^2<=H`. Hence `E|Hu-u|^2<=zeta`.
For the active coordinate, Gaussian Poincare applied to
`<u,T(G)>-<u,G>` gives squared error at most `zeta`; for the transverse
coordinate, one-dimensional Gaussian Poincare on each active Gaussian
fiber gives the same bound because its derivative is `PHu=P(Hu-u)`.
The components are orthogonal, proving the total `2zeta` estimate.

Caffarelli contraction for an extended-valued convex target is invoked via
monotone smooth convex approximation. To use this as a final-proof lemma,
one must state the convergence of the Brenier maps in `L^2(gamma)` and weak
lower semicontinuity of their Sobolev energies; this has not yet been written
out in the source note.

## 4. Midpoint and line-coherence statements

For endpoint separation `d=c_1-c_2`, the exponential-family identity is

\[
 E_{c_i}\left({d\pi_{c_m}\over d\pi_{c_i}}\right)^2
 ={Z(c_1)Z(c_2)\over Z(c_m)^2}
 \le \exp\{|d|^2/(4t)\}.                            \tag{4.1}
\]

It follows from `nabla^2 log Z=Cov<=t^{-1}I`. Cauchy--Schwarz therefore
transfers a `[0,1]`-valued bulk resolvent defect with a square-root loss.
This calculation is correct, but the resolvent integrand ignores hard
support contact and must not replace the Brenier certificate there.

For two halfspaces, whitening the at-most-two-dimensional marginal turns
their covectors into the normalized vectors `A^(1/2)u_i`. A proof of the
claimed linear angle bound can be obtained from uniform fixed-dimensional
central-quantile bounds and a uniform density lower bound on an inner ball.
The current source note only sketches those two inputs. They must either be
proved with constants or replaced by an exact cited fixed-dimensional
compactness lemma before (5.5) is audited.

Likewise, exact equality yields `Hu=u` almost everywhere for the Brenier
map and should imply Gaussian product splitting. The smooth full-support
case follows by integrating the vanishing mixed derivatives. A complete
hard-support, lower-dimensional proof still has to pass this identity
through approximation and identify the product support. Until that is
done, Section 6 is a correct rigidity blueprint rather than an audited
general theorem.
