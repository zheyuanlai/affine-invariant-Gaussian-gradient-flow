# Exact reduction layer

This file records only reductions; it does **not** contain the central
dimension-free inequality.

## 1. The ambient space may be replaced by the affine hull

Let `mu` be supported on an affine subspace `E` of `R^n`, and let `mu_E^+`
denote exterior Minkowski content computed with the Euclidean metric induced
on `E`. Then the ambient and relative Cheeger constants coincide.

Indeed, for every Borel `A subset R^n`, put `B=A intersect E`. We have
`mu(A)=mu(B)` and, for every `eps>0`,

`(B_eps^E) subset (A_eps intersect E)`.

Consequently `mu^+(A) >= mu_E^+(B)`. Conversely, if `B subset E` is regarded
as a subset of `R^n`, then for `x in E`,
`dist_{R^n}(x,B)=dist_E(x,B)`, and hence
`mu(B_eps)=mu(B_eps^E)`. Taking the two infima proves equality. This also
handles ambient Borel sets containing points off the support.

## 2. Covariance on the affine hull

Write `E=m+L`, where `m` is the barycenter and `L` is the direction space.
For `v in L`,

`v^T Cov(mu) v = int <x-m,v>^2 dmu(x)`.

If this is zero, the support is contained in the affine hyperplane
`m+(L intersect v^perp)`. Thus, when `E` is the minimal affine hull of the
support, `Cov(mu)|_L` is positive definite. On `L^perp` the covariance is
zero, so its ambient operator norm is exactly the largest eigenvalue of its
restriction to `L`. A non-point log-concave measure has `dim L >= 1` and
therefore has strictly positive covariance operator norm.

## 3. Boundary content under an invertible linear map

Let `S:E -> F` be an invertible linear map between Euclidean spaces of the
same dimension, let `nu=S_#mu`, and put `L=||S||_op`. For every Borel
`A subset E`,

`A_{eps/L} subset S^{-1}((SA)_eps)`.

It follows directly from the definition and the substitution
`delta=eps/L` that

`nu^+(SA) >= mu^+(A)/||S||_op`.

Applying this to every set and then to `S^{-1}` gives the two-sided
bi-Lipschitz comparison

`psi_mu/||S||_op <= psi_nu <= ||S^{-1}||_op psi_mu`,

or, equivalently,

`psi_mu >= psi_nu/||S^{-1}||_op`.

No boundary regularity or attainment of the infimum is used.

## 4. Isotropic form implies the full covariance form

Assume there is a universal `c>0` such that every full-dimensional isotropic
log-concave probability on every `R^k`, `k>=1`, has Cheeger constant at least
`c`. Let `mu` be any non-point log-concave probability on `R^n`. By Borell's
characterization it has a log-concave density on its minimal affine hull
`E=m+L`; by Sections 1 and 2 we work on `L` and its covariance `A` is positive
definite there. Define

`S=A^{-1/2}` and `nu=S_#(mu translated by -m)`.

Then `nu` is isotropic and log-concave on the `k`-dimensional Euclidean space
`L`. Section 3, applied in the direction from `nu` back to `mu`, yields

`psi_mu >= psi_nu/||S^{-1}||_op
        >= c/||A^{1/2}||_op
        = c/sqrt(||Cov(mu)||_op)`.

Conversely, the full covariance statement applied to an isotropic measure has
`||Cov(mu)||_op=1`, so it gives the isotropic statement. Translations and
orthogonal identifications do not change either side. This proves the claimed
equivalence, including degenerate ambient covariance and lower-dimensional
support. A point mass is precisely the zero-dimensional case and is excluded.

## 5. Poincare constants under affine maps

Let `C_P(mu)` be the optimal constant in

`Var_mu(f) <= C_P(mu) int |grad f|^2 dmu`.

For the same `S` as in Section 3, the chain rule gives

`C_P(S_#mu) <= ||S||_op^2 C_P(mu)`,

and applying the statement to `S^{-1}` gives

`C_P(mu) <= ||S^{-1}||_op^2 C_P(S_#mu)`.

Thus a dimension-free isotropic Poincare bound implies
`C_P(mu) <= C ||Cov(mu)||_op` after whitening. If a locally Lipschitz test
function has infinite Dirichlet integral, the assertion is vacuous; otherwise
its restriction to the affine hull lies in the corresponding first-order
Sobolev space. The usual density of compactly supported smooth functions in
that weighted Sobolev space supplies the full test class when a proof is first
established for smooth tests.

Under `CD(0,infinity)` (and hence for log-concave Euclidean measures after the
standard quantitative approximation), Cheeger--Maz'ya and Buser--Ledoux state,
in the normalization used here,

`C_P(mu) <= 4/psi_mu^2`, and `psi_mu >= c_BL/sqrt(C_P(mu))`.

Therefore the Cheeger and Poincare targets are equivalent up to universal
constants. The approximation step cannot be omitted if Buser--Ledoux is
invoked only in its smooth formulation.

## 6. Mean and median in the first-moment target

For any integrable real function `f` and any median `m`, the median minimizes
the `L^1` distance to constants, and

`int |f-m| dmu <= int |f-Ef| dmu
 <= 2 int |f-m| dmu`.

The second inequality follows from
`|Ef-m| <= int|f-m|dmu`. Hence the mean-centered and median-centered versions
of the first-moment target are equivalent with constants at most two. E.
Milman's theorem for metric-measure spaces satisfying the log-concave/
`CD(0,infinity)` convexity hypotheses identifies the reciprocal of the worst
such first moment with the Cheeger constant, up to universal constants.

## 7. Immediate restrictions on a bad Lipschitz witness

Let `X` be isotropic. If `ell(x)=<a,x>+b` is 1-Lipschitz, then

`E|ell(X)-E ell(X)| <= sqrt(Var ell(X))=|a|<=1`.

If `f(x)=g(|x|)` with `g` 1-Lipschitz, the dimension-free thin-shell theorem
gives

`E|f(X)-Ef(X)|
 <= 2 E|g(|X|)-g(sqrt(n))|
 <= 2 E||X|-sqrt(n)| <= C`.

Thus a sequence of witnesses with diverging centered first moment has
diverging `L^1` distance from every centered 1-Lipschitz linear or radial
witness. This is only a structural restriction, not a proof for general
Lipschitz functions.
