# Checkpoint 1

## Candidate status

No complete central proof exists yet.  The affine/lower-dimensional reduction
layer is rigorous, but every proposed central estimate below has either been
reduced to an explicit new compatibility statement or invalidated by a
dimension-tracked counterexample.  Nothing in this checkpoint is returnable as
a proof of KLS.

## Surviving mechanisms

1. **Global balanced-ray compatibility.**  A true T3 extremizer exists and is
   an exact signed-distance/Kantorovich potential.  Its transport-ray
   disintegration balances the sign cut on every ray, and a constant quotient
   mass of rays has two-sided length comparable to the bad first-moment scale.
   Thin shell makes those rays almost tangent.  The missing statement is a
   global theorem coupling long scale to the full smooth and singular
   variation of one normal congruence under a connected log-concave density.

2. **Small-eigenvalue rigidity.**  For a first eigenfunction with eigenvalue
   `lambda`, Bochner gives Hessian energy `O(lambda^2)` and gradient energy
   `lambda`.  A dimension-free almost-splitting or tail-energy theorem for the
   resulting curl-free field would close T2.  Applying Poincare separately to
   its components is forbidden because it is exactly the target.

3. **Small-perimeter dimension descent.**  Boundary stability gives a precise
   curvature/normal-alignment dichotomy, and the aligned branch can be repaired
   to a subgraph.  A valid induction still needs control of a conditional-score
   velocity term specifically in the small-perimeter stationary regime.

## New exact structural facts

- `D_1(mu)=2 sup_E p(1-p)W_1(mu_E,mu_Ec)`.
- Every extremizer is eikonal almost everywhere and is signed distance to its
  zero separator; every transport needle has the global sign proportions.
- A bad near-extremizer has constant-mass separated tails, bounded conditional
  barycenters, high-stable-rank displacement, and long nearly tangent saturated
  rays.
- Every log-concave cut has a universal gap from perfect linear prediction:
  `v^T A^{-1}v <= (1-eta)p(1-p)`.  The isotropic interval attains ratio `3/4`,
  so this slack alone cannot beat a square-root curvature loss.
- For covariance-normalized localization, set-mass quadratic variation is at
  most `g(1-g)`, covariance drift is `-A`, and accumulated curvature is
  `int A^{-1}`.
- Posterior resampling `K_T` is positive reversible and
  `<h,(I-K_T)h>=E Var_{mu_T}h`, with `I-K_T>=exp(-T)I`.

## Blocked mechanisms and exact reasons

- Projection thin-shell inequalities alone permit a `log n` quadratic-form
  obstruction (`Q(H)=(Tr(BH))^2`, eigenvalues of `B` equal `i^{-1/2`).
- Heat-semigroup and Witten-inverse estimates proposed so far have optimal
  constants exactly equal to `C_P`.
- Random-line resampling has gap `Theta(n^{-2})` on cube/simplex and
  `Theta(n^{-1})` on ball/Gaussian; its optimal conditional weight has rare
  central-chord spikes, so separate gap/weight comparison fails.
- Global average-Jacobian transport bounds do not localize target energy;
  Gaussian-to-Laplace gives an explicit tail counterexample.
- Covariance-normalized endpoint curvature cannot be transferred through full
  anisotropic perimeter: a one-dimensional exponential signal makes the
  likelihood-weighted metric grow like `sqrt(x)`.
- Scalar corrections depending on posterior mass and its linear projection
  ratio remain blind to irrelevant Laplace tail faces in a product decoration.
- The face-specific averaged `L^2` calibration norm equals
  `||L^{-1/2}(I-K_T)||`, which lies between
  `exp(-T)sqrt(C_P)` and `sqrt(C_P)`; a uniform bound is KLS itself.
- Smooth local focal estimates plus isotropy and thin shell stop at
  `n^{1/4}`.  Disconnected cylinders realize all those local constraints.
  Connected Gaussian fans show that many balanced orientations may place all
  transitions on medial endpoints, although their ray scale remains constant.
- Boundary rotation/affine second-variation identities are homogeneous in the
  normalized surface measure and contain no Cheeger height.  Projected BV has
  a conditional-score term of order `n` on a regular-simplex slice.

## Extremal-sequence constraint set

A hypothetical bad sequence may be assumed isotropic, non-product up to the
tensorization constants, and witnessed by true/near T3 extremizers.  Its
witnesses are nonlinear, nonradial, not low-dimensional ridges, and induce
balanced long tangent ray families.  Any claimed contradiction must use at
least one property absent from the labelled-cylinder countermodel: connected
global log-concavity across orientation bundles, stationarity of a genuinely
near-optimal boundary, or spectral/eigenfunction structure.

## Next audit gates

No approximation, stochastic, or clean-room audit begins until one of the
three surviving mechanisms supplies a self-contained dimension-free central
lemma.  Any such lemma will first be tested against the cube, simplex,
`l_1` ball, product exponentials, Gaussian/Laplace tail decorations, the
Gaussian fan, and the disconnected-cylinder algebraic countermodel.
