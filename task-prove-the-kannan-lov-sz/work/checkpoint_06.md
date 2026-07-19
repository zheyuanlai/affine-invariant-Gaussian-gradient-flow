# Checkpoint 06: multilevel extremality, soft curvature, and the locked-phase obstruction

## 1. Exact results added since Checkpoint 05

### 1.1 Phase capacity and its sharp diffuse obstruction

For the natural tilt law at precision `t=tau/K`, product Poincare and
Buser--Ledoux give

\[
 C_P(\nu_t)\le t(1+\tau),\qquad
 \psi_{\nu_t}\ge {c\over\sqrt{t(1+\tau)}}.             \tag{1.1}
\]

If two central posterior phase packets have tilt-law mass at least `a`,
projective angular separation `gamma`, and centroid defect at most `e`, their
tilt-space distance is at least `c_delta gamma sqrt(t)`.  Expansion of
`nu_t` therefore forces intervening mass

\[
 c_\delta a\min\{1,\gamma/\sqrt{1+\tau}\}.             \tag{1.2}
\]

If all intervening states are central and have defect at least `e`, this
mass is charged to the exact Gaussian-profile derivative.  This is a valid,
noncircular capacity theorem for separated packets.

It does not turn effective rank into a fixed-mass packet pair.  A uniform
projective-spherical phase law has arbitrarily high effective rank while all
fixed-mass packets lie within angular distance `O(rank^{-1/2})`; the Gaussian
radial model realizes this.  Conditioning the tilt law on good states is
also invalid because the conditioned law need not be log-concave.  Finally,
the audited angular modulus asks for a polynomially small defect, whereas
the exact scale budget is

\[
 \int q(\tau){d\tau\over\tau}=O(1).                   \tag{1.3}
\]

Every such polynomial threshold is integrable against `d tau/tau`.  An
explicit scalar profile obeys the exact differential identity and stays
above the required threshold at every scale.  Thus scalar scale selection,
even on additive windows, is formally blocked.

### 1.2 Covariant Hessian null fields and phase-cell countermodels

For a noncritical posterior level-set function `z`, put

\[
 u={\nabla z\over|\nabla z|},\quad P=I-uu^T,\quad H=\nabla^2z,
\]

and define

\[
 |\mathcal C[z]|^2=
 \left|{PHu\over|\nabla z|}\right|^2+
 \left\|{PHP\over|\nabla z|}
 -{\operatorname {tr}(PHP)\over(n-1)|\nabla z|}P\right\|_{HS}^2.
                                                               \tag{1.4}
\]

This tensor is invariant under every monotone reparametrization of `z`.  In
dimension at least three, `mathcal C[z]=0` on a connected noncritical region
if and only if its level sets are parallel hyperplanes or concentric spheres.
The proof is an exact Codazzi--Riccati calculation.  In dimension two the
statement is false: every signed-distance field has zero tensor.

Small covariant energy alone has no global inverse.  Polyhedral facets and
balanced Gaussian/product-exponential maximum cuts have many flat phase
cells separated through small-amplitude interfaces.  They can have diffuse
high-rank directions and vanishing normalized covariant energy without a
positive-mass concurrent packet.  Near-global Cheeger extremality, not local
second-order geometry, is the missing hypothesis.

### 1.3 Exact angular occupation calculus

The heat-flow Bernstein identity and the audited modulus yield a
dimension-free weighted occupation estimate for low-defect angular states,
and an exact dictionary between heat spatial directions and localization
phases.  The scalar profile admits a formal zero-angular allocation for every
compactly perturbed admissible counterprofile.  More strongly, the only
local scalar weights invariant under every compact profile perturbation are
endpoint/null-Lagrangian weights proportional to `1-d`, where `d` is the
profile defect.  Hence no local scalar reweighting creates a positive angular
seed.  A phase rank, incidence, or extremality input is indispensable.

### 1.4 Gradient multiplicity and the ridge no-go

For convolution `Y=X+sqrt(s)G`, every scalar test has the exact three-defect
decomposition

\[
 (K+s)E|\nabla\phi|^2-\operatorname {Var}\phi
 =\mathcal D_G+\mathcal D_\mu+K\mathcal D_J.          \tag{1.5}
\]

The Gaussian-noise defect has Hermite stability; the original-measure
Poincare defect has no factor-rigidity theorem.  Exact high-rank equality of
all wedge components is impossible, but quantitative near equality does not
force a genuine factor.  A radial Gaussian median ball exactly pays the
wedge lower bound through its angular Hessian, and irreducible radial quartic
perturbations preserve the near-saturation while destroying exact product
structure.  A valid conclusion must be approximate affine-or-radial, not
literal product splitting.

A separate ridge calculation gives a positive conditional result: if every
bad Lipschitz witness admits an `L1` ridge approximation of rank `O(A^2)`
with raw relative error below `1/2`, the known dimension-`k` KLS bound closes
a universal self-consistency inequality for its amplitude `A`.  This route
cannot be generic.  Tensor amplification shows that every fixed nonlinear
regression component survives while arbitrary rank-`O(A^2)` observations
miss a universal fraction of it.  Gradient PCA and ordinary random
projections therefore require an extremality-dependent rigidity theorem.

### 1.5 Corrected coarea cascade and a sharp multilevel theorem

All coarea identities have now passed a clean-room proof with the correct
perimeter convention.  `P_mu` inside coarea is relative relaxed weighted
`BV` perimeter on the affine hull; it is not setwise exterior Minkowski
content.  The two definitions have the same Cheeger infimum.

For a half-mass set `S`, a `[0,1]`-valued `L`-Lipschitz function `F` with
`EF=1/2`, and

\[
 U={1\over2}E|F-1_S|,
\]

the full family of thresholds gives the sharp inequality

\[
 \boxed{
 U\ge {1\over4}-{1-e^{-\kappa/2}\over2\kappa},
 \qquad \kappa={\psi_\mu\over L}.}                   \tag{1.6}
\]

The right side is `kappa/16+O(kappa^2)`, improving the optimized two-level
constant by a factor of two.  A clipped affine function under the symmetric
Laplace law attains equality, and every one of its nontrivial levels is an
exact Cheeger minimizer.  Boolean combinations measurable only through `F`
have additive boundary over disjoint levels and produce no angular cross
term.  Thus even exact nested minimality cannot improve (1.6) without
splicing different spatial patches.

The audit also corrected a model claim.  Gaussian maxima and inner cube
boxes are far from optimal, but the balanced inner box for `n` independent
one-sided exponentials has perimeter tending to `(log 2)/2`, below the
coordinate-cut perimeter `1/2`.  Its non-extremality has not been proved; it
must remain a phase/product stress test.

### 1.6 Multiway midpoint incidence

For bounded canonical packets `A_i` centered at `c_i` on the natural tilt
scale, define the Gaussian overlap weights

\[
 w_{ij}=\exp\left(-{(|c_i-c_j|+2R\sqrt t)^2\over8t}\right),
 \qquad \rho=M^{-2}\sum_{i,j}w_{ij}.                 \tag{1.7}
\]

There is an exact separator/reflection alternative: either a prescribed
high-defect set has tilt-law mass at least `a rho/4`, or on mass at least
`a rho/4` every midpoint supports at least `a rho M/4` approximate reflection
pairs, with explicit center and distance errors.  A bounded-overlap version
gives `nu_t(H)>=a rho-L/M`.

The unrestricted inverse is false.  Null-sphere inflation destroys
noncanonical Minkowski multiplicity, and even canonical Gaussian winner
cones satisfy `A_i+A_j=R^M` for every distinct pair while their orthogonal
labels have neither affine nor radial organization.  Canonical perimeter
representatives repair the null-set problem but not the unbounded-cell
problem.  Coarea extremality excludes Gaussian maxima through an order-`p`
deficit, but it supplies no pairwise multiplicity cross term.  A closing
inequality would have to lower-bound coarea deficit by an actual weighted
separator/incidence charge.

### 1.7 The soft longitudinal driver

Let, under a current posterior,

\[
 g=\mu_t(S),\quad G=g(1-g),\quad v=\operatorname {Cov}_t(1_S,X),
 \quad k={|v|^2\over G},\quad a={1\over1+k},
\]

and use

\[
 \Gamma=P_{v^\perp}+aI.                               \tag{1.8}
\]

A nonsingular regularization at `v=0`, stopped well-posedness, all moment
SDEs, and removal of stops have been written out.  With

\[
 \beta=1-a,\quad Z_T=\int_0^T a_tdt,\quad
 Q_T=\int_0^T\Gamma_tdt,
\]

one has exactly

\[
 d\langle g\rangle=G\beta dt,qquad
 dH(g)=dM-\tfrac12\beta dt,qquad Q_T\succeq Z_TI.    \tag{1.9}
\]

Consequently

\[
 \boxed{E\lambda_{\min}(Q_T)\ge T-2H(g_0).}          \tag{1.10}
\]

This is a dimension-free unweighted curvature theorem.  It does not transfer
perimeter, because curvature may arrive after the label has polarized.  In
the entropy clock `tau=int beta`, `g` is exactly the neutral Wright--Fisher
diffusion and `Z=t-tau`; an explicit bang--bang clock makes every scalar
mass-weighted curvature lower bound vanish.

There is an exact post-seed Lyapunov.  If `q_t=lambda_min(Q_t)>0`, then

\[
 \boxed{I(g_t)\sqrt{q_t}\text{ is a stopped local submartingale}.} \tag{1.11}
\]

Indeed `dot q>=a`, the strong-convexity centroid inequality gives
`q|v|^2<=I(g)^2`, and the finite-variation drift is at least

\[
 {a\{I(g)^2-q|v|^2\}\over2I(g)\sqrt q}\ge0.          \tag{1.12}
\]

Thus preservation/amplification of a positive seed is solved; uniform seed
creation is not.  If `q_1(Q_T)` is small, all other eigenvalues are at least
`T/2`, the active direction is projectively locked to one line, and almost
all time is spent with `k` large.  Determinant and active-variance Lyapunovs
do not seed the exceptional direction.  The surviving case is precisely an
initially selected, high-variance one-dimensional phase.

### 1.8 Audited fixed-scale Fisher seed

The endpoint estimate in the heat-flow profile can be made at one prescribed
scale, rather than at a scale selected by a logarithmic pigeonhole.  Let
`K=C_P(mu)`, let `S` be a balanced near-Cheeger set, and put

\[
             \alpha=10^{-10},\qquad s=\alpha K .       \tag{1.13}
\]

After choosing the relative near-minimizer error `beta<=10^{-5}`, the exact
uncertainty and Bernstein identities give a central good set `G` for the
observation law `q_s` on which

\[
\begin{gathered}
 q_s(G)\ge1.25\,10^{-6},\qquad
 1-s|\nabla z|^2\le \tau<0.098,\\
 \operatorname {tr}R_G\ge {10^{-5}\over8\pi},\qquad
 \operatorname {rank}_{\rm eff}(R_G)
       \ge {10^{-15}K\over8\pi}.                     \tag{1.14}
\end{gathered}
\]

Here

\[
 R_G=\int_Gq_s(y){v(y)v(y)^T\over g(y)(1-g(y))}\,dy,
 \quad v=\operatorname {Cov}(1_S,X\mid Y=y),          \tag{1.15}
\]

and the central cutoff is the fixed number
`delta=Phi(-sqrt(2 log(3.2*10^6)))`.  The clean-room audit checked the
scale-sensitive identity

\[
 1-{|v|\over\sqrt{s}\,I(g)}=1-\sqrt e,
 \qquad e=s|\nabla z|^2,                              \tag{1.16}
\]

the zero-defect case, the `v=0` convention, every Fisher trace/rank constant,
and the exact angular derivative

\[
 P_u\nabla F={P_uD\over s^{3/2}\sqrt{g(1-g)}}.         \tag{1.17}
\]

The seed itself is therefore complete and dimension free.  Its natural
angular closure is not: the ratio between the available angular-stability
upper scale and the rank-sensitive wedge lower scale is

\[
 {C_\delta\Omega_\delta(O(\sqrt\alpha))\over\alpha^{3/2}},
                                                               \tag{1.18}
\]

which worsens as `alpha` decreases.  Longitudinal and bad-state derivative
energies remain outside the scalar coarea charge.  The symmetric Laplace
median halfline proves this is substantive: its coarea deficit is zero at
every heat time while its longitudinal Fisher derivative energy is positive.
The affine/radial/phase-charge alternative in the fixed-scale report is
therefore recorded only as a precisely quantified sufficient lemma, not as a
proved statement.

## 2. Newly blocked statements

1. **Scalar profile occupation.**  Compact scalar perturbations and the
   null-Lagrangian classification rule out a local scalar seed.
2. **Threshold-only extremality.**  The symmetric-Laplace equality model
   saturates every level; no Boolean combination of nested levels compares
   normals.
3. **Unweighted soft curvature.**  Equation (1.10) is insufficient for
   perimeter; the Wright--Fisher clock gives a sharp scalar countermodel.
4. **Generic low-rank ridge capture.**  Tensor amplification preserves a
   nonlinear component invisible to every fixed rank-`O(A^2)` observation.
5. **Effective rank to fixed packets.**  Diffuse projective-spherical laws
   and radial Gaussian cuts refute the implication.
6. **Unrestricted midpoint inverse.**  Canonical winner cones have maximal
   Minkowski overlap with no same-ray reflection or concurrency.
7. **Covariant-Hessian globalization.**  Polyhedral phase cells and the
   two-dimensional signed-distance exception require interface/incidence
   information.
8. **Instantaneous localization-rate arguments.**  A large current angular
   or covariance rate does not control finite-time displacement; all such
   claims require a stopped martingale estimate.
9. **Unweighted random-fiber energy at scale `1/n`.**  For the isotropic
   cube, regular simplex, and shifted-exponential product, the mean
   conditional variance of a Haar line is `Theta(1/n)`, and for a linear
   test the covariance-weighted directional energy is `Theta(1/n^2)` times
   its variance.  Thus neither a `c/n` lower bound for that energy nor the
   proposed one-step fiber recurrence survives even with zero Hessian.  Any
   repair must normalize the fiber pointwise and then control its correlation
   with the test function; doing so by the unknown global spectral constant
   is circular.

## 3. Active incompatible mechanisms

1. **Fixed-scale Fisher extremality.**  At `s=alpha K`, derive central good
   Fisher mass and high effective rank with all `alpha` powers explicit, then
   combine the wedge lower bound with the coarea near-minimizer property.
   The required conclusion must allow affine, radial, and genuine product
   branches.
2. **Soft-driver winner locking.**  Use (1.11) and (1.12) to reduce failure of
   a mass-weighted seed to an initially selected locked phase.  The missing
   theorem must show that a locked phase starting from isotropy is one-
   dimensional/product-like, or that its formation supplies transverse
   curvature.  A current-state variance estimate is not enough.
3. **Locked-fiber balance.**  In the hard mass-preserving transverse flow,
   prove quantitatively that a direction locked to a deterministic line
   produces approximately half-mass conditional one-dimensional fibers.
   Integrating their sharp one-dimensional log-concave Cheeger bounds would
   give a constant perimeter because the original conditional variances have
   mean at most one.  The unresolved historical issue is early random winner
   selection: a path-dependent line may select a posterior variance that is
   not controlled by the original covariance.
4. **Extremal spatial incidence.**  Seek a finite patch-splicing or
   reflection inequality which lower-bounds coarea deficit by phase
   multiplicity.  It must be null-invariant and must distinguish a long bad
   phase from bounded-scale exponential/product winner cells.

## 4. Audit ledger

* The general angular theorem has two independent audits and retains the
  affine-hull convention.
* The coarea cascade has an independent clean-room proof.  Setwise exterior
  Minkowski content has been removed from exact coarea statements.
* The soft driver includes a smooth `v=0` regularization, compact parameter
  stops, martingale integrability, perimeter-transfer orientation, and
  approximation to arbitrary log-concave laws.
* The covariance-transverse localization audit corrected an earlier
  instantaneous-rate overclaim: an order-`log^2 n` quadratic-variation rate
  alone says nothing about finite-time winner selection.
* All new model claims were run on Gaussian halfspaces, Gaussian median
  balls, cubes, one-sided and two-sided exponential products, maximum cuts,
  simplex vertex tilts, and irreducible radial quartic perturbations.
* There is still no candidate dimension-free proof.  No item in Sections 1--3
  is being treated as a KLS seed unless its mass weighting and global
  transfer have been proved.

## 5. Next checkpoint target

The next checkpoint will accept a route only if it supplies one of the
following genuinely new statements with a full audit:

1. a mass-weighted creation theorem for the soft driver which excludes the
   locked high-variance initial layer using spatial, not scalar, information;
2. a fixed-scale extremality inequality coupling high Fisher phase rank to a
   positive coarea/finite-competitor charge, with the exponential-product
   winner model classified as a bounded/product branch; or
3. a locked-line theorem converting pathwise transverse localization into
   approximately balanced original fibers, together with a historical bound
   preventing selection of a large-variance random line.

Any conclusion whose missing step is merely renamed phase synchronization,
weighted expansion, or exceptional covariance control remains blocked.
