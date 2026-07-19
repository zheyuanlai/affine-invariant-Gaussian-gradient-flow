# KLS proof-search registry

## Scope

Goal: prove a universal Cheeger lower bound for every non-atomic-support
log-concave probability, including lower-dimensional support. Every proposed
lemma is tracked with all dependence on dimension and auxiliary parameters.

## Active mechanism families

1. **Lipschitz witnesses / extremal structure (T3).** Seek an exact
   decomposition or structural exclusion for a 1-Lipschitz function with
   unbounded first centered moment. Independent agent active.
2. **Dual elliptic / averaged transport (T2/T5).** First round completed. It
   produced exact no-go identities and a hierarchical conditional-transport
   certificate, but no general certificate bound. The slot has been reallocated.
3. **Affine and degenerate-support reductions (T1/T2).** Root is proving the
   exact transformation and approximation statements required by any central
   proof.
4. **Localization / heat-flow controls.** Not yet active; capped so that this
   family does not dominate the search. Root is stress-testing a new controlled
   variant with sequential rank-one driver
   `C_t^2 = n u_t tensor u_t`, where `u_t` is chosen from the top covariance
   eigenspace. The intended mechanism is feedback damping of the current top
   mode without the maximum of `n` simultaneous Brownian modes.
5. **Geometric minimizers / needles / angular decomposition.** Independent
   agent active; tasked with a dimension-free half-mass distance lemma or a
   counterexample to each proposed intermediate claim.
6. **Weighted boundary stability.** Independent agent active. The first exact
   variation calculation shows a structural no-go: a smooth Cheeger-ratio
   minimizer of mass strictly below `1/2` would have nonnegative second
   variation for every normal speed, while summing the translation speeds
   gives `-int_Sigma Hess(V)(N,N)`. Hence positive-Hessian regularization
   forces minimizers to the half-mass cusp (or the infimum to escape); local
   stability alone does not control the Cheeger height there.
7. **Covariance-normalized localization.** Active. With
   `C_t=A_t^{-1/2}`, a set-mass martingale has quadratic-variation density
   at most `g_t(1-g_t)` by covariance Cauchy, while
   `dA_t=(third-moment martingale)-A_t dt` and the accumulated Hessian is
   `Q_T=int_0^T A_s^{-1}ds`. The matrix inequality
   `Q_T^{-1} <= T^{-2} int_0^T A_s ds` is exact. The load-bearing issue is an
   averaged endpoint-energy estimate retaining the correlation between the
   random covariance and the random density; replacing that correlation by
   an operator norm is forbidden.
8. **Cut transport / eikonal rays.** Active. There is an exact identity
   `D_1(mu)=2 sup_E p(1-p) W_1(mu_E,mu_Ec)`. An extremizing cut has a
   Kantorovich potential with unit gradient almost everywhere and a
   nonbranching balanced-ray decomposition. Bad near-extremizers therefore
   require a positive-mass family of long, almost tangent rays; a new
   focal-curvature/Jacobian estimate is being sought.
9. **Boundary dimension descent.** A quantitative alignment-to-subgraph
   repair is proved, but projected BV contains an unavoidable conditional
   score/quantile-velocity term. A regular-simplex slice makes that term
   order `n` even for an exactly flat halfspace; small-perimeter stationarity
   must be used before any induction can work.
10. **Spectral rigidity.** A fresh route is testing whether the Bochner
    estimate for a hypothetical small first eigenvalue forces an almost
    Euclidean product factor or contradicts thin shell, without applying a
    hidden Poincare inequality to gradient components.

## Blocked claims

- A dimension-free bound on `E int_0^c ||A_t||_op dt` for stochastic
  localization is conjecture-strength unless obtained from a genuinely new
  control or functional.
- For the sequential rank-one control, the unproved load-bearing claims are:
  dimension-free control of the largest covariance eigenvalue under adaptive
  eigenvector selection, sufficient coverage of all directions so that the
  accumulated quadratic potential is uniformly positive, and dimension-free
  quadratic variation of a balanced set's mass. None may be inferred from the
  scalar third-moment bound alone.
- Uniform Poincare for arbitrary isotropic log-concave measures is exactly the
  target and cannot be used as a compactness or regularization lemma.
- Thin shell alone does not control arbitrary Lipschitz functions; the known
  implication loses a factor growing with dimension.
- Convexifying a half-mass distance witness is not a free reduction. Already
  for the isotropic interval, the union of two outer intervals has half the
  mass and positive average distance, while its convex hull is the whole
  support. An additive universal convexification error, once convex witnesses
  are controlled, is quantitatively the full first-moment KLS target.
- Random-line resampling has no universal `c/n` spectral gap. For the isotropic
  cube its linear-mode gap is `Theta(n^{-2})`, whereas for the isotropic ball
  it is `Theta(n^{-1})`. The corresponding conditional one-dimensional energy
  shrinks at the same rate, so only a joint operator comparison could exploit
  cancellation; separate gap and fiber-Poincare estimates cannot close KLS.
- The exact algebraic source of the logarithmic loss is now recorded. Thin
  shell for every marginal gives `Var(X^T P X)<=C rank(P)` for every
  orthogonal projection `P`. This does not imply the Hilbert--Schmidt
  quadratic-form bound. The positive quadratic form
  `Q(H)=(Tr(BH))^2`, with eigenvalues `lambda_i(B)=i^{-1/2}`, satisfies
  `Q(P)<=4 rank(P)` for every projection but has operator norm
  `||B||_HS^2~log n`. Any proof closing the projection-to-general-quadratic
  step must use a property beyond the collection of thin-shell marginal
  inequalities.
- A transport map whose average differential norm is bounded only through the
  target spectral gap is circular.
- Brascamp--Lieb with `(Hess V)^{-1}` gives no uniform estimate on flat regions
  and cannot handle uniform measures on convex bodies without a new term.
- The best constant in
  `int_0^infinity E(|grad P_t f|^2)dt <= K E(|grad f|^2)` is exactly
  `C_P/2`; using this as a semigroup estimate is circular.
- The norm of the Witten inverse `(A^(1))^{-1}` restricted to exact gradients
  is exactly `C_P`; a Helffer--Sjoestrand proof must add genuinely new input.
- Source-global bounds such as `E[DT DT^T] <= C Cov(mu)` do not localize the
  target Dirichlet energy. Gaussian-to-Laplace transport gives an explicit
  counterexample with an unbounded target-conditioned weight.
- Boundary second variation is homogeneous after normalizing the surface
  measure and does not contain the Cheeger height. Quantitative normal
  alignment cannot be upgraded to closeness to a halfspace without a
  transverse BV estimate; an explicit Gaussian-product construction reduces
  that upgrade to the lower-dimensional Cheeger constant.
- Scalar likelihood weights cannot repair normalized-localization perimeter
  transfer. A Laplace-tail decoration leaves the transverse posterior
  statistics `(g_t,r_t)` unchanged while multiplying an irrelevant horizontal
  face by `sqrt(L)`; the decorated cut remains asymptotically near-optimal.
- The face-specific averaged `L^2` calibration estimate is exactly KLS in
  operator form. If `K_T` is posterior resampling and `R=I-K_T`, then
  `exp(-T)sqrt(C_P) <= ||L^{-1/2}R|| <= sqrt(C_P)`. The posterior
  decomposition supplies no weakening by itself.
- Local eikonal/focal data plus isotropy and thin shell do not suffice. A
  disconnected cylinder mixture realizes long balanced flat rays, while a
  connected Gaussian fan realizes arbitrarily many exact balanced
  orientations with transition curvature concentrated at medial endpoints.
  The former violates global log-concavity; the latter has only constant ray
  scale. A closing theorem must couple scale to the full global, including
  singular, variation of the normal congruence.

## Extremal-sequence constraints (provisional)

- By Milman equivalence, a counterexample sequence yields isotropic measures
  and 1-Lipschitz witnesses with diverging centered first moment.
- A near-maximizing first-moment witness may be chosen with
  `Var(f) <= C (E|f-Ef|)^2`: apply Milman's universal equivalence
  `C_P <= C D_1^2` and use `int|grad f|^2<=1`. Hence its upper and lower
  tails at a fixed fraction of its first-moment scale both have universal
  positive mass (Paley--Zygmund). Rare-tail pathologies remain possible for
  arbitrary witnesses but not for a near-maximizer.
- Linear witnesses are uniformly controlled exactly by isotropy.
- A radial witness of the form `g(|x|)` is controlled once the dimension-free
  thin-shell theorem is invoked, since `g` is 1-Lipschitz on the radius.
- Product decompositions cannot create a smaller Cheeger scale than their
  worst factor, up to the universal constants in Bobkov--Houdre tensorization.
- These observations do not yet control functions that are nonlinear and
  nonradial; treating them as though they did is a forbidden gap.
- A cut cannot be arbitrarily well predicted by all linear functions: for a
  universal `eta>0`,
  `Cov(1_S,X)^T A^{-1}Cov(1_S,X) <= (1-eta)p(1-p)`.
  The proof is a two-point quantization lower bound for a one-dimensional
  isotropic log-concave marginal. The isotropic interval attains ratio
  `3/4`, so this slack alone is too small to beat square-root curvature.

## Audit ledger

- Workspace audit: no pre-existing files other than `work/` and `outputs/`.
- Git: current task directory is not a repository.
- Candidate proof: none yet.
- Dual/transport round: a valid general variance-decomposition certificate was
  proved for hierarchical conditional transports. Its unresolved hypothesis is
  a dimension-free pointwise bound on the conditional velocity matrix; this is
  not presently known and must not be quoted as a conclusion.

## Checkpoint 01 update (supersedes the activity labels above)

### Active

1. **Global ray compatibility / additive overlap.**  The two-sided `r`-core
   has a `1/r`-Lipschitz Gauss map even across nonsmooth charts.  Long-ray
   direction rank forces `Omega(s^2)` full absolute turning on every connected
   completion, and log-concavity forces `Omega(s^2)` midpoint multiplicity
   between separated direction packets.  The live question is an inverse
   theorem turning that multiplicity into either a parallel factor or an
   approximately concurrent/reflection structure.  Both outcomes are
   dimension-free: linear functions are controlled by isotropy, and
   `Var|X-z|<=C` for every center `z` follows from thin shell.
2. **Set-mass-preserving localization.**  The control
   `C_t=P_{v_t^perp}`, `v_t=Cov_t(1_S,X)`, preserves the mass of the witness
   exactly.  Its endpoint precision always has at least `n-1` eigenvalues at
   least `T/2`.  A balanced-needle lemma proves
   `psi>=c/(alpha^{-1}+Var<X,u>)^(1/2)` under curvature
   `Hess V>=alpha P_{u^perp}`.  The exact remaining random quantity is the
   posterior variance in the single path-dependent exceptional eigenvector.
   The product-exponential/symmetric-halfspace model is being used as an
   adversarial test for a hidden `log n` tail-selection loss.
3. **Polar/angular decomposition.**  A clean-room route is auditing whether
   polar disintegration has any usable geodesic curvature after exact radius
   normalization, with cubes, simplices, and ellipsoids as mandatory tests.

### Newly blocked

- **Natural spectral rigidity.**  If `R_mu` is the optimal constant in the
  centered-gradient Bochner/Reilly estimate, then exactly
  `C_P(mu)-1<=R_mu<=C_P(mu)`, including continuous spectral edges and convex
  Neumann domains.  Thus an almost-parallel-gradient lemma of this form is
  quantitatively KLS, not an independent input.  Nodal convexity and
  log-concavity of an eigenfunction pushforward also fail on explicit
  interval/square examples.
- **Local boundary stationarity.**  The exact aligned-branch descent is
  `h_(n-1)/8<=J_u+P sqrt(delta)`.  Quantile coordinates show that `J_u` is a
  zeroth-order connection, whereas CMC and the Jacobi form see only its
  derivative multiplied by physical slope; it cancels completely on flat
  pieces.  The generic Fisher bound is sublinear in `P` and is unstable under
  convex-body approximation.  A closing boundary argument must therefore be
  global and may be merged only with the ray/focal mechanism.
- **Smooth curvature alone on transport rays.**  The Gaussian fan has zero
  classical shape operator on every chart but carries order-`m` singular
  turning at its focal vertex.  All future curvature functionals must include
  the completed full-BV/turning charge.

### Added extremal constraint

- The thin-shell theorem controls every translated radial model, not just a
  radius about the barycenter:
  `sup_z Var_mu |X-z|<=C`.  Indeed
  `Var(|X-z|^2)<=C(n+|z|^2)` and division by
  `E|X-z|^2=n+|z|^2` gives the claim.  Hence an approximately concurrent ray
  family is no more capable of supporting a bad T3 witness than a parallel
  family.

## Checkpoint 02 update

The full smooth-plus-medial second variation is now exact on generic strata,
and the translation trace identity captures arbitrary focal sets.  Its
small-energy branch gives alignment, but an abstract exponentially large
complete-graph model saturates every current stability, covariance, curvature
and packet-entropy constraint.  The missing content is specifically
Euclidean/log-concave realizability of the medial graph.

The sharp equality list has been enlarged.  Complete cross-calibration
classifies the endpoint sets as spheres in orthogonal affine subspaces.  This
Clifford branch is dimension-free by translated thin shell on the two
orthogonal marginals plus isotropy in the separating linear coordinate.
Accordingly, the active target is a robust trichotomy: aligned, radial or
orthogonal-radial structure, unless a finite signed-distance competitor
improves the global T3 objective.

The ordinary stochastic bootstrap, the polar chart gluing route, and the
abstract junction-energy route are blocked at conjecture-strength statements
and will be reopened only if they provide a new coupling, finite competitor,
or positive-density bridge estimate.
# Checkpoint 03 update

See `checkpoint_03.md` for the full audit.  New exact assets are the
interior cross-chord slack, the null-invariant multi-Brenier entropy identity,
the all-amplitude Bregman/envelope inequality, the covariance-saturation
inverse theorem, and the self-convolution Hoeffding defect.  Covariance-only
Brenier synchronization, packet-count amplification of finite offsets,
John/inradius geometry before balance, and monotone first-moment convolution
are now blocked by explicit countermodels.  The active central lemma is the
weighted endpoint-incidence heat-bath rigidity theorem, with soft
mass-localization retained as an incompatible alternative.

## Checkpoint 04 update

See `checkpoint_04.md` for the complete constant and obstruction ledger.
The endpoint program now has a sharp projection-wise rank theorem and a
genuine smooth Clifford-cone survivor: positive family weight and vanishing
endpoint density do not force active endpoints.  The strengthened geometric
target is fixed-bulk coherence of the `sigma sqrt(k)` focal centers supplied
by rank-balanced curvature; Codazzi gives only leafwise centers and does not
yet globalize their projections.

Ordinary Gaussian-profile localization now has an exact two-deficit equality
theorem.  Near equality simultaneously forces a posterior threshold cut and
a nearly Gaussian-dilation active marginal.  The active question is whether
these scalar defects control transverse splitting and the rotation of the
active direction across overlapping posterior tilts.  Low-order Hermite
expansions, a universal bounded-degree informative seed, and a separated
random-block gap/conditional-weight recurrence are blocked by explicit
exponential, Gaussian-parity, ball, and shifted-exponential tests.

The current incompatible live mechanisms are: fixed-bulk focal-center
coherence; Gaussian-centroid temporal/transverse gluing; and non-Gaussian
convex-cell localization with exact balance.  None is presently being
treated as a proved KLS seed.

## Checkpoint 05 update

See `checkpoint_05.md` for the complete ledger.  The clean-room posterior
stability theorem has passed, and its composition with Gaussian-halfspace
pullback now proves a dimension-free general angular estimate for arbitrary
strongly log-concave posteriors:
`||P D||HS^2 <= C_delta t^-2 Omega_delta(epsilon)`, with
`Omega_delta(epsilon)=O_delta(epsilon^(1/6)log(e/epsilon)^(1/12))`.
This closes the formerly open *pointwise* nonlinear-set/nonlinear-map
intersection, but its weak modulus does not close the temporal scale
mismatch.

The Gaussian-channel/ray bridge now gives exact posterior ray alignment,
fixed original mass of long calibrated rays, and high-rank phase
orthogonality.  The exact synergy height and gradient-wedge inequality give
two possible global charges, but the former has a log-concave zero-average
countermodel and the latter is presently matched by the longitudinal
eikonal term.  Exact convex nullity forces a global Gaussian cylinder;
approximate nullity does not globalize without overlap, as shown by cusp and
polyhedral-cell countermodels.

The four active mechanisms are phase partition/capacity, gradient
multiplicity rigidity, Fisher-weighted focal incidence, and temporal
integration of angular stability.  Scalar baselines, ordinary
resampling/amplitude closure, projection-height synergy, local-cell
nullity, bounded bulk resolvents, bounded-degree seeds, and separated
slab/block recurrences are blocked.  No dimension-free KLS seed is yet being
treated as proved.

## Checkpoint 06 update

See `checkpoint_06.md` for the full ledger.  Three scalar closures are now
formally exhausted.  The angular occupation calculus admits compact
zero-angular counterprofiles; the sharp multilevel coarea bound is attained
by a clipped affine profile under symmetric Laplace, with every level exactly
Cheeger-optimal; and the soft longitudinal driver has an exact
Wright--Fisher clock which lets unweighted curvature arrive only after label
polarization.

The soft driver nevertheless gives two exact positive theorems:
`E lambda_min(Q_T)>=T-2H(g_0)`, and, after any positive seed,
`I(g_t)sqrt(lambda_min Q_t)` is a local submartingale.  Failure of a
mass-weighted seed forces one low-curvature eigenline, projective locking of
the posterior centroid direction, and normalized binary correlation of
order the reciprocal seed size.  All transverse endpoint directions already
have numerical curvature.

At the deterministic heat scale `s=alpha C_P(mu)`, a balanced near-Cheeger
set supplies central good Fisher mass at least `sqrt(alpha)/8`, Fisher trace
at least `sqrt(alpha)/(8pi)`, and effective rank at least
`alpha^(3/2) C_P(mu)/(8pi)`, with explicit fixed `alpha`.  This removes scale
pigeonholing as the source of the missing phase seed.  The wedge/angular
power comparison worsens as `alpha` decreases, and longitudinal/bad-state
energy is not a coarea charge.  An exact symmetric-Laplace halfline shows
that longitudinal energy can be positive while coarea deficit is identically
zero.

The active core is now a spatial alternative: a low-curvature path or a
high-rank fixed-scale phase family must be shown to be affine, radial, or a
bounded/product branch, unless it creates a physical finite-competitor
charge.  Plain midpoint multiplicity is blocked by canonical Gaussian winner
cones; generic ridge capture is blocked by tensor amplification; covariant
Hessian nullity is blocked by flat phase cells.  A clean-room low-mode random
fiber descent and an extremality-dependent fixed-scale phase charge are the
two current incompatible tests.  No mass-weighted seed or complete KLS proof
is being asserted.

The fixed-scale theorem has now passed a separate clean-room audit.  With the
fully frozen choice `alpha=10^-10`, its good phase set has observation mass at
least `1.25*10^-6`, Fisher trace at least `10^-5/(8pi)`, and effective rank at
least `10^-15 C_P(mu)/(8pi)`.  The correct profile defect is
`1-|v|/(sqrt(s)I(g))=1-sqrt(e)`; all zero-set and endpoint conventions were
checked independently.

The random-line low-mode test failed sharply.  On the isotropic cube and
regular simplex, the mean conditional variance of a Haar fiber is
`Theta(1/n)`, yet symmetry makes the covariance-weighted energy of every
linear mode `Theta(1/n^2)` times its variance.  The same scale holds for a
coordinate mode of the shifted-exponential product.  Thus the hoped-for
`c/n` fiber-energy lower bound is false even with identically zero Hessian.
The only live repair is a pointwise standardized fiber dynamics, whose
gradient/covariance correlation remains to be controlled without invoking
the target spectral gap.  A separate locked-fiber route is testing whether
hard mass preservation plus a deterministic exceptional line produces
balanced original fibers and hence a direct one-dimensional perimeter bound.

## Checkpoint 07 update

See `checkpoint_07.md` for the full constant and obstruction ledger.

The fixed-scale phase seed is now a literal physical coarea measure with a
globally analytic selector.  Conservatively, its trace exceeds `.005109p`
and its effective rank exceeds `19.54`; a core with selector at least
`1/2000` and alignment ratio at least `1/2` retains trace `.004489p` and
rank `17.16`.  Central restriction removes extreme-volume levels and leaves
trace `.0048495p`, rank `17.86`, and an explicit pointwise level-trace floor.
The exact same-level/between-level variance split therefore survives without
selector zeros, amplitude wells, or endpoint concentration.

The special selector has dimension-free binary Fisher and vector-direction
information bounds.  These do not close the argument.  An exact
heat-generated log-affine facet model has perfect eikonal alignment and zero
core Hessian charge while its rank is entirely cross-level; every possible
charge lies at ridges, ends, or global phase transitions.  The exact
level-transport equation isolates the remaining terms as lapse-weighted
within-trajectory rotation, orientation-dependent geometric/selector
reweighting, and a focal/contact residual.

The physical finite-splice functional is at most `6.02*10^-5p`.  A proved
ridge-graph spectral-gap inequality yields more than `1.17*10^-4p` whenever
the admissible ridge graph has gap at least `1/4`, so the numerical expander
branch is closed.  Common-slope flat translations are an exact no-go.  In
the low-incidence branch, genuine perimeter-additive log-affine components
obey

`P(A)-psi mu(A)=sum_i a_i(1-psi/kappa_i)`.

Hence every significant separated component of a near-Cheeger set must have
normal slope close to the common Cheeger slope; at heat scale `s=alpha K`
such synchronized slopes cannot create fixed softmax rotation.  The active
geometric cut lemma must turn small ridge/end conductance into actual
components with the same error budget.

For exact fully regular isoperimetric leaves, the global killed-normal-tube
theorem has passed an independent clean-room audit.  Profile near-linearity
over a fixed multiplicative volume interval implies tube length `Theta(1/p)`
and normalized flux loss `O(delta p)`, with first support contact, focal
collapse, cut collision, and transported curvature all included exactly.
Long survivors give the covariance lower bound

`Cov(mu) >= exp(-|lambda|T-h) T^3/12 int_G N tensor N`.

This closes coherent normals.  The high-rank product-exponential tail model
has collision/ridge loss only quadratic in its small volume and calibrates
the remaining almost-product inverse.  The tube theorem cannot yet be
applied to the heat packet: the latter lies on integrated near-minimizers,
not exact CMC/profile leaves.  A rank-preserving Ekeland/quasiminimal transfer
is a separate active target.

The three active teams are now: rank-preserving near-minimizer-to-CMC
transfer; low-ridge component cutting and slope synchronization; and
high-rank killed-tube collision-versus-product classification.  Gaussian
selector regularity, scalar Bernstein propagation, local quadratic
curvature, generic random fibers, and an ambient capacity estimate at scale
`K` are blocked as standalone mechanisms.

The first low-ridge component lift failed audit in its naive form.  For a
boundary phase partition the exact BV cost is twice the minimal interior
filling of the trace cut; exterior codimension-two ridge conductance does not
bound that filling.  Rectangles and balls are sharp countermodels.  The
registry now treats large filling as a separate coherence/concurrency branch,
to be attacked by max-flow/curve decomposition, long-ray covariance, and
translated thin shell.  Cheap filling retains the component-slope inverse;
large filling may not be silently called a component decomposition.

### New incompatible family: mass-preserving matrix localization

`mass_preserving_matrix_localization.md` records an exact controlled
localization law for one fixed witness set.  If `A_t` is the current
covariance and `v_t=Cov_{mu_t}(1_S,X)`, choose

`D_t=A_t^(-1/2)(I-w_t w_t^T/|w_t|^2)A_t^(-1/2)`,
`w_t=A_t^(-1/2)v_t`.

Then `D_t v_t=0`, so `mu_t(S)` is pathwise constant, while the covariance
drift is

`-A_t+v_t v_t^T/(v_t^T A_t^(-1)v_t)`

and damps all but one whitened mode.  The fixed-set weighted perimeter is a
nonnegative local martingale.  A terminal needle reduction would therefore
give `P_mu(S)>=c min(mu(S),1-mu(S)) E(1/sigma_omega)`.

The hard feedback is not well posed at `v_t=0`.  In one dimension a
Gaussian with a symmetric interval turns the convention “full control at
zero, zero control off zero” into an occupation-time contradiction, so even
a continuous weak solution need not exist.  The only valid current
formulation is a mesh-relaxed rank-`n-1` control; every terminal estimate
must hold uniformly over its relaxed limits.

The algebraic covariance disintegration supplies only
`E[sigma_omega^2 u_omega tensor u_omega]<=I`, whose trace loses `n`.
The unproved terminal estimate `E(1/sigma_omega)>=c` is now the sole
geometric obstruction in this family.  Long terminal needles of isotropically
spread directions must have nearly coincident barycenters, identifying the
same concurrent/orthogonal-radial branch as large interior filling.  No
terminal variance estimate is being assumed.

### Near-minimizer and large-filling audits after Checkpoint 07

The scalar proximal problem
`min_{mu(B)=mu(A)} P(B)+Lambda mu(A triangle B)` is dimension free:
`mu(A triangle B)<=delta/Lambda`, it is globally `Lambda`-quasiminimal,
and its regular weighted mean curvature lies in a common interval of width
`2Lambda`.  A killed-tube formula therefore loses only `exp(Lambda T)`.
This does not preserve the ranked normal packet.  Scalar divergence
calibration aligns the old normal with its field but cannot transport the
quadratic projector to the new boundary; an exact rotating divergence-free
field is the obstruction.  Nuclear-matrix fidelity really preserves rank,
but its Euler equation is anisotropic and contains `kappa grad omega`; the
fixed heat-scale derivative audit leaves a nonperturbative loss of order
`10^4`.  An isotropic asymmetric-Laplace example also shows that a nearby
exact profile minimizer may lose essentially all of a positive analytic
selected packet.  Quantitative stationarity of the physical level, not
exact replacement, is mandatory.

The physical filling error is exactly `R_cut=2 Fill`.  Continuous
max-flow/min-cut duality and corrected-flux Pythagoras are proved, and an
expanding calibrated filling graph is quantitatively projectively coherent.
Large filling supplies a lower bound on quadratic capacity, not the upper
bound needed for a cheap interpolation.  A constant-radius tube about an
arbitrarily turning embedded curve has an exact unit divergence-free
max-flow and a minimum cross-sectional cut, disproving flow-only
parallel/concurrent/radial rigidity.  The corrected live alternative is
bevel expansion, cheap physical min-cut, or a large interior flow whose
turning/Plateau excess must still be charged by near-minimality.

The terminal-needle audit proves that, conditional on an actual terminal
line disintegration, a fixed direction cap closes by covariance, exact
common barycenters are impossible in dimension at least two, and the
thin-shell theorem bounds the multiplicity of every fixed eigenvalue
threshold of `K=E[sigma^2 u tensor u]`.  Quantitatively
`N_K(alpha)<=C alpha^(-9/2)` and `E sigma^2<=C n^(7/9)`.  This does not close
the multiscale small-eigenvalue/dispersed-barycenter regime.  The exact
product-exponential maximum model permits coordinate-ray scale `log n`, so
no pointwise terminal-variance bound may be used.

`equimeasurable_matrix_replacement.md` starts a revised transfer.  The
selected central normal law is a soft thinning of the full central coarea
normal law; hence the latter has unweighted projector variance greater than
`.0045779`.  Joint minimization over functions equimeasurable with the
clipped heat function preserves the unweighted coarea normal matrix within
`Delta_F/kappa` and produces only a constant anisotropy, eliminating the
spatial `grad omega` term.  At the frozen constants the retention estimate
fails numerically by a factor about forty-three, and a smaller-deficit
fixed-scale run must be re-audited.  More importantly, global nested
minimality has not yet been converted into a.e. unconstrained levelwise CMC
geometry; this is the active transfer question.

## Checkpoint 08 update

See `checkpoint_08.md` for the complete current ledger.  The
constant-anisotropy replacement has passed two independent numerical and
geometric audits.  The normalized matrix estimate is

`||Q'-Q||_* <= (Delta_F/T)(1+1/kappa)`.

With `kappa=10^-6`, `alpha=10^-28`, and `beta=10^-14`, the central trace is
larger than `.0032827p`; matrix retention, Wulff-direction conversion, and
the explicit tube deletion leave a fixed positive angular-variance packet.
The price is an enormous but fixed rank threshold, which is harmless for a
dimension-free conclusion.  The old fixed-`kappa` power mismatch is no
longer an obstruction.

The exact smooth Wulff calculation is now complete.  For
`z=D Phi(N)` and `B=D^2 Phi(N)|_(N^perp)`, the stopped ray Jacobian is

`j_x(t)=det(I+tBS) exp(-V(x+tz)+V(x))=exp(lambda t-D_x(t))`,

and covariance is controlled by the displacement matrix `z tensor z`, not
the Euclidean normal matrix.  A two-dimensional example refutes every
two-sided Loewner comparison; only projective/nuclear stability is valid.

A near-Cheeger CMC leaf no longer needs to be a global profile minimizer.
Concavity of its own positive-lapse foliation bounds its multiplier by
secants of the Cheeger tent.  After orienting the leaf or its complement, a
short volume excursion `gamma m(v)` has length `Theta(gamma/psi)` and loses
only `O(delta+gamma)` normalized Wulff flux.  Direct Cheeger-deficit
selection deletes only `D_co/epsilon` of Euclidean surface trace, and the
heat identity

`int |mu(F_0>r)-1/2| dr <= ||F_0-1_S||_1=2U(s)`

localizes the retained trace near half volume.  The full explicit short-
tube coefficient is positive and dimension free.

Two audits prevent this smooth theorem from being promoted prematurely.
First, a BV equimeasurable minimizer has only aggregate coarea stationarity.
An `SBV` interface satisfies the divided-difference law

`H_ij=(p_j-p_i)/(a_j-a_i)`,

so one threshold may contain several curvature multipliers.  Jump matrix
cannot be bounded by profile deficit: an indicator of an isoperimetric set
has zero deficit and entirely jump-supported matrix.  Second, high flat
rank plus zero Jacobi/contact charge and negligible normalized killed flux
does not imply a global product.  The cyclic support

`{x_i>=0, x_i+x_(i+1)>=a}`

has coordinate tail leaves with normal matrix `I/m` but an irreducible
facet-normal circuit.  Global swept-flux coverage of support/Hessian
couplings is indispensable.

The newest incompatible transfer minimizes the direct Cheeger deficit
itself:

`D_psi(G)+kappa||M(G)-M(F)||_*`,

where

`D_psi(G)=TV(G)-psi int min(mu(G>r),1-mu(G>r))dr >=0`.

The distribution term is the median absolute deviation.  Comparison with
the heat function formally retains both the deficit and the full normal
matrix, while the smooth Euler multiplier is prescribed in `[-psi,psi]`.
Every `SBV` jump divided difference is also in this interval.  If the
nonsmooth audit succeeds, a variable-multiplier short Wulff tube can bypass
both equimeasurability contacts and surrounding-foliation concavity.  The
simultaneous nuclear/median subgradient, Cantor part, anisotropic
quasiminimality, and singular Wulff coverage are under clean-room audit and
are not yet treated as proved.

The active families are now: direct-deficit nonsmooth Euler/coverage;
finite layered/subgraph regularization; balanced high-rank global
saturation; and, independently, mesh-relaxed mass-preserving localization.
Rank-only product rigidity, levelwise minimization, aggregate-BV-to-CMC,
flow-only filling rigidity, fixed-`kappa` Wulff transfer, and hard feedback
at `v=0` are blocked.

## Checkpoint 09 update: global completion for one interface

- **Direct variational family:** exact common median subgradient and constant anisotropy obtained; almost every level is a global forced minimizer. Matrix losses are explicitly dimension-free. The obstacle multiplier can carry a median jump, so common calibration does not yield a foliation.
- **Diffuse-band rigidity:** stable regular bands are flat and log-affine. This excludes genuinely diffuse curved exact foliations but not a single balanced interface.
- **Polyhedral/log-affine family:** closed dimension-freely for exact balanced polyhedral minimizers by bevel elimination followed by central-cell mass completion and a one-dimensional density bound.
- **Smooth geometric family:** local stability, Jacobi, two-sided tube, and normal-matrix data are insufficient; an isotropic radial exponential median sphere is the sharp countermodel.
- **Quantile/asymmetry family:** aggregate moment and total-variation penalties cannot stop collapse to one off-center level; product exponentials give exact matched counterexamples.
- **Random-fiber family:** unweighted incidence is false at the required scale on the cube. The weighted inverse-scale quantity remains live only if coupled to a new global fiber-completion theorem.
- **Localization/transport-ray family:** per-ray inequalities are available; the unresolved step is a nontrace spectral compatibility law forced by log-concavity and nonbranching.
- **Primary live target:** measurable global flow cells for a single balanced interface, with dimension-free bounded overlap or a collision-to-perimeter-saving alternative.
- **Blocked:** covariance-process control equivalent to KLS; local-only curvature closures; unlabeled aggregate penalties; unweighted fibers; naive smooth triangulation.

## Checkpoint 10 update: contact tensor and longitudinal turning

- **Tensor Minkowski family:** the exact free-boundary identity
  \(I-Q=HX+B\) retains a Frobenius signal of size \(\sqrt d\).  Closed
  interfaces with isotropic-scale surface moment close dimension-freely;
  a small-perimeter escape must be contact-tensor or boundary-tail dominated.
- **Stability/contact family:** translated-normal tests bound support
  curvature at contact and, with a global tube defect, yield a sharp
  reciprocal-curvature lower bound.  This remains infinitesimal.
- **Quartic audit:** zero contact curvature does not imply support ruling;
  a quartically flat convex graph is the explicit counterexample.  Global
  longitudinal turning or a finite chord comparison is mandatory.
- **Weighted-fiber family:** exact Crofton normalization is complete and
  all canonical models have the desired \(d^{-1/2}\) scale.  Majority
  completion plus marginal max-flow is exactly quasiminimality and gives no
  strict saving.
- **Transport-ray family:** cross-endpoint co-Lipschitz packing is exact.
  The fourth-moment spectral improvement reaches \(d^{-6/17}\) only and does
  not sum.  A cube checkerboard proves that regular Gauss Jacobians cannot
  replace transition-stratum accounting.
- **Local contact packets:** one packet of long rays has one approximate
  projective direction and a common contact-tensor kernel.  The unresolved
  escape is a diffuse collection of spatially separated packets.
- **New live routes:** longitudinal support turning, median signed-distance
  shape variation, T3 quotient second variation, and multi-packet
  aggregation.
- **Blocked:** pointwise curvature-to-ruling; one-directional marginal
  descent; nonsummable spectral powers; regular-chart Jacobian arguments;
  local-only tube and stability estimates.

## Checkpoint 11 update: random-line uncertainty and transport rank

- **Median signed-distance family:** the exact representation of `D1` by
  half-mass average distance is complete, including level atoms and affine
  support.  A smooth maximizer bisects every normal Voronoi cell.  The full
  second variation contains a nonpositive medial switching form; the radial
  exponential law has degree-one equality.  The remaining upper bound is T3
  itself.
- **Multi-cell chord family:** exact volume cancellation removes only the
  common CMC multiplier.  Stable smooth residuals cannot become negative by
  summation.  Only genuine first-order ridges fire; this branch is blocked
  for smooth survivors.
- **Aggressive localization family:** `D=A^-1 P A^-1` preserves the chosen
  set mass and gives covariance drift `-P`.  Its stopped Itô ledger passed an
  independent audit.  It reaches only the first covariance face; direction
  rotation, divergent terminal clock, singular continuation, and survivor
  variance remain unresolved.  No terminal needle is assumed.
- **Random-line family:** the exact target
  `E_theta int min(p,1-p)/sigma >= c/sqrt(d)` would close KLS through
  one-dimensional Cheeger and Crofton.  Gaussian, cube, simplex, radial
  exponential, and product-maximum models have the sharp scale.  Direct
  Radon high moments lose dimension; the target remains unproved.
- **Transport-density family:** raywise Beckmann integration gives
  `M=int tau u tensor u=int h X tensor u`, `||M||_HS<=1`, and
  `rank M >= (tr M)^2`.  Clean-room audit validated these statements and the
  Beckmann factor.  It corrected the Bochner boundary term to a signed flux;
  a Gaussian multiple-crossing set disproves positive-part equality.  A
  general nonsmooth support/medial curvature measure remains unproved, so
  curvature is used only under explicit smooth endpoint hypotheses.
- **Live incompatible routes:** transport-density global compatibility;
  random-line uncertainty; Boolean rigidity for the normalized line form;
  and a projected high-rank ray collision test.  Normal-cell log-concavity is now
  complete locally: `E|T|` and `sd(T)` are both comparable to the reciprocal
  median density, but the aggregate remains T3-equivalent and the radial
  exponential sphere saturates directional stability.
- **Polarization family:** set-only polarization increases the distance
  objective for an already reflection-invariant measure.  Joint
  density/set polarization fails mass balance and can destroy
  log-concavity; Steiner-rebalance-whitening has both monotonicity signs.
  This route is blocked for general measures.
- **Newly blocked:** smooth chord aggregation, quotient stability alone,
  aggressive localization past the first face, and direct high-moment Radon
  comparison.

## Checkpoint 12 update: defect curvature and singular survivors

- **Defect-subspace theorem:** if \(D^2V\succeq\kappa P_E\)
  distributionally, then
  \[
  C_P(\mu)\le\kappa^{-1}
       +2\{C_P((P_F)_\#\mu)+\kappa^{-1}\}.
  \]
  The rank-one specialization is
  \(C_P(\mu)\le96(\kappa^{-1}+\operatorname {Var}\langle X,u\rangle)\).
  A clean-room audit independently rederived reinforced Prekopa, the
  exact-one-form improved Brascamp--Lieb estimate, operator domains,
  smoothing, and affine supports.
- **Ordinary bootstrap no-go:** covariance-whitened localization yields the
  raw recurrence
  \(\mathcal K_n\le C(n/m+n\mathcal K_m)\), and water-filling yields
  \(\mathcal K_n\le C((n/m)^2+n\mathcal K_m)\).  Even the counterfactual
  recurrence \(\mathcal K_n\le C(n/m+\mathcal K_m)\) cannot improve
  \(\mathcal K_m\le C\log m\).  This family is blocked absent a genuinely
  new covariance input.
- **Bounded/soft projection flow:** hard Euclidean projection preserves a
  balanced set and produces curvature at least \(1/2\) outside one terminal
  line by time one, reducing the proof to an adaptive survivor inverse
  moment.  The hard rule is ill posed at zero signal.  The continuous soft
  rule
  \[
  C_\varepsilon=I-\frac{bb^T}{|b|^2+\varepsilon}
  \]
  has set-mass quadratic variation at most \(\varepsilon/4\) and preserves
  the same rank-one-defect curvature.  It repairs existence and mass
  survival but not future-selected survivor variance.  The uniform law on
  \(\{\pm\sqrt n e_i\}\), with the positive atoms as the label, realizes
  the exact feedback and makes the survivor inverse moment
  \(O(\sqrt{\log n/n})\).  This is a non-log-concave mechanism no-go.  A
  binary quantization lemma proves
  \(\mathbb E\operatorname {Var}(Z\mid Y)\ge
  c\operatorname {Var}Z\) for every one-dimensional log-concave \(Z\),
  excluding atomic collapse but not an adaptively selected long continuous
  needle.  The active log-concavity input is an orientation-convexification
  tail bound.
- **Displacement midpoint:** exact domination, entropy, covariance, and
  Cayley-strain identities are complete.  The zero-gap branch contains
  singular piecewise translations.  Orthogonal product corridors close:
  half plateau volume gives \(\prod_j(1-\delta_j)=1/2\),
  \(\sum_j\delta_j\le\log2\), and squared translation fluctuation at most
  \((\log2)^2\).  Arbitrary noncommuting corridors also close if they are
  sequential Minkowski extrusions of a common half-volume convex core:
  \(\sum_j|v_j|\le\sqrt{24}\) and the label energy is at most \(13\).
  Simplex and crosspolytope fans close separately.  The remaining
  local-to-global step is to turn arbitrary local wall prisms into a
  laminar family of large convex cores; the wall \(L^1\) budget alone does
  not provide this.
- **Singular eikonal completion:** the canonical ray measure has regular
  weighted curvature plus linear-BV medial and support endpoint charges,
  total \(2P\).  A quadratic singular Hessian measure does not exist.  A
  fixed-mass long-cell core has charge \(O(1/D)\), bounded basepoint
  covariance, and normal tensor stable rank at least \(cD^2\).
- **Normal-congruence tests:** the endpoint-free complete smooth branch is
  a hyperplane and has signed-distance moment at most \(\sqrt2\).  The
  isotropic cube parity fan is a realized high-rank transition model with
  \(Q=I/n\), \(J=\sqrt3/(n+1)\), \(P=n/(2\sqrt3)\), and
  \(2PJ=n/(n+1)\).  Thus rank-only rigidity is blocked; the active target
  is a positive global incidence inequality retaining medial/support walls.
  On a long-cell core, the exact reach calculation gives
  \(\int\|S_y\|_{HS}^2d\eta\le C/D^2\); it controls local shape but a
  direct chord summation loses transverse trace.  Stable rank cannot be
  discretized into polynomial-mass narrow direction caps, as the radial
  exponential model demonstrates.
- **Active incompatible routes:** softened projection plus a
  past-measurable anchor; convexification/incidence of stable-rank long
  normal packets; noncommuting mixed-volume control for zero-gap
  translations; and applications of the defect theorem outside ordinary
  trace localization.

## Checkpoint 13 update: phase action and tensor normal-cone barriers

- **Tensor spectral family:** two-sided full-matrix stationarity of a tensor
  square forces the mixed stress `-lambda f tensor f` to be quadratic and
  hence forces `f` to be linear.  The quantitative factorization is
  dimension free.  The new DC audit proves that every convex
  positive/negative buffer for this stress has normalized Hessian norm at
  least `sqrt(2(1-lambda^2))`; a positive KKT multiplier pays action at
  least `1-|E Xf|^4`.  Bochner smallness therefore does not erase the
  convex-potential normal cone.
- **Directional H-minus-one family:** spectral stability near a linear
  eigenfunction would follow from the operator estimate
  `E(tau_* tau_*^T) <= C I` for the minimal Stein kernel.  The proved theorem
  controls only its trace by `Cn`; scalar moment-map control of
  `a^T tau a` does not control `|tau a|^2`.  Treating the operator estimate as
  known is circular.
- **Gaussian-observation family:** posterior probit is 1-Lipschitz.  The
  exact monotone functional `P_E(t)=sqrt(t) E I(p_t)` runs from zero to the
  perimeter, while Boolean uncertainty defines a probability action measure
  on log-SNR.  The identities do not locate this action at a universal
  scale.  The correct low-SNR chi-square divergence cancels independent
  nuisance coordinates, but exponential median labels have factorial mixed
  moments and the Hermite expansion has zero radius of convergence.  All
  fixed-degree chaos closures are blocked.
- **Displacement family:** the corrected midpoint identity is
  `Cov(nu_1/2)=I-K/4`, with a universal spectral floor `I/4800` and
  `tr K/4 <= -log det Cov <= 1200 tr K`.  A realizable interval example has
  zero entropy and local strain at every interpolation time but nonzero
  covariance loss, so physical translation gaps must be charged.
- **Projected-localization family:** a second protected scalar constraint
  leaves an exact direction-rotation coefficient.  Killing the full
  protected-variance martingale coefficient leaves a signed third-moment
  drift, and a smooth log-concave product rules out every bounded-derivative
  correction by the natural scalar energies.  Rank-two control creates
  constant curvature outside at most four dimensions, but a static ball
  posterior shows that adaptive survivor covariance is still essential.
  Conversely, an exact subcritical ball calculation bounds any eigenvector
  variance by `6/(1-T)`, eliminating the proposed argmax-sign obstruction
  before the critical time.
- **Function-specific localization:** preserving the mean of a normalized
  first eigenfunction and applying the rank-one defect theorem gives the
  exact ledger
  \`1 <= (192/T) lambda + 96 E[R_T E_T]\`.  Hence a hypothetical small gap
  forces the energy-biased adaptive survivor variance to be
  \`Omega(1/lambda)\`.  The soft driver has the complementary exact tradeoff
  \`1-E m_T^2 <= lambda/(delta T)\` when it creates full curvature
  \`delta T I\`.  This implementation cannot preserve the signal and
  regularize its selected direction simultaneously at a universal cost.
- **Zero-gap displacement family:** exact data are a cyclic firmly
  nonexpansive binary martingale dilation.  Iteration has geometric lifetime
  and quadratic variation `2 E|b|^2`; covariance splits as `I-S` and `I+S`.
  The missing theorem is `tr S <= C`.  Local singular Hessian mass and the
  canonical firm extension both fail on an explicit one-dimensional
  physical-gap example.
- **Live incompatible routes:** function-specific localization for a fixed
  eigenfunction; Bochner/PDE rigidity after tensor near-linearity;
  nonperturbative Gaussian information with the exact output denominator;
  and a cyclic physical-gap/laminar-extrusion theorem.
- **Newly blocked:** DC-buffer absorption of tensor stress; the directional
  Stein operator estimate as a consequence of its trace version; finite
  Gaussian chaos at numerical SNR; finite scalar Lyapunov corrections for
  projected localization; function-specific hard/soft projection without a
  new adaptive-survivor theorem; and local-only mixed-volume charges.

## Checkpoint 14 update: eigenfunction graph rigidity and slice curvature

- **Near-linear spectral family:** for a normalized first eigenfunction,
  the affine residual \(r=f-(EXf)\cdot x\) has exact mean gradient
  \(-(1-\lambda)EXf\), Hessian norm at most \(\lambda\), and Bochner defect
  \([B_f-\lambda\operatorname{Var}(\nabla f)]+\mathcal R_V
   =\lambda^3|EXf|^2\).  The natural full-gradient interpolation is false
  on isotropic balls by a factor \(\sqrt n\).
- **Residual mean-gradient gate:** the surviving sufficient statement is
  \(|E\nabla g|\le C(\|g\|_2+\|D^2g\|_2)\) on the affine-orthogonal
  subspace.  Its Hessian-only form holds with constant \(12\) in one
  dimension and tensorizes exactly to product laws.  The quadratic
  restriction is the directional third-moment Frobenius bound
  \(\sup_u\|E[(u\cdot X)(XX^T-I)]\|_{HS}\le C\).  Projection thin shell
  yields only \(\sqrt{\log n}\); conditional-slice differentiation leaves
  signed score or moving-wall fluxes.
- **Dependent positive classes:** paired box wedges have
  \(\sup_u\|M_u\|_{HS}\le2\), are \(4\)-Lipschitz product images, and have
  \(C_P\le192\).  One-sided wedges have \(C_P<972\).  Dirichlet laws,
  symmetric-base cones, affine box slices, and affine Gaussian covariance
  pencils also have universal directional third-moment bounds.  Axial
  log-volume curvature charges the squared covariance slopes exactly.
- **Function-specific localization:** mean preservation plus the rank-one
  defect theorem gives
  \(1\le192\lambda/T+96E[R_T\mathcal E_T]\).  Small gap forces survivor
  variance \(\Omega(1/\lambda)\) under the energy-biased path law.  A soft
  full-curvature driver loses the signal according to
  \(1-Em_T^2\le\lambda/(\delta T)\).
- **Convex-half-set distance:** for \(\mu(A)=1/2\),
  \(E d(X,A)\ge(2\log2-1)/(4\mu^+(A))\), with the sharp scalar constant.
  A universal reverse distance bound is KLS-strength; naive hull or convex
  replacement fails even for one-dimensional isotropic Laplace.
- **Cayley-potential audit:** all balanced transports have an exact
  nonconvex \(1\)-smooth certificate.  A zero-deficit interval transport
  defeats every signed convex certificate with the required orientation,
  so entropy-controlled convexification is impossible.
- **Active incompatible routes:** eigenfunction-specific mean-gradient
  rigidity; noncommuting slice-curvature compensation; low-mode use of the
  parallel coupling; and a laminar/incidence theorem for cyclic binary
  dilations.
- **Newly blocked:** unrestricted eigenfunction gradient interpolation;
  thin-shell-only proof of the directional third tensor; naive
  one-dimensional conditional slicing; and every entropy-controlled
  signed-convex replacement of the Cayley certificate.

## Checkpoint 15 update: convolution amplification and Gaussian-fiber WFI

- **Forward spectral amplification:** normalized self-convolution sends
  gaps \(a\) and \(b\) to
  \(b(1-b)\le4(b-a)\), hence
  \(b\ge(\sqrt{9+16a}-3)/2\).  This is dimension-free but purely forward;
  it gives no dimension-free lower bound on the initial gap.
- **Reverse and slice-local obstructions:** conditional expectation is not
  an \(H^1\)-contraction even for the explicit interval model, and scalar
  entropy, inverse-Fisher, and fixed-noise MMSE budgets do not reverse the
  convolution step.  Pointwise post-noise curvature--Korn control also
  fails; only an integrated Stein-weighted statement survives.
- **Verified positive classes:** the weighted-Fisher estimate is closed for
  every pure translating slice family.  For smooth full-support Gaussian
  conditional fibers with \(R\succeq I\), \(ER\preceq\Lambda I\), and the
  stated mixed-isotropy condition, the exact noncommuting calculation gives
  \[
  C_s+B_s^0=\tfrac12\operatorname {tr}
  (R^{-1}R'R^{-1}R')\le\Phi'',\qquad
  \int\rho\tau^2(C_s+B_s)\le1+C_{\rm cent}\Lambda.
  \]
  The centroid proof uses the corrected universal bounds
  \(\tau\le400(1+|s|)\) and
  \(\tau|\Phi'|\le400(1+|s|)\); the earlier constant-one pointwise score
  claim remains unverified.  This theorem is smooth and post-unit-noise
  only.  No nonsmooth, distributional, or zero-noise extension is claimed.
- **Retraction:** `scalar_shape_acceleration_lemma.md` is retracted.  Its
  midpoint expansion cancels the proposed \(U_zm''\) term and yields only
  ordinary affine-path curvature.  The precise correction is recorded in
  `scalar_shape_acceleration_lemma_retraction.md`; equations (1.2)--(1.5)
  of the original note must not be used.
- **Surviving gates:** the general integrated residual estimate (WK) and a
  mechanism applying it beyond the near-linear spectral branch remain
  load-bearing.  The smooth Gaussian-fiber subclass does not prove general
  WFI or KLS.

## Checkpoint 16 update: posterior reverse smoothing and the nonlinear output gate

- **Forward-sign correction:** the valid unequal-law convolution estimate
  remains (b(1-b)\le4(b-a)), but its correct rearrangement is
  (a\le b(3+b)/4).  It is a forward amplification theorem and gives no
  reverse lower bound on the input gap by itself.
- **Posterior reverse theorem:** for
  (S=(X+G)/\sqrt2), (a=\lambda_1(\mu)), and
  (b=\lambda_1(\mathcal L(S))), posterior unit strong log-concavity gives
  
  \[
    b\le\frac{2a}{1-a},\qquad
    a\ge\frac{b}{2+b}.
  \]
  The proof covers form-domain spectral windows, hard convex supports, and
  intrinsic lower-dimensional supports.
- **Equivalent regularized class:** the fixed-Gaussian output is positive
  analytic, isotropic, and has (0\preceq D^2V\preceq2I).  A universal gap
  on this class is quantitatively equivalent to KLS by the posterior reverse
  theorem.  Score rigidity and the audited weighted-Fisher/MMSE estimate both
  close near-affine bottom modes, but neither controls a genuinely nonlinear
  low mode.
- **Affine-interaction gate:** for
  (R=f((X+G)/\sqrt2)-E[f\mid X]-E[f\mid G]), the candidate estimate
  
  \[
    \operatorname{dist}_{L^2}(f,\mathrm{Aff})^2\le C_A ER^2
  \]
  would combine with reverse smoothing and score rigidity to close KLS.  Its
  kernel is exactly affine in full dimension and the Gaussian sharp constant
  is (2), but no universal proof is known; its quadratic sector already
  contains a dimension-free generalized quadratic-variance theorem.
- **One-form circularity audit:** for gradient fields, the best constant in
  
  \[
    E|w|^2\le C\{E\|Dw\|_{HS}^2+E[w^TD^2Vw]\}
  \]
  is exactly (C_P(\mu)) by Bochner and spectral calculus.  Likewise the
  centered-gradient Korn constant lies in
  ([C_P(\mu)-1,C_P(\mu)]).  Thus generic one-form coercivity is not an
  additional lemma; the surviving target is coherence of an almost-centered,
  slowly varying field which adaptively saturates the posterior covariance.
