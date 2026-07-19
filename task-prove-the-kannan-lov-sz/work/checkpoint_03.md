# Checkpoint 03: finite competitors, transport synchronization, and endpoint incidence

## 1. Exact results added since Checkpoint 02

### 1.1 Interior cross-ray rigidity

If two calibrated endpoint pieces of a one-Lipschitz potential retain an
outward calibrated continuation of length `h`, then

\[
 |x-y|-(f(x)-f(y))
 \ge {h(f(x)-f(y))\over 8(f(x)-f(y)+h)}|N_x-N_y|^2.
\]

On constant-mass bands of balanced scale-`s` rays this gives the uniform
cross slack

\[
 |b-c|-(f(b)-f(c))\ge c_\delta s|N_b-N_c|^2.
\]

Thus all alternative exact calibration between separated directions is
confined to medial or focal endpoints.

### 1.2 Null-invariant multi-Brenier identity

For packet conditionals `mu_i=mu(.|A_i)`, a common absolutely continuous
source `gamma`, Brenier maps `T_i`, weights `lambda_i`, and
`S=sum lambda_i T_i`, one has in the smooth nondegenerate case

\[
 D(S_\#\gamma\Vert\mu)+E\Delta_V+E\log Q
 =\sum_i\lambda_i\log {1\over\mu(A_i)},
\]

where

\[
 Q={\det(\sum_i\lambda_iDT_i)\over
          \prod_i\det(DT_i)^{\lambda_i}}\ge1.
\]

For synchronized rank-one stretches of factor `s` in a tight frame whose
direction covariance is `O(s^{-2})`, `log Q>=s-O(log s)`.  This is the
correct linear scale.

The covariance-to-differential implication is false.  Smooth Brenier maps

\[
 T_i(x)=x+L\theta_\varepsilon(\langle x,u_i\rangle)u_i
\]

have balanced scale-`L` target variance in orthogonal directions while
`E log Q -> 0`: the stretches occur on disjoint rare slabs.  Hence any use
of the determinant identity must derive synchronization from the single
global extremal ray congruence, not from packet covariance or balance.

### 1.3 Exact finite-amplitude ray-height inequality

For a true T3 extremizer and any one-Lipschitz competitor written as

\[
 g=f-h(Q)+r,
\]

extremality gives

\[
 \int B_y(h-\eta h)d\eta\le2\int|r|d\mu,
\]

where

\[
 B_y(c)=2\int_{0<t<c}(c-t)d\nu_y(t)
\]

for `c>=0`, with the reflected formula for `c<0`.  On balanced scale-`s`
rays, `B_y(c)>=c_0 min(c^2/s,|c|)`.

For a finite min envelope `F=min_i ell_i`, the switching defect is exactly

\[
 F_h=F-h_{I(x)}-\max_j
      \{h_j-h_{I(x)}-(\ell_j-F)(x)\}_+.
\]

There is no packet-entropy amplification: the total ideal gain of every
mean-zero multilevel height is at most one half of its range.  Complete or
strong vertex-expander envelopes can absorb every amplitude.  A Hall-type
rematching lemma does give explicit improvement under global same-sign
packet separation, but that separation has not been obtained for the true
bands.

### 1.4 Covariance-saturation inverse theorem

Let

\[
 A_\Omega=E[1_\Omega\sigma_Y^2U_YU_Y^T]\preceq I.
\]

If `A_Omega >= (1-epsilon)P` for an orthogonal projection `P`, then

\[
 E|PM_Y|^2\le\varepsilon\,\operatorname{rank}P,
\]

so the projected ray lines are nearly concurrent.  Distinct saturated
blocks are nearly orthogonal.  Exact smooth saturation gives orthogonal
radial cylinders.

This theorem is sharp but does not apply merely from long reach and
direction dispersion.  Regular simplices and isotropic rounded cubes give
dimension-sharp countermodels to a John/inradius route before balance and
positive quotient mass are used.

### 1.5 Convolution rigidity and its limitation

For normalized self-convolution `T mu=Law((X+Y)/sqrt(2))`, the Hoeffding
residual of `F(x,y)=f((x+y)/sqrt(2))` gives

\[
 \operatorname{Var}_{T\mu}f+{1\over4}E(\Delta_4f)^2
 \le C_P(\mu)\int|\nabla f|^2d(T\mu).
\]

Thus `C_P(T mu)<=C_P(mu)`, and attained equality forces an affine witness
and a Gaussian law.  The first-moment constant is not monotone: normalized
self-convolution increases it for variance-one Laplace and decreases it for
the uniform interval.  A forward strict contraction, even if proved, would
not by itself bound the initial constant without a reverse renormalization
estimate.

## 2. Newly blocked families

1. **Covariance-only multi-Brenier synchronization.**  The rare-slab maps
   give a smooth balanced counterexample to the proposed derivative lemma.
2. **Pairwise determinant accumulation.**  Log-determinant Jensen gaps
   telescope; an expander path does not add them without a new rank-growth
   invariant.
3. **Finite-offset entropy amplification.**  Quotient normalization bounds
   the ideal gain by the height range, independently of the number of
   packets.
4. **Inradius/John/cone geometry before balance.**  The rounded-cube and
   simplex core models have long reach and tight-frame normals but a
   dimension-dependent concurrence error.
5. **Monotone first-moment convolution.**  It fails already in dimension
   one in both directions.

These families may be reopened only with an input tying them to the same
global extremal congruence or with a genuinely new Lyapunov functional.

## 3. Active mechanisms

### 3.1 Endpoint heat-bath rigidity

Every completed ray is an edge between its positive and negative focal
endpoint fibers.  The singular charge is therefore not an arbitrary graph
form: it is a weighted two-block conditional-variance form on this incidence
space, plus the smooth turning form.  In the ideal constant-length model,
maximal gap would force both conditional expectations of
`N=(P_+-P_-)/L` to be constant.  The two reverse martingale identities then
force `N` itself to be constant.  Near-maximal gap should yield a trichotomy:

* aligned normals;
* approximately independent endpoint variables, hence the
  orthogonal-spherical cross-distance model;
* a strict finite competitor.

The load-bearing audit is whether the actual medial coarea weights give the
required conditional-variance normalization and whether the quantitative
error avoids an ambient trace loss.

### 3.2 Saturation decomposition

The exact covariance inverse theorem handles every spectral block on which
the ray covariance nearly exhausts ambient covariance.  The remaining task
is to prove that a long balanced extremizer must contain such a saturated
block, or else that the nonsaturated remainder yields a strict endpoint
competitor.  This is stronger and more precise than a generic inradius
claim.

### 3.3 Soft mass localization

A full-rank variant

\[
 C_t^2=P_{e_t^\perp}+\alpha_t^2e_te_t^T
\]

trades binary-entropy loss for curvature in the protected direction.  The
current audit asks whether an adaptive `alpha_t` can accumulate universal
minimum curvature before the balanced set mass exits a fixed interval.  A
successful estimate would bypass endpoint gluing; a product-exponential
winner-selection model is the required stress test.

## 4. Extremal constraint set

A hypothetical sequence with first-moment scale `s -> infinity` now must
satisfy all of the following simultaneously:

1. a fixed quotient mass of exactly balanced one-dimensional log-concave
   rays has scale comparable to `s`;
2. every coherent unit direction packet has mass at most `C exp(-cs)`;
3. ray covariance is bounded by the identity, while every nearly saturated
   spectral block is an orthogonal radial cylinder and therefore harmless;
4. separated-direction cross calibration is impossible throughout the
   interior core;
5. all finite ray-height perturbations are absorbed at focal/medial
   endpoints at the exact Bregman scale;
6. the endpoint incidence form is near its maximal two-block spectral
   efficiency;
7. separate packet covariance cannot be used as a proxy for synchronized
   transport derivatives.

The unresolved object is consequently a globally log-concave,
nonsaturated, endpoint-expanding incidence geometry with exponentially
diffuse directions and all-amplitude extremal stability.  The next round is
focused on proving that its weighted endpoint heat-bath form cannot exist,
or converting it into the saturated orthogonal-radial blocks above.
