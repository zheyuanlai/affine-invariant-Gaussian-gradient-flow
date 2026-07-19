# Checkpoint 10: tensor contact and longitudinal turning

Date: 2026-07-16

## Status

No complete proof has passed audit.  The exact polyhedral branch remains
closed.  In the smooth branch, a new tensor identity retains the dimension
that scalar curvature and covariance estimates lose.  Its unresolved term is
now an explicit support-contact tensor.  Pointwise support curvature does not
control this tensor; a global longitudinal or multi-packet theorem is needed.

## 1. Tensor Minkowski identity

Let \(\Sigma\) be a smooth free-boundary constant-mean-curvature interface in
a convex body \(K\), with relative perimeter \(P\), normal \(N\), contact
stratum \(\Gamma\), and support conormal \(n_K\).  For any center \(c\), set

\[
 Q=\frac1P\int_\Sigma N\otimes N,\qquad
 X=\frac1P\int_\Sigma(x-c)\otimes N,\qquad
 B=\frac1P\int_\Gamma(x-c)\otimes n_K.
\]

Surface integration by parts gives the exact identity

\[
                         \boxed{I-Q=HX+B}.                    \tag{1.1}
\]

Since \(Q\succeq0\) and \(\operatorname{tr}Q=1\),

\[
                  \|I-Q\|_F^2=d-2+\operatorname{tr}(Q^2).    \tag{1.2}
\]

Thus the left side is of order \(\sqrt d\), even for the fully diffuse law
\(Q=I/d\).  If \(B=0\), the surface-position second moment is at most \(C d\),
and \(|H|\le2P\), then \(P\ge c\).  Therefore a small-perimeter smooth escape
must either have an uncontrolled boundary-position tail or carry a contact
tensor of order \(\sqrt d\).

For a flat complete slice, \(H=0\) and
\(B=I-u\otimes u\).  The contact branch correctly retains cube and simplex
slices instead of declaring them defective.

## 2. Stability and the contact-curvature consequence

For a stable free-boundary CMC interface, test the index form with
\(u_i=N_i-m_i\), where \(m=P^{-1}\int_\Sigma N\).  If
\(\delta^2\le\|Q\|_{\rm op}\), the audited estimate is

\[
 \int_\Gamma \mathrm{II}_{\partial K}(N,N)
 \le \frac{2\delta}{(1-\delta)^2}\int_\Sigma|A|^2.            \tag{2.1}
\]

In the retained effective-rank branch, the coefficient is below \(0.150\).
Combining (2.1), a global two-sided tube-defect estimate at length \(T\),
and the contact tensor gives

\[
 \int_\Gamma
   \frac{|x-c|^2}{\mathrm{II}_{\partial K}(N,N)}
 \ge
 \frac{P d(1-\delta)^2T^2}{64\delta\varepsilon}.             \tag{2.2}
\]

At \(T=\gamma/\psi=\gamma/(2P)\), the right side has order \(d/P\).
This is the strongest closed tensor-plus-stability consequence.

## 3. The quartic support audit

The inference

\[
 \mathrm{II}_{\partial K}(N,N)=0
 \quad\Longrightarrow\quad
 \partial K\hbox{ is ruled in direction }N
\]

is false, even for \(C^\infty\) convex supports.  The graph

\[
 x_d=x_1^4+\sum_{j=2}^{d-1}x_j^2
\]

has zero second fundamental form in direction \(e_1\) along \(x_1=0\), but
every nonzero tangent segment leaves the support.  Stability sees only the
2-jet at contact.  A valid theorem must control the longitudinal quantity

\[
 g(R)=\int_0^R(R-s)g''(s)\,ds
\]

in each convex normal section, or use a finite chord competitor.  This
quartic example invalidates any pointwise-curvature-to-cell step.

## 4. Weighted random fibers

For a balanced set \(E\), direction \(\theta\), unnormalised conditional
fiber masses \(m_y^\pm\), and normalized conditional deviation
\(\sigma_y\), define

\[
 B_\theta(E)=\int_{\theta^\perp}
       \frac{\min(m_y^+,m_y^-)}{\sigma_y}\,dy.
\]

The exact Crofton and one-dimensional estimates are

\[
 J_\theta(E)\ge\frac1{2\sqrt3}B_\theta(E),\qquad
 \mathbb E_\theta J_\theta(E)=c_dP(E),\qquad
 \frac1{\sqrt{2d}}\le c_d\le\frac1{\sqrt d}.                 \tag{4.1}
\]

Hence

\[
                \mathbb E_\theta B_\theta(E)\ge c/\sqrt d   \tag{4.2}
\]

for every exact balanced minimizer would close KLS.  The scale is correct on
Gaussian halfspaces, exact cube halfspaces, simplex caps, radial exponential
spheres, and product-exponential max sets.

Exact fiberwise quantile replacement lowers only \(J_\theta\); its transverse
graph perimeter is uncontrolled.  Majority completion has mass error
\(I_\theta\), and max-flow/min-cut shows that the projected marginal
inequality recovers exactly the original quasiminimality budget.  Thus
one-directional dimension descent is tautological.  Strict saving must use
multi-directional coherence or the support-contact tensor.

## 5. Global transport-ray compatibility

If two rays of a single 1-Lipschitz eikonal potential contain
\([-R,R]\), cross-endpoint Lipschitz inequalities give

\[
                         R|u-v|\le|y-z|.                     \tag{5.1}
\]

This is a genuine global co-Lipschitz constraint on the orientation map.
Using fourth moments of projected conditional scales improves the spectral
multiplicity estimate to

\[
 N_K(\alpha)\le C\alpha^{-17/5},\qquad
 \operatorname{tr}K\le C d^{12/17},\qquad
 \mathbb E(1/\sigma)\ge c d^{-6/17}.                        \tag{5.2}
\]

The power improvement is only a search milestone and is not a returnable
result.  Dyadic summation still diverges.

A cube checkerboard gives one globally log-concave, globally eikonal,
balanced model with diffuse projective orientations and zero Gauss
derivative on every regular orientation patch.  Its transition strata make
\(\mathbb E(1/\sigma)\) large.  Thus a successful theorem must charge
orientation transitions; regular-chart Jacobians alone cannot.

At support contact, a local \(D\)-cluster of \(R\)-long rays is necessarily
one projective cap and gives an approximate common kernel of the contact
normal matrix and contact tensor.  The surviving global escape is many
spatially \(R\)-separated packets with different kernels.

## 6. New live routes

1. **Longitudinal support turning.**  Convert the convex-section identity
   \(g(R)=\int_0^R(R-s)g''(s)ds\) and normal-ray survival into a finite
   chord competitor or a bounded-reuse cell assignment.
2. **Median signed-distance extremizer.**  Replace a median-centered
   1-Lipschitz maximizer by the signed distance to its sign interface,
   compute the shape derivatives of its normal Voronoi cells, and test
   whether global maximality forces cellwise balance and transition cost.
3. **T3 quotient second variation.**  For ray-constant perturbations of a
   first-moment optimizer, construct a globally feasible second-order
   eikonal perturbation and derive the exact quotient stability form.
4. **Multi-packet aggregation.**  Charge spatially separated contact
   packets with different approximate kernels to isotropic volume, without
   falling back to a trace bound.

## 7. Blocked list

- Pointwise support curvature, including exact vanishing of one principal
  curvature at contact.
- One-directional majority-fiber descent or marginal mass repair.
- Regular Gauss-Jacobian control without transition-stratum charge.
- Dyadic spectral multiplicity with any nonsummable power law.
- Local tube, Jacobi, or normal-matrix estimates without longitudinal
  support control.
- Any smooth-to-polyhedral passage whose bevel gain is only of the same
  quadratic order as its mesh error.

## 8. Next audit gate

A closing lemma must survive both extremes:

- flat facets, where all local support curvature vanishes and the contact
  tensor is the entire signal; and
- quartically flat smooth support, where the contact 2-jet vanishes but no
  ruling exists.

It must also aggregate a diffuse family of separated packets without a
dimension-dependent covariance trace estimate.  A statement that handles
only one local packet is insufficient.
