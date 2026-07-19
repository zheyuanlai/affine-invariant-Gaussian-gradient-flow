# Physical splicing by translated patches: exact identities and a flat-phase obstruction

## 0. Conclusion

Let the notation and the fixed numerical output of
`fixed_scale_physical_splicing.md` be in force.  Thus

\[
 \operatorname {tr}M_{\rm phys}>.004p,\qquad
 {\operatorname {tr}M_{\rm phys}\over\|M_{\rm phys}\|_{op}}>17,
 \qquad \mathfrak G_{\rm fin}(F_0)\le 6.02\,10^{-5}p.       \tag{0.1}
\]

Small translations of disjoint boundary patches, followed by a compensating
patch, do **not** prove the fixed physical-splice lemma from (0.1).  The
obstruction is exact and persists for finite displacements and nonsmooth
convex potentials.

In a flat chart in which

\[
 A=\{(z,s):s<0\},\qquad e^{-V(z,s)}=w(z)e^{-\kappa s},       \tag{0.2}
\]

replace the boundary by the compactly supported graph `s=h(z)`.  If
`Delta mu` and `Delta P` denote the resulting changes of weighted volume and
perimeter, then

\[
 \boxed{\quad
 \Delta P+\kappa\Delta\mu
 =\int w(z)e^{-\kappa h(z)}
       \big(\sqrt{1+|\nabla h(z)|^2}-1\big)\,dz\ge0.
 \quad}                                                    \tag{0.3}
\]

Consequently, any finite collection of such patches with the same `kappa`,
including a same-`kappa` volume-correction patch, satisfies

\[
                    \sum_j\Delta\mu_j=0
       \quad\Longrightarrow\quad \sum_j\Delta P_j\ge0.    \tag{0.4}
\]

There is no Taylor remainder in (0.3), and the traces on the surrounding
balls agree exactly.  Thus these are legitimate finite physical splices.
Product-exponential box facets realize (0.2) with `kappa=1` in arbitrarily
high rank.  In fact, even translating their entire facets, within the class
of unequal boxes, increases perimeter at fixed mass.  Their known saving
comes from a diagonal bevel of the ridges, not from translations.

For a general convex potential, the right side of (0.3) becomes

\[
                 \text{cutoff/tilt cost}-\text{normal convexity gain}. \tag{0.5}
\]

The matrix in (0.1) controls neither term.  In the smooth infinitesimal
calculation, averaging localized translation speeds over a complete
orthonormal basis gives exactly

\[
 \sum_{i=1}^n Q(\chi N_i)
 =\int_\Sigma\big(|\nabla_\Sigma\chi|^2
             -\chi^2\nabla^2V(N,N)\big)e^{-V}d\mathcal H^{n-1}. \tag{0.6}
\]

All shape curvature cancels.  Isotropy gives no lower bound for the second
term on the selected physical patches, while (0.1) gives no upper bound for
the cutoff capacity in the first term.  Flat cells have zero normal
convexity and make (0.6) nonnegative.

Therefore translations can only contribute after adding a new, genuinely
spatial hypothesis: a fixed-scale packet/capacity bound or a ridge/focal
conductance bound.  The latter is exactly the missing incidence datum already
identified in `fixed_scale_physical_splicing.md`.  Calling the flat survivor
"product" does not solve the issue: the translation identities contain no
estimate of the Poincare constant of a remaining factor.

## 1. Exact finite graph variations for a convex density

The next lemma fixes all signs and requires no differentiability of `V`.

**Lemma 1.1 (finite flat-chart identity).**  Let `Q` be a bounded open subset
of `R^(d-1)`, let `a>0`, and suppose that in
`C=Q times (-a,a)` the initial finite-perimeter set is `A={s<0}`.  Write

\[
                         g_z(s)=e^{-V(z,s)}.                \tag{1.1}
\]

Assume that `s mapsto V(z,s)` is convex and finite on `(-a,a)` for almost
every `z`.  Suppose that one number `kappa` belongs to
`partial_s V(z,0)` for almost every `z`.  If
`h in W^{1,infinity}_0(Q)` and `|h|<a`, replace `A` in `C` by
`A_h={s<h(z)}`.  Then

\[
 \Delta\mu=\int_Q\int_0^{h(z)}g_z(s)\,ds\,dz,              \tag{1.2}
\]

\[
 \Delta P=\int_Q\left[g_z(h(z))\sqrt{1+|\nabla h(z)|^2}
                         -g_z(0)\right]dz,                 \tag{1.3}
\]

and

\[
 \boxed{\quad
 \Delta P+\kappa\Delta\mu
 =T_\kappa(h)-D_\kappa(h),\quad}                           \tag{1.4}
\]

where

\[
 T_\kappa(h)=\int_Qg_z(h(z))
               \big(\sqrt{1+|\nabla h(z)|^2}-1\big)dz\ge0 \tag{1.5}
\]

and

\[
 D_\kappa(h)=\int_Q\left[g_z(0)-g_z(h(z))
              -\kappa\int_0^{h(z)}g_z(s)ds\right]dz\ge0.  \tag{1.6}
\]

The modification has exactly the old trace on `partial C` if `h` vanishes
in a neighborhood of `partial Q`.

**Proof.**  Equations (1.2) and (1.3) are Fubini's theorem and the area
formula for a graph.  For fixed `z`, put

\[
                         W(s)=V(z,s)-\kappa s.              \tag{1.7}
\]

The subgradient hypothesis says that the convex function `W` has a minimum
at zero.  In the sense of one-dimensional distributional derivatives,

\[
 {d\over ds}g_z(s)+\kappa g_z(s)=-W'(s)g_z(s).             \tag{1.8}
\]

For `s>0` the right side is nonpositive, and for `s<0` it is nonnegative.
Integrating from zero to `h`, with the orientation of the integral retained
when `h<0`, gives

\[
                   g_z(h)-g_z(0)+\kappa\int_0^h g_z\le0.  \tag{1.9}
\]

This proves (1.6).  Adding and subtracting `g_z(h)` in (1.3) proves (1.4).
Approximation of a one-dimensional convex function by convolution proves
the same calculation at kinks and at points where only one-sided derivatives
exist.  The trace assertion follows from `h=0` near `partial Q`.  QED.

The two terms in (1.4) have an unambiguous meaning.  `T_kappa` is the area
paid for tapering a translated patch back to the unchanged trace.
`D_kappa` is the gain supplied by strict convexity of the potential in the
normal direction.  If `V(z,s)=V(z,0)+kappa s`, then `D_kappa=0`, and (1.4)
is (0.3).

A rigid ambient translation has no additional favorable first-order term.
Its tangential component only reparametrizes an infinite planar boundary;
its shape speed is the normal component.  On a bounded patch, moving the
tangential footprint creates seams at the cutoff.  Those seams are another
form of `T_kappa` (or a singular ridge cost), and the physical normal matrix
does not control them.

For disjoint charts with possibly different slopes `kappa_j`, exact mass
balance gives

\[
 \Delta P= -\sum_j\kappa_j\Delta\mu_j
             +\sum_jT_{\kappa_j}(h_j)
             -\sum_jD_{\kappa_j}(h_j),
 \qquad \sum_j\Delta\mu_j=0.                              \tag{1.10}
\]

Thus a first-order volume transfer can detect a difference between two
weighted mean curvatures.  If all slopes agree, the whole first-order term
cancels exactly.  A fixed numerical saving then requires a quantitative
lower bound for `sum D_j-sum T_j`.  Neither weighted facet area nor normal
rank supplies such a bound.

### 1.1 Exact volume correction

Let `h_1,...,h_k` be prescribed in disjoint charts and let a further chart
have a nonnegative, nonzero cutoff `phi`.  The map

\[
 t\longmapsto \int_Q\int_0^{t\phi(z)}g_z(s)dsdz             \tag{1.11}
\]

is continuous and strictly increasing while its graph remains in the chart.
Hence every sufficiently small volume error has a unique exact correction.
No implicit smoothness assumption is needed.  If the correction chart is
log-affine with the same `kappa` as the original charts, (0.4) shows that
the corrected finite competitor still cannot save perimeter.  If its slope,
area, or available height is uncontrolled, its cost is likewise uncontrolled.

### 1.2 The second-order limit

When `V` is `C^2` in the normal variable, `h=t phi`, and
`partial_sV(z,0)=kappa`, (1.4) gives

\[
 T_\kappa(t\phi)={t^2\over2}\int_Qg_z(0)|\nabla\phi|^2dz+o(t^2),
 \quad
 D_\kappa(t\phi)={t^2\over2}\int_Qg_z(0)V_{ss}(z,0)\phi^2dz+o(t^2). \tag{1.12}
\]

The little-oh is only local.  Without a uniform reach and a uniform modulus
for the density it cannot be turned into the fixed saving `10^-4 p`.
Formula (1.4), rather than (1.12), is the safe statement for a nonsmooth
log-concave measure.

## 2. Smooth shape variation and orthogonal averaging

This section is diagnostic only; Lemma 1.1 is the nonsmooth finite statement.
Let `Sigma` be a `C^3` boundary in a region where `V` is `C^2`, let `N` be
the outward normal, use `II(tau)=D_tau N`, and put

\[
 H_\mu=\operatorname {tr}\mathrm{II}-\langle\nabla V,N\rangle,
 \qquad q=|\mathrm{II}|^2+\nabla^2V(N,N).                  \tag{2.1}
\]

For normal speed `u`,

\[
 \delta\mu(A)=\int_\Sigma u\,d\sigma_\mu,
 \qquad \delta P=\int_\Sigma H_\mu u\,d\sigma_\mu.      \tag{2.2}
\]

At a weighted-CMC patch, `H_mu=lambda`, the second variation of
`P-lambda mu` is

\[
 Q(u)=\int_\Sigma\big(|\nabla_\Sigma u|^2-q u^2\big)d\sigma_\mu. \tag{2.3}
\]

Take a compactly supported cutoff `chi` and the localized translation-normal
speeds

\[
                              u_i=\chi N_i.                \tag{2.4}
\]

For every tangent vector `tau`,

\[
 \sum_i|\partial_\tau(\chi N_i)|^2
 =|\partial_\tau\chi|^2+\chi^2|D_\tau N|^2,               \tag{2.5}
\]

because the cross term is `chi partial_tau chi partial_tau |N|^2=0`.
Summing (2.3) over a full ambient orthonormal basis proves (0.6):

\[
 \boxed{\quad
 \sum_{i=1}^nQ(\chi N_i)
 =\int_\Sigma\left(|\nabla_\Sigma\chi|^2
       -\chi^2\nabla^2V(N,N)\right)d\sigma_\mu.
 \quad}                                                    \tag{2.6}
\]

The `|II|^2` terms cancel exactly.  For a flat patch, (2.6) is the
quadratic limit of Lemma 1.1.

There are three losses before (2.6) can produce a physical splice.

1. Each speed must satisfy an exact volume constraint.  Subtracting a mean
   or using a compensating patch adds a second Dirichlet/Jacobi cost which
   is not present in `M_phys`.
2. The negative term wins only if

   \[
       \int\chi^2\nabla^2V(N,N)d\sigma_\mu
             >\int|\nabla\chi|^2d\sigma_\mu.               \tag{2.7}
   \]

   Boundary area and effective normal rank imply neither side of (2.7).
3. If the `n` translation directions are treated as `n` alternative
   competitors, the splice functional can retain only the best one.  From a
   trace estimate one obtains at best

   \[
      \max_i[-Q(u_i)]_+\ge {1\over n}
                   \left[-\sum_iQ(u_i)\right]_+,           \tag{2.8}
   \]

   which loses the dimension.  Combining unit directions into one ambient
   translation requires normalization of the coefficient vector and causes
   the same loss.  Summing without loss requires disjoint spatial packets
   that can be moved simultaneously; (0.1) contains no such packetization.

Using only eighteen directions does not fix the last point.  Effective rank
larger than seventeen means only
`||M_phys||_op<tr(M_phys)/17`.  For
`M_phys=(tr M_phys/n)I`, every fixed eighteen-dimensional subspace captures
only `18/n` of the trace.

There is also a prior variational gap.  A regular heat level `A_r` is a
smooth set, but it is not known to be weighted-CMC or locally minimizing.
The bound on `G_fin` is an integrated finite deficit, not an Euler--Lagrange
identity for each level.  Obtaining (2.3) for the original levels would
require a quantitative variational-selection theorem which preserves the
weighted normal submeasure.  No such theorem follows from (0.1).  Moreover,
`M_phys` sees only the selected weight `omega`, whereas every cutoff and
volume corrector pays against the full boundary measure.

## 3. Why convexity, isotropy, and the coarea rank do not supply the missing term

### 3.1 Isotropy is a bulk condition

For a smooth full-support density, integration by parts may give bulk
identities such as

\[
 \int x_i\partial_jV\,d\mu=\delta_{ij},\qquad
 \int\partial_{ij}V\,d\mu=\int\partial_iV\partial_jV\,d\mu, \tag{3.1}
\]

under the required integrability.  These identities do not lower-bound
`nabla^2V(N,N)` on a selected boundary submeasure.  For polyhedral
potentials, the Hessian is a positive matrix-valued measure supported on the
walls between affine cells.  Boundary patches can lie wholly inside those
cells, where their normal convexity defect (1.6) is exactly zero.  Covariance
normalization does not move that curvature onto the selected patches.

An explicit full-support isotropic-up-to-scale example is

\[
 d\mu_\varepsilon(x)=Z^{-1}
       \exp\{-\|x\|_\infty-\varepsilon\|x\|_1\}\,dx.      \tag{3.2}
\]

Signed-permutation symmetry makes its covariance a scalar matrix, so a
scalar dilation makes it isotropic.  It is not a product density.  On the
core of the facet `x_i=R` of the cube `A_R={||x||_infty<R}`, with
`|x_j|<R-delta` for `j!=i`,

\[
 V_\varepsilon(z,R+s)=V_\varepsilon(z,R)+(1+\varepsilon)s. \tag{3.3}
\]

All `2n` facet cores therefore have the common slope `1+epsilon`, and
Lemma 1.1 forbids any saving by mass-balanced localized translations of
those cores.  Their normal projector is a scalar multiple of `I_n`.
The normal line through a generic point of a cube facet does not pass through
one common center, so the family is not Euclidean-concurrent.  Whatever
global competitor excludes this set must use its ridges or other global
geometry; it is not encoded by translation stability.

At `epsilon=0`, `R=||X||_infty` has the `Gamma(n,1)` law and

\[
 \operatorname {Cov}(X)={(n+1)(n+2)\over3}I_n.             \tag{3.4}
\]

On a facet at `R=n`, the core at distance one from all ridges retains the
fraction

\[
                         (1-1/n)^{n-1}\longrightarrow e^{-1} \tag{3.5}
\]

of its weighted area.  Thus even a fixed fraction of high-rank flux can be
carried by mutually separated flat cores with zero normal convexity.

### 3.2 The physical rank is aggregated over levels

The matrix in (0.1) is

\[
 M_{\rm phys}=\int_0^1M_rdr,
 \quad M_r=\int_{\partial^*A_r}\omega nn^T e^{-V}d\mathcal H^{n-1}. \tag{3.6}
\]

High rank of the integral does not imply high rank of `M_r` on a positive
measure set of levels.  Rank-one matrices with rotating directions can have
an isotropic integral.  Nesting of the sets may constrain a literal family
of global halfspaces, but it does not prevent the selected weight `omega`
from following one rotating local patch on each level.  Since mass correction
and splicing occur separately for every `r`, translations cannot exchange
volume or perimeter between different levels.  A valid proof would need a
new measurable within-level packet theorem.

Rectifiability does not provide the required scale.  It supplies tangent
planes at almost every boundary point after passing to arbitrarily small
balls.  In such a ball a cutoff has energy of inverse-square scale, whereas
the normal-convexity gain is quadratic in the displacement.  No fixed
`10^-4 p` saving follows without a lower bound on patch reach or on a
codimension-two junction capacity.

### 3.3 A transverse-factor obstruction to classification from local shape

Let `eta` be any isotropic log-concave probability and let `xi` be an
isotropic flat-phase model such as a shifted one-sided-exponential product.
Then `eta tensor xi` is isotropic, and every calculation in Lemma 1.1 on a
set and patches depending only on the `xi` coordinates is unchanged.  On the
other hand, Poincare tensorization gives

\[
 C_P(\eta\otimes\xi)=\max\{C_P(\eta),C_P(\xi)\}.           \tag{3.7}
\]

Thus no classification based only on the local translated-patch identities
can conclude a bounded Poincare constant: it would have to bound the
arbitrary transverse factor.  This is not asserted to be a counterexample to
the full heat-generated fixed-splice lemma--near-Cheeger extremality could
force the witness to live in the bad factor.  Proving that force, however,
is precisely additional global information and is absent from translation
stability.

## 4. Model tests

### 4.1 Gaussian ball

For standard Gaussian measure and the ball `B_R`,

\[
 N=x/R,\quad \mathrm{II}=R^{-1}I_T,\quad
 H_\gamma=(n-1)/R-R,\quad q=(n-1)/R^2+1.                 \tag{4.1}
\]

The translation-normal speed `u_a=a dot N` has zero surface mean.  With
surface measure normalized to a probability,

\[
 E u_a^2={|a|^2\over n},\qquad
 E|\nabla_\Sigma u_a|^2={n-1\over nR^2}|a|^2,             \tag{4.2}
\]

and hence

\[
                         Q(u_a)=-{P_\gamma(B_R)\over n}|a|^2. \tag{4.3}
\]

An exact-volume translation, with a second-order radial correction, therefore
saves

\[
                   {P_\gamma(B_R)\over2n}|a|^2t^2+o(t^2). \tag{4.4}
\]

The saving of any unit translation is `1/n` of the perimeter.  Summing the
`n` values gives a fixed trace, but they are alternative competitors, not
additive savings.  This is the dimension loss in (2.8) and the precise reason
the concurrent radial branch is necessary.

### 4.2 Uniform cube

Inside a uniform convex body, a boundary patch away from the support boundary
has `kappa=0` and `D_0=0`.  Therefore

\[
 \Delta P=\int_Q\rho
           \big(\sqrt{1+|\nabla h|^2}-1\big)dz\ge0        \tag{4.5}
\]

for every localized translation of an inner-box facet.  The equal facets of
an inner cube have normal matrix `P I_n/n`.  Moving and compensating separated
facet cores does not save; beveling their ridges does.  The actual coordinate
half-cube is the rank-one affine branch.

### 4.3 Uniform regular simplex

The same identity (4.5) holds on the interior cores of a homothetic simplex.
For a regular simplex the `n+1` equal facet normals form a tight frame, so the
normal matrix of the equal cores is again proportional to `I_n`.  Local
translations plus a flat compensator cannot save.  The gain comes from
dihedral ridge bevels or from comparison with a barycentric halfspace.  The
latter is the affine branch.  Contact with the support boundary cannot be
discarded in a smooth second-variation argument.

### 4.4 One-sided exponential products

Let

\[
 d\mu_m(x)=e^{-\sum_{i=1}^m x_i}1_{\{x_i\ge0\}}dx,
 \qquad B(q_1,\ldots,q_m)=\prod_i[0,q_i].                  \tag{4.6}
\]

Put `a_i=1-e^{-q_i}`.  Relative to the support,

\[
 \mu_m(B)=\prod_i a_i=:M,\qquad
 P_{\mu_m}(B)=M\sum_i{1-a_i\over a_i}.                    \tag{4.7}
\]

With `s_i=log a_i`, fixed mass is the constraint `sum s_i=log M`, and

\[
 {P_{\mu_m}(B)\over M}=\sum_i(e^{-s_i}-1).                 \tag{4.8}
\]

Strict convexity and Jensen show that the equal box uniquely minimizes (4.8)
among all boxes of the same mass.  Thus even whole-facet threshold
translations cannot save perimeter.

For the half-mass equal box, let

\[
 (1-e^{-q_m})^m={1\over2},\qquad u_m=e^{-q_m}=1-2^{-1/m}.  \tag{4.9}
\]

Every upper face has slope `kappa=1`, and the normal matrix is
`P I_m/m`.  The part of a face satisfying `x_j<=q_m-delta` for all `j!=i`
has relative face area

\[
 R_{m,\delta}=\left({1-e^\delta u_m\over1-u_m}\right)^{m-1}
 \longrightarrow \exp\{-(e^\delta-1)\log2\}>0.           \tag{4.10}
\]

These cores stay a fixed distance `delta` from every upper ridge and still
carry a fixed fraction of the full high-rank flux.  Lemma 1.1 says that any
finite mass-balanced collection of translated cores, with a translated core
as compensator, has nonnegative perimeter change at every admissible
amplitude.  The simultaneous diagonal bevel in
`fixed_scale_physical_splicing.md` does save a fixed amount.  Hence the model
separates the two mechanisms cleanly: translation stability sees nothing,
whereas ridge capacity closes the example.

For literal compact charts inside the open orthant, additionally impose
`x_j>=m^{-2}` on the cores and taper before that lower boundary.  For each
fixed `m` this puts the whole modification in the interior of the support,
and the extra retained-area factor tends to one as `m` tends to infinity.
Thus the hard support creates no hidden error in the saturation example.

Translation by the barycenter makes (4.6) isotropic without changing any
perimeter statement.  The density is strongly non-symmetric, so the
obstruction is not an artifact of central symmetry.

## 5. What a viable translation lemma would have to assume

The strongest statement justified by the calculations is conditional.  On
flat convex charts with a common weighted curvature `-kappa`, an exactly
balanced graph splice saves precisely when

\[
                \sum_jD_\kappa(h_j)-\sum_jT_\kappa(h_j)>0. \tag{5.1}
\]

In the smooth limit this asks for a negative-capacity inequality

\[
 \sum_j\int\chi_j^2\nabla^2V(N_j,N_j)d\sigma_\mu
 >\sum_j\int|\nabla\chi_j|^2d\sigma_\mu
       +\text{volume-correction cost}.                    \tag{5.2}
\]

To reach the numerical target in (0.1), the surplus in a finite version of
(5.2) must be at least `10^-4 p`.  The following data would suffice in
principle, but none is a consequence of the current physical trace bound:

* a measurable, within-level decomposition into disjoint patches carrying a
  fixed fraction of `M_phys`;
* a fixed lower bound on their reach and transferable volume;
* a fixed upper bound on the cutoff and compensator capacities; and
* either fixed normal convexity on those patches or a ridge/focal graph with
  a fixed conductance spectral gap.

The flat log-affine identity proves that the normal-convexity alternative can
vanish identically.  The product-exponential calculation proves that global
facet translations can remain stable even then.  What remains is therefore
the ridge/focal alternative or a genuinely global classification theorem.
Deriving a dimension-free Poincare bound merely from the words "flat",
"product-like", or "radial" would be circular: the displayed variation
identities are homogeneous in boundary measure and contain no Cheeger height
or spectral gap of a transverse factor.

Accordingly, the proposed translation/shape-stability mechanism does not
establish either a fixed saving `10^-4 p` or a bounded-`K` classification.
Its exact useful output is (1.4), which identifies the missing quantity.  A
proof of the fixed physical-splice lemma must import new cross-patch spatial
information, most naturally the admissible ridge conductance already isolated
by the finite bevel construction.
