# Large interior filling: max-flow certificates, calibrated rigidity, and two no-go theorems

## 0. Outcome

This note audits the proposed large-interior-filling branch.  It proves the
sharp statements that follow from `BV` filling, max-flow/min-cut duality,
and calibration, with all factors of two fixed.

Let `A` be a finite-perimeter set and let its exposed boundary be partitioned
into phase traces.  If `Fill` denotes the weighted area of the interior
interfaces, counted once, then the cut error in the component inverse is

\[
                         R_{\rm cut}=2\,\operatorname {Fill}.       \tag{0.1}
\]

Consequently the already proved softmax inverse applies only when

\[
 2\operatorname {Fill}+E_{\rm tail}\le1.4\,10^{-4}p.              \tag{0.2}
\]

In particular, with no tail error, the numerical cheap-filling threshold is
`7*10^-5 p`, not `1.4*10^-4 p`.

For large filling there is an exact and useful refinement.  Each cell has a
corrected exposed flux `F_i`.  If `r_i` is the area of its interior boundary,
counted on the cell side, then

\[
 r_i-|F_i|={1\over2}\int_{\Lambda_i}|N_i-u_i|^2\,d\sigma_\mu,
 \qquad u_i=-{F_i\over|F_i|}.                         \tag{0.3}
\]

Thus a large filling splits into a **calibrated part** `sum |F_i|` and a
nonnegative **Plateau excess** `sum(r_i-|F_i|)`.  In the calibrated branch,
an expanding filling-adjacency graph forces all filling normals onto one
projective line, quantitatively; a graph bottleneck is an actual aggregate
component cut.  These are Theorems 3.1 and 4.1 below.

The hoped-for unconditional last step is false at the level of the data
currently available.  Two decisive obstructions are proved.

1. Large `BV` filling gives a *lower* bound on quadratic phase capacity,
   not an upper bound from which a cheap phase interpolation can be built.
   Rectangles attain the lower bound exactly, and scaling their height makes
   the ratio between quadratic capacity and filling arbitrary.
2. A max-flow certificate need not be parallel, concurrent, or
   orthogonal-radial.  Thin tubes around arbitrary embedded curves have
   nearly maximal flows whose direction follows the curve and can undergo
   any prescribed sequence of rotations.  This remains true in arbitrarily
   high ambient dimension.  Such tubes are not near-Cheeger sets; therefore
   the example identifies, rather than evades, the missing hypothesis: one
   must use near-minimality to control the Plateau excess or the turning of
   the max flow.

Accordingly no perimeter saving of `1.4*10^-4 p` follows from “large
filling” alone.  The exact additional lemma needed for closure is stated in
Section 8.  Assuming it, even implicitly, would be a new geometric inverse
of conjecture strength.

All results are stated first for a bounded Lipschitz set and a positive
`C^1` density.  The `BV` identities extend to a convex lower-semicontinuous
potential and convex support by monotone truncation and the normal-trace
Gauss--Green theorem; no regularity limit changes the constants.

## 1. Filling and the exact factor of two

Let `Omega subset R^n` be bounded and Lipschitz and let

\[
 d\mu=\rho(x)\,dx,\qquad \rho\in C^1(\overline\Omega),
 \qquad 0<\inf\rho\le\sup\rho<\infty.                \tag{1.1}
\]

Let `A subset Omega` have finite perimeter.  A phase partition of `A` is a
Caccioppoli partition `A_1,...,A_m`: the cells are disjoint, their union is
`A` modulo null sets, and `sum_i P_mu(A_i)<infty`.  Write

\[
 \Gamma_i=\partial^*A_i\cap\partial^*A\cap\Omega,
 \qquad
 \Lambda_{ij}=A^{(1)}\cap\partial^*A_i\cap\partial^*A_j.          \tag{1.2}
\]

Put

\[
 a_i=\sigma_\mu(\Gamma_i),\quad
 \ell_{ij}=\sigma_\mu(\Lambda_{ij}),\quad
 r_i=\sum_{j\ne i}\ell_{ij},\quad
 L=\sum_{i<j}\ell_{ij}.                                \tag{1.3}
\]

Here `L` is the interior filling counted once.  The reduced-boundary
partition identity gives

\[
 \boxed{\quad
 \sum_iP_\mu(A_i)=P_\mu(A)+2L,
 \qquad \sum_i r_i=2L.\quad}                         \tag{1.4}
\]

Indeed, an exposed boundary point belongs to exactly one cell boundary,
whereas an interior interface belongs to exactly two, with opposite
measure-theoretic normals.  Multiplication by the common nonnegative density
and monotone approximation prove the weighted statement.

If the prescribed phase traces are given before the cells, define

\[
 \operatorname {Fill}_A(\Gamma_1,\ldots,\Gamma_m)
 =\inf_{\{A_i\}}\sum_{i<j}\sigma_\mu(\Lambda_{ij}),              \tag{1.5}
\]

where the infimum is over Caccioppoli partitions with the prescribed trace
on `partial^* A`.  Formula (1.4) proves (0.1), including when the infimum is
not attained, by relaxation.

**Numerical consequence.**  In Lemma 2.1 and Theorem 4.1 of
`tube_product_inverse.md`, the seam error is `R=sum_i r_i`.  Therefore its
condition `R+E_tail<=1.4*10^-4 p` is exactly (0.2).

## 2. Binary least-gradient max-flow/min-cut

The binary version admits a complete dual description.  Let
`Gamma_+` and `Gamma_-` be disjoint relatively open pieces of
`partial A cap Omega`, up to boundary-null sets.  Boundary left over by the
selected packet is free.  Define

\[
 \Phi_1(\Gamma_+,\Gamma_-)
 =\inf\left\{\int_A\rho\,|Du|:
 u\in BV(A),\ 0\le u\le1,
 \operatorname {Tr}u=1\text{ on }\Gamma_+,
 \operatorname {Tr}u=0\text{ on }\Gamma_-\right\}.              \tag{2.1}
\]

Coarea and thresholding show that the same infimum is obtained over
characteristic functions.  Thus `Phi_1` is the least weighted interior
separator for the binary trace.

**Theorem 2.1 (continuous max-flow/min-cut).**  Under (1.1),

\[
 \boxed{
 \Phi_1(\Gamma_+,\Gamma_-)=
 \sup_Z\int_{\Gamma_+}\rho Z\cdot N_A\,d\mathcal H^{n-1},}       \tag{2.2}
\]

where the supremum is over `Z in L^infty(A;R^n)` satisfying

\[
 |Z|\le1\quad\text{a.e.},\qquad
 \operatorname {div}(\rho Z)=0\quad\text{in }\mathcal D'(A),
 \qquad (\rho Z)\cdot N_A=0
 \text{ on the free boundary}.                       \tag{2.3}
\]

The normal traces on `Gamma_+` and `Gamma_-` have opposite total flux.

**Proof.**  If `u` is primal-admissible and `Z` is dual-admissible, the
Anzellotti pairing and Gauss--Green give

\[
 \int_{\Gamma_+}\rho Z\cdot N_A
 =\int_A (\rho Z,Du)
 \le\int_A\rho|Du|.                                  \tag{2.4}
\]

This is weak duality.  For strong duality, first replace the trace constraint
by the standard `BV` boundary penalty and write total variation as the
support function of the closed unit ball in `L^infty`.  Fenchel--Rockafellar
duality applies because a Lipschitz extension of the boundary datum is a
continuity point of the penalty.  The adjoint constraint is exactly
`div(rho Z)=0`, and the free-boundary part of the normal trace must vanish.
Truncation gives `0<=u<=1`; the coarea formula gives a characteristic
minimizer at almost every threshold.  Smooth approximation of `Z` proves
the displayed formulation without changing the supremum.  QED.

Removing the cyclic part of the vector measure `rho Z dx`, Smirnov's
decomposition represents it as a superposition of oriented Lipschitz arcs
from `Gamma_+` to `Gamma_-`.  Its total endpoint mass is the flux in (2.2),
and its total curve length, with multiplicity, is

\[
                         \int_A\rho|Z|\,dx\le\mu(A).              \tag{2.5}
\]

This is a useful flow certificate, but (2.5) controls neither the turning of
the arcs nor one common direction for different arcs.

## 3. Corrected flux and exact rigidity of its saturated part

Let `N_A` be the outer normal of `A`.  Let `N_i` be the outer normal of cell
`A_i` on its interior boundary `Lambda_i=union_{j ne i}Lambda_ij`.  If the
ambient support has boundary, write `C_i=partial^*A_i cap partial Omega` and
let `n_Omega` be its outer normal.  Define

\[
 F_i=\int_{\Gamma_i}N_A\,d\sigma_\mu
     +\int_{A_i}\nabla V\,d\mu
     +\int_{C_i}n_\Omega\,d\sigma_\mu,               \tag{3.1}
\]

where `rho=Z^{-1}e^{-V}`.  The last term is omitted for full support.  The
weighted Gauss--Green identity is

\[
                         \int_{\Lambda_i}N_i\,d\sigma_\mu=-F_i. \tag{3.2}
\]

In particular,

\[
                         |F_i|\le r_i.                         \tag{3.3}
\]

When `F_i ne0`, put `u_i=-F_i/|F_i|`; when it vanishes, choose any unit
vector.  Set

\[
 \varepsilon_i=r_i-|F_i|,\qquad
 \varepsilon=\sum_i\varepsilon_i
              =2L-\sum_i|F_i|.                       \tag{3.4}
\]

**Theorem 3.1 (calibration Pythagoras).**  For every cell,

\[
 \boxed{
 \varepsilon_i={1\over2}\int_{\Lambda_i}|N_i-u_i|^2d\sigma_\mu.} \tag{3.5}
\]

If `F_i=0`, the right side is `r_i` for every constant unit `u_i`, so (3.5)
still holds.  Moreover, with `U_i=u_i u_i^T`,

\[
 \boxed{
 \sum_{i<j}\ell_{ij}\|U_i-U_j\|_{HS}^2\le8\varepsilon.}        \tag{3.6}
\]

**Proof.**  Expanding the square in (3.5) and using (3.2) gives
`2r_i-2|F_i|`.  On `Lambda_ij`, `N_j=-N_i`, and for unit vectors

\[
 \|u_i u_i^T-u_j u_j^T\|_{HS}^2
 \le2|u_i+u_j|^2
 \le4\{|u_i-N_i|^2+|u_j-N_j|^2\}.                  \tag{3.7}
\]

Integrate and sum over interfaces; (3.5) yields (3.6).  QED.

Equality `epsilon=0` says that every cell-side filling normal is constant,
and adjacent constants differ only by sign.  Thus all interfaces in one
connected filling component lie in parallel hyperplanes.  Notice what the
statement does **not** say: it does not align `u_i` with the exposed phase
normal.  In a one-sided log-affine tail the drift term cancels the exposed
flux, `F_i=0`, and no interior filling is required.  In a radial cap the
exposed normal integral instead calibrates a flat interior disk.  These two
possibilities must not be conflated.

## 4. A quantitative calibrated-flow trichotomy

The filling adjacency graph has vertex weights

\[
 d_i=r_i=\sum_j\ell_{ij},\qquad D=\sum_i d_i=2L,                 \tag{4.1}
\]

and edge weights `ell_ij`.  Its normalized spectral gap `lambda` is the
largest number such that for every Hilbert-valued family `H_i`,

\[
 \sum_{i<j}\ell_{ij}|H_i-H_j|^2
 \ge\lambda\sum_i d_i|H_i-\bar H|^2,
 \qquad \bar H=D^{-1}\sum_i d_iH_i.                 \tag{4.2}
\]

Vertices with `d_i=0` are deleted; they contribute neither filling nor any
term in (4.2).

**Theorem 4.1 (large filling: cut, parallel flow, or Plateau excess).**
Fix `lambda_0 in (0,1)`.  Every phase partition satisfies the following
exhaustive alternative.

1. `2L+E_tail<=1.4*10^-4 p`; this is the cheap component branch.
2. The cheap inequality fails and `lambda<lambda_0`.  Then the normalized
   graph has a bottleneck: for some nontrivial `S`,

   \[
   \ell(S,S^c)\le\sqrt{2\lambda_0}\,
        \min\{d(S),d(S^c)\}.                         \tag{4.3}
   \]

   The union `A_S=union_{i in S}A_i` is an actual binary component cut
   whose newly inserted boundary is exactly `ell(S,S^c)`.
3. The cheap inequality fails and `lambda>=lambda_0`.  Then

   \[
   {1\over D}\sum_i d_i\|U_i-\bar U\|_{HS}^2
   \le {8\varepsilon\over\lambda_0 D}.              \tag{4.4}
   \]

   Hence, after fixing any `eta>0`, this branch further divides into
   `epsilon<=eta D`, which is projective coherence with error
   `8eta/lambda_0`, and `epsilon>eta D`, which is a fixed noncalibrated
   Plateau excess.

Here (4.3) is the standard normalized-graph Cheeger sweep inequality.
Equation (4.4) follows immediately from (3.6) and (4.2).  In particular,
with `lambda_0=1/4`, if `epsilon<=10^-4 D`, then

\[
 {1\over D}\sum_i d_i\|U_i-\bar U\|_{HS}^2\le.0032.             \tag{4.5}
\]

Thus saturation plus expansion gives a robust parallel-flow certificate.
Conversely, obtaining a cut of absolute cost `7*10^-5p` from (4.3) requires
an absolute bound on the smaller vertex volume; a small normalized
conductance alone does not supply it.  Recursive sweeping can accumulate a
logarithmic number of cuts, so it cannot be inserted into the dimension-free
constant chain without a separate bounded-reuse theorem.

Most importantly, the Plateau-excess subbranch is not a perimeter saving.  The interfaces
are auxiliary surfaces inserted inside `A`; reducing their Plateau area does
not change `P_mu(A)`.  A theorem turning `epsilon` into a saving of the
*exposed* perimeter needs an additional coupling to the near-Cheeger first
variation.  Neither Gauss--Green nor max-flow/min-cut provides that coupling.

## 5. Quadratic phase capacity has the wrong implication

For the same binary trace define

\[
 \Phi_2(\Gamma_+,\Gamma_-)
 =\inf\left\{\int_A\rho|\nabla u|^2dx:
 u\in W^{1,2}(A),\ \operatorname {Tr}u|_{\Gamma_+}=1,
 \operatorname {Tr}u|_{\Gamma_-}=0\right\}.         \tag{5.1}
\]

The remaining boundary is Neumann.  Cauchy--Schwarz and relaxation give the
sharp universal relation

\[
 \boxed{
 \Phi_2(\Gamma_+,\Gamma_-)
 \ge {\Phi_1(\Gamma_+,\Gamma_-)^2\over\mu(A)}.}      \tag{5.2}
\]

Indeed every admissible `u` satisfies

\[
 \Phi_1\le\int_A\rho|\nabla u|
 \le\sqrt{\mu(A)}\left(\int_A\rho|\nabla u|^2\right)^{1/2}.     \tag{5.3}
\]

There is no reverse bound with a universal constant or with the heat scale
inserted automatically.

**Exact rectangle audit.**  Let

\[
 A=(-L,L)\times(-h,h),\qquad \rho=1,                \tag{5.4}
\]

put value one on the top side, zero on the bottom side, and leave the two
vertical sides free.  Every horizontal line is a separator and the constant
flow `Z=e_2` calibrates it, so

\[
                         \Phi_1=2L.                            \tag{5.5}
\]

The harmonic minimizer is `u(x,y)=(y+h)/(2h)`, hence

\[
                         \Phi_2={L\over h}.                    \tag{5.6}
\]

Since `|A|=4Lh`, equality holds in (5.2).  But

\[
                         {\Phi_2\over\Phi_1}={1\over2h},       \tag{5.7}
\]

which ranges from zero to infinity under scaling.  Therefore a large
interior filling does not yield a phase interpolation with controlled
Dirichlet cost.  It asserts the opposite: every interpolation costs at least
the amount in (5.2).  To turn this into a contradiction one would need an
independent *upper* energy bound at the same selected scale and a lower bound
strong enough at the audited constants.  The heat information estimate only
gives an upper budget of order `p/s`; (5.2), at the threshold (0.2), does not
match it without a new bound relating `mu(A)`, `s`, and the selected filling.

There is also a measure mismatch which is fatal before constants are even
compared.  The proved heat energy is

\[
 \int \eta|\nabla H|^2d\mu,
 \qquad \eta=\omega|\nabla F_0|,                    \tag{5.8}
\]

whereas (5.2) concerns the unweighted interior filling measure `dmu`.
Cauchy--Schwarz can compare them only through

\[
 \int_A|\nabla u|d\mu
 \le\left(\int_A\eta|\nabla u|^2d\mu\right)^{1/2}
     \left(\int_A\eta^{-1}d\mu\right)^{1/2}.         \tag{5.9}
\]

The central trace floor bounds an integral over each level; it gives no
bound on `int_A eta^{-1}`.  The lapse `|nabla F_0|^{-1}` may be arbitrarily
large near critical points.  Thus even the direction of (5.2) cannot be
inserted into the heat estimate.

For scale, suppose unrealistically that `eta>=c p` throughout the filling
region.  At the minimal large-fill threshold and `mu(A)<=1/2`, (5.2) gives

\[
 \Phi_2\ge 2(7\,10^{-5}p)^2=9.8\,10^{-9}p^2.        \tag{5.10}
\]

The heat derivative budget at `s=alpha K`, `alpha=10^-10`, becomes, after
using the Buser--Ledoux comparison `K asymp p^{-2}`, a coefficient of order
`alpha^{-1}p^3` in the `eta`-weighted energy.  The hypothetical lower bound
from (5.10) is only of order `10^-8 p^3`.  The coefficient gap is about
`10^18`, even before selector and trace-extension losses.  Hence ordinary
quadratic capacity cannot supply the audited `10^-4 p` finite saving; one
needs a selected scale-`s` capacity or a genuinely finite competitor.

## 6. Exact stress tests

### 6.1 Long rectangles and boxes

Give the top face of `[-L,L]times[-h,h]` one phase and the remaining exposed
boundary the other.  The two exterior trace junctions can be confined to
arbitrarily short bevel neighborhoods, while any interior separator has
length at least `2L`.  The vector flux of the top face is `(0,2L)`, so the
horizontal separator saturates (3.3) and `epsilon=0`.  This is the parallel
branch, not a phase-capacity saving.

In the box `prod_k[-L_k,L_k]`, a facet has flux equal to its `(n-1)`-area,
and a parallel cross-section calibrates its binary filling.  Assigning all
`2n` facets separately gives orthogonal calibrated binary flows.  There is
no single concurrent center theorem at the level of these flows.  For the
isotropic cube `[-sqrt3,sqrt3]^n`, the ordinary facet ridge area grows like
`n` times the facet area, so the exterior ridge-bevel branch, rather than a
large-fill inverse, is the correct charge for the full facet partition.

### 6.2 Balls and hemispheres

Let `A=B_R^n`, `rho=1`, and split its boundary into the hemispheres
`x_1>0` and `x_1<0`.  The equatorial disk is calibrated by the constant flow
`e_1`, hence

\[
 \Phi_1=\kappa_{n-1}R^{n-1}.                         \tag{6.1}
\]

The ridge is the equatorial sphere, of area
`(n-1)kappa_{n-1}R^{n-2}`, so

\[
 {\Phi_1\over\mathcal H^{n-2}(\text{ridge})}
 ={R\over n-1}.                                      \tag{6.2}
\]

The exposed normals are concurrent at the center, but the max flow and the
minimal filling are parallel.  Thus concurrence is a statement about the
*exposed normal lines*, not about the flow certificate.  Max-flow alone
cannot distinguish the radial and slab models.

### 6.3 A regular simplex

For a facet `F` of a simplex, prescribe value one on the relative interior
of `F` and zero on the other facets.  The constant vector in the facet-normal
direction gives the lower bound `Phi_1>=H^{n-1}(F)`.  Parallel sections,
tapered in a vanishing neighborhood of the facet ridges, give the reverse
inequality.  Hence

\[
                         \Phi_1=\mathcal H^{n-1}(F),             \tag{6.3}
\]

although the infimum need not be attained with the discontinuous full trace.
The facet graph of a regular simplex is complete and its dihedral ridges have
positive weighted capacity.  The simultaneous bevel therefore sees the
full facet partition before the large-filling branch is invoked.

### 6.4 Product exponentials

For

\[
 d\mu=e^{-\sum_{i=1}^m x_i}1_{\{x_i\ge0\}}dx,
 \qquad A=\{\max_i x_i\ge L\},\qquad q=e^{-L},        \tag{6.4}
\]

partition `A` by the largest coordinate.  On the interface between phases
`i` and `j`, `x_i=x_j=t>=L` and every other coordinate is at most `t`.
Therefore

\[
 \ell_{ij}=\sqrt2\int_L^\infty e^{-2t}
                    (1-e^{-t})^{m-2}dt
 \le {q^2\over\sqrt2}.                              \tag{6.5}
\]

It follows that this explicit partition satisfies

\[
 2\operatorname {Fill}
 \le {m(m-1)q^2\over\sqrt2}.                         \tag{6.6}
\]

The exposed perimeter is `P=mq(1-q)^{m-1}`.  If `alpha=mq` is small, then

\[
 {2\operatorname {Fill}\over P}
 \le{(m-1)q\over\sqrt2(1-q)^{m-1}}
 ={\alpha\over\sqrt2}+O(\alpha^2+\alpha/m).          \tag{6.7}
\]

Thus the rare product tail is genuinely in the cheap-filling branch.  At
fixed `alpha`, its pairwise incidence is fixed order and is seen by the
finite simultaneous bevel.  This agrees with, and sharpens, the product
stress test in `global_tube_amplification.md`.

## 7. A decisive max-flow rigidity counterexample

The next construction rules out an inference from a large max flow to
parallel/concurrent/orthogonal-radial geometry without using near-Cheeger
minimality.

Let `gamma:[0,L]->R^n` be an embedded `C^2` unit-speed curve with reach
greater than `4epsilon`.  Its tangent may visit any prescribed finite list
of unit directions.  Let `T_epsilon(gamma)` be its radius-`epsilon` tubular
neighborhood, cut by the normal disks at the two endpoints.  Put Dirichlet
trace zero and one on the two endpoint disks and leave the lateral boundary
free.

**Proposition 7.1 (turning tube).**  With
`A_{n-1}=kappa_{n-1}epsilon^{n-1}`, one has exactly

\[
                         \Phi_1=A_{n-1}.             \tag{7.1}
\]

Moreover equality is certified by a unit divergence-free flow whose
direction at every point of the normal disk over `gamma(s)` is
`gamma'(s)`.

**Proof.**  In a parallel normal frame, the tube map has the form

\[
 \Psi(s,z)=\gamma(s)+\sum_{a=1}^{n-1}z_a e_a(s),
 \qquad |z|<\epsilon.                                \tag{7.2}
\]

Choose a relatively parallel (Bishop) normal frame.  If
`kappa_a(s)=<gamma''(s),e_a(s)>`, then

\[
 \partial_s\Psi=(1-\textstyle\sum_a\kappa_a(s)z_a)\gamma'(s),
 \qquad \partial_{z_a}\Psi=e_a(s).                  \tag{7.3}
\]

The physical vector field `Z(Psi(s,z))=gamma'(s)` is the Piola transform
of the constant axial flux.  In the coordinates (7.2),

\[
 \operatorname {div}Z
 ={1\over J}\partial_s\left(J{1\over J}\right)=0,
 \qquad J=1-\sum_a\kappa_a(s)z_a.                   \tag{7.4}
\]

It has unit norm, is tangent to the lateral boundary, and has flux
`A_{n-1}` through every normal disk.  Theorem 2.1 gives the lower bound.
Any normal disk is an admissible separator and has Euclidean area exactly
`A_{n-1}`, giving the upper bound.  QED.

By choosing a long curve with small curvature, its tangent can successively
approach `e_1,...,e_n` while (7.1) remains exact.
Neither the flow nor its Smirnov arcs are close to one line, one center, or
an orthogonal radial block decomposition.  The example can be placed inside
a region where a smooth log-concave density varies by `1+o(1)`, so weighting
does not restore the claimed rigidity.

The tube has perimeter-to-volume ratio of order `1/epsilon`, and hence is
not asserted to be near-Cheeger.  That is precisely the point: max-flow,
finite perimeter, high directional rank, and log-concavity of the ambient
density do not contain the missing rigidity.  A valid large-fill theorem
must quantitatively use the coarea near-minimality deficit or the
constant-mean-curvature/tube information.  Calling the max-flow certificate
itself “concurrent” would be false.

### 7.1 Relation to the mass-preserving terminal-needle obstruction

The flow decomposition gives one exact inequality which is superficially
similar to the terminal-needle target in
`mass_preserving_matrix_localization.md`.  Normalize the acyclic Smirnov
decomposition by its endpoint flux `Phi_1`, and let `ell_gamma` be the
length of a sampled flow arc.  From (2.5),

\[
 \mathbb E_{flow}\ell_\gamma\le{\mu(A)\over\Phi_1},
 \qquad
 \mathbb E_{flow}{1\over\ell_\gamma}
 \ge{\Phi_1\over\mu(A)}.                            \tag{7.5}
\]

If a probability supported on an interval of length `ell_gamma` is placed
on each arc, its one-dimensional standard deviation is at most
`ell_gamma/2`; hence

\[
 \mathbb E_{flow}{1\over\sigma_\gamma}
 \ge {2\Phi_1\over\mu(A)}.                          \tag{7.6}
\]

This does **not** prove the terminal-needle lemma.  The adaptive
mass-preserving localization needles are affine lines with a particular
path law; they are not the arbitrary max-flow arcs, and no coupling between
the two decompositions is known.  The turning tube calibrates the mismatch:
its flow arcs have length `L`, `mu(A)=A_{n-1}L`, and (7.5) is an equality.
Thus a long high-turning flow gives only `1/L`, precisely the bad long-needle
scale.  Isotropy and the concurrent-barycenter constraint of the adaptive
localization must supply the missing improvement; max-flow does not.

## 8. The exact missing lemma

The proved statements reduce the large-filling branch to the following
formal target.  It is written so that it can be audited without geometric
prose.

**Required large-fill inverse.**  There must be universal constants
`c_0,c_1>0` with the following property.  Let `mu` be isotropic log-concave,
let `A` be one of the central physical heat levels, and suppose

\[
 P_\mu(A)-I(\mu(A))\le d,
 \qquad \operatorname {tr}M_{sel}\ge.0048495p,
 \qquad {\operatorname {tr}M_{sel}\over\|M_{sel}\|_{op}}\ge17.86. \tag{8.1}
\]

For every exterior phase cut whose admissible bevel conductance is below
the expander threshold, one must prove at least one of:

1. an actual Caccioppoli componentization with

   \[
   2\operatorname {Fill}+E_{tail}\le1.4\,10^{-4}p;                \tag{8.2}
   \]

2. a competitor `A'` with `mu(A')=mu(A)` and

   \[
   P_\mu(A')\le P_\mu(A)-1.4\,10^{-4}p;                           \tag{8.3}
   \]

3. a quantitatively stated parallel, concurrent, or orthogonal-radial
   normal-line approximation whose error is small enough for the existing
   long-ray covariance or translated-thin-shell argument.

Theorems 2.1--4.1 show that (8.3) cannot be replaced by “the Plateau excess
is large,” and Proposition 7.1 shows that item 3 cannot be concluded from a
max-flow decomposition.  A sufficient new hypothesis would be an estimate
of the explicit form

\[
 \varepsilon
 +\int |\operatorname {turn}Z|^2\,d|Z|
 \le C\{d+\mathcal S_{bevel}\},                    \tag{8.4}
\]

for an optimal filling flow `Z`, together with a rigidity theorem for the
small right-hand side.  No present coarea, Fisher-information, or `BV`
identity proves (8.4).  Establishing (8.4), or a different coupling with the
same consequences, is the genuine large-fill problem.

The audit conclusion is therefore sharp: cheap filling is completely
quantified by (0.2); calibrated large filling is reduced to graph coherence
by (4.4); and the remaining noncalibrated/turning branch cannot be promoted
to a perimeter saving without a new near-minimality-to-flow estimate.
