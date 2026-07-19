# Audit and product inverse for the global tube argument

## 0. Verdict

The killed normal-tube identities in `global_tube_amplification.md` are
correct under the hypotheses actually stated there: a bounded smooth convex
support, a `C^2` convex potential, and an everywhere regular compact
free-boundary isoperimetric hypersurface.  In that setting the argument gives
an exact finite-distance accounting of curvature, support contact, focal
collapse, and cut-locus collisions.

There are, however, two logically separate transfer problems before that
theorem applies to the physical packet constructed from `F_0=T_s1_S`.

1. The packet is carried by an *integral of near-minimizing level sets* of
   `F_0`.  The tube theorem starts from one *exact isoperimetric region*.
   The coarea estimate supplies neither stationarity of an `F_0` level nor a
   constant weighted mean curvature.  Replacing a level by a profile
   minimizer does not preserve its normal matrix.
2. In dimensions at least eight an exact isoperimetric boundary may have a
   singular set.  The present tube proof assumes that the entire free
   boundary is `C^2`; it does not yet account for normal cones based at the
   singular set.  Smooth approximation of the density and support does not
   by itself remove this issue because minimizers and their normal measures
   also have to converge.

Thus (0.6) of that report must not presently be applied to the `F_0` levels.
An Ekeland replacement gives quasiminimality, not the constant first
variation needed for the factor `exp(lambda t)`.  This is recorded precisely
in Section 5 below.

On the other hand, the proposed high-rank inverse can be completed in the
ideal branch in which small ridge conductance has already cut the phases into
genuine one-sided log-affine components.  The result is stronger than a
generic level-law Poincare estimate and does not assume that the level law is
unimodal.  At the audited constants it proves:

\[
 \text{cross-level variance }\ge {8\over17}
 \quad\Longrightarrow\quad
 R_{\rm cut}+E_{\rm tail}>1.4\,10^{-4}p.              \tag{0.1}
\]

where `R_cut` is the perimeter needed to separate the phase components and
`E_tail` is the perimeter-unit error in the log-affine identity
`volume=boundary flux/slope`.  The numerical constant in (0.1) is proved
below with a slightly conservative value `1.4*10^-4`.  It is valid for the
central selected trace `A>=.0048495p`, scalar coarea deficit
`D<=6.02*10^-5p`, floor `omega_0=10^-5`, and heat scale
`s=10^-10 K`.

Consequently a surviving high-rank branch is not merely called “product.”
It has to pay one of two explicit geometric defects.  The remaining missing
lemma is to dominate `R_cut` by a bounded-reuse bevel saving and
`E_tail` by killed tube flux for the *same near-minimizing physical levels*.

## 1. Audit of the killed tube theorem

Let the notation of `global_tube_amplification.md` be in force.  Before the
first support, focal, or cut time, the map

\[
                         F(x,t)=x+tN_x
\]

has weighted Jacobian

\[
 j_x(t)=e^{-V(x+tN_x)+V(x)}\det(I+tS_x).
\]

The sign convention `S u=D_uN` is consistent with the outer parallel map.
Differentiating twice gives

\[
 (\log j_x)''(t)=-\operatorname {tr}
   \{S_x^2(I+tS_x)^{-2}\}-\nabla^2V(x+tN_x)[N_x,N_x]\le0.
\]

Hence `j_x(t)=exp(lambda t-D_x(t))`, with `D_x>=0`, exactly as claimed.
The area formula on the injectivity region and the coarea formula for the
distance function give

\[
 v'(t)=e^{\lambda t}R(t)
\]

for almost every `t`.  Killing a ray at its first exceptional time is the
correct convention: no Taylor remainder is hidden at a support contact,
focal point, or collision.

The profile comparison also has the stated constants.  If
`c=I(b)/b`, `s_0=min(v_0-a,b-v_0)`, and
`eta=delta b/s_0`, concavity gives

\[
                         |\lambda-c|\le\eta c.
\]

The inequalities `v'>=cv`, `R<=P_0`, and
`P_0/v_0<=(1+delta)c` then give (0.5), while

\[
 {R(T_b-)\over P_0}
 \ge {1\over1+\delta}\,\beta^{1-\lambda/c}
 \ge {1\over1+\delta}\,\beta^{-\eta}
\]

gives (0.6).  The covariance bound uses
`inf_a int_0^T(A+Bt-a)^2dt=B^2T^3/12`; its constant is also correct.

The audit therefore finds no algebraic loss in the regular theorem.  Its
limitations are hypotheses and transfer, not the finite-distance constant
chain.

## 2. Approximate separated-component charge

The following lemma makes hard-support and cutting errors explicit.  All
quantities have perimeter units except the volumes and slopes.

**Lemma 2.1 (component excess with cut and tail errors).**  Let `psi>0`.
For `i=1,...,m`, let `a_i>0`, `v_i in(0,1/2]`, `kappa_i>0`, and
`r_i,e_i>=0`.  Assume there are disjoint finite-perimeter component sets
`E_i` such that

\[
 P_\mu(E_i)\le a_i+r_i,\qquad \mu(E_i)=v_i,
 \qquad \left|v_i-{a_i\over\kappa_i}\right|
       \le {e_i\over\psi}.                            \tag{2.1}
\]

Assume also that the union has additive principal perimeter and volume,
so that

\[
 D:=\sum_i(a_i-\psi v_i)\ge0.                         \tag{2.2}
\]

Put `R=sum r_i` and `E=sum e_i`.  Then

\[
 \sum_i\left[a_i\left(1-{\psi\over\kappa_i}\right)
       \right]_+\le D+R+E,                            \tag{2.3}
\]

and, for every `tau>0`,

\[
 \sum_{\{\kappa_i\ge(1+\tau)\psi\}}a_i
 \le {1+\tau\over\tau}(D+R+E).                      \tag{2.4}
\]

Moreover,

\[
 \sum_{\{\kappa_i\le\psi/2\}}a_i\le R+E.          \tag{2.5}
\]

**Proof.**  Set

\[
 x_i=a_i(1-\psi/\kappa_i),\qquad y_i=a_i-\psi v_i.
\]

The tail estimate in (2.1) gives `|x_i-y_i|<=e_i`.  The Cheeger
inequality applied to `E_i` gives

\[
                         y_i\ge-r_i.                 \tag{2.6}
\]

Since `sum y_i=D`,

\[
 \sum_i(y_i)_+=D+\sum_i(y_i)_-\le D+R.
\]

Thus `sum(x_i)_+<=D+R+E`, proving (2.3).  On
`kappa_i>=(1+tau)psi`, one has
`x_i>=tau a_i/(1+tau)`, which proves (2.4).  Finally, on
`kappa_i<=psi/2`, `x_i<=-a_i`; (2.6) and
`|x_i-y_i|<=e_i` give `a_i<=r_i+e_i`.  Summation proves
(2.5).  QED.

In the exact branch, `r_i=e_i=0`, (2.1) says
`v_i=a_i/kappa_i`, and Cheeger itself gives `kappa_i>=psi`.
This recovers the exact component identity

\[
 D=\sum_i a_i(1-\psi/\kappa_i).                       \tag{2.7}
\]

The error `r_i` includes the seam inserted when a phase packet is cut away
from its neighbors.  The error `e_i` includes support truncation, focal
loss, collision loss, and any non-affinity of the normal density which
prevents the full one-sided tail integral from equaling `a_i/kappa_i`.
Thus (2.1) is the exact location at which hard support enters the inverse.

### 2.1 Exterior ridge conductance is not componentization cost

There is an exact `BV` distinction between a small exterior ridge and a
cheap component cut.  Let `E_1,...,E_m` be a Caccioppoli partition of a
finite-perimeter set `A`, modulo null sets, and use relative weighted
perimeter in the affine support.  The partition identity is

\[
 \boxed{
 \sum_iP_\mu(E_i)
 =P_\mu(A)+2\sum_{i<j}
   \int_{A^{(1)}\cap\partial^*E_i\cap\partial^*E_j}e^{-V},
       d\mathcal H^{n-1}.}                            \tag{2.8}
\]

It follows by applying the ordinary reduced-boundary partition identity
and then integrating the common nonnegative weight; truncation handles an
unbounded or merely lower-semicontinuous convex potential.  Denote the sum
over `i<j` in (2.8) by `Fill(A;{E_i})`.  If `a_i` is the exterior principal
boundary assigned to phase `i` and `r_i` is its newly inserted interior
boundary, then

\[
 \sum_i a_i=P_\mu(A),\qquad
 \sum_i r_i=2\operatorname {Fill}(A;\{E_i\}).         \tag{2.9}
\]

Thus the `R` in Lemma 2.1 is twice an *interior codimension-one filling*.
It is not the codimension-two measure of the ridges where exterior phase
patches meet.

This distinction is sharp.  In a planar unit disk, assign the upper and
lower semicircles to two phases.  Their exterior ridge consists of two
points and has zero one-dimensional perimeter measure, but every partition
realizing these traces contains a curve joining the two endpoints; its
length is at least the diameter two.  Equality is attained by the diameter.
Formula (2.8) therefore adds four to the sum of the two component
perimeters.  The normals in this example are exactly concurrent at the
center.  A rectangle gives the analogous phenomenon with parallel top and
bottom phases and a long interior filling segment.  Higher-dimensional
balls and boxes give radial and orthogonal-block versions.

Consequently “small ridge conductance implies cheap componentization” is
false.  The geometric inverse has at least three branches:

1. large exterior ridge expansion, which is available to a simultaneous
   bevel;
2. small interior filling, in which Lemmas 2.1 and 3.1 synchronize the
   component slopes; and
3. large minimal interior filling despite small exterior ridge measure.

The model geometries of the third branch are parallel slabs, concurrent
radial cells, and orthogonal radial blocks.  Exact concurrence is controlled
by translated thin shell, and a coherent parallel packet is controlled by
the long-ray covariance inequality.  However, (2.8) alone does not prove
that a general high-rank large-filling partition is close to either model:
a box already shows that concurrence is too narrow a conclusion.  A robust
parallel/concurrent/orthogonal-radial filling theorem is an additional
missing inverse statement.  The product conclusion cannot be asserted
before that third branch is discharged.

### 2.2 A corrected flux which forces interior filling

There is a useful exact invariant inside the large-filling branch.  Assume
temporarily that `Omega` and `V` are smooth.  Let `E_i subset A` be one cell
of a partition and decompose its relative boundary into its exposed part
`Gamma_i subset partial^*A cap Omega`, its interior filling interface
`Lambda_i subset A^(1)`, and its support contact `C_i subset partial Omega`.
Weighted Gauss--Green, applied to each constant ambient vector, gives the
vector identity

\[
 \int_{\Gamma_i}N_A\,d\sigma_\mu
 +\int_{\Lambda_i}N_i\,d\sigma_\mu
 +Z^{-1}\int_{C_i}n_\Omega e^{-V},d\mathcal H^{n-1}
 =-\int_{E_i}\nabla V\,d\mu.                          \tag{2.10}
\]

Consequently

\[
 \boxed{
 P_\mu(\Lambda_i)\ge |\mathcal F_i|,\qquad
 \mathcal F_i:=\int_{\Gamma_i}N_A\,d\sigma_\mu
 +\int_{E_i}\nabla V\,d\mu
 +Z^{-1}\int_{C_i}n_\Omega e^{-V},d\mathcal H^{n-1}.} \tag{2.11}
\]

The same formula holds for finite-perimeter cells by approximation and the
Gauss--Green theorem for divergence-measure fields.  For a nonsmooth convex
potential, (2.10) is read with the distributional gradient of `e^{-V}`;
for a nonsmooth support the last integral is the normal trace on its reduced
boundary.  Thus hard-support contact is part of the invariant rather than a
discarded error.

In a full-support log-affine region, `V(x)=<b,x>+constant`, so

\[
                    \mathcal F_i
       =\int_{\Gamma_i}N_A\,d\sigma_\mu+b\,\mu(E_i). \tag{2.12}
\]

For a genuine one-sided exponential tail, the exposed planar flux is
exactly `-b mu(E_i)`, and `mathcal F_i=0`; no interior filling is needed.
This is the product/log-affine equality.  For the upper semicircle of a disk
lying strictly inside a larger constant-density chart, the drift and contact terms vanish and the exposed normal
integral has norm two, forcing the diameter filling with its sharp length.
The same computation calibrates the long filling segment between the top
and bottom phases of a rectangle.

Formula (2.11) refines the third branch into two subcases.  If a substantial
part of the filling is forced by `sum_i|mathcal F_i|`, there is a coherent
vector-flow certificate linking each exposed phase to drift or support
contact.  Concurrent radial cells, parallel slabs, and orthogonal blocks
saturate this certificate.  If instead

\[
 \sum_iP_\mu(\Lambda_i)-\sum_i|\mathcal F_i|
\]

is large, the filling interfaces themselves have a large noncalibrated
area.  Turning that excess into a finite perimeter saving is a more focused
geometric target than inferring product structure from exterior ridges.
No estimate currently proves that either subcase is close to the listed
models; (2.11) identifies the additional invariant without assuming the
missing rigidity theorem.

There is an exact equality statement behind this refinement.  From (2.10),

\[
                         \int_{\Lambda_i}N_i\,d\sigma_\mu
                              =-\mathcal F_i.
\]

If `mathcal F_i!=0`, put `u_i=-mathcal F_i/|mathcal F_i|`; if it is zero,
choose an arbitrary unit vector.  Then

\[
 \boxed{
 P_\mu(\Lambda_i)-|\mathcal F_i|
 ={1\over2}\int_{\Lambda_i}|N_i-u_i|^2\,d\sigma_\mu.} \tag{2.13}
\]

For `mathcal F_i=0`, the right side equals `P_mu(Lambda_i)` because the
mean normal vanishes, so the same formula is valid.  Thus equality in the
corrected-flux lower bound forces the entire filling boundary of a cell to
have one normal direction; its rectifiable components lie in parallel
hyperplanes up to null sets.

The compatibility across cells is also quantitative.  On a shared
interface `Lambda_ij`, the measure-theoretic normals satisfy
`N_j=-N_i`.  Hence

\[
 \sigma_\mu(\Lambda_{ij})|u_i+u_j|^2
 \le2\int_{\Lambda_{ij}}
       \{|N_i-u_i|^2+|N_j-u_j|^2\}\,d\sigma_\mu.       \tag{2.14}
\]

Consequently a filling adjacency graph with a numerical spectral gap and
small total excess in (2.13) synchronizes all cell directions to one
projective line; this is the parallel-slab branch.  A graph bottleneck
returns a candidate cheap cluster cut.  What remains outside those two
alternatives is a fixed noncalibrated filling excess, the natural place to
seek a concurrent/orthogonal-radial finite competitor.  Equations
(2.13)--(2.14) prove the algebraic part of that refined trichotomy; they do
not yet convert the noncalibrated excess into perimeter saving.

## 3. Softmax exclusion without a level-law spectral gap

We next use the analytic selector floor.  The statement permits an
arbitrary multimodal level law.

Let `Z_i` be Hilbert-space vectors with `|Z_i|<=1`; in the application
`Z_i=n_i n_i^T` and the Hilbert norm is Hilbert--Schmidt.  Put

\[
 c_i=\sqrt{s}\,\kappa_i,
\]

let `J` be an interval (possibly a half-line or all of `R`), and let
`w:J -> [omega_0,1]`.  Define the selected joint phase law by

\[
 d\nu(i,z)={A_i1_J(z)w(z)\varphi(z-\sqrt2c_i)\,dz
                 \over A_{sel}},\qquad
 A_{sel}=\sum_iA_i\int_J w(z)\varphi(z-\sqrt2c_i)dz. \tag{3.1}
\]

Write

\[
 \bar A_i=A_i\int_J\varphi(z-\sqrt2c_i)dz.           \tag{3.1a}
\]

Thus `bar A_i` is the unselected coarea flux of phase `i` on the chosen
level interval, while its selected flux is at most `bar A_i` and at least
`omega_0 bar A_i`.

This is exactly the log-affine heat formula because

\[
 \varphi(z)e^{-c_i^2+\sqrt2c_i z}
                         =\varphi(z-\sqrt2c_i).       \tag{3.2}
\]

Let `rho` be the `z` marginal and

\[
 Q(z)=E_\nu[Z_i\mid z].                               \tag{3.3}
\]

**Lemma 3.1 (floored softmax exclusion).**  Suppose a set `G` of good
phases satisfies

\[
 c_i\in[c_-,c_+]\quad(i\in G),\qquad
 \Delta=c_+-c_-<\sqrt2.                               \tag{3.4}
\]

If `q=nu{i notin G}`, then

\[
 \boxed{
 \sqrt{\operatorname {Var}_\rho Q}
 \le {\Delta\over
       \sqrt{2\omega_0(1-\Delta^2/2)}}+2\sqrt q.}    \tag{3.5}
\]

**Proof.**  Conditional on being good, write `Q_G(z)` for the projector
mean.  Differentiating its softmax weights gives

\[
 Q_G'(z)=\sqrt2\operatorname {Cov}(c_i,Z_i\mid z,G).
\]

Cauchy--Schwarz, `Var(c_i)<=Delta^2/4`, and
`E|Z_i-EZ_i|^2<=1` imply

\[
                         |Q_G'(z)|\le{\Delta\over\sqrt2}. \tag{3.6}
\]

Remove the factor `w` temporarily and normalize the good `z` marginal on
`J`.  Its density is proportional to

\[
 h(z)=1_J(z)\sum_{i\in G}A_i\varphi(z-\sqrt2c_i).
\]

In the interior of `J`, direct differentiation gives

\[
 (\log h)''(z)=-1+2\operatorname {Var}(c_i\mid z,G)
       \le-1+{\Delta^2\over2}.                       \tag{3.7}
\]

The one-dimensional Brascamp--Lieb inequality, valid on an interval with
the natural Neumann boundary condition by approximation, therefore gives

\[
                    \operatorname {Var}_h z
       \le {1\over1-\Delta^2/2}.                     \tag{3.8}
\]

Multiplication by `w` costs only its oscillation.  If `m_h` is the mean
under the normalized density `h`, then

\[
 \operatorname {Var}(z\mid G)
 \le E[(z-m_h)^2\mid G]
 \le {1\over\omega_0(1-\Delta^2/2)}.                \tag{3.9}
\]

It follows from (3.6), using `Q_G(m_h)` as reference, that

\[
 \left(E[|Q_G(z)-Q_G(m_h)|^2\mid G]\right)^{1/2}
 \le {\Delta\over\sqrt{2\omega_0(1-\Delta^2/2)}}
 =:d.                                                 \tag{3.10}
\]

Let `beta(z)=P(i notin G\mid z)` and let
`Q_G^{sub}(z)=(1-beta(z))Q_G(z)`.  Since `|Z_i|<=1`,

\[
 \|Q-Q_G^{sub}\|_{L^2(\rho)}\le\sqrt q.             \tag{3.11}
\]

Also, comparison with the fixed vector `Q_G(m_h)` and Minkowski's
inequality give

\[
 \sqrt{\operatorname {Var}_\rho Q_G^{sub}}
 \le\|Q_G^{sub}-Q_G(m_h)\|_{L^2(\rho)}
 \le d+\sqrt q.                                      \tag{3.12}
\]

The square root of Hilbert-valued variance is a seminorm.  Combining
(3.11)--(3.12) proves (3.5).  QED.

Unlike a Poincare argument on the full `z` marginal, Lemma 3.1 remains
valid when exceptional separated slopes make that marginal strongly
multimodal.  It does not require `J` to be bounded.  The exact good-slope
Gaussian mixture is strongly log-concave because its means lie in a narrow
interval; the selector floor transfers its variance bound through (3.9).

## 4. Numerical inverse theorem

We now combine Lemmas 2.1 and 3.1.  Assume

\[
 \alpha=10^{-10},\quad s=\alpha K,\quad
 \psi\sqrt K\le2,\quad \omega_0=10^{-5},             \tag{4.1}
\]

and the audited central selected flux and scalar deficit bounds

\[
 A_{sel}\ge .0048495p,\qquad D\le .0000602p.          \tag{4.2}
\]

Suppose Lemma 2.1 applies level by level and, after integration over `J`,
its principal areas `a_i` are the unselected phase fluxes `bar A_i` in
(3.1a), with the integrated errors denoted by the same letters `D,R,E`.
This is the exact separated log-affine-core model with explicit cut and
tail errors.  Take

\[
 G=\{i:\psi/2<\kappa_i<2\psi\},\qquad J=R+E.          \tag{4.3}
\]

Equations (2.4), with `tau=1`, and (2.5) imply

\[
 \sum_{i\notin G}\bar A_i\le2D+3J.                    \tag{4.4}
\]

Since `w<=1`, the selected bad fraction satisfies

\[
                         q\le{2D+3J\over A_{sel}}.    \tag{4.5}
\]

On `G`, (4.1) gives

\[
 \Delta\le{3\over2}\sqrt{s}\psi\le3\sqrt\alpha
       =3\,10^{-5},\qquad
 |c_i|\le2\sqrt{s}\psi\le4\sqrt\alpha=4\,10^{-5}. \tag{4.6}
\]

Thus the first term in (3.5) is less than

\[
 d_0={3\,10^{-5}\over
       \sqrt{2\,10^{-5}(1-4.5\,10^{-10})}}<.006709. \tag{4.7}
\]

If `J<=1.4*10^-4p`, then

\[
 q\le{2(.0000602)+3(.00014)\over .0048495}<.11145.   \tag{4.8}
\]

Consequently

\[
 \operatorname {Var}_\rho Q
 \le\{.006709+2\sqrt{.11145}\}^2
 <.4557<{8\over17}.                                  \tag{4.9}
\]

This proves (0.1).  The exact branch `R=E=0` has the much stronger bound
`Var Q<.111` at the same constants.

The constants have useful slack.  Solving (3.5) at equality shows that the
largest permitted combined cut-and-tail error is approximately
`1.46*10^-4p`; (4.9) uses `1.4*10^-4p` to avoid reliance on rounded last
digits.

### 4.1 Product check

For independent one-sided exponentials, every coordinate tail has
`kappa_i=1` and `v_i=a_i`, so all dimensionless slopes are equal.  Hence
the cross-level softmax variance is exactly zero.  A union of coordinate
tails reduces its perimeter through overlaps; those overlaps are precisely
nonadditivity of the component perimeter and are charged to `R_cut` in the
inverse theorem.  At half volume their ridge incidence is a fixed fraction,
and the simultaneous finite bevel in
`fixed_scale_physical_splicing.md` exploits it.  At very small volume their
ridge fraction is small, but if the ambient Cheeger constant is appreciably
below one, (2.7) assigns a fixed scalar excess to every separated coordinate
tail.  Thus the product model passes every branch of the dichotomy.

The same check shows why a high-rank ancillary product factor is not enough.
For `mu=nu tensor Exp^m`, exponential facets can have arbitrary normal rank,
while `C_P(mu)=max(C_P(nu),1)`.  If `nu` is the bottleneck, those facets have
a fixed Cheeger excess and cannot carry the audited near-minimizing packet.
The missing invariant is therefore not normal rank but the componentwise
Cheeger excess together with the cut-and-tail defect `J`.

## 5. The unresolved transfer to the physical heat levels

The inverse theorem above is quantitative, but its geometric hypotheses are
not consequences of the current coarea estimate.  This is an independent
load-bearing gap.

For a regular level `A_r={F_0>r}`, the proof knows only

\[
 \int_0^1\{P_\mu(A_r)-I(\mu(A_r))\}\,dr
                         \le6.02\,10^{-5}p.           \tag{5.1}
\]

It does not know that `A_r` is stationary.  Ekeland's variational principle,
applied to perimeter under a volume constraint with symmetric-difference
metric, can produce a `Lambda`-quasiminimizer close in measure to `A_r`.
Its first variation satisfies an inequality with error `Lambda`; it does
not have one constant weighted mean curvature.  Along its normal rays the
Jacobian has base-dependent initial slope, so the factorization
`j_x(t)=exp(lambda t-D_x(t))` used in Theorem 2.1 no longer yields one
monotone normalized flux.

Nor may one replace `A_r` by an exact profile minimizer of the same volume:
the inequality `P(E_r)<=P(A_r)` contains no quantitative control of

\[
 \int_{\partial^*A_r}\omega nn^T\,d\sigma_\mu
 \quad\hbox{versus}\quad
 \int_{\partial^*E_r}nn^T\,d\sigma_\mu.              \tag{5.2}
\]

Strict `BV` convergence would preserve (5.2) in a limit by Reshetnyak's
theorem, but (5.1) is a fixed nonzero budget, not a sequence whose deficit
tends to zero at the chosen heat scale.  In addition, no dimension-free
modulus of strict `BV` stability is available here.

A complete use of the tube inverse therefore requires one of the following
new statements, with the same numerical constants.

1. A near-stationary killed-tube theorem for the actual `F_0` levels which
   permits a base-dependent mean-curvature error and charges its integral to
   (5.1), while preserving the selected normal matrix.
2. A quantitative Ekeland/profile replacement theorem which preserves at
   least the `.0048495p` selected normal trace and its cross-level variance.
3. A direct physical cut lemma producing component sets satisfying (2.1)
   and bounding `R_cut+E_tail` by the already available finite-splice and
   coarea budgets, without passing through exact minimizers.

None of these statements follows from ordinary compactness, and assuming
one would be a new conjecture-strength bridge.  The rigorous conclusion of
this note is the explicit inverse (0.1) *after* that bridge, together with
the exact identification of the error invariant the bridge must control.

### 5.1 Exact mixture assumptions and selector-floor audit

For clarity, Lemma 3.1 uses more than flat normals at one level.  The exact
Gaussian mixture (3.1) follows from a physical phase packet only when all of
the following hold on the same level interval.

1. Each phase has one persistent normal `n_i` and one persistent weighted
   normal slope `kappa_i`; phase creation, merging, and relabelling are
   absent or charged to the cut error.
2. Along the normal coordinate the density is log-affine over the Gaussian
   reach used by the heat kernel.  Departure from log-affinity and loss of
   reach are charged to `E_tail` or to a separately estimated Gaussian-tail
   error.
3. The tangential phase coefficient is independent of the profile
   coordinate, giving the fixed number `A_i` in (3.1).
4. The selector is a *common* multiplier `w(z)` for every phase.  In the
   exact log-affine core computation this is true: `r_G/R` is a universal
   function of the one-dimensional profile coordinate, independent of
   `kappa_i`, and

   \[
              w(z)=10^{-5}+(1-10^{-5}){r_G\over R}
   \]

   lies in `[10^-5,1]`.
5. The unselected phase fluxes `bar A_i` admit actual component sets whose
   volumes satisfy (2.1).  This is the interior-filling and hard-support
   hypothesis, not a consequence of the softmax formula.

The fourth item is essential.  The global analytic floor by itself gives
`10^-5<=omega<=1`, but on approximate cores it does not say that the
phase restrictions are the same function of `z`.  If arbitrary
phase-dependent factors `w_i(z)` in this range are allowed, they can rotate
the conditional phase weights even when every `c_i` is equal; the
strong-logconcavity calculation (3.7) then fails.  Analyticity alone gives
no uniform derivative bound.  A quantitative core approximation must
therefore show common-selector error is small in selected joint mass, or
charge it to the coarea/finite-splice budget.

The uses of the floor and ceiling in the proved model are otherwise exact:

* `w<=1` turns the unselected exceptional flux in (4.4) into the selected
  exceptional probability bound (4.5);
* `w>=omega_0` costs precisely `omega_0^{-1}` in (3.9), and hence
  `omega_0^{-1/2}` in (3.5);
* no derivative of `w` is taken.

If the physical packet differs from an exact common-selector mixture by
joint total mass `zeta`, the proof remains valid after replacing `q` in
(3.5) by `q+zeta`: delete the exceptional mass and use the same
`L^2` subprobability estimate as in (3.11).  Thus the numerical theorem
allows

\[
 {2D+3J\over A_{sel}}+\zeta<
 {1\over4}\left(\sqrt{8/17}-.006709\right)^2
 <.1154.                                               \tag{5.3}
\]

This is a formal error target, not a proved approximation estimate for the
physical heat levels.
