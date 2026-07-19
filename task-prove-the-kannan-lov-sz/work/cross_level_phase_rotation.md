# Cross-level phase rotation: exact rigidity, lattice straightening, and the focal obstruction

## 0. Verdict

Let the physical coarea submeasure from
`fixed_scale_physical_splicing.md` be

\[
 M=\int_0^1M_rdr,qquad
 M_r=\int_{\partial^*A_r}\omega nn^Te^{-V}d\mathcal H^{n-1},
 \quad A_r=\{F_0>r\}.                                \tag{0.1}
\]

We use the floored analytic selector
`omega=10^{-5}+(1-10^{-5})r_G/R` from Section 1.4 of that report.  It is
globally real analytic, satisfies `10^{-5}<=omega<=1`, has `trM>.0051p`,
and has effective rank larger than `18.8`; only the weaker rank `17` is used
below.

The high effective rank of `M` has an exact two-way decomposition.  Either
there is substantial angular dispersion on individual levels, or the
normalized normal projector rotates substantially between levels.  The
second case cannot be dismissed by selecting two normals on the same level.

Three rigorous facts are obtained here.

1. Exact rank-one planar patches at different levels of the entire analytic
   heat function `F_0=T_s1_S` cannot rotate.  Analytic continuation makes
   each patch a whole affine hyperplane, and two nonparallel level
   hyperplanes intersect.
2. Finite union/intersection sorting is an exact null-invariant
   straightening operation.  It never increases the sum of perimeters.  Its
   only escape is precisely a linear isoperimetric profile, the branch
   analyzed in `profile_linearity_separation.md`.
3. A fully explicit rotating pencil for two one-sided exponentials has
   rank-one normals on every level and a fixed cross-level rotation, but it
   pays a fixed scalar perimeter deficit.  Rotation occurs around the hard
   support vertex.  This is the canonical focal/contact obstruction.

The quantitative closing step is still missing.  The selected physical
weight may concentrate on small nearly planar pieces and rotate through
levels where its amplitude is negligible.  A dimension-free quantitative
unique-continuation or free-boundary rigidity theorem is needed to turn the
exact analytic statement into the required `10^{-4}p` saving.

## 1. Exact same-level/cross-level variance decomposition

Put

\[
 A=trM,qquad d\pi(r)={trM_r\over A}dr,qquad
 Q_r={M_r\over trM_r},qquad Q={M\over A}=E_\pi Q_r,              \tag{1.1}
\]

with arbitrary `Q_r` on the null set where `trM_r=0`.  Every `Q_r` is
positive semidefinite with trace one.  Therefore

\[
 \boxed{
 1-tr(Q^2)=E_\pi\{1-tr(Q_r^2)\}
              +E_\pi\|Q_r-Q\|_{HS}^2.}             \tag{1.2}
\]

Indeed `E||Q_r-Q||^2=Etr(Q_r^2)-tr(Q^2)`.

The physical transfer gives `trQ=1` and `||Q||op<1/17`, so
`tr(Q^2)<=||Q||op trQ<1/17`.  Consequently at least one of

\[
 \boxed{
 E_\pi\{1-tr(Q_r^2)\}\ge {8\over17},
 \qquad
 E_\pi\|Q_r-Q\|_{HS}^2\ge {8\over17}}              \tag{1.3}
\]

holds.  The first is the same-level branch.  The second is the cross-level
branch studied below.  Notice that the second alternative is compatible
with `Q_r=u(r)u(r)^T` for every `r`.

### 1.1 The exact one-dimensional phase energy

There is a canonical way to quantify the amount of cross-level motion,
although the present estimates do not control it.  Write

\[
 a(r)=trM_r,\qquad T(r)={1\over A}\int_0^r a(q)dq.                \tag{1.4}
\]

If `Q_r` has an absolutely continuous representative in the selected-mass
coordinate `t=T(r)`, the Hilbert-valued Neumann Poincare inequality on
`[0,1]` gives

\[
 \int_0^1\|Q(t)-E Q\|_{HS}^2dt
 \le {1\over\pi^2}\int_0^1\left\|{dQ\over dt}\right\|_{HS}^2dt.
                                                               \tag{1.5}
\]

Consequently the second alternative in (1.3) implies the exact lower bound

\[
 \boxed{
 A\int_{\{a>0\}}{\|Q'_r\|_{HS}^2\over a(r)}dr
 =\int_0^1\left\|{dQ\over dt}\right\|_{HS}^2dt
 \ge {8\pi^2\over17}.}                            \tag{1.6}
\]

If the representative is not absolutely continuous, its corresponding
Dirichlet energy is infinite and (1.6) remains true in the relaxed sense.
Thus the cross-level branch already carries a fixed inverse-amplitude phase
energy.  The difficulty is not extracting that energy; it is charging it to
physical perimeter.  Shape differentiation of `M_r` contains normal
rotation, motion and creation of selected patches, and derivatives of
`omega`.  None is bounded by the scalar coarea deficit at present.  In
particular, the factor `1/a(r)` in (1.6) makes a switch through a selected
boundary-amplitude well maximally visible analytically but potentially
almost free geometrically.  The floor `omega>=10^{-5}` rules out a well
created solely by the selector: `a(r)>=10^{-5}P_mu(A_r)`.

### 1.2 Centralization removes endpoint and trace-amplitude escapes

The floor does more than give qualitative faithfulness.  It permits an
explicit restriction away from extreme-volume levels without losing the
rank.  Write

\[
 d(r)=P_\mu(A_r)-\psi\min\{\mu(A_r),1-\mu(A_r)\},
 \qquad \int_0^1d(r)dr\le\epsilon_*p,               \tag{1.7}
\]

where `epsilon_*=6.02*10^{-5}`, and put `delta=10^{-4}`.  Let

\[
 C=\{r:\delta\le\mu(A_r)\le1-\delta\},\qquad
 M_C=\int_C M_rdr.                                  \tag{1.8}
\]

Since `a(r)=trM_r` satisfies

\[
 10^{-5}P_\mu(A_r)\le a(r)\le P_\mu(A_r),           \tag{1.9}
\]

and `psi<=2p`,

\[
 tr(M-M_C)\le\int_{C^c}P_\mu(A_r)dr
 \le\psi\delta+\mathcal D_{co}(F_0)
 \le(2\delta+\epsilon_*)p=.0002602p.               \tag{1.10}
\]

The constants in the floored transfer therefore give

\[
 \boxed{
 trM_C>.0048495p,\qquad
 {trM_C\over\|M_C\|_{op}}>17.86.}                  \tag{1.11}
\]

For clarity, if `t=b/p>=.00517`, the lower bound for the last ratio is

\[
 { (1-10^{-5})(t-.0000602)-.0002602
       \over .00001+.004t+4(.0000602)}.
\]

This is increasing in `t` and is larger than `17.863` at `t=.00517`.

There is also a pointwise trace floor on `C`.  Since
`psi>=2(1-beta)p` and `beta<=10^{-5}`,

\[
 \boxed{a(r)\ge10^{-5}\psi\delta
       >1.99998\cdot10^{-9}p\quad(r\in C).}         \tag{1.12}
\]

Moreover, for every measurable `E subset C`,

\[
 \int_Ea(r)dr\le p|E|+\epsilon_*p,
 \qquad
 \pi_C(E)\le206.21|E|+.01242,                      \tag{1.13}
\]

where `pi_C` is normalized by `trM_C`.  In particular `|C|>.004789`,
and no phase carrying more than `1.242%` of the central selected flux can be
confined to an arbitrarily short set of level values.  Applying (1.2) to
`M_C` gives the same `8/17` same-level/cross-level dichotomy entirely on
central-volume levels.  Thus neither endpoint concentration nor a vanishing
selected trace remains as an escape.

Continuity still does not convert the cross-level alternative into the
same-level one.  To see the exact logical limitation, take `N>=19`, constant
`a(r)>0`, and a real-analytic curve of unit vectors `u_epsilon(r)` which is
arbitrarily close to `e_j` on `N` successive intervals and makes its
transitions in intervals of width `epsilon`.  Such curves are obtained, for
example, by normalizing soft-max weights of `N` affine functions with
successive upper-envelope intervals.  Then

\[
 Q_r=u_\epsilon(r)u_\epsilon(r)^T,qquad
 1-tr(Q_r^2)=0,qquad
 {tr\int Q_rdr\over\|\int Q_rdr\|_{op}}\longrightarrow N.      \tag{1.14}
\]

This model obeys a uniform trace floor and analytic continuity in `r`, yet
has pure cross-level rank and transitions of arbitrarily small width.  It is
not asserted to come from nested heat levels; it proves that trace,
continuity, and the matrix decomposition alone cannot supply that spatial
realizability step.  The heat derivative identity explains the remaining
escape: although `|nabla^2F_0|op<=C/s`,

\[
             |du/dr|={|P\nabla^2F_0u|\over|\nabla F_0|^2},       \tag{1.15}
\]

and (1.12) is a surface-area bound containing no lower bound for
`|nabla F_0|`.  Hence it gives no upper bound on the transition derivative.
Controlling this lapse, or deriving a finite-splice charge directly from a
rapid transition, is the precise remaining bridge.

## 2. Exact analytic rigidity of planar level patches

The heat function is entire real analytic on the affine support: it is a
Gaussian convolution of a bounded function.  This rules out exact rotating
hyperplane patches without any isoperimetric input.

**Lemma 2.1 (no rotating analytic hyperplane patches).**  Let `F` be a
nonconstant entire real-analytic function on `R^k`, `k>=2`.  Suppose that,
for each `r` in a set `J`, the level `{F=r}` contains a relatively open
subset of an affine hyperplane `H_r`.  Then all `H_r` with distinct level
values are parallel.

**Proof.**  Restrict the analytic function `F-r` to `H_r`.  It vanishes on
a nonempty relatively open set, so the identity theorem makes it vanish on
all of the connected hyperplane `H_r`.  Two nonparallel affine hyperplanes
in `R^k` intersect.  At an intersection point `F` would equal two distinct
values.  Hence the hyperplanes are parallel.

**Corollary 2.2.**  Suppose `Q_r` is rank one for almost every contributing
level and its regular boundary pieces are exactly flat.  Then `Q_r` is
constant in `r` and the second alternative in (1.3) is impossible.

Indeed the underlying selector

\[
 r_G(x)=\int |W(y)|1_G(y)\gamma_s(y-x)dy,
 \qquad \omega(x)={r_G(x)\over R(x)}                 \tag{2.1}
\]

is strictly positive and real analytic, as is its floored version.  Since
the latter is everywhere positive, rank one of `M_r` forces the normal
projectors of the whole regular part of that level to agree almost
everywhere.  Continuity propagates this to every connected regular patch,
and Lemma 2.1 makes the resulting hyperplanes parallel across distinct
levels.  What is not quantitative is the propagation radius of an
*approximately* flat patch or a bound on the derivative of the analytic
selector.

### 2.1 The differential rotation charge

On a noncritical level write `u=nabla F/|nabla F|`, `P=I-uu^T`, and let
`H=nabla^2F`.  The flow

\[
                         {dx\over dr}={u\over|\nabla F|}              \tag{2.2}
\]

advances the level value at unit speed.  Along it,

\[
 \boxed{{du\over dr}={PHu\over|\nabla F|^2}.}        \tag{2.3}
\]

Thus genuine cross-level rotation is paid either by longitudinal covariant
curvature `PHu`, or by passage through a region where `|nabla F|` is small.
The latter has little coarea weight and is exactly the amplitude-well
obstruction already visible in the Fisher derivative decomposition.  Heat
analyticity alone supplies upper derivative bounds; it does not give a
lower bound on the coarea mass of the transition region.

### 2.2 The strongest direct Gaussian rotation charge, and its denominator

For the actual heat function there is a useful pointwise estimate which
tests whether the low-gradient escape can be removed.  Let
`H=nabla^2F_0`, `g=|m|`, `theta=m/g`, `P=I-theta theta^T`, and recall
`d=R-g`.  For every unit vector `h`, Gaussian score differentiation gives

\[
 Hh={1\over\sqrt s}E[W(x+\sqrt sG)(G\cdot h)].                 \tag{2.4}
\]

The posterior centroid bound gives `|W|<=I_0/sqrt(s)`.  Since `Pm=0`,
conditional Cauchy--Schwarz and the exact alignment identity yield

\[
 \begin{aligned}
 |PHh|^2
 &\le {1\over s}E|PW|^2\\
 &\le {I_0\over s^{3/2}}
       E[|W||u-\theta|^2]
 ={2I_0\over s^{3/2}}d(x).
 \end{aligned}                                                \tag{2.5}
\]

Thus

\[
 \boxed{\|PH\|_{op}^2\le {2I_0\over s^{3/2}}d.}              \tag{2.6}
\]

This is dimension-free.  It charges physical-arclength rotation with no
small-gradient denominator.  Indeed, on the noncritical set,

\[
 \nabla_\theta\theta={PH\theta\over g},
\]

and weighted coarea gives, for every `0<=omega<=1`,

\[
 \begin{aligned}
 \mathcal R_{phys}
 &:=\int_0^1\int_{\partial^*A_r}
       \omega|\nabla_\theta\theta|\rho\,d\mathcal H^{k-1}dr\\
 &=\int\omega|PH\theta|d\mu
 \le\left({2I_0\Delta\over s^{3/2}}\right)^{1/2}.             \tag{2.7}
 \end{aligned}
\]

At the fixed scale, `Delta<=epsilon_*p`, `s=alpha K`, and
Buser--Ledoux gives `p sqrt(K)>=c_{BL}/2`.  Normalizing (2.7) by the
central trace in (1.11) therefore gives

\[
 {\mathcal R_{phys}\over trM_C}
 \le {C_{BL}\sqrt{\epsilon_*}
          \over .0048495\,\alpha^{3/4}\sqrt K}.                \tag{2.8}
\]

In particular, within-trajectory normal change per unit physical length is
negligible in the large-`K` regime.

Equation (2.7) does **not** control the rotation relevant to conditioning on
the level value.  Along the level flow,

\[
 {d\theta\over dr}={PH\theta\over g^2},
\]

and the corresponding `L^1` quantity is

\[
 \boxed{
 \mathcal R_{lev}
 =\int_0^1\int_{\partial^*A_r}
       \omega\left|{d\theta\over dr}\right|\rho
 =\int\omega{|PH\theta|\over g}d\mu.}             \tag{2.9}
\]

The remaining factor `1/g` cannot be bounded by (1.12): that estimate is a
surface-area trace floor and contains no information about the pointwise
gradient.  Nor does the scalar coarea deficit contain `1/g`.  A normal may
therefore turn by order one over a long physical path on which `F_0` changes
through an arbitrarily short interval of level values.  Formula (2.7) is the
explicit direct charge available from the actual heat form; (2.9) is the
sharp place where the proposed argument stops.  To finish, one must either
couple phases in physical arclength/volume rather than in `r`, charge the
resulting focal low-gradient region by a finite splice, or prove a new lower
gradient (lapse) estimate.

### 2.3 Exact continuity equation for locked phase trajectories

The remaining reweighting can be isolated exactly.  Work first on an open
interval of regular values on which `g>0`, the density is smooth, and the
level flow does not contact the boundary of the support.  Let

\[
 \sigma_r=\omega\rho\,\mathcal H^{k-1}|_{\Sigma_r},\qquad
 Z=\theta\theta^T,\qquad v={\theta\over g},
 \qquad H_\Sigma=\operatorname {div}_{\Sigma_r}\theta.          \tag{2.10}
\]

For every smooth scalar or matrix-valued test field `Phi`, the hypersurface
transport formula is

\[
 {d\over dr}\int_{\Sigma_r}\Phi\,d\sigma_r
 =\int_{\Sigma_r}\{v\cdot\nabla\Phi+c\Phi\}\,d\sigma_r,        \tag{2.11}
\]

where

\[
 c={H_\Sigma-\theta\cdot\nabla V\over g}
       +{\theta\cdot\nabla\log\omega\over g}
   =:c_{geom}+c_{sel}.                                         \tag{2.12}
\]

The first term is weighted mean-curvature expansion of the trajectory
bundle; the second is selector reweighting.  Since

\[
 \dot Z:=v\cdot\nabla Z
 ={PH\theta\over g^2}\theta^T
   +\theta{(PH\theta)^T\over g^2},
 \qquad \|\dot Z\|_{HS}={\sqrt2|PH\theta|\over g^2},            \tag{2.13}
\]

(2.11) gives the exact normalized equation

\[
 \boxed{
 Q'_r={1\over a(r)}\int\dot Z\,d\sigma_r
      +{1\over a(r)}\int c\{Z-Q_r\}\,d\sigma_r.}              \tag{2.14}
\]

The first term is within-trajectory rotation.  The second is precisely the
birth/death or reweighting of orientation-locked trajectories.  In
particular, if `dot Z=0` and the logarithmic expansion rate `c` is
uncorrelated with `Z` on every level, then `Q_r` is constant.

The total action of the first term in the `r`-continuity equation is

\[
 \int a(r)\left\|{1\over a(r)}\int\dot Z\,d\sigma_r\right\|dr
 \le\sqrt2\int\omega{|PH\theta|\over g}d\mu
 =\sqrt2\mathcal R_{lev}.                                    \tag{2.15}
\]

Thus it is (2.9), not the smaller physical-arclength charge (2.7), which
controls variation of `Q_r`.  The requested bound by `R_phys` is false at
the kinematic level: synchronization by the scalar value `r` divides the
physical derivative by `g`.

The two smooth source terms have equally explicit total-variation bounds:

\[
 \begin{aligned}
 \int a\|S_{sel}(r)\|_{HS}dr
   &\le\sqrt2\int|\nabla\omega|d\mu,\\
 \int a\|S_{geom}(r)\|_{HS}dr
   &\le\sqrt2\int\omega|H_\Sigma-\theta\cdot\nabla V|d\mu.
                                                               \tag{2.16}
 \end{aligned}
\]

Here `S_sel,S_geom` denote the two covariance terms in (2.14).  The
Gaussian selector-information estimate does control the more strongly
flux-weighted quantity.  If `omega_0=r_G/R`, `B_0=E R<=p`, and
`M_0=I_0/sqrt(s)`, then

\[
 \int R|\nabla\omega_0|^2d\mu
 \le{B_0\over s}\{2+4\log(M_0/B_0)\},                         \tag{2.17}
\]

and, because `g<=R`,

\[
 \int g|\nabla\bar\omega|d\mu
 \le {B_0\over\sqrt s}
       \{2+4\log(M_0/B_0)\}^{1/2}.                            \tag{2.18}
\]

But (2.16) requires the unweighted `int|nabla omega|`, and neither the floor
nor (2.17) controls it where `R` is small.  Similarly, profile second
variation controls curvature of actual minimizers, not the weighted mean
curvature of these merely integrated near-minimizing levels.  Therefore
neither smooth source in (2.16) is presently charged.

At critical values or contact with a hard support, (2.14) acquires a
distributional residual `dB` recording creation, annihilation, or merging of
trajectory sheets.  Under a local `BV` hypothesis on `r mapsto(a(r),M_r)`,
`dB` is a matrix-valued Radon measure supported on the corresponding
critical/contact values.  A nonflat part of this residual with positive
ridge capacity is charged by the finite bevel lemma; a separated residual
with no ridge capacity is exactly the profile-linearity/log-affine branch.
No quantitative lower bound for either capacity is yet derived from
`|dB|`.  Equations (2.14)--(2.18) therefore give the requested exact
superposition decomposition and leave only three explicit obstructions:

1. the lapse factor in the within-trajectory action (2.15);
2. smooth orientation-dependent expansion/selector reweighting (2.16); and
3. the focal/contact residual `dB`.

## 3. A finite lattice-straightening operation

For finite-perimeter sets `E,F`, perimeter submodularity gives

\[
 \boxed{P(E\cap F)+P(E\cup F)\le P(E)+P(F).}         \tag{3.1}
\]

Replacing the unordered pair `(E,F)` by the nested pair
`(E union F,E intersection F)` is therefore a finite, physical,
null-invariant straightening.  It preserves the pointwise membership count
and hence

\[
 \mu(E\cup F)+\mu(E\cap F)=\mu(E)+\mu(F).           \tag{3.2}
\]

Iterating (3.1) sorts any finite family `E_1,...,E_m` into the nested order
statistics

\[
 B_j=\left\{x:\sum_{i=1}^m1_{E_i}(x)\ge j\right\},
 \quad B_{j+1}\subset B_j,                           \tag{3.3}
\]

and proves

\[
 \boxed{
 \sum_jP(B_j)\le\sum_iP(E_i),\qquad
 \sum_j\mu(B_j)=\sum_i\mu(E_i).}                   \tag{3.4}
\]

For two proxy halfspaces of masses `v_1>=v_2`, put
`c=mu(E_2 setminus E_1)`.  The sorted masses are exactly

\[
                         v_1+c,qquad v_2-c.         \tag{3.5}
\]

Concavity of the isoperimetric profile gives the exact scalar escape

\[
 \mathcal I(v_1)+\mathcal I(v_2)
 -\mathcal I(v_1+c)-\mathcal I(v_2-c)\ge0.          \tag{3.6}
\]

When the profile is strictly concave, (3.6) records the cost of crossing
the intended level order.  When `mathcal I(v)=psi v` on the relevant small
side, it vanishes identically.  Thus union/intersection cannot by itself
create an angular charge in the exponential equality model; its exact
exception is again profile linearity.

This operation becomes a candidate physical competitor as follows.  For a
finite selection of levels, replace each nearly planar level set by its
same-mass halfspace proxy, sort the proxies by (3.3), and threshold the
resulting nested family.  Equation (3.4) controls perimeter exactly.  The
unproved step is a dimension-free estimate converting concentration of the
selected normal projector into small symmetric difference and perimeter
error between the original level and the proxy halfspace.

## 4. Exact rotating-pencil stress test

Let `mu` be the product of two standard one-sided exponentials on the
positive quadrant and put

\[
                         F(x,y)={x\over x+y}.         \tag{4.1}
\]

The function is analytic in the interior.  Its superlevels are nested
halfspace sectors whose boundary lines all pass through the support vertex.
On the level `F=r`, write

\[
 x=rT,\qquad y=(1-r)T,qquad T>0.                   \tag{4.2}
\]

The normal is constant on the level and proportional to `(1-r,-r)`, so
every `Q_r` has rank one while its projective line rotates by a right angle
between the endpoint levels.

The ratio `X/(X+Y)` is uniform on `(0,1)`, hence

\[
                         \mu\{F>r\}=1-r.             \tag{4.3}
\]

The relative Euclidean perimeter is exactly

\[
 \boxed{
 P_\mu\{F>r\}
 =\int_0^\infty e^{-T}
       \sqrt{r^2+(1-r)^2}\,dT
 =\sqrt{r^2+(1-r)^2}.}                              \tag{4.4}
\]

A coordinate halfspace gives `psi_mu<=1`.  Therefore the coarea deficit of
this pencil obeys

\[
 \begin{aligned}
 \mathcal D_{co}(F)
 &\ge\int_0^1\left\{
    \sqrt{r^2+(1-r)^2}-\min(r,1-r)\right\}dr\\
 &\ge {1\over\sqrt2}-{1\over2}>0.207.              \tag{4.5}
 \end{aligned}
\]

Thus the focal rotation pays a fixed scalar charge.  Straight coordinate
cuts are the corresponding finite competitors.  After centering, the
product has covariance `I_2` and bounded Poincare constant, so it lies in the
allowed bounded-product branch.

The same geometry appears for barycentric ratios in a simplex: rotating
hyperplanes focus on a lower-dimensional support face.  A cube can realize
piecewise pencils at corners.  Both have bounded `K`, and their nonaffine
inner cuts have a fixed perimeter gap or bevel capacity.  A log-affine slab
has parallel levels and hence zero cross-level term.  For a radial Gaussian
function, the full level normal matrix already has high same-level rank; an
artificial weight selecting one point on each sphere could manufacture a
cross-level term, which is why the analytic structure of the particular
weight, rather than an arbitrary `0<=omega<=1`, matters.

## 5. Audit of the profile second-variation route

There is a genuine quantitative statement behind the suggestion that a
nearly linear profile has little curvature.  It is important, however, to
separate the curvature which it controls from the rotation of the selected
normal field.

### 5.1 What profile near-linearity actually controls

Assume for this subsection that `rho=e^{-V}` is positive and `C^2`, that
`nabla^2V>=0`, and that an isoperimetric minimizer of volume `v` has a
compact `C^2` regular boundary `Sigma_v` with constant weighted mean
curvature.  Put

\[
 Kappa(v)=\int_{\Sigma_v}
       \{\|II\|_{HS}^2+\nabla^2V(n,n)\}\rho\,d\mathcal H^{k-1}.
                                                               \tag{5.1}
\]

The outward parallel deformation gives the support inequality

\[
 \boxed{\quad \mathcal I''(v)\le
                  -{Kappa(v)\over\mathcal I(v)^2}.\quad}       \tag{5.2}
\]

Here (5.2) is in the upper-support, hence viscosity, sense.  To check the
normalization, let `q(t)` and `p(t)` be respectively the volume and weighted
perimeter of the outward parallel set.  At zero,

\[
 q'=p=\mathcal I(v),\qquad
 p'=h\mathcal I(v),\qquad
 p''=h^2\mathcal I(v)-Kappa(v),                    \tag{5.3}
\]

where `h` is the constant weighted mean curvature.  Therefore, after using
volume as parameter,

\[
 {d^2p\over dq^2}(v)=-{Kappa(v)\over\mathcal I(v)^2}.
\]

The parallel sets are competitors and are tangent to the profile at `v`,
which proves (5.2).  The same statement on the regular part, followed by the
usual cutoff around the codimension-at-least-seven singular set, is the
standard profile support inequality.  No such limiting statement is used
below without these regularity hypotheses.

The following elementary constant tracking records the useful consequence.

**Lemma 5.1 (curvature budget on a multiplicative interval).**  Let
`I:[a,b]to R` be concave and suppose, for `c>0` and `0<epsilon<1`,

\[
                       ct\le I(t)\le(1+\epsilon)ct
                       \quad(a\le t\le b).          \tag{5.4}
\]

Let `nu=-I''` be its distributional curvature measure.  If
`0<eta<1/2`, `L=b-a`, `x=a+eta L`, and `y=b-eta L`, then

\[
                 \boxed{\nu((x,y))\le
                 {2\epsilon cb\over\eta L}.}       \tag{5.5}
\]

**Proof.**  The function `e(t)=I(t)-ct` is concave and
`0<=e<=epsilon cb`.  Concavity bounds its right derivative at `x` above by
`epsilon cb/(eta L)` and its left derivative at `y` below by the negative
of the same number.  Their drop is (5.5).

If (5.2) is available measurably throughout the interval, it follows that

\[
 \int_x^y {Kappa(t)\over I(t)^2}\,dt
       \le {2\epsilon cb\over\eta L}.               \tag{5.6}
\]

In particular, on `[a,b]=[v,2v]` with `eta=1/4`, some
`t in [5v/4,7v/4]` satisfies

\[
 {Kappa(t)\over I(t)}
       \le56(1+\epsilon)\epsilon c^2.               \tag{5.7}
\]

Thus the proposed `O(epsilon psi^2)` normalized curvature conclusion is
correct, with explicit constants, whenever (5.4) holds with `c=psi` and
the stated minimizer regularity is justified.

### 5.2 Why this does not yet control cross-level phase rotation

Let an arbitrary smooth nested foliation be written in volume parameter
`v`.  On its level `Sigma_v`, its normal speed has the form `phi_v n` and

\[
                       \int_{\Sigma_v}\phi_v\rho=1.
\]

The kinematic identity for its normal is

\[
                         \partial_v n=-\nabla_{\Sigma_v}\phi_v.       \tag{5.8}
\]

For the actual levels of `F`, if

\[
 J(r)=\int_{\{F=r\}}{\rho\over|\nabla F|}\,d\mathcal H^{k-1},
\]

then, up to orientation,

\[
 \phi_v={1\over J(r)|\nabla F|},\qquad
 |\nabla_{\Sigma_v}\phi_v|
 ={ |P\nabla^2F\,n|\over J(r)|\nabla F|^2}.          \tag{5.9}
\]

The energy in (5.9), not `Kappa`, measures the rotation in (2.3).  It is the
tangential gradient of the *lapse* of the foliation.  It does not occur in
the constant-speed parallel comparison which proves (5.2).

Even for a differentiable nested family of exact minimizers, the second
variation in its own speed is

\[
 \mathcal I''(v)=\int_{\Sigma_v}
       \{ |\nabla_{\Sigma_v}\phi_v|^2
          -(\|II\|^2+\nabla^2V(n,n))\phi_v^2\}\rho. 
                                                               \tag{5.10}
\]

Thus (5.6), which controls the curvature with *constant* weight, does not
control either term in (5.10) when `phi_v` is unbounded or concentrated.
A Harnack bound `phi_v asymp 1/I(v)` would bridge this gap, but no such
dimension-free estimate is available and hard-support focal examples make
it false without additional hypotheses.

There are two further mismatches in the present application.

1. The actual sets `A_r` are only near-minimizers in an integrated sense.
   A small `int(P(A_r)-I(mu(A_r)))dr` controls no second derivative of their
   perimeter as a function of volume.  Narrow nonnegative perimeter bumps
   can have arbitrarily large second derivative.
2. The matrices `Q_r` use the selector `omega`, whereas (5.2) and (5.10)
   concern the whole boundary.  The floored selector now gives the static
   comparison `10^{-5}<=omega<=1`, so selected and whole boundary mass are
   comparable by a universal factor.  But there is still no derivative
   bound or transport identity coupling the selected mass on different
   levels.  The selector can redistribute its fixed-factor weights rapidly
   between phase packets even if no packet rotates along a normal
   trajectory.

Consequently Cauchy--Schwarz over the natural distance
`int dv/I(v)=O(1/psi)` becomes valid only after proving a lapse bound and a
dynamic transport bound for the weighted boundary packets.  Static selector
faithfulness is now available, but it does not supply either dynamic bound.
These statements, rather than the profile curvature estimate (5.6), are the
precise missing bridge.  Support contact adds a fourth issue: parallel
variation ceases at the focal boundary, as in Section 4.

## 6. Exact status of the cross-level branch

The following implication is proved:

* exact flat selected patches plus analyticity imply parallel phases;
* the floored analytic selector gives a central-level trace floor and leaves
  effective rank larger than `17.86`, so endpoint and selector-amplitude
  concentration are excluded;
* a numerical same-level ridge expander gives more than `10^{-4}p` finite
  bevel saving;
* exactly separated zero-deficit components force a linear profile and, in
  the smooth full-support case, an affine bounded-`K` branch;
* the canonical rotating hard-support pencil pays a fixed coarea charge.

What remains unproved is a quantitative theorem at the audited budget

\[
             \mathcal D_{co}(F_0)\le6.02\cdot10^{-5}p.             \tag{6.1}
\]

It must show that cross-level projector variance `8/17` forces at least one
of the following with fixed constants:

1. enough `PHu` or ridge capacity for a finite bevel/straightening saving;
2. enough low-gradient transition mass to exceed (6.1);
3. a quantitative affine exponential splitting, hence bounded `K`; or
4. a genuinely concurrent radial organization.

The continuity equation (2.14) leaves exactly three residual terms.

1. **Within-trajectory lapse.**  The actual heat form satisfies the sharp
   dimension-free bound (2.7) on physical-arclength rotation, but the
   `r`-synchronized action is (2.15) and contains one additional factor
   `1/|nabla F_0|`.  The central trace floor and the scalar coarea gap do not
   control this factor.
2. **Smooth phase reweighting.**  Geometry and selector contributions are
   exactly bounded by the two quantities in (2.16).  Selector Fisher
   information bounds only the more strongly flux-weighted expression
   (2.18), while profile second variation applies to true minimizers rather
   than the weighted mean curvature of the actual levels.  Neither available
   estimate bounds (2.16).
3. **Focal/contact birth and death.**  The distributional residual `dB`
   records critical mergers and hard-support contact.  Positive ridge
   capacity is finitely bevelled and exact zero ridge capacity leads to the
   profile-linear/log-affine branch, but there is no quantitative lower
   bound transferring `|dB|` to either charge.

Analyticity proves the zero-error planar statement but supplies no
dimension-free quantitative propagation through these residual terms.  The
rotating exponential pencil shows the focal mechanism and pays a fixed
charge; it does not provide the missing general lower bound.  Charging the
three terms above at cost below (6.1) is the remaining cross-level
phase-charge theorem.
