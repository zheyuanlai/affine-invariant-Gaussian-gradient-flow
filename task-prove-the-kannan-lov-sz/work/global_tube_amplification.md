# Finite-distance amplification for hard support

## 0. Outcome

This note gives a finite-distance version of the free-boundary profile
curvature estimate in `hard_support_profile_variation.md`.  It answers the
scaling question affirmatively, including first support contact, focal time,
and cut-locus collisions.

Let `Omega` be a smooth convex support, let

\[
                 d\mu=Z^{-1}e^{-V}1_\Omega dx,
\]

with `V` convex, and let `E` be a regular isoperimetric region of volume
`v_0`.  Write `Sigma=partial E cap Omega`, `P_0=I(v_0)`, and let `lambda` be
the constant weighted mean curvature of `Sigma`.  For `x in Sigma`, follow
the exterior normal ray `x+tN_x` only until its first support contact, focal
time, or loss of unique nearest point.  Denote this stopping time by
`tau(x)`.  No positive lower bound on `tau`, and hence no uniform reach, is
assumed.

Before it is killed, the weighted normal Jacobian is exactly

\[
 j_x(t)=\exp\{\lambda t-D_x(t)\},                         \tag{0.1}
\]

where

\[
 D_x(t)=\int_0^t(t-s)\left[
  \operatorname {tr}\big(S_x^2(I+sS_x)^{-2}\big)
  +\nabla^2V(x+sN_x)[N_x,N_x]\right]ds\ge0.               \tag{0.2}
\]

Thus the killed, exponentially normalized flux

\[
 R(t)=\int_\Sigma 1_{\{t<\tau(x)\}}e^{-D_x(t)}d\sigma_\mu(x) \tag{0.3}
\]

is nonincreasing.  Every focal/contact/collision event is charged by killing
the whole remaining base-ray weight; there is no Taylor remainder.

Here is the quantitative conclusion.  Suppose `0<a<v_0<b<=1/2`, put

\[
 c={I(b)\over b},\qquad \beta={b\over v_0},
\]

and assume

\[
 {I(a)\over a}\le(1+\delta)c,
 \qquad s:=\min(v_0-a,b-v_0)>0,
 \qquad \eta:={\delta b\over s}<1.                        \tag{0.4}
\]

Let `T_b` be the first time the exterior parallel set of `E` has volume
`b`.  Then

\[
 {1\over(1+\eta)c}
 \log\left(1+{1+\eta\over1+\delta}(\beta-1)\right)
 \le T_b\le {\log\beta\over c},                           \tag{0.5}
\]

and, with the left-limit convention at killed rays,

\[
 \boxed{\quad
 P_0-R(T_b-)\le
 \left(1-{\beta^{-\eta}\over1+\delta}\right)P_0
 \le(\delta+\eta\log\beta)P_0.\quad}                    \tag{0.6}
\]

For fixed `a/v_0` and `b/v_0`, (0.5) says `T_b asymp 1/c`.  If the relevant
volumes are fixed universal numbers and `P_0 asymp p`, then `c asymp p`, so
`T_b asymp1/p`; (0.6) is precisely the desired `O(delta p)` finite-distance
charge.  In the smooth no-cut regime, (0.2) shows directly why a transported
curvature budget `O(delta p^3)` accumulates twice over time `T asymp1/p` to
`O(delta p)`.  Formula (0.6) remains valid when this smooth picture breaks:
all missing Jacobian is charged as contact/focal/cut loss.

The lemma has a useful rigidity consequence.  Apart from a surface set of
measure at most

\[
 \left(1+{1\over1-e^{-h}}\right)
       (\delta+\eta\log\beta)P_0,                         \tag{0.7}
\]

every ray survives to `T_b` and has `D_x(T_b-)<=h`.  Long good rays give the
covariance matrix inequality

\[
 \operatorname {Cov}(\mu)\succeq
 {e^{-|\lambda|T_b-h}T_b^3\over12}
 \int_G N_x\otimes N_x\,d\sigma_\mu(x).                  \tag{0.8}
\]

Consequently an isotropic measure cannot carry surface mass comparable to
`p` in one coherent normal direction when `p` is small.  But (0.8) does not
exclude a high-rank normal distribution: it only gives

\[
 \left\|\int_GN\otimes N\,d\sigma_\mu\right\|_{op}
 \lesssim p^3.                                           \tag{0.9}
\]

An exact product-exponential model in Section 6 shows that this is a real
logical obstruction.  Near-linear flux can be split among arbitrarily many
transverse planar sheets while their pairwise collision ridge has only
quadratic mass.  A successful final argument must classify that branch as
an almost-product/multi-coordinate exponential geometry; profile curvature
and tube loss alone cannot turn it into a one-dimensional affine marginal.

There is also a separate **transfer gap**.  The theorem concerns an exact
regular isoperimetric leaf, hence a constant-weighted-mean-curvature
surface.  The ranked physical normal matrix in the surrounding argument is
extracted from levels of a near-minimizing smoothed function.  Near-minimal
perimeter by itself does not make those levels stationary and does not
preserve their normal matrix under replacement by an exact minimizer.  The
theorem must not be applied to that physical matrix without the quantitative
stability statement isolated in Section 8.

## 1. The killed normal tube formula

Assume for Sections 1--5 that `Omega subset R^n` is a bounded convex body
with `C^2` boundary, `V in C^2(overline Omega)` is convex, and `E` is taken
as a closed regular representative whose **entire** relative boundary
`Sigma=closure(partial E cap Omega)` is a compact `C^2`
hypersurface-with-boundary meeting `partial Omega` orthogonally.  Thus no
singular internal boundary has been omitted from `Sigma`.  Boundedness
guarantees that the target parallel volume is reached.  The local tube
formula also holds on an unbounded support by the sigma-finite area formula,
but extending the full theorem there additionally requires existence and
regularity of the stated minimizer and the stated profile hypotheses; no
unproved exhaustion transfer is being invoked.  In dimensions in which an
isoperimetric minimizer may have a singular set, one must either prove that
its metric projection contributes zero tube volume or include that
contribution separately; this note does not silently discard it.

Let `N` be the outer normal of `E`, and use

\[
 S_xu=D_uN,\qquad
 \lambda=\operatorname {tr}S_x-\langle\nabla V(x),N_x\rangle. \tag{1.1}
\]

For `x` in the relative interior of `Sigma`, define `tau(x)` as the supremum
of `t>0` such that, for every `0<s<t`,

1. `x+sN_x` belongs to `Omega setminus E`;
2. `x` is the unique nearest point of `E` to `x+sN_x`; and
3. `I+sS_x` is nonsingular.

The third item is redundant before a genuine focal point but makes the
Jacobian convention explicit.  The function `tau` is measurable.  The map

\[
 F(x,t)=x+tN_x                                             \tag{1.2}
\]

is one-to-one on `{(x,t):0<t<tau(x)}`.  Conversely, outside a Lebesgue-null
set, every point of `Omega setminus E` with positive distance from `E` and
whose nearest point lies on the regular part of `Sigma` has exactly this
form.  Uniqueness of the metric projection holds almost everywhere because
the distance to a closed set is Lipschitz and is differentiable only where
the nearest point is unique.

Here is the required free-boundary normal-cone argument.  If `y in Omega`
has a nearest point `x in partial Sigma`, put `w=y-x`.  The nearest-point
inequality puts `w` in the polar of the tangent cone of `E` at `x`.  Since
the two `C^2` boundary constraints meet transversely and orthogonally, that
polar is

\[
                 \{aN_x+b\nu_x:a,b\ge0\},                \tag{1.2a}
\]

where `nu` is the outer normal of `Omega`.  Convexity of `Omega` gives the
exact supporting-halfspace inequality

\[
                         \langle\nu_x,y-x\rangle\le0.     \tag{1.2b}
\]

But orthogonality and (1.2a) make the left side equal to `b`; hence `b=0`
and `y=x+aN_x`.  The image of
`partial Sigma times [0,k]` under `(x,a) mapsto x+aN_x` has dimension at
most `n-1` and therefore zero `n`-dimensional measure.  Taking the countable
union over `k` proves that projections onto `partial Sigma` contribute no
tube volume.  Applying the area formula on the injectivity domain now gives
the formulas below.  No uniform lower bound on local tubular radii is used.

**Lemma 1.1 (exact killed Jacobian).**  For `0<t<tau(x)`, the weighted
Jacobian of (1.2), relative to `d sigma_mu(x)dt`, is

\[
 j_x(t)=\exp[-V(x+tN_x)+V(x)]\det(I+tS_x).                 \tag{1.3}
\]

It satisfies (0.1)--(0.2).  In particular,

\[
 0<j_x(t)\le e^{\lambda t},\qquad
 t\longmapsto e^{-\lambda t}j_x(t)=e^{-D_x(t)}
 \quad\hbox{is nonincreasing}.                            \tag{1.4}
\]

**Proof.**  Differentiating (1.2) in a principal frame gives the Euclidean
Jacobian `det(I+tS_x)`.  Along a surviving ray this determinant is positive.
Moreover,

\[
 {d\over dt}\log j_x(t)
 =\operatorname {tr}\big(S_x(I+tS_x)^{-1}\big)
  -\partial_NV(x+tN_x),                                   \tag{1.5}
\]

and

\[
 {d^2\over dt^2}\log j_x(t)
 =-\operatorname {tr}\big(S_x^2(I+tS_x)^{-2}\big)
  -\nabla^2V(x+tN_x)[N_x,N_x]\le0.                       \tag{1.6}
\]

At zero, (1.5) equals `lambda`.  Twice integrating (1.6) gives
(0.1)--(0.2), and (1.4) follows.  QED.

Let

\[
 E_t=\{y\in\Omega:\operatorname {dist}(y,E)\le t\},
 \qquad v(t)=\mu(E_t).                                    \tag{1.7}
\]

The coarea and area formulas give, for almost every `t>0`,

\[
 v'(t)=P_\mu(E_t;\Omega)
      =\int_\Sigma1_{\{t<\tau(x)\}}j_x(t)d\sigma_\mu(x)
      =e^{\lambda t}R(t).                                 \tag{1.8}
\]

This statement includes first support contact, focal collapse, and collision
of different normal rays.  At such a time the affected ray is absent from
(1.8); hence its full remaining normalized flux is lost from `R`.

For bookkeeping at atoms of the cut-time distribution, define

\[
 R(t-)=\int_\Sigma1_{\{t\le\tau(x)\}}e^{-D_x(t)}d\sigma_\mu(x). \tag{1.9}
\]

When `tau(x)=t`, the notation in (1.9) means
`D_x(t)=lim_{s upward t}D_x(s) in [0,infinity]`; at a focal time the value
may be infinite.  Then `R(t-)` is the left limit of (0.3), and `R(0)=P_0`.

## 2. Profile amplification over a multiplicative interval

For Theorem 2.1, concavity of `I` is an explicit hypothesis.  In the smooth
setting above it is supplied by the weighted isoperimetric-profile
concavity theorem: the hypotheses are Bakry--Emery Ricci tensor
`Hess V>=0`, convex boundary `II_(partial Omega)>=0`, and relative
(Neumann) perimeter.  The same conclusion follows in viscosity form from
the constant-speed free-boundary second variation.  In particular,
`I(0)=0` and `v mapsto I(v)/v` is nonincreasing on `(0,1)`.  No Poincare or
KLS-strength input is used here.

**Theorem 2.1 (global tube amplification).**  Suppose `E` is isoperimetric,
`v_0=mu(E)`, `P_0=I(v_0)`, and `Sigma` has constant weighted mean curvature
`lambda`.  Under (0.4), formulas (0.5)--(0.6) hold.

**Proof.**  Put `F(v)=I(v)-cv`.  On `[a,b]`, monotonicity of `I(v)/v`
and (0.4) give

\[
                    0\le F(v)\le\delta cb.               \tag{2.1}
\]

At a differentiability point `v_0` of the profile, stationarity identifies
`lambda=I'(v_0)`.  Concavity and the two secants from `a` and to `b` imply

\[
 -{\delta cb\over b-v_0}
 \le\lambda-c\le {\delta cb\over v_0-a},
 \qquad |\lambda-c|\le\eta c.                             \tag{2.2}
\]

The same inequalities hold for every supporting slope at `v_0`, so the
argument does not require classical differentiability of `I`.

Until `v(t)` reaches `b`, the parallel set is an admissible competitor and

\[
 e^{\lambda t}R(t)=P_\mu(E_t;\Omega)
       \ge I(v(t))\ge cv(t)                                \tag{2.3}
\]

for almost every `t`.  Thus `v'(t)>=cv(t)` a.e., and

\[
                         T_b\le{\log\beta\over c}.        \tag{2.4}
\]

On the other hand `R(t)<=R(0)=P_0`, so

\[
 b-v_0\le P_0\int_0^{T_b}e^{\lambda t}dt.                 \tag{2.5}
\]

Because `eta<1`, (2.2) gives `lambda>0`.  Solving (2.5), using that
`lambda mapsto lambda^{-1}log(1+lambda A)` is decreasing, together with
`lambda<=(1+eta)c` and `P_0/v_0=I(v_0)/v_0<=(1+delta)c`, gives the lower
bound in (0.5).

Choose regular times `t_k upward T_b`.  From (2.3), (2.4), and
`v(t_k) to b`,

\[
 \begin{aligned}
 {R(T_b-)\over P_0}
 &\ge {cb\over P_0}e^{-\lambda T_b}\\
 &= {c\over P_0/v_0}\,\beta e^{-\lambda T_b}\\
 &\ge {1\over1+\delta}\,\beta^{1-\lambda/c}
 \ge {1\over1+\delta}\,\beta^{-\eta}.                  \tag{2.6}
 \end{aligned}
\]

Finally,

\[
 1-{\beta^{-\eta}\over1+\delta}
 ={\delta\over1+\delta}
   +{1-\beta^{-\eta}\over1+\delta}
 \le\delta+\eta\log\beta,                               \tag{2.7}
\]

which is (0.6).  QED.

### 2.1 Exact accounting for curvature and singular times

For every `T>0`, (0.3) gives the identity

\[
 \begin{aligned}
 P_0-R(T-)
 ={}&\sigma_\mu\{x:\tau(x)<T\}\\
 &+\int_{\{\tau(x)\ge T\}}(1-e^{-D_x(T-)})d\sigma_\mu(x). \tag{2.8}
 \end{aligned}
\]

Thus, if the right side of (0.6) is denoted by `epsilon P_0`,

\[
 \sigma_\mu\{\tau<T_b\}\le\epsilon P_0,                \tag{2.9}
\]

and, for every `h>0`,

\[
 \sigma_\mu\{\tau\ge T_b,\ D_x(T_b-)>h\}
 \le {\epsilon P_0\over1-e^{-h}}.                        \tag{2.10}
\]

This proves (0.7).  Also

\[
 \int_{\{\tau\ge T_b\}}\min(D_x(T_b-),1)d\sigma_\mu(x)
 \le {\epsilon P_0\over1-e^{-1}}.                        \tag{2.11}
\]

In particular, (2.9) charges every first support contact, focal time, or
collision strictly before `T_b`, with no regularity-dependent constant.  A
ray whose first event is exactly `T_b` still supplies the full open tube
`0<t<T_b` needed in Lemma 3.1; a focal endpoint is additionally recorded by
`D_x(T_b-)=infinity` in (2.10).

If there is no killing and `D_x(T)<=1`, then

\[
 P_0-R(T)\asymp
 \int_\Sigma\int_0^T(T-s)K_x(s)\,ds\,d\sigma_\mu(x),     \tag{2.12}
\]

where

\[
 K_x(s)=\operatorname {tr}(S_x^2(I+sS_x)^{-2})
       +\nabla^2V(x+sN_x)[N_x,N_x].                       \tag{2.13}
\]

Therefore a transported integral of `K` of size `delta p^3` over each
unit of the relevant surface family, sustained for `T asymp1/p`, contributes
`delta p^3 T^2=O(delta p)`.  An estimate only at `s=0` does **not** imply
(2.12) at finite `T`; Theorem 2.1 replaces that unavailable propagation
estimate by the global profile comparison.

## 3. Covariance forced by surviving rays

**Lemma 3.1 (long-ray covariance).**  Let `G subset Sigma` be measurable and
suppose, for all `x in G`,

\[
                    \tau(x)\ge T,\qquad D_x(T-)\le h.    \tag{3.1}
\]

Then (0.8) holds.

**Proof.**  Since `D_x(t)` is nondecreasing, on `0<=t<=T`

\[
                  j_x(t)\ge e^{-|\lambda|T-h}.            \tag{3.2}
\]

The stopped tube map is injective.  For `xi in R^n`, restrict the variance
integral to its image.  On the ray based at `x`, the function
`<xi,x+tN_x>` is affine in `t`.  Minimizing its squared integral over the
constant term and using (3.2) gives

\[
 \int_0^T(\langle\xi,x+tN_x\rangle-a)^2j_x(t)dt
 \ge e^{-|\lambda|T-h}{T^3\over12}
          \langle\xi,N_x\rangle^2                         \tag{3.3}
\]

for every `a in R`.  Integrate over `G` and take `a` to be the global mean.
This proves the quadratic-form inequality (0.8).  QED.

For an isotropic measure, Lemma 3.1 implies

\[
 \int_GN\otimes N\,d\sigma_\mu
 \preceq 12e^{|\lambda|T+h}T^{-3}I.                       \tag{3.4}
\]

If `T=T_b`, then (0.5), `|lambda-c|<=eta c`, and fixed `beta` make the
right side `C(delta,eta,beta)c^3I`.  Hence any packet `G` satisfying

\[
 \left\|\int_GN\otimes N\,d\sigma_\mu\right\|_{op}
 \ge\kappa\,\sigma_\mu(G),\qquad
 \sigma_\mu(G)\ge\alpha P_0                             \tag{3.5}
\]

forces, when `P_0 asymp c`, a universal lower bound
`c^2 >= c_0(alpha,kappa,delta,eta,beta)`.  This closes the coherent/parallel
normal branch.

The trace of (3.4), however, only yields

\[
              \sigma_\mu(G)\le Cc^3 n.                   \tag{3.6}
\]

Thus a packet of mass `alpha p` is compatible with isotropy whenever its
normal matrix has effective rank of order at least `1/p^2`.  The global tube
lemma does not turn a high-rank collection of sheets into one direction.

## 4. What exact intersection does and does not add

For two exact planar free-boundary sections of one convex body, projection
onto the plane spanned by their normals turns both sections into affine
diameters of the same planar convex projection.  The affine-diameter
intersection theorem therefore forces the sections to meet in the interior.
This classifies the zero-defect case: disjoint components are parallel.

At finite defect, an intersection may be hidden in a small ridge or neck.
Normal tubes based near that ridge collide early, and (2.9) charges them.
But neither convexity nor the `L^2` curvature budget supplies a lower bound
for the surface measure of the ridge neighborhood.  In dimensions at least
three, a fixed normal rotation can occur through a neck whose
`int |II|^2` tends to zero.  The exact intersection theorem therefore
upgrades (0.6) only if one already has a quantitative ridge-capacity or
conductance lower bound.  Such a bound is the missing input, not a remainder
in the tube expansion.

## 5. Model checks

1. **Log-affine cylinder.**  On a planar section normal to the affine
   direction, `S=0`, `V_NN=0`, and no ray is killed before the end cap.
   Thus `D=0` and `R=P_0` exactly.  Lemma 3.1 is the usual one-dimensional
   variance obstruction and has the correct scale.

2. **Cube.**  A coordinate section again has `D=0`; first contact with the
   opposite support face kills every ray at the common marginal width.
   On any shorter interval the formula is exact.  The uniform marginal has
   an affine profile with a nonzero intercept, so the through-origin
   near-linearity hypothesis cannot persist across arbitrarily large
   multiplicative volume scales.

3. **Euclidean ball.**  For an equatorial free-boundary disk, rays near its
   boundary hit the spherical support quickly.  Their loss is exactly the
   first term of (2.8).  A fixed multiplicative volume expansion has a
   nonnegligible killed packet, consistently with strict profile curvature.

4. **Superellipsoid.**  A coordinate section has vanishing initial curvature
   in the flat directions, so local second variation can miss the geometry.
   At distance comparable with the marginal width, support hitting produces
   the missing flux loss in (2.8).  As the superellipsoid approaches a
   product cylinder this loss can vanish, which is precisely the product
   branch.

5. **Simplex.**  A generic planar truncation does not meet every active
   support face orthogonally.  Its contact-angle defect is a first-order
   term, rather than an error in (0.2).  For a genuine free-boundary leaf,
   contact with the nonsmooth skeleton is recovered by smoothing the support
   and appears as killed flux in (2.8).

6. **Thin neck.**  Large sheets on either side may have long surviving rays
   while all orientation change is concentrated in a neck of arbitrarily
   small surface capacity.  Theorem 2.1 faithfully reports only the small
   neck flux.  It cannot propagate its orientation to the large sheets.

## 6. A sharp high-rank product model

Let `mu_m` be the product of `m` one-sided exponential laws,

\[
 d\mu_m(x)=e^{-\sum_{i=1}^m x_i}1_{\{x_i\ge0\}}dx,
\]

which is isotropic after translation (each coordinate already has variance
one).  Let `q=e^{-L}` and

\[
 A_{m,q}=\bigcup_{i=1}^m\{x_i\ge L\}.
\]

Its volume and relative weighted perimeter are exactly

\[
 v(q)=1-(1-q)^m,
 \qquad P(q)=mq(1-q)^{m-1}.                               \tag{6.1}
\]

The boundary consists of `m` planar sheets.  Their normal matrix is

\[
 M(q)={P(q)\over m}I_m,\qquad
 \operatorname {rank}_{eff}M(q)=m.                       \tag{6.2}
\]

The total weighted codimension-two incidence, counted once for each pair,
is

\[
 R_2(q)={m\choose2}q^2(1-q)^{m-2},
 \qquad {R_2(q)\over P(q)}={ (m-1)q\over2(1-q)}.          \tag{6.3}
\]

Take `alpha=mq` and let `m to infinity`, then `alpha downarrow0`.  One has

\[
 v(q)=\alpha+O(\alpha^2+\alpha/m),\quad
 P(q)=\alpha+O(\alpha^2+\alpha/m),\quad
 {R_2(q)\over P(q)}={\alpha\over2}+o(\alpha).            \tag{6.4}
\]

Thus the perimeter-through-volume slope is asymptotically constant while
the ridge charge is an arbitrarily small fraction of perimeter and the
normal rank is arbitrary.

The killed tube calculation is also exact.  On sheet `i`, the outer normal
of `A_{m,q}` is `-e_i`, `lambda=1`, and `D=0`.  A ray loses uniqueness when
it becomes equally close to a second sheet.  Hence, for `0<t<L`,

\[
 R(t)=mq(1-qe^t)^{m-1},
 \qquad v(t)=1-(1-qe^t)^m,                                \tag{6.5}
\]

and

\[
 {R(t)\over R(0)}
 =\left({1-qe^t\over1-q}\right)^{m-1}
 =\exp\{-\alpha(e^t-1)+o(\alpha)\}.                     \tag{6.6}
\]

At time zero the collision rate and pair-incidence have the exact relation

\[
 -{d\over dt}\bigg|_{t=0}\log R(t)
 ={(m-1)q\over1-q}=2{R_2(q)\over P(q)}.                  \tag{6.7}
\]

Thus the apparent factor two between the first-order flux loss and the
unoriented pair-ridge count is exactly the fact that each ridge is incident
to two sheets.

For a fixed multiplicative increase of the small volume, `t=Theta(1)` and
the total collision loss is `Theta(alpha)P`, exactly matching (0.6).  It is
distributed among `Theta(m^2)` pairwise ridges and yields no coherent
direction.

This example is an exact nested competitor family, not a claim that every
`A_{m,q}` is a global Euclidean isoperimetric minimizer.  Its role is more
precise: it disproves any attempted deduction from the normal Jacobian,
near-constant perimeter/volume slope, convex support, and pairwise sheet
intersection **alone** to a dimension-free ridge charge or a single affine
marginal.  The missing quantitative object must distinguish an almost
product family like (6.1) from a genuinely irreducible high-rank family.

## 7. Implication for the physical packet

Suppose a physical boundary packet has

\[
 \operatorname {tr}M_{phys}\ge\alpha p,
 \qquad M_{phys}=\int N\otimes N\,d\nu_{phys},            \tag{7.1}
\]

and lies on a regular isoperimetric sheet satisfying Theorem 2.1 with
`P_0 asymp p`.  If

\[
 \epsilon:=\delta+{\delta b\over s}\log\beta
 < {\alpha p\over
  2C_{dom}P_0[1+(1-e^{-1})^{-1}]},                       \tag{7.2}
\]

where domination means `d nu_phys<=C_dom d sigma_mu`, then after deleting
all cut rays and the set `D>1`, a packet of
trace at least `alpha p/2` remains.  Formula (3.4) applies to it.

For the numerical `alpha=.004`, this requires an actual global profile/tube
defect below the `10^{-3}p` scale (with the displayed domination constants),
not merely a pointwise curvature estimate.  When it holds, the coherent
normal branch contradicts small `p`.  In the high-effective-rank branch,
(3.4) remains compatible with isotropy, and the product model shows why.
The next necessary lemma is therefore a dichotomy of the form

\[
 \text{high-rank long-ray packet}
 \quad\Longrightarrow\quad
 \begin{cases}
 \text{quantitative collision/conductance saving},\\
 \text{or an almost-product decomposition with a controlled 1D factor}.
 \end{cases}                                              \tag{7.3}
\]

Neither profile concavity, the free-boundary index form, nor exact planar
intersection proves (7.3).  Formula (6.6) is the sharp calibration any such
dichotomy must survive.

### 7.1 A fixed-scale numerical chain

Take

\[
                 v_0={b\over2},\qquad a={b\over4}.
\]

Then `beta=2`, `s=b/4`, and `eta=4delta`.  For
`0<delta<=10^{-2}`, Theorem 2.1 gives

\[
 {1\over(1+4\delta)c}
 \log\left(1+{1+4\delta\over1+\delta}\right)
 \le T_b\le{\log2\over c},                               \tag{7.4}
\]

so, numerically,

\[
                         {0.680\over c}\le T_b
                         \le {0.694\over c}.              \tag{7.5}
\]

The total normalized tube defect obeys

\[
 \epsilon\le(1+4\log2)\delta<3.773\delta.                \tag{7.6}
\]

With the cutoff `h=1`, the deleted base surface (cut rays together with
surviving rays having `D>1`) has measure at most

\[
 \left(1+{1\over1-e^{-1}}\right)\epsilon P_0
 <9.742\delta P_0.                                       \tag{7.7}
\]

Moreover `|lambda|T_b<=(1+4delta)log2<0.722`, and (0.8),
(7.5) imply the fully numerical estimate

\[
 \operatorname {Cov}(\mu)
 \succeq 0.0046\,c^{-3}
       \int_GN\otimes N\,d\sigma_\mu.                   \tag{7.8}
\]

Thus, in the isotropic case,

\[
       \left\|\int_GN\otimes N\,d\sigma_\mu\right\|_{op}
       <218c^3.                                          \tag{7.9}
\]

Since `cb/2<=P_0<=(1+delta)cb/2`, a coherent surviving packet satisfying
`sigma_mu(G)>=alpha P_0` and the first inequality in (3.5) forces

\[
 P_0^2\ge {\kappa\alpha b^3\over
                 1744}.                                  \tag{7.10}
\]

For a physical trace fraction `alpha=.004`, surface domination constant
one, and `P_0=p`, (7.7) retains at least half that trace only when
`delta<.002/9.742<2.06\,10^{-4}`.  To force the tube defect itself below
`6.02\,10^{-5}P_0` using (7.6) requires
`delta<1.60\,10^{-5}`.  With domination constant `C_dom` and
`P_0<=C_Pp`, the first condition is instead

\[
 \epsilon\le {\alpha\over
  2C_{dom}C_P[1+(1-e^{-1})^{-1}]}.                       \tag{7.11}
\]

These are assumptions on a global profile deficit; the local
`L^1`-in-volume curvature estimate does not automatically provide them at
the selected physical level.

## 8. Transfer gap from near-minimizing levels

For a nonstationary smooth hypersurface the exact normal Jacobian is still

\[
 j_x(t)=\exp\{H_\mu(x)t-D_x(t)\},                         \tag{8.1}
\]

but the base exponent `H_mu(x)` now depends on `x`.  There is no single
factor `e^{lambda t}` whose removal makes total flux monotone, and the
profile comparison used in (2.6) no longer controls the spatial reweighting
`e^{H_mu(x)t}`.  A small perimeter excess does not, without an Ekeland-type
selection and a normed first-variation estimate, control the oscillation of
`H_mu` on the physical packet.

To apply Theorem 2.1 to the ranked physical matrix one needs a quantitative
lemma with all of the following conclusions.  Starting from the relevant
near-minimizing level, it must produce either an exact minimizer or a
quantitatively stationary regular competitor while preserving, up to
universal losses,

1. the boundary packet mass;
2. the operator norm and effective rank of its normal matrix;
3. its domination by weighted surface measure;
4. the hard-support contact geometry; and
5. a profile deficit small enough for (7.7).

Ordinary compactness or lower semicontinuity preserves perimeter and volume
but not items 1--3.  No such stability lemma is proved here.  Therefore the
finite-distance result settles the `delta p^3` versus `delta p` scaling and
the first-contact/focal bookkeeping for exact regular leaves, but it does
not by itself yield a saving for the physical ranked matrix.
