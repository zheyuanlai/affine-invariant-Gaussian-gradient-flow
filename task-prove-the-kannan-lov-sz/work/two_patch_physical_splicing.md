# Same-level normal dispersion, cross-level rotation, and the limits of two-patch splicing

## 0. Outcome

This note audits the proposed implication

\[
 \text{high effective rank of the physical coarea normal measure}
 \quad\Longrightarrow\quad
 \text{a fixed finite two-patch perimeter saving}.                 \tag{0.1}
\]

There is a sharp obstruction before any local geometry is used: the matrix in
`fixed_scale_physical_splicing.md` is integrated over the level parameter.
Two independent normals sampled from that matrix generally lie on boundaries
of *different* sets.  A mass-preserving splice of `A_r` must instead use two
patches of the same boundary `∂A_r`.

The exact disintegration gives a useful quantitative dichotomy.  Write `a`
for the trace of the physical matrix and `R` for its effective rank.  Then

\[
 \boxed{\mathsf W+\mathsf B=1-\operatorname{tr}(Q^2)
                    \ge 1-R^{-1},}                  \tag{0.2}
\]

where `W` is the same-level projective angular dispersion and `B` is the
between-level variance of the normalized normal matrices.  At the numerical
rank `R>17`, either

\[
 \mathsf W>{8\over17}
 \quad\hbox{or}\quad
 \mathsf B>{8\over17}.                              \tag{0.3}
\]

The first branch supplies many separated same-level pairs with exactly linear
boundary-mass accounting.  The second says that the high global rank is paid
by a fixed rotation of the normal projector distribution between levels.  It
does not supply a same-level pair.

Neither branch by itself proves a finite saving.  In the first branch, the
weighted first- and second-variation formulas for two disjoint patches contain
no cross-normal angle term.  A saving needs a physical incidence object: a
ridge, a short bridge, a focal junction, or another deformation whose cost is
linear in boundary mass.  In the second branch, one needs a new theorem which
charges the level-to-level rotation to curvature, topology change, or coarea
deficit.  No such theorem follows from nesting alone, because the scalar
subweight `omega` can select a rotating cap on an otherwise smooth nested
family.

Consequently, the fixed `10^{-4}p` splice lemma is **not proved here**.  What is
proved is the exact same-level/between-level alternative, the correct reuse
normalization, and a first/second-variation no-go theorem for any proposed
argument based only on two separated normal directions.  These results reduce
the remaining geometric statement to one of two explicit physical charges.

## 1. Physical coarea disintegration

Let `E` be a Euclidean affine space of dimension `k>=2`.  Let
`dmu=varrho dx`, where `varrho=e^{-V}` and `V` is convex on `E`.  Let
`F:E->[0,1]` be locally Lipschitz and put

\[
                         A_r=\{F>r\}.                \tag{1.1}
\]

Let `omega:E->[0,1]` be measurable.  For almost every `r`, define a finite
measure on the reduced boundary by

\[
 d\sigma_r(x)=\omega(x)\varrho(x)\,
               d\mathcal H^{k-1}\!\restriction_{\partial^*A_r}(x),
 \qquad n_r(x)=\hbox{measure-theoretic outer unit normal}.          \tag{1.2}
\]

Assume

\[
 0<a:=\int_0^1 a_r\,dr<\infty,
 \qquad a_r:=\sigma_r(\partial^*A_r).                \tag{1.3}
\]

The physical normal matrix is

\[
 M=\int_0^1M_r\,dr,
 \qquad
 M_r=\int n_rn_r^T\,d\sigma_r,
 \qquad \operatorname{tr}M=a.                       \tag{1.4}
\]

This includes (1.23) of `fixed_scale_physical_splicing.md`, with
`F=T_s1_S`.  Matrix-valued measurability follows entrywise from scalar
weighted coarea.

On `{a_r>0}`, put

\[
 d\nu_r=a_r^{-1}d\sigma_r,
 \qquad Q_r=\int nn^T\,d\nu_r(n)={M_r\over a_r}.     \tag{1.5}
\]

Set `Q=M/a` and give the level interval the probability law

\[
                         d\pi(r)={a_r\over a}\,dr.    \tag{1.6}
\]

Then

\[
 Q=\mathbb E_\pi Q_r,
 \qquad Q_r\succeq0,\quad \operatorname{tr}Q_r=1,
 \qquad Q\succeq0,\quad \operatorname{tr}Q=1.       \tag{1.7}
\]

If

\[
                         R={a\over\|M\|_{op}},       \tag{1.8}
\]

then `||Q||_op=1/R` and

\[
                         \operatorname{tr}(Q^2)
 \le\|Q\|_{op}\operatorname{tr}Q={1\over R}.        \tag{1.9}
\]

### 1.1 What the global rank really gives

Sample `(R_1,N_1)` by first choosing `R_1~pi` and then choosing
`N_1~nu_{R_1}`.  Let `(R_2,N_2)` be an independent copy.  Then

\[
 \mathbb E\langle N_1,N_2\rangle^2
 =\operatorname{tr}(Q^2)\le {1\over R}.             \tag{1.10}
\]

For `R>=17`, Markov's inequality gives

\[
 \mathbb P\{ |\langle N_1,N_2\rangle|>1/2\}
 \le {4\over17},
\quad
 \mathbb P\{ |\langle N_1,N_2\rangle|\le1/2\}
 \ge {13\over17}.                                  \tag{1.11}
\]

This is a strong projective separation statement, but it is a statement
about two independently sampled *levels*.  The event `R_1=R_2` has
probability zero whenever `pi` is nonatomic.  Therefore (1.11) cannot be
fed directly into a splice of one `A_r`.

## 2. The exact within-level/between-level identity

Define

\[
 \begin{aligned}
 \mathsf W
   &:=\mathbb E_\pi\left[1-\operatorname{tr}(Q_r^2)\right],\\
 \mathsf B
   &:=\mathbb E_\pi\|Q_r-Q\|_F^2.
 \end{aligned}                                      \tag{2.1}
\]

The letter `W` stands for within-level angular dispersion; `B` stands for
between-level matrix variation.

**Lemma 2.1 (orthogonal variance decomposition).**  Under (1.3)--(1.7),

\[
 \boxed{\mathsf W+\mathsf B=1-\operatorname{tr}(Q^2).}           \tag{2.2}
\]

In particular,

\[
 \boxed{\mathsf W+\mathsf B\ge1-{1\over R}.}         \tag{2.3}
\]

**Proof.**  Since `E_pi Q_r=Q`,

\[
 \mathsf B
 =\mathbb E_\pi\operatorname{tr}(Q_r^2)
   -2\operatorname{tr}\!\left((\mathbb E_\pi Q_r)Q\right)
   +\operatorname{tr}(Q^2)
 =\mathbb E_\pi\operatorname{tr}(Q_r^2)-\operatorname{tr}(Q^2).
\]

Adding the definition of `W` proves (2.2).  Equation (2.3) follows from
(1.9).  QED.

At the physical rank in (1.21) of the fixed-scale report,

\[
 \boxed{\max\{\mathsf W,\mathsf B\}>{8\over17}.}      \tag{2.4}
\]

This dichotomy is lossless except for the final split into two alternatives.

### 2.1 The same-level branch

Conditionally on `r`, let `N,N'` be independent with law `nu_r`.  Then

\[
 \mathbb E[1-\langle N,N'\rangle^2\mid r]
 =1-\operatorname{tr}(Q_r^2).                       \tag{2.5}
\]

Thus the same-level probability measure

\[
 d\Pi_{same}(r,x,y)
 ={1\over a}\,{d\sigma_r(x)d\sigma_r(y)\over a_r}\,dr           \tag{2.6}
\]

satisfies

\[
 \int[1-\langle n_r(x),n_r(y)\rangle^2]\,d\Pi_{same}
 =\mathsf W.                                        \tag{2.7}
\]

The division by `a_r` is essential.  It makes the total mass one and gives
both marginals

\[
 {1\over a}\,d\sigma_r(x)\,dr.                      \tag{2.8}
\]

Hence (2.6) is the unique elementary product sampling with *linear* rather
than quadratic boundary-mass accounting.  The raw measure
`d sigma_r(x)d sigma_r(y)dr` reuses a level in proportion to `a_r^2` and is
not compatible with a dimension-free perimeter budget.

If `W>=8/17`, then under `Pi_same`,

\[
 \mathbb E\langle N,N'\rangle^2\le {9\over17}.       \tag{2.9}
\]

Markov's inequality at `3/4` yields

\[
 \boxed{
 \Pi_{same}\{ |\langle N,N'\rangle|\le\sqrt3/2\}
 \ge {5\over17}.}                                  \tag{2.10}
\]

Thus at least `5/17` of the linearly normalized same-level pairs have
projective angle at least `pi/6`.

This is the strongest fixed-angle conclusion available from the half of
(2.4) assigned to `W`, up to changing numerical thresholds.  It is already
enough in constants if each unit of angular energy could be converted into a
finite saving with coefficient `1/8`: indeed, using
`a>.004p`,

\[
 {a\over8}\mathsf W
 \ge {a\over17}>2.35\cdot10^{-4}p.                  \tag{2.11}
\]

The missing issue is not the numerical constant; it is the existence of a
physical, bounded-overlap splice realizing any fixed fraction of (2.11).

### 2.2 The between-level branch

Let `R,R'` be independent with law `pi`.  From the Hilbert-space variance
identity,

\[
 \mathbb E\|Q_R-Q_{R'}\|_F^2=2\mathsf B.             \tag{2.12}
\]

For positive semidefinite trace-one matrices,
`||Q_R-Q_R'||_F^2<=2`.  If `B>=8/17`, then

\[
 \mathbb E\|Q_R-Q_{R'}\|_F^2\ge {16\over17}.        \tag{2.13}
\]

For a random variable `0<=Z<=2`,
`E Z<=t+(2-t)P{Z>=t}`.  Taking `t=1/2` in (2.13) gives

\[
 \boxed{
 (\pi\otimes\pi)\left\{
 \|Q_R-Q_{R'}\|_F^2\ge {1\over2}\right\}
 \ge {5\over17}.}                                  \tag{2.14}
\]

So the alternative to same-level angular dispersion is a macroscopic
level-to-level change of the normalized physical normal matrix.  A proof of
the fixed splice lemma must charge (2.14) to a physical event occurring on
the intervening levels.  Merely pairing the two levels violates exact
mass-preservation for each `A_r`.

## 3. Why two disjoint patches have no angular second-variation term

This section gives a local no-go theorem.  It applies even when both patches
belong to the same smooth level.

Let `Σ=∂A` be a `C^2` hypersurface in a region where `V` is `C^2`.
Use the outward unit normal `n`, the convention

\[
 H=\operatorname{div}_\Sigma n,
 \qquad H_V=H-\langle\nabla V,n\rangle,              \tag{3.1}
\]

and the normal deformation `x -> x+t f(x)n(x)`.  For compactly supported
`C^1` speeds,

\[
 {d\over dt}\bigg|_{t=0}\mu(A_t)
       =\int_\Sigma f\varrho\,d\mathcal H^{k-1},     \tag{3.2}
\]

\[
 {d\over dt}\bigg|_{t=0}P_\mu(A_t)
       =\int_\Sigma H_Vf\varrho\,d\mathcal H^{k-1}. \tag{3.3}
\]

If `H_V=lambda` on `Sigma`, the second derivative of
`P_mu-lambda mu` is

\[
 \boxed{
 \mathcal Q(f)=\int_\Sigma
 \left(|\nabla_\Sigma f|^2-
       (|\mathrm {II}|^2+\nabla^2V(n,n))f^2\right)
       \varrho\,d\mathcal H^{k-1}.}                 \tag{3.4}
\]

The formula follows by differentiating (3.3), using the standard normal
variation identities

\[
 \dot n=-\nabla_\Sigma f,qquad
 \dot H=-\Delta_\Sigma f-|\mathrm {II}|^2f,qquad
 {d\over dt}\langle\nabla V,n\rangle
 =f\nabla^2V(n,n)-\langle\nabla_\Sigma V,
                              \nabla_\Sigma f\rangle,
\]

and integrating the weighted tangential Laplacian by parts.  Compact support
removes boundary terms.  Approximation extends (3.4) to the form domain.

Now choose disjoint open patches `U_1,U_2 ⊂ Σ` and functions
`f_i` supported in `U_i`.  Set

\[
 f=f_1-cf_2,
 \qquad
 c={\int f_1\varrho\over\int f_2\varrho},            \tag{3.5}
\]

assuming the denominator is nonzero.  Then the first volume derivative is
zero.  If `H_V` is constant, the first perimeter derivative is also zero, and

\[
                         \mathcal Q(f)
 =\mathcal Q(f_1)+c^2\mathcal Q(f_2).                \tag{3.6}
\]

There is no cross term because the supports are disjoint.  In particular,
the number

\[
                         \langle n(x_1),n(x_2)\rangle             \tag{3.7}
\]

does not occur anywhere in either first or second variation.

**Proposition 3.1 (two-patch angular no-go).**  No universal implication of
the form

\[
 |\langle n_1,n_2\rangle|\le c_0<1
 \quad\Longrightarrow\quad
 \hbox{a negative first or second variation supported on the two patches}
                                                               \tag{3.8}
\]

can follow from log-concavity alone.  Any valid implication must use
additional geometry connecting the patches or a nonlocal extremality
condition.

**Proof.**  Equations (3.3) and (3.6) contain no interaction between the two
normal directions.  More concretely, on a flat patch in a log-affine region,
`II=0` and `nabla^2V=0`, so

\[
                         \mathcal Q(f_i)
 =\int|\nabla_\Sigma f_i|^2\varrho\ge0.              \tag{3.9}
\]

Two such disjoint patches can have any prescribed normal angle.  Their local
mass exchange has zero first variation when their weighted mean curvatures
agree and nonnegative second variation.  QED.

Log-concavity helps the negative potential term in (3.4), because
`nabla^2V>=0`, but this contribution vanishes identically for uniform and
log-affine densities.  Those are mandatory test cases, including the cube,
simplex, and product exponential measure.

### 3.1 Exact mass preservation does not create an angle term

Suppose `f` satisfies the first-order constraint in (3.5).  Add a scalar
second-order correction on a third regular patch, or solve for one patch's
amplitude by the implicit function theorem.  This produces a deformation
with `mu(A_t)=mu(A)` exactly for small `t`.  At a stationary surface the
coefficient of `t^2` remains (3.4); the correction contributes only through
the Lagrange multiplier already subtracted in `P_mu-lambda mu`.  Therefore
exact mass preservation does not repair (3.8).

## 4. The direction of log-concavity for a bridge or chord

Let `x,y` be points at which the density is positive.  Along their segment,
log-concavity says

\[
 \varrho((1-t)x+ty)
 \ge\varrho(x)^{1-t}\varrho(y)^t,qquad0\le t\le1.   \tag{4.1}
\]

Consequently, the weighted cost of the straight chord is bounded **below**:

\[
 \int_{[x,y]}\varrho\,d\mathcal H^1
 \ge |x-y|\int_0^1\varrho(x)^{1-t}\varrho(y)^t\,dt
 =|x-y|\,L(\varrho(x),\varrho(y)),                  \tag{4.2}
\]

where

\[
 L(a,b)=
 \begin{cases}
 (a-b)/(\log a-\log b),&a\ne b,\\
 a,&a=b
 \end{cases}                                       \tag{4.3}
\]

is the logarithmic mean.

Thus log-concavity does not make a remote connecting chord cheap.  It prevents
the density from dipping below the geometric interpolation of its endpoint
values.  A chord produces a perimeter saving only when Euclidean shortening
beats this bridge cost.  At a genuine wedge ridge the endpoints are distance
`O(ell)` apart and the exact shortening is first order, as in Section 3 of
`fixed_scale_physical_splicing.md`.  For two spatially separated patches,
the chord length can be arbitrarily larger than the amount of old boundary
allocated to the operation.  Any proof that uses (4.1) in the reverse
direction is invalid.

This also explains why midpoint overlap is a constraint rather than an
automatic saving.  Log-concavity puts mass in the bridge; a separate inverse
or incidence theorem is needed to turn that bridge mass into a competing
boundary.

## 5. A physical realization of the between-level obstruction

The matrix-theoretic possibility `Q_r` rank one for almost every `r` while
`Q=E Q_r` has rank `18` is not merely an abstract labelled-graph example.  It
can be realized by smooth nested convex level sets in a log-concave
probability space.

Take `k=18`, let `mu=gamma_k` be standard Gaussian measure, choose a positive
definite diagonal matrix

\[
 H=\operatorname{diag}(1,2,3,1,\ldots,1),            \tag{5.1}
\]

and put

\[
 F(x)=\exp(-x^THx).                                  \tag{5.2}
\]

For `0<r<1`,

\[
 A_r=\{x:x^THx<-\log r\}                            \tag{5.3}
\]

is a smooth nested ellipsoid.  For any unit vector `u` and any level `r`,
the ellipsoid has a unique point whose outward normal is `u`, namely

\[
 x(r,u)=\sqrt{-\log r}\,
 {H^{-1}u\over\sqrt{u^TH^{-1}u}}.                   \tag{5.4}
\]

Choose `epsilon>0` with `sin epsilon=10^{-2}` and the eighteen unit vectors

\[
 u_1=e_1,\qquad
 u_2=\sin\epsilon\,e_1+\cos\epsilon\,e_2,\qquad
 u_i=e_i\quad(3\le i\le18).                         \tag{5.5}
\]

Their equally weighted projector matrix is

\[
 Q_*=\frac1{18}\sum_{i=1}^{18}u_iu_i^T.             \tag{5.6}
\]

The two eigenvalues of its `(e_1,e_2)` block are
`(1±sin epsilon)/18`; all other eigenvalues are `1/18`.  Hence

\[
 {1\over\|Q_*\|_{op}}
 ={18\over1+\sin\epsilon}>17.82.                    \tag{5.7}
\]

Let `I_1,...,I_18` be disjoint compact level intervals.  For `x≠0`, write

\[
                         n(x)={Hx\over|Hx|}.         \tag{5.8}
\]

For a small aperture `eta`, take a measurable selector of the explicit form

\[
 \omega_\eta(x)=\sum_{i=1}^{18}c_i(F(x))
  \mathbf1_{I_i}(F(x))\mathbf1_{\{|n(x)-u_i|<\eta\}},
 \qquad 0\le c_i\le1.                              \tag{5.9}
\]

Thus, on levels in `I_i`, `omega_eta` is supported on a small surface cap
about `x(r,u_i)`, and it vanishes on the other intervals.  Because the Gaussian
density and the surface Jacobian are continuous and positive on each compact
cap, the cap widths and scalar multipliers in `[0,1]` may be chosen so that
all eighteen intervals carry the same total `sigma`-mass.  As the cap
apertures tend to zero,

\[
 Q={M\over\operatorname{tr}M}\longrightarrow Q_*    \tag{5.10}
\]

while, for almost every selected level,

\[
                         Q_r\longrightarrow u_iu_i^T.             \tag{5.11}
\]

Thus `R>17` and `W` can be made arbitrarily small; all the charge in (2.2)
lies in `B`.

This family is not concurrent about one center.  To see this explicitly,
write `s=sin epsilon`, `c=cos epsilon`, and choose one level from each of the
first three intervals.  Up to positive scalar factors `kappa_i`, (5.4) gives

\[
 x_1=\kappa_1 e_1,\qquad
 x_2=\kappa_2(se_1+(c/2)e_2),\qquad
 x_3=\kappa_3e_3.                                   \tag{5.12}
\]

The normal lines are

\[
 L_1=x_1+\mathbb Re_1,\qquad
 L_2=x_2+\mathbb R(se_1+ce_2),\qquad
 L_3=x_3+\mathbb Re_3.                              \tag{5.13}
\]

The first and second lines meet at

\[
                         L_1\cap L_2
 =\{(\kappa_2s/2)e_1\},                              \tag{5.14}
\]

which is nonzero, whereas `L_1 ∩ L_3={0}`.  Hence the three lines have no
common point.

The construction has `C_P(gamma_k)=1` and is not claimed to arise from the
special heat-smoothed indicator or its special `omega`.  Its role is exact:
it proves that smooth nesting, convexity of every level, log-concavity of the
ambient measure, and the numerical matrix bounds alone do not convert global
rank into same-level rank or into one-center concurrency.  A successful
argument must use the near-extremal heat-flow origin of `omega`, or must prove
a charge for the between-level alternative (2.14).

## 6. What a reuse-safe physical conversion would have to prove

The same-level law (2.6) is a coupling with correct marginals, but it is not
an admissible family of disjoint surgery tubes.  For a finite splice, one
needs a nonnegative conductance measure `Gamma_r` on pairs of boundary
patches such that:

1. both marginals of `Gamma_r` are bounded by `d sigma_r`;
2. the associated ridge, bridge, or deformation tubes can be chosen with
   bounded overlap;
3. the saving assigned to `(x,y)` is bounded below by a fixed multiple of
   `1-<n_r(x),n_r(y)>^2` times the allocated conductance; and
4. the accumulated volume error has an exactly mass-preserving correction
   whose cost is a controlled fraction of the saving.

If such a family captured the full product coupling in the sense

\[
 \int [1-\langle n(x),n(y)\rangle^2]\,d\Gamma_r\,dr
 \ge {1\over2}a\mathsf W,                            \tag{6.1}
\]

and converted angular energy with coefficient `1/8`, then the same-level
branch would save at least

\[
 {a\mathsf W\over16}\ge {a\over34}
 >1.17\cdot10^{-4}p,                                \tag{6.2}
\]

which closes the numerical gap.  Hence a factor-two loss in matching and the
`1/8` bevel coefficient are both affordable.

Proposition 3.1 shows that item 3 cannot hold for arbitrary disjoint patches.
It is precisely a lower bound on physical incidence capacity.  In the facet
model of the fixed-scale report, this capacity is the ridge conductance and
Lemma 3.2 proves the conversion once its spectral gap is at least `1/4`.

For the between-level branch, a sufficient statement of comparable strength
would be a charge functional `C_rot(F,omega)` satisfying

\[
 C_{rot}(F,\omega)\ge c\,a\mathsf B                 \tag{6.3}
\]

and producing finite same-mass splices with total saving at least
`C_rot/16`.  At `B>=8/17`, this again exceeds `10^{-4}p`.  The functional
cannot depend only on the smooth second fundamental form of the reduced
boundaries: components can be born, merge, or turn through singular focal
sets, and the scalar selector `omega` can rotate without the surface itself
rotating.  A valid `C_rot` must therefore include full `BV`/focal charge or
must control the variation of the special heat-flow selector.

## 7. Mandatory model tests

### 7.1 Radial spheres

For a rotationally invariant level and unweighted full boundary measure,

\[
                         Q_r={1\over k}I.             \tag{7.1}
\]

Hence `B=0` and `W=1-1/k`: the same-level branch is maximally active.  Yet a
Euclidean sphere is volume-preserving stable in constant density, and every
normal line meets at its center.  This proves that separated normals alone
cannot imply a saving and validates the need for a concurrent radial branch.
For the Gaussian measure, radial balls are not global isoperimetric regions,
but their instability is detected by the curvature and potential terms in
(3.4), not by a generic two-disjoint-patch angle interaction.

### 7.2 Spatially separated planar patches

On a flat patch in a uniform or log-affine region,

\[
                         \mathcal Q(f)=\int|\nabla f|^2\varrho.    \tag{7.2}
\]

Several disjoint planar patches may have arbitrarily different normals, but
their local quadratic forms add as in (3.6).  There is no saving until one
uses the place where the patches terminate, meet a ridge, hit the support, or
are globally connected.  This is the exact facet-area versus ridge-capacity
gap.

### 7.3 Product exponentials and the rounded max box

For `V(x)=sum_i x_i` on the positive orthant, `nabla^2V=0`.  Therefore a
disjoint infinitesimal cap exchange has no curvature supplied by the density.
The fixed saving for the maximum box comes instead from simultaneous finite
bevels at its codimension-two ridges.  Equations (4.3)--(4.12) of
`fixed_scale_physical_splicing.md` verify that the ridge capacity is large
enough and that all multiple overlaps can be handled exactly.  This model
supports, rather than bypasses, the incidence requirement.

### 7.4 Cube

For uniform measure on a cube, a coordinate half-cube is the affine branch.
An inner box has a high-rank normal matrix, but its coordinate facets share
codimension-two ridges of positive capacity.  Simultaneous corner beveling
gives a finite saving.  Formula (3.4) alone cannot see that saving on the
interiors of the facets because both `II` and `nabla^2V` vanish there.

### 7.5 Simplex

A barycentric halfspace supplies the affine benchmark.  A homothetic inner
simplex has high projective normal rank and a connected dihedral ridge graph;
the local wedge calculation applies.  Again the saving is charged on the
codimension-two skeleton, not on two separated facet interiors.

### 7.6 Smooth nonconcurrent nested ellipsoids

The construction in Section 5 has smooth convex levels and no singular
ridges.  Its physical submeasure has high global rank purely through
between-level selection.  It activates (2.14), not (2.10).  Because the
ambient Gaussian has `C_P=1`, it belongs to the bounded-`K` side and is not a
counterexample to the sought fixed-scale theorem.  It nevertheless rules out
any proof of that theorem which uses only (1.21), coarea, nesting, and
log-concavity.

## 8. Exact remaining geometric alternatives

For the actual fixed-scale data, `a>.004p` and `R>17`.  By Lemma 2.1, every
closing proof may be organized into the following exhaustive cases.

1. **Same-level dispersion:** `W>8/17`.  Prove a bounded-overlap physical
   incidence/matching statement of the form (6.1).  The existing local bevel
   and constants then give more than `10^{-4}p` saving.
2. **Between-level rotation:** `B>8/17`.  Prove that the special
   heat-flow-generated selector and nested levels pay a full `BV`, focal, or
   critical-value charge of the form (6.3), unless the configuration reduces
   to a bounded affine, radial, or product/cylindrical branch.

The second alternative cannot be deleted.  An argument that samples two
normals directly from `M_phys` and calls them a two-patch splice silently
pairs different levels.  The first alternative also cannot be closed by
weighted second variation on disjoint patches, because Proposition 3.1
shows that their angle is absent from the quadratic form.

No genuine isotropic log-concave example with arbitrarily large `K` and all
the fixed-scale near-extremality identities is constructed here.  Such an
example would itself be a counterexample to KLS.  The result of this audit is
instead a rigorous identification of the two load-bearing physical theorems
which the current fixed-scale seed does not yet contain.
