# The common-calibration branch: a no-switching theorem and the jump-collapse obstruction

## 0. Verdict

There is one genuine gain from having the **same** calibration on all
levels, but it does not reach the retained matrix without an additional
non-collapse statement.

* On a region on which the variation of `G` is diffuse, saturation fixes
  the calibration field pointwise.  A polyhedral change of normal then
  creates a codimension-one singular part in its divergence.  Thus a
  bounded common forcing rules out the coordinate-max mechanism of both
  the product exponential and the cyclically constrained exponential.
  Proposition 3.1 below is quantitative: the forbidden singular density is
  at least

  \[
       (1-5\kappa)|u-v|\rho
  \]

  across a switch from a unit normal `u` to a unit normal `v`.

* The retained coarea matrix does not say that any of its mass is diffuse.
  It may be carried entirely by one jump.  Proposition 2.1 shows that every
  balanced Cheeger minimizer can be encoded as a two-valued function for
  which all nontrivial calibrated levels are the **same set**.  The trace
  of its matrix can be scaled to any prescribed universal multiple of
  `psi`, while its normalized angular purity is unchanged.  Consequently
  the commonness of the calibration supplies no interlevel geometry in
  this case.

There is also an important distinction between the two direct-deficit
constructions currently in the registry.  The boxed minimizer
`0 <= G <= 1` has obstacle multipliers in its calibration, and a negative
set energy `e_*` forbids the equation

\[
                    -\operatorname {div}_\mu z=\psi q
\]

without those multipliers.  Proposition 1.1 gives an isotropic
one-dimensional example attaining the sharp value `e_*=-kappa psi/2`.
The free-law Ekeland construction on `L^1(mu)/R` can instead produce a pure
bounded-divergence calibration.  The jump-collapse obstruction applies
even to that stronger construction.

Thus the proposed global-flow argument needs at least one new datum:
either a quantitative lower bound on matrix mass carried by genuinely
separated/diffuse levels, or a dimension-free inverse theorem for a single
high-rank calibrated Cheeger interface.  The latter is exactly the
remaining high-rank geometric branch, rather than a consequence of using
one common field.

## 1. Calibration conventions and the sharp obstacle test

Let

\[
 \Phi_H(p)=|p|+\kappa {p^THp\over |p|},\qquad
 H=H^T,\quad \|H\|_{\rm op}\le1,\quad 0<\kappa<1/5.       \tag{1.1}
\]

The value at zero is defined by one-homogeneity.  We use the polar
convention

\[
 z\cdot p\le \Phi_H(p),\qquad
 (z,DG)_\mu=\Phi_H(DG).                                  \tag{1.2}
\]

Changing `z` to `-z` converts this to the sign convention used in the
boxed direct-deficit audit.  Nothing below depends on that global sign.

For the boxed problem, convex duality gives, in a bounded smooth support,

\[
 \operatorname {div}_\mu z-\psi q+\eta=0,                 \tag{1.3}
\]

where `eta` is supported on the two obstacles, is nonpositive on
`{G=0}`, and is nonnegative on `{G=1}`.  For almost every active level
`B_r={G>r}`,

\[
 P_{\Phi_H,\mu}(B_r)
  =\psi\int_{B_r}q\,d\mu-\eta(B_r),\qquad
 e_*=P_{\Phi_H,\mu}(B_r)-\psi\int_{B_r}q\,d\mu
     =-\eta(\{G=1\}).                                      \tag{1.4}
\]

The next example shows that the obstacle term cannot be deleted even at
arbitrarily small `kappa`.

**Proposition 1.1 (sharp negative energy in an isotropic interval).**
Let `mu` be uniform on `[-sqrt(3),sqrt(3)]`.  Then `mu` is isotropic and

\[
                         \psi_\mu={1\over\sqrt3}.             \tag{1.5}
\]

Take `H=-I`, so `Phi_H=(1-kappa)|.|`, put

\[
 E=(0,\sqrt3),\qquad q=1_E-1_{E^c}.                          \tag{1.6}
\]

Then

\[
 \inf_B\left\{P_{\Phi_H,\mu}(B)-\psi_\mu\int_Bq\,d\mu\right\}
 =-\frac{\kappa\psi_\mu}{2},                               \tag{1.7}
\]

and the infimum is attained by `E`.  In particular there is no field
satisfying simultaneously (1.2), a pure equation with forcing `psi q`,
and calibration of `E`.

**Proof.**  The density is `1/(2 sqrt(3))`.  A relative finite-perimeter
subset of an interval with mass at most `1/2` has at least one interior
endpoint, and hence perimeter at least `1/(2 sqrt(3))`; equality is
attained by either half interval.  This proves (1.5).

For every Borel finite-perimeter set `B`, writing
`m(B)=min(mu(B),1-mu(B))`, the mean-zero and unit bounds on `q` give

\[
             \int_Bq\,d\mu\le m(B).                          \tag{1.8}
\]

Indeed the integral is at most `mu(B)` because `q<=1`, and, since
`int q=0`, it equals `-int_{B^c}q` and is at most `mu(B^c)`.  Therefore

\[
\begin{aligned}
 P_{\Phi_H,\mu}(B)-\psi_\mu\int_Bq\,d\mu
 &=(1-\kappa)P_\mu(B)-\psi_\mu\int_Bq\,d\mu\\
 &\ge -\kappa\psi_\mu m(B)
 \ge-\frac{\kappa\psi_\mu}{2}.                            \tag{1.9}
\end{aligned}
\]

For `E`, one has `mu(E)=1/2`, `int_E q=1/2`, and

\[
 P_{\Phi_H,\mu}(E)=(1-\kappa){\psi_\mu\over2},              \tag{1.10}
\]

so equality holds.  A pure calibrated divergence equation would give by
Gauss' formula

\[
 P_{\Phi_H,\mu}(E)=\psi_\mu\int_Eq\,d\mu={\psi_\mu\over2}, \tag{1.11}
\]

contradicting (1.10).  The missing charge is exactly
`kappa psi_mu/2` and is carried by the upper obstacle in (1.4).  QED.

The lower bound `e_* >= -kappa psi/2` used above is valid in every
dimension whenever `(1-kappa)|p| <= Phi_H(p)` and `|q|<=1`, `int q=0`.
Thus (1.7) also proves that this general estimate is sharp.

## 2. A common calibration can collapse to one jump

The free-law construction avoids the obstacle in Section 1.  It does not,
however, force different thresholds to be geometrically different.

For a finite-perimeter set `E`, write

\[
 M(E)=\int_{\partial^*E}n_E\otimes n_E\,d\sigma_\mu,
 \qquad Q_E={M(E)\over P_\mu(E)}.                             \tag{2.1}
\]

**Proposition 2.1 (exact jump-collapse encoding).**  Let `mu` be a
full-dimensional log-concave probability for which a balanced Cheeger
minimizer exists:

\[
 \mu(E)=\frac12,\qquad P_\mu(E)=\frac{\psi_\mu}{2}.           \tag{2.2}
\]

Fix any `a>0` and define

\[
 G=a(1_E-1_{E^c}),\qquad q=1_E-1_{E^c},\qquad
 \xi=\psi_\mu q.                                             \tag{2.3}
\]

Use `H=0`, so that `Phi_H=|.|`.  Then:

1. every Borel finite-perimeter set `B` satisfies

   \[
      P_\mu(B)-\int_B\xi\,d\mu\ge0;                         \tag{2.4}
   \]

2. every superlevel set of `G`, except at the two values `+-a`, is one of
   `emptyset`, `E`, or the full space, and each is a global minimizer in
   (2.4);
3. `q` is a norming functional for the quotient first moment,

   \[
     R(G)=\inf_c\int|G-c|\,d\mu=a,qquad
     \int qG\,d\mu=a;                                       \tag{2.5}
   \]

4. the complete variation matrix is jump variation and obeys

   \[
      M(G)=2aM(E),\qquad
      \operatorname {tr}M(G)=a\psi_\mu,qquad
      {M(G)\over\operatorname {tr}M(G)}=Q_E.                 \tag{2.6}
   \]

Under the weighted-`BV` dual-attainment hypotheses, (2.4) is represented
by one field `z` with bounded divergence `xi` which calibrates `G` and all
of its levels.  Thus choosing `a` to be any fixed universal constant makes
the matrix trace comparable to `psi`, while an arbitrary lower bound on
`1-tr(Q_E^2)` is retained exactly.

**Proof.**  Formula (1.8), now with the `q` in (2.3), gives

\[
 \int_B\xi\,d\mu
   \le\psi_\mu\min\{\mu(B),1-\mu(B)\}
   \le P_\mu(B),                                             \tag{2.7}
\]

which is (2.4).  Equality holds for `E` by (2.2), and it holds trivially
for the empty and full sets.  This proves the levelwise claim.  Since the
two values of `G` have equal mass, zero is a median and (2.5) follows
directly.  Finally

\[
                 DG=2aD1_E,                                  \tag{2.8}
\]

and the one-homogeneity of the matrix-valued variation gives (2.6).
Setwise coarea extends (2.4) to the corresponding convex inequality for
all `BV` functions.  Hence `xi` is a subgradient of total variation at
`G`; standard weighted-`BV` dual attainment gives the asserted field in
the regular compact setting, and it is exactly the field assumed in the
question.  QED.

**Remark 2.2 (balance is not essential).**  If an attained Cheeger
minimizer has mass `v<=1/2`, put

\[
 q=1_E-{v\over1-v}1_{E^c},\qquad G=a1_E.                     \tag{2.9}
\]

Then `|q|<=1`, `int q=0`, `int_E q=v=m(v)`, and the proof above is
unchanged.  In this normalization

\[
 R(G)=av,\qquad M(G)=aM(E),\qquad
 \operatorname {tr}M(G)=a\psi_\mu v.                        \tag{2.10}
\]

Choosing `a=tau/v` makes the trace exactly `tau psi_mu` for any prescribed
`tau>0`.  Thus even the trace normalization does not exclude a one-level
collapse.  The balanced formulation was used in Proposition 2.1 only to
make the median sign and constants especially transparent.

**Consequence 2.3.**  The numerical statements

\[
 \operatorname {tr}M(G)\asymp\psi_\mu,
 \qquad 1-\operatorname {tr}
    \left({M(G)\over\operatorname {tr}M(G)}\right)^2\ge\omega
                                                                    \tag{2.11}
\]

do not force a positive-lapse foliation, distinct levels, an absolutely
continuous part of `DG`, or any motion of the normal law in the level
parameter.  They are compatible with all matrix mass sitting on one
interface.  Therefore the Jacobi identity for a common CMC foliation and
the flow-Jacobian identity for a diffuse calibration cannot be applied to
(2.11) without a separate non-collapse lemma.

This is not a counterexample to KLS.  It is a precise equivalence
obstruction: on the jump branch, a putative inverse from the assumptions
in the question must prove directly that every balanced high-purity
Cheeger interface in an isotropic log-concave measure has universal
perimeter.  That statement is already the unresolved high-rank branch.

The heat comparator does not automatically remove this possibility.  For
a smooth finite-perimeter `E`, strict `BV` convergence of heat
regularizations gives

\[
       F_s\longrightarrow1_E\quad\hbox{in }L^1,qquad
       M(F_s)\longrightarrow M(E),qquad
       \operatorname {TV}(F_s)\longrightarrow P(E)           \tag{2.12}
\]

as `s` decreases to zero (after the usual clipping).  If `E` is a Cheeger
minimizer, its coarea Cheeger deficit also tends to zero.  Hence a matrix
penalty comparing with a very small-scale heat smoothing cannot, by
itself, impose diffuse variation on the replacement.

## 3. What the same field does rule out: quantitative no-switching

The preceding obstruction should not obscure a useful new rigidity.  It
is strongest precisely on the part of the variation which fills a band.

For a unit vector `u`, direct differentiation of (1.1) gives

\[
 D\Phi_H(u)=u+\kappa B_H(u),\qquad
 B_H(u)=2Hu-(u^THu)u.                                      \tag{3.1}
\]

For unit `u,v`,

\[
       |B_H(u)-B_H(v)|\le5|u-v|.                             \tag{3.2}
\]

Indeed the `2H` term costs at most `2|u-v|`; writing
`alpha=u^THu`, `beta=v^THv`, the remaining difference is
`alpha(u-v)+(alpha-beta)v`, whose norm is at most `3|u-v|`.
Consequently

\[
 [D\Phi_H(u)-D\Phi_H(v)]\cdot(u-v)
       \ge(1-5\kappa)|u-v|^2.                               \tag{3.3}
\]

**Proposition 3.1 (a bounded common divergence forbids a polyhedral
normal switch).**  Let `U` be open, let `rho` be positive and `C^1` on
`U`, and put `dmu=rho dx`.  Let

\[
 L(x)=\max_{1\le i\le N}\{u_i\cdot x+b_i\},                 \tag{3.4}
\]

where the `u_i` are unit vectors.  Let `g` be `C^1` and satisfy `g'>0`
on an interval `J`.  Suppose that indices `i != j` have an open regular
tie sheet

\[
 S_{ij}=\{u_i\cdot x+b_i=u_j\cdot x+b_j
        >u_k\cdot x+b_k\ (k\ne i,j),\ L(x)\in J\}\subset U. \tag{3.5}
\]

If a field `z` satisfies the polar inequality in (1.2) and calibrates
`G=g(L)` almost everywhere on `{L in J}`, then on the two cells adjacent
to `S_ij`,

\[
                    z=D\Phi_H(u_i),\qquad
                    z=D\Phi_H(u_j),                           \tag{3.6}
\]

respectively.  Its weighted distributional divergence has on `S_ij` a
singular part whose total-variation density is at least

\[
                  (1-5\kappa)|u_i-u_j|\rho.                  \tag{3.7}
\]

In particular, if `div_mu z` belongs to `L^infty(mu)` on `{L in J}`,
then no such sheet with `u_i != u_j` can occur.

**Proof.**  In the open cell where `i` is uniquely active,

\[
        \nabla G=g'(L)u_i.                                   \tag{3.8}
\]

The polar inequality and equality in (1.2), together with differentiability
and strict convexity of `Phi_H`, force the unique supporting vector
`z=D Phi_H(u_i)`.  The same holds in the `j` cell.

Let `d=u_i-u_j` and orient
`nu=d/|d|` from the `j` cell into the `i` cell.  The jump formula for a
piecewise bounded vector field gives the singular part

\[
 [\operatorname {div}(\rho z)]^s
 =\rho[D\Phi_H(u_i)-D\Phi_H(u_j)]\cdot\nu\,
      \mathcal H^{n-1}\mathbin\vrule height 1.4ex depth -0.3ex
      width 0.07ex\vrule height 0.07ex depth -0.02ex width 0.8ex S_{ij}.
                                                                    \tag{3.9}
\]

Equations (3.3) and (3.9) give (3.7).  On the other hand,
`div_mu z in L^infty(mu)` means
`div(rho z)=(div_mu z)rho dx`, which has no singular part.  This is a
contradiction unless `u_i=u_j`.  QED.

The conclusion is local.  It needs neither covariance nor tube coverage,
and it remains valid inside an arbitrary convex hard support because the
sheet in (3.5) is required to lie in its interior.

There is an analogous complete rigidity statement for a smooth band.  It
is useful for locating exactly where high rank can hide.

**Proposition 3.2 (exact smooth-band rigidity).**  Let `Omega` have `C^2`
convex boundary, let `V` be `C^2` and convex, and let `Phi` be a fixed
smooth uniformly elliptic one-homogeneous anisotropy.  Suppose an open
connected band `U` is foliated by a fixed-topology family of `C^2` nested
leaves `Sigma_v`, parametrized by enclosed weighted volume, with positive
Wulff lapse `f`.  Assume the leaves satisfy the natural Young condition at
the hard wall and all have the same weighted anisotropic mean curvature
`lambda_0`.  Then

\[
 \nabla_{\Sigma_v}f=0,\qquad S_v=0,\qquad
 \nabla^2V[D\Phi(N_v),D\Phi(N_v)]=0                         \tag{3.10}
\]

on every leaf, and the support-contact curvature vanishes.  On each
connected component the leaves are parallel hyperplane pieces and their
normal matrix has rank one.

**Proof.**  In Wulff gauge, the Jacobi charge identity is

\[
\begin{aligned}
 -\lambda'(v)\int_{\Sigma_v}{1\over f}\,d\sigma_\Phi
 &=\int_{\Sigma_v}
   \langle A_\Phi\nabla\log f,\nabla\log f\rangle
       \,d\sigma_\Phi\\
 &\quad+\int_{\Sigma_v}\left\{
   \operatorname {tr}[(D^2\Phi(N)S)^2]
   +\nabla^2V[D\Phi(N),D\Phi(N)]\right\}\,d\sigma_\Phi\\
 &\quad+\int_{\partial\Sigma_v}b_\Phi\,d\tau_\Phi.       \tag{3.11}
\end{aligned}
\]

Here `A_Phi` is positive definite, `b_Phi>=0` by convexity of the support,
and `D^2 V>=0`.  The operator `D^2 Phi(N)` is positive definite on the
tangent space; hence `tr[(D^2 Phi(N)S)^2]` is nonnegative and vanishes only
when `S=0` (conjugate by the positive square root of `D^2 Phi(N)`).  Since
`lambda'(v)=0`, every nonnegative term in (3.11) vanishes.  This gives
(3.10) and `b_Phi=0`.  The normal variation formula in Wulff gauge is

\[
                         D_vN=-\Phi(N)\nabla_{\Sigma_v}f,     \tag{3.12}
\]

so the normal is also constant from leaf to leaf.  A connected hypersurface
with `S=0` is contained in a hyperplane, proving the final assertion.  QED.

In the pure free-law calibration, a band on which the median subgradient
is `q=+1` or `q=-1` has the common curvature `+psi` or `-psi`; therefore
Proposition 3.2 applies whenever that band has the stated positive-lapse
regularity.  It shows that smooth turning cannot carry the retained
high-rank matrix in the exact-forcing case.  High rank must then be in a
jump, a critical/singular/contact event, or a mixture of disconnected flat
components.  The Ekeland error `xi=psi q+e` with merely
`||e||_infty<=delta` does not satisfy the constant-curvature hypothesis:
without control of derivatives or traces of `e`, (3.11) cannot simply be
perturbed by `O(delta)`.

## 4. Product and cyclic stress tests

### 4.1 Product exponentials

Let

\[
 d\mu(x)=e^{-\sum_{i=1}^m x_i}1_{(0,\infty)^m}(x)\,dx.       \tag{4.1}
\]

Translation by `(1,...,1)` makes this measure isotropic.  On a band on
which `g'>0`, take

\[
                   G(x)=g(\max_i x_i).                        \tag{4.2}
\]

The cells have normals `e_i`.  Every tie sheet
`x_i=x_j>max_{k notin {i,j}}x_k` in the interior carries, by Proposition
3.1, singular divergence density at least

\[
                         \sqrt2(1-5\kappa)\rho.               \tag{4.3}
\]

Thus the continuous coordinate-max foliation cannot be calibrated by one
field with bounded forcing.  The absolutely continuous part
`-grad V dot z` cannot cancel (4.3).

At a single median threshold `L_m`, however, put

\[
 A_m=\{\max_i x_i\le L_m\},\qquad
 (1-e^{-L_m})^m={1\over2}.                                  \tag{4.4}
\]

Its perimeter and normal matrix are

\[
 p_m={m\over2}(2^{1/m}-1),\qquad
 M(A_m)={p_m\over m}I_m,                                    \tag{4.5}
\]

and

\[
 {\log2\over2}\le p_m\le{1\over2},\qquad
 1-\operatorname {tr}\left({M(A_m)\over p_m}\right)^2
       =1-{1\over m}.                                       \tag{4.6}
\]

For the two-valued function `1_{A_m}`, saturation constrains the trace of
`z` only on this one boundary.  It does not force `z=D Phi(e_i)` throughout
the max cells, so the interior tie-sheet argument disappears.  This is
exactly why a lower bound on coarea matrix trace and purity does not turn
Proposition 3.1 into an inverse theorem.  The set `A_m` is used here as a
stress test of the geometry, not asserted to be an exact Cheeger minimizer.

The balanced-tube calculation gives the same message: the facets have
full rank and short killed-ray loss, but the universal value of `p_m`
comes from completion of the coordinate slices, not from the covariance
current.  A common diffuse calibration would forbid the switches; a jump
calibration would still require that global completion argument.

### 4.2 Cyclically constrained exponentials

Now restrict (4.1) to

\[
 \Omega_a=\{x_i\ge0,\ x_i+x_{i+1}\ge a\quad(i\bmod m)\},
 \qquad a=m^{-6},                                           \tag{4.7}
\]

and renormalize.  The support coupling does not remove the obstruction.
For any large `L>a`, there are open tie sheets

\[
 x_i=x_j=L>\max_{k\ne i,j}x_k                                \tag{4.8}
\]

lying strictly inside `Omega_a`.  Hence the same singular density (4.3)
is present for a diffuse max foliation.  Whitening applies a uniformly
conditioned linear map; strict monotonicity of the transformed constant
anisotropy still gives a universal positive jump charge.

On the other hand, the balanced single-threshold set has perimeter between
`.33` and `.44`, a rank-`m` normal matrix, and killed loss `O(h)` before
distance `h/p`, as computed in `balanced_tube_inverse.md`.  Again this is
compatible with one jump and therefore does not test commonness across
distinct levels.

## 5. The exact missing datum

The common field becomes useful only on variation which constrains it in
an ambient band.  The present retained quantity

\[
 M(G)=\int \sigma_G\otimes\sigma_G\,d|DG|_\mu               \tag{5.1}
\]

does not distinguish the absolutely continuous, Cantor, and jump parts of
`DG`.  Proposition 2.1 shows that no inequality of the form

\[
                  |D^jG|_\mu\le C\,J(G)                      \tag{5.2}
\]

can hold: at an exact Cheeger jump the right side is zero and the left
side is the whole variation.

There are therefore two honest ways to continue this route.

1. **Non-collapse plus flow/rigidity.**  Prove that a universal fraction of
   the retained high-rank matrix is carried by a family of spatially
   separated levels on which the common field is fixed in an ambient band.
   A positive-lapse or quantitative symmetric-difference lower bound must
   be part of the conclusion.  Proposition 3.1 would then remove the
   polyhedral product branch, while a smooth turning estimate would still
   be needed.

2. **One-interface inverse.**  Prove directly that if a balanced set is a
   global minimizer of

   \[
       P_{\Phi_H,\mu}(B)-\int_B\xi\,d\mu,qquad
       \|\xi\|_\infty\le(1+o(1))\psi_\mu,                    \tag{5.3}
   \]

   and its normalized normal matrix has purity at least a universal
   `omega>0`, then `psi_mu>=c(omega)>0`.  This statement covers the jump
   encoding in Proposition 2.1.  It is the global
   completion/overlap theorem isolated in `balanced_tube_inverse.md`; the
   same calibration does not prove it automatically.

In particular, neither the pure common calibration nor the obstacle
version should be advertised as the final high-rank inverse.  What has
been proved here is a dimension-free, quantitative exclusion of abrupt
normal switching on the diffuse branch and a sharp demonstration of why
that exclusion does not yet reach the matrix supplied by the direct
replacement.
