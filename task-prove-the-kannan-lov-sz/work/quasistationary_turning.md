# Quasistationarity does not control Plateau turning

## 0. Outcome

The estimate proposed in (8.4) of `large_fill_capacity.md` is false as an
unconditional consequence of either scalar proximal quasiminimality or even
exact isoperimetric stability.  The failure is not the arbitrary turning tube
from that note.  It occurs for a fixed isotropic log-concave probability and
a sequence of **exact isoperimetric regions whose Cheeger ratios tend to the
Cheeger constant**.

The countermodel is the uniform probability on an isotropic Euclidean ball.
Its relative isoperimetric regions below half volume are spherical caps whose
free boundaries meet the support orthogonally.  Split one such free boundary
by a smooth, normal-defined trace which is not the trace of any affine
halfspace.  A least-perimeter binary filling of this trace exists.  If its
calibration/Plateau excess vanished, equality in weighted Gauss--Green would
force the filling interface to have one constant normal and hence force the
prescribed trace to be an affine-halfspace trace.  Therefore its Plateau
excess is strictly positive.  On the other hand,

\[
 P_\mu(A)-I_\mu(\mu(A))=0
 \quad\hbox{and}\quad
 \mathcal S_{\rm bevel}=0.                              \tag{0.1}
\]

The second equality is literal: the free boundary is smooth across the
artificial phase trace, and exact fixed-volume minimality rules out every
realized perimeter saving.  Thus

\[
 \varepsilon+\int |\operatorname {turn}Z|^2d|Z|
 \le C\{d+\mathcal S_{\rm bevel}\}                       \tag{0.2}
\]

fails for every finite `C`, already because `epsilon>0` and the right side is
zero.

This does not refute the large-fill trichotomy.  All normal lines of a
spherical cap meet at the center of its defining sphere, so the example is
exactly in the **concurrent escape branch** of item 3 in Section 8 of
`large_fill_capacity.md`.  It proves that this branch must be detected before
any estimate like (0.2) is invoked.  Concurrence cannot merely be added as an
error term: its defect is zero in the example while the phase-filling excess
is positive.

Two further audits are decisive.

1. A scalar `Lambda`-proximal minimizer gives an `L^infinity` first-variation
   bound, but the symmetric-difference penalty is first order in the size of
   a variation.  It gives no lower bound on the Jacobi form.
2. The turning-tube flow has an exactly computable turning energy, but a long
   tube has isoperimetric deficit comparable to its whole perimeter.  It is
   therefore not a counterexample to a statement which genuinely uses the
   profile deficit.

The exact valid conclusion is consequently negative but structural: scalar
near-stationarity cannot supply (8.4).  A usable lemma must be an alternative
which exits immediately in the parallel/concurrent/orthogonal-radial cases,
and only estimates Plateau excess and turning after those geometries have
been quantitatively excluded.

## 1. What scalar proximal minimization actually implies

Let `mu=e^{-V}dx/Z` on a smooth convex support `Omega`, and let `B` be a
regular fixed-volume minimizer of

\[
 E\longmapsto P_\mu(E)+\Lambda\mu(E\mathbin\triangle A). \tag{1.1}
\]

The same argument applies on the regular part of a `BV` minimizer.  If `u`
is the normal speed of a smooth flow and its first-order volume change is
zero, then

\[
 \left|\int_{\partial B}(H_\mu-\lambda)u\,d\sigma_\mu\right|
 \le \Lambda\int_{\partial B}|u|\,d\sigma_\mu             \tag{1.2}
\]

for a scalar `lambda`.  Equivalently,

\[
                         |H_\mu-\lambda|\le\Lambda        \tag{1.3}
\]

almost everywhere on every regular connected piece.  This is the
dimension-free conclusion proved in `rank_preserving_eikeland.md`.

There is no corresponding dimension-free second-variation conclusion.
Indeed, for a volume-preserving flow `B_t` with initial normal speed `u`,

\[
 \mu(B_t\mathbin\triangle B)
   =|t|\int_{\partial B}|u|\,d\sigma_\mu+o(|t|).           \tag{1.4}
\]

Proximal minimality gives

\[
 P_\mu(B_t)-P_\mu(B)
 \ge-\Lambda |t|\int_{\partial B}|u|\,d\sigma_\mu
       +o(|t|).                                           \tag{1.5}
\]

Using both signs of `t` proves (1.2).  If the first variation happens to
vanish and

\[
 P_\mu(B_t)-P_\mu(B)={t^2\over2}Q_B(u)+o(t^2),            \tag{1.6}
\]

then division of (1.5) by `t^2` yields only

\[
 Q_B(u)\ge-{2\Lambda\over |t|}
           \int_{\partial B}|u|\,d\sigma_\mu+o(|t|^{-1}), \tag{1.7}
\]

which is vacuous as `t` tends to zero.  The cusp in the fidelity term is
essential; it is not a regularity issue.

Nor does a scalar profile deficit produce an infinitesimal bound.  Put

\[
 d_B=P_\mu(B)-I_\mu(\mu(B)).                              \tag{1.8}
\]

For an exactly volume-preserving flow,

\[
 P_\mu(B_t)-P_\mu(B)\ge-d_B.                             \tag{1.9}
\]

After (1.6), this gives at a finite nonzero `t`

\[
 Q_B(u)\ge-{2d_B\over t^2}+o(1),                         \tag{1.10}
\]

but no limit inequality.  Optimizing (1.10) requires a uniform cubic
remainder, hence a reach, curvature, and density-modulus bound.  None is
available for the physical heat levels, and none follows from (1.1).

In particular, the scalar proximal construction cannot be used to assert
that a negative Jacobi direction costs `O(d_B)`.  It provides first-order
CMC control only.

## 2. The exact free-boundary second variation

The obstruction persists if proximal quasistationarity is replaced by exact
isoperimetric minimality.  We record the full hard-support formula to make
clear that no omitted contact term repairs it.

Let `Sigma=partial B cap int(Omega)` be `C^2`, suppose that it meets the
`C^2` boundary of the convex support orthogonally, and use

\[
 \mathrm {II}_{\partial\Omega}(Y,Y)
       =\langle D_Y\nu_\Omega,Y\rangle\ge0.               \tag{2.1}
\]

If `H_mu=lambda`, the second variation of `P_mu-lambda mu` with normal speed
`u` is

\[
 \boxed{
 Q_\Sigma(u)=
 \int_\Sigma\bigl(|\nabla_\Sigma u|^2
  -( |\mathrm {II}_\Sigma|^2+\nabla^2V(N,N))u^2\bigr)d\sigma_\mu
 -\int_{\partial\Sigma}
   \mathrm {II}_{\partial\Omega}(N,N)u^2d\tau_\mu .}    \tag{2.2}
\]

The last term includes the acceleration required to keep the contact curve
on `partial Omega`.  An exact fixed-volume minimizer satisfies

\[
 Q_\Sigma(u)\ge0
 \quad\hbox{whenever}\quad
 \int_\Sigma u\,d\sigma_\mu=0.                          \tag{2.3}
\]

Formula (2.2) concerns deformations of the **exposed boundary**.  A Plateau
interface inserted inside `B` is not a deformation speed `u` on `Sigma`.
There is no term in (2.2) containing its area, its calibration excess, or the
turning of a max-flow certificate.  The example below shows that this formal
separation is sharp.

## 3. An isotropic near-Cheeger countermodel

### 3.1 The measure and its exact isoperimetric regions

Fix `n=3` (the displayed cap formulas remain valid for every `n>=3`) and put

\[
 R=\sqrt{n+2},\qquad
 d\mu={1_{B_R}\over\kappa_nR^n}\,dx.                    \tag{3.1}
\]

The measure is log-concave, has barycenter zero, and

\[
 \int x_ix_jd\mu={R^2\over n+2}\delta_{ij}=\delta_{ij};  \tag{3.2}
\]

hence it is isotropic.  Its perimeter is relative perimeter in `B_R`, exactly
the Minkowski boundary measure from the task statement.

We use the Burago--Maz'ya relative isoperimetric theorem for a Euclidean ball
(also obtained independently by Bokowski--Sperner and Almgren):
for every `v in (0,1/2)`, an isoperimetric region of volume `v` is, up to a
rotation,

\[
 A_a=B_R\cap B_\rho(ae_1),
 \qquad a>R,\qquad \rho=\sqrt{a^2-R^2},                  \tag{3.3}
\]

with `a` uniquely chosen by `mu(A_a)=v`.  Conversely these orthogonal
spherical caps are minimizers.  The theorem is the relative version of the
spherical-cap symmetrization theorem: a perimeter minimizer in a ball has a
constant-mean-curvature free boundary and orthogonal contact, and spherical
symmetrization gives the cap (3.3).  Its hypotheses here are precisely a
Euclidean ball, constant positive density, relative perimeter, and prescribed
volume.  No statement about general convex bodies is being used.

The two spheres in (3.3) meet orthogonally.  Indeed, at an intersection
point `x`,

\[
 |x|=R,\quad |x-ae_1|^2=a^2-R^2
 \quad\Longrightarrow\quad
 (x-ae_1)\cdot x=0.                                     \tag{3.4}
\]

Write

\[
 \Sigma_a=\partial B_\rho(ae_1)\cap B_R,\qquad
 N_a(x)={x-ae_1\over\rho}.                              \tag{3.5}
\]

The normal lines of `Sigma_a` are exactly concurrent at `ae_1`.

This family contains a near-Cheeger sequence without invoking any unproved
KLS statement.  Let `I` be the relative profile and put

\[
 q(v)={I(v)\over v},\qquad 0<v\le1/2.                   \tag{3.6}
\]

By the definition of the profile,

\[
 \psi_\mu=\inf_{0<v\le1/2}q(v).                         \tag{3.7}
\]

The cap formula makes `I` continuous on `(0,1/2]`; the local relative
isoperimetric inequality gives $q(v)\to+\infty$ as $v\downarrow0$.  Choose
`v_j in (0,1/2)` with

\[
                         q(v_j)\longrightarrow\psi_\mu. \tag{3.8}
\]

The corresponding caps `A_j` are therefore exact profile minimizers and a
near-Cheeger sequence in the literal ratio sense:

\[
 P_\mu(A_j)=I(v_j),\qquad
 {P_\mu(A_j)\over\mu(A_j)}\longrightarrow\psi_\mu.       \tag{3.9}
\]

If the infimum in (3.7) is attained only at half volume, take
$v_j\uparrow1/2$
from below; the caps then converge to a half-ball but remain curved for every
finite `j`.

### 3.2 A normal-defined phase trace

Fix one cap `A=A_a`.  The Gauss map `N_a` is an affine rescaling of
`Sigma_a` onto an open subset of the unit sphere.  Choose a relatively
compact smooth domain `D`
in that Gauss image such that the smooth closed curve

\[
 C=\partial\Gamma_+,\qquad
 \Gamma_+=N_a^{-1}(D),                                  \tag{3.10}
\]

contains no nonempty open arc lying in an affine hyperplane.  For example,
in a small spherical chart take a real-analytic sinusoidally perturbed
circle; a real-analytic planar arc would force the entire curve to be planar,
which this choice is not.  Put

\[
                         \Gamma_-=\Sigma_a\setminus
                                  \overline{\Gamma_+}.   \tag{3.11}
\]

Thus the phase is a function of the exposed normal itself; it is not an
unrelated checkerboard label.

Consider binary Caccioppoli partitions `(E_+,E_-)` of `A` with traces one
and zero on `Gamma_+` and `Gamma_-`, respectively.  The support-contact part
of `partial A` is left free.  Define

\[
 L=\inf\sigma_\mu(\partial^*E_+\cap A^{(1)}).           \tag{3.12}
\]

The cap `A`, being an intersection of two balls, is a bounded Lipschitz
domain.  The relative codimension-one Plateau problem with fixed boundary
`C` has a mass minimizer by integral-current compactness.  In codimension one
and with `Z_2` coefficients it is the interior boundary of a Caccioppoli
partition and has the required trace.  Equivalently one may use the relaxed
`BV` boundary-penalty formulation from Theorem 2.1 of
`large_fill_capacity.md`.  This formulation matters: ordinary interior
`L^1` convergence alone does not preserve a prescribed `BV` trace because a
boundary layer can collapse.  The relative-current formulation retains `C`
as its boundary and covers that collapse.  The prescribed trace is
nonconstant, so `L>0`.

For a minimizing cell define the corrected flux exactly as in Section 3 of
`large_fill_capacity.md`:

\[
 F_+=\int_{\Gamma_+}N_a\,d\sigma_\mu
     +\int_{C_+}\nu_{B_R}\,d\sigma_\mu,                 \tag{3.13}
\]

where `C_+=partial^*E_+ cap partial B_R`; the potential term is zero.  The
weighted Gauss--Green formula gives

\[
 \int_{\Lambda}N_+\,d\sigma_\mu=-F_+,
 \qquad |F_+|\le L,                                     \tag{3.14}
\]

where `Lambda=partial^*E_+ cap A^{(1)}`.  The two cell fluxes are opposite,
so the total calibration excess is

\[
                         \varepsilon=2(L-|F_+|).         \tag{3.15}
\]

**Lemma 3.1 (strict Plateau excess).**  For the trace (3.10),

\[
                              \varepsilon>0.             \tag{3.16}
\]

**Proof.**  Suppose equality held in (3.14).  Since `L>0`, `F_+` is
nonzero, and equality in the triangle inequality implies

\[
 N_+(x)=-{F_+\over|F_+|}=:u
 \quad\text{for }\sigma_\mu\text{-almost every }x\in\Lambda. \tag{3.17}
\]

Equivalently, the distributional derivative of `1_{E_+}` in the interior of
the convex domain `A` is everywhere a nonnegative multiple of the fixed
vector `-u`.  Slicing in every direction orthogonal to `u` shows that
`1_{E_+}` depends, in the interior, only on the scalar `u . x`; slicing in
the `u` direction and using the fixed sign shows that this one-dimensional
`{0,1}`-valued function is monotone.  Hence, modulo null sets,

\[
                         E_+=A\cap\{u\cdot x<t\}          \tag{3.18}
\]

or its complement, for one scalar `t`.  The transition set of its trace on
`Sigma_a` lies in the single affine plane `{u . x=t}`.  This contradicts the
choice of the nonplanar curve `C` in (3.10).  Therefore `L>|F_+|`, proving
(3.16).  QED.

The argument includes every support-contact contribution through (3.13).
Omitting that term would give the wrong flux and is not needed.

### 3.3 Failure of the proposed estimate

For every cap in the near-Cheeger sequence (3.9),

\[
                         d=P_\mu(A)-I(\mu(A))=0.          \tag{3.19}
\]

The artificial phase trace `C` is not a ridge of `partial A`.  Both one-sided
surface normals equal the same smooth vector `N_a` along `C`; the wedge angle
is zero and the first-order bevel saving is zero.  More invariantly, because
`A` is a fixed-volume perimeter minimizer, no same-volume finite competitor
has positive realized saving.  Hence

\[
                         \mathcal S_{\rm bevel}=0.        \tag{3.20}
\]

Equations (3.16), (3.19), and (3.20) contradict (0.2) for every finite `C`.
This is an exact counterexample, not a limiting failure.

If a proposed graph proxy assigns a positive `S_bevel` to (3.10) by replacing
each smooth patch by one representative normal, that proxy is not a realized
bevel saving.  Exact minimality proves that no theorem lower-bounding an
actual saving by that proxy can hold on this example.

## 4. The cap passes the full stationarity audit

For the cap, `V=0`,

\[
 |\mathrm {II}_{\Sigma_a}|^2={n-1\over\rho^2},
 \qquad
 \mathrm {II}_{\partial B_R}(N_a,N_a)={1\over R}.        \tag{4.1}
\]

The second equality uses the orthogonality (3.4), which makes `N_a` tangent
to `partial B_R` at contact.  With the harmless common normalization
`(kappa_nR^n)^{-1}`, (2.2) becomes

\[
 Q_{\Sigma_a}(u)=
 \int_{\Sigma_a}\left(|\nabla_{\Sigma_a}u|^2
       -{n-1\over\rho^2}u^2\right)d\sigma_\mu
 -{1\over R}\int_{\partial\Sigma_a}u^2d\tau_\mu.        \tag{4.2}
\]

Exact isoperimetric minimality gives

\[
 Q_{\Sigma_a}(u)\ge0
 \quad\text{for all smooth }u\text{ with }
 \int_{\Sigma_a}u\,d\sigma_\mu=0.                      \tag{4.3}
\]

Thus the countermodel has constant weighted mean curvature, the correct
orthogonal contact angle, and the full nonnegative Jacobi form.  Its positive
Plateau excess is nevertheless untouched.  Any proof of (8.4) that invokes
only (1.3), (4.2), or (4.3) necessarily inserts an implication which this
example disproves.

The example also explains the geometry of the missing alternative:

\[
 (I-N_aN_a^T)(x-ae_1)=0
 \quad\text{for every }x\in\Sigma_a.                    \tag{4.4}
\]

This is exact normal-line concurrence.  It is the third outcome of the
large-fill inverse, not an error to be bounded by the scalar deficit.

## 5. Exact audit of the turning tube

Let `gamma:[0,L] -> R^n` be unit speed, let
`epsilon |gamma''| <= 1/4`, and use
the Bishop tube coordinates from Proposition 7.1 of
`large_fill_capacity.md`:

\[
 \Psi(s,z)=\gamma(s)+\sum_{a=1}^{n-1}z_ae_a(s),
 \qquad J(s,z)=1-\sum_a\kappa_a(s)z_a.                   \tag{5.1}
\]

The optimal unit flow is

\[
                         Z(\Psi(s,z))=\gamma'(s).         \tag{5.2}
\]

Since the physical vector `Z` equals `J^{-1} partial_s` in these coordinates,

\[
 (Z\cdot\nabla)Z={\gamma''(s)\over J(s,z)}.              \tag{5.3}
\]

Writing `A_{n-1}=kappa_{n-1}epsilon^{n-1}`, its natural squared turning
energy is therefore

\[
 \mathcal T(Z)
 =\int_{T_\epsilon(\gamma)}|(Z\cdot\nabla)Z|^2dx
 =\int_0^L|\gamma''(s)|^2
      \int_{|z|<\epsilon}{dz\over J(s,z)}\,ds.           \tag{5.4}
\]

Under the displayed curvature bound,

\[
 {4\over5}A_{n-1}\int_0^L|\gamma''|^2ds
 \le\mathcal T(Z)\le
 {4\over3}A_{n-1}\int_0^L|\gamma''|^2ds.                \tag{5.5}
\]

Thus the example genuinely has the advertised turning; it is not a defect
of terminology.

It is also quantitatively far from fixed-volume perimeter minimization.  By
symmetry of the normal disks, its volume and finite-perimeter boundary are
exactly

\[
 |T_\epsilon(\gamma)|=A_{n-1}L,
 \qquad
 P(T_\epsilon(\gamma))
 ={n-1\over\epsilon}A_{n-1}L+2A_{n-1}.                  \tag{5.6}
\]

A Euclidean ball of the same volume has perimeter

\[
 P_{ball}=n\kappa_n^{1/n}(A_{n-1}L)^{(n-1)/n}.           \tag{5.7}
\]

Consequently, as `L/epsilon` tends to infinity,

\[
 {P_{ball}\over P(T_\epsilon(\gamma))}
 \le C_n\left({\epsilon\over L}\right)^{1/n},
 \qquad
 {P(T_\epsilon(\gamma))-P_{ball}\over P(T_\epsilon(\gamma))}
 \longrightarrow1.                                     \tag{5.8}
\]

The same comparison works inside a convex region on which a smooth
log-concave density is nearly constant.  Hence the turning tube has profile
deficit of the order of its whole perimeter.  The term `d` in (0.2) is large
enough to pay for it.

The lateral weighted mean curvature also shows why first-order stationarity
does not rescue the tube.  In the unweighted case its principal curvatures
are `n-2` cross-sectional curvatures of size `1/epsilon` and the axial
curvature

\[
                    -{\kappa(s)\cdot\omega
                         \over1-\epsilon\kappa(s)\cdot\omega}. \tag{5.9}
\]

Endpoint smoothing introduces another curvature regime of size
`1/epsilon`.  A closed gently turning tube can make the oscillatory part in
(5.9) small, but (5.8) remains: local CMC accuracy does not imply small
global profile deficit.

## 6. Formal consequence for the large-fill program

There are three distinct pieces of information:

1. scalar proximal minimization gives (1.2)--(1.3);
2. exact profile minimization gives the exposed-boundary stability
   (2.2)--(2.3);
3. Plateau calibration gives the independent identity

   \[
   \varepsilon={1\over2}\sum_i
      \int_{\Lambda_i}|N_i-u_i|^2d\sigma_\mu.            \tag{6.1}
   \]

No implication from the first two to the third is valid.  The cap has the
strongest possible versions of 1 and 2 and still has `epsilon>0` for the
normal-defined trace (3.10).

A correct large-fill statement must therefore have the logical form

\[
 \boxed{\quad
 \text{parallel/concurrent/orthogonal-radial escape}
 \quad\text{or}\quad
 \text{a turning/Plateau estimate}.\quad}                \tag{6.2}
\]

The geometric escape must be an alternative, not a term added to the right
side of (0.2).  A convenient quantitative concurrence defect for formulating
that alternative is

\[
 \mathfrak C(\Sigma)
 =\inf_{c\in R^n}{1\over \ell^2\sigma_\mu(\Sigma)}
   \int_\Sigma|(I-NN^T)(x-c)|^2d\sigma_\mu,              \tag{6.3}
\]

with a separately specified physical scale `ell`; the cap has
`mathfrak C=0`.  The parallel defect may be written

\[
 \mathfrak P(\Sigma)
 =\inf_{|u|=1}{1\over\sigma_\mu(\Sigma)}
   \int_\Sigma(1-(N\cdot u)^2)d\sigma_\mu.               \tag{6.4}
\]

This note does not assert a gap theorem after imposing lower bounds on
`mathfrak C` and `mathfrak P`; proving such a theorem is precisely the new
geometric inverse still required.  What is proved is that the scalar
near-stationarity route cannot establish (8.4) before those escape branches
are removed.  A lemma restricted to the canonical heat-generated phase law
could still be true, but its proof would have to use a compatibility property
of that phase law which is absent from scalar quasiminimality and the Jacobi
form.  Calling that extra compatibility “near-stationarity” would be
circular.

## 7. Generality audit

1. **Hard support.**  All perimeters in Sections 3--4 are relative
   perimeters.  The support-contact flux in (3.13) and the contact term in
   (4.2) are retained.  The counterexample is therefore not produced by
   dropping the hard-wall terms.
2. **Lower-dimensional support.**  The example is full dimensional.  The
   arguments remain intrinsic on an affine support, but no reduction is
   needed here.
3. **Nonattainment.**  The cap attains the profile and the phase filling
   attains (3.12), so no limiting equality is used in the strict-excess
   argument.
4. **Regularity.**  The free surface and its phase trace are smooth; the
   support is smooth.  The two boundary pieces meet orthogonally, and the
   cap is Lipschitz at their contact.  There are no singular isoperimetric
   points in this model.
5. **Dimension tracking.**  The contradiction is `epsilon>0` versus a zero
   right side.  It does not hide a factor depending on `n`; it already works
   in dimension three and the ambient ball is exactly isotropic.
6. **Relation to the rank hypothesis.**  This is not a counterexample to
   the full three-way target (8.1)--(8.3), because its normal lines are
   exactly concurrent.  It is a counterexample to using (8.4) as an
   unconditional bridge from near-stationarity to filling rigidity.  Any
   rank-based proof must preserve the concurrent escape instead of trying
   to charge its arbitrary phase filling to the scalar deficit.

The load-bearing conclusion is therefore exact: a successful proof must
first recognize the global normal-line geometry.  Neither the scalar
Ekeland step nor the full free-boundary Jacobi inequality can do that by
controlling an auxiliary optimal filling flow.
