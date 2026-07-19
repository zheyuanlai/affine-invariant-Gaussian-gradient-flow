# Checkpoint 12: defect-subspace curvature and the singular survivor

## 1. Candidate-proof status

There is still no complete proof, and no dimension-free conclusion is
claimed at this checkpoint.  The principal positive development is a new
defect-subspace Poincare theorem which has passed an independent clean-room
reproof.  Three attempts to feed it into known localization machinery expose
one sharply isolated obstruction: a posterior may retain a single weak
direction selected from the same future noise that selected its covariance.
The displacement and signed-distance routes independently reach the same
geometric obstruction as a high-rank family of almost-translations or long
normal cells whose global incidence has not yet been controlled.

All constants below are numerical and independent of the dimension.

## 2. Audited curvature outside a weak subspace

Let \(\mathbb R^n=F\oplus E\), let \(\mu=e^{-V}dx\) be log-concave, and
assume distributionally that

\[
                         D^2V\succeq\kappa P_E .          \tag{2.1}
\]

Writing \(\bar\mu=(P_F)_\#\mu\), the new theorem is

\[
 \boxed{
 C_P(\mu)\le \kappa^{-1}
      +2\bigl(C_P(\bar\mu)+\kappa^{-1}\bigr)
 \le3\bigl(C_P(\bar\mu)+\kappa^{-1}\bigr).}             \tag{2.2}
\]

For \(F=\mathbb Ru\), the one-dimensional log-concave estimate gives

\[
 C_P(\mu)\le96\left(\kappa^{-1}
             +\operatorname {Var}_\mu\langle X,u\rangle\right). \tag{2.3}
\]

The proof is not a formal appeal to Brascamp--Lieb.  Conditional transport
velocities \(\mathcal W_z\) give the reinforced Prekopa inequality

\[
 D^2U(z)\succeq\kappa G(z),\qquad
 G(z)=\mathbb E_z\mathcal W_z^*\mathcal W_z,             \tag{2.4}
\]

and matrix Cauchy--Schwarz controls the derivative of the conditional
expectation with the weight \((I+G)^{-1}\).  The load-bearing marginal
estimate is the proved exact-one-form inequality

\[
 \operatorname {Var}_\eta h
 \le2\int\langle(D^2U+C_P(\eta)^{-1}I)^{-1}\nabla h,
                         \nabla h\rangle d\eta.          \tag{2.5}
\]

It follows from the scalar spectral gap on the closure of gradients, the
one-form intertwining identity, and the Helffer--Sjostrand variational
formula.  Gaussian convolution preserves (2.1) with
\(\kappa_\varepsilon^{-1}=\kappa^{-1}+\varepsilon\), and tensorization
gives \(C_P(\bar\mu*\gamma_\varepsilon)\le
C_P(\bar\mu)+\varepsilon\).  Thus (2.2) holds for distributional convex
potentials and intrinsic affine supports.

The independent audit in `audit_rank_defect_lichnerowicz.md` rederived the
reinforced Prekopa inequality, the operator-domain and reducing-subspace
steps in (2.5), the matrix comparison, the nonsmooth convolution passage,
and the affine-support cases.  It found no hidden dimension factor or
unproved Poincare input.

## 3. Ordinary localization does not amplify the theorem

Let

\[
 \mathcal K_d=\sup\{C_P(\nu):\nu\text{ isotropic log-concave},
                   \dim\operatorname {aff}\nu\le d\}.
\]

Inserting (2.2) into covariance-whitened stochastic localization and using
only the density-martingale and trace identities gives the exact raw
recurrence

\[
 \mathcal K_n\le C\left({n\over m}+n\mathcal K_m\right). \tag{3.1}
\]

A spectrum-adaptive water-filling stop gives

\[
 \mathcal K_n\le C\left[\left({n\over m}\right)^2
                         +n\mathcal K_m\right].          \tag{3.2}
\]

Both are weaker than the elementary trace estimate.  Even granting the
counterfactual missing assertion that the weak marginal has covariance
operator norm \(O(1)\), the resulting recurrence

\[
                         \mathcal K_n\le C(n/m+\mathcal K_m) \tag{3.3}
\]

cannot improve \(\mathcal K_m\lesssim\log m\): if \(m\le n/\log n\),
the curvature term is at least logarithmic, while otherwise
\(\log m\asymp\log n\).  Thus the theorem alone does not even meet M1
through the standard trace bootstrap.  A successful localization must
generate constant curvature outside \(n^{o(1)}\) weak directions while
also controlling the covariance of that weak marginal, not merely its
dimension.

## 4. Bounded mass-preserving localization and its precise gate

For a fixed balanced set \(S\), write

\[
 A_t=\operatorname {Cov}_{\mu_t}X,\qquad
 b_t=\operatorname {Cov}_{\mu_t}(1_S,X),\qquad
 u_t=b_t/|b_t|.
\]

On every zero-signal-free interval the bounded controller

\[
 C_t=P_t=I-u_tu_t^T                                      \tag{4.1}
\]

preserves \(\mu_t(S)\) pathwise.  Its accumulated quadratic tilt is

\[
 B_t=tI-Q_t,\qquad Q_t=\int_0^tu_su_s^Tds,\qquad
 \lambda_2(B_t)\ge t/2.                                 \tag{4.2}
\]

At \(t=1\), (2.3), Buser--Ledoux, and perimeter averaging reduce the
desired estimate to

\[
 \boxed{\mathbb E\left(1+
    \operatorname {Var}_{\mu_1}\langle X,w_1\rangle
                  \right)^{-1/2}\ge c,}                 \tag{4.3}
\]

where \(w_1\) is the weakest eigenvector of \(B_1\).  The mixture identity
\(\mathbb EA_1\preceq I\) applies only to deterministic or past-measurable
directions; it gives merely \(\mathbb Ew_1^TA_1w_1\le n\) for this
future-selected direction.  An explicit adaptive endpoint ensemble
achieves equality at scale \(n\), so covariance bookkeeping alone cannot
prove (4.3).

The hard convention at \(b=0\) is not well posed: for a one-dimensional
Gaussian and a centred interval it contradicts the occupation-density
formula.  The continuous soft controller

\[
 C_\varepsilon(b)=I-{bb^T\over |b|^2+\varepsilon}       \tag{4.4}
\]

removes this existence defect.  It obeys

\[
 |C_\varepsilon b|^2
 ={\varepsilon^2|b|^2\over(|b|^2+\varepsilon)^2}
 \le {\varepsilon\over4},                               \tag{4.5}
\]

so the balanced set mass remains in a fixed central interval with
constant probability for small fixed \(\varepsilon\).  Explicitly,
\[
 \mathbb P\{|M_t-\tfrac12|>\tfrac14\}\le4\varepsilon t,
 \qquad
 \mathbb P\{\sup_{s\le t}|M_s-\tfrac12|>\tfrac14\}
       \le16\varepsilon t.                              \tag{4.6}
\]
Moreover
\(C_\varepsilon^2=I-\alpha uu^T\), \(0\le\alpha\le1\), so all but one
eigenvalue of \(B_1\) remain at least \(1/2\).  This repairs the
zero-signal coefficient defect and mass survival on the audited stopped
flow, but it does not yet repair (4.3).

The generic-identity obstruction is now realized by an actual projected
filter, not only by a static endpoint ensemble.  Let \(\mu_n\) be uniform
on \(\{\pm\sqrt n\,e_i:1\le i\le n\}\), and let \(E\) be the positive
atoms.  The exact feedback has \(|b_t|\ge1/2\), hence never meets the
zero-signal singularity.  With failure probabilities at most
\[
 992\sqrt{\log n/n}\quad\text{and}\quad16\log n/n,
\]
the occupation line aligns with the terminal signal and its posterior
variance is at least \(3n/16\).  Consequently the left side of (4.3) is
\(O(\sqrt{\log n/n})\).  This law is isotropic but not log-concave, so it
is a rigorous mechanism no-go, not a KLS counterexample.

A static distinction supplied by log-concavity is also proved.  For every
one-dimensional log-concave \(Z\) and every binary label \(Y\),
\[
             \mathbb E\operatorname {Var}(Z\mid Y)
                  \ge c_0\operatorname {Var}Z.           \tag{4.7}
\]
Thus complementary posterior halves cannot collapse to two atoms.
Critical linear tilts of the isotropic \(\ell_1\) ball and of an
exponential marginal nevertheless produce genuine long continuous needles.
The live log-concavity-only target is the orientation-convexification tail
bound
\[
 \mathbb P\{\operatorname {Var}_{\mu_1}
        \langle X,v_1\rangle\ge R,\ \rho_1\le r\}
       \le C(r+R^{-1/2}),                                \tag{4.8}
\]
where \(\rho_1\) is the accumulated curvature in the terminal occupation
line.  It is sufficient for (4.3), but is not proved.

The remaining deterministic dichotomy is sharp.  If
\(1-\lambda_{\max}(Q_1)\) is bounded below, the terminal measure is fully
strongly log-concave.  Otherwise \(u_t\) is \(L^2(dt)\)-close to one final
line.  The unresolved step is to replace that future-selected line by a
past-measurable anchor without losing a trace factor.  The median maximum
of independent shifted exponentials has angular quadratic-variation rate
\(\Theta((\log n)^2)\) and is the principal stress test for every such
anchor argument.

## 5. Displacement midpoint: bounded strain, unbounded translation rank

For \(\mu_+=2\,1_E\mu\), \(\mu_-=2\,1_{E^c}\mu\), and the Brenier map
\(T:\mu_+\to\mu_-\), the midpoint law
\(\nu_{1/2}=((I+T)/2)_\#\mu_+\) satisfies

\[
                         \nu_{1/2}\le2\mu,              \tag{5.1}
\]

and the exact entropy identity

\[
 \log2-\operatorname {Ent}_\mu(\nu_{1/2})
 =\mathbb E\left[
 {V(X)+V(TX)\over2}-V\!\left({X+TX\over2}\right)
 +\log\det {I+DT\over2(DT)^{1/2}}\right].              \tag{5.2}
\]

The covariance identity is

\[
 \operatorname {Cov}(\nu_{1/2})=I-\tfrac14K,\qquad
 K=\mathbb E[(T-X)\otimes(T-X)]\preceq4I.              \tag{5.3}
\]

The matrix deficit gives the dimension-free Cayley strain budget

\[
 \mathbb E\left\|{DT-I\over DT+I}\right\|_{HS}^2
 \le2\log2.                                             \tag{5.4}
\]

It does not see singular Hessian mass or piecewise translations.  In the
exact zero-strain branch, \(D_aT=I\) on the source does not imply that
\(T-I\) is the gradient of a convex function; a one-dimensional two-cell
example has a downward displacement jump across the target gap.  Under the
additional convex polyhedral displacement hypothesis, translation rank
\(r\) gives \(W_2^2\le4r\), and disjoint midpoint prisms give an exact
weighted wall budget.  The missing case is a high-rank family of
small-strain translations connected through physical gaps or singular
walls.  Finally,

One high-rank corridor family is now closed.  For an isotropic box with a
coordinatewise monotone Moreau label, suppose the label is locally constant
off coordinate gaps of relative widths \(\delta_j\), and the midpoint
plateau has half the box volume.  Then

\[
 \prod_j(1-\delta_j)=\tfrac12,
 \qquad \operatorname {Var}_R b_j\le\delta_j^2,
 \qquad \sum_j\delta_j\le\log2,                         \tag{5.5}
\]

so \(\mathbb E_R|b-\mathbb Eb|^2\le(\log2)^2\), apart
from the already bounded mean-displacement term.  Thus independent
orthogonal corridors lose volume multiplicatively while their translation
energy adds quadratically.  The tensor product that would make the energy
add has source mass \(2^{-n}\); parity repair reverses some one-dimensional
branches and violates Brenier monotonicity.  The unproved extension is a
mixed-volume analogue of (5.5) for a noncommuting connected power-cell
complex.

The mixed-volume extension is now complete whenever the local corridors
share one global convex core.  If
\[
 C_0\subset\cdots\subset C_m\subset K,\qquad |C_0|=|K|/2,
 \qquad C_{j-1}+[0,v_j]\subset C_j,
\]
then isotropy alone gives the projection-thickness bound
\[
 {|A|\over|P_{u^\perp}A|}\le\sqrt{12/\alpha}
 \quad\text{for }A\subset K,\ |A|=\alpha|K|.
\]
The exact extrusion identity
\(|C+[0,v]|=|C|+|v||P_{v^\perp}C|\) therefore telescopes to
\[
                         \sum_j|v_j|\le\sqrt{24}.        \tag{5.6}
\]
The directions need not commute.  If the labels lie in the resulting
zonotope, then
\(\operatorname {tr}\operatorname {Cov}B\le12\) and
\(\mathbb E|B|^2\le13\).  Symmetric simplex and crosspolytope fans also
close by exact Dirichlet-union and homothetic-core calculations.  The
remaining premise is local-to-global: arbitrary power walls provide only
local facet prisms, not a common half-volume convex core or a laminar
sequence of global extrusions.

Finally,

\[
 \mathcal W_2(\mu)\le4\sqrt{C_P(\mu)},\qquad
 D_1(\mu)\le\tfrac12\mathcal W_2(\mu),                  \tag{5.7}
\]

so a universal bound on this displacement is itself KLS-equivalent.

## 6. Completed singular ray curvature and global incidence

For a balanced normal cell with conditional log-concave density \(q_y\),
tail transport \(W_y\), and \(\kappa_y=D^2(-\log q_y)\), the correct
completed curvature measure on the compactified ray is

\[
 d\mathcal C_y=W_y\,d\kappa_y+e_y^-\delta_{a_y}
                                      +e_y^+\delta_{b_y},
 \qquad \mathcal C_y([a_y,b_y])=2q_y(0).                \tag{6.1}
\]

After integration, \(\mathcal C=2P_\mu(E)\).  At a regular medial wall the
endpoint contribution is

\[
 e^{-V}|u^+-u^-|\,d\mathcal H^{n-1},                    \tag{6.2}
\]

which is linear BV turning; at a hard support it is
\(e^{-V}|u\cdot n_K|d\mathcal H^{n-1}\).  There is no canonical quadratic
singular Hessian measure: the transport weight vanishes at the endpoint,
and logarithmic smoothing of a Gauss jump drives the weighted quadratic
energy to zero.

For a fixed-mass core of cells of scale \(D\), the audited consequences are

\[
 \mathcal C_{\mathcal G}\le {C\over D},\qquad
 \operatorname {Cov}_{op}(Y\mid\mathcal G)\le C,
 \qquad \operatorname {tr}Q_{\mathcal G}=1,\qquad
 \|Q_{\mathcal G}\|_{HS}\le {C\over D}.               \tag{6.3}
\]

Hence \(\operatorname {rank}Q_{\mathcal G}\gtrsim D^2\).  An abstract
orthogonal star of long coordinate segments satisfies every numerical
constraint with \(D\to\infty\), but is neither log-concave nor the normal
congruence of a codimension-one interface.

There is also a new dimension-free local shape estimate.  If
\(s_y=1/q_y(0)\), the balanced cell contains the normal interval
\((-s_y/4,s_y/4)\), and every principal curvature satisfies
\(|\kappa_j(y)|\le4/s_y\).  On \(|t|\le s_y/8\), the tail transport is at
least \(1/4\) and \(1/2\le1+t\kappa_j\le3/2\).  Therefore its regular
normal-Jacobian charge obeys

\[
 \mathcal C_y^{\rm shape}\ge {s_y\over36}\|S_y\|_{HS}^2,
 \qquad
 \int_{\mathcal G}\|S_y\|_{HS}^2d\eta(y)
       \le {36A\over aD^2}.                              \tag{6.4}
\]

The associated reach inequality
\[
 |\langle N(y),z-y\rangle|
       \le {|z-y|^2\over2L_y},\qquad L_y=s_y/4,           \tag{6.5}
\]
controls only the normal component of each basepoint chord; direct
summation loses the transverse trace.  Stable rank also does not produce
\(D^2\) direction caps of mass \(D^{-2}\): the realized radial exponential
has a uniform spherical Gauss law and exponentially small narrow caps.
Prékopa--Leindler fills cross-pair midpoints, but their uncontrolled overlap
and transverse trace prevent the naive convexification sum.

A new exact harmless branch has now been proved.  If a closed \(C^2\)
interface has untruncated regular normal cells for every
\(t\in\mathbb R\), positivity of \(\det(I+tS_y)\) forces \(S_y=0\);
completeness and absence of cuts make every component a hyperplane.
For an isotropic half-mass hyperplane of offset \(a\), Cantelli gives
\(|a|\le1\), and

\[
                         \mathbb E|\langle X,u\rangle-a|\le\sqrt2. \tag{6.6}
\]

High rank itself cannot be excluded.  For the uniform isotropic cube and
the parity set, the coordinate-hyperplane fan has balanced cells,
\(Q=I/n\),

\[
 J={\sqrt3\over n+1},\qquad P={n\over2\sqrt3},
 \qquad 2PJ={n\over n+1}.                               \tag{6.7}
\]

All rigidity is carried by the transition incidence, at exactly the
linear-BV scale.  The live target is therefore the positive-incidence
inequality

\[
 \mathcal C_{\mathcal G}
      +\lambda_{\max}(Q_{\mathcal G})^{1/2}\ge c(A,\beta) \tag{6.8}
\]

under the realized log-concave normal-congruence hypotheses and the core
conditions in (6.3).  It holds in the endpoint-free branch and passes the
parity fan, radial exponential, and cube-halfspace tests; it is false for
the unrealized abstract star.  Inequality (6.8) is not yet proved and is
currently another exact formulation of the global compatibility gap.

## 7. Registry after checkpoint 12

### Active

1. **Soft bounded localization.**  Prove or refute a dimension-free
   survivor inverse moment after using (4.4), with an adapted-anchor
   argument that explicitly handles future selection.
2. **Normal-congruence incidence.**  Quantify the fact that stable-rank
   packets of long, almost-orthogonal cells must be convexified by the
   ambient log-concave density, producing medial/support incidence or
   basepoint spread.
3. **Zero-gap connectivity.**  Convert the exact wall-prism and firm
   nonexpansiveness identities into an \(L^2\) translation bound, retaining
   physical target gaps and the nonconvex displacement-potential branch.
4. **Defect theorem applications outside ordinary localization.**  Seek a
   construction whose weak marginal has controlled covariance without
   invoking a covariance-process estimate equivalent to KLS.

### Blocked

1. Ordinary stochastic-localization trace bootstraps, including
   water-filling, after the explicit recurrences (3.1)--(3.3).
2. Hard projected localization through zero signal, by the occupation-time
   contradiction.
3. Any claim that \(\mathbb EA_t\preceq I\) controls an adaptively selected
   terminal direction.
4. Quadratic singular curvature at medial or hard-support endpoints.
5. Strain-only geometric rigidity, which loses piecewise translations and
   target gaps.
6. Rank-only ray rigidity, refuted by the radial exponential law and the
   realized cube parity fan.
7. Same-interface use of Milman's equivalence to infer \(P\lesssim1/D\);
   the extremizing set and extremizing Lipschitz function need not coincide.

No item in the active list is being treated as a lemma until it has a
dimension-free proof with the approximation, affine-support, and test-model
audits requested in the task.
