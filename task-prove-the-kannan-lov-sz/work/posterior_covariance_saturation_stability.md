# Posterior covariance saturation: a dimension-free third-moment stability theorem

## 0. Result and scope

Let \(\pi\) be a centered \(1\)-strongly log-concave probability measure
on \(\mathbb R^k\), let \(Z\sim\pi\), and put

\[
 A=\operatorname{Cov}(Z),\qquad
 \delta_\pi(u)=|u|^2-u^TAu\quad (u\in\mathbb R^k).
\]

Then \(0\preceq A\preceq I\).  Define the symmetric mixed-third-moment
matrix

\[
 M_\pi(u)=E[(u\cdot Z)ZZ^T].
\]

The dimension-free matrix estimate is

\[
 \boxed{\|M_\pi(u)\|_{HS}\le4\sqrt{\delta_\pi(u)}.} \tag{FM}
\]

Consequently, for every \(u,v,w\in\mathbb R^k\),

\[
 \boxed{
 \left|E[(u\cdot Z)(v\cdot Z)(w\cdot Z)]\right|
 \le 4\sqrt{\delta_\pi(u)}\,|v|\,|w|.}             \tag{TM}
\]

By symmetry, the right side can be replaced by the minimum of the three
corresponding expressions.  The constant is independent of the dimension.
The square-root dependence is the natural local scale: an odd perturbation
of a Gaussian can have third moment of first order while its covariance
deficit is only of second order.

For the unit-Gaussian posterior

\[
 \pi_y(dx)\propto
 \exp\left[-V(x)-\frac12|x-y|^2\right]dx,
 \qquad A_y=\operatorname{Cov}_{\pi_y}(X),
\]

(FM) gives both the Hilbert--Schmidt and operator derivative estimates

\[
 \boxed{
 \left(\sum_{j=1}^k\|(D A_y[e_j])u\|^2\right)^{1/2}
 =\|M_{\pi_y}(u)\|_{HS}
 \le4\sqrt{|u|^2-u^TA_yu},}                         \tag{0.1}
\]

and hence

\[
 \|(D A_y[v])u\|
 \le 4\sqrt{|u|^2-u^TA_yu}\,|v|.                  \tag{0.2}
\]

Thus posterior saturation controls all mixed third moments containing the
saturating direction, even when that direction is chosen adaptively from
\(y\).  This is a genuine local stability theorem.  It does not by itself
prove the global adaptive-coherence estimate needed for KLS; that remaining
step is isolated in Section 5.

## 1. Caffarelli deficit ledger

Let \(G\sim\gamma_k=N(0,I_k)\).  Caffarelli's contraction theorem, in the
form for a target whose potential minus \(|x|^2/2\) is convex, supplies a
Brenier map \(T=\nabla\varphi\) with

\[
 T_\#\gamma_k=\pi,
 \qquad 0\preceq J(g):=DT(g)\preceq I
 \quad\text{for a.e. }g.                            \tag{1.1}
\]

Translations preserve (1.1), so the centering assumption gives
\(ET(G)=0\).  Gaussian Poincare applied to \(u\cdot T\) gives

\[
 u^TAu=\operatorname{Var}(u\cdot T(G))
 \le E|J(G)u|^2.                                    \tag{1.2}
\]

Consequently

\[
 E\bigl(|u|^2-|Ju|^2\bigr)\le\delta_\pi(u).         \tag{1.3}
\]

For a symmetric contraction \(J\), spectral calculus gives

\[
 (I-J)^2\preceq I-J^2.
\]

Combining this pointwise inequality with (1.3) yields

\[
 \boxed{E|(I-J)u|^2\le\delta_\pi(u).}               \tag{1.4}
\]

The scalar function

\[
 r_u(g)=u\cdot(T(g)-g)
\]

has mean zero and weak gradient \((J-I)u\).  A second use of Gaussian
Poincare therefore gives

\[
 \boxed{E|u\cdot(T(G)-G)|^2\le\delta_\pi(u).}       \tag{1.5}
\]

Equations (1.4)--(1.5) are the two quantitative equality statements used
below.  They control a fixed saturated direction without summing over an
orthonormal basis.

## 2. Quadratic-form duality without a trace loss

The target \(\pi\) itself has Poincare constant at most one.  Therefore,
for every symmetric matrix \(B\),

\[
 \operatorname{Var}_\pi(Z^TBZ)
 \le E|2BZ|^2
 =4\operatorname{tr}(B^2A)
 \le4\|B\|_{HS}^2.                                 \tag{2.1}
\]

Let

\[
 r=u\cdot(T(G)-G),\qquad e=(I-J(G))u.
\]

Equations (1.4)--(1.5) say \(Er^2\le\delta_\pi(u)\) and
\(E|e|^2\le\delta_\pi(u)\).  Since \(Er=0\), (2.1) and
Hilbert--Schmidt duality give

\[
 \left\|E[rZZ^T]\right\|_{HS}
 =\sup_{\substack{B=B^T\\\|B\|_{HS}=1}}
   |E[r(Z^TBZ-EZ^TBZ)]|
 \le2\sqrt{\delta_\pi(u)}.                         \tag{2.2}
\]

Likewise, for an arbitrary matrix \(B\) with \(\|B\|_{HS}=1\),

\[
 |E[e^TBZ]|
 \le(E|e|^2)^{1/2}(E|BZ|^2)^{1/2}
 \le\sqrt{\delta_\pi(u)},
\]

because \(A\preceq I\).  Hence

\[
 \boxed{\|E[eZ^T]\|_{HS}\le\sqrt{\delta_\pi(u)}.} \tag{2.3}
\]

This quadratic duality is the step that prevents a factor \(\sqrt k\):
one never bounds the entries of the third-moment matrix separately.

## 3. Proof of the mixed-third-moment theorem

Write \(Z=T(G)\).  Entrywise Gaussian integration by parts gives the
matrix identity

\[
 E[(u\cdot G)ZZ^T]
 =E[(Ju)Z^T+Z(Ju)^T].                               \tag{3.1}
\]

Since \(Ju=u-e\) and \(EZ=0\), the constant terms vanish.  Also
\(u\cdot Z=u\cdot G+r\).  Consequently

\[
 M_\pi(u)
 =E[rZZ^T]-E[eZ^T+Ze^T].                            \tag{3.2}
\]

Equations (2.2)--(2.3) and the triangle inequality now give

\[
 \|M_\pi(u)\|_{HS}
 \le2\sqrt{\delta_\pi(u)}
    +2\sqrt{\delta_\pi(u)}
 =4\sqrt{\delta_\pi(u)},                           \tag{3.3}
\]

which proves (FM).  Since

\[
 E[(u\cdot Z)(v\cdot Z)(w\cdot Z)]=v^TM_\pi(u)w,
\]

(TM) follows from \(\|M\|_{op}\le\|M\|_{HS}\).  The argument first
applies to a smooth positive target.  For an extended convex potential or
a hard convex support, use the standard strongly-convex approximation in
Caffarelli's theorem; (1.1) persists, Gaussian concentration gives uniform
moments, and weak convergence plus uniform integrability passes (FM) to the
limit.  Equivalently, the Lipschitz Brenier map formulation and Gaussian
Sobolev integration by parts already apply directly.  Proper affine
supports are handled intrinsically in their supporting space.

## 4. Posterior covariance and output-Hessian corollaries

For the Gaussian posterior \(\pi_y\), write

\[
 m_y=E_{\pi_y}X,qquad Z_y=X-m_y.
\]

The quadratic likelihood gives uniform Gaussian tails for each posterior,
so differentiation under the integral is valid for every polynomial
moment.  The exponential-family derivative identity gives

\[
 D A_y[v]
 =E_{\pi_y}\left[Z_yZ_y^T(v\cdot Z_y)\right].       \tag{4.1}
\]

The matrix of the linear map \(v\mapsto(D A_y[v])u\) is exactly
\(M_{\pi_y}(u)\).  Thus (FM) proves (0.1)--(0.2).  In bilinear form,

\[
 |w^T(D A_y[v])u|
 \le4\sqrt{|u|^2-u^TA_yu}\,|v|\,|w|.               \tag{4.2}
\]

Let \(U=-\log p_{X+G}\).  The posterior identity

\[
 H(y):=D^2U(y)=I-A_y
\]

turns (4.2) into the self-concordance-type bound

\[
 \boxed{
 |D^3U(y)[u,v,w]|
 \le4\sqrt{u^TH(y)u}\,|v|\,|w|.}                  \tag{4.3}
\]

In particular,

\[
 \left(\sum_{j=1}^k\|(D H(y)[e_j])u\|^2\right)^{1/2}
 \le4\sqrt{u^TH(y)u}.                              \tag{4.4}
\]

In particular,

\[
 \|(D H(y)[v])u\|
 \le4\sqrt{u^TH(y)u}\,|v|.                         \tag{4.4a}
\]

For each fixed \(u\), the function

\[
 y\longmapsto\sqrt{u^TH(y)u}
\]

is \(2|u|\)-Lipschitz.  Indeed, at positive values this follows by
differentiating and using (4.3) with the two slots \(u,u\); approximation
by \(\sqrt{u^THu+\varepsilon}\) gives the global statement.  For the
isotropic rescaling \(S=(X+G)/\sqrt2\), whose potential is
\(\widetilde U(s)=U(\sqrt2s)+\mathrm{const}\), the corresponding bound is

\[
 |D^3\widetilde U(s)[u,v,w]|
 \le8\sqrt{u^TD^2\widetilde U(s)u}\,|v|\,|w|.     \tag{4.5}
\]

## 5. What remains for adaptive saturation

Let \(W=\nabla h\) be a normalized bottom spectral-window field for the
unscaled output.  The Bochner ledger has the form

\[
 E\|DW\|_{HS}^2+E[W^THW]=O(\lambda),
 \qquad E|W|^2=1,qquad |EW|^2=O(\lambda).           \tag{5.1}
\]

Equation (4.4) now applies pointwise with the adaptive choice \(u=W(y)\):

\[
 \boxed{
 \left(\sum_{j=1}^k\|(D H(y)[e_j])W(y)\|^2\right)^{1/2}
 \le4\sqrt{W(y)^TH(y)W(y)}.}                       \tag{5.2}
\]

There is also a directly contracted derivative consequence.  Since

\[
 D(HW)[v]=(D H[v])W+H(DW[v])
\]

and \(\|H\|_{op}\le1\), (5.2) gives pointwise

\[
 \boxed{
 \|D(HW)\|_{HS}
 \le4\sqrt{W^THW}+\|DW\|_{HS}.}                    \tag{5.3}
\]

In particular, if
\(\mathcal D=E\|DW\|_{HS}^2\) and
\(\mathcal K=E[W^THW]\), then

\[
 E\|D(HW)\|_{HS}^2\le2\mathcal D+32\mathcal K,
 \qquad E|HW|^2\le\mathcal K.                       \tag{5.4}
\]

Thus a low-energy survivor makes the curvature-weighted field \(HW\)
small in the full first-order Sobolev norm, again without a trace factor.

Thus the missing field cannot exploit an arbitrary rotating kernel: the
variation of the posterior Hessian in its selected direction is controlled
by the square root of the same curvature deficit appearing in (5.1).
In particular, the derivative estimate itself has no trace loss.  This
removes the first obstruction to differentiating an adaptively selected
near-null direction.  It still does not by itself imply

\[
 E|W|^2\le C\{E\|DW\|_{HS}^2+E[W^THW]\};
\]

the latter is exactly the KLS-equivalent one-form coercivity audited in
`upper_curvature_oneform_coercivity_audit.md`.  A valid continuation must
combine (5.2) with the symmetry and exactness of \(DW\), or with the
Gaussian-fiber conditional
Poincare estimate

\[
 E|W(X+G)-E[W(X+G)\mid X]|^2\le E\|DW(X+G)\|_{HS}^2.
\]

No global coherence claim is made in this note.
