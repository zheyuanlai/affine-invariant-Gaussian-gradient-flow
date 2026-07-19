# Simple posterior eigenspaces: trace-free coherence and the remaining amplitude gate

## 0. Setting and verdict

Let \(Y=X+G\), where \(G\sim N(0,I_k)\) is independent of a log-concave
\(X\), and let

\[
 d\nu(y)=e^{-U(y)}dy,
 \qquad H(y)=D^2U(y)=I-A_y,
 \qquad A_y=\operatorname{Cov}(X\mid Y=y).
\]

Thus \(0\preceq H\preceq I\).  The Hilbert--Schmidt posterior saturation
theorem in `posterior_covariance_saturation_stability.md` says

\[
 \sum_{j=1}^k\left|(D_jH(y))u\right|^2
 \le16\,u^TH(y)u                                   \tag{0.1}
\]

for every deterministic or adaptively selected vector \(u=u(y)\); no
derivative of \(u\) occurs in (0.1).

Assume in this note that the smallest eigenvalue \(h(y)\) of \(H(y)\) is
simple and is separated uniformly from the rest of the spectrum:

\[
 \lambda_2(H(y))-h(y)\ge\kappa>0
 \qquad\text{for every }y.                          \tag{0.2}
\]

Let \(P(y)\) be its rank-one spectral projection and \(Q=I-P\).  The main
proved consequence is the trace-free eigenspace estimate

\[
 \boxed{
 \sum_{j=1}^k\|D_jP(y)\|_{HS}^2
 \le\frac{32}{\kappa^2}h(y).}                       \tag{0.3}
\]

For a low-energy exact field \(W=\nabla f\), this yields a complete
dimension-free ledger: \(W\) is close to the low line, the selected line
field varies slowly under the \(|W|^2d\nu\) law, and projecting \(W\) onto
that line preserves its norm, centering, and small deformation energy.
All constants are written below.

This closes the exact zero-defect case and proves that a simple-eigenspace
survivor cannot hide a factor of the dimension in the rotation tensor.  It
does **not** close the quantitative positive-defect case: the remaining
obstruction is amplitude-weighted global connectivity.  The probability
\(|W|^2d\nu/E|W|^2\) need not be log-concave and can place its line-field
transition in a region where \(|W|\) is tiny.  A dimension-free theorem
excluding that transition is not proved here.

## 1. Derivative of the simple spectral projection

Fix a point and choose a local unit eigenvector \(e\) with

\[
 He=he,
 \qquad P=e e^T.
\]

The sign of \(e\) is immaterial for \(P\).  Differentiating the eigenvalue
equation in the coordinate direction \(e_j\), projecting onto
\(e^\perp\), and using \(e\cdot D_je=0\) gives

\[
 D_je
 =-\left(Q(H-hI)Q\right)^{-1}Q(D_jH)e.             \tag{1.1}
\]

The inverse in (1.1) has operator norm at most \(\kappa^{-1}\) by
(0.2).  Applying (0.1) with \(u=e\) gives

\[
 \sum_j|D_je|^2
 \le\frac1{\kappa^2}
      \sum_j|(D_jH)e|^2
 \le\frac{16}{\kappa^2}h.                           \tag{1.2}
\]

Since

\[
 D_jP=(D_je)e^T+e(D_je)^T,
 \qquad \|D_jP\|_{HS}^2=2|D_je|^2,
\]

(1.2) proves (0.3).  This calculation is invariant under the local sign
choice and hence patches globally.  Analyticity of the Gaussian output and
the uniform gap make \(P\) analytic; the same calculation also holds in
weak \(C^1\) form whenever the spectral projection is defined by a
resolvent contour.

Two points are essential.  First, applying only the operator version of
posterior stability separately in each coordinate would give a factor
\(k\).  Equation (0.1) removes it exactly.  Second, the denominator is the
actual eigengap \(\kappa\); a lower bound on \(\lambda_2\) alone gives the
same result only on the region \(h\le\lambda_2-\kappa\).

## 2. Ledger for a low-energy field

Let \(W:\mathbb R^k\to\mathbb R^k\) be a \(C^1\) vector field and define

\[
 \mathcal D=E_\nu\|DW\|_{HS}^2,
 \qquad
 \mathcal K=E_\nu[W^THW],
 \qquad Z=PW.
\]

The spectral gap (0.2) implies the pointwise quadratic-form inequality

\[
 W^THW
 \ge h|W|^2+\kappa|QW|^2.                          \tag{2.1}
\]

Consequently

\[
 \boxed{E|W-Z|^2=E|QW|^2\le\frac{\mathcal K}{\kappa}.} \tag{2.2}
\]

For the derivative of the projection applied to \(W\), write locally
\(W=ae+r\), with \(r\perp e\).  If \(d_j=D_je\), then

\[
 (D_jP)W=a d_j+e(d_j\cdot r),
 \qquad |(D_jP)W|^2\le|W|^2|d_j|^2.
\]

Using (1.2) and the first term in (2.1),

\[
 \boxed{
 \sum_j|(D_jP)W|^2
 \le\frac{16}{\kappa^2}h|W|^2
 \le\frac{16}{\kappa^2}W^THW.}                    \tag{2.3}
\]

Since \(D_jZ=(D_jP)W+P D_jW\), equations (2.3) and
\(|a+b|^2\le2|a|^2+2|b|^2\) give

\[
 \boxed{
 E\|DZ\|_{HS}^2
 \le2\mathcal D+\frac{32}{\kappa^2}\mathcal K.}   \tag{2.4}
\]

If \(E|W|^2=1\), then

\[
 E|Z|^2\ge1-\frac{\mathcal K}{\kappa},             \tag{2.5}
\]

and

\[
 |EZ|\le|EW|+\sqrt{\mathcal K/\kappa}.             \tag{2.6}
\]

Thus a normalized almost-centered field with
\(\mathcal D+\mathcal K=O(\lambda)\) produces an almost-centered,
almost-normalized rank-one section \(Z=PW\) whose full derivative energy
is still \(O_\kappa(\lambda)\).  No ambient dimension occurs.

## 3. The amplitude-weighted line field

Assume \(E|W|^2=1\), and define the probability

\[
 d\sigma=|W|^2d\nu.
\]

On \(\{W\ne0\}\), let

\[
 N=\frac{WW^T}{|W|^2}
\]

be the rank-one projection onto the direction selected by \(W\).  The
distance between two rank-one projections is

\[
 \|N-P\|_{HS}^2
 =2\left|Q\frac{W}{|W|}\right|^2.
\]

Equation (2.1) therefore gives

\[
 \boxed{E_\sigma\|N-P\|_{HS}^2
 \le\frac{2}{\kappa}\mathcal K.}                   \tag{3.1}
\]

If \(n=W/|W|\), then

\[
 D_jN=(D_jn)n^T+n(D_jn)^T,
 \qquad
 |W|^2\sum_j\|D_jN\|_{HS}^2
 \le2\|DW\|_{HS}^2.
\]

Hence

\[
 \boxed{E_\sigma\sum_j\|D_jN\|_{HS}^2\le2\mathcal D.} \tag{3.2}
\]

Finally, multiplying (0.3) by \(|W|^2\), using
\(h|W|^2\le W^THW\), and averaging gives

\[
 \boxed{E_\sigma\sum_j\|D_jP\|_{HS}^2
 \le\frac{32}{\kappa^2}\mathcal K.}                \tag{3.3}
\]

Equations (3.1)--(3.3) are a dimension-free coherence theorem: under the
amplitude-biased law, the field direction and the simple posterior
eigendirection coincide up to \(O_\kappa(\sqrt{\mathcal K})\), and both
have small Dirichlet energy.  What they do not say is that either projection
is close to a single fixed projection.  Such a conclusion would require a
Poincare or connectivity estimate for \(\sigma\), not for \(\nu\).

## 4. Exact zero-defect rigidity

The local equality case is rigid.  Let \(\Omega\) be a connected open set,
let \(W=\nabla f\) be nonvanishing on \(\Omega\), and assume

\[
 HW=0,
 \qquad \ker H=\operatorname{span}\{W\}
 \quad\text{on }\Omega.                             \tag{4.1}
\]

Equation (0.1) with \(u=W(y)\) gives

\[
 (D_jH)W=0\quad\text{for every }j.                 \tag{4.2}
\]

Differentiating \(HW=0\) and using (4.2) yields

\[
 H(D_jW)=0.
\]

Thus every column of \(DW\) is parallel to \(W\).  Since \(DW=D^2f\) is
symmetric, there is a scalar function \(c\) such that

\[
 DW=c,WW^T.                                        \tag{4.3}
\]

It follows that the normalized direction \(W/|W|\) has zero derivative
on \(\Omega\).  Hence the direction is fixed, and \(W=a(t)u\) for a fixed
unit vector \(u\) and scalar coordinate \(t=u\cdot y\).  This proves the
exact version of the proposed simple-eigenspace rigidity without a
topological or dimension-dependent argument.

There is also a field-free global equality statement.  If \(h=0\)
everywhere and (0.2) holds, then (0.3) gives \(DP=0\), so
\(P=u u^T\) is fixed and

\[
 H\succeq\kappa(I-u u^T).
\]

The audited rank-one-defect theorem then gives

\[
 C_P(\nu)
 \le3\kappa^{-1}+96\operatorname{Var}_\nu(u\cdot Y). \tag{4.4}
\]

For the unscaled isotropic convolution \(\operatorname{Cov}(Y)=2I\), the
right side is at most \(3\kappa^{-1}+192\).  Thus an exact simple nullity
cannot support a vanishing gap.

## 5. Twice-Bochner decomposition: only the longitudinal scalar survives

The Hilbert--Schmidt estimate also removes the dimension loss from the
next Bochner level.  Let

\[
 L=\Delta-\nabla U\cdot\nabla,
 \qquad \mathcal A_1=-L+H
\]

be the Witten one-form operator, acting componentwise in the \(-L\) term.
For a smooth exact field \(W=\nabla f\), put

\[
 C=DW=D^2f,
 \qquad
 \mathcal T_W=(D H)[\,\cdot\,]W,
\]

where the \(j\)-th column of \(\mathcal T_W\) is \((D_jH)W\).
Both \(C\) and \(\mathcal T_W\) are symmetric.  Componentwise Bochner and
one integration by parts give the exact identity

\[
 \boxed{
 \begin{aligned}
 \|\mathcal A_1W\|_2^2
 &=E\|D C\|_{HS}^2
   +3E\operatorname{tr}(C H C)\\
 &\quad+2E\langle C,\mathcal T_W\rangle_{HS}
   +E|HW|^2.
 \end{aligned}}                                     \tag{5.1}
\]

Here \(\|DC\|_{HS}\) denotes the full third-order Euclidean tensor norm.
To verify the coefficient three, componentwise Bochner contributes one
copy of \(\operatorname{tr}(CHC)\), while

\[
 2\langle-LW,HW\rangle
 =2E\langle C,D(HW)\rangle_{HS}
\]

contributes two further copies and the cubic term.  Posterior saturation
gives pointwise

\[
 \boxed{\|\mathcal T_W\|_{HS}^2\le16W^THW.}         \tag{5.2}
\]

Thus even the fully contracted cubic tensor has no trace factor.

Under the simple eigengap (0.2), write

\[
 C=\alpha P+C_0,
 \qquad \alpha=e^TCe,
 \qquad \langle C_0,P\rangle_{HS}=0,
\]

and put \(t=e^T\mathcal T_We=D^3U[e,e,W]\).  Direct block algebra gives

\[
 \operatorname{tr}(CHC)
 \ge h\alpha^2+\frac\kappa2\|C_0\|_{HS}^2.          \tag{5.3}
\]

Splitting the cubic term into \(2\alpha t+2\langle
C_0,\mathcal T_W\rangle\), using Young's inequality on the second part,
and then (5.2), yields the dimension-free lower ledger

\[
 \boxed{
 \begin{aligned}
 \|\mathcal A_1W\|_2^2
 &\ge E\|DC\|_{HS}^2
   +3E[h\alpha^2]
   +\frac{3\kappa}{4}E\|C_0\|_{HS}^2\\
 &\quad+2E[\alpha t]+E|HW|^2
   -\frac{64}{3\kappa}E[W^THW].
 \end{aligned}}                                     \tag{5.4}
\]

For an attained normalized bottom mode,
\(\mathcal A_1W=\lambda W\) and

\[
 E\|DW\|_{HS}^2+E[W^THW]=\lambda E|W|^2.
\]

The spectral-window version follows from scalar/one-form intertwining and
closure.  Equations (5.1)--(5.4) show that every transverse component of
the deformation tensor is coercive with a dimension-free constant.  The
only unabsorbed signed quantity is

\[
 \boxed{E\left[(e^TDWe)D^3U[e,e,W]\right].}         \tag{5.5}
\]

When \(W=ae\) and the line is fixed, (5.5) is the one-dimensional
amplitude interaction \(E[a'a\,\partial_e h]\).  The coarse bound from
(5.2) is only

\[
 2|E\langle C,\mathcal T_W\rangle|
 \le8\sqrt{E\|DW\|_{HS}^2\,E[W^THW]},              \tag{5.6}
\]

which is of order \(\lambda\), whereas the left side for an eigenfield is
of order \(\lambda^2\).  It is therefore too weak to rule out a small gap.
This identifies the second-order obstruction as scalar longitudinal
amplitude, not high-dimensional eigenspace rotation.

## 6. Precise remaining gate

For a normalized bottom spectral-window field, one has

\[
 \mathcal D+\mathcal K=O(\lambda),
 \qquad |EW|^2=O(\lambda).
\]

Sections 2--3 then prove all of the following with constants depending only
on the fixed eigengap \(\kappa\):

1. the field lies in the simple low eigenspace in \(L^2(\nu)\);
2. its rank-one projection has no hidden trace-dimensional rotation cost;
3. the projected field remains normalized, almost centered, and slowly
   varying; and
4. exact zero defect forces a fixed line and falls under the rank-one-defect
   Poincare theorem.

The unresolved quantitative step is a weighted global-coherence principle
which prevents \(P\) and \(N\) from changing by order one only through a
set on which \(|W|\) is negligible.  A representative sufficient statement
would have to control the oscillation of a rank-one projection from the
three left sides of (3.1)--(3.3), uniformly for the non-log-concave tilted
law \(\sigma=|W|^2\nu\).  Ordinary Poincare for \(\nu\) would do this but is
circular; ordinary Poincare for \(\sigma\) is neither available nor implied
by log-concavity of \(\nu\).

Therefore the simple-eigenvalue branch has been reduced without dimension
loss to an amplitude/nodal connectivity problem.  Until that problem is
proved, (0.3) must not be promoted to a full simple-eigenspace KLS theorem.
