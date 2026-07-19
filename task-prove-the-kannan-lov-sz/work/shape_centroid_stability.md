# Shape-changing slices: centroid stability audit

This note records what can be proved about the conditional centroid
\\(m(s)=\\mathbb E[Z\\mid S=s]\\) when the conditional slices genuinely change
shape.  The goal is to control

\\[
  \\int \\rho(s)\\,\\tau(s)^2 |m'(s)|^2\\,ds                         \\tag{0.1}
\\]

using only joint convexity, isotropy, and the deformation budget from
weighted_fisher_prekopa.md, without inserting a conditional Poincare
inequality.  The general statement remains open.  Two nontrivial classes
are settled below: (i) conditional Gaussian slices with a changing,
possibly noncommuting covariance, and (ii) planar interval (wedge/cone)
slices.  The calculations also identify the obstruction for arbitrary
shape-changing fibers.

Throughout, \\(S\\) is centered and variance one, \\(\\rho=e^{-\\Phi}\\) is
log-concave, and \\(\\tau\\) is its canonical Stein kernel.  We use the
identities from the weighted-Fisher note:

\\[
 m'=-\\mathbb E_sF_s,\\qquad F_s=-m'(s)+F_s^0,\\qquad
 \\mathbb E_sF_s^0=0,                                      \\tag{0.2}
\\]

\\[
 m''=-\\mathbb E_s[(Z-m)(\\mathcal Q_s+\\|D_zF_s\\|_{\\rm HS}^2)],\\qquad
 \\int \\rho\\tau^2 C_s\\le 1,                                \\tag{0.3}
\\]

where \\(C_s=\\mathbb E_s\\|D_zF_s\\|_{\\rm HS}^2\\).  Isotropy in the mixed
coordinates gives
\\(\\mathbb E[\\tau(S)m'(S)]=0\\).

The bounds below use only elementary one-dimensional log-concave facts.  We
state the numerical constants generously; optimizing them is irrelevant.

## 1. One-dimensional estimates used repeatedly

We record three elementary estimates for a centered variance-one
log-concave law.  There are universal constants (M_3,M_\\tau,L_\\mathrm{mode})
such that

\\[
 |s_0|\\le L_\\mathrm{mode}\\quad (s_0\\text{ any mode}),\\qquad
 \\mathbb E|S|^3\\le M_3,\\qquad \\mathbb E\\tau^2\\le M_\\tau.       \\tag{1.1}
\\]

The first two follow from the standard one-dimensional exponential tail
bound (variance one); the last is the usual one-dimensional Stein-kernel
moment estimate.  We may take, for definiteness,
\\(L_\\mathrm{mode}=4,\\ M_3=10^4,\\ M_\\tau=400\\).

A useful pointwise estimate does not require any smoothness beyond a
one-sided derivative.  At a point where (a=|\\Phi'(s)|>0), convexity gives,
in the direction in which \\(\\Phi\\) increases,
\\( \\rho(s+u)\\le \\rho(s)e^{-au}\\).  Using the corresponding right- or
left-tail formula for \\(\\tau\\), and \\(|s+u|\\le |s|+u\\), we obtain

\\[
 \\boxed{\\ \\tau(s)|\\Phi'(s)|\\le 1+|s|.\\ }                    \\tag{1.2}
\\]

For example, on a right tail,
\\[
 \\tau(s)|\\Phi'(s)|\\le
 a\\,\\rho(s)^{-1}\\int_0^\\infty |s+u|\\rho(s+u)\\,du
 \\le |s|+1;
\\]
the left tail is identical.  At a mode the left side is zero; at a kink
use one-sided limits.  Consequently

\\[
 \\mathbb E\\!\\left[\\tau^2(1+|S|)(\\Phi')^2\\right]
 \\le \\mathbb E(1+|S|)^3\\le 8(1+M_3).                     \\tag{1.3}
\\]

We will also use the following growth fact.  If (h\\ge0) is concave on the
support interval and \\(\\mathbb Eh(S)\\le\\Lambda\\), then

\\[
 h(s)\\le L_\\rho\\Lambda(1+|s|),                            \\tag{1.4}
\\]

with a universal (L_\\rho).  To verify it, use the universal isotropic core
\\([-b,b]\\) on which \\(\\rho\\ge c_\\rho\\).  Concavity and nonnegativity imply
\\(h(t)\\ge h(0)/2\\) for \\(|t|\\le b/2\\), hence
\\(h(0)\\le 2\\Lambda/(bc_\\rho)\\).  For (s>b), the secant slope between
\\(-b\\) and (0) bounds the positive slope from (0) to (s), giving
\\(h(s)\\le h(0)(1+s/b)\\); the negative side is symmetric.  One may take
\\(L_\\rho=10^4\\) with the rough core constants (b=1/20\\) and
\\(c_\\rho=10^{-3}\\).  No dimension enters (1.4).

## 2. A Gaussian shape-changing theorem

Consider a smooth jointly log-concave density on an interval (J) times
\\(\\mathbb R^d\\) of the form

\\[
 p(s,z)=\\exp\\!\\left[-W(s)-\\frac12(z-m(s))^TQ(s)(z-m(s))\\right]/Z,
 \\qquad Q(s)\\succ0,\\quad R(s)=Q(s)^{-1}.                 \\tag{2.1}
\\]

The conditional law is \\(N(m(s),R(s))\\).  Assume

\\[
 \\mathbb E[S\\,m(S)]=0,\\qquad \\mathbb E R(S)\\preceq\\Lambda I_d. \\tag{2.2}
\\]

The first condition is the mixed-isotropy cancellation; the second holds,
for example, after unit transverse Gaussian smoothing of a globally
isotropic law, with \\(\\Lambda\\) an absolute constant.  We prove

\\[
 \\boxed{\\quad \\mathbb E[\\tau(S)^2|m'(S)|^2]\\le C_\\Lambda,\\quad} \\tag{2.3}
\\]

where (C_\\Lambda) is universal for fixed \\(\\Lambda\\), independent of (d\\),
of condition numbers, and of rotations of the eigenspaces of \\(R(s)\\).
Thus this class satisfies the requested weaker form

\\[
 \\mathbb E[\\tau^2|m'|^2]\\le C_\\Lambda\\left(1+\\int\\rho\\tau^2C_s\\right),\\tag{2.4}
\\]

since the deformation budget is nonnegative (and in the weighted-Fisher
setting is at most one).

### 2.1 Exact convexity calculation

Write (y=z-m(s)), (U=W+\\frac12y^TQy).  Its block derivatives are

\\[
\\begin{aligned}
 U_{zz}&=Q, & U_{sz}&=Q'y-Qm',\\
 U_{ss}&=W''+\\tfrac12y^TQ''y-2m'^TQ'y-m''^TQy+m'^TQm'.
\\end{aligned}                                                   \\tag{2.5}
\\]

Taking the Schur complement of (Q) gives the exact condition

\\[
 W''-m''^TQy+\\frac12y^TQ(-R'')Qy\\ge0\\quad\\text{for every }y. \\tag{2.6}
\\]

In particular (K:=-R''\\succeq0).  With (x=Qy), (2.6) is

\\[
 W''-m''\\!\\cdot x+\\frac12x^TKx\\ge0\\quad(x\\in\\mathbb R^d). \\tag{2.7}
\\]

Therefore (m''\\in\\operatorname{Ran}K) and

\\[
 (m'')^TK^\\dagger m''\\le2W''.                           \\tag{2.8}
\\]

For a unit vector (u), Cauchy--Schwarz in the (K)-metric and
\\(K=R^{1/2}(R^{-1/2}KR^{-1/2})R^{1/2}\\) imply

\\[
 |u\\!\\cdot m''|^2\\le2W''\\,u^TKu
 \\le2W''\\,\\lambda_{\\max}(R)\\operatorname{tr}(R^{-1}K). \\tag{2.9}
\\]

The marginal potential is

\\[
 \\Phi=W-\\tfrac12\\log\\det R+\\mathrm{const},\\qquad
 \\Phi''=W''+\\frac12\\operatorname{tr}(R^{-1}K)
       +\\frac12\\|R^{-1/2}R'R^{-1/2}\\|_{\\rm HS}^2.          \\tag{2.10}
\\]

Combining (2.9)--(2.10), and using (W''\\le\\Phi''), yields the pointwise
curvature-to-centroid estimate

\\[
 \\boxed{\\ |m''(s)|\\le2\\sqrt{\\lambda_{\\max}R(s)}\\,\\Phi''(s). } \\tag{2.11}
\\]

This estimate is the shape analogue of the score-range bound in the pure
translation theorem.  Notice that the rotational part of (R') does not
need to commute with (R); only the positive matrix (K=-R'') enters.

### 2.2 Covariance growth and integration

For each unit (u), (h_u(s)=u^TR(s)u) is nonnegative and concave.  By
(2.2), \\(\\mathbb Eh_u\\le\\Lambda\\).  Applying (1.4) gives

\\[
 \\lambda_{\\max}R(s)\\le L_\\rho\\Lambda(1+|s|).              \\tag{2.12}
\\]

Choose a mode (s_0) with (0\\in\\partial\\Phi(s_0)), and put (b=m'(s_0))
(in the smooth case (b=m'(s_0)); for a nonsmooth approximation take a
one-sided representative).  Since \\(\\Phi''\\) is a nonnegative measure,
(2.11)--(2.12) imply, for a.e. (s),
+
+\\[
+ |m'(s)-b|\\le
+ 2\\sqrt{L_\\rho\\Lambda}\\,\\sqrt{1+|s|+|s_0|}\\,|\\Phi'(s)|
+ \\le C_1\\sqrt{\\Lambda}\\sqrt{1+|s|}\\,|\\Phi'(s)|,          \\tag{2.13}
+\\]
+
+where \\(|s_0|\\le L_\\mathrm{mode}\\) was used in the last inequality.  Indeed,
+integrate (2.11) between (s_0) and (s), bound the covariance factor by
+the supremum on that interval, and use
+\\(\\int_{s_0}^s\\Phi''=|\\Phi'(s)|\\).
+
+The mixed-isotropy identity is
+
+\\[
+ 0=\\mathbb E[\\tau(S)m'(S)]=\\mathbb E\\tau,b+\\mathbb E[\\tau(m'-b)],\\qquad
+ \\mathbb E\\tau=1.                                         \\tag{2.14}
+\\]
+
+Thus (1.3) and (2.13) give
+
+\\[
+ |b|^2\\le C_1^2\\Lambda\\,8(1+M_3),\\qquad
+ \\mathbb E[\\tau^2|m'-b|^2]\\le C_1^2\\Lambda\\,8(1+M_3).       \\tag{2.15}
+\\]
+
+Finally, \\(|m'|^2\\le2|m'-b|^2+2|b|^2\\) and (1.1) imply
+
+\\[
+ \\mathbb E[\\tau^2|m'|^2]\\le
+ 2C_1^2\\Lambda\\,8(1+M_3)(1+M_\\tau)=:C_\\Lambda.             \\tag{2.16}
+\\]
+
+With the deliberately rough constants above, (C_\\Lambda<10^{13}\\Lambda).
+The important point is that it is independent of (d\\) and of the spectrum
+of (R).  A finite-difference/measure version follows by approximating
+(W,m,R) on compact subintervals; (2.7) is stable under this approximation.
+
+### 2.3 Relation with (F=-m'+F^0) and the deformation budget
+
+For a Gaussian conditional family the continuity equation has an affine
+gradient velocity.  In centered coordinates (y=z-m(s)),
+
+\\[
+ F_s=-m'(s)+F_s^0(y),\\qquad F_s^0(y)=K_s y,\\qquad
+ K_s=K_s^T,                                               \\tag{2.17}
+\\]
+
+where (K_s) is the unique symmetric solution of
+
+\\[
+ QK_s+K_sQ=-QR'Q.                                        \\tag{2.18}
+\\]
+
+Consequently (C_s=\\|K_s\\|_{\\rm HS}^2).  The exact budget
+\\(\\int\\rho\\tau^2 C_s\\le1\\) from (0.3)--(0.4) pays the changing-shape
+component, while (2.11)--(2.16) pay the centroid component.  No conditional
+Poincare estimate is used.  In particular, noncommuting covariance
+velocities are harmless in this Gaussian class.
+
+## 3. Interval fibers: wedges and cones
+
+We next stress-test the mechanism on a hard-support model where the slices
+are intervals.  Let
+
+\\[
+ K=\\{(s,z):s\\in J,\\ l(s)\\le z\\le u(s)\\},\\qquad
+ m=(u+l)/2,\\quad w=(u-l)/2>0.                            \\tag{3.1}
+\\]
+
+The set (K) is convex exactly when (u) is concave and (l) is convex.
+In the sense of distributions write (\\mu_+=-u''\\ge0) and
+\\(\\mu_-=l''\\ge0\\).  Then
+
+\\[
+ Dm'=\\tfrac12(-\\mu_++\\mu_-),\\qquad
+ -Dw'=\\tfrac12(\\mu_++\\mu_-),\\qquad
+ \\boxed{|Dm'|\\le-Dw'.}                                  \\tag{3.2}
+\\]
+
+The marginal density is \\(\\rho(s)\\propto w(s)\\), so
+\\(\\Phi=-\\log w+\\mathrm{const}\\), and
+
+\\[
+ D\\Phi'=(w'/w)^2\\,ds+\\frac{-Dw'}{w}.                    \\tag{3.3}
+\\]
+
+Hence \\(|Dm'|\\le w\\,D\\Phi'\\).  If the transverse variance is normalized
+so that \\(\\mathbb E_\\rho w(S)^2\\le3\\) (the interval variance is (w^2/3)),
+the concave-growth lemma (1.4), applied to \\(w\\), gives
+\\(w(s)\\le L(1+|s|)\\).  Integrating as in (2.13) now gives
+
+\\[
+ |m'(s)-b|\\le C(1+|s|)|\\Phi'(s)|,\\qquad
+ b\\text{ fixed by }\\mathbb E[\\tau m']=0.                \\tag{3.4}
+\\]
+
+Using (1.2) and the one-dimensional fourth moment bound,
+\\(\\mathbb E[\\tau^2|m'|^2]\\le C\\).  This includes triangular wedges and
+cones (piecewise affine (u,l)); curvature atoms at apices are handled by
+(3.2).  The ideal uniform fibers have a centered deformation velocity
+
+\\[
+ v_s^0(y)=(w'/w)y,\\qquad C_s^{\\rm geom}=(w'/w)^2,             \\tag{3.5}
+\\]
+
+so the same calculation explicitly charges width changes.  Passing from hard
+fibers to the unit-Gaussian-smoothed (q_s) is a quantitative step not
+provided here; one cannot silently identify (C_s^{\\rm geom}\\) with the
+post-noise Hessian charge.
+
+## 4. Rotating boxes and matrix shape changes
+
+A volume-preserving pure rotation of a conditional ellipsoid/box cannot be
+a jointly convex family.  For Gaussian ellipsoids this is immediate from
+\\(R''\\preceq0\\): if all eigenvalues of (R) were constant, then
+\\(\\operatorname{tr}R''=0\\), hence \\(R''=0\\), and constancy of
+\\(\\operatorname{tr}R^2\\) forces \\(R'\\equiv0\\).  Thus rotation must be accompanied
+by a change of eigenvalues, and the determinant curvature in (2.10) pays the
+Frobenius square of the relative velocity.
+
+For an explicit noncommuting test, take on a small interval
+\\[
+ R(s)=I+sA-s^2B,\\quad
+ A=\\begin{pmatrix}1&0\\\\0&-1\\end{pmatrix},\\quad
+ B=\\begin{pmatrix}1&\\varrho\\\\\\varrho&1\\end{pmatrix},\\quad0<|\\varrho|<1.
+\\]
+Then (R''=-2B\\preceq0), while \\([A,B]\\ne0\\).  Choose (m'') in
+\\operatorname{Ran}B and (W''\\ge\\tfrac12(m'')^TB^\\dagger m''\\); (2.7)
+then gives a jointly convex density.  Equations (2.10)--(2.16) give the
+same dimension-free centroid bound, regardless of the angle at which the
+eigenspaces rotate.
+
+## 5. What remains for arbitrary changing shapes
+
+The preceding proofs rely on a structural positive object controlling
+centroid acceleration:
+
+* for Gaussian slices, the matrix curvature (K=-R'') and the exact Schur
+  complement (2.7);
+* for interval fibers, the width curvature (-w'') and the endpoint
+  inequalities.
+
+For a general conditional potential, the exact identity from
+weighted_fisher_prekopa.md is only
+
+\\[
+ m''=-\\mathbb E_s[y(\\mathcal Q_s+\\|D F_s^0\\|_{\\rm HS}^2)]. \\tag{5.1}
+\\]
+
+The known budget controls the *mass* of the nonnegative charge,
+\\(\\int\\rho\\tau^2\\mathbb E_s(\\mathcal Q_s+\\|DF_s^0\\|^2)\\le1\\), but not its
+conditional first moment in (y).  Without an analogue of (2.7) or (3.2),
+an arbitrary nonnegative charge can concentrate in a far tail and make
+\\(|m''|\\) large while its mass is small.  Bounding (5.1) by conditional
+Cauchy--Schwarz would require a dimension-free conditional second-moment or
+Poincare estimate, which is precisely the forbidden KLS-strength step.
+
+Thus no claim of the full weighted slice-rigidity lemma is made here.  The
+Gaussian and wedge calculations are genuine shape-changing progress and
+provide stress tests for any proposed general inequality; a proof extending
+their curvature-to-centroid mechanism to arbitrary log-concave fibers is the
+remaining load-bearing problem.

