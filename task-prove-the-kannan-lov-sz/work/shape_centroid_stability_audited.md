# Shape-centroid stability: audited scope and obstruction

This is an audit of the shape-changing calculations. It deliberately separates
proved smooth subclasses from unproved approximation steps.

## A. Setup and exact missing term

Let \(p(s,z)=e^{-U(s,z)}\), \(\rho(s)=\int e^{-U(s,z)}dz=e^{-\Phi(s)}\),
and \(q_s=e^{-U}/\rho\). Assume \(S\) is centered and variance one, and
write \(m(s)=E_s Z\), \(y=z-m(s)\). For the centered Poisson field,
\[
 F_s=-m'(s)+F_s^0,\quad E_sF_s^0=0,\quad
 R_s=\mathcal Q_s+\|D_zF_s^0\|_{\rm HS}^2\ge0,
\]
\[
 m''(s)=-E_s[yR_s],\qquad
 \int\rho\tau^2 E_sR_s\le1,\qquad
 E[\tau(S)m'(S)]=0.                                    \tag{A.1}
\]
The last three identities are exact in the regular weighted-Fisher setup.
The unresolved issue is the first moment \(E_s[yR_s]\); the mass bound in
(A.1) does not by itself control it.

## B. Fully proved smooth Gaussian subclass

Assume \(J\times\mathbb R^d\) and
\[
 U(s,z)=W(s)+\tfrac12(z-m(s))^TQ(s)(z-m(s)),\qquad
 Q(s)\succ0,\quad R=Q^{-1},
\]
with \(W,m,R\in C^2\), a normalizable density, and joint convexity. Assume
\(S\) is centered/variance one, \(E[S\,m(S)]=0\), and \(ER(S)\preceq\Lambda I\).
Then there is a numerical \(C\) such that
\[
 E[\tau(S)^2|m'(S)|^2]\le C\Lambda.                       \tag{B.1}
\]

**Proof.** The Schur complement of \(Q\) is, for every \(y\),
\[
 W''-m''{}^TQy+\tfrac12y^TQ(-R'')Qy\ge0.                 \tag{B.2}
\]
Putting \(K=-R''\) and \(x=Qy\) gives
\(W''-m''\cdot x+\frac12x^TKx\ge0\). Therefore
\(m''\in\operatorname{Ran}K\) and
\(m''{}^TK^\dagger m''\le2W''\). For every unit \(u\),
\[
 |u\cdot m''|^2\le2W''u^TKu
 \le2W''\lambda_{\max}(R)\operatorname{tr}(R^{-1}K).      \tag{B.3}
\]
The marginal potential satisfies
\[
 \Phi=W-\tfrac12\log\det R,\qquad
 \Phi''=W''+\tfrac12\operatorname{tr}(R^{-1}K)
 +\tfrac12\|R^{-1/2}R'R^{-1/2}\|_{\rm HS}^2,              \tag{B.4}
\]
hence
\[
 |m''|\le2\sqrt{\lambda_{\max}R}\,\Phi''.                 \tag{B.5}
\]
For each unit \(u\), \(u^TRu\) is nonnegative concave. The elementary
one-dimensional core/growth lemma says that a nonnegative concave \(h\) with
\(Eh(S)\le\Lambda\) obeys \(h(s)\le L\Lambda(1+|s|)\); therefore
\(\lambda_{\max}R(s)\le L\Lambda(1+|s|)\).

Choose \(s_0\) with \(0\in\partial\Phi(s_0)\) and set \(b=m'(s_0)\).
The standard one-dimensional mode bound gives \(|s_0|\le4\). Integrating
(B.5) along the segment from \(s_0\) to \(s\) gives
\[
 |m'(s)-b|\le C_1\sqrt{\Lambda}\sqrt{1+|s|}\,|\Phi'(s)|.   \tag{B.6}
\]
For a centered variance-one log-concave density,
\(\tau|\Phi'|\le1+|s|\) and \(E|S|^4\le M_4\). Consequently
\(E[\tau^2|m'-b|^2]\le C_2\Lambda\). Since
\(E[\tau m']=0\), \(E\tau=1\), Cauchy--Schwarz gives
\(|b|^2\le C_2\Lambda\); using \(E\tau^2\le400\) proves (B.1).
All constants are independent of \(d\). \(\square\)

For this Gaussian family the centered Poisson field is
\(F_s=-m'(s)+A_sy\), where \(A_s=A_s^T\) uniquely solves
\(QA_s+A_sQ=-QR'Q\); hence \(C_s=\|A_s\|_{\rm HS}^2\). The deformation
budget in (A.1) applies when the family is embedded in the regular
weighted-Fisher construction. No claim is made here for nonsmooth Gaussian
limits; preserving the \(ER\) bound and the weighted form budget under such
limits needs a separate argument.

## C. Fully proved geometric interval calculation (not post-noise WFI)

Let \(K=\{(s,z):l(s)\le z\le u(s)\}\) be a planar convex body, with uniform
density, \(m=(u+l)/2\), \(w=(u-l)/2\). Assume the \(S\)-marginal is centered
variance one, \(E[S m(S)]=0\), and \(E_\rho w^2\le3\). In distributions,
\[
 |D m'|\le-Dw',\qquad
 D\Phi'=(w'/w)^2ds+(-Dw')/w,\quad \rho\propto w.          \tag{C.1}
\]
The same core/growth lemma gives \(w(s)\le L(1+|s|)\). Integrating (C.1)
from a mode and using \(\tau|\Phi'|\le1+|s|\), then the fourth moment bound,
gives \(E[\tau^2|m'|^2]\le C\), with \(C\) numerical. Piecewise-affine
wedges/cones and curvature atoms are included. The geometric centered
velocity is \(v_s^0=(w'/w)y\). This is a hard-support geometric result only:
identifying it with the post-Gaussian Poisson field is not justified.

## D. Cone counterstress (fully verified elsewhere)

The isotropic cone and exact post-noise Poisson field in
work/gaussian_curvature_korn.md, §§2--4, have \(m=0\),
\(C_s\le a^{-2}\), but \(B_s\ge(256a)^{-1}\) in a boundary layer. Thus
\(B_s/(C_s+|E_sF_s|^2)\ge a/256\to\infty\), while the full weighted
integral satisfies \(\int\rho\tau^2B_s\le3/2\). This rules out any
slice-by-slice curvature--Korn closure but does not contradict the integrated
target.

## E. Rotating boxes

For Gaussian ellipsoids, joint convexity forces \(R''\preceq0\). A
volume-preserving pure rotation would have constant eigenvalues, hence
\(R''=0\) and then \(R'=0\). The local path
\(R=I+sA-s^2B\), \(A=\operatorname{diag}(1,-1)\),
\(B=\left(\begin{smallmatrix}1&r\\r&1\end{smallmatrix}\right)\), \(0<|r|<1\),
has \(R''=-2B\preceq0\) and \([A,B]\ne0\); it is a stress test, not a
globally normalized counterexample. The Gaussian theorem above handles it
whenever the global \(ER\) and mixed-isotropy hypotheses are imposed.

## F. Status

The Gaussian and hard-interval statements are genuine dimension-free
subclass results. The general shape-changing weighted slice-rigidity lemma
remains unproved. Any argument that replaces \(E_s[yR_s]\) by conditional
Cauchy--Schwarz or a slice Poincare bound is a forbidden KLS-strength step.

