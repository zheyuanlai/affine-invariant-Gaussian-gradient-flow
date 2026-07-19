# ANOVA interaction operator: clean corrected version

This is a clean copy of the ANOVA audit; all formulas below are literal
LaTeX and the statements are self-contained.

## Operator identity

Let \(\nu=\mathcal L((X+G)/\sqrt2)\), \(H=L^2_0(\nu)\), and
\[
(Tf)(x,g)=f((x+g)/\sqrt2).
\]
Let \(P_X,P_G\) be the centered conditional-expectation projections in
\(L^2(\mu\otimes\gamma)\).  For
\[
 U_Xf=\mathbb E[Tf\mid X],\qquad U_Gf=\mathbb E[Tf\mid G],
\]
put \(Q_X=U_X^*U_X=T^*P_XT\) and
\(Q_G=U_G^*U_G=T^*P_GT\).  Independence gives \(P_XP_G=0\).
For every centered \(f\in L^2(\nu)\),
\[
\boxed{\quad
\mathbb E R^2
=\|f\|_2^2-\|U_Xf\|_2^2-\|U_Gf\|_2^2
=\langle f,(I-Q_X-Q_G)f\rangle.
\quad} \tag{1}
\]
This is an exact bounded-operator identity, so it includes all
form-domain and non-smooth functions.  The operator is positive because
\(P_XTf\) and \(P_GTf\) are orthogonal.

## Full-dimensional kernel and degenerate-support counterexample

If \(\mu\) is full-dimensional log-concave, \(R=0\) implies
\[
f((x+g)/\sqrt2)=u(x)+v(g)
\]
almost everywhere.  Differences of two \(x\)'s in the interior of the
support contain a neighborhood of zero, so all small increments of \(f\)
are a.e. constant in the base point.  Measurable Cauchy rigidity then
gives \(f\in\mathrm{Aff}\) a.e.; the converse is immediate.

If isotropy is imposed only on a lower-dimensional support but \(G\) is
ambient, the estimate is false.  Take \(X=(U,0)\in\mathbb R^2\) with
\(U\sim{\rm Unif}[-\sqrt3,\sqrt3]\), \(G=(G_1,G_2)\sim N(0,I_2)\), and
\(f(s)=s_2^2-\frac12\).  Then \(S_2=G_2/\sqrt2\),
\(\operatorname{dist}_{L^2}(f,\mathrm{Aff})^2=\operatorname{Var}(S_2^2)=1/2\),
whereas \(h_X=0\), \(h_G=f(S)\), and \(R=0\).  Thus lower-dimensional
reductions must project both \(G\) and \(\mathrm{Aff}\) to the supporting
hull.

## Gaussian model

For \(X\sim\gamma_n\), normalized Hermite degree \(m\) satisfies
\[
Q_XH_\alpha=Q_GH_\alpha=2^{-m}H_\alpha,\qquad
\mathbb E R^2=(1-2^{1-m})\|H_\alpha\|_2^2.
\]
The affine kernel is degree \(0,1\); on its orthogonal complement the
smallest defect is \(1/2\), so \(C_{\rm Gauss}=2\) is sharp.

## Quadratic sector and thin-shell obstruction

For \(A=A^T\), \(f_A(s)=s^TAs-\operatorname{tr}A\), direct expansion gives
\[
h_X(x)=\tfrac12(x^TAx-\operatorname{tr}A),\quad
h_G(g)=\tfrac12(g^TAg-\operatorname{tr}A),\quad
R=x^TAg,
\]
and hence \(\mathbb E R^2=\operatorname{tr}(A^2)\).  If
\(Q=X^TAX\) and \(d_A=\mathbb E[X(Q-\operatorname{tr}A)]\), then
\[
\operatorname{dist}(f_A,\mathrm{Aff})^2
=\tfrac14\operatorname{Var}(Q)+\tfrac32\operatorname{tr}(A^2)
-\tfrac18|d_A|^2. \tag{2}
\]
For centrally symmetric \(X\), \(d_A=0\).  A universal constant in the
proposed inequality would then imply
\[
\operatorname{Var}(X^TAX)\le(4C-6)\operatorname{tr}(A^2),
\]
and with \(A=I\), \(\operatorname{Var}|X|^2=O(n)\), the dimension-free
thin-shell variance bound.  This sector is therefore not an elementary
consequence of isotropy.

## Model values

For a centered one-sided exponential \(X=E-1\), \(E\sim{\rm Exp}(1)\), the
coordinate quadratic has \(\operatorname{dist}^2=3\) and
\(\mathbb E R^2=1\), requiring \(C\ge3\).  For symmetric Laplace the ratio
\(\mathbb E R^2/\operatorname{dist}^2=4/11\).  For the isotropic cube,
\(\operatorname{Var}|X|^2=4n/5\), and the radial quadratic ratio is
\(10/17\).  These hard-support and product checks are nondegenerate.

## Status

The operator identity, form-domain scope, full-dimensional kernel, Gaussian
spectrum, degenerate-support counterexample, and quadratic obstruction are
proved.  No full-dimensional isotropic log-concave counterexample is known
from these calculations, but a universal proof would already establish the
thin-shell variance bound (and control all higher nonlinear sectors).  The
lemma should therefore not be used as an established input without a new
KLS-strength argument.

