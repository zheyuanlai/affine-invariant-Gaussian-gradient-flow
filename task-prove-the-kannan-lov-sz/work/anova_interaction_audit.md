# ANOVA interaction operator for (S=(X+G)/\sqrt 2)

This note audits the proposed estimate

\[
 \operatorname {dist}_{L^2(\nu)}(f,\mathrm {Aff})^2
 \le C\,\mathbb E R^2,\qquad
 R=f(S)-\mathbb E[f(S)\mid X]-\mathbb E[f(S)\mid G],
 \tag{A}
\]

where (X) is centered isotropic log-concave, (G\sim N(0,I)) is
independent, and \(\nu=\mathcal L(S)\).  The conclusions below are
rigorous and do not use KLS or a second-order Poincare inequality.

## 1. Hilbert-space identity and the operator (Q_X+Q_G)

Let (H=L^2_0(\nu)), and let
\[
 (Tf)(x,g)=f((x+g)/\sqrt2),\qquad T:H\longrightarrow L^2(\mu\otimes\gamma).
\]
Because (S) has law \(\nu), (T) is an isometry.  Let (P_X,P_G) be
the conditional-expectation projections onto the centered subspaces of
functions of (X) and of (G), respectively.  Independence gives
\(P_XP_G=P_GP_X=0\) on centered functions.  If
\[
 U_Xf=\mathbb E[Tf\mid X],\qquad U_Gf=\mathbb E[Tf\mid G],
\]
then (U_X,U_G) are contractions and (Q_X=U_X^*U_X=T^*P_XT),
\(Q_G=U_G^*U_GT^*P_GT\) are positive contractions on (H).  (The
second display means (Q_G=T^*P_GT); no extra (T^*) is intended.)

For every centered (f\in H),
\[
 \boxed{\quad
 \mathbb E R^2
 =\|f\|_{L^2(\nu)}^2-\|U_Xf\|_{L^2(\mu)}^2
                    -\|U_Gf\|_{L^2(\gamma)}^2
 =\langle f,(I-Q_X-Q_G)f\rangle.\quad}       \tag{1.1}
\]
Indeed, the two conditional means are independent and centered, so their
cross inner product is zero.  Positivity of (I-Q_X-Q_G) follows directly
from orthogonality of (P_XT f) and (P_GT f).  Thus (A) is exactly a
dimension-free spectral lower bound for (I-Q_X-Q_G) on the orthogonal
complement of the affine space.

The identity is purely (L^2); it holds for every centered form-domain
function and, in fact, for every centered (L^2(\nu)) function.  No
regularity of (f) is needed to define (R).

## 2. Kernel and the lower-dimensional warning

Assume that (\mu) is full-dimensional, so its density is positive on the
interior of a convex set (K).  If the right side of (1.1) is zero, then
\[
 f((x+g)/\sqrt2)=u(x)+v(g)                         \tag{2.1}
\]
for \(\mu\otimes\gamma\)-almost every \((x,g)).  Comparing two admissible
values of (x) shows that, for every (d) in a neighborhood of zero,
\(f(y+d)-f(y)\) is a.e. independent of (y).  The resulting increment
function is measurable and additive, hence linear.  Therefore (f) is
affine a.e.  Conversely every affine (f) has (R=0).  This proves
\[
 \ker(I-Q_X-Q_G)=\mathrm {Aff}\quad\text{(full-dimensional case)}. 
 \tag{2.2}
\]
The argument is distributional and extends to arbitrary (L^2(\nu)) by
localization; no attainment or smoothness is used.

If “isotropic” is interpreted only on a lower-dimensional supporting affine
hull while (G) is taken in the ambient space, (A) is false exactly.  For
example, in (\mathbb R^2), let (X=(U,0)), where
\(U\sim\mathrm{Unif}[-\sqrt3,\sqrt3]\), and let (G=(G_1,G_2)\) be standard
Gaussian.  Then (S_2=G_2/\sqrt2).  With
\[
 f(s)=s_2^2-\tfrac12,
\]
we have (\mathbb Ef(S)=0), (\operatorname {dist}(f,\mathrm {Aff})^2
=\operatorname {Var}(S_2^2)=1/2), while
\[
 h_X=\mathbb E[f(S)\mid X]=0,\qquad h_G=f(S),\qquad R=0.
\]
Thus a lower-dimensional reduction must take the Gaussian on the supporting
hull (or project both sides to that hull) before applying (A).

## 3. Exact Gaussian spectrum (sharp model)

For (X\sim N(0,I)), also \(\nu=N(0,I)).  Let (H_\alpha) be normalized
multivariate Hermite polynomials and (m=|\alpha|\).  The Gaussian addition
formula gives
\[
 \mathbb E[H_\alpha((X+G)/\sqrt2)\mid X]
   =2^{-m/2}H_\alpha(X),
\]
and the same formula with (X,G) interchanged.  Hence
\[
 Q_XH_\alpha=Q_GH_\alpha=2^{-m}H_\alpha,
 \qquad
 \mathbb E R^2=(1-2^{1-m})\|H_\alpha\|_2^2.       \tag{3.1}
\]
The degree-one space is exactly the affine kernel.  For every (m\ge2),
the affine distance is the full norm and the smallest defect is (1/2) at
degree two.  Therefore (A) holds for the Gaussian with the optimal constant
\(C=2\).

## 4. Quadratic modes: a KLS-strength obstruction

Let (A=A^T), (T_A=\operatorname {tr}A), and
\[
 f_A(s)=s^TAs-T_A.
\]
For every centered isotropic (X), direct expansion gives
\[
 h_X(x)=\tfrac12(x^TAx-T_A),\quad
 h_G(g)=\tfrac12(g^TAg-T_A),\quad
 R=x^TAg,                                      \tag{4.1}
\]
so
\[
 \mathbb E R^2=\operatorname {tr}(A^2).          \tag{4.2}
\]
Put (Q=X^TAX) and (d_A=\mathbb E[X(Q-T_A)]).  Since
\(\mathbb E[Sf_A]=d_A/(2\sqrt2)), isotropy of (S) yields
\[
 \operatorname {dist}_{L^2(\nu)}(f_A,\mathrm {Aff})^2
 =\tfrac14\operatorname {Var}(Q)+\tfrac32\operatorname {tr}(A^2)
       -\tfrac18|d_A|^2.                            \tag{4.3}
\]
Consequently, for centrally symmetric (X) (so (d_A=0)), (A) would
imply
\[
 \operatorname {Var}(X^TAX)\le (4C-6)\operatorname {tr}(A^2). 
 \tag{4.4}
\]
Taking (A=I) gives the dimension-free thin-shell variance bound
\(\operatorname {Var}|X|^2\le (4C-6)n) for every centrally symmetric
isotropic log-concave law.  This is a known KLS-strength problem, not a
consequence of covariance normalization.  Thus a proof of (A) cannot be
claimed from elementary Gaussian/Poincare or Hessian identities; it would
already settle this substantial nonlinear sector.

For a general (possibly skew) law, (4.3) is still exact after affine
projection.  In one dimension the marginal fourth-moment bounds make the
ratio finite in all standard examples; no skew first-order obstruction
survives the affine projection.

## 5. Product-exponential, cube, and hard-support checks

For a centered one-sided exponential (X=E-1), (E\sim\mathrm{Exp}(1)),
\(\mathbb EX^3=2\), \(\mathbb EX^4=9\).  For the coordinate quadratic
\(f(s)=s_1^2-1\), (4.1)--(4.3) give
\[
 \operatorname {dist}(f,\mathrm {Aff})^2=3,
 \qquad \mathbb E R^2=1,
\]
so this mode requires (C\ge3).  For the symmetric Laplace law
\(\mathbb EX^4=6\), the corresponding ratio is
\(1/(11/4)=4/11\), hence (C\ge11/4).  For the isotropic cube
\([-√3,\u221a3]^n), (d_A=0) and, for (A=I),
\(\operatorname {Var}|X|^2=4n/5); the ratio in (A) is
\[
 \frac{\mathbb E R^2}{\operatorname {dist}(f_A,\mathrm {Aff})^2}
 =\frac{n}{(2/5)n+(3/2)n}=\frac58.           \tag{5.1}
\]
Uniform Euclidean balls and regular simplices give analogous order-one
ratios after their standard isotropic normalization.  These computations
are consistency checks, not a proof of the universal estimate.

For product measures, conditional expectation factorizes over coordinates.
Tensor-product polynomial modes therefore have defects bounded away from
zero (e.g. (f=s_i s_j) has (\mathbb E R^2=1/2) after normalization),
and no dimension loss appears in these explicit product tests.  A complete
tensorization proof for arbitrary (f) would require a separate spectral
analysis of the noncommuting one-dimensional operators (Q_X,Q_G); it is
not asserted here.

## 6. What is and is not established

* The operator identity (1.1), its continuous/form-domain scope, the full
  dimensional kernel statement (2.2), and all model calculations above are
  proved exactly.
* The Gaussian constant is (C=2), and elementary hard-support/product
  tests are uniformly nondegenerate.
* The universal full-dimensional estimate (A) is not proved by these facts.
  Its quadratic restriction already contains the dimension-free thin-shell
  variance problem.  No explicit full-dimensional isotropic log-concave law
  with a vanishing ANOVA gap is currently exhibited here.
* If ambient lower-dimensional supports are admitted without reducing the
  Gaussian and the affine space to the support, the displayed two-dimensional
  example is a rigorous counterexample.

Accordingly, (A) should not be inserted as an established lemma in a KLS
proof unless a genuinely new argument controlling the quadratic variance
sector (and the higher nonlinear sectors) is supplied.

