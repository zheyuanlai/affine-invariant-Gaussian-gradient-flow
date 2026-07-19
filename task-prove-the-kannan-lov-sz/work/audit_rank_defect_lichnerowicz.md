# Clean-room audit of the rank-defect Lichnerowicz lemmas

## 0. Scope and verdict

This audit uses only the three lemma statements supplied in the audit request. It
does not use the proofs in rank_one_defect_lichnerowicz.md.

Let \(H=F\oplus E\) be an orthogonal decomposition of a finite-dimensional
Euclidean space. Write \(P_E\) for the orthogonal projection onto \(E\), and
use the convention that the Gaussian \(\gamma_\varepsilon\) has covariance
\(\varepsilon I\).

The verdict is:

1. **(L1) is valid**, with the factor \(2\), for locally Lipschitz
   \(h\in L^2(\eta)\), with the usual convention that a divergent right-hand
   side makes the assertion vacuous. The closed gradient subspace is a
   reducing subspace of the one-form Witten Laplacian; this cannot merely be
   assumed, but follows from the scalar/one-form intertwining (or, equivalently,
   from the polar decomposition of the closed gradient operator).
2. **(L2) is valid** with the stated constant

   \[
   C_P(\mu)\le \kappa^{-1}+2\bigl(C_P(\pi_F{}_{\#}\mu)+\kappa^{-1}\bigr)
   =2C_P(\pi_F{}_{\#}\mu)+3\kappa^{-1}.
   \]

   The load-bearing estimate is a conditional kinetic-tensor bound
   \(B_x^*B_x\preceq \kappa^{-1}D^2W(x)\). Its sign and its proof are given
   explicitly below.
3. **(L3) is valid**: convolution with \(\gamma_\varepsilon\) changes
   \(\kappa\) to \(\kappa/(1+\kappa\varepsilon)\), and

   \[
   C_P(\nu*\gamma_\varepsilon)
   \le C_P(\nu)+\varepsilon.
   \]

The bare statements need one convention in lower dimension: all Hessians,
gradients, Lebesgue densities, and orthogonal decompositions must be taken in
the affine hull of the measure. An ambient Hessian in a direction normal to
the support is not defined by the measure and cannot be used in (L1) or (L2).

The proofs below are first given for positive \(C^3\) densities for which the
displayed integrations are finite. Compactly supported smooth test functions
and quadratic confinement give a direct core. Gaussian convolution, proved
in Section 3, then removes this auxiliary regularity without changing the
limiting constants; the limiting argument is written out in Section 4.

---

## 1. Audit and proof of (L1)

### 1.1 Precise statement

Let

\[
d\eta=e^{-U}\,dx/Z,
\qquad U\in C^2(\mathbb R^m),
\qquad D^2U\succeq0,
\]

and suppose \(0<C_P(\eta)<\infty\). Set

\[
\lambda=C_P(\eta)^{-1}.
\]

Then every locally Lipschitz \(h\in L^2(\eta)\) satisfies

\[
\operatorname{Var}_\eta h
\le
2\int
\left\langle (D^2U+\lambda I)^{-1}\nabla h,\nabla h\right\rangle\,d\eta.
\tag{1.1}
\]

### 1.2 Scalar and one-form operators

Let \(d\) denote the closed gradient operator

\[
d:H^1(\eta)\subset L^2(\eta)\longrightarrow
L^2(\eta;\mathbb R^m).
\]

Let

\[
A=d^*d=-\Delta+\nabla U\cdot\nabla
\]

be its nonnegative Friedrichs realization. On the mean-zero space
\(L^2_0(\eta)\), the Poincare inequality is exactly the form inequality

\[
A\succeq\lambda I.
\tag{1.2}
\]

Let

\[
\mathcal G=\overline{\operatorname{Ran}d}
\subset L^2(\eta;\mathbb R^m)
\]

be the closed gradient subspace. In fact the range is closed here: if
\(u\in H^1(\eta)\cap L^2_0(\eta)\), then

\[
\|du\|_2^2\ge\lambda\|u\|_2^2.
\]

The polar decomposition on \(L^2_0(\eta)\) is

\[
d=Q A^{1/2},
\tag{1.3}
\]

where \(Q:L^2_0(\eta)\to\mathcal G\) is unitary. Define

\[
B_{\mathcal G}=QAQ^*.
\tag{1.4}
\]

For smooth compactly supported \(u\), direct differentiation gives

\[
d(Au)=\bigl(A\otimes I+D^2U\bigr)du.
\tag{1.5}
\]

Thus \(B_{\mathcal G}\) is precisely the restriction to \(\mathcal G\) of the
one-form Witten Laplacian. Its quadratic form on smooth exact fields, and then
by closure, is

\[
q_1(v)=
\int\left(\|Dv\|_{\mathrm{HS}}^2+
\langle D^2U\,v,v\rangle\right)d\eta.
\tag{1.6}
\]

This also checks the often omitted invariance issue. Equation (1.5) shows
invariance on an operator core. Equivalently, if
\(w\in\mathcal G^\perp=\ker d^*\), then the Hilbert-complex identity gives
\(\langle Bdu,w\rangle=\langle dAu,w\rangle=0\); closure of the forms shows
that \(\mathcal G\) is reducing for the one-form resolvent and semigroup. No
projection of a multiplication operator onto \(\mathcal G\) is being silently
identified with the original multiplication operator.

### 1.3 The two form lower bounds

Unitary equivalence in (1.4) and (1.2) give

\[
q_1(v)\ge\lambda\|v\|_2^2,
\qquad v\in D(q_1)\cap\mathcal G.
\tag{1.7}
\]

The Bochner form (1.6), together with \(D^2U\succeq0\), gives independently

\[
q_1(v)\ge
\int\langle D^2U\,v,v\rangle\,d\eta.
\tag{1.8}
\]

Averaging (1.7) and (1.8), rather than adding them, is the source of the
constant \(2\):

\[
q_1(v)\ge
\frac12\int
\langle(D^2U+\lambda I)v,v\rangle\,d\eta.
\tag{1.9}
\]

### 1.4 Exact inverse identity and the variational constant

For centered \(h\in H^1(\eta)\), put \(g=dh\). Equations (1.3)--(1.4) imply

\[
\begin{aligned}
\langle B_{\mathcal G}^{-1}g,g\rangle
&=\left\langle
QA^{-1}Q^*QA^{1/2}h,QA^{1/2}h
\right\rangle\\
&=\|h\|_2^2
=\operatorname{Var}_\eta h.
\end{aligned}
\tag{1.10}
\]

The inverse quadratic form has the exact dual representation

\[
\langle B_{\mathcal G}^{-1}g,g\rangle
=\sup_{v\in D(q_1)\cap\mathcal G}
\{2\langle g,v\rangle-q_1(v)\}.
\tag{1.11}
\]

Let \(M(x)=D^2U(x)+\lambda I\). Using (1.9), then enlarging the
supremum from exact fields to all vector fields in the multiplication-form
domain,

\[
\begin{aligned}
\operatorname{Var}_\eta h
&\le
\sup_v\left\{2\langle g,v\rangle
-\frac12\int\langle Mv,v\rangle d\eta\right\}\\
&=2\int\langle M^{-1}g,g\rangle d\eta.
\end{aligned}
\tag{1.12}
\]

The maximizer in the last line is \(v=2M^{-1}g\). This verifies both the
factor \(2\) and the direction of the compression/inverse comparison. In
particular, one must not replace the inverse of the compression of \(M\) to
\(\mathcal G\) by the compression of \(M^{-1}\) as an equality; only the
variational inequality used above is valid.

For a locally Lipschitz \(h\in L^2(\eta)\), first truncate its values and
multiply by cutoffs \(\chi_R\) satisfying
\(0\le\chi_R\le1\), \(\chi_R=1\) on \(B_R\), and
\(|\nabla\chi_R|\le 2/R\). Apply (1.12) to the bounded compactly supported
approximants. Since \(M^{-1}\preceq\lambda^{-1}I\), the cutoff-gradient
error is bounded by

\[
\frac{4}{\lambda R^2}\int_{|x|\ge R}|h|^2d\eta\longrightarrow0.
\]

Value truncation and lower semicontinuity then give (1.1). If the right-hand
side of (1.1) is infinite, the assertion is automatic.

---

## 2. Audit and proof of (L2)

### 2.1 Conditional notation and the variance decomposition

Assume first that

\[
d\mu(x,y)=e^{-V(x,y)}\,dx\,dy/Z,
\qquad (x,y)\in F\oplus E,
\tag{2.1}
\]

where \(V\in C^3\), and

\[
D^2V(x,y)\succeq\kappa P_E,
\qquad\kappa>0.
\tag{2.2}
\]

Let \(\nu=\pi_F{}_{\#}\mu\), write

\[
d\nu(x)=e^{-W(x)}dx/Z_F,
\qquad
W(x)=-\log\int_E e^{-V(x,y)}dy,
\tag{2.3}
\]

and let \(\eta_x\) be the conditional law of \(y\) given \(x\). For a smooth
test function \(f\), set

\[
g(x)=\int_E f(x,y)d\eta_x(y).
\]

Since \(D^2_{yy}V\succeq\kappa I_E\), conditional Lichnerowicz gives

\[
\int_F\operatorname{Var}_{\eta_x}f(x,\cdot)d\nu(x)
\le\kappa^{-1}\int|\nabla_Ef|^2d\mu.
\tag{2.4}
\]

Thus

\[
\operatorname{Var}_\mu f
\le \kappa^{-1}\int|\nabla_Ef|^2d\mu
+\operatorname{Var}_\nu g.
\tag{2.5}
\]

The remaining issue is that differentiating \(g\) produces motion of the
conditional law. Dropping that term, or assigning it the wrong sign, would
invalidate the proof.

### 2.2 Conditional continuity equation, with its sign checked

Fix \(u\in F\), and define the centered conditional score

\[
q_{x,u}(y)=\partial_uV(x,y)
-\int\partial_uV(x,\cdot)d\eta_x.
\tag{2.6}
\]

Let

\[
A_x=-\Delta_y+\nabla_yV(x,\cdot)\cdot\nabla_y
\]

on mean-zero \(L^2(\eta_x)\), and solve

\[
A_x\phi_{x,u}=q_{x,u}.
\tag{2.7}
\]

Conditional Lichnerowicz makes \(A_x^{-1}\) bounded on the mean-zero space.
Define the conditional velocity

\[
B_xu=w_{x,u}:=-\nabla_y\phi_{x,u}
\in L^2(\eta_x;E).
\tag{2.8}
\]

Since

\[
\partial_u\eta_x=-q_{x,u}\eta_x
\]

and

\[
\operatorname{div}_y(\eta_xw_{x,u})
=\operatorname{div}_y(-\eta_x\nabla_y\phi_{x,u})
=q_{x,u}\eta_x,
\]

the correct continuity equation is

\[
\partial_u\eta_x
+\operatorname{div}_y(\eta_xw_{x,u})=0.
\tag{2.9}
\]

Consequently

\[
\partial_ug(x)
=\int\partial_uf\,d\eta_x
+\int\langle\nabla_Ef,w_{x,u}\rangle d\eta_x.
\tag{2.10}
\]

The second term equals
\(-\operatorname{Cov}_{\eta_x}(f,\partial_uV)\). Choosing
\(+\nabla\phi\) in (2.8) would reverse this sign and would not solve (2.9).

### 2.3 Kinetic tensor versus marginal Hessian

Define the conditional kinetic tensor on \(F\) by

\[
G_x=B_x^*B_x,
\qquad
\langle G_xu,u\rangle
=\int|w_{x,u}|^2d\eta_x.
\tag{2.11}
\]

We claim the dimension-free matrix inequality

\[
G_x\preceq\kappa^{-1}D^2W(x).
\tag{2.12}
\]

It suffices to check a fixed \(u\in F\). Differentiating (2.3) gives

\[
\partial_{uu}^2W
=\int V_{uu}d\eta_x-\int q_{x,u}^2d\eta_x.
\tag{2.13}
\]

The conditional Bochner identity applied to (2.7) is

\[
\int q_{x,u}^2d\eta_x
=\int\|D_y^2\phi_{x,u}\|_{\mathrm{HS}}^2d\eta_x
+\int\langle V_{yy}\nabla_y\phi_{x,u},
\nabla_y\phi_{x,u}\rangle d\eta_x.
\tag{2.14}
\]

Moreover, integration by parts gives

\[
\int\langle V_{yu},\nabla_y\phi_{x,u}\rangle d\eta_x
=\int\langle\nabla_yq_{x,u},\nabla_y\phi_{x,u}\rangle d\eta_x
=\int q_{x,u}^2d\eta_x.
\tag{2.15}
\]

Apply (2.2) pointwise to the vector
\((u,-\nabla_y\phi_{x,u})\). After integration this gives

\[
\begin{aligned}
0\le{}&\int V_{uu}d\eta_x
-2\int q_{x,u}^2d\eta_x\\
&+\int\langle(V_{yy}-\kappa I)\nabla_y\phi_{x,u},
\nabla_y\phi_{x,u}\rangle d\eta_x.
\end{aligned}
\tag{2.16}
\]

Substituting (2.14) into (2.16) yields the stronger identity-bound

\[
\partial_{uu}^2W
\ge
\int\|D_y^2\phi_{x,u}\|_{\mathrm{HS}}^2d\eta_x
+\kappa\int|\nabla_y\phi_{x,u}|^2d\eta_x.
\tag{2.17}
\]

In particular, (2.12) follows. This is the load-bearing estimate: it uses
the full block condition \(D^2V\succeq\kappa P_E\), not merely
\(V_{yy}\succeq\kappa I_E\).

### 2.4 Application of (L1) on the marginal

Let

\[
c=C_P(\nu),
\qquad\lambda_F=c^{-1},
\qquad M_x=D^2W(x)+\lambda_FI_F.
\tag{2.18}
\]

The case \(\dim F=0\) is handled separately at the end. Define

\[
a(x)=\int\nabla_Ff(x,y)d\eta_x(y),
\qquad
z_x(y)=\nabla_Ef(x,y),
\qquad
r(x)=B_x^*z_x.
\tag{2.19}
\]

Equation (2.10) says exactly

\[
\nabla_Fg=a+r.
\tag{2.20}
\]

Apply (L1) to \(\nu\):

\[
\operatorname{Var}_\nu g
\le2\int\langle M_x^{-1}(a+r),a+r\rangle d\nu(x).
\tag{2.21}
\]

Since \(M_x\succeq c^{-1}I\), Jensen's inequality gives

\[
\int\langle M_x^{-1}a,a\rangle d\nu
\le c\int|\nabla_Ff|^2d\mu.
\tag{2.22}
\]

From (2.12),

\[
B_x^*B_x\preceq\kappa^{-1}D^2W(x)
\preceq\kappa^{-1}M_x.
\]

Therefore

\[
B_xM_x^{-1}B_x^*\preceq\kappa^{-1}I
\quad\text{on }L^2(\eta_x;E),
\tag{2.23}
\]

because
\(M_x^{-1/2}B_x^*B_xM_x^{-1/2}\preceq\kappa^{-1}I\)
and \(CC^*\) and \(C^*C\) have the same nonzero operator norm. Hence

\[
\int\langle M_x^{-1}r,r\rangle d\nu
\le\kappa^{-1}\int|\nabla_Ef|^2d\mu.
\tag{2.24}
\]

Minkowski's inequality in the Hilbert space
\(L^2(\nu;F,M_x^{-1})\), followed by the two-dimensional Cauchy--Schwarz
inequality, now gives

\[
\begin{aligned}
\left(\int|a+r|_{M_x^{-1}}^2d\nu\right)^{1/2}
&\le
\sqrt{c}\left(\int|\nabla_Ff|^2d\mu\right)^{1/2}
+\frac1{\sqrt\kappa}
\left(\int|\nabla_Ef|^2d\mu\right)^{1/2}\\
&\le
\sqrt{c+\kappa^{-1}}
\left(\int|\nabla f|^2d\mu\right)^{1/2}.
\end{aligned}
\tag{2.25}
\]

Combining (2.21) and (2.25),

\[
\operatorname{Var}_\nu g
\le2(c+\kappa^{-1})\int|\nabla f|^2d\mu.
\tag{2.26}
\]

Finally, (2.5) gives

\[
\operatorname{Var}_\mu f
\le
\left(2c+3\kappa^{-1}\right)
\int|\nabla f|^2d\mu.
\tag{2.27}
\]

This proves (L2). Notice that separately applying
\(|a+r|^2\le2|a|^2+2|r|^2\) would produce the weaker constants
\(4c+5/\kappa\). The Hilbert-space Minkowski step (2.25) is necessary for
the claimed \(2c+3/\kappa\).

If \(F=\{0\}\), conditional Lichnerowicz directly yields
\(C_P(\mu)\le\kappa^{-1}\), which is stronger than (L2). If \(E=\{0\}\),
then \(\mu=\nu\) and (L2) is tautological.

---

## 3. Audit and proof of (L3)

### 3.1 Preservation of the partial curvature

Condition (2.2) is equivalent, in the classical or distributional sense, to
convexity of

\[
V_0(x)=V(x)-\frac\kappa2|P_Ex|^2.
\tag{3.1}
\]

Let \(p_\varepsilon\) be the density of
\(\mu_\varepsilon=\mu*\gamma_\varepsilon\). Completing the square in the
\(E\)-variables gives, for \(z=z_F+z_E\),

\[
\begin{aligned}
&\frac\kappa2|x_E|^2+
\frac1{2\varepsilon}|z_E-x_E|^2\\
&\qquad=
\frac{\kappa}{2(1+\kappa\varepsilon)}|z_E|^2
+\frac{1+\kappa\varepsilon}{2\varepsilon}
\left|x_E-\frac{z_E}{1+\kappa\varepsilon}\right|^2.
\end{aligned}
\tag{3.2}
\]

Consequently

\[
p_\varepsilon(z)=
\exp\left[-\frac{\kappa}{2(1+\kappa\varepsilon)}|z_E|^2\right]
L_\varepsilon(z),
\tag{3.3}
\]

where, up to a constant factor,

\[
\begin{aligned}
L_\varepsilon(z)=\int
\exp\bigg[&-V_0(x)
-\frac{|z_F-x_F|^2}{2\varepsilon}\\
&-\frac{1+\kappa\varepsilon}{2\varepsilon}
\left|x_E-\frac{z_E}{1+\kappa\varepsilon}\right|^2
\bigg]dx.
\end{aligned}
\tag{3.4}
\]

The negative logarithm of the integrand in (3.4) is jointly convex in
\((x,z)\). Prekopa's theorem therefore makes \(L_\varepsilon\) log-concave.
Writing \(p_\varepsilon=e^{-V_\varepsilon}\), (3.3) implies

\[
D^2V_\varepsilon
\succeq
\frac\kappa{1+\kappa\varepsilon}P_E.
\tag{3.5}
\]

This proof works without differentiability of \(V\): (3.1) is then understood
as convexity on the affine support, and (3.5) distributionally. The
convolution itself is smooth and positive.

### 3.2 Poincare constant of the weak marginal

Projection commutes with convolution:

\[
\pi_F{}_{\#}(\mu*\gamma_\varepsilon)
=(\pi_F{}_{\#}\mu)*\gamma_\varepsilon^F.
\tag{3.6}
\]

Let \(X\sim\nu=\pi_F{}_{\#}\mu\) and
\(G\sim N(0,\varepsilon I_F)\) be independent. For a locally Lipschitz
\(f\in L^2(\nu*\gamma_\varepsilon^F)\), conditional variance and Jensen give

\[
\begin{aligned}
\operatorname{Var}f(X+G)
&=\mathbb E\operatorname{Var}(f(X+G)\mid X)
+\operatorname{Var}(\mathbb E[f(X+G)\mid X])\\
&\le\varepsilon\mathbb E|\nabla f(X+G)|^2
+C_P(\nu)\mathbb E
|\mathbb E[\nabla f(X+G)\mid X]|^2\\
&\le\bigl(C_P(\nu)+\varepsilon\bigr)
\mathbb E|\nabla f(X+G)|^2.
\end{aligned}
\tag{3.7}
\]

Thus

\[
C_P(\pi_F{}_{\#}\mu_\varepsilon)
\le C_P(\pi_F{}_{\#}\mu)+\varepsilon.
\tag{3.8}
\]

This is the weak Poincare constant: no smoothness or full-dimensionality of
the original marginal is required. A point marginal has constant \(0\), and
(3.8) then reduces to the Gaussian constant \(\varepsilon\).

---

## 4. Approximation, domains, and lower-dimensional supports

### 4.1 Removal of smoothness in (L2)

Suppose intrinsically on \(H=F\oplus E\) that

\[
V-\frac\kappa2|P_E\cdot|^2
\quad\text{is convex},
\tag{4.1}
\]

with \(V\) possibly nonsmooth or extended-valued. Let
\(\mu_\varepsilon=\mu*\gamma_\varepsilon\), with convolution taken in \(H\).
By (L3),

\[
\kappa_\varepsilon=
\frac\kappa{1+\kappa\varepsilon},
\qquad
C_P(\pi_F{}_{\#}\mu_\varepsilon)
\le c+\varepsilon,
\tag{4.2}
\]

where \(c=C_P(\pi_F{}_{\#}\mu)\). Applying the smooth result (2.27) gives

\[
C_P(\mu_\varepsilon)
\le2(c+\varepsilon)+\frac3{\kappa_\varepsilon}
=2c+\frac3\kappa+5\varepsilon.
\tag{4.3}
\]

For \(f\in C_c^\infty(H)\), weak convergence
\(\mu_\varepsilon\Rightarrow\mu\) gives

\[
\operatorname{Var}_{\mu_\varepsilon}f
\to\operatorname{Var}_\mu f,
\qquad
\int|\nabla f|^2d\mu_\varepsilon
\to\int|\nabla f|^2d\mu.
\tag{4.4}
\]

Letting \(\varepsilon\downarrow0\) in (4.3) proves (2.27) on the smooth compact
core. Value truncation, spatial cutoffs, and convolution of the test
function inside \(H\) extend it to all locally Lipschitz
\(f\in L^2(\mu)\). For the cutoff step, choose \(|\nabla\chi_R|\le2/R\);
the extra energy is bounded by a constant times
\(R^{-2}\int_{|x|\ge R}|f|^2d\mu\), which tends to zero. Thus no compactness
limit changes the constant.

### 4.2 Intrinsic support convention

Let the affine support of a log-concave measure be \(z_0+S\). Translate by
\(-z_0\) and identify \(S\) with \(\mathbb R^{\dim S}\). Then:

* all gradients and Hessians in (L1)--(L2) are intrinsic to \(S\);
* for (L2), the relevant decomposition is an orthogonal decomposition
  \(S=F_S\oplus E_S\), and \(P_{E_S}\) replaces the ambient \(P_E\);
* densities are with respect to Lebesgue (Hausdorff) measure on \(S\);
* Gaussian convolution used for approximation is first taken inside \(S\).

With these conventions, every proof above is unchanged. If \(S\) is a
point, \(C_P=0\); if only the \(F_S\)-marginal is a point, the \(F=0\) case in
Section 2.4 applies. Normal ambient directions contribute neither variance
nor intrinsic Dirichlet energy.

Without this convention, a phrase such as
\(D^2V\succeq\kappa P_E\) for a measure supported on a proper affine subspace
has no invariant meaning: values of an ambient extension of \(V\) normal to
the support can be altered arbitrarily without changing the measure. This is
a statement-level omission, not a counterexample to the intrinsic theorem.
For an explicit illustration, let \(\mu\) be standard Gaussian measure on the
\(x\)-axis in \(\mathbb R^2\). The two ambient functions
\[
V_+(x,y)=x^2/2+\kappa y^2/2,
\qquad
V_-(x,y)=x^2/2-My^2
\]
have the same restriction \(x^2/2\) to the support and hence describe the
same intrinsic density there, while their normal Hessians have opposite
signs. Thus an ambient curvature assertion cannot be recovered from this
lower-dimensional measure; the affine-hull convention is necessary.

---

## 5. Audit checklist

* **One-form domain:** handled through the closed gradient, its polar
  decomposition, and the Friedrichs one-form form; no eigenfunction or
  attained spectral gap is assumed.
* **Closed gradient subspace:** it is reducing by the scalar/one-form
  intertwining. The proof never assumes that pointwise multiplication by
  \(D^2U+\lambda I\) preserves gradients.
* **Inverse comparison:** the exact variational supremum gives the factor
  \(2\). Compression and inversion are not interchanged as equalities.
* **Conditional sign:** \(w_{x,u}=-\nabla A_x^{-1}q_{x,u}\). This sign is
  forced by the continuity equation and produces (2.10).
* **Dimension tracking:** the only numerical constants are \(2\) from (L1),
  \(1/\kappa\) from conditional Lichnerowicz, and the two-dimensional
  Cauchy--Schwarz step. No trace, rank, or dimension occurs.
* **Lower-dimensional support:** the result is intrinsic to the affine hull;
  ambient normal Hessians are inadmissible.

No counterexample to (L1), (L2), or (L3) remains after these domain and support
conventions are imposed.
