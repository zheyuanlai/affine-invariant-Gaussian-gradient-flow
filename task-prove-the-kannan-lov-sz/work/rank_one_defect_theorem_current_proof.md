# The Rank-One-Defect Poincare Theorem

## A self-contained account of the current proof

> **Status.** This document records the strongest current proof in the
> workspace. It is a partial result toward KLS, not a proof of KLS. The proof
> is mathematically self-contained modulo the standard results listed in
> Section 2, but it has not yet undergone external peer review or a complete
> literature-novelty audit.

### Contents

1. Main theorem and its defect-subspace extension
2. Standard inputs
3. Conditional velocity and reinforced Prekopa
4. Improved marginal Brascamp--Lieb estimate
5. Assembly of the smooth proof
6. Nonsmooth potentials, limits, and affine support
7. Alternative scalar weighted-Hardy proof
8. Adversarial examples and scaling checks
9. Exact relationship to localization and KLS
10. Constant ledger and publication audit

Throughout, all Euclidean decompositions are orthogonal. For a probability
measure \(\mu\), its Poincare constant is

\[
 C_P(\mu)=\inf\left\{C\ge 0:
 \operatorname{Var}_\mu f\le C\int |\nabla f|^2\,d\mu
 \text{ for every locally Lipschitz }f\in L^2(\mu)\right\}.
\]

The convention is \(C_P(\delta_x)=0\).

---

## 1. Main statements

### 1.1 Rank-one-defect theorem

Let \(u\in S^{n-1}\), and let

\[
 P_{u^\perp}=I-u\otimes u.
\]

Let \(\mu\) be a log-concave probability measure on \(\mathbb R^n\) with
full-dimensional density

\[
 d\mu(x)=Z^{-1}e^{-V(x)}\,dx,
\]

where \(V:\mathbb R^n\to(-\infty,+\infty]\) is proper, lower
semicontinuous, and convex. Suppose, in the sense of distributions, that

\[
 D^2V\succeq \kappa P_{u^\perp}                         \tag{1.1}
\]

for some \(\kappa>0\). Equivalently,

\[
 x\longmapsto V(x)-\frac{\kappa}{2}|P_{u^\perp}x|^2
\]

is convex. Then

\[
 \boxed{
 C_P(\mu)
 \le 3\kappa^{-1}
      +96\operatorname{Var}_\mu\langle X,u\rangle
 \le 96\left(
      \kappa^{-1}+\operatorname{Var}_\mu\langle X,u\rangle
      \right).}                                      \tag{1.2}
\]

The constant is independent of \(n\), \(\mu\), \(u\), and \(\kappa\).
No centering assumption is required. In particular, if
\(\operatorname{Cov}(\mu)\preceq I\), then

\[
 C_P(\mu)\le 96(1+\kappa^{-1}).                       \tag{1.3}
\]

The theorem says that uniform convexity in all but one direction is enough:
the uncontrolled direction costs only its actual one-dimensional variance.

### 1.2 Strong-conditionals/weak-marginal theorem

The rank-one result is a corollary of a more general statement.

Let

\[
 \mathbb R^n=F\oplus E,
\]

let \(P_E\) be the projection onto \(E\), and let

\[
 \bar\mu=(P_F)_\#\mu
\]

be the marginal on \(F\). If

\[
 D^2V\succeq \kappa P_E                              \tag{1.4}
\]

distributionally, then

\[
 \boxed{
 C_P(\mu)
 \le \kappa^{-1}
      +2\bigl(C_P(\bar\mu)+\kappa^{-1}\bigr)
 =2C_P(\bar\mu)+3\kappa^{-1}.}                       \tag{1.5}
\]

No dimension-free estimate for the weak marginal is hidden in (1.5): it
uses its **actual** Poincare constant. When \(\dim F=1\), the standard
one-dimensional log-concave estimate

\[
 C_P(\bar\mu)le
 48\operatorname{Var}_{\bar\mu}(T)                   \tag{1.6}
\]

turns (1.5) into (1.2).

### 1.3 Intrinsic affine-support formulation

If \(\mu\) is supported on a proper affine subspace \(x_0+H\), translate
to \(H\) and formulate (1.4) intrinsically there:

\[
 H=F_H\oplus E_H,
 \qquad D_H^2V\succeq\kappa P_{E_H}.
\]

Then (1.5) holds with intrinsic gradients, covariance, convolution, and
Lebesgue measure on \(H\). For the rank-one statement, the weak unit vector
should be chosen in \(H\). This avoids ambiguity caused by restricting an
ambient projection \(I-u\otimes u\) when \(u\notin H\). A point mass is
handled by the convention \(C_P=0\).

---

## 2. Standard ingredients used

The proof uses the following established facts.

1. **Prekopa's theorem.** A marginal of a log-concave density is
   log-concave.
2. **Lichnerowicz/Bakry--Emery.** If \(D^2W\succeq\kappa I\), then
   \(C_P(e^{-W})\le\kappa^{-1}\).
3. **Bochner identity.** For
   \(L=\Delta-\nabla W\cdot\nabla\),
   \[
   \int(L\phi)^2d\nu
   =\int\|D^2\phi\|_{\mathrm{HS}}^2d\nu
    +\int\langle D^2W\nabla\phi,\nabla\phi\rangle d\nu.
   \]
4. **Lax--Milgram.** A centered conditional score has a unique centered
   weak Poisson solution under a conditional spectral gap.
5. **Total variance.** For a disintegration \(\mu(dz,dy)=\bar\mu(dz)
   \nu_z(dy)\),
   \[
   \operatorname{Var}_\mu f
   =\int\operatorname{Var}_{\nu_z}f\,d\bar\mu(z)
    +\operatorname{Var}_{\bar\mu}(\mathbb E_zf).
   \]
6. **One-dimensional log-concave Poincare inequality.** The deliberately
   nonoptimized constant (1.6) is sufficient.
7. **Weighted one-dimensional Muckenhoupt/Hardy criterion.** This is used
   only in the alternative rank-one proof in Section 7.
8. **Brascamp--Lieb for linear functions.** It is used in the Gaussian
   regularization step.

All other estimates needed for the theorem are proved below.

---

## 3. Smooth proof of the defect-subspace theorem

We first assume that \(V\in C^\infty(\mathbb R^n)\), all derivatives used
below are integrable, and a small global quadratic confinement has been
added. Section 6 removes these assumptions with uniform constants.

Write

\[
 x=z+y,\qquad z\in F,\quad y\in E.
\]

In these coordinates,

\[
 D^2V(z,y)=
 \begin{pmatrix}
  A&B^*\\
  B&C
 \end{pmatrix}.                                      \tag{3.1}
\]

Hypothesis (1.4) is exactly

\[
 \begin{pmatrix}
  A&B^*\\
  B&C-\kappa I_E
 \end{pmatrix}\succeq0.                              \tag{3.2}
\]

In particular,

\[
 C\succeq\kappa I_E.                                 \tag{3.3}
\]

### 3.1 Disintegration

Define the weak marginal potential \(U\) and the conditional measures by

\[
 e^{-U(z)}=\int_Ee^{-V(z,y)}\,dy,                     \tag{3.4}
\]

\[
 d\nu_z(y)=e^{U(z)-V(z,y)}\,dy.                       \tag{3.5}
\]

Thus

\[
 d\mu(z,y)=d\bar\mu(z)d\nu_z(y),
 \qquad d\bar\mu(z)=e^{-U(z)}dz.
\]

By Prekopa, \(U\) is convex. By (3.3), every \(\nu_z\) is
\(\kappa\)-strongly log-concave and

\[
 C_P(\nu_z)\le\kappa^{-1}.                            \tag{3.6}
\]

### 3.2 Conditional scores and canonical velocity fields

Fix \(\xi\in F\). Define the centered conditional score

\[
 s_\xi(y)
 =\langle\nabla_zV(z,y),\xi\rangle
  -\mathbb E_z\langle\nabla_zV(z,\cdot),\xi\rangle. \tag{3.7}
\]

Let

\[
 L_z=\Delta_y-\nabla_yV(z,\cdot)\cdot\nabla_y          \tag{3.8}
\]

be the conditional generator. Since \(s_\xi\) is centered and (3.6)
holds, there is a unique centered weak solution \(\phi_\xi\) of

\[
 L_z\phi_\xi=s_\xi.                                  \tag{3.9}
\]

Set

\[
 \mathcal W_z(y)\xi=\nabla_y\phi_\xi(y).              \tag{3.10}
\]

The map \(\xi\mapsto\mathcal W_z(y)\xi\) is linear. Define the positive
operator \(G(z):F\to F\) by

\[
 \langle G(z)\xi,\xi\rangle
 =\mathbb E_z|\mathcal W_z\xi|^2.                     \tag{3.11}
\]

Differentiating (3.5) in direction \(\xi\) gives

\[
 \partial_\xi\nu_z=-s_\xi\nu_z.                     \tag{3.12}
\]

Equation (3.9) gives

\[
 \operatorname{div}_y(\nu_z\mathcal W_z\xi)
 =\nu_zL_z\phi_\xi=s_\xi\nu_z.                      \tag{3.13}
\]

Therefore

\[
 \boxed{
 \partial_\xi\nu_z
 +\operatorname{div}_y(\nu_z\mathcal W_z\xi)=0.}     \tag{3.14}
\]

This is the canonical continuity equation for the moving conditional
law. It is the device that retains, rather than discards, the derivative of
the conditional measure.

### 3.3 Reinforced Prekopa inequality

The central estimate is

\[
 \boxed{D^2U(z)\succeq\kappa G(z).}                   \tag{3.15}
\]

We prove it direction by direction. Fix \(\xi\in F\), and abbreviate

\[
 a=\langle A\xi,\xi\rangle,
 \qquad b=B\xi,
 \qquad w=\mathcal W_z\xi,
 \qquad s=s_\xi.
\]

Marginal differentiation gives

\[
 \langle D^2U\xi,\xi\rangle
 =\mathbb E_za-\operatorname{Var}_z
   \langle\nabla_zV,\xi\rangle
 =\mathbb E_za-\mathbb E_zs^2.                       \tag{3.16}
\]

Bochner applied to \(L_z\phi_\xi=s\) yields

\[
 \mathbb E_zs^2
 =\mathbb E_z\|D_y^2\phi_\xi\|_{\mathrm{HS}}^2
  +\mathbb E_z\langle Cw,w\rangle.                   \tag{3.17}
\]

Since \(\nabla_ys=b\), integration by parts gives

\[
 \mathbb E_z\langle b,w\rangle
 =\mathbb E_z\langle\nabla_ys,\nabla_y\phi_\xi\rangle
 =-\mathbb E_zsL_z\phi_\xi
 =-\mathbb E_zs^2.                                   \tag{3.18}
\]

Equations (3.16) and (3.18) imply

\[
 \langle D^2U\xi,\xi\rangle
 =\mathbb E_z\bigl[a+\langle b,w\rangle\bigr].       \tag{3.19}
\]

Apply the full block inequality (3.2) pointwise to \((\xi,w)\):

\[
 a+2\langle b,w\rangle+\langle Cw,w\rangle
 \ge\kappa|w|^2.                                     \tag{3.20}
\]

From (3.17)--(3.18),

\[
 \mathbb E_z\langle Cw,w\rangle
 \le\mathbb E_zs^2
 =-\mathbb E_z\langle b,w\rangle.                   \tag{3.21}
\]

Consequently, the expectation of the left side of (3.20) is at most
the right side of (3.19). Hence

\[
 \langle D^2U\xi,\xi\rangle
 \ge\kappa\mathbb E_z|w|^2
 =\kappa\langle G\xi,\xi\rangle.
\]

Since this holds for every \(\xi\), (3.15) follows.

Two points are load-bearing:

* the cross Hessian block \(B\) is used in (3.20);
* merely knowing \(C\succeq\kappa I_E\) would not prove (3.15).

### 3.4 Differentiating the conditional expectation

Let \(f\in C_c^\infty(\mathbb R^n)\), and define

\[
 g(z)=\mathbb E_z f(z,\cdot).                          \tag{3.22}
\]

Using (3.12)--(3.14), for every \(\xi\in F\),

\[
\begin{aligned}
 \langle\nabla g,\xi\rangle
 &=\mathbb E_z\langle\nabla_zf,\xi\rangle
   -\mathbb E_z(fs_\xi)\\
 &=\mathbb E_z\left[
   \langle\nabla_zf,\xi\rangle
   +\langle\nabla_yf,\mathcal W_z\xi\rangle
   \right].
\end{aligned}                                        \tag{3.23}
\]

Define

\[
 \mathcal B_z(y):F\to F\oplus E,
 \qquad
 \mathcal B_z(y)\xi=(\xi,\mathcal W_z(y)\xi).        \tag{3.24}
\]

If \(q=(\nabla_zf,\nabla_yf)\), then

\[
 \nabla g=\mathbb E_z\mathcal B_z^*q,
 \qquad
 \mathbb E_z\mathcal B_z^*\mathcal B_z=I_F+G(z).     \tag{3.25}
\]

Matrix Cauchy--Schwarz gives

\[
 \boxed{
 \langle(I_F+G)^{-1}\nabla g,\nabla g\rangle
 \le\mathbb E_z|\nabla f|^2.}                        \tag{3.26}
\]

For completeness, put \(\alpha=(I_F+G)^{-1}\nabla g\). Then

\[
 \langle\alpha,\nabla g\rangle
 =\mathbb E_z\langle\mathcal B_z\alpha,q\rangle.
\]

Cauchy--Schwarz implies

\[
 \langle\alpha,\nabla g\rangle^2
 \le
 \langle\alpha,(I_F+G)\alpha\rangle
 \mathbb E_z|q|^2.
\]

Since
\(\langle\alpha,\nabla g\rangle
=\langle\alpha,(I_F+G)\alpha\rangle\), cancellation proves (3.26).

---

## 4. The improved marginal Brascamp--Lieb estimate

We need an estimate that combines the local Hessian of the weak marginal
with its actual global spectral gap.

### Lemma 4.1

Let

\[
 d\eta(z)=e^{-U(z)}dz
\]

be a smooth log-concave probability on \(\mathbb R^k\). Put

\[
 H=D^2U\succeq0,
 \qquad C_0=C_P(\eta),
 \qquad \lambda=C_0^{-1}.
\]

Assume \(0<C_0<\infty\). Then every smooth \(h\in L^2(\eta)\) satisfies

\[
 \boxed{
 \operatorname{Var}_\eta h
 \le2\int
 \langle(H+\lambda I)^{-1}\nabla h,\nabla h\rangle
 \,d\eta.}                                           \tag{4.1}
\]

### Proof

Let

\[
 L=\Delta-\nabla U\cdot\nabla,
 \qquad A_0=-L                                      \tag{4.2}
\]

on centered scalar functions. Let \(\mathcal G\) be the
\(L^2(\eta;\mathbb R^k)\)-closure of gradients. On one-forms define

\[
 A_1=-L\otimes I+H.                                  \tag{4.3}
\]

On a smooth core, and then by closure,

\[
 \nabla A_0=A_1\nabla.                               \tag{4.4}
\]

For \(r=\nabla a\in\mathcal G\), the scalar spectral gap gives

\[
\begin{aligned}
 \langle r,A_1r\rangle
 &=\|A_0a\|_2^2\\
 &\ge\lambda\langle a,A_0a\rangle\\
 &=\lambda\|r\|_2^2.
\end{aligned}                                        \tag{4.5}
\]

The one-form energy identity also gives

\[
 \langle r,A_1r\rangle
 =\int\left(\|\nabla r\|_{\mathrm{HS}}^2
             +\langle Hr,r\rangle\right)d\eta
 \ge\int\langle Hr,r\rangle d\eta.                 \tag{4.6}
\]

Combining (4.5) and (4.6), as quadratic forms on \(\mathcal G\),

\[
 A_1\succeq
 \frac12P_{\mathcal G}(H+\lambda I)P_{\mathcal G}.  \tag{4.7}
\]

For centered \(h\), set \(a=A_0^{-1}h\). The intertwining relation gives

\[
 A_1\nabla a=\nabla h.                               \tag{4.8}
\]

Integration by parts yields the Helffer--Sjostrand identity

\[
 \operatorname{Var}_\eta h
 =\langle\nabla h,A_1^{-1}\nabla h\rangle_{\mathcal G}. \tag{4.9}
\]

Using the variational formula for the inverse and enlarging the supremum
from \(\mathcal G\) to all vector fields,

\[
\begin{aligned}
 \operatorname{Var}_\eta h
 &=\sup_{r\in\mathcal G}
   \left\{2\langle\nabla h,r\rangle
          -\langle r,A_1r\rangle\right\}\\
 &\le\sup_{r\in\mathcal G}
   \left\{2\langle\nabla h,r\rangle
          -\frac12\langle r,(H+\lambda I)r\rangle\right\}\\
 &\le2\int
   \langle(H+\lambda I)^{-1}\nabla h,\nabla h\rangle
   d\eta.
\end{aligned}                                        \tag{4.10}
\]

This proves (4.1). The factor \(2\) comes solely from averaging the two
lower bounds (4.5) and (4.6).

### Why ordinary Brascamp--Lieb is insufficient

Ordinary Brascamp--Lieb would use the weight \(H^{-1}\), which is
singular on flat directions of \(U\). Lemma 4.1 inserts the genuine
global gap \(\lambda\) and replaces that weight by
\((H+\lambda I)^{-1}\). No dimension-free gap for \(\eta\) is assumed.

---

## 5. Completion of the smooth proof

Apply Lemma 4.1 to \(\eta=\bar\mu\) and \(h=g\). Thus

\[
 H=D^2U,
 \qquad C_0=C_P(\bar\mu),
 \qquad \lambda=C_0^{-1}.                             \tag{5.1}
\]

From (3.15),

\[
 H\succeq\kappa G.                                   \tag{5.2}
\]

Therefore

\[
 (I_F+G)^{-1}
 \succeq\kappa(\kappa I_F+H)^{-1}.                   \tag{5.3}
\]

For every scalar \(r\ge0\),

\[
 \frac1{r+\lambda}
 \le
 \left(\frac1\lambda+\frac1\kappa\right)
 \frac\kappa{\kappa+r}.                              \tag{5.4}
\]

Indeed, after multiplying by the positive denominators, (5.4) is
equivalent to

\[
 \lambda(\kappa+r)
 \le(\kappa+\lambda)(r+\lambda),
\]

whose right side minus left side is \(\kappa r+\lambda^2\ge0\).
Functional calculus and (5.3) give

\[
 (H+\lambda I)^{-1}
 \preceq
 (C_0+\kappa^{-1})(I_F+G)^{-1}.                       \tag{5.5}
\]

Lemma 4.1, (5.5), and (3.26) imply

\[
\begin{aligned}
 \operatorname{Var}_{\bar\mu}g
 &\le2(C_0+\kappa^{-1})
 \int\langle(I_F+G)^{-1}\nabla g,\nabla g\rangle
 \,d\bar\mu\\
 &\le2(C_0+\kappa^{-1})
 \int|\nabla f|^2d\mu.                               \tag{5.6}
\end{aligned}
\]

For the conditional term, (3.6) gives

\[
 \int_F\operatorname{Var}_{\nu_z}f\,d\bar\mu(z)
 \le\kappa^{-1}\int|\nabla_yf|^2d\mu
 \le\kappa^{-1}\int|\nabla f|^2d\mu.               \tag{5.7}
\]

The total-variance identity now yields

\[
\begin{aligned}
 \operatorname{Var}_\mu f
 &\le
 \left[\kappa^{-1}
 +2(C_0+\kappa^{-1})\right]
 \int|\nabla f|^2d\mu\\
 &=\left(2C_0+3\kappa^{-1}\right)
 \int|\nabla f|^2d\mu.                               \tag{5.8}
\end{aligned}
\]

This proves (1.5) in the smooth setting.

If \(\dim F=1\), (1.6) gives

\[
 C_P(\mu)
 \le96\operatorname{Var}_\mu\langle X,u\rangle
      +3\kappa^{-1}
 \le96\left(
 \operatorname{Var}_\mu\langle X,u\rangle
 +\kappa^{-1}\right),                                \tag{5.9}
\]

which is the rank-one theorem.

---

## 6. Nonsmooth potentials, unbounded support, and test functions

We now remove smoothness and confinement without changing the constants.

### 6.1 Gaussian regularization and retained curvature

Assume that \(V\) is merely proper, lower semicontinuous, convex, and

\[
 V(x)-\frac\kappa2|P_Ex|^2
\]

is convex. Let

\[
 \gamma_\varepsilon=N(0,\varepsilon I),
 \qquad
 \mu_\varepsilon=\mu*\gamma_\varepsilon
 =e^{-V_\varepsilon}dx.                               \tag{6.1}
\]

Then \(V_\varepsilon\) is smooth and

\[
 \boxed{
 D^2V_\varepsilon
 \succeq\kappa_\varepsilon P_E,
 \qquad
 \kappa_\varepsilon
 =\frac\kappa{1+\kappa\varepsilon}.}                 \tag{6.2}
\]

To prove this, let \(X\sim\mu\), \(G\sim N(0,\varepsilon I)\), and
condition on \(X+G=z\). The posterior potential is

\[
 V(x)+\frac{|z-x|^2}{2\varepsilon},                   \tag{6.3}
\]

whose distributional Hessian is at least

\[
 \kappa P_E+\varepsilon^{-1}I.
\]

Brascamp--Lieb for linear functions gives

\[
 \operatorname{Cov}(X\mid X+G=z)
 \preceq(\kappa P_E+\varepsilon^{-1}I)^{-1}.          \tag{6.4}
\]

The Gaussian-channel Hessian identity is

\[
 D^2V_\varepsilon(z)
 =\varepsilon^{-1}I
  -\varepsilon^{-2}
   \operatorname{Cov}(X\mid X+G=z).                   \tag{6.5}
\]

On \(F\), (6.5) gives the nonnegative lower bound zero. On \(E\), it
gives

\[
 \varepsilon^{-1}
 -\varepsilon^{-2}(\kappa+\varepsilon^{-1})^{-1}
 =\frac\kappa{1+\kappa\varepsilon},
\]

proving (6.2).

### 6.2 The weak marginal under convolution

The weak marginal becomes

\[
 \bar\mu_\varepsilon
 =\bar\mu*N(0,\varepsilon I_F).                       \tag{6.6}
\]

Tensorization with the independent Gaussian and Jensen's inequality give

\[
 C_P(\bar\mu_\varepsilon)
 \le C_P(\bar\mu)+\varepsilon.                        \tag{6.7}
\]

Indeed, apply Poincare first in \(Z\sim\bar\mu\) and then in
\(G_F\sim N(0,\varepsilon I_F)\) to \(h(Z+G_F)\). Both resulting
gradient terms are bounded by the Dirichlet integral under the convolution.

Applying the smooth theorem to \(\mu_\varepsilon\) gives

\[
 C_P(\mu_\varepsilon)
 \le2\bigl(C_P(\bar\mu)+\varepsilon\bigr)
    +3\bigl(\kappa^{-1}+\varepsilon\bigr).            \tag{6.8}
\]

Thus the error is \(5\varepsilon\), independent of dimension.

### 6.3 Temporary confinement

If additional decay is needed to justify differentiation and integration
by parts, replace \(V_\varepsilon\) by

\[
 V_{\varepsilon,\delta}(x)
 =V_\varepsilon(x)+\frac\delta2|x|^2.                 \tag{6.9}
\]

The curvature lower bound is not weakened. Apply the smooth proof and
then let \(\delta\downarrow0\). For compactly supported test functions,
the variance, weak marginal moments, and Dirichlet integral converge by
dominated convergence.

### 6.4 Passage to general locally Lipschitz test functions

Let \(f\in L^2(\mu)\) be locally Lipschitz with
\(\int|\nabla f|^2d\mu<\infty\).

First truncate values:

\[
 T_M(r)=\max(-M,\min(r,M)).                            \tag{6.10}
\]

Then

\[
 T_M(f)\to f\quad\text{in }L^2(\mu),
 \qquad
 |\nabla T_M(f)|
 =1_{\{|f|<M\}}|\nabla f|\quad\mu\text{-a.e.}       \tag{6.11}
\]

Next choose \(\chi_R\) equal to one on \(B_R\), zero outside
\(B_{2R}\), with \(|\nabla\chi_R|\le2/R\). For fixed \(M\),

\[
 \chi_RT_M(f)\to T_M(f)\quad\text{in }L^2(\mu),      \tag{6.12}
\]

and

\[
 \int|\nabla(\chi_RT_M(f))|^2d\mu
 \longrightarrow
 \int|\nabla T_M(f)|^2d\mu.                          \tag{6.13}
\]

The cutoff-gradient term has \(L^2\)-norm at most \(2M/R\); the main
term converges by dominated convergence.

Finally, mollify the compactly supported Lipschitz function. Uniform
convergence and almost-everywhere convergence of gradients give convergence
of both variance and energy because a full-dimensional log-concave measure
is absolutely continuous.

The logically safe order is the following.

1. For a fixed smooth compactly supported test function, first remove the
   temporary confinement \(\delta\downarrow0\), and then let the Gaussian
   regularization scale \(\varepsilon\downarrow0\).
2. Having obtained the inequality for \(C_c^\infty\) functions under the
   original measure, extend the test class in the order
   \[
   \text{mollification scale}\downarrow0,
   \qquad R\uparrow\infty,
   \qquad M\uparrow\infty.                             \tag{6.14}
   \]

In particular, no convergence assertion for the numerical constants
\(C_P(\mu_\varepsilon)\) themselves is needed: one passes to the limit in
the inequality for each fixed test function. This proves (1.5) for all
locally Lipschitz \(L^2\) test functions. If the Dirichlet integral is
infinite, the desired inequality is automatic.

### 6.5 Lower-dimensional support

All steps above are performed intrinsically on the affine hull:

* Lebesgue measure is the Lebesgue measure on the supporting subspace;
* gradients and Hessians are intrinsic;
* Gaussian convolution uses a Gaussian on that subspace;
* the orthogonal decomposition is a decomposition of its direction space.

Hence no artificial ambient variance or covariance is introduced.

---

## 7. Alternative direct rank-one proof

The proof above is the cleanest route and gives constant \(96\). There is
also a direct scalar argument that gives constant \(201\). It is useful
because it makes the compensation mechanism completely explicit.

Use coordinates

\[
 x=tu+y,
 \qquad t\in\mathbb R,
 \quad y\in u^\perp.
\]

Let

\[
 e^{-U(t)}=\int_{u^\perp}e^{-V(t,y)}dy,
 \qquad d\nu_t(y)=e^{U(t)-V(t,y)}dy.                  \tag{7.1}
\]

The scalar version of Section 3 constructs

\[
 w_t=\nabla_y\phi_t,
 \qquad L_t\phi_t=V_t-\mathbb E_tV_t,                \tag{7.2}
\]

and proves

\[
 \partial_t\nu_t+\operatorname{div}_y(\nu_tw_t)=0,   \tag{7.3}
\]

\[
 \boxed{
 \kappa W(t)\le U''(t),
 \qquad W(t)=\mathbb E_t|w_t|^2.}                    \tag{7.4}
\]

For

\[
 g(t)=\mathbb E_tf(t,\cdot),                          \tag{7.5}
\]

one has

\[
 g'(t)=\mathbb E_t
 \bigl(\partial_tf+\nabla_yf\cdot w_t\bigr),         \tag{7.6}
\]

and therefore

\[
 \boxed{
 \frac{g'(t)^2}{1+W(t)}
 \le\mathbb E_t|\nabla f|^2.}                        \tag{7.7}
\]

The remaining input is the following weighted one-dimensional inequality.

### Lemma 7.1

Let

\[
 d\nu(t)=\rho(t)dt=Z^{-1}e^{-U(t)}1_{(\alpha,\beta)}dt
\]

be a one-dimensional log-concave probability of variance \(s^2\), with
\(U\in C^2\). Suppose \(W\ge0\) and

\[
 \kappa W\le U''.
\]

Then

\[
 \operatorname{Var}_\nu h
 \le200(s^2+\kappa^{-1})
 \int\frac{h'^2}{1+W}d\nu.                            \tag{7.8}
\]

### Proof of Lemma 7.1

Set

\[
 a_0(t)=\frac\kappa{\kappa+U''(t)}.                   \tag{7.9}
\]

Since \(\kappa W\le U''\),

\[
 \frac1{1+W}\ge a_0.                                 \tag{7.10}
\]

Let \(m\) be a median. The weighted Muckenhoupt criterion bounds the
optimal constant in

\[
 \operatorname{Var}_\nu h\le C_a\int a_0h'^2d\nu
\]

by \(4\max(B_+,B_-)\), where

\[
 B_+
 =\sup_{m<x<\beta}\nu([x,\beta))
   \int_m^x\frac{dt}{a_0(t)\rho(t)},                  \tag{7.11}
\]

\[
 B_-
 =\sup_{\alpha<x<m}\nu((\alpha,x])
   \int_x^m\frac{dt}{a_0(t)\rho(t)}.                 \tag{7.12}
\]

The identity

\[
 \frac1{a_0\rho}
 =\frac1\rho+\frac1\kappa\frac{U''}{\rho}            \tag{7.13}
\]

separates the ordinary and curvature contributions. The ordinary part is
at most

\[
 C_P(\nu)\le48s^2.                                    \tag{7.14}
\]

For the right curvature term, put \(Q(x)=\nu([x,\beta))\). Integration
by parts gives

\[
 \int_m^x\frac{U''(t)}{\rho(t)}dt
 \le\frac{U'(x)}{\rho(x)}-rac{U'(m)}{\rho(m)}.      \tag{7.15}
\]

If \(U'(x)>0\), convexity gives

\[
 Q(x)\le\frac{\rho(x)}{U'(x)}.                        \tag{7.16}
\]

If \(U'(x)\le0\), the first term in (7.15) is
nonpositive. At the median, if \(U'(m)<0\), the left tangent estimate
gives

\[
 \frac{|U'(m)|}{\rho(m)}\le2,                         \tag{7.17}
\]

and \(Q(x)\le1/2\). If \(U'(m)\ge0\), its contribution has the favorable
sign. Consequently

\[
 \sup_{x>m}Q(x)
 \int_m^x\frac{U''(t)}{\rho(t)}dt\le2.                \tag{7.18}
\]

The reflected argument gives the same bound on the left. Hence

\[
 B_\pm\le48s^2+2\kappa^{-1}.                          \tag{7.19}
\]

The Muckenhoupt factor \(4\), followed by harmless rounding, proves
(7.8).

Finite endpoints cause no problem: tangent-line bounds are first integrated
only to the actual endpoint and then bounded by the corresponding infinite
exponential integral.

### Completion of the scalar route

Total variance gives

\[
 \operatorname{Var}_\mu f
 =\int\operatorname{Var}_{\nu_t}f\,d\bar\mu(t)
  +\operatorname{Var}_{\bar\mu}g.                     \tag{7.20}
\]

Conditional Lichnerowicz gives

\[
 \int\operatorname{Var}_{\nu_t}f\,d\bar\mu(t)
 \le\kappa^{-1}\int|\nabla_yf|^2d\mu.               \tag{7.21}
\]

Lemma 7.1, (7.4), and (7.7) give

\[
 \operatorname{Var}_{\bar\mu}g
 \le200(s^2+\kappa^{-1})
 \int|\nabla f|^2d\mu.                               \tag{7.22}
\]

Thus

\[
 C_P(\mu)\le201(s^2+\kappa^{-1}).                    \tag{7.23}
\]

The general defect-subspace proof improves \(201\) to \(96\) in the
rank-one case.

---

## 8. Sanity checks and adversarial models

### 8.1 Sheared Gaussian conditionals

Consider

\[
 V(t,y)=\frac c2t^2+\frac L2(y-mt)^2,
 \qquad L>\kappa.                                    \tag{8.1}
\]

The Hessian blocks are

\[
 a=c+Lm^2,
 \qquad b=-Lm,
 \qquad C=L.                                         \tag{8.2}
\]

The defect-curvature condition is equivalent to

\[
 c\ge\frac{\kappa L}{L-\kappa}m^2.                   \tag{8.3}
\]

The conditional law is \(N(mt,L^{-1})\), and the canonical velocity is

\[
 w_t=m.                                               \tag{8.4}
\]

The reinforced Prekopa estimate becomes

\[
 \kappa m^2\le c,                                    \tag{8.5}
\]

which is implied by (8.3). Thus a large shear cannot be inserted freely:
it necessarily creates curvature in the weak marginal.

For an affine test function \(f(t,y)=\alpha t+\beta y\), (7.7) becomes

\[
 \frac{(\alpha+\beta m)^2}{1+m^2}
 \le\alpha^2+\beta^2,                                \tag{8.6}
\]

and equality occurs when \((\alpha,\beta)\) is parallel to \((1,m)\).
This confirms the normalization and signs in the continuity equation.

### 8.2 Rapidly moving conditional centers

Let

\[
 m(t)=\int y\,d\nu_t(y).
\]

The continuity equation gives

\[
 m'(t)=\int w_t\,d\nu_t.                              \tag{8.7}
\]

Therefore

\[
 \kappa|m'(t)|^2
 \le\kappa\int|w_t|^2d\nu_t
 \le U''(t).                                         \tag{8.8}
\]

Any rapid conditional translation produces a curvature spike in the
marginal. Lemma 7.1 shows that this spike costs only \(O(\kappa^{-1})\)
in the final Poincare estimate.

### 8.3 Products

If

\[
 V(z,y)=U(z)+W(y),
 \qquad D^2W\succeq\kappa I_E,
\]

then all conditional velocities vanish, \(G=0\), and the theorem reduces
to a nonoptimized form of ordinary tensorization. The interest of the
proof is that it also handles nonproduct dependence and moving conditional
centers.

### 8.4 Scaling

Under \(x\mapsto ax\), both \(C_P\) and the weak variance scale by
\(a^2\), while \(\kappa^{-1}\) scales by \(a^2\). Thus (1.2) has the
correct homogeneity.

### 8.5 Extreme dimensions

* If \(F=\{0\}\), ordinary Lichnerowicz gives the sharper
  \(C_P(\mu)\le\kappa^{-1}\); (1.5) gives the valid but nonsharp
  \(3\kappa^{-1}\).
* If \(E=\{0\}\), (1.5) is tautological and makes no assertion beyond
  the actual marginal gap.
* For \(\dim F=1\), the dimension-free one-dimensional estimate is exactly
  what makes the rank-one result universal.

---

## 9. Relation to stochastic localization and KLS

Suppose a localization argument produces a posterior measure satisfying

\[
 D^2V\succeq\kappa P_{u^\perp}
\]

and

\[
 \operatorname{Var}_\mu\langle X,u\rangle\le R^2.
\]

The theorem immediately gives

\[
 C_P(\mu)\le96(\kappa^{-1}+R^2).                      \tag{9.1}
\]

Thus a localization process may leave one direction unregularized without
losing dimension-free control, provided its surviving variance is bounded.
The conditional mean may rotate or move rapidly; the reinforced Prekopa
estimate already charges that motion to marginal curvature.

What the theorem does **not** prove is that every isotropic log-concave
measure can be transformed into such a posterior while retaining a
universal bound on \(R\), or that the resulting posterior estimate transfers
back with no loss. That adaptive survivor-variance problem is the remaining
KLS-strength obstruction in this route.

Consequently:

* the theorem is a genuine dimension-free result for a broad curvature
  class;
* it is compatible with a localization proof of KLS;
* it is not, by itself, a new general dimension-dependent KLS bound;
* it must not be cited as a proof of KLS.

---

## 10. Constant ledger

The final constants arise as follows.

1. Conditional Lichnerowicz contributes \(\kappa^{-1}\).
2. Lemma 4.1 contributes the factor \(2\).
3. Weight comparison contributes \(C_P(\bar\mu)+\kappa^{-1}\).
4. Hence
   \[
   C_P(\mu)
   \le\kappa^{-1}
      +2(C_P(\bar\mu)+\kappa^{-1}).
   \]
5. In one dimension,
   \[
   C_P(\bar\mu)\le48\operatorname{Var}(T).
   \]
6. Therefore
   \[
   C_P(\mu)
   \le3\kappa^{-1}+96\operatorname{Var}(T)
   \le96(\kappa^{-1}+\operatorname{Var}(T)).
   \]

The scalar weighted-Hardy route instead gives

\[
 C_P(\mu)\le201(\kappa^{-1}+\operatorname{Var}(T)).
\]

No constant depends on the ambient dimension or on a smoothing,
confinement, truncation, or cutoff parameter.

---

## 11. Audit checklist before publication

The following points should be checked independently before treating this
as a finished preprint.

1. **Novelty.** Compare Theorem (1.5) against existing two-scale Poincare,
   asymmetric Brascamp--Lieb, reinforced Prekopa, and improved
   Lichnerowicz inequalities.
2. **Operator domains.** Recheck the closure of \(A_1\) on
   \(\mathcal G\), the intertwining identity, and the variational inverse
   comparison in Lemma 4.1.
3. **Weak Poisson solutions.** Verify measurability in \(z\) and the
   differentiations leading to (3.23) under the chosen approximation.
4. **Nonsmooth curvature.** Verify the posterior Brascamp--Lieb step for
   distributional convexity using a fully specified monotone approximation.
5. **Limits.** Check the exact order of confinement, mollification, cutoff,
   truncation, and Gaussian-regularization limits.
6. **Affine support.** State the theorem intrinsically; do not silently
   restrict an ambient rank-one projector.
7. **Constants.** Distinguish the direct constant \(201\) from the improved
   final constant \(96\).
8. **No circularity.** In the general theorem, only the actual
   \(C_P(\bar\mu)\) is used. Dimension-free control enters only when the
   weak marginal is one-dimensional.
9. **Boundary models.** If convex hard supports are included through
   \(V=+\infty\), verify that intrinsic Gaussian regularization and the
   limiting argument account for boundary terms.
10. **Test functions.** Confirm that the final statement covers every
    locally Lipschitz \(f\in L^2(\mu)\), not merely smooth compactly
    supported functions.

---

## 12. Condensed proof chain

For quick reference, the proof is

\[
 D^2V\succeq\kappa P_E
 \Longrightarrow
 \begin{cases}
 C_P(\nu_z)\le\kappa^{-1},\\[2mm]
 D^2U\succeq\kappa G,
 \end{cases}                                         \tag{12.1}
\]

\[
 \langle(I+G)^{-1}\nabla g,\nabla g\rangle
 \le\mathbb E_z|\nabla f|^2,                          \tag{12.2}
\]

\[
 \operatorname{Var}_{\bar\mu}g
 \le2\int\langle(D^2U+C_P(\bar\mu)^{-1}I)^{-1}
 \nabla g,\nabla g\rangle d\bar\mu,                \tag{12.3}
\]

\[
 (D^2U+C_P(\bar\mu)^{-1}I)^{-1}
 \preceq
 (C_P(\bar\mu)+\kappa^{-1})(I+G)^{-1},               \tag{12.4}
\]

and therefore

\[
 \boxed{
 C_P(\mu)
 \le\kappa^{-1}
 +2(C_P(\bar\mu)+\kappa^{-1}).}                      \tag{12.5}
\]

For a one-dimensional weak marginal,

\[
 C_P(\bar\mu)\le48\operatorname{Var}(T),             \tag{12.6}
\]

which yields the rank-one-defect theorem with universal constant \(96\).
