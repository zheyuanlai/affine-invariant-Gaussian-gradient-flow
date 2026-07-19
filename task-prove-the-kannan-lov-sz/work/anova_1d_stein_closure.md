# One-dimensional closure of the affine--orthogonal Gaussian ANOVA gate

This note proves a genuine dimension-free subcase of the affine--orthogonal
interaction estimate. It is a calibrated lemma; it does **not** extend to
arbitrary irreducible multidimensional log-concave laws.

## 1. Statement

Let \(X\) have a centered, variance-one, one-dimensional log-concave law
\(\mu\), let \(G\sim N(0,1)\) be independent, and set
\[
 S=(X+G)/\sqrt2,\qquad \nu=\mathcal L(S).
\]
For \(f\in L^2_0(\nu)\), define
\[
 F(x,g)=f((x+g)/\sqrt2),\quad
 h(x)=E_G F(x,G),\quad v(g)=E_XF(X,g),
\]
\[
 R=F-h(X)-v(G).
\]
There is a numerical constant \(C_{\rm 1d}\) (one may take
\(C_{\rm 1d}=900\)) such that
\[
 \operatorname {dist}_{L^2(\nu)}(f,\operatorname {Aff})^2
 \le C_{\rm 1d} E R^2. \tag{1.1}
\]
The assertion is intrinsic and remains valid when the original ambient law
is supported on a one-dimensional affine line (the Gaussian is also taken in
that line).

## 2. Conditional Hermite ledger

Write \(H_k\) for the probabilists' Hermite polynomials and
\(\phi_k=H_k/\sqrt{k!}\). For bounded smooth \(f\), Gaussian integration by
parts gives
\[
 c_k(x):=E_G[F(x,G)\phi_k(G)]={h^{(k)}(x)\over\sqrt{k!}}. \tag{2.1}
\]
Conditional Parseval, followed by subtracting the \(X\)- and \(G\)-conditional
projections, yields
\[
 E R^2=\sum_{k\ge1}{\operatorname {Var}_\mu(h^{(k)})\over k!}. \tag{2.2}
\]
Put
\[
 V_k=\operatorname {Var}_\mu(h^{(k)}),\qquad
 M_k=(E_\mu h^{(k)})^2.
\]
The same Parseval identity gives
\[
 \|f\|_{L^2(\nu)}^2=V_0+\sum_{k\ge1}{V_k+M_k\over k!}, \tag{2.3}
\]
because \(Eh=Ef=0\). Since \(\nu\) is centered and variance one,
\[
 E_\nu[Sf(S)]={E_\mu[Xh]+E_\mu h'\over\sqrt2}, \tag{2.4}
\]
and therefore
\[
 \operatorname {dist}_{L^2(\nu)}(f,\operatorname {Aff})^2
 =V_0+\sum_{k\ge1}{V_k+M_k\over k!}
  -{(E[Xh]+Eh')^2\over2}. \tag{2.5}
\]
For the upper bound we may discard the final non-positive term.

## 3. A reverse derivative estimate from the one-dimensional Stein kernel

Let \(\rho\) be a centered variance-one log-concave density, including finite
interval support by taking one-sided limits. Its canonical Stein kernel is
\[
 \tau(x)={1\over\rho(x)}\int_x^\infty t\rho(t)\,dt,
\]
with the usual zero extension at endpoints. It satisfies, for every locally
absolutely continuous \(u\) in the Stein form domain,
\[
 E[Xu(X)]=E[\tau(X)u'(X)],\qquad E\tau=1. \tag{3.1}
\]
The median/hazard argument in
upper_curvature_wfi_median_addendum.md gives the uniform bound
\[
 E\tau^2\le400. \tag{3.2}
\]
Set \(K^2=E(\tau-1)^2\le E\tau^2\le400\). For any \(u\) with
\(u,u'\in L^2(\mu)\), write \(\bar u=u-Eu\) and \(m=Eu'\). From (3.1),
\[
 E[X\bar u]=m+E[(\tau-1)(u'-m)]. \tag{3.3}
\]
Since \(EX^2=1\) and \(E(\tau-1)=0\), Cauchy--Schwarz gives the key estimate
\[
 |Eu'|^2\le 2\operatorname {Var}_\mu(u)+2K^2\operatorname {Var}_\mu(u'). \tag{3.4}
\]
This is a reverse Poincare estimate for the *mean derivative*; it is not a
spectral-gap assertion.
We also use the elementary one-dimensional log-concave Poincare bound
\[
 \operatorname {Var}_\mu(u)\le 12 E(u')^2. \tag{3.5}
\]
The numerical value 12 is valid for every centered variance-one log-concave
law, with the closed-form interpretation on an interval.
Apply (3.4) to \(u=h^{(k-1)}\) for \(k\ge1\):
\[
 M_k\le 2V_{k-1}+800V_k. \tag{3.6}
\]

## 4. Summation and constant

Using (3.5) with \(u=h\), \(V_0\le12V_1\). Hence, with
\[
 \mathcal R:=\sum_{k\ge1}{V_k\over k!}=ER^2,
\]
\[
\begin{aligned}
 \sum_{k\ge1}{M_k\over k!}
 &\le 2\sum_{k\ge1}{V_{k-1}\over k!}+800\mathcal R \\
 &\le 2V_0+2\mathcal R+800\mathcal R
 \le 24V_1+802\mathcal R
 \le826\mathcal R. \tag{4.1}
\end{aligned}
\]
Combining (2.3), \(V_0\le12V_1\), (2.5), and (4.1),
\[
 \operatorname {dist}(f,\operatorname {Aff})^2
 \le 12V_1+\mathcal R+826\mathcal R
 \le 839\mathcal R.
\]
Thus \(C_{\rm 1d}=900\) is a safe rounded constant.

## 5. Form-domain and hard-support passage

For bounded smooth \(f\), all preceding identities follow by ordinary
Gaussian integration by parts. For general \(f\in L^2(\nu)\), choose bounded
smooth \(f_j\to f\) in \(L^2(\nu)\). The lifted functions
\[
 F_j(X,G)=f_j((X+G)/\sqrt2)
\]
converge in \(L^2(\mu\otimes\gamma)\). Conditional expectation is an
\(L^2\)-contraction, so \(R_j\to R\) in \(L^2\). The coefficient maps
\[
 f\longmapsto E_G[F(X,G)\phi_k(G)]
\]
are also contractions from \(L^2(\nu)\) to \(L^2(\mu)\). Consequently (2.2)
passes to the limit by Parseval and monotone convergence of the partial sums.
The Stein estimate (3.4) is first applied to finite Hermite jets of \(h_j\),
then passed by weak lower semicontinuity. Finally (2.5) is an \(L^2\) identity,
so the stated constant survives the limit.
If \(\mu\) is uniform on a bounded interval, use the canonical Stein kernel on
the interior and approximate \(u\) by restrictions of smooth functions on a
slightly larger interval. The boundary terms vanish because the Stein
kernel is zero at each endpoint. Thus no unverified smooth-density or
strict-convexity assumption enters (1.1).

## 6. Why this does not close the multidimensional target

The step (3.4) uses the scalar Stein kernel and the scalar identity
\(E[Xu]=E[\tau u']\). In dimensions \(n>1\), the analogous estimate with a
uniform operator bound on a Stein kernel is not known; its trace bound loses
\(\sqrt n\), and the resulting assertion would contain the unresolved
generalized quadratic-variance/KLS gate. Thus (1.1) is a valid
one-dimensional calibration, not a proof of the conjecture.

