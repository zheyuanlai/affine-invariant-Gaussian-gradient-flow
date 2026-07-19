# Audited product ANOVA closure (final ledger)

Let \(\mu=\bigotimes_{i=1}^n\mu_i\), with each factor centered,
variance-one, and log-concave. Let \(G\sim N(0,I_n)\) be independent,
\(S=(X+G)/\sqrt2\), \(\nu=\mathcal L(S)\). For centered
\(f\in L^2(\nu)\), define
\[
 h(X)=E[f(S)\mid X],\quad v(G)=E[f(S)\mid G],\quad
 R=f(S)-h(X)-v(G).
\]
Then
\[
 \operatorname {dist}_{L^2(\nu)}(f,\operatorname {Aff})^2
 \le 900\,E R^2. \tag{1}
\]

## 1. Reduction and jet identities

Subtract the (linear, since \(Ef=0\)) \(L^2(\nu)\)-projection onto the affine
space. Its residual \(R\) is unchanged because
\(E[a\cdot S\mid X]+E[a\cdot S\mid G]=a\cdot S\). Thus assume
\(f\perp\operatorname{Aff}\). Let
\[
 c_\alpha(x)=E_G\!\left[f((x+G)/\sqrt2)
              H_\alpha(G)/\sqrt{\alpha!}\right]
            =\partial^\alpha h(x)/\sqrt{\alpha!}.
\]
Conditional Parseval gives
\[
 ER^2=\sum_{|\alpha|\ge1}
       {\operatorname {Var}_\mu(\partial^\alpha h)\over\alpha!}. \tag{2}
\]
Set
\[
 V_k=\sum_{|\alpha|=k}{\operatorname {Var}(\partial^\alpha h)\over\alpha!},
 \qquad
 M_k=\sum_{|\alpha|=k}{(E\partial^\alpha h)^2\over\alpha!}.
\]
Then
\[
 \|f\|_2^2=V_0+\sum_{k\ge1}(V_k+M_k),\qquad ER^2=\sum_{k\ge1}V_k, \tag{3}
\]
where \(V_0=Eh^2\), and affine orthogonality gives
\(E[Xh]+E\nabla h=0\).

Put \(m=E\nabla h\), \(u=h-m\cdot x\). Then \(Eu=0\),
\(E[Xu]=-2m\), \(E|\nabla u|^2=V_1\), and hence, by tensorized
one-dimensional Poincare with constant 12,
\[
 V_0+3|m|^2=Eu^2\le12V_1. \tag{4}
\]

## 2. Higher mean jets

Let \(\tau_i\) be the canonical Stein kernel of \(\mu_i\). The scalar
hazard bound gives \(E\tau_i=1\), \(E\tau_i^2\le400\); set
\(K^2=\sup_iE(\tau_i-1)^2\le400\).
For ordered derivative tensors \(H_k=D^kh\), use
\[
 V_k={E|H_k-EH_k|_{\rm ord}^2\over k!},\qquad
 M_k={|EH_k|_{\rm ord}^2\over k!}. \tag{5}
\]
For \(k\ge2\), coordinatewise Stein integration and slot symmetrization give
\[
 EH_k=\operatorname{Sym}B_k-\operatorname{Sym}C_k,
\]
with
\[
 (B_k)_{i_1\ldots i_k}
 =E[X_{i_1}(H_{k-1}-EH_{k-1})_{i_2\ldots i_k}],
\quad
 (C_k)_{i_1\ldots i_k}
 =E[(\tau_{i_1}-1)(H_k-EH_k)_{i_1\ldots i_k}].
\]
Bessel's inequality for the orthonormal family \(X_1,\ldots,X_n\) gives
\[
 |B_k|_{\rm ord}^2\le E|H_{k-1}-EH_{k-1}|_{\rm ord}^2,
\]
and Cauchy--Schwarz gives
\[
 |C_k|_{\rm ord}^2\le K^2 E|H_k-EH_k|_{\rm ord}^2.
\]
Since symmetrization is contractive,
\[
 M_k\le {2\over k}V_{k-1}+2K^2V_k
      \le {2\over k}V_{k-1}+800V_k.                    \tag{6}
\]

## 3. Summation

Because \(1/k\le1\),
\[
 \sum_{k\ge2}M_k
 \le2\sum_{k\ge2}V_{k-1}+800\sum_{k\ge2}V_k
 \le802\,ER^2.                                        \tag{7}
\]
Using (3), (4), and (7),
\[
 \|f\|_2^2
 \le12V_1+ER^2+802ER^2
 \le815\,ER^2.                                        \tag{8}
\]
Here \(M_1=|m|^2\) is included in (4), and the explicit \(V_1\) is bounded by
\(ER^2\). Since \(f\perp\operatorname{Aff}\), (8) is the desired distance
bound, with the rounded constant 900.

## 4. Form-domain details

For bounded smooth \(f\), all identities follow from Gaussian and
one-dimensional Stein integration by parts. For general \(L^2(\nu)\) functions,
approximate by bounded smooth functions. The lift, conditional expectations,
and Hermite coefficient maps are \(L^2\)-contractions. Apply the finite-jet
bounds, pass at each finite order, and then use monotone convergence. The
one-dimensional Stein identity and \(E\tau_i^2\le400\) hold in the closed form
on finite intervals by the canonical endpoint limits. Thus the constant is
uniform over all product factors and all dimensions.

