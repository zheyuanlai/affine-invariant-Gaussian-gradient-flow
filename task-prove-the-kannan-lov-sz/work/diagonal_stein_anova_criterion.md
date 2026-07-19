# A diagonal-Stein criterion for the Gaussian ANOVA gate

Let \(\mu\) be centered and isotropic on \(\mathbb R^n\). Assume:

1. (Diagonal Stein identities) There are measurable \(\tau_i\) such that for every
   smooth compactly supported \(u\),
   \[
   E_\mu[X_i u(X)]=E_\mu[\tau_i(X)\partial_i u(X)],\qquad i=1,\ldots,n.
   \tag{1}
   \]
   Testing \(u=X_i\) gives \(E\tau_i=1\).
2. (Uniform row budget) \(K^2:=\sup_i E(\tau_i-1)^2<\infty\).
3. (Poincare) \(\operatorname {Var}_\mu u\le C_0E|\nabla u|^2\).

Let \(G\sim N(0,I_n)\), \(S=(X+G)/\sqrt2\), and \(\nu=\mathcal L(S)\). Then
\[
 \operatorname {dist}_{L^2(\nu)}(f,\operatorname {Aff})^2
 \le (2C_0+2+2K^2+1)\, E R^2, \tag{2}
\]
for every centered \(f\in L^2(\nu)\), where
\(R=f(S)-E[f(S)\mid X]-E[f(S)\mid G]\). In particular, the constant is
dimension-free whenever \(C_0,K\) are.

## Proof

Subtract the linear \(L^2(\nu)\)-projection of \(f\) onto the affine space; this
does not change \(R\). Thus \(f\perp\operatorname{Aff}\). Let
\(h(x)=E_G f((x+G)/\sqrt2)\). The Gaussian Hermite expansion gives, with
\[
 V_k=\sum_{|\alpha|=k}{\operatorname {Var}(\partial^\alpha h)\over\alpha!},
 \qquad
 M_k=\sum_{|\alpha|=k}{(E\partial^\alpha h)^2\over\alpha!},
\]
the exact identities
\[
 ER^2=\sum_{k\ge1}V_k,\qquad
 \|f\|_2^2=V_0+\sum_{k\ge1}(V_k+M_k). \tag{3}
\]
Affine orthogonality says \(E[Xh]+E\nabla h=0\). Put
\(m=E\nabla h\), \(u=h-m\cdot x\). Then
\[
 V_0+3|m|^2=Eu^2\le C_0V_1. \tag{4}
\]

For ordered tensors \(H_k=D^kh\), define
\[
 V_k={E|H_k-EH_k|_{\rm ord}^2\over k!},\qquad
 M_k={|EH_k|_{\rm ord}^2\over k!}. \tag{5}
\]
For \(k\ge2\), diagonal Stein integration by parts, averaged over the symmetric
slots, gives \(EH_k=\operatorname{Sym}(B_k-C_k)\), with
\[
 (B_k)_{i_1\ldots i_k}
 =E[X_{i_1}(H_{k-1}-EH_{k-1})_{i_2\ldots i_k}],
\]
\[
 (C_k)_{i_1\ldots i_k}
 =E[(\tau_{i_1}-1)(H_k-EH_k)_{i_1\ldots i_k}].
\]
Bessel's inequality for the orthonormal family \(X_i\) and Cauchy--Schwarz
give
\[
 |B_k|_{\rm ord}^2\le E|H_{k-1}-EH_{k-1}|_{\rm ord}^2,\qquad
 |C_k|_{\rm ord}^2\le K^2E|H_k-EH_k|_{\rm ord}^2.
\]
Therefore
\[
 M_k\le {2\over k}V_{k-1}+2K^2V_k. \tag{6}
\]
Summing (6) for \(k\ge2\) and using \(1/k\le1\),
\[
 \sum_{k\ge2}M_k\le 2(1+K^2)ER^2. \tag{7}
\]
From (3), (4), and \(V_1\le ER^2\),
\[
 \|f\|_2^2
 \le C_0V_1+ER^2+2(1+K^2)ER^2
 \le (C_0+3+2K^2)ER^2. \tag{8}
\]
The displayed bound (2) is a looser rounded form. Form-domain passage follows
by bounded smooth approximation and closedness of the conditional Hermite
and Stein forms.

