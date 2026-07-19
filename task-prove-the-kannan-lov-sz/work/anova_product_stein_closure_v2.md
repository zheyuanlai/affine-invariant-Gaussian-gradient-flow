# Affine--orthogonal Gaussian ANOVA for product log-concave inputs (audited wording)

## Scope

Let \(\mu=\bigotimes_{i=1}^n\mu_i\), where each \(\mu_i\) is centered,
variance one, and log-concave. Let \(G\sim N(0,I_n)\) be independent and
\(S=(X+G)/\sqrt2\), with law \(\nu\). For \(f\in L^2_0(\nu)\), put
\[
 F=f(S),\quad h(X)=E[F\mid X],\quad v(G)=E[F\mid G],\quad
 R=F-h(X)-v(G).
\]
Then
\[
 \operatorname {dist}_{L^2(\nu)}(f,\operatorname {Aff})^2
 \le 900\,E R^2. \tag{1}
\]
The constant is independent of \(n\) and of the factors. This is a
tensorized calibration of the unresolved irreducible multidimensional gate.

## 1. Affine reduction and Gaussian jets

Since \(Ef=0\) and \(\nu\) is centered, the \(L^2(\nu)\)-projection of \(f\)
onto the affine space has no constant term; write it \(a\cdot S\). For every
linear \(\ell(S)=a\cdot S\),
\[
 E[\ell(S)\mid X]+E[\ell(S)\mid G]=\ell(S).
\]
Thus subtracting the affine projection leaves \(R\) unchanged. We may assume
\(f\perp\operatorname{Aff}\), so
\[
 Ef=0,\qquad E[Sf(S)]=0. \tag{2}
\]

For a multi-index \(\alpha\), let \(H_\alpha/\sqrt{\alpha!}\) be the
orthonormal Gaussian Hermite basis. Gaussian integration by parts gives
\[
 c_\alpha(x):=E_G\!\left[f((x+G)/\sqrt2)
       {H_\alpha(G)\over\sqrt{\alpha!}}\right]
 ={\partial^\alpha h(x)\over\sqrt{\alpha!}}. \tag{3}
\]
Consequently
\[
 ER^2=\sum_{|\alpha|\ge1}
 {\operatorname {Var}_\mu(\partial^\alpha h)\over\alpha!}. \tag{4}
\]
Define
\[
 V_k=\sum_{|\alpha|=k}{\operatorname {Var}_\mu(\partial^\alpha h)\over\alpha!},
 \qquad
 M_k=\sum_{|\alpha|=k}{(E_\mu\partial^\alpha h)^2\over\alpha!}.
\]
Then Parseval gives
\[
 \|f\|_2^2=V_0+\sum_{k\ge1}(V_k+M_k),\qquad ER^2=\sum_{k\ge1}V_k. \tag{5}
\]
Here \(V_0=\operatorname {Var}_\mu(h)=E h^2\), since \(Eh=Ef=0\).
The affine orthogonality in (2) and
\(E[Sf]=(E[Xh]+E\nabla h)/\sqrt2\) imply
\[
 E[Xh]+E\nabla h=0. \tag{6}
\]

## 2. Zeroth and first jets

Let \(m=E\nabla h\) and \(u=h-m\cdot x\). Since \(EX=0\), \(Eu=0\),
and (6) gives \(E[Xh]=-m\), hence \(E[Xu]=-2m\). Also
\(E|\nabla u|^2=V_1\). Tensorized one-dimensional Poincare (constant 12)
therefore gives
\[
 V_0+3|m|^2=E u^2\le12V_1. \tag{7}
\]
Indeed \(Eh^2=E|u+m\cdot X|^2=Eu^2-3|m|^2\).

## 3. Mean higher jets from diagonal Stein kernels

For each factor let \(\tau_i\) be its canonical Stein kernel. The scalar
median/hazard estimate gives
\[
 E_{\mu_i}\tau_i^2\le400,\qquad E_{\mu_i}\tau_i=1. \tag{8}
\]
Set \(K^2=\sup_iE(\tau_i-1)^2\le400\).

Use ordered derivative tensors
\[
 H_k(x)=(\partial_{i_1}\cdots\partial_{i_k}h(x))_{i_1,\ldots,i_k},
\]
with ordered Hilbert--Schmidt norm. Then
\[
 V_k={E|H_k-EH_k|_{\rm ord}^2\over k!},\qquad
 M_k={|EH_k|_{\rm ord}^2\over k!}. \tag{9}
\]
For \(k\ge2\), integrate by parts in the first index and average over the
\(k\) symmetric slots:
\[
 EH_k=\operatorname {Sym}B_k-\operatorname {Sym}C_k,
\]
where
\[
 (B_k)_{i_1,\ldots,i_k}
 =E[X_{i_1}(H_{k-1}-EH_{k-1})_{i_2,\ldots,i_k}],
\]
\[
 (C_k)_{i_1,\ldots,i_k}
 =E[(\tau_{i_1}-1)(H_k-EH_k)_{i_1,\ldots,i_k}].
\]
The coordinate functions \(X_i\) are an orthonormal family in
\(L^2(\mu)\). Bessel's inequality, applied to each component of
\(H_{k-1}-EH_{k-1}\), gives
\[
 |B_k|_{\rm ord}^2\le E|H_{k-1}-EH_{k-1}|_{\rm ord}^2. \tag{10}
\]
Cauchy--Schwarz and (8) give
\[
 |C_k|_{\rm ord}^2\le K^2 E|H_k-EH_k|_{\rm ord}^2. \tag{11}
\]
Since symmetrization is an orthogonal contraction, (9)--(11) imply
\[
 M_k\le {2\over k}V_{k-1}+2K^2V_k
 \le {2\over k}V_{k-1}+800V_k,\qquad k\ge2. \tag{12}
\]

## 4. Summation

From (12),
\[
 \sum_{k\ge2}M_k
 \le2\sum_{k\ge2}{V_{k-1}\over k}+800\sum_{k\ge2}V_k
 \le2V_1+800\,ER^2. \tag{13}
\]
Using (5), (7), and (13), and absorbing the explicit \(V_1\) into
\(ER^2\),
\[
 \|f\|_2^2
 \le12V_1+ER^2+2V_1+800ER^2
 \le815ER^2. \tag{14}
\]
Because \(f\perp\operatorname{Aff}\), this is exactly the desired squared
distance. Thus the rounded constant 900 is valid.

## 5. Form-domain passage

The calculation is immediate for bounded smooth \(f\). For general
\(f\in L^2(\nu)\), approximate in \(L^2(\nu)\) by bounded smooth functions.
The lift \(f\mapsto f((X+G)/\sqrt2)\), both conditional expectations, and
each Hermite coefficient map are \(L^2\)-contractions. Apply the finite-jet
inequality to the approximants and pass first at fixed jet order, then by
monotone convergence of the nonnegative sums. The one-dimensional Stein
identity is valid in the closed Sobolev form, including interval endpoints;
the same truncation argument used for the median/hazard estimate gives the
factor bound (8). Thus (1) holds for the full \(L^2\) form domain.

