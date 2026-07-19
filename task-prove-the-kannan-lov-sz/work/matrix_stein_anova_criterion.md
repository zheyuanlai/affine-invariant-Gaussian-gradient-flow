# Frame-free matrix-Stein criterion and its sharp obstruction

## 1. Criterion

Let \(\mu\) be centered and isotropic on \(\mathbb R^n\). Suppose it admits a
(square-integrable) matrix Stein kernel \(\tau(x)\) in the sense that, for
every smooth compactly supported scalar \(u\),
\[
 E[Xu(X)]=E[\tau(X)\nabla u(X)],\qquad E\tau=I. \tag{1}
\]
No symmetry is needed below; replace \((\tau-I)^2\) by
\((\tau-I)(\tau-I)^T\). Assume
\[
 K^2:=\left\|E[(\tau-I)(\tau-I)^T]\right\|_{\rm op}<\infty. \tag{2}
\]
If, in addition, \(\mu\) has Poincare constant \(C_0\), then for
\(S=(X+G)/\sqrt2\) and
\(R=f(S)-E[f(S)\mid X]-E[f(S)\mid G]\),
\[
 \operatorname {dist}_{L^2(\mathcal L(S))}(f,\operatorname {Aff})^2
 \le (C_0+3+2K^2)\,E R^2. \tag{3}
\]
Thus (2) is a frame-free replacement for the diagonal row budget in the
product proof.

## 2. Ordered-tensor proof

Subtract the linear affine projection of a centered \(f\); \(R\) is unchanged,
and assume \(f\perp\operatorname{Aff}\). Let
\(h(x)=E_G f((x+G)/\sqrt2)\). The Gaussian Hermite expansion gives
\[
 ER^2=\sum_{k\ge1}V_k,\qquad
 \|f\|_2^2=V_0+\sum_{k\ge1}(V_k+M_k), \tag{4}
\]
where
\[
 V_k={E|D^kh-E D^kh|_{\rm ord}^2\over k!},\qquad
 M_k={|E D^kh|_{\rm ord}^2\over k!}. \tag{5}
\]
Affine orthogonality gives \(E[Xh]+E\nabla h=0\). If
\(m=E\nabla h\) and \(u=h-m\cdot x\), then
\[
 V_0+3|m|^2=Eu^2\le C_0V_1. \tag{6}
\]

For \(k\ge2\), write \(H_k=D^kh\) as an ordered symmetric tensor. Integrate
the Stein identity in one distinguished slot and average over the \(k\) slots:
\[
 E H_k=\operatorname{Sym}(B_k-C_k),
\]
with
\[
 (B_k)_{i_1\ldots i_k}
 =E[X_{i_1}(H_{k-1}-EH_{k-1})_{i_2\ldots i_k}],
\]
\[
 (C_k)_{i_1\ldots i_k}
 =E[((\tau-I)_{i_1j})(H_k-EH_k)_{j\,i_2\ldots i_k}].
\]
Bessel's inequality for the orthonormal family \(X_1,\ldots,X_n\) gives
\[
 |B_k|_{\rm ord}^2\le E|H_{k-1}-EH_{k-1}|_{\rm ord}^2. \tag{7}
\]
For each fixed \((i_2,\ldots,i_k)\), put
\(Y=(H_k-EH_k)_{\bullet\,i_2\ldots i_k}\). Then
\[
 \left|E[(\tau-I)Y]\right|^2
 \le \left\|E[(\tau-I)(\tau-I)^T]\right\|_{\rm op}\,E|Y|^2
 \le K^2E|Y|^2, \tag{8}
\]
by duality and Cauchy--Schwarz. Summing and using contractivity of slot
symmetrization gives
\[
 M_k\le {2\over k}V_{k-1}+2K^2V_k,\qquad k\ge2. \tag{9}
\]
Consequently
\[
 \sum_{k\ge2}M_k\le2(1+K^2)ER^2. \tag{10}
\]
Combining (4), (6), and (10), with \(V_1\le ER^2\), proves (3).

## 3. Why the criterion is not currently available for general log-concavity

Known isotropic log-concave Stein kernels have a Hilbert--Schmidt trace
budget of order \(n\),
\[
 E\|\tau-I\|_{\rm HS}^2\lesssim n,
\]
but this only gives
\(\operatorname {tr}E[(\tau-I)(\tau-I)^T]\lesssim n\). It does not imply (2).
The ordered-tensor error at \(k=2\) already contains
\[
 E[(\tau-I)Y],\qquad Y=D^2h-E D^2h,
\]
whose norm is controlled by the operator, not trace, norm in (2).

The algebraic gap is sharp: a deterministic rank-one matrix
\(A=\sqrt n\,e_1e_1^T\) has \(\|A\|_{\rm HS}^2=n\) but
\(\|A^2\|_{\rm op}=n\). Thus a trace budget permits an error of order
\(\sqrt n\) in one direction. No rearrangement of the Hilbert--Schmidt
bound removes this possibility.

Moreover, a universal bound (2) for a suitable minimal Stein kernel would be
a strong directional statement. Testing the Stein identity on quadratic
functions gives the third-moment map
\[
 u\longmapsto E[(u\cdot X)(XX^T-I)].
\]
Dimension-free Hilbert--Schmidt control of this map is implied by (2), while
the converse spectral stability estimate is known to be KLS-strength. Hence
(2) cannot be inserted as a routine regularization fact.

The criterion is therefore a genuine conditional extension (covering, for
example, products with diagonal kernels), not a proof for arbitrary
log-concave measures.

