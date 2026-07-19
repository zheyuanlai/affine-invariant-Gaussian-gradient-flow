# Conditional-coordinate Stein obstruction beyond product laws (corrected)

This note records the exact point at which the product ANOVA proof fails for a
general isotropic log-concave law. It is an obstruction identity, not an
assumption.

Let \(X\) be centered and isotropic in \(\mathbb R^n\). Disintegrate with
respect to \(X_{-i}\), and write
\[
 m_i(x_{-i})=E[X_i\mid X_{-i}=x_{-i}],\qquad
 \sigma_i^2(x_{-i})=\operatorname {Var}(X_i\mid X_{-i}=x_{-i}).
\]
The conditional centered law \(Z_i=X_i-m_i(X_{-i})\) is one-dimensional
log-concave. Let \(\tau_i^{\rm c}(x_i;x_{-i})\) be its canonical Stein kernel,
so \(E[\tau_i^{\rm c}\mid X_{-i}]=\sigma_i^2\). For smooth \(u\), conditional
integration by parts gives
\[
 E[X_i u(X)]
 =E[m_i(X_{-i})u(X)]
  +E[\tau_i^{\rm c}(X)\,\partial_i u(X)]. \tag{1}
\]
Write \(\widetilde\tau_i=\tau_i^{\rm c}-\sigma_i^2\), which has conditional
mean zero.

For a smooth \(h\), put \(H_k=D^kh\) in ordered tensor notation. Applying (1)
to each slot and averaging over the symmetric slots gives, for \(k\ge1\),
\[
 EH_k=\operatorname{Sym}(B_k-L_k-D_k-C_k), \tag{2}
\]
where
\[
 (B_k)_{i_1\ldots i_k}
 =E[X_{i_1}(H_{k-1}-EH_{k-1})_{i_2\ldots i_k}],
\]
\[
 (L_k)_{i_1\ldots i_k}
 =E[m_{i_1}(X_{-i_1})
       (H_{k-1}-EH_{k-1})_{i_2\ldots i_k}],
\]
\[
 (D_k)_{i_1\ldots i_k}
 =E[(\sigma_{i_1}^2-1)
       (H_k-EH_k)_{i_1\ldots i_k}],
\]
\[
 (C_k)_{i_1\ldots i_k}
 =E[\widetilde\tau_{i_1}
       (H_k-EH_k)_{i_1\ldots i_k}].
\]
Indeed \(E[X_iH_{k-1}]=E[m_iH_{k-1}]
+E[\sigma_i^2H_k]+E[\widetilde\tau_iH_k]\), and \(EX_i=0\).

The \(B_k\)-term has the dimension-free Bessel estimate
\[
 |B_k|_{\rm ord}^2
 \le E|H_{k-1}-EH_{k-1}|_{\rm ord}^2, \tag{3}
\]
because isotropy makes \(X_1,\ldots,X_n\) an orthonormal family. Conditional
one-dimensional moment comparison controls \(C_k\) only through a quantity
involving \(E\sigma_i^4\); no dimension-free sum over \(i\) follows from that
scalar fact. The variance-deviation term \(D_k\) has the same issue.

The regression term admits only the crude bound
\[
 |L_k|_{\rm ord}^2
 \le \Big(\sum_i E[m_i(X_{-i})^2]\Big)
      E|H_{k-1}-EH_{k-1}|_{\rm ord}^2, \tag{4}
\]
which can be \(O(n)\). A sharper bound would require a dimension-free
operator estimate for the nonlinear conditional-regression map
\[
 g\longmapsto
 \big(E[m_i(X_{-i})g(X)\mid X_{-i}]\big)_{i=1}^n.
\]
Likewise, controlling \(D_k\) dimension-freely would require a bound on the
operator-valued conditional variance deviations
\(\sigma_i^2-1\). Neither estimate follows from isotropy.

For product laws, \(m_i\equiv0\) and \(\sigma_i^2\equiv1\), so \(L_k=D_k=0\)
and the proof in anova_product_stein_closure_v3.md applies. Unconditional
symmetry forces \(m_i=0\) but generally does **not** force
\(\sigma_i^2=1\); the \(D_k\) term remains and the product proof does not
automatically extend.

Thus the failure of the product argument is explicit: it is not legitimate
to replace the diagonal factor Stein kernels by conditional ones without
paying both the nonlinear regression \(L_k\) and conditional-variance \(D_k\)
terms.

