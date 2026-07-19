# Conditional-coordinate Stein obstruction beyond product laws

This note records the exact point at which the product ANOVA proof fails for a
general isotropic log-concave law. It is an obstruction identity, not an
assumption.

Let \(X\) be centered and isotropic in \(\mathbb R^n\). Disintegrate \(\mu\)
with respect to \(X_{-i}\), and write
\[
 m_i(x_{-i})=E[X_i\mid X_{-i}=x_{-i}],
 \qquad Z_i=X_i-m_i(X_{-i}).
\]
The conditional law of \(Z_i\) is centered one-dimensional log-concave.
Let \(\tau_i^{\rm c}(x_i;x_{-i})\) be its canonical conditional Stein kernel.
For smooth \(u\), conditional integration by parts gives the exact identity
\[
 E[X_i u(X)]
 =E[m_i(X_{-i})u(X)]
  +E[\tau_i^{\rm c}(X)\,\partial_i u(X)]. \tag{1}
\]

For a smooth \(h\), put \(H_k=D^kh\) in ordered tensor notation. The
coordinatewise step used in the product proof becomes, for \(k\ge1\),
\[
 E H_k=\operatorname{Sym}(B_k+L_k-C_k), \tag{2}
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
and
\[
 (C_k)_{i_1\ldots i_k}
 =E[(\tau_{i_1}^{\rm c}-E[\tau_{i_1}^{\rm c}\mid X_{-i_1}])
      (H_k-EH_k)_{i_1\ldots i_k}].
\]
(The centering in \(C_k\) is harmless because the conditional centered Stein
kernel has zero conditional mean after subtracting its conditional variance;
one may instead retain the uncentered kernel and absorb its mean into \(L_k\).)

The \(B_k\)-term has the dimension-free Bessel estimate
\[
 |B_k|_{\rm ord}^2
 \le E|H_{k-1}-EH_{k-1}|_{\rm ord}^2, \tag{3}
\]
because isotropy makes \(X_1,\ldots,X_n\) an orthonormal family. The conditional
Stein term is bounded by the conditional fourth-moment budget
\[
 |C_k|_{\rm ord}^2
 \le \Big(\sup_i E[(\tau_i^{\rm c}-E[\tau_i^{\rm c}\mid X_{-i}])^2]\Big)
      E|H_k-EH_k|_{\rm ord}^2, \tag{4}
\]
whenever the displayed supremum is finite.

The unresolved term is \(L_k\). A direct Hilbert bound is
\[
 |L_k|_{\rm ord}^2
 \le \sum_i E[m_i(X_{-i})^2]\,
      E\!\left|H_{k-1}-EH_{k-1}\right|_{\rm ord}^2, \tag{5}
\]
which is \(O(n)\) in general and therefore useless. A sharper operator
version would require a dimension-free bound on the regression map
\[
 T:g\longmapsto
 \big(E[m_i(X_{-i})g(X)\mid X_{-i}]\big)_{i=1}^n
\]
from centered scalar/tensor \(L^2(\mu)\) to \(L^2(\mu;\mathbb R^n)\). Such a
bound is a nonlinear conditional-correlation statement; no universal bound
is presently known, and inserting one would be a KLS-strength step.

Two classes eliminate \(L_k\) exactly:

* Product laws: \(m_i\equiv0\).
* Unconditional laws: each conditional fiber is symmetric, so \(m_i\equiv0\).

For a general law, even isotropy only says \(E[X_iX_j]=0\) for \(i\ne j\); it
does not force \(m_i=0\) or control the nonlinear regression operator in (5).
Thus the product proof cannot be promoted by simply replacing scalar Stein
kernels with conditional ones.

