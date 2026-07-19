# Nonzero conditional means: Gaussian weighted-Fisher addendum

The conditional-Gaussian proof in section 5 of dmmse_audit_current.md is
written for \(q_t=N(0,R(t))\). The same output-channel conclusion holds for
smooth Gaussian slices with a varying mean:
\[
q_t=N(m(t),R(t)),\qquad \widetilde R(t):=R(t)+I\succeq I.
\]

Assume the joint \((T,Y)\) law is isotropic, so
\[
E[T\,m(T)]=0,\qquad E\widetilde R(T)\preceq\Lambda I.
\]
For the noisy output \(Z=Y+G\), the conditional law is
\(N(m(t),\widetilde R(t))\). Its scalar-parameter Fisher information is
\[
I_Z(t)
=m'(t)^T\widetilde R(t)^{-1}m'(t)
+\frac12\operatorname{tr}\!\left(
\widetilde R^{-1}\widetilde R'\widetilde R^{-1}\widetilde R'\right).
\]

For the centered conditional Poisson field, the covariance/shape component
has the exact noncommuting identity
\[
C_t+B_t^0
=\frac12\operatorname{tr}\!\left(
\widetilde R^{-1}\widetilde R'\widetilde R^{-1}\widetilde R'\right).
\]
Indeed, in an eigenbasis of \(\widetilde R\), the Lyapunov solution is
\(A_{ij}=-\widetilde R'_{ij}/(\tilde r_i+\tilde r_j)\), and pairing the
ordered \(C\) and \(B^0\) terms gives the covariance Fisher coefficient
exactly.

Joint convexity gives
\[
\Phi''=W''+\frac12\operatorname{tr}\bigl(
\widetilde R^{-1}(-\widetilde R'')\bigr)
+\frac12\operatorname{tr}\!\left(
\widetilde R^{-1}\widetilde R'\widetilde R^{-1}\widetilde R'\right),
\]
so the shape component is bounded by the marginal curvature. The scalar
centroid lemma for matrix-concave \(\widetilde R\), using
\(E[Tm(T)]=0\) and \(E\widetilde R\preceq\Lambda I\), gives
\[
E[\tau(T)^2|m'(T)|^2]\le C\Lambda.
\]
Since \(\widetilde R^{-1}\preceq I\),
\[
E[\tau(T)^2 I_Z(T)]\le 1+C\Lambda.
\]
Thus the weighted Fisher gate used in (0.4) is dimension-free for this
nonzero-mean Gaussian subclass, including noncommuting covariance rotation.
The assertion is \(C^2\), full-support, post-noise only; no zero-noise limit
or general-slice extension is implied.
