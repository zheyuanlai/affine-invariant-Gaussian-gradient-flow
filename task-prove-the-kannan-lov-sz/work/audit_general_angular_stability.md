# Independent audit of general angular stability

## Verdict

The argument in general_angular_stability.md is correct.  In particular,
the threshold pullback and the Hilbert--Schmidt quadratic perturbation are
genuinely dimension free.  No hidden trace, ambient fourth moment, or
\(\sqrt n\) loss occurs.

Two presentation improvements are advisable but do not change the result.

1. The threshold-smoothing proof can be made shorter and completely explicit
   by comparing both events to the same threshold \(a\); no preliminary
   Euclidean estimate on \(|a-b|\) is needed.
2. For a lower-dimensional affine hull, every Gaussian, matrix, and
   Hilbert--Schmidt norm should be read on the parallel linear space \(E\),
   not on the ambient \(\mathbb R^n\).

With
\[
 \zeta_\delta(\varepsilon)
 =C_\delta\varepsilon\sqrt{\log(e/\varepsilon)}
\]
and
\[
 \rho=C_\delta\{\sqrt\varepsilon+\zeta_\delta(\varepsilon)^{1/3}\},
\]
the composition proves
\[
 \|PD\|_{\mathrm{HS}}^2
 \le C_\delta t^{-2}\{\varepsilon+\sqrt\rho\}.
\]
Thus the displayed modulus
\(\Omega_\delta(\varepsilon)=
\varepsilon+\sqrt{\rho_\delta(\varepsilon)}\)
is valid and behaves as
\(O_\delta(\varepsilon^{1/6}\log(e/\varepsilon)^{1/12})\).

---

## 1. Threshold-source error

After standardization, put
\[
 Y=\langle u,T(G)\rangle,\qquad Z=\langle u,G\rangle.
\]
The clean-room Brenier estimate gives
\[
 \mathbb E|Y-Z|^2\le\zeta.
                                                                    \tag{1.1}
\]
Let
\[
 E=\{Y\ge a\},\qquad F=\{Z\ge b\},
\]
where both events have probability \(g\).  For an arbitrary \(h>0\), compare
\(E\) first with \(F_a=\{Z\ge a\}\).  On
\(\{|Y-Z|\le h\}\), the events \(E\) and \(F_a\) can differ only when
\(|Z-a|\le h\).  Hence
\[
 \mathbb P(E\mathbin{\triangle}F_a)
 \le \frac{\zeta}{h^2}
     \mathbb P(|Z-a|\le h)
 \le \frac{\zeta}{h^2}+\frac{2h}{\sqrt{2\pi}}.                  \tag{1.2}
\]

The two Gaussian threshold events \(F_a\) and \(F\) are nested.  Since
\(\mathbb P(E)=\mathbb P(F)=g\),
\[
\begin{aligned}
 \mathbb P(F_a\mathbin{\triangle}F)
 &=|\mathbb P(F_a)-\mathbb P(F)|\\
 &=|\mathbb P(F_a)-\mathbb P(E)|\\
 &\le\mathbb P(E\mathbin{\triangle}F_a).                        \tag{1.3}
\end{aligned}
\]
Combining (1.2)--(1.3),
\[
 \mathbb P(E\mathbin{\triangle}F)
 \le2\left(\frac{\zeta}{h^2}+\frac{2h}{\sqrt{2\pi}}\right).
                                                                    \tag{1.4}
\]
Taking \(h=\zeta^{1/3}\) proves
\[
 \mathbb P(E\mathbin{\triangle}F)\le C\zeta^{1/3}.              \tag{1.5}
\]
This estimate is universal; centrality is only needed elsewhere to keep
the active centroid bounded away from zero.

Now let
\[
 A=T^{-1}(S),\qquad B=\{\langle G,u\rangle\ge b\},
\]
and let \(H=\{\langle X,u\rangle\ge a\}\) be the target halfspace.  Then
\[
\begin{aligned}
 \gamma_E(A\mathbin{\triangle}B)
 &\le\gamma_E(A\mathbin{\triangle}T^{-1}H)
   +\gamma_E(T^{-1}H\mathbin{\triangle}B)\\
 &=\pi(S\mathbin{\triangle}H)
   +\mathbb P(E\mathbin{\triangle}F)\\
 &\le C_\delta\sqrt\varepsilon+
       C\zeta_\delta(\varepsilon)^{1/3}.                        \tag{1.6}
\end{aligned}
\]
There is no dimension-dependent anti-concentration term: only the density
bound of the one-dimensional standard Gaussian appears.

---

## 2. First-moment perturbation

Write
\[
 h=\mathbf1_A-\mathbf1_B,\qquad
 \rho=\gamma_E(A\mathbin{\triangle}B).
\]
The equal masses give \(\mathbb Eh=0\), while
\(\mathbb Eh^2=\rho\).  If \(X=T(G)\) is centered and
\(\operatorname{Cov}X\preceq I_E\), then, for every unit \(\theta\),
\[
\begin{aligned}
 |\langle\theta,v_A-v_B\rangle|
 &=|\mathbb E[h\langle\theta,X\rangle]|\\
 &\le\sqrt{\rho}\,
      \sqrt{\operatorname{Var}\langle\theta,X\rangle}
 \le\sqrt\rho.                                                  \tag{2.1}
\end{aligned}
\]
Taking the supremum in \(\theta\) proves
\[
 |v_A-v_B|\le\sqrt\rho.                                        \tag{2.2}
\]

---

## 3. Hilbert--Schmidt quadratic perturbation

Let \(M=M^{\mathsf T}\) with \(\|M\|_{\mathrm{HS}}=1\).  The target is
one-strongly log-concave, so its closed Poincare form gives
\[
\begin{aligned}
 \operatorname{Var}(X^{\mathsf T}MX)
 &\le\mathbb E|\nabla(X^{\mathsf T}MX)|^2\\
 &=4\mathbb E|MX|^2\\
 &=4\operatorname{tr}(M^2\operatorname{Cov}X)
 \le4\operatorname{tr}M^2=4.                                  \tag{3.1}
\end{aligned}
\]
This is the crucial dimension-free step: the covariance bound is paired
with \(\operatorname{tr}M^2=\|M\|_{\mathrm{HS}}^2\), not with
\(\operatorname{tr}\operatorname{Cov}X\).

Because \(\mathbb Eh=0\),
\[
\begin{aligned}
 |\langle M,D_A-D_B\rangle_{\mathrm{HS}}|
 &=|\mathbb E[hX^{\mathsf T}MX]|\\
 &=|\mathbb E[h\{X^{\mathsf T}MX-\mathbb E(X^{\mathsf T}MX)\}]|\\
 &\le\sqrt\rho\,
       \sqrt{\operatorname{Var}(X^{\mathsf T}MX)}
 \le2\sqrt\rho.                                                 \tag{3.2}
\end{aligned}
\]
The matrix \(D_A-D_B\) is symmetric, so Hilbert--Schmidt duality over
symmetric \(M\) gives
\[
 \|D_A-D_B\|_{\mathrm{HS}}\le2\sqrt\rho.                        \tag{3.3}
\]

For the operator bound on \(D_B\), take
\(M=\theta\theta^{\mathsf T}\).  Now
\[
 \mathbb E(\mathbf1_B-g)^2=g(1-g)\le\frac14.
\]
Consequently (3.1) gives
\[
 |\theta^{\mathsf T}D_B\theta|
 \le\sqrt{g(1-g)}\,
       \sqrt{\operatorname{Var}\langle\theta,X\rangle^2}
 \le1,                                                         \tag{3.4}
\]
and hence
\[
 \|D_B\|_{\mathrm{op}}\le1.                                    \tag{3.5}
\]

The indicator \(B\) need not be measurable with respect to \(X\).  This
causes no problem: Poincare is applied only to the quadratic function of
\(X\), while Cauchy--Schwarz is taken in the joint Gaussian-source
probability space.

For a hard support, truncate the quadratic, use the closed Poincare form,
and pass by the subgaussian moment bound.  No boundary term or ambient trace
appears.

---

## 4. Composition with the halfspace-pullback theorem

The halfspace-pullback theorem applies to the source halfspace \(B\), the
same centered positive Brenier contraction \(T\), and
\[
 v_B=\mathbb E[(\mathbf1_B-g)X],\quad
 D_B=\mathbb E[(\mathbf1_B-g)XX^{\mathsf T}].
\]
It supplies
\[
 \|P_BD_B\|_{\mathrm{HS}}^2
 \le C_\delta\{I(g)-|v_B|\}.                                   \tag{4.1}
\]
Since \(v_A=v=(1-\varepsilon)I(g)u\), (2.2) gives
\[
 I(g)-|v_B|
 \le I(g)\varepsilon+\sqrt\rho.                                \tag{4.2}
\]
For sufficiently small \(\varepsilon\), centrality and (2.2) also give
\[
 |u-u_B|\le C_\delta\sqrt\rho.                                 \tag{4.3}
\]

For rank-one projectors,
\[
 \|P-P_B\|_{\mathrm{HS}}
 =\|uu^{\mathsf T}-u_Bu_B^{\mathsf T}\|_{\mathrm{HS}}
 \le\sqrt2\,|u-u_B|.                                           \tag{4.4}
\]
Therefore (3.3), (3.5), and (4.1)--(4.4) imply
\[
\begin{aligned}
 \|PD_A\|_{\mathrm{HS}}
 &\le\|D_A-D_B\|_{\mathrm{HS}}
 +\|(P-P_B)D_B\|_{\mathrm{HS}}
 +\|P_BD_B\|_{\mathrm{HS}}\\
 &\le C_\delta\left\{
 \sqrt\rho+(\varepsilon+\sqrt\rho)^{1/2}\right\}.               \tag{4.5}
\end{aligned}
\]
After squaring and using \(\rho\le\sqrt\rho\) in the small-defect regime,
\[
 \|PD_A\|_{\mathrm{HS}}^2
 \le C_\delta\{\varepsilon+\sqrt\rho\}.                         \tag{4.6}
\]

For large \(\varepsilon\), direct Hilbert--Schmidt duality with
\(\mathbf1_A-g\) and (3.1) gives \(\|D_A\|_{\mathrm{HS}}\le1\), so enlarging
\(C_\delta\) closes the full interval.

Finally, replacing \(X\) by \(\widetilde X=\sqrt t(X-m)\) gives
\(\widetilde D=tD\) without changing \(P\).  Thus (4.6) restores the factor
\(t^{-2}\), exactly as claimed.

---

## 5. Audit conclusion

Both potential dimension leaks are absent:

* threshold disagreement is one-dimensional Gaussian anti-concentration;
* the quadratic perturbation uses Poincare plus Hilbert--Schmidt duality,
  pairing \(\operatorname{Cov}X\preceq I\) with
  \(\operatorname{tr}M^2=1\).

The rank-two projector perturbation is also dimension free.  Subject to the
already proved halfspace-pullback theorem and the clean-room spatial
stability lemmas, the general angular-stability result is load-bearing.
