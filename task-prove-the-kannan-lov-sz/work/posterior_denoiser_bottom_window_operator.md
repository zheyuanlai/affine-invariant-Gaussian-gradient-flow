# Denoiser operators on a bottom spectral window

## 0. Purpose and verdict

This note records a dimension-free operator calculation for a fixed Gaussian
output. It uses the posterior covariance estimate
\[
 A_y=\operatorname {Cov}(X\mid Y=y)\preceq I,\qquad Dm(y)=A_y,
\]
and the posterior Hilbert--Schmidt third-moment estimate
\[
 \left(\sum_j\|(DA_y[e_j])u\|^2\right)^{1/2}
 \le 4\sqrt{u^T(I-A_y)u}.                         \tag{0.1}
\]

The map \(W\mapsto A_yW\) has no trace loss: if \(W\) is the gradient of a
bottom spectral-window vector, then its derivative remains small after this
map. A bounded inverse surrogate
\[
 \mathcal C_y=2(I+A_y)^{-1}
\]
has the same property and gives a canonical posterior conditional lift.
These estimates remove \(C_0\) from a *backward transfer* of a low output
mode to the input, but they do not give the affine--orthogonal ANOVA
coercivity. The exact minimal Stein kernel of the output still contains an
inverse-generator term, whose directional row covariance is precisely the
unresolved bottom-spectrum quantity. Thus the calculation does not prove
KLS; it isolates why the HS saturation theorem alone cannot remove \(C_0\).

All statements below are intrinsic to the supporting affine space. We use
the unscaled channel
\[
 Y=X+G,\qquad X\sim\mu,\quad G\sim N(0,I),
\]
so that \(\operatorname {Cov}(Y)=2I\). The normalized output is
\(S=Y/\sqrt2\); its spectral gap is twice the unscaled one.

## 1. Posterior notation and a directional row estimate

Write \(m(y)=E[X\mid Y=y]\), \(A(y)=\operatorname {Cov}(X\mid Y=y)\), and
\[
 H(y)=I-A(y)=D^2U(y),\qquad d\nu(y)=e^{-U(y)}dy.
\]
Gaussian convolution gives a positive analytic density. Posterior
Brascamp--Lieb gives \(0\preceq A\preceq I\), and differentiation of the
exponential family gives \(Dm=A\). Let \(W:\mathbb R^n\to\mathbb R^n\) be
locally Sobolev. The product rule and (0.1) imply
\[
 D(AW)[v]=(DA[v])W+A\,DW[v].                      \tag{1.1}
\]
Therefore, with the Hilbert--Schmidt norm taken over the \(v\)-slot,
\[
 \boxed{\ \|D(AW)\|_{\rm HS}
 \le \|DW\|_{\rm HS}+4\sqrt{W^THW}.\ }              \tag{1.2}
\]
Indeed \(\|A DW\|_{\rm HS}\le\|DW\|_{\rm HS}\), while (0.1) controls the
other summand. In particular,
\[
 E\|D(AW)\|_{\rm HS}^2
 \le 17\,E\bigl(\|DW\|_{\rm HS}^2+W^THW\bigr).       \tag{1.3}
\]
The numerical constant follows from \((a+4b)^2\le17(a^2+b^2)\). Also
\[
 E|AW-W|^2=E|HW|^2\le E[W^THW].                    \tag{1.4}
\]
Equations (1.2)--(1.4) are the promised directional row-covariance
control. They are valid for an adaptively selected \(W(y)\); no summation
over directions and no \(n\)-dependent trace estimate occurs.

### 1.1 A bounded inverse surrogate

Put \(B=(I+A)^{-1}\) and \(C=2BW\). Since \(A,H,B\) commute pointwise,
\[
 W-AC=HBW,\qquad
 DC[v]=2B\,DW[v]-2B\,(DA[v])BW.                    \tag{1.5}
\]
Applying (0.1) to \(BW\), and using \(0\preceq B\preceq I\), gives
\[
 \boxed{\ \|DC\|_{\rm HS}
 \le 2\|DW\|_{\rm HS}+8\sqrt{W^THW}.\ }              \tag{1.6}
\]
and hence
\[
 E\|DC\|_{\rm HS}^2
 \le68\,E\bigl(\|DW\|_{\rm HS}^2+W^THW\bigr).       \tag{1.7}
\]
The point of \(C=2(I+A)^{-1}W\) is that it is uniformly bounded even when
\(A\) has small eigenvalues, while \(AC\approx W\) on a direction in which
\(H=I-A\) is small. Formula (1.6) is a dimension-free row estimate for
this regularized posterior inverse.

## 2. Exact bottom-window bookkeeping

Let \(\mathcal A=-L_\nu\), \(L_\nu=\Delta-\nabla U\cdot\nabla\), and let
\(f\) be centered, \(\|f\|_{L^2(\nu)}=1\). Set
\[
 q=E|\nabla f|^2,\quad W=\nabla f,\quad K=D^2f,\quad
 B_f=E\|K\|_{\rm HS}^2,\quad C_f=E[W^THW].
\]
For an exact eigenfunction \(\mathcal Af=\lambda f\), the weighted Bochner
identity gives
\[
 q=\lambda,\qquad B_f+C_f=\lambda^2.              \tag{2.1}
\]
For a unit vector in the spectral window \([b,b+\varepsilon]\), the same
identities hold with
\[
 q\le b+\varepsilon,\qquad B_f+C_f
 =\|\mathcal Af\|_2^2\le(b+\varepsilon)^2,           \tag{2.2}
\]
after the standard closed-form approximation. Since the coordinate test has
Rayleigh quotient \(1/2\) for \(Y\), one always has \(b\le1/2\).

For the normalized field \(\widehat W=W/\sqrt q\), (1.3) and (2.2) read
\[
 E\|D(A\widehat W)\|_{\rm HS}^2
 \le 17\,\frac{(b+\varepsilon)^2}{q}
 \le17\,\frac{(b+\varepsilon)^2}{b}
\]
when \(b>0\), and
\[
 E|A\widehat W-\widehat W|^2
 \le\frac{(b+\varepsilon)^2}{q}.
                                                               \tag{2.3}
\]
In the usual bottom-window choice \(q\ge b\) and
\(\varepsilon=o(b)\), these are \(O(b)\). Thus denoising preserves the
slowly varying, nearly saturated structure of a bottom field with a
universal constant. This remains valid when the saturating direction
rotates with \(y\).

## 3. A posterior conditional lift (no input Poincare constant)

Let \(Z=X-m(Y)\), so \(E[Z\mid Y]=0\) and
\(E[ZZ^T\mid Y]=A(Y)\). Define the first-order lift
\[
 \Phi(X,G)=f(Y)+W(Y)\cdot Z.                         \tag{3.1}
\]
For smooth \(f\), direct differentiation gives
\[
 \nabla_G\Phi=H W+KZ,\qquad
 \nabla_X\Phi=(I+H)W+KZ.                             \tag{3.2}
\]
Conditioning on \(Y\), the cross terms with \(Z\) vanish. Since
\(H^2\preceq H\), \(A\preceq I\), and \((I+H)^2\preceq I+3H\),
\[
 E|\nabla_G\Phi|^2\le C_f+B_f,\qquad
 E|\nabla_X\Phi|^2\le q+3C_f+B_f.                  \tag{3.3}
\]
Moreover \(E[\Phi\mid Y]=f(Y)\) and
\[
 \|\Phi\|_2^2=1+q-C_f.                              \tag{3.4}
\]
Let \(u(X)=E_G[\Phi(X,G)\mid X]\). Conditional Gaussian Poincare and
(3.3) give
\[
 \|\Phi-u\|_2^2\le B_f+C_f,\qquad
 \|u\|_2^2\ge1+q-2(B_f+C_f),\qquad
 E|\nabla u|^2\le q+3(B_f+C_f).                    \tag{3.5}
\]
Consequently, whenever the denominator is positive,
\[
 \boxed{\ C_P(\mu)^{-1}\le
 {q+3(B_f+C_f)\over1+q-2(B_f+C_f)}.\ }              \tag{3.6}
\]
For an exact bottom eigenfunction this is
\[
 C_P(\mu)^{-1}\le
 {\,\lambda+3\lambda^2\,\over1+\lambda-2\lambda^2}. \tag{3.7}
\]
The inequality is an *upper* comparison for the input gap (a low output
mode produces a comparably low input mode). It uses no input Poincare
constant and is valid on the closed form domain by approximation.

### 3.1 HS-enhanced lift

To make explicit where the posterior HS theorem enters, replace \(W\) in
(3.1) by \(C=2(I+A)^{-1}W\):
\[
 \Phi_C=f(Y)+C(Y)\cdot Z.                            \tag{3.8}
\]
Set \(e=W-AC=HBW\). From (1.6), conditioning on \(Y\), and Gaussian
Poincare,
\[
 E|\nabla_G\Phi_C|^2
 \le E|e|^2+E\|DC\|_{\rm HS}^2
 \le69(B_f+C_f),                                    \tag{3.9}
\]
\[
 \|\Phi_C\|_2^2=1+E[C^TAC]\ge1+q-C_f.               \tag{3.10}
\]
The last inequality uses the scalar bound
\(4a/(1+a)^2\ge a\) for \(0\le a\le1\). Its \(X\)-gradient obeys
\[
 E|\nabla_X\Phi_C|^2\le q+76(B_f+C_f),               \tag{3.11}
\]
because \(W+HC=\frac{2I+H}{2I-H}W\) and
\((2+h)^2/(2-h)^2\le1+8h\) for \(0\le h\le1\).
Thus \(u_C=E_G[\Phi_C\mid X]\) satisfies the explicit estimate
\[
 C_P(\mu)^{-1}\le
 {q+76(B_f+C_f)\over1+q-70(B_f+C_f)}.               \tag{3.12}
\]
The constants are worse than (3.6), but (3.12) is the version whose
conditional operator is controlled directly by the audited HS saturation
estimate. For a bottom window, \(B_f+C_f=O(b^2)\), so it still gives
\(C_P(\mu)^{-1}\le b/2+O(b^2)\) after the \(S=Y/\sqrt2\) rescaling.

## 4. The canonical output Stein kernel and the surviving obstruction

The preceding lifts should not be confused with a Stein kernel. Let
\(m_a=a\cdot m\), and let \(u_a\) be the minimal weak solution of
\[
 \mathcal A u_a=m_a,\qquad E u_a=0.
\]
Define the matrix \(K_\nu\) by its rows \(\nabla u_{e_i}\). Then
\[
 E[m_a\,\varphi(Y)]
 =E[K_\nu^Ta\cdot\nabla\varphi(Y)]                \tag{4.1}
\]
for every form-domain test \(\varphi\). Since
\[
 E[Y\varphi]=E[\nabla\varphi]+E[m\varphi],       \tag{4.2}
\]
the matrix
\[
 \tau_Y=I+K_\nu                                      \tag{4.3}
\]
is a Stein kernel for \(Y\) (its expectation is \(2I\)). This is the
canonical minimal kernel; it is not a pointwise function of \(A_y\).

There is an exact spectral formula for its directional row covariance. If
\(\sigma_a\) is the spectral measure of the coordinate \(Y_a=a\cdot Y\)
for \(\mathcal A\), then
\[
 \|Y_a\|_2^2=2|a|^2,\qquad
 \langle Y_a,\mathcal A Y_a\rangle=|a|^2,
\]
and \(m_a=(I-\mathcal A)Y_a\). Hence
\[
 E|K_\nu^Ta|^2
 =\int {(1-t)^2\over t}\,d\sigma_a(t)
 =\langle Y_a,\mathcal A^{-1}Y_a\rangle-3|a|^2,   \tag{4.4}
\]
\[
 \boxed{\ E|(K_\nu-I)^Ta|^2
 =\langle Y_a,\mathcal A^{-1}Y_a\rangle-4|a|^2.\ }   \tag{4.5}
\]
The right side is the directional \(H^{-1}(\nu)\) norm of the coordinate,
minus its Gaussian value. It is controlled by
\(C_P(\nu)=1/\lambda_1(\nu)\), not by the pointwise bounds
\(0\preceq A\preceq I\) or by (0.1). In particular, the HS saturation
estimate controls \(D(AW)\) in a selected direction, whereas (4.4)
involves the inverse of the full scalar generator at the bottom of its
spectrum.

For an exact eigenfunction \(\mathcal Af=\lambda f\), writing
\(\ell_a=E[Y_a f]\), (4.1) gives
\[
 E[K_\nu^Ta\cdot\nabla f]=(1-\lambda)\ell_a,\qquad
 E[a\cdot\nabla f]=\lambda\ell_a,                  \tag{4.6}
\]
so
\[
 E[(K_\nu-I)^Ta\cdot\nabla f]=(1-2\lambda)\ell_a.   \tag{4.7}
\]
On an affine-orthogonal mode this pairing vanishes, but (4.5) still has no
dimension-free bound. On a bottom spectral window the error terms in
(4.6)--(4.7) are \(O(\varepsilon)\), while the row covariance remains an
inverse-generator quantity. Therefore the posterior operator estimates
above do not remove the \(C_0\) term in the matrix-Stein ANOVA criterion on
the nonlinear survivor; they only provide a sharp, dimension-free backward
transfer and a controlled denoiser derivative.

## 5. Intrinsic approximation

All displayed identities are first proved for smooth \(f\) and a smooth
positive input density. Gaussian convolution makes the output and all
posterior moments analytic. Truncation in \(X\) and \(Y\), followed by
closed-form convergence, gives (2.1)--(3.12) for spectral-window vectors.
If the original measure is supported on a proper affine subspace, perform
the construction inside that subspace and take the Gaussian noise there;
no ambient null direction is included. The constants do not depend on the
intrinsic dimension.

