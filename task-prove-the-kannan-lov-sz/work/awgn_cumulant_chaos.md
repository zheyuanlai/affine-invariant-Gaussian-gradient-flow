# Low-SNR AWGN chaos with the unconditional output denominator

## 0. Verdict

Let \(X\) be isotropic log-concave in \(\mathbb R^n\), let
\(Z\in\{-1,1\}\) be a balanced Borel function of \(X\), and observe

\[
Y_t=\sqrt t\,X+G,\qquad G\sim N(0,I_n).
\]

The correct binary chi-square quantity is

\[
\mathcal D(t)
=\mathbb E\bigl[\mathbb E(Z\mid Y_t)^2\bigr]
=\int_{\mathbb R^n}\gamma(y)\frac{N_t(y)^2}{L_t(y)}\,dy,
\tag{0.1}
\]

where

\[
L_t(y)=\mathbb E e^{\sqrt t\langle y,X\rangle-t|X|^2/2},
\qquad
N_t(y)=\mathbb E Z
e^{\sqrt t\langle y,X\rangle-t|X|^2/2}.
\tag{0.2}
\]

The unconditional likelihood \(L_t\), rather than the standard Gaussian
density alone, is essential.  If \(X=(X_{\mathrm{sig}},X_{\mathrm{nuis}})\)
and \(X_{\mathrm{nuis}}\) is independent of \((Z,X_{\mathrm{sig}})\), then
its likelihood factor cancels exactly from (0.1).  Thus (0.1) has no false
dimension amplification from independent nuisance coordinates.

Put

\[
a=\mathbb E[ZX],\qquad B=\mathbb E[ZXX^T],
\qquad C=\mathbb E[ZX^{\otimes3}],\qquad R=\mathbb EX^{\otimes3},
\]

and use the averaged symmetrization convention

\[
\operatorname{Sym}(a\otimes B)_{ijk}
=\frac13(a_iB_{jk}+a_jB_{ik}+a_kB_{ij}).
\]

Set

\[
K_3=C-3\operatorname{Sym}(a\otimes I).
\tag{0.3}
\]

The exact expansion through third order is

\[
\boxed{
\begin{aligned}
\mathcal D(t)
={}&t|a|^2
+t^2\left(\frac12\|B\|_{\mathrm{HS}}^2-|a|^2\right)\\
&+t^3\left(
|a|^2-\|B\|_{\mathrm{HS}}^2
+\frac16\|K_3\|_{\mathrm{HS}}^2
-\langle\operatorname{Sym}(a\otimes B),R\rangle
\right)
+O(t^4).
\end{aligned}}
\tag{0.4}
\]

All odd powers of \(\sqrt t\) vanish.  Section 3 derives (0.4) directly
from Hermite orthogonality.

The degree-one tensor is universally bounded, \(|a|\le1\).  Under the
general quadratic-form variance theorem

\[
\operatorname{Var}(X^TQX)\le C\|Q\|_{\mathrm{HS}}^2
\tag{0.5}
\]

for isotropic log-concave \(X\), duality gives
\(\|B\|_{\mathrm{HS}}\le C\).  Thus the first two orders are
dimension-free.

There are two rigorous obstructions at the next stage.

1. **Individual class-conditional degree-three cumulants have no
   dimension-free Hilbert--Schmidt bound.**  For the isotropic regular
   simplex in dimension \(n\), the third cumulant tensor has norm

   \[
   \frac{2\sqrt{n(n-1)(n+2)}}{n+3}\asymp\sqrt n.
   \tag{0.6}
   \]

   Taking a product with one independent Gaussian coordinate and letting
   \(Z\) be the sign of that coordinate gives a smoothable isotropic
   log-concave counterexample: both conditional laws contain the same
   simplex nuisance tensor of size \(\sqrt n\).  The tensor cancels from
   the *contrast* and from (0.1).  Consequently, bounding the two
   class-conditional cumulants separately is mathematically false and
   recreates precisely the nuisance-coordinate error that (0.1) avoids.

2. **The Hermite expansion has no positive uniform analytic radius under
   log-concavity.**  Already in one dimension, for a centered exponential
   variable and its balanced median-threshold label, the signed moments
   satisfy \(|\mathbb E[ZX^k]|\ge c\,k!\).  The \(k\)-th Hermite term of
   \(N_t\) has squared \(L^2(\gamma)\)-norm at least \(c\,k!\).
   Therefore the standard-Gaussian \(L^2\) chaos series has radius zero.
   The denominator in (0.1) may resum these terms, but no proof at a
   numerical \(t\) can come from absolute termwise chaos bounds.

The only potentially meaningful degree-three target is a bound on the
**relative contrast** \(K_3\), together with its coupling to the
unconditional skew tensor \(R\) in (0.4).  Neither a bound on the
individual conditional cumulants nor any finite truncation gives the
desired fixed-SNR contraction.  A successful low-SNR proof would require
a nonperturbative resummation of (0.1), uniform over exponential tails;
the polynomial-chaos route by itself stops here.

---

## 1. Binary testing and the correct divergence

Let \(p_\pm^t\) be the density of \(Y_t\) conditional on \(Z=\pm1\), and
let

\[
p_t=\frac12(p_+^t+p_-^t).
\]

The posterior magnetization is

\[
m_t(y)=\mathbb E[Z\mid Y_t=y]
=\frac{p_+^t(y)-p_-^t(y)}{p_+^t(y)+p_-^t(y)}.
\]

Define

\[
\mathcal D(t)=\int p_t(y)m_t(y)^2\,dy
=\frac12\int
\frac{(p_+^t-p_-^t)^2}{p_+^t+p_-^t}\,dy.
\tag{1.1}
\]

This is half the triangular discrimination.  It is also the explained
variance of the balanced binary variable.  Hence

\[
\operatorname{mmse}(Z\mid Y_t)=1-\mathcal D(t).
\tag{1.2}
\]

If

\[
\operatorname{err}(t)
=\mathbb E\min\{\mathbb P(Z=1\mid Y_t),
\mathbb P(Z=-1\mid Y_t)\},
\]

then the elementary inequalities

\[
p(1-p)\le\min(p,1-p)\le2p(1-p)
\]

give

\[
\frac{1-\mathcal D(t)}4
\le\operatorname{err}(t)
\le\frac{1-\mathcal D(t)}2.
\tag{1.3}
\]

Thus a numerical upper bound \(\mathcal D(t)\le1-c\) is sufficient for a
universal Bayes-error lower bound.

Dividing the Gaussian channel densities by \(\gamma\) gives (0.2), and

\[
p_t(y)=\gamma(y)L_t(y),\qquad
\frac12(p_+^t-p_-^t)=\gamma(y)N_t(y).
\]

Substitution in (1.1) gives (0.1).

### Exact nuisance cancellation

Suppose \(X=(S,U)\), where \(U\) is independent of \((Z,S)\).  Then

\[
L_t(y_S,y_U)=L_t^S(y_S)L_t^U(y_U),\qquad
N_t(y_S,y_U)=N_t^S(y_S)L_t^U(y_U).
\]

Therefore

\[
\begin{aligned}
\mathcal D(t)
&=\int\gamma_S\gamma_U
\frac{(N_t^S)^2(L_t^U)^2}{L_t^SL_t^U}\\
&=\left(\int\gamma_S\frac{(N_t^S)^2}{L_t^S}\right)
\left(\int\gamma_U L_t^U\right)
=\mathcal D_S(t),
\end{aligned}
\tag{1.4}
\]

because \(\int\gamma_U L_t^U=1\).  This identity is the audit that every
tensor expansion must pass.

---

## 2. Hermite tensors

Let \(H_k(y)\) be the probabilists' symmetric Hermite tensor, normalized by

\[
e^{\varepsilon\langle y,x\rangle-\varepsilon^2|x|^2/2}
=\sum_{k=0}^{\infty}
\frac{\varepsilon^k}{k!}H_k(y):x^{\otimes k}.
\tag{2.1}
\]

For symmetric \(k\)-tensors \(S,T\),

\[
\int\gamma(y)(H_k(y):S)(H_\ell(y):T)\,dy
=\mathbf 1_{\{k=\ell\}}\,k!\,\langle S,T\rangle.
\tag{2.2}
\]

Put \(\varepsilon=\sqrt t\), and define

\[
n_k(y)=\frac1{k!}H_k(y):\mathbb E[ZX^{\otimes k}],
\]

\[
\ell_k(y)=\frac1{k!}H_k(y):\mathbb E[X^{\otimes k}].
\]

Since \(\mathbb EZ=0\), \(\mathbb EX=0\), and
\(\mathbb EXX^T=I\),

\[
N_t=\varepsilon n_1+\varepsilon^2n_2+\varepsilon^3n_3+\cdots,
\tag{2.3}
\]

\[
L_t=1+\varepsilon^2\ell_2+\varepsilon^3\ell_3+\cdots,
\qquad
\ell_2(y)=\frac12(|y|^2-n).
\tag{2.4}
\]

Although (2.3)--(2.4) always give finite-order derivatives at zero, they
need not converge for any nonzero \(\varepsilon\); Section 6 gives an
explicit example.

---

## 3. Coefficients through order \(t^3\)

Expanding \(N_t^2/L_t\) through \(\varepsilon^6\), using
\(\ell_1=0\), gives

\[
\begin{aligned}
[t]\mathcal D&=\int n_1^2\,d\gamma,\\
[t^2]\mathcal D&=\int(n_2^2-\ell_2n_1^2)\,d\gamma,\\
[t^3]\mathcal D&=\int\bigl(
n_3^2-\ell_2n_2^2-2\ell_2n_1n_3
-2\ell_3n_1n_2+\ell_2^2n_1^2
\bigr)\,d\gamma.
\end{aligned}
\tag{3.1}
\]

Terms \(2n_2n_4\), \(2n_1n_5\), and
\(-\ell_4n_1^2\) integrate to zero by Hermite orthogonality.

The first two lines immediately give

\[
[t]\mathcal D=|a|^2,
\tag{3.2}
\]

\[
[t^2]\mathcal D=\frac12\|B\|_{\mathrm{HS}}^2-|a|^2.
\tag{3.3}
\]

For the third line, standard Gaussian contractions give

\[
\int n_3^2\,d\gamma=\frac16\|C\|_{\mathrm{HS}}^2,
\]

\[
\int\ell_2n_2^2\,d\gamma=\|B\|_{\mathrm{HS}}^2,
\]

\[
2\int\ell_2n_1n_3\,d\gamma
=\langle a,\operatorname{tr}_{12}C\rangle,
\]

\[
2\int\ell_3n_1n_2\,d\gamma
=\langle\operatorname{Sym}(a\otimes B),R\rangle,
\]

\[
\int\ell_2^2n_1^2\,d\gamma
=\left(\frac n2+2\right)|a|^2.
\tag{3.4}
\]

Using

\[
\left\|\operatorname{Sym}(a\otimes I)\right\|_{\mathrm{HS}}^2
=\frac{n+2}{3}|a|^2,
\qquad
\operatorname{tr}_{12}\operatorname{Sym}(a\otimes I)
=\frac{n+2}{3}a,
\]

the dimension-dependent pieces in (3.4) cancel exactly, yielding the
third line of (0.4).

This cancellation is the first nontrivial appearance of the unconditional
denominator.  Dropping \(L_t\) would leave an artificial term of order
\(n|a|^2\).

### Relation to class-conditional cumulants

Let \(\kappa_3^\pm\) denote the third central cumulant tensor of
\(X\mid Z=\pm1\), and put

\[
\Delta\kappa_3=\frac12(\kappa_3^+-\kappa_3^-).
\]

Because the conditional means are \(a\) and \(-a\), a direct expansion
gives

\[
\Delta\kappa_3=K_3+2a^{\otimes3}.
\tag{3.5}
\]

Thus \(K_3=\Delta\kappa_3-2a^{\otimes3}\) is the relative third-order
contrast that actually enters (0.4).  Individual values of
\(\kappa_3^\pm\) contain common nuisance skewness and need not be bounded.

---

## 4. Dimension-free control of degrees one and two

For every unit vector \(\theta\),

\[
|\langle a,\theta\rangle|
=|\mathbb E[Z\langle X,\theta\rangle]|
\le\sqrt{\mathbb E\langle X,\theta\rangle^2}=1.
\]

Taking the supremum gives

\[
|a|\le1.
\tag{4.1}
\]

For a symmetric matrix \(Q\), balancedness and isotropy give

\[
\langle B,Q\rangle
=\mathbb E\left[
Z(X^TQX-\operatorname{tr}Q)\right].
\]

Hence (0.5) and Cauchy--Schwarz imply

\[
|\langle B,Q\rangle|\le C\|Q\|_{\mathrm{HS}}.
\]

Duality gives

\[
\|B\|_{\mathrm{HS}}\le C.
\tag{4.2}
\]

Equations (4.1)--(4.2) make (3.2)--(3.3) dimension-free.  They do not
control the relative cubic tensor in (3.5).

---

## 5. A degree-three class-cumulant counterexample with exact nuisance cancellation

Let \(w_1,\ldots,w_{n+1}\in\mathbb R^n\) be the vertices of a regular
simplex normalized by

\[
\langle w_i,w_j\rangle=
\begin{cases}
n,&i=j,\\
-1,&i\ne j.
\end{cases}
\]

Let \(P=(P_1,\ldots,P_{n+1})\) be uniform on the probability simplex,
that is, Dirichlet with all parameters equal to one.  Then

\[
X_{\mathrm{simp}}
=\sqrt{n+2}\sum_{i=1}^{n+1}P_iw_i
\tag{5.1}
\]

is the isotropic uniform measure on a regular simplex.

Put

\[
S_3=\sum_{i=1}^{n+1}w_i^{\otimes3}.
\]

The Dirichlet third-moment formula and \(\sum_iw_i=0\) give

\[
\kappa_3(X_{\mathrm{simp}})
=\mathbb EX_{\mathrm{simp}}^{\otimes3}
=\frac{2\sqrt{n+2}}{(n+1)(n+3)}S_3.
\tag{5.2}
\]

Moreover,

\[
\begin{aligned}
\|S_3\|_{\mathrm{HS}}^2
&=\sum_{i,j}\langle w_i,w_j\rangle^3\\
&=(n+1)n^3-(n+1)n\\
&=n(n-1)(n+1)^2.
\end{aligned}
\tag{5.3}
\]

Combining (5.2)--(5.3) gives the exact formula

\[
\left\|\kappa_3(X_{\mathrm{simp}})\right\|_{\mathrm{HS}}
=\frac{2\sqrt{n(n-1)(n+2)}}{n+3}
\asymp\sqrt n.
\tag{5.4}
\]

Now take the product of (5.1) with an independent standard Gaussian
\(W\), and let

\[
Z=\operatorname{sign}W.
\]

The product is isotropic log-concave and \(Z\) is balanced.  Conditional
on either sign of \(W\), the simplex factor is unchanged.  Consequently
both class-conditional third cumulant tensors contain the tensor (5.2),
whose norm diverges as \(\sqrt n\).  If a positive \(C^\infty\) density is
desired, convolve the simplex factor with an arbitrarily small Gaussian
and whiten it; (5.4) persists up to a fixed factor by moment continuity.

This is an explicit counterexample to every claim that
\(\|\kappa_3^\pm\|_{\mathrm{HS}}\) is universally bounded.  It is not a
counterexample to AWGN contraction.  Indeed, (1.4) cancels the simplex
factor exactly, and

\[
\Delta\kappa_3
\]

has zero simplex component.  The example pinpoints why the expansion must
be organized in relative contrasts such as (3.5), never in separate
class tensors.

---

## 6. Zero analytic radius for exponential tails

Let

\[
X=E-1,\qquad E\sim\operatorname{Exp}(1).
\]

Then \(X\) is centered, has variance one, and has a one-dimensional
log-concave density.  Let

\[
Z=
\begin{cases}
+1,&E\ge\log2,\\
-1,&E<\log2.
\end{cases}
\tag{6.1}
\]

This is balanced.  For every sufficiently large integer \(k\), the
positive exponential tail dominates the bounded interval
\(-1\le X\le\log2-1\), and elementary gamma integration gives

\[
|\mathbb E[ZX^k]|\ge c\,k!.
\tag{6.2}
\]

The \(k\)-th term in the Hermite expansion of \(N_t\) is

\[
\frac{t^{k/2}}{k!}\mathbb E[ZX^k]H_k(y).
\]

Since \(\|H_k\|_{L^2(\gamma)}^2=k!\), its squared norm after removing the
factor \(t^{k/2}\) is

\[
\frac{|\mathbb E[ZX^k]|^2}{k!}\ge c\,k!.
\tag{6.3}
\]

For every \(t>0\), the norms of the individual terms
\(t^{k/2}\sqrt{k!}\) fail to tend to zero.  Hence the Hermite series for
\(N_t\) has radius zero in \(L^2(\gamma)\).

The exact ratio \(N_t^2/L_t\) in (0.1) is finite for every \(t>0\);
the Gaussian factor \(\exp(-tX^2/2)\) performs a nonperturbative
resummation.  Formula (6.3) proves that this resummation cannot be
replaced by absolute control of polynomial-chaos degrees.  In particular,
dimension-free estimates of any fixed collection of degrees, even if
extended through degrees three and four, do not provide a finite-radius
bound at a numerical signal-to-noise ratio.

---

## 7. Final assessment

The unconditional-output binary chi-square (0.1) is the right low-SNR
object.  It passes exact nuisance cancellation and its first two
coefficients are controlled by isotropy and the quadratic-form variance
theorem.  At third order it exposes the relative tensor \(K_3\), not the
two class cumulants separately.

The regular-simplex product proves that separate degree-three tensor
bounds are false by a factor \(\sqrt n\), while their common component is
irrelevant to Bayes error.  The centered-exponential example proves that
the full Hermite expansion has no positive analytic radius even in one
dimension.  Therefore a finite-order polynomial-chaos bootstrap cannot
establish a universal AWGN contraction at fixed \(t\).

A viable continuation would have to estimate the ratio in (0.1)
nonperturbatively--for example by a comparison principle for the
posterior likelihood ratio--while preserving the exact product
cancellation (1.4).  Merely adding degree-three or degree-four moment
bounds cannot close the argument.
