# Angular stability for a posterior that splits in the active direction

## 1. Statement

Fix `t>0`.  Let `rho` be a centered `t`-strongly log-concave
probability on `R`, let `nu` be a centered `t`-strongly log-concave
probability on `R^m`, and let

\[
 \pi=\rho\otimes\nu
\]

on `R times R^m`.  Write a point as `X=(Y,W)`.  Let `S` be Borel and
put

\[
 g=\pi(S),\qquad v=\mathbb E[(\mathbf1_S-g)X],\qquad
 D=\mathbb E[(\mathbf1_S-g)XX^T].                    \tag{1.1}
\]

Assume that `v=(|v|,0)` with `|v|>0`; in other words, the product
coordinate is the active centroid direction.  Let `P` be projection onto
the `W` coordinates and let

\[
 \Delta={\mathcal I(g)\over\sqrt t}-|v|.             \tag{1.2}
\]

**Theorem 1 (split-posterior angular stability).**  For every
`delta in (0,1/2)` there is a finite numerical `C_delta` such that, if

\[
 \delta\le g\le1-\delta,
\]

then

\[
 \boxed{\|PD\|_{HS}^2\le C_\delta t^{-3/2}\Delta.}  \tag{1.3}
\]

No Gaussian hypothesis is imposed on either factor.  The result therefore
identifies the missing issue in the nonsplit case precisely: it is not
non-Gaussianity of the active marginal or of the transverse law, but
dependence between them.

## 2. The two scalar defects

Let `c` be the upper `g`-quantile of `rho` and put

\[
 H(y)=\mathbf1_{\{y\ge c\}},\qquad
 \sigma(y,w)=H(y)-\mathbf1_S(y,w).                   \tag{2.1}
\]

For `w in R^m`, define

\[
 m_0(w)=\int\sigma(y,w)d\rho(y),\qquad
 m_1(w)=\int(y-c)\sigma(y,w)d\rho(y).                \tag{2.2}
\]

The sign of `sigma` is the sign of `y-c`, so `m_1>=0`.  Equality of the
masses and the assumption that `v` has no transverse component give

\[
 \mathbb E_\nu m_0=0,\qquad
 \mathbb E_\nu[Wm_0(W)]=0.                           \tag{2.3}
\]

Let `T` be the increasing transport from a standard Gaussian `Z` to
`rho`, put `L=t^{-1/2}`, `z_g=Phi^{-1}(1-g)`, and

\[
 R(z)=T(z)-Lz.
\]

One-dimensional Caffarelli contraction says that `T` is `L`-Lipschitz,
hence `R` is nonincreasing.  After changing `T` on a null set we may take
`T(z_g)=c`.  The exact centroid decomposition is

\[
 \Delta=D_{\rm cut}+D_{\rm map},                    \tag{2.4}
\]

where

\[
 D_{\rm cut}=\mathbb E_\nu m_1(W),\qquad
 D_{\rm map}=g(1-g)
 \{\mathbb E[R\mid Z<z_g]-\mathbb E[R\mid Z\ge z_g]\}. \tag{2.5}
\]

Both terms are nonnegative.

## 3. Transport-weighted flip cost

The point that permits a non-Gaussian active marginal is the following
estimate.

**Lemma 2.**  Under the hypotheses above,

\[
 \boxed{\mathbb E_\nu m_0(W)^2
 \le C_\delta\sqrt t\,(D_{\rm cut}+D_{\rm map}).}    \tag{3.1}
\]

**Proof.**  For fixed `w`, put

\[
 a_w(z)=|\sigma(T(z),w)|,qquad A_w=\mathbb E a_w(Z).
\]

Then `0<=a_w<=1` and `A_w>=|m_0(w)|`.  Since the standard Gaussian
density is bounded by `M=(2pi)^{-1/2}`, the elementary rearrangement
bound around `z_g` gives

\[
 \mathbb E[|Z-z_g|a_w(Z)]\ge {A_w^2\over4M}.         \tag{3.2}
\]

Monotonicity of `T` and antitonicity of `R` give the pointwise identity

\[
 L|z-z_g|=|T(z)-c|+|R(z)-R(z_g)|.                   \tag{3.3}
\]

The first term in (3.3), integrated against `a_w dgamma_1`, is exactly
`m_1(w)`.  Therefore

\[
 {L\over4M}m_0(w)^2
 \le m_1(w)+\mathbb E\{|R(Z)-R(z_g)|a_w(Z)\}.       \tag{3.4}
\]

The last expectation is at most

\[
 K:=\mathbb E|R(Z)-R(z_g)|.
\]

Writing

\[
 A_-=\mathbb E[R\mid Z<z_g]-R(z_g),\qquad
 A_+=R(z_g)-\mathbb E[R\mid Z\ge z_g],
\]

we have `A_-,A_+>=0` and

\[
 K=(1-g)A_-+gA_+,qquad
 D_{\rm map}=g(1-g)(A_-+A_+).                       \tag{3.5}
\]

Thus `K<=delta^{-1}D_map`.  Averaging (3.4) in `w`, using
`L=t^{-1/2}`, proves (3.1).  `square`

The contraction of `T` also implies the pointwise upper bound

\[
 0\le m_1(w)\le\int|y-c|d\rho(y)
 \le L\mathbb E|Z-z_g|\le {C_\delta\over\sqrt t}.   \tag{3.6}
\]

Consequently

\[
 \mathbb E_\nu m_1(W)^2
 \le {C_\delta\over\sqrt t}D_{\rm cut}.             \tag{3.7}
\]

## 4. Transverse quadratic forms

Let `Sigma=Cov(nu)`.  Strong log-concavity and Lichnerowicz give

\[
 C_P(\nu)\le t^{-1},\qquad \Sigma\preceq t^{-1}I.   \tag{4.1}
\]

For every symmetric matrix `A` with `||A||HS=1`, Poincare applied to
`w mapsto w^TAw` yields

\[
 \begin{aligned}
 \operatorname {Var}_\nu(W^TAW)
 &\le {1\over t}\mathbb E|2AW|^2\\
 &={4\over t}\operatorname {tr}(A^2\Sigma)
 \le {4\over t^2}.                                  \tag{4.2}
 \end{aligned}
\]

It follows by Hilbert--Schmidt duality that

\[
 \left\|\mathbb E[m_0(W)(WW^T-\Sigma)]\right\|_{HS}^2
 \le {4\over t^2}\mathbb E m_0(W)^2.               \tag{4.3}
\]

Similarly, by scalar Cauchy--Schwarz and `Sigma preceq t^{-1}I`,

\[
 \left|\mathbb E[Wm_1(W)]\right|^2
 \le {1\over t}\operatorname {Var}_\nu(m_1(W)).    \tag{4.4}
\]

Here (4.4) means: take the supremum over unit vectors `a` in
`|E<a,W>(m_1-Em_1)|` and use
`Var<a,W><=1/t`.

## 5. Assembly

The threshold halfspace has zero transverse row in its correlated second
moment, by independence and centering.  Hence the transverse-transverse
block of `PD` is, up to sign,

\[
 B=\mathbb E_\nu[m_0(W)WW^T]
  =\mathbb E_\nu[m_0(W)(WW^T-\Sigma)],              \tag{5.1}
\]

where (2.3) was used.  Its transverse-active column is, again up to sign,

\[
 b=\mathbb E[W\{cm_0(W)+m_1(W)\}]
  =\mathbb E[Wm_1(W)].                               \tag{5.2}
\]

The two blocks are orthogonal in Hilbert--Schmidt norm.  Equations
(3.1), (3.7), (4.3), and (4.4) give

\[
 \begin{aligned}
 \|PD\|_{HS}^2
 &=\|B\|_{HS}^2+|b|^2\\
 &\le {4\over t^2}C_\delta\sqrt t\,\Delta
     +{1\over t}{C_\delta\over\sqrt t}D_{\rm cut}\\
 &\le C_\delta t^{-3/2}\Delta.
 \end{aligned}
\]

This proves Theorem 1.

## 6. Boundary of the result

The proof uses product splitting in exactly two places:

1. the threshold `H(Y)` has no transverse correlated first or second
   moment; and
2. the same transverse law `nu` appears in every active slice, allowing
   the quadratic Poincare estimate (4.2) to be applied to `m_0(W)`.

For a general `t`-strongly log-concave posterior, the active marginal
still obeys the exact transport decomposition (2.4), but the conditional
transverse law may vary with `Y`.  A complete angular gluing lemma must
quantify this failure of splitting by the same map defect, or find an
integrated heat-flow identity in which the conditional-law terms cancel.
The theorem above neither asserts nor assumes such a stability result.
