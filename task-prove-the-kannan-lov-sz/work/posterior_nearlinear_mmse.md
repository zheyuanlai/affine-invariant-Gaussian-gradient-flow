# Gaussian posterior identities for a near-linear low mode

## 1. Setting

Let \(\mu\) be centered, isotropic, and log-concave on \(\mathbb R^n\).
Let \(X\sim\mu\), \(G\sim N(0,I_n)\) be independent, and put

\[
                         Y=X+G.
\]

For every observation \(y\), denote the posterior law of \(X\) by
\(\mu_y\), its mean by \(m(y)\), and its covariance by

\[
                         A(y)=\operatorname {Cov}(X\mid Y=y).
\]

The posterior potential is the sum of the original convex potential and
\(|x-y|^2/2\).  Thus \(\mu_y\) is 1-strongly log-concave on the same
convex support and

\[
                         0\preceq A(y)\preceq I.       \tag{1.1}
\]

All statements below first follow for a smooth positive density and then
for an arbitrary log-concave probability by intrinsic Gaussian
regularization on its affine support.  The conditional-expectation
identities themselves require only finite second moments.

Let \(f\in W^{1,2}(\mu)\) satisfy

\[
 \mathbb Ef=0,\qquad \mathbb Ef^2=1,
 \qquad \mathbb E|\nabla f|^2=\lambda.
\]

Set

\[
 a=\mathbb E[Xf],\qquad f=a\cdot x+r,
 \qquad \delta=\|r\|_2.
\]

If \(r\) is the orthogonal residual from the affine projection, then

\[
 \mathbb Er=0,\qquad \mathbb E[Xr]=0,
 \qquad |a|^2=1-\delta^2.                             \tag{1.2}
\]

## 2. Posterior gradient and the exact affine-residual split

Define the posterior transform

\[
 F(y)=\mathbb E[f(X)\mid Y=y].
\]

Differentiating the exponential-family posterior gives

\[
 \nabla F(y)=\operatorname {Cov}(X,f(X)\mid Y=y)=:c(y). \tag{2.1}
\]

The affine decomposition of \(f\) yields the pointwise identity

\[
 \boxed{\quad c(y)=A(y)a+e(y),\qquad
 e(y)=\operatorname {Cov}(X,r(X)\mid Y=y).\quad}       \tag{2.2}
\]

The posterior Poincare inequality and the law of total variance give two
different useful estimates.  First,

\[
 |c(y)|^2
 \le \operatorname {Var}(f\mid Y=y)
 \le \mathbb E[|\nabla f|^2\mid Y=y],
\]

and hence

\[
                         \mathbb E|c(Y)|^2\le\lambda. \tag{2.3}
\]

Second, Cauchy--Schwarz and (1.1) give

\[
 |e(y)|^2
 \le \operatorname {Var}(r\mid Y=y),
\]

so that

\[
                         \mathbb E|e(Y)|^2\le\delta^2.\tag{2.4}
\]

Combining (2.2)--(2.4),

\[
 \boxed{\quad
 \mathbb E|A(Y)a|^2\le2(\lambda+\delta^2).
 \quad}                                               \tag{2.5}
\]

No spectral-gap assertion is used in (2.5); only posterior strong
log-concavity is used.

## 3. A directional MMSE floor would close the near-linear branch

For a unit vector \(b\), define the directional unit-noise MMSE

\[
 \operatorname {mmse}_\mu(b)
 =\mathbb E\operatorname {Var}(b\cdot X\mid X+G)
 =\mathbb E[b^TA(Y)b].                                \tag{3.1}
\]

Because \(0\preceq A\preceq I\),

\[
 b^TAb\le |Ab|,
 \qquad
 \operatorname {mmse}_\mu(b)
 \le \big(\mathbb E|A(Y)b|^2\big)^{1/2}.             \tag{3.2}
\]

Consequently, the putative dimension-free statement

\[
 \boxed{\qquad
 \inf_{\mu,b}\operatorname {mmse}_\mu(b)\ge c_0>0,
 \qquad}                                              \tag{DMMSE}
\]

where the infimum runs over all centered isotropic log-concave laws and
all unit vectors, would imply from (2.5), with \(b=a/|a|\), that

\[
 c_0^2(1-\delta^2)
 \le2(\lambda+\delta^2).                              \tag{3.3}
\]

In fact the Hilbert-space triangle inequality gives the sharper version

\[
 \sqrt\lambda\ge\|c(Y)\|_{L^2}
 \ge |a|\,\|A(Y)b\|_{L^2}-\|e(Y)\|_{L^2}
 \ge c_0\sqrt{1-\delta^2}-\delta.                     \tag{3.3a}
\]

Thus, for \(\delta\le c_0/4\), one obtains a numerical lower bound on
\(\lambda\).  This is a clean replacement for the terminal
parallel-coupling alignment lemma in the very-near-linear regime.

The statement (DMMSE) must not be treated as elementary.  Let
\(q\) be the density of \(Y\).  Tweedie's identity gives

\[
 m(y)=y+\nabla\log q(y),\qquad \nabla m(y)=A(y).
\]

For \(h_b(y)=b\cdot m(y)\),

\[
 \operatorname {Var}_q(h_b)=1-\operatorname {mmse}_\mu(b),
 \qquad
 \int|\nabla h_b|^2dq=\mathbb E|A(Y)b|^2.             \tag{3.4}
\]

Hence failure of (DMMSE) produces a restricted, monotone posterior-mean
test function with variance tending to one and Dirichlet energy tending
to zero for the log-concave Gaussian convolution \(q\), whose covariance
is \(2I\).  Proving (DMMSE) therefore requires a genuine nonlinear
correlation or posterior-covariance input; the scalar MMSE bound for the
one-dimensional marginal \(b\cdot X\) is insufficient because the
orthogonal noisy observations may carry nonlinear information about
\(b\cdot X\).

## 4. Covariance and flow martingales

In standard stochastic localization, let \(A_t\) be the posterior
covariance and let \(M_t\) solve

\[
 M_0=I,\qquad \dot M_t=A_tM_t.
\]

The covariance SDE has the form

\[
 dA_t=H_t\,dB_t-A_t^2dt,
\]

where \(H_t\,dB_t\) denotes the centered third-moment matrix noise.
It follows exactly that

\[
 \boxed{\quad J_t=A_tM_t\text{ is a matrix martingale},
 \qquad J_0=I.\quad}                                  \tag{4.1}
\]

Indeed, the two drift terms in \(d(A_tM_t)\) cancel.  Its transpose
\(M_t^TA_t\) is the corresponding transpose martingale.  For the full
posterior covariance \(c_t=\operatorname {Cov}_t(X,f)\), the same
calculation gives

\[
                         M_t^Tc_t\text{ is a martingale}. \tag{4.2}
\]

The affine and residual parts in (2.2) are separately compatible with
(4.2):

\[
 \mathbb E[M_1^TA_1a]=a,
 \qquad
 \mathbb E[M_1^Te_1]=0.                              \tag{4.3}
\]

These identities explain the cancellation a hypothetical near-linear
small mode would require.  They do not by themselves upper-bound the
directional amplification: a matrix martingale with mean \(I\) can have a
large one-dimensional second moment.  Any use of (4.1) therefore still
needs either (DMMSE), a terminal normal-alignment estimate, or a new
matrix-martingale inequality that uses the posterior-gradient structure.

## 5. Scope

The deductions above close no general KLS branch without an additional
theorem.  Their value is to separate three distinct statements:

1. the proved posterior split and energy estimate (2.5);
2. the precise directional posterior-covariance target (DMMSE); and
3. the stronger pathwise terminal-alignment target from parallel
   coupling.

Trace MMSE bounds, entropy bounds of order \(n\), and scalar observation
bounds do not imply (DMMSE), because they do not control the least
eigenvalue of the average posterior covariance.
