# Affine stability for the strongly log-concave Poincare inequality

## 0. The dimension-free deficit theorem

Let \(\pi\) be a centered \(1\)-strongly log-concave probability on its
intrinsic Euclidean support, let \(Z\sim\pi\), and let
\(f\in W^{1,2}(\pi)\).  Put

\[
 D_\pi(f)=E_\pi|\nabla f|^2-\operatorname {Var}_\pi(f)\ge0.
\]

Then there is a vector \(\ell\) in the intrinsic support such that

\[
 \boxed{
 E_\pi\bigl|f-E_\pi f-\ell\cdot Z\bigr|^2
 \le22D_\pi(f).}                                      \tag{AS}
\]

The constant is universal and independent of the dimension.  Thus equality
in the Poincare-one inequality forces an affine function along a Gaussian
factor, while a small deficit gives an averaged quantitative affine
approximation.  The statement is local to a strongly log-concave law; it
does not assert synchronization of the affine directions for a family of
Gaussian posteriors.

## 1. Proof by Caffarelli contraction and Gaussian chaos

Let \(T=\nabla\varphi\) be the Brenier map taking
\(G\sim N(0,I)\) to \(Z\), and write \(J=DT\).  Caffarelli contraction gives

\[
 0\preceq J\preceq I\quad\text{a.e.}
\]

First take \(f\) smooth and put

\[
 g=f\circ T,\qquad p=(\nabla f)\circ T,qquad
 \ell=E_\gamma\nabla g=E_\gamma[Jp].
\]

The chain rule and \(T_\#\gamma=\pi\) split the deficit into two
nonnegative pieces:

\[
 \begin{aligned}
 D_\pi(f)
 &=D_1+D_2,\\
 D_1&=E_\gamma\bigl(|p|^2-|Jp|^2\bigr),\\
 D_2&=E_\gamma|\nabla g|^2-\operatorname {Var}_\gamma(g).
 \end{aligned}                                      \tag{1.1}
\]

Because \((I-J)^2\preceq I-J^2\),

\[
 E|(I-J)p|^2\le D_1.                                \tag{1.2}
\]

The vector \(\ell\) is the degree-one Gaussian-chaos coefficient of
\(g\).  Gaussian spectral decomposition therefore gives

\[
 \operatorname {Var}_\gamma(g-\ell\cdot G)\le D_2,
 \qquad
 E|\nabla g-\ell|^2\le2D_2.                         \tag{1.3}
\]

Since \(\nabla g=Jp\), equations (1.2)--(1.3) imply

\[
 E|p-\ell|^2
 \le2E|(I-J)p|^2+2E|Jp-\ell|^2
 \le2D_1+4D_2\le4D_\pi(f).                         \tag{1.4}
\]

Using \(\|I-J\|_{op}\le1\) once more,

\[
 \begin{aligned}
 E|(I-J)\ell|^2
 &\le2E|(I-J)(\ell-p)|^2+2E|(I-J)p|^2\\
 &\le 10D_\pi(f).                                  \tag{1.5}
 \end{aligned}
\]

The scalar \(r(G)=\ell\cdot(T(G)-G)\) is centered because both \(T(G)\)
and \(G\) are centered, and \(\nabla r=(J-I)\ell\).  Gaussian Poincare
and (1.5) yield

\[
 E r^2\le10D_\pi(f).                                \tag{1.6}
\]

Finally,

\[
 f(T)-E f-\ell\cdot T
 =\bigl(g-Eg-\ell\cdot G\bigr)+\ell\cdot(G-T).
\]

The squared triangle inequality, (1.3), and (1.6) give

\[
 E|f(T)-Ef-\ell\cdot T|^2
 \le2D_2+20D_\pi(f)\le22D_\pi(f),
\]

which is (AS).

For a general Sobolev function, approximate in \(W^{1,2}(\pi)\), use the
Lipschitz Sobolev chain rule for \(T\), and take a weakly convergent
subsequence of the finite-dimensional vectors \(\ell\).  Extended convex
potentials and hard convex supports are covered directly by the generalized
Caffarelli contraction or by centered strong-convex approximation.  Proper
affine supports are treated intrinsically.  No spectral attainment is used.

## 2. Consequence for a hypothetical small-gap sequence

Let \(X\sim\mu\) be centered and isotropic, let \(Y=X+G\), and write

\[
 a=\lambda_1(\mu),\qquad
 \lambda=\lambda_1(\mathcal L(Y)).
\]

The posterior reverse comparison and the first-order posterior lift give,
for \(a\le1/2\),

\[
 a-3a^2\le\lambda\le a+2a^2.                       \tag{2.1}
\]

Choose a centered normalized input spectral-window vector \(f\) with
energy \(q\le a+\eta\), and put

\[
 d=E_Y\operatorname {Var}(f(X)\mid Y).
\]

The reverse-transfer test gives

\[
 d\le q,
 \qquad
 \lambda\le\frac d{1-d},
 \qquad
 d\ge\frac\lambda{1+\lambda}.                      \tag{2.2}
\]

Hence, as \(a\downarrow0\) and \(\eta=o(a^2)\),

\[
 q-d=O(a^2).                                        \tag{2.3}
\]

For each posterior

\[
 \pi_y(dx)\propto e^{-V(x)-|x-y|^2/2}\,dx,
\]

apply (AS) intrinsically to \(f\).  Since

\[
 E_YD_{\pi_Y}(f)
 =q-d,
\]

one obtains the dimension-free averaged posterior-affine estimate

\[
 \boxed{
 E_Y\inf_{c\in\mathbb R,\,\ell\in\mathbb R^n}
 E_{\pi_Y}|f(X)-c-\ell\cdot(X-E_{\pi_Y}X)|^2
 \le22(q-d)=O(a^2).}                                \tag{2.4}
\]

Thus a hypothetical counterexample sequence has input bottom modes which
are affine on almost every strongly log-concave posterior up to total
squared error \(O(a^2)\), while their total posterior conditional variance
is \(a+O(a^2)\).  The relative error is \(O(a)\).

This does not yet prove KLS.  The affine coefficient in (2.4) depends on
the posterior state \(y\); (2.4) permits switching on fibers carrying little
signal.  The audited Hilbert--Schmidt posterior theorem controls the local
variation of a nearly saturated direction without a trace loss, but a
global synchronization/amplitude-connectivity theorem is still required.
