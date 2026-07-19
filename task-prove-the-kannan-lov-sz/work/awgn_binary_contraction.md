# A fixed AWGN channel and balanced Boolean contraction

## Verdict

Let \(X\sim\mu\) be isotropic and log-concave in \(\mathbb R^n\), let
\(G\sim N(0,I_n)\) be independent, put \(Y=X+G\), and let
\(S=1_E(X)\) for a set \(E\) with \(\mu(E)=1/2\).  Write

\[
  \operatorname{err}_\mu(E)
  :=\mathbb E\min\{\mathbb P(S=1\mid Y),\mathbb P(S=0\mid Y)\}.
\]

I found neither a proof nor a counterexample to

\[
  \inf_{n,\mu,E:\,\mu(E)=1/2}\operatorname{err}_\mu(E)>0. \tag{AWGN}
\]

There is a precise reason.  A universal form of (AWGN) implies KLS by an
exact posterior-perimeter averaging identity, while a dimension-free
Poincare/KLS bound implies (AWGN), with explicit constants.  Thus (AWGN)
is quantitatively equivalent to KLS, rather than an easier consequence of
posterior strong log-concavity.

There is, however, a genuinely stronger pointwise fact which seems useful:
if

\[
  p(y):=\mathbb P(X\in E\mid Y=y),
\]

then

\[
  \boxed{\quad \|\nabla \Phi^{-1}(p(y))\|_2\le 1 \quad} \tag{P}
\]

at every differentiability point.  Here \(\Phi\) is the standard normal
cdf.  This is sharp in the class of 1-strongly log-concave posteriors and is
strictly stronger than the elementary estimate
\(\|\nabla p\|^2\le p(1-p)\).  It says that high-confidence posterior phases
must be macroscopically separated.  The exact remaining structural target
is therefore a **transition-mass lemma** for posterior-derived probit
functions, not a generic restatement of noise contraction:

> Prove that every posterior probit \(F=\Phi^{-1}(p)\) arising above
> satisfies \((\mu*\gamma)(|F|\le1)\ge c\) for a universal \(c>0\).

This would imply \(\operatorname{err}_\mu(E)\ge c\Phi(-1)\).  It is strictly
stronger than the desired averaged error bound: a positive average can be
carried entirely by points with \(|F|>1\), whereas this target forces mass
at the unit noise scale.  The special posterior origin of \(F\), beyond
merely being 1-Lipschitz, has to be used; proving such a statement for
arbitrary Lipschitz functions under arbitrary isotropic log-concave laws
runs directly into KLS-type concentration.

All explicit tests below are consistent with a universal constant.  The
smallest exact canonical value found is
\(0.2205682480\ldots\), for a coordinate sign of the isotropic cube.  A
balanced OR of rare exponential-tail events, which is the most plausible
coding-style obstruction, has Bayes error at least
\(0.1514348782\ldots+o(1)\).

---

## 1. Exact channel identities

Let \(\gamma(x)=(2\pi)^{-n/2}e^{-\|x\|^2/2}\), let
\(\nu=\mu*\gamma\) be the law of \(Y\), and let

\[
 q_y(dx)=\mathcal L(X\mid Y=y)
 =\frac{\gamma(y-x)\,\mu(dx)}{\nu(y)}.
\]

If \(\mu(dx)=e^{-V(x)}dx\), with \(V\) convex in the extended-valued
sense, then

\[
 q_y(dx)\propto
 \exp\left[-V(x)-\frac12\|x-y\|^2\right]dx. \tag{1.1}
\]

Thus every \(q_y\) is 1-strongly log-concave.

Set \(p(y)=q_y(E)\), \(Z=2S-1\), and
\(m(y)=\mathbb E[Z\mid Y=y]=2p(y)-1\).  If

\[
 a(y)=\int_E\gamma(y-x)\,\mu(dx),\qquad
 b(y)=\int_{E^c}\gamma(y-x)\,\mu(dx),
\]

then \(a+b=\nu\), and equal-prior binary testing gives

\[
\begin{aligned}
 \operatorname{err}_\mu(E)
 &=\int \min\{a(y),b(y)\}\,dy\\
 &=\mathbb E_\nu\min\{p(Y),1-p(Y)\}\\
 &=\frac12\left(1-
 \|\mathcal L(Y\mid S=1)-\mathcal L(Y\mid S=0)\|_{\rm TV}\right).
                                                        \tag{1.2}
\end{aligned}
\]

Let

\[
 B(E):=\mathbb E_\nu p(Y)(1-p(Y)).
\]

Since, for \(0\le p\le1\),

\[
 p(1-p)\le \min(p,1-p)\le 2p(1-p),
\]

we have the exact comparison

\[
 B(E)\le \operatorname{err}_\mu(E)\le 2B(E). \tag{1.3}
\]

If, conditionally on \(Y\), \(X'\) is an independent draw from \(q_Y\),
then the reversible two-step kernel \(X\to Y\to X'\) satisfies

\[
 D(E):=\mathbb P(1_E(X)\ne1_E(X'))=2B(E), \tag{1.4}
\]

and hence

\[
 \frac12D(E)\le\operatorname{err}_\mu(E)\le D(E). \tag{1.5}
\]

Equivalently,

\[
 \operatorname{mmse}(Z\mid Y)=1-\mathbb E m(Y)^2=4B(E), \tag{1.6}
\]

so

\[
 2\operatorname{err}_\mu(E)
 \le \operatorname{mmse}(Z\mid Y)
 \le 4\operatorname{err}_\mu(E).                  \tag{1.7}
\]

These formulas are the exact one-step noise-sensitivity, total-variation,
and Boolean-MMSE versions of the question.

---

## 2. Perimeter averaging: (AWGN) implies KLS

Let \(I(t)=\phi(\Phi^{-1}(t))\) be the Gaussian isoperimetric profile.
The Gaussian isoperimetric comparison for a 1-strongly log-concave law
(equivalently, Caffarelli contraction followed by Gaussian isoperimetry)
gives, for every posterior,

\[
 q_y^+(E)\ge I(p(y)). \tag{2.1}
\]

Moreover,

\[
 I(t)\ge \sqrt{\frac2\pi}\min(t,1-t),\qquad 0\le t\le1, \tag{2.2}
\]

with equality in (2.2) at \(t=1/2\).

For a finite-perimeter set, multiplying its posterior boundary density by
\(\nu(y)\) and integrating in \(y\) gives

\[
\begin{aligned}
 \mathbb E_\nu q_Y^+(E)
 &=\int_{\mathbb R^n}\int_{\partial^*E}
      \gamma(y-x)e^{-V(x)}\,d\mathcal H_{n-1}(x)\,dy\\
 &=\mu^+(E).                                      \tag{2.3}
\end{aligned}
\]

Approximation gives the same statement for the Minkowski perimeter.  From
(2.1)--(2.3),

\[
 \boxed{\quad
 \mu^+(E)\ge \sqrt{\frac2\pi}\,
 \operatorname{err}_\mu(E).
 \quad}                                           \tag{2.4}
\]

Consequently, a universal lower bound
\(\operatorname{err}_\mu(E)\ge\varepsilon_0\) for every balanced \(E\)
would give the balanced KLS expansion bound

\[
 \mu^+(E)\ge \sqrt{\frac2\pi}\,\varepsilon_0.   \tag{2.5}
\]

Nothing in posterior strong log-concavity supplies the missing lower bound
on the left side of (1.2); it only supplies the upper comparison of error
by the original perimeter.

---

## 3. The converse: KLS/Poincare implies AWGN contraction

Let \(C_P(\mu)\) be the optimal Poincare constant.  Convolution and variance
decomposition give

\[
 C_P(\nu)\le C_P(\mu)+1.                          \tag{3.1}
\]

Indeed, condition first on \(G\), use the Poincare inequality of \(\mu\),
and then use the Gaussian Poincare inequality; Jensen bounds the gradient
of the conditional expectation.

Differentiating the exponential-family posterior gives

\[
 \nabla p(y)=\operatorname{Cov}_{q_y}(1_E(X),X).  \tag{3.2}
\]

Brascamp--Lieb gives \(\operatorname{Cov}_{q_y}(X)\preceq I\), and
Cauchy--Schwarz therefore gives

\[
 \|\nabla p(y)\|^2\le p(y)(1-p(y)).               \tag{3.3}
\]

Since \(\mathbb E_\nu p=1/2\), the total-variance identity and (3.1)--(3.3)
give

\[
\begin{aligned}
 \frac14
 &=\operatorname{Var}_\mu(1_E)\\
 &=\operatorname{Var}_\nu(p)+B(E)\\
 &\le (C_P(\mu)+1)B(E)+B(E).
\end{aligned}
\]

Thus

\[
 \boxed{\quad
 \operatorname{err}_\mu(E)\ge B(E)
 \ge \frac{1}{4(C_P(\mu)+2)}.
 \quad}                                           \tag{3.4}
\]

For example, a Cheeger bound \(h(\mu)\ge h_0\), together with
\(C_P(\mu)\le4/h_0^2\), gives

\[
 \operatorname{err}_\mu(E)
 \ge \frac{h_0^2}{16+8h_0^2}.                    \tag{3.5}
\]

Equations (2.4) and (3.4) show the quantitative KLS equivalence.

### Full maximal correlation

Let \(Tf(y)=\mathbb E[f(X)\mid Y=y]\).  For arbitrary square-integrable
\(f\), the same calculation gives

\[
 \|\nabla Tf(y)\|^2\le\operatorname{Var}_{q_y}(f). \tag{3.6}
\]

Writing

\[
 A_f=\operatorname{Var}_\nu(Tf),\qquad
 B_f=\mathbb E_\nu\operatorname{Var}_{q_Y}(f),
\]

we have \(\operatorname{Var}_\mu f=A_f+B_f\) and

\[
 A_f\le(C_P(\mu)+1)B_f.                           \tag{3.7}
\]

Therefore the Hirschfeld--Gebelein--Renyi maximal correlation obeys

\[
 \rho_{\rm HGR}(X,Y)^2
 =\sup_f\frac{A_f}{A_f+B_f}
 \le\frac{C_P(\mu)+1}{C_P(\mu)+2}.               \tag{3.8}
\]

Conversely, if \(\rho_{\rm HGR}^2\le1-\delta\), then
\(B_f\ge\delta\operatorname{Var}_\mu f\).  Posterior Poincare gives

\[
 B_f\le \int\|\nabla f\|^2d\mu,
\]

so \(C_P(\mu)\le1/\delta\).  Thus a universal full maximal-correlation
gap is itself quantitatively equivalent to KLS.  It cannot be inserted as
an independent strong-data-processing premise.

For a balanced Boolean \(Z\),

\[
 \rho_Z^2:=\mathbb E(\mathbb E[Z\mid Y]^2),
 \qquad 1-\rho_Z^2=4B(E),                         \tag{3.9}
\]

and (1.3) becomes

\[
 \frac{1-\rho_Z^2}{4}
 \le\operatorname{err}_\mu(E)
 \le\frac{1-\rho_Z^2}{2}.                       \tag{3.10}
\]

This is the exact sense in which the desired result is a maximal-correlation
gap restricted to balanced Boolean functions.

### Entropy

Let \(h_2\) denote binary entropy in nats.  Then

\[
 H(S\mid Y)=\mathbb E h_2(p(Y)).                  \tag{3.11}
\]

For \(0\le t\le1/2\), concavity gives
\(h_2(t)\ge2(\log2)t\), while binary Fano gives

\[
 2(\log2)\operatorname{err}_\mu(E)
 \le H(S\mid Y)
 \le h_2(\operatorname{err}_\mu(E)).             \tag{3.12}
\]

Hence a balanced-binary entropy SDPI is equivalent, up to explicit scalar
functions, to (AWGN), and therefore also has KLS content.  A full KL-SDPI
for arbitrary likelihood ratios would be closer to a log-Sobolev statement
and is stronger than what KLS provides.  The vector Gaussian I-MMSE formula
for \(I(S;\sqrt tX+G)\) involves the difference between the estimators of
\(X\) with and without \(S\); it does not directly yield a lower bound on
the Boolean MMSE in (1.6).

---

## 4. A sharp structural lemma: posterior probit is 1-Lipschitz

The following conditional-centroid inequality is the main positive lemma.

### Lemma 4.1 (strongly log-concave conditional centroid)

Let \(Q\) be 1-strongly log-concave on \(\mathbb R^n\), let \(A\) be
measurable with \(Q(A)=t\), and let \(I(t)=\phi(\Phi^{-1}(t))\).  Then

\[
 \left\|\operatorname{Cov}_Q(1_A(X),X)\right\|_2\le I(t). \tag{4.1}
\]

#### Proof

Fix a unit vector \(u\) and put \(W=u\cdot X\).  The one-dimensional
marginal of a 1-strongly log-concave law is again 1-strongly log-concave:
after extracting the factor \(e^{-w^2/2}\), this is exactly Prekopa's
theorem.  Its increasing transport \(T\) from a standard normal variable
\(N\) is 1-Lipschitz (the one-dimensional Caffarelli contraction theorem).

Let \(r(w)=Q(A\mid W=w)\).  Among all \(0\le r\le1\) with
\(\mathbb E r(W)=t\), the bathtub principle says that
\(\operatorname{Cov}(r(W),W)\) is maximized by the indicator of the upper
\(t\)-tail of \(W\).  Write that tail as \(N\ge a\), where
\(\mathbb P(N\ge a)=t\).  If \(N_+\) and \(N_-\) are independent normals
conditioned on \(N\ge a\) and \(N<a\), respectively, then

\[
 T(N_+)-T(N_-)\le N_+-N_-
\]

pointwise, because \(T\) is increasing and 1-Lipschitz.  Consequently the
upper-tail covariance for \(W=T(N)\) is no larger than the corresponding
Gaussian covariance, which is

\[
 \operatorname{Cov}(1_{\{N\ge a\}},N)=\phi(a)=I(t).
\]

Applying the same argument to the lower tail controls the absolute value.
Taking the supremum over \(u\) proves (4.1). \(\square\)

Applying Lemma 4.1 to (3.2) yields

\[
 \|\nabla p(y)\|\le I(p(y)).                      \tag{4.2}
\]

Because \((\Phi^{-1})'(t)=1/I(t)\), (4.2) proves (P):

\[
 F(y):=\Phi^{-1}(p(y))\quad\hbox{is 1-Lipschitz}. \tag{4.3}
\]

The weaker Brascamp--Lieb estimate (3.3) only says that
\(\arcsin(2p-1)\) is 1-Lipschitz.  Thus (4.3) contains substantially more
tail information.

### A quantitative two-phase certificate

Since

\[
 \operatorname{err}_\mu(E)=\mathbb E_\nu\Phi(-|F(Y)|), \tag{4.4}
\]

small error forces two separated posterior phases.  Precisely, suppose
\(\operatorname{err}_\mu(E)=\varepsilon\), fix \(a>0\), and put
\(q=\Phi(-a)\).  Define

\[
 A_+=\{F\ge a\},\quad A_-=\{F\le-a\},\quad R=\{|F|<a\}.
\]

Then

\[
 \nu(R)\le\frac{\varepsilon}{q},                 \tag{4.5}
\]

and the balance identity \(\mathbb E_\nu p=1/2\) gives

\[
 \nu(A_\pm)\ge
 \frac{1/2-q}{1-q}-\frac{\varepsilon}{q}.        \tag{4.6}
\]

Finally, (4.3) gives

\[
 \operatorname{dist}(A_+,A_-)\ge2a.              \tag{4.7}
\]

For example, taking \(q=\sqrt\varepsilon\) shows that error tending to zero
would create two sets whose masses tend to \(1/2\), separated by

\[
 2\Phi^{-1}(1-\sqrt\varepsilon)
 \asymp 2\sqrt{\log(1/\varepsilon)}.              \tag{4.8}
\]

This is a useful inverse formulation: a counterexample to (AWGN) must
produce a Gaussian-smoothed isotropic log-concave law with two large,
widely separated confidence phases generated by one common posterior set.
Generic second-moment control only rules out separation of order
\(\sqrt n\); excluding (4.8) at a universal scale is exactly where KLS-type
concentration enters.

---

## 5. Why Gaussian rearrangement and joint log-concavity stop short

For \(\mu=\gamma_n\), \(Y/\sqrt2\) is standard Gaussian and has correlation
\(1/\sqrt2\) with \(X\).  Borell's two-set Gaussian rearrangement theorem
then solves the balanced problem exactly: for any balanced \(E\) and any
decision set \(B\) in observation space,

\[
 \mathbb P(X\in E,Y/\sqrt2\in B)
 +\mathbb P(X\notin E,Y/\sqrt2\notin B)\le\frac34. \tag{5.1}
\]

Thus every balanced Gaussian label has Bayes error at least \(1/4\), with
equality for halfspaces.

For a general log-concave prior, each fixed posterior admits a Caffarelli
contraction from a Gaussian.  Lemma 4.1 is the part of Gaussian
rearrangement which survives consistently: it controls every posterior
centroid and hence the derivative of \(p\).  What does **not** survive is a
single rearrangement of \(E\) valid for all \(y\).  The transport map to
\(q_y\) depends on \(y\); rearranging \(E\) separately in each posterior
destroys both its common origin in \(x\)-space and the global balance under
\(\mu\).  This is why posterior-by-posterior Gaussian isoperimetry gives
(2.4) but not a lower bound on (1.2).

The joint density

\[
 (x,y)\longmapsto e^{-V(x)-\|y-x\|^2/2}           \tag{5.2}
\]

is log-concave.  It is not uniformly strongly log-concave: when \(V\) is
locally affine, the Hessian of \(\|y-x\|^2/2\) has zero directions
\((u,u)\).  After an affine normalization, applying a dimension-free
isoperimetric theorem to this joint log-concave law would again be a KLS
input.  There is also no useful closure statement for the two testing
subdensities: \(1_Ee^{-V}\) is log-concave only when \(E\) is convex, and
\(1_{E^c}e^{-V}\) generally is not; for arbitrary Boolean labels neither
output subdensity is log-concave.  Joint log-concavity alone therefore does
not control their overlap.

---

## 6. Explicit stress tests

All numerical constants below use noise variance one.

### 6.1 Gaussian halfspaces and radial sets

For \(X\sim N(0,I_n)\) and \(E=\{u\cdot X\ge0\}\), the Bayes rule is
\(u\cdot Y\ge0\).  The correlation of \(u\cdot X\) and
\(u\cdot Y/\sqrt2\) is \(1/\sqrt2\), so

\[
 \operatorname{err}(E)
 =\frac1\pi\arccos\frac1{\sqrt2}=\frac14.        \tag{6.1}
\]

By (5.1), this is the minimum over all balanced Gaussian sets.

Now let \(E_n=\{\|X\|^2\ge m_n\}\), where \(m_n\) is the median of
\(\chi_n^2\).  The posterior is radial and the Bayes decision is a radial
threshold in \(Y\).  The joint central limit theorem for

\[
 \frac{\|X\|^2-n}{\sqrt{2n}},\qquad
 \frac{\|Y/\sqrt2\|^2-n}{\sqrt{2n}}
\]

has correlation \((1/\sqrt2)^2=1/2\).  Hence

\[
 \operatorname{err}(E_n)\longrightarrow
 \frac1\pi\arccos\frac12=\frac13.               \tag{6.2}
\]

Radial coding is therefore less stable than a halfspace in the Gaussian
model.

### 6.2 The isotropic cube: coordinate, parity, and majority

Let \(X_i\) be iid uniform on \([-a,a]\), with \(a=\sqrt3\).  For the
coordinate sign \(E=\{X_1\ge0\}\), monotone likelihood ratio gives the
Bayes rule \(Y_1\ge0\), and

\[
\begin{aligned}
 e_{\square}
 &=\frac1a\int_0^a\Phi(-x)\,dx\\
 &=\Phi(-a)+\frac{\phi(0)-\phi(a)}a\\
 &=0.22056824804474245\ldots .                    \tag{6.3}
\end{aligned}
\]

Let \(Z_i=\operatorname{sign}(X_i)\) and
\(m_i(Y_i)=\mathbb E[Z_i\mid Y_i]\).  For parity
\(Z=\prod_{i=1}^nZ_i\), independence gives

\[
 \operatorname{err}_{\rm parity}(n)
 =\frac12\left[1-(\mathbb E|m_1(Y_1)|)^n\right]
 =\frac12\left[1-(1-2e_{\square})^n\right]
 \longrightarrow\frac12.                        \tag{6.4}
\]

Thus parity is destroyed, not protected, by redundant dimension.

For odd \(n\), let \(Z=\operatorname{sign}(\sum_i Z_i)\).  Put

\[
 r_{\square}^2=\mathbb E[m_1(Y_1)^2].             \tag{6.5}
\]

If

\[
 a_+(y)=\frac{\Phi(y)-\Phi(y-a)}{2a},\qquad
 a_-(y)=\frac{\Phi(y+a)-\Phi(y)}{2a},
\]

then

\[
 r_{\square}^2
 =\int_{\mathbb R}\frac{(a_+(y)-a_-(y))^2}
                         {a_+(y)+a_-(y)}\,dy
 =0.3968031373376871\ldots .                     \tag{6.6}
\]

A joint and conditional CLT shows that the asymptotic Bayes rule is the
sign of \(\sum_i m_i(Y_i)\).  Therefore

\[
 \operatorname{err}_{\rm majority}(n)
 \longrightarrow
 \frac1\pi\arccos(r_{\square})
 =0.2830863849427284\ldots .                     \tag{6.7}
\]

There is also an all-label check.  The sharp Poincare constant of
\(\mathrm{Unif}[-\sqrt3,\sqrt3]\) is \(12/\pi^2\), and it tensorizes.
Thus (3.4) proves, for **every** balanced measurable set in every product
cube,

\[
 \operatorname{err}(E)
 \ge\frac{1}{4(2+12/\pi^2)}
 =0.07773984271791194\ldots .                    \tag{6.8}
\]

So no continuous checkerboard or error-correcting-code partition of the
cube can make the error vanish.

### 6.3 A simplex barycentric cut

Let \(U=(U_1,\ldots,U_{n+1})\) be uniform on the standard \(n\)-simplex,
so \(U\sim\mathrm{Dirichlet}(1,\ldots,1)\).  In its affine span the
isotropic version is

\[
 X=\sqrt{(n+1)(n+2)}\left(U-\frac1{n+1}{\bf1}\right). \tag{6.9}
\]

Let \(v\) be the unit direction toward the first vertex and set

\[
 T_n=v\cdot X
 =(n+1)\sqrt{\frac{n+2}{n}}
   \left(U_1-\frac1{n+1}\right).                 \tag{6.10}
\]

Then \(T_n\Rightarrow U-1\), where \(U\sim\mathrm{Exp}(1)\).  The balanced
barycentric cut is \(E_n=\{U_1\ge1-2^{-1/n}\}\), whose threshold in the
\(T_n\)-coordinate tends to \(\log2-1\).

One must check that the other \(n-1\) observed coordinates do not secretly
decode this cut.  The Dirichlet factorization gives the exact decomposition

\[
 X=T_nv+a_n(T_n)Z_n,\qquad
 a_n(t)=\frac{n}{n+1}\sqrt{\frac{n+2}{n}}-\frac{t}{n+1}, \tag{6.11}
\]

where \(Z_n\), independent of \(T_n\), is isotropic in \(v^\perp\) and has
dimension \(n-1\).  Let \(K_t\) be the law of
\(a_n(t)Z_n+G_\perp\).  Coupling with the same \(Z_n\), and using convexity
of relative entropy under Gaussian convolution, gives

\[
 D_{\rm KL}(K_t\|K_s)
 \le\frac12\mathbb E\|(a_n(t)-a_n(s))Z_n\|^2
 =\frac{(n-1)(t-s)^2}{2(n+1)^2}.                 \tag{6.12}
\]

Conditioning also on \(T_n+G_v\), the pairwise-KL bound yields

\[
 I(T_n;Y_\perp\mid T_n+G_v)
 \le\frac{n-1}{(n+1)^2}.                         \tag{6.13}
\]

Pinsker's inequality shows that the improvement in binary Bayes risk from
seeing \(Y_\perp\), in addition to the scalar observation \(T_n+G_v\), is
at most

\[
 \sqrt{\frac{n-1}{2(n+1)^2}}=O(n^{-1/2}).        \tag{6.14}
\]

Thus the full Bayes error converges to the scalar shifted-exponential
error calculated below:

\[
 \operatorname{err}_{\rm simplex}(E_n)
 \longrightarrow 0.2861485156834012\ldots .      \tag{6.15}
\]

Corners of the simplex therefore do not create an asymptotically noiseless
balanced barycentric bit.

### 6.4 The isotropic \(\ell_1\) ball

For the uniform law on an isotropically scaled \(\ell_1^n\) ball, the
coordinate signs are independent fair signs, independent of the vector of
absolute values.  For \(E=\{X_1\ge0\}\), the likelihood ratio has the sign
of \(Y_1\), even after conditioning on all other observed coordinates:
for every possible magnitude \(r>0\),

\[
 \frac{\phi(y_1-r)}{\phi(y_1+r)}=e^{2r y_1}.
\]

A positive mixture preserves this sign.  Hence the Bayes rule is exactly
\(Y_1\ge0\) and

\[
 \operatorname{err}(E)=\mathbb E\Phi(-|X_1|).    \tag{6.16}
\]

The isotropic coordinate converges to the variance-one Laplace law, for
which \(|X_1|\sim\mathrm{Exp}(\sqrt2)\).  Therefore

\[
\begin{aligned}
 \operatorname{err}_{\ell_1}(E)
 &\longrightarrow
 \int_0^\infty\sqrt2e^{-\sqrt2x}\Phi(-x)\,dx\\
 &=\frac12-e\,\Phi(-\sqrt2)\\
 &=0.28620821192209656\ldots .                   \tag{6.17}
\end{aligned}
\]

The dependence among coordinate magnitudes does not yield extra sign
information.

### 6.5 Product exponentials, including a balanced rare-event OR

Let \(U\sim\mathrm{Exp}(1)\), so \(X=U-1\) is isotropic.  The balanced
one-coordinate label is \(S=1_{\{U\ge\log2\}}\).  For the translated
observation \(W=U+G\), the two joint subdensities are

\[
\begin{aligned}
 A(w)&=e^{-w+1/2}\Phi(w-1-\log2),\\
 B(w)&=e^{-w+1/2}
       [\Phi(1+\log2-w)-\Phi(1-w)].              \tag{6.18}
\end{aligned}
\]

The likelihood ratio is monotone.  The decision threshold is
\(w_0=1+z_0\), where

\[
 2\Phi(z_0-\log2)=\Phi(z_0),\qquad
 z_0=0.04924142985630897\ldots .                 \tag{6.19}
\]

Direct integration gives

\[
 e_{\exp}
 =\int_0^{\log2}e^{-u}\Phi(u-w_0)\,du
  +\int_{\log2}^\infty e^{-u}\Phi(w_0-u)\,du
 =0.2861485156834012\ldots .                    \tag{6.20}
\]

The sharp Poincare constant of a rate-one exponential is \(4\), and it
tensorizes.  Consequently, for **every** balanced label under every product
of shifted rate-one exponentials, (3.4) gives

\[
 \operatorname{err}(E)\ge\frac1{24}.             \tag{6.21}
\]

This already rules out parity, majority, and coding counterexamples in the
whole product family.  It is still instructive to compute the most dangerous
rare-event construction exactly.

Let \(p_n=1-2^{-1/n}\), let \(t_n=-\log p_n\), and define

\[
 S_n=1_{\{\max_{1\le i\le n}U_i\ge t_n\}}.
\]

Then \(\mathbb P(S_n=1)=1/2\).  For one coordinate let
\(q_n(W)=\mathbb P(U\ge t_n\mid W)\) and

\[
 b_n=\mathbb E[q_n(W)(1-q_n(W))].
\]

After writing \(W=t_n+s\) and completing the square, the two local
subdensities, divided by \(p_n\), converge to

\[
 A_0(s)=e^{-s+1/2}\Phi(s-1),\qquad
 B_0(s)=e^{-s+1/2}\Phi(1-s).                     \tag{6.22}
\]

Since \(A_0+B_0=e^{-s+1/2}\), dominated convergence gives the exact limit

\[
 \frac{b_n}{p_n}\longrightarrow
 \kappa:=\int_{\mathbb R}e^{-s+1/2}
             \Phi(s-1)\Phi(1-s)\,ds
 =2\Phi(1/\sqrt2)-1
 =0.5204998778130465\ldots .                    \tag{6.23}
\]

The corresponding one-coordinate rare-event Bayes error divided by its
prior mass tends to

\[
 \int e^{-s+1/2}\min\{\Phi(s-1),\Phi(1-s)\}\,ds
 =2\Phi(1)-1
 =0.6826894921370859\ldots .                    \tag{6.24}
\]

Thus even normalized rare-tail maximal correlation stays away from one in
this example.

For the balanced OR, conditional independence gives

\[
 Q_n:=\mathbb P(S_n=0\mid W_1,\ldots,W_n)
      =\prod_{i=1}^n(1-q_n(W_i)).                \tag{6.25}
\]

Now \(\mathbb E Q_n=1/2\), while

\[
 \mathbb E Q_n^2=(1-p_n-b_n)^n
 \longrightarrow 2^{-(1+\kappa)}.               \tag{6.26}
\]

Therefore

\[
\begin{aligned}
 \liminf_{n\to\infty}\operatorname{err}(S_n)
 &\ge\lim_{n\to\infty}\mathbb E[Q_n(1-Q_n)]\\
 &=\frac12-2^{-(1+\kappa)}\\
 &=0.15143487828269592\ldots .                  \tag{6.27}
\end{aligned}
\]

The absolute height \(t_n\sim\log n\) of a rare exponential coordinate
does not create an almost perfectly observable balanced OR: the
order-one overshoot and the order-one Gaussian noise survive in the limit.

---

## 7. Rare-event maximal-correlation counterexamples

Rare events do invalidate naive maximal-correlation heuristics for general
priors.  For example, let a scalar input have mass \(1-\varepsilon\) near
one point and mass \(\varepsilon\) near a second point at distance
\(\asymp\varepsilon^{-1/2}\), and center and rescale it to variance one.
Unit Gaussian noise identifies the rare component with error
\(e^{-\Theta(1/\varepsilon)}\).  The normalized rare indicator then has
correlation with the observation tending to one.  Smooth two-Gaussian
mixtures have the same property.

This does not refute the present statement for two separate reasons.

1. Such separated mixtures are not log-concave.
2. The almost preserved indicator has mass \(\varepsilon\), not \(1/2\).
   Adding an arbitrary split of the bulk to balance the label makes the
   rare component contribute only \(O(\varepsilon)\) to the unnormalized
   Bayes advantage.  Vanishing balanced error would still require the bulk
   split itself to be nearly observable.

It would be circular to claim that log-concavity universally excludes every
rare-event HGR obstruction: by Section 3, a full universal HGR gap is
equivalent to KLS.  What can be said rigorously is that the standard
log-concave suspects do not exhibit the pathology.  One-dimensional
Gaussian tails become less observable at extreme thresholds; compact
tails are blurred by the noise; and the extremal exponential-type tail has
the explicit nonzero residual constants (6.23)--(6.24).  Tensorization then
settles product log-concave examples once a one-dimensional Poincare bound
is available.

---

## 8. What remains and what would count as progress

The following routes do **not** close the problem on their own:

* posterior 1-strong log-concavity, because it controls the conditional
  boundary from below but not how often \(p(Y)\) lies in transition;
* a universal HGR or chi-square SDPI, because it is quantitatively
  equivalent to a universal Poincare/KLS bound;
* a balanced binary entropy SDPI, because (3.12) turns it back into the same
  Bayes-error question;
* Gaussian rearrangement separately inside each posterior, because the
  rearrangements depend on \(y\) and are incompatible with one fixed set
  \(E\);
* log-concavity of the joint \((X,Y)\) density, because it has flat diagonal
  directions and arbitrary Boolean cuts destroy log-concavity of the two
  testing subdensities.

The concrete new reduction is the posterior-probit geometry:

\[
 \operatorname{err}_\mu(E)=\mathbb E_\nu\Phi(-|F|),
 \qquad F=\Phi^{-1}(q_Y(E)),qquad \|\nabla F\|\le1, \tag{8.1}
\]

together with the fact that \(F\) comes from one common set through the
exponential-family posteriors (1.1).  A useful next lemma must exploit that
last compatibility.  The transition-mass statement in the Verdict is a
strictly stronger, falsifiable target.  Conversely, any proposed
counterexample must, by (4.5)--(4.8), exhibit two almost-half-mass posterior
confidence regions separated by a distance diverging like
\(\sqrt{\log(1/\operatorname{err})}\), while retaining isotropy and global
log-concavity.  None of the Gaussian, cube, simplex, crosspolytope, or
product-exponential tests does so.
