# Multiscale posterior phases in the Gaussian channel

## Executive verdict

Let \(X\sim\mu\) be isotropic and log-concave, let \(S=1_E(X)\) with
\(\mu(E)=1/2\), and observe the signal-to-noise-ratio (SNR) channel

\[
 Y_t=\sqrt t\,X+G,\qquad G\sim N(0,I_n),\qquad t>0.       \tag{0.1}
\]

Write

\[
 p_t(y)=\mathbb P(S=1\mid Y_t=y),\qquad
 F_t(y)=\Phi^{-1}(p_t(y)),                              \tag{0.2}
\]

and let \(I(a)=\phi(\Phi^{-1}(a))\) be the Gaussian isoperimetric
profile.  The fixed-time probit contraction strengthens at all scales to

\[
 \boxed{\ \|\nabla_yF_t(y)\|\le1\quad\text{for every }t>0.\ } \tag{0.3}
\]

In the equivalent heat channel

\[
 Z_s=X+\sqrt s\,G,\qquad s=t^{-1},                     \tag{0.4}
\]

the same statement is

\[
 \|\nabla_zF_s(z)\|\le s^{-1/2}.                       \tag{0.5}
\]

The main new multiscale functional is

\[
 \boxed{\quad
 \mathscr P_E(t):=\sqrt t\,\mathbb E I(p_t(Y_t)).
 \quad}                                                \tag{0.6}
\]

It is nondecreasing, and its derivative is the exact probit-slope deficit

\[
 \boxed{\quad
 \mathscr P_E'(t)=\frac1{2\sqrt t}\,
 \mathbb E\!\left[I(p_t)(1-\|\nabla F_t\|^2)\right]\ge0.
 \quad}                                                \tag{0.7}
\]

For every finite-perimeter \(E\),

\[
 \mathscr P_E(0+)=0,
 \qquad
 \mathscr P_E(\infty)=\mu^+(E).                       \tag{0.8}
\]

Thus perimeter has the exact scale decomposition

\[
 \boxed{\quad
 \mu^+(E)=\int_0^\infty\frac1{2\sqrt t}
 \mathbb E\!\left[I(p_t)(1-\|\nabla F_t\|^2)\right]dt.
 \quad}                                                \tag{0.9}
\]

There is a second exact multiscale object.  Let

\[
 B_E(t):=\mathbb E[p_t(Y_t)(1-p_t(Y_t))]
 =\frac14\operatorname{mmse}(2S-1\mid Y_t).            \tag{0.10}
\]

Then

\[
 -B_E'(t)=\frac1t\mathbb E\|\nabla p_t\|^2
 =\frac1t\mathbb E[I(p_t)^2\|\nabla F_t\|^2],         \tag{0.11}
\]

and, because \(B_E(0+)=1/4\) and \(B_E(\infty)=0\),

\[
 \boxed{\quad
 \frac14=\int_0^\infty\frac1t
 \mathbb E[I(p_t)^2\|\nabla F_t\|^2]dt.
 \quad}                                                \tag{0.12}
\]

Equivalently,

\[
 d\omega_E(t):=4\,\frac1t
 \mathbb E[I(p_t)^2\|\nabla F_t\|^2]dt               \tag{0.13}
\]

is a probability measure on log-SNR scale, and

\[
 B_E(t_0)=\frac14\omega_E([t_0,\infty)).               \tag{0.14}
\]

This identifies the precise obstruction.  Integrating over all scales
fixes the **total** Boolean uncertainty dissipation, but gives no control
on where its log-SNR measure \(\omega_E\) is located.  A universal Bayes
lower bound at \(t=1\) is exactly a universal lower bound on the tail
\(\omega_E([1,\infty))\).  Neither (0.3), (0.7), nor the posterior
covariance constraints below prevent the action from moving to earlier
scales.  Proving that isotropy and log-concavity do prevent this is KLS
content.

This loss is realized exactly by a log-concave model before isotropic
normalization: for \(X_L\sim\mathrm{Unif}[-L,L]\) and
\(E=\{X_L\ge0\}\),

\[
 \operatorname{err}_L(t)
 =\Phi(-L\sqrt t)+
   \frac{\phi(0)-\phi(L\sqrt t)}{L\sqrt t}
 \sim\frac{\phi(0)}{L\sqrt t}.                         \tag{0.15}
\]

Its entire relaxation profile is translated by \(-2\log L\) on the
log-SNR axis.  The posterior is strongly log-concave and obeys every local
identity in this report, while its fixed-time error tends to zero.
Whitening forces \(L=\sqrt3\), so this is not an isotropic counterexample;
an isotropic construction with the same delay would disprove KLS.  The
example pinpoints the missing global input: local posterior curvature
controls transition width, but not transition mass or relaxation scale.

---

## 1. Two equivalent parameterizations

Let \(f=e^{-V}\) be the density of \(\mu\), where \(V\) is convex in the
extended-valued sense.  Denote by \(\nu_t\) the density of \(Y_t\).  The
posterior is

\[
 q_{t,y}(dx)
 =\frac{f(x)\phi(y-\sqrt t\,x)dx}{\nu_t(y)}
 \propto
 \exp\left[-V(x)-\frac t2\|x\|^2+\sqrt t\,y\cdot x\right]dx. \tag{1.1}
\]

It is \(t\)-strongly log-concave.

For evolution calculations it is cleaner to use

\[
 Z_s=X+\sqrt s\,G,\qquad v_s=f*\gamma_s,              \tag{1.2}
\]

where \(\gamma_s\) is the \(N(0,sI)\) density.  Its posterior is

\[
 q_{s,z}(dx)
 \propto\exp\left[-V(x)-\frac1{2s}\|x-z\|^2\right]dx, \tag{1.3}
\]

which is \(s^{-1}\)-strongly log-concave.  The parameterizations are
related by

\[
 s=t^{-1},\qquad z=y/\sqrt t,\qquad
 p_t(y)=p_s(z),\qquad F_t(y)=F_s(z).                    \tag{1.4}
\]

Both label subdensities and the total output density solve the same heat
equation in \(s\):

\[
 \partial_sv_s=\frac12\Delta v_s.                     \tag{1.5}
\]

In SNR coordinates the corresponding anti-Ornstein--Uhlenbeck equation is

\[
 \partial_t\nu_t
 =-\frac1{2t}\left(\Delta\nu_t+\nabla\cdot(y\nu_t)\right). \tag{1.6}
\]

Equation (1.6) follows directly from
\(\nu_t(y)=t^{-n/2}v_{1/t}(y/\sqrt t)\).

---

## 2. Exact posterior moment calculus

All expectations in this section are under \(q_{t,y}\).  Put

\[
 m=\mathbb E[X],\qquad C=\operatorname{Cov}(X),
 \qquad p=\mathbb E[S],
 \qquad a=\operatorname{Cov}(S,X).                    \tag{2.1}
\]

Let \(m_i=\mathbb E[X\mid S=i,Y_t=y]\) and
\(C_i=\operatorname{Cov}(X\mid S=i,Y_t=y)\).  Then

\[
 a=p(1-p)(m_1-m_0),                                  \tag{2.2}
\]

and total covariance decomposes as

\[
 C=pC_1+(1-p)C_0+p(1-p)(m_1-m_0)(m_1-m_0)^T.         \tag{2.3}
\]

### 2.1 Output score and posterior covariance

Differentiating the exponential family gives the exact Tweedie identities

\[
 \nabla_y\log\nu_t(y)=-y+\sqrt t\,m,                 \tag{2.4}
\]

\[
 \nabla_y m=\sqrt t\,C,                              \tag{2.5}
\]

and

\[
 \nabla_y^2\log\nu_t=-I+tC.                         \tag{2.6}
\]

Brascamp--Lieb and log-concavity of convolution give

\[
 0\preceq C\preceq t^{-1}I,
 \qquad
 -I\preceq\nabla^2\log\nu_t\preceq0.               \tag{2.7}
\]

Thus the posterior mean map has Lipschitz constant at most \(t^{-1/2}\).
In heat coordinates the same identities are

\[
 m_s(z)=z+s\nabla\log v_s(z),                        \tag{2.8}
\]

\[
 C_s=sI+s^2\nabla^2\log v_s,\qquad
 0\preceq C_s\preceq sI,\qquad
 -s^{-1}I\preceq\nabla^2\log v_s\preceq0.           \tag{2.9}
\]

Put \(\widetilde X=X-m\).  The next spatial derivative is the posterior
third cumulant:

\[
 \partial_{y_k}C_{ij}=\sqrt t\,
 \mathbb E[\widetilde X_i\widetilde X_j\widetilde X_k]. \tag{2.10}
\]

There is no sign constraint on (2.10).

### 2.2 Gradient and Hessian of the posterior label

Differentiation gives

\[
 \boxed{\quad \nabla_y p_t=\sqrt t\,a.\quad}         \tag{2.11}
\]

If

\[
 M:=\mathbb E[(S-p)(X-m)(X-m)^T],                    \tag{2.12}
\]

then

\[
 \boxed{\quad \nabla_y^2p_t=tM.\quad}               \tag{2.13}
\]

There is a useful equivalent class-conditional form.  Let

\[
 L_t(y)=\log\frac{p_t(y)}{1-p_t(y)}.                 \tag{2.14}
\]

Then

\[
 \boxed{\quad
 \nabla L_t=\sqrt t\,(m_1-m_0),
 \qquad
 \nabla^2L_t=t(C_1-C_0).
 \quad}                                              \tag{2.15}
\]

Consequently,

\[
 \nabla^2p_t
 =tp(1-p)\left[C_1-C_0+(1-2p)(m_1-m_0)(m_1-m_0)^T\right]. \tag{2.16}
\]

The conditioned posteriors \(q_{t,y}(\cdot\mid S=i)\) need not be
log-concave because \(E\) is arbitrary.  Nevertheless, (2.3) and (2.7)
give the weighted constraints

\[
 C_1\preceq\frac1{tp}I,\qquad
 C_0\preceq\frac1{t(1-p)}I,                         \tag{2.17}
\]

and

\[
 t\,p(1-p)\|m_1-m_0\|^2\le1.                       \tag{2.18}
\]

In particular,

\[
 \|\nabla^2L_t\|_{\rm op}
 \le\frac1{\min(p,1-p)}.                            \tag{2.19}
\]

The blow-up in (2.19) near a rare posterior phase is real; posterior
strong log-concavity does not pass to an arbitrary Boolean conditioning.

### 2.3 Probit derivatives and central-band Hessian control

Since \(p=\Phi(F)\),

\[
 \nabla p=I(p)\nabla F,                              \tag{2.20}
\]

and

\[
 \nabla^2p=I(p)\left(\nabla^2F-F\nabla F\nabla F^T\right). \tag{2.21}
\]

Combining (2.13) and (2.21),

\[
 \boxed{\quad
 \nabla^2F_t=\frac{tM}{I(p)}+F_t\nabla F_t\nabla F_t^T.
 \quad}                                              \tag{2.22}
\]

Equivalently, using (2.16),

\[
\begin{aligned}
 \nabla^2F_t
 &=\frac{tp(1-p)}{I(p)}
   \left[C_1-C_0+(1-2p)(m_1-m_0)(m_1-m_0)^T\right]\\
 &\qquad+F_t\nabla F_t\nabla F_t^T.                \tag{2.23}
\end{aligned}
\]

Once the gradient bound in Section 3 is used, (2.17)--(2.18) give the
explicit central-band estimate

\[
 \|\nabla^2F_t\|_{\rm op}
 \le
 \frac{\max(p,1-p)+|1-2p|}{I(p)}+|\Phi^{-1}(p)|.     \tag{2.24}
\]

Thus \(F_t\) has a dimension- and time-independent Hessian bound on every
fixed band \(\delta\le p\le1-\delta\).  In heat coordinates the right side
of (2.24) is multiplied by \(s^{-1}\).

### 2.4 Exact time derivatives of posterior moments

At fixed \(y\), define the likelihood-time score

\[
 \ell_t(x;y)=-\frac12\|x\|^2+
              \frac{y\cdot x}{2\sqrt t}.             \tag{2.25}
\]

Then

\[
 \partial_t\log\nu_t(y)=\mathbb E\ell_t(X;y),         \tag{2.26}
\]

and, for every posterior observable \(H(X)\) independent of \(t\),

\[
 \partial_t\mathbb E[H(X)]
 =\operatorname{Cov}(H(X),\ell_t(X;y)).               \tag{2.27}
\]

In particular,

\[
 \partial_tp_t(y)=\operatorname{Cov}(S,\ell_t),
 \qquad
 \partial_tF_t(y)=\frac{\operatorname{Cov}(S,\ell_t)}{I(p_t)}, \tag{2.28}
\]

\[
 \partial_tm=\operatorname{Cov}(X,\ell_t),\qquad
 \partial_tC=\operatorname{Cov}((X-m)(X-m)^T,\ell_t). \tag{2.29}
\]

These identities are exact but have no favorable sign; the score contains
both a radial quadratic term and a moving linear term.

---

## 3. Scale-sharp probit contraction and phase geometry

### Lemma 3.1 (conditional centroid bound at curvature \(t\))

If \(Q\) is \(t\)-strongly log-concave and \(Q(A)=p\), then

\[
 \left\|\operatorname{Cov}_Q(1_A(X),X)\right\|
 \le\frac{I(p)}{\sqrt t}.                            \tag{3.1}
\]

#### Proof

After scaling by \(\sqrt t\), it is enough to take \(t=1\).  For a unit
direction \(u\), the marginal \(u\cdot X\) is 1-strongly log-concave.
Its monotone transport from a standard normal is 1-Lipschitz.  Among all
conditional weights in \([0,1]\) of mean \(p\), the covariance with this
marginal is maximized by an upper-tail indicator.  Coupling the upper and
lower Gaussian tails through the monotone contraction bounds their mean
gap by the Gaussian one.  The resulting covariance is
\(\phi(\Phi^{-1}(p))=I(p)\).  Taking the supremum over \(u\) proves (3.1).
\(\square\)

Applying (3.1) to (2.11) gives

\[
 \|\nabla p_t\|\le I(p_t),\qquad
 \boxed{\|\nabla F_t\|\le1}.                       \tag{3.2}
\]

In heat coordinates,

\[
 \|\nabla_zp_s\|\le\frac{I(p_s)}{\sqrt s},
 \qquad
 \boxed{\|\nabla_zF_s\|\le s^{-1/2}}.             \tag{3.3}
\]

The class-centroid form is

\[
 p(1-p)\|m_1-m_0\|
 \le\frac{I(p)}{\sqrt t},                          \tag{3.4}
\]

which strictly improves (2.18), especially in posterior tails.

### Phase separation

For \(0<\alpha<\beta<1\), define

\[
 A_\alpha(t)=\{y:p_t(y)\le\alpha\},\qquad
 B_\beta(t)=\{y:p_t(y)\ge\beta\}.                 \tag{3.5}
\]

Then (3.2) gives the exact Gaussian-quantile separation

\[
 \boxed{\quad
 \operatorname{dist}(A_\alpha(t),B_\beta(t))
 \ge\Phi^{-1}(\beta)-\Phi^{-1}(\alpha).
 \quad}                                             \tag{3.6}
\]

In heat coordinates the right side is multiplied by \(\sqrt s\).

Let

\[
 e_t=\mathbb E\min(p_t(Y_t),1-p_t(Y_t)).             \tag{3.7}
\]

Fix \(a>0\), put \(q=\Phi(-a)\), and define

\[
 H_+(t)=\{F_t\ge a\},\quad
 H_-(t)=\{F_t\le-a\},\quad
 R_a(t)=\{|F_t|<a\}.                                \tag{3.8}
\]

Then

\[
 \nu_t(R_a(t))\le\frac{e_t}{q},                     \tag{3.9}
\]

and balance gives

\[
 \nu_t(H_\pm(t))
 \ge\frac{1/2-q}{1-q}-\frac{e_t}{q}.               \tag{3.10}
\]

Moreover,

\[
 \operatorname{dist}(H_+(t),H_-(t))\ge2a.           \tag{3.11}
\]

Taking \(q=\sqrt{e_t}\) shows that vanishing Bayes error creates two
almost-half-mass output phases separated by

\[
 2\Phi^{-1}(1-\sqrt{e_t})
 \sim2\sqrt{\log(1/e_t)}.                           \tag{3.12}
\]

The multiscale estimates fix the transition **width** at every SNR.  They
do not lower-bound the output mass in that transition.

---

## 4. Evolution equations

### 4.1 Heat evolution of \(p_s\) and \(F_s\)

Let \(u_s=(1_Ef)*\gamma_s\), so \(p_s=u_s/v_s\).  Since both \(u_s\) and
\(v_s\) solve the heat equation,

\[
 \boxed{\quad
 \partial_sp_s
 =\frac12\Delta p_s+\nabla\log v_s\cdot\nabla p_s.
 \quad}                                              \tag{4.1}
\]

Using \(p_s=\Phi(F_s)\),

\[
 \boxed{\quad
 \partial_sF_s
 =\frac12\Delta F_s+\nabla\log v_s\cdot\nabla F_s
  -\frac12F_s\|\nabla F_s\|^2.
 \quad}                                              \tag{4.2}
\]

Thus probit converts the quotient heat equation into a viscous transport
equation with the cubic Hamilton--Jacobi term
\(-F\|\nabla F\|^2/2\).

Let

\[
 \mathcal L_s=\frac12\Delta+\nabla\log v_s\cdot\nabla,
 \qquad g=\nabla F_s,\qquad H=\nabla^2F_s.            \tag{4.3}
\]

A direct Bochner calculation gives

\[
\begin{aligned}
 (\partial_s-\mathcal L_s)\|g\|^2
 &=2g^T(\nabla^2\log v_s)g-\|H\|_{\rm HS}^2
   -\|g\|^4-2F_sg^THg\\
 &=2g^T(\nabla^2\log v_s)g
   -\|H+F_sgg^T\|_{\rm HS}^2
   -(1-F_s^2)\|g\|^4.                              \tag{4.4}
\end{aligned}
\]

Because \(\nabla^2\log v_s\preceq0\), (4.4) implies the central-phase
subsolution property

\[
 (\partial_s-\mathcal L_s)\|\nabla F_s\|^2\le0
 \qquad\text{on }\{|F_s|\le1\}.                    \tag{4.5}
\]

This is stronger than a static Lipschitz estimate: in the interior of the
unit probit band, forward heat flow cannot create gradient-square maxima.
The term
\((F_s^2-1)\|g\|^4\) changes sign in the posterior tails, so (4.5) does not
extend globally.  Most importantly, it still says nothing about how much
\(v_s\)-mass enters \(|F_s|\le1\).

### 4.2 SNR evolution

Transforming (4.1) using (1.4), or taking a quotient in (1.6), gives

\[
\boxed{\quad
 \partial_tp_t
 =-\frac1{2t}\Delta p_t
  -\frac1t\left(\nabla\log\nu_t+\frac y2\right)\cdot\nabla p_t.
\quad}                                               \tag{4.6}
\]

The probit equation is

\[
\boxed{\quad
 \partial_tF_t
 =-\frac1{2t}\Delta F_t
  -\frac1t\left(\nabla\log\nu_t+\frac y2\right)\cdot\nabla F_t
  +\frac{F_t}{2t}\|\nabla F_t\|^2.
\quad}                                               \tag{4.7}
\]

The backward-parabolic sign in SNR time is expected: increasing \(t\)
reveals rather than blurs the signal.

---

## 5. A master dissipation identity

Let \(\Psi:(0,1)\to\mathbb R\) be \(C^2\) with sufficient integrability,
and define

\[
 \mathcal M_\Psi(s)=\int v_s(z)\Psi(p_s(z))dz.       \tag{5.1}
\]

Using (1.5), (4.1), and integration by parts gives

\[
 \boxed{\quad
 \frac d{ds}\mathcal M_\Psi(s)
 =-\frac12\mathbb E_{v_s}
   [\Psi''(p_s)\|\nabla_zp_s\|^2].
 \quad}                                              \tag{5.2}
\]

In SNR time this becomes

\[
 \boxed{\quad
 \frac d{dt}\mathbb E_{\nu_t}\Psi(p_t)
 =\frac1{2t}\mathbb E_{\nu_t}
   [\Psi''(p_t)\|\nabla_yp_t\|^2].
 \quad}                                              \tag{5.3}
\]

This one identity yields all of the useful multiscale laws.

### 5.1 Boolean MMSE/action

For \(\Psi(p)=p(1-p)\), (5.3) gives

\[
 B_E'(t)=-\frac1t\mathbb E\|\nabla p_t\|^2.          \tag{5.4}
\]

Since \(\nabla p_t=I(p_t)\nabla F_t\), this is (0.11).  Integrating from
\(t_0\) to infinity gives

\[
 B_E(t_0)=\int_{t_0}^\infty\frac1t
 \mathbb E[I(p_t)^2\|\nabla F_t\|^2]dt.             \tag{5.5}
\]

The Bayes error satisfies

\[
 B_E(t)\le e_t\le2B_E(t),                           \tag{5.6}
\]

so (5.5) is an exact multiscale representation of the desired fixed-time
quantity up to a factor two.

### 5.2 Conditional entropy

For binary entropy \(h_2\), \(h_2''(p)=-1/[p(1-p)]\).  Therefore

\[
 \frac d{dt}H(S\mid Y_t)
 =-\frac1{2t}\mathbb E
   \frac{\|\nabla p_t\|^2}{p_t(1-p_t)}.              \tag{5.7}
\]

This is the exact binary entropy dissipation.  It has the same scale
localization problem as (5.4).

### 5.3 Bayes flux through the decision surface

For \(\Psi(p)=\min(p,1-p)\), interpreted distributionally,
\(\Psi''=-2\delta_{1/2}\).  Hence

\[
 \boxed{\quad
 e_t'=-\frac1t\int_{\{p_t=1/2\}}
          \nu_t(y)\|\nabla p_t(y)\|\,d\mathcal H_{n-1}(y).
 \quad}                                              \tag{5.8}
\]

At the decision surface, \(\|\nabla p_t\|=\phi(0)\|\nabla F_t\|\).
Thus Bayes error is lost exactly through posterior probability flux across
the moving zero-probit surface.

### 5.4 Gaussian-profile uncertainty and perimeter

The Gaussian profile obeys

\[
 I''(p)=-\frac1{I(p)}.                               \tag{5.9}
\]

Let

\[
 A_E(t)=\mathbb E_{\nu_t}I(p_t).                     \tag{5.10}
\]

Then (5.3) gives

\[
 A_E'(t)=-\frac1{2t}\mathbb E[I(p_t)\|\nabla F_t\|^2]. \tag{5.11}
\]

Consequently,

\[
\begin{aligned}
 \frac d{dt}[\sqrt t\,A_E(t)]
 &=\frac1{2\sqrt t}\mathbb E
   [I(p_t)(1-\|\nabla F_t\|^2)]\\
 &\ge0,                                             \tag{5.12}
\end{aligned}
\]

which proves (0.7).

For a regular finite-perimeter set, the local halfspace blow-up as
\(s\downarrow0\) gives

\[
 \lim_{s\downarrow0}\frac1{\sqrt s}
 \mathbb E_{v_s}I(p_s)=\mu^+(E).                    \tag{5.13}
\]

Indeed, in signed normal coordinate \(r\sqrt s\), the posterior tends to
\(\Phi(r)\), and

\[
 \int_{\mathbb R}I(\Phi(r))dr
 =\int_{\mathbb R}\phi(r)dr=1.                     \tag{5.14}
\]

Approximation extends (5.13) to finite-perimeter sets.  Equations
(5.12)--(5.14) give (0.8)--(0.9).  In heat time the equivalent monotonicity
is

\[
 s\longmapsto\frac1{\sqrt s}\mathbb E_{v_s}I(p_s)
 \quad\text{is nonincreasing}.                      \tag{5.15}
\]

For comparison, the small-noise boundary asymptotics are

\[
 B_E(s)\sim\sqrt{\frac s\pi}\,\mu^+(E),\qquad
 e_E(s)\sim\sqrt{\frac{2s}\pi}\,\mu^+(E),           \tag{5.16}
\]

using

\[
 \int_{\mathbb R}\Phi(r)\Phi(-r)dr=\frac1{\sqrt\pi}. \tag{5.17}
\]

The functional \(\mathscr P_E\) is therefore a true multiscale
renormalization of perimeter, not merely another form of Bayes error.

---

## 6. Why integrating over scales does not close the fixed-time bound

First, every positive numerical SNR is already KLS-hard.  Posterior
Gaussian isoperimetry at curvature \(t\), followed by exact perimeter
averaging, gives

\[
 \mu^+(E)
 =\mathbb E_{\nu_t}q_{t,Y_t}^+(E)
 \ge\sqrt t\,\mathbb E I(p_t)
 \ge\sqrt{\frac{2t}{\pi}}\,e_t.                    \tag{6.1}
\]

Thus a universal lower bound on \(e_{t_0}\) at any fixed \(t_0>0\) would
prove KLS.  Conversely, if \(C_P(\mu)=C\), then
\(C_P(\nu_t)\le tC+1\).  The variance decomposition and
\(\|\nabla p_t\|^2\le p_t(1-p_t)\) give

\[
 e_t\ge B_E(t)\ge\frac1{4(tC+2)}.                  \tag{6.2}
\]

So a dimension-free KLS/Poincare bound supplies a fixed-time bound at
every numerical SNR.  Varying \(t\) does not weaken the logical barrier;
the only possible gain would be a new theorem localizing the multiscale
action using isotropy.

Equation (0.12) says that every balanced Boolean label has total action
\(1/4\), but the desired uncertainty at SNR one is only the tail

\[
 B_E(1)=\int_1^\infty\frac1t
 \mathbb E[I(p_t)^2\|\nabla F_t\|^2]dt.              \tag{6.3}
\]

There is always a cut-dependent median relaxation scale \(t_*(E)\) defined
by \(B_E(t_*)=1/8\), or equivalently
\(\omega_E([t_*,\infty))=1/2\).  The action identity proves existence and
nothing more: it contains no dimension-free lower or upper bound on
\(t_*(E)\).  The delay model in Section 7 has
\(t_*(E_L)=L^{-2}t_*(E_1)\) exactly.

The local constraints allow three independent losses.

1. **Scale location.**  The probability measure \(\omega_E\) may, a
   priori, concentrate on \(t\ll1\).  Total action does not control its
   tail at one.
2. **Transition mass.**  The bound \(\|\nabla F_t\|\le1\) controls how fast
   posterior confidence can change in space, but \(\nu_t\) may put little
   mass in every fixed probit band.
3. **Tail Hessians.**  The central-band evolution (4.5) is favorable, but
   class conditioning destroys strong log-concavity in the tails, as
   quantified by (2.19).  There is no global Hessian maximum principle
   which pushes mass back into the central band.

The profile monotonicity does not repair the first loss.  It decomposes
the original perimeter into the complementary slope deficit
\(1-\|\nabla F_t\|^2\), while Boolean action uses
\(\|\nabla F_t\|^2\):

\[
\begin{array}{ll}
 \text{perimeter density:}&
 \displaystyle \frac{I(p_t)}{2\sqrt t}
 (1-\|\nabla F_t\|^2),\\[6pt]
 \text{MMSE-action density:}&
 \displaystyle \frac{I(p_t)^2}{t}\|\nabla F_t\|^2.
\end{array}                                          \tag{6.4}
\]

A nearly Gaussian transition has \(\|\nabla F_t\|\approx1\): it spends
Boolean action but almost no perimeter-defect action at that scale.  A flat
posterior does the reverse.  Without a lower bound on transition mass,
(6.4) gives no scale pinning.

### What isotropy gives for free

For isotropic \(X\),

\[
 \operatorname{Cov}(Y_t)=(1+t)I.                    \tag{6.5}
\]

This yields only a dimension-dependent low-SNR guarantee.  Gaussian
capacity gives

\[
 I(S;Y_t)\le I(X;Y_t)
 \le\frac n2\log(1+t).                              \tag{6.6}
\]

For equal priors, Pinsker and the Jensen--Shannon representation imply

\[
 I(S;Y_t)\ge\frac12(1-2e_t)^2.                      \tag{6.7}
\]

At \(t=1/(2n)\), (6.6)--(6.7) give the explicit bound

\[
 e_{1/(2n)}
 \ge\frac12\left(1-\frac1{\sqrt2}\right)
 =0.1464466094\ldots .                              \tag{6.8}
\]

This uses no log-concavity, but the SNR scale is \(1/n\), not a positive
universal numerical time.  Upgrading (6.8) to \(t=1\) for every isotropic
log-concave law is exactly the unresolved geometric step.

---

## 7. A sharp delay model showing the precise loss

Let

\[
 X_L\sim\mathrm{Unif}[-L,L],\qquad E_L=\{X_L\ge0\}. \tag{7.1}
\]

This is a balanced log-concave labeled model.  By symmetry and the monotone
likelihood-ratio property of the Gaussian kernel, the Bayes rule for
\(Y_t=\sqrt tX_L+G\) is \(Y_t\ge0\).  Therefore

\[
\begin{aligned}
 e_L(t)
 &=\mathbb E\Phi(-\sqrt t|X_L|)\\
 &=\frac1L\int_0^L\Phi(-\sqrt t x)dx\\
 &=\Phi(-r)+\frac{\phi(0)-\phi(r)}r,\qquad
 r=L\sqrt t.                                       \tag{7.2}
\end{aligned}
\]

Hence

\[
 e_L(t)\sim\frac{\phi(0)}{L\sqrt t}
 \quad\text{when }L\sqrt t\to\infty,              \tag{7.3}
\]

while \(e_L(t)\) is order one only when \(t=O(L^{-2})\).

The scale translation is exact, not just asymptotic.  If
\(X_L=LX_1\), with the corresponding scaled label, then every posterior
functional satisfies

\[
 e_L(t)=e_1(L^2t),\qquad
 B_L(t)=B_1(L^2t).                                  \tag{7.4}
\]

Consequently, on logarithmic SNR \(r=\log t\),

\[
 -\frac d{dr}B_L(e^r)
 =-\frac d{dr}B_1(e^{r+2\log L}).                   \tag{7.5}
\]

Thus \(\omega_{E_L}\) is precisely the translate of \(\omega_{E_1}\) by
\(-2\log L\).  All of the fixed total action in (0.12) moves below SNR one
as \(L\to\infty\).

The posterior in this example is \(t\)-strongly log-concave on its convex
support, (0.3) holds, and for \(L\sqrt t\gg1\) the transition far from the
endpoints is locally the flat-prior Gaussian transition

\[
 p_t(y)\approx\Phi(y),\qquad F_t(y)\approx y,        \tag{7.6}
\]

which nearly saturates \(\|\nabla F_t\|\le1\).  The failure is solely that
the output density in this unit-width transition is
\(\asymp(L\sqrt t)^{-1}\).

The variance is \(L^2/3\).  Whitening sends

\[
 \widetilde X_L=\frac{\sqrt3}{L}X_L
 \sim\mathrm{Unif}[-\sqrt3,\sqrt3],                \tag{7.7}
\]

and translates the relaxation scale back to order one.  Thus this model
does not contradict the isotropic conjecture.  It proves something more
limited but exact: posterior curvature, probit contraction, Hessian control,
and all-scale conservation laws by themselves do not pin the dissipation
scale.  Covariance normalization pins linear dilation; ruling out a
nonlinear, high-dimensional analogue of (7.5) is precisely the missing KLS
geometry.  No isotropic log-concave model with \(e_1\to0\) is known, since
such a model would disprove KLS through posterior perimeter averaging.

---

## 8. Gaussian halfspace check

Let \(X\sim N(0,I_n)\) and \(E=\{X_1\ge0\}\).  Then

\[
 p_t(y)=\Phi\left(\sqrt{\frac t{1+t}}\,y_1\right),
 \qquad
 F_t(y)=\sqrt{\frac t{1+t}}\,y_1.                  \tag{8.1}
\]

Thus

\[
 \|\nabla F_t\|^2=\frac t{1+t},                    \tag{8.2}
\]

and

\[
 e_t=\frac1\pi\arccos\sqrt{\frac t{1+t}}.          \tag{8.3}
\]

The profile functional is explicit:

\[
 A_E(t)=\mathbb E I(p_t)=\frac{\phi(0)}{\sqrt{1+t}}, \tag{8.4}
\]

\[
 \mathscr P_E(t)=\phi(0)\sqrt{\frac t{1+t}}
 \nearrow\phi(0)=\gamma_1^+(\mathbb R_+).          \tag{8.5}
\]

Equation (0.7) is an equality after substituting (8.2)--(8.4).  The
Boolean-MMSE functional is

\[
 B_E(t)=\frac14\left[1-\frac2\pi
 \arcsin\left(\frac t{1+t}\right)\right].           \tag{8.6}
\]

This example confirms every endpoint and normalization above.  It also
shows that even the Gaussian extremizer spreads action over a continuum of
SNR scales; the multiscale laws are not a disguised one-time inequality.

---

## 9. Final assessment and a genuinely stronger next target

The multiscale calculation produces two new exact measures:

1. the perimeter-defect measure from (0.9), which records where posterior
   probit slopes fail to saturate the Gaussian centroid inequality; and
2. the Boolean-action probability measure \(\omega_E\) from (0.13), which
   records where the label's MMSE is dissipated on log-SNR scale.

They are complementary through the factors
\(1-\|\nabla F_t\|^2\) and \(\|\nabla F_t\|^2\), but neither controls the
common transition mass \(I(p_t)\).  All-scale integration therefore does
not yield a universal lower bound at \(t=1\) without a new geometric input.

A structural target stronger than fixed-time Bayes contraction is a
**scale-localization theorem**: find universal \(c,C>0\) such that, for
every isotropic log-concave \(\mu\) and balanced \(E\),

\[
 \omega_E([C^{-1},C])\ge c                           \tag{9.1}
\]

and simultaneously

\[
 \inf_{t\in[C^{-1},C]}
 \nu_t\{|F_t|\le1\}\ge c.                           \tag{9.2}
\]

Either statement contains more information than the scalar error at one:
(9.1) localizes uncertainty loss on an entire scale window, while (9.2)
forces actual unit-probit transition mass throughout that window.  The
uniform-interval model shows exactly how both fail without isotropic scale
pinning.  The central-band Bochner inequality (4.5), the profile
monotonicity (0.7), and the class-covariance formulas (2.15)--(2.24) are the
available rigid structures for attacking (9.1)--(9.2).  What remains is a
global argument converting isotropy of \(\nu_t\) into mass in the central
posterior phase; none of the exact evolution identities supplies that
conversion automatically.
