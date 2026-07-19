# Gaussian-profile stochastic localization: the exact monotonicity and the seed barrier

## 1. Executive conclusion

Run ordinary stochastic localization, so that the accumulated quadratic
potential is \(Q_t=tI\).  For a fixed Borel set \(S\), put

\[
 g_t=\mu_t(S),\qquad
 v_t=\operatorname{Cov}_{\mu_t}(\mathbf 1_S,X).
\]

Let

\[
 \mathcal I(s)=\varphi(\Phi^{-1}(s)),\qquad 0<s<1,             \tag{1.1}
\]

be the standard Gaussian isoperimetric profile.  The sharp facts are

\[
 \boxed{\sqrt t\,|v_t|\leq \mathcal I(g_t)}                    \tag{1.2}
\]

and

\[
 \boxed{
 d\big(\sqrt t\,\mathcal I(g_t)\big)
 =\sqrt t\,\mathcal I'(g_t)v_t^T dW_t
 +\frac{\mathcal I(g_t)^2-t|v_t|^2}
        {2\sqrt t\,\mathcal I(g_t)}\,dt .}                   \tag{1.3}
\]

Thus

\[
 F(t):=\sqrt t\,\mathbb E\mathcal I(g_t)                     \tag{1.4}
\]

is nondecreasing.  The power \(1/2\) is critical: for
\(t^\alpha\mathcal I(g_t)\), the worst-case drift is nonnegative exactly
when \(\alpha\geq1/2\).  Gaussian halfspaces attain equality after initial
Gaussian curvature is included.

For a regular set and density, ordinary localization has exact perimeter
transfer and

\[
 F(t)\leq P_\mu(S),\qquad
 \lim_{t\to\infty}F(t)=P_\mu(S).                                \tag{1.5}
\]

Consequently a universal positive lower bound on \(F(t)\), at a deterministic
time or at a bounded randomized stopping time, for every isotropic
log-concave half-set is already a dimension-free Cheeger theorem.  The
profile monotonicity starts from

\[
 \lim_{t\downarrow0}F(t)=0,                                   \tag{1.6}
\]

so it supplies no seed by itself.

Three proposed ways of manufacturing a seed were audited.

* Randomized stopping changes no inequality: for every bounded stopping time
  \(\tau\),

  \[
  \mathbb E[\sqrt\tau\,\mathcal I(g_\tau)]\leq P_\mu(S).       \tag{1.7}
  \]

  A universal lower bound for the left side is therefore precisely a KLS
  seed, not a preliminary estimate.
* Adding deterministic centered Gaussian curvature can destroy the posterior
  mass of a half-set exponentially in dimension.  Randomizing the center in
  the canonical mixture-preserving way is exactly localization to a positive
  time, whose seed is the unknown quantity \(F(t)\).  Removing the curvature
  has an exponentially bad density ratio on the same example.
* A scalar correction using posterior covariance is not closed under Ito
  calculus: third central moments enter its generator.  The natural
  one-dimensional Gaussian-calibrated candidate
  \(\mathcal I(g_t)/\sqrt{A_t}\) has negative drift on explicit exponential
  posterior states, and its expectation cannot be monotone by the isotropic
  uniform-interval example.

No universal seed theorem is proved here.  What is proved is the sharp
profile identity, its exact endpoint meaning, a dimension-dependent trace
seed of order \(n^{-1/2}\), and rigorous obstructions to the three scalar
repairs above.  The missing step is a set-specific occupation or rigidity
estimate preventing the process from entering the critical regime
\(t|v_t|^2\simeq\mathcal I(g_t)^2\) at arbitrarily small physical time.

## 2. Ordinary localization and filtering identities

Initially take a positive smooth log-concave density \(p\) with enough
moments, and stop all natural parameters and moments on compact sets.  Define

\[
 p_t(x)=\frac1{Z_t}
  \exp\left(c_t\cdot x-\frac t2|x|^2\right)p(x),
 \qquad dc_t=dW_t+a_tdt,qquad a_t=\mathbb E_tX.                \tag{2.1}
\]

The usual Ito calculation gives, pointwise in \(x\),

\[
 dp_t(x)=p_t(x)\langle x-a_t,dW_t\rangle.                     \tag{2.2}
\]

Consequently, for

\[
 A_t=\operatorname{Cov}_{\mu_t}(X),\qquad
 v_t=\mathbb E_t[(\mathbf1_S-g_t)(X-a_t)],                    \tag{2.3}
\]

one has

\[
 dg_t=v_t^TdW_t,qquad d\langle g\rangle_t=|v_t|^2dt.         \tag{2.4}
\]

The filtering realization is exact.  If \(X\sim\mu\) and \(B_t\) is an
independent Brownian motion, then

\[
 c_t\ \stackrel{d}=\ tX+B_t,qquad
 \mu_t=\operatorname{Law}(X\mid c_t).                         \tag{2.5}
\]

In particular \(p_t(x)/p(x)\) is the likelihood martingale and

\[
 \mathbb E p_t(x)=p(x).                                       \tag{2.6}
\]

All formulas below are first proved with the compact parameter stop.  They
pass to general full-dimensional log-concave laws by convex truncation,
smoothing, and monotone removal of the stops.  When an unbounded stopping
time is used, only the Fatou direction justified by nonnegative likelihoods
is asserted.

## 3. The sharp strongly-log-concave centroid bound

### Theorem 3.1

Let \(\pi\) be a \(t\)-strongly log-concave probability on
\(\mathbb R^n\), \(t>0\), let \(S\) be Borel, and set

\[
 g=\pi(S),\qquad v=\operatorname{Cov}_\pi(\mathbf1_S,X).
\]

Then

\[
 \sqrt t\,|v|\leq\mathcal I(g).                               \tag{3.1}
\]

The constant is sharp.  Equality is attained by a Gaussian of covariance
\(t^{-1}I\) and a halfspace.

### Proof

Translate so that \(\mathbb E_\pi X=0\).  If \(v=0\), there is nothing to
prove.  Otherwise put \(u=v/|v|\) and \(Y=\langle X,u\rangle\).  The marginal
law \(\rho\) of \(Y\) is \(t\)-strongly log-concave.  Indeed, write the density
of \(\pi\) as

\[
 e^{-t|x|^2/2}h(x)
\]

with \(h\) log-concave.  Integrating the transverse variables and applying
Prékopa's theorem writes the density of \(Y\) as
\(e^{-ty^2/2}\widetilde h(y)\), with \(\widetilde h\) log-concave.

Let \(q(y)=\mathbb P_\pi(S\mid Y=y)\).  Then \(0\leq q\leq1\),
\(\int q\,d\rho=g\), and

\[
 |v|=\int yq(y)\,d\rho(y).                                    \tag{3.2}
\]

The bathtub principle shows that the last integral is maximized, among all
such \(q\), by the upper \(g\)-tail of \(\rho\).  Indeed, if
\(H_Y=\mathbf1_{\{Y\geq y_g\}}\) has \(\rho\)-mass \(g\), then
\((y-y_g)(q(y)-H_Y(y))\leq0\) pointwise and
\(\int(q-H_Y)\,d\rho=0\).  Hence
\(\int yq\,d\rho\leq\int yH_Y\,d\rho\).

Let \(Z\) be standard Gaussian and let \(T\) be the increasing transport
from \(Z\) to \(Y\).  The one-dimensional Caffarelli contraction theorem gives

\[
 0\leq T'(z)\leq L:=t^{-1/2}\quad\text{a.e.}                  \tag{3.3}
\]

Here is the one-dimensional contraction calculation.  In the smooth case
write the density of \(Y\) as \(Z_\rho^{-1}e^{-V(y)}\), where
\(V''\geq t\).  The transport equation
\(\varphi(z)=Z_\rho^{-1}e^{-V(T(z))}T'(z)\) gives
\[
 (\log T')'=-z+V'(T)T',
 \qquad
 (\log T')''=-1+V''(T)(T')^2+V'(T)T''.
\]
At an interior maximum of \(T'\), one has \(T''=0\) and
\((\log T')''\leq0\), so \(t(T')^2\leq1\).  A vanishing quadratic barrier
handles a supremum attained only at infinity.  Approximation of the convex
part of \(V\), followed by convergence of the quantile maps, proves (3.3)
without smoothness.  Put \(R(z)=T(z)-Lz\).  Then \(R\) is nonincreasing.
If \(z_g=\Phi^{-1}(1-g)\) and \(H(z)=\mathbf1_{\{z\geq z_g\}}\), the
independent-copy covariance identity gives

\[
 \operatorname{Cov}(R(Z),H(Z))
 =\frac12\mathbb E[(R(Z)-R(Z'))(H(Z)-H(Z'))]\leq0.            \tag{3.4}
\]

Since \(\mathbb EY=0\), (3.2)--(3.4) imply

\[
 |v|\leq L\operatorname{Cov}(Z,H(Z))
 =L\varphi(z_g)=\frac{\mathcal I(g)}{\sqrt t}.
\]

This proves (3.1).  Every inequality is an equality for the stated Gaussian
halfspace model.  QED.

Applying Theorem 3.1 to \(\mu_t\), whose potential has Hessian at least
\(tI\), proves (1.2).

## 4. The Gaussian profile and the exact Ito formula

Let \(z=\Phi^{-1}(s)\).  Since \(\Phi'(z)=\varphi(z)\) and
\(\varphi'(z)=-z\varphi(z)\), the chain rule gives
\(\mathcal I'(s)=\varphi'(z)/\Phi'(z)=-z\), and a second differentiation
gives

\[
 \mathcal I'(s)=-\Phi^{-1}(s),
 \qquad
 \boxed{\mathcal I''(s)=-\frac1{\mathcal I(s)}}.              \tag{4.1}
\]

Stopping first while \(g_t\in[\varepsilon,1-\varepsilon]\), Ito's formula
and (2.4) give

\[
 d\mathcal I(g_t)
 =\mathcal I'(g_t)v_t^TdW_t
  -\frac{|v_t|^2}{2\mathcal I(g_t)}dt.                         \tag{4.2}
\]

Multiplication by \(\sqrt t\) gives exactly (1.3).  The drift is
nonnegative by (1.2).  After bounded stopping, the stochastic integral is a
true martingale; letting \(\varepsilon\downarrow0\) gives the corresponding
local-submartingale.  On each bounded time interval,
\(0\leq\sqrt t\,\mathcal I(g_t)\leq\sqrt T/\sqrt{2\pi}\), so it is of class
\(D\) and therefore a true submartingale.  This proves expectation
monotonicity.
For \(0<s<t\), one may equivalently write

\[
 F(t)-F(s)
 =\frac12\mathbb E\int_s^t
 \frac{\mathcal I(g_r)^2-r|v_r|^2}
      {\sqrt r\,\mathcal I(g_r)}dr\geq0.                       \tag{4.3}
\]

### 4.1 Criticality of the square root

For an arbitrary exponent \(\alpha\), the finite-variation part of
\(t^\alpha\mathcal I(g_t)\) is

\[
 t^{\alpha-1}\left[
 \alpha\mathcal I(g_t)
 -\frac{t|v_t|^2}{2\mathcal I(g_t)}\right]dt.                  \tag{4.4}
\]

Under only the sharp bound \(t|v_t|^2\leq\mathcal I(g_t)^2\), this is
nonnegative for every state if and only if \(\alpha\geq1/2\).  Thus
\(\alpha=1/2\) is the least, and critical, exponent.  This algebraic
criticality is realized by probability laws.  For
\(\mu=N(0,\kappa^{-1}I)\) and a centered halfspace,
\[
 \mathbb E\mathcal I(g_t)
 =\frac{\sqrt\kappa}{\sqrt{2\pi(\kappa+t)}}.
\]
Hence \(t^\alpha\mathbb E\mathcal I(g_t)\) has logarithmic derivative
\(\alpha/t-1/[2(\kappa+t)]\), which is negative for all sufficiently large
\(t\) whenever \(\alpha<1/2\).

More generally, if the initial law is \(\kappa\)-strongly log-concave, then
the posterior is \((\kappa+t)\)-strongly log-concave and

\[
 d\big(\sqrt{\kappa+t}\,\mathcal I(g_t)\big)
 =\sqrt{\kappa+t}\,\mathcal I'(g_t)v_t^TdW_t
 +\frac{\mathcal I(g_t)^2-(\kappa+t)|v_t|^2}
 {2\sqrt{\kappa+t}\,\mathcal I(g_t)}dt.                       \tag{4.5}
\]

For a Gaussian of covariance \(\kappa^{-1}I\) and a halfspace, the drift in
(4.5) vanishes identically.  This is an honest probability model showing
sharpness, not merely a formal equality case.

## 5. Exact perimeter meaning of the profile

### 5.1 Perimeter transfer

For a smooth set \(S\) and density \(p\), define

\[
 P_\mu(S)=\int_{\partial S}p(x)\,d\mathcal H^{n-1}(x).
\]

Equation (2.6) and Tonelli give

\[
 \mathbb E P_{\mu_t}(S)=P_\mu(S).                              \tag{5.1}
\]

Every \(\mu_t\) is \(t\)-strongly log-concave.  Sharp Gaussian
isoperimetry for strongly log-concave measures, obtained from Caffarelli
contraction, says

\[
 P_{\mu_t}(S)\geq\sqrt t\,\mathcal I(g_t).                    \tag{5.2}
\]

Taking expectations proves the first part of (1.5).  For finite-perimeter
sets, (5.1) follows by integrating \(p_t\) on the reduced boundary.  For
exterior Minkowski content, Fatou gives the direction
\(\mathbb E\mu_t^+(S)\leq\mu^+(S)\), which is all that is needed for the
seed implications below.

The same argument works at a bounded stopping time \(\tau\).  Optional
sampling for the stopped likelihood gives

\[
 P_\mu(S)=\mathbb E P_{\mu_\tau}(S)
 \geq\mathbb E[\sqrt\tau\,\mathcal I(g_\tau)].                 \tag{5.3}
\]

For unbounded \(\tau\), first stop at \(\tau\wedge m\) and use Fatou.

### 5.2 Heat-kernel representation and the endpoint limit

Assume now that \(S\) has compact \(C^2\) boundary, that \(p\) is continuous
and positive on a neighborhood of \(\partial S\), and that no support
boundary meets \(\partial S\).  Put \(s=1/t\), and let \(P_s\) denote
convolution with the centered Gaussian of covariance \(sI\).  From (2.5),
\(Y_t=c_t/t=X+\sqrt sZ\).  Bayes' formula
therefore gives

\[
 \mathbb E\mathcal I(g_t)
 =\int_{\mathbb R^n}P_sp(y)\,
 \mathcal I\left(\frac{P_s(p\mathbf1_S)(y)}{P_sp(y)}\right)dy.
                                                               \tag{5.4}
\]

Consequently

\[
 F(1/s)=\frac1{\sqrt s}\int P_sp\,
 \mathcal I\left(\frac{P_s(p\mathbf1_S)}{P_sp}\right).        \tag{5.5}
\]

In tubular coordinates \(y=x+\sqrt s\,rN(x)\) about the boundary,

\[
 \frac{P_s(p\mathbf1_S)(y)}{P_sp(y)}\longrightarrow\Phi(-r),
 \qquad P_sp(y)\longrightarrow p(x).                          \tag{5.6}
\]

Since \(\mathcal I(\Phi(-r))=\varphi(r)\) and
\(\int_{\mathbb R}\varphi(r)\,dr=1\), dominated localization to the
\(O(\sqrt s)\) boundary tube yields

\[
 \lim_{s\downarrow0}F(1/s)
 =\int_{\partial S}p\,d\mathcal H^{n-1}=P_\mu(S).              \tag{5.7}
\]

This proves the second part of (1.5) under the stated regularity.  The same
formula extends to weighted finite perimeter by the standard blow-up of a BV
set at almost every reduced-boundary point; this extension is not needed for
any counterexample below.

Thus \(F(t)\) is an increasing boundary-layer approximation to the original
perimeter, starting at zero.  A positive lower bound on it is not weaker than
the desired isoperimetry.

## 6. A baseline seed using only trace covariance

It is useful to record what follows without an operator norm.  Suppose
\(\mu\) is isotropic and \(g_0=1/2\).  The mixture identity
\(\mathbb E\mu_t=\mu\) and total covariance give

\[
 \mathbb EA_t+\operatorname{Cov}(a_t)=I,
 \qquad \mathbb E\operatorname{tr}A_t\leq n.                  \tag{6.1}
\]

Covariance Cauchy--Schwarz yields

\[
 |v_t|^2\leq g_t(1-g_t)\|A_t\|_{\rm op}
 \leq\frac14\operatorname{tr}A_t.                             \tag{6.2}
\]

Therefore

\[
 \mathbb E\langle g\rangle_T
 =\int_0^T\mathbb E|v_t|^2dt\leq\frac{nT}{4}.                 \tag{6.3}
\]

Doob's inequality gives

\[
 \mathbb P\left\{\sup_{t\leq T}|g_t-1/2|\geq1/4\right\}
 \leq4nT.                                                     \tag{6.4}
\]

Taking \(T=(8n)^{-1}\), with probability at least \(1/2\) one has
\(g_T\in[1/4,3/4]\).  Hence

\[
 F(T)\geq\frac{\mathcal I(1/4)}{2\sqrt{8n}},
 \qquad P_\mu(S)\geq\frac{c}{\sqrt n}.                        \tag{6.5}
\]

This is a genuine seed theorem without an operator norm, but only at the
classical \(n^{-1/2}\) scale.  Replacing its physical time \(1/n\) by a
universal time is exactly the unresolved issue.

## 7. Why randomized stopping does not create a free seed

The process \(Y_t=\sqrt t\,\mathcal I(g_t)\) is a true submartingale on
every bounded horizon by (1.3) and the bound in Section 4.  For every
bounded stopping time \(\tau\leq T\),
optional sampling and (5.2) give

\[
 \mathbb EY_\tau\leq\mathbb EY_T=F(T)\leq P_\mu(S).            \tag{7.1}
\]

Independent randomization of the stopping rule is covered by enlarging the
initial sigma field.  Thus a claim

\[
 \mathbb E[\sqrt\tau\,\mathcal I(g_\tau)]\geq c>0             \tag{7.2}
\]

for every isotropic log-concave half-set proves the dimension-free Cheeger
inequality immediately.

For example, fix \(0<\delta<1/2\) and let

\[
 \tau_\delta=\inf\{t:g_t\notin(\delta,1-\delta)\}.
\]

At \(\tau=T\wedge\tau_\delta\), retaining only paths with
\(\tau_\delta>T\) in (5.3) gives

\[
 P_\mu(S)\geq\sqrt T\,\mathcal I(\delta)
          \mathbb P\{\tau_\delta>T\}.                         \tag{7.3}
\]

A universal \(T>0\) and survival probability in (7.3) would solve KLS.  The
profile identities do not prove that probability bound; they merely express
the exact reward if it is proved.  Randomizing \(T\) does not change this
logical status.

## 8. Adding or removing initial Gaussian curvature

### 8.1 What works for an already strongly log-concave law

If \(\mu\) is \(\kappa\)-strongly log-concave, (4.5) starts with the positive
seed

\[
 \sqrt\kappa\,\mathcal I(g_0).                                \tag{8.1}
\]

Together with perimeter transfer this reproves the sharp strongly
log-concave isoperimetric inequality.  This is the legitimate use of initial
curvature.

For general \(\mu\), a posterior with curvature \(\kappa I\) and random
linear tilt \(c_\kappa\) is exactly \(\mu_\kappa\) in (2.1).  Averaging its
putative initial seed gives

\[
 \sqrt\kappa\,\mathbb E\mathcal I(g_\kappa)=F(\kappa),        \tag{8.2}
\]

the quantity that needed a lower bound in the first place.  Random centering
therefore preserves the mixture but is circular.

### 8.2 Centered curvature loses half-set mass: an exact Gaussian model

Let \(\mu=\gamma_n=N(0,I_n)\), choose a median \(m_n\) of \(\chi_n^2\), and
put

\[
 S_n=\{x:|x|^2\geq m_n\}.
\]

Then \(\gamma_n(S_n)=1/2\).  Deterministically adding centered curvature
\(\kappa>0\) changes the probability to

\[
 d\gamma_{n,\kappa}(x)
 =(1+\kappa)^{n/2}e^{-\kappa|x|^2/2}d\gamma_n(x)
 =N(0,(1+\kappa)^{-1}I_n),                                   \tag{8.3}
\]

and hence

\[
 \gamma_{n,\kappa}(S_n)
 =\mathbb P\{\chi_n^2\geq(1+\kappa)m_n\}
 \leq e^{-c_\kappa n}.                                      \tag{8.4}
\]

Here \(m_n/n\to1\), and (8.4) follows from the elementary chi-square
Chernoff bound.  Thus the centered-curvature profile
\(\sqrt\kappa\,\mathcal I(\gamma_{n,\kappa}(S_n))\) tends to zero
exponentially, although the original Gaussian perimeter of \(S_n\) stays
bounded above and below by numerical constants (in fact it tends to
\(1/\sqrt\pi\)).
Indeed, if \(R_n=|G|\), then
\[
 P_{\gamma_n}(S_n)=f_{R_n}(\sqrt{m_n}),\qquad
 f_{R_n}(r)=\frac{2^{1-n/2}}{\Gamma(n/2)}r^{n-1}e^{-r^2/2}.
\]
The central limit theorem implies
\((m_n-n)/\sqrt{2n}\to0\) for any choice of median, and Stirling's formula
in the last display gives \(f_{R_n}(\sqrt{m_n})\to1/\sqrt\pi\).

Nor can the curvature be removed with a dimension-free density comparison.
On the boundary \(|x|^2=m_n\),

\[
 \frac{d\gamma_n}{d\gamma_{n,\kappa}}(x)
 =(1+\kappa)^{-n/2}e^{\kappa m_n/2}
 =\exp\left(\frac n2[\kappa-\log(1+\kappa)]+o(n)\right),       \tag{8.5}
\]

which is exponentially large.  This model rules out a deterministic
add-then-remove curvature seed based only on the original half-mass
condition.

## 9. Nonlinear profiles involving posterior covariance

### 9.1 The scalar state is not closed

Let \(Y=X-a_t\), and define the third central moment tensor

\[
 T_{ijk}=\mathbb E_t[Y_iY_jY_k].
\]

Ordinary localization gives

\[
 dA_{ij}=\sum_kT_{ijk}dW_k-(A_t^2)_{ij}dt.                    \tag{9.1}
\]

For a \(C^2\) scalar profile \(\Psi(t,g,A)\), its drift is

\[
\begin{aligned}
 \mathcal L\Psi={}&\partial_t\Psi
 -\langle D_A\Psi,A^2\rangle
 +\frac12\Psi_{gg}|v|^2 \\
 &+\sum_{i,j,k}\Psi_{g,A_{ij}}v_kT_{ijk}
 +\frac12\sum_{i,j,k,\ell,m}
   \Psi_{A_{ij},A_{\ell m}}T_{ijk}T_{\ell mk}.                \tag{9.2}
\end{aligned}
\]

Thus neither \(g\) nor \((g,A)\) is a closed Markov state.  A covariance
multiplier creates mixed third-moment terms with no sign.  Discarding them is
not a valid scalar comparison.

### 9.2 The natural one-dimensional correction has negative drift

In one dimension write

\[
 a=A_t,\qquad \tau_3=\mathbb E_t[(X-a_t)^3],qquad
 J_t=\frac{\mathcal I(g_t)}{\sqrt{a_t}}.
\]

From (2.4), (9.1), and \(d[g,a]_t=v_t\tau_3dt\), direct Ito calculus gives

\[
\begin{aligned}
 \operatorname{drift}(J_t)={}&
 \frac{\mathcal I(g)\sqrt a}{2}
 -\frac{v^2}{2\mathcal I(g)\sqrt a}
 -\frac{\mathcal I'(g)v\tau_3}{2a^{3/2}}
 +\frac{3\mathcal I(g)\tau_3^2}{8a^{5/2}}.                  \tag{9.3}
\end{aligned}
\]

This drift can be negative on an explicit log-concave state.  Let \(Z\) have
the exponential density \(re^{-rz}\mathbf1_{\{z\geq0\}}\), take \(X=Z-1\),
and use the fixed set \(S=\{X\geq L-1\}=\{Z\geq L\}\), where \(L=\log2\).
Then

\[
 g=e^{-rL},\qquad a=r^{-2},\qquad \tau_3=2r^{-3},
 \qquad v=gL.                                                 \tag{9.4}
\]

Writing \(z_g=\Phi^{-1}(g)\), substitution in (9.3) yields

\[
 b_r=\frac{2\mathcal I(g)}r
 -\frac{rg^2L^2}{2\mathcal I(g)}+z_ggL.                       \tag{9.5}
\]

As \(r\to\infty\), the Gaussian Mills asymptotics give

\[
 |z_g|\sim\sqrt{2rL},\qquad
 \mathcal I(g)\sim g\sqrt{2rL},
\]

so the positive term in (9.5) is \(O(g/\sqrt r)\), whereas each negative
term has order \(g\sqrt r\).  Thus \(b_r<0\) for every sufficiently large
\(r\).  These are not artificial
states: \(X=Z-1\) for a rate-one exponential \(Z\) is isotropic, and the
displayed set has initial mass \(1/2\).  Its posteriors with small quadratic
time and linear tilt
near \(1-r\) approximate (9.4), and that open set of natural parameters has
positive probability because of the Brownian tilt.  Hence \(J_t\) is not a
local submartingale.

There is also a global obstruction.  Let \(\mu\) be uniform on
\([-\sqrt3,\sqrt3]\), which is isotropic, and let \(S=[0,\sqrt3]\).  Initially

\[
 J_0=\mathcal I(1/2)=\frac1{\sqrt{2\pi}}.                      \tag{9.6}
\]

Write \(s=1/t\) and condition on \(Y=X+\sqrt sZ\).  In the boundary layer
\(Y=\sqrt s\,r\) about the sole relative boundary point \(0\), the posterior
is, uniformly for bounded \(r\), a Gaussian of variance \(s(1+o(1))\).
Writing \(h_s\) for the density of \(Y\), one has
\[
 h_s(\sqrt s\,r)\longrightarrow\frac1{2\sqrt3},\qquad
 g_t(\sqrt s\,r)\longrightarrow\Phi(r),\qquad
 \frac{A_t(\sqrt s\,r)}s\longrightarrow1.
\]
After the change of variables \(dy=\sqrt s\,dr\), the local integrand
therefore tends to \(\varphi(r)\,dr/(2\sqrt3)\).
Splitting at \(|Y|\leq M\sqrt s\), then sending first \(s\downarrow0\) and
then \(M\to\infty\), proves by Gaussian tail bounds that the complement
contributes zero.  Consequently

\[
 \lim_{t\to\infty}\mathbb EJ_t
 =P_\mu(S)=\frac1{2\sqrt3}
 <\frac1{\sqrt{2\pi}}.                                       \tag{9.7}
\]

Thus the most direct covariance replacement of the missing curvature is not
monotone even in dimension one.

### 9.3 The general scalar-profile barrier

The preceding examples do not prove that no ingenious nonlinear scalar
functional can prove KLS.  Such a statement would be unjustified.  What can
be stated exactly is the following.

Suppose a nonnegative posterior functional \(H_t\), scalar or otherwise,
satisfies for every regular log-concave pair

\[
 \mathbb EH_\tau\leq P_\mu(S)                                 \tag{9.8}
\]

at the stopping rule used, and suppose isotropy and \(\mu(S)=1/2\) imply

\[
 \mathbb EH_\tau\geq c>0.                                    \tag{9.9}
\]

Then (9.8)--(9.9) are already the dimension-free Cheeger inequality.  If
\(H_t\) is designed to converge to the boundary layer in (5.7), a monotonicity
from a positive isotropic initial value is likewise a complete KLS proof.
Posterior covariance does not make the seed (9.9) an easier theorem; its Ito
equation (9.2) identifies the new third-moment information which would have
to prove it.

## 10. Gaussian calibration

Take \(\mu=N(0,I_n)\) and \(S=\{x_1\geq0\}\).  The posterior is

\[
 \mu_t=N\left(\frac{c_t}{1+t},\frac1{1+t}I\right),
 \qquad
 g_t=\Phi\left(\frac{c_{t,1}}{\sqrt{1+t}}\right).              \tag{10.1}
\]

Under the original path law,

\[
 Z_t:=\Phi^{-1}(g_t)\sim N(0,t),qquad
 |v_t|=\frac{\mathcal I(g_t)}{\sqrt{1+t}}.                    \tag{10.2}
\]

Therefore

\[
 \mathbb E\mathcal I(g_t)
 =\frac1{\sqrt{2\pi(1+t)}},
 \qquad
 F(t)=\frac1{\sqrt{2\pi}}\sqrt{\frac t{1+t}}.                \tag{10.3}
\]

The function \(F\) rises from zero to the Gaussian halfspace perimeter.
Meanwhile

\[
 \sqrt{1+t}\,\mathbb E\mathcal I(g_t)=\frac1{\sqrt{2\pi}}    \tag{10.4}
\]

is constant, exactly as (4.5) predicts.  This model verifies every constant
and shows that no improvement of the centroid coefficient or of the critical
power is available.

## 11. Status and the precise missing statement

The Gaussian-profile route removes one genuine source of slack: it replaces
the crude estimate on \(g_t(1-g_t)\) by the sharp pointwise inequality
\(\sqrt t|v_t|\leq\mathcal I(g_t)\).  It also avoids an anisotropic endpoint
transfer, because ordinary perimeter transfers exactly.

It does not create the initial positive scale.  A sufficient new statement
would be any one of the following, with numerical constants and every
isotropic log-concave half-set:

\[
 \sqrt{t_0}\,\mathbb E\mathcal I(g_{t_0})\geq c;              \tag{11.1}
\]

or, for a bounded stopping time,

\[
 \mathbb E[\sqrt\tau\,\mathcal I(g_\tau)]\geq c;              \tag{11.2}
\]

or a covariance-corrected functional satisfying a positive isotropic seed
and an endpoint comparison with perimeter.  By Sections 5, 7, and 9.3, each
is already KLS-strength.

In terms of the exact drift, the unresolved phenomenon is occupation of the
critical set

\[
 \mathcal R_t:=\frac{t|v_t|^2}{\mathcal I(g_t)^2}\simeq1.     \tag{11.3}
\]

The scalar inequality only says \(0\leq\mathcal R_t\leq1\).  A successful
proof would need to use the history and geometry of posteriors reachable from
one isotropic initial law to show that (11.3) cannot persist from arbitrarily
small times, or else show that near equality forces a directly controlled
Gaussian-halfspace structure.  No such occupation or rigidity theorem is
proved here.  Treating it as a consequence of scalar profile monotonicity
would be the missing KLS step.

## 12. Reference for the contraction input

The only named geometric input beyond standard localization calculus is
Caffarelli's contraction theorem: the Brenier map from a Gaussian measure to
a log-concave perturbation of that Gaussian is a contraction.  A convenient
review is A. Saumard and J. Wellner, *Log-concavity and strong log-concavity:
a review*, Statistics Surveys 8 (2014),
[doi:10.1214/14-SS107](https://doi.org/10.1214/14-SS107).  The one-dimensional
argument used in Theorem 3.1 is the monotone specialization of that theorem.
