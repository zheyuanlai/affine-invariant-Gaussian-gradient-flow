# The soft longitudinal driver: an unweighted curvature theorem and the mass-weighted seed obstruction

## 0. Verdict

Fix a set `S` and, at the current posterior, write

\[
 g=\mu_t(S),\qquad G=g(1-g),\qquad
 v=\operatorname {Cov}_t({\bf1}_S,X),\qquad \ell=|v|,
 \qquad u={v\over\ell}.
\]

The letter `\ell` is used for `|v|`; below, `k=\ell^2/G`.  This avoids the
otherwise dangerous ambiguity between `|v|` and `|v|^2`.

Consider the literal driver covariance proposed in the question,

\[
 a={G\over G+\ell^2}={1\over1+k},\qquad
 \Gamma=P_{u^\perp}+aI=(1+a)I-uu^T.                 \tag{0.1}
\]

Its eigenvalue in the `u` direction is `a`, and every transverse eigenvalue
is `1+a`.  Thus transverse Euclidean curvature accumulates at rate at least
one.  The control has the following exact information--curvature identity.
Put

\[
 \beta={\ell^2\over G+\ell^2}={k\over1+k}=1-a,
 \qquad Z_t=\int_0^t a_s\,ds.                         \tag{0.2}
\]

Then

\[
 d\langle g\rangle_t=G_t\beta_t\,dt,
 \qquad Z_t=t-\int_0^t\beta_s\,ds,                   \tag{0.3}
\]

and, for binary Shannon entropy
`H(x)=-x\log x-(1-x)\log(1-x)`,

\[
 dH(g_t)=dM_t-\frac12\beta_t\,dt.                    \tag{0.4}
\]

Moreover the accumulated quadratic potential satisfies, pathwise,

\[
 Q_t=Z_tI+\int_0^tP_{u_s^\perp}\,ds\succeq Z_tI.     \tag{0.5}
\]

Consequently, for every deterministic `T` for which the stopped construction
has been removed,

\[
 \boxed{\quad
 \mathbb E\lambda_{\min}(Q_T)\ge \mathbb EZ_T
 =T-2\{H(g_0)-\mathbb EH(g_T)\}
 \ge T-2H(g_0).
 \quad}                                               \tag{0.6}
\]

In particular, for a half-mass set the right side is positive as soon as
`T>2\log2`.  This completely answers the **unweighted expectation** question
in the affirmative.  It does not prove a Cheeger bound.  Endpoint transfer
requires

\[
 \mathbb E\left[\mathcal I(g_T)
                 \sqrt{\lambda_{\min}(Q_T)}\right],   \tag{0.7}
\]

or the weaker expression with `min(g_T,1-g_T)`, rather than the unweighted
expectation in (0.6).  Curvature in (0.6) may be accumulated on paths on
which `g_T` has already polarized.

There is nevertheless an exact mass-weighted Lyapunov.  If
`q_t=\lambda_{\min}(Q_t)>0`, then

\[
 \boxed{\quad \mathcal I(g_t)\sqrt{q_t}
 \quad\hbox{is a stopped local submartingale}.\quad}   \tag{0.8}
\]

Thus the driver preserves and amplifies every positive curvature seed.  It
does not create a uniformly positive seed by any argument below.  The exact
remaining event is a short initial layer in which `a` could fall from its
isotropic initial value `a_0\ge1/2` to a value of order `q_t`, while `u_t`
locks onto one line.  Excluding that event is already sufficient for KLS.

The added transverse term `aP_{u^\perp}` in (0.1), compared with the minimal
smooth entropy-capped driver, does not touch a coherent exceptional line.
For product states in which `S` depends on one factor, the direction is
fixed and the longitudinal dynamics is exactly one-dimensional.  A balanced
large-variance exponential factor gives an explicit state at which
`a\asymp1/\operatorname {Var}`.  This is a statewise obstruction, not a
counterexample starting from an isotropic law along the same controlled
trajectory.

## 1. Stopped construction and a rigorous zero convention

### 1.1 General controlled posterior

First let `p_0=e^{-V_0}` be smooth and positive on a compact convex set in
its minimal affine hull.  For a predictable positive-semidefinite covariance
`\Gamma_t`, let `C_t=\Gamma_t^{1/2}` and set

\[
 p_t(x)=Z_t^{-1}\exp\left(c_t\cdot x-{1\over2}x^TQ_tx\right)p_0(x), \tag{1.1}
\]

where

\[
 dc_t=C_t\,dW_t+\Gamma_tm_t\,dt,\qquad
 dQ_t=\Gamma_t\,dt,\qquad m_t=\mathbb E_tX.           \tag{1.2}
\]

Stop before `(c,Q)` leaves a compact natural-parameter set, before the
moments used below exceed a fixed constant, and, when convenient, before
`g` leaves `[\delta,1-\delta]`.  Direct Ito differentiation gives

\[
 dp_t(x)=p_t(x)\langle x-m_t,C_t\,dW_t\rangle.         \tag{1.3}
\]

On such a stop all coefficients and posterior moments are bounded.

### 1.2 Why the literal formula is singular at `v=0`

Formula (0.1) has no direction-independent limit at `v=0`:

\[
 \Gamma\longrightarrow 2I-uu^T
\]

as `v` tends to zero along the line `u`.  Merely declaring an arbitrary
`u` at a zero does not prove strong existence or pathwise uniqueness.

A convenient rigorous regularization is

\[
 \boxed{
 \Gamma_\varepsilon
 =(1+a)I-{vv^T\over\ell^2+\varepsilon G},
 \qquad \varepsilon\in(0,1].}                         \tag{1.4}
\]

At `v=0`, `\Gamma_\varepsilon=2I`.  Its transverse eigenvalue is `1+a`
and its active eigenvalue is

\[
 a_\varepsilon=a+d_\varepsilon,
 \qquad d_\varepsilon={\varepsilon G\over
                              \ell^2+\varepsilon G}.  \tag{1.5}
\]

On every central, moment-bounded parameter stop this is a smooth function
of the natural parameters and is uniformly positive definite.  Its
principal square root is locally Lipschitz, so (1.2) has a unique strong
solution up to the stop.  Also

\[
 v^T\Gamma_\varepsilon v
 =G\{\beta+e_\varepsilon\},\qquad
 0\le e_\varepsilon
 ={\varepsilon\ell^2\over\ell^2+\varepsilon G}
 \le\varepsilon,                                      \tag{1.6}
\]

and

\[
 \Gamma_\varepsilon\succeq aI.                       \tag{1.7}
\]

Thus every lower estimate below which uses only `Q\succeq ZI` and
`\mathbb E\int\beta\le2H(g_0)` is uniform in `\varepsilon`.  Tightness on
the stopped compact parameter set gives a weak limit as
`\varepsilon\downarrow0`.  Equivalently, one may retain (1.4) throughout;
the extra information rate is at most `\varepsilon`, and all universal
conclusions survive after sending `\varepsilon` to zero.  This is the zero
regularization used here.

For the exact differential identities in Sections 2--5, we write the
unregularized coefficient on intervals where `v\ne0`.  The identities
which enter (0.6) hold for the regularization with the indicated favorable
error and then pass to the limit.

## 2. Complete moment equations

Put

\[
 Y=X-m_t,\qquad \xi={\bf1}_S-g_t,
\]

\[
 A=\mathbb E_tYY^T,qquad
 v=\mathbb E_t[\xi Y],\qquad
 D=\mathbb E_t[\xi YY^T].                             \tag{2.1}
\]

For a vector `z`, define

\[
 \mathcal T(z)=\mathbb E_t[YY^T\langle Y,z\rangle],
 \qquad
 \mathcal K(z)=\mathbb E_t[\xi YY^T\langle Y,z\rangle]. \tag{2.2}
\]

Posterior differentiation gives the exact stopped SDEs

\[
 \boxed{
 \begin{aligned}
 dm&=AC\,dW,\\
 dg&=v^TC\,dW,\\
 dv&=DC\,dW-A\Gamma v\,dt,\\
 dA&=\mathcal T(C\,dW)-A\Gamma A\,dt,\\
 dQ&=\Gamma\,dt.
 \end{aligned}}                                      \tag{2.3}
\]

For reference, the signed second moment satisfies

\[
\boxed{
\begin{aligned}
dD={}&\mathcal K(C\,dW)-A\langle v,C\,dW\rangle
      -v(AC\,dW)^T-(AC\,dW)v^T\\
&-\{\mathcal T(\Gamma v)+A\Gamma D+D\Gamma A\}\,dt.
\end{aligned}}                                       \tag{2.4}
\]

No derivatives of the feedback matrix occur in these equations: `C_t` is
predictable and enters only as the integrand in the likelihood SDE.

For the literal control, `\Gamma u=au`.  Let

\[
 P=I-uu^T,qquad K=D\Gamma D.
\]

On an interval on which `\ell>0`, normalization of the vector semimartingale
in (2.3) yields

\[
\boxed{
\begin{aligned}
d\ell={}&u^TDC\,dW-a\ell\,u^TAu\,dt
              +{\operatorname {tr}(PK)\over2\ell}\,dt,\\
du={}&{PDC\over\ell}\,dW-aPAu\,dt
      -{PKu\over\ell^2}\,dt
      -{\operatorname {tr}(PK)\over2\ell^2}u\,dt,\\
d[u]_t={}&{PD\Gamma DP\over\ell^2}\,dt.
\end{aligned}}                                       \tag{2.5}
\]

In particular, the angular quadratic-variation rate is

\[
\boxed{
 \omega_t={\operatorname {tr}(PD\Gamma DP)\over\ell^2}
 ={(1+a)\|PDP\|_{HS}^2+a|PDu|^2\over\ell^2}.}        \tag{2.6}
\]

The formula includes both the transverse block and the active column of
`D`.  Replacing it by `\|PDP\|_{HS}^2/\ell^2` would omit the second term.

### 2.1 The active variance has no closed scalar drift

Let `s=u^TAu`, `b=PAu`, and let `c_j=Ce_j` be the columns of `C`.  Put
`T_j=\mathcal T(c_j)` and

\[
 \tau_u=\mathbb E_t[(u^TY)^2Y].
\]

An explicit expansion of `d(u^TAu)` gives

\[
 ds=\left(C\tau_u+{2\over\ell}CDP Au\right)^TdW
       +\mathcal B_s\,dt,                              \tag{2.7}
\]

where

\[
\boxed{
\begin{aligned}
\mathcal B_s={}&-(Au)^T\Gamma(Au)-2a\,u^TAPA u
 -{2\over\ell^2}u^TAPKu
 -{s\over\ell^2}\operatorname {tr}(PK)\\
&+{1\over\ell^2}\operatorname {tr}(APK P)
 +{2\over\ell}\sum_j(PDc_j)^TT_ju.
\end{aligned}}                                       \tag{2.8}
\]

The final term is the covariation of the covariance martingale and the
direction martingale.  It has no sign.  Thus the tempting scalar assertion
that `a(1+s)\ge1` forces `s` to decrease is not a consequence of the SDE.
Controlling (2.8) uniformly would require new information about the full
third-moment tensor and the adaptive direction.

## 3. Exact information--curvature complementarity

From `\Gamma u=au`,

\[
 d\langle g\rangle_t=v^T\Gamma v\,dt
 =a\ell^2\,dt
 ={G\ell^2\over G+\ell^2}\,dt
 =G\beta\,dt.                                         \tag{3.1}
\]

The covariance Cauchy inequality gives

\[
 \ell^2\le G\,u^TAu=Gs,\qquad
 a={1\over1+k}\ge {1\over1+s}.                       \tag{3.2}
\]

For log-concave posteriors one may use the strict two-point quantization
gap `v^TA^{-1}v\le(1-\eta)G`, with numerical `\eta>0`, to improve the last
bound to `k\le(1-\eta)s`.  This does not remove the zero-seed degeneracy.

Since `H''(g)=-1/G`, Ito's formula gives (0.4), namely

\[
 dH(g)=H'(g)\sqrt{G\beta}\,dB-\frac12\beta\,dt.       \tag{3.3}
\]

On the other hand, `a+\beta=1`.  Thus, pathwise,

\[
 Z_t+\mathcal E_t=t,
 \qquad
 \mathcal E_t:=\int_0^t\beta_s\,ds.                  \tag{3.4}
\]

For every bounded stopped interval,

\[
 \mathbb E\mathcal E_T
 =2\{H(g_0)-\mathbb EH(g_T)\}\le2H(g_0).             \tag{3.5}
\]

Equations (0.5)--(0.6) follow.  A quantitative unweighted probability
consequence is: for `0<\rho<T-2H(g_0)`,

\[
 \mathbb P\{\lambda_{\min}(Q_T)\ge\rho\}
 \ge {T-2H(g_0)-\rho\over T-\rho}.                   \tag{3.6}
\]

Indeed `0\le Z_T\le T` and `\lambda_{\min}(Q_T)\ge Z_T`.
There is no mass weight in (3.6).

### 3.1 Exact Wright--Fisher time change

Let `\tau(t)=\mathcal E_t`.  Dambis--Dubins--Schwarz applied to (3.1)
shows that, in its intrinsic time, the set mass is the neutral
Wright--Fisher diffusion

\[
 d\widehat g_\tau
 =\sqrt{\widehat g_\tau(1-\widehat g_\tau)}\,dB_\tau,
 \qquad g_t=\widehat g_{\tau(t)},                     \tag{3.7}
\]

up to the usual endpoint convention.  Simultaneously,

\[
 Z_t=t-\tau(t).                                       \tag{3.8}
\]

This is an exact adversarial formulation of the scalar issue.  If one
forgets the spatial origin of `\beta`, the admissible scalar system is

\[
 dg=\sqrt{G\beta}\,dB,\qquad dZ=(1-\beta)dt,
 \qquad 0\le\beta\le1.                               \tag{3.9}
\]

It has no positive mass-weighted curvature lower bound.  To see this while
even respecting the isotropic point constraint `\beta_0\le1/2`, fix
`\varepsilon,\delta>0` and an integer `N`.  Use `\beta=1/2` on
`[0,\varepsilon]`, then `\beta=1-1/N` until the diffusion exits
`(\delta,1-\delta)` or time `T`, and set `\beta=0` after exit.  At time `T`,

\[
 \min(g_T,1-g_T)\sqrt{Z_T}
 \le {1\over2}\sqrt{\varepsilon+T/N}+\delta\sqrt T.  \tag{3.10}
\]

First send `N` to infinity, then `\varepsilon` and `\delta` to zero.  Smooth
time and state interpolations give the same conclusion.  Thus entropy,
the pointwise initial inequality `a_0\ge1/2`, and continuity alone cannot
prove the needed mass-weighted estimate.  The coefficient path (3.10) is
not asserted to arise from a log-concave localization trajectory; its
precise role is to identify the additional geometric input required.

## 4. Endpoint transfer and the exact profile Lyapunov

For a compact stopped construction, (1.3) makes the likelihood at every
fixed `x` a true martingale.  Hence

\[
 \mathbb E\mu_T(B)=\mu_0(B)                           \tag{4.1}
\]

for every Borel set `B`.  Fatou applied to exterior enlargements gives the
correctly oriented perimeter transfer

\[
 \mathbb E\mu_T^+(S)\le\mu_0^+(S).                   \tag{4.2}
\]

If `q_T=\lambda_{\min}(Q_T)>0`, then the endpoint posterior is
`q_T`-strongly log-concave.  Gaussian isoperimetry for strongly
log-concave measures gives

\[
 \mu_T^+(S)\ge\sqrt{q_T}\,\mathcal I(g_T).            \tag{4.3}
\]

Thus (0.7) is exactly the quantity produced by endpoint transfer.

There is a matching centroid inequality.  At every time with `q_t>0`,

\[
 \ell_t\le{\mathcal I(g_t)\over\sqrt{q_t}}.           \tag{4.4}
\]

This is the sharp Gaussian centroid inequality for a `q_t`-strongly
log-concave measure, obtained by scaling the unit-curvature statement.

The matrix `Q_t` is absolutely continuous.  Its least eigenvalue `q_t` is
also absolutely continuous, and for almost every `t`, directional
perturbation theory gives

\[
 \dot q_t\ge\lambda_{\min}(\Gamma_t)=a_t.             \tag{4.5}
\]

Using `\mathcal I''(g)=-1/\mathcal I(g)`, (3.1), (4.4), and (4.5), Ito's
formula gives, on a stop where `q\ge q_0>0`,

\[
\begin{aligned}
d\{\mathcal I(g)\sqrt q\}_{FV}
&\ge {a\mathcal I(g)\over2\sqrt q}\,dt
 -{G\beta\sqrt q\over2\mathcal I(g)}\,dt\\
&={a\{\mathcal I(g)^2-q\ell^2\}
       \over2\mathcal I(g)\sqrt q}\,dt\ge0.         \tag{4.6}
\end{aligned}
\]

Here `\beta=ak` and `Gk=\ell^2`.  This proves (0.8).  More generally,
`\mathcal I(g)q^p` has nonnegative drift for every `p\ge1/2`; the exponent
`1/2` is the endpoint scale and the sharp one.

The calculation cannot be started with a positive value at `q=0`:

\[
 \mathcal I(g_0)\sqrt{q_0}=0.                         \tag{4.7}
\]

If `\tau_\rho=\inf\{t:q_t\ge\rho\}` occurs before
`g` leaves `[\delta,1-\delta]`, (4.6) preserves at least
`\mathcal I(\delta)\sqrt\rho` thereafter.  Conversely, proving a universal
positive probability for this seed event and using (4.2)--(4.3) is already
a dimension-free Cheeger theorem.  Hence the seed is load-bearing, not a
routine continuity lemma.

For comparison, `w(g)=\min(g,1-g)` is less well adapted.  Away from
`g=1/2`, the finite-variation part of `w(g)\sqrt q` coming from `q` is
nonnegative, but Tanaka's formula adds

\[
 -\sqrt q\,dL_t^{1/2}(g),                             \tag{4.8}
\]

where `L^{1/2}(g)` is local time.  The smooth Gaussian profile cancels this
median-crossing loss through the sharp centroid inequality.

## 5. What small endpoint curvature forces geometrically

The literal extra transverse exposure has an exact spectral description:

\[
 Q_T=(T+Z_T)I-M_T,qquad
 M_T=\int_0^Tu_tu_t^T\,dt,qquad \operatorname {tr}M_T=T. \tag{5.1}
\]

Let `q_1\le q_2\le\cdots` be the eigenvalues of `Q_T`, and let `\theta`
be a unit eigenvector for `q_1`.  Then

\[
 \boxed{
 q_1=Z_T+\int_0^T|P_{\theta^\perp}u_t|^2\,dt.}        \tag{5.2}
\]

Also, if `\lambda_1\ge\lambda_2` are the two largest eigenvalues of `M_T`,
then `\lambda_1+\lambda_2\le T`.  Hence

\[
 \boxed{q_2\ge {T+2Z_T\over2}.}                      \tag{5.3}
\]

Thus `q_1\le\rho\ll T` forces all three properties

\[
 Z_T\le\rho,\qquad
 \int_0^T|P_{\theta^\perp}u_t|^2dt\le\rho,
 \qquad q_2\ge T/2.                                  \tag{5.4}
\]

The active direction is locked in projective `L^2` to one line, while the
posterior has order-`T` curvature in every orthogonal endpoint direction.

There is also a precise high-signal occupation statement.  Since
`a=1/(1+k)`, for every `K>0`,

\[
 \left|\{t\le T:k_t\le K\}\right|
 \le(1+K)Z_T.                                         \tag{5.5}
\]

Consequently a low-curvature path spends almost all its time with both a
locked direction and large normalized binary correlation `k=\ell^2/G`.
Covariance Cauchy says its active variance is at least `k`.

At an intermediate time for which `q_t>0`, (4.4) gives

\[
 k_tq_t\le{\mathcal I(g_t)^2\over G_t}.               \tag{5.6}
\]

On a central mass band the right side is numerical.  Thus the critical bad
profile is exactly

\[
 k_t\asymp q_t^{-1},                                  \tag{5.7}
\]

created in an arbitrarily short initial layer and maintained along one
line.  Brascamp--Lieb or (5.6) only yields

\[
 \dot Z={1\over1+k}\gtrsim Z,                         \tag{5.8}
\]

which amplifies an existing seed exponentially but has the zero solution.
Over a fixed information budget it amplifies a dimension-dependent seed by
only a fixed factor.

### 5.1 Angular stability does not by itself exclude the bad phase

The exact angular clock is (2.6).  The available general angular-stability
theorem for a `q`-strongly log-concave central posterior states, with

\[
 \epsilon=1- {\sqrt q\,\ell\over\mathcal I(g)},
\]

that

\[
 \|PD\|_{HS}^2\le C_\delta q^{-2}\Omega_\delta(\epsilon). \tag{5.9}
\]

Since `\|\Gamma\|_{op}\le2`, this gives, near centroid saturation,

\[
 \omega_t\le
 {2\|PD\|_{HS}^2\over\ell^2}
 \le {C_\delta\over q_t}\Omega_\delta(\epsilon_t).   \tag{5.10}
\]

This is an **upper** bound on angular quadratic variation and is compatible
with a locked, nearly one-dimensional phase.  If the centroid defect is not
small, (5.9) supplies no lower angular motion.  Moreover the drift in
`du`, including the `PKu/\ell^2` term, can adaptively select a direction, so
large quadratic variation alone cannot be converted into (5.2) without a
drift estimate.  The general angular theorem therefore narrows the missing
case to a coherent high-signal factor but does not exclude it.

## 6. Determinant and covariance Lyapunovs

Assume `A\succ0` on a bounded stop and let `c_j=Ce_j`.  From (2.3),

\[
\boxed{
\begin{aligned}
d\log\det A={}&
 \sum_j\operatorname {tr}(A^{-1}\mathcal T(c_j))\,dW_j
 -\operatorname {tr}(\Gamma A)\,dt\\
&-\frac12\sum_j
 \operatorname {tr}\{A^{-1}\mathcal T(c_j)
                       A^{-1}\mathcal T(c_j)\}\,dt.
\end{aligned}}                                       \tag{6.1}
\]

The last term is nonpositive.  If `s=u^TAu`, covariance Cauchy implies

\[
 \operatorname {tr}(\Gamma A)
 \ge a s\ge{s\over1+s}.                               \tag{6.2}
\]

Thus a large active variance spends determinant at a numerical rate.  This
does not bound its occupation: transverse exposure can drive
`\log\det A` to `-\infty`, and the determinant has no dimension-free lower
barrier.  In particular it can pay indefinitely in the already strongly
contracted `n-1` directions while one selected variance remains large.

Likewise, `d\log\det(Q+\varepsilon I)
=\operatorname {tr}((Q+\varepsilon I)^{-1}\Gamma)dt` is exact but controls
a product of eigenvalues, not the exceptional eigenvalue.  Formula (2.8)
shows why the active variance itself has no signed scalar drift.  These
three natural Lyapunovs do not improve (5.8).

## 7. Explicit stress tests

### 7.1 Gaussian halfspace: the control succeeds with room to spare

Let `\mu_0=N(0,I)` and `S=\{x_1\ge0\}`.  Product structure keeps
`u=e_1`.  If `q` is the accumulated curvature in this line, its posterior
standard deviation is `\sigma=(1+q)^{-1/2}` and

\[
 \ell=\sigma\mathcal I(g),\qquad
 k={\sigma^2\mathcal I(g)^2\over G}.                  \tag{7.1}
\]

The elementary Gaussian bound
`\mathcal I(g)^2/G\le2/\pi` gives

\[
 a\ge{1\over1+2/\pi}={\pi\over\pi+2}.                \tag{7.2}
\]

Thus the exceptional curvature grows linearly and the direction never
rotates.  This model does not exhibit the seed obstruction.

### 7.2 A radial Gaussian cut: no high-dimensional initial collapse

Let `G_n` be standard Gaussian and let `S` be a centered ball of half
Gaussian mass.  Initially `v=0` and rotational symmetry gives `D=dI`.
More generally, for every Gaussian set,

\[
 \|D\|_{HS}^2\le2g(1-g).                              \tag{7.3}
\]

Indeed `M\mapsto\operatorname {Cov}({\bf1}_S,
G_n^TMG_n-\operatorname {tr}M)` has representing matrix `D`, while the
Gaussian second-chaos identity gives variance `2\|M\|_{HS}^2`; Cauchy--
Schwarz and duality prove (7.3).  Under the regularization
`\Gamma_\varepsilon=2I` at `v=0`, the initial quadratic-variation rate of
`v` is at most `4g(1-g)`, independent of `n`.  The radial model therefore
does not make `a` collapse on a dimension-dependent infinitesimal scale.

### 7.3 Product exponentials

Let `Z_i` be independent rate-one exponentials and set `X_i=Z_i-1`.
For the coordinate cut `S=\{Z_1\ge\log2\}`,

\[
 g={1\over2},\qquad
 v={\log2\over2}e_1,
 \qquad k=(\log2)^2,
 \qquad a={1\over1+(\log2)^2}.                        \tag{7.4}
\]

The posterior remains a product, `u=e_1` for the whole stopped evolution,
and the extra transverse rate `aP` cannot contribute to its exceptional
curvature.  Starting isotropically, however, (7.4) supplies a universal
seed.

For the balanced maximum cut

\[
 S_n=\{\max_i Z_i\ge L_n\},\qquad
 (1-e^{-L_n})^n={1\over2},                            \tag{7.5}
\]

exchangeability gives

\[
 k_0={nL_n^2e^{-2L_n}\over(1-e^{-L_n})^2}
 \sim{(\log2)^2(\log n)^2\over n}\longrightarrow0.   \tag{7.6}
\]

Thus `a_0\to1`; the controller initially exposes every direction at a
universal rate.  The maximum cut has a fast angular winner mechanism, but
the literal extra transverse factor is between one and two and changes no
critical dimension scale.  Existing one-coordinate Laplace estimates show
that full exposure at a fixed `1/\log n` scale does not give one coordinate
order-one posterior tail mass.  Hence this is a stress test for the
initial-layer coherence theorem, not a proved counterexample to the driver.

### 7.4 Isotropic simplex: an order-one tilt creates a small-`a` state

For the uniform law on the standard `n`-simplex, let `R=nX_1`.  Its density
is

\[
 f_n(r)=(1-r/n)^{n-1}{\bf1}_{[0,n]}(r),               \tag{7.7}
\]

and `\operatorname {Var}R\to1`; this is a vertex direction in isotropic
units up to a universal factor.  Under the order-one linear tilt `e^R`,
put `r=\sqrt n\,y`.  Uniformly on compact `y` intervals,

\[
 r+(n-1)\log(1-r/n)=-{y^2\over2}+o(1).                \tag{7.8}
\]

Thus `R/\sqrt n` under the tilted law converges to a half-normal random
variable and its variance is `\Theta(n)`.  A median halfspace for this
tilted one-dimensional marginal has

\[
 k=\Theta(n),\qquad a=\Theta(n^{-1}).                 \tag{7.9}
\]

This proves that neither log-concavity, balance at the current posterior,
nor an order-one natural parameter prevents a small instantaneous active
rate.  The median halfspace in (7.9) has very small mass under the original
simplex, so (7.9) is not an isotropic half-set trajectory obstruction.
It is a valid test against state-local determinant or covariance claims.

### 7.5 Rare exponential tilts, with the same set balanced at both ends

The statewise issue can be made sharper.  Let `\pi_\lambda` be the
rate-`\lambda` exponential law on `[0,\infty)`, where `0<\lambda<1`.
There are unique numbers `0<A_\lambda<B_\lambda` such that

\[
 S_\lambda=[0,A_\lambda]\cup[B_\lambda,\infty)
\]

has probability `1/2` under both `\pi_1` and `\pi_\lambda`.  To see this,
put `x=e^{-A}`, `y=e^{-B}`.  The two equations are

\[
 x-y={1\over2},\qquad x^\lambda-y^\lambda={1\over2}. \tag{7.10}
\]

For `x\in(1/2,1)`, the second difference decreases continuously from
`2^{-\lambda}` to `1-2^{-\lambda}`, so there is exactly one solution.
As `\lambda\downarrow0`,

\[
 A_\lambda\to\log2,qquad
 B_\lambda\sim{\log2\over\lambda}.                  \tag{7.11}
\]

At the rare tilted state `\pi_\lambda`, direct integration gives

\[
 \left|\operatorname {Cov}_{\pi_\lambda}
       ({\bf1}_{S_\lambda},X)\right|
 ={\log2+o(1)\over2\lambda}.                          \tag{7.12}
\]

Therefore

\[
 k_\lambda={ (\log2)^2+o(1)\over\lambda^2},
 \qquad a_\lambda\asymp\lambda^2.                   \tag{7.13}
\]

The original law `\pi_1`, after translation, is isotropic in one dimension,
and the **same** set has half mass there.  The tilted law is obtained from
it by the linear natural parameter `c=1-\lambda`.  However, the controlled
process also accumulates a quadratic parameter `Q`, so (7.10)--(7.13) do
not show that this state is reached with appreciable probability along the
soft trajectory.  Moreover the original exterior perimeter is

\[
 e^{-A_\lambda}+e^{-B_\lambda}\longrightarrow{1\over2}, \tag{7.14}
\]

so the set is not a small-boundary counterexample.  The example precisely
rules out any argument based only on current balance, current
log-concavity, the original balance of the same set, and the size of the
linear tilt.

## 8. Approximation and scope

All stochastic identities above are first statements on bounded stopping
times for smooth compactly supported posteriors.  For an arbitrary
log-concave law, restrict to its minimal affine hull, truncate by increasing
compact convex sets, smooth inside that hull, recenter, and whiten.  On a
fixed natural-parameter and moment stop, the posterior moment maps converge
locally uniformly.  The regularized controls (1.4) are bounded by `2I`, so
the stopped solution laws are tight and every displayed martingale identity
passes to a weak limit.  The entropy estimate uses bounded `H` and Fatou;
the perimeter transfer already has the correct lower-semicontinuous
orientation.  One then removes parameter stops, sends
`\varepsilon\downarrow0`, and finally removes truncation.  A quantitative
seed theorem, if proved uniformly before these limits, would therefore
transfer with the same universal constant.

What has been proved is:

1. exact stopped SDEs for the driver and a nonsingular approximation at
   `v=0`;
2. the pathwise information--curvature complementarity and the universal
   unweighted expectation (0.6);
3. the exact post-seed Lyapunov (0.8);
4. the locked high-signal structure (5.4)--(5.7) of every low-curvature
   path;
5. explicit statewise obstructions to a scalar active-variance or
   determinant proof.

The missing assertion is a uniform lower bound, for one numerical
`\rho,\delta,p>0`, on

\[
 \mathbb P\{\lambda_{\min}(Q_t)\hbox{ reaches }\rho
       \hbox{ before }g_t\notin[\delta,1-\delta]\}\ge p. \tag{8.1}
\]

By (4.2)--(4.3), (8.1) already implies the half-mass KLS inequality.  The
scalar Wright--Fisher reduction permits its failure, while the spatial
equations show that failure must be an initially selected, projectively
locked, high-variance one-dimensional phase.  Neither instantaneous rates,
the determinant, Brascamp--Lieb, quantization, nor the currently available
angular-stability estimate excludes that phase.
