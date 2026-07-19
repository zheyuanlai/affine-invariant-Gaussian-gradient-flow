# Aggressive mass-preserving stochastic localization

## 0. Verdict

Let \(\mu_t\) be a stochastic-localization posterior, let

\[
 A_t=\operatorname {Cov}_{\mu_t}X,\qquad
 v_t=\operatorname {Cov}_{\mu_t}({\bf 1}_S,X),
 \qquad q_t=A_t^{-1}v_t,
\]

and, while \(A_t\succ0\) and \(q_t\ne0\), put

\[
 u_t={q_t\over |q_t|},\qquad
 P_t=I-u_tu_t^T,\qquad
 D_t=A_t^{-1}P_tA_t^{-1}.                       \tag{0.1}
\]

Here \(D_t\) is the covariance of the Brownian coefficient in the
density equation; a square root \(C_t\) satisfies \(C_tC_t^T=D_t\).
The attraction of (0.1) is genuine:

\[
 D_tv_t=0,\qquad A_tD_tA_t=P_t.                 \tag{0.2}
\]

Consequently the chosen set mass is pathwise constant and the covariance
drift is exactly \(-P_t\), rather than a covariance-weighted projection.
All Itô identities requested in the problem are derived below.

The terminal conclusion, however, does not follow.

1.  Up to a legitimate bounded stopping time,

    \[
      d\,\operatorname {tr}A_t
      =dM_t-(n-1)\,dt.                            \tag{0.3}
    \]

    Thus an isotropic process reaches the first singular face in expected
    time at most \(n/(n-1)\), provided no earlier coefficient stop occurs.
    A first singular face normally has rank \(n-1\), not rank one.

2.  The feedback direction \(u_t=A_t^{-1}v_t/|A_t^{-1}v_t|\) can rotate.
    Even for a Gaussian starting law and a smooth fixed set, this rotation
    is real.  For an off-centre Euclidean ball its initial angular
    quadratic-variation rate is exactly

    \[
                         {n-1\over d^2},          \tag{0.4}
    \]

    where \(d\) is the distance between the Gaussian mean and the centre
    of the ball.  For a Gaussian posterior

    \[
      A_t=(1-t)I+\int_0^t u_su_s^T\,ds.           \tag{0.5}
    \]

    If \(u_s\) is not constant, then \(A_1\) need not have rank one.  The
    covariance drift protects a moving line, and the time average of
    moving lines can have any rank.

3.  At a transverse eigenvalue \(\lambda\), the driver has size
    \(D\asymp\lambda^{-2}\).  In the exact Gaussian-halfspace solution,

    \[
      \int_0^{1-\varepsilon}\operatorname {tr}(D_tA_t)\,dt
        =(n-1)\log {1\over\varepsilon},\qquad
      \int_0^{1-\varepsilon}\operatorname {tr}D_t\,dt
        =(n-1)(\varepsilon^{-1}-1).                \tag{0.6}
    \]

    Thus the posterior has a perfectly good weak rank-one limit while
    the natural-parameter stochastic clock diverges.  The stopped SDE
    cannot simply be evaluated at the terminal time.

4.  At \(v=0\) the hard direction selector is undefined.  A smooth
    relaxation gives approximate, not exact, mass conservation.  Under
    an additional stopped tightness-and-identification argument, a
    vanishing-mesh subsequence would produce a relaxed predictable
    projection, but that compactness argument is only outlined here.
    Such a limit is noncanonical and still stops at the first covariance
    singularity.  The one-dimensional Gaussian symmetric-interval
    example rules out the naive full-driver restart as a continuous weak
    SDE.

5.  Neither (0.3) nor the determinant identity controls the last
    variance.  Product exponential posteriors contain protected
    one-dimensional factors of arbitrarily large variance.  The
    product-maximum test has initial angular quadratic variation
    \(\Theta((\log n)^2)\), so the aggressive normalization gives no
    leading-order suppression of the winner-selection mechanism.

The positive conclusion is therefore limited but exact: (0.1) gives a
mass-preserving, constant-trace-speed localization up to bounded stops.
It proves universal expected time to the first covariance boundary.  It
does not prove universal expected time to a rank-one posterior, nor a
universal bound on the surviving variance.  Any use of it for KLS needs a
new singular-continuation theorem and a new protected-variance estimate;
neither can be replaced by a terminal-needle assumption.

## 1. Stopped construction

Let \(p_0\) be a positive smooth log-concave probability density on
\(\mathbb R^n\), initially with all moments required below.  Introduce
natural parameters

\[
 p_t(x)={1\over Z_t}
 \exp\left\{c_t\cdot x-\frac12x^TQ_tx\right\}p_0(x).
                                                               \tag{1.1}
\]

For a predictable positive-semidefinite \(D_t\), choose a predictable
square root \(C_t\) with

\[
                         C_tC_t^T=D_t,              \tag{1.2}
\]

and set

\[
 dc_t=C_t\,dW_t+D_ta_t\,dt,\qquad
 dQ_t=D_t\,dt,\qquad
 a_t=\mathbb E_tX.                                   \tag{1.3}
\]

Fix \(\varepsilon,\delta,R>0\).  Initially stop before

\[
 \lambda_{\min}(A_t)\le\varepsilon,\qquad
 |q_t|\le\delta,\qquad
 \|A_t\|_{\rm op}+|a_t|+\|c_t\|+\|Q_t\|>R,            \tag{1.4}
\]

and before the moments appearing below exceed \(R\).  On this stopped
region the posterior moments are smooth functions of \((c,Q)\), the
matrix in (0.1) has constant rank, and its principal square root is
locally Lipschitz.  Standard finite-dimensional SDE theory therefore
gives a unique strong solution up to the stop.  Every identity in
Sections 1--4 is first asserted in this setting.

The moment and natural-parameter stops matter.  In particular,
\(\lambda_{\min}(A)\ge\varepsilon\) bounds \(D\) by
\(\varepsilon^{-2}I\), but it does not by itself bound all third moments
uniformly over an unbounded posterior family.  Log-concavity supplies
moment bounds after covariance and barycentre control; the explicit
stops avoid using such bounds circularly.

## 2. Density equation, including the normalization correction

Let

\[
 r_t(x)=
 \exp\left\{c_t\cdot x-\frac12x^TQ_tx\right\}p_0(x).
\]

The quadratic variation of \(c_t\cdot x\) is
\(x^TD_tx\,dt\), which cancels the Itô correction from
\(-\frac12x^TQ_tx\).  From (1.3),

\[
 {dr_t(x)\over r_t(x)}
   =x^TC_t\,dW_t+x^TD_ta_t\,dt.                      \tag{2.1}
\]

Differentiating \(Z_t=\int r_t\), and then applying the quotient rule
with its cross-variation term, gives the drift-free normalized equation

\[
 \boxed{\quad
 dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle .
 \quad}                                                \tag{2.2}
\]

Equivalently, at each fixed \(x\),

\[
 d\log p_t(x)
 =\langle x-a_t,C_t\,dW_t\rangle
 -{1\over2}(x-a_t)^TD_t(x-a_t)\,dt.                   \tag{2.3}
\]

Formula (2.3) is the full pointwise Itô correction.  There is no hidden
finite-variation term in (2.2).

For every fixed integrable scalar or vector function \(F\),

\[
 d\,\mathbb E_tF(X)
 =\operatorname {Cov}_t(F(X),X)^TC_t\,dW_t,           \tag{2.4}
\]

with the evident matrix interpretation for vector-valued \(F\).

## 3. Exact moment SDEs

Write

\[
 Y=X-a_t,\qquad h={\bf1}_S-g_t,\qquad
 g_t=\mathbb E_t{\bf1}_S,
\]

\[
 A_t=\mathbb E_tYY^T,\qquad
 v_t=\mathbb E_t[hY],\qquad
 H_t=\mathbb E_t[hYY^T],                              \tag{3.1}
\]

and define the centred third-moment map

\[
 {\cal T}_t(z)
   =\mathbb E_t[YY^T\langle Y,z\rangle].               \tag{3.2}
\]

If \(c_{t,k}=C_te_k\), abbreviate

\[
                         T_{t,k}={\cal T}_t(c_{t,k}).
\]

Applying (2.4) and retaining every product covariation yields

\[
 \boxed{
 \begin{aligned}
 da_t&=A_tC_t\,dW_t,\\
 dg_t&=v_t^TC_t\,dW_t,\\
 dv_t&=H_tC_t\,dW_t-A_tD_tv_t\,dt,\\
 dA_t&=\sum_kT_{t,k}\,dW_{t,k}-A_tD_tA_t\,dt.
 \end{aligned}}                                       \tag{3.3}
\]

Here is the drift check in the third line.  Put
\(m_t=\mathbb E_t[{\bf1}_SX]\), so \(v=m-ga\).  The martingale
coefficient of \(m\) is

\[
 \mathbb E_t[{\bf1}_SXY^T]C
 =\{a v^T+\mathbb E_t[{\bf1}_SYY^T]\}C.
\]

Subtracting the coefficients of \(g\,a\) leaves \(HC\), while

\[
                         d[g,a]_t=A_tD_tv_t\,dt.       \tag{3.4}
\]

This is the finite-variation term in \(dv\).  It vanishes under (0.1),
but it must not be omitted for a softened or frozen-mesh controller.

The covariance drift in (3.3) includes the term
\(-d[a,a]_t=-A_tD_tA_t\,dt\).  No fourth-moment drift is missing.

## 4. The aggressive hard controller

Assume \(A\succ0\) and \(q=A^{-1}v\ne0\), and define \(u,P,D\) by
(0.1).  Since \(v=Aq\),

\[
 Dv=A^{-1}PA^{-1}Aq=A^{-1}Pq=0.                       \tag{4.1}
\]

Because \(D=C C^T\),

\[
                         |C^Tv|^2=v^TDv=0,
\]

and hence \(C^Tv=0\).  It follows from (3.3) that

\[
 \boxed{
 g_t=g_0,\qquad
 dv_t=H_tC_t\,dW_t,\qquad
 dA_t=\sum_kT_{t,k}\,dW_{t,k}-P_t\,dt .
 }                                                      \tag{4.2}
\]

The last equality uses

\[
 A D A=A A^{-1}P A^{-1}A=P.                            \tag{4.3}
\]

This is the exact advantage over the covariance-normalized driver
\(A^{-1/2}PA^{-1/2}\): the ordinary Euclidean covariance loses one unit
per unit time in every instantaneous direction orthogonal to \(q\).

The mean equation also has an especially simple quadratic variation:

\[
 da=AC\,dW,\qquad
                         d[a]_t=ADA\,dt=P_t\,dt.        \tag{4.4}
\]

Thus the mean diffuses with unit covariance in the same \(n-1\)
directions in which the covariance has drift \(-I\).

### 4.1 The feedback-vector SDE

The controller depends on \(q=A^{-1}v\), not merely on \(v\).  Its Itô
equation contains terms which are invisible in (4.2).  Let \(B=A^{-1}\)
and

\[
                         h_k=Hc_k.
\]

Matrix Itô calculus gives

\[
 dB=-\sum_kBT_kB\,dW_k+BPB\,dt
       +\sum_kBT_kBT_kB\,dt.                           \tag{4.5}
\]

Multiplying \(q=Bv\), including \(d[B,v]\), and using \(Pq=0\), gives

\[
 \boxed{
 \begin{aligned}
 dq={}&\sum_k B(h_k-T_kq)\,dW_k\\
 &+\sum_kBT_kB(T_kq-h_k)\,dt.
 \end{aligned}}                                       \tag{4.6}
\]

In particular, if

\[
 K_k=B(h_k-T_kq),                                      \tag{4.7}
\]

then the angular quadratic variation of \(u=q/|q|\) is

\[
 d[u]_t={1\over |q|^2}
 P\left(\sum_kK_kK_k^T\right)P\,dt.                    \tag{4.8}
\]

Both singularities are explicit: \(B=A^{-1}\) blows up at a covariance
face, and the normalization blows up at \(q=0\).  A calculation based
only on \(dv=HC\,dW\) misses the covariance-noise correction
\(-BT_kq\) in (4.6).

## 5. Trace, determinant, and perimeter

Define

\[
                         m_{3,t}=\mathbb E_t[|Y|^2Y].
\]

Taking traces in (4.2) yields

\[
 \boxed{
 d\,\operatorname {tr}A_t
   =m_{3,t}^TC_t\,dW_t-(n-1)\,dt.
 }                                                      \tag{5.1}
\]

Consequently, for every bounded stopping time \(\tau\) before the
coefficient stops,

\[
 \mathbb E\operatorname {tr}A_\tau
 =\operatorname {tr}A_0-(n-1)\mathbb E\tau.             \tag{5.2}
\]

This equality requires the stopped stochastic integral to be a true
martingale; that is why the moment stop was included in (1.4).

For the determinant, matrix Itô calculus gives the cleanest formula for
the logarithm:

\[
 \boxed{
 \begin{aligned}
 d\log\det A_t
 ={}&\sum_k\operatorname {tr}(A_t^{-1}T_{t,k})\,dW_{t,k}\\
 &-\operatorname {tr}(A_t^{-1}P_t)\,dt\\
 &-\frac12\sum_k
 \operatorname {tr}
 \left(A_t^{-1}T_{t,k}A_t^{-1}T_{t,k}\right)dt .
 \end{aligned}}                                       \tag{5.3}
\]

The last line is nonpositive because
\(A^{-1/2}T_kA^{-1/2}\) is symmetric.  Notice that the deterministic
drift is not \(-(n-1)\); it is the sum of inverse variances in the
instantaneous transverse space.

For completeness, the determinant itself satisfies

\[
 \boxed{\begin{aligned}
 {d\det A\over\det A}
 ={}&\sum_k\operatorname {tr}(A^{-1}T_k)\,dW_k
 -\operatorname {tr}(A^{-1}P)\,dt\\
 &+{1\over2}\sum_k\left[
       \{\operatorname {tr}(A^{-1}T_k)\}^2
       -\operatorname {tr}(A^{-1}T_kA^{-1}T_k)
                    \right]dt.
 \end{aligned}}                                       \tag{5.4}
\]

Unlike the correction in (5.3), the last line of (5.4) has no fixed
sign.

Now suppose \(S\) has compact \(C^1\) boundary and define its weighted
perimeter

\[
 \Pi_t(S)=\int_{\partial S}p_t(x)\,
                         d{\cal H}^{n-1}(x).            \tag{5.5}
\]

Let

\[
 b_{\partial,t}
 ={1\over\Pi_t(S)}
   \int_{\partial S}x\,p_t(x)\,d{\cal H}^{n-1}(x)
\]

when \(\Pi_t(S)>0\).  Integrating (2.2) over the fixed boundary gives

\[
 \boxed{
 d\Pi_t(S)=\Pi_t(S)
 \langle b_{\partial,t}-a_t,C_t\,dW_t\rangle .
 }                                                      \tag{5.6}
\]

Thus

\[
 \boxed{
 d\log\Pi_t(S)
 =\langle b_{\partial,t}-a_t,C_t\,dW_t\rangle
 -{1\over2}(b_{\partial,t}-a_t)^TD_t
                 (b_{\partial,t}-a_t)\,dt .
 }                                                      \tag{5.7}
\]

Stopped perimeter is a martingale.  Unstopped perimeter is a nonnegative
local martingale and therefore a supermartingale.  Since \(g_t=g_0\), the
same statement holds for the stopped isoperimetric ratio
\(\Pi_t(S)/\min(g_0,1-g_0)\).

For a finite-perimeter set, (5.6) first holds on the reduced boundary
after smooth approximation on a bounded coefficient stop.  Passing to an
unbounded or singular terminal time additionally requires uniform
integrability or a lower-semicontinuity argument; (5.6) alone does not
supply either one.

## 6. The zero-signal state and mesh relaxation

The formula \(P=I-qq^T/|q|^2\) is not defined at \(q=0\).  This is not a
removable notational defect.

### 6.1 A smooth approximate controller

For \(\eta>0\), define

\[
 P^{(\eta)}
   =I-{qq^T\over |q|^2+\eta^2},\qquad
 D^{(\eta)}=A^{-1}P^{(\eta)}A^{-1}.                    \tag{6.1}
\]

On \(A\succeq\varepsilon I\), this is a continuous positive-definite
controller, including at \(q=0\).  Its covariance drift is
\(-P^{(\eta)}\), but it only approximately preserves mass:

\[
 \boxed{\begin{aligned}
 d[g]_t
 &=v^TD^{(\eta)}v\,dt\\
 &=q^TP^{(\eta)}q\,dt
 ={\,|q|^2\eta^2\over |q|^2+\eta^2}\,dt
 \le\eta^2dt.
 \end{aligned}}                                       \tag{6.2}
\]

For every deterministic \(T\), the stopped process therefore satisfies

\[
 \mathbb E\sup_{t\le T}|g_t-g_0|^2\le4\eta^2T.          \tag{6.3}
\]

The equations for \(v\) and the trace now contain

\[
 -AD^{(\eta)}v
 =-P^{(\eta)}q
 =-{\eta^2q\over |q|^2+\eta^2},                        \tag{6.4}
\]

and

\[
 d\,\operatorname {tr}A
 =dM-\left(n-{|q|^2\over |q|^2+\eta^2}\right)dt.        \tag{6.5}
\]

Thus the softened process has a classical stopped SDE and its mass error
vanishes in \(L^2\) as \(\eta\downarrow0\).  This is the simplest precise
meaning of mesh relaxation.

### 6.2 Frozen meshes and relaxed limits

An exactly projection-valued alternative is to choose a partition
\(0=t_0<t_1<\cdots\), freeze on
\([t_j,t_{j+1})\) a rank-\((n-1)\) projection annihilating \(q_{t_j}\),
and use

\[
 D_t=A_t^{-1}P_{t_j}A_t^{-1}.
\]

On a bounded coefficient stop, the increment
\(q_t-q_{t_j}\) is \(O_{L^2}(\sqrt{t-t_j})\).
Consequently the total set-mass quadratic variation caused by the
lagging projection tends to zero with the mesh.  The bounded stopped
coefficients give the expected tightness estimates.  If one additionally
proves joint convergence of the posterior moments and identifies the
limiting martingale problem, any resulting relaxed matrix \(\bar P_t\)
satisfies

\[
 0\preceq\bar P_t\preceq I,\qquad
 \bar P_tq_t=0                                         \tag{6.6}
\]

for almost every \(t\), and the limiting mass is constant.  Away from
\(q=0\), necessarily
\(\bar P=I-qq^T/|q|^2\).  At \(q=0\), different mesh selectors can give
different relaxed limits.  A limit of (6.1) can also have trace between
\(n-1\) and \(n\) on its zero-signal occupation set.

This paragraph is a compactness scheme, not a complete weak-convergence
theorem: the topology, joint convergence, and martingale-problem
identification would have to be written out for a final construction.  It
does not prove uniqueness of the relaxed law.  It also does not pass
through \(A=0\): the coefficient \(A^{-1}\bar P A^{-1}\) is still
singular at the first covariance face.

### 6.3 The naive full restart can fail to exist

The obstruction is already one-dimensional.  Let
\(\mu_0=N(0,1)\) and \(S=[-r,r]\).  A Gaussian posterior is
\(N(a,A)\), and

\[
                         v=A\,\partial_a\mu_{a,A}(S).
\]

Locally, \(v=0\) exactly at \(a=0\).  In dimension one the hard
projection is \(P=0\) whenever \(v\ne0\).  If one declares \(P=1\) at
\(v=0\), the natural parameter would have to satisfy

\[
                         d[c]_t=A_t^{-2}
                         {\bf1}_{\{c_t=0\}}\,dt.        \tag{6.7}
\]

The occupation-density formula for a continuous semimartingale gives

\[
 \int_0^T{\bf1}_{\{c_t=0\}}\,d[c]_t=0.                 \tag{6.8}
\]

Equations (6.7)--(6.8) force \([c]_T=0\), hence \(c\equiv0\), but then
(6.7) forces strictly positive quadratic variation.  This contradiction
shows that the pointwise full-restart convention has no continuous weak
solution.  Smoothing or a genuine relaxed mesh is essential.

## 7. What the trace drift proves

Let

\[
 \tau_\varepsilon
 =\inf\{t:\lambda_{\min}(A_t)\le\varepsilon\}
\]

and retain the auxiliary coefficient stops.  If \(A_0=I\), optional
stopping in (5.1) gives, in the ideal case where no auxiliary stop occurs,

\[
 \mathbb E(\tau_\varepsilon\wedge T)
 ={n-\mathbb E\operatorname {tr}
       A_{\tau_\varepsilon\wedge T}\over n-1}
 \le {n\over n-1}.                                    \tag{7.1}
\]

If \(\tau_\varepsilon<\infty\) almost surely and the stopped martingales
are uniformly integrable, then

\[
                 \mathbb E\tau_\varepsilon
 \le {n(1-\varepsilon)\over n-1}.                      \tag{7.2}
\]

This is a dimension-free expected-time estimate for the **first**
covariance boundary.

At \(\tau_\varepsilon\), only one eigenvalue is known to equal
\(\varepsilon\).  The other \(n-1\) eigenvalues can be macroscopic.  The
trace equation contains no count of small eigenvalues.  The determinant
equation is even more biased toward the first small eigenvalue, because
\(\operatorname {tr}(A^{-1}P)\) becomes large as soon as one transverse
variance becomes small.

### 7.1 A moving protected line

When the posterior is Gaussian, every centred third moment vanishes, so
the covariance martingale in (4.2) vanishes identically.  Starting from
\(A_0=I\),

\[
 \boxed{
 A_t=I-\int_0^tP_s\,ds
 =(1-t)I+\int_0^tu_su_s^T\,ds.
 }                                                      \tag{7.3}
\]

If \(u_s\equiv u\), then \(A_1=uu^T\), exactly as desired.  If \(u_s\)
visits two noncollinear directions on sets of positive time, then
\(A_1\) has rank at least two.  If the occupation measure of \(u_s\) is
approximately isotropic, then

\[
 A_t\approx\left(1-{n-1\over n}t\right)I,              \tag{7.4}
\]

which contracts all directions together.  A relaxed control with
\(\int_0^tu_su_s^Tds=(t/n)I\) reaches the zero matrix at
\(t=n/(n-1)\) and never has a rank-one stage.

Equation (7.3) is a deterministic algebraic counterexample to the
inference

\[
 \text{drift }-P_t\quad\Longrightarrow\quad
 \text{one fixed surviving line}.
\]

A terminal-line result must control the occupation measure of \(u_t\);
the trace and determinant do not do so.

### 7.2 Gaussian off-centre balls rotate the line

The rotation in the preceding paragraph occurs for an actual fixed set.
Let the current posterior be \(N(a,I)\) and let

\[
                         S=B(b,R),\qquad d=|b-a|>0.
\]

Write

\[
                         g(a)=G(|b-a|).
\]

Gaussian differentiation gives

\[
 q=A^{-1}v=\nabla_ag,\qquad H=\nabla_a^2g,\qquad
 {\cal T}=0.                                          \tag{7.5}
\]

The tangential eigenvalue of \(H\) is \(G'(d)/d\), while
\(|q|=|G'(d)|\).  Since \(C=P\) at \(A=I\), (4.8) gives the exact initial
angular rate

\[
 \boxed{
 {d\over dt}\operatorname {tr}[u]_0={n-1\over d^2}.
 }                                                      \tag{7.6}
\]

For every \(b\ne a\), choose \(R\) so that the ball has any prescribed
mass in \((0,1)\), for example \(1/2\).  The radial derivative is
nonzero, and hence (7.6) is nondegenerate.  When \(d\) is universal, the
protected line turns on the \(1/n\) time scale, much faster than the
order-one covariance-collapse scale.

Moreover, (7.3) implies \(A_t\succeq(1-t)I\) for \(t<1\).  Thus covariance
singularity cannot prevent this initial rotation.  On every path for
which \(u_s\) explores more than one direction before time \(1\),
\(A_1\) is not rank one.  If subsequently a single eigenvalue is the
first one to vanish, the first singular state has rank \(n-1\); no claim
that this generic-looking event always occurs is needed.  In either case,
the stopped controller supplies no continuation theorem.

This example is not a Gaussian isoperimetric minimizer; halfspaces behave
better.  It nevertheless disproves any terminal-rank theorem for (0.1)
which uses only log-concavity, mass preservation, and the covariance
drift.

## 8. Noise blow-up at a covariance face

Suppose \(Ae=\lambda e\) and \(e\perp q\).  Then

\[
                         e^TDe=\lambda^{-2}.           \tag{8.1}
\]

The covariance drift in that direction has size one, but the density
driver has size \(1/\lambda\).  Its posterior-averaged pointwise
quadratic variation is

\[
 \mathbb E_t|C_t^TY|^2
 =\operatorname {tr}(D_tA_t),                         \tag{8.2}
\]

which contains \(1/\lambda\).  The natural-parameter quadratic variation
contains \(\operatorname {tr}D_t\), hence \(1/\lambda^2\).

The covariance martingale has coefficients
\({\cal T}(C e_k)\), and the feedback martingale in (4.6) contains one
additional \(A^{-1}\).  Log-concavity bounds standardized third moments,
but it does not cancel these inverse powers.  Therefore bounded physical
time does not imply bounded stochastic clock or uniform integrability at
the terminal face.

This is not merely a pessimistic estimate.  The next Gaussian test attains
both divergences in (0.6) exactly.

## 9. Exact and model stress tests

### 9.1 Gaussian halfspace: exact success and exact blow-up

Let \(\mu_0=N(0,I_n)\) and

\[
                         S=\{x_1\ge0\}.
\]

Then \(g=1/2\),

\[
 v={1\over\sqrt{2\pi}}e_1,\qquad q=v,\qquad
 P=I-e_1e_1^T
\]

for the whole evolution.  For \(0\le t<1\),

\[
 \boxed{
 A_t=e_1e_1^T+(1-t)P,\qquad
 D_t=(1-t)^{-2}P.
 }                                                      \tag{9.1}
\]

One may choose \(C_t=(1-t)^{-1}P\).  The posterior mean satisfies

\[
                         da_t=P\,dW_t,                 \tag{9.2}
\]

and the posterior converges weakly at \(t=1\) to a Gaussian line parallel
to \(e_1\).  The first-coordinate marginal is untouched, so the weighted
perimeter is constantly

\[
                         \Pi_t(S)=(2\pi)^{-1/2}.        \tag{9.3}
\]

On the other hand,

\[
 \boxed{\begin{aligned}
 \int_0^{1-\varepsilon}\operatorname {tr}(D_tA_t)\,dt
   &=(n-1)\log(1/\varepsilon),\\
 \int_0^{1-\varepsilon}\operatorname {tr}D_t\,dt
   &=(n-1)(\varepsilon^{-1}-1).
 \end{aligned}}                                       \tag{9.4}
\]

Thus even the best-aligned example has no finite natural-parameter clock
at its rank-one limit.

### 9.2 Gaussian symmetric intervals: the zero-signal pathology

Let

\[
                         S=\{|x_1|\le r\}.
\]

At the standard Gaussian state, \(v=q=0\), whereas

\[
 H_{11}
 =\mathbb E[({\bf1}_{\{|Z|\le r\}}-g)Z^2]<0.           \tag{9.5}
\]

A full driver at the zero-signal instant tries to create \(v\) in the
\(e_1\) direction; as soon as \(v\ne0\), the hard controller deletes
that same driver direction.  In dimension one this is exactly the
nonexistence argument of Section 6.3.  In higher dimension the transverse
coordinates can continue to localize, but the \(e_1\) restart remains
mesh-dependent.  The frozen convention \(D=0\) at \(q=0\) is a valid
exact solution and performs no localization at all.

### 9.3 Cube halfspace

Let

\[
 \mu_0=\operatorname {Unif}[-\sqrt3,\sqrt3]^n,\qquad
 S=\{x_1\ge0\}.
\]

The product structure is preserved, \(q\) remains parallel to \(e_1\),
and the first marginal is never changed.  For each transverse coordinate
\(j\ge2\), write its posterior variance and centred third moment as
\(A_j\) and \(m_{3,j}\).  Then

\[
 \boxed{
 da_j=dW_j,\qquad
 dA_j={m_{3,j}\over A_j}\,dW_j-dt.
 }                                                      \tag{9.6}
\]

Initially \(m_{3,j}=0\), but a random tilt immediately destroys the
symmetry.  The first marginal retains variance one and

\[
                         \Pi_t(S)={1\over2\sqrt3}.      \tag{9.7}
\]

Thus the aligned cube test has the correct one-dimensional survivor.
It also shows why the constant trace drift is not a deterministic
eigenvalue drift: after tilting, each transverse variance has a
nontrivial martingale coefficient.

### 9.4 Coordinate caps for product exponentials

Let \(X_j=E_j-1\), where the \(E_j\) are independent
\(\operatorname {Exp}(1)\) variables.  For a set depending only on
\(X_1\), the vector \(q\) remains parallel to \(e_1\), the first marginal
is frozen, and every transverse factor obeys

\[
                         dA_j={m_{3,j}\over A_j}\,dW_j-dt.
                                                               \tag{9.8}
\]

At time zero, \(A_j=1\) and \(m_{3,j}=2\), so

\[
                         dA_j=2\,dW_j-dt
\]

to first order.  The covariance noise is therefore order one even at the
isotropic state.

There is no statewise universal bound on the protected variance.  A
linear exponential tilt can change a rate-one exponential factor to rate
\(\lambda>0\), giving variance \(\lambda^{-2}\).  If the set depends only
on that factor, (0.1) never acts on it, so this arbitrarily large variance
survives forever.  This does not by itself prove that a pure aggressive
flow started at the isotropic product produces a large survivor with
large probability.  It does prove that a protected-variance theorem
cannot follow only from the current posterior's mass-preservation
identities.

### 9.5 Product-exponential maximum: exact fast angular motion

The nonlinear product test displays the winner mechanism.  Let
\(X_i=E_i-1\) as above and

\[
 S=\{\max_i E_i\le r\},\qquad
 a=e^{-r},\qquad F=1-a,\qquad g=F^n.                  \tag{9.9}
\]

At the isotropic initial state,

\[
 v_i=-F^{n-1}ar,\qquad
 |q|=|v|=\sqrt n\,F^{n-1}ar,                          \tag{9.10}
\]

and \(u=-n^{-1/2}(1,\ldots,1)\).  For every unit
\(z\perp(1,\ldots,1)\), direct product integration gives

\[
 Hz=-F^{n-2}ar^2z,                                    \tag{9.11}
\]

and, using the centred exponential third moment \(2\),

\[
 {\cal T}(z)q=-2F^{n-1}ar\,z.                         \tag{9.12}
\]

At \(A=I\), \(C=P\), so (4.6) and (4.8) yield the exact initial rate

\[
 \boxed{
 {d\over dt}\operatorname {tr}[u]_0
 ={n-1\over n}\,{(r-2F)^2\over F^2}.
 }                                                      \tag{9.13}
\]

Choose the median maximum, \(F^n=1/2\).  Then

\[
 F=2^{-1/n},\qquad
 r=\log n-\log\log2+o(1),
\]

and (9.13) is \(\Theta((\log n)^2)\).  The protected direction moves by
order one on the \(1/\log^2n\) scale.  Because \(A=I\) at time zero, the
aggressive driver initially equals the ordinary Euclidean transverse
driver.  Its two inverse covariances give no leading-order improvement
against this angular instability.

### 9.6 Isotropic simplex cap: an exact benign cancellation

Let \(X=(X_1,\ldots,X_{n+1})\) be uniform on the barycentric simplex,
that is, Dirichlet\((1,\ldots,1)\) on

\[
 X_i\ge0,\qquad \sum_{i=1}^{n+1}X_i=1.
\]

After multiplying the centred affine tangent space by
\(\sqrt{(n+1)(n+2)}\), the law is isotropic.  Put

\[
 S=\{X_1\ge s\},\qquad
 a=1-s,\qquad g=a^n,
\]

and choose the median cap \(g=1/2\), so \(a=2^{-1/n}\).
Let \(u\) be the unit cap-normal direction in the isotropic tangent
space.  Then

\[
                         q=v=r_0u,\qquad
 r_0=g\,s\sqrt{n(n+2)}.                                \tag{9.14}
\]

For a unit tangent contrast \(z\perp u\), conditional Dirichlet moments
give

\[
                         Hz=g(a^2-1)z.                 \tag{9.15}
\]

The unconditional third moment is

\[
 \mathbb E[(z\cdot Y)^2(u\cdot Y)]
 =-{2\over n+3}\sqrt{\,{n+2\over n}\,}.                \tag{9.16}
\]

Consequently

\[
 {\cal T}(z)q
 =-{2g s(n+2)\over n+3}\,z
\]

and

\[
 (H-{\cal T}(\,\cdot\,)q)z
 ={g s\{(n+1)-(n+3)a\}\over n+3}\,z.                  \tag{9.17}
\]

Equations (4.6)--(4.8) now give the exact initial angular rate

\[
 \boxed{
 {d\over dt}\operatorname {tr}[u]_0
 ={n-1\over n(n+2)}
 \left\{{(n+1)-(n+3)a\over n+3}\right\}^2.
 }                                                      \tag{9.18}
\]

Since \(a=2^{-1/n}\), this is \(O(n^{-3})\).  The covariance-noise term
\(-{\cal T}(\,\cdot\,)q\), which would be missed by differentiating only
\(v\), nearly cancels the signed second-moment term.  Thus simplex caps
are benign for this controller at the initial state; they do not repair
the Gaussian-ball or product-maximum obstructions.

## 10. Surviving variance and the exact missing theorem

Suppose, without justification, that one could continue the process
through successive covariance faces until

\[
                         A_\infty=\sigma_\omega^2
                         u_\omega u_\omega^T.
\]

The trace identity would give only

\[
 \mathbb E\sigma_\omega^2
 =n-\mathbb E\int_0^\infty
       \operatorname {rank}(P_t)\,dt,                 \tag{10.1}
\]

where the rank is taken in the current affine support.  This is an
identity, not a universal upper bound.  Early covariance-martingale
fluctuations, rotation of the protected line, and noncanonical choices at
\(q=0\) all affect the integral.

The determinant identity does not improve (10.1): it controls a product
of active eigenvalues and becomes dominated by the smallest one.  The
regression inequality

\[
                         v^TA^{-1}v\le g(1-g)          \tag{10.2}
\]

controls the explained label variance, but it places no upper bound on
the Rayleigh variance in the direction \(A^{-1}v\).

For a terminal-needle Cheeger argument one would need, at minimum,

\[
 \mathbb E\sigma_\omega^2\le C
 \quad\hbox{or directly}\quad
 \mathbb E{1\over\sigma_\omega}\ge c.                 \tag{10.3}
\]

Neither (10.3) follows from (4.2), (5.1), (5.3), or (10.2).  The
arbitrarily tilted exponential factor in Section 9.4 disproves every
pathwise version of the first assertion.  The product-maximum calculation
shows that adaptive winner selection acts before substantial covariance
drift.  Controlling its probability-weighted contribution is a genuinely
new problem, not an Itô correction.

There are therefore two separate load-bearing statements still missing:

1. a singular-continuation theorem producing a canonical rank-one
   disintegration from the stopped or mesh-relaxed aggressive process,
   with perimeter lower semicontinuity; and
2. a universal estimate of the surviving one-dimensional variance for
   that disintegration.

Either statement would go substantially beyond the identities proved in
this report.  Assuming a terminal needle, or assuming (10.3), would insert
the KLS-strength step rather than prove it.

## 11. Final audit table

| Claim | Status | Reason |
|---|---:|---|
| Density SDE | proved on bounded stops | (2.1)--(2.3), including normalization and log Itô term |
| Set-mass conservation | proved away from \(q=0\); approximate for the soft controller; conditional for identified relaxed limits | \(D v=0\), (4.1)--(4.2), (6.2)--(6.6) |
| Covariance drift \(-P\) | proved | \(ADA=P\), (4.3) |
| \(v\)-SDE | proved | cross-variation retained in (3.3)--(3.4) |
| Feedback \(q=A^{-1}v\)-SDE | proved | inverse-matrix and cross corrections in (4.5)--(4.8) |
| Trace and determinant SDEs | proved | (5.1), (5.3), (5.4) |
| Perimeter local-martingale SDE | proved on smooth bounded stops | (5.6)--(5.7) |
| Universal time to first covariance face | conditional stopped theorem | (7.1)--(7.2) |
| Canonical passage through \(v=0\) | false for naive restart; relaxed only | Section 6 |
| Finite stochastic clock at the first face | false | Gaussian halfspace, (9.4) |
| Rank-one collapse from trace drift | false as an inference | moving-line identity (7.3) |
| Universal surviving variance | not proved; pathwise version false | Sections 9.4 and 10 |
| KLS consequence | not obtained | requires both missing statements in Section 10 |
