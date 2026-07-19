# Terminal-needle inverse for the covariance-normalized hard controller

## 0. Audit verdict

Let \(\mu_0\) be log-concave and let \(S\) be a fixed Borel set.  The
proposed controller deletes, in covariance-normalized coordinates, the
single mode which changes \(\mu_t(S)\).  Away from zero label covariance
its algebra is correct: the density equation, pathwise mass conservation,
covariance drift, and stopped perimeter-martingale identity all hold.

There are, however, two logically prior gaps in the proposed terminal
argument.

1. The convention which applies full whitening when the label covariance
   is zero is not merely nonsmooth.  It need not admit a continuous weak
   solution.  Section 3 gives a one-dimensional isotropic Gaussian example
   for which the defining SDE has no solution on any positive time
   interval.
2. Exponential decay of expected log covariance volume does not prove
   convergence to a one-dimensional posterior.  Even conditional on such
   convergence, passage of ambient perimeter to the perimeter of the trace
   requires a separate lower-semicontinuity and uniform-integrability
   theorem.  These facts are isolated precisely in Section 5.

Thus the terminal-needle lemma is not currently a theorem about the stated
controller, because for legitimate sets the controller itself may be
undefined.  Conditional on a valid terminal line disintegration, this note
does prove three inverse branches:

* a fixed projective direction, or a positive-probability fixed cap, gives
  a universal lower bound immediately from covariance disintegration;
* exactly concurrent terminal barycenters are impossible in affine
  dimension at least two;
* a high-dimensional family which carries a fixed covariance fraction in
  every direction is impossible by the thin-shell theorem.  Quantitatively,
  the number of eigenvalues at least \(\alpha\) of the internal covariance
  matrix is bounded by a number depending only on \(\alpha\).

The remaining case is a multiscale family: terminal directions spread over
more and more dimensions while the covariance carried by each spectral
scale tends to zero, and terminal barycenters remain sufficiently dispersed
to avoid the concurrent-center theorem.  None of the arguments below
controls that case, and no KLS-strength assertion is assumed.

Throughout Sections 1--5 the initial density is positive, smooth, and
log-concave on \(\mathbb R^k\), with all moments.  All stochastic identities
are asserted only before explicit coefficient and integrability stops.
Sections 6--9 concern an abstract terminal mixture and need no smoothness.

## 1. Natural-parameter process and exact density equation

Let

\[
 p_t(x)=Z_t^{-1}\exp\left\{c_t\cdot x-
             \frac12 x^TQ_tx\right\}p_0(x),
 \qquad Q_t\succeq0.                                      \tag{1.1}
\]

For a predictable positive-semidefinite matrix \(D_t\), choose a
predictable square root \(C_t\) satisfying

\[
                         C_tC_t^T=D_t                       \tag{1.2}
\]

and set

\[
 dc_t=C_t\,dW_t+D_ta_t\,dt,
 \qquad dQ_t=D_t\,dt,
 \qquad a_t=\int x p_t(x)\,dx.                            \tag{1.3}
\]

On a stop on which the natural parameters remain in a compact subset of
the moment-generating domain and the required moments are bounded, direct
Itô differentiation gives

\[
 dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle.             \tag{1.4}
\]

For completeness, the quadratic variation of \(c_t\cdot x\) contributes
\(\frac12x^TD_tx\,dt\), which cancels the differential of
\(-\frac12x^TQ_tx\).  The drift \(D_ta_t\) and the normalization \(Z_t\)
then cancel all remaining finite-variation terms.  This verifies (1.4)
without an isoperimetric or spectral hypothesis.

Define, with \(Y=X-a_t\),

\[
 A_t=\mathbb E_tYY^T,
 \quad g_t=\mathbb E_t{\bf1}_S,
 \quad v_t=\mathbb E_t[({\bf1}_S-g_t)Y],                  \tag{1.5}
\]

and

\[
 H_t=\mathbb E_t[({\bf1}_S-g_t)YY^T],
 \qquad
 \mathcal T_t(z)=\mathbb E_t[YY^T\langle Y,z\rangle].    \tag{1.6}
\]

Stochastic differentiation of these moments gives

\[
 \boxed{
 \begin{aligned}
 da_t&=A_tC_t\,dW_t,\\
 dg_t&=v_t^TC_t\,dW_t,\\
 dv_t&=H_tC_t\,dW_t-A_tD_tv_t\,dt,\\
 dA_t&=\mathcal T_t(C_t\,dW_t)-A_tD_tA_t\,dt.
 \end{aligned}}                                           \tag{1.7}
\]

The drift in the third line is sometimes omitted.  It is the cross
variation in differentiating
\(v_t=\mathbb E_t[{\mathbf{1}}_SX]-g_ta_t\).  It vanishes for the hard controller,
but it is present for a general driver.

## 2. The hard controller away from its singular state

Assume \(A_t\succ0\) and \(v_t\ne0\).  Put

\[
 q_t=v_t^TA_t^{-1}v_t,
 \qquad
 w_t={A_t^{-1/2}v_t\over\sqrt{q_t}},
 \qquad
 P_t=I-w_tw_t^T,                                        \tag{2.1}
\]

and

\[
 D_t=A_t^{-1/2}P_tA_t^{-1/2}
 =A_t^{-1}-{A_t^{-1}v_tv_t^TA_t^{-1}\over q_t}.          \tag{2.2}
\]

Then

\[
 D_tv_t=0,
 \qquad C_t^Tv_t=0,                                    \tag{2.3}
\]

the second identity following from
\(|C_t^Tv_t|^2=v_t^TD_tv_t=0\).  Equations (1.7) therefore imply

\[
                  g_t=g_0\quad\hbox{pathwise}            \tag{2.4}
\]

throughout every interval on which the solution exists and (2.2) is used.
The covariance drift is

\[
 -A_tD_tA_t=-A_t+{v_tv_t^T\over q_t}.                   \tag{2.5}
\]

In particular,

\[
 \operatorname {tr}(A_tD_t)=k-1.                        \tag{2.6}
\]

If \(T_{t,j}=\mathcal T_t(C_te_j)\), matrix Itô calculus yields

\[
 \begin{aligned}
 d\log\det A_t
 &=\sum_j\operatorname {tr}(A_t^{-1}T_{t,j})\,dW_{t,j}
 -(k-1)\,dt\\
 &\quad-{1\over2}\sum_j
 \left\|A_t^{-1/2}T_{t,j}A_t^{-1/2}\right\|_{HS}^2dt.
 \end{aligned}                                               \tag{2.7}
\]

Thus the advertised log-volume drift is correct.  It controls a product of
eigenvalues, not the number of surviving eigenvalues or the largest
surviving variance.

Fix \(R<\infty\) and stop while

\[
 R^{-1}I\preceq A_t\preceq RI,
 \qquad |v_t|\ge R^{-1},                                \tag{2.8}
\]

and while the natural parameters and moments needed in (1.7) remain in a
compact set.  On this region (2.2) has constant rank and its principal
square root is a smooth function of the natural parameters.  Standard
finite-dimensional SDE theory then supplies a unique strong solution up to
the stop.  This is the maximal well-posedness conclusion justified by the
formula.  Bounded covariance alone does not remove the singularity at
\(v=0\).

## 3. The full-whitening restart at \(v=0\) has no weak solution

The proposed convention is

\[
 D=A^{-1}\quad\hbox{when }v=0,
 \qquad
 D=A^{-1/2}(I-ww^T)A^{-1/2}\quad\hbox{when }v\ne0.
                                                               \tag{3.1}
\]

The following example shows that (3.1) is not a stochastic differential
equation in the usual weak sense.

### Theorem 3.1 (one-dimensional zero-signal nonexistence)

Let \(k=1\), let \(\mu_0=N(0,1)\), and let

\[
                         S=[-r,r],\qquad r>0.           \tag{3.2}
\]

There is no continuous weak solution of (1.1)--(1.3) with the convention
(3.1) on any deterministic interval \([0,T]\) with \(T>0\).

#### Proof

Every posterior in (1.1) is Gaussian \(N(a,A)\), where

\[
                         A=(1+Q)^{-1},\qquad a=Ac.      \tag{3.3}
\]

Write \(G(a,A)=\mathbb P_{N(a,A)}\{|X|\le r\}\).  Gaussian score
differentiation gives

\[
 v=\mathbb E[({\bf1}_S-G)(X-a)]=A\,\partial_aG(a,A).   \tag{3.4}
\]

Moreover,

\[
 \partial_aG(a,A)={1\over\sqrt A}
 \left[
 \varphi\left({-r-a\over\sqrt A}\right)
 -\varphi\left({r-a\over\sqrt A}\right)
 \right].                                               \tag{3.5}
\]

For \(a>0\), the first Gaussian density in brackets is strictly smaller
than the second; for \(a<0\) it is strictly larger.  Hence

\[
                         v=0\quad\Longleftrightarrow
                         a=0\quad\Longleftrightarrow c=0. \tag{3.6}
\]

In one dimension the projection in (3.1) is zero whenever \(v\ne0\).
Consequently the proposed equations reduce exactly to

\[
 \begin{aligned}
 dc_t&=\sqrt{1+Q_t}\,{\bf1}_{\{c_t=0\}}\,dW_t,\\
 dQ_t&=(1+Q_t){\bf1}_{\{c_t=0\}}\,dt,
 \qquad c_0=Q_0=0.
 \end{aligned}                                         \tag{3.7}
\]

The drift \(D_ta_t\) in \(dc_t\) is identically zero: if \(c_t=0\), then
\(a_t=0\), while if \(c_t\ne0\), then \(D_t=0\).

Suppose a continuous weak solution existed.  Its quadratic variation would
be

\[
 [c]_t=\int_0^t(1+Q_s){\bf1}_{\{c_s=0\}}\,ds.          \tag{3.8}
\]

The occupation-density formula for a continuous semimartingale gives

\[
 \int_0^t{\bf1}_{\{c_s=0\}}\,d[c]_s
 =\int_{\mathbb R}{\bf1}_{\{0\}}(x)L_t^x(c)\,dx=0.    \tag{3.9}
\]

But (3.8) is supported on \(\{c=0\}\), so (3.9) says \([c]_t=0\).
The local martingale \(c\) is therefore identically zero.  Substitution
back into (3.7) gives \(Q_t=e^t-1\) and hence

\[
                         [c]_t=\int_0^te^sds=e^t-1>0,
                                                               \tag{3.10}
\]

a contradiction.  \(\square\)

This is an isotropic log-concave starting measure and a nontrivial Borel
set.  Stopping all covariances and moments in a compact range does not help,
because the contradiction occurs on every positive interval before such a
stop.

There is also a deterministic reason the singularity cannot be repaired by
a nonzero continuous exact controller.  Suppose \(D(v)\) is continuous at
zero and \(D(v)v=0\) for all \(v\ne0\).  For every unit \(e\), letting
\(v=se\to0\) gives \(D(0)e=0\).  Thus

\[
                         D(0)=0.                       \tag{3.11}
\]

The only continuous exact mass-preserving extension freezes a zero-signal
state.  Such a frozen posterior does not converge to a needle.  A soft
controller or a separately chosen protected direction may be well posed,
but either is a different construction and requires a new analysis.

## 4. Audit of the perimeter martingale

Let \(S\) have compact \(C^1\) boundary.  Its weighted perimeter is

\[
 P_t(S)=\int_{\partial S}p_t(x)\,d\mathcal H^{k-1}(x). \tag{4.1}
\]

On a stop for which stochastic Fubini is justified, integration of (1.4)
over the fixed boundary gives

\[
 dP_t(S)=\left\langle C_t^T
 \int_{\partial S}(x-a_t)p_t(x)\,d\mathcal H^{k-1}(x),
 dW_t\right\rangle.                                    \tag{4.2}
\]

Hence stopped perimeter is a local martingale.  If

\[
 \mathbb E\int_0^T\left|C_t^T
 \int_{\partial S}(x-a_t)p_t(x)\,d\mathcal H^{k-1}(x)
 \right|^2dt<\infty,                                   \tag{4.3}
\]

it is a true martingale through time \(T\).  Along an increasing
localizing sequence, nonnegativity makes the unrestricted process a
supermartingale and gives

\[
                         \mathbb EP_t(S)\le P_0(S).     \tag{4.4}
\]

For a fixed finite-perimeter set, (4.2) continues to hold with
\(\partial S\) replaced by its reduced boundary, provided the corresponding
surface integrals satisfy (4.3).  This follows directly by stochastic
Fubini for the fixed rectifiable measure
\(\mathcal H^{k-1}\mathbin\llcorner\partial^*S\); it does not require
changing the set.

An approximation argument in which \(S_j\to S\) is more delicate than the
usual weighted-perimeter approximation, because the controller itself
depends on \(S_j\) through \(v_t\).  Stability of the controlled SDE under
that approximation is not a consequence of (4.2), especially near
\(v=0\).  Thus (4.4) is fully verified for regular sets under stopped
integrability, but the blanket extension to every Borel set has not been
proved for this controller.

## 5. What is required for the terminal reduction

Assume, solely for this section, that a valid process exists globally and
that the following four assertions have been proved.

1. Almost surely \(p_tdx\) converges weakly to a nondegenerate
   one-dimensional log-concave probability \(\nu_\omega\) supported on an
   affine line \(L_\omega\) of projective direction \(u_\omega\).
2. The set mass passes to the limit:
   \(\nu_\omega(S)=g_0\).
3. The stopped perimeters are uniformly integrable in the direction needed
   below and

   \[
   \liminf_{t\to\infty}P_t(S)
   \ge \nu_\omega^+(S\cap L_\omega).                  \tag{5.1}
   \]
4. The posterior interpretation holds at terminal time, so that
   \(\nu_\omega\) is a regular conditional law of the original random
   vector given the terminal observation sigma-field.

For a smooth transverse crossing, (5.1) has the local form

\[
 \lim P_t(S)=\sum_{x\in\partial S\cap L_\omega}
 {\rho_\omega(x)\over|n_S(x)\cdot u_\omega|}
 \ge\sum_{x\in\partial(S\cap L_\omega)}\rho_\omega(x).
                                                               \tag{5.2}
\]

Tangencies and nonsmooth traces require lower semicontinuity.  Formula
(5.2) explains the proposed inequality but is not, by itself, a proof of
(5.1) for a random collapsing sequence.

Every one-dimensional log-concave probability of standard deviation
\(\sigma\) has Cheeger constant at least \(c_1/\sigma\), with a universal
\(c_1>0\).  This follows from one-dimensional half-line isoperimetry and
the standard bound of the density at a median by a universal multiple of
\(1/\sigma\).  Therefore (2.4) and (5.1) would imply

\[
 P_0(S)\ge c_1\min(g_0,1-g_0)
                  \mathbb E{1\over\sigma_\omega}.      \tag{5.3}
\]

Here Fatou is used in the valid direction:
\(\mathbb E\liminf P_t\le\liminf\mathbb EP_t\le P_0\).
The factor \(g_0(1-g_0)\) can replace the minimum after reducing the
constant, since \(\min(g,1-g)\ge g(1-g)\).

None of the four terminal assertions follows from (2.7).  In particular,
decay of \(\det A_t\) allows rotating small eigenspaces, collapse to a
point on the measure-theoretic boundary of \(S\), or loss of moment and
surface uniform integrability.  A point limit is not ruled out merely by
\(0<g_0<1\): shrinking centered measures can retain a fixed fraction of a
set when their centers lie on its boundary.

Under assertion 4, let

\[
 b_\omega=\int x\,d\nu_\omega(x),
 \qquad
 \operatorname {Cov}(\nu_\omega)
 =\sigma_\omega^2u_\omega u_\omega^T.                 \tag{5.4}
\]

The Hilbert-space law of total covariance gives, for isotropic \(\mu_0\),

\[
 I=\operatorname {Cov}(b_\omega)
   +K,
 \qquad
 K:=\mathbb E[\sigma_\omega^2u_\omega u_\omega^T]
 \preceq I.                                            \tag{5.5}
\]

This identity is exact if the terminal laws are actual conditional laws.
Deriving it only from weak convergence would require uniform integrability
of second moments.

## 6. The fixed-direction branch

The first inverse branch is purely algebraic and loses no dimension.

### Proposition 6.1 (a projective cap of terminal directions)

Assume (5.5).  Let \(e\) be a deterministic unit vector and let
\(E\) be an event of probability \(q>0\) such that

\[
                         |u_\omega\cdot e|\ge\delta>0
                         \quad\hbox{on }E.              \tag{6.1}
\]

Then

\[
 \mathbb E\left[{\bf1}_E{1\over\sigma_\omega}\right]
 \ge \delta q^{3/2}.                                  \tag{6.2}
\]

In particular, if all terminal directions equal one projective direction,
then

\[
                         \mathbb E{1\over\sigma_\omega}\ge1. \tag{6.3}
\]

#### Proof

Taking the quadratic form of (5.5) in direction \(e\) gives

\[
 \delta^2\mathbb E[{\mathbf{1}}_E\sigma_\omega^2]
 \le\mathbb E[\sigma_\omega^2(u_\omega\cdot e)^2]
 \le1.                                                  \tag{6.4}
\]

Cauchy--Schwarz and then Cauchy applied to \(\sigma\) and \(1/\sigma\)
give

\[
 \mathbb E[{\mathbf{1}}_E\sigma]
 \le\sqrt{q\mathbb E[{\mathbf{1}}_E\sigma^2]}
 \le{\sqrt q\over\delta},
 \qquad
 q^2\le
 \mathbb E[{\mathbf{1}}_E\sigma]\,
 \mathbb E[{\mathbf{1}}_E/\sigma].                           \tag{6.5}
\]

Combining them proves (6.2).  \(\square\)

Thus a bad terminal family cannot put a fixed positive probability in a
fixed projective cap of fixed aperture.  The adaptively selected direction
must escape every deterministic cap in the weighted sense of (6.2).

## 7. Exact concurrence is impossible

The next result uses the one-dimensional log-concavity of the terminal
laws, rather than only covariance disintegration.

### Proposition 7.1 (no common terminal barycenter)

Let \(\mu\) be a full-dimensional log-concave probability on
\(\mathbb R^k\), \(k\ge2\).  There is no probability mixture

\[
                         \mu=\int\nu_\omega\,d\pi(\omega) \tag{7.1}
\]

such that every \(\nu_\omega\) is a nondegenerate one-dimensional
log-concave probability supported on a line through one fixed point \(b\)
and has barycenter exactly \(b\).

#### Proof

Choose a unit orientation on each line and write points on that line as
\(b+tu_\omega\).  The arclength density \(f_\omega(t)\) is log-concave,
nondegenerate, and has mean zero.  Hence zero lies in the interior of its
support and its upper-semicontinuous log-concave version satisfies

\[
                         f_\omega(0)>0.                 \tag{7.2}
\]

Consequently

\[
 \lim_{\varepsilon\downarrow0}{1\over\varepsilon}
 \nu_\omega(B(b,\varepsilon))=2f_\omega(0)>0.          \tag{7.3}
\]

Fatou's lemma applied to (7.1) yields

\[
 \liminf_{\varepsilon\downarrow0}{\mu(B(b,\varepsilon))
                                      \over\varepsilon}
 \ge2\int f_\omega(0)\,d\pi(\omega)>0.                \tag{7.4}
\]

On the other hand, every integrable full-dimensional log-concave density
is locally bounded.  Thus for some finite \(M\) and all small
\(\varepsilon\),

\[
                         \mu(B(b,\varepsilon))
                         \le M\,\kappa_k\varepsilon^k, \tag{7.5}
\]

which makes the left side of (7.4) zero when \(k\ge2\).  This is a
contradiction.  \(\square\)

If a terminal law is degenerate, \(1/\sigma=+\infty\) and the desired
inverse estimate is already true.  Thus degeneracies do not weaken the
proposition for the terminal-needle application.

The proof identifies the geometric reason exact concurrence fails.  A
one-dimensional mean-zero log-concave law places order-\(\varepsilon\)
mass in an \(\varepsilon\)-ball about its barycenter.  A
full-dimensional log-concave density places only order
\(\varepsilon^k\) mass there.  Quantitative near-concurrence requires an
annular, rather than pointwise, version of this comparison; the next
section supplies one high-rank form.

## 8. A high-rank inverse from thin shell

We use two dimension-free facts.

**One-dimensional central mass.**  For every \(a>0\), there is a number
\(\beta(a)>0\) such that every mean-zero, variance-one log-concave
probability \(\lambda\) on \(\mathbb R\) satisfies

\[
                         \lambda([-a,a])\ge\beta(a).    \tag{8.1}
\]

One self-contained proof is by compactness.  Standardized one-dimensional
log-concave laws have a uniform exponential tail, hence form a tight family
whose first two moments are uniformly integrable.  A weak limit is again
mean-zero, variance-one, and log-concave.  Its support contains zero in its
interior, so it has positive mass in \((-a,a)\).  A sequence violating
(8.1) would contradict the open-set part of the Portmanteau theorem.

**Thin shell.**  There is a universal \(C_{TS}\) such that every isotropic
log-concave \(Y\) in \(\mathbb R^r\) satisfies

\[
                         \mathbb E(|Y|-\sqrt r)^2\le C_{TS}. \tag{8.2}
\]

This is the dimension-free thin-shell theorem, used here only for a
log-concave marginal.

### Theorem 8.1 (bounded spectral multiplicity of terminal covariance)

Let an isotropic log-concave \(\mu\) admit a terminal mixture satisfying
(5.5).  For every \(\alpha\in(0,1]\), the number of eigenvalues of

\[
                         K=\mathbb E[\sigma^2uu^T]      \tag{8.3}
\]

which are at least \(\alpha\) is at most a finite number
\(N(\alpha)\) depending only on \(\alpha\), not on the ambient dimension.
One admissible explicit definition is

\[
 \begin{gathered}
 s_\alpha=\sqrt{1-\alpha/2},
 \qquad d_\alpha=1-s_\alpha,
 \qquad q_\alpha={\alpha\over2-\alpha},
 \qquad M_\alpha=\sqrt{2/q_\alpha},\\
 a_\alpha={d_\alpha\over2M_\alpha},\\
 N(\alpha)=
 \left\lceil{8C_{TS}\over
        d_\alpha^2\beta(a_\alpha)q_\alpha}\right\rceil.
 \end{gathered}                                             \tag{8.4}
\]

#### Proof

Let \(E\) be the span of eigenvectors of \(K\) with eigenvalue at least
\(\alpha\), and put \(r=\dim E\).  Write

\[
                         Y=P_EX,
 \quad y_\omega=P_Eb_\omega,
 \quad \tau_\omega=\sigma_\omega|P_Eu_\omega|.        \tag{8.5}
\]

The marginal law of \(Y\) is isotropic and log-concave on \(E\).  From
(5.5),

\[
 \mathbb E|y_\omega|^2
 =\operatorname {tr}\big(P_E(I-K)P_E\big)
 \le(1-\alpha)r,                                      \tag{8.6}
\]

while

\[
                         \mathbb E\tau_\omega^2
                         =\operatorname {tr}(P_EK)
                         \le r.                        \tag{8.7}
\]

Markov's inequality and the definitions in (8.4) give

\[
 \begin{aligned}
 \mathbb P\{|y_\omega|\le s_\alpha\sqrt r\}
 &\ge1-{1-\alpha\over s_\alpha^2}=q_\alpha,\\
 \mathbb P\{\tau_\omega>M_\alpha\sqrt r\}
 &\le M_\alpha^{-2}=q_\alpha/2.
 \end{aligned}                                         \tag{8.8}
\]

Hence the event \(G\) on which both desired bounds hold has probability at
least \(q_\alpha/2\).  Conditional on \(\omega\), the projected law is a
possibly degenerate one-dimensional log-concave law with barycenter
\(y_\omega\) and standard deviation \(\tau_\omega\).  By (8.1), on \(G\)
it puts at least \(\beta(a_\alpha)\) mass in

\[
 |Y-y_\omega|\le a_\alpha\tau_\omega.
\]

The choice of \(a_\alpha\) makes this set lie inside

\[
 |Y|\le(s_\alpha+d_\alpha/2)\sqrt r
       =(1-d_\alpha/2)\sqrt r.                         \tag{8.9}
\]

Therefore

\[
 \mathbb P\{|Y|\le(1-d_\alpha/2)\sqrt r\}
 \ge{q_\alpha\over2}\beta(a_\alpha).                 \tag{8.10}
\]

On the other hand, (8.2) and Chebyshev give

\[
 \mathbb P\{|Y|\le(1-d_\alpha/2)\sqrt r\}
 \le{4C_{TS}\over d_\alpha^2r}.                       \tag{8.11}
\]

Combining (8.10)--(8.11) gives \(r\le N(\alpha)\).
\(\square\)

### Corollary 8.2 (uniformly high-rank directions close the inverse)

Suppose all terminal directions lie in a deterministic subspace \(E\),
and for one fixed universal \(\alpha>0\),

\[
                         K\succeq\alpha P_E.           \tag{8.12}
\]

Then

\[
 \dim E\le N(\alpha),
 \qquad
 \mathbb E{1\over\sigma}
 \ge{1\over\sqrt{\mathbb E\sigma^2}}
 ={1\over\sqrt{\operatorname {tr}K}}
 \ge{1\over\sqrt{N(\alpha)}}.                        \tag{8.13}
\]

Thus a covariance-saturating family of long needles cannot evade the
inverse estimate merely by distributing its directions over a high-rank
frame.  For example, the algebraic model
\(\sigma^2=k\) with \(u\) isotropic on the sphere has \(K=I\); Theorem 8.1
shows that such a model cannot be a terminal log-concave conditional
mixture in arbitrarily high dimension.

More generally, if an event \(H\) satisfies

\[
 \mathbb E[\sigma^2uu^T{\mathbf{1}}_H]\succeq\alpha P_E,     \tag{8.14}
\]

then \(\dim E\le N(\alpha)\).  This is the precise high-rank obstruction.
It does not treat a multiscale spectrum for which every fixed spectral
threshold has bounded multiplicity while the number of progressively
smaller eigenvalues diverges.

The spectral restriction can be made polynomial.  The elementary
one-dimensional density bounds for a standardized log-concave law give a
universal \(c_0>0\) such that

\[
                         \beta(a)\ge c_0a,
                         \qquad 0<a\le1.                \tag{8.15}
\]

Indeed, a mean-zero variance-one log-concave density is bounded below by a
universal positive constant on a universal interval about zero; integrate
over \([-a,a]\), reducing \(c_0\) for the remaining values of \(a\le1\).
Since

\[
 {\alpha\over4}\le d_\alpha\le{\alpha\over2},
 \qquad {\alpha\over2}\le q_\alpha\le\alpha,
 \qquad a_\alpha\ge{\alpha^{3/2}\over16},              \tag{8.16}
\]

Theorem 8.1 gives

\[
                         N_K(\alpha)le C\alpha^{-9/2}, \tag{8.17}
\]

where \(N_K(\alpha)\) is the eigenvalue-counting function of \(K\).
Integrating the counting function and using \(N_K(\alpha)\le k\) gives

\[
 \mathbb E\sigma^2=\operatorname {tr}K
 =\int_0^1N_K(\alpha)\,d\alpha
 \le Ck^{7/9},
 \qquad
 \mathbb E{1\over\sigma}\ge c k^{-7/18}.              \tag{8.18}
\]

This is dimension dependent and therefore does not approach the required
terminal lemma.  Its purpose is diagnostic: log-concavity and thin shell do
exclude the algebraic \(\operatorname {tr}K=k\) example by an unbounded
factor, but leave a genuine small-eigenvalue accumulation problem.

The proof also quantifies approximate concurrence.  On the spectral
subspace \(E\), (8.6) says the terminal barycenters have mean-square radius
at most \(\sqrt{1-\alpha}\sqrt r\).  A central piece of every projected
needle would then put fixed mass strictly inside the thin shell, which is
impossible for large \(r\).

## 9. A product-exponential terminal stress test

The following exact calculation explains why no pointwise bound on a
terminal needle can be true, even for an isotropic product measure.

Let \(Z_1,\ldots,Z_k\) be independent rate-one exponentials and put
\(X_i=Z_i-1\).  Then \(X\) is isotropic and log-concave.  Choose \(L_k\) by

\[
                         (1-e^{-L_k})^k={1\over2}       \tag{9.1}
\]

and set

\[
                         S=\{\max_iZ_i\ge L_k\}.       \tag{9.2}
\]

Consider any terminal posterior which is supported on the \(i\)-th
coordinate ray, whose other coordinates are below \(L_k\), and whose
surviving coordinate has exponential rate \(r\).  Its trace of \(S\) is the
half-line \(\{Z_i\ge L_k\}\).  Exact terminal mass preservation forces

\[
                         e^{-rL_k}={1\over2},
 \qquad r={\log2\over L_k},                            \tag{9.3}
\]

so its standard deviation is

\[
                         \sigma={L_k\over\log2}
                         \asymp\log k.                 \tag{9.4}
\]

This is not a counterexample to the terminal-needle lemma, because it has
not been proved that the adaptive controller terminates in these
coordinate-ray states with substantial probability.  It is an exact
stress test: any successful inverse theorem must use the path law and the
distribution of terminal barycenters and directions.  A line-by-line
claim \(\sigma\le C\) is false.

## 10. Final status of this route

The controlled identities that survive audit are:

\[
 \begin{gathered}
 dp_t=p_t\langle x-a_t,C_t\,dW_t\rangle,\\
 D_tv_t=0\Longrightarrow \mu_t(S)=\mu_0(S),\\
 -A_tD_tA_t=-A_t+v_tv_t^T/(v_t^TA_t^{-1}v_t),\\
 dP_t(S)\text{ is a stopped local martingale}.
 \end{gathered}                                         \tag{10.1}
\]

They are valid on stopped intervals with \(A\succ0\) and \(v\ne0\).  The
specific full-whitening restart at \(v=0\) is invalid by Theorem 3.1, and
the continuous exact restart freezes by (3.11).  Consequently a global
terminal process has not been constructed in the stated generality.

Even after postulating a valid terminal line disintegration, the inverse
estimate is proved here only in the following regimes:

1. a fixed projective cap of terminal directions with fixed positive
   probability (Proposition 6.1);
2. exact concurrence, which is impossible in dimension at least two
   (Proposition 7.1);
3. a direction family carrying a fixed covariance fraction on a
   high-dimensional subspace, which is reduced to bounded dimension by
   Theorem 8.1 and then closed by (8.13).

The unresolved regime has no fixed cap, no common barycenter, and no fixed
spectral threshold of unbounded multiplicity.  Proving a dimension-free
lower bound for \(\mathbb E(1/\sigma)\) there would require a new
multiscale incidence or historical estimate.  Covariance disintegration,
thin shell, Paouris tails, and the stopped martingale identities above do
not by themselves provide it.
