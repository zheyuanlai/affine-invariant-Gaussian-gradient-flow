# Normalized and bounded mass-preserving localization

## 0. Verdict

Fix a Borel set \(E\), and for a stochastic-localization posterior
\(\mu_t\) write

\[
 a_t={\mathbb E}_tX,\qquad A_t=\operatorname {Cov}_tX,
 \qquad b_t=\operatorname {Cov}_t({\bf1}_E,X).
\]

There are two natural mass-preserving controls.  The first is the
covariance-normalized control proposed in the assignment:

\[
 q_t=A_t^{-1/2}b_t,\qquad
 P_t=I-{q_tq_t^T\over |q_t|^2},\qquad
 D_t=A_t^{-1/2}P_tA_t^{-1/2}.                         \tag{0.1}
\]

Away from \(q_t=0\), (0.1) gives a rigorous stopped strong solution and
the exact identities

\[
 D_tb_t=0,\qquad A_tD_tA_t=A_t^{1/2}P_tA_t^{1/2},
 \qquad {d\mathsf B_t\over dt}=D_t,                   \tag{0.2}
\]

where \(\mathsf B_t\) is the accumulated quadratic tilt.  In particular,
\(\mu_t(E)\) is pathwise constant.  This normalization genuinely removes
the finite-time clock obstruction of the more aggressive control
\(A^{-1}PA^{-1}\):

\[
 \operatorname {tr}(D_tA_t)=n-1.                     \tag{0.3}
\]

Thus the posterior-averaged quadratic variation of the log-density is
constant.  In the Gaussian halfspace test the transverse covariance is
\(e^{-t}I\), rather than \((1-t)I\), and every natural parameter is finite
at every finite time.  The quadratic tilt still diverges as
\(t\to\infty\), which is harmless for finite-time well-posedness.

Away from a zero of \(q_t\), standardized third-moment bounds also rule
out a covariance face at every finite time.  Thus the normalized process
extends until its first zero-signal time, with no finite covariance-clock
stop.

The normalized control nevertheless does not complete the argument.  If
\(\beta_1(t)\le\beta_2(t)\le\cdots\) are the eigenvalues of
\(\mathsf B_t\), then

\[
 \boxed{\quad
 \beta_2(t)\ge {1\over2}\int_0^t{ds\over\|A_s\|_{\rm op}}.
 \quad}                                               \tag{0.4}
\]

Consequently a path on which \(\sup_s\|A_s\|_{\rm op}<\infty\) has at
most one weak-curvature direction in the limit.  But (0.4) gives no
dimension-free curvature at a fixed universal time, and it gives no bound
on the variance in the remaining direction.

There is an even simpler bounded control.  Put

\[
 u_t={b_t\over|b_t|},\qquad P_t=I-u_tu_t^T,qquad
 C_t=P_t,\quad D_t=C_tC_t^T=P_t.                      \tag{0.5}
\]

It also preserves \(\mu_t(E)\) on every interval on which \(b_t\ne0\).
For any such path that exists through time \(t\),

\[
 \mathsf B_t=\int_0^tP_sds=tI-Q_t,qquad
 Q_t=\int_0^tu_su_s^Tds,qquad \operatorname {tr}Q_t=t,
\]

and hence

\[
 \boxed{\qquad \lambda_2(\mathsf B_t)\ge t/2.\qquad} \tag{0.6}
\]

The deterministic spectral statement is unconditional for every
measurable unit-direction path, but the hard SDE need not supply such a
path through a zero of \(b_t\).  Conditional on existence through
\(t=1\), the rank-one-defect Lichnerowicz lemma gives

\[
 C_P(\mu_1)\le C\bigl(1+
        \operatorname {Var}_{\mu_1}\langle X,w_1\rangle\bigr), \tag{0.7}
\]

where \(w_1\) is the weakest eigenvector of \(\mathsf B_1\).  This is a
clean reduction, but its remaining average is not controlled by the
usual mixture identity.  That identity only gives
\({\mathbb E}A_1\preceq I\) in the isotropic case; for an adaptively chosen
\(w_1\) it yields merely \({\mathbb E}w_1^TA_1w_1\le n\).

The obstruction is visible, rather than formal.  For the median maximum
of \(n\) independent shifted exponentials, the bounded controller has
initial angular quadratic-variation rate \(\Theta((\log n)^2)\), and

\[
 {d\over dt}{\mathbb E}\,[u_t^TA_tu_t]\big|_{t=0}
       ={4(n-1)R\over nF}=\Theta(\log n),              \tag{0.8}
\]

with the notation of Section 8.  The signal-weighted quantity
\(b_t^TA_tb_t\) has relative initial growth \(\Theta((\log n)^2)\).
Thus adaptive selection and covariance noise correlate in exactly the
uncontrolled direction.

Finally, the convention demanded at zero signal, namely \(P=I\) when
\(q=0\) (or \(b=0\)), has no continuous weak solution in general.  A
one-dimensional Gaussian symmetric interval gives an occupation-time
contradiction.  Hence neither hard control is a globally defined process
under that convention.  A rank-\((n-1)\) selector, a soft control, or a
mesh relaxation is required at zero signal, and each is noncanonical or
loses exact mass preservation.

The precise surviving gate is

\[
 \boxed{
 \inf_{t=1}\ {mathbb E}
 \left(1+\operatorname {Var}_{\mu_1}
                   \langle X,w_1\rangle\right)^{-1/2}\ge c.
 }                                                      \tag{0.9}
\]

No estimate in this note proves (0.9).  Section 7 gives exact statewise
and adaptive-mixture countermodels showing that (0.9) cannot follow from
curvature, mass preservation, and total covariance alone.

## 1. Common stochastic-localization identities

Let \(p_0=e^{-V_0}\) be a positive smooth log-concave probability density
on \(\mathbb R^n\), initially with the moments used below.  For a
predictable positive-semidefinite matrix \(D_t\), choose a predictable
square root \(C_t\) satisfying \(C_tC_t^T=D_t\), and define

\[
 p_t(x)={1\over Z_t}
 \exp\left\{c_t\cdot x-\frac12x^T\mathsf B_tx\right\}p_0(x), \tag{1.1}
\]

\[
 dc_t=C_t\,dW_t+D_ta_t\,dt,qquad
 d\mathsf B_t=D_t\,dt.                                \tag{1.2}
\]

The normalization correction in (1.2) gives

\[
 dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle.         \tag{1.3}
\]

Thus, for every fixed integrable scalar or vector function \(F\),

\[
 d{\mathbb E}_tF(X)
   =\operatorname {Cov}_t(F(X),X)^TC_t\,dW_t.          \tag{1.4}
\]

Put \(g_t={\mathbb E}_t{\bf1}_E\), \(Y=X-a_t\),
\(h={\bf1}_E-g_t\), and

\[
 H_t={\mathbb E}_t[hYY^T],\qquad
 \mathcal T_t(z)={\mathbb E}_t[YY^T\langle Y,z\rangle]. \tag{1.5}
\]

Writing \(T_{t,k}=\mathcal T_t(C_te_k)\), direct Itô calculation gives

\[
 \boxed{
 \begin{aligned}
 da_t&=A_tC_t\,dW_t,\\
 dg_t&=b_t^TC_t\,dW_t,\\
 db_t&=H_tC_t\,dW_t-A_tD_tb_t\,dt,\\
 dA_t&=\sum_kT_{t,k}\,dW_{t,k}-A_tD_tA_t\,dt.
 \end{aligned}}                                      \tag{1.6}
\]

The drift in \(db_t\) is the cross variation between \(g_t\) and
\(a_t\); it cannot be discarded for a control that does not annihilate
\(b_t\).

Every assertion below is first made before a bounded stopping time on
which \(A_t\succ0\), the relevant signal is separated from zero, the
natural parameters and barycentre are bounded, and all moments in (1.6)
are bounded.  On such a region posterior moments are smooth functions of
\((c,\mathsf B)\).  The matrix square roots used below have constant rank
and a spectral gap at zero, so their principal square roots are locally
Lipschitz.  Standard finite-dimensional SDE theory then gives a unique
strong solution up to the stop.  This is the exact level of
well-posedness used in the identities.

If the stopped coefficient is frozen after the stopping time, (1.3)
remains a true martingale after the usual integrability localization.
For isotropic \(p_0\), the law-of-total-covariance identity is

\[
 \boxed{\qquad
 {\mathbb E}A_t+\operatorname {Cov}(a_t)=I,qquad
 {\mathbb E}A_t\preceq I.
 \qquad}                                               \tag{1.7}
\]

The direction in (1.7) must be deterministic.  It gives no bound on
\({\mathbb E}\sup_{|v|=1}v^TA_tv\), nor on a direction selected from the
same posterior.

## 2. The covariance-normalized control

Assume \(A\succ0\) and \(q=A^{-1/2}b\ne0\), and use (0.1).  Since
\(Pq=0\),

\[
 D b=A^{-1/2}Pq=0,
\]

and therefore \(C^Tb=0\).  Equations (1.6) become

\[
 \boxed{
 \begin{aligned}
 dg_t&=0,\\
 db_t&=H_tC_t\,dW_t,\\
 dA_t&=\sum_kT_{t,k}\,dW_{t,k}
       -A_t^{1/2}P_tA_t^{1/2}\,dt,\\
 d\mathsf B_t&=A_t^{-1/2}P_tA_t^{-1/2}\,dt.
 \end{aligned}}                                      \tag{2.1}
\]

This verifies pathwise set-mass preservation, the requested covariance
drift, and the quadratic-tilt equation.

### 2.1 Trace, effective rank, and clocks

Let

\[
 K_t=A_t^{1/2}P_tA_t^{1/2}.
\]

Then \(K_t\) has rank \(n-1\) and

\[
 d\operatorname {tr}A_t
 =m_{3,t}^TC_t\,dW_t
  -\bigl(\operatorname {tr}A_t-u_t^TA_tu_t\bigr)dt,   \tag{2.2}
\]

where \(u_t=q_t/|q_t|\) and
\(m_{3,t}={\mathbb E}_t[|Y|^2Y]\).  In particular,

\[
 {\mathbb E}\int_0^t\operatorname {tr}K_sds
 \le\operatorname {tr}A_0.                            \tag{2.3}
\]

The Euclidean effective rank of the covariance drift satisfies

\[
 {\operatorname {tr}(K_t)^2\over\operatorname {tr}(K_t^2)}
 \ge {\operatorname {tr}K_t\over\|A_t\|_{\rm op}},   \tag{2.4}
\]

and can be as small as one in an anisotropic state.  In the current
covariance metric, however,

\[
 A_t^{-1/2}K_tA_t^{-1/2}=P_t,                          \tag{2.5}
\]

so its normalized rank and effective rank are exactly \(n-1\).

The mean and log-density clocks are

\[
 d[a]_t=K_tdt,
 \qquad
 {d\over dt}{\mathbb E}_t
 \bigl[\langle Y,C_t\,dW_t\rangle^2\bigr]/dt
 =\operatorname {tr}(D_tA_t)=\operatorname {tr}P_t=n-1. \tag{2.6}
\]

Thus the normalized control has no inverse-covariance singularity in the
posterior-standardized likelihood clock.  Its natural-parameter clock is

\[
 \operatorname {tr}D_t
 =\operatorname {tr}(P_tA_t^{-1})
 =\operatorname {tr}A_t^{-1}-u_t^TA_t^{-1}u_t,        \tag{2.7}
\]

which can diverge only as covariance directions collapse.

There is also an exact logarithmic-volume identity.  Let

\[
 R_{t,k}=A_t^{-1/2}T_{t,k}A_t^{-1/2}.
\]

Matrix Itô calculus gives

\[
\begin{aligned}
 d\log\det A_t
  ={}&\sum_k\operatorname {tr}(R_{t,k})\,dW_{t,k}
       -(n-1)dt\\
 &-{1\over2}\sum_k\operatorname {tr}(R_{t,k}^2)dt.
\end{aligned}                                        \tag{2.8}
\]

The fixed drift \(-(n-1)\) follows from
\(\operatorname {tr}(D_tA_t)=n-1\).  The relevant standardized bounds
are uniform over the posterior state.  Put \(Z=A_t^{-1/2}Y\) and
\(r_{t,k}=A_t^{1/2}C_te_k\).  Then \(Z\) is isotropic log-concave and

\[
 \sum_kr_{t,k}r_{t,k}^T=P_t,\qquad
 R_{t,k}={\mathbb E}_t[ZZ^T\langle Z,r_{t,k}\rangle].  \tag{2.9}
\]

Every one-dimensional isotropic log-concave marginal has universally
bounded third absolute moment.  Hölder's inequality therefore gives,
with a numerical \(C\),

\[
 \sum_k\|R_{t,k}\|_{HS}^2\le Cn^3,\qquad
 \sum_k\operatorname {tr}(R_{t,k})^2\le Cn^3.          \tag{2.10}
\]

For the first inequality, write the third-moment tensor as
\(M_{ij\ell}={\mathbb E}Z_iZ_jZ_\ell\) and use

\[
 \sum_{i,j,k}\langle M_{ij,\cdot},r_k\rangle^2
 \le\sum_{i,j,\ell}M_{ij\ell}^2\le Cn^3.
\]

For the second, use
\(\sum_k\operatorname {tr}(R_k)^2
\le|{\mathbb E}(|Z|^2Z)|^2\le Cn^3\).
Consequently the drift and quadratic variation in (2.8) are bounded on
every finite time interval by constants depending only on \(n\), not on
the posterior state.  Its right side cannot tend to \(-\infty\) at a
finite time.

This rules out a finite covariance face.  Here is the continuation
argument explicitly.  The trace is a nonnegative local supermartingale
by (2.2), so its stopped maximal inequality rules out
\(\sup_{s\le T}\operatorname {tr}A_s=\infty\).  Also,

\[
 b_tb_t^T\preceq g(1-g)A_t,\qquad |q_t|^2\le g(1-g),   \tag{2.11}
\]

by Cauchy--Schwarz.  On an interval on which \(|q_t|\ge\delta>0\),
(2.8)--(2.10) keep \(\lambda_{\min}(A_t)>0\), the trace bound keeps
\(\lambda_{\max}(A_t)<\infty\), and the mean martingale has finite
quadratic variation because \(d[a]_t=K_tdt\preceq A_tdt\).
Log-concavity then bounds every required centred moment in terms of
\(A_t\), while Section 4.1 rules out loss of normalization through an
exponential tail.  The natural-parameter coefficients remain locally
bounded and the stopped solutions concatenate.  Hence, for a smooth
full-dimensional initial law, the normalized process has a unique strong
continuation through every finite time before \(q_t\) hits zero.  The
only hard finite-time obstruction left is the zero-signal singularity.

### 2.2 Spectrum of the accumulated quadratic tilt

Assume \(n\ge2\); in dimension one the projector is zero whenever the
signal is nonzero, so there is only the survivor direction.

For eigenvalues in increasing order, rank-one interlacing applied to

\[
 D_t=A_t^{-1}-(A_t^{-1/2}u_t)(A_t^{-1/2}u_t)^T
\]

gives

\[
 \lambda_1(D_t)=0,qquad
 \lambda_2(D_t)\ge\lambda_1(A_t^{-1})
                 ={1\over\|A_t\|_{\rm op}}.           \tag{2.12}
\]

The sum of the two smallest eigenvalues is superadditive on the
positive-semidefinite cone.  Equivalently, by the Ky Fan variational
formula,

\[
\begin{aligned}
 \beta_1(t)+\beta_2(t)
 &=\min_{U^TU=I_2}\operatorname {tr}(U^T\mathsf B_tU)\\
 &\ge\int_0^t\bigl(\lambda_1(D_s)+\lambda_2(D_s)\bigr)ds\\
 &\ge\int_0^t{ds\over\|A_s\|_{\rm op}}.
\end{aligned}                                        \tag{2.13}
\]

Since \(\beta_2\ge(\beta_1+\beta_2)/2\), (0.4) follows.

Suppose on a path that \(\sup_{s\ge0}\|A_s\|_{\rm op}\le R<\infty\).
Then \(\beta_2(t)\ge t/(2R)\).  The posterior potential satisfies

\[
 D^2V_t\succeq\mathsf B_t
 \succeq\beta_2(t)(I-w_tw_t^T),                        \tag{2.14}
\]

where \(w_t\) is a weakest eigenvector.  Brascamp--Lieb, with an
arbitrarily small regularization in the weak direction, implies

\[
 \operatorname {Var}_t\langle X,z\rangle
 \le {1\over\beta_2(t)}\qquad(z\perp w_t,\ |z|=1).    \tag{2.15}
\]

For fixed \(n\), any tight subsequential weak limit is therefore
supported on an affine line.  If the weak direction rotates without a
limit, every direction may collapse; there still cannot be two
uncollapsed directions along a convergent subsequence.  This is the
precise meaning of “at most one weak direction.”

The hypothesis \(\sup_s\|A_s\|<\infty\) cannot simply be deleted.  From
(1.7),

\[
 {\mathbb E}\|A_t\|_{\rm op}\le n,
 \qquad
 {\mathbb E}{1\over\|A_t\|_{\rm op}}\ge {1\over n},   \tag{2.16}
\]

which only suggests time of order \(n\).  Likewise the nonnegative
trace supermartingale gives
\(\mathbb P(\sup_s\operatorname {tr}A_s\ge Rn)\le R^{-1}\),
again with an \(n\)-scale survivor bound.

### 2.3 Gaussian halfspace

Let \(p_0=N(0,I_n)\) and \(E=\{x_1\ge0\}\).  The protected direction is
fixed at \(e_1\), third centred moments vanish, and

\[
 \boxed{
 \begin{aligned}
 A_t&=e_1e_1^T+e^{-t}(I-e_1e_1^T),\\
 D_t&=e^t(I-e_1e_1^T),\\
 \mathsf B_t&=(e^t-1)(I-e_1e_1^T).
 \end{aligned}}                                      \tag{2.17}
\]

The transverse mean coefficient is \(e^{-t/2}\), and

\[
 \int_0^T\operatorname {tr}(D_tA_t)dt=(n-1)T,
 \qquad
 \int_0^T\operatorname {tr}D_tdt=(n-1)(e^T-1).       \tag{2.18}
\]

Both are finite for every finite \(T\); the posterior converges weakly
to the surviving Gaussian line only as \(T\to\infty\).  Its survivor
variance is exactly one.  Thus the normalized driver fixes the earlier
finite-time singular-clock defect in the best-aligned model.

## 3. The zero-signal convention is not well posed

The declaration \(P=I\) at \(q=0\) is not a removable convention.
Consider \(n=1\), \(p_0=N(0,1)\), and

\[
 E=[-r,r],\qquad 0<\mu_0(E)<1.
\]

For the one-dimensional Gaussian posterior near the symmetric state,
\(b=0\) exactly when its linear natural parameter \(c=0\), while

\[
 H={\mathbb E}[({\bf1}_E-g)X^2]\ne0.                  \tag{3.1}
\]

In dimension one, \(P=0\) whenever \(q\ne0\), and the proposed
convention gives \(P=1\) when \(q=0\).  For the normalized driver a
continuous weak solution would therefore satisfy, locally,

\[
 d[c]_t=A_t^{-1}{\bf1}_{\{c_t=0\}}dt.                 \tag{3.2}
\]

On a bounded stop, \(A_t^{-1}\) is bounded above and below by positive
constants.  The occupation-density formula for a continuous
semimartingale gives

\[
 \int_0^T{\bf1}_{\{c_t=0\}}d[c]_t=0.                  \tag{3.3}
\]

Equations (3.2)--(3.3) force \([c]_T=0\).  Hence \(c\equiv0\), but then
(3.2) forces strictly positive quadratic variation.  This is a
contradiction.  The bounded control of Section 4 gives the identical
contradiction with \(A_t^{-1}\) replaced by one.

A centred radial set under any rotationally invariant law has the same
zero-signal issue: \(b_0=0\), while for a nontrivial centred ball the
matrix \(H_0\) is generally a nonzero scalar multiple of \(I\).  The
full driver attempts to create a direction, and the hard projection
deletes that direction immediately.  In higher dimension a relaxed
weak solution may exist, but it is selector-dependent; the one-
dimensional example already disproves a general hard-convention theorem.

One may instead freeze a rank-\((n-1)\) projection at zero, soften the
projector, or use a vanishing mesh.  A frozen selector is noncanonical,
a soft projector has \(D_tb_t\ne0\) and hence loses exact mass
preservation, and a mesh limit requires a separate martingale-problem
identification.  None of these choices supplies the survivor estimate.

## 4. The bounded Euclidean projection driver

Assume \(b_t\ne0\) and use (0.5).  Because \(P_tb_t=0\), (1.6) gives

\[
 \boxed{
 \begin{aligned}
 dg_t&=0,\\
 db_t&=H_tP_t\,dW_t,\\
 dA_t&=\sum_k\mathcal T_t(P_te_k)dW_{t,k}
        -A_tP_tA_tdt,\\
 d\mathsf B_t&=P_tdt.
 \end{aligned}}                                      \tag{4.1}
\]

The coefficients \(C_t=P_t\) and \(D_t=P_t\) are bounded.  There is no
inverse covariance and no covariance-face clock.  The natural-parameter
quadratic variation is exactly

\[
 \operatorname {tr}D_t=n-1.                           \tag{4.2}
\]

The posterior likelihood and mean clocks are respectively

\[
 \operatorname {tr}(P_tA_t),
 \qquad \operatorname {tr}(A_tP_tA_t).                \tag{4.3}
\]

They may be large if an adaptively chosen transverse covariance is
large, but the SDE coefficient itself remains bounded in natural
coordinates.

### 4.1 Exponential tails do not create a hidden finite-time boundary

There is a useful deterministic integrability fact.  For any continuous
path of projectors \(P_s\),

\[
 \mathsf B_t=\int_0^tP_sds\succeq0.                   \tag{4.4}
\]

If \(\mathsf B_t\succ0\), the factor
\(\exp(c_t\cdot x-x^T\mathsf B_tx/2)\) is bounded above after completing
the square, so it is integrable against every probability \(p_0\), even
one with merely exponential tails.

If \(\mathsf B_t\) is singular, a vector \(w\in\ker\mathsf B_t\) obeys

\[
 0=w^T\mathsf B_tw=\int_0^t|P_sw|^2ds.                 \tag{4.5}
\]

Away from zero signal, \(P_s\) is continuous.  Hence \(P_sw=0\) for all
\(s\le t\), so all protected directions are the same line.  That line is
fixed by the deterministic initial signal.  Both terms in

\[
 c_t=\int_0^tP_s\,dW_s+\int_0^tP_sa_sds               \tag{4.6}
\]

are orthogonal to \(w\).  Thus \(c_t\in\operatorname {Ran}\mathsf B_t\),
and completing the square on that range again bounds the likelihood
factor above.  Therefore rank deficiency of the quadratic tilt does not
drive a one-sided exponential marginal past its integrability boundary.
The same argument applies to the normalized driver with \(P_s\) replaced
by \(D_s\), before its coefficient stops.

This proves absence of a tail-normalization obstruction on every
zero-signal-free stopped interval.  It does not repair the discontinuity
at \(b=0\), nor does it replace moment/nonexplosion estimates needed to
remove all auxiliary stops.

### 4.2 Deterministic rank-one-defect curvature

Again assume \(n\ge2\).

This subsection is a deterministic matrix theorem for an arbitrary
measurable unit-vector path \(u_s\).  It does not assert that the hard
stochastic controller continues through \(b_s=0\); Section 3 proves that
the stipulated full-driver continuation can fail to exist.

Put

\[
 Q_t=\int_0^tu_su_s^Tds,
 \qquad \mathsf B_t=tI-Q_t.                            \tag{4.7}
\]

Let \(\alpha_1\ge\alpha_2\ge\cdots\ge0\) be the eigenvalues of \(Q_t\).
Since \(\operatorname {tr}Q_t=t\),

\[
 \alpha_2\le {t\over2}.                               \tag{4.8}
\]

The two smallest eigenvalues of \(\mathsf B_t\) are
\(t-\alpha_1\) and \(t-\alpha_2\).  Therefore

\[
 \lambda_2(\mathsf B_t)=t-\alpha_2\ge {t\over2}.      \tag{4.9}
\]

If \(w_t\) is a top eigenvector of \(Q_t\), then

\[
 D^2V_t\succeq {t\over2}(I-w_tw_t^T).                 \tag{4.10}
\]

The rank-one-defect Lichnerowicz lemma, in its distributional form,
states that a log-concave probability satisfying
\(D^2V\succeq\kappa(I-ww^T)\) obeys

\[
 C_P(\mu)\le C\left(\kappa^{-1}
       +\operatorname {Var}_\mu\langle X,w\rangle\right). \tag{4.11}
\]

Applying it to (4.10) gives

\[
 C_P(\mu_t)\le C\left({2\over t}
       +V_t\right),\qquad
 V_t=\operatorname {Var}_{\mu_t}\langle X,w_t\rangle. \tag{4.12}
\]

If the zero-signal-free SDE exists through \(t=1\), all curvature and
clock terms are universal.  Only \(V_1\) remains.

### 4.3 Exact perimeter averaging and the survivor gate

Let the stopped process be defined through time \(t\), and suppose
\(g=\mu_0(E)\in(0,1)\).  For every \(\varepsilon>0\), the density
martingale gives

\[
 {\mathbb E}\mu_t(E_\varepsilon)=\mu_0(E_\varepsilon),
 \qquad \mu_t(E)=g\quad\hbox{pathwise}.                \tag{4.13}
\]

Fatou's lemma, applied to the nonnegative exterior difference quotients,
therefore yields

\[
 \mu_0^+(E)\ge {\mathbb E}\mu_t^+(E).                 \tag{4.14}
\]

Buser--Ledoux for log-concave measures and (4.12) give

\[
 \boxed{
 \mu_0^+(E)\ge c\min(g,1-g)\,
 {\mathbb E}\left({2\over t}+V_t\right)^{-1/2}.
 }                                                      \tag{4.15}
\]

Since \(x\mapsto(2/t+x)^{-1/2}\) is convex, a universal estimate
\({\mathbb E}V_1\le C\) would close this route.  The total covariance
identity only gives

\[
 {\mathbb E}V_1
 ={\mathbb E}\,w_1^TA_1w_1
 \le {\mathbb E}\operatorname {tr}A_1\le n.          \tag{4.16}
\]

Replacing the last \(n\) by a universal constant is the adaptive-
survivor problem.

## 5. Aligned model audits

### 5.1 Gaussian halfspace

For \(N(0,I_n)\) and \(E=\{x_1\ge0\}\), the bounded control has fixed
\(u=e_1\).  The exact solution is

\[
 \mathsf B_t=t(I-e_1e_1^T),\qquad
 A_t=e_1e_1^T+{1\over1+t}(I-e_1e_1^T).                 \tag{5.1}
\]

The survivor variance is one.  At \(t=1\), (4.12) is universal.  The
normalized solution was recorded in (2.14).  Both controls preserve the
Gaussian halfspace perimeter exactly.

### 5.2 Product asymmetric exponentials: coordinate cap

Let \(X_j=Y_j-1\), where the \(Y_j\) are independent
\(\operatorname {Exp}(1)\), and take a half-mass set depending only on
\(X_1\).  For both controls the protected direction remains \(e_1\), the
first marginal is unchanged, and its variance remains one.  Under the
bounded control each transverse factor satisfies

\[
 dA_j=m_{3,j}\,dW_j-A_j^2dt,
 \qquad dA_j\big|_{t=0}=2\,dW_j-dt,                    \tag{5.2}
\]

while under the normalized control

\[
 dA_j={m_{3,j}\over\sqrt{A_j}}dW_j-A_jdt,
 \qquad dA_j|_{t=0}=2\,dW_j-dt.                       \tag{5.3}
\]

Thus even the aligned exponential model has order-one covariance noise,
but its protected variance is benign because the controller never tilts
that factor.

### 5.3 Off-centre and centred radial sets

At a standard Gaussian state, let \(E=B(z,R)\) with \(d=|z|>0\).  If
\(g(a)=G(|z-a|)\) denotes its mass under \(N(a,I)\), then

\[
 b=\nabla_ag,qquad H=\nabla_a^2g,qquad \mathcal T=0.
\]

The tangential eigenvalue of \(H\) is \(G'(d)/d\), and \(|b|=|G'(d)|\).
At \(A=I\), both controls therefore have initial angular rate

\[
 {d\over dt}\operatorname {tr}[u]_0={n-1\over d^2}.  \tag{5.4}
\]

The protected line genuinely rotates.  Gaussian covariance nevertheless
remains bounded by \(I\), so survivor variance is not the obstruction in
this model.  For a centred radial ball, \(b=0\), and the hard process is
instead blocked by Section 3.

## 6. Statewise protected-variance countermodel

The rank-one curvature lemma cannot control the survivor variance from
the current state alone.  This is already exact for a two-dimensional
product exponential.

Let \(Y_1\sim\operatorname {Exp}(\lambda)\), with \(0<\lambda<1\), and
let the second coordinate be any one-dimensional log-concave law whose
potential has second derivative at least \(\kappa>0\).  Their product is
a natural-parameter posterior of the isotropic shifted-exponential
product: relative to \(\operatorname {Exp}(1)\), the first coordinate has
linear tilt \(1-\lambda\), and the second can be obtained by a quadratic
tilt.  Put

\[
 E_\lambda=\left\{Y_1\ge{\log2\over\lambda}\right\}.  \tag{6.1}
\]

Then \(\mu(E_\lambda)=1/2\),

\[
 b={\log2\over2\lambda}e_1,
 \qquad \operatorname {Var}(Y_1)=\lambda^{-2}.         \tag{6.2}
\]

Both the bounded and normalized controllers protect \(e_1\) at this
state.  The posterior has transverse curvature \(\kappa\), but

\[
 C_P(\mu)\ge\operatorname {Var}(Y_1)=\lambda^{-2}.     \tag{6.3}
\]

Letting \(\lambda\downarrow0\) proves that no statewise combination of
mass preservation, transverse curvature, and log-concavity bounds the
survivor.  This state is admissible in the natural-parameter family.  The
calculation does not claim that the controlled flow from the isotropic
initial state reaches it with enough probability.  That probability is
exactly the missing mixture-average question.

## 7. Sharp adaptive-mixture no-go for covariance bookkeeping

There is a finite-dimensional endpoint model showing why (1.7) cannot
control an adaptive survivor.  Let \(I\) be uniform on
\(\{1,\ldots,n\}\), and conditionally on \(I=i\) let \(\nu_i\) be a
one-dimensional log-concave probability on the line \(\mathbb Re_i\),
centred and with variance \(n\).  Then

\[
 A_i=ne_ie_i^T,qquad
 {\mathbb E}_I A_I=I.                                 \tag{7.1}
\]

Give component \(i\) the quadratic curvature matrix

\[
 \mathsf B_i=I-e_ie_i^T.                              \tag{7.2}
\]

It has curvature one on the complement of its adaptive weak direction
\(w_i=e_i\), yet

\[
 w_i^TA_iw_i=n,qquad
 {\mathbb E}_I(1+w_i^TA_iw_i)^{-1/2}=(n+1)^{-1/2}.    \tag{7.3}
\]

Every component in this endpoint ensemble is log-concave on its affine
support, and (7.1) saturates the total-covariance information used in the
localization argument.  The mixture \(n^{-1}\sum_i\nu_i\) is not
log-concave, so (7.1)--(7.3) are not a counterexample to KLS and are not
asserted to be the endpoint law of the proposed SDE.  They are a sharp
countermodel to the inference

\[
 {\mathbb E}A\preceq I
 \quad+\quad\hbox{rank-one-defect curvature}
 \quad\Longrightarrow\quad
 {\mathbb E}w^TAw=O(1).                               \tag{7.4}
\]

A successful proof must use compatibility of the posterior ensemble with
one common initial log-concave density and with the mass-preserving
martingale, not merely covariance bookkeeping.

## 8. Product-exponential maximum: exact adaptive growth

Let \(Y_1,\ldots,Y_n\) be independent \(\operatorname {Exp}(1)\) variables,
put \(X_i=Y_i-1\), and consider

\[
 E=\{\max_iY_i\le R\},\qquad
 a=e^{-R},\qquad F=1-a,qquad g=F^n.                  \tag{8.1}
\]

At the isotropic initial state, set

\[
 s=F^{n-1}aR.
\]

Then

\[
 b=-s(1,\ldots,1),\qquad |b|=\sqrt n\,s,qquad
 u=-{1\over\sqrt n}(1,\ldots,1).                     \tag{8.2}
\]

For every \(z\perp(1,\ldots,1)\), direct product integration gives

\[
 Hz=-{sR\over F}z,
 \qquad \mathcal T(z)b=-2s z.                         \tag{8.3}
\]

### 8.1 Angular motion

For the bounded controller, \(db=HP\,dW\) at time zero.  Therefore

\[
 \boxed{
 {d\over dt}\operatorname {tr}[u]_0
 ={1\over|b|^2}\operatorname {tr}(PH^2P)
 ={n-1\over n}\left({R\over F}\right)^2.
 }                                                      \tag{8.4}
\]

For the normalized controller, the martingale differential of
\(q=A^{-1/2}b\) at \(A=I\) in a transverse direction is

\[
 H z-{1\over2}\mathcal T(z)b
 =-s\left({R\over F}-1\right)z.
\]

Hence its whitened protected direction has initial angular rate

\[
 {n-1\over n}\left({R\over F}-1\right)^2.             \tag{8.5}
\]

For the median maximum, \(F^n=1/2\), so

\[
 F=2^{-1/n},\qquad
 R=\log n-\log\log2+o(1),                              \tag{8.6}
\]

and both (8.4) and (8.5) are \(\Theta((\log n)^2)\).

### 8.2 Growth of the protected variance

Let

\[
 U_t=u_t^TA_tu_t,
 \qquad u_t={b_t\over|b_t|}
\]

for the bounded controller.  At \(A_0=I\), all Itô terms involving only
the normalization of \(u_t\) cancel because \(|u_t|=1\).  The covariance
drift contributes \(-u^TPu=0\).  The only expectation drift is the cross
variation between the martingale parts of \(u_t\) and \(A_t\):

\[
 {d\over dt}{\mathbb E}U_t\big|_{t=0}
 ={2\over|b|}\sum_k
 \left\langle PHPe_k,\mathcal T(Pe_k)u\right\rangle.  \tag{8.7}
\]

Using (8.2)--(8.3),

\[
 \boxed{
 {d\over dt}{\mathbb E}U_t\big|_{t=0}
 ={4(n-1)R\over nF}=\Theta(\log n).
 }                                                      \tag{8.8}
\]

Thus the adaptive protected variance has an unbounded positive initial
drift even though \({\mathbb E}A_t\preceq I\) in every deterministic
direction.

For the literal signal-weighted energy \(G_t=b_t^TA_tb_t\), Itô's formula
instead gives

\[
 \boxed{
 {d\over dt}{\mathbb E}G_t\big|_{t=0}
 =(n-1)s^2\left[\left({R\over F}\right)^2
                    +{4R\over F}\right].
 }                                                      \tag{8.9}
\]

Since \(G_0=ns^2\), its relative initial growth rate is
\(\Theta((\log n)^2)\).  The absolute derivative in (8.9) tends to zero
for the median maximum, so (8.9) alone is not a KLS counterexample.  The
unit-direction quantity (8.8) is the one appearing in the survivor
variance gate.

These formulas do not prove that \({\mathbb E}U_1\) diverges: the angular
time scale is \(1/\log^2n\), on which the change predicted by (8.8) is
only \(O(1/\log n)\).  They do prove that a universal differential bound
or a fixed-direction use of (1.7) cannot control the process.  Any
successful survivor estimate needs a compensating angular/covariance
identity beyond the present localization equations.

## 9. What is proved and what remains

The covariance-normalized control proves, on zero-signal-free stopped
intervals:

1. exact set-mass preservation;
2. covariance drift \(-A^{1/2}PA^{1/2}\);
3. quadratic tilt derivative \(\mathsf B'=A^{-1/2}PA^{-1/2}\);
4. constant standardized likelihood clock \(n-1\);
5. exponential Gaussian collapse with no finite-time covariance clock;
6. the curvature estimate (0.4), and hence at most one uncollapsed
   direction along every bounded-covariance tight limit.

The bounded Euclidean projection driver proves, on the same kind of
stopped interval:

1. exact set-mass preservation with bounded natural coefficients;
2. no exponential-tail integrability obstruction away from zero signal;
3. deterministic transverse curvature \(t/2\) by time \(t\);
4. the exact perimeter estimate (4.15).

Neither control proves a dimension-free survivor variance.  The hard
zero-signal convention is not even weakly well posed, the product maximum
has rapid adaptive rotation and positive protected-variance drift, and
the statewise and endpoint models in Sections 6--7 show why local
curvature plus total covariance cannot close the gate.  A new theorem
must control the compatibility of the random weak direction with the
posterior mixture; assuming that theorem would be a KLS-strength step.
