# Function-specific stochastic localization: exact ledgers and the adaptive survivor

## Executive conclusion

Let \(\mu\) be isotropic and log-concave, and let \(f\) be a normalized
first eigenfunction:
\[
 \int f\,d\mu=0,\qquad \int f^2\,d\mu=1,
 \qquad \int|\nabla f|^2\,d\mu=\lambda.                  \tag{0.1}
\]
For a stochastic-localization posterior \(\mu_t\), put
\[
 m_t=\int f\,d\mu_t,qquad
 b_t=\operatorname{Cov}_{\mu_t}(f,X),qquad
 \mathcal E_t=\int|\nabla f|^2\,d\mu_t.                 \tag{0.2}
\]
The proposed hard driver
\[
 C_t=P_{b_t^\perp}                                       \tag{0.3}
\]
has a precise benefit and a precise cost.

* Since \(C_tb_t=0\), \(m_t=0\) pathwise.
* Consequently it does not dissipate any of the variance of \(f\):
  \[
   \mathbb E\operatorname{Var}_{\mu_T}(f)=1.
  \]
* The localized Dirichlet energy is also conserved in expectation:
  \[
   \mathbb E\mathcal E_T=\lambda.
  \]
* The accumulated quadratic curvature
  \(B_T=\int_0^TC_t^2dt\) has at most one weak direction.  There is a
  path-dependent unit vector \(v_T\) such that
  \[
   B_T\succeq\frac T2P_{v_T^\perp}.
  \]

Let
\[
 R_T=\operatorname{Var}_{\mu_T}\langle X,v_T\rangle.
\]
The rank-one-defect theorem with its audited constant \(96\) gives the
pathwise estimate
\[
 \operatorname{Var}_{\mu_T}(f)
 \le96\left(\frac2T+R_T\right)\mathcal E_T.              \tag{0.4}
\]
Averaging and using the exact ledgers yields
\[
 \boxed{
 1\le\frac{192}{T}\lambda
       +96\,\mathbb E[R_T\mathcal E_T].
 }                                                        \tag{0.5}
\]
In particular,
\[
 \boxed{
 \mathbb E[R_T\mathcal E_T]
 \ge\frac1{96}-\frac{2\lambda}{T}.
 }                                                        \tag{0.6}
\]
If \(\lambda\le T/384\), the right side is at least \(1/192\).
Equivalently, under the energy-biased path law
\[
 \frac{d\mathbb Q_T}{d\mathbb P}=\frac{\mathcal E_T}{\lambda},
\]
one necessarily has
\[
 \mathbb E_{\mathbb Q_T}R_T
 \ge\frac1{96\lambda}-\frac2T
 \ge\frac1{192\lambda}
 \quad\text{if }\lambda\le T/384.                       \tag{0.7}
\]

Thus the function-specific control does not remove the adaptive covariance
obstruction.  It makes the obstruction quantitatively unavoidable: a
hypothetical small-gap first eigenfunction forces the surviving posterior
direction to have variance of order \(1/\lambda\) under exactly the path
bias relevant to its energy.

For example, any estimate
\[
 \mathbb E[R_T\mathcal E_T]\le K\lambda                 \tag{0.8}
\]
with numerical \(K\) would imply
\[
 \lambda\ge\frac1{96K+192/T},                            \tag{0.9}
\]
and hence KLS.  The missing estimate is therefore at least KLS-strength.
It cannot be supplied merely by the rank/trace curvature calculation.

The soft control
\[
 C_t=I-\alpha P_{b_t},\qquad 0\le\alpha<1,               \tag{0.10}
\]
does not evade the obstruction.  Put
\(\delta=(1-\alpha)^2\).  It gains full curvature
\(B_T\succeq\delta TI\), but the posterior mean of \(f\) acquires
quadratic variation.  If
\[
 Q_T:=\mathbb E m_T^2,
\]
then the exact variance ledger and Brascamp--Lieb give
\[
 1=Q_T+\mathbb E\operatorname{Var}_{\mu_T}(f),
 \qquad
 1-Q_T\le\frac{\lambda}{\delta T}.                       \tag{0.11}
\]
Therefore keeping \(Q_T\le\eta<1\) forces
\[
 \delta T\le\frac{\lambda}{1-\eta}.                     \tag{0.12}
\]
Universal full curvature and approximate preservation of the fixed
eigenfunction are incompatible when \(\lambda\) is small.

The report below proves all of these claims, treats the zero-signal
well-posedness issue, gives exact Gaussian and product tests, and constructs
a smooth isotropic critical-tail posterior showing that rank-one terminal
curvature alone permits arbitrarily large survivor variance.  No genuine
small-gap isotropic eigenfunction counterexample is claimed: producing one
would itself disprove KLS.

## 1. Controlled localization and well-posedness

Let \(p_0\) be a smooth log-concave density on its affine support.  For a
bounded predictable symmetric positive-semidefinite control \(C_t\), define
the posterior
\[
 \frac{dp_t}{dp_0}(x)
 =Z_t^{-1}\exp\left(\langle c_t,x\rangle
              -\frac12\langle B_tx,x\rangle\right),      \tag{1.1}
\]
where
\[
 dc_t=C_t\,dW_t+C_t^2a_t\,dt,
 \qquad dB_t=C_t^2dt,
 \qquad a_t=\int x\,dp_t.                                \tag{1.2}
\]
After the usual moment and parameter stopping, Itô's formula gives
\[
 dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle.           \tag{1.3}
\]
Consequently, for every fixed integrable test \(g\),
\[
 d\mathbb E_tg
 =\operatorname{Cov}_t(g,X)^TC_t\,dW_t.                 \tag{1.4}
\]
The stopped density is a true measure-valued martingale:
\[
 \mathbb E\mu_t=\mu.                                     \tag{1.5}
\]
All identities below are first stated under bounded stopping and then pass
to any nonexplosive or relaxed limit for which the indicated quantities
are uniformly integrable.

### 1.1 The hard feedback at \(b=0\)

When \(b\ne0\), the orthogonal projection
\[
 C(b)=I-\frac{bb^T}{|b|^2}                               \tag{1.6}
\]
is smooth and annihilates \(b\).  It is discontinuous at \(b=0\).  A
Borel convention at zero does not by itself give a strong or weak solution;
the same occupation-time obstruction as in set-mass-preserving localization
appears already in dimension one.

There are two honest formulations.

1. Stop before \(|b_t|\) reaches a prescribed positive threshold.
2. Use mesh-relaxed limits.  On each mesh interval freeze a rank-
   \((n-1)\) projection annihilating the preceding \(b\); at \(b=0\),
   choose any rank-\((n-1)\) projection.  A stopped weak limit has a
   predictable symmetric control satisfying, almost everywhere,
   \[
    C_tb_t=0,\qquad 0\preceq C_t\preceq I,
    \qquad\operatorname{tr}(I-C_t^2)\le1.                \tag{1.7}
   \]

The conclusions below use only (1.7), so they apply uniformly to every
such relaxed limit.  They do not assert a canonical continuation through a
zero of \(b_t\).

A globally smooth alternative is
\[
 C_{\alpha,\varepsilon}(b)
 =I-\alpha\frac{bb^T}{|b|^2+\varepsilon},
 \qquad 0\le\alpha\le1.                                 \tag{1.8}
\]
Section 5 records its quantitative errors.

## 2. Exact fixed-function ledgers

For the fixed initial function \(f\), define
\[
 m_t=\mathbb E_tf,qquad
 V_t=\operatorname{Var}_t(f),qquad
 b_t=\operatorname{Cov}_t(f,X),qquad
 \mathcal E_t=\mathbb E_t|\nabla f|^2.                  \tag{2.1}
\]
Also put
\[
 r_t=\operatorname{Cov}_t((f-m_t)^2,X),
 \qquad
 d_t=\operatorname{Cov}_t(|\nabla f|^2,X).               \tag{2.2}
\]
Equation (1.4) and Itô's formula give
\[
 dm_t=b_t^TC_t\,dW_t,                                    \tag{2.3}
\]
\[
 dV_t=r_t^TC_t\,dW_t-|C_tb_t|^2dt,                      \tag{2.4}
\]
and
\[
 d\mathcal E_t=d_t^TC_t\,dW_t.                          \tag{2.5}
\]
The coefficient in (2.4) follows from
\[
 \operatorname{Cov}_t(f^2,X)-2m_tb_t
 =\operatorname{Cov}_t((f-m_t)^2,X).
\]

Optional sampling and (1.5) give the exact total-variance and energy
ledgers
\[
 \boxed{
 1=\mathbb E V_T+\mathbb E[m_T^2]
  =\mathbb E V_T+
    \mathbb E\int_0^T|C_tb_t|^2dt,
 }                                                        \tag{2.6}
\]
\[
 \boxed{
 \mathbb E\mathcal E_T=\lambda.
 }                                                        \tag{2.7}
\]

Under the hard condition \(C_tb_t=0\), these become
\[
 m_t=0\quad\text{pathwise},
 \qquad \mathbb EV_T=1,
 \qquad \mathbb E\mathcal E_T=\lambda.                 \tag{2.8}
\]
Thus preservation of the posterior mean is exactly preservation of all
unexplained variance in expectation.  In the filtering coupling,
\(m_t=\mathbb E[f(X)\mid\mathcal F_t]\); (2.8) says that the adaptive
observations acquire no conditional-mean information about \(f\).

### 2.1 Evolution of the protected direction

Let
\[
 A_t=\operatorname{Cov}_t(X),
 \qquad
 H_t=\mathbb E_t[(f-m_t)(X-a_t)(X-a_t)^T].                \tag{2.9}
\]
A direct product-rule calculation gives
\[
 db_t=H_tC_t\,dW_t-A_tC_t^2b_t\,dt.                     \tag{2.10}
\]
For the hard control, \(C_tb_t=0\) implies \(C_t^2b_t=0\), so
\[
 db_t=H_tC_t\,dW_t.                                      \tag{2.11}
\]
Away from \(b_t=0\), the angular quadratic variation of
\(u_t=b_t/|b_t|\) is therefore
\[
 d[u]_t^{\mathrm{angular}}
 =\frac{\|P_{u_t^\perp}H_tC_t\|_{\mathrm{HS}}^2}
        {|b_t|^2}\,dt.                                  \tag{2.12}
\]
Any proof controlling the final adaptive survivor from the special nature
of \(b_t\) would have to use information of this kind.  The variance and
energy ledgers alone do not control (2.12).

### 2.2 What remains of the eigenfunction equation

At time zero, testing the weak eigenfunction equation with the coordinate
functions gives
\[
 \int\nabla f\,d\mu=\lambda b_0,
 \qquad |b_0|\le1.                                      \tag{2.13}
\]
After localization, \(f\) is no longer an eigenfunction of the posterior
generator.  If
\[
 \ell_t(x)=\langle c_t,x\rangle-\frac12\langle B_tx,x\rangle
\]
is the log-likelihood up to normalization, then
\[
 \mathcal E_t
 =\lambda\,\mathbb E_t f^2
   +\mathbb E_t\!\left[f(B_tx-c_t)\cdot\nabla f\right]. \tag{2.14}
\]
For the hard flow \(m_t=0\), the first term is \(\lambda V_t\).  The
second term has no fixed sign, and its path average is zero by
(2.7)--(2.8).  Thus the initial eigenfunction equation supplies no
pathwise replacement for the missing survivor estimate.

## 3. Finite-time curvature and the exact rank-one reduction

Assume the hard relaxed conditions (1.7), and put
\[
 D_t=I-C_t^2\succeq0,
 \qquad U_T=\int_0^TD_tdt.
\]
Then
\[
 B_T=TI-U_T,
 \qquad\operatorname{tr}U_T\le T.                       \tag{3.1}
\]
Let \(v_T\) be a top eigenvector of \(U_T\).  Its second eigenvalue is at
most half its trace:
\[
 \lambda_2(U_T)\le\frac12\operatorname{tr}U_T\le\frac T2.
\]
Therefore
\[
 \boxed{
 B_T\succeq\frac T2P_{v_T^\perp}.
 }                                                        \tag{3.2}
\]
The posterior potential is the initial convex potential plus the affine
tilt and \(\langle B_Tx,x\rangle/2\), so distributionally
\[
 D^2V_T\succeq\frac T2P_{v_T^\perp}.                     \tag{3.3}
\]

Write
\[
 R_T=v_T^TA_Tv_T
 =\operatorname{Var}_{\mu_T}\langle X,v_T\rangle.        \tag{3.4}
\]
The audited rank-one-defect Lichnerowicz theorem says that a log-concave
law satisfying \(D^2V\succeq\kappa P_{v^\perp}\) obeys
\[
 C_P\le96(\kappa^{-1}+\operatorname{Var}\langle X,v\rangle).
                                                                    \tag{3.5}
\]
With \(\kappa=T/2\), this yields
\[
 C_P(\mu_T)\le96\left(\frac2T+R_T\right).                \tag{3.6}
\]
Applying the posterior Poincare inequality to the fixed \(f\) gives
\[
 V_T\le96\left(\frac2T+R_T\right)\mathcal E_T.           \tag{3.7}
\]
Taking expectations and inserting (2.8) proves (0.5)--(0.7).

### 3.1 Two useful biased formulations

Besides the energy-biased law in (0.7), define the variance-biased path
law
\[
 \frac{d\widehat{\mathbb Q}_T}{d\mathbb P}=V_T.
\]
This is a probability law because \(\mathbb EV_T=1\).  On \(V_T>0\),
(3.7) gives
\[
 \frac{\mathcal E_T}{V_T}
 \ge\frac1{96(2/T+R_T)}.
\]
Since
\[
 \mathbb E_{\widehat{\mathbb Q}_T}
 \frac{\mathcal E_T}{V_T}=\mathbb E\mathcal E_T=\lambda,
\]
one obtains
\[
 \boxed{
 \mathbb E_{\widehat{\mathbb Q}_T}
 \frac1{2/T+R_T}\le96\lambda.
 }                                                        \tag{3.8}
\]
Thus a small-gap eigenfunction forces large survivor variance both under
energy bias, through (0.7), and under variance bias, through the inverse
moment (3.8).

## 4. Why no covariance-free closure follows

### 4.1 The exact missing estimate

For a deterministic unit vector \(v\), the measure-martingale identity
does give
\[
 \mathbb E[v^TA_Tv]\le v^TA_0v=1.                        \tag{4.1}
\]
It does not apply to \(v_T\), which is selected from the same path through
\(U_T\).  Even a bound on \(\mathbb ER_T\) would not directly control
\(\mathbb E[R_T\mathcal E_T]\), because the localized energy is correlated
with the selected endpoint.

The exact sufficient estimate is
\[
 \mathbb E[R_T\mathcal E_T]\le K\lambda,                 \tag{4.2}
\]
or equivalently
\[
 \mathbb E_{\mathbb Q_T}R_T\le K.                        \tag{4.3}
\]
Equations (0.5) and (0.9) show that either statement, with numerical
\(K\) at one numerical time \(T\), proves KLS.  Conversely, (0.6) shows
that every hypothetical small-gap eigenfunction must violate (4.2) by a
factor of order \(1/\lambda\).  This is an exact reduction, not merely an
analogy with the usual covariance-process obstruction.

A trace estimate does not help.  Pathwise \(R_T\le\operatorname{tr}A_T\),
and \(\mathbb E\operatorname{tr}A_T\le n\), but neither controls the
energy-weighted product dimension freely.  Stopping when the trace is
bounded also stops curvature time; proving that a numerical amount of time
elapses before this stop is the standard covariance-process problem in a
different form.

### 4.2 The survivor term is deterministically necessary

There cannot be a terminal theorem which deletes \(R_T\) from (3.7).  For
\(R>0\) and \(\kappa>0\), consider the smooth log-concave Gaussian
\[
 \nu_{R,\kappa}
 =N(0,R)\otimes N(0,\kappa^{-1}I_{n-1}).                 \tag{4.4}
\]
It satisfies
\[
 D^2V\succeq\kappa P_{e_1^\perp}.
\]
For \(g(x)=x_1/\sqrt R\),
\[
 \operatorname{Var}_{\nu_{R,\kappa}}g=1,
 \qquad
 \int|\nabla g|^2d\nu_{R,\kappa}=\frac1R.              \tag{4.5}
\]
Thus no bound depending only on \(\kappa^{-1}\) can control the fixed
function variance.  The weak marginal variance is not a proof artifact.

The example (4.4) is a terminal geometric model, not an isotropic initial
state.  The next test shows that arbitrarily long weak marginals also occur
inside the exponential-quadratic posterior family of a smooth isotropic
log-concave initial law.

### 4.3 Smooth critical-tail endpoint test

Choose a smooth even convex potential \(V\) on \(\mathbb R\) such that the
probability \(d\nu\propto e^{-V(x)}dx\) is isotropic and
\[
 V'(x)\longrightarrow a>0\quad(x\to+\infty).
\]
A scaled version of \(V(x)=\sqrt{1+x^2}\) is a concrete choice.  For
\(0<\eta<a/2\), tilt by \(a-\eta\):
\[
 d\nu_\eta(x)
 =Z_\eta^{-1}e^{(a-\eta)x}\,d\nu(x).                    \tag{4.6}
\]
On the positive tail the linear part of \(V\) cancels, leaving decay
\(e^{-\eta x+o(1)}\).  Direct one-dimensional tail integration gives
numerical constants depending only on the fixed smoothing such that
\[
 c\eta^{-2}\le\operatorname{Var}_{\nu_\eta}X
 \le C\eta^{-2}.                                        \tag{4.7}
\]

Now start from the smooth isotropic product
\[
 \mu=\nu\otimes\gamma_{n-1}.
\]
The exponential-quadratic posterior with
\[
 c=(a-\eta)e_1,
 \qquad B=T P_{e_1^\perp}                                \tag{4.8}
\]
is
\[
 \mu_{c,B}=\nu_\eta\otimes
 N\left(0,(1+T)^{-1}I_{n-1}\right).                     \tag{4.9}
\]
It is smooth and log-concave, has the exact all-but-one curvature profile
produced by the trace calculation, and has survivor variance
\(R\asymp\eta^{-2}\), arbitrarily large.  The normalized weak coordinate
has variance one and energy \(\asymp1/R\), exactly as in (4.5).

This is a genuine smooth isotropic-initial terminal posterior model.  It is
not claimed to be the endpoint of the hard \(b_t\)-annihilating feedback:
if the protected direction were fixed at \(e_1\), its driver could not
generate the tilt in (4.8).  Its role is precise.  It proves that terminal
log-concavity, finite-time rank/trace curvature, and the initial isotropic
normalization do not exclude a long survivor.  Any exclusion must use the
coupled evolution of \(b_t\), especially its rotation (2.12).

Producing the same endpoint together with a normalized small-gap initial
first eigenfunction and the exact hard feedback would be a counterexample
to KLS itself.  No such claim is made.

## 5. The soft tradeoff \(C=I-\alpha P_b\)

Assume first that \(b_t\ne0\), take a constant \(0\le\alpha<1\), and set
\[
 u_t=\frac{b_t}{|b_t|},
 \qquad P_t=u_tu_t^T,
 \qquad C_t=I-\alpha P_t.                               \tag{5.1}
\]
Define
\[
 \delta=(1-\alpha)^2,
 \qquad \beta=2\alpha-\alpha^2=1-\delta.                \tag{5.2}
\]
Then
\[
 C_tb_t=(1-\alpha)b_t,
 \qquad C_t^2=I-\beta P_t.                              \tag{5.3}
\]
The exact variance ledger is
\[
 Q_T:=\mathbb E m_T^2
 =\delta\,\mathbb E\int_0^T|b_t|^2dt,
 \qquad
 \mathbb EV_T=1-Q_T.                                    \tag{5.4}
\]
The accumulated curvature is
\[
 B_T=TI-\beta\int_0^TP_tdt.                             \tag{5.5}
\]
Since \(\int_0^TP_tdt\preceq TI\),
\[
 \boxed{B_T\succeq\delta TI.}                          \tag{5.6}
\]
Brascamp--Lieb and (2.7) therefore give
\[
 1-Q_T=\mathbb EV_T
 \le\frac1{\delta T}\mathbb E\mathcal E_T
 =\frac\lambda{\delta T}.                              \tag{5.7}
\]
This proves (0.11)--(0.12).

There is simultaneously a slightly stronger transverse rank-one bound.
The second eigenvalue of \(\beta\int_0^TP_tdt\) is at most \(\beta T/2\),
so for a path-dependent \(v_T\),
\[
 B_T\succeq \frac{T(1+\delta)}2P_{v_T^\perp}.            \tag{5.8}
\]
This improves the transverse constant but does not remove the weak
direction; (5.6) controls that direction only by \(\delta T\).

The tradeoff is exact in the sense relevant to the proposal:

* If \(\delta T\) is numerical and \(\lambda\to0\), then (5.7) forces
  \(Q_T\to1\).  The localization observations learn almost all of \(f\)
  through the posterior mean.
* If one requires \(Q_T\le\eta<1\), then (0.12) forces the gained weak
  curvature to be only \(O(\lambda)\).  Using that curvature in
  Brascamp--Lieb is then tautological.

### 5.1 Smooth regularization at zero signal

For (1.8), the eigenvalue of \(C_{\alpha,\varepsilon}\) in the \(b\)
direction is at least \(1-\alpha\).  Hence for \(\alpha<1\),
\[
 B_T\succeq(1-\alpha)^2TI,                               \tag{5.9}
\]
and the same implication (5.7) holds with
\[
 Q_T=\mathbb E\int_0^T
       |C_{\alpha,\varepsilon}(b_t)b_t|^2dt.
\]
Near \(b=0\) the regularization leaks more information, not less, so it
cannot improve the preservation--curvature tradeoff.

For \(\alpha=1\),
\[
 C_{1,\varepsilon}(b)b
 =\frac{\varepsilon}{|b|^2+\varepsilon}b,
\]
and
\[
 |C_{1,\varepsilon}(b)b|^2
 \le\frac\varepsilon4.                                  \tag{5.10}
\]
Therefore
\[
 \mathbb Em_T^2\le\frac{\varepsilon T}{4}.              \tag{5.11}
\]
Also \(\operatorname{tr}(I-C_{1,\varepsilon}^2)\le1\), so
the rank-one curvature conclusion (3.2) survives.  Letting
\(\varepsilon\downarrow0\) repairs well-posedness but returns exactly to
the hard adaptive-survivor term in (0.5).

## 6. Canonical tests

### 6.1 Gaussian linear mode

Let \(\mu=\gamma_n\) and \(f(x)=x_1\).  Then \(\lambda=1\),
\(b_0=e_1\), and the hard feedback remains
\[
 C_t=P_{e_1^\perp}.
\]
The posterior covariance is
\[
 A_t=1\oplus(1+t)^{-1}I_{n-1},
\]
while
\[
 m_t=0,qquad V_t=1,qquad\mathcal E_t=1,
 \qquad v_T=e_1,qquad R_T=1.                            \tag{6.1}
\]
The hard ledger is exact and the survivor is the eigenfunction direction.

For the soft control, \(C\) is constant with eigenvalue \(1-\alpha\) on
\(e_1\).  Writing \(\delta=(1-\alpha)^2\), one obtains exactly
\[
 \mathbb EV_T=\frac1{1+\delta T},
 \qquad Q_T=\frac{\delta T}{1+\delta T}.                 \tag{6.2}
\]
This is the elementary Gaussian realization of the information--curvature
tradeoff.

### 6.2 A product first eigenfunction

Let \(\nu\) be a smooth isotropic one-dimensional log-concave law with an
attained normalized first eigenfunction \(\phi\), and assume
\(a=\operatorname{Cov}_\nu(\phi,Z)\ne0\).  Put
\[
 \mu=\nu\otimes\gamma_{n-1},
 \qquad f(x)=\phi(x_1).
\]
If \(\lambda_1(\nu)\le1\), this is a first eigenfunction of \(\mu\).
The hard control is again \(P_{e_1^\perp}\).  It localizes only the
Gaussian factors and leaves \(\nu\), \(f\), and their Rayleigh quotient
unchanged:
\[
 V_t=1,qquad\mathcal E_t=\lambda_1(\nu),
 \qquad R_T=\operatorname{Var}_\nu Z=1.                  \tag{6.3}
\]
This test shows exactly what the hard flow does: it removes directions
which carry conditional-mean information and preserves the entire
one-dimensional spectral problem.

### 6.3 Zero signal is substantive

Nothing in the identities proves \(b_t\ne0\).  The relation
\(\int\nabla f=\lambda b_0\) does not give a lower bound on \(|b_0|\).
Thus any proof using the literal projector must either establish a new
nonvanishing theorem for first eigenfunctions or use the relaxed/smoothed
formulations of Section 1.  This issue is independent of the adaptive
survivor estimate; resolving it does not improve (0.5).

## 7. Audited status and the remaining possible route

The following statements are dimension-free and complete, subject only to
the standard stopped/relaxed construction:

1. the fixed-function mean, variance, energy, and \(b_t\) equations;
2. pathwise preservation of \(m_t\) under \(C_tb_t=0\);
3. exact conservation of expected posterior variance and energy;
4. all-but-one curvature at time \(T\);
5. the rank-one-defect estimate with constant \(96\);
6. the lower ledger (0.6) and the energy-biased bound (0.7);
7. the soft information--curvature tradeoff (0.11)--(0.12).

These identities do not prove a dimension-free spectral gap.  They show
why the proposed shortcut cannot avoid a covariance estimate.  The precise
remaining quantity is not the ordinary posterior operator norm but the
more specialized energy-biased adaptive survivor
\[
 \mathbb E_{\mathbb Q_T}
 \operatorname{Var}_{\mu_T}\langle X,v_T\rangle.
\]
Bounding it numerically would prove KLS by (0.9), while a hypothetical
small-gap eigenfunction forces it to be \(\Omega(1/\lambda)\) by (0.7).

The only plausible function-specific escape left by this audit is to use
the dynamics of the protected signal itself: the tensor \(H_t\) and the
angular quadratic variation (2.12).  One would need a theorem saying that
creation of a long survivor along \(v_T\) forces enough rotation or enough
localized energy loss to contradict (2.8).  Such a theorem would be new;
it is not a consequence of finite-time trace curvature, log-concavity of
the endpoint, or the first-eigenfunction identity at time zero.
