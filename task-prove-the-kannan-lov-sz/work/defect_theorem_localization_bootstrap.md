# What the defect theorem gives under stochastic localization

## 0. Outcome

For \(d\ge1\), define the monotone worst-case Poincare constant

\[
 \mathcal K_d=\sup\left\{C_P(\nu):
  \begin{array}{l}
   \nu\text{ is isotropic and log-concave on an affine space}\\
   \text{of dimension at most }d
  \end{array}\right\}.                                  \tag{0.1}
\]

The weak-subspace defect theorem says that if
\(\mathbb R^n=F\oplus E\) and

\[
                         D^2V\succeq\kappa P_E,          \tag{0.2}
\]

then, writing \(\bar\nu=(P_F)_\#\nu\),

\[
 C_P(\nu)\le\kappa^{-1}
     +2\bigl(C_P(\bar\nu)+\kappa^{-1}\bigr)
 \le3\bigl(\kappa^{-1}+C_P(\bar\nu)\bigr).               \tag{0.3}
\]

This note inserts (0.3) into two localization schemes.  Because
\(\mathcal K_d\) is a Poincare constant, its square root is the KLS length
scale \(1/\psi\), up to universal constants.  The exact raw **actual
trace-based dimension recurrences** obtainable from only the standard
measure-martingale and covariance identities, when one retains the
dimension-\(m\) affine estimate for the weak marginal, are

\[
 \boxed{\mathcal K_n\le C\left({n\over m}
                         +n\mathcal K_m\right)}
 \qquad(1\le m\le n)                                    \tag{0.4}
\]

for a full covariance-whitened localization, and

\[
 \boxed{\mathcal K_n\le C\left[\left({n\over m}\right)^2
                         +n\mathcal K_m\right]}          \tag{0.5}
\]

for a spectrum-adaptive water-filling stopping rule.  All constants in
these two displays are numerical.  There is also a standard trace
Poincare estimate \(C_P(\nu)\le C\operatorname {Tr}\operatorname {Cov}\nu\).
Using it on the weak marginal gives the sharper *nonrecursive* trace-only
conclusion \(\mathcal K_n\le Cn\).  Thus (0.4)--(0.5) are the exact raw
dimension recurrences, but they are already dominated by the elementary
trace conclusion.  We keep them because they expose the term a hypothetical
covariance improvement would have to change.

Neither actual recurrence improves the known KLS-length bound
\(\sqrt{\mathcal K_n}\lesssim\sqrt{\log n}\), equivalently
\(\mathcal K_n\lesssim\log n\).  For a separate, explicitly
**counterfactual audit**, grant the currently unavailable estimate that
the weak covariance is \(O(1)\).  Only under that hypothetical assumption
does the first recurrence become

\[
                         \mathcal K_n\le C(n/m+\mathcal K_m), \tag{0.6}
\]

and inserting \(\mathcal K_m\lesssim\log m\) cannot produce
\(o(\log n)\), equivalently an unbounded-factor improvement of the
\(\sqrt{\log n}\) KLS-length bound.  The obstruction is therefore not
merely a poor constant: the curvature cost \(n/m\), and in the actual
trace-based argument the additional weak-covariance factor \(n\), remain.

The exact new estimate that would be required is identified in Section 5.
No estimate on \(\int\|A_t\|_{\mathrm{op}}dt\), uniform covariance-process
bound, or KLS-strength statement is assumed below.

## 1. Localization identities and the transfer back to the initial measure

We first work with a smooth, full-dimensional log-concave probability
\(p_0=e^{-V_0}\) having all moments.  The usual truncation and convolution
limits are uniform in every parameter used below.  More explicitly, for
an arbitrary isotropic log-concave \(\mu\), first condition on \(B_R\),
convolve intrinsically with a Gaussian of variance \(\varepsilon\), and
apply the affine map that makes the resulting measure isotropic.  Let
\(R\uparrow\infty\) and then \(\varepsilon\downarrow0\).  Log-concave
exponential tails give convergence of the barycenters and covariance
matrices; hence the isotropizing maps converge to the identity on the
affine support.  Once the estimate below is proved for the approximants,
apply it to \(C_c^\infty\) test functions and pass to the limit in
variance and Dirichlet energy.  Density by value truncation, spatial
cutoff, and mollification gives the complete locally Lipschitz Poincare
class.  Thus it suffices to prove the estimates in the smooth setting.
The localization regularization \(\delta\) below keeps every driver
bounded, and all resulting constants are independent of
\(R,\varepsilon,\delta\).

Let \(C_t\) be a bounded predictable symmetric positive semidefinite
matrix.  Eldan's measure-valued martingale is

\[
 dp_t(x)=p_t(x)\langle C_t(x-a_t),dW_t\rangle,\qquad
 a_t=\int x\,p_t(x)dx.                                  \tag{1.1}
\]

It has the explicit form

\[
 p_t(x)=Z_t^{-1}
 \exp\left(\langle c_t,x\rangle-\tfrac12\langle Q_tx,x\rangle\right)
 p_0(x),
 \qquad
 Q_t=\int_0^t C_s^2ds.                                  \tag{1.2}
\]

In particular, if \(V_t=-\log p_t\), then

\[
                         D^2V_t\succeq Q_t.              \tag{1.3}
\]

For an integrable function \(h\), put
\(m_t(h)=\int h\,dp_t\), and let \(A_t=\operatorname {Cov}(p_t)\).  Direct
integration of (1.1) gives

\[
 dm_t(h)=\operatorname {Cov}_{p_t}(h,X)^TC_t\,dW_t.      \tag{1.4}
\]

The covariance satisfies

\[
 dA_t=\mathcal T_t(C_t\,dW_t)-A_tC_t^2A_t\,dt,           \tag{1.5}
\]

where \(\mathcal T_t\) is the centered third-moment tensor.  We will never
bound its martingale part.  The facts needed below follow without such a
bound:

\[
 \mathbb E A_t=I-\mathbb E(a_ta_t^T)\preceq I,
 \qquad
 \mathbb E\operatorname {Tr}A_t\le n,                   \tag{1.6}
\]

and \(\operatorname {Tr}A_t\) is a nonnegative supermartingale.  The first
identity follows by integrating \(xx^T\) in the measure martingale.  The
second assertion follows from (1.5), whose trace drift is
\(-\operatorname {Tr}(A_tC_t^2A_t)\le0\).  Consequently, for every bounded
stopping time \(\tau\),

\[
 \mathbb P\left\{\sup_{s\le\tau}\operatorname {Tr}A_s\ge R\right\}
 \le {n\over R}.                                        \tag{1.7}
\]

### 1.1 Set-mass survival

Let \(S\) satisfy \(p_0(S)=1/2\), and set \(M_t=p_t(S)\).  If
\(b_t=\operatorname {Cov}_{p_t}(1_S,X)\), then

\[
 d[M]_t=|C_tb_t|^2dt
 \le M_t(1-M_t)\|C_tA_tC_t\|_{\mathrm{op}}dt.            \tag{1.8}
\]

Indeed, for every vector \(v\),

\[
 |\langle b_t,v\rangle|^2
 \le M_t(1-M_t)\langle A_tv,v\rangle,
\]

and (1.8) follows by taking the supremum after replacing \(v\) by
\(C_tv\).  Therefore, if

\[
 \int_0^\tau\|C_tA_tC_t\|_{\mathrm{op}}dt\le T
 \quad\text{almost surely},                             \tag{1.9}
\]

then

\[
 \mathbb P\{|M_\tau-\tfrac12|\ge\tfrac14\}
 \le {\mathbb E\,[M]_\tau\over(1/4)^2}
 \le4T.                                                  \tag{1.10}
\]

Thus \(T\) is the total set-mass survival budget.

### 1.2 Boundary and Poincare transfer

Optional sampling gives \(\mathbb E p_\tau(B)=p_0(B)\) for every Borel
set \(B\).  Fatou's lemma applied to
\([p_\tau(S_\varepsilon)-p_\tau(S)]/\varepsilon\) gives

\[
                         p_0^+(S)\ge\mathbb E p_\tau^+(S). \tag{1.11}
\]

Suppose an event \(G\) has probability at least \(q_0>0\), on \(G\)

\[
 C_P(p_\tau)\le D,
 \qquad
 M_\tau\in[1/4,3/4].                                    \tag{1.12}
\]

The Buser--Ledoux inequality for log-concave probabilities, including
lower-dimensional affine supports, gives
\(\psi_{p_\tau}\ge cD^{-1/2}\).  Hence (1.11) yields

\[
                         p_0^+(S)\ge {cq_0\over4\sqrt D}. \tag{1.13}
\]

The standard localization/concavity theorem for the isoperimetric profile
of a log-concave probability implies

\[
 \psi_\nu\ge c\inf\{\nu^+(S):\nu(S)=1/2\}.               \tag{1.14}
\]

Indeed, if
\(I_\nu(s)=\inf_{\nu(A)=s}\nu^+(A)\), the log-concave isoperimetric
profile is symmetric and concave after the standard lower-semicontinuous
regularization.  Hence \(I_\nu(s)/s\ge2I_\nu(1/2)\) for \(s\le1/2\).
The same conclusion, with a universal constant, follows directly from
one-dimensional needle localization and is stable under affine-support
approximation.  Combining (1.13), (1.14), and Cheeger's inequality gives

\[
                         C_P(p_0)\le C(q_0)D.             \tag{1.15}
\]

All later applications have a fixed numerical \(q_0\).

Finally, if \(\nu\) is log-concave of dimension \(k\), with covariance
\(B\), affine whitening and reduction to the support give

\[
                         C_P(\nu)\le\mathcal K_k\|B\|_{\mathrm{op}}. \tag{1.16}
\]

This is the only use of a worst-dimensional constant.

We will also audit against the standard trace estimate

\[
                         C_P(\nu)\le C\operatorname {Tr}B. \tag{1.17}
\]

For completeness, if \(X,Y\) are independent with law \(\nu\), then every
one-Lipschitz \(f\) satisfies
\[
 \mathbb E|f-\operatorname {med}f|
 \le\mathbb E|f-\mathbb Ef|
 \le\mathbb E|f(X)-f(Y)|
 \le\sqrt{2\operatorname {Tr}B}.
\]
E. Milman's first-moment/Poincare equivalence for log-concave
probabilities gives (1.17), with a universal constant and with the same
statement on affine supports.

## 2. Full covariance-whitened localization

Fix a small \(\delta>0\), and use

\[
                         C_t=(A_t+\delta I)^{-1/2}.       \tag{2.1}
\]

This is a bounded predictable driver and

\[
 C_tA_tC_t=A_t(A_t+\delta I)^{-1}\preceq I.             \tag{2.2}
\]

Thus a deterministic time \(T\) uses at most survival budget \(T\).
Let

\[
 Q_T=\int_0^T(A_t+\delta I)^{-1}dt.                      \tag{2.3}
\]

For a threshold \(\kappa>0\), let \(F\) be the spectral subspace of
\(Q_T\) corresponding to eigenvalues below \(\kappa\), and put
\(k=\dim F\).  On \(E=F^\perp\), (1.3) gives

\[
                         D^2V_T\succeq\kappa P_E.        \tag{2.4}
\]

### Lemma 2.1 (number of weak directions)

If \(\delta\le T/(2\kappa)\), then

\[
                         \mathbb E k\le {2\kappa n\over T}. \tag{2.5}
\]

#### Proof

Let \(v_1,\ldots,v_k\) be an orthonormal basis of the terminal subspace
\(F\).  For each \(i\),

\[
 \int_0^T\langle(A_t+\delta I)^{-1}v_i,v_i\rangle dt
 =\langle Q_Tv_i,v_i\rangle<\kappa.                     \tag{2.6}
\]

For every positive definite matrix \(B\) and unit vector \(v\),
\(\langle Bv,v\rangle\langle B^{-1}v,v\rangle\ge1\).
Cauchy--Schwarz in time and (2.6) therefore give

\[
 T^2
 <\kappa\int_0^T\bigl(\langle A_tv_i,v_i\rangle+\delta\bigr)dt. \tag{2.7}
\]

Summing over \(i\) yields

\[
 k(T^2/\kappa-\delta T)
 \le\int_0^T\operatorname {Tr}(P_FA_t)dt
 \le\int_0^T\operatorname {Tr}A_tdt.                    \tag{2.8}
\]

Take expectations, use (1.6), and use the assumed bound on \(\delta\).
\(\square\)

Let

\[
                         R_F=\|P_FA_TP_F\|_{\mathrm{op}}. \tag{2.9}
\]

The weak marginal of \(p_T\) has covariance \(P_FA_TP_F\).  Hence, on
\(\{k\le m,\ \operatorname {Tr}A_T\le R\}\), (0.3), (1.16), and
(1.17) give

\[
 C_P(p_T)\le3\left(\kappa^{-1}
       +\min\{R\mathcal K_m,CR\}\right).                 \tag{2.10}
\]

The three failure probabilities are explicitly

\[
\begin{aligned}
 \mathbb P\{|M_T-\tfrac12|\ge\tfrac14\}&\le4T,\\
 \mathbb P\{k>m\}&\le {2\kappa n\over Tm},\\
 \mathbb P\{\operatorname {Tr}A_T>R\}&\le {n\over R}.   \tag{2.11}
\end{aligned}
\]

Consequently, whenever

\[
                         4T+{2\kappa n\over Tm}+{n\over R}<1, \tag{2.12}
\]

the localization transfer proves the parameterized recurrence

\[
 \boxed{\mathcal K_n
 \le C\bigl(\kappa^{-1}+R\mathcal K_m\bigr),}            \tag{2.13}
\]

where \(C\) depends only on the numerical margin in (2.12).
Independently, the trace alternative in (2.10) gives

\[
                         \mathcal K_n\le C(\kappa^{-1}+R). \tag{2.13a}
\]

For a concrete choice, take

\[
 T={1\over64},\qquad
 \kappa={Tm\over16n},\qquad
 R=8n,\qquad
 0<\delta\le {T\over2\kappa}.                            \tag{2.14}
\]

The total failure probability in (2.11) is at most
\(1/16+1/8+1/8=5/16\).  Equations (2.13)--(2.14) give

\[
                         \mathcal K_n
 \le C\left({n\over m}+n\mathcal K_m\right),             \tag{2.15}
\]

which is (0.4).  Every parameter is fixed before the smoothing limit, and
all estimates are uniform as \(\delta\downarrow0\).

With the same parameters, (2.13a) gives
\[
                         \mathcal K_n\le C(n/m+n).
\]
Its best trace-only consequence is the familiar \(O(n)\) bound.  Hence
(2.15) is not advertised as a numerically sharper bound; it is the
dimension-recursive form whose fixed-point behavior is under audit.

### 2.1 Why the weak-covariance factor cannot be dropped from these inputs

The only available estimate for \(R_F\) above is
\(\mathbb E\operatorname {Tr}A_T\le n\).  It forces a constant-probability
threshold \(R\asymp n\).  This loss is not an artifact of replacing an
operator norm by a trace after taking expectation.  The matrix-valued
random variable

\[
                         A=n\,vv^T,                      \tag{2.16}
\]

where \(v\) is uniform on \(S^{n-1}\), satisfies
\(\mathbb EA=I\), while \(\|A\|_{\mathrm{op}}=n\) almost surely.  If the
weak subspace is \(\operatorname {span}(v)\), it also has dimension one.
Thus the moment identity (1.6), even combined with a very small weak
dimension, logically permits \(R_F=n\).  Formula (2.16) is an audit of the
information used here, not a claim that this abstract matrix law is a
localization trajectory.  Ruling it out for actual log-concave
localization is precisely an additional covariance-process theorem.

For comparison, on the isotropic Gaussian the obstruction is absent.
With \(\delta=0\), Gaussian posterior calculus gives
\[
 A_t=(I+Q_t)^{-1},\qquad
 \dot Q_t=A_t^{-1}=I+Q_t,
\]
and therefore \(Q_t=(e^t-1)I\), \(A_t=e^{-t}I\).  Every direction acquires
constant curvature in constant survival time.  Thus the losses in
(2.5) and (2.11) measure possible non-Gaussian covariance fluctuations;
they are not caused by a normalization error in the driver.

## 3. Fixed-point audit of the first recurrence

### 3.1 Counterfactual bounded-weak-covariance audit

There are two independent losses in (2.15).

First, even if one grants the unproved ideal estimate \(R_F\le C\) with
constant probability, (2.5) forces

\[
                         \kappa\lesssim {m\over n}        \tag{3.1}
\]

when one wants \(k\le m\) with constant probability.  The defect theorem
then costs \(\kappa^{-1}\gtrsim n/m\), giving (0.6).

Insert the known input
\(\mathcal K_j\le C_0\log(e+j)\).  For every \(2\le m\le n\),

\[
 {n\over m}+\log(e+m)
 \ge c\log(e+n).                                        \tag{3.2}
\]

Indeed, if \(m\le n/\log(e+n)\), the first term has the required
size.  Otherwise
\(\log(e+m)\ge\log n-\log\log(e+n)\), which is at least a fixed
multiple of \(\log(e+n)\) outside a finite range.  Therefore the map

\[
 B_n\longmapsto\inf_{m\le n}\{n/m+B_m\}                 \tag{3.3}
\]

does not turn the \(\log n\) Poincare input into \(o(\log n)\).  After
taking square roots, it cannot improve the \(\sqrt{\log n}\) KLS-length
bound by an unbounded factor, so it cannot reach milestone M1.

### 3.2 Return to the actual trace-based recurrence

The actual covariance estimate gives \(R\asymp n\), so (2.15) is
much weaker still.  Since \(\mathcal K_m\ge1\) (test a linear function in
an isotropic measure), optimizing (2.15) yields no better than \(O(n)\).

Choosing \(m=\rho n\) and iterating does not create a fixed point.
Even under the ideal covariance bound, one only gets

\[
                         \mathcal K_n\le C(1+\mathcal K_{\rho n}). \tag{3.4}
\]

The coefficient is not a contraction, and the iteration has
\(\Theta(\log n)\) levels.  Rescaling the driver to spend only
\(1/\log n\) survival budget per level rescales \(Q_t\) by the same factor,
reintroducing the reciprocal curvature loss.  Thus a multilevel
restatement of (2.13) does not remove the dimension dependence.

For comparison, the isotropic normalized driver

\[
                         C_t={I\over\sqrt{\|A_t\|_{\mathrm{op}}}} \tag{3.5}
\]

also has \(\|C_tA_tC_t\|_{\mathrm{op}}=1\), but its accumulated tilt is
\[
 Q_T=s_TI,\qquad
 s_T=\int_0^T{dt\over\|A_t\|_{\mathrm{op}}}.             \tag{3.6}
\]

Cauchy--Schwarz gives

\[
 s_T\ge {T^2\over\int_0^T\|A_t\|_{\mathrm{op}}dt}.       \tag{3.7}
\]

Using only \(\mathbb E\operatorname {Tr}A_t\le n\) gives, with constant
probability, \(s_T\gtrsim T/n\), again only an \(O(n)\) Poincare bound.
A constant lower bound on \(s_T\) is exactly a dimension-free bound on
\(\int_0^T\|A_t\|_{\mathrm{op}}dt\), the prohibited M3/KLS-strength
estimate.  The covariance-whitened driver in Section 2 avoids assuming
that estimate, but pays for the surviving high-covariance directions in
(2.10).

## 4. A spectrum-adaptive water-filling stop

We next change the driver, rather than merely changing the threshold in
Section 2.  The process observes only directions whose accumulated
curvature has not yet reached \(\kappa\).

Let

\[
 P_t=1_{[0,\kappa)}(Q_t),\qquad r_t=\operatorname {rank}P_t,\qquad
 B_t=P_tA_tP_t\big|_{\operatorname {Ran}P_t}.            \tag{4.1}
\]

With a regularization \(\delta>0\), use

\[
 C_t=P_t(B_t+\delta I_{\operatorname {Ran}P_t})^{-1/2}P_t. \tag{4.2}
\]

At spectral crossing times one uses a smooth cutoff and then takes a
monotone limit.  Before the trace and time stops below, (4.2) is bounded;
the standard well-posedness theorem for predictable bounded localization
drivers applies.  Moreover,

\[
 C_tA_tC_t=B_t(B_t+\delta I)^{-1}\preceq P_t,            \tag{4.3}
\]

so the set-mass budget is again at most elapsed time.

Stop when \(r_t\le m\), or at time \(T\), whichever comes first.  Also
consider the good trace event

\[
                         H_R=\{\sup_{t\le T}\operatorname {Tr}A_t<R\}. \tag{4.4}
\]

By (1.7), \(\mathbb P(H_R^c)\le n/R\).

Define the remaining curvature deficit

\[
                         \Phi_t=\operatorname {Tr}(\kappa I-Q_t)_+. \tag{4.5}
\]

For almost every time before the rank stop, spectral calculus gives

\[
\begin{aligned}
 -{d\Phi_t\over dt}
 &=\operatorname {Tr}\bigl(P_tC_t^2\bigr)
   =\operatorname {Tr}(B_t+\delta I)^{-1}\\
 &\ge {r_t^2\over\operatorname {Tr}B_t+\delta r_t}.      \tag{4.6}
\end{aligned}
\]

The last line is the arithmetic-harmonic mean inequality.  On \(H_R\), if
the rank has not reached \(m\), then

\[
                         -\Phi_t'\ge {m^2\over R+\delta n}. \tag{4.7}
\]

Since \(\Phi_0=n\kappa\) and \(\Phi_t\ge0\), the rank stop must occur
before \(T\) on \(H_R\) whenever

\[
                         n\kappa<{Tm^2\over R+\delta n}. \tag{4.8}
\]

At that stop, the subspace \(F=\operatorname {Ran}P_\tau\) has dimension
at most \(m\), \(Q_\tau\succeq\kappa P_{F^\perp}\), and on \(H_R\)

\[
                         \|P_FA_\tau P_F\|_{\mathrm{op}}<R. \tag{4.9}
\]

Take, for example,

\[
 T={1\over64},\qquad R=16n,\qquad
 0<\delta\le {R\over n},\qquad
 \kappa={Tm^2\over64n^2}.                               \tag{4.10}
\]

Then (4.8) holds strictly.  The trace failure probability is at most
\(1/16\), and the set-mass survival failure probability is at most
\(4T=1/16\).  Thus the posterior estimate (0.3), boundary transfer, and
(4.9) give

\[
 \mathcal K_n
 \le C\left(\kappa^{-1}+R\mathcal K_m\right)
 \le C\left[\left({n\over m}\right)^2+n\mathcal K_m\right], \tag{4.11}
\]

which is (0.5).

Using (1.17) instead gives the sharper trace-only but nonrecursive bound
\[
                         \mathcal K_n
 \le C\left[\left({n\over m}\right)^2+n\right].          \tag{4.11a}
\]

This rule genuinely adapts to the spectrum: completed directions are no
longer driven, and the stopping time is defined by the number of
unfinished eigenvalues.  Nevertheless it is worse than (2.15).  The
deficit calculation exposes why: under only a trace ceiling, filling
\(n-m\) curvature levels forces the cost
\(\kappa^{-1}\gtrsim(n/m)^2\).  The same weak-covariance threshold
\(R\asymp n\) remains.

There is also an instantaneous optimality check.  Any symmetric driver
satisfying the unit survival constraint
\[
                         C A C\preceq I                  \tag{4.12}
\]
obeys
\[
                         C^2\preceq A^{-1}.              \tag{4.13}
\]
Indeed, \(CAC\preceq I\) implies
\(\|CA^{1/2}\|_{\mathrm{op}}\le1\), hence
\(A^{1/2}C^2A^{1/2}\preceq I\).  Full whitening attains equality in
(4.13).  A different spectrum allocation can reduce covariance
fluctuations by observing fewer directions, but it cannot accumulate more
quadratic curvature at the current covariance without spending more
set-mass budget.  Estimate (4.11) quantifies this tradeoff for the
water-filling allocation.

## 5. Exact missing estimates

The defect theorem converts a localization posterior into a useful
global estimate once the following three quantities are controlled:

\[
 \kappa,\qquad
 k=\dim 1_{[0,\kappa)}(Q_\tau),\qquad
 R_F=\|P_FA_\tau P_F\|_{\mathrm{op}}.                    \tag{5.1}
\]

The standard identities prove only

\[
 \mathbb Ek\lesssim{\kappa n\over T},
 \qquad
 \mathbb P\{R_F>R\}\le {n\over R}.                       \tag{5.2}
\]

To improve the \(\sqrt{\log n}\) KLS-length bound by an unbounded factor
through (2.13), one would need, at constant survival budget, a statement
such as

\[
\begin{split}
 &\mathbb P\{\lambda_{m+1}(Q_T)\ge\kappa_0,\ 
       \|P_FA_TP_F\|_{\mathrm{op}}\le R_0\}\ge c_0,\\
 &\kappa_0,R_0,c_0^{-1}=O(1),
 \qquad \log m=o(\log n).                               \tag{5.3}
\end{split}
\]

Here \(\lambda_1(Q_T)\le\cdots\le\lambda_n(Q_T)\), so the first event in
(5.3) says that at most \(m\) curvature eigenvalues are below
\(\kappa_0\).

Indeed, (5.3) and the known bound in dimension \(m\) would give
\(\mathcal K_n\lesssim1+\log m=o(\log n)\), and hence
\(\sqrt{\mathcal K_n}=o(\sqrt{\log n})\).  For a full dimension-free
conclusion one needs
\(m=O(1)\), or a genuinely contractive recurrence replacing the additive
defect estimate.

Neither clause of (5.3) follows from (1.6):

* the first requires constant curvature in all but \(n^{o(1)}\)
  directions, whereas (2.5) gives only a constant-fraction reduction at
  constant \(\kappa\);
* the second is a covariance-aligned tail bound, and the random-needle
  matrix audit (2.16) shows that averaged isotropy does not imply it.

For the isotropic driver (3.5), the analogous missing statement is exactly
\[
                         \mathbb E\int_0^T
                         \|A_t\|_{\mathrm{op}}dt\le C,   \tag{5.4}
\]
which is the prohibited M3 estimate.  For the whitened and water-filling
drivers, (5.3) is the corresponding spectrum-resolved form of the same
obstruction.  Proving it would be a materially new covariance-process
theorem, not a consequence of the defect lemma or of standard martingale
identities.

Thus the weak-subspace theorem is useful and exact once a localization
leaves only controlled weak directions, but the two standard
spectrum-selection mechanisms above do not themselves create such a
posterior strongly enough to yield M1, let alone a dimension-free fixed
point.
