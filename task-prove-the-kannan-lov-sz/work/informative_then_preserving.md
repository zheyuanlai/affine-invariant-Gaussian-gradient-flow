# Informative-then-preserving localization: exact control and the switch obstruction

## 1. Outcome

Let

\[
 g_t=\mu_t(S),\qquad
 v_t=\operatorname{Cov}_{\mu_t}(\mathbf1_S,X),\qquad
 e_t=\frac{v_t}{|v_t|}.
\]

Consider the proposed two-phase control:

\[
 \text{Phase A: } C_t=e_te_t^T,
 \qquad
 \text{Phase B: } C_t=P_t:=I-e_te_t^T.
 \tag{1.1}
\]

The following facts are proved below.

1. Phase A is a genuinely scalar localization.  The exact equations are driven by one Brownian motion, and the Hilbert-valued natural-parameter martingale has quadratic variation exactly \(t\), independent of dimension.

2. While \(\|c_t\|\), \(\operatorname{tr}Q_t\), and \(g_t\) remain in fixed numerical ranges, the entire posterior covariance is bounded by a numerical multiple of \(I\).  This follows from a two-dimensional log-concave mgf estimate, not from KLS.

3. If \(|v_0|\ge\eta>0\), then for time \(c\eta^2\) there is universal positive probability that \(e_t\) stays aligned with \(e_0\), the set mass stays central, and
   \[
   e_0^TQ_Ae_0\ge c\eta^2.
   \]
   This is a directional Rayleigh exposure, not a Loewner rank-one lower bound.  It is the strongest seed furnished by the scalar martingale control and degenerates quadratically as \(|v_0|\downarrow0\).

4. After the switch,
   \[
   Q_A+Q_B
   =Q_A+T_BI-\int_{\mathrm B}e_te_t^T\,dt.
   \]
   If this matrix has a small eigenvalue, then Phase B must lock for almost all its duration onto a line which received almost no Phase-A exposure.  This is an exact deterministic alternative.

5. There is no dimension-free lower bound on the signal needed by the Phase-A lemma.  For Gaussian sign parity in dimension \(n\ge3\),
   \[
   v_0=0,\qquad D_0:=\mathbb E[(\mathbf1_S-\tfrac12)XX^T]=0,
   \]
   and in fact every label-correlated polynomial moment of degree \(<n\) vanishes.  Thus no dimension-independent bounded-degree SDE or small-parameter expansion can create a universal informative direction at \(v=0\).

The hybrid therefore narrows the missing theorem but does not prove a universal full-rank seed.  The unresolved event is sharply specified: Phase B would have to select, in a short initial layer, a high-variance line which was missed by Phase A and then protect it coherently.  Radial Gaussian laws cannot realize the high-variance part; the product-exponential sum does not realize the missed-line part; and the product-exponential maximum remains a sharp winner-selection stress test rather than a proved counterexample.

No KLS estimate, covariance-operator survival theorem, or second-order Poincare inequality is used.

## 2. Stopped localization and exact SDEs

Write

\[
 p_t(x)=Z_t^{-1}
 \exp\left(c_t\cdot x-\frac12x^TQ_tx\right)p_0(x),
\]

with

\[
 dc_t=C_t\,dW_t+C_t^2a_t\,dt,\qquad
 dQ_t=C_t^2dt,\qquad a_t=\mathbb E_tX.
 \tag{2.1}
\]

Stop before \(g_t\notin(\delta,1-\delta)\) or the displayed natural parameters leave fixed compact ranges.  Put

\[
 Y=X-a_t,\qquad h=\mathbf1_S-g_t,
\]

\[
 A_t=\mathbb E_tYY^T,\qquad
 D_t=\mathbb E_t[hYY^T],\qquad
 \mathcal T_t(z)=\mathbb E_t[YY^T\langle Y,z\rangle].
 \tag{2.2}
\]

The density martingale equation is

\[
 dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle.
 \tag{2.3}
\]

For any symmetric predictable \(C_t\),

\[
 \boxed{
 \begin{aligned}
 da_t&=A_tC_t\,dW_t,\\
 dg_t&=v_t^TC_t\,dW_t,\\
 dv_t&=D_tC_t\,dW_t-A_tC_t^2v_t\,dt,\\
 dA_t&=\mathcal T_t(C_t\,dW_t)-A_tC_t^2A_t\,dt.
 \end{aligned}}
 \tag{2.4}
\]

The drift in \(dv_t\) is the quadratic covariation between \(g_t\) and \(a_t\); omitting it gives an incorrect Phase-A equation.

### 2.1 Phase A

Let \(\ell_t=|v_t|>0\), \(e_t=v_t/\ell_t\), \(P_t=I-e_te_t^T\), and

\[
 b_t=D_te_t.
\]

Since \(C_t=e_te_t^T\), all noise is the scalar

\[
 d\beta_t=e_t^TdW_t.
\]

Equation (2.4) becomes

\[
 \boxed{
 \begin{aligned}
 da_t&=A_te_t\,d\beta_t,\\
 dg_t&=\ell_t\,d\beta_t,\\
 dv_t&=b_t\,d\beta_t-\ell_tA_te_t\,dt,\\
 dA_t&=\mathcal T_t(e_t)\,d\beta_t
       -A_te_te_t^TA_t\,dt,\\
 dQ_t&=e_te_t^Tdt.
 \end{aligned}}
 \tag{2.5}
\]

Ito normalization of \(v_t\) gives

\[
 \boxed{
 \begin{aligned}
 d\ell_t
 &=e_t^Tb_t\,d\beta_t
   -\ell_te_t^TA_te_t\,dt
   +\frac{|P_tb_t|^2}{2\ell_t}\,dt,\\
 de_t
 &=\frac{P_tb_t}{\ell_t}\,d\beta_t
   -P_tA_te_t\,dt
   -\frac{P_tb_t\,b_t^Te_t}{\ell_t^2}\,dt
   -\frac{|P_tb_t|^2}{2\ell_t^2}e_t\,dt,\\
 d[e]_t
 &=\frac{P_tb_tb_t^TP_t}{\ell_t^2}\,dt.
 \end{aligned}}
 \tag{2.6}
\]

The singular factor \(|P_tb_t|/\ell_t\) is the exact small-\(v\) obstruction.

The set-mass information expenditure is

\[
 d\langle g\rangle_t=\ell_t^2dt.
 \tag{2.7}
\]

The Phase-A precision is

\[
 Q_A(t)=\int_0^t e_se_s^T\,ds,\qquad
 \operatorname{tr}Q_A(t)=t.
 \tag{2.8}
\]

### 2.2 Phase B

With \(C_t=P_t\), equation (2.4) gives

\[
 \boxed{
 \begin{aligned}
 dg_t&=0,\\
 da_t&=A_tP_t\,dW_t,\\
 dv_t&=D_tP_t\,dW_t,\\
 dA_t&=\mathcal T_t(P_t\,dW_t)-A_tP_tA_t\,dt,\\
 dQ_t&=P_t\,dt.
 \end{aligned}}
 \tag{2.9}
\]

Thus the set mass is preserved exactly after the switch.  If Phase B has duration \(T_B\),

\[
 Q_B=T_BI-\int_{\mathrm B}e_te_t^T\,dt.
 \tag{2.10}
\]

At a zero of \(v_t\), neither (2.6) nor the literal controls in (1.1) are defined.  A rigorous formulation must specify a convention or use a relaxed rank-one control.  In Phase A this means a predictable \(0\preceq R_t\preceq I\) with \(\operatorname{tr}R_t=1\) and \(R_tv_t=v_t\) when \(v_t\ne0\); in Phase B it means \(\operatorname{tr}R_t=n-1\) and \(R_tv_t=0\).  Mesh approximation gives stopped weak limits, but it does not manufacture an informative direction at \(v=0\).

## 3. Dimension-free scalar estimates

Two elementary inequalities are central.

### Lemma 3.1 (dimension-free bound on the Phase-A vector noise)

For \(b=D e\) and \(s=e^TAe\),

\[
 \boxed{b^TA^{-1}b\le
 \mathbb E[h^2(e^TY)^2]\le s.}
 \tag{3.1}
\]

Consequently, if \(A\preceq KI\), then

\[
 |b|\le K.
 \tag{3.2}
\]

#### Proof

For every \(z\),

\[
 z^Tb=\mathbb E[h(z^TY)(e^TY)].
\]

Cauchy--Schwarz gives

\[
 |z^Tb|^2
 \le (z^TAz)\,\mathbb E[h^2(e^TY)^2].
\]

Optimize over \(z\), regularizing \(A\) if necessary.  Since \(|h|\le1\), the last expectation is at most \(s\).  If \(A\preceq KI\), then

\[
 |b|^2\le K\,b^TA^{-1}b\le K s\le K^2.
\]

\(\square\)

### Lemma 3.2 (small natural parameters bound every direction)

There are numerical \(c_0,q_0,C<\infty\) with the following property.  If \(\mu\) is centered isotropic log-concave and

\[
 d\mu_{c,Q}(x)
 =Z^{-1}e^{c\cdot x-x^TQx/2}\,d\mu(x),
\qquad Q\succeq0,
\]

with

\[
 |c|\le c_0,\qquad \operatorname{tr}Q\le q_0,
\]

then

\[
 \boxed{\operatorname{Cov}(\mu_{c,Q})\preceq CI.}
 \tag{3.3}
\]

#### Proof

Jensen gives

\[
 Z\ge
 \exp\left(\mathbb E[c\cdot X-\tfrac12X^TQX]\right)
 =e^{-\operatorname{tr}Q/2}.
 \tag{3.4}
\]

For a unit \(u\),

\[
 \mathbb E_\mu[(u\cdot X)^2
 e^{c\cdot X-X^TQX/2}]
 \le
 \mathbb E_\mu[(u\cdot X)^2e^{c\cdot X}].
 \tag{3.5}
\]

The right side depends only on the marginal on the at-most-two-dimensional span of \(u\) and \(c\).  Every isotropic log-concave law in fixed dimension has a numerical exponential tail.  Hence, for sufficiently small numerical \(c_0\),

\[
 \sup_{|c|\le c_0,\ |u|=1}
 \mathbb E[(u\cdot X)^2e^{c\cdot X}]
 \le C_0.
 \tag{3.6}
\]

Combine (3.4)--(3.6) and subtract the posterior mean.  \(\square\)

This is a fixed-dimensional marginal argument, not a covariance-survival theorem.

### Lemma 3.3 (Hilbert martingale and Phase-A survival)

Let

\[
 M_t=\int_0^t e_s\,d\beta_s.
\]

For every \(p\ge2\),

\[
 \left(\mathbb E\sup_{s\le T}|M_s|^p\right)^{1/p}
 \le C\sqrt{pT}.
 \tag{3.7}
\]

In particular, for numerical \(\theta>0\),

\[
 \mathbb E\exp\left(
 \theta\frac{\sup_{s\le T}|M_s|^2}{T}
 \right)\le C.
 \tag{3.8}
\]

These are the Hilbert-space Burkholder--Davis--Gundy inequality and its moment-series consequence; the constants do not depend on \(n\).

Under Phase A,

\[
 c_t=M_t+\int_0^t e_se_s^Ta_s\,ds.
 \tag{3.9}
\]

Lemmas 3.2--3.3 give a bootstrap: for a sufficiently small numerical \(T_0\), with probability at least a numerical \(p_0>0\),

\[
 \sup_{t\le T_0}|c_t|\le c_0,\qquad
 \operatorname{tr}Q_t=t\le q_0,\qquad
 A_t\preceq CI.
 \tag{3.10}
\]

Indeed, on the region in Lemma 3.2, every one-dimensional posterior second moment is bounded, so

\[
 |e_t^Ta_t|\le C.
\]

Thus the drift in (3.9) is \(O(T_0)\), while (3.8) controls its martingale part.  Also

\[
 \ell_t^2\le g_t(1-g_t)e_t^TA_te_t\le C/4.
 \tag{3.11}
\]

Consequently, if \(g_0=1/2\) and \(a_\delta=1/2-\delta\), Doob's inequality gives

\[
 \mathbb P\left\{
 \sup_{t\le T_0}|g_t-\tfrac12|\ge a_\delta
 \right\}
 \le \frac{CT_0}{4a_\delta^2}.
 \tag{3.12}
\]

Taking \(T_0\) smaller if needed proves simultaneous natural-parameter and mass survival with universal positive probability.
In fact, the same estimates make this probability arbitrarily close to one as \(T_0\downarrow0\).

This part of the hybrid works even when \(v_0\) is small.  What it does not control is the direction of the rank-one trace in (2.8).

## 4. What can be seeded when the initial signal is nonzero

### Proposition 4.1 (quantitative Phase-A alignment)

Suppose \(A_t\preceq KI\) on a stopped Phase-A interval, with \(K\ge1\), and

\[
 |v_0|\ge\eta,\qquad 0<\eta\le1.
\]

There is a numerical \(c>0\) such that, for

\[
 T_A=c\frac{\eta^2}{K^2},
 \tag{4.1}
\]

with probability at least \(3/4\),

\[
 \sup_{t\le T_A}|v_t-v_0|\le\frac{\eta}{3}.
 \tag{4.2}
\]

On this event,

\[
 |\langle e_t,e_0\rangle|^2\ge\frac12
 \quad(0\le t\le T_A)
 \tag{4.3}
\]

and hence

\[
 \boxed{e_0^TQ_A(T_A)e_0\ge
 \frac{T_A}{2}.}
 \tag{4.4}
\]

Taking the constants so that the survival event in Lemma 3.3 also has probability at least \(3/4\), their intersection has probability at least \(1/2\).

#### Proof

From (2.5),

\[
 v_t-v_0=\int_0^t b_s\,d\beta_s
 -\int_0^tA_sv_s\,ds.
 \]

Lemma 3.1 gives \(|b_s|\le K\).  The Hilbert-valued martingale maximal inequality gives

\[
 \mathbb P\left\{
 \sup_{t\le T_A}
 \left|\int_0^tb_s\,d\beta_s\right|>\frac{\eta}{4}
 \right\}
 \le \frac{CK^2T_A}{\eta^2}.
 \tag{4.5}
\]

Moreover, the covariance inequality

\[
 |v_t|^2\le g_t(1-g_t)e_t^TA_te_t\le K/4
\]

shows that the drift has norm at most \(K^{3/2}/2\).  Its contribution by time \(T_A\) is therefore at most \(c\eta^2/(2\sqrt K)\le c\eta/2\).  Choosing \(c\) small proves (4.2).  Normalization gives

\[
 |e_t-e_0|\le
 \frac{2|v_t-v_0|}{|v_0|}\le\frac23,
\]

which implies (4.3).  Integrating the scalar identity

\[
 e_0^TQ_A(T_A)e_0
 =\int_0^{T_A}|\langle e_t,e_0\rangle|^2dt
\]

proves (4.4).  Notice that (4.3) alone does not imply the Loewner inequality \(Q_A\succeq (T_A/2)e_0e_0^T\).  \(\square\)

The scale \(\eta^2\) cannot be removed by this argument: the martingale perturbation of \(v\) is \(O(\sqrt t)\).  The examples below show that isotropicity gives no positive universal \(\eta\).

## 5. Exact algebra at the switch

Let \(\tau\) be the Phase-A stopping time, and run Phase B for deterministic active duration \(T_B\).  Put

\[
 Q_A=\int_0^\tau e_te_t^Tdt,\qquad
 M_B=\int_\tau^{\tau+T_B}e_te_t^Tdt.
\]

Then

\[
 Q_{\mathrm{tot}}=Q_A+T_BI-M_B.
 \tag{5.1}
\]

### Lemma 5.1 (missed-line/locked-line alternative)

If

\[
 \lambda_{\min}(Q_{\mathrm{tot}})<\varepsilon,
\]

then there is a unit vector \(u\) such that

\[
 \boxed{
 u^TQ_Au<\varepsilon,
 \qquad
 \int_\tau^{\tau+T_B}
 \left(1-|\langle u,e_t\rangle|^2\right)dt<\varepsilon.
 }
 \tag{5.2}
\]

Equivalently,

\[
 \int_\tau^{\tau+T_B}
 \|e_te_t^T-uu^T\|_F^2dt<2\varepsilon.
 \tag{5.3}
\]

#### Proof

Choose a minimizing unit vector \(u\).  Both terms in

\[
 u^TQ_{\mathrm{tot}}u
 =u^TQ_Au+
 \int_\tau^{\tau+T_B}
 \left(1-|\langle u,e_t\rangle|^2\right)dt
\]

are nonnegative, proving (5.2).  Formula (5.3) uses

\[
 \|ee^T-uu^T\|_F^2=2(1-|\langle e,u\rangle|^2).
\]

\(\square\)

Thus failure of a full-rank seed is much more specific than failure of Phase-A alignment: Phase B must persist near a line almost unexposed in Phase A.

As in the pure mass-preserving construction, \(Q_B\) alone has all but at most one eigenvalue at least \(T_B/2\).  On a low-curvature path, the endpoint is strongly log-concave transverse to the \(u\) in (5.2), and the anisotropic one-flat-direction Cheeger lemma reduces the endpoint estimate to

\[
 \operatorname{Var}_{\mu_{\tau+T_B}}\langle u,X\rangle.
 \tag{5.4}
\]

Lemma 3.2 bounds this variance at the switch on the good Phase-A event.  What is not controlled is its adaptive creation during the initial layer of Phase B.  Controlling (5.4) for the path-selected, Phase-A-missed line \(u\) is precisely the remaining theorem; ordinary Loewner survival for deterministic directions does not apply.

## 6. The sharp \(v_0=0\) obstruction

Let \(X\sim N(0,I_n)\), \(n\ge3\), and define

\[
 S=\left\{\prod_{i=1}^n\operatorname{sgn}(X_i)=1\right\}.
 \tag{6.1}
\]

Then \(\gamma_n(S)=1/2\) and

\[
 \mathbf1_S-\frac12
 =\frac12\prod_{i=1}^n\operatorname{sgn}(X_i)
 \quad\text{a.s.}
 \tag{6.2}
\]

Independence and symmetry give

\[
 \boxed{v_0=0,\qquad D_0=0.}
 \tag{6.3}
\]

Indeed, every entry of \(v_0\) or \(D_0\) leaves at least one independent sign with mean zero.  More generally,

\[
 \mathbb E\left[
 \left(\mathbf1_S-\frac12\right)P(X)
 \right]=0
 \tag{6.4}
\]

for every polynomial \(P\) of total degree \(<n\): each monomial misses at least one coordinate, whose sign remains unpaired.

Therefore every label-sensitive Taylor coefficient whose associated polynomial has total degree \(<n\) vanishes at the isotropic starting point.  Equivalently, assign weight one to a derivative in \(c\) and weight two to a derivative in \(Q\): every mixed derivative of weighted order \(<n\) vanishes.  In particular, (2.5)--(2.6) provide no direction at time zero and no dimension-free bounded-order coefficient from which to start Proposition 4.1.

The obstruction is not an approximation artifact.  If the convention at \(v=0\) chooses a coordinate direction, rank-one localization in that coordinate preserves product structure and leaves at least one other fair sign; then

\[
 g_t=\frac12,\qquad v_t=0
\]

throughout that coordinate Phase A.  A rotated parity construction gives the same conclusion for any predetermined initial line.

This example does **not** have a large endpoint variance: for a Gaussian base,

\[
 A_t=(I+Q_t)^{-1}\preceq I
 \tag{6.5}
\]

pathwise.  It therefore defeats a universal *informative-direction or full-rank-seed deduction* from scalar data, but not the alternative endpoint exceptional-variance bound.  Any complete hybrid theorem must explicitly use such an alternative at \(v=0\); it cannot claim that isotropicity creates a positive informative seed.

The same construction with independent median bits of product log-concave coordinates shows that the vanishing of low-order label moments is not specifically Gaussian.  In that case, however, controlling the eventual adaptively selected tail variance is again the missing Phase-B problem.

## 7. Model tests

### 7.1 Radial Gaussian

For every set \(S\), every posterior generated by either phase from \(N(0,I)\) is Gaussian with

\[
 A_t=(I+Q_t)^{-1}\preceq I.
 \tag{7.1}
\]

Thus radial half-sets with \(v_0=0\), Gaussian fans, and sign parity all satisfy the endpoint exceptional-variance alternative with constant one.  They show that full-rank curvature is not necessary and that \(v_0=0\) must be handled through variance, not through a fictitious direction.

### 7.2 Product exponential sum

Let \(X_i=Z_i-1\) with iid rate-one exponentials, and take

\[
 S=\left\{\sum_iZ_i\ge m_n\right\},
\]

where \(m_n\) is a median of \(\operatorname{Gamma}(n,1)\).  Exchangeability gives

\[
 e_0=\frac{\mathbf1}{\sqrt n},
\qquad
 |v_0|\longrightarrow\frac1{\sqrt{2\pi}}.
 \tag{7.2}
\]

Phase A along this line preserves permutation symmetry, so \(e_t=e_0\) and

\[
 Q_A=T_Ae_0e_0^T.
\]

This model satisfies the nondegenerate hypothesis of Proposition 4.1 and does not defeat the hybrid.  Under Phase B its initial transverse angular quadratic variation is only order one, as follows from the exact Dirichlet conditional-moment calculation for the gamma sum.

### 7.3 Product exponential maximum

For

\[
 S=\{\max_iZ_i\ge L\},
\qquad (1-e^{-L})^n=\frac12,
\]

put \(q=e^{-L}=1-2^{-1/n}\).  The exact initial values are

\[
 |v_0|=\frac{\sqrt n\,Lq}{2(1-q)}
 \asymp\frac{\log n}{\sqrt n},
\qquad
 e_0=\frac{\mathbf1}{\sqrt n}.
 \tag{7.3}
\]

Again Phase A preserves permutation symmetry, so \(e_t=e_0\) despite the small signal and \(Q_A=T_Ae_0e_0^T\).  Thus the \(\eta^2\) loss in Proposition 4.1 is not always sharp; symmetry can stabilize a small informative direction.

At the initial isotropic state, the exact transverse eigenvalue of \(D_0\) gives

\[
 \left.\frac d{dt}\operatorname{tr}[e]_t\right|_{0}
 \asymp(\log n)^2.
 \tag{7.4}
\]

This is the sharp switch stress test: the protected direction can leave the Phase-A-exposed diagonal line on a \(1/\log^2n\) scale.  To defeat the hybrid it must then lock onto a transverse line of large posterior variance before that line receives order-one Phase-B curvature.

The one-coordinate large-deviation calculation for full exposure does not produce such locking: on the \(1/\log n\) scale, the largest tilt still gives exponentially small posterior tail mass in every coordinate.  Hence the maximum set is not an established counterexample.  It identifies exactly the winner-selection estimate required beyond Lemma 5.1.

### 7.4 Simplex

For a regular isotropic simplex and a half-set defined by thresholding the vertex-direction coordinate at its median, the stabilizer of that vertex fixes a unique line.  Phase A preserves this symmetry and accumulates precision on that line.  Thus the basic simplex halfspace behaves like the product-sum model, not like a missed-line counterexample.

More complicated simplex cuts can have \(v_0=0\) by symmetry, but then the same dichotomy as in Section 6 applies: no informative seed follows, while a counterexample must still produce a large endpoint variance along the Phase-B locked line.  Compactness alone gives only dimension-dependent diameter control and does not settle that adaptive variance.

## 8. What the hybrid has and has not achieved

The hybrid supplies three rigorous dimension-free ingredients:

\[
 \begin{array}{c}
 \text{scalar Phase-A natural-parameter control},\\[2mm]
 \text{a }c|v_0|^2\text{ aligned seed when }v_0\ne0,\\[2mm]
 \text{the exact missed-in-A/locked-in-B alternative (5.2).}
 \end{array}
\]

It does not supply a universal seed when \(v_0\approx0\).  Gaussian parity proves that informative moments of arbitrarily high prescribed bounded degree may vanish (by taking \(n\) larger), so this is a real structural obstruction rather than a loose estimate.  Nor does the switch algebra control the exceptional variance of a line selected from the same Phase-B path.

A complete theorem would need the following new, sharply restricted statement:

> Starting from a covariance-bounded posterior, a mass-preserving Phase B cannot, with large probability, create a high-variance line \(u\) in a short initial layer, keep \(e_t\) locked near \(u\), and simultaneously have \(u^TQ_Au\ll1\).

This is materially narrower than ordinary covariance-operator survival, but it is not proved by scalar martingale bounds, small-parameter mgfs, or isotropic initial data.  The product maximum is the relevant stress test; radial Gaussian, product sum, and simplex halfspaces do not violate it.  No explicit isotropic log-concave half-set realizing all three bad properties is currently obtained.
