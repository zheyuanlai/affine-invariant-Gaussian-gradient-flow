# Covariance-normalized mass-preserving localization

## Executive verdict

Let

\[
 A=\operatorname {Cov}_tX,\qquad
 v=\operatorname {Cov}_t({\bf1}_S,X),\qquad
 q=v^TA^{-1}v.
\]

On the stopped region where \(A\) is positive definite and \(q>0\), use
the **driver covariance**

\[
 B=A^{-1}-{A^{-1}vv^TA^{-1}\over q},\qquad
 \Gamma=B+\alpha P_{v^\perp},\quad \alpha\ge0.           \tag{0.1}
\]

Thus the Brownian coefficient in the density equation is
\(\Gamma^{1/2}\), and the accumulated quadratic curvature is
\(\int\Gamma\,dt\).  This convention is important: using \(B\) itself as
the Brownian coefficient would instead accumulate \(B^2\).

The controlled process has several exact dimension-free features.

1. \(\Gamma v=0\), so the set mass is preserved pathwise:
   \(g_t=\mu_t(S)=g_0\).
2. The covariance drift is

\[
 dA=dA^{\rm mart}
 -\left(A-{vv^T\over q}+\alpha AP_{v^\perp}A\right)dt.  \tag{0.2}
\]

Thus the normalized part removes all covariance except the single
whitened direction \(A^{-1}v\).
3. The log-volume has deterministic drift exactly

\[
 -(n-1)-\alpha\{\operatorname {tr}A-e^TAe\},
 \qquad e={v\over|v|},
\]

plus an additional nonpositive Ito correction.
4. The explained label variance \(q=v^TA^{-1}v\) is a bounded
submartingale.  Its drift is the square

\[
 \sum_k\left\|A^{-1/2}
 \{D c_k-\mathcal T(c_k)A^{-1}v\}\right\|^2,
 \qquad c_k=\Gamma^{1/2}e_k,                            \tag{0.3}
\]

and \(q\le g_0(1-g_0)\).  This is a new finite information budget.
5. On a stop where \(A\preceq MI\),

\[
 \Gamma\succeq(M^{-1}+\alpha)P_{v^\perp}.              \tag{0.4}
\]

Consequently any small eigenvalue of the accumulated endpoint curvature
forces \(v_t/|v_t|\) to lock near one line for almost the whole active
time.

These identities do **not** bound the exceptional posterior variance.
If a product law has one factor of variance \(L\), the set depends only on
that factor, and \(v\) points along it, then the controller never acts on
that factor and its variance remains \(L\).  Such states occur as
one-coordinate exponential tilts with \(L\) arbitrarily large.  Thus no
dimension-free exceptional-variance theorem follows from mass preservation,
\(q\le g(1-g)\), or the determinant drift.

Nor does the normalized term suppress the winner mechanism at the
infinitesimal isotropic state.  There,

\[
                         B=P_{v^\perp};
\]

more precisely, for \(A=I+H\),

\[
 B=P_{v^\perp}-P_{v^\perp}HP_{v^\perp}+O(\|H\|^2).     \tag{0.5}
\]

The product-exponential maximum therefore has exactly the same initial
\(\Theta((\log n)^2)\) angular quadratic-variation rate as under the
ordinary transverse controller.  This proves that covariance
normalization gives no leading infinitesimal suppression.  It does not by
itself prove an order-one displacement on the \(1/\log^2n\) scale, because
the coefficients evolve on that scale.  In the frozen-coefficient model,
such a displacement would consume only order \(q\) of the complete-square
budget.

The report is rigorous up to explicitly stated bounded stopping times.  At
\(v=0\) the adopted convention is \(\Gamma=0\), which gives a well-defined
frozen solution.  A nontrivial restart at \(v=0\) has a discontinuous
direction selector and is not asserted to possess a strong nonexplosive
solution.

## 1. Stopped construction and the meaning of the driver

Write the posterior in natural parameters as

\[
 p_t(x)=Z_t^{-1}
 \exp\left\{c_t\cdot x-\frac12x^TQ_tx\right\}p_0(x).
                                                               \tag{1.1}
\]

For a predictable positive semidefinite driver covariance \(\Gamma_t\),
put \(C_t=\Gamma_t^{1/2}\) and solve

\[
 dc_t=C_t\,dW_t+\Gamma_ta_t\,dt,\qquad
 dQ_t=\Gamma_tdt,\qquad a_t=\mathbb E_tX.              \tag{1.2}
\]

The density then satisfies

\[
                         dp_t(x)
 =p_t(x)\langle x-a_t,C_t\,dW_t\rangle.                \tag{1.3}
\]

Assume first that \(A_0\succ0\) and \(q_0>0\).  Fix \(R<\infty\) and stop
before any of the following occurs:

\[
 R^{-1}I\not\preceq A_t\not\preceq RI,\qquad
 q_t\notin[R^{-1},R],                                  \tag{1.4}
\]

or before \((c_t,Q_t)\) leaves a compact subset of the interior of the
posterior moment-generating domain.  On this stopped set all posterior
moments appearing below are bounded.  They are smooth functions of the
natural parameters, and the matrix in (0.1) has constant rank with its
nonzero eigenvalues bounded away from zero.  Its principal square root is
smooth there.  Standard finite-dimensional SDE theory therefore gives a
unique strong solution up to this stopping time.  This is the only
existence assertion made here.

At \(v=0\), set

\[
                         \Gamma=0.                     \tag{1.5}
\]

Starting from such a state, (1.2) is the frozen solution.  One may instead
try \(\Gamma=A^{-1}\) at the instant \(v=0\), because it still has
\(\Gamma v=0\), but the coefficient jumps to a direction-dependent
rank-\((n-1)\) matrix as soon as \(v\ne0\).  Neither pathwise uniqueness nor
nonexplosion for that restart is used or claimed.

For \(v\ne0\), define

\[
 H=A^{-1},\qquad
 q=v^THv,\qquad e={v\over|v|},\qquad
 P=I-ee^T,
\]

and use

\[
 B=H-{Hvv^TH\over q},\qquad
 \Gamma=B+\alpha P.                                    \tag{1.6}
\]

If

\[
                         w={A^{-1/2}v\over\sqrt q},
\]

then \(|w|=1\) and

\[
 B=A^{-1/2}(I-ww^T)A^{-1/2}\succeq0,\qquad Bv=0.       \tag{1.7}
\]

Thus \(\Gamma v=0\).

## 2. Exact moment SDEs and mass preservation

Put

\[
 Y=X-a_t,\qquad h={\bf1}_S-g_t,
\]

\[
 A=\mathbb EYY^T,\qquad
 v=\mathbb E[hY],\qquad
 D=\mathbb E[hYY^T],\qquad
 \mathcal T(z)=\mathbb E[YY^T\langle Y,z\rangle].      \tag{2.1}
\]

For the columns \(c_k=Ce_k\), write

\[
                         T_k=\mathcal T(c_k).
\]

The standard posterior differentiation, retained here with every drift,
gives

\[
 \boxed{
 \begin{aligned}
 da&=AC\,dW,\\
 dg&=v^TC\,dW,\\
 dv&=DC\,dW-A\Gamma v\,dt,\\
 dA&=\mathcal T(C\,dW)-A\Gamma A\,dt,\\
 dQ&=\Gamma\,dt.
 \end{aligned}}                                         \tag{2.2}
\]

Since \(\Gamma v=0\),

\[
 \boxed{
 dg=0,\qquad
 dv=DC\,dW,\qquad
 g_t=g_0
 }                                                       \tag{2.3}
\]

up to the stopping time.  Thus a central initial set survives with
probability one; there is no separate set-mass stopping event.

The covariance drift simplifies because

\[
                         ABA=A-{vv^T\over q}.           \tag{2.4}
\]

Hence

\[
 \boxed{
 dA=\mathcal T(C\,dW)
 -\left(A-{vv^T\over q}+\alpha APA\right)dt.
 }                                                       \tag{2.5}
\]

The normalized deterministic drift is positive semidefinite:

\[
 A-{vv^T\over q}
 =A^{1/2}(I-ww^T)A^{1/2}\succeq0.                     \tag{2.6}
\]

Its null line is \(\operatorname {span}\{A^{-1}v\}\), because

\[
 \left(A-{vv^T\over q}\right)A^{-1}v=0.                \tag{2.7}
\]

This differs from the instantaneous null line of the accumulated curvature,
which is \(\operatorname {span}\{v\}\).

Taking traces in (2.5) gives

\[
 \boxed{
 \begin{aligned}
 d\,\operatorname {tr}A
 &=\operatorname {tr}\mathcal T(C\,dW)\\
 &\quad-\left\{\operatorname {tr}A-{|v|^2\over q}
 +\alpha(\operatorname {tr}A^2-e^TA^2e)\right\}dt.
 \end{aligned}}                                         \tag{2.8}
\]

The term

\[
                         {|v|^2\over q}
 ={1\over e^TA^{-1}e}                                  \tag{2.9}
\]

is the harmonic Rayleigh variance in the active line.  In particular,
the drift of the ordinary Rayleigh variance \(e^TAe\), with \(e\) frozen,
contains

\[
                         -\left(e^TAe-
                         {1\over e^TA^{-1}e}\right)\le0. \tag{2.10}
\]

The covariance normalization removes correlation-induced excess variance,
but it leaves an exact exception when \(e\) is an eigenvector of \(A\).

## 3. Exact log-volume collapse

Matrix Ito applied to \(\log\det A\) yields

\[
 \begin{aligned}
 d\log\det A
 &=\sum_k\operatorname {tr}(HT_k)dW_k
 -\operatorname {tr}(A\Gamma)dt\\
 &\quad-\frac12\sum_k
 \operatorname {tr}(HT_kHT_k)dt.                      \tag{3.1}
 \end{aligned}
\]

Now

\[
 \operatorname {tr}(AB)
 =\operatorname {tr}\left(I-{vv^TH\over q}\right)
 =n-1,                                                 \tag{3.2}
\]

and

\[
 \operatorname {tr}(AP)=\operatorname {tr}A-e^TAe.
\]

Therefore

\[
 \boxed{
 \begin{aligned}
 d\log\det A
 &=\sum_k\operatorname {tr}(HT_k)dW_k\\
 &\quad-\left[n-1+\alpha(\operatorname {tr}A-e^TAe)\right]dt\\
 &\quad-\frac12\sum_k
 \|A^{-1/2}T_kA^{-1/2}\|_{HS}^2dt.
 \end{aligned}}                                         \tag{3.3}
\]

After localization, the martingale is a true martingale under the bounded
stop.  More precisely, if \(\sigma\) is any bounded stopping time on which
the preceding integrability hypotheses hold, then

\[
 \mathbb E\log\det A_{T\wedge\sigma}
 \le \log\det A_0-(n-1)\mathbb E(T\wedge\sigma).       \tag{3.4}
\]

Thus the expected log-volume falls at rate at least \(n-1\) per unit of
active stopped time.
This controls the product of covariance eigenvalues, not the largest
exceptional eigenvalue.

## 4. A complete-square SDE for the explained label variance

The covariance inequality

\[
                         {vv^T\over g(1-g)}\preceq A
                                                               \tag{4.1}
\]

implies

\[
                         0\le q=v^TA^{-1}v\le g(1-g).   \tag{4.2}
\]

Since \(g\) is constant, \(q\) has a fixed numerical upper bound for a
balanced set.

To derive its SDE, matrix Ito gives

\[
 \begin{aligned}
 dH
 &=-H\,dA\,H+H\,dA\,H\,dA\,H\\
 &=-\sum_kHT_kH\,dW_k
 +\left\{\Gamma+\sum_kHT_kHT_kH\right\}dt.             \tag{4.3}
\end{aligned}
\]

For each \(k\), put

\[
 n_k=Dc_k,\qquad
 L_k=n_k-T_kHv.                                        \tag{4.4}
\]

Applying Ito to \(q=v^THv\), including the cross-variation between \(v\)
and \(H\), gives

\[
 \boxed{
 \begin{aligned}
 dq
 &=\sum_k\left\{
 2v^THn_k-v^THT_kHv\right\}dW_k\\
 &\quad+\sum_k L_k^THL_k\,dt.
 \end{aligned}}                                         \tag{4.5}
\]

The otherwise present drift \(v^T\Gamma v\) is zero.  Thus \(q\) is a
bounded submartingale.  For any bounded stopping time \(\sigma\) before
the coefficient stop,

\[
 \boxed{
 \mathbb E\int_0^\sigma
 \sum_k\|A^{-1/2}
 \{Dc_k-\mathcal T(c_k)A^{-1}v\}\|^2dt
 =\mathbb E q_\sigma-q_0
 \le g_0(1-g_0)-q_0.
 }                                                       \tag{4.6}
\]

This identity is dimension free and exact.  It is the strongest new
information functional supplied by the normalized controller.

It does not control a small-signal initial rotation.  At \(A=I\), the
angular quadratic variation of \(e=v/|v|\) is

\[
 d[e]_t={PDC\,C^TDP\over|v|^2}dt
 ={PD\Gamma DP\over|v|^2}dt.                           \tag{4.7}
\]

For \(\alpha=0\), \(\Gamma=P\), so an order-one angular movement costs
only order \(|v|^2=q\) in (4.6).  When the initial Fisher signal is tiny,
the fixed budget \(g(1-g)-q_0\) does not prevent such a movement.

## 5. Accumulated curvature and the locked-line alternative

The normalized covariance has the variational identity

\[
 \theta^TB\theta
 =\min_{a\in\mathbb R}(\theta-av)^TA^{-1}(\theta-av).
                                                               \tag{5.1}
\]

Consequently, on a stopped interval with \(A\preceq MI\),

\[
 \boxed{
 B\succeq {1\over M}P_{v^\perp},\qquad
 \Gamma\succeq\left({1\over M}+\alpha\right)P_{v^\perp}.
 }                                                       \tag{5.2}
\]

Let

\[
 Q_T-Q_0=\int_0^T\Gamma_tdt,\qquad
 \beta=M^{-1}+\alpha.
\]

For every unit vector \(u\),

\[
 u^T(Q_T-Q_0)u
 \ge\beta\int_0^T
 \{1-|\langle u,e_t\rangle|^2\}dt.                     \tag{5.3}
\]

Thus if

\[
                         \lambda_{\min}(Q_T-Q_0)<\varepsilon,
\]

there is a unit \(u\) satisfying

\[
 \boxed{
 \int_0^T\|e_te_t^T-uu^T\|_{HS}^2dt
 <{2\varepsilon\over\beta}.
 }                                                       \tag{5.4}
\]

The driver either builds full curvature or locks its unique missing line.
This is the covariance-normalized analogue of the ordinary transverse
localization alternative.

At the endpoint,

\[
 d\mu_T(x)\propto
 \exp\left\{c_T\cdot x-\frac12x^TQ_Tx\right\}d\mu_0(x).
                                                               \tag{5.5}
\]

If \(d\mu_0=e^{-V_0}dx\), its potential is

\[
                         V_T(x)=V_0(x)-c_T\cdot x+
                         \frac12x^TQ_Tx.                \tag{5.6}
\]

Hence, in the distributional convexity sense,

\[
                         \nabla^2V_T\succeq Q_T.        \tag{5.7}
\]

If (5.4) selects an exceptional line \(u\), the standard anisotropic
endpoint argument controls all transverse directions by \(Q_T\) and leaves
only

\[
                         \operatorname {Var}_{\mu_T}
                         \langle u,X\rangle             \tag{5.8}
\]

to be bounded.  None of (2.8), (3.3), or (4.6) supplies that upper bound.

## 6. Exact one-factor obstruction to exceptional-variance control

Let

\[
                         \mu=\nu_L\otimes\rho
\]

on \(\mathbb R\times\mathbb R^{n-1}\), where \(\nu_L\) is a
one-dimensional log-concave probability of variance \(L\), \(\rho\) is a
centered product log-concave law, and

\[
                         S=\{x_1\ge m_L\}
\]

for a median \(m_L\) of \(\nu_L\).  Assume the covariance is diagonal.
Then

\[
                         v=\ell e_1,\qquad
                         A=\operatorname {diag}(L,A_\rho).
                                                               \tag{6.1}
\]

The normalized driver has block form

\[
                         B=
 \begin{pmatrix}
 0&0\\
 0&A_\rho^{-1}
 \end{pmatrix},
\qquad
 P_{v^\perp}=
 \begin{pmatrix}
 0&0\\0&I
 \end{pmatrix}.                                        \tag{6.2}
\]

The density martingale depends only on the transverse coordinates.
Product structure is preserved, the first marginal remains exactly
\(\nu_L\), and

\[
                         \operatorname {Var}_tX_1=L
                                                               \tag{6.3}
\]

for the entire stopped evolution.  Meanwhile all identities above hold,
including \(q=\ell^2/L\le1/4\) and the rate-\((n-1)\) determinant collapse.
Since \(L\) is arbitrary, these identities cannot bound (5.8).

This obstruction can occur inside the exponential family generated from an
isotropic product exponential.  If \(Z\) has rate one and is tilted by
\(e^{cZ}\), \(c<1\), then the tilted law is exponential of rate \(1-c\)
and variance

\[
                         L=(1-c)^{-2}.                 \tag{6.4}
\]

Thus a preceding informative or uncontrolled layer can create an
arbitrarily large one-coordinate variance.  Once the set centroid locks to
that coordinate, the mass-preserving normalized controller never acts on
it.  The example is not an isotropic starting counterexample: at the
original isotropic state \(L=1\).  It proves that a closing theorem must use
the history before locking, not merely the current covariance, mass
preservation, or the \(q\)-budget.

## 7. Initial-layer expansion and the winner test

Let \(A=I+H\), \(v=\ell e\), and \(P=I-ee^T\).  A direct first-order
expansion gives

\[
 A^{-1}=I-H+O(\|H\|^2),\qquad
 q=\ell^2(1-e^THe)+O(\ell^2\|H\|^2),
\]

and

\[
 \boxed{
 B=P-PHP+O(\|H\|^2).
 }                                                       \tag{7.1}
\]

In particular, at an isotropic state,

\[
                         \Gamma=(1+\alpha)P.            \tag{7.2}
\]

Therefore the initial SDEs are exactly those of ordinary
mass-preserving transverse localization, with time multiplied by
\(1+\alpha\):

\[
                         dv=D\sqrt{1+\alpha}\,P\,dW,
\]

\[
 {d\over dt}\operatorname {tr}[e]_0
 ={1+\alpha\over|v|^2}\|PDP\|_{HS}^2.                 \tag{7.3}
\]

### 7.1 Product exponential maximum

For iid rate-one exponentials \(Z_i\), centered and scaled to variance one,
take the balanced maximum event

\[
 S=\{\max_iZ_i\ge L_n\},\qquad
 (1-e^{-L_n})^n={1\over2}.
\]

At the isotropic initial state,

\[
 |v_0|\asymp{\log n\over\sqrt n},\qquad
 {1\over|v_0|^2}\|P_0D_0P_0\|_{HS}^2
 \asymp(\log n)^2.                                    \tag{7.4}
\]

Equations (7.2)--(7.3) show

\[
 \left.{d\over dt}\operatorname {tr}[e]_t\right|_0
 \asymp(1+\alpha)(\log n)^2.                           \tag{7.5}
\]

If the coefficients in (7.5) were frozen, its angular
quadratic-variation clock would become order one on the
\(1/\{(1+\alpha)\log^2n\}\) scale.  Equation (7.1) therefore proves an
infinitesimal no-go statement: at \(A=I\), covariance normalization does
not reduce the old controller's initial angular rate.  It does **not**
prove an order-one displacement or winner selection on that scale;
\(D_t\), \(v_t\), and \(A_t\) all evolve, and in particular the covariance
has martingale fluctuations rather than a deterministic \(O(t)\) change.
A finite-time conclusion requires a separate historical SDE estimate.

The complete-square budget is consistent with the frozen-coefficient
clock.  Since \(q_0\asymp(\log n)^2/n\), its candidate order-one angular
layer would cost only order \(q_0\), while (4.6) has a numerical total
budget.  This is not a proof that the actual process realizes that layer.
It shows only that the budget alone cannot exclude it.

### 7.2 Product sum and simplex

For the product-exponential sum half-set, permutation symmetry keeps
\(v_t\) on the diagonal line.  The normalized controller localizes all
transverse factors at the same relative covariance rate and does not create
a new exceptional direction.

For a regular isotropic simplex and a median cut in one vertex direction,
the vertex stabilizer likewise fixes one line.  At \(A_0=I\), the controller
again begins as \(P_{v^\perp}\); symmetry prevents winner switching.
Neither model refutes the process, but neither supplies a bound for a
general adaptively selected exceptional line.

More symmetric simplex cuts may have \(v_0=0\).  Under convention (1.5) the
process then freezes.  A nontrivial convention would require a separate
well-posed restart theorem.

### 7.3 Radial Gaussian and parity

For a Gaussian base, every posterior has

\[
                         A_t=(A_0^{-1}+Q_t)^{-1}
\]

and hence \(A_t\preceq A_0\).  Starting isotropically gives the universal
exceptional-variance bound one, independently of the controller.  Radial
sets and sign parity can have \(v_0=0\); with (1.5) they simply freeze at
the harmless isotropic covariance.

These examples show why a zero signal need not be forced into an arbitrary
direction.  They do not address non-Gaussian posterior variance creation.

## 8. What the normalized controller proves and what remains

The covariance-normalized process improves the algebra of ordinary
mass-preserving localization:

\[
 \begin{array}{c}
 ABA=A-vv^T/q,\\[1mm]
 d\,\mathbb E\log\det A/dt\le-(n-1),\\[1mm]
 q=v^TA^{-1}v\text{ has the complete-square budget (4.6)},\\[1mm]
 \Gamma\succeq(M^{-1}+\alpha)P_{v^\perp}
 \text{ on }A\preceq MI.
 \end{array}
\]

It rigorously reduces low endpoint curvature to one locked line.  It does
not bound the variance in that line, and its infinitesimal behavior at an
isotropic state is identical to the old transverse process.  The
product-factor obstruction shows that no unrestricted current-state
inequality can fill the gap.

A genuinely closing lemma would have to be historical:

> Starting from an isotropic law, the process cannot first create a
> variance \(L\gg1\) in a direction receiving normalized transverse
> curvature, then rotate \(v_t\) into that direction within the
> small-signal layer, and finally keep that line locked while preserving
> the set mass.

The complete-square \(q\)-budget is a plausible bookkeeping device for such
a theorem, but it is too weak at small \(q\): the frozen-coefficient
product-maximum angular clock can spend only order \(q_0\) on its candidate
initial layer.  No dimension-free historical inequality excluding this
sequence is proved here.  Consequently the proposed driver does not
presently close KLS; it fails to remove the winner obstruction at the
infinitesimal level, while the corresponding finite-time obstruction
remains unproved.
