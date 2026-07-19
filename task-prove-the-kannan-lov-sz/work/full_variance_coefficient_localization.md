# Projecting the full variance coefficient: exact drift and a scalar-correction no-go

## 0. Conclusion

Use the notation

\[
M=\mu_t(E),\quad a=\mathbb E_tX,\quad Z=X-a,\quad
A=\mathbb E_tZZ^T,
\]

\[
b=\operatorname{Cov}_t(\mathbf 1_E,X),\quad
\rho=|b|,\quad u=b/\rho,\quad P=I-uu^T,
\]

\[
D=\mathbb E_t[(\mathbf 1_E-M)ZZ^T],\quad
T=\mathbb E_tZ^{\otimes3},\quad q=T[u,u,\cdot].
\]

The exact martingale coefficient of

\[
r=u^TAu
\]

before applying the driver is the vector

\[
\xi=q+\frac2\rho DPAu.
\tag{0.1}
\]

Thus the strengthened driver is

\[
C=P_{\operatorname{span}\{b,\xi\}^{\perp}}.
\tag{0.2}
\]

It makes both \(M\) and \(r\) drift processes: \(Cb=C\xi=0\).  It still
has defect rank at most two, so the accumulated-curvature conclusion from
the previous report remains

\[
B_t\succeq \frac t3P_{F_t^\perp}
\]

for a future-selected two-plane \(F_t\).

The complete drift of \(r\) under (0.2) is

\[
\begin{aligned}
\Gamma_r={}&-|CAu|^2
-\frac2{\rho^2}\langle PAu,DCDu\rangle
-\frac r{\rho^2}\operatorname{tr}(PDCD)\\
&+\frac1{\rho^2}\operatorname{tr}(APDCDP)
+\frac2\rho\sum_j T[PDCe_j,u,Ce_j].
\end{aligned}
\tag{0.3}
\]

There is no sign.  In particular, at an isotropic state \(A=I\), (0.3)
reduces to

\[
\Gamma_r=\frac2\rho\sum_jT[PDCe_j,u,Ce_j],
\tag{0.4}
\]

which can be strictly positive.

Three natural scalar corrections have exact, informative drifts:

* \(\rho^2=|b|^2\) is a submartingale with drift
  \(\operatorname{tr}(DCD)\).
* The explained variance

  \[
  e=b^TA^{-1}b
  \]

  is a submartingale, with the exact square drift

  \[
  \Gamma_e=\sum_j
  \left\|A^{-1/2}\bigl(DCe_j-T(Ce_j)A^{-1}b\bigr)\right\|^2.
  \tag{0.5}
  \]

  Moreover \(0\le e\le M(1-M)\le1/4\).
* The covariance volume satisfies

  \[
  \Gamma_{\log\det A}
  =-\operatorname{tr}(CA)
  -\frac12\sum_j
  \left\|A^{-1/2}T(Ce_j)A^{-1/2}\right\|_{\mathrm{HS}}^2\le0.
  \tag{0.6}
  \]

The helpful signs would be \(r-\alpha\rho^2\), \(r-\beta e\), or
\(r+\gamma\log\det A\).  The first two fail quantitatively.  Section 6
constructs an explicit family of smooth isotropic log-concave states and
balanced Borel sets for which

\[
\Gamma_r=c_0>0,\qquad
\Gamma_{\rho^2}+\Gamma_e+d[\rho^2]/dt+d[e]/dt=O(\varepsilon^2).
\tag{0.7}
\]

Consequently no universal correction

\[
r+\phi(\rho^2,e)
\tag{0.8}
\]

with \(\phi\) having universally bounded first two derivatives can be a
supermartingale.  This includes every fixed linear combination of
\(|b|^2\) and \(b^TA^{-1}b\).

The favorable sign \(+\gamma\log\det A\) in (0.6) does not control \(r\):
\(\log\det A\) has no universal lower bound along log-concave posterior
states.  A supermartingale \(r+\gamma\log\det A\) therefore gives no upper
bound on \(r\).  Singular corrections such as \(-\gamma\log e\) are
coercive near \(b=0\), but their initial value is arbitrarily large even
for smooth isotropic measures and balanced sets; the same construction has
\(e=\Theta(\varepsilon^2)\).  Such a correction necessarily pays
\(\log(1/\varepsilon)\), so it cannot yield a universal KLS bound.

Thus projecting the full coefficient repairs the Ito error in the
\(q\)-driver but does not produce a dimension-free terminal Poincare
estimate.  The obstruction has moved from martingale noise to the unsigned
third-moment drift (0.3).

---

## 1. Exact moment equations

For a symmetric predictable contraction \(C_t\), stochastic localization
has density equation

\[
dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle.
\tag{1.1}
\]

When \(Cb=0\), the label mass is constant and the basic moment equations
are

\[
dM=0,\qquad da=AC\,dW,
\tag{1.2}
\]

\[
db=DC\,dW,
\tag{1.3}
\]

\[
dA=T[\cdot,\cdot,C\,dW]-AC^2A\,dt.
\tag{1.4}
\]

For an exact orthogonal projection, \(C^2=C\).  All formulas below are
first asserted before the stopping time \(\inf\{|b|\le\delta\}\).  At
\(b=0\), one may use

\[
u_\eta=\frac{b}{\sqrt{|b|^2+\eta}}
\]

and the ridge projector

\[
C_\eta=I-U_\eta(U_\eta^TU_\eta+\eta I_2)^{-1}U_\eta^T,
\qquad U_\eta=[b,\xi_\eta].
\tag{1.5}
\]

This is globally smooth on compact support.  As in the two-constraint
report,

\[
\operatorname{rank}(I-C_\eta^2)\le2,\qquad
\operatorname{tr}(I-C_\eta^2)\le2,
\tag{1.6}
\]

so every curvature estimate below is uniform in \(\eta\).  The exact
formulas are recovered on stopped constant-rank regions as
\(\eta\downarrow0\).

From (1.3), with

\[
Q=DCD
\]

for an exact projection, Ito's formula gives

\[
du=L\,dW+h\,dt,
\qquad L=\frac1\rho PDC,
\tag{1.7}
\]

\[
h=-\frac1{\rho^2}PQu
-\frac1{2\rho^2}u\,\operatorname{tr}(PQ).
\tag{1.8}
\]

---

## 2. The full coefficient and complete drift

For each Brownian coordinate, put

\[
L_j=Le_j,\qquad H_j=T[\cdot,\cdot,Ce_j].
\]

Expanding \(d(u^TAu)\) gives

\[
\begin{aligned}
dr={}&u^T(dA)u+2(Au)^Tdu
+du^TA\,du+2du^T(dA)u\\
={}&(C\xi)^T\,dW+\Gamma_r\,dt,
\end{aligned}
\tag{2.1}
\]

where

\[
\xi=q+\frac2\rho DPAu
\tag{2.2}
\]

and

\[
\Gamma_r=-|CAu|^2
+2\langle Au,h\rangle
+\sum_j\langle AL_j,L_j\rangle
+2\sum_j\langle L_j,H_ju\rangle.
\tag{2.3}
\]

Substituting (1.7)--(1.8) into (2.3), and using \(C^2=C\), yields the
fully expanded formula

\[
\boxed{
\begin{aligned}
\Gamma_r={}&-|CAu|^2
-\frac2{\rho^2}\langle PAu,DCDu\rangle
-\frac r{\rho^2}\operatorname{tr}(PDCD)\\
&+\frac1{\rho^2}\operatorname{tr}(APDCDP)
+\frac2\rho\sum_jT[PDCe_j,u,Ce_j].
\end{aligned}}
\tag{2.4}
\]

The condition \(C\xi=0\) removes the martingale term in (2.1).  It creates
no further algebraic identity in (2.4).

At an isotropic state, \(A=I\), one has \(r=1\), \(PAu=0\), and \(Cu=0\)
because \(Cb=0\).  The third and fourth terms of (2.4) cancel:

\[
-\rho^{-2}\operatorname{tr}(PDCD)
+\rho^{-2}\operatorname{tr}(PDCDP)=0.
\]

Thus

\[
\Gamma_r=\frac2\rho\sum_jT[PDCe_j,u,Ce_j].
\tag{2.5}
\]

This term has either sign.

---

## 3. The signal norm

Equation (1.3) gives

\[
d\rho^2=2b^TDC\,dW+\operatorname{tr}(DCD)\,dt.
\tag{3.1}
\]

Thus \(\rho^2\) is a submartingale.  On \(\rho>0\),

\[
d\log\rho^2
=\frac2\rho u^TDC\,dW
+\frac{\operatorname{tr}(DCD)-2u^TDCDu}{\rho^2}\,dt.
\tag{3.2}
\]

The logarithmic drift has no fixed sign.  The negative logarithm is
coercive as \(b\to0\), but its value is not uniformly bounded at the
initial time: a balanced set may have \(b=0\), or arbitrarily small
nonzero \(b\), even under a smooth isotropic log-concave measure.

---

## 4. Explained variance

Assume \(A\) is positive definite and put

\[
H=A^{-1},\qquad w=Hb,\qquad e=b^THb.
\]

For \(R_j=T[\cdot,\cdot,Ce_j]\), the inverse covariance equation obtained
from (1.4) is

\[
dH=-\sum_jHR_jH\,dW_j
+\left(C+\sum_jHR_jHR_jH\right)dt.
\tag{4.1}
\]

Applying Ito's formula to \(e=b^THb\) gives

\[
de=\zeta^T\,dW+\Gamma_e\,dt,
\tag{4.2}
\]

where

\[
\zeta=C(2Dw-q_w),\qquad q_w=T[w,w,\cdot],
\tag{4.3}
\]

and

\[
\boxed{
\Gamma_e=\sum_j
\left\|H^{1/2}(DCe_j-R_jw)\right\|^2\ge0.}
\tag{4.4}
\]

To verify (4.4), expand its square.  The three terms are respectively

\[
\operatorname{tr}(HDCD),\qquad
\sum_j\langle R_jw,HR_jw\rangle,\qquad
-2\sum_j\langle DCe_j,HR_jw\rangle,
\]

which are exactly the three Ito drift terms in \(d(b^THb)\); the apparent
term \(b^TCb\) vanishes because \(Cb=0\).

There is also a universal state bound.  The covariance matrix of the
joint random vector \((\mathbf 1_E,X)\) is positive semidefinite.  Taking
its Schur complement gives

\[
0\le e=b^TA^{-1}b\le\operatorname{Var}(\mathbf 1_E)
=M(1-M)\le\frac14.
\tag{4.5}
\]

Hence \(-\alpha e\) is a bounded-below correction and is the most plausible
linear scalar candidate.  Section 6 shows that its drift can be
arbitrarily too small to dominate (2.5).

For completeness,

\[
\Gamma_{\log e}
=\frac{\Gamma_e}{e}-\frac{|\zeta|^2}{2e^2}
\tag{4.6}
\]

whenever \(e>0\).  It has no universal sign.

---

## 5. Covariance volume

Applying Ito's formula to \(\log\det A\) in (1.4) gives

\[
d\log\det A
=\sum_j\operatorname{tr}(A^{-1}R_j)\,dW_j
+\Gamma_{\det}\,dt,
\tag{5.1}
\]

where

\[
\boxed{
\Gamma_{\det}
=-\operatorname{tr}(CA)
-\frac12\sum_j
\left\|A^{-1/2}R_jA^{-1/2}\right\|_{\mathrm{HS}}^2
\le0.}
\tag{5.2}
\]

The sign is attractive only for \(r+\gamma\log\det A\), \(\gamma>0\).
However, this scalar does not dominate \(r\): log-concave
exponential-quadratic posterior states can have
\(\log\det A\to-\infty\).  Therefore a supermartingale estimate for
\(r+\gamma\log\det A\) allows \(r\) to be as large as
\(-\gamma\log\det A\), with no dimension-free bound.  With the opposite
sign, \(-\log\det A\) is coercive but its drift is nonnegative and cannot
cancel a positive \(\Gamma_r\).

---

## 6. Smooth isotropic small-signal states

This section gives the promised no-go family for every bounded scalar
correction built from \(\rho^2\) and \(e\).

### 6.1 A smooth skew isotropic product

Let \(E_1,E_2\) be independent mean-one exponential variables, let
\(G_1,G_2\) be independent standard Gaussians, and fix a sufficiently
small \(\sigma>0\).  Define

\[
X_i=\frac{E_i-1+\sigma G_i}{\sqrt{1+\sigma^2}},
\qquad i=1,2.
\tag{6.1}
\]

Each \(X_i\) has a positive \(C^\infty\) log-concave density: it is a
Gaussian convolution of a log-concave exponential density.  It has mean
zero, variance one, and third central moment

\[
\tau=\mathbb E X_i^3=\frac2{(1+\sigma^2)^{3/2}}>0.
\tag{6.2}
\]

Let \(V\) be an independent standard Gaussian.  The law

\[
\mu=\mathcal L(X_1,X_2,V)
\tag{6.3}
\]

is a smooth, full-dimensional, isotropic log-concave probability measure
on \(\mathbb R^3\).

Put

\[
u=\frac{e_1+e_2}{\sqrt2},\qquad
y=\frac{e_1-e_2}{\sqrt2}.
\]

Let \(m_\sigma\) be the median of \(|X_1-X_2|\), and define the balanced
two-coordinate event

\[
H=\{|X_1-X_2|\ge m_\sigma\}.
\tag{6.4}
\]

It is invariant under interchanging \(X_1,X_2\).  At \(\sigma=0\),
\(|X_1-X_2|\) is exponential with mean one and median \(\log2\), while

\[
\mathbb E[X_1+X_2\mid |X_1-X_2|=r]=r-1.
\]

It follows directly that, for \(g_H=\mathbf 1_H-1/2\),

\[
\rho_H:=\mathbb E[g_H\,u\cdot X]
=\frac{\log2}{2\sqrt2}>0
\tag{6.5}
\]

at \(\sigma=0\), and

\[
k_H:=\mathbb E[g_H(y\cdot X)^2]
=\frac{(\log2)^2+2\log2}{4}>0.
\tag{6.6}
\]

Both quantities remain positive for all sufficiently small fixed
\(\sigma>0\) by dominated convergence.  Swap symmetry gives

\[
\mathbb E[g_HX]=\rho_Hu,\qquad
\mathbb E[g_HXX^T]\,y=k_Hy,
\tag{6.7}
\]

and the \(u\)-\(y\) matrix entry in the second identity is zero.

### 6.2 Diluting the label without changing the measure

We now construct a deterministic balanced Borel set whose degree-one and
degree-two correlations are an arbitrarily small multiple of those in
(6.7).

Choose an even Borel set \(A_0\subset\mathbb R\) for the Gaussian \(V\)
such that

\[
\mathbb P(V\in A_0)=\frac12,\qquad
\mathbb E[V^2\mathbf 1_{A_0}]=\frac12.
\tag{6.8}
\]

Such a set may be taken as a union of a central interval and two tails.
Existence follows continuously by varying their two endpoints; one
endpoint fixes the mass and the other moves the conditional second moment
from below one to above one.

For every sufficiently small \(\varepsilon>0\), nonatomicity permits
even sets

\[
S_+\subset A_0,\qquad S_-\subset A_0^c
\]

such that

\[
\mathbb P(V\in S_+)=\mathbb P(V\in S_-)=\frac\varepsilon2,
\qquad
\mathbb E[V^2\mathbf 1_{S_+}]
=\mathbb E[V^2\mathbf 1_{S_-}].
\tag{6.9}
\]

For an explicit construction, take two short symmetric intervals at two
radii in \(A_0\), two at bracketing radii in \(A_0^c\), and choose their
four lengths from the two linear equations in (6.9).  The central-plus-tail
form of \(A_0\) makes the required second-moment intervals overlap.

Define \(E_\varepsilon\subset\mathbb R^3\) as follows:

* if \(V\notin S_+\cup S_-\), use the baseline label
  \(\mathbf 1_{\{V\in A_0\}}\);
* if \(V\in S_+\cup S_-\), use the label \(\mathbf 1_H(X_1,X_2)\).

The equalities (6.8)--(6.9), the evenness of all four \(V\)-sets, and the
balance of \(H\) imply exactly

\[
\mu(E_\varepsilon)=\frac12,
\tag{6.10}
\]

\[
b_\varepsilon=\varepsilon\rho_Hu,
\tag{6.11}
\]

and, for

\[
D_\varepsilon
=\mathbb E[(\mathbf 1_{E_\varepsilon}-1/2)XX^T],
\]

\[
D_\varepsilon y=\varepsilon k_Hy,\qquad
D_\varepsilon e_3=0,\qquad
\langle u,D_\varepsilon y\rangle=0.
\tag{6.12}
\]

All unwanted \(V\)- and cross-moments cancel by (6.9) and evenness.

### 6.3 Drift under the full-coefficient driver

The measure, not merely its covariance, is fixed as \(\varepsilon\)
varies.  Its third tensor has only

\[
T_{111}=T_{222}=\tau
\]

among the first two coordinates.  Therefore

\[
q=T[u,u,\cdot]=\frac{\tau}{\sqrt2}u.
\tag{6.13}
\]

Since \(A=I\), one has \(PAu=0\).  The full coefficient (0.1) is just
\(\xi=q\), which is collinear with \(b_\varepsilon\).  Hence

\[
C=P_{u^\perp}.
\tag{6.14}
\]

The \(e_3\) direction makes no contribution to (2.5), and

\[
T[y,u,y]=\frac{\tau}{\sqrt2}.
\]

Equations (6.11)--(6.14) give the exact positive drift

\[
\boxed{
\Gamma_r
=\frac{2}{\varepsilon\rho_H}
T[\varepsilon k_Hy,u,y]
=\frac{\sqrt2\,\tau k_H}{\rho_H}
=:c_0>0.}
\tag{6.15}
\]

The constant is independent of \(\varepsilon\).

On the other hand, (3.1) gives

\[
\Gamma_{\rho^2}=O(\varepsilon^2).
\tag{6.16}
\]

At this isotropic state \(e=\rho^2=\varepsilon^2\rho_H^2\).  Formula
(4.4) gives

\[
\Gamma_e=O(\varepsilon^2),
\tag{6.17}
\]

and the martingale coefficients and quadratic variations of both
\(\rho^2\) and \(e\) are \(O(\varepsilon^2)\) and \(O(\varepsilon^4)\),
respectively.  All constants in (6.16)--(6.17) depend only on the fixed
smooth measure and the fixed event \(H\), not on \(\varepsilon\).

Let \(\phi\) be any \(C^2\) function on a neighborhood of
\([0,1/4]^2\) with finite derivative bound

\[
\|\nabla\phi\|_\infty+\|\nabla^2\phi\|_\infty=L<\infty.
\]

Ito's formula and (6.16)--(6.17) yield

\[
\left|\Gamma_{\phi(\rho^2,e)}\right|\le C L\varepsilon^2.
\tag{6.18}
\]

For every proposed universal \(\phi\), choose
\(\varepsilon<\sqrt{c_0/(2CL)}\).  Then

\[
\Gamma_{r+\phi(\rho^2,e)}\ge\frac{c_0}{2}>0.
\tag{6.19}
\]

Thus no member of this clearly specified family is a universal
supermartingale.  In particular, no fixed coefficient choice in

\[
r-\alpha|b|^2-\beta b^TA^{-1}b
\]

works for all smooth isotropic log-concave states and balanced Borel sets.

Singular choices evade the bounded-derivative hypothesis only by losing
the desired universal control.  For example,

\[
-\log e_\varepsilon
=2\log(1/\varepsilon)-\log\rho_H^2
\]

already diverges at the initial state.  A supermartingale containing this
term can at best return a bound depending on the arbitrary dilution
parameter \(\varepsilon\), not a KLS constant.

---

## 7. Curvature conclusion and final obstruction

The full-coefficient driver (0.2) has exactly the same spectral bookkeeping
as the \(q\)-driver.  If

\[
K_t=\int_0^tP_{\operatorname{span}\{b_s,\xi_s\}}\,ds,
\]

then

\[
\operatorname{tr}K_t\le2t.
\]

For its top-two eigenspace \(F_t\),

\[
B_t=tI-K_t\succeq\frac t3P_{F_t^\perp}.
\]

The generalized defect theorem therefore stops at

\[
C_P(\mu_t)
\le C\left(t^{-1}
+\|P_{F_t}A_tP_{F_t}\|_{\mathrm{op}}\right).
\tag{7.1}
\]

Projecting \(\xi\) makes the *current* scalar \(r_t\) a drift process, but
the unsigned drift (2.4) can be positive by a fixed amount even when the
signal \(b\) is arbitrarily small.  It also does not control the second
direction in the history-selected plane \(F_t\).

The exact identities (0.3), (0.5), and (0.6) exhaust the elementary scalar
corrections based on signal size, explained variance, and covariance
volume:

* bounded corrections from \(|b|^2\) and \(e\) are defeated by the
  small-signal family;
* logarithmic signal corrections begin with an arbitrarily large value;
* the useful sign of \(\log\det A\) destroys coercivity and gives no upper
  bound on \(r\).

Accordingly, the full-coefficient projection does not yield a
dimension-free terminal Poincare bound without a new estimate controlling
the third-moment drift and the adaptive two-plane covariance.
