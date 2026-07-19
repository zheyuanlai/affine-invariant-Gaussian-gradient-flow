# Two-constraint projected localization: exact Ito calculus and the remaining two-dimensional defect

## 0. Executive conclusion

Let

\[
M_t=\mu_t(E),\qquad
b_t=\operatorname{Cov}_{\mu_t}(\mathbf 1_E,X),\qquad
u_t=\frac{b_t}{|b_t|},
\]

and, writing \(a_t=\mathbb E_tX\) and \(Z=X-a_t\), let

\[
q_t=\operatorname{Cov}_{\mu_t}((u_t\cdot Z)^2,X)
=\mathbb E_t[(u_t\cdot Z)^2Z].
\tag{0.1}
\]

The proposed exact driver is the orthogonal projection

\[
C_t=P_{\operatorname{span}\{b_t,q_t\}^{\perp}}.
\tag{0.2}
\]

It has one valid cancellation and one invalid advertised cancellation.

1. Since \(C_tb_t=0\), the label mass is exactly preserved:
   \(M_t=M_0\).
2. Since \(C_tq_t=0\), the stochastic coefficient of
   \(u^TAu\) is killed **only while \(u\) is frozen**.  The direction
   \(u=b/|b|\) is itself a semimartingale.  Its Ito motion contributes the
   additional martingale coefficient

   \[
   \frac{2}{|b|}\,C D P_{u^\perp}Au,
   \qquad
   D=\mathbb E[(\mathbf 1_E-M)ZZ^T].
   \tag{0.3}
   \]

   There is no identity forcing (0.3) to vanish.  Section 4 gives an
   explicit three-dimensional Gaussian state and a balanced sign partition
   for which \(q=0\) but (0.3) is nonzero.  Thus the claimed instantaneous
   variance preservation is false before any terminal estimate is attempted.

The curvature bookkeeping is nevertheless exact.  Since the kernel in
(0.2) has dimension at most two,

\[
B_t=\int_0^tC_s^2ds=tI-K_t,qquad
\operatorname{tr}K_t\le2t.
\]

If \(F_t\) is the span of the top two eigenvectors of \(K_t\), then

\[
B_t\succeq\frac t3P_{F_t^\perp}.
\tag{0.4}
\]

The audited defect-subspace theorem consequently gives

\[
C_P(\mu_t)
\le 2C_P((P_{F_t})_\#\mu_t)+\frac9t
\le C\left(t^{-1}+\|P_{F_t}A_tP_{F_t}\|_{\mathrm{op}}\right).
\tag{0.5}
\]

The second inequality uses only the fixed two-dimensional log-concave
Poincare bound.  Thus the one-dimensional adaptive-survivor problem has
been replaced by an adaptive **two-plane** covariance problem.  Neither
\(u_t\) nor the instantaneous plane \(\operatorname{span}\{b_t,q_t\}\)
need equal the history-selected plane \(F_t\).

Section 7 gives a log-concave terminal-state obstruction to every proposed
pointwise completion of (0.5).  For the isotropic ball and a balanced
halfspace label there are exponential-quadratic posterior states satisfying

\[
\operatorname{span}\{b,q\}=F,qquad C=P_{F^\perp},qquad
B=P_{F^\perp},qquad M=\frac12,
\]

but

\[
\|P_FA P_F\|_{\mathrm{op}}\ge c\sqrt n,qquad
C_P(\mu_t)\ge c\sqrt n.
\tag{0.6}
\]

The direct \(q\)-coefficient and, in this symmetric state, the rotation
coefficient both vanish; the protected variance is already large.  Hence
the two constraints do not yield a dimension-free terminal Poincare bound.
An additional probabilistic theorem excluding states such as (0.6) with
sufficient probability would still be needed.  Such a theorem is not a
consequence of the two cancellations or of the eigenvalue count.

---

## 1. Rigorous exact and regularized drivers

### 1.1 Exponential-quadratic localization

Let \(\mu\) be a log-concave probability measure on its affine support and
let

\[
\frac{d\mu_t}{d\mu}(x)
=Z_t^{-1}\exp\{\langle c_t,x\rangle-\tfrac12\langle B_tx,x\rangle\}.
\]

For a predictable symmetric contraction \(C_t\), use

\[
dc_t=C_t\,dW_t+C_t^2a_t\,dt,
\qquad dB_t=C_t^2dt.
\tag{1.1}
\]

After the usual moment stopping, the density satisfies

\[
dp_t(x)=p_t(x)\langle x-a_t,C_t\,dW_t\rangle.
\tag{1.2}
\]

All calculations below are first made before a stopping time on which the
displayed moments are finite.  For compactly supported measures they hold
globally for every smooth regularized driver.

### 1.2 Exact projection and its rank singularities

Set \(U=[b,q]\), regarded as an \(n\times2\) matrix.  Formally,

\[
P_U=U(U^TU)^\dagger U^T,\qquad C=I-P_U,
\tag{1.3}
\]

where \(\dagger\) denotes the Moore--Penrose inverse.  Formula (1.3) is a
Borel function of \((b,q)\), but it is discontinuous when the rank of
\(U\) changes.  In particular, specifying (1.3) does not establish strong
existence or uniqueness through \(b=0\), through \(q\in\mathbb Rb\), or
through \(q=0\).  The exact Ito identities below are valid for every
semimartingale solution and, classically, on stopped constant-rank regions
such as

\[
|b|\ge\delta,qquad
|P_{b^\perp}q|\ge\delta
\tag{1.4}
\]

in the rank-two case.  A separate rank-one chart applies when
\(q\in\mathbb Rb\) throughout an open time interval.  It is not legitimate
to glue these charts across a rank change without an existence theorem.

### 1.3 A global smooth regularization

For \(\eta>0\), define

\[
u_\eta=\frac{b}{\sqrt{|b|^2+\eta}},
\qquad
q_\eta=\mathbb E[(u_\eta\cdot Z)^2Z],
\qquad
U_\eta=[b,q_\eta],
\tag{1.5}
\]

and

\[
C_\eta=I-U_\eta(U_\eta^TU_\eta+\eta I_2)^{-1}U_\eta^T.
\tag{1.6}
\]

This is smooth in all posterior moments, including at \(b=0\).  It is
symmetric and \(0\preceq C_\eta\preceq I\).  If \(\sigma\) is a singular
value of \(U_\eta\), the corresponding singular value of
\(C_\eta U_\eta\) is

\[
\frac{\eta\sigma}{\sigma^2+\eta}\le\frac{\sqrt\eta}{2}.
\tag{1.7}
\]

Consequently

\[
|C_\eta b|\le\frac{\sqrt\eta}{2},
\qquad |C_\eta q_\eta|\le\frac{\sqrt\eta}{2}.
\tag{1.8}
\]

For example, the mass martingale has quadratic variation at most
\(\eta t/4\), exactly as in the one-constraint soft projection.  On compact
support (1.6) gives a unique nonexplosive strong solution.  On unbounded
support one first stops the moments and then passes to weak limits.  The
eigenvalue count in Section 6 is uniform in \(\eta\).

The regularization does not remove the rotation term (0.3).  It merely
replaces the singular factors \(|b|^{-1}\) and \(|b|^{-2}\) by smooth
expressions for the nonunit vector \(u_\eta\).  Since the exact formulas are
most transparent for a unit vector, Sections 2--5 work on \(|b|>0\); (1.5)
is the rigorous approximation at zero.

---

## 2. Moment SDEs

Suppress the time subscript.  Put

\[
Z=X-a,\qquad g=\mathbf 1_E-M,
\]

and define

\[
A=\mathbb E[ZZ^T],qquad
D=\mathbb E[gZZ^T],qquad
T=\mathbb E[Z^{\otimes3}].
\tag{2.1}
\]

For \(y\in\mathbb R^n\), write

\[
\mathcal T(y)=T[\cdot,\cdot,y]
=\mathbb E[ZZ^T\langle Z,y\rangle].
\]

Equation (1.2) gives

\[
dM=b^TC\,dW,
\qquad
da=AC\,dW,
\tag{2.2}
\]

and the standard covariance equation

\[
dA=\mathcal T(C\,dW)-AC^2A\,dt.
\tag{2.3}
\]

For the exact driver \(Cb=0\), (2.2) gives \(dM=0\).  Since

\[
b=\mathbb E[\mathbf 1_EX]-Ma,
\]

and \(M\) is constant, direct application of (1.2) gives the particularly
simple equation

\[
db=DC\,dW.
\tag{2.4}
\]

There is no drift in (2.4).  The cancellation of the apparently present
term \(ab^TC\,dW\) uses \(Cb=0\).

For reference, let \(\mathcal K_4\) be the fourth cumulant tensor

\[
(\mathcal K_4)_{ijkl}
=\mathbb E[Z_iZ_jZ_kZ_l]
-(A_{ij}A_{kl}+A_{ik}A_{jl}+A_{il}A_{jk}).
\]

Then the exact third-cumulant equation is

\[
dT=\mathcal K_4[\cdot,\cdot,\cdot,C\,dW]-\mathcal J\,dt,
\tag{2.5}
\]

where

\[
\mathcal J_{ijk}
=\sum_\ell\bigl((AC^2)_{i\ell}T_{\ell jk}
+(AC^2)_{j\ell}T_{\ell ik}
+(AC^2)_{k\ell}T_{\ell ij}\bigr).
\tag{2.6}
\]

One quick verification of (2.5) is to expand the conditional cumulant
generating function.  If \(K_t(\theta)=\log\mathbb E_t e^{\langle
\theta,X\rangle}\), then

\[
dK_t(\theta)
=(a_t(\theta)-a_t)^TC\,dW
-\frac12|C(a_t(\theta)-a_t)|^2dt.
\]

The cubic coefficient is exactly (2.5)--(2.6).

---

## 3. Exact Ito motion of the normalized signal

Let

\[
s=|b|,\qquad u=b/s,\qquad P=I-uu^T,\qquad G=DC,
\qquad Q=GG^T=DC^2D.
\tag{3.1}
\]

On \(s>0\), applying Ito's formula to \(b\mapsto b/|b|\) and using
(2.4) gives

\[
du=L\,dW+h_u\,dt,
\tag{3.2}
\]

where

\[
L=\frac1sPG=\frac1sPDC,
\tag{3.3}
\]

and

\[
h_u=-\frac1{s^2}PQu
-\frac1{2s^2}u\,\operatorname{tr}(PQ).
\tag{3.4}
\]

Formula (3.4) includes both tangential Ito drift and the radial drift that
keeps \(|u|=1\).  In particular, setting \(Cb=0\) does not make \(u\)
constant.  It would be constant only under the additional tensor identity
\(PDC=0\).

---

## 4. The protected variance SDE and an explicit rotation counterexample

Define

\[
r=u^TAu=\operatorname{Var}_{\mu_t}\langle X,u\rangle,
\qquad
q=T[u,u,\cdot].
\]

For each coordinate Brownian direction let

\[
L_j=Le_j,qquad H_j=\mathcal T(Ce_j).
\]

Expanding \(d(u^TAu)\), including both quadratic covariations
\(du^TAdu\) and \(2du^TdA\,u\), gives

\[
dr=\left(Cq+\frac2sCDPAu\right)^T dW+\Gamma\,dt,
\tag{4.1}
\]

with the exact drift

\[
\Gamma
=-|CAu|^2+2\langle Au,h_u\rangle
+\sum_j\langle AL_j,L_j\rangle
+2\sum_j\langle L_j,H_ju\rangle.
\tag{4.2}
\]

The first term in (4.1) is the frozen-direction coefficient.  The proposed
second constraint makes \(Cq=0\).  The second term is caused by the Ito
motion of \(u\), and it is not killed by \(Cb=Cq=0\).

### A concrete log-concave state where the rotation term is nonzero

Let \(X_1,X_2,X_3\) be independent centered Gaussians with standard
deviations \((1,2,3)\).  This is a smooth full-dimensional log-concave
state.  Let \(\epsilon_i=\operatorname{sign}X_i\), and let \(E\) consist of
the following four of the eight sign patterns:

\[
(-,-,-),\quad(-,-,+),\quad(-,+,-),\quad(+,-,+).
\tag{4.3}
\]

Every sign pattern has probability \(1/8\), so \(M=1/2\).  If
\(\sigma=2\mathbf 1_E-1\), direct averaging over the eight signs gives

\[
\mathbb E(\sigma\epsilon_1,\sigma\epsilon_2,
\sigma\epsilon_3)=(-1/2,-1/2,0),
\tag{4.4}
\]

and the only nonzero pair Fourier coefficients are

\[
\mathbb E[\sigma\epsilon_1\epsilon_3]=1/2,
\qquad
\mathbb E[\sigma\epsilon_2\epsilon_3]=-1/2.
\tag{4.5}
\]

Since magnitudes and signs are independent for centered Gaussians, (4.4)
implies

\[
b\parallel(-1,-2,0),
\qquad u=\frac{(-1,-2,0)}{\sqrt5}.
\]

All third central moments of a Gaussian vanish, so \(q=0\).  Put
\(m=\sqrt{2/\pi}\).  Equations (4.5) give

\[
D_{13}=D_{31}=\frac{3m^2}{4},
\qquad
D_{23}=D_{32}=-\frac{3m^2}{2},
\tag{4.6}
\]

with all other entries zero.  Also

\[
A=\operatorname{diag}(1,4,9),
\qquad
PAu=\frac1{5\sqrt5}(12,-6,0).
\]

Therefore

\[
CDPAu=\frac{18m^2}{5\sqrt5}e_3\ne0,
\tag{4.7}
\]

where \(C=P_{b^\perp}\) because \(q=0\).  Thus (4.1) has a nonzero
martingale coefficient even though \(Cq=0\).  This is an exact
log-concave counterexample to the frozen-eigenvector calculation.

Killing (4.7) would require adding the third kernel vector
\(DPAu\).  That would increase the curvature defect to dimension three.
Moreover, differentiating the new vector introduces fourth and fifth
moments, beginning an unclosed constraint hierarchy rather than preserving
a fixed one-dimensional mode.

---

## 5. The SDE of the second constraint

Although the failure is already visible in (4.1), the motion of \(q\) is
useful for checking well-posedness.  Let

\[
\Theta_j=\mathcal K_4[\cdot,\cdot,\cdot,Ce_j]
\]

be the noise coefficient of \(T\) in Brownian direction \(j\).  Applying
Ito's formula to \(q=T[u,u,\cdot]\), using (2.5) and (3.2), gives

\[
\begin{aligned}
dq={}&\left(\mathcal K_4[u,u,\cdot,C\,dW]
+2T[u,L\,dW,\cdot]\right)\\
&+\Bigl(-\mathcal J[u,u,\cdot]
+2T[u,h_u,\cdot]
+\sum_jT[L_j,L_j,\cdot]
+2\sum_j\Theta_j[u,L_j,\cdot]\Bigr)dt.
\end{aligned}
\tag{5.1}
\]

Thus even the plane \(\operatorname{span}\{b,q\}\) is driven by fourth
cumulants and by covariations with the rotating direction.  Formula (5.1)
also exposes the rank-change problem in the exact projection: a vector
\(q\) of arbitrarily small nonzero magnitude removes an entire direction
from the noise, whereas \(q=0\) removes none.  The ridge driver (1.6) is a
genuine SDE regularization; merely assigning a value to the projection at
\(q=0\) is not.

---

## 6. Accumulated curvature and the generalized defect theorem

For the exact projection, \(C^2=C\) and

\[
B_t=\int_0^tC_sds
=tI-K_t,qquad
K_t=\int_0^tP_{\operatorname{span}\{b_s,q_s\}}ds.
\tag{6.1}
\]

Every integrand in \(K_t\) has rank at most two and operator norm at most
one.  Hence

\[
K_t\succeq0,qquad \operatorname{tr}K_t\le2t.
\tag{6.2}
\]

Let \(\lambda_1\ge\lambda_2\ge\lambda_3\ge\cdots\) be its eigenvalues.
Then

\[
\lambda_3\le\frac{\operatorname{tr}K_t}{3}\le\frac{2t}{3}.
\tag{6.3}
\]

Let \(F_t\) be a measurable top-two eigenspace, enlarged arbitrarily to
dimension two if \(K_t\) has rank less than two.  For \(x\perp F_t\),

\[
\langle B_tx,x\rangle
=t|x|^2-\langle K_tx,x\rangle
\ge\frac t3|x|^2.
\tag{6.4}
\]

For the ridge regularization, set

\[
K_t^\eta=\int_0^t(I-C_{\eta,s}^2)ds.
\]

The range of \(I-C_\eta^2\) lies in \(\operatorname{span}U_\eta\), its
rank is at most two, and its eigenvalues lie in \([0,1]\).  Thus
\(\operatorname{tr}K_t^\eta\le2t\), and (6.3)--(6.4) remain unchanged.

The posterior potential satisfies, distributionally,

\[
D^2V_t\succeq B_t\succeq\frac t3P_{F_t^\perp}.
\]

The audited defect-subspace theorem, with \(\kappa=t/3\), gives

\[
C_P(\mu_t)
\le2C_P((P_{F_t})_\#\mu_t)+\frac9t.
\tag{6.5}
\]

Every log-concave probability law in dimension at most two satisfies

\[
C_P(\nu)\le C\|\operatorname{Cov}(\nu)\|_{\mathrm{op}}
\tag{6.6}
\]

with a numerical constant; this is the fixed-dimensional log-concave
Cheeger/Poincare theorem.  Applying (6.6) to the marginal in (6.5) yields

\[
C_P(\mu_t)
\le C\left(t^{-1}
+\|P_{F_t}A_tP_{F_t}\|_{\mathrm{op}}\right).
\tag{6.7}
\]

This is the exact generalized defect conclusion.  The plane \(F_t\) is
selected from the full future history of the kernel planes.  The
law-of-total-variance identity controls a fixed deterministic plane, not
\(F_t\).  Equations (4.1) and (5.1) do not control it either.

---

## 7. A log-concave terminal-state obstruction

We now construct states satisfying all instantaneous two-constraint
geometry but having an unbounded weak-plane covariance.

Let \(n\ge3\), let

\[
R=\sqrt{n+2},\qquad F=\operatorname{span}\{e_1,e_2\},
\qquad B=P_{F^\perp},
\]

and, for \(h\in\mathbb R\), let

\[
d\nu_{n,h}(x)
=Z_{n,h}^{-1}\exp\left\{hx_2-\frac12|P_{F^\perp}x|^2\right\}
\mathbf 1_{\{|x|\le R\}}dx.
\tag{7.1}
\]

This is a full-dimensional log-concave probability measure and is an
exponential-quadratic posterior of the isotropic uniform measure on
\(RB_2^n\).  Take the balanced label

\[
E=\{x_1\ge0\}.
\tag{7.2}
\]

Reflection in \(x_1\) gives \(\nu_{n,h}(E)=1/2\), and all remaining
coordinate symmetries give

\[
b=b_1e_1,\quad b_1>0,\qquad u=e_1,\qquad
q=\operatorname{Cov}_{\nu_{n,h}}(x_1^2,X)=q_2(h)e_2.
\tag{7.3}
\]

We can choose an arbitrarily small \(h>0\) for which \(q_2(h)\ne0\).
Indeed, put

\[
m_n(h)=\mathbb E_{\nu_{n,h}}x_1^2.
\]

The compact support makes \(m_n\) real analytic, and

\[
m_n'(h)=\operatorname{Cov}_{\nu_{n,h}}(x_1^2,x_2)=q_2(h).
\tag{7.4}
\]

The function is not constant: as \(h\to+\infty\), the measure concentrates
at the exposed point \(Re_2\), so \(m_n(h)\to0\), while \(m_n(0)>0\).
Therefore \(m_n'\) cannot vanish on any interval.  Choose

\[
0<h_n<\frac1{100R}
\tag{7.5}
\]

outside its discrete zero set.  Equations (7.3)--(7.5) imply

\[
\operatorname{span}\{b,q\}=F,qquad C=P_{F^\perp}.
\tag{7.6}
\]

It remains to estimate the covariance.  At \(h=0\), write
\(r^2=x_1^2+x_2^2\) and \(d=n-2\).  Integrating the transverse coordinates
gives the radial density in \(F\)

\[
r\longmapsto
r\,\mathbb P\{Q_d\le d+4-r^2\},
\qquad 0\le r\le\sqrt{d+4},
\tag{7.7}
\]

where \(Q_d\) is chi-square with \(d\) degrees of freedom.  The chi-square
lower-tail bound

\[
\mathbb P\{Q_d\le d-s\}\le e^{-s^2/(4d)}
\]

and Berry--Esseen on intervals \(r\asymp d^{1/4}\) give

\[
c\sqrt d\le\mathbb E_{\nu_{n,0}}r^2\le C\sqrt d.
\tag{7.8}
\]

Rotational symmetry inside \(F\) yields

\[
\operatorname{Var}_{\nu_{n,0}}x_1
=\frac12\mathbb E r^2\asymp\sqrt n.
\tag{7.9}
\]

Condition (7.5) makes the density ratio between \(\nu_{n,h_n}\) and
\(\nu_{n,0}\) lie between two fixed constants: \(|h_nx_2|\le1/100\), and
the normalizers obey the same comparison.  Reflection still makes
\(\mathbb E x_1=0\).  Hence

\[
\operatorname{Var}_{\nu_{n,h_n}}x_1\ge c\sqrt n.
\tag{7.10}
\]

At the parameter state (7.1),

\[
B=I-P_F,qquad K=I-B=P_F,
\]

so \(F\) is exactly the history-defect plane appearing in (6.7).  Moreover,
the driver plane determined by the two constraints is also exactly \(F\)
by (7.6).  Both frozen stochastic coefficients are zero.  Symmetry makes
\(Au\in\mathbb Ru\), so the rotation coefficient (0.3) vanishes at this
particular state as well.  Nevertheless (7.10) holds.  Testing the
Poincare inequality with \(f(x)=x_1\) also gives

\[
C_P(\nu_{n,h_n})\ge\operatorname{Var}_{\nu_{n,h_n}}x_1
\ge c\sqrt n.
\tag{7.11}
\]

Thus no pointwise terminal estimate with a universal constant follows from
mass preservation, frozen-variance preservation, and the rank-two
curvature count.  The example is a terminal **state** obstruction, not a
claim that a specified exact rank-changing feedback reaches this state with
positive probability.  In fact, (5.1) shows why rank changes require a
regularization before that dynamic question is even well posed.  For the
ridge driver, choosing \(\eta\ll|q_2(h_n)|^2\) makes the driver and curvature
arbitrarily close to (7.6), so the state obstruction persists uniformly as
the regularization is removed.

---

## 8. Final assessment

The two-constraint idea produces a correct deterministic improvement in the
curvature geometry: all but a two-dimensional subspace acquire curvature
at least \(t/3\).  It does not produce a dimension-free Poincare bound.
There are three independent gaps.

1. The normalized direction moves, and the exact missing martingale term is
   (0.3).  The Gaussian sign-pattern example makes it nonzero.
2. Even if a third constraint kills (0.3), the accumulated defect plane is
   selected from the entire path and need not coincide with the current
   protected directions.
3. Even in a state where the instantaneous plane, the accumulated defect
   plane, and both killed coefficients agree exactly, the log-concave ball
   posterior (7.1) has weak-plane variance and Poincare constant of order
   \(\sqrt n\).

Therefore the scheme stops rigorously at (6.7).  Completing it would require
a new dimension-free probability estimate showing that large-covariance
states of the type in Section 7 have sufficiently small weight along the
regularized stochastic flow.  Assuming such an estimate would simply
reintroduce the adaptive covariance obstruction in a two-dimensional form.
