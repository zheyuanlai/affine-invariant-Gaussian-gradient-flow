# Bounded projection localization: rigorous identities and an adaptive-survivor obstruction

## 0. Executive conclusion

Let \(\mu_t\) be a stochastic-localization posterior and let

\[
b_t=\operatorname{Cov}_{\mu_t}(\mathbf 1_E,X).
\]

The formal projection choice

\[
u_t=\frac{b_t}{|b_t|},\qquad C_t=I-u_tu_t^T
\tag{0.1}
\]

has two exact attractive properties as long as \(b_t\ne0\):

* the set mass \(\mu_t(E)\) is constant;
* the accumulated quadratic potential is

  \[
  B_t=tI-\int_0^t u_su_s^Tds.
  \]

The discontinuity at \(b=0\) is a genuine well-posedness issue, not a
notational issue. A Borel rule at zero makes the control adapted but does not
by itself supply existence or uniqueness for the degenerate feedback SDE.
The smooth replacement

\[
C_\varepsilon(b)=I-\frac{bb^T}{|b|^2+\varepsilon}
\tag{0.2}
\]

does resolve it. It preserves the set mass up to quadratic variation at most
\(\varepsilon t/4\), and retains the same one-dimensional spectral defect.
The stopped SDE and all limiting statements are proved below.

For either the exact process away from \(b=0\), or the soft process, a
deterministic spectral argument supplies a unit vector \(v_t\) such that

\[
B_t\succeq \frac t2 P_{v_t^\perp}.
\tag{0.3}
\]

The audited rank-defect Lichnerowicz lemma then gives, for a log-concave
posterior,

\[
C_P(\mu_t)
\le C\left(t^{-1}
+\operatorname{Var}_{\mu_t}\langle X,v_t\rangle\right).
\tag{0.4}
\]

Thus the remaining gate is precisely

\[
\mathbb E\left[
\left(t^{-1}+\operatorname{Var}_{\mu_t}
\langle X,v_t\rangle\right)^{-1/2}\right]\ge c.
\tag{0.5}
\]

There is no proof of (0.5) from isotropy and the generic localization
identities alone. Section 6 gives an explicit isotropic \(2n\)-point model,
with a balanced set and the exact feedback (0.1), for which

\[
\mathbb E\left[
\left(1+\operatorname{Var}_{\mu_1}
\langle X,v_1\rangle\right)^{-1/2}\right]
\le C\sqrt{\frac{\log n}{n}}\longrightarrow0.
\tag{0.6}
\]

This atomic measure is not log-concave. Therefore (0.6) is a rigorous
**mechanism/identity obstruction**, not a counterexample to KLS and not yet a
counterexample to the projection driver on the log-concave class. Section 7
identifies one exact feature that log-concavity forbids and explains why a
long-needle version of the obstruction remains the unresolved case.

---

## 1. Parameterized localization and stopped well-posedness

### 1.1 Exponential-quadratic posterior

Work intrinsically on the affine support \(H\) of \(\mu\). For
\(c\in H\) and a symmetric positive-semidefinite operator \(B\), set

\[
\frac{d\mu_{c,B}}{d\mu}(x)
=\frac{\exp(\langle c,x\rangle-\frac12\langle Bx,x\rangle)}
{Z(c,B)}.
\tag{1.1}
\]

Whenever the right-hand side is defined, write

\[
a(c,B)=\mathbb E_{c,B}X,\qquad
A(c,B)=\operatorname{Cov}_{c,B}(X),
\]

\[
M(c,B)=\mu_{c,B}(E),\qquad
b(c,B)=\operatorname{Cov}_{c,B}(\mathbf1_E,X).
\tag{1.2}
\]

For a predictable symmetric control \(C_t\), the parameter equations are

\[
dc_t=C_t\,dW_t+C_t^2a_t\,dt,\qquad
dB_t=C_t^2dt,\qquad (c_0,B_0)=(0,0).
\tag{1.3}
\]

Here and below all operators act on \(H\).

For the soft feedback, take

\[
C_t=C_\varepsilon(b_t)
=I-\frac{b_tb_t^T}{|b_t|^2+\varepsilon}.
\tag{1.4}
\]

The map in (1.4) is smooth, bounded, and has operator norm at most \(1\).
On any parameter region on which the moments in (1.2) are bounded, the maps
\((c,B)\mapsto a,A,M,b\) are locally Lipschitz, by differentiation under
the integral in (1.1). Standard finite-dimensional SDE theory therefore
gives a unique strong solution up to the exit time

\[
\tau_R=\inf\left\{s:\,
|c_s|+\|B_s\|+\mathbb E_s|X|^4>R\right\}\wedge R.
\tag{1.5}
\]

For a log-concave initial law, the moment generating function is finite in a
neighborhood of zero. At time zero the smallest eigenvalue of
\(C_\varepsilon(b_0)^2\) is strictly positive. By continuity,
\(B_s\succeq \delta s I\) for some random-free \(\delta>0\) on a sufficiently
small stopped interval. This positive quadratic term makes (1.1) integrable
for every finite \(c_s\). Thus the solution does not encounter an
instantaneous integrability boundary. The statements below are first made
for \(t\wedge\tau_R\); one may then let \(R\uparrow\infty\) whenever
non-explosion is available.

For the exact feedback (0.1), fix a deterministic unit vector \(e_0\) and use
the nonanticipating Borel rule

\[
\widehat u(b)=
\begin{cases}
b/|b|,&b\ne0,\\
e_0,&b=0.
\end{cases}
\tag{1.6}
\]

This specifies an adapted coefficient, but it is discontinuous at zero.
Strong well-posedness follows on the stopped regions
\(\{|b|\ge\delta\}\), where the coefficient is locally Lipschitz. It is not
legitimate to claim global strong well-posedness through \(b=0\) from (1.6)
alone. The soft flow is the rigorous repair. A stopped weak limit as
\(\varepsilon\downarrow0\) has exactly preserved set mass and agrees with
the exact projection whenever its limiting \(b\) is nonzero; on the zero set
its limiting covariance may be any soft limit of the form described in
Section 3. This relaxed zero-signal behavior is sufficient for the spectral
argument.

### 1.2 Density SDE

Let \(p_t=d\mu_t/d\mu\), with \(\mu_t=\mu_{c_t,B_t}\).
Ito's formula applied to (1.1)--(1.3) gives, up to every stopping time
\(\tau_R\),

\[
dp_t(x)
=p_t(x)\langle x-a_t,C_t\,dW_t\rangle.
\tag{1.7}
\]

There is no drift in (1.7). The drift \(C_t^2a_tdt\) in \(dc_t\) and the
quadratic-potential increment \(dB_t=C_t^2dt\) are exactly what cancel the
Ito and normalization drifts.

Consequently, for every integrable test function \(f\),

\[
d\mathbb E_t f
=\operatorname{Cov}_{\mu_t}(f,X)^TC_t\,dW_t.
\tag{1.8}
\]

After stopping, the nonnegative density martingale is a true mixture:

\[
\mathbb E\,\mu_{t\wedge\tau_R}=\mu.
\tag{1.9}
\]

Indeed each \(p_{t\wedge\tau_R}(x)\) is a nonnegative local martingale and
hence has expectation at most \(1\), while
\(\int\mathbb E p_{t\wedge\tau_R}\,d\mu=1\); equality follows almost
everywhere.

---

## 2. Exact and approximate preservation of the balanced mass

Let

\[
M_t=\mu_t(E),\qquad
b_t=\operatorname{Cov}_{\mu_t}(\mathbf1_E,X).
\]

Equation (1.8) gives

\[
dM_t=\langle b_t,C_t\,dW_t\rangle.
\tag{2.1}
\]

### 2.1 Exact projection

If \(b_t\ne0\) and \(C_t=I-u_tu_t^T\) with \(u_t=b_t/|b_t|\), then

\[
C_tb_t=0.
\]

Therefore, on every interval on which the exact feedback is well posed,

\[
M_t=M_0.
\tag{2.2}
\]

In particular, a balanced set stays exactly balanced.

### 2.2 Soft projection

For (1.4),

\[
C_\varepsilon(b)b
=\frac{\varepsilon}{|b|^2+\varepsilon}b,
\]

and hence

\[
|C_\varepsilon(b)b|^2
=\frac{\varepsilon^2|b|^2}{(|b|^2+\varepsilon)^2}
\le\frac\varepsilon4.
\tag{2.3}
\]

Thus, for \(M_0=1/2\),

\[
\mathbb E\left(M_{t\wedge\tau_R}-\frac12\right)^2
=\mathbb E\langle M\rangle_{t\wedge\tau_R}
\le\frac{\varepsilon t}{4}.
\tag{2.4}
\]

The all-times survival event is also quantitative. Doob's \(L^2\)
inequality gives

\[
\mathbb E\sup_{s\le t\wedge\tau_R}
\left|M_s-\frac12\right|^2
\le4\,\mathbb E\langle M\rangle_{t\wedge\tau_R}
\le\varepsilon t.
\tag{2.5}
\]

Consequently

\[
\mathbb P\left\{
\inf_{s\le t\wedge\tau_R}M_s<\frac14
\ \text{or}\
\sup_{s\le t\wedge\tau_R}M_s>\frac34
\right\}
\le16\varepsilon t.
\tag{2.6}
\]

At terminal time alone, (2.4) gives the sharper bound

\[
\mathbb P\left\{
\left|M_{t\wedge\tau_R}-\frac12\right|>\frac14
\right\}
\le4\varepsilon t.
\tag{2.7}
\]

More generally, Burkholder--Davis--Gundy gives, for every \(p\ge2\),

\[
\mathbb E\sup_{s\le t\wedge\tau_R}
\left|M_s-\frac12\right|^p
\le C_p(\varepsilon t)^{p/2}.
\tag{2.8}
\]

Thus a fixed sufficiently small \(\varepsilon\) gives a balanced survival
event with universal positive probability; taking
\(\varepsilon\downarrow0\) gives exact mass preservation in every stopped
weak limit.

---

## 3. The accumulated quadratic potential and its spectral defect

For \(b\ne0\), write \(u=b/|b|\) and

\[
\beta_\varepsilon(b)=\frac{|b|^2}{|b|^2+\varepsilon}.
\]

Then

\[
C_\varepsilon(b)=I-\beta_\varepsilon(b)uu^T
\]

and

\[
C_\varepsilon(b)^2
=I-\alpha_\varepsilon(b)uu^T,
\tag{3.1}
\]

where

\[
\alpha_\varepsilon(b)
=2\beta_\varepsilon(b)-\beta_\varepsilon(b)^2
=\frac{|b|^2(|b|^2+2\varepsilon)}
{(|b|^2+\varepsilon)^2}
\in[0,1).
\tag{3.2}
\]

At \(b=0\), set \(\alpha_\varepsilon=0\) and choose any predictable unit
\(u\); the product \(\alpha_\varepsilon uu^T\) is then zero. Integrating
\(dB_t=C_t^2dt\) gives

\[
B_t=tI-\int_0^t\alpha_su_su_s^Tds.
\tag{3.3}
\]

For the exact projection away from zero, \(\alpha_s=1\), so

\[
B_t=tI-\int_0^t u_su_s^Tds.
\tag{3.4}
\]

Define

\[
K_t=\int_0^t\alpha_su_su_s^Tds.
\]

It is positive semidefinite and

\[
\operatorname{tr}K_t=\int_0^t\alpha_sds\le t.
\tag{3.5}
\]

Let its eigenvalues be
\(\lambda_1\ge\lambda_2\ge\cdots\ge0\), and choose a measurable unit top
eigenvector \(v_t\). Since

\[
\lambda_1+\lambda_2\le\operatorname{tr}K_t\le t,
\]

one has

\[
\lambda_2\le\frac t2.
\tag{3.6}
\]

For every \(z\perp v_t\),

\[
\langle B_tz,z\rangle
=t|z|^2-\langle K_tz,z\rangle
\ge(t-\lambda_2)|z|^2
\ge\frac t2|z|^2.
\]

Thus

\[
B_t\succeq\frac t2P_{v_t^\perp}.
\tag{3.7}
\]

This is completely deterministic and contains no dimension-dependent
constant.

The protected eigenvalue is

\[
\rho_t:=\langle B_tv_t,v_t\rangle=t-\lambda_1.
\tag{3.8}
\]

For the soft flow,

\[
\rho_t=\int_0^t
\left(1-\alpha_s\langle u_s,v_t\rangle^2\right)ds
\]

\[
=\int_0^t(1-\alpha_s)ds
+\int_0^t\alpha_s|P_{v_t^\perp}u_s|^2ds.
\tag{3.9}
\]

Thus a small protected eigenvalue simultaneously requires little soft
leakage and \(L^2\)-angular concentration of the entire history around the
terminal line \(\mathbb Rv_t\). For every fixed \(\varepsilon>0\), the first
integrand is positive whenever \(b_s\) is finite, so \(B_t\) is positive
definite for \(t>0\), although (3.9) supplies no dimension-free lower bound
on \(\rho_t\).

---

## 4. Rank-one-defect reduction and boundary transfer

Assume now that \(\mu\) is log-concave. The terminal potential is

\[
V_t(x)=V(x)-\langle c_t,x\rangle+\frac12\langle B_tx,x\rangle.
\tag{4.1}
\]

In the distributional sense,

\[
D^2V_t\succeq B_t
\succeq\frac t2P_{v_t^\perp}.
\tag{4.2}
\]

Apply the audited rank-defect lemma with

\[
F=\mathbb Rv_t,\qquad E=v_t^\perp,\qquad\kappa=t/2.
\]

It gives

\[
C_P(\mu_t)
\le 2C_P\bigl((\langle v_t,\cdot\rangle)_{\#}\mu_t\bigr)
+\frac6t.
\tag{4.3}
\]

Every one-dimensional log-concave probability law \(\nu\) satisfies

\[
C_P(\nu)\le C_1\operatorname{Var}(\nu)
\tag{4.4}
\]

with a universal numerical constant. This follows, for example, from the
one-dimensional log-concave Cheeger bound and Cheeger's inequality.
Therefore

\[
C_P(\mu_t)
\le C\left(
t^{-1}+\operatorname{Var}_{\mu_t}\langle X,v_t\rangle
\right).
\tag{4.5}
\]

For a stopped mixture, (1.9) and Fatou's lemma give, for every Borel set \(E\),

\[
\mu^+(E)\ge\mathbb E\,\mu_{t\wedge\tau_R}^+(E).
\tag{4.6}
\]

The Buser--Ledoux implication for a log-concave posterior and (4.5) yield

\[
\mu_t^+(E)
\ge c\,
\min(M_t,1-M_t)
\left(t^{-1}+\operatorname{Var}_{\mu_t}
\langle X,v_t\rangle\right)^{-1/2}.
\tag{4.7}
\]

For the exact balanced flow this becomes

\[
\mu^+(E)\ge
c\,\mathbb E\left[
\left(t^{-1}+\operatorname{Var}_{\mu_t}
\langle X,v_t\rangle\right)^{-1/2}\right].
\tag{4.8}
\]

For the soft flow, put

\[
G_t=\{|M_t-1/2|\le1/4\}.
\]

Since the inverse-square-root factor is at most \(\sqrt t\), (2.7) gives

\[
\begin{aligned}
\mu^+(E)
&\ge\frac c4\,\mathbb E\left[
\mathbf1_{G_t}
\left(t^{-1}+\operatorname{Var}_{\mu_t}
\langle X,v_t\rangle\right)^{-1/2}\right]\\
&\ge\frac c4\left(
\mathbb E\left[
\left(t^{-1}+\operatorname{Var}_{\mu_t}
\langle X,v_t\rangle\right)^{-1/2}\right]
-4\varepsilon t^{3/2}\right).
\end{aligned}
\tag{4.9}
\]

Thus the soft mass error is harmless once the survivor inverse moment has a
universal lower bound. It does not help prove that lower bound.

---

## 5. What the protected-eigenvalue dichotomy does and does not prove

Let \(\rho_t\) be (3.8). If

\[
\rho_t\ge t/4,
\]

then (3.7) implies \(B_t\succeq(t/4)I\). Brascamp--Lieb gives

\[
\operatorname{Var}_{\mu_t}\langle X,v_t\rangle
\le \frac4t.
\tag{5.1}
\]

Hence the survivor gate is automatic on this event.

For the exact projection, a small \(\rho_t\) means

\[
\int_0^t|P_{v_t^\perp}u_s|^2ds=\rho_t.
\tag{5.2}
\]

For a fixed deterministic vector \(v\), the mixture identity gives the exact
law-of-total-variance bound

\[
\mathbb E\operatorname{Var}_{\mu_t}\langle X,v\rangle
\le\operatorname{Var}_{\mu}\langle X,v\rangle.
\tag{5.3}
\]

For isotropic \(\mu\), the right-hand side is \(1\). The obstruction is that
\(v_t\) is future-measurable and selected from the same localization path.
Replacing \(v_t\) by a fixed vector in (5.3) is invalid. A stopping-time
anchor would need a stability estimate paying quantitatively for changes of
the protected line. The next section shows that no such estimate follows
from the generic posterior identities alone.

---

## 6. A rigorous finite-state adaptive-survivor obstruction

### 6.1 Initial isotropic measure and balanced set

Let

\[
\Omega_n=\{x_i^+,x_i^-:1\le i\le n\},
\qquad
x_i^\pm=\pm\sqrt n\,e_i,
\]

and let \(\mu_n\) be uniform on these \(2n\) points. Then

\[
\mathbb E_{\mu_n}X=0,\qquad
\mathbb E_{\mu_n}XX^T=I.
\tag{6.1}
\]

Thus \(\mu_n\) is isotropic. Let

\[
E_n=\{x_i^+:1\le i\le n\}.
\tag{6.2}
\]

It is balanced.

Run the exact projection feedback. Let \(p_i^\pm(t)\) be the posterior
weights. Exact mass preservation gives

\[
\sum_i p_i^+(t)=\sum_i p_i^-(t)=\frac12.
\]

Define the conditional probability vectors

\[
q_i^\pm(t)=2p_i^\pm(t),\qquad
\sum_iq_i^\pm(t)=1.
\tag{6.3}
\]

Their conditional means are

\[
m_t^+=\sqrt n\,q_t^+,\qquad
m_t^-=-\sqrt n\,q_t^-,
\]

where a probability vector is identified with its coordinate vector.
Consequently

\[
b_t=\frac14(m_t^+-m_t^-)
=\frac{\sqrt n}{4}(q_t^++q_t^-).
\tag{6.4}
\]

Since the vector \(q_t^++q_t^-\) is nonnegative and has coordinate sum \(2\),

\[
|b_t|
\ge\frac{\sqrt n}{4}\frac2{\sqrt n}
=\frac12.
\tag{6.5}
\]

Thus this example never encounters \(b=0\). The exact feedback is smooth
along the entire finite-time path, and the finite support makes all
coefficients bounded. Hence the SDE is globally well posed on every finite
time interval. Moreover

\[
u_t=\frac{q_t^++q_t^-}{|q_t^++q_t^-|}
\tag{6.6}
\]

has nonnegative coordinates.

### 6.2 Conditional posterior SDE

Let \(P_t=I-u_tu_t^T\). Equation (1.7), together with
\(P_t(m_t^+-a_t)=P_t(2b_t)=0\), gives

\[
dq_i^+
=\sqrt n\,q_i^+
\langle P_t(e_i-q_t^+),dW_t\rangle.
\tag{6.7}
\]

Similarly,

\[
dq_i^-
=-\sqrt n\,q_i^-
\langle P_t(e_i-q_t^-),dW_t\rangle.
\tag{6.8}
\]

Each coordinate is a bounded martingale and each probability vector remains
in the simplex.

### 6.3 Entropy collapse

For a probability vector \(q\), write

\[
\mathcal H(q)=-\sum_iq_i\log q_i,
\qquad
R(q)=1-\sum_iq_i^2.
\tag{6.9}
\]

Let

\[
S(q)=\operatorname{diag}(q)-qq^T.
\]

Ito's formula applied to (6.7) or (6.8) yields

\[
\frac d{dt}\mathbb E\mathcal H(q_t)
=-\frac n2\,
\mathbb E\operatorname{tr}(P_tS(q_t)).
\tag{6.10}
\]

For any nonnegative unit vector \(u\),

\[
u^TS(q)u
=\sum_{i<j}q_iq_j(u_i-u_j)^2
\le\sum_{i<j}q_iq_j
=\frac12R(q),
\tag{6.11}
\]

because \((u_i-u_j)^2\le u_i^2+u_j^2\le1\). Since
\(\operatorname{tr}S(q)=R(q)\), equations (6.10)--(6.11) give

\[
\frac d{dt}\mathbb E\mathcal H(q_t)
\le-\frac n4\mathbb ER(q_t).
\tag{6.12}
\]

Initially \(q_0\) is uniform, so \(\mathcal H(q_0)=\log n\). Hence

\[
\int_0^t\mathbb ER(q_s)ds
\le\frac{4\log n}{n}.
\tag{6.13}
\]

Because \(\sum_iq_i(t)^2\) is a submartingale,
\(\mathbb ER(q_t)\) is nonincreasing. Therefore, for every \(t>0\),

\[
\mathbb ER(q_t)
\le\frac{4\log n}{nt}.
\tag{6.14}
\]

This holds separately for \(q^+\) and \(q^-\).

For \(0<s\le t\), the vector-martingale identity gives

\[
\mathbb E|q_t-q_s|^2
=\mathbb E|q_t|^2-\mathbb E|q_s|^2
\le\mathbb ER(q_s)
\le\frac{4\log n}{ns}.
\tag{6.15}
\]

### 6.4 Stabilization of the protected direction

Put

\[
r_s=q_s^++q_s^-,
\qquad u_s=r_s/|r_s|.
\]

If

\[
R(q_s^+)+R(q_s^-)\le\frac12,
\tag{6.16}
\]

then, since \(q_s^+\cdot q_s^-\ge0\),

\[
|r_s|^2
\ge|q_s^+|^2+|q_s^-|^2
\ge\frac32.
\tag{6.17}
\]

On pairs of times for which (6.16) holds, normalization is uniformly
Lipschitz. On its complement, use \(|u_s-u_1|\le2\) and Markov's inequality.
Equations (6.14)--(6.15) therefore imply, for a universal \(C\),

\[
\mathbb E|u_s-u_1|^2
\le C\min\left\{1,\frac{\log n}{ns}\right\},
\qquad0<s\le1.
\tag{6.18}
\]

For explicit bookkeeping, one may take \(C=192\). Indeed the probability
that (6.16) fails at either \(s\) or \(1\) is at most
\(32\log n/(ns)\); on the complementary event, squared Lipschitz
normalization and (6.15) contribute at most
\(64\log n/(ns)\), while the exceptional event contributes at most
\(128\log n/(ns)\).

For unit vectors,

\[
\|uu^T-vv^T\|_{\mathrm{op}}\le2|u-v|.
\]

Let

\[
K_1=\int_0^1u_su_s^Tds.
\]

Using (6.18), Cauchy--Schwarz, and splitting the integral at
\((\log n)/n\), one obtains

\[
\begin{aligned}
\mathbb E\|K_1-u_1u_1^T\|_{\mathrm{op}}
&\le2\int_0^1
\left(\mathbb E|u_s-u_1|^2\right)^{1/2}ds\\
&\le C\sqrt{\frac{\log n}{n}}.
\end{aligned}
\tag{6.19}
\]

The preceding explicit value permits \(C=248\) in (6.19).

Let \(v_1\) be a measurable top eigenvector of \(K_1\). On the event

\[
\|K_1-u_1u_1^T\|_{\mathrm{op}}\le\frac14,
\tag{6.20}
\]

the top eigenspace is one-dimensional and

\[
|\langle v_1,u_1\rangle|^2\ge\frac12.
\tag{6.21}
\]

Indeed the top eigenvalue is at least \(1-\delta\), while it is at most
\(|\langle v_1,u_1\rangle|^2+\delta\), with
\(\delta=\|K_1-u_1u_1^T\|_{\mathrm{op}}\). By (6.19), the failure
probability in (6.20) is at most

\[
C\sqrt{\frac{\log n}{n}}.
\tag{6.22}
\]

More explicitly, it is at most
\(992\sqrt{\log n/n}\).

### 6.5 Terminal survivor variance

The covariance decomposition across the two balanced halves is

\[
A_1
=\frac12A_1^+
+\frac12A_1^-
+\frac14(m_1^+-m_1^-)(m_1^+-m_1^-)^T
\succeq4b_1b_1^T.
\tag{6.23}
\]

By (6.4),

\[
|b_1|^2
=\frac n{16}|q_1^++q_1^-|^2
\ge\frac n{16}
\left(2-R(q_1^+)-R(q_1^-)\right).
\tag{6.24}
\]

Equations (6.14) and Markov's inequality show that

\[
R(q_1^+)+R(q_1^-)\le\frac12
\tag{6.25}
\]

except on an event of probability at most \(16\log n/n\). On the intersection
of (6.20) and (6.25),

\[
\begin{aligned}
\operatorname{Var}_{\mu_1}\langle X,v_1\rangle
&=v_1^TA_1v_1\\
&\ge4|\langle v_1,b_1\rangle|^2\\
&=4|b_1|^2|\langle v_1,u_1\rangle|^2\\
&\ge\frac{3n}{16}.
\end{aligned}
\tag{6.26}
\]

The inverse-square-root factor is at most \(1\). Combining
(6.22), (6.25), and (6.26) gives

\[
\mathbb E\left[
\left(1+\operatorname{Var}_{\mu_1}
\langle X,v_1\rangle\right)^{-1/2}\right]
\le
C\sqrt{\frac{\log n}{n}}
+\left(1+\frac{3n}{16}\right)^{-1/2}.
\tag{6.27}
\]

Before absorbing constants, the two exceptional probabilities in this
display are at most
\[
992\sqrt{\frac{\log n}{n}}
\quad\text{and}\quad
\frac{16\log n}{n},
\]
respectively.

After changing \(C\),

\[
\mathbb E\left[
\left(1+\operatorname{Var}_{\mu_1}
\langle X,v_1\rangle\right)^{-1/2}\right]
\le C\sqrt{\frac{\log n}{n}}\longrightarrow0.
\tag{6.28}
\]

This proves the announced identity-level obstruction.

For each fixed \(n\), the finite-state coefficients depend continuously on
the soft parameter \(\varepsilon\), and (6.5) keeps the exact trajectory away
from the singular set. Hence the same conclusion, with an arbitrarily small
error, holds for all sufficiently small \(\varepsilon=\varepsilon_n>0\).
Thus no lower bound uniform simultaneously in dimension and in the
zero-signal regularization parameter can be obtained from the generic
posterior identities.

Again, \(\mu_n\) is not log-concave. Nothing in Section 6 disproves (0.5) on
the log-concave class.

---

## 7. The minimal log-concavity input: binary convexification

The atomic failure concentrates each half of the posterior near a single
point. A one-dimensional marginal of a log-concave measure cannot do this
for two complementary labels. The following quantitative statement makes
that distinction exact.

### 7.1 Binary quantization lemma

**Lemma.** There is a universal \(c_0>0\) with the following property. Let
\(Z\) have a one-dimensional log-concave law with
\(\operatorname{Var}Z=\sigma^2>0\), and let \(Y\in\{0,1\}\) be any binary
random variable coupled with \(Z\). Then

\[
\mathbb E\operatorname{Var}(Z\mid Y)
\ge c_0\sigma^2.
\tag{7.1}
\]

No balance assumption on \(Y\) is needed.

**Proof.** A standard one-dimensional log-concavity estimate gives

\[
\|f_Z\|_\infty\le \frac{C_0}{\sigma}
\tag{7.2}
\]

for a universal \(C_0\). For completeness, normalize
\(\|f_Z\|_\infty=1\) and choose a mode \(m\). The superlevel interval
\(\{f_Z\ge1/2\}\) has length at most \(2\). Concavity of \(\log f_Z\)
then bounds either tail beyond this interval by an exponential with rate at
least \((\log2)/2\). Hence
\(\mathbb E(Z-m)^2\le C_0^2\) for a numerical \(C_0\). Rescaling gives
(7.2). Put
\(a_j=\mathbb E[Z\mid Y=j]\), omitting a null class if necessary. Then

\[
\mathbb E\operatorname{Var}(Z\mid Y)
=\mathbb E(Z-a_Y)^2
\ge\mathbb E\min_{j=0,1}(Z-a_j)^2.
\tag{7.3}
\]

If the last expectation were smaller than
\(\delta^2\sigma^2\), Markov's inequality would put at least \(3/4\) of the
law of \(Z\) inside the union of the two intervals

\[
[a_j-2\delta\sigma,a_j+2\delta\sigma],\qquad j=0,1.
\]

Their total length is at most \(8\delta\sigma\), so (7.2) bounds their mass
by \(8C_0\delta\). Choosing
\(\delta=(16C_0)^{-1}\) gives a contradiction. Thus (7.1) holds with
\(c_0=(16C_0)^{-2}\). \(\square\)

Apply the lemma to

\[
Z=\langle X,u\rangle,\qquad Y=\mathbf1_E
\]

under any log-concave posterior. When \(M=\mu(E)=1/2\), the covariance
decomposition gives

\[
\operatorname{Var}Z
=\frac12\operatorname{Var}(Z\mid E)
+\frac12\operatorname{Var}(Z\mid E^c)
+4\langle b,u\rangle^2.
\tag{7.4}
\]

The lemma implies

\[
\frac12\operatorname{Var}(Z\mid E)
+\frac12\operatorname{Var}(Z\mid E^c)
\ge c_0\operatorname{Var}Z.
\tag{7.5}
\]

Thus the two conditional halves of a log-concave posterior cannot become
two atoms carrying essentially all of the variance. This is exactly the
step at which the example in Section 6 violates log-concavity.

### 7.2 Tests on the isotropic cross-polytope and the exponential law

Let \(K_R=\{x:\|x\|_1\le R\}\). A coordinate of the uniform law on \(K_R\)
has density

\[
f_n(x)=\frac n{2R}\left(1-\frac{|x|}{R}\right)^{n-1},
\qquad |x|\le R,
\tag{7.6}
\]

and variance

\[
\frac{2R^2}{(n+1)(n+2)}.
\]

Thus isotropy corresponds to

\[
R=\sqrt{\frac{(n+1)(n+2)}2}.
\]

For the sign split \(E=\{x_1\ge0\}\), the ratio of the within-half variance
in (7.5) to the total variance is exactly

\[
\frac{n}{2(n+1)}\longrightarrow\frac12.
\tag{7.7}
\]

For a centered one-sided exponential coordinate, split at its median. The
upper conditional variance is \(1\), the lower conditional variance is
\(1-2(\log2)^2\), and their average is

\[
1-(\log2)^2>0.51.
\tag{7.8}
\]

Thus the binary convexification lemma has the correct qualitative behavior
for both canonical stress tests.

It does **not** bound the total weak-direction variance. Indeed tilt the
coordinate density (7.6) by

\[
\exp(c_nx),\qquad c_n=\frac{n-1}{R}\asymp1.
\]

On the positive side, the linear terms in

\[
c_nx+(n-1)\log(1-x/R)
\]

cancel. Elementary Laplace estimates then give coordinate fluctuations on
the scale \(R/\sqrt n\), and hence variance of order \(n\). Similarly, for a
one-sided exponential density, the critical tilt \(c=1\) together with a
small quadratic curvature \(\rho x^2/2\) produces variance of order
\(\rho^{-1}\).

Therefore pointwise log-concavity and (7.5) rule out atomic collapse but do
not rule out a long continuous needle. For the product-exponential maximum
set, the dangerous scenario is that one of the \(n\) transverse Brownian
tilts approaches the critical value \(1\), after which the set-covariance
direction locks onto that coordinate before much quadratic curvature has
accumulated. For the cross-polytope, the analogous critical coordinate tilt
is \(c_n\asymp1\). A maximum among \(n\) transverse Brownian coordinates
reaches such a scale rapidly, so any successful proof must use the coupled
evolution of \(b_t\), not a static log-concavity bound.

### 7.3 A precise remaining convexification gate

Let

\[
Q_t=\mathbb E_t\left[
(\mathbf1_E-M_t)(X-a_t)(X-a_t)^T\right].
\tag{7.9}
\]

When \(M_t\) is exactly constant, differentiation of \(b_t\) gives

\[
db_t=Q_tC_t\,dW_t.
\tag{7.10}
\]

The angular quadratic variation of \(u_t=b_t/|b_t|\), away from \(b_t=0\),
is controlled by

\[
\frac{\|P_{u_t^\perp}Q_tC_t\|_{\mathrm{HS}}^2}{|b_t|^2}\,dt.
\tag{7.11}
\]

The atomic two-point posterior makes (7.11) vanish while keeping a huge
protected variance. The binary quantization lemma forces longitudinal
within-half spread for a log-concave posterior, but it does not by itself
force the transverse tensor in (7.11) to be large.

A sufficient genuinely new lemma would be a pathwise
**orientation-convexification inequality** saying that creation of a
log-concave survivor with variance \(R\gg1\), starting from an isotropic law,
must pay either

\[
\rho_t=\int_0^t
\left(1-\alpha_s\langle u_s,v_t\rangle^2\right)ds
\]

of universal size, or universal angular quadratic variation in (7.11)
before the direction stabilizes. One concrete form to test is

\[
\mathbb P\left\{
\operatorname{Var}_{\mu_t}\langle X,v_t\rangle\ge R,\
\rho_t\le r
\right\}
\le C\left(r+\frac1{\sqrt R}\right),
\tag{7.12}
\]

for \(t=1\) and \(0<r<1/4\). Equation (7.12), combined with
Brascamp--Lieb on \(\{\rho_t>r\}\), would imply the required inverse moment.

Neither (7.5) nor a static entropy estimate proves (7.12). The critical
tilts of the cross-polytope and product exponential show why: a log-concave
posterior may have a genuine long needle of variance \(R\), rather than two
atoms. Controlling how rapidly the adaptive maximum-set direction can select
such a needle is the surviving, specifically log-concave problem.

---

## 8. Final audit classification

The following parts are complete and dimension free:

1. the stopped density and moment identities;
2. exact mass preservation away from the zero-signal singularity;
3. soft mass preservation with the exact Doob/BDG survival bounds;
4. the formula for \(B_t\);
5. the deterministic rank-one spectral defect;
6. reduction of terminal Poincare to
   \(t^{-1}+\operatorname{Var}_{\mu_t}\langle X,v_t\rangle\);
7. the finite-state obstruction (6.28);
8. the binary log-concave convexification lemma.

What is not proved is a dimension-free survivor inverse moment for
log-concave initial laws. The atomic construction proves that this estimate
cannot come from stochastic-localization identities, isotropy, and
mass preservation alone. The remaining gap is exactly the long,
adaptively-selected log-concave needle described in Sections 7.2--7.3. This
is a failure of the proposed generic mechanism, not a KLS counterexample.
