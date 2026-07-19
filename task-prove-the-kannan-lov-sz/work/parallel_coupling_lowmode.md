# Parallel coupling at a near-linear low mode

This note audits the directional content of
`/tmp/parallel_src/thin-shell-arXiv-v2.tex`.

## 0. Verdict

Let \(\mu\) be isotropic and log-concave, let \(f\) be a normalized first
eigenfunction,

\[
 Af=\lambda f,\qquad \mathbb Ef=0,\qquad \mathbb Ef^2=1,
\]

and write

\[
 a=\mathbb E[Xf],\qquad f=a\cdot x+r,\qquad
 \delta=\|r\|_2,\qquad |a|^2=1-\delta^2.
 \tag{0.1}
\]

The parallel-coupling theorem gives, for every unit vector \(v\),

\[
 H_\mu(v)\le \mathbb E|M_1v|^2,
 \tag{0.2}
\]

where \(M_t\) solves \(M'_t=A_tM_t\) along stochastic localization and
\(A_t=\operatorname {Cov}(\mu_t)\).  The paper proves the trace estimate

\[
 \mathbb E\operatorname {Tr}(M_1^TM_1)\le Cn.
 \tag{0.3}
\]

The trace proof does not become directional merely because \(v\) is the
linear projection of a near-linear eigenfunction.  In fact the exact
low-mode evolution goes in the opposite direction.  Put

\[
 c_t=\operatorname {Cov}_{\mu_t}(X,f).
\]

Then

\[
 \boxed{\quad Z_t=M_t^Tc_t\text{ is a martingale},\qquad
 \mathbb EZ_t=a.\quad}
 \tag{0.4}
\]

At time one, the posterior is \(1\)-strongly log-concave, and hence

\[
 \mathbb E|c_1|^2\le\lambda.
 \tag{0.5}
\]

If \(b=a/|a|\), (0.4)--(0.5) give the fully explicit inequality

\[
 \boxed{
 \mathbb E|M_1b|^2\ge {1-\delta^2\over\lambda}.}
 \tag{0.6}
\]

Thus a hypothetical near-linear small-gap mode necessarily selects an
expensive direction of the parallel coupling.  Guan's eigenvalue-rank
estimate says that there cannot be many such directions; it does not rule
out one direction chosen adversarially by \(f\).

The weakest natural final-time alignment estimate that closes this argument
is not the full bound \(\mathbb E|M_1b|^2\le C\).  It is only

\[
 \boxed{
 \mathbb E\left[
 {\langle M_1b,c_1\rangle^2\over |c_1|^2}
 \mathbf1_{\{c_1\ne0\}}
 \right]\le C.}
 \tag{FA}
\]

Indeed, (FA) and (0.4)--(0.5) imply

\[
 1-\delta^2\le C\lambda.
 \tag{0.7}
\]

(FA) is strictly weaker than directional control of \(M_1b\): it sees only
the random component used by the martingale pairing.  It holds in all the
requested canonical tests--Gaussian, cube, product exponential, regular
simplex, and isotropic ball--because those models have either product
factorization or enough symmetry to upgrade the trace theorem.  It is not a
consequence of Guan's rank estimates: an explicit one-spike matrix path
below obeys their rank/trace scale while the left side of (FA) grows like
\((\log n)^{16}\).

The exact Bochner defect remains useful but does not supply the missing
alignment.  After localization it controls posterior Hessian energy and
posterior variation of \(\nabla f\), whereas (FA) is weighted by the
path-ordered amplification \(M_1\).  No estimate in the parallel-coupling
proof controls that correlation.  The formal low-mode lemma is therefore
(FA), or the stronger high-covariance occupation estimate in Section 4;
neither follows from the trace theorem presently proved.

## 1. What the parallel-coupling theorem actually supplies

For compactly supported \(\mu\), define

\[
 d\mu_{t,\theta}(x)
 ={\exp(\theta\cdot x-t|x|^2/2)\over
   \int\exp(\theta\cdot y-t|y|^2/2)d\mu(y)}d\mu(x).
 \tag{1.1}
\]

Along the Brownian flow of the paper, write

\[
 \mu_t=\mu_{t,\theta_t},\qquad
 a_t=\mathbb E_tX,\qquad
 A_t=\operatorname {Cov}_{\mu_t}(X).
 \tag{1.2}
\]

The differential of the flow with respect to its initial exponential tilt
is the product integral

\[
 M_0=I,\qquad {d\over dt}M_t=A_tM_t.
 \tag{1.3}
\]

For every \(v\), the parallel-coupling corollary and the infinitesimal
Wasserstein argument in the source give, at a general time \(T>0\),

\[
 H_\mu(v)\le {1\over T^2}\mathbb E|M_Tv|^2.
 \tag{1.4}
\]

We use \(T=1\).  The covariance process satisfies

\[
 dA_t=\sum_{k=1}^n H_{k,t}\,dB_{k,t}-A_t^2dt,
 \qquad A_t\preceq t^{-1}I,
 \tag{1.5}
\]

where the \(H_{k,t}\) are centered third-moment matrices.

Let

\[
 \lambda_1(t)\ge\cdots\ge\lambda_n(t)>0
\]

be the eigenvalues of \(A_t\), and set

\[
 \tau_k=\inf\{t>0:\lambda_k(t)\ge3\}.
\]

The Guan-type rank theorem proved in the source gives

\[
 \mathbb P\{\tau_k\le t\}
 \le C{n\over k}e^{-t^{-1/8}},
 \qquad
 \mathbb E\tau_k^{-2}
 \le C\left(1+\log{n\over k}\right)^{16}.
 \tag{1.6}
\]

The deterministic product-integral inequality is

\[
 \|M_1\|_{\mathrm{HS}}^2
 \le\sum_{k=1}^n
 \exp\left(2\int_0^1\lambda_k(t)dt\right).
 \tag{1.7}
\]

Combining (1.6)--(1.7) over all ranks yields (0.3).  For one fixed vector,
however, the argument can use only the top rank without further alignment:

\[
\begin{aligned}
 |M_1v|^2
 &\le\|M_1\|_{\mathrm{op}}^2
 \le\exp\left(2\int_0^1\lambda_1(t)dt\right),\\
 \mathbb E|M_1v|^2
 &\le C(1+\log n)^{16}.
\end{aligned}
 \tag{1.8}
\]

This polylogarithmic directional consequence is genuine, but it is not
dimension free.  The averaging over \(k\) which turns (1.8) into (0.3)
has no analogue for a single prescribed \(v\).

## 2. Exact low-mode martingale

### 2.1 Evolution of the localized covariance with \(f\)

For any fixed integrable function \(h\), stochastic localization gives

\[
 d\mathbb E_th
 =\operatorname {Cov}_{\mu_t}(h,X)\cdot dB_t.
 \tag{2.1}
\]

Put

\[
 m_t=\mathbb E_tf,
 \qquad
 c_t=\mathbb E_t[(X-a_t)(f-m_t)],
 \tag{2.2}
\]

and define the symmetric matrix

\[
 D_t=\mathbb E_t[(f-m_t)(X-a_t)(X-a_t)^T].
 \tag{2.3}
\]

Equation (2.1) gives

\[
 da_t=A_t,dB_t,
 \qquad dm_t=c_t\cdot dB_t.
 \tag{2.4}
\]

To differentiate \(c_t\), first write

\[
 c_t=\mathbb E_t[Xf]-a_tm_t.
\]

The martingale coefficient of \(\mathbb E_t[Xf]\) is

\[
 \operatorname {Cov}_t(Xf,X)
 =D_t+m_tA_t+a_tc_t^T.
 \tag{2.5}
\]

On the other hand, Itô's product rule gives

\[
 d(a_tm_t)
 =m_tA_t,dB_t+a_t(c_t\cdot dB_t)+A_tc_tdt.
 \tag{2.6}
\]

Subtracting (2.6) from (2.5) yields the exact vector SDE

\[
 \boxed{dc_t=D_t,dB_t-A_tc_tdt.}
 \tag{2.7}
\]

Since \(dM_t^T=M_t^TA_tdt\), there is no quadratic-covariation term and

\[
\begin{aligned}
 d(M_t^Tc_t)
 &=M_t^TA_tc_tdt
   +M_t^T(D_t,dB_t-A_tc_tdt)\\
 &=M_t^TD_t,dB_t.
\end{aligned}
 \tag{2.8}
\]

Thus \(Z_t=M_t^Tc_t\) is a martingale.  At time zero,

\[
 c_0=\operatorname {Cov}_\mu(X,f)=a,
 \qquad M_0=I,
\]

and hence

\[
 \boxed{a=\mathbb E[M_1^Tc_1].}
 \tag{2.9}
\]

Compact support makes all martingales above integrable.  General
log-concave measures follow by the same truncation and lower-semicontinuity
scheme used in the source.

### 2.2 The endpoint energy of \(c_1\)

The posterior \(\mu_1\) is \(1\)-strongly log-concave on a convex support.
Consequently

\[
 A_1\preceq I,
 \qquad C_P(\mu_1)\le1.
 \tag{2.10}
\]

For every unit vector \(u\),

\[
 |u\cdot c_1|^2
 \le\operatorname {Var}_1(u\cdot X)\operatorname {Var}_1(f)
 \le\mathbb E_1|\nabla f|^2.
\]

Taking the supremum over \(u\) gives the pointwise posterior estimate

\[
 |c_1|^2\le\mathbb E_1|\nabla f|^2.
 \tag{2.11}
\]

The localization densities average back to \(\mu\).  Therefore

\[
 \boxed{
 \mathbb E_B|c_1|^2
 \le\mathbb E_B\mathbb E_1|\nabla f|^2
 =\mathbb E_\mu|\nabla f|^2=\lambda.}
 \tag{2.12}
\]

### 2.3 The Cauchy--Schwarz step, with every normalization visible

Assume \(a\ne0\), and set

\[
 b={a\over|a|}.
\]

Take the scalar product of (2.9) with \(b\).  Since \(M_1b\) and \(c_1\)
live on the same Brownian probability space,

\[
\begin{aligned}
 |a|
 &=b\cdot a\\
 &=\mathbb E\,[b\cdot M_1^Tc_1]\\
 &=\mathbb E\,\langle M_1b,c_1\rangle.
\end{aligned}
 \tag{2.13}
\]

Cauchy--Schwarz first in \(\mathbb R^n\), then in probability, gives

\[
\begin{aligned}
 |a|
 &\le\mathbb E[|M_1b|\,|c_1|]\\
 &\le
 \big(\mathbb E|M_1b|^2\big)^{1/2}
 \big(\mathbb E|c_1|^2\big)^{1/2}.
\end{aligned}
 \tag{2.14}
\]

Squaring and using (2.12) and \(|a|^2=1-\delta^2\) gives

\[
 1-\delta^2
 \le\lambda\,\mathbb E|M_1b|^2,
 \tag{2.15}
\]

which is (0.6).  Notice that this has the same direction as the elementary
test of \(H_\mu(b)\) with \(f\).  The coupling does not dilute the bad mode;
its derivative must amplify enough to carry the fixed pairing (2.13) while
\(c_1\) has energy at most \(\lambda\).

## 3. What near-linearity and Bochner add after localization

The affine residual gives an endpoint relation which is stronger than
arbitrary low-mode information.  Since

\[
 f=a\cdot x+r,
\]

we have at every time

\[
 c_t=A_ta+e_t,
 \qquad e_t=\operatorname {Cov}_{\mu_t}(X,r).
 \tag{3.1}
\]

For \(t>0\), \(A_t\preceq t^{-1}I\), so conditional Cauchy--Schwarz gives

\[
 |e_t|^2
 \le {1\over t}\operatorname {Var}_t(r).
\]

Averaging and using \(\mathbb E_B\mathbb E_tr^2=\delta^2\),

\[
 \boxed{\mathbb E|e_t|^2\le{\delta^2\over t}.}
 \tag{3.2}
\]

Similarly, \(C_P(\mu_t)\le t^{-1}\) gives

\[
 \boxed{\mathbb E|c_t|^2\le {\lambda\over t^2}.}
 \tag{3.3}
\]

Combining (3.1)--(3.3),

\[
 \mathbb E|A_ta|^2
 \le {2\lambda\over t^2}+{2\delta^2\over t}.
 \tag{3.4}
\]

In particular, a hypothetical regime \(\lambda,\delta\to0\) forces the
time-one posterior covariance to collapse in direction \(a\).  This does
not bound \(M_1a\): the product integral remembers covariance encountered
over the whole path, and covariance contraction at the endpoint is naturally
paired with expansion of the derivative flow.

The exact Bochner defect is

\[
 \left(B_f-\lambda\operatorname {Var}_\mu(\nabla f)\right)
 +\mathcal R_V(f)
 =\lambda^3|a|^2,
 \qquad
 B_f=\mathbb E_\mu\|D^2f\|_{\mathrm{HS}}^2.
 \tag{3.5}
\]

In particular,

\[
 \lambda^2(1-\lambda|a|^2)\le B_f\le\lambda^2.
 \tag{3.6}
\]

This information also averages through localization:

\[
 \mathbb E_B\mathbb E_t\|D^2f\|_{\mathrm{HS}}^2=B_f.
 \tag{3.7}
\]

Applying the posterior Poincare inequality componentwise gives

\[
 \mathbb E_B\operatorname {Var}_{\mu_t}(\nabla f)
 \le {B_f\over t}.
 \tag{3.8}
\]

At time one, the localized gradient is therefore almost deterministic when
\(\lambda\) is small.  Equations (3.2), (3.4), and (3.8) are the full direct
benefit of near-linearity and Bochner for this route.

They are not weighted by \(M_t\).  The missing estimates would have to
control quantities such as \(M_t^Te_t\), or the correlation of
\(M_tb\) with high-covariance eigenspaces.  Cauchy--Schwarz with (3.2) and
the trace theorem loses a factor of order \(\sqrt n\); using the top
singular-value estimate loses a polylogarithm.  The scalar defect (3.5)
does not constrain this pathwise alignment.

## 4. The precise low-mode alignment lemmas

### 4.1 Weakest final-time form

Define

\[
 Y_1=
 \begin{cases}
 \langle M_1b,c_1\rangle/|c_1|,&c_1\ne0,\\
 0,&c_1=0.
 \end{cases}
 \tag{4.1}
\]

Then (2.13) is simply

\[
 |a|=\mathbb E[|c_1|Y_1].
 \tag{4.2}
\]

Consequently the following is a formal low-mode lemma.

**Final-alignment lemma.**  If, for every near-linear normalized first
eigenfunction under consideration,

\[
 \mathbb EY_1^2\le C_A,
 \tag{4.3}
\]

then

\[
 \boxed{\lambda\ge {1-\delta^2\over C_A}.}
 \tag{4.4}
\]

Indeed, (4.2), Cauchy--Schwarz, and (2.12) give

\[
 |a|^2
 \le\mathbb E|c_1|^2\,\mathbb EY_1^2
 \le C_A\lambda.
\]

This is exactly (FA).  It is weaker than
\(\mathbb E|M_1b|^2\le C_A\), because pointwise

\[
 Y_1^2\le|M_1b|^2.
\]

It is also the weakest natural quadratic estimate furnished by the pairing
(2.13): components of \(M_1b\) orthogonal to \(c_1\) are irrelevant.

### 4.2 A stronger pathwise occupation form tied to Guan's theorem

Put

\[
 w_t=M_tb,
 \qquad u_b(t)=\mathbb E|w_t|^2.
\]

From (1.3),

\[
 u_b'(t)=2\mathbb E\langle w_t,A_tw_t\rangle.
 \tag{4.5}
\]

Let

\[
 Q_t=(A_t-3I)_+.
\]

Since \(A_t\preceq3I+Q_t\),

\[
 {d\over dt}\log u_b(t)
 \le6+2{\mathbb E\langle w_t,Q_tw_t\rangle\over u_b(t)}.
 \tag{4.6}
\]

Hence the high-covariance occupation estimate

\[
 \boxed{
 \mathcal A_\mu(b):=
 \int_0^1
 {\mathbb E\langle M_tb,(A_t-3I)_+M_tb\rangle
  \over\mathbb E|M_tb|^2},dt
 \le C_A}
 \tag{OA}
\]

implies

\[
 \mathbb E|M_1b|^2\le e^{6+2C_A},
 \tag{4.7}
\]

and therefore implies (FA).  This formulation requires only bounded
*amplification-weighted occupation* of the exceptional covariance space.
It does not ask that \(b\) be uniformly spread among all \(n\) directions.

If \(P_t=\mathbf1_{[3,\infty)}(A_t)\), then

\[
 Q_t\preceq t^{-1}P_t.
 \tag{4.8}
\]

Guan's estimate controls the unweighted rank

\[
 \mathbb E\operatorname {rank}P_t
 =\sum_i\mathbb P\{\lambda_i(t)\ge3\},
 \tag{4.9}
\]

whereas (OA) contains the adapted, already-amplified vector \(M_tb\).
The flow itself favors high covariance.  If
\(q_t=w_t/|w_t|\), then pathwise

\[
 q_t'=(A_t-\langle q_t,A_tq_t\rangle I)q_t,
 \tag{4.10}
\]

so in a commuting eigenbasis the weight on an eigenvector of eigenvalue
\(\lambda_i\) grows relative to that on \(\lambda_j\) at rate
\(\lambda_i-\lambda_j\).  There is no automatic delocalization; the
product integral selects exceptional eigenspaces.

## 5. Why rank and trace do not imply alignment

### 5.1 Final singular-space formulation

Let

\[
 S=M_1^TM_1.
\]

For every unit \(b\), layer cake gives

\[
 \mathbb E[b^TSb]
 =\int_0^\infty
 \mathbb E\|\mathbf1_{(s,\infty)}(S)b\|^2ds,
 \tag{5.1}
\]

while

\[
 \mathbb E\operatorname {Tr}S
 =\int_0^\infty
 \mathbb E\operatorname {rank}\mathbf1_{(s,\infty)}(S)ds.
 \tag{5.2}
\]

Thus the exact missing trace-to-direction statement would be an alignment
bound such as

\[
 \mathbb E\|\mathbf1_{(s,\infty)}(S)b\|^2
 \le {C\over n}
 \mathbb E\operatorname {rank}\mathbf1_{(s,\infty)}(S).
 \tag{5.3}
\]

Near-linearity identifies \(b\), but supplies no randomness of \(b\)
relative to these singular spaces.  The weaker (FA) avoids asking for all
of (5.3), yet it still needs information absent from ranks.

### 5.2 A one-spike process permitted by the rank scale

The following is not asserted to be an actual stochastic-localization
covariance process.  It is an information-theoretic countermodel showing
that the deterministic product-integral bounds, the Guan rank scale, and
the trace conclusion do not imply either (OA) or (FA).

Fix \(b=e_1\), and let

\[
 \varepsilon_n\asymp(1+\log n)^{-8}.
\]

Consider a commuting covariance path whose exceptional eigenvalue starts
at one, rises continuously to \(1/\varepsilon_n\) by time
\(\varepsilon_n\), follows \(1/t\) until a fixed time, and then decreases
below the threshold three.  It can be chosen throughout so that

\[
 A_0=I,
 \qquad0\preceq A_t\preceq t^{-1}I,
 \qquad\operatorname {rank}\mathbf1_{[3,\infty)}(A_t)=1
\]

whenever the exceptional eigenvalue is above three.  One explicit rising
segment is

\[
 \alpha(t)=
 {1\over1-(1-\varepsilon_n)t/\varepsilon_n},
 \qquad0\le t\le\varepsilon_n,
 \tag{5.4}
\]

followed by \(\alpha(t)=1/t\).  All other eigenvalues may be
\((1+t)^{-1}\).

The exceptional product-integral factor obeys

\[
 |M_1b|^2\asymp\varepsilon_n^{-2}
 \asymp(1+\log n)^{16},
 \tag{5.5}
\]

up to a universal factor if the path is made to decay after a fixed time.
At the same time,

\[
 \operatorname {Tr}(M_1^TM_1)
 \le4(n-1)+C(1+\log n)^{16}=O(n).
 \tag{5.6}
\]

The first high-eigenvalue hitting time has inverse-square scale
\((1+\log n)^{16}\), exactly the \(k=1\) allowance in (1.6), while no
second eigenvalue is exceptional.  Thus the rank theorem is designed to
permit this single spike.

The occupation in (OA) is

\[
 \mathcal A_\mu(b)\asymp\log(1/\varepsilon_n)
 \asymp\log\log n,
 \tag{5.7}
\]

and is unbounded.  At the endpoint, choose algebraic vectors
\(c_1=M_1^{-T}a\) with \(a\parallel b\).  Then

\[
 M_1^Tc_1=a,
 \qquad
 {\langle M_1b,c_1\rangle^2\over|c_1|^2}
 =|M_1b|^2,
 \tag{5.8}
\]

so (FA) fails by the same factor.  By letting the exceptional covariance
decay to \(|M_1b|^{-1}\) at time one, one may also arrange
\(c_1=A_1a\), matching the endpoint affine relation (3.1) with zero error.

Finally, the scalar Bochner data can be overlaid without contradiction:
take \(\lambda\asymp|M_1b|^{-2}\), choose \(|a|^2=1-\delta^2\), and set
\(B_f=\lambda^2\).  Then

\[
 B_f-\lambda\operatorname {Var}(\nabla f)
 =\lambda^3|a|^2,
\]

with zero curvature term.  This is only a formal package, not a
log-concave counterexample, but it proves that the scalar Bochner identity
cannot exclude the one-spike scenario without an additional coupling to
the evolving eigenspaces.

## 6. Tests of the alignment lemma

None of the standard models disproves (FA) or (OA).

### 6.1 Gaussian

For \(\gamma_n\),

\[
 A_t=(1+t)^{-1}I,
 \qquad M_t=(1+t)I.
\]

For \(f=b\cdot x\), \(c_1=b/2\), and hence

\[
 {\langle M_1b,c_1\rangle^2\over|c_1|^2}=4.
\]

Moreover \((A_t-3I)_+=0\), so (OA) holds with zero exceptional
occupation.

### 6.2 Isotropic cube

For \([ -\sqrt3,\sqrt3]^n\), localization preserves the product structure,
so \(A_t\) and \(M_t\) are diagonal in the coordinate basis.  Every
one-dimensional posterior is supported on an interval of length
\(2\sqrt3\).  Popoviciu's inequality gives

\[
 \operatorname {Var}_t(X_i)\le3.
\]

Thus \((A_t-3I)_+=0\) pathwise and (OA) is immediate.  Consequently (FA)
holds for every first-mode combination, including the near-linear
interval modes tested previously.

### 6.3 Product exponential laws

For a product of centered isotropic one-dimensional exponential or Laplace
laws, localization again factorizes.  Write

\[
 M_t=\operatorname {diag}(m_1(t),\ldots,m_n(t)).
\]

The dimension-one instance of the parallel-coupling theorem gives a
universal bound

\[
 \mathbb Em_i(1)^2\le C
\]

for every factor.  Hence, for any deterministic unit \(b\),

\[
 \mathbb E|M_1b|^2
 =\sum_i b_i^2\mathbb Em_i(1)^2\le C.
 \tag{6.1}
\]

This proves (FA).  It also proves (OA), since its integrand is bounded by
the full directional logarithmic growth and

\[
 {1\over2}\log\mathbb E|M_1b|^2
 =\int_0^1
 {\mathbb E\langle M_tb,A_tM_tb\rangle
  \over\mathbb E|M_tb|^2}dt.
 \tag{6.2}
\]

### 6.4 Isotropic Euclidean ball

For the isotropic ball, the law of the parallel flow is orthogonally
equivariant.  Therefore

\[
 K:=\mathbb E[M_1^TM_1]
\]

commutes with every orthogonal transformation and hence \(K=\kappa I\).
The trace theorem gives \(\kappa\le C\), so

\[
 \mathbb E|M_1b|^2\le C
\]

for every \(b\).  In particular, (FA) holds for the degree-one Neumann
modes whose affine error is \(\delta\sim n^{-1}\).  Thus the ball
counterexample to the multiplicative Hessian interpolation is not a
counterexample to parallel-coupling alignment.

### 6.5 Regular simplex

The symmetry group of a regular simplex acts irreducibly on its centered
\(n\)-dimensional representation.  Equivariance again implies that
\(K=\mathbb E[M_1^TM_1]\) is a scalar on this representation.  From the
trace theorem,

\[
 K\preceq CI.
\]

Hence (FA) and (OA) hold for every direction.  This conclusion does not
require an explicit formula for the simplex's first Neumann eigenfunction.

### 6.6 General symmetry lemma

The last two tests are instances of a useful formal result.  Suppose a
compact orthogonal group \(G\) preserves \(\mu\), and the \(G\)-orbit of
\(b\) spans an irreducible subspace \(E\) of dimension \(d\).  Equivariance
and irreducibility give

\[
 \mathbb E|M_1b|^2
 ={1\over d}\operatorname {Tr}
 \left(P_E\mathbb E[M_1^TM_1]P_E\right)
 \le {Cn\over d}.
 \tag{6.3}
\]

Thus the parallel trace theorem is dimension-free on any low mode whose
symmetry orbit has dimension comparable to \(n\).  A single near-linear
mode for a general anisotropic isotropic measure can have \(d=1\), in which
case (6.3) reduces to the useless trace bound.

## 7. A separate high-rank low-mode lemma

There is one other rigorous way in which rank helps.  Suppose the first
eigenspace contains orthonormal functions \(f_1,\ldots,f_m\) whose linear
projections have singular values at least \(\sqrt{1-\delta^2}\).  Choose
corresponding orthonormal directions \(b_1,\ldots,b_m\).  Testing
\(H_\mu(b_j)\) with the associated singular vectors and summing gives

\[
 m(1-\delta^2)
 \le\lambda\sum_{j=1}^mH_\mu(b_j)
 \le\lambda\sum_{i=1}^nH_\mu(e_i)
 \le C\lambda n.
 \tag{7.1}
\]

Therefore

\[
 \boxed{\lambda\ge {m\over n}{1-\delta^2\over C}.}
 \tag{7.2}
\]

This closes the low-mode problem when \(m\ge cn\).  The exact Bochner
defect says that the components of the centered gradient occupy a bottom
spectral cluster, but it does not give a lower bound on the effective rank
of that cluster.  Its covariance may be rank one.  Hence (7.2) does not
upgrade a single near-linear eigenfunction automatically.

## 8. Precise remaining target

The full parallel-coupling theorem has now been used down to its
directional endpoint.  The source proves ranks and trace; the low-mode
calculation produces the exact martingale pairing

\[
 |a|=\mathbb E\langle M_1b,c_1\rangle,
 \qquad\mathbb E|c_1|^2\le\lambda.
\]

The minimal new statement needed for Target B is (FA): a dimension-free
bound on the component of the amplified direction along the localized
low-mode covariance.  A more geometric sufficient statement is (OA): a
dimension-free bound on amplification-weighted occupation of covariance
eigenvalues above three.

Both statements survive Gaussian, cube, product exponential, simplex, and
ball tests.  Neither follows from Guan's eigenvalue-rank estimate, and the
one-spike process shows exactly why: rank one is cheap in the trace budget
and the product integral actively aligns with it.  The Bochner defect
controls unweighted posterior derivatives but does not couple them to this
adapted spike.  Therefore the parallel-coupling route does not yet close
Target B; it reduces it to the explicit final-alignment estimate (FA).
