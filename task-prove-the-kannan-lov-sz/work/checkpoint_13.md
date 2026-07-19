# Checkpoint 13: all-scale phase action and tensor normal-cone barriers

## 1. Candidate-proof status

There is still no complete dimension-free proof, and no KLS conclusion is
claimed at this checkpoint.  Four independent mechanisms have nevertheless
been reduced to exact, dimension-free statements.

1. Tensoring a spectral extremizer creates an off-diagonal stress equal to
   \(-\lambda f\otimes f\).  If both signs of every moment-preserving
   potential variation were available, this would force \(f\) to be affine.
   The convex-potential normal cone is now shown to carry order-one, rather
   than perturbative, normalized action.
2. Gaussian observation of a balanced label has an exact monotone perimeter
   functional over all signal-to-noise scales.  Its derivative is a
   nonnegative rigidity defect, but the associated unit action measure can
   migrate arbitrarily far along the logarithmic scale.
3. Bounded projected stochastic localization creates constant curvature
   outside finitely many weak directions.  Killing successively more
   martingale coefficients does not close the argument: rotation and a
   signed third-moment drift survive, and explicit smooth log-concave tests
   rule out the most natural finite-dimensional Lyapunov corrections.
4. The singular zero-strain displacement branch is exactly a cyclic binary
   martingale dilation.  The desired estimate becomes a dimension-free
   quadratic-variation bound for that dilation; local Hessian and entropy
   charges provably miss physical gaps.

All constants below are numerical and independent of the dimension.

## 2. Tensor-square spectral rigidity

Let \(\mu\) be isotropic and let \(f\) be a normalized first eigenfunction,

\[
 -L_\mu f=\lambda f,\qquad \int f\,d\mu=0,\qquad
 \int f^2\,d\mu=1.
\]

On \(\mu\otimes\mu\), the first eigenspace contains
\(F_1(x,y)=f(x)\) and \(F_2(x,y)=f(y)\).  Their mixed stress is exactly

\[
 \nabla F_1\mathbin\cdot\nabla F_2-\lambda F_1F_2
       =-\lambda f(x)f(y).                              \tag{2.1}
\]

At a two-sided interior minimizer, the full matrix first-variation equation
says that every eigenspace stress is a polynomial of total degree at most
two, modulo the barycenter and covariance constraints.  Double centering
then gives the dimension-free factorization lemma

\[
 \operatorname {dist}_{L^2(\mu\otimes\mu)}
       (f\otimes f,\mathcal P_2)\le\varepsilon
 \quad\Longrightarrow\quad
 \operatorname {dist}_{L^2(\mu)}(f,U)^2
       \le1-\sqrt{1-\varepsilon^2},                    \tag{2.2}
\]

where \(U\) is the space of linear functions.  Exact stationarity therefore
makes \(f\) affine.  On full space the eigenvalue equation then gives a
standard Gaussian factor and \(\lambda=1\); on a bounded Neumann domain an
affine first eigenfunction is impossible.

Convexity prevents promoting this argument directly.  If
\(h=f\otimes f\), \(P\in\mathcal P_2\), and

\[
                    h-P=k_+-k_-                        \tag{2.3}
\]

with convex \(k_\pm\), the new difference-of-convex audit proves

\[
 {1\over\lambda}
 \left\|D^2k_++D^2k_-\right\|_{L^2(\mu\otimes\mu;\mathrm{HS})}
 \ge \sqrt{2\left(1-\lambda^2|a|^4\right)}
 \ge \sqrt{2(1-\lambda^2)},                            \tag{2.4}
\]

where \(a=\mathbb E[Xf(X)]\).  The proof isolates the mixed Hessian block
\(\nabla f(x)\nabla f(y)^T\); its best constant approximation is
\(mm^T\), with \(m=\mathbb E\nabla f=\lambda a\).  Thus along a hypothetical
small-gap sequence the normalized convex buffer costs asymptotically
\(\sqrt2\), not \(o(1)\).

There is also a dual action barrier.  If a positive-semidefinite KKT measure
\(N\) carries the normalized residual modulo \(\mathcal P_2\), then

\[
 \int (D^2k_++D^2k_-):dN
       \ge \|h-\Pi_{\mathcal P_2}h\|_2^2
       =1-|a|^4.                                       \tag{2.5}
\]

The result holds for distributional Hessian measures and survives honest
compact truncation.  Therefore Bochner smallness of \(D^2f\) cannot make the
normalized normal-cone defect disappear.  Eliminating (2.5) would itself
require new information of the strength of the missing full-matrix
stationarity theorem.

## 3. The exact Target B reformulation

For \(a\in\mathbb R^n\), define the directional negative Sobolev norm

\[
 H_\mu(a)=\sup_g
 {\left(\int (a\cdot x)g\,d\mu\right)^2
       \over \int|\nabla g|^2\,d\mu}.                  \tag{3.1}
\]

It has the exact divergence formulation

\[
 H_\mu(a)=\min\left\{
       \int|F|^2d\mu:
       -\operatorname {div}(\mu F)=(a\cdot x)\mu
                         \right\}.                    \tag{3.2}
\]

The minimizing rows assemble the minimal Stein kernel \(\tau_*\), and

\[
                 H_\mu(a)=a^T\mathbb E(\tau_*\tau_*^T)a. \tag{3.3}
\]

If a normalized first eigenfunction satisfies
\(\operatorname {dist}(f,U)=\delta<1\), its linear projection gives

\[
                 1-\delta^2\le \lambda H_\mu(a)        \tag{3.4}
\]

for a suitable unit direction \(a\).  Consequently the operator estimate

\[
                 \mathbb E(\tau_*\tau_*^T)\preceq CI   \tag{3.5}
\]

would close the spectral stability target immediately.  What is presently
proved is only the trace estimate

\[
                 \operatorname {tr}\mathbb E
                    (\tau_*\tau_*^T)\le Cn.            \tag{3.6}
\]

Moment-map estimates for \(a^T\tau a\) do not bound
\(|\tau a|^2=a^T\tau^2a\); a positive rank-one matrix already separates the
two quantities.  Thus (3.5) is not a consequence of the known trace theorem
and, because \(H_\mu(a)\le C_P(\mu)\), a growing genuine log-concave example
would itself be KLS-relevant.  The live replacement is to combine the
near-linearity in (2.2) with the eigenfunction Bochner identity rather than
assume (3.5).

## 4. Fixed-noise and all-scale Gaussian observation

Let \(E\) be balanced, put \(U=2\mathbf1_E(X)-1\), and observe

\[
                  Y_t=\sqrt t\,X+G,\qquad G\sim N(0,I).
\]

Write \(p_t(y)=\mathbb P(E\mid Y_t=y)\),
\(F_t=\Phi^{-1}(p_t)\), and
\(I(s)=\varphi(\Phi^{-1}s)\).  Every posterior law is
\(t\)-strongly log-concave in the original variable.  The sharp posterior
covariance estimate gives

\[
                         |\nabla F_t|\le1.              \tag{4.1}
\]

At unit noise the Bayes error
\(\operatorname {err}(E)=\mathbb E\min(p_1,1-p_1)\) obeys

\[
 \mu^+(E)\ge\sqrt{2/\pi}\,\operatorname {err}(E),
 \qquad
 \operatorname {err}(E)\ge {1\over4(C_P(\mu)+2)}.      \tag{4.2}
\]

Thus a universal fixed-noise transition bound is quantitatively equivalent
to KLS, not a weaker smoothing statement.

The multiscale calculation yields the monotone perimeter functional

\[
 \mathscr P_E(t)=\sqrt t\,\mathbb E I(p_t),
 \qquad
 \mathscr P_E'(t)={1\over2\sqrt t}\,
     \mathbb E\left[I(p_t)(1-|\nabla F_t|^2)\right]\ge0, \tag{4.3}
\]

with endpoints

\[
                 \mathscr P_E(0)=0,\qquad
                 \mathscr P_E(\infty)=\mu^+(E).        \tag{4.4}
\]

For the posterior Bernoulli variance
\(B(t)=\mathbb E[p_t(1-p_t)]\), one has

\[
 -B'(t)=t^{-1}\mathbb E
          [I(p_t)^2|\nabla F_t|^2].                    \tag{4.5}
\]

Hence

\[
 d\omega_E(t)=4t^{-1}\mathbb E
          [I(p_t)^2|\nabla F_t|^2]dt                   \tag{4.6}
\]

is a probability measure on logarithmic signal-to-noise scale, and

\[
                  B(t_0)={1\over4}\omega_E([t_0,\infty)). \tag{4.7}
\]

This is an exact unit-action decomposition of Boolean uncertainty.  It
does not locate the action.  A uniform-interval model translates the whole
profile by \(-2\log L\) under dilation; whitening deletes that elementary
delay, but no proved identity pins a fixed positive fraction of
\(\omega_E\) to a universal scale for arbitrary isotropic log-concave
measures.  Such pinning would already imply (4.2).

The nuisance-cancelled low-signal expansion was also completed.  With
\(Z=2\mathbf1_E-1\), let \(N_t\) be the signed output density relative to
the standard Gaussian and \(L_t\) the unconditional output density.  The
correct chi-square information is

\[
 D(t)=\int {N_t(y)^2\over L_t(y)}\,d\gamma(y)
     =\mathbb E\left[\mathbb E(Z\mid Y_t)^2\right].     \tag{4.8}
\]

Writing \(a=\mathbb E[ZX]\), \(B=\mathbb E[ZXX^T]\),
\(R=\mathbb E[X^{\otimes3}]\), and
\(K_3=\mathbb E[ZX^{\otimes3}]-3\operatorname {Sym}(a\otimes I)\), the
exact expansion begins

\[
\begin{aligned}
 D(t)={}&t|a|^2+t^2\left({1\over2}\|B\|_{\mathrm{HS}}^2-|a|^2\right)\\
 &+t^3\left(|a|^2-\|B\|_{\mathrm{HS}}^2
       +{1\over6}\|K_3\|_{\mathrm{HS}}^2
       -\langle\operatorname {Sym}(a\otimes B),R\rangle\right)
       +O(t^4).
\end{aligned}                                           \tag{4.9}
\]

Formula (4.8) cancels independent asymmetric nuisance factors exactly.
Separate conditional tensor estimates do not: an isotropic regular simplex
has third-cumulant norm of order \(\sqrt n\).  More decisively, for a
centered one-dimensional exponential law with its median label,
\(|\mathbb E[ZX^k]|\) grows factorially.  The corresponding Hermite term
has squared norm of order \(k!\), and the \(L^2(\gamma)\) chaos series has
zero radius of convergence for every \(t>0\).  Thus no fixed-degree analytic
estimate can reach a numerical signal-to-noise ratio; the denominator in
(4.8) would need a nonperturbative resummation.

## 5. Displacement midpoint: corrected covariance and determinant ledger

For the two normalized halves \(\mu_+=2\mathbf1_E\mu\) and
\(\mu_-=2\mathbf1_{E^c}\mu\), let \(T\) be the Brenier map and
\(D=T(X)-X\).  If \(K=\mathbb EDD^T\), then the midpoint law satisfies

\[
                         \nu_{1/2}\le2\mu              \tag{5.1}
\]

and, after centering the half-measure displacement correctly,

\[
                 \operatorname {Cov}(\nu_{1/2})=I-{1\over4}K. \tag{5.2}
\]

The previously tempting all-time affine covariance formula is false; the
correct identity is

\[
 \operatorname {Cov}(Z_t)=(1-t)C_++tC_- -t(1-t)K_c.    \tag{5.3}
\]

Using only domination, log-concavity, and one-dimensional moment bounds,
the midpoint covariance has the universal spectral floor

\[
                  {1\over4800}I\preceq
                  \operatorname {Cov}(\nu_{1/2})\preceq I. \tag{5.4}
\]

Consequently

\[
 {1\over4}\operatorname {tr}K
 \le-\log\det\operatorname {Cov}(\nu_{1/2})
 \le1200\operatorname {tr}K.                           \tag{5.5}
\]

This is a useful absolute equivalence, but a determinant lower bound for
every balanced transport remains KLS-equivalent.  Moreover, covariance or
log-determinant concavity along the geodesic is false.  An explicit
one-dimensional union-of-intervals example has zero entropy, AM--GM, and
strain deficits for every \(t\), while its midpoint variance is
\((9\sqrt5-5)/16<1\).  Translation gaps, not local strain, carry the cost.

## 6. Projected localization: what extra constraints do and do not buy

For a protected balanced label, let

\[
 b=\operatorname {Cov}(\mathbf1_E,X),\qquad
 u=b/|b|,\qquad
 q=\operatorname {Cov}((u\cdot Z)^2,X).
\]

Killing both \(b\) and \(q\) in the driving matrix removes the direct
martingale term in the protected directional variance, but the direction
itself rotates.  The exact remaining stochastic coefficient contains

\[
                 {2\over|b|},C D P A u,               \tag{6.1}
\]

where \(D=\mathbb E[(\mathbf1_E-p)ZZ^T]\) and
\(P=I-uu^T\).  An anisotropic Gaussian with a balanced sign partition has
\(q=0\) while (6.1) is nonzero.  Thus a second scalar constraint does not
freeze the weak line.

A rank-two kernel does give a deterministic curvature statement: after
time \(t\), the accumulated quadratic tilt is at least \(t/3\) outside a
subspace of dimension at most four.  A static ball posterior nevertheless
shows why this alone is insufficient.  With curvature one outside a fixed
two-plane and a nearly balanced sign label in that plane, the surviving
variance and Poincare constant can both be of order \(\sqrt n\).

The exact uniform-ball calculation also rules out a misleading dynamic
counterexample.  If

\[
 d\nu\propto e^{\langle c,x\rangle-\langle Bx,x\rangle/2}
       \mathbf1_{|x|\le\sqrt{n+2}}dx,\qquad 0\preceq B\preceq TI,
\]

and \(v\) is an eigenvector of \(B\), then

\[
                  \operatorname {Var}_\nu(v\cdot X)
                         \le {6\over1-T}.               \tag{6.2}
\]

It follows from a log-concave transverse section and a one-dimensional
divergence estimate.  Therefore the apparent argmax-sign obstruction is
critical only at \(T=1\); it cannot defeat a fixed subcritical projected
flow.

Killing the full martingale coefficient of the protected variance leaves a
signed third-moment drift.  Even at \(A=I\), it need not have a favorable
sign.  A smooth isotropic log-concave product of Gaussian-smoothed centered
exponentials, with a diluted balanced label, makes this drift strictly
positive while the drifts and quadratic variations of the natural
corrections \(|b|^2\) and \(b^TA^{-1}b\) are \(O(\varepsilon^2)\).
Consequently every bounded-derivative Lyapunov correction
\(r+\phi(|b|^2,b^TA^{-1}b)\), including fixed linear combinations, fails
uniformly as \(\varepsilon\downarrow0\).

The function-specific version has now been closed to an exact ledger.  Fix a
normalized first eigenfunction \(f\), let
\(\mathcal E_t=\int|\nabla f|^2d\mu_t\), and choose \(C_tb_t=0\) for
\(b_t=\operatorname {Cov}_{\mu_t}(f,X)\).  Then

\[
 \mathbb E\operatorname {Var}_{\mu_T}f=1,\qquad
 \mathbb E\mathcal E_T=\lambda.                         \tag{6.3}
\]

For every relaxed hard driver with
\(\operatorname {tr}(I-C_t^2)\le1\), the accumulated tilt has curvature
\(T/2\) outside one spectral direction \(v_T\).  If
\(R_T=\operatorname {Var}_{\mu_T}(v_T\cdot X)\), the audited rank-one
defect theorem gives

\[
             1\le {192\over T}\lambda
                    +96\,\mathbb E[R_T\mathcal E_T].    \tag{6.4}
\]

Under the energy-biased path law
\(d\mathbb Q=(\mathcal E_T/\lambda)d\mathbb P\), a small gap therefore
forces \(\mathbb E_{\mathbb Q}R_T\ge c/\lambda\).  Conversely, the soft
driver \(C=I-\alpha P_b\), with
\(\delta=(1-\alpha)^2\), has full curvature
\(B_T\succeq\delta TI\) but obeys the exact mean/variance tradeoff

\[
             1-\mathbb E m_T^2\le{\lambda\over\delta T}. \tag{6.5}
\]

Thus a universal bound on the energy-biased adaptive survivor in (6.4)
would prove KLS, while enough soft curvature to bypass it necessarily loses
the protected mean when \(\lambda\) is small.  The calculation does not
provide a weaker bootstrap.

## 7. The zero-strain branch as a cyclic binary dilation

For a uniform isotropic convex body \(K\), the exact zero-gap midpoint data
consist of a half-volume set \(R\) and a locally constant cyclic firmly
nonexpansive map \(b:R\to\mathbb R^n\) such that

\[
                  K=(I-b)(R)\mathbin{\dot\cup}(I+b)(R). \tag{7.1}
\]

If \(Z\) is uniform on \(R\) and \(\varepsilon\) is a fair sign, then

\[
                         U=Z+\varepsilon b(Z)           \tag{7.2}
\]

is uniform on \(K\) and is a martingale dilation.  Iterating (7.2) until
the chain first exits \(R\) gives

\[
 \mathbb P(\tau>m)=2^{-m},\qquad
 \mathbb E|Z_\tau-Z_0|^2=2\mathbb E_R|b|^2.            \tag{7.3}
\]

Writing \(S=\mathbb E_R[bb^T]\), isotropy gives the exact covariance split

\[
 \operatorname {Cov}(R)=I-S,\qquad
 \operatorname {Cov}(K\setminus R)=I+S.               \tag{7.4}
\]

The desired absolute theorem is therefore

\[
                         \operatorname {tr}S\le C.      \tag{7.5}
\]

Local singular-Hessian mass cannot prove (7.5).  The realizable
one-dimensional example

\[
 K=[0,6],\qquad E=[0,1]\cup[2,4]                       \tag{7.6}
\]

uses translations \(+1\) and \(+2\), has zero source-wall functional, and
after isotropic normalization satisfies

\[
 \mathbb E_R|b|^2={1\over4},\qquad
 \operatorname {Var}(R)={3\over4}.                    \tag{7.7}
\]

The canonical firm extension sends part of the complementary set outside
\(K\), invalidating the natural ambient mixed-volume containment.  A valid
proof of (7.5) must charge the physical gaps between plateau components and
use cyclic realizability, not merely pairwise firm nonexpansiveness.

## 8. Updated mechanism registry

The following routes are now blocked in their tested forms.

* A small difference-of-convex buffer derived only from the Bochner bound;
  its normalized curvature is bounded below by (2.4).
* Replacing the directional Stein-kernel operator estimate by the known
  trace estimate or by scalar moment-map control.
* Fixed-noise Gaussian smoothing without a new transition-mass theorem, and
  multiscale action without a scale-location principle.  Fixed-degree
  low-signal chaos is also blocked by factorial mixed moments.
* Two protected scalar moments, or the full protected-variance martingale
  coefficient, without a new control of rotation and third moments.
* Static finite-defect curvature without control of the covariance in the
  adaptively selected defect space.  The function-specific ledger (6.4)
  shows that its energy-biased form grows like \(1/\lambda\) along every
  hypothetical small-gap sequence.
* Entropy, Cayley strain, source-wall Hessian mass, or a canonical firm
  extension as the sole charge in the zero-gap displacement branch.

The live incompatible routes are:

1. a function-specific projected localization which preserves the mean of
   a first eigenfunction and estimates only its terminal variance;
2. a Bochner/PDE rigidity theorem for a first eigenfunction already close
   to a linear function;
3. a nonperturbative Gaussian-information estimate using the exact
   denominator in (4.8); and
4. a global physical-gap or laminar-extrusion theorem for cyclic binary
   dilations.

No live route is counted as a proof until it supplies a dimension-free
estimate for every log-concave measure and every admissible test function.

## 9. Audit state

The rank-one defect theorem from Checkpoint 12 remains clean-room verified.
The identities in Sections 2--7 have been checked on the Gaussian, cube,
isotropic simplex, isotropic crosspolytope, products of centered exponential
laws, a radial exponential law, and the asymmetric interval construction
where applicable.  These checks expose no hidden dimension factor, but they
also do not bridge any of the boxed missing global statements.

The next audit gate is deliberately binary.  A candidate must either:

* turn the one-direction trace defect in function-specific localization
  into an averaged estimate with a numerical constant; or
* turn tensor near-linearity plus the Bochner identity into a numerical
  lower bound for \(\lambda\), without invoking the unresolved directional
  \(H^{-1}\) estimate.

Only after one of those gates closes will the approximation, affine-support,
all-test-function, and canonical-instantiation audits be restarted for a
candidate full proof.
