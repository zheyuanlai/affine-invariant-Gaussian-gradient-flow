# Tensor extremizers on slow-growth dimensions: exact plateaus and a buffer no-go

## 0. Verdict

For \(n\geq1\), let

\[
 \mathcal C_n=
 \sup\{C_P(\mu):\mu\text{ isotropic and log-concave on }\mathbb R^n\}.
 \tag{0.1}
\]

Assume hypothetically that \(\mathcal C_n\to\infty\), and use the known
bound

\[
 1\leq\mathcal C_n\leq K\log(en)                         \tag{0.2}
\]

with a numerical \(K\).  The conclusions of this audit are as follows.

1. There are dimensions \(n_j\to\infty\) such that

   \[
      {\mathcal C_{2n_j}\over\mathcal C_{n_j}}\longrightarrow1.
      \tag{0.3}
   \]

   In fact there are integers \(r_j\to\infty\) for which

   \[
      {\mathcal C_{r_jn_j}\over\mathcal C_{n_j}}\longrightarrow1.
      \tag{0.4}
   \]

   Thus an \(r_j\)-fold tensor power of an almost worst \(n_j\)-dimensional
   law is a relative near-minimizer of the spectral gap in dimension
   \(r_jn_j\), while its first-eigenvalue multiplicity tends to infinity.

2. This growing multiplicity does not amplify the mixed eigenvalue splitting
   relative to the convex buffer.  If \(f\) is a normalized first
   eigenfunction with gap \(\lambda\), \(a=\mathbb E[Xf]\),
   \(\beta=|a|^2\), and \(A\) is any symmetric zero-diagonal \(r\times r\)
   matrix, the moment-tangent mixed perturbation has first-cluster branch
   matrix

   \[
       \lambda(1-\beta^2)A.                             \tag{0.5}
   \]

   Its Hessian, and therefore every difference-of-convex buffer realizing
   it, has the exact lower bound

   \[
   \|D^2Q_A\|_2
   \ge
   \sqrt{\lambda^2+(1-2\lambda^2)\beta^2}\,\|A\|_{\rm F}.
                                                               \tag{0.6}
   \]

   Consequently

   \[
   {\text{mixed spectral splitting}\over
    \text{least possible }L^2\text{ convex-buffer curvature}}
   \le {\|A\|_{\rm op}\over\|A\|_{\rm F}}\le1.           \tag{0.7}
   \]

   This is independent of \(r\).  Complete graphs attain only a constant
   ratio; disjoint-pair or random-sign constructions make the ratio worse.

3. Strong-convexity regularization plus relative near-minimality yields
   two-sided approximate stationarity only under the additional scale
   condition

   \[
                    \Delta=o(\rho),                    \tag{0.8}
   \]

   where \(\Delta\) is the excess in the \(\log\lambda_1\) objective and
   \(\rho\) is the Hessian margin of the regularized potential in the
   perturbation norm.  Mere convergence of strongly convex approximations
   gives no such rate.  A first-order regularization loss
   \(\Delta\asymp\rho\) is consistent with all available information and
   with the buffer calculation (0.5)--(0.7).  Moreover, an Ekeland move in
   the full potential space generically splits the tensor eigenspace, while
   an Ekeland move restricted to product potentials does not test mixed
   directions.

4. Guan's ordered-eigenvalue stopping theorem supplies no missing selector.
   Tensoring \(r\) copies creates \(r\) block low-mode directions in ambient
   dimension \(rn\).  A rank-\(r\) exceptional localization subspace may
   contain all of them, and Guan's bound then reads

   \[
    \mathbb E\tau_r^{-2}
    \le C\bigl(1+\log((rn)/r)\bigr)^{16}
    =C(1+\log n)^{16}.                                  \tag{0.9}
   \]

   Hence tensor multiplicity does not move the selected low mode out of the
   small-rank tail.  The exact low-mode martingale actually gives the reverse
   pressure

   \[
      |a|^2\le\lambda\,\mathbb E|M_1b|^2,
      \qquad b=a/|a|.                                   \tag{0.10}
   \]

   Thus a near-linear small-gap eigenfunction must select an expensive
   localization direction.

The slow-growth selection is therefore valid and stronger than doubling,
but neither tensor powers, strong-convexity regularization, nor Guan's rank
bound turns it into the two-sided mixed stationarity needed by the tensor
rigidity proof.

## 1. The dimension profile and slow-growth selection

### 1.1 Monotonicity and the reciprocal gap profile

Every isotropic law has \(C_P\ge1\), by testing the Poincare quotient on a
linear coordinate.  If \(\gamma_m\) is standard Gaussian in
\(\mathbb R^m\), tensorization gives

\[
 C_P(\mu\otimes\gamma_m)=\max(C_P(\mu),1)=C_P(\mu).
\]

Taking suprema proves

\[
 \mathcal C_{n+m}\ge\mathcal C_n.                       \tag{1.1}
\]

Thus \((\mathcal C_n)\) is nondecreasing.  If

\[
 \alpha_n=\inf\{\lambda_1(\mu):\mu\text{ isotropic log-concave in }
 \mathbb R^n\},
\]

then, with infima and suprema interpreted by approximating sequences,

\[
                         \alpha_n=\mathcal C_n^{-1}.     \tag{1.2}
\]

### 1.2 Doubling subsequence

**Lemma 1.1.**  Let \((c_n)\) be nondecreasing, unbounded, and satisfy
\(1\le c_n\le K\log(en)\).  Then there are \(n_j\to\infty\) such that

\[
                         {c_{2n_j}\over c_{n_j}}\to1.   \tag{1.3}
\]

**Proof.**  If (1.3) failed, there would be \(\varepsilon>0\) and
\(N\) such that

\[
 c_{2n}\ge(1+\varepsilon)c_n\qquad(n\ge N).
\]

Iteration gives

\[
 c_{2^kN}\ge(1+\varepsilon)^kc_N,
\]

whereas (0.2) gives

\[
 c_{2^kN}\le K(1+\log N+k\log2).
\]

The first expression grows exponentially in \(k\) and the second linearly,
a contradiction.  Since every ratio is at least one, a diagonal choice
with ratio at most \(1+1/j\) proves the claim.  Unboundedness and
monotonicity also give \(c_{n_j}\to\infty\). \(\square\)

### 1.3 Arbitrarily long dyadic plateaus

The same envelope gives substantially more than Lemma 1.1.

**Lemma 1.2.**  Under the assumptions of Lemma 1.1, there are integers
\(k_j,L_j\to\infty\) such that

\[
 {c_{2^{k_j+L_j}}\over c_{2^{k_j}}}\to1.               \tag{1.4}
\]

One may arrange

\[
 k_j\in[R_j,2R_j],\qquad
 L_j=\left\lfloor\sqrt{R_j/\log(R_j+1)}\right\rfloor    \tag{1.5}
\]

for any sufficiently rapidly increasing integers \(R_j\).

**Proof.**  Put \(b_k=\log c_{2^k}\).  Then \(b_k\) is nondecreasing and

\[
 0\le b_k\le \log K+\log(1+k\log2)\le C_0+\log(k+1).   \tag{1.6}
\]

Fix a large \(R\) and let

\[
 L=\left\lfloor\sqrt{R/\log(R+1)}\right\rfloor.
\]

Partition \([R,2R]\cap\mathbb Z\) into consecutive blocks of length \(L\),
discarding at most two incomplete blocks.  There are at least \(R/(2L)\)
complete blocks.  The sum of their increments is at most

\[
 b_{2R+L}-b_R\le C_0+\log(2R+L+1).
\]

Hence one complete block \([k,k+L]\) satisfies

\[
 0\le b_{k+L}-b_k
 \le {2L\over R}\,[C_0+\log(3R+1)]
 =O\!\left(\sqrt{\log R\over R}\right).                \tag{1.7}
\]

This tends to zero.  Exponentiating and taking \(R=R_j\to\infty\)
proves (1.4). \(\square\)

Set

\[
 n_j=2^{k_j},\qquad r_j=2^{L_j}.                        \tag{1.8}
\]

Then \(n_j,r_j\to\infty\), (1.4) is exactly (0.4), and monotonicity also
implies (0.3) along the same sequence.

## 2. Tensor powers are relative near-minimizers

Choose \(\varepsilon_j\downarrow0\) and isotropic log-concave
\(\mu_j\) in dimension \(n_j\) such that

\[
 C_P(\mu_j)\ge e^{-\varepsilon_j}\mathcal C_{n_j}.      \tag{2.1}
\]

Let

\[
 \nu_j=\mu_j^{\otimes r_j},\qquad N_j=r_jn_j.
\]

Tensorization gives

\[
 C_P(\nu_j)=C_P(\mu_j),
 \qquad
 \lambda_1(\nu_j)=C_P(\mu_j)^{-1}.                    \tag{2.2}
\]

The excess of \(\nu_j\) over the globally least spectral gap in dimension
\(N_j\), measured in the correct relative objective, is

\[
\begin{aligned}
 \Delta_j
 &:=
 \log\lambda_1(\nu_j)-\log\alpha_{N_j}\\
 &=\log{\mathcal C_{N_j}\over C_P(\mu_j)}\\
 &\le
 \log{\mathcal C_{r_jn_j}\over\mathcal C_{n_j}}
 +\varepsilon_j\longrightarrow0.                     \tag{2.3}
\end{aligned}
\]

This is the exact gain furnished by slow growth.  If the first gap of
\(\mu_j\) is attained by a normalized eigenfunction \(f_j\), the product
first eigenspace contains

\[
 F_i(x_1,\ldots,x_{r_j})=f_j(x_i),\qquad1\le i\le r_j. \tag{2.4}
\]

There is no attainment issue after a law is approximated by a smooth
strongly convex density on a bounded smooth convex domain: the weighted
Neumann form then has compact resolvent.  The \(r\)-fold product domain is
bounded and Lipschitz, so its product Neumann form also has compact
resolvent and the eigenfunctions in (2.4) remain exact, despite the product
corners.  Such approximations can retain (2.1) to arbitrary prescribed
accuracy by first choosing a compactly supported smooth Poincare witness
and passing its Rayleigh quotient through convolution, a small quadratic
potential, truncation, and affine isotropization.  What is not controlled
is the size of the resulting strong-convexity parameter; Section 5 isolates
that loss.

## 3. Exact \(r\)-fold mixed branch calculation

Fix one smooth isotropic factor \(\mu\), a normalized first eigenfunction
\(f\), and its first gap \(\lambda\):

\[
 \mathbb Ef=0,\qquad\mathbb Ef^2=1,\qquad
 \mathbb E|\nabla f|^2=\lambda.                         \tag{3.1}
\]

Put

\[
 a=\mathbb E[Xf],\qquad \beta=|a|^2\le1,
 \qquad \ell(x)=a\cdot x.                              \tag{3.2}
\]

Testing the eigenfunction equation against the coordinates gives

\[
                         m:=\mathbb E\nabla f=\lambda a.\tag{3.3}
\]

On \(\mu^{\otimes r}\), write \(f_i=f(x_i)\), \(\ell_i=\ell(x_i)\), and
let \(A\in\mathbb S^r\) have zero diagonal.

### 3.1 The unconstrained mixed perturbation

Define

\[
 W_A={1\over2}\sum_{i,j=1}^rA_{ij}f_if_j
     =\sum_{i<j}A_{ij}f_if_j.                           \tag{3.4}
\]

For the product first-eigenspace basis \(F_i=f_i\), the stress entries are

\[
 g_{ij}=\nabla F_i\cdot\nabla F_j-\lambda F_iF_j
 =-\lambda f_if_j\quad(i\ne j).                        \tag{3.5}
\]

The potential-variation branch matrix is

\[
 H(W)_{ij}=-\mathbb E[Wg_{ij}].                         \tag{3.6}
\]

Independence and centering give, on the displayed \(r\)-dimensional
principal block,

\[
                         H(W_A)=\lambda A.              \tag{3.7}
\]

All diagonal entries vanish.  Therefore the greater descent available
from the two signs is

\[
 \max\{-\lambda_{\min}H(W_A),
       -\lambda_{\min}H(-W_A)\}
 =\lambda\|A\|_{\rm op}.                               \tag{3.8}
\]

This calculation is valid even if the full first eigenspace is larger:
the bottom derivative of the full branch matrix is no larger than the
bottom derivative of a principal block.

### 3.2 The exact moment-tangent correction

The double-centered projection of \(f_i f_j\) onto total-degree-two
polynomials is \(\ell_i\ell_j\).  Hence

\[
 q_{ij}=f_if_j-\ell_i\ell_j                             \tag{3.9}
\]

is orthogonal to constants, all coordinates, and all centered quadratic
moment statistics in the \(rn\) variables.  In particular it is an exact
first-order tangent direction to the barycenter and covariance constraints.
Set

\[
 Q_A=\sum_{i<j}A_{ij}q_{ij}.                            \tag{3.10}
\]

Since \(\|\ell\|_2^2=\beta\) and
\(\langle f,\ell\rangle=\beta\),

\[
 \|q_{ij}\|_2^2=1-\beta^2.                             \tag{3.11}
\]

The same branch calculation now gives

\[
                  \boxed{H(Q_A)=\lambda(1-\beta^2)A.} \tag{3.12}
\]

Thus the moment correction removes exactly the degree-two part which the
tensor-factorization lemma is meant to detect.

For a full-dimensional smooth density, the Gram matrix of the coordinate
and centered quadratic statistics is positive definite.  The moment map is
therefore a submersion.  By the finite-dimensional implicit function
theorem, the tangent \(Q_A\) can be completed by an \(O(t^2)\) linear and
quadratic correction to an exactly centered and isotropic path.  This
correction does not alter (3.12).  If the base potential has a positive
Hessian margin and \(Q_A\) has bounded Hessian, both signs exist for a
nonzero interval whose length is controlled by that margin.

### 3.3 Exact mixed-Hessian cost

For \(i\ne j\), the \(x_i x_j\) Hessian block of \(q_{ij}\) is

\[
 R(x_i,x_j)=\nabla f(x_i)\nabla f(x_j)^T-aa^T.          \tag{3.13}
\]

Using (3.1), (3.3), and independence,

\[
\begin{aligned}
 \mathbb E\|R\|_{\rm HS}^2
 &=\lambda^2-2\lambda^2\beta^2+\beta^2\\
 &=\lambda^2+(1-2\lambda^2)\beta^2.                   \tag{3.14}
\end{aligned}
\]

The ordered off-diagonal Hessian blocks are orthogonal in the ambient
Hilbert--Schmidt sum.  Since
\(\sum_{i\ne j}A_{ij}^2=\|A\|_{\rm F}^2\), (3.13)--(3.14) imply

\[
 \boxed{
 \|D^2Q_A\|_{L^2(\mu^{\otimes r};\mathrm{HS})}
 \ge
 \sqrt{\lambda^2+(1-2\lambda^2)\beta^2}\,
 \|A\|_{\rm F}.}                                      \tag{3.15}
\]

Only mixed blocks were retained, so diagonal blocks can only increase the
left side.

If \(Q_A=k_+-k_-\) with convex \(k_\pm\), then pointwise

\[
 \|D^2k_++D^2k_-\|_{\rm HS}
 \ge\|D^2k_+-D^2k_-\|_{\rm HS}=\|D^2Q_A\|_{\rm HS},   \tag{3.16}
\]

because the two Hessians are positive semidefinite.  Thus (3.15) is also a
lower bound for every \(L^2\) difference-of-convex buffer.  Distributional
singular Hessian parts occur with the same sign in the total buffer and do
not reduce this bound.

Finally,

\[
 \lambda^2(1-\beta^2)^2
 \le \lambda^2+(1-2\lambda^2)\beta^2,                  \tag{3.17}
\]

because the difference between the right and left sides is
\(\beta^2(1-\lambda^2\beta^2)\ge0\).  Combining
(3.12), (3.15), and (3.17) yields the dimension-free no-amplification
inequality

\[
 \boxed{
 {\lambda(1-\beta^2)\|A\|_{\rm op}
  \over
  \inf\|D^2k_++D^2k_-\|_2}
 \le {\|A\|_{\rm op}\over\|A\|_{\rm F}}\le1.}          \tag{3.18}
\]

The infimum is over all convex decompositions of \(Q_A\), and is \(+\infty\)
if none exists.  Equation (3.18) is generous to the proposed route: it uses
only an \(L^2\) curvature cost, whereas actual global convexity requires
pointwise or measure-valued control.

The only zero-denominator case in (3.18) is
\(\lambda=\beta=1\), when \(f=\ell\), \(Q_A=0\), and the assertion is
vacuous.  In every nonzero mixed-residual case the displayed quotient is
well-defined.

## 4. The buffer has first-order spectral size

The same obstruction can be expressed entirely in the first-eigenspace
branch matrix.  Suppose an exact one-sided local minimum has a flat
potential, and write a difference-of-convex decomposition

\[
                        Q_A=k_+-k_-                    \tag{4.1}
\]

with \(k_\pm\) convex.  Each \(k_\pm\) is then a feasible one-sided
potential direction.  Local minimality gives

\[
 B_+:=H(k_+)\succeq0,\qquad B_-:=H(k_-)\succeq0.       \tag{4.2}
\]

By linearity and (3.12),

\[
 B_+-B_-=\lambda(1-\beta^2)A.                          \tag{4.3}
\]

The trace norm triangle inequality and positivity give

\[
\boxed{
 \operatorname {Tr}B_++\operatorname {Tr}B_-
 \ge\lambda(1-\beta^2)\|A\|_*.}                       \tag{4.4}
\]

Thus the total buffer first variation is of the same order as the mixed
splitting for every nonzero \(A\).  After division by \(\lambda\), as
required for stationarity of \(\log\lambda_1\), both are order one.

Equivalently, if a common convex buffer \(k\) makes both \(k+Q_A\) and
\(k-Q_A\) convex, exact one-sided stationarity gives

\[
 B\pm\lambda(1-\beta^2)A\succeq0,
 \qquad B=H(k),                                        \tag{4.5}
\]

and hence

\[
 \operatorname {Tr}B
 \ge\lambda(1-\beta^2)\|A\|_* .                       \tag{4.6}
\]

For (4.6), apply the trace-norm triangle inequality to the two PSD matrices
in (4.5), whose difference is \(2\lambda(1-\beta^2)A\).

Multiplicity does not improve these inequalities.  For example:

* one swap block has \(\|A\|_{\rm op}=1\),
  \(\|A\|_{\rm F}=\sqrt2\), and \(\|A\|_*=2\);
* \(r/2\) disjoint swap blocks still have operator norm one, but Frobenius
  norm \(\sqrt r\) and nuclear norm \(r\);
* for the complete graph \(A=(J-I)/(r-1)\), the operator norm is one and
  the Frobenius and nuclear norms remain of numerical order one;
* a typical sign matrix has
  \(\|A\|_{\rm op}/\|A\|_{\rm F}=O(r^{-1/2})\).

No choice makes the left side of (3.18) grow with \(r\).

## 5. Strong convexity and near-minimality: the exact scale

The phrase “regularize, then use both signs” hides a necessary comparison
of two small parameters.  The following abstract lemma records it without
assuming a second derivative of the spectral objective.

### Lemma 5.1 (Ekeland interior scale)

Let \(\mathcal X\) be a complete Banach chart of moment-normalized
potentials with norm controlling Hessians:

\[
                         \|D^2W\|_\infty\le\|W\|_\mathcal X.\tag{5.1}
\]

Let \(\mathcal K=\{V:D^2V\succeq0\}\) be closed in this chart, let
\(F(V)=\log\lambda_1(V)\) be lower semicontinuous and bounded below, and
suppose

\[
 D^2V\succeq\rho I,\qquad
 F(V)\le\inf_{\mathcal K}F+\Delta.                    \tag{5.2}
\]

For every \(\eta>0\), Ekeland's principle produces \(\widetilde V\in
\mathcal K\) such that

\[
 \|\widetilde V-V\|_\mathcal X\le{\Delta\over\eta},\qquad
 D^+F(\widetilde V)[W]\ge-\eta\|W\|_\mathcal X         \tag{5.3}
\]

for every feasible direction \(W\).  The point \(\widetilde V\) is still a
two-sided interior point with Hessian margin at least \(\rho/2\) provided

\[
                         {\Delta\over\eta}\le{\rho\over2}.\tag{5.4}
\]

There exists a choice with both \(\eta\to0\) and (5.4) precisely when

\[
                         {\Delta\over\rho}\to0.        \tag{5.5}
\]

For example, if (5.5) holds, take
\(\eta=\sqrt{\Delta/\rho}\) after fixing compatible units in the chart.

**Proof.**  Equation (5.3) is the standard Ekeland variational inequality.
The Hessian bound (5.1) gives

\[
 D^2\widetilde V\succeq
 \left(\rho-\|\widetilde V-V\|_\mathcal X\right)I,
\]

which proves (5.4).  Conditions \(\eta\to0\) and
\(\Delta/\eta=o(\rho)\) can hold simultaneously if and only if
\(\Delta/\rho\to0\). \(\square\)

At an interior point with a first-eigenvalue cluster, applying (5.3) to
both \(W\) and \(-W\) gives the full branch-matrix estimate

\[
 \left\|{H(W)\over\lambda}\right\|_{\rm op}
 \le\eta\|W\|_\mathcal X.                              \tag{5.6}
\]

Indeed the two one-sided derivatives of the bottom cluster are
\(\lambda_{\min}(H(W)/\lambda)\) and
\(-\lambda_{\max}(H(W)/\lambda)\).

For the mixed tangent \(Q_A\), (5.6) can overcome the normalized splitting
in (3.12) only if

\[
 \eta\|Q_A\|_{\mathcal X}
 =o\bigl((1-\beta^2)\|A\|_{\rm op}\bigr).              \tag{5.6a}
\]

Together with the interior requirement
\(\eta\ge2\Delta/\rho\), this forces the more precise necessary scale

\[
 {\Delta\over\rho}\,
 {\|Q_A\|_{\mathcal X}\over(1-\beta^2)\|A\|_{\rm op}}
 \longrightarrow0.                                    \tag{5.6b}
\]

Even under the optimistic assumption that the second factor is numerical,
one still needs \(\Delta=o(\rho)\).  A pointwise Hessian norm can make that
factor grow; the \(L^2\) calculation of Section 3 cannot bound it from
above.

Lemma 5.1 identifies two independent gaps in the proposed use of (2.3).

1. Smooth strongly convex approximation proves only that the restricted
   infimum tends to the unrestricted infimum as \(\rho\downarrow0\).  It
   gives no rate \(\Delta=o(\rho)\).  The scalar law

   \[
          \inf_{D^2V\succeq\rho I}F(V)
          =\inf_{D^2V\succeq0}F(V)+c\rho              \tag{5.7}
   \]

   is fully consistent with such convergence and violates (5.5).  The
   first-order buffer estimates (3.18) and (4.4) show why a linear loss is
   the natural, rather than an artificially pessimistic, scale.

2. Ekeland's point \(\widetilde V\), obtained in the full potential space,
   need not remain a tensor product.  A generic arbitrarily small
   perturbation splits the \(r\)-fold first eigenvalue, so (5.6) no longer
   applies to the original \(r\times r\) mixed stress matrix.  Restricting
   Ekeland to product potentials preserves the cluster but removes all
   mixed directions.  Restricting only to permutation-invariant potentials
   retains a small representation-theoretic part of the cluster, not the
   full collection of off-diagonal entries required by tensor
   factorization.

Near-minimality by itself gives no derivative estimate at the original
product.  The scalar cusp family

\[
 F_\Delta(t)=F_\Delta(0)-\min\{|t|,\Delta\}             \tag{5.8}
\]

has objective excess \(\Delta\to0\) but order-one descending directional
derivative in both signs at zero.  A uniform Taylor radius or the Ekeland
interior scale is additional information, not a consequence of (2.3).

## 6. Guan stopping subspaces do not select the low-mode direction

This section uses the exact low-mode notation from the parallel-coupling
route.  Along stochastic localization let

\[
 M_0=I,\qquad M_t'=A_tM_t,
\]

where \(A_t\) is the posterior covariance.  If \(f\) is a normalized first
eigenfunction and

\[
 a=\mathbb E[Xf],\qquad b={a\over|a|},
\]

the localized covariance \(c_t=\operatorname {Cov}_{\mu_t}(X,f)\) obeys

\[
 M_t^Tc_t\text{ is a martingale},\qquad
 \mathbb E|c_1|^2\le\lambda.                            \tag{6.1}
\]

Therefore

\[
 |a|=\mathbb E\langle M_1b,c_1\rangle,
 \qquad
 \boxed{|a|^2\le\lambda\mathbb E|M_1b|^2.}             \tag{6.2}
\]

This is a lower, not an upper, bound on the cost of the selected direction.

Let the eigenvalues of \(A_t\) be ordered decreasingly and put

\[
 \tau_k=\inf\{t:\lambda_k(A_t)\ge3\}.
\]

Guan's stopping estimate has the form

\[
 \mathbb E\tau_k^{-2}
 \le C\left(1+\log{n\over k}\right)^{16}.              \tag{6.3}
\]

It controls the number of exceptional covariance directions, but contains
no spectral weights of \(b\) or of the adapted vector \(M_tb\).  A
directional conclusion would require a new inequality such as

\[
 \mathbb E\|P_k(t)M_tb\|^2
 \lesssim {k\over n}\,\mathbb E|M_tb|^2,               \tag{6.4}
\]

where \(P_k(t)\) projects onto the top \(k\) covariance eigenspace.  Neither
the eigenfunction equation, (6.1), nor tensor stationarity implies (6.4).
In fact, if \(w_t=M_tb\) and \(q_t=w_t/|w_t|\), then wherever \(w_t\ne0\),

\[
 q_t'=\bigl(A_t-\langle q_t,A_tq_t\rangle I\bigr)q_t.  \tag{6.4a}
\]

Thus the product integral increases, rather than suppresses, the relative
weight of a currently high-covariance eigendirection.

### 6.1 Tensor powers preserve the bad rank fraction

For an \(r\)-fold product, the product localization posterior remains a
product and

\[
 A_t=\operatorname {diag}(A_t^{(1)},\ldots,A_t^{(r)}),
 \qquad
 M_t=\operatorname {diag}(M_t^{(1)},\ldots,M_t^{(r)}). \tag{6.5}
\]

The \(r\) low-mode directions are

\[
 b_i=(0,\ldots,0,b,0,\ldots,0),\qquad1\le i\le r.      \tag{6.6}
\]

They span a subspace of dimension \(r\) inside ambient dimension \(rn\),
so their rank fraction is exactly

\[
                              {r\over rn}={1\over n}.   \tag{6.7}
\]

More invariantly, identify the ambient coordinate space with
\(\mathbb R^r\otimes\mathbb R^n\).  Every linear combination of the tensor
low modes has physical linear projection in

\[
 \mathcal S_b=\mathbb R^r\otimes\operatorname {span}\{b\},
 \qquad \dim\mathcal S_b=r.                            \tag{6.7a}
\]

The mixed branch matrix \(A\) from Section 3 acts only on the
\(\mathbb R^r\) factor.  Varying or extremizing over \(A\) changes the
coefficient vector in \(\mathbb R^r\), but cannot rotate the physical
factor \(b\).  Hence the rank-\(r\) projector
\(I_r\otimes bb^T\) may contain the linear projection of every tensor
low mode simultaneously, regardless of how the mixed extremizer is
selected.

If each factor develops one covariance spike aligned with its \(b_i\), the
global exceptional subspace has rank \(r\) and contains every selected
low-mode direction.  Guan's estimate at that rank is exactly (0.9); it is
independent of \(r\).  For the permutation-symmetric direction

\[
 \bar b=r^{-1/2}\sum_i b_i,
\]

block diagonality gives

\[
 \mathbb E|M_1\bar b|^2
 ={1\over r}\sum_i\mathbb E|M_1^{(i)}b|^2
 =\mathbb E|M_1b|^2.                                   \tag{6.8}
\]

Thus averaging across tensor factors does not dilute the expensive
direction.

The one-spike countermodel can be tensorized verbatim.  If one factor has
an exceptional rank-one path with

\[
 |M_1b|^2\asymp(1+\log n)^{16},
\]

then \(r\) copies have rank \(r\), trace contribution
\(r(1+\log n)^{16}\), and the same directional amplification for every
\(b_i\) and for \(\bar b\).  This is compatible with the global trace scale
\(O(rn)\) whenever the one-factor model is, and it saturates the right side
of (0.9).  It is an information-theoretic covariance-path model, not an
assertion that a KLS counterexample exists; it proves that the stopping-time
data plus tensor symmetry do not imply a selector.

### 6.2 Exact rank bookkeeping for a low eigenspace

Suppose more generally that a first eigenspace has \(m\) orthonormal modes
whose linear projections give \(m\) orthonormal directions
\(b_1,\ldots,b_m\), each with squared projection at least \(1-\delta^2\).
The trace bound for the parallel coupling yields

\[
 m(1-\delta^2)
 \le C\lambda n,                                      \tag{6.9}
\]

and therefore

\[
 \lambda\ge {m\over n}{1-\delta^2\over C}.            \tag{6.10}
\]

For an \(r\)-fold tensor of a single \(n\)-dimensional mode, \(m=r\) and
the ambient dimension is \(rn\).  Formula (6.10) remains only

\[
                         \lambda\ge {1-\delta^2\over Cn},\tag{6.11}
\]

with no improvement as \(r\to\infty\).

### 6.3 Abstract selector criterion

Rank information alone cannot choose an extremizing direction.  If a
linear operator \(T:\mathcal H\to\mathbb R^n\) sends an extremizer \(B\)
to its low-mode coordinate \(a=TB/|TB|\), then an exact extremizer aligned
with a subspace \(E\) exists only if

\[
 \|P_ET\|=\|T\|                                       \tag{6.12}
\]

and the top right-singular space contains a vector whose \(T\)-image lies
in \(E\).  Quantitatively, if

\[
 \|TB\|\ge(1-\eta)\|T\|\|B\|,
 \qquad \|P_Ea\|^2\ge1-\zeta^2,
\]

then necessarily

\[
 \|P_ET\|\ge(1-\eta)\sqrt{1-\zeta^2}\,\|T\|.           \tag{6.13}
\]

A rank bound on \(E\) says nothing about (6.12) or (6.13).  The algebraic
counterexample \(T(B)=\langle B,B_0\rangle e_n\) and
\(E=\operatorname {span}\{e_1,\ldots,e_k\}\), \(k<n\), has complete
misalignment.  The directional third-moment operator of an actual isotropic
product law realizes the same pattern: take \(n-1\) Gaussian factors and
one centered variance-one exponential factor.  Then

\[
 T(B)_i=\mathbb E X_i^3\,B_{ii},
\]

so its unique nonzero left-singular direction is the exponential
coordinate.  A separate rank estimate on a localization covariance
subspace cannot align it without a new correlation theorem.

## 7. Definitive conclusion for this route

The \(O(\log n)\) KLS bound does give exactly the slow-growth extremizing
dimensions required by the tensor idea, and even permits tensor powers of
unbounded multiplicity with vanishing relative optimization defect.  This
part of the route is complete.

The remaining inference fails at a quantified scale:

\[
\begin{array}{c}
 \text{mixed branch splitting}\sim
 \lambda(1-\beta^2)\|A\|_{\rm op},\\[2mm]
 \text{least convex-buffer curvature}\gtrsim
 \sqrt{\lambda^2+(1-2\lambda^2)\beta^2}\|A\|_{\rm F},\\[2mm]
 \text{buffer branch action}\gtrsim
 \lambda(1-\beta^2)\|A\|_*.
\end{array}                                             \tag{7.1}
\]

Every ratio is dimension-free and none improves with multiplicity.  To
turn strong-convexity regularization into full two-sided stationarity would
additionally require \(\Delta=o(\rho)\) and preservation of the tensor
eigenvalue cluster.  Slow growth supplies neither.  Guan's theorem controls
the rank of exceptional localization directions, but \(r\) tensor low modes
may all lie in a rank-\(r\) exceptional subspace, leaving the critical ratio
\(r/(rn)=1/n\) unchanged.

Therefore tensor powers do not amplify the mixed splitting past the convex
buffer or the Guan small-rank tail.  Any successful continuation needs a
new tensor-stable variational principle that controls the normalized
normal-cone action itself, or a new correlation theorem coupling the
first-eigenfunction direction \(b\) to the evolving covariance eigenspaces.
