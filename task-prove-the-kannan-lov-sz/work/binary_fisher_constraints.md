# Binary Fisher constraints for a hypothetical bad separator

## 1. Gaussian observation notation

Let `X` have a centered isotropic log-concave law `mu` on `R^n`, let
`B=1_S` with `P(B=1)=p in (0,1)`, and, for `t>0`, observe

\[
 C_t=tX+\sqrt t\,G,\qquad G\sim N(0,I)                 \tag{1.1}
\]

independently.  Conditional on `C_t=c`, the law of `X` is the ordinary
stochastic-localization posterior

\[
 d\mu_{t,c}(x)\propto
   \exp\{c\cdot x-t|x|^2/2\}\,d\mu(x).              \tag{1.2}
\]

Write

\[
 \begin{aligned}
 g&=P(B=1\mid C_t),&m&=E[X\mid C_t],\\
 A&=Cov(X\mid C_t),&v&=Cov(B,X\mid C_t).
 \end{aligned}                                      \tag{1.3}
\]

When `v ne0`, put `u=v/|v|` and

\[
 \eta={\sqrt t|v|\over I(g)}\in[0,1].               \tag{1.4}
\]

The upper bound is the sharp centroid inequality for the `t`-strongly
log-concave posterior.

## 2. Exact matrix conditional-covariance bound

For every posterior state,

\[
 \boxed{{vv^T\over g(1-g)}\preceq A.}                \tag{2.1}
\]

Indeed, for every vector `theta`, scalar covariance Cauchy--Schwarz gives

\[
 \langle v,\theta\rangle^2
 \le Var(B\mid C_t)Var(\langle X,\theta\rangle\mid C_t)
 =g(1-g)\theta^TA\theta.
\]

The law of total covariance and isotropy give

\[
 E A+Cov(m)=I,
\]

and hence, after averaging (2.1),

\[
 \boxed{
 R_t:=E\left[{\eta^2 I(g)^2\over g(1-g)}uu^T\right]
 =tE\left[{vv^T\over g(1-g)}\right]\preceq tI.}     \tag{2.2}
\]

This is an operator inequality, not merely a trace estimate.  It remains
valid without log-concavity; log-concavity enters only through `eta<=1`.

It follows that

\[
 \boxed{
 \operatorname {rank}_{eff}(R_t):={tr R_t\over\|R_t\|_{op}}
 \ge {1\over t}\,tr R_t.}                           \tag{2.3}
\]

Thus any mechanism forcing `tr R_t` to be a fixed positive constant at
time `t=K^{-1}` simultaneously forces at least order `K` mutually
dispersed active directions.

## 3. Exact binary I--MMSE identity

Let

\[
 \mathsf I(t)=I(B;C_t)
 =h(p)-E h(g),
\]

where `h(a)=-a log a-(1-a)log(1-a)`.  Then

\[
 \boxed{
 \mathsf I'(t)
 ={1\over2}E{|v|^2\over g(1-g)}
 ={1\over2t}\operatorname {tr}R_t.}                 \tag{3.1}
\]

One proof differentiates the two Gaussian channel densities.  Equivalently,
the vector I--MMSE identity for the pair `(B,X)` gives

\[
 {d\over dt}I(B;C_t)
 ={1\over2}\{E|X-E[X\mid C_t]|^2
 -E|X-E[X\mid B,C_t]|^2\},
\]

and the conditional law of total variance turns the expression in braces
into `E|v|^2/[g(1-g)]`.  Approximation by bounded `X`, followed by monotone
convergence of mutual information and `L^2` convergence of conditional
means, justifies the identity under the finite second moment assumed here.

Since `0<=mathsf I(t)<=h(p)<=log 2` for a balanced label,

\[
 \boxed{\int_0^\infty {tr R_t\over2t}\,dt\le\log2.} 
                                                               \tag{3.2}
\]

This is exactly a one-bit budget.

## 4. Relation to Gaussian-profile localization

Put

\[
 J(t)=E I(g_t),\qquad F(t)=\sqrt t\,J(t).
\]

The ordinary localization profile identity says

\[
 F'(t)={1\over2\sqrt t}
 E\left[I(g_t)(1-\eta_t^2)\right].                  \tag{4.1}
\]

For all `a in (0,1)`,

\[
 {I(a)\over a(1-a)}\ge c_0>0,                       \tag{4.2}
\]

with a numerical constant `c_0`.  Consequently

\[
 tr R_t
 =E\left[{\eta_t^2I(g_t)^2\over g_t(1-g_t)}\right]
 \ge c_0E[\eta_t^2 I(g_t)].                         \tag{4.3}
\]

Suppose on some time `t` that

\[
 J(t)\ge j_0,
 \qquad E[I(g_t)(1-\eta_t^2)]\le\varepsilon j_0.   \tag{4.4}
\]

Then (4.3) gives

\[
 tr R_t\ge c_0(1-\varepsilon)j_0,                  \tag{4.5}
\]

and (2.3) gives the explicit effective-rank constraint

\[
 \boxed{
 \operatorname {rank}_{eff}(R_t)
 \ge {c_0(1-\varepsilon)j_0\over t}.}               \tag{4.6}
\]

For a hypothetical Poincare constant `K>>1`, the natural bad-separator
plateau has `t=K^{-1}`, `J(t)` of constant order, and small relative profile
drift.  Equation (4.6) then forces effective rank `Omega(K)`.  A single
coherent direction would instead make `R_t` nearly rank one; (2.2) would
immediately contradict `tr R_t` of constant order when `t=K^{-1}`.

## 5. Why scalar information alone does not close the argument

On an idealized bad profile plateau,

\[
 J(t)\asymp {1\over\sqrt{Kt}},\qquad \eta_t\asymp1,
 \qquad K^{-1}\lesssim t\lesssim1.                 \tag{5.1}
\]

Equations (3.1), (4.2), and the matching upper profile estimates give the
scale

\[
 \mathsf I'(t)\asymp {1\over\sqrt K\,t^{3/2}}.      \tag{5.2}
\]

Its integral over `[K^{-1},1]` is of constant order:

\[
 \int_{1/K}^1{dt\over\sqrt K\,t^{3/2}}
 =2(1-K^{-1/2}).                                    \tag{5.3}
\]

Thus the hypothetical plateau consumes exactly one bit, not more than one
bit.  The scalar entropy budget (3.2) can recover the usual equivalence
scale but cannot yield a power improvement in `K`.

The matrix constraint is strictly sharper: it says that the one bit must be
assembled from `Omega(K)` orthogonal weak directional channels.  A complete
argument must now use global log-concavity to prove one of the following:

1. such a high-effective-rank posterior-label Fisher matrix has an
   additional profile/heat-Hessian cost; or
2. it forces an approximate product decomposition, after which
   Bobkov--Houdre tensorization reduces to a lower-dimensional factor.

Neither conclusion follows from (2.2)--(3.2) alone.  Gaussian parity gives
the diagnostic high-rank one-bit pattern, but its two label classes are not
log-concave and it is not a counterexample within the admissible class.
