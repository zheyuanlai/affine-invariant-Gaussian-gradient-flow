# Gaussian-selector information and angular Sobolev control

## 0. Outcome

This note records a dimension-free regularity statement for the special
physical selector in `fixed_scale_physical_splicing.md`.  It uses only the
Gaussian channel and is valid for an arbitrary bounded selector; no
log-concavity of the selector is assumed.

Let `X` have an arbitrary probability law `mu`, let `G` be standard Gaussian,
and put `Y=X+sqrt(s)G`.  If `0<=h<=M`, define

\[
 r(x)=E_Gh(x+\sqrt sG),\qquad b=\int r\,d\mu>0.
\]

If `u:R^n->R^d` is measurable with `|u|<=1`, put

\[
 m(x)=E_G[h(Y)u(Y)\mid X=x],\qquad v(x)=m(x)/r(x).
\]

Then

\[
 \boxed{\int {|\nabla r|^2\over r}\,d\mu
 \le {2b\over s}\log {M\over b},}                    \tag{0.1}
\]

and

\[
 \boxed{\int r\|\nabla v\|_{HS}^2d\mu
 \le {b\over s}\left(2+4\log {M\over b}\right).}    \tag{0.2}
\]

Consequently the ordinary `L^2(mu)` feature

\[
                         \Phi=\sqrt r\,v={m\over\sqrt r}
\]

obeys

\[
 \boxed{\int\|\nabla\Phi\|_{HS}^2d\mu
 \le {b\over s}\left(5+10\log {M\over b}\right).}   \tag{0.3}
\]

For the heat selector `h=1_G|W|`, one has `M<=I(1/2)/sqrt(s)`.
The good-direction Fisher matrix survives conditional averaging: if its
trace is `b`, its effective rank is `R`, and the directional alignment loss
is `Delta_G`, then

\[
 \operatorname {tr}\int rvv^T\,d\mu\ge b-2\Delta_G,
 \qquad
 \left\|\int rvv^T\,d\mu\right\|_{op}\le {b\over R}. \tag{0.4}
\]

Thus the fixed-scale construction produces a genuinely high-rank smooth
vector field on the original space.  This does not close KLS: applying the
only available Poincare inequality for `mu` to (0.3) multiplies the right
side by `K/s=1/alpha`.  The channel estimate therefore sharpens the selector
regularity issue but, by itself, reproduces the fixed-scale power
obstruction.

## 1. Conditional entropy identities

For a fixed `x`, write

\[
 h_x(g)=h(x+\sqrt s g),\qquad
 dP_x(g)={h_x(g)\over r(x)}d\gamma(g),
 \qquad D_x=D(P_x\|\gamma).
\]

Gaussian score differentiation, first for bounded smooth `h` and then by
Gaussian mollification, gives

\[
 \nabla r={r\over\sqrt s}E_{P_x}G,
 \qquad
 \nabla v={1\over\sqrt s}\operatorname {Cov}_{P_x}(u,G). \tag{1.1}
\]

The average conditional entropy has the exact expression

\[
 \int rD_xd\mu
 =E[h(Y)\log h(Y)]-\int r\log r\,d\mu.              \tag{1.2}
\]

Since `h<=M`, the first term is at most `b log M`; convexity of
`t log t` and `int r=b` make the second term at least `b log b`.  Hence

\[
 \boxed{\int rD_xd\mu\le b\log(M/b).}               \tag{1.3}
\]

This remains valid when `h` vanishes, with `0 log 0=0`.  Since `b<=M`, the
right side is nonnegative.

## 2. Scalar Fisher bound

The entropy variational formula, tested with linear Gaussian exponentials,
gives

\[
 D_x\ge {1\over2}|E_{P_x}G|^2.                       \tag{2.1}
\]

Indeed

\[
 D(P\|\gamma)\ge a\cdot E_PG-\log E_\gamma e^{a\cdot G}
 =a\cdot E_PG-|a|^2/2,
\]

and optimization over `a` proves (2.1).  Equations (1.1), (1.3), and
(2.1) prove (0.1).

## 3. Vector angular bound

Put

\[
 C_x=\operatorname {Cov}_{P_x}(u,G),\qquad
 \Sigma_x=\operatorname {Cov}_{P_x}(G).
\]

Covariance Cauchy--Schwarz gives the matrix inequality

\[
 C_xC_x^T\preceq\lambda_{max}(\Sigma_x)
                    \operatorname {Cov}_{P_x}(u).
\]

Because `|u|<=1`,

\[
 \|C_x\|_{HS}^2
 \le\lambda_{max}(\Sigma_x)\operatorname {tr}
          \operatorname {Cov}_{P_x}(u)
 \le\lambda_{max}(\Sigma_x).                       \tag{3.1}
\]

For any unit `a`, data processing under `g mapsto a dot g`, followed by
the fact that a Gaussian maximizes entropy at fixed mean and variance, gives

\[
 D_x\ge {1\over2}\{m_a^2+\sigma_a^2-1-\log\sigma_a^2\}, \tag{3.2}
\]

where `m_a=E_P<a,G>` and
`sigma_a^2=Var_P<a,G>`.  If `sigma_a^2<=2` there is nothing to prove.  If
`sigma_a^2>2`, then `log sigma_a^2<=sigma_a^2/2`, so (3.2) yields

\[
                         \sigma_a^2\le4D_x+2.
\]

Taking the supremum in `a`, (3.1) gives

\[
                         \|C_x\|_{HS}^2\le4D_x+2.   \tag{3.3}
\]

Equations (1.1), (1.3), and (3.3) prove (0.2).  The bound is independent of
both the input dimension and the output dimension.

For (0.3), differentiate `Phi=sqrt(r)v` and use

\[
 \|\nabla\Phi\|_{HS}^2
 \le2r\|\nabla v\|_{HS}^2+{|\nabla r|^2\over2r}|v|^2,
 \qquad |v|\le1.
\]

The slightly enlarged integer constants in (0.3) follow from (0.1)--(0.2).

## 4. Rank survives the conditional mean

Let

\[
 B=E[h(Y)u(Y)u(Y)^T],\qquad \operatorname {tr}B=b,
 \qquad R={b\over\|B\|_{op}}.
\]

Suppose a unit field `theta(X)` satisfies

\[
 E[h(Y)|u(Y)-\theta(X)|^2]\le2\Delta_G.             \tag{4.1}
\]

Conditionally on `X`, `v=E_{P_x}u`.  Therefore

\[
 B-\int rvv^Td\mu
 =\int r\operatorname {Cov}_{P_x}(u)d\mu\succeq0.   \tag{4.2}
\]

Its trace is

\[
 E[h|u-v|^2]
 \le E[h|u-\theta|^2]\le2\Delta_G,                 \tag{4.3}
\]

because a conditional mean minimizes conditional squared error.  Equations
(4.2)--(4.3) prove the trace inequality in (0.4), while the PSD order in
(4.2) gives the operator inequality.

At the fixed numerical scale, `b>=b_0p` with `b_0>.00517`,
`Delta_G<=Delta<6.02*10^-5p`, and `R>=500`.  Hence

\[
 {\operatorname {tr}\int rvv^Td\mu\over
       \|\int rvv^Td\mu\|_{op}}
 \ge R\left(1-{2\Delta\over b}\right)>488.          \tag{4.4}
\]

Thus neither the physical transfer nor conditional angular averaging loses
the high-rank seed.

## 5. Fixed-scale audit and remaining obstruction

In the selector application,

\[
                         M={I(1/2)\over\sqrt s}.
\]

Moreover `b>=b_0p`, and Buser--Ledoux plus
`p>=psi_mu/2` gives `p sqrt(K)>=c` with a universal positive `c`.  Since
`s=alpha K`,

\[
                         \log(M/b)\le C+{1\over2}\log(1/alpha). \tag{5.1}
\]

Apply the Poincare inequality of `mu` componentwise to `Phi`.  Even after
subtracting its vector mean, (0.3) gives only

\[
 \int|\Phi-E\Phi|^2d\mu
 \le {CbK\over s}\{1+\log(1/alpha)\}
 ={Cb\over\alpha}\{1+\log(1/alpha)\}.              \tag{5.2}
\]

The high-rank covariance in (4.4) has trace of order `b`, so (5.2) is much
larger than the scale which would contradict it when `alpha` is small.
Effective rank does not improve a componentwise Poincare trace estimate.

The same issue prevents (0.2) from controlling the cross-level matrices.
The flux probability `r dmu/b` is generally not log-concave and has no
known dimension-free Poincare inequality.  Its pushforward by the scalar
level `F_0` may put mass on separated level intervals, allowing its
conditional normal projector to switch in a low-flux gap.  A lower bound on
that one-dimensional occupation density, a full `BV` focal charge, or a
selector-faithfulness theorem would be additional input, not a consequence
of (0.1)--(0.4).

The proved content is therefore the dimension-free channel regularity and
rank preservation, not a phase-charge theorem.
