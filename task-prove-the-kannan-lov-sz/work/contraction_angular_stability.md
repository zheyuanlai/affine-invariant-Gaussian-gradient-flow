# Angular stability under a Brenier contraction: the Gaussian-halfspace case

## 1. Normalization and result

The proposed contraction theorem is homogeneous.  This note first proves the
following exact special case, which is the part of the theorem in which the
pullback of the witnessing set is already a Gaussian halfspace.

Let `G` be standard Gaussian in `R^n`.  Let

\[
 T=\nabla\Phi:\mathbb R^n\longrightarrow\mathbb R^n,
 \qquad 0\preceq \nabla T\preceq I\quad\hbox{a.e.},                 \tag{1.1}
\]

and put `X=T(G)-E T(G)`.  Fix a unit vector `e`, a number
`g in [delta,1-delta]`, and the upper Gaussian `g`-quantile `c`, so that

\[
 B=\{z:\langle z,e\rangle\ge c\},\qquad \gamma_n(B)=g.
\]

Write

\[
 h=\mathbf1_B-g,\quad v=\mathbb E[hX],\quad
 u={v\over |v|},\quad P=I-u u^T,quad
 D=\mathbb E[hXX^T].                                  \tag{1.2}
\]

**Theorem 1 (halfspace-pullback contraction theorem).**  If `v ne 0`, then

\[
 \boxed{\|PD\|_{HS}^2\le C_\delta\{\mathcal I(g)-|v|\}.}          \tag{1.3}
\]

If `T` has Lipschitz constant `L` instead of one, scaling `T/L` gives

\[
 \boxed{\|PD\|_{HS}^2
 \le C_\delta L^3\{L\mathcal I(g)-|v|\}.}                         \tag{1.4}
\]

This includes every nonlinear contraction when the target set has a
Gaussian-halfspace pullback.  Section 7 separately handles every set under a
linear positive-semidefinite contraction.  The present theorem does not
assume that the pushforward is absolutely continuous.

## 2. The Gaussian threshold flux

Rotate the source so that `e=e_1`, write `G=(S,Z)`, and let `varphi` and
`Phi` denote the one-dimensional Gaussian density and distribution
function.  Define

\[
 k(s)=\begin{cases}
  g\,\Phi(s)/\varphi(s),&s<c,\\
  (1-g)\{1-\Phi(s)\}/\varphi(s),&s\ge c.
 \end{cases}                                                       \tag{2.1}
\]

Then `k` is positive and continuous and

\[
 (k\varphi)'=-h\varphi,
 \qquad \mathbb E k(S)=\mathbb E[hS]=\varphi(c)=\mathcal I(g).    \tag{2.2}
\]

Consequently, for every locally Lipschitz `F` for which the two sides are
integrable,

\[
 \mathbb E[hF(S,Z)]=\mathbb E[k(S)\partial_sF(S,Z)].               \tag{2.3}
\]

For central `g`, Mills' bounds imply

\[
 0<k(s)\le C_\delta.                                                \tag{2.4}
\]

We also need a weighted covariance inequality.

**Lemma 2 (threshold-flux covariance).**  Put `I_0=I(g)`.  For every
locally absolutely continuous `f:R -> R^m`,

\[
 \left|\mathbb E[(k(S)-I_0)f(S)]\right|^2
 \le C_\delta\mathbb E[k(S)|f'(S)|^2].                             \tag{2.5}
\]

**Proof.**  Define `a` by

\[
 a(s)\varphi(s)=\int_s^\infty(k(r)-I_0)\varphi(r)dr.                \tag{2.6}
\]

The mean-zero identity in (2.2) also gives the corresponding integral from
`-infinity` to `s`.  Integration by parts and Cauchy--Schwarz give

\[
 \begin{aligned}
 \mathbb E[(k-I_0)f]
 &=\mathbb E[a f'],\\
 \left|\mathbb E[a f']\right|^2
 &\le \mathbb E[a^2/k],\mathbb E[k|f'|^2].                      \tag{2.7}
 \end{aligned}
\]

It remains only to bound the first factor.  Since `g` is central, `|c|` is
bounded in terms of `delta`.  On a fixed interval containing `c`, `k` is
bounded below and `a` is bounded, uniformly in `g in [delta,1-delta]`.
Outside that interval, the two-sided Mills bounds give

\[
 {c_\delta\over1+|s|}\le k(s)\le {C_\delta\over1+|s|},
 \qquad |a(s)|\le {C_\delta\over1+|s|}.                            \tag{2.8}
\]

Thus `a(s)^2/k(s)<=C_delta/(1+|s|)` in the tails and
`E[a^2/k]<=C_delta`.  The vector-valued assertion follows componentwise,
or directly by Hilbert-space Cauchy--Schwarz.  This proves (2.5).

## 3. Matrix contraction rigidity

Let

\[
 H(z)=\nabla T(z),\qquad q(z)=H(z)e.                                \tag{3.1}
\]

The assumptions imply `0 preceq H preceq I`, hence `|q|<=1`.  Applying
(2.3) componentwise gives

\[
 v=\mathbb E[kq].                                                    \tag{3.2}
\]

It follows that the centroid deficit has the exact nonnegative form

\[
 \epsilon:=\mathcal I(g)-|v|
 =\mathbb E\{k(1-u\cdot q)\}.                                      \tag{3.3}
\]

Two pointwise consequences are

\[
 |Pq|^2\le 2(1-u\cdot q),
 \qquad |q-u|^2\le2(1-u\cdot q).                                  \tag{3.4}
\]

Indeed both follow from `|q|<=1`.  Notice that this is stronger than merely
using Lipschitzness: symmetry and positivity of `H` identify the same source
and target direction in the equality case.

The covariance of `X` is bounded by the identity.  In fact, for every unit
`theta`, Gaussian Poincare and (1.1) give

\[
 \operatorname {Var}\langle\theta,T(G)\rangle
 \le\mathbb E|H\theta|^2\le1.                                     \tag{3.5}
\]

We use the following elementary Bessel consequence of (3.5): for every
square-integrable random vector `F`,

\[
 \|\mathbb E[F X^T]\|_{HS}^2\le\mathbb E|F|^2.                   \tag{3.6}
\]

To verify it, apply scalar Cauchy--Schwarz to each row: if `F_j` is one
component, then

\[
 |\mathbb E[F_jX]|^2
 =\sup_{|theta|=1}|\operatorname {Cov}(F_j,\langle theta,X\rangle)|^2
 \le\operatorname {Var}(F_j),
\]

and sum in `j`.

## 4. Second moment and proof of Theorem 1

Apply (2.3) to the matrix-valued function `XX^T`.  Since
`partial_s X=q`,

\[
 D=\mathbb E[k(qX^T+Xq^T)].                                       \tag{4.1}
\]

After projecting its first row space,

\[
 PD=A+B_1+B_2,                                                      \tag{4.2}
\]

where

\[
 A=\mathbb E[k(Pq)X^T],\qquad
 B_1=\mathbb E[k(PX)(q-u)^T],\qquad
 B_2=\mathbb E[kPX]u^T.                                            \tag{4.3}
\]

By (3.6), (2.4), (3.4), and (3.3),

\[
 \begin{aligned}
 \|A\|_{HS}^2
 &\le\mathbb E[k^2|Pq|^2]\le C_\delta\epsilon,\\
 \|B_1\|_{HS}^2
 &=\|\mathbb E[k(q-u)(PX)^T]\|_{HS}^2\\
 &\le\mathbb E[k^2|q-u|^2]\le C_\delta\epsilon.                 \tag{4.4}
 \end{aligned}
\]

For the remaining vector, use `E PX=0`, condition on `Z`, and apply Lemma 2
to `f_Z(s)=PX(s,Z)`.  Since `partial_s(PX)=Pq`, Jensen gives

\[
 \begin{aligned}
 |\mathbb E[kPX]|^2
 &=\left|\mathbb E_Z\mathbb E_S[(k-I_0)PX\mid Z]\right|^2\\
 &\le C_\delta\mathbb E[k|Pq|^2]
 \le C_\delta\epsilon.                                            \tag{4.5}
 \end{aligned}
\]

Thus `||B_2||HS^2<=C_delta epsilon`.  Squaring the triangle inequality in
(4.2) proves (1.3).

For (1.4), replace `T` by `T/L`.  The normalized correlated second moment is
`D/L^2` and its centroid is `v/L`; multiplying the normalized estimate by
`L^4` gives exactly (1.4).

## 5. Equality and what the proof says about splitting

If `epsilon=0`, then (3.4)--(4.5) give `PD=0` directly.  More structurally,
(3.3) forces `q=u` almost everywhere because `k>0`.  Thus `H e=u` a.e.
For a symmetric positive contraction this can happen only when `e=u` and
`He=e`: indeed the maximum of `a^THb`, over `0 preceq H preceq I`, is
`(1+a\cdot b)/2`, so `u\cdot He=1` forces `u=e`.  Symmetry then gives
`PHe=0` and the map splits off the active Gaussian coordinate.

The quantitative proof shows precisely why a varying transverse conditional
law is harmless in this halfspace-pullback case.  The two terms containing
`PHe` are paid for by the pointwise contraction defect.  The only apparently
unpaid term is `E[kPX]`; Lemma 2 and the identity
`partial_e(PT)=PHe` pay for it by the same defect.

## 6. Boundary of the theorem

For a general target set `A`, its pullback `T^{-1}(A)` is not a Gaussian
halfspace.  The Gaussian integration-by-parts field solving

\[
 \mathbb E[(\mathbf1_{T^{-1}A}-g)F]
 =\mathbb E\langle W,\nabla F\rangle                              \tag{6.1}
\]

is vector-valued and need not have a fixed direction.  Formula (4.1) becomes

\[
 D=\mathbb E[(HW)X^T+X(HW)^T].                                     \tag{6.2}
\]

The scalar column corresponding to the varying magnitude and direction of
`W` is the exact generalization of (4.5).  Neither ordinary Gaussian
Poincare nor the bound `Cov(X) preceq I` controls it by
`I(g)-|E HW|`: both leave a nonvanishing term in the halfspace equality
case.  Thus extending Theorem 1 requires a quantitative rigidity theorem for
the Gaussian divergence flux of a nearly extremal Lipschitz-set pair.  No
such theorem is assumed here.

This identifies a strictly narrower obstruction than arbitrary
``transverse dependence'': dependence itself is already allowed in Theorem
1; what remains is rotation/branching of the Gaussian transport flux of the
pullback set.

## 7. A second exact case: arbitrary sets under a linear contraction

There is another useful closure which is complementary to Theorem 1.  Let
`M` be a symmetric positive-semidefinite matrix with `||M||op<=1`, put
`X=MG`, and let `B` now be an arbitrary Borel set of central Gaussian mass
`g`.  Define

\[
 a=\mathbb E[(\mathbf1_B-g)G],\qquad
 D_G=\mathbb E[(\mathbf1_B-g)GG^T].                   \tag{7.1}
\]

Then `v=Ma` and `D=MD_GM`.  If `a ne0`, set `e=a/|a|`.  Since
`v=|a|Me`, its target direction `u` satisfies

\[
 PMe=0,
 \qquad PM=PMP_e,\quad P_e=I-ee^T.                  \tag{7.2}
\]

The exact Gaussian angular theorem gives

\[
 \|P_eD_G\|_{HS}^2\le C_\delta\{\mathcal I(g)-|a|\}. \tag{7.3}
\]

Consequently

\[
 \begin{aligned}
 \|PD\|_{HS}
 &=\|PMP_eD_GM\|_{HS}\le\|P_eD_G\|_{HS},\\
 \|PD\|_{HS}^2
 &\le C_\delta\{\mathcal I(g)-|a|\}
 \le C_\delta\{\mathcal I(g)-|Ma|\}.                \tag{7.4}
 \end{aligned}
\]

If `a=0` then `v=0`, so the active projection is not defined and the
near-equality regime is absent.  Scaling gives the same `L^3` formula as in
(1.4).  Thus the proposed pointwise theorem is proved for:

1. every set and every anisotropic Gaussian target; and
2. every nonlinear Brenier contraction when the pullback set is a Gaussian
   halfspace.

The genuinely open intersection is simultaneous nonlinearity of the map and
non-halfspace geometry of the pullback.

## 8. Exact integrated substitute for the remaining intersection

Although the last static intersection is not resolved above, ordinary
Gaussian localization gives an exact integrated-over-tilts statement in
precisely the required generality.  Let `mu_0` be `t_0`-strongly
log-concave, let

\[
 d\mu_s(x)\propto
 \exp\{c_s\cdot x-s|x|^2/2\}\,d\mu_0(x),
 \qquad dc_s=dW_s+a_sds,                              \tag{8.1}
\]

and fix a Borel set `A`.  With

\[
 \begin{gathered}
 g_s=\mu_s(A),\quad
 v_s=\mathbb E_s[(\mathbf1_A-g_s)(X-a_s)],\quad
 r_s=|v_s|,\\
 u_s=v_s/r_s,\quad P_s=I-u_su_s^T,\quad
 D_s=\mathbb E_s[(\mathbf1_A-g_s)(X-a_s)(X-a_s)^T],\\
 \lambda_s=t_0+s,qquad
 \Delta_s={\mathcal I(g_s)\over\sqrt{\lambda_s}}-r_s,
 \end{gathered}                                                       \tag{8.2}
\]

Ito calculus gives, up to bounded moment and `r_s>0` stops,

\[
 \boxed{
 d\Delta_s=dM_s-left\{
 {\|P_sD_s\|_{HS}^2\over2r_s}
 +r_su_s^T(\lambda_s^{-1}I-A_s)u_s
 +{\mathcal I(g_s)(1-\eta_s)^2\over2\lambda_s^{3/2}}
 \right\}ds.}                                                        \tag{8.3}
\]

where `A_s=Cov(mu_s)` and
`eta_s=sqrt(lambda_s)r_s/I(g_s)`.  Brascamp--Lieb gives
`A_s preceq lambda_s^{-1}I`, so all three
drifts are nonnegative.  Therefore, for every bounded stopping time `tau`,

\[
 \boxed{
 \mathbb E\int_0^\tau {\|P_sD_s\|_{HS}^2\over2r_s}\,ds
 \le\Delta_0.}                                                       \tag{8.4}
\]

The law of the natural parameter at a deterministic time is the Gaussian
channel mixture

\[
 c_s\ \stackrel d=\ sX+\sqrt s\,Z,
 \qquad X\sim\mu_0,\quad Z\sim\gamma_n\quad\hbox{independent}.       \tag{8.5}
\]

Thus (8.4) is literally an average over Gaussian natural tilts, not merely a
pathwise heuristic.  On a multiplicative curvature window, stopped while
`g_s` is central and `eta_s>=alpha>0`, it implies

\[
 \mathbb E\int {\|P_sD_s\|_{HS}^2\over r_s^2}\,ds
 \le C_{\delta,\alpha,\Lambda}(1-\eta_0).             \tag{8.6}
\]

This integrated identity has been independently checked from
`dg=v^TdW`, `dv=DdW-Avds`, and the Ito formula for `|v|`.  It avoids the
static flux-branching term in (6.2) by following the actual posterior path;
the angular quadratic variation is itself one of the dissipations.

## 9. Exact conditional splitting matrix from Prekopa

There is a static nonlinear splitting fact which is useful in analyzing the
remaining intersection.  Normalize the curvature to one and write a smooth
joint density in active/transverse coordinates as

\[
 p(y,z)={1\over Z}\exp\left{-{y^2+|z|^2\over2}-W(y,z)\right\},
 \qquad \nabla^2W\succeq0.                             \tag{9.1}
\]

Let `nu_y` be the conditional law of `Z` given `Y=y`, and put

\[
 m(y)=\mathbb E_{\nu_y}Z,qquad
 C(y)=\operatorname {Cov}_{\nu_y}(Z).                 \tag{9.2}
\]

Write the active marginal as

\[
 p_Y(y)={1\over Z_1}\exp\{-y^2/2-\overline W(y)\}.    \tag{9.3}
\]

**Lemma 9 (conditional splitting matrix).**  At every point of twice
differentiability,

\[
 \boxed{
 \begin{pmatrix}
  \overline W''(y)&-m'(y)^T\\
  -m'(y)&I-C(y)
 \end{pmatrix}\succeq0.}                              \tag{9.4}
\]

In particular, for every transverse vector `theta`,

\[
 \boxed{
 |\langle m'(y),\theta\rangle|^2
 \le\overline W''(y)\,
       \langle(I-C(y))\theta,\theta\rangle.}          \tag{9.5}
\]

**Proof.**  Define, up to an irrelevant additive constant,

\[
 F(y,b)=-\log\int_{\mathbb R^{n-1}}
 \exp\left\{-{|z-b|^2\over2}-W(y,z)\right\}dz.       \tag{9.6}
\]

The function

\[
 (y,z,b)\longmapsto {|z-b|^2\over2}+W(y,z)
\]

is convex.  Prekopa's theorem therefore says that `F` is convex jointly in
`(y,b)`.  At `b=0`, direct differentiation gives

\[
 F_{yy}=\overline W'',\qquad
 \nabla_bF=-m,qquad
 \nabla_{yb}^2F=-m',\qquad
 \nabla_{bb}^2F=I-C.                                  \tag{9.7}
\]

Thus (9.4) is exactly `nabla^2F(y,0) succeq0`, and (9.5) is its Schur
complement.  Approximation gives the distributional/Alexandrov version for
an arbitrary convex `W`.

The lemma is sharp.  For a quadratic interaction

\[
 W(y,z)={1\over2}(ay^2+2yb^Tz+z^TCz),
 \qquad \begin{pmatrix}a&b^T\\b&C\end{pmatrix}\succeq0,
\]

one has

\[
 \overline W''=a-b^T(I+C)^{-1}b,quad
 m'=-(I+C)^{-1}b,quad
 \operatorname {Cov}(Z\mid Y)= (I+C)^{-1},           \tag{9.8}
\]

and (9.4) reduces to the ordinary Schur-complement inequality for the
quadratic block matrix.

After curvature scaling, (9.4) quantifies a fact that was previously only
heuristic: variation of the conditional transverse mean cannot occur for
free when the active marginal is close to its Gaussian curvature equality
case.  The cost is the active marginal's excess convexity, paired with the
transverse covariance deficit.

The transverse--transverse part of `PD` also involves variation of
`C(y)`.  This is a third derivative `-partial_y F_{bb}` of the convex
function in (9.6).  Convexity of `F` proves (9.4) but, by itself, does not
give a pointwise bound on that third derivative by the Hessian slack at the
same point.  Heat-time integration controls precisely this term through the
Bernstein identity in Section 8; a static pointwise self-concordance estimate
for (9.6), strong enough in Hilbert--Schmidt norm, is the residual question.
