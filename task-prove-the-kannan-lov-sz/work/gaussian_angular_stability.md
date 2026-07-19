# Gaussian centroid deficit controls transverse angular noise

## 1. The theorem

Let `X` be standard Gaussian in `R^n`, let `S` be Borel, and put

\[
 g=\gamma_n(S),\qquad
 v=\mathbb E[(\mathbf1_S-g)X].                         \tag{1.1}
\]

Assume `v ne 0`, put `u=v/|v|`, `P=I-u u^T`, and define

\[
 D=\mathbb E[(\mathbf1_S-g)XX^T].                      \tag{1.2}
\]

Let `I(g)=varphi(Phi^{-1}(g))` be the Gaussian isoperimetric profile.

**Theorem 1 (dimension-free angular stability).**  For every
`delta in (0,1/2)` there is a finite `C_delta` such that, whenever

\[
 \delta\le g\le1-\delta,
\]

one has

\[
 \boxed{
 \|PD\|_{HS}^2\le C_\delta\{\mathcal I(g)-|v|\}.}     \tag{1.3}
\]

The right side is nonnegative by the sharp Gaussian centroid inequality.
It vanishes for every Gaussian halfspace.  Both the projection on the left
and the subtraction on the right are essential: a noncentral halfspace has
a nonzero `u-u` entry in `D`, while rotating a halfspace changes its
unprojected first and second moments without creating a centroid deficit.

After scaling, if `X` is Gaussian with covariance `t^{-1}I`, then

\[
 \boxed{
 \|PD\|_{HS}^2
 \le {C_\delta\over t^{3/2}}
 \left\{{\mathcal I(g)\over\sqrt t}-|v|\right\}.}      \tag{1.4}
\]

Indeed `sqrt(t)X` is standard Gaussian, its centroid is `sqrt(t)v`, and
its second correlated moment is `tD`.

## 2. Slice notation and the exact flip cost

Rotate so that `u=e_1`, and write `X=(U,Z)` with
`U~gamma_1`, `Z~gamma_{n-1}` independent.  Let

\[
 c=\Phi^{-1}(1-g),\qquad H=\mathbf1_{\{U\ge c\}},
 \qquad \sigma=H-\mathbf1_S.                          \tag{2.1}
\]

For each transverse point `z`, define

\[
 \begin{aligned}
 m_0(z)&=\int\sigma(r,z)\varphi(r)dr,\\
 m_1(z)&=\int(r-c)\sigma(r,z)\varphi(r)dr.
 \end{aligned}                                        \tag{2.2}
\]

The sign pattern in (2.1) gives

\[
 (r-c)\sigma(r,z)\ge0,
 \qquad m_1(z)\ge0.                                   \tag{2.3}
\]

The threshold and `S` have the same mass, and `u` is the direction of the
centroid of `S`.  Hence

\[
 E m_0(Z)=0,qquad E[Z m_0(Z)]=0.                      \tag{2.4}
\]

Moreover,

\[
 \begin{aligned}
 E m_1(Z)
 &=E[(U-c)(H-1_S)]\\
 &=E[UH]-E[U1_S]
 =\mathcal I(g)-|v|.                                  \tag{2.5}
 \end{aligned}
\]

The `c` term disappears by equality of the masses.

We need one elementary one-dimensional estimate.

**Lemma 2 (net flip mass costs squared distance).**  There is a numerical
`c_0>0` such that, for every `z`,

\[
 m_1(z)\ge c_0m_0(z)^2.                               \tag{2.6}
\]

**Proof.**  Put `a(r)=|sigma(r,z)|`.  Then `0<=a<=1` and

\[
 \int a\,d\gamma_1\ge\left|\int\sigma\,d\gamma_1\right|
 =|m_0(z)|.                                           \tag{2.7}
\]

Also, by (2.3),

\[
 m_1(z)=\int|r-c|a(r)d\gamma_1(r).                    \tag{2.8}
\]

Since the standard Gaussian density is at most
`M=(2pi)^(-1/2)`, the interval `[c-s,c+s]` has Gaussian mass at most
`2Ms`.  If `A=int a dgamma`, the layer-cake/rearrangement minimum of
(2.8), among `0<=a<=1` of mass `A`, is bounded below by the cost of
placing density `M` on the two-sided interval of total length `A/M` about
`c`.  Thus

\[
 \int|r-c|a(r)d\gamma_1(r)\ge {A^2\over4M}.           \tag{2.9}
\]

Equations (2.7)--(2.9) prove (2.6), with `c_0=1/(4M)`.
`square`

For later use, centrality of `g` also gives

\[
 0\le m_1(z)\le\int|r-c|\varphi(r)dr\le C_\delta,     \tag{2.10}
\]

and therefore

\[
 E m_0(Z)^2\le C E m_1(Z),
 \qquad E m_1(Z)^2\le C_\delta E m_1(Z).              \tag{2.11}
\]

## 3. First- and second-chaos Bessel bounds

The halfspace correlated second moment

\[
 D_H=E[(H-g)XX^T]
\]

has no row in `u^perp`: independence and centering give `P D_H=0`.
Consequently

\[
 PD=-P E[\sigma XX^T].                                \tag{3.1}
\]

Its transverse-transverse block is

\[
 B=E[m_0(Z)ZZ^T]
  =E[m_0(Z)(ZZ^T-I)],                                  \tag{3.2}
\]

where the second equality uses (2.4).  Gaussian second-chaos orthogonality
gives

\[
 \boxed{\|B\|_{HS}^2\le2E m_0(Z)^2.}                  \tag{3.3}
\]

Indeed, for every symmetric matrix `A` with `||A||HS=1`,

\[
 E\{Z^TAZ-\operatorname {tr}A\}^2=2,
\]

and duality followed by Cauchy--Schwarz proves (3.3).

The transverse-`u` column is

\[
 \begin{aligned}
 b&=E[Z\,U\sigma(U,Z)]\\
  &=E[Z\{c m_0(Z)+m_1(Z)\}]
  =E[Z m_1(Z)],                                        \tag{3.4}
 \end{aligned}
\]

where (2.4) removes the first term.  Gaussian first-chaos orthogonality
and (2.10) give

\[
 |b|^2\le E m_1(Z)^2\le C_\delta E m_1(Z).            \tag{3.5}
\]

The two blocks in (3.2) and (3.4) are orthogonal in Hilbert--Schmidt norm.
Using (2.11), (2.5), (3.3), and (3.5),

\[
 \begin{aligned}
 \|PD\|_{HS}^2
 &=\|B\|_{HS}^2+|b|^2\\
 &\le C_\delta E m_1(Z)
 =C_\delta\{\mathcal I(g)-|v|\}.
 \end{aligned}
\]

This proves Theorem 1.

## 4. Relation to posterior direction dynamics

For a natural exponential tilt of a probability, differentiation with
respect to the tilt parameter gives

\[
 \nabla_c g=v,\qquad \nabla_c v=D.                    \tag{4.1}
\]

Thus, wherever `v ne0`,

\[
 \nabla_c\left({v\over|v|}\right)={PD\over|v|}.       \tag{4.2}
\]

For a Gaussian posterior of covariance `t^{-1}I`, (1.4) therefore gives

\[
 \left\|\nabla_c{v\over|v|}\right\|_{HS}^2
 \le {C_\delta\over t}
 {\mathcal I(g)-\sqrt t|v|\over t|v|^2}.              \tag{4.3}
\]

On central near-equality states, where `sqrt(t)|v|` is a fixed fraction of
`I(g)`, this is, up to a numerical factor,

\[
 {C_\delta\over t}\left(1-
 {t|v|^2\over\mathcal I(g)^2}\right).                 \tag{4.4}
\]

The same factor is the drift defect in Gaussian-profile localization.
This proves the desired scalar-to-angular implication for an *exact
Gaussian posterior*.

It does not yet prove it for a general `t`-strongly log-concave posterior.
The exact centroid decomposition only controls the active one-dimensional
Caffarelli map.  Transferring (1.3) requires quantitative transverse
splitting of the full contraction map; that step remains unproved and is the
current gluing obstruction.
