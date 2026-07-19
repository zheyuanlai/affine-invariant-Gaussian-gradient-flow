# Covariant Hessian rigidity and the phase-cell no-go

## 0. Result

There is a natural second-order tensor for the heat posterior field which
removes all one-dimensional profile amplitude and treats the Gaussian radial
ball as an exact zero-energy model.  Its exact local null fields are parallel
or concentric radial fields in dimensions at least three.  Nevertheless the
desired quantitative global conclusion is false under only

\[
                    C_P(q_s)\le K+s.                            \tag{0.1}
\]

The failure is not merely qualitative.  Balanced maximum cuts of Gaussian
products give a completely explicit sequence, and centered exponential
products give the same mandatory phase-cell test, for which

* the normalized covariant energy tends to zero;
* the Fisher effective rank tends to infinity;
* every individual flat phase has Fisher mass `1/n`;
* no region of any prescribed positive Fisher mass has approximately
  concurrent normals; and
* `C_P(q_s)<=4+s=4(1+alpha)` with `alpha=s/4` throughout.

Thus a small covariant Hessian plus high Fisher effective rank does **not**
force a positive-mass concurrent radial packet.  The missing input is a
quantitative lower bound on the mass/overlap of phase interfaces.  Such a
bound is not supplied by (0.1).  A possible theorem with an additional fixed
lower bound on `tr R` is not refuted: in the counterexample
`tr R=Theta(sqrt(alpha))` tends to zero.  This parameter loss is explicit.

The Gaussian halfspace, Gaussian median ball, inner square in the cube,
product-exponential maximum, and radial quartic perturbation are all computed
below.

---

## 1. Heat posterior notation and the Fisher law

Let `mu` be isotropic log-concave with `C_P(mu)=K`, let `S` be a Borel set,
and, for `s>0`, write

\[
 q_s=P_s\mu,\qquad
 g_s=\frac{P_s(\mathbf1_S\mu)}{P_s\mu},\qquad
 z=\Phi^{-1}(g_s),\qquad h=2\arcsin\sqrt{g_s}.                 \tag{1.1}
\]

Put

\[
 A(z)=\frac{\varphi(z)}{\sqrt{\Phi(z)(1-\Phi(z))}},\qquad
 F=\sqrt{s}\,\nabla h=\sqrt{s}\,A(z)\nabla z,                \tag{1.2}
\]

and

\[
 R=\mathbb E_{q_s}[FF^{\mathsf T}],\qquad
 \tau=\operatorname {tr}R.                                   \tag{1.3}
\]

When `tau>0`, define the Fisher probability and its direction covariance by

\[
 d\nu_F(y)=\frac{|F(y)|^2q_s(y)}{\tau}\,dy,\qquad
 u(y)=\frac{\nabla z(y)}{|\nabla z(y)|},\qquad
 M_F=\int uu^{\mathsf T}d\nu_F=\frac R{\operatorname {tr}R}.  \tag{1.4}
\]

The values assigned to `u` on `{nabla z=0}` are irrelevant.  The operator
effective rank is

\[
 r_{\rm eff}(R)=\frac{\operatorname {tr}R}{\|R\|_{\rm op}}
                =\frac1{\|M_F\|_{\rm op}}.                   \tag{1.5}
\]

The convolution variance decomposition, and no KLS estimate, gives

\[
 C_P(q_s)\le K+s=K(1+\alpha),\qquad \alpha=\frac sK.          \tag{1.6}
\]

---

## 2. The reparametrization-covariant Hessian

At a point where `nabla z` is nonzero, set

\[
 H=\nabla^2z,\qquad P=I-uu^{\mathsf T},\qquad a=|\nabla z|.   \tag{2.1}
\]

The derivative of the direction field is

\[
             \nabla u=\frac{PH}{a}.                            \tag{2.2}
\]

Decompose it into the derivative along the normal and the shape operator of
the level set:

\[
 b=\frac{PHu}{a},\qquad
 \mathcal S=\frac{PHP}{a},\qquad
 \kappa=\frac{\operatorname {tr}\mathcal S}{n-1}.             \tag{2.3}
\]

For `n>=2`, define the covariant Hessian defect

\[
\boxed{
 |\mathcal C[z]|^2
 =|b|^2+\|\mathcal S-\kappa P\|_{\rm HS}^2.}                  \tag{2.4}
\]

This definition removes exactly the two pieces which should not be charged:

* `H` in the `u-u` entry changes only the speed of a one-dimensional
  profile;
* the scalar transverse curvature `kappa P` is the curvature of a sphere.

It is genuinely covariant under profile changes.  If `w=f(z)` and `f'` does
not vanish, then

\[
 \nabla^2w=f'\nabla^2z+f''\nabla z\nabla z^{\mathsf T}.        \tag{2.5}
\]

The second term is annihilated by `P`, while division by `|nabla w|`
cancels `|f'|`.  Hence

\[
                   |\mathcal C[f\circ z]|=|\mathcal C[z]|.    \tag{2.6}
\]

In particular (2.4) removes the nonlinear amplitude of every flat
one-dimensional Gaussian profile, not only an affine `z` profile.

The dimensionless normalized energy at heat scale `s` is

\[
\boxed{
 \mathscr E_{\rm cov}(s)
 =s\int|\mathcal C[z](y)|^2d\nu_F(y).}                        \tag{2.7}
\]

Equivalently,

\[
 \mathscr E_{\rm cov}(s)
 =\frac{s^2}{\tau}\int q_sA(z)^2
 \left\{|PHu|^2+\left\|PHP-
       \frac{\operatorname {tr}(PHP)}{n-1}P\right\|_{\rm HS}^2
 \right\}dy.                                                 \tag{2.8}
\]

The full Hessian in the heat Bernstein identity dominates the braces in
(2.8), but the weights `q A(z)^2` and `q varphi(z)` are not uniformly
comparable in the posterior tails.  The phase-cell counterexample below lies
in the central boundary layer and works for either weight; it is therefore
not an artifact of this distinction.

---

## 3. Exact local affine-or-radial classification

The local zero set of (2.4) has the intended geometry, subject to an
important two-dimensional exception.

**Lemma 3.1 (exact local null classification).**  Let `n>=3`, let `U` be a
connected open subset of `R^n`, and let `z in C^3(U)` have no critical
point.  Suppose `mathcal C[z]=0` on `U`.

On each connected component on which `kappa=0`, the unit field `u` is
constant and the level sets of `z` are parallel hyperplanes.  On each
connected component on which `kappa` is nonzero, there is a single point
`z_0` such that

\[
                  u(y)=\epsilon\frac{y-z_0}{|y-z_0|},
 \qquad \epsilon\in\{-1,1\},                                 \tag{3.1}
\]

and `z` is a one-dimensional reparametrization of `|y-z_0|`.

**Proof.**  Equations (2.2)--(2.4) give

\[
 D_u u=0,\qquad D_vu=\kappa v\quad(v\perp u).                 \tag{3.2}
\]

The Euclidean Codazzi identity, applied to two independent tangent
directions, gives

\[
                         P\nabla\kappa=0.                      \tag{3.3}
\]

Transport a tangent vector `v` so that `[u,v]=0`.  Since the connection is
torsion free, `D_uv=D_vu=kappa v`.  Differentiating (3.2) along `u` gives

\[
                         D_u\kappa=-\kappa^2.                  \tag{3.4}
\]

If `kappa=0`, (3.2)--(3.3) give `Du=0`.  If `kappa` is nonzero, define

\[
                         c(y)=y-\frac{u(y)}{\kappa(y)}.         \tag{3.5}
\]

Equations (3.2)--(3.4) show `D_vc=0` for tangent `v` and `D_uc=0`.
Thus `c=z_0` is constant.  Formula (3.1) follows, and the fact that
`nabla z` is parallel to `y-z_0` makes `z` a radial reparametrization. QED.

The restriction `n>=3` is essential for this second-order tensor.  In two
dimensions the traceless shape term is identically zero.  The signed
distance to any `C^3` plane curve satisfies `D_u u=0`, and hence
`mathcal C=0`, even when the normals are neither parallel nor concurrent.
An ellipse already defeats a two-dimensional pointwise classification.
Adding a derivative of `kappa` would repair this local defect, but it would
be a third-order quantity not controlled by the Hessian term in the
Bernstein identity.

Lemma 3.1 is local.  It does not glue components separated by a set on which
`|F|` is small.  The next two sections show that this is a sharp failure.

---

## 4. Polyhedral boundary-layer lemma

The following standard heat-content asymptotic makes the phase-cell
obstruction formal.

**Lemma 4.1 (flat-facet Fisher limit).**  Let a log-concave measure have a
`C^2`, positive density in a neighborhood of the relative interiors of the
finitely many facets of a polyhedral set `S`.  Assume its weighted facet
area and weighted codimension-two skeleton area are finite.  Let `n_j` be
the oriented unit normal of facet `F_j`.  Then, as `s downarrow0`,

\[
\begin{aligned}
 R_s
 &=c_*\sqrt{s}\sum_j
    \left(\int_{F_j}\rho_\mu\,d\mathcal H^{n-1}\right)
       n_jn_j^{\mathsf T}+o(\sqrt{s}),\\
 \operatorname {tr}R_s
 &=c_*\sqrt{s}\,\operatorname {Per}_\mu(S)+o(\sqrt{s}),       \tag{4.1}
\end{aligned}
\]

where

\[
                     c_*=\int_{\mathbb R}A(t)^2dt\in(0,\infty).
                                                                    \tag{4.2}
\]

Moreover, `nu_F` converges weakly to normalized weighted surface measure
with direction `u=n_j` on `F_j`, and

\[
                         \mathscr E_{\rm cov}(s)=O(\sqrt{s}). \tag{4.3}
\]

The constant in (4.3) may depend on the fixed polyhedron.  The same result
holds when the ambient log-concave density has a polyhedral support, provided
the intersections with the support boundary are included in the skeleton.

**Proof.**  At a point a fixed positive distance from the skeleton, write
`y=x+sqrt(s)t n_j` in normal coordinates.  Gaussian convolution and one
Taylor expansion of the log density give

\[
 g_s(y)=\Phi(-t)+O(\sqrt{s}),\qquad
 z_s(y)=-t+O(\sqrt{s}),\qquad
 \sqrt{s}\nabla z_s=-n_j+O(\sqrt{s}).                         \tag{4.4}
\]

The errors and their first two scaled derivatives are uniform on compact
facet pieces.  Substitution into (1.2)--(1.4), followed by dominated
convergence in `t`, proves (4.1) and the weak limit.

On those facet pieces `mathcal C=O(sqrt(s))`.  In a `C sqrt(s)`
neighborhood of the skeleton, scaled cone coordinates give

\[
 |\nabla z|=O(s^{-1/2}),\qquad |\nabla^2z|=O(s^{-1}),
 \qquad |\mathcal C|=O(s^{-1/2}).                             \tag{4.5}
\]

The skeleton neighborhood has Fisher mass `O(sqrt(s))` relative to the
facet layer.  Hence `s|mathcal C|^2=O(1)` there and its contribution to
(2.7) is `O(sqrt(s))`.  Gaussian tail truncation handles the complement of
the fixed scaled neighborhood.  This proves (4.3). QED.

The lemma may equally be proved by the explicit Gaussian cone integrals.  It
uses no isoperimetric minimality and no Poincare inequality.

---

## 5. First no-go: the inner square in the isotropic cube

Let `mu` be uniform on the isotropic square

\[
                  [-\sqrt3,\sqrt3]^2,
\]

and let

\[
                  S=[-\sqrt{3/2},\sqrt{3/2}]^2.                \tag{5.1}
\]

Then `mu(S)=1/2` and

\[
                         K=C_P(\mu)=\frac{12}{\pi^2}.          \tag{5.2}
\]

The inner boundary is a fixed positive distance from the support boundary,
so Lemma 4.1 applies without a contact correction.  By the four symmetries,

\[
 \frac{R_s}{\operatorname {tr}R_s}\longrightarrow\frac12I_2,
 \qquad r_{\rm eff}(R_s)\longrightarrow2,                     \tag{5.3}
\]

while

\[
 \operatorname {tr}R_s=\Theta(\sqrt{s})=\Theta(\sqrt\alpha),
 \qquad
 \mathscr E_{\rm cov}(s)=O(\sqrt{s})=O(\sqrt\alpha),         \tag{5.4}
\]

where `s=alpha K`.  Of course

\[
                         C_P(q_s)\le K(1+\alpha).              \tag{5.5}
\]

The limiting Fisher law is uniform on the four faces with its corresponding
normal line.  Define the projective concurrence error

\[
 \Gamma(\nu,u)=\inf_{z_0\in\mathbb R^2}
  \int\left[1-\left\langle u(y),
       \frac{y-z_0}{|y-z_0|}\right\rangle^2\right]d\nu(y).    \tag{5.6}
\]

There is a numerical `c_square>0` such that

\[
                         \liminf_{s\downarrow0}\Gamma(\nu_F,u)
                         \ge c_\square.                        \tag{5.7}
\]

Indeed compactify the possible centers by adding directions at infinity.
The limiting functional is continuous on that compactification.  It cannot
vanish at a finite center, because the normal lines along even one open face
of a square do not pass through one point.  At infinity the candidate radial
line becomes one constant line, whose error is `1/2` against
`M_F=I_2/2`.  Compactness proves (5.7).

Thus even the whole Fisher field need not be affine or concurrent when the
normalized covariant energy is small.  This example alone does not rule out
selecting one face: each face has limiting Fisher mass `1/4`.  The next
example removes that escape.

---

## 5A. A fully quantitative fixed-mass no-go: the Gaussian maximum

Before treating the required exponential model, it is useful to record a
version whose small-ball step is a one-line calculation.  Let

\[
 \mu_n=\gamma_n,\qquad
 S_n^G=\{\max_iX_i\ge L_n^G\},\qquad
 \Phi(L_n^G)^n=\frac12.                                      \tag{5A.1}
\]

Then `K=1`, `L_n^G=sqrt(2log n)+o(sqrt(log n))`, and the cut is balanced.
The limiting Fisher boundary law is obtained by choosing `I` uniformly,
putting `u=e_I` and `Y_I=L_n^G`, and taking the other coordinates to be
independent standard Gaussians conditioned to lie below `L_n^G`.  Call this
law `beta_n^G`.

**Lemma 5A.1 (uniform Gaussian phase nonconcurrence).**  For every `c>0`
there is `eta_G(c)>0` such that, for all sufficiently large `n`, every
`A` with `beta_n^G(A)>=c` and every `z_0 in R^n` obey

\[
 \frac1{\beta_n^G(A)}\int_A
 \left[1-\left\langle e_I,
       \frac{Y-z_0}{|Y-z_0|}\right\rangle^2\right]d\beta_n^G
 \ge\eta_G(c).                                                \tag{5A.2}
\]

**Proof.**  We first record the uniform noncentral Gaussian small-ball
bound.  If `G~N(0,I_m)`, then for numerical `c_0,c_1>0` and every
`z in R^m`,

\[
 \mathbb P\left\{|G-z|^2\le c_0(m+|z|^2)\right\}
 \le e^{-c_1(m+|z|^2)}.                                      \tag{5A.3}
\]

Indeed, for `t>0`,

\[
 \mathbb E e^{-t|G-z|^2}
 =(1+2t)^{-m/2}
   \exp\left(-\frac{t}{1+2t}|z|^2\right).                    \tag{5A.4}
\]

Markov's inequality with `t=1/2` and a sufficiently small `c_0` proves
(5A.3).  Conditioning all coordinates to be below `L_n^G` multiplies a
probability by at most

\[
 \Phi(L_n^G)^{-(n-1)}\le3.                                   \tag{5A.5}
\]

Let `m_n` be the mean of the truncated Gaussian coordinate, put

\[
 d_n=L_n^G-m_n=O(\sqrt{\log n}),\qquad
 b_j=m_n-z_{0,j},\qquad B^2=\sum_jb_j^2.                      \tag{5A.6}
\]

On facet `i`, squared radial correlation at least `1-delta` requires

\[
 \sum_{j\ne i}(Y_j-z_{0,j})^2
 \le C_\delta(d_n+b_i)^2.                                    \tag{5A.7}
\]

By (5A.3)--(5A.5), its conditional probability is exponentially small
unless

\[
 (d_n+b_i)^2\ge c_\delta(n+B^2-b_i^2).                       \tag{5A.8}
\]

Since `d_n^2=o(n)` and
`(d_n+b_i)^2<=2d_n^2+2b_i^2`, every index satisfying (5A.8) also satisfies

\[
                         b_i^2\ge c'_\delta(n+B^2).           \tag{5A.9}
\]

There are at most `C_delta` such indices, because their squared `b_i`
sum to at most `B^2`.  Averaging over the uniform active index gives

\[
 \sup_{z_0}\beta_n^G\left\{
 \left\langle e_I,\frac{Y-z_0}{|Y-z_0|}\right\rangle^2
       \ge1-\delta\right\}
 \le\frac{C_\delta}{n}+3e^{-c_\delta n}=o(1).                \tag{5A.10}
\]

For `beta_n^G(A)>=c`, at least half of `A` is outside the event in
(5A.10), for large `n`.  Taking `eta_G(c)=delta/2` proves (5A.2). QED.

Lemma 4.1 and the permutation symmetry give, at fixed `n` as `s downarrow0`,

\[
 \frac{R_{n,s}}{\operatorname {tr}R_{n,s}}\longrightarrow\frac1nI,
 \qquad \mathscr E_{\rm cov}(n,s)\longrightarrow0.           \tag{5A.11}
\]

Choose a diagonal `s_n downarrow0` so that, outside Fisher mass `o(1)`, the
field is coupled to `beta_n^G` with spatial displacement `o(1)` and angular
error `o(1)`, the energy is at most `1/n`, and the effective rank is at least
`n/2`.  The estimates (5A.7)--(5A.10) have fixed slack and therefore remain
valid under this coupling; the singular case `Y=z_0` has zero boundary mass.
Since `q_s=N(0,(1+s)I)`, this sequence satisfies the exact allowed spectral
statement

\[
                         C_P(q_{s_n})=1+s_n=1+\alpha_n,
 \qquad \alpha_n=s_n.                                        \tag{5A.12}
\]

Lemma 5A.1 then proves the fixed-mass no-go with no unproved expansion or
small-ball input.  Its Fisher trace is

\[
 \operatorname {tr}R_{n,s}
 =\Theta\left(\sqrt{s}\,n\varphi(L_n^G)\Phi(L_n^G)^{n-1}\right)
 =\Theta\left(\sqrt{s\log n}\right),                         \tag{5A.13}
\]

so it too leaves open a theorem with a fixed trace threshold.

---

## 6. Fixed-mass no-go: the product-exponential maximum

Let `Z_1,...,Z_n` be independent rate-one exponentials and put

\[
 X_i=Z_i-1.
\]

This product measure is isotropic, log-concave, and has

\[
                         C_P(\mu_n)=4.                         \tag{6.1}
\]

Choose `L_n` by

\[
                 (1-e^{-L_n})^n=\frac12,
 \qquad L_n=\log n+O(1),                                     \tag{6.2}
\]

and set

\[
                 S_n=\{\max_iZ_i\ge L_n\}.                   \tag{6.3}
\]

This is a balanced polyhedral cut.  Its `n` principal facets have equal
weighted area, and the total weighted perimeter is bounded above and below
by numerical constants.  Lemma 4.1 gives, for every fixed `n`,

\[
\begin{aligned}
 \operatorname {tr}R_{n,s}&=\Theta(\sqrt{s}),\\
 \frac{R_{n,s}}{\operatorname {tr}R_{n,s}}&\longrightarrow
                         \frac1nI_n,\\
 r_{\rm eff}(R_{n,s})&\longrightarrow n,\\
 \mathscr E_{\rm cov}(n,s)&\longrightarrow0.                 \tag{6.4}
\end{aligned}
\]

All constants in the first line are numerical; the rate in the last line
may depend polynomially on `n` because of support contacts.  We are free to
choose a diagonal sequence `s_n downarrow0` after fixing `n`.

The limiting Fisher boundary law has an exact description.  Choose
`I` uniformly from `{1,...,n}`, put `u=e_I`, set `Y_I=L_n`, and, conditional
on `I`, let the other coordinates be independent exponentials conditioned
to be below `L_n`.  Denote this law by `beta_n`.

**Lemma 6.1 (no positive-mass concurrent packet).**  For every `c>0` there
is `eta(c)>0` such that, for all sufficiently large `n`, every Borel set
`A` with `beta_n(A)>=c` and every `z_0 in R^n` satisfy

\[
 \frac1{\beta_n(A)}\int_A
 \left[1-\left\langle e_I,
       \frac{Y-z_0}{|Y-z_0|}\right\rangle^2\right]d\beta_n
 \ge\eta(c).                                                  \tag{6.5}
\]

The same assertion includes the affine limit `|z_0|=infinity`.

**Proof.**  Let `xi_n` be an exponential conditioned below `L_n`, let
`m_n=E xi_n`, and put

\[
 d_n=L_n-m_n=\log n+O(1),\qquad b_j=m_n-z_{0,j},qquad
 B^2=\sum_jb_j^2.                                             \tag{6.6}
\]

For a point on facet `i`, the squared radial correlation is

\[
 C_i=\frac{(d_n+b_i)^2}
 {(d_n+b_i)^2+\sum_{j\ne i}(\xi_{n,j}-z_{0,j})^2}.            \tag{6.7}
\]

An elementary product small-ball estimate, uniform in all real `b_j`, says
that outside a set of probability `o_n(1)`,

\[
 \sum_{j\ne i}(\xi_{n,j}-z_{0,j})^2
 \ge c_0\left(n+\sum_{j\ne i}b_j^2\right),                   \tag{6.8}
\]

except possibly for `O(1)` indices carrying a fixed fraction of `B^2`.
To verify (6.8), split the indices into `|b_j|<=4` and `|b_j|>4`.
For the first group, a fixed-probability deviation of `xi_n` and Chernoff's
inequality give a contribution proportional to its cardinality.  For the
second group, `|xi_n-z_{0,j}|>=|b_j|/2` except on an interval whose
exponential mass is summable; dyadically grouping the `|b_j|` gives the
weighted assertion.  A coordinate carrying a fixed fraction of `B^2` can
create only one exceptional active facet.

Fix `delta in (0,1/4)`.  Equations (6.7)--(6.8) show that `C_i>=1-delta`
can have nonnegligible conditional probability only if

\[
 (d_n+b_i)^2\ge c_\delta
       \left(n+B^2-b_i^2\right).                               \tag{6.9}
\]

Since `(d_n+b_i)^2<=2d_n^2+2b_i^2` and `d_n^2=o(n)`, (6.9) implies

\[
                         b_i^2\ge c'_\delta(n+B^2).            \tag{6.10}
\]

Summing (6.10) over all such indices shows that their number is at most
`C_delta`, independently of `n`.  Since `I` is uniform,

\[
             \sup_{z_0}\beta_n\{C_I\ge1-\delta\}=o_n(1).      \tag{6.11}
\]

If `beta_n(A)>=c`, at least half of `A` lies outside this good-correlation
set for large `n`.  Taking `eta(c)=delta/2` proves (6.5). QED.

Choose `s_n>0` so small that, off Fisher mass `o(1)`, the heat field admits
the analogous facet coupling to `beta_n` with `o(1)` spatial and angular
error, and

\[
 \mathscr E_{\rm cov}(n,s_n)\le\frac1n,
 \qquad r_{\rm eff}(R_{n,s_n})\ge\frac n2.                    \tag{6.12}
\]

Set

\[
                         \alpha_n=\frac{s_n}{4}.               \tag{6.13}
\]

Then

\[
 C_P(q_{s_n})\le4+s_n=4(1+\alpha_n),                          \tag{6.14}
\]

while (6.5) persists for every Fisher region of mass at least `c`.  This is
the promised formal no-go counterexample.

Its only degenerating parameter is visible:

\[
                 \operatorname {tr}R_{n,s_n}
                 =\Theta(\sqrt{s_n})
                 =\Theta(\sqrt{\alpha_n}).                    \tag{6.15}
\]

Therefore an inverse theorem which additionally assumes
`tr R>=r_0>0` may have a modulus depending on `r_0`; the present example
does not refute it.  No conclusion depending only on effective rank and
normalized covariant energy can survive.

---

## 7. The two exact positive models

### 7.1 Gaussian halfspace

Let `mu=N(0,I)`, so `K=1`, and let `S={x_1<=a}`.  With
`beta=1+s`, the conditional law is

\[
 X\mid Y=y\sim N\left(\frac y\beta,\frac s\beta I\right).
\]

Consequently

\[
 z(y)=a\sqrt{\frac\beta s}-\frac{y_1}{\sqrt{s\beta}}.         \tag{7.1}
\]

Thus `H=0`, `mathcal C=0`, and `R` has rank one.  This is the affine
branch.  The only spectral input remains

\[
                         C_P(q_s)=1+s=1+\alpha.                \tag{7.2}
\]

### 7.2 Gaussian median ball

Let `mu=N(0,I)` and let `S` be a centered median ball.  Both `q_s` and the
posterior probability are radial, so `z(y)=f_s(r)` with `r=|y|`.  Wherever
`f_s'` is nonzero,

\[
 u=\epsilon\frac y{|y|},\qquad
 H=f_s''uu^{\mathsf T}+\frac{f_s'}rP.                         \tag{7.3}
\]

It follows exactly, at every `n` and `s`, that

\[
                  b=0,\qquad \mathcal S=\frac\epsilon rP,
 \qquad \mathcal C=0.                                        \tag{7.4}
\]

Rotational symmetry gives `R=(tr R/n)I`, so the effective rank is `n`, and
the normal lines are concurrent at the origin.  This is precisely why the
scalar transverse curvature was removed in (2.4).  The explicit trace and
wedge saturation calculation is given in
`gradient_multiplicity_rigidity.md`; none of it is evidence for a product
factor.

---

## 8. Radial quartic perturbation

Let

\[
 d\mu_\epsilon(x)=Z_\epsilon^{-1}
 \exp\left(-\frac{|x|^2}{2}-\frac\epsilon{4n}|x|^4\right)dx,  \tag{8.1}
\]

followed by the scalar dilation making the measure isotropic, and let `S`
be its centered median ball.  This is smooth, full-dimensional, strongly
log-concave, and has no nontrivial affine product decomposition when
`epsilon>0`.

Gaussian convolution preserves rotational symmetry.  Hence `q_s`, `g_s`,
and `z_s` are radial for every `s>0`.  Equations (7.3)--(7.4) apply without
change:

\[
 \mathcal C[z_s]=0,\qquad
 R_s=\frac{\operatorname {tr}R_s}{n}I,qquad
 r_{\rm eff}(R_s)=n.                                         \tag{8.2}
\]

The field is exactly concurrent at the origin although the measure has no
Gaussian or product factor.  Thus the radial conclusion is the correct
positive branch; replacing it by a splitting conclusion would be false.
Writing `K_epsilon=C_P(mu_epsilon)` and `s=alpha K_epsilon`, the only
available general estimate is still

\[
 C_P(q_s)\le K_\epsilon+s=K_\epsilon(1+\alpha).               \tag{8.3}
\]

No quotient Poincare inequality has been used.

---

## 9. Sharp scope of the no-go

The proposed implication

\[
 \mathscr E_{\rm cov}\ll1
 \quad+\quad r_{\rm eff}(R)\gg1
 \quad\Longrightarrow\quad
 u(y)\approx\frac{P(y-z_0)}{|P(y-z_0)|}
 \text{ on fixed Fisher mass}                                \tag{9.1}
\]

is false even if `P` is the projection onto `ran R`: in Section 6 that
projection is the identity.  The counterexample also rules out the affine
limit, because each affine phase has mass `1/n`.

Both fixed-mass counterexamples use `alpha_n downarrow0`.  They do not
refute a statement at a fixed lower scale ratio `alpha>=alpha_0>0`, nor a
statement whose required energy threshold is an explicit function of
`alpha` and `tr R`.  The inner-square calculation identifies the natural
loss: both its trace and its normalized interface error are of order
`sqrt(alpha)`.

There are three precise ways to strengthen the hypotheses which are not
disposed of here.

1. Impose `tr R>=r_0>0` and retain the full dependence on `r_0`.  The
   phase-cell construction pays `tr R=Theta(sqrt(alpha))`.
2. Charge the codimension-two interface by a scale-independent junction
   measure.  The bulk Hessian energy weights its `sqrt(s)` neighborhood and
   therefore loses `sqrt(s)`.
3. Assume that a positive fraction of the Fisher law lies in one connected
   overlap component.  Deriving that assertion from log-concavity and
   (0.1), however, is an expansion problem for a reweighted law; ordinary
   Poincare for `q_s` does not provide it.

If "fixed mass" in (9.1) means `q_s`-mass rather than Fisher mass, an
additional trace lower bound is indispensable even to state the conclusion:
in every small-heat boundary-layer example the entire Fisher field occupies
only `Theta(sqrt(alpha))` ordinary mass.

The formal output is therefore a local exact lemma (Lemma 3.1) together with
an explicit global no-go sequence (6.12)--(6.15).  A valid affine-or-radial
theorem must add an interface/overlap hypothesis or a fixed Fisher-trace
threshold; high effective rank and small covariant Hessian alone are
insufficient.
