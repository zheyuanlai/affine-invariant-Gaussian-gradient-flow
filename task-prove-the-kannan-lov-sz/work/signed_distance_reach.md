# Signed distance, reach, and the two branches that can be closed

## 0. Verdict

There is a clean dimension-free theorem for two of the three possible
geometries of a long balanced signed-distance packet.

* A **coherent** packet is excluded directly by the covariance matrix.
* A **concurrent-radial** packet is excluded by the translated thin-shell
  estimate.
* Any packet not covered by those conclusions has effective normal rank
  of order at least the square of its length.  Its two-sided core has a
  globally Lipschitz Gauss map and its smooth Hilbert--Schmidt curvature is
  of order at most the reciprocal length.  The missing variation must occur
  at separated components, ridges, focal endpoints, or singular junctions.

The last, high-rank branch is not closed here.  Closing it for the ray packet
of a true T3 extremizer would prove the desired T3 bound and hence KLS; it
cannot be inserted as an auxiliary tube lemma.  Section 5 gives an exact
product-exponential counterexample to the tempting claim that high normal
rank by itself forces a fixed amount of ridge collision.  Thus the result of
this note is a rigorous branchwise exclusion and a precise residual case, not
a purported proof of KLS.

## 1. The balanced core and its exact reach

Let \(\mu\) be an isotropic log-concave probability on \(\mathbb R^n\), and
let \(f\) be a 1-Lipschitz function with a transport-ray disintegration

\[
 \mu=\int \nu_y\,d\eta(y),\qquad
 x=z_y+tN_y,\qquad f(z_y+tN_y)=t.                         \tag{1.1}
\]

Here \(|N_y|=1\), \(f(z_y)=0\), and each \(\nu_y\) is regarded as a
probability law for the coordinate \(T=t\).  For an extremizer of the
mean-centered first-moment problem, this is the nonbranching disintegration
from `t3_extremal_report.md`; the measures \(\nu_y\) are one-dimensional
log-concave and put the same sign mass on almost every ray.

Fix \(r>0\), \(0<\beta\leq 1/2\), and a measurable quotient packet \(G\)
of mass

\[
 \eta(G)=:\alpha_G\geq\alpha>0.                          \tag{1.2}
\]

The only quantitative balance hypothesis used below is

\[
 \nu_y\{T\geq r\}\geq\beta,
 \qquad
 \nu_y\{T\leq-r\}\geq\beta
 \quad (y\in G).                                         \tag{1.3}
\]

In particular, every \((z_y,N_y)\) in the packet belongs to the two-sided
\(r\)-core of \(f\):

\[
 f(z_y+tN_y)=t\qquad (|t|\leq r).                        \tag{1.4}
\]

Indeed, any calibrated point in each of the two tails and equality in the
Lipschitz inequality force equality on the segment between them.

### Lemma 1.1 (sharp global reach inequality)

For any two members \((z,N)\), \((z',N')\) of the two-sided \(r\)-core,
put \(d=z-z'\).  Then

\[
 \boxed{
 |d|^2-r^2|N-N'|^2
 \geq 2r\,|\langle d,N+N'\rangle|.}                     \tag{1.5}
\]

Consequently

\[
 |z-z'|\geq r|N-N'|,                                    \tag{1.6}
\]

so the oriented normal is unique at a base point and the Gauss map on the
whole core, including different smooth charts, is \(1/r\)-Lipschitz.  For
\(|s|<r\), the normal map \(F_s(z)=z+sN(z)\) is injective and

\[
 (1-|s|/r)|z-z'|
 \leq |F_s(z)-F_s(z')|
 \leq(1+|s|/r)|z-z'|.                                   \tag{1.7}
\]

#### Proof

Crossing the positive endpoint of either ray with the negative endpoint of
the other and using global 1-Lipschitzness gives

\[
 |d+r(N+N')|\geq2r,
 \qquad |d-r(N+N')|\geq2r.                              \tag{1.8}
\]

Square both inequalities and use
\(|N+N'|^2=4-|N-N'|^2\).  The two resulting inequalities are (1.5)
with the two possible signs of the inner product.  Dropping the right side
gives (1.6), and the triangle inequality applied to
\(d+s(N-N')\) gives (1.7).  \(\square\)

This is a genuine reach statement, not a pointwise curvature calculation.
Direction changes of order one either require base separation of order
\(r\), or occur outside the two-sided core at a cut, focal, or singular set.

### Lemma 1.2 (smooth tube Jacobian and flatness)

Suppose near one base point the zero set is a \(C^2\) hypersurface with
shape operator \(S_y=D_\Sigma N_y\), and write \(d\mu=e^{-V}dx/Z\) with
\(V\) convex and \(C^2\).  On the core,

\[
 J_y(t)=\det(I+tS_y),                                    \tag{1.9}
\]

and the conditional density has the form

\[
 q_y(t)=c_y e^{-V(z_y+tN_y)}\det(I+tS_y).                \tag{1.10}
\]

It is log-concave, since

\[
 -\frac{d^2}{dt^2}\log q_y(t)
 =\nabla^2V(z_y+tN_y)[N_y,N_y]
  +\operatorname{tr}\!\left(S_y^2(I+tS_y)^{-2}\right)
 \geq0.                                                  \tag{1.11}
\]

Under (1.3),

\[
 \|S_y\|_{\mathrm{op}}\leq\frac1r,
 \qquad
 \|S_y\|_{\mathrm{HS}}^2\leq\frac{C_\beta}{r^2}.      \tag{1.12}
\]

#### Proof

The eikonal equation gives \(D^2f\,N=0\), so its characteristics are the
straight normal lines in (1.1).  Differentiating the normal map in principal
directions proves (1.9), and then (1.10)--(1.11) follow by the area formula.
No factor \(1+t\kappa_i\) can vanish for \(|t|<r\); applying this at both
ends gives \(|\kappa_i|\leq1/r\).

For completeness, the Hilbert--Schmidt estimate does not hide a dimension
factor.  If \(\sigma_y^2=\operatorname{Var}_{\nu_y}T\), (1.3) and an
independent copy \(T'\) give

\[
 \sigma_y^2=\tfrac12\mathbb E(T-T')^2\geq4\beta^2r^2.   \tag{1.13}
\]

The standard one-dimensional log-concave density bound gives
\(\|q_y\|_\infty\leq C/\sigma_y\leq C/(\beta r)\).  Put
\(W=-\log q_y\).  If \(W'_+(r)>0\), convexity gives

\[
 \beta\leq\int_r^\infty q_y(t)dt
 \leq\frac{q_y(r)}{W'_+(r)},
\]

and otherwise the required upper bound on \(W'_+(r)\) is automatic.
The left tail is identical.  Hence

\[
 W'_+(r)-W'_-(-r)\leq\frac{C}{\beta^2r}.                \tag{1.14}
\]

On \([-r/2,r/2]\), (1.12)'s operator bound makes
\((1+t\kappa_i)^{-2}\geq4/9\).  Integrating (1.11) on that interval and
using (1.14) yields

\[
 \frac{4r}{9}\|S_y\|_{\mathrm{HS}}^2
 \leq \int_{-r/2}^{r/2}W''(t)dt
 \leq\frac{C}{\beta^2r},
\]

which proves (1.12), with for example \(C_\beta=C\beta^{-2}\).
The same proof works with convex second-derivative measures.  \(\square\)

Thus a long balanced smooth ray is flat in the full Hilbert--Schmidt sense.
This does not control how many far-separated flat charts occur.

## 2. Covariance closes the coherent branch

Define the unnormalized normal matrix of the packet by

\[
 M_G=\int_G N_yN_y^T\,d\eta(y),
 \qquad \operatorname{tr}M_G=\alpha_G.                  \tag{2.1}
\]

### Theorem 2.1 (coherent exclusion or forced high rank)

Under (1.1)--(1.3),

\[
 \boxed{4\beta^2r^2M_G\preceq I.}                       \tag{2.2}
\]

In particular, if the packet is \(\kappa\)-coherent in the sense that

\[
 \|M_G\|_{\mathrm{op}}\geq\kappa\alpha_G,               \tag{2.3}
\]

then

\[
 \boxed{r\leq\frac{1}{2\beta\sqrt{\kappa\alpha}}.}     \tag{2.4}
\]

Without (2.3), its effective rank obeys the exact lower bound

\[
 \boxed{
 \operatorname{rank}_{\mathrm{eff}}(M_G)
 :=\frac{\operatorname{tr}M_G}{\|M_G\|_{\mathrm{op}}}
 \geq4\beta^2\alpha_G r^2
 \geq4\beta^2\alpha r^2.}                              \tag{2.5}
\]

#### Proof

Equation (1.13) holds on every ray in \(G\).  Conditional covariance on a
ray is

\[
 \operatorname{Cov}(X\mid y)=\sigma_y^2N_yN_y^T.
\]

The law of total covariance and isotropy give

\[
 I\succeq\int_G\sigma_y^2N_yN_y^T\,d\eta(y)
 \succeq4\beta^2r^2M_G,
\]

which is (2.2).  Taking the operator norm and using (2.3) gives (2.4);
dividing the trace by the same operator-norm bound gives (2.5).  \(\square\)

There is also an exact topological price for the last alternative.  Let
\(\nu_G=(N_\#\eta|_G)/\alpha_G\).  From (2.2),

\[
 \int uu^T\,d\nu_G(u)
 \preceq\frac{1}{4\beta^2\alpha_Gr^2}I.                 \tag{2.6}
\]

Every unit chordal ball \(\{u:|u-v|\leq1\}\) therefore has
\(\nu_G\)-mass at most \((\beta^2\alpha_Gr^2)^{-1}\), because
\(u\cdot v\geq1/2\) there.  A maximal 1-separated subset of the direction
support has at least \(\beta^2\alpha_Gr^2\) points.  Any connected
rectifiable set in the sphere containing that support consequently has

\[
 \mathcal H^1(K)
 \geq\frac12\bigl(\beta^2\alpha_Gr^2-1\bigr).            \tag{2.7}
\]

Indeed, disjoint half-unit balls about the separated points each require at
least one half-unit of connecting length, apart from one ball.  Formula
(2.7) is unweighted.  The unresolved issue is to prevent the connecting
turning or focal graph from lying where the log-concave density is
arbitrarily small.

## 3. Translated thin shell closes the concurrent-radial branch

We first record the exact translated form of the radial estimate.  Assume the
available thin-shell bound in quadratic-radius form

\[
 \operatorname{Var}_\mu|X|^2\leq C_{\mathrm{TS}}n.       \tag{3.1}
\]

This is not a Poincare or KLS inequality for arbitrary Lipschitz functions.

### Lemma 3.1 (thin shell about every center)

For every deterministic \(c\in\mathbb R^n\),

\[
 \operatorname{Var}_\mu|X-c|
 \leq C_{\mathrm{rad}}:=\max\{2C_{\mathrm{TS}},8\}.      \tag{3.2}
\]

#### Proof

Put \(Z=|X-c|^2\).  Isotropy gives \(\mathbb EZ=n+|c|^2\), while

\[
 \operatorname{Var}Z
 \leq2\operatorname{Var}|X|^2+8\operatorname{Var}\langle c,X\rangle
 \leq2C_{\mathrm{TS}}n+8|c|^2
 \leq C_{\mathrm{rad}}\,\mathbb EZ.                     \tag{3.3}
\]

For any nonnegative \(Z\),

\[
 \operatorname{Var}\sqrt Z
 \leq\mathbb E(\sqrt Z-\sqrt{\mathbb EZ})^2
 \leq\frac{\operatorname{Var}Z}{\mathbb EZ}.            \tag{3.4}
\]

Equations (3.3)--(3.4) prove the claim.  \(\square\)

### Theorem 3.2 (concurrent-ray exclusion)

Assume, in addition to (1.1)--(1.3), that there is a point \(c\) such that
for every \(y\in G\) the line \(z_y+\mathbb RN_y\) passes through \(c\),
and the interior of the conditional ray does not cross \(c\).  Then

\[
 \boxed{r\leq
 \frac{\sqrt{C_{\mathrm{rad}}}}{2\beta\sqrt\alpha}.}     \tag{3.5}
\]

#### Proof

Write \(c=z_y+s_yN_y\).  Since the conditional interval stays on one side
of \(s_y\), on that ray

\[
 |z_y+tN_y-c|=|t-s_y|=a_y+\varepsilon_y t,
 \qquad \varepsilon_y\in\{-1,1\}.                       \tag{3.6}
\]

Thus
\(\operatorname{Var}_{\nu_y}|X-c|=\operatorname{Var}_{\nu_y}T
=\sigma_y^2\).  Total variance, (1.13), and Lemma 3.1 give

\[
 C_{\mathrm{rad}}
 \geq\operatorname{Var}_\mu|X-c|
 \geq\int_G\sigma_y^2d\eta(y)
 \geq4\beta^2r^2\alpha,
\]

which is (3.5).  \(\square\)

For a nonbranching congruence, exact concurrence normally puts the common
point at a ray endpoint: two distinct ray interiors cannot meet.  Hence the
one-sidedness in the theorem is the natural concurrent/focal geometry, not
an extra radial smoothness assumption.  A connected smooth hypersurface
whose normal lines all pass through \(c\) is a piece of a sphere, since
\(D_v|z-c|^2=2\langle v,z-c\rangle=0\) for every tangent vector \(v\).

## 4. Consequence for a true T3 extremizer

Let

\[
 A=\sup_{\operatorname{Lip}(f)\leq1}
       \int|f-\mu f|\,d\mu,                              \tag{4.1}
\]

and take a maximizer.  The exact cut--transport argument in
`t3_extremal_report.md` shows that, after centering, it is \(\mu\)-almost
everywhere the signed distance to its zero set.  Its optimal rays are
one-dimensional log-concave needles with the same positive-sign mass on
almost every ray.  The constant-tail lemma in that report supplies universal
constants \(c_0,\alpha_0,\beta_0>0\) and a packet satisfying (1.2)--(1.3)
with

\[
 r=c_0A,qquad \alpha=\alpha_0,qquad\beta=\beta_0.       \tag{4.2}
\]

Theorems 2.1 and 3.2 therefore give the following rigorous trichotomy.

1. If \(\|M_G\|_{\mathrm{op}}\geq\kappa\alpha_0\) for a universal
   \(\kappa>0\), then \(A=O_\kappa(1)\).
2. If the long packet is concurrent-radial, then \(A=O(1)\).
3. Otherwise

   \[
   \operatorname{rank}_{\mathrm{eff}}(M_G)\geq cA^2,    \tag{4.3}
   \]

   different core normals obey the sharp separation (1.5), every smooth
   chart has \(\|S\|_{\mathrm{HS}}^2\leq C/A^2\), and any connected
   completion of the direction support has total unweighted turning at
   least \(cA^2\).

The third alternative is the exact high-rank/nonconcurrent obstruction.  A
nonlinear extremizer must realize its direction changes outside the long
core, in focal or singular geometry, or split into many almost-product
components.  Proving that global log-concavity assigns enough weight to that
completion is precisely the missing high-rank tube inverse.  None of
(1.5), (1.11), (2.2), or translated thin shell proves it.

## 5. A decisive counterexample to rank-only ridge amplification

The following tempting sublemma is false.

> **False rank-only sublemma.**  There is a universal \(c>0\) such that a
> flat normal packet for an isotropic log-concave measure with
> \(\|M/\operatorname{tr}M\|_{\mathrm{op}}\leq1/2\) must lose at least a
> \(c\)-fraction of its normalized tube flux by distance one.

Let \(\mu_m\) be the product of \(m\) one-sided exponential laws,

\[
 d\mu_m(x)=e^{-\sum_{i=1}^m x_i}
             \mathbf1_{\{x_i\geq0\}}dx.                 \tag{5.1}
\]

Translation by the all-ones vector makes this measure isotropic.  Put
\(q=e^{-L}\) and

\[
 E_{m,q}=\bigcup_{i=1}^m\{x_i\geq L\}.                  \tag{5.2}
\]

Its volume, weighted perimeter, and normal matrix are exactly

\[
 v(q)=1-(1-q)^m,
 \qquad P(q)=mq(1-q)^{m-1},
 \qquad M(q)=\frac{P(q)}m I_m.                           \tag{5.3}
\]

Every open facet is flat.  The killed normalized flux at distance
\(0<t<L\) is

\[
 R(t)=mq(1-qe^t)^{m-1},
 \qquad
 \frac{R(t)}{R(0)}
 =\left(\frac{1-qe^t}{1-q}\right)^{m-1}.                \tag{5.4}
\]

To see (5.4), a point on facet \(i\) survives its inward coordinate-normal
ray to time \(t\) exactly when none of the other \(m-1\) coordinates has
entered the competing level \(L-t\).  Independence gives the displayed
factor.

Take \(q=\gamma/m\) and first let \(m\to\infty\).  For fixed \(t\),

\[
 \frac{R(t)}{R(0)}
 \longrightarrow\exp\{-\gamma(e^t-1)\}.                 \tag{5.5}
\]

Now let \(\gamma\downarrow0\).  The normalized normal matrix is \(I_m/m\),
of arbitrarily high effective rank, while \(R(1)/R(0)\to1\).  This
contradicts the proposed sublemma.

The example does **not** contradict Theorems 2.1 or 3.2.  Its coordinate
exponential scale is one and it does not provide fixed two-sided tail mass at
distance \(r\gg1\) on every ray.  Its role is narrower and decisive: high
normal rank, flat Jacobians, convex support, and pairwise ridge incidence do
not by themselves yield a dimension-free collision charge.  Any theorem
used on the T3 packet must retain the balanced long-ray hypothesis or add a
genuine global extremality/product inverse.

## 6. Model audit

| Model | Normal geometry | What stops a long balanced packet |
|---|---|---|
| Sphere | For a rotational packet, \(M/\operatorname{tr}M=I/n\); all normal lines meet at the center. | The coherent estimate is intentionally weak, but Theorem 3.2 applies.  The inward focal time is exactly the radius, and translated thin shell bounds the conditional radial variance. |
| Isotropic cube | A coordinate midplane is coherent.  The boundary of an inner cube has normals \(\{\pm e_i\}\) and is product/high-rank. | The midplane needles have length at most the fixed isotropic half-width \(\sqrt3\).  On the inner boundary, facet rays hit support, a medial sheet, or a facet ridge at a fixed scale; the nonsmooth endpoint charge cannot be replaced by classical curvature. |
| Isotropic regular simplex | The \(n+1\) facet normals form a tight frame: their normalized second moment is \(I/n\). | If the inradius is \(\rho\), uniform covariance is \(n\rho^2/(n+2)I\); isotropy gives \(\rho=\sqrt{(n+2)/n}\).  A facet's intrinsic inradius is \(\rho\sqrt{(n+1)/(n-1)}=O(1)\), so the high-rank normal change is met at ridges on a fixed scale. |
| Product exponentials | Coordinate facets have normalized normal matrix \(I/m\), zero smooth curvature, and many pair ridges. | The one-dimensional normal scale is one.  Equations (5.3)--(5.5) show both why long balanced needles are absent and why rank alone does not amplify ridge loss. |
| Wavy high-reach graph | It may have high rank after subtracting the mean normal, but it retains one coherent uncentered direction. | The coherent covariance estimate, not curvature, is the correct test.  The calculation below shows that arbitrarily large reach alone is harmless. |

Here is the promised wavy calculation.  For \(d\geq1\), \(R>0\), and
\(0<a\leq1\), let

\[
 h(x)=\frac{aR}{\sqrt d}\sum_{i=1}^d\cos(x_i/R),
 \qquad \Sigma=\{(x,h(x)):x\in\mathbb R^d\}.             \tag{6.1}
\]

Then

\[
 |\nabla h|\leq a,
 \qquad \|D^2h\|_{\mathrm{op}}
 \leq\frac{a}{R\sqrt d}.                                \tag{6.2}
\]

The graph has reach at least \(R\sqrt d/a\).  One quick proof uses the
global chord criterion.  If \(p=(x,h(x))\), \(q=(x',h(x'))\), and \(N_p\)
is either unit normal, Taylor's inequality gives

\[
 |\langle q-p,N_p\rangle|
 \leq\frac{a}{2R\sqrt d}|x-x'|^2
 \leq\frac{a}{2R\sqrt d}|p-q|^2,                        \tag{6.3}
\]

which is the positive-reach chord inequality at radius \(R\sqrt d/a\).
On the other hand, for any probability weight \(\theta\) on the graph,

\[
 N=\frac{(-\nabla h,1)}{\sqrt{1+|\nabla h|^2}},
 \qquad
 e_{d+1}^T\left(\int NN^T\,d\theta\right)e_{d+1}
 \geq\frac1{1+a^2}.                                     \tag{6.4}
\]

Thus the uncentered normal matrix is universally coherent.  If \(\theta\)
is uniform on one period cell, the horizontal block of the **centered**
normal covariance has rank \(d\) by symmetry.  This is a useful warning:
centered angular rank can be arbitrarily large on a very high-reach surface,
yet it is irrelevant to (2.2), which correctly uses the uncentered matrix
\(\int NN^T\).  Geometry alone permits (6.1); what it cannot support, by
Theorem 2.1, is a fixed-mass isotropic packet with balanced tails at a length
comparable to that reach.

## 7. Exact remaining statement

For the long packet extracted from a hypothetical bad T3 extremizer, the
coherent and concurrent-radial possibilities are now excluded with universal
constants.  The residual data are simultaneously

\[
 \operatorname{rank}_{\mathrm{eff}}M_G\gtrsim A^2,
 \qquad |z-z'|\gtrsim A|N-N'|,
 \qquad \|S\|_{\mathrm{HS}}^2\lesssim A^{-2}             \tag{7.1}
\]

on a fixed quotient mass, together with exact one-dimensional balance and
one global signed-distance potential.  The direction completion has
unweighted turning \(\gtrsim A^2\), but current arguments do not give that
completion positive \(\mu\)-weight.  The sphere shows that high rank may be
concurrent; product exponentials show that high rank may be product-like and
have very small ridge loss; the wavy graph shows that high centered rank and
large reach may still be coherent.

Accordingly, a valid closing lemma must use **all** of global log-concavity,
joint T3 extremality, fixed two-sided balance, and the singular/focal
completion.  A rank-only collision bound, a smooth-curvature bound, or an
unqualified assertion that large reach aligns the normals is false or
insufficient.
