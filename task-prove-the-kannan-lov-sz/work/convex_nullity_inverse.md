# Convex nullity inverse: exact rigidity and the quantitative obstruction

## Verdict

Let
\[
 \pi_c(dx)=Z(c)^{-1}
 \exp\{\langle c,x\rangle-t|x|^2/2-W(x)\}\,dx,
 \qquad t>0,                                                   \tag{0.1}
\]
where \(W\) is proper lower-semicontinuous convex and may take the value
\(+\infty\).  There are two sharply different conclusions.

1. **Exact zero for a genuine posterior certificate is completely rigid.**
   If the variance deficit or the standardized Brenier-Jacobian defect
   vanishes in a direction \(u\) at even one tilt, then the domain of \(W\)
   is a cylinder in direction \(u\), \(W\) is globally affine in that
   direction, and every posterior has an independent \(N(\cdot,t^{-1})\)
   factor in that direction.  For a weighted family, the global Gaussian
   factor contains the whole range of its direction matrix \(R\).  This
   theorem includes nonsmooth potentials and hard support.
2. **The proposed dimension-free quantitative global inverse is false from
   the stated inputs alone.**  A bounded bulk-resolvent defect can miss
   singular curvature completely.  Even the genuine Brenier defect is
   posterior-local: widely separated tilts may sit deep inside unrelated
   affine cells of one polyhedral convex function.  Every local defect can
   tend to zero and the direction matrix can have effective rank \(n\), while
   the global lineality space is zero and no fixed-mass family of local cells
   has a coherent focal center.

The safe dimension-free quantitative output is therefore a **separate local
Gaussian product approximation at each tilt**.  Converting those separate
products into a common global factor or a fixed-mass coherent focal family
requires an additional overlap/incidence hypothesis.  Log-concavity of the
tilt law alone does not supply that hypothesis without a dimension-free
expansion statement, which is excluded here.

Strictly speaking, the proposed approximate dichotomy has no truth value
until “approximately global,” “near-null leaf,” “focal,” and “fixed-mass
coherence” are assigned a metric, scale, and reference measure.  The no-go
results below apply to the natural strong readings: global epigraph/support
control, or posterior-independent affine directions and focal centers with
an error tending to zero with the stated defect.

The sharp obstructions are:

* singular curvature on a set of vanishing Lebesgue thickness, which the
  bounded resolvent does not see;
* exponential loss of posterior overlap at tilt distance
  \(|c-c'|/\sqrt t\); and
* absence of a universal convex error bound converting small Bregman or
  transport defect into distance to a complete null leaf or a common focus.

---

## 1. Definitions which distinguish the valid and invalid statements

Let \(E\) be the linear space parallel to the affine hull of
\(\operatorname{dom}W\).  All statements below are relative to \(E\); a
direction orthogonal to \(E\) cannot have saturated variance.  Translating
the affine hull, if necessary, causes no change.

### 1.1 Global affine directions

Define the affine lineality space
\[
\begin{aligned}
 \mathcal L(W)=\{u\in E:\ &\text{there is a number }\ell(u)\text{ such that}\\
 &W(x+su)=W(x)+s\ell(u)
 \quad\text{for every }x\in\operatorname{ri}(\operatorname{dom}W),
 \ s\in\mathbb R\}.
                                                                  \tag{1.1}
\end{aligned}
\]
The equality includes the assertion that \(x+\mathbb Ru\) remains in the
domain.  The set \(\mathcal L(W)\) is a vector space, and \(\ell\) is a linear
functional on it.  If \(L=\mathcal L(W)\) and \(P=P_{L^\perp}\), then
\[
 \overline{\operatorname{dom}W}=K_0+L,\qquad
 W(y+z)=W_0(y)+\ell(z),\quad y\in K_0,\ z\in L.                   \tag{1.2}
\]
Thus (1.1), unlike an almost-everywhere Hessian condition, records both the
potential and its hard support.

### 1.2 The genuine posterior defects

Let \(m_c=\mathbb E_cX\).  Let \(T_c\) be the standardized centered Brenier
map
\[
 T_c:\gamma_E\longrightarrow
 \operatorname{Law}_{\pi_c}\big(\sqrt t(X-m_c)\big).              \tag{1.3}
\]
Caffarelli contraction, interpreted by convex approximation when necessary,
gives an almost-everywhere symmetric Jacobian
\[
 0\preceq H_c=DT_c\preceq I_E.                                  \tag{1.4}
\]
For a unit \(u\in E\), set
\[
 \alpha_c(u)=1-t\,\operatorname{Var}_{\pi_c}\langle X,u\rangle,
 \qquad
 \beta_c(u)=\mathbb E_{\gamma_E}|H_cu-u|^2.                      \tag{1.5}
\]
Both defects include contact with nonsmooth faces and support boundaries.

### 1.3 The bulk resolvent

When \(W\) is finite and \(C^2\), define
\[
 K_t(x)=t(tI+\nabla^2W(x))^{-1},\qquad
 \rho_c(u)=\mathbb E_c\langle u,(I-K_t(X))u\rangle.              \tag{1.6}
\]
Here \(0\preceq K_t\preceq I\), and Brascamp--Lieb gives
\[
 t\operatorname{Var}_c\langle X,u\rangle
 \le \mathbb E_c\langle u,K_t(X)u\rangle
 =1-\rho_c(u).                                                  \tag{1.7}
\]
Thus \(\alpha_c(u)\ge\rho_c(u)\).  Small \(\alpha\) forces small
\(\rho\), but small \(\rho\) does not force small \(\alpha\).

For nonsmooth \(W\), inserting the Alexandrov Hessian into (1.6) is not a
closed operation.  It discards the singular part of the distributional
Hessian and all normal-cone curvature of a hard support.  Section 4 gives an
explicit failure even under smooth approximation.

---

## 2. Complete exact-zero theorem, including hard support

### Theorem 2.1 (one direction)

For a unit vector \(u\in E\), the following are equivalent.

1. \(\alpha_c(u)=0\) for some tilt \(c\).
2. \(\beta_c(u)=0\) for some tilt \(c\).
3. \(u\in\mathcal L(W)\).
4. For every tilt \(c'\), in the decomposition
   \(E=\mathbb Ru\oplus u^\perp\),
   \[
   \pi_{c'}=
   N\left(\frac{\langle c',u\rangle-\ell(u)}t,\frac1t\right)
   \otimes\nu_{P c'}                                               \tag{2.1}
   \]
   for a log-concave transverse law \(\nu_{P c'}\).

The statement is relative to the affine hull and remains valid when
\(W=+\infty\) off a closed convex support.

#### Proof

Let
\[
 f(G)=\langle T_c(G),u\rangle .
\]
Gaussian Poincare and (1.4) give
\[
\begin{aligned}
 t\operatorname{Var}_c\langle X,u\rangle
 &=\operatorname{Var}f\\
 &\le\mathbb E|H_cu|^2
 \le\mathbb E\langle u,H_cu\rangle
 \le1.                                                         \tag{2.2}
\end{aligned}
\]
If \(\alpha_c(u)=0\), equality holds throughout.  Since
\(|H_cu|^2\le\langle u,H_cu\rangle\le1\),
\[
\begin{aligned}
\beta_c(u)
&=1-2\mathbb E\langle u,H_cu\rangle+\mathbb E|H_cu|^2\\
&\le1-\mathbb E\langle u,H_cu\rangle=0.                        \tag{2.3}
\end{aligned}
\]
Thus (1) implies (2).

The same calculation without assuming equality records the useful general
comparison
\[
 0\le\beta_c(u)\le\alpha_c(u).                                 \tag{2.3a}
\]

Conversely, \(\beta_c(u)=0\) gives \(H_cu=u\) almost everywhere.  Write
\(G=zu+w\).  In the weak Sobolev sense,
\[
 \partial_z\{T_c(zu+w)-zu\}=0.
\]
Hence
\[
 T_c(zu+w)=zu+A(w).                                             \tag{2.4}
\]
Symmetry of \(H_c\) also gives
\[
 \nabla_w\langle T_c(zu+w),u\rangle=P H_cu=0.
\]
The \(u\)-component of \(A(w)\) is therefore constant.  After absorbing that
constant into the mean already removed in (1.3), (2.4) has the product form
\[
 T_c(zu+w)=zu+A_\perp(w),\qquad A_\perp(w)\in u^\perp.          \tag{2.5}
\]
Its pushforward is a standard Gaussian in direction \(u\), independent of
the transverse pushforward.  Undoing the standardization proves that
\(\pi_c\) is a product with variance \(1/t\) in direction \(u\).

The topological support of a product with a nondegenerate Gaussian factor is
a cylinder.  Since the density of \(\pi_c\) is positive on the relative
interior of \(\operatorname{dom}W\), this gives
\[
 \overline{\operatorname{dom}W}=K_0+\mathbb Ru.                 \tag{2.6}
\]
Compare the logarithm of the product density with (0.1) on the relative
interior.  The quadratic term already supplies \(-ts^2/2\), so the remaining
coefficient of \(s\) must be constant:
\[
 W(y+su)=W_0(y)+s\ell(u).                                      \tag{2.7}
\]
Thus (2) implies (3).  This argument uses neither differentiability of
\(W\) nor full support.

Substituting (2.7) into (0.1) and completing the square proves (4) for every
tilt.  Finally (4) plainly implies (1).  This closes the cycle. \(\square\)

The only approximation hidden above is the standard extended-valued form of
Caffarelli contraction.  One may take finite smooth convex \(W_j\) which
epiconverge to \(W\).  The target laws converge in
second moment, their Brenier maps converge in \(L^2(\gamma_E)\), and the
limit is a Lipschitz cyclically monotone map satisfying (1.4) almost
everywhere.  The proof then uses only weak derivatives of that limit.  The
support cylinder is recovered from the support of the product law, rather
than from an interior Hessian.

### Corollary 2.2 (a weighted family)

Let \(\omega\) be a finite positive measure on tilts, let \(u(c)\) be
measurable and unit, and put
\[
 R=\int u(c)u(c)^{\mathsf T}\,d\omega(c).                       \tag{2.8}
\]
If
\[
 \int\beta_c(u(c))\,d\omega(c)=0
 \quad\text{or}\quad
 \int\alpha_c(u(c))\,d\omega(c)=0,                              \tag{2.9}
\]
then
\[
 \operatorname{ran}R\subset\mathcal L(W).                      \tag{2.10}
\]
Consequently, if either common definition
\[
 r_{\mathrm{op}}(R)=\frac{\operatorname{tr}R}{\|R\|_{\mathrm{op}}},
 \qquad
 r_{\mathrm{st}}(R)=\frac{(\operatorname{tr}R)^2}{\operatorname{tr}R^2}
                                                                        \tag{2.11}
\]
is at least \(r\), then the global Gaussian factor has dimension at least
\(r\).

Indeed, the nonnegative integrand in (2.9) vanishes almost everywhere, so
Theorem 2.1 puts almost every \(u(c)\) in \(\mathcal L(W)\).  Equation (2.10)
then follows from
\[
 a^{\mathsf T}Ra=\int\langle a,u(c)\rangle^2\,d\omega(c).
\]

There is no separate concurrent or focal exact-zero branch.  Such a branch
can only describe incomplete, approximate, posterior-local leaves.

---

## 3. Exact local null leaves of a convex function

For \(p\in E\), define the subgradient fiber
\[
 F_p=(\partial W)^{-1}(p)
 =\operatorname*{argmin}_{x}\{W(x)-\langle p,x\rangle\}.        \tag{3.1}
\]
Every \(F_p\) is closed and convex, and
\[
 W(x)-W(y)=\langle p,x-y\rangle,\qquad x,y\in F_p.              \tag{3.2}
\]
These are the canonical nonsmooth affine leaves.  They include:

* full-dimensional affine cells of a polyhedral maximum;
* rays of a strictly convex norm or perspective cone;
* faces created by a hard support; and
* affine fibers of a smooth constant-rank gradient map.

If \(W\) is \(C^2\) on an open set and \(\nabla^2W\) has locally constant
rank, the constant-rank theorem makes the fibers of \(\nabla W\) smooth
submanifolds tangent to \(\ker\nabla^2W\).  Convexity of (3.1) then shows
that each connected leaf is an open subset of an affine plane, and \(W\)
is affine on it.

This local fact does not provide global incidence.  Different fibers may be
parallel, concurrent, terminate on a focal set, or form an arbitrary regular
polyhedral subdivision.  Global lineality is the much stronger assertion
that one fixed subspace is contained in the direction space of every
relevant fiber and that the leaves are complete in those directions.

For a quantitative version, the natural scalar is the Bregman slack
\[
 D_W(x,y;p)=W(x)-W(y)-\langle p,x-y\rangle,\qquad p\in\partial W(y).
                                                                    \tag{3.3}
\]
It is nonnegative and vanishes exactly when \(x\in F_p\).  There is no
universal convex error bound of the form
\[
 \operatorname{dist}(x,F_p)^2\le C\,D_W(x,y;p).                 \tag{3.4}
\]
Flat powers, hinges with remote breakpoints, and polyhedral cells violate
every dimension-free choice of \(C\).  Hence small analytic nullity does not
by itself give a metric leaf, leaf length, or focal radius.

---

## 4. The bulk-resolvent blind spot is sharp

For finite \(C^2\) potentials, exact zero of (1.6) does imply a global
factor.  Indeed, the integrand is nonnegative.  If \(\rho_c(u)=0\), then
\[
 \nabla^2W(x)u=0
\]
for every \(x\), by positivity of the posterior density and continuity.
Symmetry makes all mixed derivatives vanish as well, and integration gives
\(u\in\mathcal L(W)\).

This statement fails completely for nonsmooth limits.

### 4.1 The one-dimensional cusp

Take \(E=\mathbb R\), \(c=0\), and
\[
 W_0(x)=\lambda|x|,\qquad \lambda>0.                            \tag{4.1}
\]
The Alexandrov Hessian is zero almost everywhere.  Thus the naive
pointwise resolvent reports
\[
 \rho_0(1)=0,                                                   \tag{4.2}
\]
although \(W_0\) has no lineality and the posterior is not Gaussian.

The failure persists under smooth convex approximation.  Let
\[
 W_\varepsilon(x)=\lambda\sqrt{x^2+\varepsilon^2},\qquad
 W_\varepsilon''(x)=
 \frac{\lambda\varepsilon^2}{(x^2+\varepsilon^2)^{3/2}}.       \tag{4.3}
\]
The scalar resolvent defect is
\[
 \rho_\varepsilon
 =\mathbb E_{\pi_\varepsilon}
 \frac{W_\varepsilon''(X)}{t+W_\varepsilon''(X)}.               \tag{4.4}
\]
Put \(a=(\lambda\varepsilon^2/t)^{1/3}\).  The posterior densities are
uniformly bounded near zero.  On \(|x|\le a\), bound the integrand by one;
on \(|x|>a\), bound it by \(W_\varepsilon''/t\).  This gives
\[
 \rho_\varepsilon
 \le C_{\lambda,t}a+
 \frac{C}{t}\int_a^\infty
       \frac{\lambda\varepsilon^2}{x^3}\,dx
 \le C_{\lambda,t}\varepsilon^{2/3}\longrightarrow0.           \tag{4.5}
\]

On the other hand, \(\pi_\varepsilon\) converges in moments to
\[
 \pi_0(dx)\propto e^{-t x^2/2-\lambda|x|}\,dx.
\]
Integration by parts against its score gives the exact identity
\[
 t\,\mathbb E_{\pi_0}X^2+\lambda\mathbb E_{\pi_0}|X|=1.
                                                                    \tag{4.6}
\]
Therefore, if \(\alpha_{\varepsilon,0}\) denotes the variance defect for
\(W_\varepsilon\),
\[
 \lim_{\varepsilon\downarrow0}\alpha_{\varepsilon,0}(1)
 =\lambda\mathbb E_{\pi_0}|X|>0.                               \tag{4.7}
\]
The bounded resolvent goes to zero while the genuine variance defect stays
strictly positive.

Thus the quantitative failure already occurs within the class of smooth
finite convex potentials \(W_\varepsilon\); what degenerates is the
thickness of the curvature layer.  Smoothness without a uniform curvature
regularity or thickness modulus does not repair the inverse theorem.

The mechanism is precise: curvature of size \(1/\varepsilon\) is compressed
onto a layer whose posterior mass goes to zero.  The bounded function
\(H(tI+H)^{-1}\) pays at most the mass of that layer.  A polyhedral seam or a
hard facet has exactly the same pathology.

Thus no inverse theorem covering nonsmooth support may use the bulk
resolvent as its only defect.  It must retain the Brenier/variance
certificate or an unbounded measure-valued curvature/contact term.

---

## 5. What the Brenier defect does prove quantitatively

The exact argument has a stable posterior-local form.  Decompose
\(G=Zu+PG\), and define
\[
 A_{c,u}(w)=\mathbb E_Z[P T_c(w+Zu)]\in u^\perp.                \tag{5.1}
\]
Gaussian Poincare, first for the active scalar component and then
conditionally on each \(u\)-fiber for the transverse component, gives
\[
\mathbb E\left|
T_c(G)-\{u\langle G,u\rangle+A_{c,u}(PG)\}
\right|^2
\le2\beta_c(u).                                                 \tag{5.2}
\]
Consequently there is a transverse law \(\nu_{c,u}\) such that
\[
 W_2^2\left(
 \pi_c,\,
 N(\langle m_c,u\rangle,t^{-1})\otimes\nu_{c,u}
 \right)
 \le\frac{2\beta_c(u)}t.                                       \tag{5.3}
\]
This is dimension free, includes hard support, and is the strongest
unconditional inverse supplied by one posterior defect.

It is intentionally local in \(c\).  It neither says that the couplings for
different tilts use the same transverse map nor that \(W\) is close to an
affine function outside the mass seen by \(\pi_c\).  A convex function may be
changed beyond \(L\) posterior standard deviations at a cost
\(\exp(-\Theta(L^2))\) in (5.3).

---

## 6. Posterior overlap is the sharp globalization cost

Let \(b(c)=\log Z(c)\).  Brascamp--Lieb gives
\[
 \nabla^2b(c)=\operatorname{Cov}_{\pi_c}(X)\preceq t^{-1}I.    \tag{6.1}
\]
For \(c_0,c_1\) and their midpoint \(c_m=(c_0+c_1)/2\),
\[
\mathbb E_{c_i}
\left(\frac{d\pi_{c_m}}{d\pi_{c_i}}\right)^2
=\frac{Z(c_0)Z(c_1)}{Z(c_m)^2}
\le\exp\left\{\frac{|c_1-c_0|^2}{4t}\right\}.                  \tag{6.2}
\]
If \(0\le F\le1\), Cauchy--Schwarz yields
\[
 \mathbb E_{c_m}F
 \le
 \exp\left\{\frac{|c_1-c_0|^2}{8t}\right\}
 \sqrt{\mathbb E_{c_i}F}.                                     \tag{6.3}
\]
Thus a defect of size \(\eta\) can be transferred across a midpoint only
while
\[
 \frac{|c_1-c_0|^2}{t}\lesssim\log(1/\eta).                    \tag{6.4}
\]
Gaussian tail events show that this exponential scale is sharp.

The high effective rank of
\[
 R=\int u(c)u(c)^{\mathsf T}\,d\omega(c)
\]
contains no information about the overlap graph of the tilts.  The mass of
\(\omega\) may be divided among arbitrarily many components separated by
more than the scale in (6.4).  Convexity of \(W\) does not reconnect the
separate Brenier maps.

If \(\omega\) is the natural same-time tilt law, it is log-concave but need
not be uniformly strongly log-concave.  Showing that a positive amount of
its low-defect mass lies in one controlled-overlap component is a
dimension-free conductance/expansion assertion.  It cannot be inserted here
without assuming the very type of estimate ruled out in the question.

---

## 7. Mandatory model tests

### 7.1 Norms: the coherent radial survivor

For
\[
 W(x)=\lambda|x|
\]
on \(\mathbb R^n\), away from the origin,
\[
 \nabla^2W(x)=\frac{\lambda}{|x|}
\left(I-\frac{xx^{\mathsf T}}{|x|^2}\right),\qquad
 \nabla^2W(x)x=0.                                              \tag{7.1}
\]
The exact pointwise null leaves are rays, all concurrent at the origin.
A posterior whose mode \(q\) satisfies \(\sqrt t\,|q|\gg\sqrt n\) sees only
a narrow tube around one ray.  Its active radial direction
\(u=q/|q|\) has a small Brenier/variance defect, although \(W\) has no global
line.  Sampling \(q/|q|\) over many angles gives a high-rank direction
matrix.  This is a legitimate quantitative radial branch, with the coherent
center \(0\).

In one dimension this can be seen explicitly.  For \(c>\lambda\),
\[
 \pi_c(dx)\propto
 \begin{cases}
 e^{-t x^2/2+(c-\lambda)x}\,dx,&x\ge0,\\
 e^{-t x^2/2+(c+\lambda)x}\,dx,&x<0.
 \end{cases}                                                   \tag{7.2}
\]
As \((c-\lambda)/\sqrt t\to\infty\), this converges in weighted moments to
\(N((c-\lambda)/t,t^{-1})\).  Hence
\(\alpha_c(1)\to0\), and (2.2)--(2.3) give
\(\beta_c(1)\le\alpha_c(1)\to0\).  The missing half-line terminates at the
focal kink \(0\).

The \(\ell_1\) norm is an additional warning: it is affine on every orthant
and has Alexandrov Hessian zero almost everywhere, so the raw resolvent sees
no curvature at all.  Its genuine posterior defect still detects crossings
of the coordinate seams.

### 7.2 Perspective cones: exact local rays and a hard apex

On \(\mathbb R^k\times\mathbb R\), let
\[
 W(z,y)=
 \begin{cases}
 |z|^2/(2y),&y>0,\\
 0,&(z,y)=(0,0),\\
 +\infty,&\text{otherwise}.
 \end{cases}                                                   \tag{7.3}
\]
This is the perspective of a quadratic and is convex and one-homogeneous.
Its Hessian is
\[
 \nabla^2W(z,y)=
 \begin{pmatrix}
 y^{-1}I&-y^{-2}z\\
 -y^{-2}z^{\mathsf T}&y^{-3}|z|^2
 \end{pmatrix},
 \qquad
 \nabla^2W(z,y)\binom zy=0.                                   \tag{7.4}
\]
The null leaves are rays from the apex.  Tilts concentrating far from the
apex and the hard boundary have small defect in the local radial direction;
directions can have effective rank \(k+1\), but the leaves have the exact
coherent center \(0\).  Translating the perspective translates the focus.

This model validates a radial/focal alternative only when a fixed positive
amount of the tilt family samples the same perspective piece.  A maximum of
several translated perspective pieces is still convex and permits different
packets to see different apices, separated by tie regions.  No common focus
follows without overlap of those packets.

### 7.3 Polyhedral maxima: the global inverse countermodel

The strongest obstruction is that local affine cells contain no canonical
null direction at all.

Set \(t=1\) for clarity.  Fix \(L\gg1\), and let
\[
 w(s)=\lambda\sum_{k=-M}^{M}|s-4Lk|,\qquad
 W(x)=\sum_{j=1}^n w(x_j).                                    \tag{7.5}
\]
This is a finite convex polyhedral function.  Since a sum of finite maxima
is again a finite maximum over all choices of active affine pieces, \(W\) is
a polyhedral maximum; moreover \(e^{-W}\) is integrable.  It is affine on
each open box between consecutive breakpoint hyperplanes,
\[
 B_\kappa=\prod_{j=1}^n(4L\kappa_j,4L(\kappa_j+1)).             \tag{7.6}
\]
Its global lineality space is zero because the slope changes in every
coordinate.  Its Alexandrov Hessian is nevertheless zero almost everywhere,
so a naive pointwise bulk resolvent declares every direction exactly flat at
every tilt.  Moreover, for every nonzero \(u\), translating far enough along
a coordinate on which \(u\) is nonzero crosses a slope jump of fixed size.
Thus \(W(x+su)-W(x)\) cannot be uniformly approximated by one linear
function of \(s\), even on a one-dimensional nonzero global factor.

Choose \(n\) boxes with centers \(q_i\) pairwise separated by \(DL\), where
\(D\) may be arbitrarily large, and let \(a_i\) be the constant gradient of
\(W\) on the \(i\)-th box.  Set
\[
 c_i=q_i+a_i.                                                  \tag{7.7}
\]
The affine function
\[
 x\longmapsto W(q_i)+\langle a_i,x-q_i\rangle
\]
supports \(W\) globally and agrees with it on the box.  Relative to the
Gaussian \(\gamma_i=N(q_i,I)\), the posterior \(\pi_{c_i}\) has density
proportional to
\[
 h_i(x)=
 \exp\{-W(x)+W(q_i)+\langle a_i,x-q_i\rangle\},
 \qquad 0<h_i\le1,\quad h_i=1\text{ on }B_i.                   \tag{7.8}
\]
If the Gaussian distance from \(q_i\) to every face of \(B_i\) is at least
\(L\), then
\[
 \mathbb E_{\gamma_i}\big[(1+|X-q_i|^2)\mathbf1_{B_i^c}\big]
 \le Cn(1+n+L^2)e^{-L^2/2}.                                   \tag{7.9}
\]
Here is the moment bookkeeping behind the next step.  Put
\[
 d_i=\mathbb E_{\gamma_i}(1-h_i),\qquad
 A_i(u)=\mathbb E_{\gamma_i}
 \big[(1-h_i)\langle u,X-q_i\rangle^2\big].
\]
The normalizer is \(1-d_i\ge\gamma_i(B_i)\), and for
\(Y=\langle u,X-q_i\rangle\),
\[
 \mathbb E_{\pi_{c_i}}Y^2=\frac{1-A_i(u)}{1-d_i},\qquad
 |\mathbb E_{\pi_{c_i}}Y|^2
 \le\frac{d_iA_i(u)}{(1-d_i)^2}.                               \tag{7.9a}
\]
For \(L\) large, \(d_i\le1/2\), so (7.9)--(7.9a) imply, uniformly
for every unit \(u\),
\[
 0\le
 1-\operatorname{Var}_{\pi_{c_i}}\langle X,u\rangle
 \le Cn(1+n+L^2)e^{-L^2/2}.                                   \tag{7.10}
\]
The lower bound is strong log-concavity.  By (2.3),
\[
 \beta_{c_i}(u)
 \le Cn(1+n+L^2)e^{-L^2/2}.                                   \tag{7.11}
\]

Assign \(u(c_i)=e_i\) and give the \(n\) tilts equal weight.  Then
\[
 R=\frac1nI,\qquad r_{\mathrm{op}}(R)=r_{\mathrm{st}}(R)=n,    \tag{7.12}
\]
while every defect in (7.11) tends to zero as \(L\to\infty\) and
\(\mathcal L(W)=\{0\}\).

The directions can also arise from one fixed measurable set, rather than
being assigned artificially.  Because the selected boxes are disjoint,
define
\[
 S\cap B_i=
 \{x\in B_i:\langle e_i,x-q_i\rangle\ge0\},                    \tag{7.13}
\]
and define \(S\) arbitrarily off their union.  Under
\(\gamma_i=N(q_i,I)\), this is the centered \(e_i\)-halfspace except on
\(B_i^c\).  Equations (7.8)--(7.9) therefore give
\[
 \pi_{c_i}(S)=\frac12+o_L(1),\qquad
 \operatorname{Cov}_{\pi_{c_i}}(\mathbf1_S,X)
 =I(1/2)e_i+o_L(1),                                           \tag{7.14}
\]
uniformly over the finite selected family after taking \(L\) large relative
to \(n\).  Thus the genuine active directions are
\(e_i+o_L(1)\), their Gaussian-profile deficits tend to zero, and their
direction matrix still has effective rank \(n-o_L(n)\).

There is also no forced focal conclusion.  Inside an affine box every
direction is null, so the assignment \(u(c_i)=e_i\) carries no information
about the box's singular skeleton.  The boxes and their nearest vertices can
be separated by \(DL\) with \(D\to\infty\).  Any neighborhood on the local
focal scale \(O(L)\) then captures at most one weight \(1/n\), not a fixed
positive fraction.  Thickening each \(c_i\) slightly gives a positive-volume
family with the same properties.

This construction refutes any dimension-free theorem whose hypotheses are
only:

* positive total family weight;
* high effective rank of \(R\); and
* small separate posterior Brenier or bulk-resolvent defects.

To exclude it one must assume a fixed-mass connected overlap component,
control the amount of defect on the tie regions, or provide an independent
incidence principle for the singular skeleton.

If the intended \(\omega\) is not an arbitrary positive weighting but a
submeasure of a prescribed natural tilt law, with density at most one and
mass bounded below independently of dimension, that requirement must be
stated explicitly.  The construction above is then a local obstruction, not
by itself a counterexample to that stronger formulation.  The unresolved
step is precisely to prove that such fixed natural mass cannot be dispersed
among many low-overlap affine cells.  Log-concavity of the tilt law does not
prove this without an additional expansion or incidence argument.

---

## 8. What a viable quantitative theorem would have to assume

A well-posed replacement must specify all of the following.

1. **Use a closed defect.**  Hard support and singular seams require
   \(\alpha_c\), \(\beta_c\), or a measure-valued subgradient/contact energy.
   The bounded bulk resolvent alone is insufficient.
2. **State the topology of approximation.**  Global epigraph or support
   closeness is impossible from posterior data because remote convex
   modifications are exponentially invisible.  Posterior \(W_2\), relative
   entropy, or a common weighted Bregman norm is meaningful.
3. **Supply common mass and overlap.**  There must be a subfamily of weight at
   least \(\kappa>0\) connected by edges
   \(|c-c'|\le C\sqrt t\), with a controlled path or likelihood-ratio budget.
   High effective rank does not imply this.
4. **Condition the direction frame.**  To extract a \(k\)-plane from
   directional defects, one needs a lower singular-value or effective-rank
   condition on the directions within the same overlap component, not merely
   after summing disconnected packets.
5. **Give a leaf error bound and completeness scale.**  Small Bregman slack
   must control distance to a leaf, and leaves must persist for a specified
   two-sided length.  Otherwise hinges and perspective apices survive.
6. **For a focal branch, assume nontrivial transverse curvature.**  A
   full-dimensional affine polyhedral cell has every direction null and no
   canonical focus.  A focal center becomes meaningful only after a
   rank-balanced nonzero curvature block or an equivalent normal-cone
   incidence condition is present.

Under these additions, the plausible exact geometry is:

* complete parallel leaves give a cylinder/global Gaussian factor;
* incomplete leaves with a common singular termination give a cone,
  perspective, or radial focus; and
* a polyhedral family requires an explicit fixed-mass incidence theorem for
  its faces and vertices.

The last bullet is not a consequence of convexity plus high effective rank.
It is precisely the missing fixed-mass center-coherence input.

---

## 9. Final conclusion

The exact-zero problem has a clean and complete answer:
\[
\boxed{
\text{genuine zero posterior defect}
\Longleftrightarrow
\text{global affine lineality of }W\text{ and its support}
\Longleftrightarrow
\text{a Gaussian product factor}.
}                                                               \tag{9.1}
\]
For a weighted family, the factor contains all of \(\operatorname{ran}R\).

The analogous approximate global statement does not follow dimension
freely.  Norms and perspectives exhibit the coherent radial survivor, while
generic polyhedral maxima exhibit disconnected affine-cell packets with
arbitrary null directions and no common focus.  The sharp analytic losses
are the singular-curvature blind spot and the exponential overlap factor
\(\exp(|c-c'|^2/(8t))\).

Therefore the defensible inverse theorem at present is the local product
estimate (5.3).  A global cylinder-versus-focal dichotomy needs a new
fixed-mass overlap/incidence hypothesis; it cannot be obtained from the
listed assumptions without importing a dimension-free expansion principle.
