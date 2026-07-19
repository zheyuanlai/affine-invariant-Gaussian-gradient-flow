# Median extremizers as signed-distance interfaces

## 0. Verdict and exact scope

Let \(F=\operatorname{aff}(\operatorname{supp}\mu)\), equipped with its
induced Euclidean metric, and assume that \(\mu\) is not a point mass.
Define

\[
 D_1(\mu)=\sup_{\operatorname{Lip}(f)\le1}
       \inf_{a\in\mathbb R}\int_F|f-a|\,d\mu .             \tag{0.1}
\]

The infimum over \(a\) is attained precisely at the medians of \(f\).
There is an exact set formulation:

\[
 \boxed{
 D_1(\mu)=
 \sup_{\substack{E\subset F\ {\rm open}\\ \mu(E)=1/2}}
       J_\mu(E),\qquad
 J_\mu(E)=\int_F d(x,\partial_FE)\,d\mu(x).}               \tag{0.2}
\]

No regularity of \(E\) is asserted in (0.2).  The proof uses absolute
continuity on the affine support and a finite first moment.  A non-point
log-concave probability has both properties on \(F\), by Borell's
characterization.  Atoms of
\(f_\#\mu\), including a positive-mass zero plateau, are removed by a
generic linear perturbation; Section 2 gives the complete argument.

For a hypothetical smooth maximizer, its nearest-boundary normal cells
obey an exact Euler equation.  If \(a(y)\) and \(b(y)\) are the \(\mu\)-mass
densities, per unit boundary area at \(y\), on the two halves of the normal
Voronoi cell based at \(y\), then

\[
                              \boxed{a(y)=b(y)}
                                                               \tag{0.3}
\]

for almost every \(y\).  The half-mass constraint is important: the
Lagrange equation initially says \(b-a=c\rho|_{\partial E}\); integrating
over all cells gives \(c=0\).

For a normal variation with speed \(h\), and in the absence of a
codimension-one switching part of the medial axis, the second variation is

\[
 J''(h)=2\int_{\partial E}\rho h^2\,dA
 -\int_{\partial E}\langle G_y\nabla_\Sigma h,
                              \nabla_\Sigma h\rangle\,dA, \tag{0.4}
\]

where \(G_y\) is the reciprocal-Jacobian first-moment tensor of the normal
cell; it is defined in (4.8).  A genuine medial switching hypersurface
contributes an additional nonpositive quadratic form.  Omitting it is an
incorrect second-variation argument.

The cube halfspace, Gaussian halfspace, and median sphere for the isotropic
radial exponential law all satisfy (0.3) and the exact stability
inequality.  For the radial exponential sphere the degree-one stability
condition is an equality, by an integration-by-parts identity using the
radial median.  A half-volume simplex cap fails (0.3), including on normal
cells based outside the support.

Finally, (0.2) shows that a universal bound \(J_\mu(E)\le C\) for isotropic
log-concave \(\mu\) is exactly the median-centered T3 form of KLS.  Global
maximality and isotropy alone do not prove it: an isotropic atomic regular
simplex, which is not log-concave, has a globally maximizing
signed-distance interface with \(J\simeq\sqrt n\).  Thus global
log-concavity across the normal cells remains indispensable.

## 1. Signed distance dominates every median-centered function

All boundaries and distances in this section are relative to \(F\).

### Lemma 1.1 (signed distance to the zero interface)

Let \(f:F\to\mathbb R\) be 1-Lipschitz and suppose that \(0\) is a median:

\[
                      \mu\{f\ge0\}\ge\tfrac12,\qquad
                      \mu\{f\le0\}\ge\tfrac12.             \tag{1.1}
\]

Put \(Z=\{f=0\}\) and define

\[
 g_Z(x)=
 \begin{cases}
  d(x,Z),&f(x)>0,\\
  0,&f(x)=0,\\
  -d(x,Z),&f(x)<0.
\end{cases}                                               \tag{1.2}
\]

The set \(Z\) is nonempty: (1.1) supplies points of the support on both
sides unless one of them already belongs to \(Z\), and continuity on the
convex affine support gives an intermediate zero.

Then \(g_Z\) is 1-Lipschitz, \(0\) is a median of \(g_Z\), and

\[
                              |g_Z(x)|\ge|f(x)|
                              \qquad(x\in F).              \tag{1.3}
\]

#### Proof

For \(z\in Z\), the Lipschitz inequality gives
\(|f(x)|\le|x-z|\).  Taking the infimum over \(z\) proves (1.3).

On each of \(\{f>0\}\) and \(\{f<0\}\), the assertion follows from the
1-Lipschitz property of distance to a fixed closed set.  If \(f(x)>0\)
and \(f(y)<0\), continuity on the segment \([x,y]\subset F\) supplies
\(z\in[x,y]\cap Z\).  Hence

\[
 |g_Z(x)-g_Z(y)|
 =d(x,Z)+d(y,Z)
 \le|x-z|+|z-y|=|x-y|.                                  \tag{1.4}
\]

The cases where one point belongs to \(Z\) follow by the same estimate.
The signs of \(g_Z\) are the signs of \(f\), so (1.1) proves the median
claim.  \(\square\)

The lemma includes a zero plateau of arbitrary \(\mu\)-mass.  What it does
not yet provide is a positive set of mass exactly one half.

### Lemma 1.2 (signed distance of an open set)

If \(E\subset F\) is nonempty and open, define

\[
 s_E(x)=d(x,E^c)-d(x,E).                                  \tag{1.5}
\]

At each point one of the two distances is zero.  Consequently

\[
 |s_E(x)|=d(x,\partial_FE),                               \tag{1.6}
\]

and \(s_E\) is 1-Lipschitz.

#### Proof

On \(E\), \(s_E=d(\,\cdot\,,E^c)\); on \(E^c\), it equals
\(-d(\,\cdot\,,E)\).  The same-side estimate is the Lipschitz property of
distance.  If \(x\in E\) and \(y\in E^c\), the segment first exits \(E\)
at some \(z\in\partial_FE\), and

\[
 s_E(x)-s_E(y)
 \le|x-z|+|z-y|=|x-y|.                                   \tag{1.7}
\]

This also proves (1.6).  \(\square\)

If \(\mu(E)=1/2\), then \(0\) is a median of \(s_E\), even when
\(\mu(\partial_FE)>0\).  Therefore

\[
                              J_\mu(E)\le D_1(\mu).        \tag{1.8}
\]

## 2. Zero plateaus, level atoms, and the proof of (0.2)

We first remove atoms of the one-dimensional distribution \(f_\#\mu\).

### Lemma 2.1 (generic tilt has no level atoms)

Let \(\mu\) be absolutely continuous on the \(k\)-dimensional affine space
\(F\), and let \(f:F\to\mathbb R\) be Lipschitz.  For every
\(\varepsilon>0\), there is a vector \(v\in F-F\) with \(0<|v|\le1\)
such that

\[
                         h(x)=f(x)+\varepsilon\langle v,x\rangle
                                                                  \tag{2.1}
\]

satisfies

\[
                              \mu\{h=t\}=0
                              \qquad(t\in\mathbb R).       \tag{2.2}
\]

#### Proof

By Rademacher's theorem, \(\nabla f\) exists almost everywhere.  The
probability distribution of \(\nabla f\) has at most countably many atoms.
Choose \(v\) in the unit ball so that \(-\varepsilon v\) is not one of
them.  This also covers \(k=1\), where restricting \(v\) to the unit
sphere would leave only two choices.

If \(\{h=t\}\) had positive \(\mu\)-mass, it would have positive
\(k\)-dimensional Lebesgue measure.  At almost every density point of this
level set where \(h\) is differentiable, \(\nabla h=0\).  This standard
fact follows directly from differentiability by approaching the density
point through the level set in \(k\) independent limiting directions.
Thus \(\nabla f=-\varepsilon v\) on a set of positive \(\mu\)-mass,
contrary to the choice of \(v\).  \(\square\)

### Theorem 2.2 (exact half-mass set representation)

Let \(\mu\) be absolutely continuous on \(F\) and have finite first
moment.  Then (0.2) holds.

#### Proof

Only the inequality \(D_1\le\sup_EJ_\mu(E)\) remains after (1.8).  Fix a
1-Lipschitz \(f\), translated so that \(0\) is a median.  Let
\(\varepsilon_j\downarrow0\), and use Lemma 2.1 to choose \(v_j\) for

\[
                         h_j=f+\varepsilon_j\langle v_j,x\rangle .
                                                                  \tag{2.3}
\]

Its law is atomless.  Continuity of its distribution function gives a
median \(m_j\) satisfying

\[
                   \mu\{h_j>m_j\}=\mu\{h_j<m_j\}=\tfrac12. \tag{2.4}
\]

Set \(E_j=\{h_j>m_j\}\).  This is open and has mass one half.  Since
\(\partial_FE_j\subset\{h_j=m_j\}\) and
\(\operatorname{Lip}(h_j)\le1+\varepsilon_j\),

\[
 J_\mu(E_j)
 \ge {1\over1+\varepsilon_j}
       \int|h_j-m_j|\,d\mu .                              \tag{2.5}
\]

For an integrable function \(u\), write

\[
                         \Phi_\mu(u)=\inf_a\int|u-a|\,d\mu .
\]

The triangle inequality gives

\[
                    |\Phi_\mu(u)-\Phi_\mu(w)|
                    \le\int|u-w|\,d\mu.                   \tag{2.6}
\]

Here

\[
 \int|h_j-f|\,d\mu
 \le\varepsilon_j\int|x|\,d\mu\longrightarrow0.           \tag{2.7}
\]

Since medians minimize absolute deviation, (2.5)--(2.7) imply

\[
 \liminf_jJ_\mu(E_j)\ge\Phi_\mu(f).
\]

Taking the supremum over \(f\) proves (0.2).  \(\square\)

### 2.1 What “atoms” means here

A non-point log-concave probability is absolutely continuous on its affine
support and hence has no spatial atoms.  The approximation above is needed
for atoms of \(f_\#\mu\), not atoms of \(\mu\).

For a genuinely atomic probability, the literal right side of (0.2) may
be empty: an atom of mass greater than one half prevents the existence of
any set of mass one half.  One would have to split atoms in an enlarged
probability space.  Thus (0.2) is not asserted for arbitrary atomic
measures.

### 2.2 Lower-dimensional support

If \(\dim F=k<n\), all preceding arguments are carried out in \(F\).
Relative distance from \(x\in F\) to a relative boundary
\(\partial_FE\subset F\) is the same Euclidean distance computed in
\(\mathbb R^n\).  Every 1-Lipschitz function on \(F\) has a
1-Lipschitz McShane extension to \(\mathbb R^n\).  Hence both sides of
(0.2) agree with their ambient formulations.  No nonsingular covariance
in \(\mathbb R^n\) is being assumed.

## 3. Normal Voronoi cells and the first shape derivative

This section assumes \(F=\mathbb R^k\), \(k\ge2\),

\[
                         d\mu(x)=\rho(x)\,dx,              \tag{3.1}
\]

where \(\rho>0\) is \(C^1\), and \(E\) has a \(C^2\) embedded boundary
\(\Sigma=\partial E\) of finite weighted area.  Noncompact boundaries are
allowed, but variations are compactly supported and every displayed
quantity is assumed integrable.

Orient the unit normal \(n\) from \(\Sigma\) into \(E\).  For almost every
point of \(\mathbb R^k\), the nearest point on \(\Sigma\) is unique.  If
that point is \(y\), write

\[
                           x=y+t n(y),                    \tag{3.2}
\]

where \(t>0\) on \(E\).  Let

\[
 I_y=(-\ell_-(y),\ell_+(y))                              \tag{3.3}
\]

be the maximal interval on which \(y\) remains the unique nearest point
and no focal value is crossed.  Put

\[
 S_y=D_\Sigma n(y),\qquad
 j_y(t)=\det(I+tS_y)>0,\qquad
 w_y(t)=\rho(y+tn(y))j_y(t).                              \tag{3.4}
\]

The normal area formula, after deleting the Lebesgue-null medial axis,
gives

\[
\begin{aligned}
 \mu(E)&=\int_\Sigma a(y)\,dA(y),&
 a(y)&=\int_0^{\ell_+(y)}w_y(t)\,dt,\\
 \mu(E^c)&=\int_\Sigma b(y)\,dA(y),&
 b(y)&=\int_{-\ell_-(y)}^0w_y(t)\,dt,                    \tag{3.5}\\
 J_\mu(E)&=\int_\Sigma\int_{I_y}|t|w_y(t)\,dt\,dA(y).
\end{aligned}
\]

These are densities per unit surface area, not normalized conditional
probabilities.  The normalized normal-cell conditional is
\(w_y(t)dt/(a(y)+b(y))\).

### Theorem 3.1 (first variation)

Let \(h\in C_c^2(\Sigma)\), and form the normal graph

\[
 \Sigma_\varepsilon
 =\{y+\varepsilon h(y)n(y):y\in\Sigma\}.                  \tag{3.6}
\]

Choose \(E_\varepsilon\) so that locally its signed-distance coordinate
is \(t-\varepsilon h(y)+o(\varepsilon)\).  Then

\[
\boxed{
\begin{aligned}
 {d\over d\varepsilon}\bigg|_0\mu(E_\varepsilon)
   &=-\int_\Sigma \rho(y)h(y)\,dA(y),\\
 {d\over d\varepsilon}\bigg|_0J_\mu(E_\varepsilon)
   &=\int_\Sigma [b(y)-a(y)]h(y)\,dA(y).
\end{aligned}}                                             \tag{3.7}
\]

#### Proof

The volume formula is the standard transport formula with \(n\) oriented
into \(E\): moving the boundary in the \(+n\) direction removes the
normal slab of thickness \(\varepsilon h\).

Fix \(x=y+tn(y)\) having a unique nearest point.  Differentiating the
nearest-point equations gives

\[
 {d\over d\varepsilon}\bigg|_0s_{E_\varepsilon}(x)
 =-h(y).                                                   \tag{3.8}
\]

Therefore the derivative of its absolute value is \(-h(y)\) for \(t>0\)
and \(+h(y)\) for \(t<0\).  The cut locus has Lebesgue measure zero.
Moreover a normal graph satisfies

\[
 |d(x,\Sigma_\varepsilon)-d(x,\Sigma)|
 \le |\varepsilon|\|h\|_\infty+O(\varepsilon^2),
\]

locally uniformly.  This supplies domination for the difference
quotients.  Integrating (3.8) with the normal area formula (3.5) proves
the second line of (3.7).  \(\square\)

### Corollary 3.2 (every normal cell is bisected)

Suppose \(\mu(E)=1/2\) and \(E\) is stationary for \(J_\mu\) under all
smooth mass-preserving normal variations.  Then

\[
                               a(y)=b(y)                   \tag{3.9}
\]

for surface-almost every \(y\).

#### Proof

Stationarity and (3.7) say

\[
 \int_\Sigma(b-a)h\,dA=0
 \quad\hbox{whenever}\quad
 \int_\Sigma\rho h\,dA=0.                                \tag{3.10}
\]

The fundamental lemma, including variations transferring weighted area
between different connected components, gives one constant \(c\) such
that

\[
                              b(y)-a(y)=c\rho(y).          \tag{3.11}
\]

On the other hand, (3.5) and the half-mass constraint give

\[
 \int_\Sigma(b-a)\,dA=\mu(E^c)-\mu(E)=0.                  \tag{3.12}
\]

Since \(0<\int_\Sigma\rho\,dA<\infty\), equations
(3.11)--(3.12) imply \(c=0\).  \(\square\)

Thus the half-mass constraint does force equality of the two masses on
almost every normal Voronoi cell.  It is not an extra conjectural
localization assertion.

### 3.1 Vanishing density and boundary portions outside the support

If \(\rho\) vanishes on part of \(\Sigma\), a variation supported there
has zero first-order volume cost.  Stationarity then directly forces
\(a=b\) there.  On \(\{\rho>0\}\), the preceding multiplier argument
still gives \(b-a=c\rho\), and global half mass again gives \(c=0\).
This observation is important for a convex-body cap: a complete
hypersurface may have basepoints outside the convex support while their
normal cells still carry mass.

For a discontinuous density such as the uniform law on a convex body,
(3.7) is first justified for variations transverse to the regular part of
the support, or after log-concave smoothing.  At a genuine free boundary,
lower-dimensional contact Voronoi cells must be included separately;
Section 6.5 records the distinction.

## 4. Second variation and the medial-axis correction

Assume now that \(\rho\in C^2\), \(\Sigma\in C^3\), and the normal graph
variation in (3.6) has zero acceleration.  We first calculate at a point
\(x=y+tn(y)\) whose closest boundary branch remains unique.

### Lemma 4.1 (second derivative of the signed distance)

At such a point,

\[
\begin{aligned}
 \dot t&=-h(y),\\
 \ddot t&=-t\,
 \left\langle (I+tS_y)^{-1}\nabla_\Sigma h(y),
                         \nabla_\Sigma h(y)\right\rangle. \tag{4.1}
\end{aligned}
\]

#### Proof

Use local orthonormal coordinates \(q\) on \(\Sigma\), and minimize

\[
 D_\varepsilon(q)=
 |x-[X(q)+\varepsilon h(q)n(q)]|^2.                       \tag{4.2}
\]

At the minimizer \(q=0\),

\[
\begin{aligned}
 D_{\varepsilon\varepsilon}&=2h^2,\\
 D_{qq}&=2(I+tS_y),\\
 D_{q\varepsilon}&=-2t\nabla_\Sigma h.                   \tag{4.3}
\end{aligned}
\]

The envelope Hessian is therefore

\[
 \ddot D
 =2h^2-2t^2
   \langle(I+tS_y)^{-1}\nabla_\Sigma h,\nabla_\Sigma h\rangle.
                                                                  \tag{4.4}
\]

Since \(D=t^2\) and \(\dot t=-h\), comparison of
\(\ddot D=2\dot t^2+2t\ddot t\) with (4.4) proves (4.1).
\(\square\)

The distributional second derivative of absolute value adds the crossing
term

\[
 {d^2\over d\varepsilon^2}|t_\varepsilon|\bigg|_0
 =\operatorname{sgn}(t)\ddot t+2h(y)^2\delta_0(t).        \tag{4.5}
\]

Define the positive-semidefinite tangent tensor

\[
 G_y=\int_{I_y}|t|(I+tS_y)^{-1}w_y(t)\,dt.                \tag{4.6}
\]

### Theorem 4.2 (second variation without switching)

Suppose there is a closed exceptional set \(C\) such that, locally,
\(|C_r|=O(r^2)\), the nearest-boundary branch is unique and \(C^2\) on
\(\mathbb R^k\setminus C\) for all sufficiently small variations, and the
second difference quotients are uniformly integrable off \(C_r\).
Equivalently for the calculation, assume that the moving cut/focal locus
contributes \(o(\varepsilon^2)\) to \(J_\mu(E_\varepsilon)\).  Then

\[
 \boxed{
 J_\mu''(h)=
 2\int_\Sigma\rho h^2\,dA
 -\int_\Sigma
       \langle G_y\nabla_\Sigma h,\nabla_\Sigma h\rangle\,dA.}
                                                                  \tag{4.7}
\]

If \(E\) is a stationary local maximizer at mass one half, then for all
\(h\) satisfying

\[
                           \int_\Sigma\rho h\,dA=0,        \tag{4.8}
\]

one has

\[
 2\int_\Sigma\rho h^2\,dA
 \le\int_\Sigma
       \langle G_y\nabla_\Sigma h,\nabla_\Sigma h\rangle\,dA.
                                                                  \tag{4.9}
\]

#### Proof

Insert (4.1) into (4.5), integrate by (3.5), and use
\(w_y(0)=\rho(y)\).  This proves (4.7).

By Corollary 3.2, the first derivative of \(J_\mu\) vanishes for every
normal speed, not only for speeds satisfying (4.8).  Hence the normal
acceleration needed to make a path volume-preserving to second order does
not contribute to \(J_\mu''\).  Local maximality gives (4.9).
\(\square\)

### 4.3 A codimension-one medial axis adds a negative term

Suppose on an open region the unsigned distance is the minimum of two
smooth nearest-branch distances:

\[
 d_\varepsilon(x)
 =\min\{d_{1,\varepsilon}(x),d_{2,\varepsilon}(x)\}.      \tag{4.10}
\]

Let \(C=\{d_{1,0}=d_{2,0}\}\) be a regular hypersurface and put
\(\alpha_i=\dot d_{i,0}\).  Since

\[
 \min(r_1,r_2)={r_1+r_2-|r_1-r_2|\over2},
\]

the distributional identity \(d^2|q+\varepsilon\beta|/d\varepsilon^2
=2\beta^2\delta_0(q)\) gives the additional term

\[
 -\int_C {(\alpha_1-\alpha_2)^2\over
             |\nabla(d_{1,0}-d_{2,0})|}\,
             \rho\,d\mathcal H^{k-1}.                    \tag{4.11}
\]

It is nonpositive.  Generic multiple-branch strata give the analogous
sum after stratification.  Consequently the general second variation has
the form

\[
 J_\mu''(h)=
 2\int_\Sigma\rho h^2\,dA
 -\int_\Sigma\langle G_y\nabla h,\nabla h\rangle\,dA
 -\mathcal M(h),\qquad \mathcal M(h)\ge0.                 \tag{4.12}
\]

At a local maximum, (4.12) does **not** imply (4.9), because the medial
term itself can supply stability.

### 4.4 Explicit audit in one dimension

Let \(\rho=c\) on \([-L,L]\), and let the boundary consist of two points
\(p<q\).  Direct integration gives

\[
 {1\over c}J(p,q)
 ={(p+L)^2\over2}+{(q-p)^2\over4}+{(L-q)^2\over2}.        \tag{4.13}
\]

For velocities \(h_1,h_2\), the naive sum of the two boundary crossing
terms would give \(2c(h_1^2+h_2^2)\).  The Voronoi switch occurs at
\((p+q)/2\).  Formula (4.11) subtracts

\[
                         {c\over2}(h_1+h_2)^2,             \tag{4.14}
\]

and the resulting Hessian is exactly

\[
 c\left({3\over2}h_1^2-h_1h_2+{3\over2}h_2^2\right),     \tag{4.15}
\]

as obtained by differentiating (4.13).  Thus the medial correction is a
real second-order contribution, despite the medial set itself having
Lebesgue measure zero.

## 5. What first and second variation do—and do not—buy

For a smooth half-mass maximizer, (3.9) says that the normal-cell
disintegration is balanced:

\[
 d\mu(x)
 =w_y(t)\,dt\,dA(y),\qquad
 \int_{t>0}w_y(t)\,dt=\int_{t<0}w_y(t)\,dt.               \tag{5.1}
\]

This is the exact signed-distance version of balanced transport rays.  If

\[
 m_y={\int t w_y(t)\,dt\over\int w_y(t)\,dt},\qquad
 \sigma_y^2={\int(t-m_y)^2w_y(t)\,dt\over\int w_y(t)\,dt},
                                                                  \tag{5.2}
\]

then isotropy gives the same conditional covariance constraint as in the
transport-ray formulation:

\[
 \int \sigma_y^2 n(y)n(y)^T\,\eta(dy)\preceq I,           \tag{5.3}
\]

where \(\eta(dy)=(a(y)+b(y))dA(y)\).  Its trace permits
\(\int\sigma_y^2d\eta\) of order \(k\).  Thin shell gives only a fourth
moment of order \(k\).  Thus (3.9), isotropy, and thin shell reproduce the
known diffuse-direction obstruction.

When \(\mathcal M=0\), stability adds the weighted surface Poincaré
inequality (4.9).  The coefficient \(G_y\) itself contains the first
moment of the cell length.  Summing coordinate test functions introduces
a factor of the dimension, and no dimension-free estimate for
\(J_\mu(E)\) follows formally.  When \(\mathcal M\ne0\), even (4.9) need
not hold.

Therefore the variational calculation is a genuine extra necessary
condition, but it is not a completed inverse.  Any use of (4.9) must
either rule out codimension-one medial switching for a global maximizer,
or retain and exploit the full negative form \(\mathcal M\).

## 6. Exact stress tests

### 6.1 The isotropic cube and a coordinate halfspace

Let

\[
 \mu=\operatorname{Unif}[-a,a]^k,\qquad a=\sqrt3,\qquad
 E=\{x_1>0\}.                                             \tag{6.1}
\]

Write \(y=(0,x_2,\ldots,x_k)\) and \(n=e_1\).  For
\(y\in[-a,a]^{k-1}\), the normal cell is \([-a,a]\), \(S=0\), and
\(\rho=(2a)^{-k}\).  Hence

\[
 a(y)=b(y)=\rho a,\qquad
 J_\mu(E)=\mathbb E|X_1|={a\over2}={\sqrt3\over2},\qquad
 G_y=\rho a^2I_{T_y\Sigma}.                              \tag{6.2}
\]

For a smooth graph variation supported in the cross-section,

\[
 J_\mu''(h)
 =\rho\int_{[-a,a]^{k-1}}
       [\,2h^2-a^2|\nabla h|^2\,]\,dy.                   \tag{6.3}
\]

The first-order volume constraint is \(\int h=0\).  The sharp Neumann
Poincaré inequality on the \((k-1)\)-cube is

\[
 \int h^2\le {4a^2\over\pi^2}\int|\nabla h|^2.
\]

Since \(8/\pi^2<1\), (6.3) is strictly negative for every nonzero
admissible \(h\).  Thus the coordinate halfspace passes both variation
tests with dimension-free slack.

### 6.2 Gaussian halfspace

For \(\mu=\gamma_k=N(0,I)\) and \(E=\{x_1>0\}\),

\[
 a(y)=b(y)={1\over2}\varphi_{k-1}(y),\qquad
 J_{\gamma_k}(E)=\sqrt{2\over\pi}.                        \tag{6.4}
\]

Since

\[
 \int_{\mathbb R}|t|\varphi(t)\,dt
 =\sqrt{2\over\pi}=2\varphi(0),
\]

one has \(G_y=2\rho(y)I_{T_y\Sigma}\), and

\[
 J_{\gamma_k}''(h)
 =2\varphi(0)\int_{\mathbb R^{k-1}}
       [h^2-|\nabla h|^2]\,d\gamma_{k-1}.                 \tag{6.5}
\]

The Gaussian Poincaré inequality makes (6.5) nonpositive for
\(\int h\,d\gamma_{k-1}=0\).  Equality occurs for linear \(h\), the
infinitesimal rotations of a halfspace.

### 6.3 The isotropic radial exponential sphere

Let

\[
 d\mu(x)=c_k e^{-\lambda|x|}\,dx,\qquad
 \lambda=\sqrt{k+1}.                                      \tag{6.6}
\]

If \(R=|X|\), then \(R\) has the gamma density

\[
 p(r)={\lambda^k\over\Gamma(k)}r^{k-1}e^{-\lambda r},
 \qquad \mathbb E R^2={k(k+1)\over\lambda^2}=k,           \tag{6.7}
\]

so \(\mu\) is isotropic.  Let \(r_0\) be the median of \(R\), and take
\(E=\{|x|<r_0\}\).  Then

\[
 J_\mu(E)=\mathbb E|R-r_0|
 \le\mathbb E|R-\mathbb ER|
 \le\sqrt{\operatorname{Var}R}
 ={ \sqrt k\over\sqrt{k+1}}<1.                           \tag{6.8}
\]

At \(y=r_0\theta\), orient \(n=-\theta\), so
\(t=r_0-R\), \(S=-r_0^{-1}I\), and

\[
 w_y(t)=c_ke^{-\lambda R}\left({R\over r_0}\right)^{k-1}.
                                                                  \tag{6.9}
\]

Rotational symmetry and the median identity give \(a(y)=b(y)\).  The
tensor in (4.6) is \(G_y=gI_{T_y\Sigma}\), where

\[
 g=c_k\int_0^\infty |r-r_0|e^{-\lambda r}
                      \left({r\over r_0}\right)^{k-2}dr. \tag{6.10}
\]

There is an exact identity

\[
                         {k-1\over r_0^2}g
                         =2c_ke^{-\lambda r_0}
                         =2\rho(r_0).                     \tag{6.11}
\]

To prove it, set

\[
 A=\int_0^\infty|r-r_0|r^{k-2}e^{-\lambda r}dr.
\]

The two integrals of \(r^{k-1}e^{-\lambda r}\) on the sides of \(r_0\)
are equal.  Expanding the absolute value therefore gives

\[
 A=r_0\left[
   \int_0^{r_0}r^{k-2}e^{-\lambda r}dr
  -\int_{r_0}^\infty r^{k-2}e^{-\lambda r}dr\right].      \tag{6.12}
\]

Integrating

\[
 (r^{k-1}e^{-\lambda r})'
 =(k-1)r^{k-2}e^{-\lambda r}
   -\lambda r^{k-1}e^{-\lambda r}
\]

on the two sides of \(r_0\), and again using the median equality, yields

\[
                         A={2r_0^k e^{-\lambda r_0}\over k-1},
\]

which is (6.11).

The first nonzero eigenvalue of the sphere of radius \(r_0\) is
\((k-1)/r_0^2\).  Equations (4.7) and (6.11) therefore show

\[
 J_\mu''(h)\le0\qquad\left(\int_{\mathbb S^{k-1}}h=0\right), \tag{6.13}
\]

with equality precisely in the degree-one spherical harmonics.  The
radial exponential sphere is thus an exact semistable test, not merely an
order-of-magnitude example.

### 6.4 A half-volume cap in an isotropic regular simplex

Represent a regular simplex using a facet as base:

\[
 K=\{(y,z):0\le z\le H,\quad y\in(1-z/H)B\},              \tag{6.14}
\]

where \(B\subset\mathbb R^{k-1}\) is the base facet, centered at the
origin.  Translation to the barycenter and a scalar dilation make the
regular simplex isotropic, so the following normal-cell calculation is
unchanged up to scale.

The cap \(E_c=\{z>c\}\) has relative volume

\[
                         \mu(E_c)=(1-c/H)^k.
\]

Thus its half-volume level is

\[
                         c=H(1-2^{-1/k}).                 \tag{6.15}
\]

The boundary is the complete plane \(z=c\).  On a vertical normal line
based at \(y\), the simplex occupies

\[
             0\le z\le z_+(y),\qquad
             z_+(y)=H(1-\|y\|_B),                        \tag{6.16}
\]

whenever \(y\in B\).  For \(y\in(1-c/H)B\), let \(\rho\) denote the
constant volume density.  The two cell masses are

\[
                {a(y)\over\rho}=z_+(y)-c,\qquad
                {b(y)\over\rho}=c.                       \tag{6.17}
\]

They are unequal except on one level set of the gauge \(\|y\|_B\).
For

\[
 y\in B\setminus(1-c/H)B,
\]

the basepoint \((y,c)\) is outside the support, but its normal cell
contains only mass below the plane: \(a(y)=0<b(y)\).  Hence the simplex
cap fails the first Euler equation (3.9), including on its
support-contact completion cells.  Half volume by itself does not imply
cellwise balance; stationarity does.

### 6.5 Relative interfaces versus complete ambient boundaries

For \(J_\mu(E)\), the relevant boundary is the complete relative boundary
of the open set \(E\subset F\).  A hyperplane cut of a convex body is a
complete hyperplane, including the part outside the support.  Normal cells
based on that exterior part can carry positive \(\mu\)-mass, as in
(6.16)--(6.17).

If one instead truncates the interface to
\(\Sigma\cap\operatorname{supp}\mu\), its boundary on the support creates
lower-dimensional Voronoi cells of positive ambient volume.  Those are
contact cells and contribute separate first- and second-variation terms.
Discarding them changes the signed-distance functional.  The two
formulations agree only after the exterior continuation or the contact
cells have been specified and included.

## 7. Global maximality plus isotropy is insufficient without log-concavity

The following exact example separates the metric variational facts from
global log-concavity.

Let \(N=k+1\) be even, and choose the vertices
\(v_1,\ldots,v_N\in\mathbb R^k\) of a centered isotropic regular simplex:

\[
 |v_i|^2=k,\qquad
 \langle v_i,v_j\rangle=-1\quad(i\ne j),\qquad
 {1\over N}\sum_i v_iv_i^T=I.                            \tag{7.1}
\]

Every pair of distinct vertices has distance

\[
                              D=\sqrt{2(k+1)}.            \tag{7.2}
\]

Partition the vertices into two classes \(P,N\) of equal size, and let
\(\mu\) be the uniform probability on the vertices.  Define the open
Voronoi-dominance set

\[
 E=\{x:d(x,P)<d(x,N)\}.                                   \tag{7.3}
\]

For \(p\in P\), the 2-Lipschitz function
\(d(\,\cdot\,,N)-d(\,\cdot\,,P)\) has value \(D\) at \(p\).  Hence
\(d(p,\partial E)\ge D/2\).  If \(q\in N\), the midpoint
\((p+q)/2\) is equidistant from \(p\) and \(q\).  Every other simplex
vertex is at distance \(\sqrt3D/2>D/2\) from that midpoint.  Thus

\[
 d(v_i,\partial E)={D\over2}\qquad(1\le i\le N),          \tag{7.4}
\]

and

\[
                         J_\mu(E)={D\over2}
                         =\sqrt{k+1\over2}.               \tag{7.5}
\]

This is globally optimal.  Indeed, the values of any 1-Lipschitz function
on the vertices lie in an interval of length at most \(D\).  For the
uniform law on an even number of points, the mean absolute deviation from
a median is at most half that range.  Equality is obtained by assigning
\(+D/2\) on \(P\) and \(-D/2\) on \(N\).  Therefore

\[
                              D_1(\mu)=J_\mu(E)=D/2.       \tag{7.6}
\]

The measure is isotropic and the interface is a globally maximizing signed
distance interface, yet its value grows like \(\sqrt k\).  The failed
hypothesis is log-concavity.  Replacing each vertex atom by a sufficiently
small equal-radius ball and applying one scalar covariance normalization
gives an absolutely continuous, still non-log-concave example with

\[
                              J_\mu(E)\ge c\sqrt k.        \tag{7.7}
\]

The bisector interface has a large stratified medial structure; its
second variation cannot be audited by (4.7) while omitting
\(\mathcal M\).

## 8. Exact equivalence with T3 and the remaining obstruction

For a non-point log-concave probability, Theorem 2.2 gives

\[
 \sup_{\operatorname{Lip}(f)\le1}
      \int|f-\operatorname{med}f|\,d\mu
 =
 \sup_{\substack{E\ {\rm open}\\\mu(E)=1/2}}
      \int d(x,\partial_FE)\,d\mu(x).                    \tag{8.1}
\]

Consequently the assertion

\[
 \sup_{\mu(E)=1/2}J_\mu(E)\le C
 \quad\hbox{for every isotropic log-concave }\mu          \tag{8.2}
\]

is neither weaker nor stronger than the median-centered first-moment
concentration target: it is exactly T3.  By E. Milman's equivalence theorem
for log-concave probabilities, (8.2) is equivalent up to universal
constants to a dimension-free Cheeger lower bound.

The first Euler equation converts a smooth maximizer into globally
compatible balanced normal cells.  The covariance and thin-shell
constraints on those cells are the same ones already present in the
balanced transport-ray formulation.  The second variation supplies either
the interface inequality (4.9) or the fuller inequality (4.12), but no
proved argument converts either one into (8.2).  Such a conversion would
itself complete KLS.

Thus the remaining statement can be isolated as follows.

> **Signed-distance inverse needed.**  Let \(\mu\) be isotropic and
> log-concave, and let \(\Sigma\) be a globally maximizing half-mass
> signed-distance interface.  Use global log-concavity, the cellwise
> balance (3.9), and the full stability form including medial/contact
> terms to prove
> \[
>                 \int d(x,\Sigma)\,d\mu(x)\le C.
> \]

Neither cellwise balance alone nor global maximality plus isotropy without
log-concavity suffices, by Section 7.

## 9. Approximation and regularity audit

### 9.1 Stability of the functional under measure approximation

For probabilities \(\mu,\nu\) with finite first moments,

\[
                              |D_1(\mu)-D_1(\nu)|
                              \le W_1(\mu,\nu).            \tag{9.1}
\]

Indeed, for any 1-Lipschitz \(f\), any \(a\), and any coupling
\(\pi\) of \(\mu,\nu\),

\[
 \left|\int|f-a|\,d\mu-\int|f-a|\,d\nu\right|
 \le\int|x-y|\,d\pi(x,y).
\]

Take the infimum over \(a\), then the supremum over \(f\), and optimize the
coupling.

Convolution of a log-concave probability on \(F\) with a small Gaussian
on \(F\) produces a positive smooth log-concave density and converges in
\(W_1\).  Second moments also converge; recentering and covariance
normalization tend to the identity when the original law is isotropic.
Thus a universal estimate proved for all smooth positive log-concave
measures, without an extra regularity-dependent constant, transfers to
arbitrary log-concave measures.

### 9.2 Approximation by smooth half-mass sets

The proof of Theorem 2.2 already gives half-mass level sets with values
arbitrarily close to \(D_1\), but their boundaries need not be smooth.
One may mollify the defining Lipschitz function, choose nearby regular
levels, and correct the small mass error by a small normal bump on a
regular patch.  On a large ball this changes the boundary in Hausdorff
distance by \(o(1)\), hence changes its distance function uniformly by
\(o(1)\); outside the ball, the finite first moment controls the error.
This proves that smooth half-mass sets approximate the **value** of the
supremum.

It does not produce a smooth maximizer.  A maximizing sequence may develop
corners, topology changes, boundary escaping through a zero-density
region, or a codimension-one medial switching set.  Near-maximality also
does not imply the exact Euler equation (3.9) or an exact second-variation
inequality without a quantitative Ekeland-type argument and compactness in
a topology controlling normal cells.

### 9.3 What may and may not be assumed

The following chain is rigorous:

\[
\text{arbitrary log-concave }\mu
\ \xrightarrow[\text{Theorem 2.2}]{}
\ \text{half-mass signed-distance maximizing sequence}.  \tag{9.2}
\]

The following stronger chain is not currently justified:

\[
\text{maximizing sequence}
\ \Longrightarrow\
\text{smooth globally maximizing interface}
\ \Longrightarrow\
\text{formula (4.9) with no medial/contact term}.         \tag{9.3}
\]

Accordingly, the Euler and Hessian formulas above are exact necessary
conditions under their stated smoothness and normal-cell hypotheses, not
a hidden regularity reduction.  Any complete proof through this route must
either establish the missing compactness and regularity with uniform
constants, or work directly with nonsmooth stratified interfaces and keep
all medial and support-contact contributions.
