# Singular eikonal curvature: the exact ray measure and the failure of a quadratic ambient measure

## 0. Verdict

Let \(f\) be the signed distance of a half-mass interface and assume that
its normal cells are individually balanced.  On a smooth normal chart, put
\(u=\nabla f\), let \(\tau\) be the half-transport density, and write
\(d\mu=e^{-V}dx\).  The familiar pointwise identity is

\[
 -\partial_u Lf
 =\|\nabla ^2f\|_{HS}^{2}+\nabla ^2V(u,u).                 \tag{0.1}
\]

There is an exact nonsmooth completion of the **integrated** identity, but
it is not of the form \(\tau K\) for a locally finite ambient curvature
measure \(K\).  The correct object lives first on the compactified normal-ray
graph.  If \(q_y(t)dt\) is the (not necessarily normalized) density on a
balanced normal cell, \(W_y\) is its two-sided tail transport, and

\[
 \kappa_y=D^2(-\log q_y)
\]

is the nonnegative one-dimensional Hessian measure, then the completed
curvature measure is

\[
 d\mathcal C_y=W_y\,d\kappa_y
       +e_y^-\delta_{a_y}+e_y^+\delta_{b_y}.              \tag{0.2}
\]

The endpoint coefficients are explicit nonnegative traces.  The exact
one-dimensional Stieltjes identity proved below is

\[
 \boxed{\mathcal C_y([a_y,b_y])=2q_y(0).}                 \tag{0.3}
\]

Consequently, after integration over the interface quotient,

\[
 \boxed{\mathcal C(\widehat{\mathcal R})=2P_\mu(E).}      \tag{0.4}
\]

In a smooth interior chart, (0.2) is precisely

\[
 W_y\,d\kappa_y
 =\tau\bigl(\|\nabla ^2f\|_{HS}^{2}
                  +\nabla ^2V(u,u)\bigr)d\mu.           \tag{0.5}
\]

At a regular two-sheet medial hypersurface \(M\), however, the two endpoint
atoms push forward to

\[
 \boxed{e_1\,dA_1+e_2\,dA_2
       =e^{-V}|u^+-u^-|\,d\mathcal H^{d-1}\lfloor M.}     \tag{0.6}
\]

Thus the singular term is **linear** BV turning, not quadratic turning.
It is also exactly the trace, over translation modes, of the medial
switching form in the signed-distance second variation.  At a hard convex
support the analogous endpoint term is the normal escape flux

\[
 e^{-V}|u\cdot n_K|\,d\mathcal H^{d-1}\lfloor\partial K. \tag{0.7}
\]

Equations (0.5)--(0.7) have the same scaling as perimeter and give a
rigorous mixed smooth/singular budget.  They do **not** produce the
lower-semicontinuous quadratic measure requested in the proposed curvature
route:

* the transport trace satisfies \(\tau=0\) at every medial or hard-support
  endpoint, so a positive endpoint atom cannot equal \(\tau K\) for any
  locally finite Radon measure \(K\);
* a Gauss jump has zero relaxed capacity for the weighted quadratic energy.
  An explicit logarithmic smoothing below converges uniformly through
  1-Lipschitz functions while
  \(\int |s|\,\|\nabla ^2f_j\|^2\to0\).  A single-scale smoothing of the
  same kink has a positive limit, and another smoothing makes the energy
  diverge.  Hence quadratic singular turning is not a measure and is not
  regularization independent.

The literal assertion that some nonnegative \(K\) obeys only
\(\int\tau\,dK\le2P\) is vacuous (one may discard all singular turning, or
take \(K=0\)).  The substantive assertion needed by the global-compatibility
argument is that \(K\) also retain medial and support turning.  That
assertion is false.  The route remains open only in the modified form using
the mixed measure \(\mathcal C\), or equivalently the full medial switching
form, and not an ambient quadratic Hessian measure.

## 1. What distributional Hessians do and do not provide

Let \(U\subset\mathbb R^d\) be open.  If \(f\) is locally semiconcave on
\(U\), then \(u=\nabla f\) belongs locally to \(BV\), and

\[
 D u=D^2f
\]

is a symmetric matrix-valued Radon measure.  On the rectifiable jump set
\(J_u\), its jump part is

\[
 D^j u=(u^+-u^-)\otimes n_u\,
             \mathcal H^{d-1}\lfloor J_u.                \tag{1.1}
\]

Symmetry forces \(u^+-u^-\) to be parallel to \(n_u\).  A general signed
distance is semiconcave on one side of the interface and semiconvex on the
other; on regions containing both kinds of medial ridge it is more safely
treated as a locally \(DC\) function.  Its gradient is still \(BV\) in the
polyhedral and positive-reach models used below.  Formula (1.1) is the
appropriate first-order Gauss-jump measure in those models.

There is no operation

\[
                         |D u|^2                         \tag{1.2}
\]

on a singular Radon measure.  If a jump of size \(A\) is smoothed through a
layer of width \(\varepsilon\), the unweighted integral of
\(|\nabla u|^2\) is ordinarily of order \(A^2/\varepsilon\).  At a medial
endpoint the transport density vanishes linearly, and the weighted energy
has critical weight \(|s|\).  Section 4 shows that this critical energy has
zero point capacity, so even its relaxed value does not retain the jump.

For a convex potential \(V\), \(D^2V\) is a positive-semidefinite
matrix-valued measure.  The expression \(D^2V(u,u)\) is canonical when \(u\)
has a continuous trace on the support of the singular part of \(D^2V\).
It need not be canonical when a potential ridge and a Gauss jump coincide.
The ray formulation in Section 2 avoids this product: it uses directly the
scalar convex function \(-\log q_y\), whose second derivative is an
unambiguous nonnegative measure.

In a \(C^3\) eikonal chart with \(V\in C^2\), the weighted Bochner formula is

\[
 \frac12L|\nabla f|^2
 =\|\nabla^2f\|_{HS}^2
  +\langle\nabla f,\nabla Lf\rangle
  +\nabla^2V(\nabla f,\nabla f).                          \tag{1.3}
\]

Since \(|\nabla f|=1\), this gives (0.1).  Replacing the smooth Hessians in
(1.3) by their distributional Hessian measures does not give a
distributional identity: the first term is quadratic in \(D^2f\), and the
product \(u\cdot\nabla Lf\) is undefined at a Gauss jump.  The exact
distributional substitute is the following raywise Stieltjes formula.

## 2. The exact one-dimensional Stieltjes identity

### 2.1 Statement

Let \(-\infty\le a<0<b\le+\infty\), and let

\[
 q:(a,b)\longrightarrow(0,+\infty)
\]

be an integrable log-concave function.  Assume that \(q(0)<\infty\) and

\[
 \int_a^0q(t)dt=\int_0^bq(t)dt.                          \tag{2.1}
\]

Put \(w=-\log q\), let

\[
 \kappa=Dw'=D^2w                                        \tag{2.2}
\]

be its nonnegative distributional second derivative on \((a,b)\), and
define

\[
 W(t)=
 \begin{cases}
  \displaystyle\int_a^tq(s)ds,&a<t<0,\\[5pt]
  \displaystyle\int_t^bq(s)ds,&0<t<b.
 \end{cases}                                             \tag{2.3}
\]

The two values at zero agree by (2.1).  Define the endpoint traces

\[
 \begin{split}
 e^-&=\lim_{t\downarrow a}
       \{q(t)+W(t)w'_+(t)\},\\
 e^+&=\lim_{t\uparrow b}
       \{q(t)-W(t)w'_-(t)\}.
 \end{split}                                             \tag{2.4}
\]

Then both limits exist in \([0,+\infty)\), and

\[
 \boxed{
   \int_{(a,b)}W(t)\,d\kappa(t)+e^-+e^+=2q(0).}          \tag{2.5}
\]

Atoms of \(\kappa\), including a kink of \(V\) at \(t=0\), are included
with their full mass in the Stieltjes integral.

### 2.2 Proof, including endpoint signs

Fix \(a<\alpha<0<\beta<b\) at which the relevant one-sided derivatives are
finite.  Stieltjes integration by parts gives

\[
 \begin{aligned}
 \int_{(\alpha,\beta)}W\,d\kappa
 ={}&W(\beta)w'_-(\beta)-W(\alpha)w'_+(\alpha)\\
   &-\int_\alpha^0q(t)w'(t)dt
     +\int_0^\beta q(t)w'(t)dt .
 \end{aligned}                                           \tag{2.6}
\]

The identity \(q'=-qw'\), valid almost everywhere, turns this into

\[
 \begin{aligned}
 \int_{(\alpha,\beta)}W\,d\kappa
 =2q(0)
  &-\{q(\alpha)+W(\alpha)w'_+(\alpha)\}\\
  &-\{q(\beta)-W(\beta)w'_-(\beta)\}.
 \end{aligned}                                           \tag{2.7}
\]

On the left half-line the function

\[
 A(t)=q(t)+W(t)w'_+(t)
\]

has distributional derivative \(W\kappa\ge0\).  On the right half-line

\[
 B(t)=q(t)-W(t)w'_-(t)
\]

has distributional derivative \(-W\kappa\le0\).  Hence the limits in
(2.4) exist.  They are nonnegative.  For example, if \(w'(t)<0\), convexity
gives, for \(s<t\),

\[
 q(s)\le q(t)\exp[-(-w'(t))(t-s)],
\]

and therefore \(W(t)(-w'(t))\le q(t)\), which proves \(A(t)\ge0\).
The other three sign cases are identical.  Finally let
\(\alpha\downarrow a\) and \(\beta\uparrow b\) in (2.7).  Monotone
convergence on the left proves (2.5).

### 2.3 Endpoint taxonomy

The formula distinguishes three geometries.

1. If an endpoint is regular and \(q\) has a positive finite trace there,
   while (w') has a finite trace, then \(W\to0\) and
   \[
                         e=q(\text{endpoint}).           \tag{2.8}
   \]
   This is the case at a regular medial sheet and at a transverse hard
   support.

2. If a focal endpoint is at (s=0) and
   \(q(s)=c s^m(1+o(1))\) with \(m\ge1\), then
   \[
    W(s)={c\over m+1}s^{m+1}(1+o(1)),\qquad
    w'(s)=-{m\over s}+O(1),
   \]
   and (e=0).  Moreover
   \[
        W(s)\,d\kappa(s)=O(s^{m-1})ds,                  \tag{2.9}
   \]
   which is integrable.  Focal curvature is therefore already contained in
   the interior term and creates no extra atom.

3. At an infinite endpoint of an integrable log-concave density, the trace
   in (2.4) is zero.  This follows either from the preceding monotonicity or
   from the exponential-tail bound for one-dimensional log-concave
   functions.

## 3. From the Stieltjes identity to a global mixed curvature measure

Let \(\Sigma=\partial E\) be a regular part of a signed-distance interface.
Use the inward orientation, so \(t=f\), and disintegrate its normal Voronoi
cells as

\[
 d\mu(F(y,t))=q_y(t)dt\,d\eta(y),\qquad
 F(y,t)=y+t u(y).                                        \tag{3.1}
\]

Here \(q_y\) is allowed to be unnormalized.  Assume the cellwise Euler
equation

\[
             \int_{a_y}^0q_y(t)dt=\int_0^{b_y}q_y(t)dt  \tag{3.2}
\]

for almost every \(y\).  This is the balanced-normal-cell condition forced
by a smooth half-mass signed-distance maximizer.  Set

\[
 \tau(F(y,t))={W_y(t)\over q_y(t)}.                      \tag{3.3}
\]

The compactified ray space is the disjoint union

\[
 \widehat{\mathcal R}=\{(y,t):a_y\le t\le b_y\}.
\]

On it define

\[
 \begin{aligned}
 \mathcal C(A)=\int\bigg[&\int_{(a_y,b_y)}
    1_A(y,t)W_y(t)d\kappa_y(t)\\
   &+e_y^-1_A(y,a_y)+e_y^+1_A(y,b_y)\bigg]d\eta(y).
 \end{aligned}                                           \tag{3.4}
\]

The usual measurable-kernel approximation defines (3.4) if the ray
disintegration is only measurable.  Applying (2.5) cell by cell gives

\[
 \mathcal C(\widehat{\mathcal R})
 =2\int q_y(0)d\eta(y)=2P_\mu(E).                        \tag{3.5}
\]

The last equality is the area formula at \(t=0\).  Singular strata of the
interface itself are of codimension at least two in the standard regularity
setting and do not alter this perimeter identity.

### 3.1 Agreement with the smooth Bochner integrand

Suppose a chart is \(C^3\), \(V\in C^2\), and

\[
 q_y(t)=e^{-V(y+t u(y))}\det(I+tS_y),                    \tag{3.6}
\]

with the quotient factor independent of \(t\).  Before the first focal
time, all factors in the determinant are positive.  Direct differentiation
gives

\[
 {d^2\over dt^2}[-\log q_y(t)]
 =\nabla^2V(u,u)
  +\operatorname{tr}\bigl[S_y^2(I+tS_y)^{-2}\bigr].     \tag{3.7}
\]

The second term is \(\|\nabla^2f\|_{HS}^2\) at \(F(y,t)\).  Since
\(d\mu=q_y(t)dt\,d\eta(y)\), equations (3.3) and (3.7) prove (0.5).
If \(V\) is merely convex along the ray, its one-dimensional Hessian
measure is included in \(\kappa_y\) without choosing a representative of
\(u\) on a potential ridge.

### 3.2 Regular medial sheets give linear BV turning

Let \(M\) be a \(C^1\) piece of the medial set at which exactly two regular
ray charts terminate.  Denote their endpoint maps by

\[
 G_i(y)=y+b_i(y)u_i(y)\in M,\qquad i=1,2.                \tag{3.8}
\]

Assume \(q_i\) and the ray logarithmic slopes have finite positive traces.
Then (2.8) gives \(e_i=q_i(b_i)\).  If

\[
 J_i(y,t)=\det(I+tS_i(y)),
\]

the area formula for (3.8) is

\[
 J_i(y,b_i(y))dA_y
 =|u_i\cdot n_M|d\mathcal H^{d-1}_M.                    \tag{3.9}
\]

Indeed, wedging \(DG_i(v_1),\ldots,DG_i(v_{d-1})\) with
\(u_i\) kills all \(db_i\otimes u_i\) terms.  The resulting determinant is
both \(J_i\) and \(|u_i\cdot n_M|\) times the area Jacobian of \(G_i\).

The two smooth eikonal branches have the same trace of \(f\) on \(M\).
Their tangential derivatives therefore agree:

\[
 P_{T M}u_1=P_{T M}u_2.                                 \tag{3.10}
\]

Because both vectors have length one and approach \(M\) from opposite ray
cells, their normal components have opposite signs and equal magnitudes.
Consequently

\[
 |u_1\cdot n_M|+|u_2\cdot n_M|=|u_1-u_2|.               \tag{3.11}
\]

Multiplying (3.9) by the density trace and summing the two incidences proves

\[
 (G_1)_*(e_1dA_1)+(G_2)_*(e_2dA_2)
 =e^{-V}|u_1-u_2|d\mathcal H^{d-1}_M.                   \tag{3.12}
\]

This is precisely the weighted total variation of the jump part (1.1).
It is linear in the turning angle.  Higher-multiplicity junctions have
codimension at least two generically; a non-generic junction is represented
without ambiguity on the incidence graph by the sum of its endpoint
traces.

### 3.3 The same term is the trace of the switching Hessian

Suppose locally the signed-distance envelope switches between two smooth
branches \(g_1,g_2\), and normal deformations have heights \(h_1,h_2\).
The exact envelope calculation gives the nonnegative switching energy

\[
 \mathcal Q_M(h)=
 \int_M{(h_1-h_2)^2e^{-V}\over|\nabla g_1-\nabla g_2|}
       d\mathcal H^{d-1}.                                \tag{3.13}
\]

For an ambient translation vector \(a\), the normal heights are

\[
 h_i^a=a\cdot u_i.
\]

If \(e_1,\ldots,e_d\) is an orthonormal basis, then pointwise on \(M\),

\[
 \sum_{k=1}^d
 {\bigl((e_k\cdot u_1)-(e_k\cdot u_2)\bigr)^2
       \over|u_1-u_2|}
 =|u_1-u_2|.                                             \tag{3.14}
\]

Therefore

\[
 \boxed{\sum_{k=1}^d\mathcal Q_M(h^{e_k})
       =\int_Me^{-V}|u_1-u_2|d\mathcal H^{d-1}.}          \tag{3.15}
\]

The endpoint curvature charge is exactly the translation-mode trace of
the medial second variation.  It is not an ad hoc correction.

### 3.4 Hard support gives normal escape, not Gauss variation

If a regular ray terminates on a \(C^1\) portion of the boundary of a convex
support \(K\), the same wedge-determinant proof gives

\[
 J(y,b(y))dA_y=|u\cdot n_K|d\mathcal H^{d-1}_{\partial K}. \tag{3.16}
\]

At a finite positive density trace the endpoint contribution is therefore

\[
 e^{-V}|u\cdot n_K|d\mathcal H^{d-1}_{\partial K}.       \tag{3.17}
\]

This is a normal escape flux.  It is sensitive to first-order incidence,
not to the second fundamental form of the support.  The quartic example in
Section 6.3 shows why support curvature cannot replace it.

### 3.5 Scaling

Under the dilation \(x\mapsto\lambda x\), perimeter scales by
\(\lambda^{-1}\), the smooth term in (0.5) scales by
\(\lambda\cdot\lambda^{-2}=\lambda^{-1}\), the medial term in
(3.12) scales as
\(\lambda^{-d}\lambda^{d-1}=\lambda^{-1}\), and the support term has the
same scaling.  The mixed budget is therefore scale invariant in exactly the
same sense as the Cheeger constant.

## 4. Exact no-go for a lower-semicontinuous quadratic jump measure

### 4.1 No factorization of endpoint atoms

On every balanced cell, (W\(a\)=W\(b\)=0).  Hence the canonical ray trace of

\[
                         \tau={W\over q}                 \tag{4.1}
\]

is zero at any endpoint with a positive finite density trace.  Let \(S\) be
a medial or support stratum carrying a positive endpoint measure.  If \(K\)
is a locally finite nonnegative Radon measure and \(\tau\) is given its
canonical zero trace on \(S\), then

\[
                         (\tau K)(S)=0.                  \tag{4.2}
\]

But (3.12) or (3.17) gives \(\mathcal C(S)>0\).  Thus

\[
                         \mathcal C\ne\tau K             \tag{4.3}
\]

for every locally finite \(K\).  Redefining \(\tau\) on the
\(\mu\)-null set \(S\) would make the answer depend on an arbitrary
representative.  Allowing \(K(S)=+\infty\) instead makes
\(0\cdot\infty\) undefined and leaves the Radon category.

The same observation applies to a weighted linear measure
\(\tau|D u|\): since the jump part of \(D u\) is supported where
\(\tau=0\), this product also discards all medial turning.  The correct
linear term is the **unmultiplied** density-weighted jump measure in
(3.12).

### 4.2 A completely explicit Gauss kink

Fix \(0<A<1\), set \(C=\sqrt{1-A^2}\), and on

\[
 Q=(-1,1)_s\times(-1,1)_z
\]

define

\[
 f_0(s,z)=Cz-A|s|.                                      \tag{4.4}
\]

Then \(f_0\) is concave and 1-Lipschitz,

\[
 |\nabla f_0|=1\quad(s\ne0),
\]

and its Gauss field has the nonzero jump

\[
 u^-=(A,C),\qquad u^+=(-A,C),\qquad |u^+-u^-|=2A.        \tag{4.5}
\]

The function is the local signed-distance branch inside a wedge.  The
critical endpoint behavior of the transport density is \(\tau(s,z)\sim
|s|\).

Choose

\[
 R_j=e^{-j},\qquad \delta_j=e^{-j^2},\qquad
 L_j=\log(R_j/\delta_j)=j^2-j.                           \tag{4.6}
\]

For \(s\ge0\), define

\[
 p_j(s)=
 \begin{cases}
  0,&0\le s\le\delta_j,\\[2pt]
  \displaystyle {A\log(s/\delta_j)\over L_j},
       &\delta_j<s<R_j,\\[6pt]
  A,&s\ge R_j,
 \end{cases}                                             \tag{4.7}
\]

and extend \(p_j\) oddly.  Smooth the two corners inside intervals of
relative width \(o(1)\), preserving \(|p_j|\le A\); the corner contributions
to the estimate below are \(O(A^2/L_j^2)\).  Let

\[
 g_j(s)=\int_0^sp_j(r)dr,\qquad f_j(s,z)=Cz-g_j(s).       \tag{4.8}
\]

Then \(g_j\) is even after adding an irrelevant constant,

\[
 \|f_j-f_0\|_{L^\infty(Q)}\le2AR_j\longrightarrow0,
 \qquad |\nabla f_j|^2=C^2+p_j^2\le1.                   \tag{4.9}
\]

On the logarithmic transition interval,

\[
                         p_j'(s)={A\over L_js}.
\]

Consequently

\[
 \begin{aligned}
 \int_Q|s|\,\|\nabla^2f_j\|_{HS}^2\,ds\,dz
 &=4\int_0^1s|p_j'(s)|^2ds\\
 &={4A^2\over L_j}+O(A^2/L_j^2)\longrightarrow0.
 \end{aligned}                                           \tag{4.10}
\]

The harmless factor (4) is the \(z\)-length times the two \(s\)-sides.
In contrast,

\[
 \int_Q\|\nabla^2f_j\|_{HS}^2\,ds\,dz
 \ge {4A^2\over L_j^2}
       \left({1\over\delta_j}-{1\over R_j}\right)
 \longrightarrow+\infty.                               \tag{4.11}
\]

A single-scale transition \(p_j(s)=As/R_j\) on \(0<s<R_j\), followed by
\(p_j=A\), has

\[
 \int_Q|s|\,\|\nabla^2f_j\|^2\,ds\,dz\longrightarrow2A^2, \tag{4.12}
\]

up to an arbitrarily small corner-smoothing error.  Concentrating the
transition in a layer of width \(o(R_j)\) centered at \(R_j\) makes the same
energy diverge.  Thus the three possible limiting values include zero, a
positive constant, and \(+\infty\), for smooth 1-Lipschitz approximations of
the identical eikonal kink.

The construction persists with a smooth positive density and with any
weight \(\tau(s,z)=c(z)|s|+o(|s|)\), \(c(z)>0\).  It proves:

> **Quadratic-jump no-go.**  The lower-semicontinuous relaxation, under
> locally uniform convergence of 1-Lipschitz potentials, of the smooth
> energy \(\int\tau\|\nabla^2f\|^2d\mu\) assigns zero energy to a medial
> Gauss jump.  Hence no lower-semicontinuous extension which agrees with
> the smooth quadratic integrand can also assign the positive endpoint
> charge (3.12).

Exact smooth eikonal functions are not a closed regularization class across
a ray collision: in one dimension a \(C^2\) solution of \(|f'|=1\) has
(f''=0).  Restricting approximations to globally smooth exact eikonal
functions therefore does not repair the problem; it removes every sequence
capable of resolving the singularity and supplies no definition of the
limit charge.

### 4.3 Relation to linear turning and to switching

The unweighted linear BV energy of the limit kink is

\[
 \int_{\{s=0\}}|u^+-u^-|\,dz=4A.                        \tag{4.13}
\]

It is lower semicontinuous as part of the full BV total variation.  But the
smooth budget (0.5) contains the quadratic weighted energy, not the full
linear variation of smooth charts.  Adding \(\int|\nabla u|\) in smooth
regions would violate smooth agreement and is not bounded by (0.4).

By contrast, the medial switching form automatically changes homogeneity at
the endpoint.  Its trace is linear by (3.15), and it is exactly the endpoint
part of the ray identity.  Therefore the only scale-correct retained
singular quantity found here is

\[
 e^{-V}|D^j u|,
\]

or the equivalent switching Hessian on the completed ray graph.  It cannot
be written as \(\tau\) times a locally finite curvature measure.

## 5. Convex supports and smooth convex-potential approximation

Let \(K\subset\mathbb R^d\) be a bounded convex body.  A standard soft-wall
approximation is

\[
 V_k(x)=k\,d(x,K)^2,\qquad
 d\mu_k=Z_k^{-1}e^{-V_k(x)}dx.                            \tag{5.1}
\]

The squared distance to a closed convex set is convex and \(C^{1,1}\).
Convolution with a symmetric smooth mollifier gives smooth convex
\(V_{k,\varepsilon}\); first let \(\varepsilon\downarrow0\), then
\(k\to\infty\).  The densities converge in \(L^1\) to
\(1_K/|K|\).  Their second moments converge as well: the exterior tubular
volume has polynomial growth, while \(e^{-k d^2}\) gives a Gaussian bound
in the distance to \(K\).

Suppose \(K\) and the mollifier are invariant under \(x_1\mapsto-x_1\), and
take \(E=\{x_1>0\}\), \(f=x_1\).  The half-mass condition is exact for every
approximation.  The smooth identity reads

\[
 \int\tau_{k,\varepsilon}\,
       \partial_{11}V_{k,\varepsilon}\,d\mu_{k,\varepsilon}
 =2P_{\mu_{k,\varepsilon}}(E).                           \tag{5.2}
\]

After the two limits, the right side converges to the perimeter of the
hard-support cut.  The left side converges to the support endpoint flux in
(3.17), not to \(\tau K\) for a finite weak limit \(K\).

This failure is already exact in one dimension.  For

\[
 V_k(t)=k(|t|-a)_+^2,\qquad q_k(t)=Z_k^{-1}e^{-V_k(t)},   \tag{5.3}
\]

one has

\[
 \int\tau_kV_k''q_kdt=2q_k(0),                          \tag{5.4}
\]

whereas

\[
 \int V_k''q_kdt
 =2\sqrt{\pi k}\,q_k(0)\longrightarrow+\infty.          \tag{5.5}
\]

At the wall, \(\tau_k\asymp k^{-1/2}\).  Thus a diverging unweighted
curvature layer times a vanishing transport density has a finite endpoint
limit.  No subsequence of the unweighted curvature measures is locally
finite near the wall.  Smoothing the two corners of (5.3) leaves (5.4) and
the asymptotic (5.5) unchanged.

The one-dimensional Stieltjes identity also gives the precise convergence
mechanism.  On compact subsets of an open limiting ray, convex convergence
of \(-\log q_k\) implies weak convergence of its Hessian measures and
uniform convergence of \(W_k\).  The mass not retained in the open ray is,
by (2.5), exactly the endpoint defect (2.4).  This proves convergence of the
**weighted completed measures** on a fixed compactified ray.  It does not
prove convergence of an unweighted ambient \(K\), which (5.5) rules out.

## 6. Mandatory exact models

### 6.1 A one-dimensional multi-component balanced set

Let \(\mu\) be uniform on ([-4,4]) and

\[
 E=(-3,-1)\cup(1,3).                                    \tag{6.1}
\]

Then \(\mu(E)=1/2\).  Its four boundary points have normal Voronoi cells

\[
 [-4,-2],\quad[-2,0],\quad[0,2],\quad[2,4].              \tag{6.2}
\]

Orient every cell inward to \(E\).  In signed normal coordinate each cell
is ([-1,1]),

\[
 q(t)={1\over8},\qquad W(t)={1-|t|\over8},\qquad
 \tau(t)=1-|t|.                                         \tag{6.3}
\]

Thus every cell has exactly one sign crossing and is balanced.  The regular
curvature measure \(\kappa\) is zero.  Each of the eight cell endpoints has
charge (1/8), so

\[
 \mathcal C=1=2P_\mu(E),\qquad P_\mu(E)={4\over8}={1\over2}. \tag{6.4}
\]

At the three internal medial points (-2,0,2), the two incident endpoint
charges sum to (1/4).  The Gauss field jumps by (2), and

\[
 {1\over8}|u^+-u^-|={1\over4},                           \tag{6.5}
\]

in agreement with (3.12).  The remaining two charges occur at the support
endpoints (-4,4).

This example verifies the monotone-branch restriction.  A global oriented
line can cross \(E\) several times and its cumulative flux need not remain
nonnegative.  The signed-distance transport instead decomposes the line
into the four normal cells; on each one \(h=\operatorname{sign}t\) and the
transport is monotone.  Multiple connected components are allowed, but
multiple sign crossings on one signed-distance cell are not.

### 6.2 The isotropic cube halfspace

Let \(K=[-a,a]^d\), \(a=\sqrt3\), let \(\mu\) be uniform on \(K\), and take
\(E=\{x_1>0\}\).  For every basepoint
\(y\in[-a,a]^{d-1}\),

\[
 q_y(t)=(2a)^{-d},\quad -a<t<a,\qquad
 W_y(t)=q_y(a-|t|).                                     \tag{6.6}
\]

The interior curvature is zero.  The two endpoint traces are both \(q_y\),
and hence

\[
 \mathcal C=2(2a)^{-d}(2a)^{d-1}={1\over a}=2P_\mu(E).  \tag{6.7}
\]

All curvature is support escape.  The soft-wall computation (5.3)--(5.5)
shows rigorously that the weighted measures converge to these wall atoms
while the unweighted Hessian measures diverge.

### 6.3 Quartic contact

Let

\[
 K=\{(x_1,x'):\ x_1^4+|x'|^2\le1\},                    \tag{6.8}
\]

let \(\mu\) be uniform on \(K\), and take \(E=\{x_1>0\}\).
For \(y=x'\in B_2^{d-1}\), put

\[
 a(y)=(1-|y|^2)^{1/4}.
\]

The normal cell is \((-a(y),a(y))\), \(q_y\) is the constant volume
density, and both endpoint traces equal \(q_y\).  Therefore

\[
 \mathcal C=2q\,|B_2^{d-1}|=2P_\mu(E).                  \tag{6.9}
\]

For the defining function \(F=x_1^4+|x'|^2-1\), at the upper endpoint
\(x_1=a(y)\),

\[
 n_K\cdot e_1={4a(y)^3\over
      \sqrt{16a(y)^6+4|y|^2}},\qquad
 II_{\partial K}(e_1,e_1)={12a(y)^2\over
      \sqrt{16a(y)^6+4|y|^2}}.                          \tag{6.10}
\]

Both quantities tend to zero at the quartic contact
\(a(y)\downarrow0\).  Nevertheless the projected endpoint charge is the
constant measure \(q\,dy\), because

\[
 dy=|n_K\cdot e_1|d\mathcal H^{d-1}_{\partial K}.        \tag{6.11}
\]

Thus zero support curvature in the transported direction does not remove
the endpoint budget.  The correct term is first-order escape flux, not a
support second fundamental form.  Soft potentials \(k\,d(x,K)^2\), followed
by symmetric mollification, preserve the balanced halfspace and converge to
(6.9) by (5.2).  On every transverse subpatch the associated unweighted
Hessian mass diverges as in (5.5).

### 6.4 Polyhedral checkerboard and ridges

On \(K=[-a,a]^d\), define

\[
 E=\{x:\prod_{i=1}^d x_i>0\},\qquad
 f(x)=\operatorname{sgn}\!\left(\prod_i x_i\right)
                      \min_i|x_i|.                      \tag{6.12}
\]

Then \(\mu(E)=1/2\), and \(f\) is its signed distance.  Away from ties, a
normal cell is indexed by a coordinate \(i\) and a basepoint
\(y\in\{x_i=0\}\).  If

\[
 m(y)=\min_{j\ne i}|y_j|,
\]

then the cell is \((-m(y),m(y))\), its density is constant, and it is
balanced.  The classical Hessian vanishes on every open cell.  All internal
curvature is on the tie ridges \(|x_i|=|x_j|\), where the Gauss field jumps.

For each of the \(d\) coordinate hyperplanes,

\[
 \int q_y(0)dA_y={1\over2a}.
\]

Thus

\[
 P_\mu(E)={d\over2a},\qquad \mathcal C={d\over a}=2P_\mu(E). \tag{6.13}
\]

At a generic tie ridge exactly two affine branches meet.  Formula (3.12)
identifies the two endpoint traces with
\(q|u^+-u^-|d\mathcal H^{d-1}\).  The trace of the switching form gives the
same quantity by (3.15).  Since \(\tau=0\) on every tie ridge,
\(\tau|D^ju|=0\); this is an exact counterexample to weighting the linear
turning measure by \(\tau\).

### 6.5 Radial exponential and the focal center

For \(d\ge2\), let

\[
 d\mu=c_de^{-\lambda|x|}dx
\]

and let \(r_0\) be a median radius.  Take
\(f(x)=|x|-r_0\).  With angular quotient measure, the radial conditional
density is

\[
 q(r)=c_dr^{d-1}e^{-\lambda r},\qquad 0<r<\infty.        \tag{6.14}
\]

Its logarithmic curvature is

\[
 D^2(-\log q)={d-1\over r^2}dr.                          \tag{6.15}
\]

At the focal endpoint \(r=0\),

\[
 q(r)\sim c_dr^{d-1},\quad
 W(r)\sim{c_d\over d}r^d,\quad
 q(r)+W(r)\left(\lambda-{d-1\over r}\right)
 \sim {c_d\over d}r^{d-1}\to0.                         \tag{6.16}
\]

There is no center atom.  Formula (2.5) gives

\[
 \int\tau{d-1\over r^2}d\mu=2P_\mu(E).                 \tag{6.17}
\]

Thus a codimension-\(d\) focal collapse is accounted for by integrable
smooth curvature and must not be assigned a hypersurface Gauss-jump charge.
For (d=1), the center separates two one-dimensional incidences with
positive density trace and is a medial endpoint instead; the endpoint term
then applies.

### 6.6 A smooth-density fan

Let \(\rho(x)=\varphi(|x|)>0\) be a smooth radial density on
\(\mathbb R^2\), and let the interface be the union of \(m\ge2\) lines
through the origin at angles \(j\pi/m\).  Color the (2m) sectors
alternately and let \(f\) be the signed distance.  Reflection across each
interface line exchanges the two colors and preserves \(\rho\), so every
normal cell is balanced.

Write \(\alpha=\pi/m\), parametrize one interface line by
(y=r e), and let \(n\) be its normal.  Its normal cell is

\[
       -\ell(r)<t<\ell(r),\qquad
       \ell(r)=|r|\tan(\alpha/2),                        \tag{6.18}
\]

with

\[
 q_r(t)=\varphi(\sqrt{r^2+t^2}).                         \tag{6.19}
\]

For the standard Gaussian, \(D^2(-\log q_r)=dt\).  Hence the exact cell
identity is

\[
 \int_{-\ell(r)}^{\ell(r)}W_r(t)dt
       +2q_r(\ell(r))=2q_r(0).                           \tag{6.20}
\]

The two endpoint terms occur on the adjacent angular bisectors, where the
Gauss field jumps.  Their pushforward is

\[
 \rho|u^+-u^-|d\mathcal H^1,                             \tag{6.21}
\]

and equals the translation trace of the fan switching form.  The density is
smooth; the singularity belongs entirely to the eikonal foliation.  This
model rules out the hope that smoothing \(V\) alone creates an ambient
quadratic curvature measure.

## 7. Consequence for the curvature route

For balanced signed-distance cells with regular codimension-one endpoints,
the exact budget can be summarized schematically as

\[
 \begin{aligned}
 &\int_{\mathrm{regular}}
   \tau\bigl(\|\nabla^2f\|_{HS}^2+\nabla^2V(u,u)\bigr)d\mu\\
 &\quad+
   \int_{J_u}e^{-V}|u^+-u^-|d\mathcal H^{d-1}
   +\int_{\partial K}e^{-V}|u\cdot n_K|d\mathcal H^{d-1}
 \le 2P_\mu(E),
 \end{aligned}                                           \tag{7.1}
\]

with equality after all ray incidences, focal contributions, and singular
potential terms are represented on the compactified ray graph.  The first
line is quadratic and transport weighted.  The second line is linear and
unweighted by \(\tau\).  The change of homogeneity is forced by the exact
Stieltjes identity.

The following statements are therefore established.

1. A finite, nonnegative completed **weighted** curvature measure exists on
   the normal-ray graph and has total mass \(2P_\mu(E)\).
2. Its regular part agrees with the smooth Bochner integrand.
3. Its generic medial part is the linear BV Gauss-jump measure and is the
   trace of the signed-distance switching Hessian.
4. Its hard-support part is normal escape flux.  Focal endpoints with
   vanishing Jacobian need no atom.
5. It cannot be factored as \(\tau K\) with \(K\) locally finite, and no
   lower-semicontinuous quadratic relaxation retains its Gauss jumps.

Thus the proposed nonsmooth quadratic \(K\) is refuted.  A proof strategy
may still use (7.1), but it must exploit a mixed quadratic/BV functional or
the full switching graph.  Any argument that writes “singular turning” as a
quadratic Hessian measure, or weights the BV jump measure by \(\tau\), loses
exactly the polyhedral and fan geometries that the singular completion was
supposed to control.

## 8. Quantified consequence on a long-cell core, and the remaining rigidity theorem

### 8.1 What follows immediately

Use the normalization of the transport-density tensor:

\[
 D=\int\tau\,d\mu,\qquad
 M=\int\tau\,u\otimes u\,d\mu,\qquad Q={M\over D}.        \tag{8.1}
\]

The exact tensor identity gives

\[
 \operatorname{tr}M=D,\qquad \|M\|_{HS}\le1,
\]

and hence

\[
 \operatorname{tr}Q=1,\qquad \|Q\|_{HS}\le {1\over D},
 \qquad \operatorname{rank}Q\ge D^2.                    \tag{8.2}
\]

Write the three parts of the completed curvature budget as

\[
 \begin{aligned}
 \mathfrak I&=\int_{\rm reg}
   \tau\bigl(\|\nabla^2f\|_{HS}^2+\nabla^2V(u,u)\bigr)d\mu,\\
 \mathfrak J&=\int_{J_u}e^{-V}|u^+-u^-|d\mathcal H^{d-1},\\
 \mathfrak S&=\int_{\partial K}e^{-V}|u\cdot n_K|
                                      d\mathcal H^{d-1}.
 \end{aligned}                                           \tag{8.3}
\]

with incidence multiplicity when needed.  Equations (3.5) and (7.1) give

\[
                         \mathfrak I+\mathfrak J+\mathfrak S=2P. \tag{8.4}
\]

One must not infer \(P\le A/D\) for this same interface from Milman's
equivalence.  Milman's theorem compares the optimal first-moment and optimal
isoperimetric scales, which may be witnessed by different functions and
sets.  Indeed, for a cell mixture, \(D\asymp\mathbb E s_y\) and
\(P\asymp\mathbb E s_y^{-1}\), which gives the lower bound \(PD\gtrsim1\)
but no upper bound.  A same-interface estimate \(P\lesssim D^{-1}\) is an
additional conjecture-strength hypothesis unless it is proved from
extremality or stability.

There is nevertheless a rigorous fixed-mass long-cell extraction.  Normalize
the quotient so that \(\eta\) is a probability and every \(q_y\) is a
conditional probability density.  Put

\[
 d_y=\int|t|q_y(t)dt,\qquad s_y={1\over q_y(0)},\qquad
 D=\int d_y\,d\eta(y).                                   \tag{8.5}
\]

For a balanced one-dimensional log-concave density, the sharp elementary
cell comparisons give

\[
 c_0s_y\le d_y\le C_0s_y,\qquad
 c_0s_y^2\le\int t^2q_y(t)dt\le C_0s_y^2.                \tag{8.6}
\]

Assume that the signed distance has value at least a fixed fraction of
\(D_1(\mu)\).  Milman's Poincare comparison, applied to this 1-Lipschitz
function, gives

\[
 \int s_y^2d\eta(y)\le C D^2,\qquad
 \int s_y d\eta(y)\asymp D.                              \tag{8.7}
\]

Here the possible nonzero mean of \(f\) costs at most
\((\int|f|d\mu)^2=D^2\).  Choose universal \(0<a<b<\infty\) so that

\[
 \mathcal G=\{y:aD\le s_y\le bD\}
\]

satisfies

\[
 \eta(\mathcal G)\ge\beta,\qquad
 D_{\mathcal G}:=\int_{\mathcal G}d_y\,d\eta(y)\ge\beta D \tag{8.8}
\]

for a universal \(\beta>0\).  To verify this, the lower tail contributes at
most \(aD\), while Cauchy--Schwarz and (8.7) make the expectation above
\(bD\) at most \(CD/b\); then choose \(a\) small and \(b\) large.

Apply the cell identity (2.5) only on \(\mathcal G\).  Since
\(q_y(0)=s_y^{-1}\),

\[
 \boxed{\mathcal C_{\mathcal G}
 =2\int_{\mathcal G}{d\eta(y)\over s_y}\le {2\over aD}.} \tag{8.9}
\]

This is the valid replacement for the unjustified whole-interface bound.
It includes regular turning, medial switching, and support escape on the
selected rays.  Put

\[
 M_{\mathcal G}=\int_{\{y\in\mathcal G\}}\tau u\otimes u\,d\mu,
 \qquad Q_{\mathcal G}={M_{\mathcal G}\over D_{\mathcal G}}.
\]

Since \(0\preceq M_{\mathcal G}\preceq M\),

\[
 \operatorname{tr}Q_{\mathcal G}=1,\qquad
 \|Q_{\mathcal G}\|_{HS}\le {C\over D},\qquad
 \operatorname{rank}Q_{\mathcal G}\ge cD^2.             \tag{8.10}
\]

With \(d\nu_{\mathcal G}=\tau1_{\mathcal G}d\mu/D_{\mathcal G}\),
the normalized regular turning is \(O(D^{-2})\), while the medial jump
charge is \(O(D^{-1})\).  Finally, conditioning the basepoint second-moment
bound on an event of quotient probability at least \(\beta\) gives

\[
 \operatorname{Cov}_{op}(Y\mid\mathcal G)\le {4\over\beta}. \tag{8.11}
\]

No power of the dimension occurs in (8.8)--(8.11).  Short cells outside
\(\mathcal G\) may carry arbitrarily large interface perimeter; no claim
about their total curvature is being made.

### 8.2 These numerical constraints do not force \(D=O(1)\)

There is a sharp algebraic countermodel.  For an integer \(r\), let the
quotient choose \(i\in\{1,\ldots,r\}\) uniformly, take \(u_i=e_i\), put
\(Y_i=0\), and let \(T_i\) be uniform on

\[
                         [-L,L],\qquad L=\sqrt{3r}.       \tag{8.12}
\]

The resulting abstract ray law has covariance

\[
 {1\over r}\sum_{i=1}^r\mathbb E[T_i^2]e_i\otimes e_i=I.
                                                               \tag{8.13}
\]

Each cell is balanced at zero.  Its signed-distance moment and perimeter
scale are

\[
 D=\mathbb E|T|={L\over2}={\sqrt{3r}\over2},\qquad
 P={1\over2L}.                                           \tag{8.14}
\]

Moreover

\[
 Q={I_r\over r},\qquad
 \|Q\|_{HS}=r^{-1/2},\qquad
 \|M\|_{HS}={D\over\sqrt r}={\sqrt3\over2},              \tag{8.15}
\]

the basepoint covariance is zero, the interior and medial energies vanish,
and the two endpoint escape charges give

\[
 \mathfrak S=2P={1\over L}=O(D^{-1}).                    \tag{8.16}
\]

Thus every core-scale numerical line in (8.2) and (8.8)--(8.11) holds while
\(D\to\infty\).  This ray law is not the disintegration of a log-concave
probability: its support is a union of coordinate segments and is not
convex.  Nor is it realized by the normal congruence of one codimension-one
interface: the abstract rays are superposed without satisfying global
nonbranching/eikonal compatibility.  That failure is the point.  The tensor, curvature,
and basepoint-covariance estimates by themselves do not encode global
log-concavity or compatibility between spatially distinct direction
packets.

The same countermodel shows that retaining the support term does not repair
the algebraic implication.  A proof must use a genuinely global
log-concave incidence theorem, not another summation of the existing
budgets.

### 8.3 A sharp formal statement

The missing assertion can be isolated without referring to a quadratic
singular measure.  It is an explicit **unproved target**, not a consequence
of Sections 2--7.

> **Mixed-BV ray-congruence rigidity (MBR).**  For every
> \(A<\infty\) and \(\beta>0\) there is a universal
> \(C(A,\beta)<\infty\) with the following property.  Let \(\mu\)
> be isotropic and log-concave, and let a half-mass signed-distance
> interface admit a balanced normal-cell disintegration.  Suppose a
> measurable ray core \(\mathcal G\) satisfies
> \[
> \begin{gathered}
> D_{\mathcal G}\ge\beta D,\qquad
> \operatorname{tr}M_{\mathcal G}=D_{\mathcal G},\qquad
> \|M_{\mathcal G}\|_{HS}\le1,\\
> \operatorname{Cov}_{op}(Y\mid\mathcal G)\le A,\qquad
> \mathcal C_{\mathcal G}\le {A\over D}.
> \end{gathered}                                         \tag{8.17}
> \]
> Then \(D\le C(A,\beta)\).

The support term may instead be handled by a separate contact-tensor
alternative, but it cannot simply be deleted: the cube halfspace has
\(\mathfrak I=\mathfrak J=0\).

For a fixed-factor signed-distance maximizer, (8.17) is automatic from the
long-cell extraction (8.5)--(8.11), without any same-interface estimate on
\(P\).  Therefore MBR proves the first-moment form of KLS.  Conversely,
the first-moment form of KLS makes the conclusion of MBR immediate.
Subject only to the usual attainment/approximation step for the
signed-distance maximizer, MBR is consequently equivalent in strength to
KLS.  It is a clean statement of the remaining global compatibility gap,
not an intermediate theorem supplied by the curvature calculation.

An equivalent quantitative form, useful for a proof attempt, would be

\[
 D\le C_\beta\left[
  1+\|M_{\mathcal G}\|_{HS}
   +D\mathcal C_{\mathcal G}
   +\|\operatorname{Cov}(Y\mid\mathcal G)\|_{op}^{1/2}\right].
                                                               \tag{8.18}
\]

Every term on the right is universally bounded in the extremal branch.
The abstract star (8.12) disproves (8.18) without the hypothesis that the
ray congruence comes from one log-concave ambient measure.

### 8.4 Stress tests for MBR

**Radial exponential.**  For the median sphere,

\[
 Q={I_d\over d},\qquad \|Q\|_{HS}=d^{-1/2}.
\]

There is no medial or support charge, while
\(\mathfrak I=2P\) by (6.17).  Both \(D\) and \(P\) are of universal order.
This is a high-rank equality model.  MBR cannot assert that low
\(\|Q\|_{HS}\) forces near-affinity, nor can it require a strict curvature
gap; it must allow the translated-radial alternative.

**Cube halfspace.**  Here \(Q=e_1\otimes e_1\),
\(\mathfrak I=\mathfrak J=0\), and \(\mathfrak S=2P\).  The moment is
\(D=\sqrt3/2\).  Any MBR proof based only on interior or medial turning
fails this model; the contact/escape branch is mandatory.

**Product maximum.**  Let \(X_i\) be independent symmetric exponentials
with variance one, so \(Y_i=|X_i|\) are exponential of rate
\(\lambda=\sqrt2\).  If \(r\) is the median of
\(R=\max_iY_i\), the Euclidean signed distance to the inner box
\(\{Y_i<r\ \forall i\}\) is

\[
 |f(X)|=
 \begin{cases}
  r-R,&R\le r,\\
  \bigl(\sum_i(Y_i-r)_+^2\bigr)^{1/2},&R>r.
 \end{cases}                                             \tag{8.19}
\]

Since

\[
 (1-e^{-\lambda r})^n={1\over2},\qquad
 ne^{-\lambda r}\asymp1,                                \tag{8.20}
\]

one has

\[
 \mathbb E\!\left[\sum_i(Y_i-r)_+^2\right]
 ={2n e^{-\lambda r}\over\lambda^2}=O(1).               \tag{8.21}
\]

Also

\[
 \mathbb P\{R\le r-t\}
 \le\exp[-c e^{\lambda t}],
\]

so \(\mathbb E(r-R)_+=O(1)\).  Hence

\[
                              \mathbb E|f(X)|=O(1).       \tag{8.22}
\]

Hyperoctahedral symmetry forces its normalized direction tensor, wherever
the signed-distance transport is defined, to be a scalar multiple of the
identity.  Thus diffuse coordinate directions are compatible with bounded
first moment.  The product maximum is a stress test against any rigidity
claim that treats high rank alone as exceptional.  This median box is not
asserted to satisfy the cellwise Euler equation of a global maximizer.

**Polyhedral fan.**  In the radial fan of Section 6.6,

\[
 |f(x)|\le |x|\sin{\pi\over2m},
\]

so \(D=O(m^{-1})\) for an isotropic radial law with
\(\mathbb E|X|=O(1)\) in dimension two.  Its interface perimeter is
\(P\asymp m\), hence \(PD\asymp1\).  The direction tensor is
\(Q=I_2/2\), and the linear medial jump term is of order \(P\), not of order
the square of the sector angle.  The fan therefore saturates the
scale-invariant switching budget and confirms that replacing
\(\mathfrak J\) by an angle-squared atom would give the wrong scale.

These models leave exactly one viable use of (8.6)--(8.9): prove MBR by a
new global theorem that couples direction packets through log-concavity.
The mixed curvature identity supplies the correct functional and removes
the nonsmooth bookkeeping gap, but it does not supply that global theorem.
