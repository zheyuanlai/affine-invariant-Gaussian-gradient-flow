# Quotient second variation for a \(D_1\)-maximizing potential

## 0. Outcome

For a probability \(\mu\) with finite first moment, put

\[
 {\cal D}_\mu(f)=\int|f-\mu f|\,d\mu,\qquad
 D_1(\mu)=\sup_{\operatorname {Lip}(f)\le1}{\cal D}_\mu(f). \tag{0.1}
\]

Suppose a maximizer \(f\), centered by \(\mu f=0\), admits a regular
transport-ray foliation

\[
 x=\Phi(y,t)=y+t\,u(y),\qquad f(\Phi(y,t))=t,                \tag{0.2}
\]

and a probability disintegration

\[
 d\mu(\Phi(y,t))=q_y(t)\,dt\,d\nu(y).                       \tag{0.3}
\]

Write

\[
 \beta=\mu(\operatorname {sign}f),\qquad
 s_y=\int\operatorname {sign}(t)q_y(t)\,dt.                 \tag{0.4}
\]

A globally feasible normal-graph perturbation by a ray-constant function
\(h(y)\) gives the first-order balance condition

\[
                              s_y=\beta                     \tag{0.5}
\]

for every quotient point seen by the admissible perturbation class.
Thus “balanced rays” means \(\beta=0\) and each conditional measure puts
mass \(1/2\) on each side of \(t=0\).

The genuine second-order eikonal correction is

\[
 k(y,t)=-{1\over2}\int_0^t
              |\nabla H(\Phi(y,r))|^2\,dr,\qquad H(\Phi(y,t))=h(y).
                                                                  \tag{0.6}
\]

Define the longitudinal tail weight

\[
 W_y(r)=
 \begin{cases}
 (1-\beta)\displaystyle\int_r^\infty q_y(t)\,dt,&r>0,\\[5pt]
 (1+\beta)\displaystyle\int_{-\infty}^r q_y(t)\,dt,&r<0.
 \end{cases}                                                   \tag{0.7}
\]

Under the explicit no-switching hypotheses of Theorem 3.1 below,
maximality gives the quotient stability inequality

\[
 \boxed{\quad
 \int q_y(0)(h-\nu h)^2\,d\nu(y)
 \le {1\over2}\int\!\!\int W_y(t)
       |\nabla H(\Phi(y,t))|^2\,dt\,d\nu(y).
 \quad}                                                       \tag{0.8}
\]

For parallel rays, \(|\nabla H|\) is independent of \(t\).  If the
conditional ray coordinate has mean zero, (0.8) becomes

\[
 \int q_y(0)(h-\nu h)^2\,d\nu
 \le {1\over2}\int
       \big(\mathbb E_y|T|-\beta\mathbb E_yT\big)
       |\nabla h|^2\,d\nu.                                  \tag{0.9}
\]

In the balanced centered case this is simply

\[
 \int q_y(0)(h-\nu h)^2\,d\nu
 \le {1\over2}\int\mathbb E_y|T|\,|\nabla h|^2\,d\nu.       \tag{0.10}
\]

This is the exact appearance of the \(q_y(0)h^2\) charge and the
longitudinal Dirichlet cost.

For an arbitrary possibly branching ray foliation, the formal correction
(0.6) is not by itself admissible.  A signed-distance perturbation is
globally \(1\)-Lipschitz but may create a nonnegative medial/contact
switching cost on the right side of (0.8).  Without controlling that cost,
the universally rigorous fallback is

\[
 \int q_y(0)(h-\nu h)^2\,d\nu
 \le {D_1(\mu)\over2}\operatorname {Lip}(H)^2.              \tag{0.11}
\]

It follows by normalizing \(f+\varepsilon H\) by
\(\sqrt{1+\varepsilon^2\operatorname {Lip}(H)^2}\); no formal eikonal
expansion is used.

The matrix tests \(h=a\cdot u\) and \(h=a\cdot y\) yield valid weighted
matrix inequalities, recorded in Section 5.  They do **not** by themselves
give \(\mathbb E_\nu(1/\sigma_y)\ge c\).  Isotropy only gives the
unweighted conditional covariance constraint

\[
 \int\sigma_y^2u(y)\otimes u(y)\,d\nu(y)\preceq I,          \tag{0.12}
\]

whereas the left sides of the stability inequalities contain
\(q_y(0)\asymp1/\sigma_y\) multiplied by centered position or normal
moments.  The needed correlation is absent.  The right sides contain
exactly the curvature/metric energy that would have to be bounded.

The model audit is sharp:

* a Gaussian halfspace makes (0.10) the Gaussian Poincaré inequality,
  with equality for \(h=a\cdot y\);
* the half-cube makes it a quotient-cube Poincaré inequality, and
  checkerboard modes pay their full frequency;
* for a radial sphere, \(h=a\cdot u\) gives the exact necessary condition
  \(2q(0)\le(n-1)\int W(t)/(r_0+t)^2\,dt\).  Gaussian radial distance
  violates it asymptotically by a factor \(2\), whereas the isotropic
  radial-exponential model is at exact equality in every \(n\ge2\).

Thus the second variation is informative—it excludes radial distance for
the Gaussian—but summing its matrix tests is sharp rather than
self-improving.  A KLS-closing use would require a new bound on the
longitudinal curvature energy or a proof that medial/contact switching
supplies a finite competitor.

## 1. Maximizers and the balanced ray condition

### 1.1 Existence

Fix a point \(x_0\) in the support and normalize candidate functions by
\(f(x_0)=0\).  Then \(|f(x)|\le|x-x_0|\).  A maximizing sequence is
locally uniformly precompact by Arzelà--Ascoli.  Its locally uniform limit
is \(1\)-Lipschitz, and uniform integrability follows from the common
linear bound and the finite first moment.  Hence the limit attains
\(D_1(\mu)\).  Log-concave probabilities have moments of every order, so
all later fixed-measure Taylor remainders are integrable.

Subtracting \(\mu f\) leaves both the Lipschitz constant and
\({\cal D}_\mu(f)\) unchanged.  We henceforth take

\[
                            \mu f=0.                         \tag{1.1}
\]

### 1.2 Saturated rays

For a \(1\)-Lipschitz \(f\), a transport ray is a maximal line segment on
which

\[
                        f(x)-f(z)=|x-z|.                    \tag{1.2}
\]

Orient it in the direction of increasing \(f\).  On every nondegenerate
ray one may use \(t=f\) as arclength, giving (0.2).  The report assumes
that the part of \(\mu\) under consideration is covered, up to a null set,
by rays that cross \(S=\{f=0\}\), and that the quotient disintegration
(0.3) is normalized by \(\int q_y=1\).

If \(x=\Phi(y,t)\) lies on such a ray, then

\[
 \operatorname {dist}(x,S)=|t|.                             \tag{1.3}
\]

Indeed the point \(y\in S\) on the same ray gives
\(\operatorname {dist}(x,S)\le|t|\), while the Lipschitz inequality
\(|f(x)-f(z)|\le|x-z|\) for every \(z\in S\) gives the reverse inequality.
Thus \(f\) is the signed distance to its zero set on the ray-covered
region.

### 1.3 First variation

Let \(H(\Phi(y,t))=h(y)\) be a globally feasible ray-constant first
variation.  Center \(h\) by \(\nu h=0\).  If \(q_y\) has no atom at zero,

\[
 {d\over d\varepsilon}\bigg|_{\varepsilon=0}
 {\cal D}_\mu(f+\varepsilon H)
 =\int s_yh(y)\,d\nu(y).                                    \tag{1.4}
\]

Centering the perturbed function subtracts its mean.  Equivalently, for
an arbitrary \(h\), the right side is

\[
                    \int(s_y-\beta)(h-\nu h)\,d\nu.         \tag{1.5}
\]

If both signs of \(\varepsilon\) are feasible and the admissible quotient
functions are dense, maximality forces (0.5).  This deduction is only as
strong as the global perturbation class; no balance is inferred on
quotient components that cannot be varied independently.

## 2. The exact scalar Taylor expansion

Let

\[
 f_\varepsilon=f+\varepsilon H+\varepsilon^2K+o(\varepsilon^2)
                                                                  \tag{2.1}
\]

in \(L^1(\mu)\), and assume the remainder is uniformly integrable after
division by \(\varepsilon^2\).  Put

\[
 h_0=h-\nu h,\qquad K_0=K-\mu K.                            \tag{2.2}
\]

For a continuous one-dimensional density \(q\) at zero,

\[
\begin{split}
 \int|t+\varepsilon a+\varepsilon^2b(t)|q(t)\,dt
 ={}&\int|t|q(t)\,dt
  +\varepsilon a\int\operatorname {sign}(t)q(t)\,dt\\
 &+\varepsilon^2\left(
      q(0)a^2+\int\operatorname {sign}(t)b(t)q(t)\,dt\right)
  +o(\varepsilon^2).
\end{split}                                                   \tag{2.3}
\]

To prove (2.3), split the integral at
\(|t|\le C|\varepsilon|\).  Away from zero, \(|\cdot|\) is affine and
the expansion is exact to second order.  In the shrinking interval,
substitute \(t=\varepsilon r\); continuity of \(q\) gives the term
\(q(0)a^2\).  Truncation and dominated convergence handle an unbounded
\(b\).

Applying (2.3) on the rays and then centering gives

\[
\begin{split}
 {\cal D}_\mu(f_\varepsilon)
 ={}&{\cal D}_\mu(f)
 +\varepsilon\int(s_y-\beta)h_0\,d\nu\\
 &+\varepsilon^2\left\{
   \int q_y(0)h_0(y)^2\,d\nu(y)
   +\int(\operatorname {sign}(t)-\beta)
          K(y,t)q_y(t)\,dt\,d\nu(y)\right\}
 +o(\varepsilon^2).
\end{split}                                                   \tag{2.4}
\]

The additive ray-constant part of \(K\) disappears from the second line
when \(s_y=\beta\).

## 3. Global feasibility

### 3.1 A rigorous normal-graph theorem

Assume in this subsection:

1. \(S=\{f=0\}\) is a \(C^3\) oriented hypersurface;
2. its normal map
   \(\Phi(y,t)=y+tN(y)\) is one-to-one on the ray intervals
   \(I_y=(a_y,b_y)\) and covers the support up to a null set;
3. the ray endpoints and the support boundary create no
   \(O(\varepsilon^2)\) switching mass under the perturbation below;
4. \(h\in C^2(S)\) and its first two derivatives are bounded so that the
   normal graphs have uniform local reach.

Let \({\mathsf A}_y=-D_SN(y)\) be the shape operator, with the convention

\[
 D\Phi(y,t)|_{T_yS}=I-t{\mathsf A}_y.                       \tag{3.1}
\]

Define the perturbed zero set

\[
 S_\varepsilon
 =\{\,y-\varepsilon h(y)N(y):y\in S\,\},                    \tag{3.2}
\]

and let \(f_\varepsilon\) be its signed distance, with the original
orientation.  Signed distance is globally \(1\)-Lipschitz: on the same
side this is the ordinary distance inequality, and for points on opposite
sides the segment joining them crosses \(S_\varepsilon\).

On compact subsets before focal endpoints, the nearest-point equations
and the implicit function theorem give

\[
 f_\varepsilon(\Phi(y,t))
 =t+\varepsilon h(y)+\varepsilon^2k(y,t)+O(\varepsilon^3),  \tag{3.3}
\]

where

\[
\begin{split}
 \Gamma_h(y,t)
 &=|\nabla H(\Phi(y,t))|^2\\
 &=\left|(I-t{\mathsf A}_y)^{-1}\nabla_Sh(y)\right|^2,\\
 k(y,t)&=-{1\over2}\int_0^t\Gamma_h(y,r)\,dr.
\end{split}                                                   \tag{3.4}
\]

One may also derive (3.4) without coordinates.  Expanding the exact
eikonal equation \(|\nabla f_\varepsilon|^2=1\) gives

\[
 \nabla f\cdot\nabla H=0,\qquad
 2\nabla f\cdot\nabla k+|\nabla H|^2=0.                     \tag{3.5}
\]

The first identity says that \(H\) is ray-constant; integrating the second
along a ray, with \(k(y,0)=0\), gives (3.4).

The signed-distance construction, rather than (3.5) alone, proves global
admissibility.  Under hypotheses 1--4, the fixed-measure moment bounds and
the no-switching assumption justify inserting (3.3) into (2.4).

**Theorem 3.1 (regular quotient stability).**  If \(f\) is a global
maximizer and \(s_y=\beta\), then every \(h\) satisfying hypotheses 1--4
obeys

\[
 \int q_y(0)(h-\nu h)^2\,d\nu
 \le {1\over2}\int\!\!\int_{I_y}W_y(t)
       \left|(I-t{\mathsf A}_y)^{-1}\nabla_Sh(y)\right|^2
       \,dt\,d\nu(y).                                       \tag{3.6}
\]

**Proof.**  Insert (3.4) into the \(K\)-term in (2.4).  Fubini gives

\[
\begin{split}
 &\int(\operatorname {sign}(t)-\beta)k(y,t)q_y(t)\,dt\\
 &\hspace{35mm}
 =-{1\over2}\int_{I_y}W_y(r)\Gamma_h(y,r)\,dr.
\end{split}                                                   \tag{3.7}
\]

The first-order term vanishes by balance.  The coefficient of
\(\varepsilon^2\) is nonpositive by global maximality, proving (3.6).
\(\square\)

### 3.2 Cut loci, branching, and contact

If the normal-ray cells meet at a medial set, the signed-distance
perturbation (3.2) remains globally feasible, but (3.3) need not hold to
second order after integration.  Near a tie between two nearest-point
branches, distance is the minimum of their branch distances.  Moving the
branches changes the winning cell in a region of thickness
\(O(|\varepsilon|)\), producing an \(O(\varepsilon^2)\) term.

On a smooth two-cell stratum this term is a nonnegative quadratic
*cost* on the right side of (3.6), proportional to the square of the
difference of the two traces of \(h\), divided by the angle at which the
distance branches meet.  Support endpoints similarly create a contact
trace.  Thus the honest schematic inequality is

\[
 \int q_y(0)(h-\nu h)^2\,d\nu
 \le {1\over2}\int\!\!\int W_y\Gamma_h
       +{\cal C}_{\rm medial}(h)+{\cal C}_{\rm contact}(h),  \tag{3.8}
\]

with nonnegative terms that must be constructed from the actual ray-cell
stratification.  Dropping them is not justified.  A complete formula for
arbitrary nonsmooth transport rays would require a second-variation
theorem for the global signed-distance envelope; it is not assumed here.

This is the quotient analogue of the contact tensor in the companion
tensor-Minkowski report: ray endpoints can carry a
dimension-sized part of a globally feasible variation even though they
are invisible in the interior eikonal equation.

### 3.3 A globally rigorous normalized fallback

There is a weaker construction that does not use a cut-locus expansion.
Let \(H\) be globally Lipschitz on the convex support, ray-constant, and
assume

\[
             \nabla f\cdot\nabla H=0\quad\hbox{a.e.}         \tag{3.9}
\]

Put \(L=\operatorname {Lip}(H)\) and

\[
              g_\varepsilon
       ={f+\varepsilon(H-\mu H)\over\sqrt{1+\varepsilon^2L^2}}.
                                                                  \tag{3.10}
\]

Rademacher's theorem and (3.9) give

\[
 |\nabla(f+\varepsilon H)|^2
 \le1+\varepsilon^2L^2\quad\hbox{a.e.}                       \tag{3.11}
\]

on the convex support.  Hence \(g_\varepsilon\) is \(1\)-Lipschitz there,
and it has a \(1\)-Lipschitz extension to \(\mathbb R^n\).

Using (2.3), balance, and
\((1+\varepsilon^2L^2)^{-1/2}
 =1-\varepsilon^2L^2/2+o(\varepsilon^2)\), maximality yields

\[
 \boxed{\quad
 \int q_y(0)(h-\nu h)^2\,d\nu
       \le {D_1(\mu)\over2}L^2.
 \quad}                                                       \tag{3.12}
\]

Unlike (3.6), (3.12) is valid without a formal second-order ray
correction whenever a global ray-constant \(H\) satisfying (3.9) exists.

## 4. Parallel rays and conditional scale

Suppose \(u(y)=u_0\) is constant, \(S\subset u_0^\perp\), and
\({\mathsf A}=0\).  Then \(\Gamma_h(y,t)=|\nabla_Sh(y)|^2\).  Fubini gives

\[
 \int_{I_y}W_y(t)\,dt
 =(1-\beta)\mathbb E_yT_+
  +(1+\beta)\mathbb E_yT_-
 =\mathbb E_y|T|-\beta\mathbb E_yT.                         \tag{4.1}
\]

This proves (0.9).  If \(\mathbb E_yT=0\), the last expression is
\(\mathbb E_y|T|\), even when \(\beta\ne0\).

If the ray laws are log-concave, balanced at zero, and have standard
deviation \(\sigma_y\), the one-dimensional median-density estimates give

\[
 q_y(0)\ge{1\over4\sqrt3\,\sigma_y},\qquad
 c\,\sigma_y\le\mathbb E_y|T|\le\sigma_y.                   \tag{4.2}
\]

The first bound follows from
\(\|q_y\|_\infty\le2q_y(0)\) and
\(\sigma_y^2\ge1/(12\|q_y\|_\infty^2)\).  The upper second bound is
Cauchy--Schwarz; the lower bound is the standard one-dimensional
log-concave first-moment comparison.  Thus (0.10) has the scale

\[
 \int{(h-\nu h)^2\over\sigma_y}\,d\nu
 \le C\int\sigma_y|\nabla h|^2\,d\nu.                       \tag{4.3}
\]

This is a weighted quotient Poincaré inequality.  It is not a bound on
\(\int1/\sigma_y\,d\nu\) unless the quotient admits a test function with a
uniformly nontrivial centered value and controlled gradient.

## 5. Matrix-valued quotient tests

Return to the regular normal foliation of Theorem 3.1.  Put

\[
 \bar u=\int u\,d\nu,\qquad \bar y=\int y\,d\nu.             \tag{5.1}
\]

### 5.1 Normal tests

For \(h_a(y)=a\cdot u(y)\),

\[
 \nabla_Sh_a=-{\mathsf A}_yP_{T_yS}a.                       \tag{5.2}
\]

Applying (3.6) for every \(a\) gives the Loewner inequality

\[
\begin{split}
 &\int q_y(0)(u-\bar u)\otimes(u-\bar u)\,d\nu\\
 &\preceq {1\over2}\int\!\!\int W_y(t)\,
 P_T{\mathsf A}_y(I-t{\mathsf A}_y)^{-2}
        {\mathsf A}_yP_T\,dt\,d\nu.
\end{split}                                                   \tag{5.3}
\]

Taking traces,

\[
 \int q_y(0)|u-\bar u|^2\,d\nu
 \le {1\over2}\int\!\!\int W_y(t)
       \left|{\mathsf A}_y(I-t{\mathsf A}_y)^{-1}\right|_F^2
       \,dt\,d\nu.                                          \tag{5.4}
\]

### 5.2 Position tests

For \(h_a(y)=a\cdot y\),

\[
 \nabla_Sh_a=P_{T_yS}a.
\]

Hence

\[
\begin{split}
 &\int q_y(0)(y-\bar y)\otimes(y-\bar y)\,d\nu\\
 &\preceq {1\over2}\int\!\!\int W_y(t)\,
 P_T(I-t{\mathsf A}_y)^{-2}P_T\,dt\,d\nu,
\end{split}                                                   \tag{5.5}
\]

and

\[
 \int q_y(0)|y-\bar y|^2\,d\nu
 \le {1\over2}\int\!\!\int W_y(t)
       \operatorname {tr}_{T_yS}(I-t{\mathsf A}_y)^{-2}
       \,dt\,d\nu.                                          \tag{5.6}
\]

### 5.3 What isotropy does and does not give

Let

\[
 m_y=\mathbb E_yT,\qquad \sigma_y^2=\operatorname {Var}_y(T).
\]

Conditional variance and isotropy imply, for every \(a\),

\[
 \int\sigma_y^2(a\cdot u(y))^2\,d\nu(y)
 \le\operatorname {Var}_\mu(a\cdot X)=|a|^2.                \tag{5.7}
\]

Equivalently, (0.12) holds.  Its trace is only

\[
                         \int\sigma_y^2\,d\nu\le n,          \tag{5.8}
\]

which by Jensen gives

\[
                         \int{d\nu\over\sigma_y}\ge n^{-1/2}.
                                                                  \tag{5.9}
\]

To improve (5.9), one would need to convert (5.3) or (5.5) into an
uncharged lower bound for \(\int q_y(0)d\nu\).  Their left sides instead
contain the factors \(|u-\bar u|^2\) and \(|y-\bar y|^2\), and their right
sides contain uncontrolled longitudinal curvature/metric energy.
Isotropy is unweighted and does not prevent \(q_y(0)\) from correlating
with quotient points where both centered factors are small.

Thus summing the matrix tests does not yield
\(\int1/\sigma_y\,d\nu\ge c\).  A proof would need an additional
coercive dichotomy:

* coherent normals must be converted through (5.7) into a fixed-direction
  variance bound; and
* dispersed normals must come with an independent upper bound on the
  right side of (5.4), including medial/contact terms.

The second bullet is not supplied by maximality or isotropy.

There is a useful algebraic sharpness model.  Let \(u\) be uniform on
\(S^{n-1}\), let \(T\sim N(0,n)\) be independent, and set \(X=Tu\).
Then

\[
 \mathbb E X\otimes X
 =\mathbb ET^2\,\mathbb Eu\otimes u=I,\qquad
 \mathbb E{1\over\sigma_y}={1\over\sqrt n}.                 \tag{5.10}
\]

This is not a log-concave measure on \(\mathbb R^n\), and the putative
zero set \(y=0\) is a singular common endpoint of every ray.  In
particular \(H(x)=a\cdot x/|x|\) is not globally Lipschitz at the origin,
so the normal matrix test is inadmissible unless one pays the medial
term.  The example is not a KLS counterexample; it proves that
(5.7), isotropy, and formal raywise matrix summation alone cannot improve
\(n^{-1/2}\).  Global log-concavity and the endpoint/medial geometry must
enter essentially.

## 6. Model audit

### 6.1 Gaussian halfspace

Let \(\mu=\gamma_n\) and \(f(x)=x_1\).  The rays are parallel to \(e_1\),
\(\nu=\gamma_{n-1}\), and

\[
 q(0)={1\over\sqrt{2\pi}},\qquad
 \mathbb E|T|=\sqrt{2\over\pi}=2q(0).                       \tag{6.1}
\]

Equation (0.10) becomes

\[
 \operatorname {Var}_{\gamma_{n-1}}(h)
 \le\int|\nabla h|^2\,d\gamma_{n-1}.                        \tag{6.2}
\]

For \(h=a\cdot y\), (6.2) is equality.  The normal test
\(h=a\cdot u\) is constant and disappears after centering.  Thus even the
cleanest exact model uses all the available quotient stability.

### 6.2 Half-cube and checkerboards

Let \(\mu\) be uniform on \([-\sqrt3,\sqrt3]^n\) and \(f(x)=x_1\).
Writing \(a=\sqrt3\),

\[
                  q(0)={1\over2a},\qquad
                  \mathbb E|T|={a\over2}.                   \tag{6.3}
\]

Equation (0.10) is

\[
 \operatorname {Var}_\nu(h)
 \le {a^2\over2}\int|\nabla h|^2\,d\nu
 ={3\over2}\int|\nabla h|^2\,d\nu.                          \tag{6.4}
\]

The sharp Neumann Poincaré constant of the quotient cube is
\(4a^2/\pi^2=12/\pi^2<3/2\), so (6.4) is valid with slack.

For a \(k\)-coordinate checkerboard mode, take

\[
 h_S(y)=2^{k/2}\prod_{j\in S}
       \sin\left({\pi y_j\over2a}\right),\qquad |S|=k.       \tag{6.5}
\]

Then \(\nu h_S=0\),

\[
 \int h_S^2\,d\nu=1,\qquad
 \int|\nabla h_S|^2\,d\nu={k\pi^2\over4a^2}
 ={k\pi^2\over12}.                                         \tag{6.6}
\]

The stability condition reads \(1\le k\pi^2/8\), already true for
\(k=1\) and increasingly wasteful for a high-frequency checkerboard.
A discontinuous checkerboard has infinite limiting Dirichlet cost.  The
quotient inequality therefore excludes oscillation, but it produces no
dimension-free gain beyond the known quotient Poincaré scale.

### 6.3 Radial distance

Let \(\mu\) be rotationally invariant and isotropic, let \(R=|X|\),
\(r_0=\mathbb ER\), and

\[
                              f(x)=|x|-r_0.                  \tag{6.7}
\]

The quotient is the unit sphere, \(u\) is uniform, and every ray has the
same density \(q(t)\), the law of \(R-r_0\).  The zero set is the sphere
of radius \(r_0\).  For \(h_a(u)=a\cdot u\),

\[
 |\nabla H((r_0+t)u)|^2
 ={\,|a|^2-(a\cdot u)^2\over(r_0+t)^2}.                     \tag{6.8}
\]

The exact necessary stability condition is

\[
 \boxed{\quad
 2q(0)\le(n-1)\int_{-r_0}^\infty
                 {W(t)\over(r_0+t)^2}\,dt.
 \quad}                                                       \tag{6.9}
\]

The position test \(h=a\cdot y=r_0a\cdot u\) is the same inequality
multiplied by \(r_0^2\).

If \(R-r_0\) is concentrated on an \(O(1)\) scale and \(r_0^2/n\to1\),
then

\[
 (n-1)\int {W(t)\over(r_0+t)^2}\,dt
 =\mathbb E|R-r_0|+o(1),                                   \tag{6.10}
\]

because \(\mathbb E(R-r_0)=0\).  Hence a radial-distance maximizer must
satisfy asymptotically

\[
                         2q(0)\le\mathbb E|R-r_0|.          \tag{6.11}
\]

For the standard Gaussian, \(R-r_0\) converges to
\(N(0,1/2)\).  Thus

\[
 q(0)\longrightarrow{1\over\sqrt\pi},\qquad
 \mathbb E|R-r_0|\longrightarrow{1\over\sqrt\pi},           \tag{6.12}
\]

and (6.11) fails by a factor \(2\).  The globally feasible perturbation is
geometrically a translation/ellipsoidal first-order deformation of the
sphere, so this is a genuine exclusion of Gaussian radial distance as a
maximizer.

For the isotropic radial-exponential law

\[
 d\mu_n(x)\propto e^{-\sqrt{n+1}|x|}\,dx,
\]

the radial variable is
\({\rm Gamma}(n,\alpha)\), \(\alpha=\sqrt{n+1}\), and
\(r_0=n/\alpha\).  In fact (6.9) is an equality for every \(n\ge2\).
To see this, write \(F_0=F_R(r_0)\), \(S_0=1-F_0\).  Since
\(\beta=S_0-F_0\),

\[
\begin{split}
 {1\over2}\int_{-r_0}^{\infty}{W(t)\over(r_0+t)^2}\,dt
 &=
 S_0\int_0^{r_0}{F_R(r)\over r^2}\,dr
 +F_0\int_{r_0}^{\infty}{1-F_R(r)\over r^2}\,dr\\
 &=S_0\int_0^{r_0}{q_R(r)\over r}\,dr
   -F_0\int_{r_0}^{\infty}{q_R(r)\over r}\,dr.
\end{split}                                                   \tag{6.13}
\]

For the gamma density,

\[
 {q_R(r)\over r}={q_R'(r)+\alpha q_R(r)\over n-1}.          \tag{6.14}
\]

Substitution in (6.13) cancels the two \(\alpha F_0S_0\) terms and gives

\[
 (n-1)\int_{-r_0}^{\infty}{W(t)\over(r_0+t)^2}\,dt
 =2q_R(r_0)=2q(0).                                          \tag{6.15}
\]

Moreover \(R-r_0\) converges to \(N(0,1)\), so

\[
 q(0)\longrightarrow{1\over\sqrt{2\pi}},\qquad
 \mathbb E|R-r_0|\longrightarrow\sqrt{2\over\pi}=2q(0).     \tag{6.16}
\]

Thus the radial-exponential translation mode saturates the quotient
second variation exactly, not just in the Gaussian limit.  More generally,
because the longitudinal coefficient is independent of \(u\), (3.6)
reduces for every smooth quotient function to

\[
 \operatorname {Var}_{S^{n-1}}(h)
 \le {1\over n-1}\int_{S^{n-1}}|\nabla_{S^{n-1}}h|^2\,d\nu.
                                                                  \tag{6.17}
\]

This is precisely the sharp spherical Poincaré inequality.  Its
degree-one modes \(h=a\cdot u\) are the equality cases.  Hence the full
quotient stability, not merely its matrix trace, is tautologically sharp
on this non-Gaussian radial model.

## 7. Circularity and feasibility audit

1. **Global admissibility.**  Equation (3.5) alone is not an admissible
   perturbation.  The sharp inequality (3.6) is asserted only for a
   signed-distance normal graph with the no-switching hypotheses of
   Theorem 3.1.  Equation (3.12) is the general rigorous fallback.
2. **Medial and support terms.**  Branch changes occur on
   \(O(|\varepsilon|)\)-thick regions and contribute at order
   \(\varepsilon^2\).  They may not be discarded.  In compact support
   they are the quotient counterpart of the support-contact tensor.
3. **Centering.**  The factor
   \(\operatorname {sign}(t)-\beta\), rather than merely
   \(\operatorname {sign}(t)\), is forced by subtracting
   \(\mu f_\varepsilon\).  Conditional mean zero is used only to simplify
   (4.1).
4. **Atoms at zero.**  If a conditional law has an atom at \(t=0\), the
   objective has a first-order \(|\varepsilon|\) term instead of
   \(q_y(0)h^2\).  Such rays require a separate nonsmooth analysis.
5. **No hidden Poincaré input.**  The Gaussian and cube Poincaré
   inequalities are model checks, not steps in the derivation.
6. **Inverse scale.**  The estimate
   \(q_y(0)\gtrsim1/\sigma_y\) uses balanced one-dimensional
   log-concavity.  Turning the matrix-weighted left sides of (5.3) and
   (5.5) into \(\int1/\sigma_y\) is precisely the missing correlation
   mechanism.
7. **Dimension tracking.**  Isotropy alone yields (5.9), with the forbidden
   \(n^{-1/2}\) loss.  Neither the normal nor position test removes it;
   the radial-exponential model makes the curvature balance asymptotically
   exact.
8. **Degenerate support.**  All gradients, rays, and dimensions are
   relative to the affine support.  A point mass is excluded.

The quotient second variation is therefore rigorous in its stated
feasibility class and has an exact normalized fallback.  It supplies a
useful stationarity and stability structure, but its matrix summation is
sharp/tautological rather than a dimension-free inverse-scale theorem.
