# Second variation of a signed-distance extremizer

This note studies the following restricted but exact setting.  Let
\(d\mu=\rho\,dx\), where \(\rho>0\) is smooth, and let \(f\) be a
mean-zero maximizer of

\[
  \mathcal D(\mu)=\sup_{\operatorname{Lip}(g)\leq 1}
       \int |g-\mu g|\,d\mu .
\]

Suppose first that the zero set \(\Sigma=\{f=0\}\) is an oriented
\(C^3\) hypersurface and that, away from a negligible cut locus, \(f\)
is its signed distance.  Write

\[
  F(y,t)=y+tN_y,\qquad N=\nabla f|_\Sigma,\qquad
  S_y=D_\Sigma N(y).
\]

On a regular ray, \(t\in I_y=(a_y,b_y)\), and

\[
 \begin{split}
  J_y(t)&=\det(I+tS_y),\\
  Z_y&=\int_{I_y}\rho(F(y,t))J_y(t)\,dt,\\
  q_y(t)&=Z_y^{-1}\rho(F(y,t))J_y(t),\\
  \eta(dy)&=Z_y\,d\mathcal H^{n-1}(y).
 \end{split}                                                    \tag{1}
\]

Thus \(\eta\) is a probability measure and

\[
  q_y(0)\eta(dy)=\rho(y)d\mathcal H^{n-1}(y).                    \tag{2}
\]

Put \(p=\mu(f>0)\), \(q=1-p\).  Optimal transport along the
calibrated rays gives the exact balance

\[
  \int_0^{b_y}q_y(t)dt=p,\qquad
  \int_{a_y}^0q_y(t)dt=q                                      \tag{3}
\]

for \(\eta\)-almost every ray.

The calculations below are classical differentiations under the integral
under, for example, the following sufficient hypotheses: the support of
\(h\) is contained in finitely many regular charts, \(h\in C_c^2(\Sigma)\),
the normal exponential maps in those charts have a uniform positive reach
on the part of the support of \(\mu\) under consideration, and all displayed
terms are integrable.  Truncation gives the same formulas whenever the two
sides are uniformly integrable.  Medial interfaces are treated separately
in Section 5.

## 1. Pointwise normal-graph variation

For small \(\varepsilon\), deform the zero surface by

\[
  \Sigma_\varepsilon
  =\{y+\varepsilon h(y)N_y:y\in\Sigma\},                         \tag{4}
\]

and let \(f_\varepsilon\) be its oriented signed distance.  At a fixed
point \(x=F(y,t)\), set

\[
  u(x)=\left.\partial_\varepsilon f_\varepsilon(x)\right|_0,
  \qquad
  w(x)=\left.\partial_{\varepsilon\varepsilon}
                    f_\varepsilon(x)\right|_0.
\]

Then

\[
 \boxed{
  u(F(y,t))=-h(y),\qquad
  w(F(y,t))=-\int_0^t
     \left|(I+sS_y)^{-1}\nabla_\Sigma h(y)\right|^2ds .}        \tag{5}
\]

Indeed, differentiating \(|\nabla f_\varepsilon|^2=1\) once gives

\[
  \partial_tu=\nabla f\mathbin\cdot\nabla u=0.
\]

Differentiating the boundary identity
\(f_\varepsilon(y+\varepsilon hN)=0\) once gives \(u(y)=-h(y)\).
This proves the first formula.  A second differentiation of the eikonal
equation gives

\[
  \partial_tw=-|\nabla u|^2.                                    \tag{6}
\]

Since

\[
  \nabla u(F(y,s))=-(I+sS_y)^{-1}\nabla_\Sigma h(y),             \tag{7}
\]

it remains only to check \(w(y)=0\).  Differentiating the boundary
identity twice gives

\[
  w(y)+2h\partial_Nu(y)+h^2\nabla^2f(y)[N,N]=0.
\]

Both of the last terms vanish: \(u\) is constant on normal rays and the
eikonal equation implies \(\nabla^2f\,N=0\).  Integration of (6) now proves
(5).  Notice that there is no missing factor two.

## 2. Exact centering term

Let

\[
  m_\varepsilon=\mu f_\varepsilon,
  \qquad g_\varepsilon=f_\varepsilon-m_\varepsilon,
  \qquad \alpha=\mu(\operatorname{sgn}f)=p-q.
\]

Since \(|\nabla f|=1\) on \(\Sigma\), the distributional identity
\((|r|)''=2\delta_0\), followed by coarea, gives

\[
 \begin{split}
  \left.\frac{d^2}{d\varepsilon^2}
       \int|g_\varepsilon|d\mu\right|_0
   &=\int \operatorname{sgn}(f)(w-\mu w)d\mu\\
   &\quad+2\int_\Sigma(u-\mu u)^2\rho\,d\mathcal H^{n-1}\\
   &=\int(\operatorname{sgn}f-\alpha)w\,d\mu\\
   &\quad+2\int_\Sigma(u-\mu u)^2\rho\,d\mathcal H^{n-1}.       \tag{8}
 \end{split}
\]

The mean of the first variation is

\[
  \mu u=-\int h\,d\eta.                                         \tag{9}
\]

Consequently, if \(\bar h=\int h\,d\eta\), the positive boundary term is

\[
  2\int q_y(0)(h(y)-\bar h)^2\eta(dy).                           \tag{10}
\]

The centering in (10) is against the ray-quotient probability \(\eta\),
not against normalized weighted area on \(\Sigma\).  Replacing it by the
latter is in general incorrect.

## 3. The smooth nonlocal stability inequality

Define the tail weight

\[
 \beta_y(s)=
 \begin{cases}
   q\displaystyle\int_s^{b_y}q_y(t)dt,&s>0,\\[6pt]
   p\displaystyle\int_{a_y}^{s}q_y(t)dt,&s<0.
 \end{cases}                                                     \tag{11}
\]

On the positive half-ray, \(\operatorname{sgn}f-\alpha=2q\), and on the
negative half-ray it equals \(-2p\).  Substitution of (5), followed by
Fubini, gives the exact identity

\[
  \int(\operatorname{sgn}f-\alpha)w\,d\mu
  =-2\int\!\int_{I_y}\beta_y(s)
     |(I+sS_y)^{-1}\nabla_\Sigma h(y)|^2ds\,\eta(dy).          \tag{12}
\]

Every \(f_\varepsilon\) is 1-Lipschitz.  Maximality of \(f\) therefore
implies that (8) is nonpositive, and hence

\[
 \boxed{
  \int q_y(0)(h-\bar h)^2d\eta
  \leq \int\!\int_{I_y}\beta_y(s)
       |(I+sS_y)^{-1}\nabla_\Sigma h|^2ds\,d\eta .}          \tag{13}
\]

This is a necessary stability inequality, not an independently available
Poincare inequality on the quotient.  In particular, invoking a quotient
Poincare estimate to bound its left side would be circular.

## 4. The normal test and the exact smooth budget

Take \(h_a(y)=a\mathbin\cdot N_y\) and sum (13) over an orthonormal basis
\((a_j)_{j=1}^n\).  With

\[
  m_N=\int N_y\eta(dy),                                          \tag{14}
\]

one has

\[
 \sum_j(h_{a_j}-\overline{h_{a_j}})^2=|N-m_N|^2                 \tag{15}
\]

and, because \(\nabla_\Sigma(a\cdot N)=S P_{T\Sigma}a\),

\[
 \sum_j|(I+sS)^{-1}\nabla_\Sigma h_{a_j}|^2
 =\operatorname{tr}\big(S^2(I+sS)^{-2}\big).                   \tag{16}
\]

Thus the smooth normal test reads

\[
 \int q_y(0)|N_y-m_N|^2d\eta
 \leq \int R_y\eta(dy),                                        \tag{17}
\]

where

\[
  R_y=\int_{I_y}\beta_y(s)
       \operatorname{tr}\big(S_y^2(I+sS_y)^{-2}\big)ds.         \tag{18}
\]

There is an exact dimension-free upper bound on (18).  Put

\[
  W_y(s)=-\log q_y(s).
\]

Log-concavity along the ray and the normal Jacobian formula give

\[
 W_y''(s)=N_y^T\nabla^2V(F(y,s))N_y
       +\operatorname{tr}\big(S_y^2(I+sS_y)^{-2}\big).          \tag{19}
\]

Let \(q_y(a_y+)\) and \(q_y(b_y-)\) denote the endpoint densities, with
value zero at an infinite endpoint.  Integration by parts on the two
half-rays gives

\[
 \begin{split}
  \int_0^{b_y}\!\left(\int_s^{b_y}q_y(t)dt\right)W_y''(s)ds
    &=-pW_y'(0)+q_y(0)-q_y(b_y-),\\
  \int_{a_y}^0\!\left(\int_{a_y}^{s}q_y(t)dt\right)W_y''(s)ds
    &= qW_y'(0)+q_y(0)-q_y(a_y+).                                \tag{20}
 \end{split}
\]

The \(W_y'(0)\) terms cancel after multiplying the first line by \(q\)
and the second by \(p\).  Therefore

\[
 \boxed{
  R_y\leq q_y(0)-q\,q_y(b_y-)-p\,q_y(a_y+)\leq q_y(0).}          \tag{21}
\]

The same statement follows by approximation when \(W_y''\) is a convex
second-derivative measure.  Formula (21) also identifies the exact omitted
nonnegative term:

\[
 \begin{split}
  q_y(0)-R_y\geq{}&q\,q_y(b_y-)+p\,q_y(a_y+)\\
  &+\int_{I_y}\beta_y(s)
          N_y^T\nabla^2V(F(y,s))N_y\,ds.                         \tag{22}
 \end{split}
\]

For \(p\in[\delta,1-\delta]\), the polynomial-Jacobian estimate

\[
  \sigma_y^2\|S_y\|_{HS}^2\leq C_\delta                         \tag{23}
\]

and the one-dimensional quantile estimate \(q_y(0)\asymp_\delta
\sigma_y^{-1}\) show that both sides of the smooth budget have the same
scale \(\sigma_y^{-1}\).  On a central interval this follows directly
from (23); the tails are controlled without a focal denominator loss by
(20).  There is no factor tending to zero as \(\sigma_y\to\infty\).

Consequently (17), even together with (23), does not force normal
alignment or concurrence.  If \(m_N=0\), its left side is exactly the
weighted perimeter \(P=\int q_y(0)d\eta\), while (21) gives at most \(P\).
Stability then forces equality in all smooth deficits in (22), but does not
by itself force \(S=0\) or a common normal.  Quantitatively, the smooth
second variation has exactly the wrong homogeneity to control the length
of a bad witness.

## 5. Exact medial-interface charge

Smooth ray formulas omit switching between competing distance charts.
For two smooth charts, the missing term can be computed without an
approximation.  Let

\[
  \phi_\varepsilon(x)
    =\min\{a(x)-\varepsilon h_i(x),
            b(x)-\varepsilon h_j(x)\},\qquad D=a-b.               \tag{24}
\]

Then

\[
 \phi_\varepsilon
 =\frac{a+b-\varepsilon(h_i+h_j)
       -|D-\varepsilon(h_i-h_j)|}{2}.                             \tag{25}
\]

Since \(\partial_{\varepsilon\varepsilon}|D-\varepsilon H|_{\varepsilon=0}
=2H^2\delta_0(D)\),

\[
 \left.\partial_{\varepsilon\varepsilon}\phi_\varepsilon
 \right|_0
 =-(h_i-h_j)^2\delta_0(D).                                       \tag{26}
\]

If \(0\) is a regular value of \(D\), coarea gives the exact coefficient

\[
 \boxed{
  \left.\frac{d^2}{d\varepsilon^2}
       \int\phi_\varepsilon\rho\,dx\right|_0
  =-\int_{\{a=b\}}
       \frac{(h_i-h_j)^2}{|\nabla a-\nabla b|}\rho
       \,d\mathcal H^{n-1}.}                                    \tag{27}
\]

For distance charts meeting on a medial interface \(M_{ij}\),
\(\nabla a=N_i\) and \(\nabla b=N_j\) on the positive side.  On the
negative side both gradients change sign, leaving the denominator
\(|N_i-N_j|\) unchanged.  Combining (27) with the factors
\(\operatorname{sgn}f-\alpha=2q\) and \(-2p\), respectively, yields the
graph energies

\[
 \begin{split}
  \mathcal C_+(h)
    &=\sum_{i<j}\int_{M_{ij}^+}
       \frac{(h_i-h_j)^2}{|N_i-N_j|}\rho\,d\mathcal H^{n-1},\\
  \mathcal C_-(h)
    &=\sum_{i<j}\int_{M_{ij}^-}
       \frac{(h_i-h_j)^2}{|N_i-N_j|}\rho\,d\mathcal H^{n-1}.     \tag{28}
 \end{split}
\]

Each generic two-chart interface is counted once.  The complete
piecewise-smooth stability inequality is

\[
 \boxed{
  \int q_y(0)(h-\bar h)^2d\eta
  \leq \int R_y(h)d\eta+q\mathcal C_+(h)+p\mathcal C_-(h).}       \tag{29}
\]

Thus the omitted object is a weighted graph Dirichlet form on the medial
adjacency graph, not a harmless null-set error.

For \(h_a=a\cdot N\), summation over an orthonormal basis gives, on every
interface,

\[
 \sum_j\frac{(a_j\cdot(N_i-N_j))^2}{|N_i-N_j|}
 =|N_i-N_j|.                                                      \tag{30}
\]

Hence the summed singular term is precisely the weighted total normal
jump

\[
 q\sum_{i<j}\int_{M_{ij}^+}\rho|N_i-N_j|d\mathcal H^{n-1}
 +p\sum_{i<j}\int_{M_{ij}^-}\rho|N_i-N_j|d\mathcal H^{n-1}.      \tag{31}
\]

At a non-generic multiway junction, strata of codimension at least two do
not contribute directly to (27).  Focal strata at which two footpoints
coalesce are limits of (28); the denominator can vanish and must not be
discarded before the numerator.  A fully general statement for an arbitrary
cut locus requires a second epi-derivative theorem for the distance
envelope.  Formula (28) proves the needed extension for finite polyhedral
charts and generic stratified medial axes, but such a theorem has not been
established here for an arbitrary singular focal set.  Treating the latter
as measure zero is invalid.

## 6. Flat two-chart wedge

Let a planar wedge have opening angle \(\theta\in(0,\pi)\).  Its unsigned
distance is the minimum of the two affine distance charts with inward unit
normals \(N_1,N_2\).  Their angle is \(\pi-\theta\), so

\[
  |N_1-N_2|=2\cos(\theta/2).                                     \tag{32}
\]

The smooth curvature term vanishes.  If the two faces are displaced by
constants \(h_1,h_2\), (27) gives

\[
  F''(0)=-\frac{(h_1-h_2)^2}{2\cos(\theta/2)}
          \int_M\rho\,d\mathcal H^1,                             \tag{33}
\]

where \(M\) is the angular bisector.  For the vector family
\(h_a=a\cdot N\), the sum of the positive charges is exactly

\[
  2\cos(\theta/2)\int_M\rho\,d\mathcal H^1.                     \tag{34}
\]

This example verifies both the denominator and the sign in (28).

## 7. Gaussian angular fan

Let \(\gamma_2\) be standard Gaussian measure and fix \(m\geq2\).  Put
\(\delta=\pi/m\), let the zero set be the union of the \(2m\) rays with
angles \(k\delta\), and give the intervening wedges alternating signs.
The resulting signed distance \(f_m\) has \(p=q=1/2\).  Every zero branch
and every distance chart is flat.  Consecutive oriented normals satisfy

\[
  |N_k-N_{k+1}|=2\cos(\delta/2).                                 \tag{35}
\]

There is one medial ray in each wedge.  Set

\[
 c_G=\int_0^\infty(2\pi)^{-1}e^{-r^2/2}dr
     =\frac1{2\sqrt{2\pi}}.                                     \tag{36}
\]

The Gaussian perimeter of the zero fan is

\[
  P=2mc_G=\frac m{\sqrt{2\pi}}.                                 \tag{37}
\]

Rotational symmetry gives \(m_N=0\).  More explicitly,
\(f_m(R_\delta x)=-f_m(x)\), hence
\(v=\int\nabla f_m\,d\gamma_2\) satisfies \(v=-R_\delta v\), and
\(v=0\) for \(m\geq2\).  Therefore the summed left side of (29) is
exactly \(P\).

There is no smooth term.  Each of the \(2m\) medial rays has Gaussian
weighted length \(c_G\); its sign weight is \(1/2\).  Formula (31) gives

\[
  \frac12(2m)\,2\cos(\delta/2)c_G
   =P\cos\!\left(\frac\pi{2m}\right).                            \tag{38}
\]

Thus (29) would require

\[
  P\leq P\cos\!\left(\frac\pi{2m}\right),                      \tag{39}
\]

which is false.  The angular fan is therefore excluded as an exact
first-moment extremizer by the singular second variation, even though it
has balanced calibrated rays.  This is a genuine gain over the smooth
normal-flow calculation.

The quantitative deficit is weak:

\[
  P\left(1-\cos\frac\pi{2m}\right)\asymp \frac1m.                \tag{40}
\]

It tends to zero as the number of sectors grows.  Consequently this test
does not by itself give a scale-free exclusion of fan-like approximate
extremizers.

## 8. A global translation identity: all focal charge at once

For the special tests \(h_a=a\cdot N\), the arbitrary focal-set issue in
Section 5 can be bypassed completely.  Define

\[
  f_{\varepsilon,a}(x)=f(x-\varepsilon a).                       \tag{41}
\]

This is 1-Lipschitz for every \(\varepsilon\), and its zero set is the
translate \(\Sigma+\varepsilon a\).  At every regular boundary point its
normal velocity is \(a\cdot N\).  More importantly, (41) moves the entire
cut locus correctly, without choosing distance charts.

Assume in this subsection that \(V=-\log\rho\) is \(C^1\), with enough
integrability to justify the displayed integrations by parts.  The result
extends by smoothing and truncation whenever the boundary terms converge.
Let

\[
  m_N=\int\nabla f\,d\mu=\int N_y\eta(dy),\qquad
  B_N=\int_\Sigma N\rho\,d\mathcal H^{n-1},                      \tag{42}
\]

and set

\[
  K=\int(\operatorname{sgn}f-\alpha)
        \langle\nabla V,\nabla f\rangle d\mu.                  \tag{43}
\]

The trace, over an orthonormal basis of translation directions, of the
centered second derivative is

\[
 \boxed{
  \sum_{j=1}^n\left.\frac{d^2}{d\varepsilon^2}
       \int|f(x-\varepsilon a_j)-\mu f(\cdot-\varepsilon a_j)|
          d\mu(x)\right|_0
  =2L-2P+K,}                                                     \tag{44}
\]

where

\[
  L=\int_\Sigma|N-m_N|^2\rho\,d\mathcal H^{n-1},qquad
  P=\int_\Sigma\rho\,d\mathcal H^{n-1}.                         \tag{45}
\]

Here is a proof which never differentiates \(f\) twice and therefore
remains valid across a stratified cut locus.  After the change of variables
\(z=x-\varepsilon a\), put

\[
 \rho_{\varepsilon,a}(z)=\rho(z+\varepsilon a),\qquad
 m_{\varepsilon,a}=\int f\rho_{\varepsilon,a}dz.
\]

At zero,

\[
  m'_{0,a}=\int f\,\partial_a\rho=-\int\partial_af\,d\mu
            =-a\cdot m_N.                                      \tag{46}
\]

Direct differentiation of
\(\int|f-m_{\varepsilon,a}|\rho_{\varepsilon,a}\) gives

\[
 \begin{split}
 J_a''(0)={}&\int(|f|-\alpha f)\partial_{aa}\rho\,dx
 -2m'_{0,a}\int\operatorname{sgn}(f)\partial_a\rho\,dx\\
 &+2P(m'_{0,a})^2.                                               \tag{46a}
 \end{split}
\]

The last term follows from \((|r|)''=2\delta_0\) and coarea.  On the
reduced zero boundary,

\[
  D(\operatorname{sgn}f)
   =2N\,\mathcal H^{n-1}\!\restriction\Sigma,
\]

and hence

\[
  \int\operatorname{sgn}(f)\partial_a\rho\,dx=-2a\cdot B_N.     \tag{46b}
\]

Finally, the weak gradient of the Lipschitz function
\(|f|-\alpha f\) is \((\operatorname{sgn}f-\alpha)\nabla f\).
Summing the first term of (46a) over a basis and integrating by parts gives

\[
 \int(|f|-\alpha f)\Delta\rho\,dx
 =-\int(\operatorname{sgn}f-\alpha)
       \langle\nabla f,\nabla\rho\rangle dx=K.                \tag{46c}
\]

Thus the sum of (46a) is
\(K-4m_N\cdot B_N+2P|m_N|^2=2L-2P+K\), proving (44).

The new term has a sign.  On every balanced ray,

\[
 \begin{split}
 K_y
 &=\int(\operatorname{sgn}t-\alpha)V_y'(t)q_y(t)dt\\
 &=2pq\left(\mathbb E[V_y'(T)\mid T>0]
             -\mathbb E[V_y'(T)\mid T<0]\right)\geq0,            \tag{47}
 \end{split}
\]

because \(V_y(t)=V(y+tN_y)\) is convex and hence \(V_y'\) is
nondecreasing.  Thus \(K=\int K_y d\eta\geq0\).

Equivalently, integration by parts against the tail kernel (11) gives the
curvature representation

\[
  \frac12K_y=\int_{I_y}\beta_y(s)V_y''(s)ds,                    \tag{47a}
\]

with convex second derivatives interpreted as measures.  This is exactly
the density-curvature deficit that appeared in (22).

Maximality makes every directional second derivative in (44) nonpositive.
Consequently

\[
 \boxed{
   L+\frac12K\leq P.}                                           \tag{48}
\]

Formula (48) is the exact full translation-stability inequality, including
all smooth, medial, multiway, and focal contributions.  In the notation of
Sections 4--5, the total translation charge is therefore

\[
  \boxed{\mathcal C_{\rm full}=P-\frac12K\leq P.}                \tag{49}
\]

For a finite generic chart system, (49) equals the sum of (18) and (31).
For an arbitrary cut locus, (49) gives the aggregate charge even though a
stratum-by-stratum epi-derivative has not been constructed.

Expanding \(L\) gives a useful rigidity form:

\[
  L=P(1+|m_N|^2)-2m_N\cdot B_N.                                 \tag{50}
\]

Since \(|B_N|\leq P\), (48) implies

\[
  \frac K{2P}\leq 2|m_N|-|m_N|^2.                               \tag{51}
\]

In particular, a zero quotient-mean normal forces \(K=0\) and equality in
the entire charge budget.  Conditional equality in (47) says that
\(V_y'\) is constant across the two sign halves of almost every active ray.
This is a genuine rigidity condition, but it still allows the geometry or
the support boundary to provide the one-dimensional log-concavity.

For the Gaussian fan, \(V(x)=|x|^2/2\), \(\alpha=0\), and Euler's identity
for the one-homogeneous \(f_m\) gives \(x\cdot\nabla f_m=f_m\).  Hence

\[
  K=\int|f_m|d\gamma_2
   =2P\left(1-\cos\frac\pi{2m}\right).                           \tag{52}
\]

Equations (49) and (52) reproduce exactly the polyhedral calculation
\(\mathcal C_{\rm full}=P\cos(\pi/(2m))\).  Thus (48) also verifies that no
unaccounted point charge at the common fan vertex is missing.

There is a further scale-free estimate when the separator perimeter itself
is small.  For \(p\in[\delta,1-\delta]\), one-dimensional log-concavity
gives \(q_y(0)\geq c_\delta/\sigma_y\), and isotropy gives

\[
  \int\sigma_y^2N_yN_y^T\eta(dy)\preceq I.                      \tag{53}
\]

For any \(r>0\), split the quotient into \(\{\sigma_y\geq r\}\) and its
complement.  In the direction \(m_N/|m_N|\), Cauchy--Schwarz and (53) bound
the contribution of the first set by \(r^{-1}\), while the mass of the
second is at most

\[
  r\int\sigma_y^{-1}d\eta\leq C_\delta rP.
\]

Optimizing at \(r\asymp P^{-1/2}\) yields

\[
  \boxed{|m_N|\leq C_\delta\sqrt P.}                             \tag{54}
\]

Combining (51) and (54) gives \(K/P\leq C_\delta\sqrt P\).  This says that
a small-perimeter extremal separator must be almost translation-invariant
in the convexity sense (47).  It does not, by itself, show that the zero cut
of a first-moment extremizer has small perimeter, nor does it classify the
near-equality cases of (47); those are the two missing steps in turning
(48) into KLS.

## 9. Product cylinders and the remaining obstruction

For a flat product separator \(\Sigma=\{0\}\times\mathbb R^{n-1}\) and
\(f(t,z)=t\), one has \(S=0\), no medial interfaces, and \(N\) is constant.
For \(h_a=a\cdot N\), centering makes the left side of (29) zero, so the
normal test is an identity and gives no information.

More generally, if \(f(x,z)=f_0(x)\) and
\(\mu=\mu_0\otimes\nu\), then the regular and medial terms, the boundary
term, and the quotient centering in (29) all integrate out the \(z\)
factor and reduce exactly to the corresponding quantities for
\((\mu_0,f_0)\).  The second variation therefore cannot rule out product
cylinders or produce an induction gain; it only sees the irreducible core.

An artificial disjoint collection of flat cylinders with dispersed normals
would have positive left side and zero graph/curvature energy, and hence
could not be an extremizer if its components could be varied independently.
Such a disconnected construction is not the support of a nondegenerate
log-concave measure.  Connected fan-like configurations replace the missing
energy by the medial graph charge (31).

## 10. Verdict for this route

The second variation supplies an exact new necessary condition, (29), and
it detects the Gaussian fan through a singular normal-jump charge that is
invisible in smooth ray coordinates.  It does not yet prove a dimension-free
bound:

1. the smooth charge has the same \(1/\sigma_y\) homogeneity as weighted
   perimeter, even after the polynomial-Jacobian estimate;
2. the medial charge can be arbitrarily close to saturating the boundary
   term, as the many-sector Gaussian fan shows;
3. the whole inequality tensorizes exactly on product cylinders; and
4. a general focal-set extension for arbitrary deformations needs a rigorous
   second epi-derivative measure for the distance envelope, although the
   translation identity (48) already captures its full aggregate for
   \(h=a\cdot N\).

A successful continuation would need a strict, dimension-free deficit in
the combined smooth-plus-medial charge for a *global maximizer*, using
log-concavity and the large ray scale.  Neither (21), (23), nor the graph
energy (31) alone provides such a deficit.
