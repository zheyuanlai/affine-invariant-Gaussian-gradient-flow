# Variance-normalized random lines for a median halfspace

## 0. Outcome

Let \(\mu\) be an isotropic log-concave probability on \(\mathbb R^n\),
let \(u\in S^{n-1}\), let \(h\) be a median of \(u\cdot X\), and put

\[
 E=\{x:u\cdot x\le h\}.
\]

For a direction \(\theta\), disintegrate over the lines parallel to
\(\theta\), and let \(p_{\theta,y}\) and \(\sigma_{\theta,y}\) be the
conditional mass of \(E\) and the conditional standard deviation of the
line coordinate.  The desired special-case estimate is

\[
 \mathbb E_\theta\int
 {\min(p_{\theta,y},1-p_{\theta,y})\over\sigma_{\theta,y}}
 \,d(\pi_{\theta^\perp}\mu)(y)\ge {c\over\sqrt n}.       \tag{0.1}
\]

This report obtains three exact pieces of progress, but not a proof of
(0.1).

1.  For every fixed \(\theta\), (0.1) is reduced without loss to a
    conditional-threshold inequality for a log-concave pair \((Y,T)\):

    \[
      \mathbb E\left[{\min(\mathbb P(T\le0\mid Y),
                    \mathbb P(T>0\mid Y))
                    \over\sqrt{\operatorname {Var}(T\mid Y)}}\right]
       \ge {c\over\sqrt{\operatorname {Var}T}},            \tag{CT}
    \]

    where zero is a median of \(T\).  If (CT) holds, then the left side
    for a fixed direction is at least \(c|u\cdot\theta|\), and averaging
    gives (0.1).

2.  A one-constraint needle decomposition in the base variable reduces
    (CT) in every dimension to (CT) for a log-concave probability on
    \(\mathbb R^2\).  The reduction loses only a factor \(\sqrt2\), not a
    dimension-dependent factor.  Thus the exact remaining assertion is a
    bivariate conditional-threshold theorem.

3.  For a uniform convex body, the two reach functions on the cut section
    are concave.  If their section means are comparable, Berwald's
    inequality and a Brunn--Minkowski midpoint argument give a completely
    explicit chord-balance bound.  Global half-mass, however, does **not**
    imply comparable reach means: fibers lying wholly in one phase are
    absent from the cut section.  This is the precise pure-fiber gap.

The bivariate theorem passes a sharply asymmetric tilted exponential wedge:
its uncertainty is asymptotic to \((\log2)/(2\sqrt{\operatorname {Var}T})\).
A standard epigraph lift to uniform chords does not prove it.  An explicit
plateau--exponential line law loses a factor \(\log(1/q)\), where \(q\) is
the conditional minority mass.  Consequently a level-set or lifting proof
must control how these logarithmic losses aggregate; discarding them is
invalid.

Gaussian halfspaces, cube cuts, regular-simplex caps, right circular cones,
and product exponential measures all have the asserted \(n^{-1/2}\) scale.
No counterexample is obtained.  The unresolved statement is (CT) already
in dimension two.

## 1. Exact fixed-direction formula

Assume first that \(\mu\) is full-dimensional with density
\(\rho=e^{-V}\), normalized so that \(\int\rho=1\).  Put

\[
 H=\{x:u\cdot x=h\},\qquad P_H=\int_H\rho\,d\mathcal H^{n-1}.
\]

The one-dimensional marginal \(u\cdot X\) is log-concave, has variance
one, and has median \(h\).  The standard one-dimensional median-density
bound gives

\[
                         P_H\ge {1\over4\sqrt3}.            \tag{1.1}
\]

Fix \(\theta\in S^{n-1}\), and reverse its sign if necessary so that

\[
                         a=u\cdot\theta>0.                 \tag{1.2}
\]

For \(x\in H\), define the unnormalized line density and its two masses

\[
 g_x(t)=\rho(x+t\theta),\qquad
 M_\pm(x)=\int_0^\infty g_x(\pm t)\,dt,qquad
 Z(x)=M_+(x)+M_-(x).                                      \tag{1.3}
\]

Let \(\sigma(x)\) be the standard deviation of the probability density
\(g_x/Z(x)\).  The orthogonal projection

\[
                         x\longmapsto \pi_{\theta^\perp}x
\]

maps \(H\) bijectively onto \(\theta^\perp\), and its
\((n-1)\)-dimensional Jacobian is \(a\).  Therefore

\[
\boxed{
 \mathcal U_\theta(E)
 =a\int_H{\min(M_+(x),M_-(x))\over\sigma(x)}
           \,d\mathcal H^{n-1}(x).}                       \tag{1.4}
\]

Equivalently, with the boundary probability

\[
 d\nu_H(x)={\rho(x)\over P_H}\,d\mathcal H^{n-1}(x)
\]

and the one-dimensional normalized tail lengths

\[
 \ell_\pm(x)={M_\pm(x)\over\rho(x)},\qquad
 r(x,\theta)={\min(\ell_+(x),\ell_-(x))\over\sigma(x)},  \tag{1.5}
\]

formula (1.4) is

\[
             \boxed{\mathcal U_\theta(E)
                     =aP_H\,\mathbb E_{\nu_H}r(X,\theta).} \tag{1.6}
\]

Thus a universal lower bound on the boundary average of \(r\) would prove
the fixed-direction estimate \(\mathcal U_\theta(E)\ge ca\).

There is also a coordinate-free probabilistic formulation.  Define the
affine line coordinate and base point

\[
 T={u\cdot X-h\over a},\qquad Y=X-T\theta\in H.            \tag{1.7}
\]

The map \(X\mapsto(Y,T)\) is affine and invertible, so its law is
log-concave.  Conditional on \(Y=x\), the law of \(T\) is exactly
\(g_x/Z(x)\), and \(T=0\) is the cutting threshold.  Moreover,

\[
 \operatorname {Var}T={1\over a^2},
 \qquad \mathbb P(T\le0)=\mathbb P(T\ge0)={1\over2}.      \tag{1.8}
\]

Consequently (CT) implies

\[
                         \mathcal U_\theta(E)\ge ca.       \tag{1.9}
\]

Finally,

\[
 \mathbb E_{\theta\sim S^{n-1}}|u\cdot\theta|
 ={\Gamma(n/2)\over\sqrt\pi\,\Gamma((n+1)/2)}
 \ge {1\over\sqrt{2n}}.                                  \tag{1.10}
\]

Hence (CT), with one universal constant in every base dimension, implies
(0.1) with constant \(c/\sqrt2\).

## 2. Reduction of (CT) to a bivariate theorem

This section isolates the remaining dimensional issue.  Let \((Y,T)\)
have a log-concave density \(G(y,t)\) on \(\mathbb R^d\times\mathbb R\),
and suppose that zero is a median of \(T\).  Put

\[
\begin{aligned}
 Z(y)&=\int_{\mathbb R}G(y,t)\,dt,\\
 M_-(y)&=\int_{-\infty}^0G(y,t)\,dt,\qquad
 M_+(y)=\int_0^\infty G(y,t)\,dt,\\
 b(y)&={M_-(y)-M_+(y)\over Z(y)}.
\end{aligned}                                             \tag{2.1}
\]

The base law \(Z(y)dy\) is log-concave by Pr\'ekopa, and

\[
                         \int b(y)Z(y)dy=0.                \tag{2.2}
\]

Use the one-function convex localization theorem in the following exact
form.  A log-concave probability \(\nu\) on \(\mathbb R^d\) and an
integrable function \(b\) with \(\int b\,d\nu=0\) admit a measurable
needle disintegration

\[
                         \nu=\int\nu_\alpha\,d\lambda(\alpha), \tag{2.3}
\]

such that, for \(\lambda\)-almost every \(\alpha\), \(\nu_\alpha\) is
supported on a segment or line, has a log-concave one-dimensional density
there, and

\[
                         \int b\,d\nu_\alpha=0.            \tag{2.4}
\]

For unbounded support and merely integrable \(b\), this follows by first
restricting to a large convex compact set, replacing \(b\) by a bounded
truncation plus a constant correction, applying the compact needle theorem,
and taking a tight subsequential limit.  The conditional quantities below
are nonnegative, so Fatou handles the uncertainty functional.  The first
and second moments of \(T\) are uniformly integrable under the same
truncation.

On a needle parametrized by \(s\mapsto y_\alpha(s)\), the needle density
has the usual form proportional to
\(e^{\ell_\alpha(s)}Z(y_\alpha(s))\), with \(\ell_\alpha\) affine.  Restoring
the conditional law of \(T\) gives the joint needle density

\[
 G_\alpha(s,t)=c_\alpha e^{\ell_\alpha(s)}
                         G(y_\alpha(s),t).                  \tag{2.5}
\]

It is log-concave on \(\mathbb R^2\).  Equation (2.4) says exactly that
zero is a median of its \(T\)-marginal.  Multiplication by the factor in
(2.5) does not change the conditional law of \(T\) at a fixed \(s\).
Therefore the original uncertainty disintegrates exactly:

\[
                         \mathcal U=\int\mathcal U_\alpha
                                      \,d\lambda(\alpha). \tag{2.6}
\]

Assume now that (CT) has been proved for every log-concave probability on
\(\mathbb R^2\), with constant \(c_2\).  If

\[
                         v_\alpha=\operatorname {Var}_\alpha T,
\]

then

\[
                         \mathcal U_\alpha\ge {c_2\over\sqrt{v_\alpha}}.
                                                                    \tag{2.7}
\]

Since each needle has median zero,

\[
 v_\alpha\le\mathbb E_\alpha T^2,
 \qquad
 \int v_\alpha\,d\lambda(\alpha)\le\mathbb ET^2.        \tag{2.8}
\]

For any real random variable with median zero, Cantelli's inequality gives

\[
                         |\mathbb ET|\le\sqrt{\operatorname {Var}T}.
                                                                    \tag{2.9}
\]

Indeed, if \(m=\mathbb ET>0\), then
\(1/2\le\mathbb P(T-m\le-m)\le
\operatorname {Var}T/(\operatorname {Var}T+m^2)\), and reflection handles
\(m<0\).  Hence

\[
                         \mathbb ET^2\le2\operatorname {Var}T.      \tag{2.10}
\]

Jensen's inequality for \(x\mapsto x^{-1/2}\), followed by
(2.6)--(2.10), gives

\[
 \mathcal U
 \ge c_2\int v_\alpha^{-1/2}d\lambda(\alpha)
 \ge {c_2\over\sqrt{\int v_\alpha d\lambda(\alpha)}}
 \ge {c_2\over\sqrt{2\operatorname {Var}T}}.             \tag{2.11}
\]

Degenerate needles with \(v_\alpha=0\) either have both phases of zero
mass, contrary to (2.4), or are interpreted by approximation; they cannot
decrease the left side.  We have proved:

> **Exact reduction.**  The bivariate version of (CT) implies (CT) in
> every dimension, with only the factor \(\sqrt2\).

No proof of the bivariate assertion is supplied below.  This is the single
remaining lemma in this route.

## 3. What convexity gives for uniform bodies

Let \(K\subset\mathbb R^n\) be a convex body, let \(\mu\) be uniform on
\(K\), and retain \(H,u,h,\theta,a\) from Section 1.  Put

\[
                         \Sigma=K\cap H.
\]

For \(x\in\Sigma\), define the two reaches

\[
 \tau_\pm(x)=\sup\{t\ge0:x\pm t\theta\in K\}.             \tag{3.1}
\]

The conditional law on a mixed chord is uniform, so

\[
 \sigma(x)={\tau_+(x)+\tau_-(x)\over\sqrt{12}},\qquad
 q(x)={\min(\tau_+(x),\tau_-(x))\over\tau_+(x)+\tau_-(x)}. \tag{3.2}
\]

The projection Jacobian from \(H\) to \(\theta^\perp\) is \(a\).  Thus

\[
\boxed{
 \mathcal U_\theta(E)=
 \sqrt{12}\,a\,{\mathcal H^{n-1}(\Sigma)\over |K|}
 \mathbb E_{x\sim\mathrm{Unif}(\Sigma)}
 {\min(\tau_+(x),\tau_-(x))\over\tau_+(x)+\tau_-(x)}.}    \tag{3.3}
\]

### 3.1 Concavity of the reaches

Both \(\tau_+\) and \(\tau_-\) are nonnegative concave functions on
\(\Sigma\).  For example, if \(t_i<\tau_+(x_i)\), convexity of \(K\) gives

\[
 \lambda(x_1+t_1\theta)+(1-\lambda)(x_2+t_2\theta)
 \in K,
\]

and letting \(t_i\uparrow\tau_+(x_i)\) proves

\[
 \tau_+(\lambda x_1+(1-\lambda)x_2)
 \ge\lambda\tau_+(x_1)+(1-\lambda)\tau_+(x_2).            \tag{3.4}
\]

The same argument applies to \(\tau_-\).

### 3.2 A quantitative overlap lemma

The following lemma is useful whenever a separate argument balances the
two reach means.

**Lemma 3.1.**  Let \(D\subset\mathbb R^d\) be a convex body and let
\(f,g:D\to[0,\infty)\) be concave.  Suppose

\[
                         \mathbb E_D f=\mathbb E_D g=m>0.  \tag{3.5}
\]

Then

\[
 \mathbb E_D{\min(f,g)\over f+g}\ge {1\over2048}.         \tag{3.6}
\]

Here the ratio is set to zero where \(f+g=0\).

**Proof.**  Berwald's inequality, with \(p=1,q=2\), gives

\[
 \mathbb E_Df^2\le {2(d+1)\over d+2}(\mathbb E_Df)^2
 \le2m^2,                                                 \tag{3.7}
\]

and the same for \(g\).  Paley--Zygmund therefore gives

\[
 \mathbb P_D(f\ge m/2)\ge {1\over8},\qquad
 \mathbb P_D(g\ge m/2)\ge {1\over8}.                    \tag{3.8}
\]

Let \(A=\{f\ge m/2\}\) and \(B=\{g\ge m/2\}\).  By
Brunn--Minkowski,

\[
                         |(A+B)/2|\ge |D|/8.               \tag{3.9}
\]

Concavity and nonnegativity show that on \((A+B)/2\), both \(f\) and
\(g\) are at least \(m/4\).  Markov's inequality gives

\[
                         \mathbb P_D(f+g>32m)\le {1\over16}.\tag{3.10}
\]

Thus a set of normalized volume at least \(1/16\) satisfies

\[
                         f,g\ge m/4,qquad f+g\le32m.      \tag{3.11}
\]

On this set the ratio in (3.6) is at least \(1/128\), proving (3.6).
\(\square\)

The constants can be improved, but their size is irrelevant.

### 3.3 The pure-fiber obstruction

It is tempting, but wrong, to infer (3.5) from \(|K\cap E|=|K|/2\).
The coordinate map

\[
                         (x,t)\longmapsto x+t\theta       \tag{3.12}
\]

has Jacobian \(a\), but

\[
 a\int_\Sigma\tau_+(x)dx
\]

is only the volume of the positive portions of those \(\theta\)-fibers
which actually cross \(H\) inside \(K\).  Fibers whose complete chord lies
in \(E^c\) do not meet \(\Sigma\), and their positive volume is missing.
The analogous statement holds on the negative side.  Therefore

\[
 a\int_\Sigma\tau_+\le |K\cap E^c|={|K|\over2},\qquad
 a\int_\Sigma\tau_-\le |K\cap E|={|K|\over2},             \tag{3.13}
\]

but equality need not hold in either inequality and the two deficits need
not agree.

This is not a null-set issue.  In a sheared cylinder or a cone, a positive
volume family of lines can lie completely on one side of the cutting
hyperplane.  Equations (3.3)--(3.4) and Lemma 3.1 are exact; the missing
statement is a global comparison of the two reach means, or a way of
charging their imbalance to the pure fibers.  Simply writing half-volume
as \(a\int_\Sigma\tau_\pm\) is an invalid proof.

## 4. The epigraph lift and its unavoidable logarithm

The standard epigraph representation makes the geometry in Section 3
available for a general log-concave density, but it does not preserve the
functional with a universal constant.

Write \(\rho=e^{-V}\), and consider

\[
 \mathcal C=\{(x,r):r\ge V(x)\},\qquad
 d\widetilde\mu(x,r)=e^{-r}\mathbf1_{\mathcal C}(x,r)\,dx\,dr. \tag{4.1}
\]

The set \(\mathcal C\) is convex, \(\widetilde\mu\) is log-concave, and
its \(x\)-marginal is \(\mu\).  Conditional on \(r\) and on the projection
orthogonal to \((\theta,0)\), the line law is uniform on a sublevel chord
of \(V\).  Thus the refined lift sees uniform chords.  The original
functional conditions only on the \(x\)-projection and mixes all these
levels.

The loss under this coarsening is unbounded.  Let \(0<\varepsilon<1/4\)
and consider the one-dimensional log-concave probability density

\[
 g_\varepsilon(t)={1\over1+\varepsilon}
 \begin{cases}
  1,&-\varepsilon\le t\le0,\\
  e^{-t},&t\ge0,\\
  0,&t<-\varepsilon.
 \end{cases}                                               \tag{4.2}
\]

At the threshold zero its minority mass is

\[
                         q_\varepsilon={\varepsilon\over1+\varepsilon},\tag{4.3}
\]

and direct integration gives

\[
 \mathbb ET={1-\varepsilon^2/2\over1+\varepsilon},\qquad
 \mathbb ET^2={2+\varepsilon^3/3\over1+\varepsilon},
 \qquad \sigma_\varepsilon\longrightarrow1.              \tag{4.4}
\]

Hence the original normalized uncertainty is

\[
                         {q_\varepsilon\over\sigma_\varepsilon}
                         \asymp\varepsilon.                \tag{4.5}
\]

For an unnormalized version of (4.2), the superlevel chord at height
\(s=e^{-r}\), \(r\ge0\), has reaches \(\varepsilon\) and \(r\).  The
refined uniform-chord contribution, divided by \(\sqrt{12}\), is

\[
 L_\varepsilon
 =\int_0^\infty e^{-r}{\min(\varepsilon,r)\over\varepsilon+r}\,dr.
                                                                    \tag{4.6}
\]

Splitting at \(r=\varepsilon\) gives

\[
\begin{aligned}
 L_\varepsilon
 &=\int_0^\varepsilon e^{-r}{r\over\varepsilon+r}\,dr
   +\varepsilon\int_\varepsilon^\infty {e^{-r}\over\varepsilon+r}\,dr\\
 &=\varepsilon\log(1/\varepsilon)+O(\varepsilon).        \tag{4.7}
\end{aligned}
\]

After division by the total mass \(1+\varepsilon\), the lifted refined
uncertainty is therefore of order
\(\varepsilon\log(1/\varepsilon)\).  The ratio between (4.7) and (4.5)
tends to infinity.  This example is a legitimate conditional line family:
conditional thresholds in a globally balanced joint law need not themselves
be medians.  A proof through (4.1) must exploit global compatibility of the
levels; a line-by-line constant comparison is false.

## 5. A sharply tilted exponential wedge

The bivariate theorem survives the most asymmetric elementary cone model.
Let \(A\) and \(B\) be independent exponentials of rates
\(\varepsilon\) and \(1\), respectively, where \(0<\varepsilon<1\).  Put

\[
                         T=A+B,\qquad Y=A-B.               \tag{5.1}
\]

The law of \((Y,T)\) is log-affine on the wedge
\(\{t\ge|y|\}\).  Let \(h_\varepsilon\) be the median of \(T\), and cut
at \(T=h_\varepsilon\).  The survival function of \(T\) is

\[
 \mathbb P(T>t)={e^{-\varepsilon t}-\varepsilon e^{-t}\over1-\varepsilon},
                                                                    \tag{5.2}
\]

so

\[
                         \varepsilon h_\varepsilon\longrightarrow\log2.
                                                                    \tag{5.3}
\]

The density of \(Y\) is

\[
 f_Y(y)={\varepsilon\over1+\varepsilon}
 \begin{cases}
  e^{-\varepsilon y},&y\ge0,\\
  e^y,&y<0.
 \end{cases}                                               \tag{5.4}
\]

Conditional on \(Y=y\ge0\), one has

\[
 B\mid Y=y\sim\operatorname {Exp}(1+\varepsilon),qquad
 T=y+2B.                                                    \tag{5.5}
\]

For \(y<0\), the analogous formula is
\(A\mid Y=y\sim\operatorname {Exp}(1+\varepsilon)\) and
\(T=-y+2A\).  Thus every conditional standard deviation is

\[
                         \sigma_Y={2\over1+\varepsilon}.   \tag{5.6}
\]

For \(|y|<h_\varepsilon\), put \(s=h_\varepsilon-|y|\).  The conditional
lower-tail mass is

\[
 p_y=1-e^{-(1+\varepsilon)s/2},qquad
 q_y=\min(p_y,1-p_y),                                      \tag{5.7}
\]

and it is zero for \(|y|\ge h_\varepsilon\).  The positive-\(y\)
contribution to the uncertainty is exactly

\[
 \mathcal U_+(\varepsilon)
 ={\varepsilon\over2}e^{-\varepsilon h_\varepsilon}
   \int_0^{h_\varepsilon}e^{\varepsilon s}
   \min\bigl(1-e^{-(1+\varepsilon)s/2},
                    e^{-(1+\varepsilon)s/2}\bigr)\,ds.    \tag{5.8}
\]

The negative-\(y\) contribution is nonnegative and is in fact
exponentially smaller.  Dominated convergence, (5.3), and

\[
 \int_0^\infty\min(1-e^{-s/2},e^{-s/2})\,ds=2\log2        \tag{5.9}
\]

give

\[
                         {\mathcal U(\varepsilon)\over\varepsilon}
                         \longrightarrow {\log2\over2}.   \tag{5.10}
\]

Finally,

\[
                         \operatorname {Var}T={1\over\varepsilon^2}+1,
                                                                    \tag{5.11}
\]

and hence

\[
 \boxed{
 \mathcal U(\varepsilon)\sqrt{\operatorname {Var}T}
 \longrightarrow {\log2\over2}.}                         \tag{5.12}
\]

Thus arbitrarily asymmetric conditional families do occur, but the narrow
transition layer in the base has exactly the reciprocal width needed by
(CT).  This model is sharp in scale and is not a counterexample.

## 6. Canonical stress tests

### 6.1 Gaussian

For \(\gamma_n\) and \(E=\{x_1\le0\}\), put
\(a=|\theta_1|\).  The line conditional standard deviation is one and the
exact calculation is

\[
                         \mathcal U_\theta(E)={\arcsin a\over\pi}.   \tag{6.1}
\]

Therefore

\[
 {a\over\pi}\le\mathcal U_\theta(E)\le {a\over2},
 \qquad
 \mathcal U(E)\asymp n^{-1/2}.                             \tag{6.2}
\]

### 6.2 Cube

For the isotropic cube \([ -\sqrt3,\sqrt3]^n\) and its central coordinate
cut, the exact reach calculation in `random_line_uncertainty.md` gives,
for every \(\theta\),

\[
                         {|\theta_1|\over48}
 \le\mathcal U_\theta(E)\le {|\theta_1|\over2}.           \tag{6.3}
\]

This example is important for normalization: a typical random chord has
standard deviation of order \(n^{-1/2}\), and division by that standard
deviation restores the boundary scale.

### 6.3 Regular simplex and tilted cuts

Every simplex becomes the regular simplex after isotropic normalization.
For the balanced cap parallel to a facet, barycentric exponential races
give

\[
                         c\,\mathbb E|\theta_1|
 \le\mathcal U(E)\le C\,\mathbb E|\theta_1|.              \tag{6.4}
\]

The proof, including the deterministic truncation caused by the cap
coordinate, is in Section 7 of `random_line_uncertainty.md`.  The same race
calculation is stable under a fixed tilt of the cutting hyperplane: the
derivative vector in barycentric coordinates is merely replaced by its
orthogonal projection onto the zero-sum subspace.  When one coefficient is
large, the remaining exponential clocks have total rate of the same order;
when no coefficient is large, the two sign rates are comparable.  Thus no
tilted-simplex degeneration is visible.  A uniform proof over every tilt
would, however, be another instance of (CT), and is not asserted here.

### 6.4 Right circular cone

For a right circular cone with axial coordinate \(S\), isotropic scaling
has height \(\Theta(n)\) and base radius \(\Theta(\sqrt n)\).  The median
axial section is at

\[
                         S_0=H,2^{-1/n}=H-\Theta(1).      \tag{6.5}
\]

A point uniform on this section lies at radial deficit
\(\Theta(n^{-1/2})\), and its chord in a random tangent direction has
two reaches of order one with a fixed probability.  For a random ambient
direction, the axial component is \(\Theta(n^{-1/2})\); over such an
order-one tangent chord the axial displacement is \(O(n^{-1/2})\), while
the base truncation is at axial distance \(\Theta(1)\).  Hence the two
reaches remain comparable on a fixed-probability event, and (3.3) gives

\[
                         \mathcal U(E)\asymp n^{-1/2}.     \tag{6.6}
\]

This also shows why a cone does not realize the pure-fiber obstruction at
the wrong scale: isotropic normalization makes its side slope exactly
comparable to a typical random axial component.

### 6.5 Product exponentials

Let \(X_i=E_i-1\), where the \(E_i\) are independent unit exponentials,
and cut at the median of \(X_1\).  On a line through a boundary point,
the density is log-affine and the support endpoints are exponential races:

\[
 \tau_+=\min_{\theta_i<0}{E_i\over|\theta_i|},\qquad
 \tau_- =\min_{\theta_i>0}{E_i\over\theta_i},              \tag{6.7}
\]

up to the deterministic, much more distant, cap from coordinate one.
For a uniform random \(\theta\), with fixed positive probability,

\[
 \sum_{i\ge2,\theta_i>0}\theta_i\asymp\sqrt n,qquad
 \sum_{i\ge2,\theta_i<0}|\theta_i|\asymp\sqrt n,qquad
 \left|\sum_i\theta_i\right|=O(1).                       \tag{6.8}
\]

On that event the two independent races have comparable rate
\(\Theta(\sqrt n)\), while the log-affine tilt changes the density by only
a universal factor over the resulting chord.  Thus the boundary-normalized
ratio in (1.5) is bounded below with fixed probability.  Since the median
density of \(X_1\) is \(1/2\), averaging the factor \(|\theta_1|\) gives

\[
                         \mathcal U(E)\asymp n^{-1/2}.     \tag{6.9}
\]

The two-dimensional extreme-rate limit of this calculation is exactly the
tilted wedge in Section 5.

## 7. Lower-dimensional support

If \(\mu\) is supported on an affine subspace \(L\) of dimension \(k\ge1\),
identify \(L\) isometrically with \(\mathbb R^k\).  The covariance on
\(L\) is positive definite after quotienting by any further null direction;
whitening is performed intrinsically in \(L\).  The halfspace is replaced
by its intersection with \(L\), and random directions are sampled from
\(S(L)\), not from the ambient sphere.  Formulas (1.3)--(1.10) then hold
with \(n\) replaced by \(k\).  A one-dimensional nonpoint law has only the
two line directions and its uncertainty at a median is

\[
                         {1\over2\sqrt{\operatorname {Var}T}},       \tag{7.1}
\]

so the base case is immediate.  A point mass is excluded by convention.

## 8. Concentrated normal matrices

Suppose a finite-perimeter interface has normal matrix

\[
 Q={1\over P_\mu(E)}\int_{\partial^*E}N_E\otimes N_E\,dP_\mu
\]

and

\[
                         \operatorname {tr}((I-u\otimes u)Q)\le\delta^2.
                                                                    \tag{8.1}
\]

Then all but an \(O(\delta^2/\eta^2)\) fraction of its perimeter has normal
within projective angle \(\eta\) of \(u\).  This is only a local statement.
It does not imply that the set is close in measure to one halfspace.  A
union of many parallel slabs has \(Q=u\otimes u\) exactly while its
symmetric difference from every halfspace can be bounded away from zero.

Such oscillatory examples have *more* random-line uncertainty, so they do
not disprove a concentrated-normal extension.  They do show that no bound
of the form

\[
 \mu(E\triangle H)\le C\delta^\alpha                         \tag{8.2}
\]

can be used to transfer the halfspace result.  A valid extension would need
a dichotomy: either the parallel sheets create extra line mixing, or a
single-sheet component carries nearly all the phase mass and can be compared
to a halfspace.  Neither conclusion follows from (8.1) alone in the present
audit.  Consequently the halfspace theorem, even if (CT) is proved, does not
yet extend quantitatively to concentrated normal matrices without an
additional global sheet-count or phase-ordering lemma.

## 9. Precise remaining statement

All dimension dependence has disappeared from the following assertion.

> **Bivariate conditional-threshold lemma.**  There is a numerical
> \(c_2>0\) such that, for every nondegenerate log-concave probability
> \(\eta\) on \(\mathbb R^2\), every linear coordinate \(T\) whose median
> is zero, and the complementary coordinate \(Y\),
> \[
>  \int {\min(\eta(T\le0\mid Y=y),\eta(T>0\mid Y=y))
>             \over\sqrt{\operatorname {Var}(T\mid Y=y)}}
>       \,d\eta_Y(y)
>  \ge {c_2\over\sqrt{\operatorname {Var}_\eta T}}.
> \]

The tilted wedge proves that both sides can be of the same arbitrarily small
order before variance normalization.  The plateau--exponential example
proves that uniform-chord epigraph refinement can exceed the original left
side by an unbounded logarithm.  The pure-fiber calculation proves that
half-mass alone does not balance the two concave reach functions.  Any proof
of the lemma must overcome both effects simultaneously; neither may be
discarded as a technical approximation term.
