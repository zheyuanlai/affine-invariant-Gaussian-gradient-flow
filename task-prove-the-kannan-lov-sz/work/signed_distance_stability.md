# Second variation of an extremal signed-distance potential

This note verifies the normal-height second variation suggested by the ray
picture.  The statement is first given in a setting where every
differentiation is classical.  It then identifies exactly which term is
missing at focal or polyhedral junctions.

## 1. Smooth tubular setting

Let `Sigma` be an oriented `C^3` hypersurface with unit normal `N` and
shape operator `S=D_Sigma N`.  Assume that a neighborhood containing the
support under consideration is represented injectively by

\[
 F(y,t)=y+tN(y),\qquad f(F(y,t))=t,                       \tag{1}
\]

and that every `I+tS_y` is invertible there.  Let

\[
 d\mu=\int q_y(t)dt\,d\eta(y)                            \tag{2}
\]

be the ray disintegration.  Assume

\[
 \int_{t>0}q_y(t)dt=p,qquad
 \int_{t<0}q_y(t)dt=1-p=:\bar p                         \tag{3}
\]

for `eta`-almost every `y`, with `p in (0,1)`.  Shift `f` by a
constant so that `mu f=0`; then its zero surface is `Sigma`.

For `h in C_c^2(Sigma)`, deform the zero surface by

\[
 \Sigma_\varepsilon=\{y+\varepsilon h(y)N(y):y\in\Sigma\}, \tag{4}
\]

and let `f_epsilon` be its oriented signed distance.  For small
`epsilon`, assume the relevant tube remains free of a cut locus.  Dots below
denote derivatives at zero evaluated at a fixed ambient point.

**Lemma 1 (shape derivatives).**  At `x=F(y,t)`,

\[
 \dot f(x)=-h(y),                                        \tag{5}
\]

and

\[
 \ddot f(x)=-\int_0^t
  \left|(I+sS_y)^{-1}\nabla_\Sigma h(y)\right|^2ds.    \tag{6}
\]

**Proof.**  Differentiate the eikonal equation
`|grad f_epsilon|^2=1`.  Its first derivative gives

\[
 \langle\nabla f,\nabla\dot f\rangle=0,                  \tag{7}
\]

so `dot f` is constant on every normal ray.  Differentiating
`f_epsilon(y+epsilon hN)=0` at zero gives `dot f(y)=-h(y)`,
proving (5).

The second eikonal derivative is

\[
 \langle\nabla f,\nabla\ddot f\rangle=-|\nabla\dot f|^2. \tag{8}
\]

Since `D_yF=I+tS_y` on the tangent space, (5) gives

\[
 \nabla\dot f(F(y,t))
 =-(I+tS_y)^{-1}\nabla_\Sigma h(y).                   \tag{9}
\]

The second derivative of the boundary identity in the preceding paragraph
gives `ddot f(y)=0`: the other terms vanish because
`partial_N dot f=0` and `nabla^2f(N,N)=0`.  Integrating (8) from zero to
`t`, and using (9), proves (6).  QED.

## 2. Stability inequality from T3 extremality

Suppose now that `f` maximizes

\[
 \mathcal J(g)=\int|g-\mu g|d\mu                         \tag{10}
\]

over all one-Lipschitz `g`.  Every `f_epsilon` is one-Lipschitz, so
`J(epsilon):=mathcal J(f_epsilon)` has a local maximum at zero.

Set

\[
 \bar h=\int h\,d\eta,
 \qquad
 A_h(y,s)=|(I+sS_y)^{-1}\nabla_\Sigma h(y)|^2.          \tag{11}
\]

**Proposition 2 (signed-distance stability).**  Under the smooth tubular
hypotheses,

\[
 \boxed{\begin{aligned}
 &\int q_y(0)(h(y)-\bar h)^2d\eta(y)\\
 &\quad\le
 \int\!\left[
  \bar p\int_{0}^{\infty}q_y(t)\int_0^t A_h(y,s)dsdt
 +p\int_{-\infty}^{0}q_y(t)\int_t^0 A_h(y,s)dsdt
 \right]d\eta(y).
 \end{aligned}}                                               \tag{12}
\]

All integrals are restricted to the actual ray intervals.  The same formula
holds by approximation whenever the displayed quantities are integrable and
the deformed signed distances have no positive-`mu` cut-locus error.

**Proof.**  Centering gives

\[
 {d\over d\varepsilon}\bigl(f_\varepsilon-\mu f_\varepsilon\bigr)
 \big|_0=-(h-\bar h).                                     \tag{13}
\]

For a continuous density `q` and constants `a,b`, the distributional
identity `(|r|)''=2delta_0` gives

\[
 {d^2\over d\varepsilon^2}\int
 |t+\varepsilon a+\tfrac12\varepsilon^2b+o(\varepsilon^2)|q(t)dt
 \big|_0
 =2q(0)a^2+\int\operatorname{sgn}(t)bq(t)dt.             \tag{14}
\]

Apply (14) conditionally and then integrate in `eta`.  The term from
centering `ddot f` is

\[
 -(2p-1)\int\ddot f\,d\mu.                               \tag{15}
\]

Hence, using (3),

\[
 \begin{aligned}
 J''(0)
 &=2\int q_y(0)(h-\bar h)^2d\eta\\
 &\quad+2\int\left[
  \bar p\int_{t>0}\ddot f\,q_y(t)dt
  -p\int_{t<0}\ddot f\,q_y(t)dt
 \right]d\eta.
 \end{aligned}                                             \tag{16}
\]

Formula (6) makes the second line equal to minus twice the right-hand side
of (12).  Since `J''(0)<=0`, (12) follows.  QED.

The boundary term in (12) is strictly positive for every nonconstant
height.  It is the gain obtained by moving the zero of the absolute value.
The right side is the exact second-order loss required by the eikonal
constraint.

## 3. Translation-normal test and an exact trace cancellation

Take `h_a(y)=<a,N_y>`.  Then

\[
 \nabla_\Sigma h_a=S_ya_T,
 \qquad
 A_{h_a}(y,s)=|(I+sS_y)^{-1}S_ya_T|^2.                  \tag{17}
\]

Sum (12) over an orthonormal ambient basis `a=e_1,...,e_n`.  Put

\[
 m=\int N_y\,d\eta(y),qquad
 C_y(s)=\operatorname{tr}\big[S_y^2(I+sS_y)^{-2}\big].   \tag{18}
\]

The result is

\[
 \int q_y(0)|N_y-m|^2d\eta(y)
 \le\int R_y\,d\eta(y),                                  \tag{19}
\]

where

\[
 \begin{aligned}
 R_y={}&\bar p\int_0^\infty
     C_y(s)\,\mathbb P_y(T\ge s)ds\\
 &+p\int_{-\infty}^0
     C_y(s)\,\mathbb P_y(T\le s)ds.
 \end{aligned}                                             \tag{20}
\]

For a normal chart of `e^{-V}`, write `q_y=e^{-W_y}`.  In the sense of
second-derivative measures,

\[
 W_y''=V_{NN}(y+sN_y)+C_y(s).                             \tag{21}
\]

If `W_y` is differentiable at zero, there is an exact one-dimensional
identity

\[
 \begin{aligned}
 &\bar p\int_0^\infty W_y''(s)\mathbb P_y(T\ge s)ds\\
 &\quad+p\int_{-\infty}^0W_y''(s)\mathbb P_y(T\le s)ds
 =q_y(0).                                                  \tag{22}
\end{aligned}
\]

Indeed, integration by parts gives

\[
 \int_0^\infty W''(s)\mathbb P(T\ge s)ds
 =q(0)-pW'_+(0),                                          \tag{23}
\]

and

\[
 \int_{-\infty}^0 W''(s)\mathbb P(T\le s)ds
 =q(0)+\bar pW'_-(0);                                     \tag{24}
\]

the two slope terms cancel after the weights are applied.  If `W` has a
kink of size `a=W'_+(0)-W'_-(0)`, the two open-half-line integrals instead
sum, with the displayed weights, to `q(0)-p\bar p a`; this is even smaller.
Thus (21)--(22) imply

\[
 R_y\le q_y(0),                                           \tag{25}
\]

with the quantitative gap

\[
 \begin{aligned}
 q_y(0)-R_y={}&\bar p\int_0^\infty
 V_{NN}(y+sN_y)\mathbb P_y(T\ge s)ds\\
 &+p\int_{-\infty}^0
 V_{NN}(y+sN_y)\mathbb P_y(T\le s)ds                    \tag{26}
\end{aligned}
\]

when `V` is `C^2`.  For nonsmooth convex `V`, use its second-derivative
measure on the two open half-lines and add the nonnegative kink contribution
`p\bar p(V'_+(0)-V'_-(0))` when present.

Equation (25) shows both the value and the limitation of the smooth
translation-normal test.  After summing all ambient directions, its entire
shape-curvature cost is at most the boundary density, with no factor depending
on ray length or dimension.  Strict coercivity can come from the ambient
convexity term (26), but it vanishes for uniform convex bodies.  Therefore
the smooth trace estimate alone cannot close KLS; it is consistent with the
known improved-Lichnerowicz mechanism.

## 4. Singular junctions are exactly the omitted term

If the smooth pieces are flat (`S=0`), the right side of (12) vanishes for a
height which is constant on each piece.  Thus a globally admissible independent
translation of those pieces would contradict extremality unless all constants
agree.  A fan evades the classical calculation because its flat pieces meet
at medial or polyhedral junctions.  Under unequal piecewise-constant heights,
the cut locus sweeps a strip of thickness `Theta(epsilon)` while the two
envelopes differ there by `Theta(epsilon)`.  The resulting
`Theta(epsilon^2)` integral contributes to `J''(0)` but is absent from (6),
which was derived inside one fixed tubular chart.

Accordingly a general version of (12) must have the form

\[
 \text{boundary gain}
 \le \text{smooth ray-metric cost}+\text{singular junction cost}. \tag{27}
\]

The missing cost is a full-BV/focal measure on the completed ray graph.  It
must be nonnegative, must vanish for a global parallel foliation, and must
charge the vertex of the Gaussian fan.  Establishing a dimension-free lower
bound for that singular term in terms of orientation-packet variation would
provide precisely the weighted endpoint estimate missing from the midpoint
approach.  Formula (12) verifies the smooth part and fixes all coefficients;
it does not by itself construct or bound the singular completion term.

## 5. Exact envelope charge and the Gaussian-fan test

There is a local formula for the simplest medial junction.  Let `g_1,g_2`
be `C^1`, suppose `Gamma={g_1=g_2}` is regular, and let `h_1,h_2`
be constants (the same formula holds for smooth heights).  Put

\[
 G_\varepsilon(x)=\min\{g_1(x)-\varepsilon h_1,
                         g_2(x)-\varepsilon h_2\}.        \tag{28}
\]

Then, distributionally in `epsilon` and after integration against a smooth
weight `rho`,

\[
 {d^2\over d\varepsilon^2}\int G_\varepsilon\rho dx\big|_0
 =-\int_\Gamma {(h_1-h_2)^2\rho\over
                    |\nabla g_1-\nabla g_2|}\,d\mathcal H^{n-1}. \tag{29}
\]

Indeed,

\[
 \min(A,B)={A+B-|A-B|\over2},                            \tag{30}
\]

and `d^2|r-epsilon a|/d epsilon^2|_0=2a^2 delta_0(r)`.
The coarea formula turns `delta_0(g_1-g_2)` into (29).  Formula (29) is
negative, so after moving it to the right side of the maximality inequality
it is exactly a positive junction energy.

Apply this to the two-dimensional fan.  Index its `2m` oriented zero rays
cyclically.  In the sector between rays `j` and `j+1`,

\[
 |f|=\min\{\ell_j,\ell_{j+1}\},                           \tag{31}
\]

where the `ell`'s are Euclidean distances to the two boundary rays.  If the
sector angle is `alpha=pi/m`, their inward unit gradients have angle
`pi-alpha`, and therefore

\[
 |\nabla\ell_j-\nabla\ell_{j+1}|=2\cos(\alpha/2).         \tag{32}
\]

For a radial density, every medial bisector ray has the same weighted
one-dimensional measure `w`.  Piecewise-constant normal heights
`h_1,...,h_{2m}` consequently produce the singular cost

\[
 {w\over2\cos(\pi/(2m))}
   \sum_{j=1}^{2m}(h_j-h_{j+1})^2.                       \tag{33}
\]

The smooth cost is zero because every open fan chart is flat.  The boundary
gain is a comparable constant times
`sum_j(h_j-bar h)^2`.  Thus the singular stability operator is the cycle
Laplacian.  Its first nonzero eigenvalue is

\[
 4\sin^2(\pi/(2m))\asymp m^{-2}.                         \tag{34}
\]

For large `m`, (33) cannot dominate the boundary gain for a slowly varying
height.  This recovers, by a genuine second-variation calculation, the fact
that the Gaussian fan is not a global first-moment extremizer.  It also shows
that merely counting total turning at the common vertex is not the right
coercive quantity: the relevant object is the spectral gap of the entire
weighted medial-junction graph.

For a hypothetical bad extremizer made from many nearly flat patches, the
completed ray graph must therefore satisfy a weighted Poincare inequality for
normal heights at precisely the scale in (12).  Establishing that the bridge
overlap from global log-concavity forces such expansion would close the
singular branch.  The indicator midpoint count does not yet imply this graph
spectral estimate.
