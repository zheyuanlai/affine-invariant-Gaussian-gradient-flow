# Audit of the equimeasurable matrix replacement

## 0. Verdict

The soft-thinning argument and the direct-method construction in the source
note are sound after two qualifications:

1. weighted variation has to be interpreted as relative weighted variation
   on the relative interior of the convex support; and
2. at a nonsmooth nuclear-norm term one cannot simply choose an arbitrary
   subgradient.  In a smooth first-variation problem there **exists** one
   compatible constant subgradient, by a minimax argument.

The numerical situation is better than stated there.  If

\[
 r={\Delta_F\over T},\qquad T=\operatorname {tr}M(F),
\]

then the normalized matrices satisfy

\[
 \|Q'-Q\|_*\le r(1+\kappa^{-1}),
 \qquad
 \operatorname {tr}(Q'^2)-\operatorname {tr}(Q^2)
 \le e+{e^2\over2},
 \quad e=r(1+\kappa^{-1}).                         \tag{0.1}
\]

This removes a factor four from the first-order purity estimate in the
original note.  Moreover the fixed-scale construction can be retuned,
without a dimension-dependent choice, to

\[
 \alpha=10^{-14},\qquad \beta\le10^{-8}.
\]

For `K>=2*10^26` this gives

\[
 {\Delta_F\over T}<9.525\cdot10^{-5}<2.86\cdot10^{-4},
 \qquad 1-\operatorname {tr}Q^2>.0053349>.0045.    \tag{0.2}
\]

Thus the requested simultaneous retuning is possible.  The complementary
case `K<2*10^26` is already a universal Poincare bound, so the very large
threshold is harmless for a dimension-free contradiction argument.

There is also a genuine geometric advance.  For a smooth minimizer with a
non-atomic regular value law, exact quantile restoration proves

\[
 -\operatorname {div}_\mu D\Phi_H(\nabla G)=\lambda(G)          \tag{0.3}
\]

for one constant matrix `H`, `||H||op<=1`.  Hence almost every smooth level
has constant weighted anisotropic mean curvature.  Wulff-normal rays
`x+tD Phi_H(n)` obey the same log-concave killed-flux identity as Euclidean
normal rays.

What does **not** follow is that almost every level is a global
fixed-volume anisotropic perimeter minimizer.  Section 8 gives an explicit
uniform-rectangle counterexample.  The unresolved step is therefore
narrower but still load bearing: extend (0.3) and the Wulff killed-tube
formula, with usable curvature multipliers, from smooth regular minimizers
to the actual `BV` minimizer in arbitrary dimension and at hard support.

## 1. Soft thinning and the constant

Let `P=Lambda/T`, let `0<=omega<=1`, and put

\[
 t=\int\omega\,d\Lambda,
 \qquad P_s={\omega\Lambda\over t}.
\]

Introduce a Bernoulli variable `I` with
`P(I=1|x)=omega(x)`.  For every Hilbert-valued `Z`, total variance gives

\[
 \operatorname {Var}_P Z
 \ge {t\over T}\operatorname {Var}_{P_s}Z.          \tag{1.1}
\]

This proof requires neither a positive floor nor a derivative of `omega`.
For `Z=theta theta^T`, `|Z|_{HS}=1`, and therefore

\[
 \operatorname {Var}_P Z=1-\operatorname {tr}Q^2.
\]

The central selected matrix in the frozen run obeys

\[
 t>.0048495p,\qquad \|Q_s\|_{op}<1/17.86,
 \qquad T\le p.
\]

Since `Q_s` is positive semidefinite with trace one,
`tr(Q_s^2)<=||Q_s||op`.  Hence

\[
 1-\operatorname {tr}Q^2
 >.0048495\left(1-{1\over17.86}\right)
 =.00457797144\ldots .                              \tag{1.2}
\]

The decimal `.0045779` in the source note is therefore correct.

## 2. Weighted `BV` lower semicontinuity and existence

Let `E` be the affine hull of the measure, let `K` be its convex support in
`E`, and write `rho=e^{-V}` on `K^o=ri(K)`.  Weighted perimeter means
relative perimeter on `K^o`; a portion of a reduced boundary lying on
`partial K` has zero cost, exactly as required by exterior Minkowski
content for a measure supported on `K`.

For `G in BV_loc(K^o)` define

\[
 |DG|_\mu=\rho |DG|,
 \quad DG=\sigma_G|DG|,
 \quad M(G)=\int_{K^o}\sigma_G\sigma_G^T\,d|DG|_\mu.            \tag{2.1}
\]

For a symmetric `H` with `||H||op<=1`, put

\[
 \Phi_H(\xi)=|\xi|+\kappa{\xi^TH\xi\over|\xi|},
 \qquad \Phi_H(0)=0.                               \tag{2.2}
\]

At a unit vector `n` and for `h perpendicular n`,

\[
 D^2\Phi_H(n)[h,h]
 =|h|^2+\kappa\{2h^THh-(n^THn)|h|^2\}.             \tag{2.3}
\]

Homogeneity gives `D^2 Phi_H(n)n=0`; thus there are no omitted radial
cross-terms.  Consequently

\[
 (1-3\kappa)|h|^2
 \le D^2\Phi_H(n)[h,h]
 \le(1+3\kappa)|h|^2.                              \tag{2.4}
\]

For `0<kappa<1/3`, `Phi_H` is convex, positive, one-homogeneous, and
uniformly elliptic on tangent spaces.  Nuclear/operator duality gives the
exact identity

\[
 \begin{split}
 \operatorname {TV}_\mu(G)
  +\kappa\|M(G)-M(F)\|_*
 =\sup_{\substack{H=H^T\\\|H\|_{op}\le1}}
 \left\{\int_{K^o}\rho\,\Phi_H(dDG)
       -\kappa\operatorname {tr}(HM(F))\right\}.
 \end{split}                                        \tag{2.5}
\]

On every compact subset of `K^o`, `rho` is continuous and bounded above
and below by positive constants.  Reshetnyak lower semicontinuity applies
to each integral in (2.5).  Exhaustion by compact subsets makes the global
relative integral lower semicontinuous, and a supremum of lower
semicontinuous functionals is lower semicontinuous.

For completeness, the direct-method compactness on an unbounded support is
as follows.  A minimizing sequence `G_j`, truncated to `[0,1]`, has bounded
weighted variation.  On each compact subset of `K^o` this bounds ordinary
variation, so a diagonal subsequence converges locally in `L^1`.  Since
`|G_j-G|<=1` and `mu` is tight,

\[
 \int_{K\setminus L}|G_j-G|d\mu\le\mu(K\setminus L)
\]

upgrades the convergence to `L^1(mu)`.  Convergence in probability and
equality of all the pushforward laws imply that the limit is still
equimeasurable with `F`.  Equation (2.5) supplies lower semicontinuity.
Thus a minimizer exists.  This argument also handles an extended-valued
convex potential directly on `ri(dom V)`; no assertion of a positive
global lower bound for the density is needed.

Weighted coarea applies componentwise to the even one-homogeneous matrix
integrand and gives

\[
 M(G)=\int_{\mathbb R}M(\{G>r\})dr,
 \qquad \operatorname {tr}M(G)=\operatorname {TV}_\mu(G).       \tag{2.6}
\]

## 3. The exact central deficit

Let `F_0=T_s1_S`, and let `C=[r_-,r_+]` be the interval on which

\[
 \delta\le\mu(F_0>r)\le1-\delta,
 \qquad \delta=10^{-4}.
\]

Let `F=h_C(F_0)` be the monotone clipping with derivative `1_C`.  For
almost every `r in C`, `{F>r}={F_0>r}`; outside `C` its superlevel is either
the whole support or the empty set.  Since `I_mu(0)=I_mu(1)=0`, weighted
coarea gives the **identity**

\[
 \begin{split}
 \Delta_F
 &=\operatorname {TV}_\mu(F)
   -\int_{\mathbb R}I_\mu(\mu(F>r))dr\\
 &=\int_C\left[P_\mu(F_0>r)
       -I_\mu(\mu(F_0>r))\right]dr.
 \end{split}                                        \tag{3.1}
\]

By the definition of the Cheeger constant,

\[
 I_\mu(v)\ge\psi\min(v,1-v).
\]

Thus the exact identity (3.1), rather than the finite-splice functional,
implies

\[
 0\le\Delta_F\le\mathcal D_{co}(F_0)
 \le p\left\{\beta+{4\sqrt\alpha\over
              c_G(1-\beta)}\right\}.              \tag{3.2}
\]

The source note's numerical bound is correct, but the logical chain is
`Delta_F <= D_co`; finite physical splicing is not needed for (3.2).

## 4. Matrix retention and the sharper purity estimate

Write

\[
 M=M(F),\quad M'=M(G_\kappa),\quad
 T=\operatorname {tr}M,\quad T'=\operatorname {tr}M',\\
 Q={M\over T},\quad Q'={M'\over T'}.
\]

Comparison with `F` and the profile lower bound give

\[
 T'\le T,
 \qquad T-T'\le\Delta_F,
 \qquad \|M'-M\|_*\le{\Delta_F\over\kappa}.        \tag{4.1}
\]

Assume `Delta_F<T`, which is the only case used below.  The normalization
estimate can be written with denominator `T`, not `T'`:

\[
 \begin{split}
 \|Q'-Q\|_*
 &\le {\|M'-M\|_*\over T}
   +\left\|M'\left({1\over T'}-{1\over T}\right)\right\|_*\\
 &={\|M'-M\|_*\over T}+{|T-T'|\over T}
 \le {\Delta_F\over T}(1+\kappa^{-1}).
 \end{split}                                        \tag{4.2}
\]

This proves the first part of (0.1) without the factor two and without the
assumption `Delta_F<T/2`.

There is a second improvement.  Put `X=Q'-Q`; then `tr X=0`.  If
`e=||X||_*`, the positive and negative parts of `X` both have trace `e/2`,
so

\[
 \operatorname {tr}X^2\le {e^2\over2}.             \tag{4.3}
\]

Also, subtracting the midpoint of the extreme eigenvalues of `Q` and using
`0<=Q<=I` gives

\[
 2\operatorname {tr}(QX)\le e.                    \tag{4.4}
\]

Therefore

\[
 \operatorname {tr}(Q'^2)-\operatorname {tr}(Q^2)
 =2\operatorname {tr}(QX)+\operatorname {tr}X^2
 \le e+{e^2\over2}.                               \tag{4.5}
\]

Equations (4.2)--(4.5) are dimension free.

## 5. A fixed retuning which passes both tests

Take

\[
 \boxed{\alpha=10^{-14},\qquad\beta\le10^{-8},
        \qquad\eta_{floor}=10^{-5},\qquad\delta=10^{-4}.}       \tag{5.1}
\]

The symbolic hypotheses of the fixed-scale report hold because
`alpha<alpha_*` and `beta<sqrt(alpha)`.  Its constant

\[
 C_*={8\over c_0}\left(1+{8\sqrt2\over c_G}\right)
 =1217.591\ldots
\]

gives

\[
 \tau\le8C_*\sqrt\alpha<.000975.                  \tag{5.2}
\]

The posterior cutoff is

\[
 M_\alpha=\sqrt{2\log(32/\sqrt\alpha)}
 =6.258407\ldots .                                  \tag{5.3}
\]

We next recheck the constants `a_-` and `a_+`; the old decimal cannot just
be reused after changing `alpha`.  For `r=Phi(-x)<=.1`, Mills' inequality
and `x>=Phi^{-1}(.9)>1` give

\[
 {r(1-r)\over I(r)}
 \ge .9{x\over1+x^2}
 \ge .9{M_\alpha\over1+M_\alpha^2}
 >.14022.
\]

For `.1<=r<=.5`, the ratio is at least `.09/I_0>.225`; symmetry handles
the other half.  Hence

\[
 a_->.14,
 \qquad a_+\le1/c_G=\sqrt{\pi/2}<1.26.             \tag{5.4}
\]

The physical good-flux matrix has effective rank at least

\[
 {a_-\sqrt{1-\tau}\over a_+}
 {\alpha^{3/2}K\over8\pi}.
\]

For `K>=2*10^26`, this is larger than `883`, and in particular larger than
the conservative value `500` used below.  Furthermore

\[
 b\ge b_0p,qquad
 b_0={.14(1-10^{-8})\over8\pi}>.0055704.            \tag{5.5}
\]

Both the coarea deficit and the alignment deficit are at most

\[
 \epsilon_\alpha p,qquad
 \epsilon_\alpha
 :=10^{-8}+{4\cdot10^{-7}\over c_G(1-10^{-8})}
 <5.114\cdot10^{-7}.                               \tag{5.6}
\]

Apply the floored analytic selector and then delete the volume-endpoint
levels.  The endpoint deletion has trace at most
`(2delta+epsilon_alpha)p`.  The resulting central selected matrix
`M_s` satisfies, with `q=b/p>=b_0`,

\[
 \begin{split}
 {\operatorname {tr}M_s\over p}
 &\ge(1-10^{-5})(q-\epsilon_\alpha)
       -(2\delta+\epsilon_\alpha),\\
 {\|M_s\|_{op}\over p}
 &\le10^{-5}+(1-10^{-5})
       \left({2q\over500}+4\epsilon_\alpha\right).
 \end{split}                                        \tag{5.7}
\]

The ratio in (5.7) is worst at the lower endpoint `b=b_0`: the numerator
and denominator have the form `Aq-C` and `Dq+E`, with all four constants
positive.  Its derivative is `(AE+DC)/(Dq+E)^2>0`.  Thus

\[
 {\operatorname {tr}M_s\over p}>.0053693,
 \qquad
 {\operatorname {tr}M_s\over\|M_s\|_{op}}>156.4.  \tag{5.8}
\]

Let `T` be the unweighted central coarea mass.  Heat contraction and
thinning give

\[
 \operatorname {tr}M_s\le T\le p.
\]

Equations (1.1), (5.7), and (5.8) now yield

\[
 1-\operatorname {tr}Q^2
 >.0053693(1-1/156.4)>.0053349.                    \tag{5.9}
\]

On the other hand, (3.2), (5.6), and `T>=tr M_s` give

\[
 {\Delta_F\over T}
 <{5.114\cdot10^{-7}\over.0053693}
 <9.525\cdot10^{-5}.                               \tag{5.10}
\]

This proves (0.2).  As one concrete choice for the replacement functional,
take `kappa=1/4`.  Then (4.2) gives `e<4.763*10^{-4}`, and (4.5) shows that
the output normal law still has angular variance larger than `.004858`.
The anisotropy remains uniformly elliptic, with ellipticity ratio bounded
using `1/4<=D^2 Phi_H<=7/4` on tangent spaces.

## 6. The compatible subgradient and exact distribution restoration

The sentence "choose a nuclear-norm subgradient" in the source note is too
strong.  An arbitrary subgradient need not give stationarity.  The correct
smooth statement is the following.

**Proposition 6.1 (smooth equimeasurable Euler equation).**  Suppose `rho`
is `C^2` and positive on an open set, `G` is a smooth local minimizer, and
on a value interval `J` one has `|grad G|>0` and the law of `G` has a
positive `C^1` density.  Then there is one matrix

\[
 H\in\partial\|M(G)-M(F)\|_*,\qquad \|H\|_{op}\le1,             \tag{6.1}
\]

and a measurable scalar function `lambda` such that

\[
 -\operatorname {div}_\mu D\Phi_H(\nabla G)=\lambda(G)          \tag{6.2}
\]

on `{G in J}` in the distributional sense.  In particular, every regular
level in `J` has constant weighted `Phi_H`-anisotropic mean curvature.

**Proof.**  First construct exact admissible curves.  Let `u` be smooth and
compactly supported in `{G in J}`, put `H_t=G+tu`, and let `C_t` be the CDF
of `H_t`.  If `C_0^{-1}` denotes the quantile map, then

\[
                  G_t=C_0^{-1}(C_t(H_t))            \tag{6.3}
\]

has exactly the same law as `G` for every sufficiently small `t` (with the
usual generalized-inverse convention).  If

\[
 \bar u(r)=E[u\mid G=r],
\]

coarea differentiation gives

\[
 \left.{dG_t\over dt}\right|_{t=0}=u-\bar u(G).    \tag{6.4}
\]

Indeed, if `w` is the density of the law of `G`, then
`partial_t C_t(r)|_0=-w(r)bar u(r)`; differentiating (6.3) proves (6.4).
Thus the tangent space contains every smooth `u` satisfying
`E[u|G]=0`.

Let `X=M(G)-M(F)`, let `S=partial||X||_*`, let `A(u)` be the first
variation of `M`, and let `L_0(u)` be the first variation of weighted total
variation.  Minimality in both signs gives

\[
 L_0(u)+\kappa\sup_{H\in S}\operatorname {tr}(HA(u))\ge0
\]

for every tangent `u`.  Apply the elementary minimax theorem to the
bilinear function

\[
 (H,u)\longmapsto \kappa\operatorname {tr}(HA(u))+L_0(u)
\]

on the compact convex set `S` and the linear tangent space.  Since `u=0`
is available and the tangent space contains `u` and `-u`, there is a
single `H in S` for which

\[
 L_0(u)+\kappa\operatorname {tr}(HA(u))=0           \tag{6.5}
\]

for every tangent `u`.  This is the required compatible, rather than
arbitrary, subgradient.

For this `H`, the left side is the first variation of

\[
 \int\rho\Phi_H(\nabla G)dx.
\]

Integration by parts and (6.4) show that
`-div_mu D Phi_H(grad G)` is orthogonal to all functions with zero
conditional mean given `G`.  It is therefore a function of `G`, which is
(6.2).  Homogeneity makes `D Phi_H(grad G)=D Phi_H(n)`, so (6.2) is exactly
the constant anisotropic mean-curvature equation on each regular level.
QED.

The same calculation can be obtained from the differentiable nuclear
regularization

\[
 X\longmapsto\operatorname {tr}\sqrt{X^2+\varepsilon^2I}-n\varepsilon.
\]

Its gradient has operator norm at most one, and its comparison estimate
loses only `n epsilon`; taking `epsilon` to zero at fixed `n` recovers a
compatible subgradient.  This is useful as a construction, but it does not
by itself supply uniform regularity of the limiting `BV` minimizer.

## 7. Wulff-normal killed flux

The constant anisotropy in Proposition 6.1 is precisely rigid enough for a
killed-tube calculation.  Let `Sigma` be a `C^2` level, let `n` be its
chosen unit normal, and put

\[
 z=D\Phi_H(n),\qquad B=D^2\Phi_H(n)|_{T\Sigma}.
\]

Let `S=D_\Sigma n` be the Euclidean shape operator and define the Wulff ray

\[
                         X_t(x)=x+t z(x).            \tag{7.1}
\]

Kill the ray at its first contact with the support, its first focal time,
or the first loss of injectivity of the Wulff normal map.  Before that
time,

\[
 D_\Sigma X_t=I+tBS,
 \qquad
 \det[D_\Sigma X_t,z]
 =\Phi_H(n)\det(I+tBS).                             \tag{7.2}
\]

The second identity uses Euler's formula `z dot n=Phi_H(n)`.  Relative to
the initial anisotropic surface measure

\[
 d\sigma_{\Phi,\mu}=\Phi_H(n)e^{-V(x)}d\mathcal H^{n-1}(x),
\]

the flux Jacobian is

\[
 j_x(t)=\det(I+tBS)
       \exp\{-V(x+tz)+V(x)\}.                       \tag{7.3}
\]

Although `BS` need not be symmetric, it is similar to the symmetric matrix
`B^{1/2}SB^{1/2}`.  Consequently

\[
 {d^2\over dt^2}\log j_x(t)
 =-\operatorname {tr}\left[((I+tBS)^{-1}BS)^2\right]
  -z^T\nabla^2V(x+tz)z\le0.                        \tag{7.4}
\]

The initial logarithmic slope is

\[
 h_{\Phi,\mu}
 =\operatorname {tr}(BS)-\nabla V\mathbin\cdot z,
\]

the weighted anisotropic mean curvature.  If it equals the constant
`lambda` on the level, then the exact Taylor formula reads

\[
 j_x(t)=\exp\{\lambda t-D_x(t)\},                  \tag{7.5}
\]

where

\[
 D_x(t)=\int_0^t(t-s)\left\{
 \operatorname {tr}[((I+sBS)^{-1}BS)^2]
 +z^T\nabla^2V(x+sz)z\right\}ds\ge0.              \tag{7.6}
\]

Thus

\[
 R(t)=\int_\Sigma1_{\{t<\tau(x)\}}e^{-D_x(t)}
                 d\sigma_{\Phi,\mu}(x)             \tag{7.7}
\]

is nonincreasing.  This is the exact anisotropic analogue of the killed
Euclidean normal-flux formula.  It uses stationarity, not global
isoperimetric minimality.

For `kappa=1/4`, one also loses no qualitative rank information by changing
from `n` to the ray direction `z/|z|`.  Indeed

\[
 {1\over4}I\le B\le{7\over4}I,
 \quad {3\over4}\le\Phi_H(n)\le{5\over4},
 \quad {3\over4}\le|z|\le{7\over4}.               \tag{7.8}
\]

The odd map `n -> z/|z|` is a diffeomorphism of the sphere and its inverse
has differential norm at most `49/3`.  One way to see the lower differential
bound is that, for `h perpendicular n`, `Bh perpendicular n` and

\[
 \left|D\left({z\over|z|}\right)h\right|
 ={ |P_{z^\perp}Bh|\over|z|}
 \ge {3\over49}|h|.                                \tag{7.9}
\]

It descends to a bi-Lipschitz projective map.  Since projector variance is
one half the mean squared projective chord distance, and the reweighting
`Phi_H(n)` lies between `3/4` and `5/4`, a fixed positive fraction of the
normal-projector variance survives in the Wulff-ray directions.  All of
these constants are dimension free.

## 8. Levelwise global minimality is false

The nesting qualification is real, not merely technical.  Consider the
uniform probability on the rectangle

\[
                         R_L=[0,L]\times[0,1],qquad L>2,
\]

with relative Euclidean perimeter.  For physical areas `a<1/2`, the planar
relative-isoperimetric classification gives

\[
 I(a)={1\over L}\min\{\sqrt{\pi a},1\}.            \tag{8.1}
\]

Below `a=1/pi` the minimizers are quarter disks at a corner; above it they
are vertical end strips.  This also follows directly from the standard
first-variation classification: a minimizing free curve is a circular arc
meeting the sides orthogonally or a straight cross-cut, after which the
quarter circle and the vertical segment are the only competitive curves
in this area range.

Take

\[
                         a_1=.31<1/\pi<a_2=.32.
\]

A minimizing quarter disk of area `a_1` has radius
`2sqrt(a_1/pi)>.62`, whereas a minimizing end strip of area `a_2` has width
`.32`.  No choice among the four quarter disks and two end strips makes the
first set a subset of the second.

Now prescribe the three-atom law of a function taking value `1` on a set of
area `a_1`, value `1/2` on an additional area `a_2-a_1`, and value zero
elsewhere.  Equimeasurable functions are exactly nested pairs

\[
                         E_1\subset E_2,qquad |E_i|=a_i,
\]

and coarea gives

\[
 \operatorname {TV}_\mu(G)
 ={1\over2}\{P_\mu(E_1)+P_\mu(E_2)\}.              \tag{8.2}
\]

The direct method gives a minimizing nested pair.  If both of its levels
were individual isoperimetric minimizers, they would be one of the
incompatible pairs just classified, a contradiction.  Thus an
equimeasurable total-variation minimizer need not have globally minimizing
levels, even for a bounded log-concave probability and even when the
anisotropy is Euclidean (`H=0`).  One level fails on an interval of level
values of length `1/2`, so this is not an exceptional-level issue.

The example does not contradict Proposition 6.1: a compromise foliation
can have constant-mean-curvature levels without making each one globally
isoperimetric.  It shows exactly why independent level replacement is not
available and why nested phase incidence must remain part of the argument.

## 9. Remaining load-bearing statements

The audit leaves the following precise list.

1. Prove a `BV` version of Proposition 6.1 for the actual minimizer,
   including a single compatible constant `H`, generalized anisotropic
   mean curvature on almost every reduced boundary, hard-support contact,
   and a singular set negligible for the tube calculation.
2. Control the multiplier `lambda(r)` on levels carrying the retained
   matrix.  In a smooth foliation, differentiation gives
   `dP_Phi/dv=lambda` up to orientation, but the integrated scalar deficit
   alone does not control this derivative.
3. Pass the Wulff killed-flux formula through the singular and unbounded
   approximations with the same constants.
4. Use the retained cross-level incidence.  The single function `G` keeps
   the family nested and Proposition 6.1 uses one anisotropy for all levels,
   but neither fact alone converts integrated angular variance into a
   single long-surviving ranked tube packet.

The fixed-scale ratio and the selector derivative are no longer numerical
obstructions.  The unresolved obstruction is now geometric regularity and
the conversion of nested anisotropic CMC incidence into the ranked
killed-tube alternative.
