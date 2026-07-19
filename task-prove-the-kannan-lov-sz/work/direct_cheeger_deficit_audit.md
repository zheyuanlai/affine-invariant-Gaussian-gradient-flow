# Clean-room audit of the direct Cheeger-deficit replacement

## 0. Verdict

The direct-deficit functional has a substantially stronger Euler theorem
than the preliminary discussion in Section 10 of
`cmc_near_cheeger_amplification.md` suggests.  The correct argument is not
to choose a nuclear-norm subgradient formally.  It has three convex-analytic
steps:

1. the augmented variation

   \[
   A_\kappa(G)={\rm TV}_\mu(G)
       +\kappa\|M(G)-M(F)\|_*
   \]

   is convex and lower semicontinuous when \(0<\kappa<1/3\);
2. the median deviation

   \[
   W(G)=\inf_{c\in\mathbb R}\int|G-c|\,d\mu
   \]

   is a continuous convex quotient norm; and
3. a convex interpolation followed by a genuine minimax argument produces
   one median subgradient \(q\) and one compatible, spatially constant
   matrix \(H\), \(\|H\|_{\rm op}\le1\), such that the direct-deficit
   minimizer also minimizes

   \[
       {\rm TV}_{\Phi_H,\mu}(G)-\psi\int qG\,d\mu,
   \qquad
       \Phi_H(\xi)=|\xi|+
          \kappa{\xi^TH\xi\over|\xi|}.                 \tag{0.1}
   \]

The words "one compatible" here refer to \(q\) and \(H\), not to an
unconstrained divergence calibration.  If a single calibration field
\(z\) is introduced, the box constraint necessarily contributes obstacle
multipliers on \(\{G=0\}\) and \(\{G=1\}\).  In particular one must not
write a global equation equating the calibration divergence to \(\psi q\)
without an obstacle term.  Section 5.1 gives
the exact calibrated formula and explains how a negative value of (0.2) is
carried by the upper obstacle.

Coarea then implies that almost every superlevel set is a **global**
minimizer of the same forced set functional

\[
          P_{\Phi_H,\mu}(E)-\psi\int_Eq\,d\mu.          \tag{0.2}
\]

This removes both the equimeasurable nesting obstruction and the need for
secants in the volume parameter.  Every nontrivial such level obeys

\[
 (1-\kappa)\psi m(\mu(E))
 \le P_{\Phi_H,\mu}(E)
 \le \psi m(\mu(E)),                                   \tag{0.3}
\]

is a \(\psi\)-quasiminimizer, and has generalized weighted anisotropic
mean curvature bounded in absolute value by \(\psi\) on every regular
reduced-boundary patch.  The bound remains true on an `SBV` jump patch in
the divided-difference sense.  At a jump crossing a median, equality with
one prescribed sign is generally unavailable; the interval
\([-\psi,\psi]\) is the strongest conclusion.

A common sign is not needed.  Section 8 proves a variable-multiplier killed
Wulff-tube lemma using only \(|\lambda(x)|\le\psi\).  At the frozen scales

\[
 \kappa=10^{-6},\qquad \alpha=10^{-28},\qquad
 \beta=10^{-14},                                        \tag{0.4}
\]

the direct replacement loses less than \(3.664\cdot10^{-5}\) of normalized
projector variance, and the complete smooth short-tube deletion leaves more
than \(.00293\) angular variance.

This is not yet a complete KLS proof.  For smooth density and a smooth hard
support, the remaining geometric issue is a global area/coverage theorem for
the Wulff rays at free-boundary singularities.  For an arbitrary nonsmooth
convex support, approximation must additionally preserve the forced
minimizers, their surface matrices, and their killed-ray loss.  The final
high-rank inverse also remains separate.  None of these gaps affects the
existence, retention, common-\(H\), levelwise minimality, or curvature
statements proved below.

## 1. Precise setting and the median identity

Work on the affine support \(E\simeq\mathbb R^k\) of a non-atomic
log-concave probability, \(k\ge1\).  All derivatives and perimeters below
are relative to this support.  Write

\[
 d\mu=Z^{-1}e^{-V}1_\Omega\,dx,
\]

where \(\Omega\) is the relative interior of the convex support and \(V\)
is convex.  Let \(P_\mu\) and \({\rm TV}_\mu\) denote the relaxed weighted
relative perimeter and variation.  Smooth approximation in compact
subsets of \(\Omega\) gives

\[
                   P_\mu(B)\ge\psi\,m(\mu(B))          \tag{1.1}
\]

for every finite-perimeter set \(B\), where \(\psi=\psi_\mu\) and
\(m(v)=\min(v,1-v)\).  Coarea then gives the corresponding inequality for
every \([0,1]\)-valued `BV` function.

For such a function put

\[
 v_G(r)=\mu(G>r),\qquad
 W(G)=\int_0^1m(v_G(r))\,dr.                          \tag{1.2}
\]

**Lemma 1.1 (atoms and nonunique medians included).**

\[
 \boxed{W(G)=\inf_{c\in\mathbb R}\int|G-c|\,d\mu.}    \tag{1.3}
\]

Moreover, \(W\) is convex and

\[
                    |W(G)-W(U)|\le\|G-U\|_{L^1(\mu)}. \tag{1.4}
\]

**Proof.**  Let \(c\) be any median, so
\(\mu(G<c)\le1/2\) and \(\mu(G>c)\le1/2\).  For \(r<c\),
\(v_G(r)\ge1/2\), while for \(r>c\), \(v_G(r)\le1/2\).  Therefore

\[
\begin{aligned}
 W(G)&=\int_0^c\mu(G\le r)\,dr
       +\int_c^1\mu(G>r)\,dr\\
 &=\int(c-G)_+\,d\mu+\int(G-c)_+\,d\mu.
\end{aligned}
\]

The value at the single level \(r=c\) is irrelevant.  A median minimizes
the last absolute-deviation functional, proving (1.3), including when the
median is nonunique or is an atom.  Finally, \(W\) is the distance in
\(L^1(\mu)\) to the closed subspace of constants.  A distance to a linear
subspace is the quotient norm, hence is convex and one-Lipschitz. \(\square\)

The direct Cheeger deficit is

\[
 D_\psi(G)={\rm TV}_\mu(G)-\psi W(G)\ge0.             \tag{1.5}
\]

The nonnegativity follows by applying (1.1) in the coarea formula.

## 2. The matrix functional and existence

Write the polar decomposition of the weighted derivative as

\[
 DG=\sigma_G|DG|_\mu,qquad |\sigma_G|=1
       \quad |DG|_\mu\hbox{-a.e.},
\]

and set

\[
 M(G)=\int\sigma_G\sigma_G^T\,d|DG|_\mu,
 \qquad {\rm tr}\,M(G)={\rm TV}_\mu(G).              \tag{2.1}
\]

Fix \(0<\kappa<1/3\) and minimize on
\(\mathcal C=\{G\in BV(\mu):0\le G\le1\}\)

\[
 \mathcal K_\kappa(G)=D_\psi(G)
          +\kappa\|M(G)-M(F)\|_*.                    \tag{2.2}
\]

The lower semicontinuity of the matrix term must not be asserted by
itself.  The correct lower-semicontinuous object is its sum with total
variation.  Indeed

\[
 A_\kappa(G):={\rm TV}_\mu(G)
       +\kappa\|M(G)-M(F)\|_*
 =\sup_{\substack{H=H^T\\\|H\|_{\rm op}\le1}}
 \left\{\int\Phi_H(dDG)-\kappa{\rm tr}(HM(F))\right\}, \tag{2.3}
\]

where \(\Phi_H\) is given by (0.1).  If \(n\) is a unit vector and
\(h\perp n\), its spherical Hessian satisfies

\[
 D^2\Phi_H(n)[h,h]
 =|h|^2+\kappa\{2h^THh-(n^THn)|h|^2\}
 \ge(1-3\kappa)|h|^2.                                \tag{2.4}
\]

Thus every \(\Phi_H\) is convex, one-homogeneous, and uniformly elliptic.
Each integral in (2.3) is convex and lower semicontinuous in local
\(L^1(\mu)\), and so is their supremum.  In particular,
\(A_\kappa\) is convex and lower semicontinuous.

**Theorem 2.1 (existence in the full log-concave support).**  If
\(F\in\mathcal C\) has finite variation, (2.2) has a minimizer.

**Proof.**  A minimizing sequence may be chosen with
\(\mathcal K_\kappa(G_j)\le D_\psi(F)+1\).  Since both terms in (2.2) are
nonnegative and \(W(G_j)\le1/2\),

\[
 \operatorname {TV}_\mu(G_j)
 =D_\psi(G_j)+\psi W(G_j)
 \le D_\psi(F)+1+\psi/2.                              \tag{2.5}
\]

On every compact subset of the relative interior, the log-concave density
is continuous and bounded above and below by positive constants.  Ordinary
local `BV` compactness and a diagonal argument give convergence in local
\(L^1(\mu)\).  Tightness of the probability and \(0\le G_j\le1\) upgrade
this to global \(L^1(\mu)\) convergence.  The box constraint is closed.

To make the hard-wall point explicit, choose an exhaustion
\(K_i\Subset\Omega\).  Lower semicontinuity on each \(K_i\) is ordinary
weighted `BV` lower semicontinuity because the density is continuous and
strictly positive there.  Relative variation is the supremum of these
local variations.  No trace term on \(\partial\Omega\) is inserted: this is
the relaxed **relative** perimeter, consistent with exterior Minkowski
content because mass outside the support is zero.  Taking the supremum over
\(i\) proves global lower semicontinuity even when the density drops
discontinuously to zero at a hard wall.  By (2.3), (1.4), and

\[
 \mathcal K_\kappa=A_\kappa-\psi W,
\]

the objective is lower semicontinuous.  The limit is a minimizer.
The argument is relative to the affine support and therefore also covers
degenerate covariance.  A point mass is the excluded \(k=0\) case.
\(\square\)

This proof includes unbounded support and an extended-valued convex
potential.  It proves existence in the relative `BV` space; it does not
claim smoothness at a nonsmooth hard wall.

## 3. Exact deficit and matrix retention

Let \(G\) be any minimizer of (2.2), and abbreviate
\(\Delta=D_\psi(F)\).  Comparison with \(F\), together with (1.5), gives

\[
 \boxed{
 D_\psi(G)\le\Delta,
 \qquad
 \|M(G)-M(F)\|_*\le{\Delta\over\kappa}.}             \tag{3.1}
\]

No equimeasurability, profile minimizer, or regularity assertion is used.
In particular, if \(T={\rm tr}M(F)>0\), \(T'={\rm tr}M(G)\),
\(Q=M(F)/T\), and \(Q'=M(G)/T'\), then with
\(E=\Delta/\kappa<T\),

\[
 T'\ge T-E,
 \qquad
 \|Q'-Q\|_*le {2E\over T}.                           \tag{3.2}
\]

Indeed \(|T'-T|\le\|M(G)-M(F)\|_*\le E\), and

\[
 Q'-Q={M(G)-M(F)\over T}
       +\left({1\over T'}-{1\over T}\right)M(G).
\]

Both normalized matrices are positive semidefinite contractions of trace
one.  Since \({\rm tr}(Q'-Q)=0\),

\[
 |{\rm tr}(Q'^2)-{\rm tr}(Q^2)|
 =|{\rm tr}[(Q'-Q)(Q'+Q-I)]|
 \le\|Q'-Q\|_*.                                      \tag{3.3}
\]

Thus (3.2) is also the exact projector-variance loss bound used below.

## 4. Median subgradients, including atoms

The dual description of the quotient norm is

\[
 W(G)=\sup\left\{\int qG\,d\mu:
       q\in L^\infty(\mu),\ |q|\le1,\ \int q\,d\mu=0\right\}.
                                                               \tag{4.1}
\]

Choose a median \(c\) of \(G\).  A concrete norming functional is

\[
 q=1\quad\hbox{on }\{G>c\},\qquad
 q=-1\quad\hbox{on }\{G<c\},                                \tag{4.2}
\]

and on \(\{G=c\}\) choose a constant in \([-1,1]\) which makes
\(\int q=0\).  Such a constant exists by the two median inequalities.
Then

\[
 \int qG\,d\mu=\int |G-c|\,d\mu=W(G),                         \tag{4.3}
\]

so \(q\in\partial W(G)\).  This explicitly treats a median atom.  If
there is an interval of medians, all choices differ only on a null value
gap or on a median atom, and the construction remains valid.

**Lemma 4.1 (DC linearization at a minimizer).**  For every
\(q\in\partial W(G)\), the same \(G\) minimizes over \(\mathcal C\)

\[
                         A_\kappa(U)-\psi\int qU\,d\mu.       \tag{4.4}
\]

**Proof.**  Fix \(U\in\mathcal C\) and put
\(G_t=(1-t)G+tU\).  Minimality of \(A_\kappa-\psi W\) and the subgradient
inequality give

\[
 A_\kappa(G_t)-A_\kappa(G)
 \ge\psi[W(G_t)-W(G)]
 \ge\psi t\int q(U-G)\,d\mu.                         \tag{4.5}
\]

Convexity of \(A_\kappa\) bounds the left side above by
\(t[A_\kappa(U)-A_\kappa(G)]\).  Divide by \(t\). \(\square\)

This elementary argument is the point at which the sign of the
direct-deficit term matters.  It does not assert that the nonconvex
functional (2.2) itself has a convex Euler equation.

## 5. One compatible constant matrix: the saddle argument

A formal choice of a nuclear-norm subgradient is insufficient: a minimizer
of a supremum need not minimize an arbitrary active summand.  Here a true
saddle exists.

Let \(\mathbb B=\{H=H^T:\|H\|_{\rm op}\le1\}\) and

\[
 L(U,H)=\int\Phi_H(dDU)-\kappa{\rm tr}(HM(F))
             -\psi\int qU\,d\mu.                              \tag{5.1}
\]

For fixed \(H\), this is convex and lower semicontinuous in \(U\); for
fixed \(U\), it is affine and continuous in the finite-dimensional matrix
\(H\).  Since \(\Phi_H\ge(1-\kappa)|\cdot|\), comparison with \(U=0\)
shows that every minimizer in \(U\), uniformly over \(H\), lies in a
fixed total-variation sublevel.  The intersection of that sublevel with
\(\mathcal C\) is compact in \(L^1(\mu)\) by the argument of Theorem 2.1.
Here is the explicit uniform restriction.  If \(U_H\) minimizes
\(L(\cdot,H)\), comparison with zero and
\(\Phi_H\ge(1-\kappa)|\cdot|\) give

\[
 (1-\kappa)\operatorname {TV}_\mu(U_H)
 \le \psi+2\kappa\operatorname {TV}_\mu(F).                    \tag{5.2}
\]

Likewise, comparison of the minimizer \(G\) of the supremum with zero gives
\(\operatorname {TV}_\mu(G)\le\psi+\kappa\operatorname {TV}_\mu(F)\).
Choose \(R\) strictly larger than both bounds and restrict \(U\) to

\[
 \mathcal C_R=\{0\le U\le1:\operatorname {TV}_\mu(U)\le R\}.
\]

This set is convex and compact in \(L^1(\mu)\), contains \(G\), and contains
a minimizer of every \(L(\cdot,H)\).  Thus neither side of the minimax value
changes under the restriction.  The compact convex minimax theorem gives

\[
 \min_{U\in\mathcal C_R}\sup_{H\in\mathbb B}L(U,H)
 =\sup_{H\in\mathbb B}\min_{U\in\mathcal C_R}L(U,H).          \tag{5.3}
\]

The left side is attained at \(G\) by Lemma 4.1.  Let \(H_*\) maximize
the right side.  If the common value is \(a\), then

\[
 a=\min_U L(U,H_*)\le L(G,H_*)
 \le\sup_HL(G,H)=a.                                           \tag{5.4}
\]

Thus equality holds throughout: \(G\) minimizes

\[
 \boxed{{\rm TV}_{\Phi_{H_*},\mu}(U)-\psi\int qU\,d\mu}       \tag{5.5}
\]

over \(0\le U\le1\), and \(H_*\) is active at \(G\).  Equivalently,

\[
 \|H_*\|_{\rm op}\le1,qquad
 {\rm tr}\{H_*[M(G)-M(F)]\}=\|M(G)-M(F)\|_*.                  \tag{5.6}
\]

This proves the existence of one compatible constant anisotropy for all
levels.  No exchange of an infimum and supremum is being assumed without
compactness.

### 5.1 A single calibration and the indispensable obstacle multipliers

Common \(q\) and \(H_*\) do not by themselves say that an unconstrained
vector field calibrates all levels.  In a bounded Lipschitz support with
positive smooth density, the relative-`BV` dual-attainment theorem does give
one field for (5.5).  With the convention

\[
 \operatorname {TV}_{\Phi,\mu}(U)=
 \sup_{\Phi^\circ(z)\le1}\int U\,d(\operatorname {div}_\mu z),
                                                               \tag{5.7}
\]

where the admissible fields have zero relative normal trace on the hard
wall,
there are a bounded divergence-measure field \(z\) and a finite signed
obstacle measure \(\eta\) satisfying

\[
\begin{aligned}
 &\Phi^\circ(z)\le1,
 \qquad -z\cdot\sigma_G=\Phi(\sigma_G)
                   \quad |DG|_\mu\hbox{-a.e.},\\
 &\operatorname {div}_\mu z-\psi q\,\mu+\eta=0,\\
 &\eta\le0\ \hbox{on }\{G=0\},\qquad
   \eta=0\ \hbox{on }\{0<G<1\},\qquad
   \eta\ge0\ \hbox{on }\{G=1\}.
\end{aligned}                                                 \tag{5.8}
\]

The signs in the first line follow from convention (5.7).  Equality in
coarea shows that the same \(z\) calibrates almost every
\(B_r=\{G>r\}\).  Gauss' formula and the second line give, for almost every
\(0<r<1\),

\[
\begin{aligned}
 P_{\Phi,\mu}(B_r)
 &=(\operatorname {div}_\mu z)(B_r)\\
 &=\psi\int_{B_r}q\,d\mu-\eta(B_r).
\end{aligned}                                                 \tag{5.9}
\]

Because \(\eta=0\) where \(0<G<1\), every active threshold contains the
same positive-obstacle charge and excludes the negative-obstacle charge.
Consequently

\[
 \mathcal E_q(B_r)=-\eta(B_r)
 =-\eta(\{G=1\})=e_*\le0.                                   \tag{5.10}
\]

This resolves the apparent contradiction when \(e_*<0\).  Dropping
\(\eta\) would force \(\mathcal E_q(B_r)=0\), and is wrong.  The upper
obstacle carries the common negative energy; after integrating (5.8), the
lower obstacle balances it.  On an unbounded support or a nonsmooth hard
wall, exhaustion gives local divergence-measure calibrations, but a global
normal trace and dual attainment require an additional compactness theorem.
None of Sections 6--8 uses a common calibration: their curvature bound
follows directly from set quasiminimality.

There is a one-dimensional exact stress test.  Let \(\mu\) be uniform on
\([-1,1]\), so \(\psi=1\), take
\(\Phi(\xi)=(1-\kappa)|\xi|\), and let
\(B=(0,1)\), \(q=1_B-1_{B^c}\).  Then

\[
 P_{\Phi,\mu}(B)={1-\kappa\over2},
 \qquad \int_Bq\,d\mu={1\over2},
 \qquad \mathcal E_q(B)=-{\kappa\over2}.                     \tag{5.11}
\]

For every set \(C\), the Cheeger inequality and (6.4) below give
\(\mathcal E_q(C)\ge-\kappa m(\mu(C))\ge-\kappa/2\), so \(B\)
is a global minimizer and \(e_*=-\kappa/2<0\).  This explicitly disproves
the obstacle-free divergence equation even in a smooth one-dimensional
model.

## 6. Exact levelwise minimality and near-Cheeger bounds

For a finite-perimeter set \(B\subset\Omega\), put

\[
 \mathcal E_q(B)=P_{\Phi_{H_*},\mu}(B)-\psi\int_Bq\,d\mu,
 \qquad e_*:=\inf_B\mathcal E_q(B).                           \tag{6.1}
\]

The empty and full sets both have value zero because \(\int q=0\), hence
\(e_*\le0\).  Coarea and layer cake give, for every \(0\le U\le1\),

\[
 {\rm TV}_{\Phi_{H_*},\mu}(U)-\psi\int qU\,d\mu
 =\int_0^1\mathcal E_q(\{U>r\})\,dr\ge e_*.                  \tag{6.2}
\]

Conversely, indicators show that the infimum of the left side is exactly
\(e_*\).  Applying (6.2) to the minimizer \(G\) proves

\[
 \boxed{\mathcal E_q(\{G>r\})=e_*
        \quad\hbox{for a.e. }r\in(0,1).}                      \tag{6.3}
\]

Thus almost every level is a global forced-perimeter minimizer; no
independent level replacement and no foliation concavity is required.

For \(|q|\le1\), \(\int q=0\), and any set \(B\),

\[
 \int_Bq\,d\mu\le\mu(B),
 \qquad
 \int_Bq\,d\mu=-\int_{B^c}q\,d\mu\le1-\mu(B).                \tag{6.4}
\]

If \(B=\{G>r\}\) is a nontrivial level satisfying (6.3), then

\[
\begin{aligned}
 P_{\Phi_{H_*},\mu}(B)
 &=e_*+\psi\int_Bq\,d\mu
 \le\psi m(\mu(B)),\\
 P_{\Phi_{H_*},\mu}(B)
 &\ge(1-\kappa)P_\mu(B)
 \ge(1-\kappa)\psi m(\mu(B)).
\end{aligned}                                                \tag{6.5}
\]

This is (0.3).  Notice that (6.5) is stronger than selecting good levels
from an integrated deficit.

Comparison in (6.3) with an arbitrary finite-perimeter set \(C\) gives

\[
 \boxed{P_{\Phi_{H_*},\mu}(B)
 \le P_{\Phi_{H_*},\mu}(C)+\psi\,\mu(B\mathbin\triangle C).} \tag{6.6}
\]

Thus every relevant level is a global \(\psi\)-quasiminimizer.

## 7. Euler theorem, jumps, Cantor variation, and hard walls

Assume first that \(V\in C^2\) in \(\Omega\) and that the surface tension
is the uniformly elliptic \(\Phi_{H_*}\).  Let \(B\) be a level satisfying
(6.3).  If \(X\) is a compactly supported \(C^1\) vector field in
\(\Omega\) and \(\varphi_t\) is its flow, (6.6) applied to
\(C=\varphi_t(B)\) yields

\[
 |\delta P_{\Phi_{H_*},\mu}(B)[X]|
 \le\psi\int_{\partial^*B}|X\cdot n_B|\,d\sigma_\mu.         \tag{7.1}
\]

Here one uses the exact transport limit

\[
 \lim_{t\to0}{\mu(B\mathbin\triangle\varphi_tB)\over|t|}
 =\int_{\partial^*B}|X\cdot n_B|\,d\sigma_\mu.                \tag{7.2}
\]

The anisotropic first variation is reparametrization invariant and hence
depends only on \(X\cdot n_B\).  The \(L^1\)-duality in (7.1) produces a
measurable generalized weighted anisotropic mean curvature
\(\lambda_B\) such that

\[
 \delta P_{\Phi_{H_*},\mu}(B)[X]
 =\int_{\partial^*B}\lambda_B(X\cdot n_B)\,d\sigma_\mu,
 \qquad
 \boxed{|\lambda_B|\le\psi\quad\sigma_\mu\hbox{-a.e.}}       \tag{7.3}
\]

Equivalently, for a Wulff normal variation \(X=fD\Phi_{H_*}(n_B)\),

\[
 \delta P_{\Phi_{H_*},\mu}(B)[X]
 =\int\lambda_B f\,d\sigma_{\Phi_{H_*},\mu}.                 \tag{7.4}
\]

Thus the same bound \(|\lambda_B|\le\psi\) holds in exactly the gauge
used by the Wulff tube; no factor \((1-\kappa)^{-1}\) is lost.

The precise regularity input used here is the interior regularity theorem
for parametric elliptic \(\Lambda\)-minimizers: if the density is positive
and \(C^2\) on \(U\Subset\Omega\), the one-homogeneous integrand is \(C^2\)
off zero and has tangent ellipticity bounded below, and (6.6) holds for all
competitors compactly supported in \(U\), then
\(\partial B\cap U\) is the disjoint union of a relatively open
\(C^{1,\eta}\) reduced boundary and a closed singular set of
\(\mathcal H^{k-1}\)-measure zero.  This is Tamanini's
elliptic-quasiminimizer theorem in its weighted parametric form; the weight
is absorbed into the \(C^2\) position-dependent integrand.  Its constants
depend locally only on ellipticity, the \(C^2\) density bounds, and
\(\Lambda=\psi\), not on the particular level.

On a regular graph, (7.3) is a uniformly elliptic prescribed-mean-curvature
equation with \(L^\infty\) right side.  Interior difference-quotient and
Calderon--Zygmund regularity give local \(W^{2,s}\) for every finite \(s\).
In particular the second fundamental form and (7.3) are defined almost
everywhere.  These two explicitly stated regularity theorems are the only
`BV`-to-surface inputs in the Jacobian calculation.  They do not assert
absence of high-codimension singularities or regularity at a nonsmooth
wall.

### 7.1 `SBV` jumps and median crossings

Suppose \(x\in J_G\) has approximate traces \(a<b\).  For almost every
\(r\in(a,b)\), the same oriented jump patch is contained in
\(\partial^*\{G>r\}\).  Equations (6.6)--(7.3) apply to every such regular
threshold patch.  If one writes a jump multiplier by averaging in the
level variable,

\[
 \bar\lambda(x)={1\over b-a}\int_a^b\lambda_r(x)\,dr,
\]

then

\[
                         |\bar\lambda(x)|\le\psi.              \tag{7.5}
\]

If the jump does not cross the chosen median and the forcing has the same
trace on both sides, the classical equation has the corresponding constant
sign.  If it crosses the median, \(q\) can jump from \(-1\) to \(+1\).
Outward and inward variations then give only the two inequalities
\(-\psi\le\lambda\le\psi\).  Hence the pointwise equality proposed in
(10.4) of the source note is not a valid general `SBV` statement.  The
bounded interval (7.5) is both valid and sufficient for Section 8.

### 7.2 Cantor part

Nothing in the argument proves \(D^cG=0\), and such a claim is unnecessary.
The matrix coarea identity

\[
 M(G)=\int_0^1M(\{G>r\})\,dr                              \tag{7.6}
\]

includes absolutely continuous, jump, and Cantor variation.  By (6.3),
almost every surface in this identity is a forced-perimeter minimizer; by
regularity, its singular stratum has zero surface measure.  Fubini therefore
places all of the retained matrix mass on regular threshold patches even
when the function itself has Cantor variation.  There is no need to assign
a pointwise Euler multiplier directly to \(D^cG\).

### 7.3 Hard support

For a \(C^2\) convex support, the same comparison gives (7.1) for flows
tangent to the wall and the weak anisotropic Young condition on regular
free-boundary points.  Interior regular patches retain the exact bound
(7.3).  At a nonsmooth convex wall, (6.6) remains valid in relative `BV`,
but neither a single classical Young normal nor boundary \(W^{2,s}\)
regularity follows.  Any use of a classical contact formula there requires
a separate smooth-support approximation theorem.

## 8. A variable-multiplier short killed-Wulff tube

This section proves that a common constant multiplier is not required.
Let \(B\) be a smooth regular level from Section 6, put

\[
 v_0=\mu(B),\qquad m_0=m(v_0),\qquad P_0=P_{\Phi,\mu}(B),
 \qquad c=(1-\kappa)\psi,                                    \tag{8.1}
\]

and assume the regular killed-ray coverage theorem.  Along the outward
Wulff ray \(x+t z_x\), \(z_x=D\Phi(n_x)\), let \(\tau(x)\) be the first
cut, focal, or support-contact time.  The weighted Jacobian has the exact
form

\[
 j_x(t)=\exp\{\lambda(x)t-D_x(t)\},
 \qquad D_x(0)=0,quad D_x\ge0,quad D_x'\ge0,                 \tag{8.2}
\]

because the second derivative of \(\log j_x\) is nonpositive.  Equation
(7.3) gives \(|\lambda(x)|\le\psi\).

Define

\[
 R(t)=\int_{\partial B}1_{\{t<\tau(x)\}}e^{-D_x(t)}
                         \,d\sigma_{\Phi,\mu}(x).             \tag{8.3}
\]

Then \(R\) is nonincreasing, \(R(0)=P_0\), and, if
\(v(t)=\mu(B_t)\),

\[
 e^{-\psi t}R(t)\le v'(t)\le e^{\psi t}R(t)                  \tag{8.4}
\]

for almost every regular time.  Fix \(0<\gamma<1\), and let \(T\) be the
first time at which \(v(T)=v_0+\gamma m_0\).  From (6.5),

\[
                  cm_0\le P_0\le\psi m_0.                    \tag{8.5}
\]

The upper inequality in (8.4) and \(R(t)\le P_0\) give

\[
 \gamma m_0\le P_0{e^{\psi T}-1\over\psi}
 \le m_0(e^{\psi T}-1).
\]

The global anisotropic Cheeger lower bound, together with
\(m(v(t))\ge(1-\gamma)m_0\), gives the reverse time bound.  Hence

\[
 \boxed{{\log(1+\gamma)\over\psi}\le T
 \le{\gamma\over(1-\kappa)(1-\gamma)\psi}.}                  \tag{8.6}
\]

At regular times, (8.4) and the same Cheeger lower bound imply

\[
 R(t)\ge e^{-\psi t}(1-\kappa)\psi(1-\gamma)m_0.
\]

Letting \(t\uparrow T\), dividing by \(P_0\le\psi m_0\), and using the
upper bound in (8.6) yields

\[
 \boxed{{R(T-)\over P_0}\ge a_{\kappa,\gamma}:=
 (1-\kappa)(1-\gamma)
 \exp\left\{-{\gamma\over(1-\kappa)(1-\gamma)}\right\}.}     \tag{8.7}
\]

Put \(\xi=1-a_{\kappa,\gamma}\).  For \(h>0\), delete rays killed before
\(T\) and surviving rays with \(D_x(T-)>h\).  The exact defect identity
shows that their total base surface is at most

\[
 \boxed{{\xi\over1-e^{-h}}P_0.}                               \tag{8.8}
\]

This slightly sharpens the bound obtained by estimating the killed and
high-\(D\) sets separately.  On every remaining ray and every \(0<t<T\),

\[
                         j_x(t)\ge e^{-\psi T-h}.              \tag{8.9}
\]

Consequently, for the good base set \(\Gamma_h\),

\[
 \boxed{
 {\rm Cov}(\mu)\succeq
 {A_{\rm var}(\kappa,\gamma,h)\over\psi^3}
 \int_{\Gamma_h}z_x\otimes z_x\,d\sigma_{\Phi,\mu}(x),}    \tag{8.10}
\]

where

\[
 A_{\rm var}(\kappa,\gamma,h)
 ={[\log(1+\gamma)]^3\over12}
 \exp\left\{-h-{\gamma\over(1-\kappa)(1-\gamma)}\right\}.  \tag{8.11}
\]

To prove (8.10), restrict the variance integral to the injective good-ray
image and use

\[
 \int_0^T(\langle\theta,x+tz_x\rangle-a)^2dt
 \ge{T^3\over12}\langle\theta,z_x\rangle^2
\]

for every constant \(a\), followed by (8.6) and (8.9).

The proof of Section 8 uses only a two-sided curvature bound.  In
particular, median-crossing jumps do not need to be oriented to have one
common sign.

## 9. Singular origins and coverage: exact remaining scope

The Euler theorem does not by itself imply the killed-ray coverage used in
(8.4).  There are two separate issues.

First, at an interior singular point which is the closest point of an
exterior Wulff ball, blow-up of the \(\psi\)-quasiminimizer gives an
anisotropic perimeter-minimizing cone contained in a halfspace.  The strong
maximum principle for elliptic minimizing cones forces that cone to be a
halfspace; excess regularity then makes the original point regular.  Thus,
conditional only on these standard two regularity theorems, an interior
singular point cannot be the origin of a positive-length closest Wulff ray.
Together with almost-everywhere uniqueness of the anisotropic distance
projection, this gives interior coverage by regular rays.

Second, the corresponding statement at a free-boundary singularity in a
hard support requires classification of relative minimizing cones with the
Young condition.  It is not proved by the interior argument.  At a regular
free-boundary edge in a \(C^2\) convex support, the Young condition gives
the single-ray formula and the edge has zero tube volume.  At a nonsmooth
wall the normal cone can have several generators.  Therefore the strongest
fully justified scope of (8.10) from the calculations in this report is:

* all regular rays in a smooth density, stopped at first contact, cut, or
  focal time;
* global coverage in the interior once the minimizing-cone maximum
  principle and excess theorem are invoked; and
* global hard-wall coverage only under the explicit regular-coverage
  hypothesis.

This limitation is geometric, not an Euler or matrix-fidelity gap.

## 10. Explicit constants at \(\kappa=10^{-6}\),
\(\alpha=10^{-28}\)

Use the audited heat constants

\[
 {\Delta\over\psi}\le6.0132566\cdot10^{-14},\qquad
 {T\over\psi}>.0032827359,qquad
 1-{\rm tr}Q^2>.00326880.                                    \tag{10.1}
\]

At \(\kappa=10^{-6}\), (3.1)--(3.3) give

\[
\begin{aligned}
 {\|M(G)-M(F)\|_*\over\psi}
 &\le6.0132566\cdot10^{-8},\\
 {T'\over\psi}&>.0032826757,\\
 |{\rm tr}(Q'^2)-{\rm tr}(Q^2)|
 &<3.663565\cdot10^{-5}.
\end{aligned}                                                \tag{10.2}
\]

Thus the Euclidean projector variance after direct replacement is greater
than \(.00323216\).  The audited pointwise conversion to the normalized
Wulff-displacement law costs at most \(42\kappa=4.2\cdot10^{-5}\), leaving
more than \(.0031901\).

Take \(\gamma=2\cdot10^{-5}\) and \(h=1\).  Equations (8.7)--(8.8) give

\[
\begin{aligned}
 {\gamma\over(1-\kappa)(1-\gamma)}
 &=2.0000421\cdot10^{-5},\\
 \xi&<4.0999781\cdot10^{-5},\\
 {\xi\over1-e^{-1}}&<6.486070\cdot10^{-5}.
\end{aligned}                                                \tag{10.3}
\]

After changing from anisotropic base surface to the
\(|z|^2d\sigma_\Phi\) trace law and applying the exact deletion estimate
\(|\Delta(1-{\rm tr}Q^2)|\le4\varepsilon\), the surviving long-ray packet
has angular variance greater than

\[
                              \boxed{.00293}.                  \tag{10.4}
\]

Finally,

\[
 \boxed{A_{\rm var}(10^{-6},2\cdot10^{-5},1)
       >2.4524\cdot10^{-16}.}                                 \tag{10.5}
\]

All these numbers are universal.  The enormous heat-rank threshold
\(k\ge1.81\cdot10^{47}\) remains universal as well; below it, the existing
finite-dimensional Buser--Ledoux branch supplies a (very small) universal
bound.  The arithmetic therefore has no hidden dependence on dimension,
smoothing order, or bootstrap depth.

## 11. What is proved and what remains

The direct functional rigorously supplies:

1. a minimizer for every non-point-mass log-concave measure on its affine
   support;
2. exact direct-deficit and normal-matrix retention;
3. one median subgradient \(q\), including atoms and nonunique medians;
4. one compatible constant nuclear matrix \(H_*\), obtained by a saddle
   theorem rather than a formal subgradient choice;
5. global forced minimality of almost every threshold;
6. the near-Cheeger estimate (0.3), \(\psi\)-quasiminimality, and the sharp
   generalized curvature bound \(|\lambda|\le\psi\);
7. the same bound on `SBV` jump divided differences, with the Cantor part
   handled through coarea; and
8. a dimension-free variable-multiplier short-tube theorem in the regular
   coverage setting.

It does not by itself prove hard-wall singular coverage, stability of that
coverage under arbitrary convex approximation, or the final inverse which
excludes a high-rank surviving packet.  Those are the remaining
load-bearing statements; the former equimeasurable nesting/concavity gap is
no longer among them.
