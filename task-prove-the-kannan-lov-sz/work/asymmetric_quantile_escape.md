# Asymmetric quantiles do not remove the single-jump branch

## 0. Verdict

There is a quantitatively consistent asymmetric perturbation of the median
quotient, but it does **not** give the missing conclusion.

For a cusp at \(\tau=1/2-d<1/2\), the largest asymmetric tent which is
dominated by the Cheeger tent is

\[
 B_\tau(v)=\min\left\{v,{\tau\over1-\tau}(1-v)\right\}.
 \tag{0.1}
\]

Its relative loss at volume \(1/2\) is

\[
 1-{B_\tau(1/2)\over1/2}
 ={4d\over1+2d}.                                      \tag{0.2}
\]

Consequently the heat deficit can indeed be retuned so that the normal
matrix is retained: if the matrix penalty is \(\kappa\) and the normalized
projector variance to be retained is \(\nu\), it is enough, at the level of
the comparison estimate, to take \(d=O(\kappa\nu)\).  One may then tune the
heat error so that the comparator's active volumes lie on one side of the
cusp.

The fatal point is that the minimizer is free to move its *surface matrix*
to the cusp.  Total nuclear-matrix fidelity does not remember the level
label.  This is not merely a logical possibility.  Section 4 gives an
explicit isotropic log-concave example in every dimension.  For the product
of one-sided exponentials there are sets \(A_v\) with

\[
 M(A_v)={P_n(v)\over n}I_n.
\]

A binary central comparator and a binary cusp competitor can be scaled to
have exactly the same full normal matrix, whose normalized projector
variance is \(1-1/n\).  Nevertheless the cusp competitor has strictly
smaller asymmetric direct deficit.  The conclusion is unchanged by the
nuclear penalty, by an exact mean constraint, or by any regularizer which
depends only on total variation or the aggregate normal matrix.

There is also a general concavity obstruction.  For every concave coarea
profile \(B\le m\), an off-centre breakpoint approximates the Cheeger tent
at least as well, in relative terms, as the centre does.  Thus lowering the
tent cannot make a cusp expensive without paying at least the same relative
cost on the central heat matrix.

Exact equimeasurability does remove a positive-height cusp plateau when the
prescribed value law has none, but restores only volume-constrained
stability.  Smooth regularization of the tent replaces the cusp by a
positive rank-one term in the second variation; it likewise does not give
unconstrained stability.  These failures are quantified below.  Therefore
an asymmetric quantile, a finite collection of scalar moment constraints,
or an infinitesimal lexicographic tilt does not eliminate the
median-crossing branch while retaining one spatially constant matrix
anisotropy.

## 1. The exact asymmetric quotient

Let \(a,b>0\), and for \(G\in L^1(\mu)\) define

\[
 W_{a,b}(G)=\inf_{c\in\mathbb R}
 \left\{a\int(G-c)_+\,d\mu+b\int(c-G)_+\,d\mu\right\}.
 \tag{1.1}
\]

Put

\[
 \tau={b\over a+b},\qquad
 B_{a,b}(v)=\min\{av,b(1-v)\}.                       \tag{1.2}
\]

The minimizers in (1.1) are the \(\tau\)-quantiles: when the value law has
no atom at \(c\), the Euler equation is

\[
                     \mu(G>c)=\tau.                 \tag{1.3}
\]

Atoms are handled by the usual interval of one-sided derivatives.  The
dual formula, valid without any atomlessness assumption, is

\[
 W_{a,b}(G)=\sup\left\{\int qG\,d\mu:
       -b\le q\le a,\quad\int q\,d\mu=0\right\}.      \tag{1.4}
\]

Indeed, subtracting a constant from \(G\) does not affect either side.  For
a minimizing quantile \(c\), take \(q=a\) on \(\{G>c\}\), \(q=-b\) on
\(\{G<c\}\), and choose \(q\in[-b,a]\) on \(\{G=c\}\) so that its mean
is zero.  This attains (1.4); the reverse inequality follows pointwise from

\[
 q(G-c)\le a(G-c)_++b(c-G)_+.
\]

For every \([0,1]\)-valued \(G\), with \(v_G(r)=\mu(G>r)\), layer cake
gives the exact identity

\[
 \boxed{W_{a,b}(G)=\int_0^1 B_{a,b}(v_G(r))\,dr.}     \tag{1.5}
\]

One proof is to split the layer-cake integral at a minimizing quantile.
Equivalently, both sides of (1.5) are translation invariant, positively
homogeneous, and agree on every indicator:

\[
 W_{a,b}(\mathbf1_A)=\min\{a\mu(A),b(1-\mu(A))\}.
 \tag{1.6}
\]

The convenient normalized parametrization is

\[
 B_{\tau,s}(v)=s\,m(\tau)
 \min\left\{{v\over\tau},{1-v\over1-\tau}\right\},
 \qquad 0<s\le1.                                    \tag{1.7}
\]

Thus the cusp is at \(\tau\), its height is \(s m(\tau)\), and

\[
                         0\le B_{\tau,s}(v)\le m(v)  \tag{1.8}
\]

for every \(v\in[0,1]\).  Conversely, among all single tents with cusp
\(\tau\) which are bounded above by \(m\), (1.7) with \(s=1\) is pointwise
maximal.  This follows because the cusp height cannot exceed
\(m(\tau)\), and a tent is linear from the cusp to each endpoint.

Let \(\psi=\psi_\mu\).  The asymmetric direct deficit is

\[
 D_{\tau,s}(G)={\rm TV}_\mu(G)
       -\psi\int_0^1B_{\tau,s}(v_G(r))\,dr\ge0.       \tag{1.9}
\]

The inequality is exactly coarea plus
\(P_\mu(E)\ge\psi m(\mu(E))\) and (1.8).  Hence this modification does not
assume any strengthening of Cheeger.

## 2. The exact asymmetry--matrix tradeoff

Write

\[
 D_m(F)={\rm TV}_\mu(F)-\psi\int_0^1m(v_F(r))\,dr,
 \qquad T={\rm TV}_\mu(F)={\rm tr}\,M(F).            \tag{2.1}
\]

Then

\[
 D_{\tau,s}(F)=D_m(F)+\psi\int_0^1
       [m(v_F(r))-B_{\tau,s}(v_F(r))],dr.            \tag{2.2}
\]

Assume first that \(\tau=1/2-d<1/2\).  At the cusp and at the centre the
relative losses are respectively

\[
 \delta_{\rm cusp}=1-s,
 \qquad
 \delta_0=1-s{\tau\over1-\tau}
          =1-s{1-2d\over1+2d}.                       \tag{2.3}
\]

In particular

\[
 \boxed{\delta_{\rm cusp}\le\delta_0,}               \tag{2.4}
\]

with strict inequality whenever \(s>0\) and \(d>0\).  For the maximal
normalization \(s=1\),

\[
                 \delta_0={4d\over1+2d}.             \tag{2.5}
\]

More precisely, if \(0<h<d\) and
\(v\in[1/2-h,1/2+h]\), all of these volumes lie on the right arm of the
asymmetric tent.  For \(v\ge1/2\), the relative loss equals \(\delta_0\).
For \(v\in[1/2-h,1/2]\), it lies between \(\delta_0\) and

\[
 \delta_-=
 1-s{(1/2-d)(1/2+h)\over(1/2+d)(1/2-h)}.             \tag{2.6}
\]

When \(s=1\), this endpoint value simplifies to

\[
 \delta_-={d-h\over(1/2+d)(1/2-h)}.                  \tag{2.7}
\]

The formulas for \(\tau>1/2\) follow by replacing \(v\) with \(1-v\).

Now minimize

\[
 \mathcal J_{\tau,s,\kappa}(G)=D_{\tau,s}(G)
       +\kappa\|M(G)-M(F)\|_*,\qquad 0<\kappa<1/3,   \tag{2.8}
\]

over \(0\le G\le1\).  The same convexification and saddle argument as for
the median quotient applies: (1.4) supplies the scalar subgradient, while
nuclear/operator duality supplies one spatially constant matrix
\(H\), \(\|H\|_{\rm op}\le1\).  Comparison with \(F\) gives, for every
minimizer,

\[
 D_{\tau,s}(G)\le D_{\tau,s}(F),\qquad
 \|M(G)-M(F)\|_*\le{D_{\tau,s}(F)\over\kappa}.       \tag{2.9}
\]

Suppose first that all comparator levels lie in the central band above.
Since \(\psi\int m(v_F)\le T\), (2.2) gives

\[
 {D_{\tau,s}(F)\over T}
 \le {D_m(F)\over T}+\delta_0.                       \tag{2.10}
\]

More generally, if the levels outside that band carry total surface trace
\(T_{\rm out}\), the completely quantitative replacement is

\[
 {D_{\tau,s}(F)\over T}
 \le {D_m(F)\over T}+\delta_0+{T_{\rm out}\over T}.  \tag{2.10a}
\]

Indeed the relative asymmetric loss is at most \(\delta_0\) on the central
band and at most one elsewhere, while \(\psi m(v)\le
P_\mu(\{F>r\})\).  The clipped heat construction can make
\(T_{\rm out}/T\) arbitrarily small by restricting its active value band
and retuning the heat error.  It does not affect the relocation issue
below.

Let

\[
 Q={M(F)\over T},\qquad
 \nu=1-\operatorname{tr}Q^2.                        \tag{2.11}
\]

If \(E=D_{\tau,s}(F)/\kappa<T\), then the normalized matrix \(Q_G\) of
the minimizer satisfies

\[
 \|Q_G-Q\|_*\le{2E\over T},\qquad
 1-\operatorname{tr}Q_G^2
 \ge\nu-{2E\over T}.                                \tag{2.12}
\]

Thus retaining at least \(\nu/2\) by this comparison argument requires

\[
 {D_m(F)\over T}+\delta_0\le{\kappa\nu\over4}.       \tag{2.13}
\]

For \(s=1\), if

\[
 \varepsilon={\kappa\nu\over4}-{D_m(F)\over T}>0,
\]

then (2.5) and (2.13) force

\[
 \boxed{d\le{\varepsilon\over4-2\varepsilon}.}      \tag{2.14}
\]

For orientation, with \(\kappa=10^{-6}\), \(\nu=3\cdot10^{-3}\), and
\(D_m(F)/T\le2\cdot10^{-11}\), the right side of (2.14) is less than
\(1.83\cdot10^{-10}\).  The heat and clipping parameters can in principle
be made small enough that the comparator's surface matrix lies in a band
\(|v-1/2|<h<d\).  This calculation is not the obstruction.

The obstruction is that (2.9) controls only the aggregate matrix.  It does
not imply that the minimizer's matrix remains on levels in that central
band.

## 3. Why lowering or smoothing the tent cannot charge the cusp

The inequality (2.4) is an instance of a general fact.

**Lemma 3.1 (relative-contact monotonicity).**  Let
\(B:[0,1]\to[0,\infty)\) be concave with \(B(0)=B(1)=0\).  Then

\[
 {B(u)\over u}\ge {B(v)\over v}
       \quad(0<u<v<1),                               \tag{3.1}
\]

and

\[
 {B(u)\over1-u}\le {B(v)\over1-v}
       \quad(0<u<v<1).                               \tag{3.2}
\]

**Proof.**  Concavity above the chord joining \((0,0)\) to
\((v,B(v))\) gives (3.1).  Concavity above the chord joining
\((u,B(u))\) to \((1,0)\) gives (3.2). \(\square\)

Assume also \(B\le m\).  If \(\tau<1/2\), (3.1) gives

\[
 1-{B(\tau)\over m(\tau)}
 \le1-{B(1/2)\over m(1/2)}.                         \tag{3.3}
\]

If \(\tau>1/2\), (3.2) gives the same statement.  Hence every off-centre
breakpoint of a concave coarea profile has relative Cheeger gap no larger
than the central gap.  In particular, making an off-centre cusp strictly
costly by replacing \(B\) with \(sB\) costs at least as much, relatively,
on a central comparator.

This covers the usual broad convex class, not just the single tent.  If
\(B\) is Lipschitz and concave with zero endpoint values, then

\[
 \mathcal W_B(G)=\int_{\mathbb R}B(\mu(G>r))\,dr     \tag{3.4}
\]

is a translation-invariant convex homogeneous functional.  Indeed, if
\(G^*\) is the decreasing rearrangement of \(G\), integration by parts
gives

\[
 \mathcal W_B(G)=\int_0^1G^*(u)B'_+(u)\,du.          \tag{3.5}
\]

The right side is the supremum of \(\int Gq\,d\mu\) over rearrangements
of the decreasing mean-zero function \(B'_+\), by the
Hardy--Littlewood inequality; taking the closed convex hull of those
rearrangements gives an explicit support-function representation.  Thus
(3.4) is convex.  Lemma 3.1 applies to every piecewise-affine profile in
this class.

There is also a direct deficit interpretation.  At a level of volume
\(v\),

\[
 P_\mu(E)-\psi B(v)
 \ge \left(1-{B(v)\over m(v)}\right)P_\mu(E),        \tag{3.6}
\]

because \(P_\mu(E)\ge\psi m(v)\).  Therefore a scalar deficit can bound
the cusp surface by \(D/\delta_{\rm cusp}\) only if the cusp has relative
gap \(\delta_{\rm cusp}>0\).  By (3.3), a central comparator pays at least
that relative gap.  Even in the ideal case in which all comparator levels
saturate Cheeger, the resulting estimate is only
\({\rm tr}\,M_{\rm cusp}\le T\); it deletes no positive universal fraction.

## 4. A high-rank, exact-matrix relocation counterexample

The preceding loss comparison is realized by actual isotropic log-concave
measures.

Let \(X_1,\ldots,X_n\) be independent exponential random variables of
mean one.  The law of

\[
                         X-(1,\ldots,1)              \tag{4.1}
\]

is isotropic and log-concave on \([-1,\infty)^n\).  Translation does not
change any perimeter or normal matrix, so calculations may be made in the
positive orthant.

For \(0<v<1\), define \(t_v>0\) by

\[
                    \mu(A_v)=v,qquad
 A_v=\{x:\max_i x_i\ge t_v\}.                       \tag{4.2}
\]

Put \(q_v=e^{-t_v}=1-(1-v)^{1/n}\).  The reduced boundary of \(A_v\)
consists, up to a codimension-two null set, of the \(n\) facets

\[
 x_i=t_v,qquad x_j<t_v\quad(j\ne i).
\]

Each facet has weighted area \(q_v(1-q_v)^{n-1}\), and its normal
projector is \(e_ie_i^T\).  Consequently

\[
 \boxed{
 P_n(v)=n[1-(1-v)^{1/n}](1-v)^{(n-1)/n},\qquad
 M(A_v)={P_n(v)\over n}I_n.}                         \tag{4.3}
\]

In particular the normalized normal matrix is \(I_n/n\), and its
projector variance is

\[
                         1-{1\over n}.               \tag{4.4}
\]

Fix \(0<\tau<1/2\), \(0<s\le1\), and use the tent (1.7).  Let
\(P_0=P_n(1/2)\) and \(P_\tau=P_n(\tau)\).  Choose \(a_0>0\) sufficiently
small and put

\[
 F=a_0\mathbf1_{A_{1/2}},\qquad
 G=a_1\mathbf1_{A_\tau},\qquad
 a_1=a_0{P_0\over P_\tau}.                          \tag{4.5}
\]

Taking \(a_0\le\min\{1,P_\tau/P_0\}\) ensures that both functions take
values in \([0,1]\).  By (4.3),

\[
 M(F)=M(G)={a_0P_0\over n}I_n.                      \tag{4.6}
\]

Thus every aggregate nuclear-fidelity penalty is exactly zero for both
functions, for every value of \(\kappa\).

The relevant strict inequality is elementary.  Formula (4.3) gives

\[
 {P_n(v)\over1-v}
 =n\left((1-v)^{-1/n}-1\right),                     \tag{4.7}
\]

which is strictly increasing in \(v\).  Therefore

\[
                  P_\tau<2(1-\tau)P_0.              \tag{4.8}
\]

On the other hand,

\[
 B_{\tau,s}(\tau)=s\tau,qquad
 B_{\tau,s}(1/2)={s\tau\over2(1-\tau)}.             \tag{4.9}
\]

Combining (4.8)--(4.9) yields

\[
 {B_{\tau,s}(\tau)\over P_\tau}
 >{B_{\tau,s}(1/2)\over P_0}.                       \tag{4.10}
\]

Let \(\psi\) be the actual Euclidean Cheeger constant of this product
measure; its numerical value is not needed.  Since both functions in
(4.5) have total variation \(T=a_0P_0\),

\[
\begin{aligned}
 D_{\tau,s}(G)-D_{\tau,s}(F)
 &=-\psi T\left[
 {B_{\tau,s}(\tau)\over P_\tau}
 -{B_{\tau,s}(1/2)\over P_0}\right]\\
 &<0.                                                     \tag{4.11}
\end{aligned}
\]

The strict improvement can be written without hidden constants as

\[
 \psi T,{s\tau\over P_\tau}
 \left[1-
 {n((1-\tau)^{-1/n}-1)\over n(2^{1/n}-1)}\right]>0.   \tag{4.12}
\]

Thus the asymmetric functional actively prefers moving the entire
high-rank normal matrix from the central volume to its own cusp.

### 4.1 An exact mean constraint does not repair the example

The same example can be made to preserve the mean of the value function.
First note that \(P_n(v)/v\) is strictly decreasing.  To see this, write
\(r=(1-v)^{1/n}\).  Then

\[
 {P_n(v)\over v}
 ={nr^{n-1}\over1+r+\cdots+r^{n-1}}.                 \tag{4.13}
\]

The right side is strictly increasing in \(r\), because the logarithmic
derivative of the numerator is \((n-1)/r\), while that of the denominator
is a weighted average of \(0/r,1/r,\ldots,(n-1)/r\) which is strictly
smaller.  Since \(r\) decreases with \(v\), the assertion follows.
Consequently

\[
                       {P_\tau\over\tau}>2P_0.       \tag{4.14}
\]

Keep \(F\) as in (4.5), and replace \(G\) by

\[
 \widetilde G=c+a_1\mathbf1_{A_\tau},qquad
 c={a_0\over2}-a_1\tau.                             \tag{4.15}
\]

Equation (4.14) gives \(c>0\), and

\[
                   \int\widetilde G\,d\mu
                   ={a_0\over2}=\int F\,d\mu.       \tag{4.16}
\]

Translation invariance of \(W_{a,b}\), total variation, and the normal
matrix leaves (4.6) and (4.11) unchanged.  Finally,

\[
 \|\widetilde G\|_\infty
 =a_0\left[{1\over2}+{(1-\tau)P_0\over P_\tau}\right],
\]

so choosing

\[
 a_0\le
 \left[{1\over2}+{(1-\tau)P_0\over P_\tau}\right]^{-1}
 \tag{4.17}
\]

keeps \(0\le\widetilde G\le1\).  Hence even an exact mean constraint does
not prevent cusp relocation.  A prescribed-median constraint alone also
does not help: before the harmless shifts, zero is a median of both binary
functions (for the balanced function it is one of the medians).

### 4.2 Robustness

The inequality (4.11) is strict.  Therefore it persists if the central
binary comparator is replaced by any sequence \(F_j\) converging to \(F\)
strictly in weighted `BV` with

\[
 \|M(F_j)-M(F)\|_*\longrightarrow0,qquad
 D_{\tau,s}(F_j)\longrightarrow D_{\tau,s}(F).       \tag{4.18}
\]

Indeed, the fixed cusp competitor then pays only an \(o(1)\) nuclear
penalty, while the right side of (4.11) stays negative.  Thus the
phenomenon is not an artefact of literal characteristic functions.  This
observation does not assert that the particular sets (4.2) are the
near-Cheeger heat comparators used elsewhere; its role is narrower but
exact: aggregate matrix rank, even with zero matrix-fidelity error and a
constant nuclear anisotropy, contains no information about the level label
of that matrix.  Any theorem specialized to the small-deficit heat
comparator would need an additional anti-relocation input not present in
(2.9).

## 5. Audits of the proposed repairs

### 5.1 Exact value-law constraints

Exact equimeasurability can remove the literal cusp plateau.  If
\(v_F(r)=\mu(F>r)\) satisfies

\[
             \mathcal L^1\{r:v_F(r)=\tau\}=0,        \tag{5.1}
\]

then every equimeasurable \(G\) satisfies

\[
 \int_{\{r:v_G(r)=\tau\}}M(\{G>r\})\,dr=0.          \tag{5.2}
\]

This is immediate because the trace of the left side is an integral over a
Lebesgue-null set.  However, equimeasurable competitors preserve the volume
of each level.  They give only volume-constrained stability, not stability
of

\[
                        P_\Phi(E)-\lambda\mu(E)      \tag{5.3}
\]

under all compactly supported variations.

The distinction is genuine even in the smoothest model.  A Gaussian
halfspace is a global perimeter minimizer at its prescribed Gaussian
volume.  Its weighted Jacobi form for a constant normal variation is

\[
                         Q(1)=-P_\gamma(E)<0,         \tag{5.4}
\]

because its second fundamental form vanishes and
\(\nabla^2V[n,n]=1\).  It is stable under zero-volume first variations but
not under the constant variation used in the flatness argument.  Thus
(5.2) does not supply exact unconstrained stability.

There is a second contact issue.  If one physical jump has traces \(a<b\)
and the volumes \(v_G(r)\), \(a<r<b\), pass from one arm of the tent to the
other, moving that interface changes all of those levels simultaneously.
The one cusp level has zero coarea measure, but the kink produces a
one-sided contact term in the coherent interface variation.  Merely saying
that almost every individual level is different from \(\tau\) does not
remove this contact.  One must keep the entire level interval carried by
each jump on one affine arm.

A hard quantile constraint only moves the same problem.  For example,

\[
       \mu(G=1)\ge\tau+\eta                             \tag{5.4a}
\]

would force every nontrivial superlevel volume onto the right affine arm.
The admissible class in (5.4a) is closed under \([0,1]\)-valued strong
\(L^1\) convergence but is highly nonconvex (average two functions whose
unit plateaux have small intersection).  It therefore also invalidates the
convex interpolation/minimax argument which produced one compatible
constant matrix \(H\).  More importantly, if equality
holds at a minimizer,
the top jump admits only one-sided changes of its enclosed volume.  It may
carry the entire matrix and has no unconstrained second variation.  Replacing
\(\ge\) by equality gives exactly volume-constrained stability; replacing it
by a strict inequality destroys compactness of the minimization problem.
The product construction in Section 4 can place its full matrix on a set of
any prescribed volume in \((0,1/2]\), so this is a genuine contact branch,
not a vacuous endpoint case.

### 5.2 Smooth concave regularization

Let \(B_\epsilon\le m\) be \(C^2\) and concave, with the cusp rounded on a
scale \(\epsilon\).  Suppose a smooth set \(E\), of volume \(v\), locally
minimizes

\[
                 P_{\Phi,\mu}(E)-\psi B_\epsilon(\mu(E)).
 \tag{5.5}
\]

The first variation gives multiplier
\(\lambda=\psi B_\epsilon'(v)\).  If \(u\) is the Wulff-normal speed and
\(Q_\lambda(u)\) is the Jacobi form of
\(P_{\Phi,\mu}-\lambda\mu\), the exact second variation is

\[
 Q_\lambda(u)-\psi B_\epsilon''(v)
       \left(\int_{\partial E}u\,d\sigma_{\Phi,\mu}\right)^2
 \ge0.                                                       \tag{5.6}
\]

Because \(B_\epsilon''\le0\), the second term in (5.6) is nonnegative.  It
can conceal a negative value of \(Q_\lambda(1)\); (5.6) therefore does not
imply unconstrained stability.  More precisely, if \(B_\epsilon=m\)
outside \([1/2-\epsilon,1/2+\epsilon]\), its derivative changes from
\(1\) to \(-1\) across an interval of length \(2\epsilon\), so

\[
                 \|B_\epsilon''\|_\infty\ge {1\over\epsilon}
 \tag{5.7}
\]

for every such rounding.  Although its scalar approximation error is
\(O(\epsilon)\), its second-variation error is not infinitesimal.  Exact
unconstrained stability is recovered only on an affine arm, where
\(B''=0\).

A bulk term must be distinguished from a gradient regularizer.  A term
\(\epsilon\int\phi(G)d\mu\) has the layer-cake representation

\[
 \epsilon\int\phi(G)d\mu
 =\epsilon\phi(0)+\epsilon\int_0^1
       \phi'(r)\mu(G>r)\,dr.                         \tag{5.7a}
\]

It therefore adds a level-dependent *linear* volume multiplier and leaves
the matrix anisotropy spatially constant.  But adding a linear function of
\(v\) does not remove the derivative jump of \(B_{\tau,s}(v)\) at
\(v=\tau\).  At a fixed cusp, an infinitesimal bulk term is also defeated
quantitatively by (5.8) below.  A nonlinear penalty on the value law which
really forces \(v_G(r)\) toward \(v_F(r)\) contributes its own second
derivative in \(v\), hence a rank-one term of exactly the form (5.6); a
piecewise-linear law penalty merely introduces additional contact cusps.

A Dirichlet or superlinear gradient regularizer is different: it is not a
one-homogeneous coarea integral.  Its Euler equation contributes a
gradient-dependent surface term, so the nuclear subgradient is no longer
the sole anisotropy in the flatness calculation.  Sending either kind of
regularization to zero does not, by itself, prevent the limiting derivative
matrix from reconcentrating on a cusp.

### 5.3 Lowering the cusp, total-variation penalties, and lexicographic limits

Replacing \(B_{\tau,1}\) by \(B_{\tau,s}\), \(s<1\), gives cusp charge
\(1-s\), but (2.3) shows that the central charge is larger.  The product
example of Section 4 still strictly prefers the cusp, since the factor
\(s\) cancels from (4.10).  Adding \(\epsilon\,{\rm TV}(G)\), or any
functional of \(M(G)\) alone, also has no effect on that comparison:
\(F\) and \(G\) have identical total variation and identical full matrix.

For a bounded auxiliary regularizer \(R\), the strict improvement in
(4.12) persists whenever

\[
 \epsilon|R(G)-R(F)|
 <\psi T\,{s\tau\over P_\tau}
 \left[1-{(1-\tau)^{-1/n}-1\over2^{1/n}-1}\right].   \tag{5.8}
\]

Thus a genuinely infinitesimal tie-breaker cannot remove the branch at a
fixed asymmetric cusp.

Finally, let \(\tau_j\to1/2\) and minimize lexicographically or pass to a
sequence of asymmetric minimizers.  Section 4 supplies, for every
\(\tau_j\), a competitor whose entire matrix is at the \(\tau_j\)-cusp.
The distance of that matrix from the central volume tends to zero, but its
fraction does not: it remains one.  The limit is precisely the balanced
single-jump obstruction.  Hence no uniform separation survives the
lexicographic limit.

## 6. What would actually be needed

The asymmetric construction proves only the following conditional
statement: if one separately knows that every jump interval carrying a
fixed positive fraction of \(M(G)\) stays inside one affine arm of
\(B_{\tau,s}\), then that fraction has exact unconstrained stability and
the flatness/Jacobi argument applies to it.  Neither (2.9) nor any scalar
moment constraint proved above supplies that premise.

There are three possible ways to add it, none of which is furnished by the
asymmetric quotient itself.

1. A label-resolved matrix penalty, for example fidelity of
   \(\int\chi(r)M(\{G>r\})\,dr\) for enough functions \(\chi\), can forbid
   relocation.  Its active matrix is generally \(H(r)\), not one spatially
   constant \(H\), and its derivatives reintroduce the selector terms that
   the unweighted construction was designed to remove.
2. Exact value-law constraints forbid a cusp plateau under (5.1), but a new
   theorem must upgrade volume-constrained stability to unconstrained
   stability.  The Gaussian calculation (5.4) shows that no such upgrade is
   formal.
3. One may accept and analyze the high-rank cusp branch directly.  That is
   the balanced-jump inverse problem, not an asymmetric-quantile escape
   from it.

Accordingly, within the broad convex coarea class with aggregate nuclear
matrix fidelity and constant anisotropy, asymmetric quantiles do not put a
universal fraction of the retained normal matrix strictly away from every
tent cusp.  The obstruction is exact and high-rank, rather than a loss in
the heat-parameter bookkeeping.
