# Clean-room audit of the coarea transition cascade

## 0. Verdict and required corrections

The main conclusions of `coarea_transition_cascade.md` are correct in the
log-concave Euclidean setting, but not with the unrestricted opening phrase
"let `mu` be a probability with Cheeger constant `psi`".  Three conventions
are essential.

1. In the coarea formulas, `P_mu` must be the relaxed weighted relative
   `BV` perimeter on the affine hull of `mu`.  It is not the exterior
   Minkowski content of each individual superlevel set.  The two notions
   have the same *Cheeger infimum* for a log-concave Euclidean probability,
   but they need not agree set by set.
2. The hypothesis `||grad F||_infinity <= L_F` in the transition lemma must
   mean that the chosen representative of `F` is globally `L_F`-Lipschitz
   on the convex support.  A bounded weak gradient supplies such a
   representative here, but this implication is false on a general metric
   support without further hypotheses.
3. The two radius-`d/2` closed neighborhoods used in the proposed proof can
   meet when their mutual distance is exactly `d`.  The proof is repaired by
   using radii strictly smaller than `d/2` and then taking a limit.

With these corrections, (1.2)--(1.10), (2.1)--(2.8), and (3.2)--(3.3) are
valid.  The phrase following (3.3) should be read as a flux-weighted joint
average statement; it is not a pointwise assertion on every level surface.

There is also one unsupported model claim.  Gaussian maximum cuts and inner
boxes in a uniform cube are quantitatively far from the known halfspace or
coordinate competitors.  The same conclusion does **not** follow for the
maximum cut of a product of one-sided exponential laws: its balanced
perimeter tends to `(log 2)/2`, below the coordinate-cut perimeter `1/2`.
Without the exact isoperimetric profile of that product, near-extremality of
this example cannot be ruled out by the argument in the audited file.

The audit proves below a new sharp all-threshold strengthening of Lemma 2.1:
if `F` is `[0,1]`-valued, `L`-Lipschitz, `E F=1/2`, and `mu(S)=1/2`, then,
with `kappa=psi_mu/L`,

\[
 \boxed{
 U={1\over2}\int|F-\mathbf1_S|\,d\mu
 \ \ge\ {1\over4}-{1-e^{-\kappa/2}\over2\kappa}.}       \tag{0.1}
\]

The right side is interpreted as zero at `kappa=0` and is
`kappa/16+O(kappa^2)` at zero.  It is sharp for a clipped affine function
under a one-dimensional two-sided exponential law.  Thus the complete
nested cascade gives a stronger dimension-free transition obstruction, but
even exact near-minimality of every level does not by itself create a cross
term between normals.

## 1. Perimeter convention and exterior Minkowski content

Let `E` be the minimal affine hull of a non-point log-concave probability
`mu`, let `k=dim E>=1`, and identify `E` with `R^k`.  On the relative
interior `Omega` of its convex support,

\[
                         d\mu=w(x)\,dx,
 \qquad w=e^{-V},                                      \tag{1.1}
\]

where `V` is convex and `w` is positive and continuous on `Omega`.
The boundary of `Omega` is Lebesgue-null, hence `mu`-null.  For a set of
locally finite perimeter define the relative weighted perimeter

\[
 P_\mu(A)=\int_{\partial^*A\cap\Omega}w(x)\,
                    d\mathcal H^{k-1}(x).             \tag{1.2}
\]

Equivalently, this is the relaxation on `Omega` of
`int |grad h| dmu` over Lipschitz approximations to `1_A`.

In this setting the following equality of infima is the exact fact needed
below:

\[
 \psi_\mu
 =\inf_{0<\mu(A)<1}{P_\mu(A)\over
                  \min(\mu(A),1-\mu(A))}.             \tag{1.3}
\]

Here the left side is the exterior-Minkowski definition in the task.
To see the two directions, distance cutoffs give relaxed perimeter no
larger than lower exterior Minkowski content.  Conversely, relative strict
`BV` approximation gives smooth sets `A_j` with

\[
 \mu(A_j)\longrightarrow\mu(A),\qquad
 P_\mu(A_j)\longrightarrow P_\mu(A).
\]

For a smooth relative set, differentiation of its outer tubes gives
`mu^+(A_j)=P_mu(A_j)`.  Exhausting `Omega` before taking the strict
approximation handles portions approaching the boundary of the support.
This proves (1.3).  It also shows why one must not replace `P_mu(A_r)` by
`mu^+(A_r)` inside an exact coarea identity.

The distinction is not cosmetic.  For uniform measure on `[0,1]`, the set

\[
 (0,1/2)\cup(\mathbb Q\cap(1/2,1))
\]

is equal almost everywhere to `(0,1/2)` and has relaxed relative perimeter
one, but its closure is `[0,1]`; its exterior Minkowski content is infinite.
Thus exterior content is sensitive to the chosen Borel representative,
whereas the coarea perimeter is a `BV` equivalence-class quantity.

In particular, (1.3) implies, for every finite-perimeter set,

\[
 P_\mu(A)\ge\psi_\mu\min(\mu(A),1-\mu(A)).            \tag{1.4}
\]

The half-mass near-minimizer used later requires one further log-concave
fact; it does not follow from the definition of `psi` for an arbitrary
probability.  Define

\[
 \mathcal I_\mu(v)=\inf\{P_\mu(A):\mu(A)=v\}.
\]

For a log-concave Euclidean probability, the lower-semicontinuous
isoperimetric profile is concave on `(0,1)` and symmetric under
`v mapsto 1-v`.  In the smooth case this is the `CD(0,infinity)` profile
concavity theorem (in this normalization, E. Milman, *Inventiones
Mathematicae* 177 (2009), Theorem 1.8 and its approximation convention);
convolution/convex-potential approximation, exhaustion
of the convex support, and lower semicontinuity give exactly the present
nonsmooth relative version.  Concavity and `I_mu(0)=0` imply

\[
 {\mathcal I_\mu(v)\over\min(v,1-v)}
 \ge2\mathcal I_\mu(1/2).
\]

Together with (1.3), this gives

\[
 \psi_\mu=2\mathcal I_\mu(1/2).                       \tag{1.5}
\]

Consequently, for every `epsilon>0` there is a half-mass finite-perimeter
set with `P_mu(S)<=psi_mu/2+epsilon`; no minimizer is asserted to exist.

All gradients and distances in what follows are relative to `E`.  The
ambient and relative exterior-Minkowski Cheeger constants agree because
distance between points of `E` is unchanged.  Thus lower-dimensional
support introduces no additional term.

## 2. Reproof of the exact coarea cascade

Let `S` be a finite-perimeter Borel set with `mu(S)=1/2`.  Let
`F:Omega->[0,1]` be locally Lipschitz, with

\[
                              \int F\,d\mu={1\over2},
 \qquad
 U={1\over2}\int|F-\mathbf1_S|\,d\mu,                \tag{2.1}
\]

and put `A_r={F>r}`.

### 2.1 Equality of the two one-sided errors

Write

\[
 U_-:=\int_S(1-F)d\mu,\qquad
 U_+:=\int_{S^c}F d\mu.
\]

The equality of the means gives

\[
 U_--U_+=\mu(S)-\int Fd\mu=0.
\]

Their sum is `int|F-1_S|dmu=2U`, so

\[
                         U_-=U_+=U.                  \tag{2.2}
\]

### 2.2 Weighted coarea

The relative Euclidean `BV` coarea theorem, multiplied by the nonnegative
continuous weight `w`, gives the extended-real identity

\[
 \int|\nabla F|d\mu
     =\int_0^1P_\mu(A_r)dr.                           \tag{2.3}
\]

In particular `A_r` has finite relative weighted perimeter for almost every
`r` whenever the left side is finite.

### 2.3 Distance to constants

Let `m` be any median of `F`.  For `r<m`,
`mu(F<=r)<=1/2`; for `r>m`, `mu(F>r)<=1/2`.  Layer cake therefore gives

\[
\begin{aligned}
 \int_0^1\min\{\mu(A_r),1-\mu(A_r)\}\,dr
 &=\int_0^m\mu(F\le r)dr+\int_m^1\mu(F>r)dr\\
 &=\int|F-m|d\mu.
\end{aligned}                                        \tag{2.4}
\]

A median minimizes absolute deviation, so the last expression equals
`inf_c int|F-c|dmu`.  An atom at `m` changes only one threshold and hence
does not alter the integral.

Combining (2.3), (2.4), and (1.4) proves the exact nonnegative
decomposition

\[
\begin{aligned}
 &\int|\nabla F|d\mu-
   \psi_\mu\inf_c\int|F-c|d\mu\\
 &\quad=\int_0^1\left[P_\mu(A_r)-
  \psi_\mu\min\{\mu(A_r),1-\mu(A_r)\}\right]dr.      \tag{2.5}
\end{aligned}
\]

### 2.4 Symmetric-difference cascade

For every `x`, direct integration in `r` gives

\[
 \int_0^1|\mathbf1_{\{F(x)>r\}}-\mathbf1_S(x)|dr
 =|F(x)-\mathbf1_S(x)|.
\]

Tonelli's theorem and (2.1) therefore give

\[
               \int_0^1\mu(A_r\mathbin\triangle S)dr=2U.       \tag{2.6}
\]

For every constant `c`, the triangle inequality gives

\[
 \int|F-c|d\mu
 \ge\int|\mathbf1_S-c|d\mu-
        \int|F-\mathbf1_S|d\mu.
\]

The first term has infimum `1/2`, hence

\[
              \inf_c\int|F-c|d\mu\ge{1\over2}-2U.   \tag{2.7}
\]

Suppose

\[
 P_\mu(S)=p\le {\psi_\mu\over2}+\varepsilon,
 \qquad \int|\nabla F|d\mu\le p.                    \tag{2.8}
\]

Equations (2.5), (2.7), and (2.8) yield

\[
 \int_0^1\left[P_\mu(A_r)-
  \psi_\mu\min(\mu(A_r),1-\mu(A_r))\right]dr
 \le\varepsilon+2\psi_\mu U.                        \tag{2.9}
\]

Both integrands in (2.6) and (2.9) are nonnegative.  Markov's inequality
and the union bound therefore prove

\[
 \left|\left\{r:\ \mu(A_r\mathbin\triangle S)>\lambda_1
 \ \hbox{or}\
 P_\mu(A_r)-\psi_\mu\min(\mu(A_r),1-\mu(A_r))>\lambda_2
 \right\}\right|
 \le {2U\over\lambda_1}
 +{\varepsilon+2\psi_\mu U\over\lambda_2}.          \tag{2.10}
\]

This reproves (1.2)--(1.10) of the audited file with all perimeter
conventions explicit.

## 3. Reproof of the transition lower bound

We first isolate the only metric-measure input.

**Neighborhood-growth lemma.**  Let a probability on Euclidean space (or,
more generally, on a length space with the neighborhood-semigroup property)
have exterior-Minkowski Cheeger constant `psi`.  If `B` is closed and
`b=mu(B)<=1/2`, then

\[
 \mu(B_t)\ge b e^{\psi t}                              \tag{3.1}
\]

for every `t` before the left side first reaches `1/2`.  Once it reaches
`1/2`, it is strictly larger than `1/2` at every later positive radius,
unless it has already reached mass one.

Indeed, `h(t)=mu(B_t)` is right-continuous and
`(B_t)_r=B_{t+r}`.  At every `t` with `h(t)<1/2`, the definition of the
Cheeger constant gives

\[
 \liminf_{r\downarrow0}{h(t+r)-h(t)\over r}\ge\psi h(t).
\]

The lower-Dini-derivative comparison applied to `e^{-psi t}h(t)` proves
(3.1).  At a point where `h=1/2`, the same lower derivative is at least
`psi/2`, which proves the strict statement.  If `psi=0`, all conclusions
used below are immediate by continuity in `psi`.

Now assume that the representative of `F` is globally `L_F`-Lipschitz.
Fix `a in (0,1/2)` and set

\[
 L_a=\{F\le a\},\qquad H_a=\{F\ge1-a\},
 \qquad m={1\over2}-{U\over a}.                       \tag{3.2}
\]

From (2.2) and Markov's inequality,

\[
 \mu(L_a)\ge m,\qquad \mu(H_a)\ge m.                 \tag{3.3}
\]

For example, on `S^c setminus L_a` one has `F>a`, whose mass is at most
`U/a`; the estimate for `H_a` is identical on `S`.

If `m<=0`, then `U>=a/2` and the desired inequality is automatic.  If
`m>0`, Lipschitzness gives

\[
 \operatorname {dist}(L_a,H_a)\ge d:={1-2a\over L_F}.           \tag{3.4}
\]

If `m exp(psi_mu d/2)>1/2`, the neighborhood-growth lemma says that, for
some radius `r<d/2` sufficiently close to `d/2`, both `(L_a)_r` and
`(H_a)_r` have mass strictly larger than `1/2`.  They are disjoint because
`2r<d`, a contradiction.  Hence

\[
                         m e^{\psi_\mu d/2}\le{1\over2}.         \tag{3.5}
\]

Substitution of (3.2) and (3.4) yields

\[
 \boxed{
 U\ge {a\over2}\left(1-
       e^{-\psi_\mu(1-2a)/(2L_F)}\right).}            \tag{3.6}
\]

This is (2.1).  Notice that the use of `r<d/2` removes the closed-neighborhood
gap in the original proof.

The middle layer splits between the two error sets:

\[
 \{a<F<1-a\}\subset
 \bigl(S\cap\{F<1-a\}\bigr)
 \cup\bigl(S^c\cap\{F>a\}\bigr).
\]

Each part has mass at most `U/a`; consequently

\[
                         \mu(a<F<1-a)\le {2U\over a}.            \tag{3.7}
\]

## 4. Heat regularization and the sharp Lipschitz constant

This section states explicitly the object denoted by `F` in Sections 2--3
of the audited file.  Let `f=1_S`, let `Z_s` be a centered Gaussian of
covariance `s I_E`, independent of `X~mu`, and put `Y=X+Z_s`.  Define

\[
 g_s(y)=\mathbb E[f(X)\mid Y=y],\qquad
 F_s(x)=\mathbb E[g_s(x+Z_s)],                        \tag{4.1}
\]

where the second expectation is in a fresh channel noise.  Thus
`F_s=T_s f` for the posterior-resampling Markov operator.

The posterior law of `X` given `Y=y` is `s^{-1}`-strongly log-concave on
`E`.  If `g=g_s(y)` and

\[
 v=\operatorname {Cov}(f(X),X\mid Y=y),
\]

differentiating the Gaussian likelihood gives

\[
                         \nabla g_s(y)={v\over s}.     \tag{4.2}
\]

The sharp centroid inequality for a `t`-strongly log-concave probability is

\[
 |\operatorname {Cov}(\mathbf1_B,X)|
       \le {I(\pi(B))\over\sqrt t},                  \tag{4.3}
\]

where `I(r)=phi(Phi^{-1}(r))`.  Here is a complete reduction.  If the
covariance vector is nonzero, project the centered probability onto its
direction and call the resulting real variable `Q`.  Strong log-concavity
is preserved under this marginal.  If
`q(z)=P(B|Q=z)`, then the bathtub principle gives

\[
 |\operatorname {Cov}(\mathbf1_B,X)|
 =\mathbb E[Qq(Q)]
 \le\mathbb E[Q\mathbf1_{\{Q\ge c\}}],               \tag{4.3a}
\]

where the upper tail has mass `pi(B)`.  Let `T` be the increasing transport
from a standard Gaussian `Z` to `Q`.  The one-dimensional Caffarelli
contraction gives `T'<=1/sqrt(t)`.  Hence
`R(z)=T(z)-z/sqrt(t)` is nonincreasing.  Since `E Q=E Z=0`, one has
`E R(Z)=0`; negative association of the decreasing `R` and the increasing
upper-tail indicator gives

\[
 \mathbb E[Q\mathbf1_{\{Q\ge c\}}]
 \le {1\over\sqrt t}\mathbb E[Z\mathbf1_{\{Z\ge z_g\}}]
 ={I(\pi(B))\over\sqrt t}.                            \tag{4.3b}
\]

Here `z_g=Phi^{-1}(1-pi(B))`, so the Gaussian upper tail has the same
mass as `B` and `phi(z_g)=I(pi(B))`.

This proves (4.3); the zero covariance-vector case is immediate.
Nonsmooth convex potentials follow by monotone convex approximation and
weak convergence; all quantities in (4.3) are bounded by second moments.
Applying (4.3) with `t=1/s` in (4.2) gives

\[
 |\nabla g_s(y)|\le {I(g_s(y))\over\sqrt s}
                \le {I(1/2)\over\sqrt s}.            \tag{4.4}
\]

Dominated differentiation in (4.1) and conditional Jensen now prove the
claimed heat-smoothing bound, with its exact constant:

\[
 \boxed{\|\nabla F_s\|_\infty
            \le {I(1/2)\over\sqrt s}.}               \tag{4.5}
\]

The Markov operator preserves the mean, so `E F_s=1/2`.  Self-adjointness
of posterior resampling gives

\[
\begin{aligned}
 \mathbb E|f-F_s|
 &=1-2\langle f,T_sf\rangle\\
 &=1-2\mathbb E g_s^2
 =2\mathbb E[g_s(1-g_s)].
\end{aligned}                                        \tag{4.6}
\]

Thus the `U` of (2.1) is exactly

\[
                         U(s)=\mathbb E[g_s(1-g_s)].   \tag{4.7}
\]

Substitution of (4.5) in (3.6) gives (2.7).  If
`psi_mu sqrt(s)<=1`, the elementary inequality
`1-e^{-x}>=c x` on the resulting fixed compact interval gives (2.8), with
a positive constant depending only on `a`.

For completeness, the heat-profile upper bound can be checked without any
boundary smoothness.  Write

\[
 q_s=P_s\mu,\qquad h_s=P_s(f\mu),\qquad
 z_s=\Phi^{-1}(g_s),\qquad
 J(s)=\int q_s I(g_s),\qquad H(s)={J(s)\over\sqrt s}.
\]

The two convolutions solve the heat equation with generator `Delta/2`.
A direct quotient differentiation gives

\[
 \partial_s(q_sI(g_s))
 ={1\over2}\Delta(q_sI(g_s))
  +{1\over2}q_sI(g_s)|\nabla z_s|^2.                 \tag{4.8a}
\]

The centroid inequality (4.3) is exactly
`s|grad z_s|^2<=1`.  Integration of (4.8a), first with spatial cutoffs and
then by monotone domination, yields

\[
 H'(s)=-{1\over2s^{3/2}}\int q_sI(g_s)
 (1-s|\nabla z_s|^2)\le0.        \tag{4.8b}
\]

One rigorous cutoff procedure is to replace `f` by
`epsilon+(1-2epsilon)f`; then `Phi^{-1}(g_s)` stays in a compact interval.
Gaussian convolution supplies smoothness, and the score estimate
`int q_s|grad log q_s|^2<=k/s` removes spatial cutoffs.  Sending
`epsilon downarrow0` uses `0<=I(g)<=I(1/2)` and
`0<=s|grad z|^2<=1`.  Thus (4.8b) remains valid for the original Borel
label and lower-dimensional support.

Weighted `BV` blow-up, or strict approximation first for a smooth cut,
gives `lim_{s downarrow0}H(s)=P_mu(S)=p`.  Finally, concavity and symmetry
of the Gaussian profile give
`I(r)>=sqrt(2/pi)r(1-r)`.  Equations (4.7)--(4.8b) therefore prove

\[
 U(s)\le {\sqrt s\,p\over\sqrt{2/\pi}}.              \tag{4.8}
\]

The same pointwise centroid bound gives the estimate used below:

\[
 \mathbb E|W|=\int q_s I(g_s)|\nabla z_s|
       \le {J(s)\over\sqrt s}=H(s)\le p.             \tag{4.8c}
\]

Consequently, if `p<=psi_mu/2+epsilon` with `epsilon=o(psi_mu)`, the lower
and upper transition bounds have the same `psi_mu sqrt(s)` order.  The
phrase `p=psi/2+o(psi)` in the audited file is valid only after explicitly
choosing such a near-minimizing sequence.

### 4.1 Jensen direction identity

Let

\[
 W=\nabla g_s(Y),\qquad m(X)=\mathbb E[W\mid X]
                              =\nabla F_s(X).          \tag{4.9}
\]

Set `theta=m/|m|` where `m!=0`, choosing it arbitrarily where `m=0`, and
set `u=W/|W|` where `W!=0`.  Then

\[
\begin{aligned}
 \mathbb E\{|W|\langle u,\theta\rangle\mid X\}
 &=\langle\mathbb E[W\mid X],\theta\rangle\\
 &=|m(X)|.
\end{aligned}
\]

This remains true on `{m=0}` because the conditional vector mean vanishes.
Therefore

\[
\begin{aligned}
 \delta_J
 &:=\mathbb E|W|-\mathbb E|m(X)|\\
 &=\mathbb E\{|W|(1-\langle u,\theta\rangle)\}
 ={1\over2}\mathbb E\{|W|\,|u-\theta|^2\}.           \tag{4.10}
\end{aligned}
\]

Weighted coarea and the `L^1` Cheeger inequality give

\[
 \mathbb E|m(X)|=\int|\nabla F_s|d\mu
 \ge\psi_\mu(1/2-2U(s)).                             \tag{4.11}
\]

The heat-profile upper bound gives `E|W|<=p`.  Since
`p>=psi_mu/2` by (1.4),

\[
\begin{aligned}
 0\le\delta_J
 &\le p-\psi_\mu(1/2-2U(s))\\
 &\le\varepsilon+4pU(s).                             \tag{4.12}
\end{aligned}
\]

This proves (3.2)--(3.3).  It aligns directions in the joint law of
`(X,Y)`, weighted by `|W|`.  A pointwise or uniform assertion on an
individual level surface requires an additional disintegration or Markov
estimate and is not contained in (4.10).

All formulas above are relative to `E`.  If one instead convolves in the
ambient space, the normal Gaussian factor cancels from the quotient
`g_s`, so `g_s` is constant in directions perpendicular to `E`; (4.2)--
(4.12) are unchanged.  Approximation of a finite-perimeter `S` by smooth
relative sets, followed by lower semicontinuity, supplies the only limit
needed.  No regularity of the original boundary is assumed.

## 5. A sharp all-threshold consequence

Lemma 2.1 uses only two thresholds.  Using every threshold gives a strictly
stronger and sharp result.

**Theorem 5.1 (sharp multilevel transition bound).**  Let a metric
probability have Cheeger constant `psi>0`.  Let `F` be a `[0,1]`-valued,
globally `L`-Lipschitz function on its support, with
`E F=1/2`.  For every set `S` of mass `1/2`, put

\[
 U={1\over2}\mathbb E|F-\mathbf1_S|,
 \qquad \kappa={\psi\over L}.
\]

Then

\[
 \boxed{
 U\ge R(\kappa):={1\over4}-
              {1-e^{-\kappa/2}\over2\kappa}.}        \tag{5.1}
\]

At `kappa=0`, define `R(0)=0`.  In particular

\[
 R(\kappa)={\kappa\over16}+O(\kappa^2)
 \quad(\kappa\downarrow0).                           \tag{5.2}
\]

**Proof.**  First, Lipschitz contraction gives

\[
                         \psi_{F_\#\mu}\ge {\psi_\mu\over L}.   \tag{5.3}
\]

Indeed, for every Borel `B subset R`,

\[
 (F^{-1}B)_\delta\subset F^{-1}(B_{L\delta}).
\]

Divide the corresponding mass increment by `L delta`, take a lower limit,
and apply the Cheeger inequality to `F^{-1}B`.

Let `nu=F_#mu` and let

\[
                         G(r)=\nu(({-\infty},r]),\qquad0<r<1.
\]

Choose a median `m in [0,1]`.  The neighborhood-growth lemma for `nu` and
(5.3) imply

\[
\begin{aligned}
 G(r)&\le {1\over2}e^{-\kappa(m-r)},&&0<r<m,\\
 1-G(r)&\le {1\over2}e^{-\kappa(r-m)},&&m<r<1.
\end{aligned}                                        \tag{5.4}
\]

For example, expand the half-line `(-infinity,r]` up to radii tending to
`m-r`.  It remains of mass at most `1/2` before the median, so (3.1)
forces the first inequality.  Apply the same argument to the upper tail
for the second.  This proof is unaffected by an atom at the median.

Set

\[
 A_-:=\int_0^m(1/2-G(r))dr,
 \qquad A_+:=\int_m^1(G(r)-1/2)dr.                   \tag{5.5}
\]

Since `E F=int_0^1(1-G(r))dr=1/2`, one has
`A_-=A_+`.  From (5.4),

\[
\begin{aligned}
 A_-&\ge h_\kappa(m),\\
 A_+&\ge h_\kappa(1-m),\\
 h_\kappa(t)&:={t\over2}-{1-e^{-\kappa t}\over2\kappa}.
\end{aligned}                                        \tag{5.6}
\]

For every threshold,

\[
 \mu(A_r\mathbin\triangle S)
 \ge|\mu(A_r)-1/2|=|G(r)-1/2|.                       \tag{5.7}
\]

Integrating and using (2.6) yields

\[
 2U\ge A_-+A_+=2A_-,
 \qquad U\ge\max\{h_\kappa(m),h_\kappa(1-m)\}.      \tag{5.8}
\]

The function `h_kappa` is increasing.  The minimum over `m` of the last
maximum occurs at `m=1/2` and equals the right side of (5.1).  This proves
the theorem.

### 5.1 Sharpness and equality model

Let

\[
 d\mu_\kappa(x)={\kappa\over2}e^{-\kappa|x|}dx
\]

be the two-sided exponential law.  Its one-dimensional isoperimetric
profile is

\[
                   \mathcal I_{\mu_\kappa}(v)
       =\kappa\min(v,1-v),                            \tag{5.9}
\]

so `psi_{mu_kappa}=kappa`.  Put

\[
 S=[0,\infty),\qquad
 F(x)=\min\{1,\max\{0,1/2+x\}\}.                    \tag{5.10}
\]

Then `F` is one-Lipschitz, `E F=1/2`, and every nontrivial superlevel set
of `F` is a half-line attaining equality in (5.9).  Direct integration
gives

\[
\begin{aligned}
 {1\over2}\mathbb E|F-\mathbf1_S|
 &= {1\over2}\left\{
 {1\over2}-{1-e^{-\kappa/2}\over\kappa}\right\}\\
 &=R(\kappa).                                        \tag{5.11}
\end{aligned}
\]

Thus (5.1) is sharp.  More importantly for the proposed proof strategy,
this example has

\[
 P_\mu(A_r)-\psi_\mu\min(\mu(A_r),1-\mu(A_r))=0
\]

for every `r in (0,1)`, not merely almost every `r`.  Exact nested
near-minimality alone therefore cannot yield a strict improvement of the
transition estimate.

### 5.2 Finite extraction from the near-minimal cascade

There is nevertheless an exact finite consequence of (2.6) and (2.9).
Put

\[
 D=\int_0^1\left[P_\mu(A_r)-
      \psi_\mu\min(\mu(A_r),1-\mu(A_r))\right]dr.
\]

For `lambda,eta>0`, let `G_{lambda,eta}` be the set of thresholds satisfying

\[
 \mu(A_r\mathbin\triangle S)\le\lambda,
 \qquad P_\mu(A_r)-\psi_\mu\min(\mu(A_r),1-\mu(A_r))\le\eta.
\]

Then

\[
 |G_{\lambda,\eta}|
 \ge1-{2U\over\lambda}-{D\over\eta}.                \tag{5.12}
\]

A maximal `h`-separated subset of `G_{lambda,eta}` has cardinality at
least

\[
 {1\over2h}\left(1-{2U\over\lambda}-{D\over\eta}\right)_+.     \tag{5.13}
\]

The corresponding sets are nested, all are `lambda`-close to the same
half-mass set, and every consecutive pair of levels separated by `h`
satisfies

\[
 \operatorname {dist}(\{F\ge r_{i+1}\},\{F\le r_i\})
                       \ge {h\over L}.               \tag{5.14}
\]

Thus a near-minimizing heat regularization really does supply a large
finite family of separated, nested near-minimizers.  Summing the Cheeger
growth constraints across this family is exactly the multilevel argument
in Theorem 5.1.

### 5.3 Why threshold-only finite competitors have no angular cross term

There is a formal obstruction to extracting more from the nested family
without using ambient geometry.  Suppose, temporarily, that `w` and `F`
are smooth and `t_1<...<t_N` are regular values.  If `B subset (0,1)` is a
finite union of intervals whose boundary is contained in `{t_i}`, then

\[
 P_\mu(\{F\in B\})
 =\sum_{t_i\in\partial B}\int_{\{F=t_i\}}w\,
                         d\mathcal H^{k-1}
 =\sum_{t_i\in\partial B}P_\mu(A_{t_i}).             \tag{5.15}
\]

The reduced-boundary pieces are disjoint level surfaces, so no inner
product between their normals occurs.  Unions and intersections of the
nested upper sets themselves reduce to one of the selected sets; more
general Boolean combinations only produce the level bands in (5.15).
Approximation extends the inequality form needed for nonsmooth levels.

Therefore any finite competitor measurable solely with respect to `F`
uses only the scalar law of `F` and the scalar perimeter of each level.  It
cannot distinguish a constant normal field from a rotating one.  A genuine
normal-incidence inequality must splice different *spatial patches* of
level surfaces, use a second variation, or use a log-concave midpoint
overlap mechanism.  This is precisely the missing cross term identified in
the audited file.

## 6. Model tests

### 6.1 Gaussian halfspace

Let `mu=gamma_n` and `S={x_1>=0}`.  A direct Gaussian regression
calculation gives

\[
 g_s(y)=\Phi\left({y_1\over\sqrt{s(1+s)}}\right),
 \qquad
 F_s(x)=\Phi\left({x_1\over\sqrt{s(s+2)}}\right).     \tag{6.1}
\]

Thus the actual Lipschitz constant is

\[
 L_s={I(1/2)\over\sqrt{s(s+2)}},
\]

which is below the universal bound (4.5), and every level normal is the
same vector.  Since `psi_gamma=2I(1/2)`, the parameter in Theorem 5.1 is
`kappa=2sqrt{s(s+2)}`.  Moreover, the two posterior-resampled Gaussian
copies have correlation `1/(1+s)`, so

\[
 U(s)={1\over4}-{1\over2\pi}
             \arcsin\left({1\over1+s}\right).        \tag{6.2}
\]

Both (6.2) and (5.1) have the `sqrt(s)` small-time scale, and both tend to
`1/4` at large time.  No constant degenerates.

### 6.2 Gaussian median balls

For a centered Gaussian median ball, every smoothed level is concentric
and the normal field is radial.  Gaussian isoperimetry says that a
half-mass ball has perimeter strictly larger than `I(1/2)`, the halfspace
perimeter.  Thus it does not satisfy vanishing `epsilon` in (2.8), but it
does show that any ambient rigidity theorem must retain a radial branch.
In high dimension the radial density at its median tends to `1/sqrt(pi)`,
whereas `I(1/2)=1/sqrt(2pi)`, so the excess stays of constant order.

### 6.3 Uniform cubes

For the uniform law on `[-a,a]^n`, let

\[
 S=[-a,t]^n,
 \qquad b={t+a\over2a}=2^{-1/n}.
\]

Then `mu(S)=1/2` and its relative perimeter is

\[
 p_{\rm box}={n\over2a}b^{n-1}={n\over4ab}.          \tag{6.3}
\]

A coordinate half-cube has perimeter `p_coord=1/(2a)`, so the half-mass
isoperimetric profile is at most `p_coord`, while

\[
 {p_{\rm box}\over p_{\rm coord}}={n\over2b}.         \tag{6.4}
\]

The many-facet inner box is therefore quantitatively non-extremal (by a
factor of order `n`), whereas coordinate slabs have a constant normal.

### 6.4 Gaussian and exponential products

For a Gaussian maximum cut `S=(-infinity,t]^n` with `Phi(t)^n=1/2`,

\[
 p=n\varphi(t)\Phi(t)^{n-1}
   ={n\varphi(t)\over2\Phi(t)}\asymp\sqrt{\log n}.   \tag{6.5}
\]

The Gaussian halfspace competitor has perimeter `I(1/2)`, so this phase
example is far from globally minimizing.

For independent one-sided `Exp(1)` coordinates, the analogous inner box
has `(1-e^{-t})^n=1/2` and

\[
 p=n e^{-t}(1-e^{-t})^{n-1}
   ={n e^{-t}\over2(1-e^{-t})}
   \longrightarrow {\log2\over2}.                   \tag{6.6}
\]

A coordinate median cut has perimeter `1/2`.  Hence (6.6), unlike (6.3)
and (6.5), does not prove non-extremality; it is even the better of these
two explicit competitors.  The claim that exponential maximum cuts are
automatically excluded by (2.9) is therefore not established.

### 6.5 One-dimensional exponential laws

For the one-sided exponential density `e^{-x}1_{x>=0}`, the Cheeger
constant is one.  Every upper tail of mass at most `1/2` has perimeter
equal to its mass and is an exact Cheeger set.  Thus there is already a
continuum of nested exact levels on one side of the median.

The symmetric exponential example (5.9)--(5.11) is stronger: *every*
nontrivial level of the clipped affine `F` is an exact Cheeger set, the
coarea deficit (2.5) vanishes identically in the threshold variable, and
the all-threshold lower bound is attained.  This is the exact equality
test for the cascade.

## 7. Consequence for the proof search

The corrected cascade gives two rigorous pieces of information:

* near-global minimality supplies many nested levels that are simultaneously
  close to the same half-mass set and close to the Cheeger profile;
* the entire nested family forces the sharp transition cost (5.1), improving
  the two-threshold constant from `kappa/32` (optimized at `a=1/4`) to
  `kappa/16` at small `kappa`.

For the heat regularization, however,

\[
 \kappa={\psi_\mu\sqrt s\over I(1/2)},
 \qquad U(s)\le {\sqrt s\,p\over\sqrt{2/\pi}},        \tag{7.1}
\]

and, when `p` is close to `psi_mu/2`, the lower and upper sides of (7.1)
remain compatible by universal constants.  The equality example proves
that no threshold-only finite competitor can improve this comparison.
The unresolved load-bearing step is therefore genuinely spatial: it must
penalize incompatible normals on different boundary patches.  None of the
audited identities supplies that penalty, and no KLS-equivalent assertion
has been inserted in this report.
