# Near-Cheeger CMC leaves: secant control and orientation-free Wulff amplification

## 0. Verdict

There is a valid repair in the smooth positive-lapse setting.  It does not
require a leaf to minimize perimeter at its volume, and it does not use
concavity of the global anisotropic isoperimetric profile.

Let

\[
                    m(v)=\min(v,1-v),\qquad B(v)=c\,m(v),
\]

where `c` is any global anisotropic Cheeger lower bound.  Suppose a nested
constant-anisotropy CMC foliation has perimeter `P(v)`, positive Wulff
lapse, and therefore concave `P`.  If one leaf satisfies

\[
                 B(v_0)\le P(v_0)\le(1+\delta)B(v_0),       \tag{0.1}
\]

and the foliation exists on

\[
 [v_0-\rho m(v_0),v_0+\rho m(v_0)],\qquad 0<\rho<1,         \tag{0.2}
\]

then its CMC multiplier `lambda` obeys

\[
                         |\lambda|\le
                \left(1+{\delta\over\rho}\right)c.         \tag{0.3}
\]

The central cusp at volume `1/2` causes no failure in (0.3).  To amplify the
leaf, orient its boundary as the boundary of either the leaf or its
complement so that the outward multiplier is nonpositive.  Increase the
chosen side's volume by

\[
                         \gamma m(v_0),\qquad0<\gamma<1.     \tag{0.4}
\]

The killed Wulff tube then has length comparable to `1/c`, loses at most
`(delta+gamma)/(1+delta)` of its exponentially normalized base flux, and
gives an explicit covariance inequality.  This handles volumes on both
arms and volumes arbitrarily close to `1/2` in one statement.  The small
number `gamma` is fixed once and for all; thus a length `gamma/c` is still
a dimension-free multiple of `1/c`.

The direct, baseline-relative **Euclidean** coarea deficit of the heat
construction is exactly the right input for selecting such leaves.  The
joint replacement cannot increase Euclidean total variation, so Markov's
inequality deletes at most `D/epsilon` of Euclidean surface mass.  The
constant anisotropy is converted only after this deletion; its relative
loss is `1+O(kappa)`.  No `O(kappa)/epsilon` integrated loss occurs.

Two gaps remain.  For the actual `BV` equimeasurable minimizer one still
needs a singular theorem giving concavity of the **chosen nested
foliation's** perimeter across a fixed surrounding volume interval, not
concavity of the global isoperimetric profile.  One also needs singular
Wulff-tube coverage, including hard-support contact.  Smooth a.e. CMC
stationarity alone does not supply either theorem.

## 1. Smooth anisotropic setup

Let

\[
 d\mu=Z^{-1}e^{-V}1_\Omega dx,
\]

where `Omega` is a smooth convex domain and `V` is smooth and convex.  Let
`Phi` be positive, even, convex, one-homogeneous, smooth off zero, and
uniformly elliptic on tangent spaces.  For a unit normal `N`, write

\[
 z(N)=D\Phi(N),\qquad d\sigma_\Phi=\Phi(N)d\sigma_\mu.         \tag{1.1}
\]

Assume that `{E_v}_{v in J}` is a nested smooth family with
`mu(E_v)=v`.  Its boundaries have constant weighted anisotropic mean
curvature `lambda(v)` and are generated in Wulff gauge by

\[
                         \partial_vX=fz(N),\qquad f>0,
 \qquad \int_{\partial E_v}f\,d\sigma_\Phi=1.                  \tag{1.2}
\]

At a hard wall the natural anisotropic Young condition is assumed.  Put

\[
                         P(v)=P_\Phi(E_v).                      \tag{1.3}
\]

The first variation and the positive-lapse Jacobi identity give

\[
                         P'(v)=\lambda(v),\qquad
                         \lambda'(v)\le0.                       \tag{1.4}
\]

For completeness, the sign in (1.4) follows from

\[
 -\lambda'(v)\int_{\partial E_v}{1\over f}\,d\sigma_\Phi
 =\int\langle A_\Phi\nabla\log f,\nabla\log f\rangle
       d\sigma_\Phi
  +\int q_\Phi\,d\sigma_\Phi
  +\int_{\partial(\partial E_v)}b_\Phi\,d\tau_\Phi,           \tag{1.5}
\]

where every term on the right is nonnegative under convexity and the Young
condition.  Thus `P` is concave on every interval on which the positive
lapse foliation is smooth.  No minimality of any `E_v` is used here.

Let `c>0` satisfy the global inequality

\[
                         P_\Phi(A)\ge c\,m(\mu(A))              \tag{1.6}
\]

for every finite-perimeter set `A`.  It is permitted, and most natural, to
take `c` to be the anisotropic Cheeger constant.  Equation (1.6), rather
than any regularity or concavity of the global profile, is the only global
isoperimetric input below.

## 2. A single near-baseline leaf controls its multiplier

The following scalar lemma is the key observation.

**Lemma 2.1 (secants above a lower baseline).**  Let `P` be concave on
`[a,b]`, differentiable at `v_0 in (a,b)`, and suppose

\[
                 P(v)\ge B(v)\quad(a\le v\le b),\qquad
                 0\le e_0:=P(v_0)-B(v_0).                       \tag{2.1}
\]

Then

\[
 {B(b)-B(v_0)-e_0\over b-v_0}
 \le P'(v_0)\le
 {B(v_0)+e_0-B(a)\over v_0-a}.                                \tag{2.2}
\]

**Proof.**  Concavity puts the derivative between the right and left
secants:

\[
 {P(b)-P(v_0)\over b-v_0}\le P'(v_0)
 \le {P(v_0)-P(a)\over v_0-a}.
\]

Use `P(a)>=B(a)`, `P(b)>=B(b)`, and `P(v_0)=B(v_0)+e_0`.
QED.

The tent `B(v)=c m(v)` is `c`-Lipschitz, including across its cusp.  Take

\[
 a=v_0-\rho m_0,\qquad b=v_0+\rho m_0,
 \qquad m_0=m(v_0),                                    \tag{2.3}
\]

and suppose `e_0<=delta c m_0`.  Lemma 2.1 gives

\[
                    -c-{e_0\over\rho m_0}
       \le\lambda(v_0)\le
                     c+{e_0\over\rho m_0},                     \tag{2.4}
\]

hence

\[
                  \boxed{|\lambda(v_0)|\le Lc,\qquad
                         L=1+{\delta\over\rho}.}                \tag{2.5}
\]

If `[a,b]` lies in the left arm, (2.2) gives the sharper

\[
                         |\lambda-c|\le{\delta c\over\rho};     \tag{2.6}
\]

on the right arm it gives `|lambda+c|<=delta c/rho`.  These
same-arm estimates are not available at the central cusp, but the absolute
bound (2.5) is all that the orientation-free tube argument needs.

The constants are sharp at the scalar level.  On one arm, prescribe the
three values `P(a)=B(a)`, `P(v_0)=B(v_0)+e_0`, and `P(b)=B(b)` and join them
linearly.  The two inequalities in (2.2) become equalities and the slopes
jump downward at `v_0`.

## 3. Orientation-free small-excursion amplification

Fix the leaf `E=E_{v_0}` and write `P_0=P(v_0)`.  Because `Phi` is even,
the complement has the same anisotropic perimeter, the opposite Wulff
normal, and CMC multiplier `-lambda`.  Choose

\[
 A=E\quad\hbox{if }\lambda\le0,
 \qquad A=\Omega\setminus E\quad\hbox{if }\lambda>0.            \tag{3.1}
\]

Let

\[
 q_0=\mu(A),\qquad \ell=-|\lambda|\le0.                         \tag{3.2}
\]

Then `m(q_0)=m_0`, `P_Phi(A)=P_0`, and `ell` is the outward CMC
multiplier of `A`.

Follow each outward Wulff ray

\[
                         x+t z(N_x)                              \tag{3.3}
\]

until its first anisotropic cut, focal, or support-contact time `tau(x)`.
Assume in this section the smooth killed-tube coverage theorem: outside a
null set, the anisotropic parallel layer is covered exactly once by the
surviving regular rays.  Relative to `d sigma_Phi(x)dt`, its weighted
Jacobian is

\[
                         j_x(t)=e^{\ell t-D_x(t)},
 \qquad D_x(t)\ge0,quad D_x(t)\text{ nondecreasing}.             \tag{3.4}
\]

Define the killed normalized flux

\[
 R(t)=\int_{\partial A}1_{\{t<\tau(x)\}}
                         e^{-D_x(t)}d\sigma_\Phi(x).             \tag{3.5}
\]

It is nonincreasing and `R(0)=P_0`.  If `A_t` is the anisotropic outer
parallel set and `q(t)=mu(A_t)`, then for almost every `t`,

\[
                         q'(t)=e^{\ell t}R(t)
                              =P_\Phi(A_t).                       \tag{3.6}
\]

Choose `0<gamma<1`, let

\[
                         q_1=q_0+\gamma m_0,                     \tag{3.7}
\]

and denote by `T` the first time `q(T)=q_1`.  The target is always below
one because `1-q_0>=m_0`.

**Theorem 3.1 (small-excursion killed-tube amplification).**  Under
(0.1)--(0.2) and the smooth coverage hypotheses above,

\[
 \boxed{
 {\gamma\over(1+\delta)c}\le T
 \le{\gamma\over(1-\gamma)c}.}                              \tag{3.8}
\]

Moreover,

\[
 \boxed{
 {R(T-)\over P_0}\ge{1-\gamma\over1+\delta},\qquad
 {P_0-R(T-)\over P_0}
 \le\xi:={\delta+\gamma\over1+\delta}.}                    \tag{3.9}
\]

**Proof.**  The tent function `m` is one-Lipschitz.  Consequently, for
`q_0<=q<=q_1`,

\[
                         m(q)\ge m_0-(q-q_0)
                                  \ge(1-\gamma)m_0.              \tag{3.10}
\]

Equations (1.6) and (3.6) imply

\[
                         q'(t)\ge c(1-\gamma)m_0.                \tag{3.11}
\]

This proves the upper bound in (3.8).  Since `ell<=0` and
`R(t)<=P_0`, (3.6) also gives `q'(t)<=P_0`.  Thus

\[
 T\ge{q_1-q_0\over P_0}
   ={\gamma m_0\over P_0}
   \ge{\gamma\over(1+\delta)c},                                \tag{3.12}
\]

which is the lower bound.

Take regular times increasing to `T`.  From (1.6), (3.6), and (3.10),

\[
 R(T-)\ge e^{-\ell T}c(1-\gamma)m_0
        \ge c(1-\gamma)m_0.                                    \tag{3.13}
\]

Divide by `P_0<=(1+delta)cm_0` to obtain (3.9).  QED.

Several points are worth isolating.

1. The proof works whether `q_0` is below, above, or equal to `1/2`.
   If the excursion crosses `1/2`, (3.10) is still exact.
2. The orientation choice (3.1) is what makes `q'<=P_0`, hence the lower
   tube-length bound.  No guess about which side is the smaller side is
   needed.
3. A tiny fixed `gamma`, such as `10^-5`, makes the lost normalized flux
   tiny while retaining a tube of length at least a fixed multiple of
   `1/c`.
4. Global profile concavity never appears in the proof.  Concavity of the
   chosen CMC foliation is used only once, through (2.5).

## 4. Explicit long-ray covariance

Use the left-limit convention at a ray killed exactly at `T`.  The exact
flux identity is

\[
\begin{split}
 P_0-R(T-)={}&\sigma_\Phi\{\tau<T\}\\
 &+\int_{\{\tau\ge T\}}
       (1-e^{-D_x(T-)})d\sigma_\Phi(x).                       \tag{4.1}
\end{split}
\]

For `h>0`, put

\[
 G_h=\{x:\tau(x)\ge T,\ D_x(T-)\le h\}.                       \tag{4.2}
\]

Equations (3.9) and (4.1) give

\[
 \boxed{
 \sigma_\Phi((\partial A)\setminus G_h)
 \le\left(1+{1\over1-e^{-h}}\right)\xi P_0.}                 \tag{4.3}
\]

On each ray based in `G_h`, (2.5), (3.8), and (3.4) imply, for
`0<=t<=T`,

\[
 j_x(t)\ge
 \exp\left[-h-{L\gamma\over1-\gamma}\right].                  \tag{4.4}
\]

Indeed `ell>=-Lc`, `t<=T<=gamma/((1-gamma)c)`, and
`D_x(t)<=D_x(T-)<=h`.

**Corollary 4.1 (covariance amplification).**  Under the hypotheses of
Theorem 3.1,

\[
 \boxed{
 \operatorname {Cov}(\mu)\succeq
 {A(\rho,\delta,\gamma,h)\over c^3}
 \int_{G_h}z(N_x)\otimes z(N_x)d\sigma_\Phi(x),}               \tag{4.5}
\]

where

\[
 A(\rho,\delta,\gamma,h)
 ={\gamma^3\over12(1+\delta)^3}
   \exp\left[-h-
     {\left(1+\delta/\rho\right)\gamma\over1-\gamma}\right]. \tag{4.6}
\]

**Proof.**  For `theta in R^n`, restrict the variance integral to the
injective image of the good rays.  For every constant `a`,

\[
 \int_0^T
  (\langle\theta,x+t z(N_x)\rangle-a)^2dt
 \ge {T^3\over12}\langle\theta,z(N_x)\rangle^2.                \tag{4.7}
\]

Multiply by the lower density bound (4.4), integrate over `G_h`, and take
`a` to be the global mean of the linear functional.  Finally use the lower
bound for `T` in (3.8).  QED.

In the isotropic case, (4.5) is equivalently

\[
 \int_{G_h}z\otimes z\,d\sigma_\Phi
 \preceq A(\rho,\delta,\gamma,h)^{-1}c^3I.                      \tag{4.8}
\]

As in every killed-tube argument, (4.8) immediately excludes a coherent
normal packet when `c` is small.  It does not by itself exclude a
high-effective-rank packet; the product-exponential obstruction remains.

## 5. The equimeasurable construction supplies near-baseline leaves

This section distinguishes two deficits that must not be conflated.  Let
`p` be the Euclidean Cheeger constant and let `F:[0,1]` be the clipped heat
function.  Put

\[
 v(r)=\mu(F>r),\qquad
 L_0=\int_0^1m(v(r))dr,                               \tag{5.1}
\]

and define the **direct Euclidean Cheeger deficit**

\[
 D_E(F)=\int_0^1
   [P_E(\{F>r\})-p\,m(v(r))]dr
       =\operatorname {TV}_E(F)-pL_0.                         \tag{5.2}
\]

The fixed-scale heat estimate controls (5.2), not merely the smaller
deficit to the isoperimetric profile.  This extra fact is essential here.

Let `G` minimize the matrix-fidelity functional

\[
 \operatorname {TV}_E(G)
       +\kappa\|M(G)-M(F)\|_*,qquad 0<\kappa<1/3,              \tag{5.3}
\]

among functions equimeasurable with `F`.  Choose a compatible subgradient
`H`, with `||H||op<=1`, such that

\[
 \operatorname {tr}H(M(G)-M(F))
                  =\|M(G)-M(F)\|_*.                            \tag{5.4}
\]

Its constant surface tension is

\[
 \Phi_H(\xi)=|\xi|+
       \kappa{\xi^TH\xi\over|\xi|}.                            \tag{5.5}
\]

Comparison of (5.3) with `F` gives two inequalities.  The first one is the
important one for selecting levels:

\[
 \boxed{
 \operatorname {TV}_E(G)\le\operatorname {TV}_E(F),\qquad
 D_E(G)\le D_E(F).}                                            \tag{5.6}
\]

The second assertion follows from the first because `G` and `F` have the
same value law and hence the same `L_0`.  The compatible subgradient also
gives

\[
 \operatorname {TV}_{\Phi_H}(G)
 \le\operatorname {TV}_{\Phi_H}(F),                             \tag{5.7}
\]

but (5.7) is not needed for the Markov deletion.

For almost every level of `G`, let

\[
 d_E(r)=P_E(\{G>r\})-p\,m(v(r))\ge0.                            \tag{5.8}
\]

Fix `0<epsilon<1` and call a level bad if

\[
                         d_E(r)>\epsilon P_E(\{G>r\}).          \tag{5.9}
\]

Equations (5.6) and (5.8) give the exact surface-trace deletion

\[
 \boxed{
 \int_{\rm bad}P_E(\{G>r\})dr
 \le{D_E(F)\over\epsilon}.}                                  \tag{5.10}
\]

Every retained level satisfies

\[
                         P_E(\{G>r\})
 \le{p\,m(v(r))\over1-\epsilon}.                              \tag{5.11}
\]

Now, and only now, convert to the anisotropy.  Put

\[
                         c=(1-\kappa)p.                         \tag{5.12}
\]

Since `(1-kappa)P_E<=P_{Phi_H}<=(1+kappa)P_E`, the number `c`
is a valid global anisotropic Cheeger lower bound, and every retained leaf
obeys

\[
 \boxed{
 P_{\Phi_H}(\{G>r\})\le(1+\delta_\Phi)c\,m(v(r)),
 \qquad
 \delta_\Phi={1+\kappa\over(1-\kappa)(1-\epsilon)}-1.}        \tag{5.13}
\]

Thus exactly the Euclidean surface mass `D_E(F)/epsilon` is discarded,
while the anisotropy enters only through the pointwise relative error
`delta_Phi=epsilon+2kappa+O(epsilon^2+kappa^2)`.  This is stronger than
first forming an integrated anisotropic deficit, which would introduce an
unnecessary `O(kappa)/epsilon` loss.

### 5.1 The heat value law concentrates the retained levels near one half

Let `F_0=T_s1_S`, where `mu(S)=1/2`, before central clipping.  For every
`r in (0,1)`,

\[
 |\mu(F_0>r)-1/2|
 \le\mu(\{F_0>r\}\mathbin\triangle S).                         \tag{5.14}
\]

Layer cake therefore gives the exact estimate

\[
 \boxed{
 \int_0^1|\mu(F_0>r)-1/2|dr
 \le\|F_0-1_S\|_{L^1(\mu)}=2U(s).}                            \tag{5.15}
\]

Clipping only deletes level values, and `G` is equimeasurable with the
clipped function, so the same upper bound holds on its active value band.
For `h_0>0`, Markov's inequality gives

\[
 |\{r:|v(r)-1/2|>h_0\}|\le{2U(s)\over h_0}.                    \tag{5.16}
\]

Using `P_E=p m+d_E`, the Euclidean surface trace on these noncentral levels
is bounded by

\[
 \boxed{
 \int_{\{|v-1/2|>h_0\}}P_E(\{G>r\})dr
 \le {pU(s)\over h_0}+D_E(F).}                                \tag{5.17}
\]

This is the missing link between value measure and surface measure.  It is
valid because the direct Cheeger deficit, rather than merely deficit to the
isoperimetric profile, controls the error term.

Choose, for example, `h_0=10^-2`.  Every retained central volume then lies
in `[.49,.51]`.  If the chosen CMC foliation exists on `[1/4,3/4]`, then

\[
 [v-.48m(v),v+.48m(v)]\subset[1/4,3/4]                         \tag{5.18}
\]

for all such `v`.  Hence Theorem 3.1 applies with the uniform value
`rho=.48`, including at and across the half-volume cusp.  The trace deleted
by (5.17) is tiny at the fixed heat scale.

### 5.2 Euclidean normals versus Wulff displacements

For (5.5) and a unit vector `n`,

\[
 z(n)=D\Phi_H(n)
 =n+\kappa\{2Hn-(n^THn)n\}.                                   \tag{5.19}
\]

Consequently

\[
 |z-n|\le3\kappa,qquad
 1-\kappa\le\Phi_H(n)\le1+\kappa.                             \tag{5.20}
\]

For `kappa<=10^-2`, nuclear-norm comparison gives the explicit pointwise
bound

\[
 \|\Phi_H(n)z\otimes z-n\otimes n\|_*
 \le8\kappa.                                                   \tag{5.21}
\]

To verify it, use

\[
\begin{split}
 \|\Phi z\otimes z-n\otimes n\|_*
 &\le |\Phi-1||z|^2
      +|z-n|(|z|+1)\\
 &\le\kappa(1+3\kappa)^2+3\kappa(2+3\kappa)<8\kappa.
\end{split}
\]

Thus making `kappa` small controls the pointwise perimeter conversion
(5.13) and conversion of the retained Euclidean normal matrix to the
Wulff-displacement matrix used in (4.5).  No Loewner comparison of
individual rank-one matrices is asserted.

## 6. One explicit compatible hierarchy

The already-audited small-anisotropy scale is sufficient once selection is
performed in the Euclidean perimeter before conversion.  Freeze

\[
 \kappa=10^{-6},\qquad \alpha=10^{-28},\qquad
 \beta=10^{-14},\qquad \epsilon=10^{-5},
 \qquad h_0=10^{-2},\qquad h=1.                               \tag{6.1}
\]

The heat estimates give

\[
 {D_E(F)\over p}<6.014\cdot10^{-14},
 \qquad U(s)<1.254\cdot10^{-14}.                              \tag{6.2}
\]

Hence the relative-gap deletion (5.10) costs less than
`6.014*10^-9 p`, and the noncentral-volume deletion (5.17) costs less
than `1.315*10^-12 p`.  For every remaining leaf, (5.13) gives, with
`c=(1-10^-6)p`,

\[
 \delta_\Phi
 ={1+10^{-6}\over(1-10^{-6})(1-10^{-5})}-1
 =1.2000123\cdot10^{-5}.                                     \tag{6.3}
\]

All remaining volumes lie in `[.49,.51]`; (5.18) permits `rho=.48`, and

\[
                         L=1+{\delta_\Phi\over.48}
                         <1.000026.                             \tag{6.4}
\]

Increase the chosen oriented side by the fixed volume

\[
                         H_0=10^{-5}.                           \tag{6.5}
\]

This is Theorem 3.1 with `gamma=H_0/m(v)`, for which

\[
                         2\cdot10^{-5}\le\gamma
                         \le2.040817\cdot10^{-5}.               \tag{6.6}
\]

Uniformly over all retained levels,

\[
 {1.99997\cdot10^{-5}\over c}<T
 <{2.04086\cdot10^{-5}\over c}.                               \tag{6.7}
\]

The normalized killed-flux loss is below `3.24080*10^-5`, and after also
deleting rays with `D_x(T-)>1`, the deleted base surface fraction is below

\[
                         8.3677\cdot10^{-5}.                    \tag{6.8}
\]

The uniform covariance coefficient obtained from (4.6), using the lower
value `gamma=2*10^-5` in the cubic factor and the upper value in the
exponential factor, is

\[
                         A>2.4523\cdot10^{-16}.                 \tag{6.9}
\]

The previous matrix audit at this same `alpha,kappa` scale retained Wulff
projective angular variance greater than `.0031535` before the new level
and tube deletions.  The two level deletions above are negligible compared
with its trace.  Formula (7.13) of `wulff_tube_audit.md`, applied to (6.8),
costs less than `3.348*10^-4` of angular variance.  Thus a fixed positive
aggregate angular-variance packet still survives the tubes.  This does not
close its high-rank branch, but it verifies that the constant chain for the
present transfer is compatible.

Symbolically the new hierarchy is

\[
 {D_E\over\epsilon P_*}\ll\mathcal V,
 \qquad {U(s)\over h_0}\ll\mathcal V,
 \qquad \epsilon+\kappa+{H_0\over m_{\min}}\ll\mathcal V,
 \qquad {D_I\over\kappa tP_*}\ll\mathcal V.                  \tag{6.10}
\]

Here the last condition is only the old matrix-retention condition, with
`D_I<=D_E`; it is not used in selecting near-Cheeger leaves.  In
particular, `kappa=10^-6` passes all displayed requirements.

<!-- Superseded first-pass hierarchy retained only in source history.

\[
 M_\alpha=10.18038,qquad
 {3M_\alpha\over4(1+M_\alpha^2)}>0.072967,
\]

so the raw physical trace lower bound remains larger than
`0.002903 P_*` before the already-audited endpoint and alignment deletions.
The rank threshold becomes enormous but remains universal.  A final proof
would rerun all fixed-scale constants in one place; the calculation here
shows that there is no incompatible power law.

The symbolic hierarchy is clearer than the particular numbers.  If
`mathcal V` is the fixed angular-variance seed and `t` is the retained
trace fraction, it is enough to choose

\[
 \epsilon,\gamma\ll\mathcal V,qquad
 \kappa\ll\epsilon\mathcal V,qquad
 {D_E\over P_*}\ll
       \min\{\kappa t\mathcal V,\epsilon\mathcal V\}.          \tag{6.11}
\]

This first-pass hierarchy is superseded by the Euclidean-good-level
selection above; its `O(kappa)/epsilon` requirement was unnecessary.
-->

## 7. Central volumes and both orientations

The cusp of `m(v)` must be handled explicitly.

### 7.1 What fails if one insists on one arm

At `v_0=1/2`, a near-baseline leaf need not have multiplier close to `c` or
to `-c`.  A smooth symmetric log-concave approximation to the Laplace
density has a central half-line leaf with

\[
 P(1/2)=c/2,
 \qquad \lambda(1/2)=0.                                      \tag{7.1}
\]

For the exact Laplace density the profile is the tent itself.  Expanding a
half-line from mass `1/2` to `3/4` reduces the normalized flux by exactly a
factor two.  Thus the assertion that one central near-Cheeger leaf forces a
near-lossless tube over a fixed macroscopic volume change is false.

### 7.2 Why the small-excursion theorem survives

At the same Laplace central leaf, expanding only by
`gamma m_0=gamma/2` retains exactly the fraction `1-gamma` of normalized
flux.  This shows that (3.9) has the correct order and that the `gamma`
term cannot be removed.

For a general leaf, the two possible outward multipliers are `lambda` and
`-lambda`.  Choosing the nonpositive one is legitimate only because the
anisotropy is even.  The smaller-side mass during the chosen excursion may
increase, decrease, or pass through `1/2`; the single Lipschitz estimate
(3.10) covers all three cases.  The surrounding secants cross the cusp
without difficulty because they are used only to prove the absolute bound
(2.5).

Thus central levels do not need to be discarded.  They do, however, need
the same fixed surrounding **foliation** interval as every other retained
level.  Direct coarea deficit alone does not provide that interval.

## 8. Model and failure-mode audit

1. **Symmetric Laplace and smooth approximations.**  These refute
   same-arm multiplier control at `1/2` and macroscopic near-lossless
   central amplification.  They saturate the order `gamma` in (3.9), but
   satisfy Theorem 3.1.

2. **One-sided exponential.**  A small tail has outward multiplier `c`.
   Rule (3.1) chooses the complement, whose multiplier is `-c`.  The tube
   moves into the exponential tail and has the required `Theta(1/c)`
   length.  This verifies the complement orientation and the density
   factor in (4.4).

3. **Gaussian.**  The median halfspace has multiplier zero and exact
   near-baseline perimeter.  A small central excursion loses less than the
   general `gamma` bound.  No curvature or log-Sobolev estimate is being
   smuggled into the proof.

4. **Cube.**  For a coordinate cut, `lambda=0`, the Wulff rays are straight,
   and normalized flux stays constant until support contact.  Equations
   (3.8)--(4.5) remain valid with the contact recorded as killing.

5. **Ball caps.**  Moving relative-isoperimetric caps can rotate through
   support contact even with zero deficit to the global profile.  Theorem
   3.1 does not claim otherwise.  It says that a near-Cheeger cap can lose
   only `O(delta+gamma)` normalized flux during the small excursion; the
   contact charge is included in (4.1).

6. **Sharp scalar concave model.**  The piecewise-linear construction after
   Lemma 2.1 shows that the term `delta/rho` cannot be improved using only
   the hypotheses stated here.

7. **No surrounding interval.**  Let the left and right secant distances
   tend to zero in that scalar model.  The multiplier becomes arbitrarily
   large while the pointwise relative deficit stays fixed.  A single
   isolated CMC equation therefore cannot replace (0.2).

8. **Non-even anisotropy.**  If `Phi(-n) != Phi(n)`, passage to the
   complement need not preserve perimeter and the two multipliers need not
   be opposites.  The orientation trick can then fail.  The matrix
   anisotropy (5.5) is even, so this is not an obstruction to the present
   route.

9. **High-rank product geometry.**  Corollary 4.1 bounds each surviving
   displacement matrix by `C c^3 I` in the isotropic case.  It does not
   turn a matrix with trace of order `c` and effective rank of order
   `1/c^2` into a coherent direction.  This is the already identified
   product/high-rank escape, not a flaw in the amplification theorem.

## 9. The exact remaining singular theorem

The smooth argument reduces the transfer problem to the following
formalizable statement.

> **Singular nested-CMC concavity and Wulff-coverage theorem.**  Let `G` be
> an equimeasurable matrix-fidelity minimizer and let `Phi_H` be one
> compatible constant anisotropy.  Outside a set of levels whose total
> `Phi_H`-perimeter is controlled by the direct deficit, every retained
> level `r_0` must have all of the following properties:
>
> 1. on the volume interval
>    `[v_0-rho m(v_0),v_0+rho m(v_0)]`, the neighboring reduced
>    boundaries form one nested generalized positive-lapse foliation;
> 2. their perimeter has a finite concave representative, or equivalently
>    the CMC multiplier has a distributional derivative which is a
>    nonpositive Radon measure including all critical, topology-change,
>    singular, and hard-contact charges;
> 3. the retained leaf has an exact killed Wulff area formula: almost every
>    point in its outer layer with unique anisotropic projection comes from
>    the regular reduced boundary, and first support contact, focal loss,
>    and cut collision contribute only through the killed-flux term;
> 4. these conclusions are stable under the smooth convex approximation
>    used to treat extended-valued potentials and nonsmooth supports, with
>    the surface and displacement matrices converging with the same
>    constants.

The formal smooth Euler equation supplies constant anisotropic mean
curvature on regular levels, and the Jacobi identity supplies concavity
between critical values.  Neither proves the theorem above across the
exceptional set.  In particular, a `BV` singular stratum can in principle
emit positive Wulff tube volume even if it has zero perimeter, and a
topology change can prevent the two secants in Lemma 2.1 from belonging to
one positive-lapse component.

This is narrower than concavity of the global anisotropic isoperimetric
profile: no leaf is required to be globally minimizing.  It is nevertheless
global in the volume parameter around each retained leaf, and it remains
load bearing.  Once it is proved, the direct-deficit selection
(5.6)--(5.13), central concentration (5.15)--(5.18), multiplier bound
(2.5), and covariance theorem (4.5) compose without any profile-minimizer
replacement.

## 10. A possible further simplification: minimize the direct deficit

There is a promising way to remove even the surrounding-foliation secant
hypothesis.  Instead of imposing equimeasurability, minimize over
`0<=G<=1`

\[
 \mathcal K(G)=
 \underbrace{\left\{\operatorname {TV}_E(G)
 -p\int_0^1m(\mu(G>r))dr\right\}}_{D_p(G)\ge0}
 +\kappa\|M(G)-M(F)\|_*.                                    \tag{10.1}
\]

Comparison with `F` gives immediately

\[
 D_p(G)\le D_p(F),\qquad
 \|M(G)-M(F)\|_*\le{D_p(F)\over\kappa}.                       \tag{10.2}
\]

There is no prescribed value law, but none is needed for these two
conclusions.  Existence looks sound: `D_p(G)>=0`, the range constraint and
`TV(G)<=D_p(F)+p/2` give `BV` compactness, the distribution term is
continuous under `L^1` convergence by convergence in measure and dominated
convergence in the level variable, and the matrix term is lower
semicontinuous as before.

For a smooth minimizer, let `v(r)=mu(G>r)`.  The first variation of the
distribution term is

\[
 {d\over dt}\bigg|_{t=0}
 \int_0^1m(\mu(G+th>r))dr
 =\int h(x)m'(v(G(x)))d\mu(x).                               \tag{10.3}
\]

Consequently a compatible matrix subgradient gives the prescribed-CMC
equation

\[
 -\operatorname {div}_\mu D\Phi_H(\nabla G)
 =p\,m'(v(G))\in[-p,p].                                      \tag{10.4}
\]

Away from half volume the multiplier is exactly `p` or `-p`; at the kink
it belongs to the displayed interval.  For an `SBV` jump of height
`b-a`, displacement of the jump interface averages the right side over
`r in (a,b)`, so its divided-difference multiplier also lies in
`[-p,p]`.  Thus (10.4), if established for the actual minimizer, supplies
the multiplier bound needed in Section 4 directly and makes Lemma 2.1
unnecessary.

This is a sanity-checked candidate, not yet a completed replacement.  The
load-bearing analytic issue becomes a nonsmooth Euler theorem for (10.1):
one must justify (10.3) at atoms and critical values, prove a single
compatible `H`, describe the Cantor part of `DG`, and show that almost every
surface packet retained by (10.2) has the singular killed-Wulff coverage
needed for (4.5).  The proposal appears to remove the global concavity
problem, but it does not remove the singular coverage problem or the final
high-rank inverse.
