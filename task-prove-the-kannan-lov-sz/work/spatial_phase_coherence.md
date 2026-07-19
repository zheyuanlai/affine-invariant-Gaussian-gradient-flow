# Spatial phase coherence across Gaussian natural tilts

## 1. Purpose and conclusion

Fix `t>0` and a common convex potential `W`.  Consider the natural-tilt
family

\[
 d\pi_c(x)={1\over Z(c)}
 \,\exp\{c\cdot x-t|x|^2/2-W(x)\}\,dx,              \tag{1.1}
\]

with the usual approximation interpretation for an extended-valued convex
`W`.  Let `S` be one fixed Borel set.  At a tilt `c`, put

\[
 g_c=\pi_c(S),\qquad
 v_c=\mathbb E_c[(\mathbf1_S-g_c)(X-\mathbb E_cX)],
 \qquad u_c={v_c\over|v_c|},                         \tag{1.2}
\]

and define the relative sharp-centroid defect

\[
 \epsilon_c=1-{sqrt t|v_c|\over\mathcal I(g_c)}.    \tag{1.3}
\]

This note proves three dimension-free facts.

1. Near equality makes `S` close in probability to its active halfspace,
   with error `O_delta(sqrt(epsilon_c))`.
2. Near equality forces a dimension-free approximate Gaussian product
   splitting of the Brenier map in the active direction.  It also gives a
   weaker common bulk-resolvent certificate for the excess potential.
3. At two tilts separated by `O(sqrt(t))`, either the active lines cohere or
   the midpoint posterior has a bulk-resolvent-flat two-dimensional plane;
   each endpoint separately has an approximate Gaussian factor.  Exact
   equality at two tilts forces the same active line.

These statements are local in natural-tilt distance.  A global claim that
all near-equality phases have one direction is false: a Gaussian parity set
has asymptotically saturated orthogonal phases at tilts whose separation
tends to infinity.  Turning the local dichotomy into a single global phase
would require a dimension-free bridge estimate for a log-concave tilt law;
without extra structure that bridge estimate is a KLS-strength input.

## 2. From centroid defect to total threshold error

We first strengthen the weighted threshold error to an unweighted one.

**Lemma 2.1 (total threshold approximation).**  Fix
`delta in (0,1/2)`.  Let `pi` be `t`-strongly log-concave, let
`g=pi(S) in [delta,1-delta]`, and use the notation in (1.2)--(1.3).  Let

\[
 Y=\langle X-\mathbb EX,u\rangle,
 \qquad H=\mathbf1_{\{Y\ge a\}},\qquad \pi(H)=g.      \tag{2.1}
\]

Then

\[
 \boxed{\pi(S\mathbin\triangle H)
       \le C_\delta\sqrt{\epsilon}.}                \tag{2.2}
\]

**Proof.**  Let `rho` be the law of `Y`, and let
`q(y)=pi(S|Y=y)`.  The exact centroid decomposition gives

\[
 \Delta={\mathcal I(g)\over\sqrt t}-|v|
        =D_{cut}+D_{map},                             \tag{2.3}
\]

where

\[
 D_{cut}=\int(y-a)(H-q)d\rho,
 \qquad
 D_{map}\ge0.                                        \tag{2.4}
\]

Set `d=sqrt(t)D_map`.  Let `T` be the increasing contraction from a
standard Gaussian to `sqrt(t)Y`, put `z_0=Phi^{-1}(1-g)`, and write
`R(z)=T(z)-z`.  Thus `R` is nonincreasing and

\[
 d=g(1-g)\mathbb E[R(Z_-)-R(Z_+)],                   \tag{2.5}
\]

where `Z_-` and `Z_+` are Gaussian conditioned below and above `z_0`.

For `b>0`, let

\[
 p=\rho\{|Y-a|\le b/\sqrt t\}.                      \tag{2.6}
\]

Suppose first that at least `p/2` of this mass is above `a`.  Split the
corresponding Gaussian quantile interval into two pieces of equal Gaussian
mass.  The first piece has Euclidean length at least `c p`, because the
Gaussian density is bounded above.  Its image under `T` has length at most
`b`.  Hence `R` drops by at least `cp-b` between its endpoints.  In (2.5),
pair the farther half of the upper interval, of unconditional mass at least
`p/4`, with the entire lower central side, whose mass is at least `delta`.
Monotonicity of `R` gives

\[
 d\ge c_\delta p(cp-b).                              \tag{2.7}
\]

If the lower part of the strip has mass at least `p/2`, the same argument is
applied with the two sides interchanged.  In either case, (2.7), with the
trivial alternative `cp<=2b`, implies

\[
 p\le C_\delta(b+\sqrt d).                            \tag{2.8}
\]

Outside the strip, (2.4) gives

\[
 \int_{|y-a|>b/\sqrt t}|H-q|d\rho
 \le {\sqrt tD_{cut}\over b}.                        \tag{2.9}
\]

Since `sqrt(t)Delta=I(g)epsilon<=C epsilon`, combine (2.8)--(2.9) and
choose `b=sqrt(epsilon)`.  This proves (2.2).  All quantile statements can
be made with generalized inverses; strong log-concavity in fact gives a
continuous one-dimensional density.

## 3. A one-quantile defect still forces almost maximal variance

The logarithmic modulus in the next lemma is necessary for this argument:
a one-dimensional contraction may differ from the identity only in a far
Gaussian tail.

**Lemma 3.1 (variance from one central centroid).**  Under the hypotheses of
Lemma 2.1,

\[
 \boxed{
 1-t\operatorname {Var}_\pi\langle X,u\rangle
 \le C_\delta\epsilon
       \sqrt{\log(e/\epsilon)}.}                     \tag{3.1}
\]

The right side is interpreted as zero at `epsilon=0` and enlarged to a
constant when `epsilon>=1/2`.

**Proof.**  Again scale to `t=1` and write the monotone contraction as
`T(Z)=Z+R(Z)`, with `ER=0`, `-1<=R'<=0`.  Put `q=-R'`.  The threshold flux

\[
 k_g(z)=\begin{cases}
 g\Phi(z)/\varphi(z),&z<z_0,\\
 (1-g)(1-\Phi(z))/\varphi(z),&z\ge z_0
 \end{cases}                                         \tag{3.2}
\]

satisfies

\[
 D_{map}=\mathbb E[k_g(Z)q(Z)].                       \tag{3.3}
\]

For central `g`, Mills' bounds give

\[
 k_g(z)\ge {c_\delta\over1+|z|}.                     \tag{3.4}
\]

For every `L>=1`, split at `|Z|=L` and use `0<=q<=1`:

\[
 \mathbb E q
 \le C_\delta(1+L)D_{map}+\mathbb P\{|Z|>L\}.       \tag{3.5}
\]

Choose `L=sqrt(2log(e/D_map))`, with the evident convention if the defect is
not small.  Then

\[
 \mathbb E q\le C_\delta D_{map}
                    \sqrt{\log(e/D_{map})}.          \tag{3.6}
\]

Finally Gaussian integration by parts gives

\[
 \begin{aligned}
 1-\operatorname {Var}(T(Z))
 &=-2\mathbb E[ZR]-\mathbb E R^2\\
 &=2\mathbb E q-\mathbb E R^2\le2\mathbb E q.       \tag{3.7}
 \end{aligned}
\]

Since `D_map<=Delta` and `sqrt(t)Delta=I(g)epsilon`, scaling proves (3.1).

## 4. Approximate Gaussian splitting and common bulk flatness

### 4.1 A genuine Brenier product certificate

The cleanest certificate uses the Brenier map and includes nonsmooth support
boundaries.  Let `T_c` be the Brenier contraction from
`N(0,t^{-1}I)` to `pi_c`, and define the standardized centered map

\[
 \widetilde T_c(G)=\sqrt t\{T_c(G/\sqrt t)-\mathbb E_cX\},
 \qquad G\sim N(0,I).                                  \tag{4.0a}
\]

At differentiability points its Jacobian `H_c` is symmetric and satisfies
`0 preceq H_c preceq I`.  Lemma 3.1 and Gaussian Poincare give

\[
 1-\zeta_\delta(\epsilon_c)
 \le\operatorname {Var}\langle u_c,\widetilde T_c(G)\rangle
 \le\mathbb E|H_c(G)u_c|^2\le1,                       \tag{4.0b}
\]

where `zeta_delta(r)=C_delta r sqrt(log(e/r))`.  For a positive
contraction,

\[
 |Hu-u|^2\le1-|Hu|^2.                                  \tag{4.0c}
\]

Consequently

\[
 \boxed{\mathbb E|H_cu_c-u_c|^2
       \le\zeta_\delta(\epsilon_c).}                  \tag{4.0d}
\]

This derivative statement integrates to an approximate product coupling.
Write `G=su_c+Z`, let `P_c=I-u_cu_c^T`, and put

\[
 \overline R_c(Z)=\mathbb E_s[P_c\widetilde T_c(su_c+Z)].          \tag{4.0e}
\]

Gaussian Poincare, first globally for the active component and then on each
one-dimensional `s`-fiber for the transverse component, gives

\[
 \boxed{
 \mathbb E\left|\widetilde T_c(G)-
       \{u_c\langle G,u_c\rangle+\overline R_c(P_cG)\}
       \right|^2
 \le2\zeta_\delta(\epsilon_c).}                       \tag{4.0f}
\]

Thus every near-equality phase has an actual approximate Gaussian factor,
not merely small bulk Hessian.  This remains meaningful for a polytope or a
one-sided exponential: boundary/contact contraction is encoded in the
global Brenier map.

### 4.2 A common bulk-potential certificate

Assume first that `W` is finite and smooth.  Put

\[
 H_W(x)=\nabla^2W(x),
 \qquad
 K_t(x)=t(tI+H_W(x))^{-1},                            \tag{4.1}
\]

so `0 preceq K_t preceq I`.  Brascamp--Lieb gives

\[
 \operatorname {Var}_c\langle X,u_c\rangle
 \le {1\over t}\mathbb E_c[u_c^TK_t(X)u_c].         \tag{4.2}
\]

Combining with Lemma 3.1 yields the resolvent-flatness certificate

\[
 \boxed{
 \mathbb E_c[u_c^T(I-K_t(X))u_c]
 \le \zeta_\delta(\epsilon_c),
 \quad
 \zeta_\delta(r)=C_\delta r\sqrt{\log(e/r)}.}       \tag{4.3}
\]

The matrix in the integrand is

\[
 I-K_t=H_W(tI+H_W)^{-1}.                              \tag{4.4}
\]

It vanishes exactly on bulk-flat directions of the smooth common potential.
For an extended-valued `W`, this certificate must be interpreted only after
a specified smooth approximation.  It can miss contact with a hard facet:
the interior Hessian of a polytope is zero even though its boundary changes
the law.  The Brenier certificate (4.0d)--(4.0f), not (4.3), is the genuine
approximate-product statement in that case.

We now transfer two such certificates to one common posterior.  Let
`c_1,c_2` be two tilts, put `c_m=(c_1+c_2)/2`, and assume

\[
 |c_1-c_2|\le R\sqrt t.                              \tag{4.5}
\]

The log-partition Hessian satisfies
`nabla^2 log Z(c)=Cov_c(X) preceq t^{-1}I`.  Therefore

\[
 \mathbb E_{c_i}
 \left({d\pi_{c_m}\over d\pi_{c_i}}\right)^2
 \le \exp(R^2/4).                                    \tag{4.6}
\]

Apply Cauchy--Schwarz to the `[0,1]`-valued integrand in (4.3).  If both
endpoint defects are at most `epsilon`, then, with

\[
 \kappa=C_R\sqrt{\zeta_\delta(\epsilon)},
 \qquad B=\mathbb E_{c_m}K_t(X),                     \tag{4.7}
\]

one has

\[
 \boxed{u_i^T(I-B)u_i\le\kappa,qquad i=1,2.}        \tag{4.8}
\]

This gives the promised dichotomy.  If `theta` is the unoriented angle
between `u_1,u_2`, then either `sin(theta)<=alpha`, or, on their span `E`,

\[
 \boxed{
 \operatorname {Tr}_E(I-B)
 \le {5\kappa\over\alpha^2}.}                        \tag{4.9}
\]

Indeed, take `e_1=u_1` and
`e_2=(u_2-<u_1,u_2>u_1)/sin(theta)`.  The seminorm
`x mapsto |(I-B)^{1/2}x|` and (4.8) give

\[
 e_2^T(I-B)e_2\le {4\kappa\over\sin^2(theta)},
\]

which proves (4.9).  Thus separated direction phases force a common
bulk-resolvent-flat two-plane, while (4.0f) supplies a genuine approximate
Gaussian factor at each endpoint.  More phases give the analogous bulk
statement on every well-conditioned subspace of their span; the loss is the
inverse square of the smallest singular value of the direction frame, which
is unavoidable.  Synchronizing the separate endpoint Brenier couplings into
one common multi-factor coupling is an additional problem; common-source
optimal maps do not provide that synchronization automatically.

## 5. Overlap also glues the threshold lines

The profile coordinate

\[
 z(c)=\Phi^{-1}(g_c)
\]

is `t^{-1/2}`-Lipschitz because

\[
 \nabla_cz={v_c\over\mathcal I(g_c)},qquad
 |v_c|\le {\mathcal I(g_c)\over\sqrt t}.             \tag{5.1}
\]

Hence, under (4.5), centrality at the endpoints implies
`g_{c_m} in [delta_R,1-delta_R]` for a positive `delta_R`.

Let `H_i` be the active halfspace from Lemma 2.1 at `c_i`.  Equations
(2.2) and (4.6) give

\[
 \pi_{c_m}(S\mathbin\triangle H_i)
 \le C_{\delta,R}\epsilon_i^{1/4}.                  \tag{5.2}
\]

Consequently

\[
 \boxed{
 \pi_{c_m}(H_1\mathbin\triangle H_2)
 \le C_{\delta,R}(\epsilon_1^{1/4}+epsilon_2^{1/4}).} \tag{5.3}
\]

There is a useful purely two-dimensional interpretation.  Let `A_m` be the
covariance of `pi_{c_m}` and define the covariance-normalized covectors

\[
 n_i={A_m^{1/2}u_i\over|A_m^{1/2}u_i|}.              \tag{5.4}
\]

Whenever the denominators are nonzero, fixed-dimensional log-concave
halfspace geometry gives

\[
 \boxed{
 |n_1-n_2|
 \le C_{\delta,R}
       (\epsilon_1^{1/4}+\epsilon_2^{1/4}).}          \tag{5.5}
\]

For completeness, project `pi_{c_m}` onto the span of the two linear forms
and whiten that marginal.  It is an isotropic log-concave law in dimension
at most two.  In fixed dimension, its density is bounded below on a ball
whose radius and lower bound depend only on the centrality parameter, and
all central one-dimensional quantiles lie in a ball of the same type.  The
Euclidean area inside that ball of the symmetric difference of two central
halfplanes is at least a constant times the angle between their normals.
Multiplying by the density lower bound proves (5.5).  This proof is entirely
two-dimensional and introduces no ambient-dimensional constant.

Equation (5.5) is the correct unconditional line statement.  Euclidean
normals cannot in general be recovered from covariance-normalized normals
without a lower covariance bound on their span.  The separate resolvent
dichotomy (4.9) supplies only a common bulk-flatness certificate when
Euclidean phases remain separated; the true endpoint product certificates
are (4.0d)--(4.0f).

## 6. Exact equality rigidity

If `epsilon_c=0`, the active marginal is exactly
`N(m,t^{-1})`, and `S` is, modulo `pi_c`-null sets, its active halfspace.
Equality in the one-dimensional Lichnerowicz inequality for a
`t`-strongly log-concave law forces the full measure to split:

\[
 \pi_c=N(m,t^{-1})\otimes\nu
\]

in the active/transverse decomposition.  Equivalently, `W` is affine in the
active direction and that affine part is absorbed into the natural tilt.

All members of (1.1) are mutually absolutely continuous on the relative
interior of their common convex support.  Therefore, if exact equality holds
at two tilts, the same fixed set `S` agrees almost everywhere with two affine
halfspaces.  Two nontrivial affine halfspaces which agree on an open convex
set have the same oriented normal.  Hence

\[
 \boxed{\epsilon_{c_1}=\epsilon_{c_2}=0
        \quad\Longrightarrow\quad u_{c_1}=u_{c_2}.}  \tag{6.1}
\]

The approximate two-plane alternative in (4.9) is thus a quantitative
certificate, not a new exact equality case for the fixed set.

## 7. Model audit and the global obstruction

### 7.1 Gaussian parity: a sharp no-go to unrestricted global coherence

For `X~N((a,b),I)` and

\[
 S=\{x_1x_2\ge0\},
\]

put `p=Phi(a)`, `q=Phi(b)`.  Along `b=0`,

\[
 g={1\over2},\qquad
 v=\varphi(0)(2\Phi(a)-1)e_2,
\]

and therefore

\[
 \mathcal I(1/2)-|v|=2\varphi(0)\Phi(-|a|).          \tag{7.1}
\]

Along the orthogonal arm `(0,a)`, the active direction is `e_1` with the
same defect.  Thus the two phases become arbitrarily close to equality while
their directions stay orthogonal.  Their natural-tilt distance is
`sqrt(2)|a|`, which tends to infinity, and every path between them crosses a
region with non-small defect (indeed `v=0` at the origin).  Here `W=0`, so
the alternative (4.9) correctly reports a full Gaussian product factor.
This example rules out a single direction valid at all near-equality tilts.

### 7.2 Radial phases

For a radial cut under a Gaussian posterior, `u(c)=c/|c|`.  Near equality
can occur only on far tilt spheres.  Two points on a sphere of radius `R`
whose tilt distance is `O(sqrt(t))` have angular separation
`O(sqrt(t)/R)`, exactly as (5.5) predicts.  Globally the sphere contains a
continuum of phases, but it surrounds the zero-centroid state at the origin;
again the target is a full Gaussian factor and there is no global line.

### 7.3 Product exponentials and polyhedral cells

For the balanced maximum cut of `n` independent shifted exponentials at
curvature zero,

\[
 |v|=\Theta((\log n)/\sqrt n).
\]

At a small positive localization curvature `t`, its normalized centroid
ratio is `O(sqrt(t)log(n)/sqrt(n))`, so the coordinate phases are far from
the hypothesis of this note.  A large tilt which isolates one coordinate
can create a locally Gaussian phase only after the posterior core is far
from the one-sided boundary; in that regime the excess potential is locally
flat in precisely the selected direction.  Different coordinate cells are
separated by tie regions where the defect is not small.  The same description
applies to cube, simplex, and general polyhedral facet phases: local
saturation is a Gaussian-core phenomenon, while changes of facet cross a
non-saturated contact region.

### 7.4 Why the remaining global bridge is KLS-strength

At fixed localization time, the natural parameter has density

\[
 \nu_t(dc)\propto e^{-|c|^2/(2t)}Z(c)\,dc,
 \qquad
 \nabla^2\log\nu_t(c)=-t^{-1}I+\operatorname {Cov}_{\pi_c}(X)
 \preceq0.                                             \tag{7.2}
\]

Thus the tilt law is log-concave, but not uniformly strongly log-concave.
The local overlap lemmas above say that each connected `O(sqrt(t))` cluster
of low-defect states has a coherent covariance line, unless it carries an
approximate Gaussian factor.  To show that two positive-mass clusters cannot
be separated by a small bridge of high-defect states would require a
dimension-free expansion estimate for the log-concave law `nu_t`.  Such an
estimate is a Cheeger/KLS statement for `nu_t` itself.  Log-concavity alone
does not supply it, and inserting it would be circular.

Therefore the exact output of the Brenier/overlap analysis is the local
dichotomy (4.9) plus covariance-line coherence (5.5).  A noncircular global
closure must add a new phase-count, entropy, or bridge mechanism; ordinary
overlap and equality rigidity do not by themselves provide one.
