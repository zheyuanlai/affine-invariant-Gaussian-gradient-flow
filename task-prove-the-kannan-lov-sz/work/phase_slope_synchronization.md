# Phase-slope synchronization in the log-affine core model

## 0. Purpose and scope

The exact heat-generated model in `heat_selector_variation.md` has phase
weights

\[
 \pi_i(z)={A_i e^{-c_i^2+\sqrt2c_i z}\over
                 \sum_jA_j e^{-c_j^2+\sqrt2c_j z}},
 \qquad c_i=\sqrt s\,\kappa_i,                         \tag{0.1}
\]

and phase projectors `Z_i=n_i n_i^T`.  This note proves two finite-dimensional
facts which isolate the exact geometric input still required.

1. A fixed between-level projector variance forces a fixed variance of the
   dimensionless normal slopes `c_i`.
2. If facet-height variations have cutoff cost at most `C/s`, that slope
   variance yields a fixed physical perimeter saving.  Equal slopes are the
   only zero-charge branch; in that branch the softmax proportions are
   independent of the level.

The second statement is conditional on an explicit capacity hypothesis.  It
does not assert that arbitrary physical phase packets have that capacity.
The product-exponential model has equal slopes and correctly falls into the
zero-charge branch.

## 1. Softmax rotation is paid by slope variance

Let `H` be a finite-dimensional Hilbert space, let `Z_1,...,Z_m in H` satisfy
`|Z_i|<=1`, and let `A_i>0`, `c_i in R`.  Define `pi_i(z)` by (0.1) and

\[
 Q(z)=\sum_i\pi_i(z)Z_i,
 \qquad \bar c(z)=\sum_i\pi_i(z)c_i.                    \tag{1.1}
\]

Let `rho` be a probability on an interval and suppose its Poincare constant
for the ordinary derivative is at most `C_z`.

**Lemma 1.1 (softmax slope charge).**  One has

\[
 \boxed{
 \operatorname {Var}_\rho Q
 \le2C_z\int\operatorname {Var}_{\pi(z)}(c)\,d\rho(z).}
                                                               \tag{1.2}
\]

Here `Var_rho Q=int|Q-E_rho Q|^2d rho` and
`Var_pi(c)=sum_i pi_i(c_i-bar c)^2`.

**Proof.**  Direct differentiation gives

\[
 \pi_i'(z)=\sqrt2\,\pi_i(z)(c_i-\bar c(z))
\]

and therefore

\[
 Q'(z)=\sqrt2\sum_i\pi_i(z)(c_i-\bar c(z))(Z_i-Q(z)). \tag{1.3}
\]

Hilbert-space Cauchy--Schwarz yields

\[
 |Q'(z)|^2
 \le2\operatorname {Var}_{\pi(z)}(c)
       \sum_i\pi_i(z)|Z_i-Q(z)|^2.                    \tag{1.4}
\]

The last sum equals
`sum_i pi_i|Z_i|^2-|Q|^2` and is at most one.  Apply the
Poincare inequality for `rho`.  This proves (1.2).  QED.

For normal projectors `Z_i=n_i n_i^T`, the Hilbert norm is the
Hilbert--Schmidt norm and `|Z_i|=1`.  Thus a between-level variance `B`
implies

\[
 \int\operatorname {Var}_{\pi(z)}(c)d\rho(z)
 \ge {B\over2C_z}.                                      \tag{1.5}
\]

If `rho` is supported on an interval of length `L_z`, the one-dimensional
Neumann inequality gives `C_z<=L_z^2/pi^2`.  No dimension enters.

**Corollary 1.2 (exact zero branch).**  If the right side of (1.2) vanishes,
then all slopes with positive phase weight agree for `rho`-almost every
connected phase cluster.  If all `c_i` agree, every `pi_i(z)` and hence
`Q(z)` is constant in `z`.

The last assertion is immediate from (0.1).  Notice that equal slopes do
not force the normals to agree; the product-exponential facets have many
orthogonal normals and the common slope one.

## 2. A capacity-to-saving lemma

Fix one physical level and let its disjoint log-affine facet cores have
weighted areas `a_i>0`, normal slopes `kappa_i`, and total area
`a=sum_i a_i`.  Put

\[
 \bar\kappa={1\over a}\sum_i a_i\kappa_i.             \tag{2.1}
\]

The exact graph identity in `physical_translation_splicing.md` shows that a
height `h_i` on the `i`-th core has first-order term
`-a_i kappa_i h_i`.  The following proposition records precisely the
cutoff estimate which would make that term useful.

**Proposition 2.1 (slope synchronization from splice capacity).**  Suppose
that for every real vector `h=(h_i)` with

\[
                         \sum_i a_i h_i=0              \tag{2.2}
\]

and `sum_i a_i h_i^2` sufficiently small, there is an exactly
mass-preserving finite physical splice satisfying

\[
 P_\mu(B_h)-P_\mu(A)
 \le-\sum_i a_i\kappa_i h_i
       +{C_{cap}\over2s}\sum_i a_i h_i^2.             \tag{2.3}
\]

Then the supremal finite-splice saving obeys

\[
 \boxed{
 \mathfrak G_{fin}(A)
 \ge {s\over2C_{cap}}
       \sum_i a_i(\kappa_i-\bar\kappa)^2}             \tag{2.4}
\]

whenever the optimizer below lies in the stated small-height range.  If
only a factor `theta in(0,1]` of that optimizer is admissible, the right
side is replaced by
`theta(2-theta)s/(2C_cap)` times the same variance.

**Proof.**  Set

\[
 h_i={s\over C_{cap}}(\kappa_i-\bar\kappa).            \tag{2.5}
\]

It satisfies (2.2).  Substitution into (2.3) gives

\[
 P_\mu(B_h)-P_\mu(A)
 \le-{s\over2C_{cap}}
       \sum_i a_i(\kappa_i-\bar\kappa)^2.             \tag{2.6}
\]

Replacing `h` by `theta h` gives the final assertion.  QED.

Since `c_i=sqrt(s)kappa_i`, formula (2.4) is

\[
 \mathfrak G_{fin}(A)
 \ge {1\over2C_{cap}}\sum_i a_i(c_i-\bar c)^2.         \tag{2.7}
\]

Thus, in the exact log-affine softmax model, Lemma 1.1 and Proposition 2.1
convert a fixed between-level variance into a fixed fraction of physical
facet area, provided the facet-height capacity is `O(1/s)`.

## 3. Constants at the physical budget

Suppose a family of levels carries total selected area `A>=.004p`, its
between-level projector variance is at least `8/17`, the profile-coordinate
law has `C_z<=4`, and (2.3) holds after integration over the levels with
`C_cap<=10`.  Then (1.5) and (2.7) give the model saving

\[
 {A\over 2C_{cap}}{8/17\over 2C_z}
 \ge {.004p\over20}{1\over17}
 >1.17\,10^{-5}p.                                    \tag{3.1}
\]

This does not yet exceed the audited `6.02*10^-5p` budget.  The constants
therefore matter: one needs either `C_z C_cap<1.95`, more of the analytic
trace `.005109p`, a direct total-variation version of Lemma 1.1, or a
simultaneous use of the same-level charge.  The purpose of (3.1) is to make
the numerical requirement explicit; no favorable constant is being hidden.

The ideal full-core log-affine variation has no taper cost, corresponding
to `C_cap` arbitrarily small over a sufficiently small admissible height.
At the opposite extreme, a core of tangential radius much smaller than
`sqrt(s)` has cutoff cost much larger than `1/s`.  Such a core lies within
Gaussian scale of a ridge or end.  Proving that the resulting incidence
capacity supplies the missing saving is exactly the unresolved geometric
alternative.

### 3.1 Separated log-affine components are charged without a capacity
inequality

There is a stronger estimate when low incidence has already produced actual
set components.  Let `A` be the disjoint union of finite-perimeter sets
`A_i`, assume `mu(A)<=1/2` and `mu(A_i)<=1/2`, and suppose perimeter is
additive:

\[
 v_i=\mu(A_i),\qquad a_i=P_\mu(A_i),\qquad
 P_\mu(A)=\sum_i a_i.                                  \tag{3.2}
\]

If the `i`-th component is an exact log-affine tail with
`a_i=kappa_i v_i`, the definition of the Cheeger constant gives
`kappa_i>=psi`.  More precisely,

\[
 \boxed{
 P_\mu(A)-\psi\mu(A)
 =\sum_i a_i\left(1-{\psi\over\kappa_i}\right).}       \tag{3.3}
\]

Consequently, for every `tau>0`,

\[
 \boxed{
 \sum_{\{i:\kappa_i\ge(1+\tau)\psi\}}a_i
 \le {1+\tau\over\tau}
       \{P_\mu(A)-\psi\mu(A)\}.}                      \tag{3.4}
\]

The proof is substitution in (3.2); on the displayed bad set each summand
in (3.3) is at least `tau a_i/(1+tau)`.  Thus a separated component carrying
a fixed fraction of the boundary flux cannot have a normal slope separated
from `psi` when the Cheeger deficit is much smaller than that flux.

At the audited scales, deficit `6.02*10^-5p` and selected flux `.004p`
give a relative budget about `.0151`.  For example, (3.4) says that slopes
at least `1.1 psi` can carry at most `.166` of the selected flux even if the
entire scalar deficit is assigned to them.  Iterating smaller slope bins
gives concentration of all but an explicitly controlled flux near the
common Cheeger slope.

Since `c_i=sqrt(s)kappa_i` and `s=alpha K`, while Buser--Ledoux makes
`psi sqrt(K)` bounded below and Cheeger's inequality bounds it above,
common `kappa_i=O(psi)` implies `c_i=O(sqrt(alpha))`.  Such phases cannot
form the widely separated softmax plateaus of Section 7 of
`heat_selector_variation.md` when `alpha` is frozen sufficiently small.

Equations (3.3)--(3.4) are exact only after a low-incidence packet group has
been converted into genuine perimeter-additive components.  Establishing
that conversion with an error controlled by the same ridge/end conductance
is the remaining geometric cut lemma.  This is preferable to assuming a
Poincare bound for the level law: separated softmax modes can make that law
arbitrarily multimodal, so its Poincare constant need not be universal.

### 3.2 A quantitative softmax exclusion after component synchronization

The preceding component charge has enough numerical room to rule out the
`8/17` cross-level branch in the ideal separated-facet model.  The following
elementary lemma tracks the loss from exceptional components.

**Lemma 3.1 (narrow slopes plus a small exceptional flux).**  Let `z` range
over an interval of length `L_z`.  Let a joint phase law have projectors
`Z_i` with `|Z_i|_{HS}=1`.  Suppose a set of good phases has conditional
weights proportional to

\[
                         A_i e^{-c_i^2+\sqrt2c_i z}       \tag{3.5}
\]

and `max_good c_i-min_good c_i<=Delta_c`.  Let `q` be the total joint mass
of all bad phases.  Then, for the conditional projector mean `Q(z)`,

\[
 \boxed{
 \operatorname {Var}Q(z)
 \le2\{e^{\sqrt2L_z\Delta_c}-1\}^2+8q.}              \tag{3.6}
\]

**Proof.**  Let `Q_G(z)` be the projector mean after conditioning on the
good phases.  For two values `z,z'`, the logarithm of the ratio between any
two good phase likelihoods changes by at most
`sqrt(2) Delta_c |z-z'|`.  After normalization, every likelihood ratio
between the two good conditional laws lies between `e^{-ell}` and `e^ell`,
where `ell=sqrt(2)L_z Delta_c`.  Hence their `l^1` distance is at most
`e^ell-1`.  Since `|Z_i|<=1`,

\[
 |Q_G(z)-Q_G(z')|_{HS}
 \le e^{\sqrt2L_z\Delta_c}-1.                         \tag{3.7}
\]

Thus `Var Q_G` is at most the square of the right side.  If `beta(z)` is
the conditional bad-phase probability, then
`|Q(z)-Q_G(z)|<=2beta(z)` and `E beta=q`.  The Hilbert
inequality `Var X<=2 Var Y+2E|X-Y|^2` gives

\[
 \operatorname {Var}Q
 \le2\operatorname {Var}Q_G+8E\beta^2
 \le2\operatorname {Var}Q_G+8q,                       \tag{3.8}
\]

which proves (3.6).  QED.

Now take `tau=1` in (3.4).  If total selected flux is
`A>=.004489p` and the whole scalar deficit is at most
`D=6.02*10^-5p`, the bad phases `kappa_i>=2psi` have joint fraction

\[
 q\le{2D\over A}<.02683.                               \tag{3.9}
\]

Every good exact tail component has `psi<=kappa_i<2psi`.  Cheeger's
inequality gives `psi sqrt(K)<=2`, so with `s=alpha K`,

\[
 \Delta_c\le\sqrt{s}\,\psi\le2\sqrt\alpha=2\,10^{-5}. \tag{3.10}
\]

If one additionally truncates the ideal profile coordinate to an interval
of length less than `10.96`, then

\[
 2\{e^{\sqrt2L_z\Delta_c}-1\}^2<2.0\,10^{-7},
 \qquad 8q<.21464<{8\over17}.                          \tag{3.11}
\]

Thus this *bounded-profile* ideal component model cannot realize the required
cross-level variance.  Centrality of `mu(A_r)` does not itself bound the
profile coordinate `z`, and the globally floored analytic selector has
support at every `z`.  Consequently (3.11) must not be applied directly to
the physical selector.

The unbounded-coordinate repair is proved in `tube_product_inverse.md`.
For the exact Gaussian-mixture phase law and selector floor `omega_0`, it
uses the dimension-free second moment of every good shifted Gaussian and
the pointwise derivative bound `|Q_G'|<=Delta_c/sqrt(2)`.  It gives

\[
 \sqrt{\operatorname {Var}Q}
 \le {\Delta_c\over\sqrt{2\omega_0}}
       \sqrt{1+2C_c^2}+2\sqrt q.                     \tag{3.12}
\]

With cut-plus-tail error at most `1.4*10^-4p`, the audited physical
constants make the square of the right side smaller than `.4557<8/17`.
This repaired estimate requires no bounded `z` interval and no Poincare
inequality for the multimodal level law.

The applicability condition remains load-bearing: the low-incidence
physical packets must be cut into actual components whose volume and
perimeter errors are small enough that (3.3), the log-affine tail relation,
and the heat softmax formula hold after integration.  Lemma 3.1 proves the
complete algebra and constants once that geometric transfer is supplied.

## 4. Product audit and remaining theorem

For independent one-sided exponentials, all coordinate facets have
`kappa_i=1`.  Hence the slope variance in (2.4) is zero and their proportions
do not rotate through the log-affine heat profile.  Their non-affine max-box
cuts are instead improved by the explicit simultaneous ridge bevel in
`fixed_scale_physical_splicing.md`.  The affine coordinate cut is the
one-dimensional product branch and has a universal Cheeger constant.

The finite-dimensional calculation therefore points to the following
geometric theorem, still unproved.

> On the selected physical levels, either phase cores have enough
> `sqrt(s)`-scale cutoff capacity for (2.3), or their failed capacity is
> carried by a bounded-reuse ridge/end/contact measure; after removing that
> charge, every surviving low-capacity component has a common weighted
> normal slope.  Common-slope components have no between-level softmax
> rotation and form the log-affine/product branch.

This is stronger than a generic boundary Poincare inequality and weaker
than the original fixed physical-splice lemma: it only tests the actual
normal-slope field.  Equations (1.2) and (2.4) prove the algebraic and
variational core.  The missing assertion is the physical capacity/incidence
decomposition, together with transfer from the integrated near-minimizer
levels of `F_0` to the log-affine core approximation.
