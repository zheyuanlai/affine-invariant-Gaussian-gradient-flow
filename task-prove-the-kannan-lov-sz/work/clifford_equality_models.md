# Cross-distance equality and the orthogonal-radial branch

This note records a third harmless model, in addition to parallel and
concurrent normal rays.  It arises in the sharp equality case of the
cross-ray Lipschitz inequalities and must be included in any global inverse
theorem.

## 1. Complete bipartite equality is orthogonal-spherical

**Lemma 1 (cross-distance classification).**  Let `B,C` be nonempty subsets
of a Euclidean space and let `d>0`.  Suppose

\[
 |b-c|=d\qquad(b\in B,\ c\in C).                       \tag{1}
\]

Then there are mutually orthogonal linear subspaces `U,V,W`, points
`b_0,c_0` with `w=b_0-c_0\in W`, and numbers `r_B,r_C>=0` such that

\[
 B\subset b_0+r_B S(U),\qquad
 C\subset c_0+r_C S(V),\qquad
 r_B^2+r_C^2+|w|^2=d^2.                                \tag{2}
\]

Here `S(U)` denotes the unit sphere of `U`, with the evident singleton
interpretation when the radius or the subspace is zero.  Conversely, every
configuration in (2) satisfies (1).

**Proof.**  Put

\[
 U=\operatorname {span}(B-B),\qquad
 V=\operatorname {span}(C-C).
\]

For `b,b' in B` and `c,c' in C`, subtracting the four instances of
(1) after squaring gives

\[
 \langle b-b',c-c'\rangle=0.                           \tag{3}
\]

Thus `U` and `V` are orthogonal.  The affine hulls `aff B` and `aff C`
therefore have unique closest affine fibers after quotienting their common
orthogonal complement.  Choose `b_0 in aff B` and `c_0 in aff C` so that

\[
 B-b_0\subset U,\quad C-c_0\subset V,\quad
 w=b_0-c_0\perp U\oplus V.                             \tag{4}
\]

Writing `b=b_0+u` and `c=c_0+v`, orthogonality gives

\[
 d^2=|b-c|^2=|w|^2+|u|^2+|v|^2.                       \tag{5}
\]

Fixing `v` shows that `|u|` is constant on `B-b_0`; fixing `u` shows
the analogous fact on `C-c_0`.  Call the two radii `r_B,r_C`.  Equation
(5) proves (2).  The converse follows by the same orthogonal decomposition.
QED.

The basic Clifford realization takes orthogonal subspaces `U,V` of the same
dimension, `w=0`, and equal radii `d/sqrt(2)`.  Every positive endpoint on
the first sphere is then at distance exactly `d` from every negative endpoint
on the second sphere.  Consequently a large family of calibrated pairs can
be permuted arbitrarily without changing cost, even though the originally
labelled rays need be neither parallel nor concurrent.

There is also an exact distributional stability identity.  It is a useful
target for any bridge argument because it avoids choosing point-set
representatives.

**Lemma 1.1 (ANOVA identity for cross distances).**  Let `B,C` be independent
square-integrable random vectors.  Put

\[
 b=EB,\quad c=EC,\quad U=B-b,\quad V=C-c,
 \quad a=b-c,
\]

and let `Sigma_B,Sigma_C` be their covariance matrices.  Then

\[
\boxed{\begin{aligned}
 \operatorname {Var}|B-C|^2
 ={}&\operatorname {Var}(|U|^2+2\langle a,U\rangle)\\
 &+\operatorname {Var}(|V|^2-2\langle a,V\rangle)
   +4\operatorname {Tr}(\Sigma_B\Sigma_C).
\end{aligned}}                                         \tag{5a}
\]

**Proof.**  Expand

\[
 |B-C|^2=|a|^2+
 (|U|^2+2\langle a,U\rangle)+
 (|V|^2-2\langle a,V\rangle)-2\langle U,V\rangle.
\]

The first two random summands are independent.  Their covariances with
`<U,V>` vanish after conditioning, because `EU=EV=0`.  Finally

\[
 E\langle U,V\rangle^2=\operatorname {Tr}(\Sigma_B\Sigma_C).
\]

Adding the four orthogonal `L^2` components proves (5a).  QED.

Consequently, if the product distribution of positive and negative endpoint
packets has cross-distance squared almost constant, then both endpoint laws
are concentrated near (possibly off-center) spheres and their covariance
subspaces have small Hilbert--Schmidt overlap.  Exact zero variance recovers
Lemma 1 at the measure level.  The unresolved implication is the preceding
one: current midpoint-set inequalities do not bound the product-law variance
in (5a).

## 2. Orthogonal-radial witnesses are dimension-free

We use the dimension-free translated thin-shell consequence

\[
 \sup_z \operatorname {Var}|Y-z|\le C_{TS}             \tag{6}
\]

for every isotropic log-concave random vector `Y` in every dimension.  It
applies intrinsically to every orthogonal marginal of an isotropic
log-concave law.

**Lemma 2 (finitely many radial coordinates).**  Let `X` be isotropic and
log-concave in `R^n`.  Let `P_1,...,P_k` be mutually orthogonal projections,
let `z_i in ran P_i`, and put

\[
 R_i=|P_iX-z_i|.
\]

If `Phi:R^k->R` is one-Lipschitz for the Euclidean norm, then

\[
 \mathbb E|\Phi(R)-\mathbb E\Phi(R)|
 \le 2\sqrt{kC_{TS}}.                                  \tag{7}
\]

In particular the constant is universal for the two-block orthogonal-radial
models furnished by Lemma 1.

**Proof.**  Each marginal `P_iX` is isotropic and log-concave on its range,
so (6) gives `Var R_i<=C_TS`.  With `m_i=ER_i`, Lipschitzness and
Cauchy--Schwarz give

\[
 \begin{aligned}
 E|\Phi(R)-E\Phi(R)|
 &\le 2E|\Phi(R)-\Phi(m)|\\
 &\le2E|R-m|_2
 \le2\left(\sum_i\operatorname {Var}R_i\right)^{1/2}.
 \end{aligned}                                         \tag{8}
\]

This is (7).  QED.

**Corollary 3 (constant-mass Clifford tails have bounded separation).**
Let `B,C` satisfy Lemma 1, and let an isotropic log-concave probability
`mu` obey

\[
 \mu(B)\ge\delta,\qquad\mu(C)\ge\delta.                \tag{9}
\]

Then

\[
 d\le C/\delta,                                         \tag{10}
\]

where `C` is universal.  The same conclusion holds if `B,C` are replaced by
sets on which the three radial/linear coordinates below are constant with
the values specified by (2).

**Proof.**  Use the decomposition in Lemma 1.  If `w!=0`, put
`e=w/|w|`; omit the third coordinate when `w=0`.  Define the contraction

\[
 T(x)=\left(
 |P_U(x-b_0)|,
 |P_V(x-c_0)|,
 \langle x-c_0,e\rangle
 \right).                                               \tag{11}
\]

The three underlying projections have orthogonal ranges, and radial
projection is a contraction.  Hence `T` is one-Lipschitz from Euclidean
space to Euclidean three-space.  On `B` and `C`, respectively,

\[
 T(B)=(r_B,0,|w|),\qquad T(C)=(0,r_C,0),                \tag{12}
\]

and the distance between these two points is `d`.  Let `ell` be the unit
linear functional in their difference direction and put `g=ell o T`.
Then `g` is one-Lipschitz and its two constant values on `B,C` differ by
`d`.

On the other hand, translated thin shell bounds the variances of the first
two coordinates of `T`, while isotropy gives variance one for the last.
Therefore

\[
 E|g-Eg|\le\sqrt{\operatorname {Var}g}\le C.            \tag{13}
\]

If `m=Eg`, the triangle inequality and (9) give

\[
 E|g-m|\ge
 \delta\big(|g(B)-m|+|g(C)-m|\big)\ge\delta d.          \tag{14}
\]

Equations (13)--(14) prove (10).  QED.

The positive-mass version needed in applications is the following immediate
robust form.  Let `t_B,t_C` denote the two feature points in (12), and assume
instead that

\[
 E[|T(X)-t_B|\mid B]+E[|T(X)-t_C|\mid C]\le\varepsilon d. \tag{15}
\]

For the same `g`, conditional Jensen and the triangle inequality replace
(14) by

\[
 E|g-Eg|\ge\delta(1-\varepsilon)d.                     \tag{16}
\]

Thus `d<=C/[delta(1-epsilon)]` whenever `epsilon<1`.  Unlike the exact
spherical statement, (15) is compatible with full-dimensional tail bands.
It is the appropriate endpoint of an approximate inverse theorem.

Thus an inverse theorem which concludes only “parallel or concurrent” is
false at the exact algebraic level.  The correct harmless list also includes
a two-block orthogonal-radial branch (and degenerate subspace versions).
What remains unproved is a robust theorem converting high weighted medial
conductance and global log-concave bridge mass into proximity to one of these
classified equality configurations.
