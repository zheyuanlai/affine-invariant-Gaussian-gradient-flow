# Ridge dimension reduction for a T3 witness

## 0. Verdict

A sufficiently accurate ridge approximation would close KLS by a short and
fully quantitative self-consistency argument.  There are two corrections to
the proposed formulation.

1. If the hypothesis is the raw error

   \[
          \mathbb E|f-\phi(PX)|\le cA,
   \]

   then the argument requires `c<1/2`, not merely `c<1`.  A centered
   approximation with error `cA` requires only `c<1`.
2. Near extremality, variance comparability, isotropy, and tensorization do
   not themselves imply such a ridge theorem.  A tensor power of any
   genuinely nonlinear hypothetical bad block remains a constant-factor T3
   near-extremizer at the same `A` scale, while every rank `O(A^2)` linear
   observation asymptotically misses its nonlinear regression residual.

The second statement is formalized below.  It is a decisive no-go for
gradient-covariance truncation, conditional expectation, and random
projection as *generic* proofs of the ridge assertion.  It says that a valid
ridge theorem must first prove a new rigidity result: every bad T3 extremizer
has almost all of its variance in one linear observation, or belongs to a
separately controlled convex/translated-radial class.

There cannot be an unconditional explicit log-concave example with
`A->infinity` unless KLS is false.  Instead, Proposition 4.1 gives the sharp
countermodel inside the contradiction framework: from any putative
nonlinear bad sequence with `A->infinity`, it constructs an isotropic
log-concave tensor sequence with `A->infinity` that violates every prescribed
`O(A^2)` ridge assertion.  This is exactly the setting in which the proposed
lemma would be invoked and does not assume a marginal-to-whole KLS transfer.

---

## 1. The exact self-consistency theorem

Let

\[
 D_1(\mu)=\sup_{\operatorname {Lip}(f)\le1}
              \mathbb E_\mu|f-\mathbb E_\mu f|.                \tag{1.1}
\]

Use the known dimension-dependent KLS estimate in the following explicit
form: there is a numerical `c_K>0` such that every isotropic log-concave
probability `nu` in dimension `k` satisfies

\[
 \psi_\nu\ge\frac{c_K}{\sqrt{\log(e+k)}}.                       \tag{1.2}
\]

Cheeger's inequality then gives

\[
 C_P(\nu)\le\frac{4}{c_K^2}\log(e+k).                          \tag{1.3}
\]

**Theorem 1.1 (ridge approximation implies a universal T3 bound).**
Suppose there are constants `C_r<infinity` and `eta<1/2` such that the
following holds whenever `mu` is isotropic log-concave and `f` is a
1-Lipschitz witness in the class under consideration.  With

\[
                       A=\mathbb E|f-\mathbb Ef|,               \tag{1.4}
\]

there are an orthogonal projection `P` of rank

\[
                       k\le C_rA^2                              \tag{1.5}
\]

and a 1-Lipschitz `phi` on `ran P` such that

\[
                       \mathbb E|f-\phi(PX)|\le\eta A.         \tag{1.6}
\]

Then every such witness satisfies

\[
 A\le b+\sqrt{b\log(e+C_r)},\qquad
 b=\frac{4}{c_K^2(1-2\eta)^2}.                                \tag{1.7}
\]

In particular, if every uncontrolled near-extremizer satisfies
(1.5)--(1.6), while the exceptional convex/translated-radial branches have a
universal T3 bound, then KLS follows through Milman's equivalence theorem.

**Proof.**  Put `Y=PX` and `h=phi(Y)`.  The marginal law of `Y` is isotropic
and log-concave on the `k`-dimensional range of `P`.  By (1.3),

\[
 \mathbb E|h-\mathbb Eh|
 \le\sqrt{\operatorname {Var}(h)}
 \le\sqrt{C_P(P_\#\mu)}
 \le\frac2{c_K}\sqrt{\log(e+k)}.                              \tag{1.8}
\]

For arbitrary integrable random variables,

\[
 \mathbb E|f-\mathbb Ef|
 \le2\mathbb E|f-h|+\mathbb E|h-\mathbb Eh|.                  \tag{1.9}
\]

Therefore (1.5)--(1.8) give

\[
 (1-2\eta)A
 \le\frac2{c_K}\sqrt{\log(e+C_rA^2)}.                         \tag{1.10}
\]

After squaring, put `b=4/[c_K^2(1-2eta)^2]`.  Since

\[
 e+C_rA^2\le(e+C_r)(1+A^2),\qquad
 \log(1+A^2)\le A,                                           \tag{1.11}
\]

one has

\[
                         A^2\le b\{\log(e+C_r)+A\}.            \tag{1.12}
\]

Solving this quadratic gives the slightly sharper bound

\[
 A\le\frac{b+\sqrt{b^2+4b\log(e+C_r)}}2,
\]

which implies (1.7). QED.

If instead the ridge theorem supplies the centered estimate

\[
 \mathbb E\left|(f-\mathbb Ef)-
   \{\phi(PX)-\mathbb E\phi(PX)\}\right|\le\eta A,            \tag{1.13}
\]

then (1.9) has coefficient one.  The same proof works for every `eta<1`,
with

\[
                         b=\frac4{c_K^2(1-\eta)^2}.             \tag{1.14}
\]

The factor two in (1.9) cannot be silently removed from the raw formulation:
`|Ef-Eh|<=E|f-h|` is a second, independent centering cost.

No marginal-to-whole KLS statement appears in this proof.  Only KLS in the
actual marginal dimension `k` is used through the already known estimate
(1.2).

---

## 2. What rank `O(A^2)` means

Two elementary facts calibrate the proposed rank.

First, for every rank-`k` projection and every 1-Lipschitz `phi`,

\[
 \mathbb E|f(X)-\phi(PX)|\ge\frac12(A-\sqrt k)_+.              \tag{2.1}
\]

Indeed `Var(phi(PX))<=k`, and (1.9) gives (2.1).  Thus rank proportional to
`A^2` is the first scale at which a fixed relative approximation is even
possible.

Second, if `n<=C A^2`, taking `P=I` is already an exact ridge
representation.  The assertion has content only in the regime

\[
                              n\gg A^2.                         \tag{2.2}
\]

For a true T3 extremizer, the known structural theorem gives
`Var(f)asymp A^2` and constant-mass positive and negative tails separated by
`Theta(A)`.  Independent samples from those tails have a difference
covariance of stable rank `Omega(A^2)`.  This is a *lower* rank statement;
it does not say that the separation is contained in an `O(A^2)`-dimensional
linear subspace.  Turning stable rank into such an upper-dimensional ridge
is the unproved step.

---

## 3. Linear regression and the nonlinear residual

Let `nu` be any isotropic probability on `R^d`, let `f` be centered and
square-integrable, and define

\[
 a=\mathbb E[Xf(X)]\in\mathbb R^d,
 \qquad r(X)=f(X)-\langle a,X\rangle.                           \tag{3.1}
\]

Isotropy gives the exact orthogonal decomposition

\[
 \mathbb E[Xr(X)]=0,
 \qquad
 \sigma^2:=\mathbb Er^2=\operatorname {Var}(f)-|a|^2\ge0.    \tag{3.2}
\]

The scalar `sigma` is the part of the witness invisible to all linear
observations at the covariance level.  If

\[
                         \sigma\ge\theta A                      \tag{3.3}
\]

for a numerical `theta>0`, call the block genuinely nonlinear.  If (3.3)
fails with a very small `theta`, then the witness has an almost rank-one
linear-regression description in `L^2`; proving that this is a controlled
1-Lipschitz linear branch is a separate inverse problem, but it is already
far stronger structure than near extremality alone.

The tensor obstruction preserves exactly the residual (3.2).

---

## 4. Tensor amplification theorem

Let `X_1,...,X_m` be independent with law `nu`, and put

\[
 \mu_m=\nu^{\otimes m},\qquad
 F_m(X_1,\ldots,X_m)=\frac1{\sqrt m}\sum_{j=1}^mf(X_j).        \tag{4.1}
\]

If `f` is 1-Lipschitz, then `F_m` is 1-Lipschitz because

\[
             \sum_{j=1}^m|\nabla_jF_m|^2
             =\frac1m\sum_{j=1}^m|\nabla f(X_j)|^2\le1.       \tag{4.2}
\]

Also

\[
 \operatorname {Var}(F_m)=V:=\operatorname {Var}_\nu(f),
 \qquad
 \mathbb E|F_m|\longrightarrow\sqrt{\frac{2V}{\pi}}.         \tag{4.3}
\]

The second statement is the scalar central limit theorem plus uniform
integrability.

**Proposition 4.1 (fixed-rank linear observations miss the nonlinear
residual).**  In the setup (3.1)--(4.1), for every fixed integer `k`,

\[
\boxed{
 \liminf_{m\to\infty}
 \inf_{\substack{P:\,\operatorname {rank}P\le k\\
                   \operatorname {Lip}(\phi)\le1}}
 \mathbb E\left|F_m-\phi(P(X_1,\ldots,X_m))\right|
 \ge\sqrt{\frac2\pi}\,\sigma.}                              \tag{4.4}
\]

Allowing an arbitrary measurable predictor in the limiting Gaussian
experiment gives the same lower bound.

**Proof.**  Choose orthonormal coordinates on `ran P`.  Then the observation
has the form

\[
 W=BX=\sum_{j=1}^mB_jX_j\in\mathbb R^k,
 \qquad BB^{\mathsf T}=I_k,
 \qquad \sum_j\|B_j\|_{\rm HS}^2=k.                            \tag{4.5}
\]

Take a sequence `delta_m downarrow0` so slowly that
`k/delta_m^2=o(m)`, and call a block high leverage if
`||B_j||HS>delta_m`.  There are at most `k/delta_m^2=o(m)` high-leverage
blocks.  Their contribution to `F_m` tends to zero in `L^2`.

Give the predictor the high- and low-leverage observations separately.
This only makes prediction easier.  The high-leverage observation is
independent of all low blocks.  On the low blocks, the triangular-array
Lindeberg theorem applies jointly to

\[
 \left(\frac1{\sqrt m}\sum_{j\in L}f(X_j),
       \sum_{j\in L}B_jX_j\right),                             \tag{4.6}
\]

because the largest block coefficient tends to zero and an isotropic
log-concave law has finite moments of every fixed order.  After passage to a
subsequence, (4.6) converges in `W_1` to a centered jointly Gaussian pair
`(Z,W_L)`.

Its covariance satisfies

\[
 \operatorname {Var}(Z)=V,
 \qquad
 \operatorname {Var}(Z\mid W_L)\ge V-|a|^2=\sigma^2.         \tag{4.7}
\]

To verify the second inequality, put

\[
 A_m=\frac1{\sqrt m}(a,\ldots,a)
\]

on the low block space.  The variance explained by `W_L` is the squared
norm of the orthogonal projection of `A_m` onto the row space of the low
observation matrix.  It is at most `|A_m|^2<=|a|^2`.

The high-leverage limiting observation is independent of `(Z,W_L)` and
therefore cannot reduce (4.7).  Conditional on all limiting observations,
the remaining Gaussian has standard deviation at least `sigma`; its least
absolute-deviation predictor is its conditional mean/median and has error at
least `sqrt(2/pi)sigma`.

Finally, the original predictors are 1-Lipschitz.  The loss

\[
                         (z,w)\longmapsto|z-\phi(w)|            \tag{4.8}
\]

has a universal Lipschitz constant, uniformly in `phi`.  Thus `W_1`
convergence passes the lower bound to every sequence of admissible
predictors.  Tightness permits extraction of a locally uniform subsequence
after harmlessly bounding `phi(0)` by comparison with the zero predictor.
This proves (4.4). QED.

The same proof works for `k=k_j` varying along an outer sequence: for each
fixed outer index choose the tensor power `m` only after `k_j` is fixed.

### 4.1 A countersequence inside the KLS contradiction

Suppose there is a hypothetical sequence of isotropic log-concave measures
`nu_j` with 1-Lipschitz T3 near-extremizers `f_j` such that

\[
 A_j=\mathbb E|f_j-\mathbb Ef_j|\longrightarrow\infty,
 \qquad V_j\asymp A_j^2,
 \qquad \sigma_j\ge\theta A_j.                               \tag{4.9}
\]

Fix any proposed rank constant `C_r` and put

\[
                         k_j=\lceil C_rA_j^2\rceil.            \tag{4.10}
\]

Choose `m_j` sufficiently large in Proposition 4.1.  Then

\[
 \widetilde\mu_j=\nu_j^{\otimes m_j},\qquad
 \widetilde f_j=m_j^{-1/2}\sum_{l=1}^{m_j}f_j(X_l)            \tag{4.11}
\]

are isotropic log-concave and 1-Lipschitz, and

\[
 \mathbb E|\widetilde f_j-\mathbb E\widetilde f_j|
 \asymp A_j,\qquad
 \operatorname {Var}(\widetilde f_j)\asymp A_j^2.             \tag{4.12}
\]

For every rank-`k_j` projection and every 1-Lipschitz ridge,

\[
 \mathbb E|\widetilde f_j-\phi(P\widetilde X_j)|
 \ge\left(\sqrt{\frac2\pi}\theta-o(1)\right)A_j.            \tag{4.13}
\]

If `f_j` is a constant-factor T3 near-extremizer, so is `tilde f_j`, with
a changed universal factor.  Indeed Bobkov--Houdre tensorization gives

\[
 \psi_{\nu_j^{\otimes m}}\asymp\psi_{\nu_j},                  \tag{4.14}
\]

and Milman's equivalence gives

\[
 D_1(\nu_j^{\otimes m})\asymp D_1(\nu_j).                     \tag{4.15}
\]

Equations (4.3), (4.9), and near extremality of `f_j` then prove the claim.

The construction does not manufacture a negative solution of KLS from
nothing; it transforms any putative bad sequence, which is exactly the
starting assumption of the affirmative contradiction argument.  It proves
that near extremality plus `V asymp A^2` cannot by itself yield the desired
ridge approximation whenever a fixed fraction of the variance is genuinely
nonlinear.

If an exact branch distinction is desired, note that restriction of
`tilde f_j` to one block is a translate and scalar multiple of `f_j`.
Therefore `tilde f_j` cannot be convex if `f_j` is not convex.  A separable
sum over at least two blocks cannot be a nonconstant translated-radial
function unless its block summands are quadratic; global Lipschitzness then
forces the quadratic coefficient to vanish.  Thus exact nonconvex,
nonradial bad blocks remain outside those two exceptional branches.

---

## 5. Why the four proposed generic mechanisms stop

### 5.1 Gradient covariance

Let

\[
                         G_f=\mathbb E[\nabla f\nabla f^{\mathsf T}].
\]

For the tensor witness (4.1),

\[
 G_{F_m}=\frac1m\operatorname {diag}(G_f,\ldots,G_f).          \tag{5.1}
\]

Every eigenvalue is divided by `m` and repeated `m` times.  A rank `k`
top-eigenspace captures at most

\[
                         \frac{k}{m}\|G_f\|_{\rm op}          \tag{5.2}
\]

of the total gradient energy.  For fixed `k=O(A^2)` and arbitrarily large
`m`, this tends to zero.  Hence a gradient-PCA proof is incompatible with
tensorization unless it first proves that the nonlinear tensor witness is
not near-extremal.

### 5.2 Conditional expectation

The variance identity

\[
 \operatorname {Var}(f)
 =\mathbb E\operatorname {Var}(f\mid PX)
  +\operatorname {Var}(\mathbb E[f\mid PX])                   \tag{5.3}
\]

does not bound the first term by discarded gradient energy.  Such a bound
would require uniform Poincare inequalities for the affine fibers.  They do
not follow from isotropy or log-concavity.

Nor is the conditional mean automatically 1-Lipschitz.  For a smooth
density `e^{-V(y,z)}`, differentiation gives

\[
 \nabla_y\mathbb E[f(y,Z)\mid y]
 =\mathbb E[\nabla_yf\mid y]
  -\operatorname {Cov}(f,\nabla_yV\mid y),                    \tag{5.4}
\]

with an additional boundary flux for moving hard supports.  Controlling the
score covariance in (5.4) is another curvature/fiber hypothesis, not a free
consequence of (global) Poincare.  Proposition 4.1 is stronger: even the
best Lipschitz predictor of an arbitrary rank-`k` linear observation misses
the nonlinear residual after tensor amplification.

### 5.3 Random projections

A random rank-`k` orthogonal projection preserves only a fraction `k/n` of
a typical squared distance.  Taking `k asymp A^2` makes its typical projected
distance of order `A`, but does not preserve the pairs on which `f` changes
by order `A`.  In the tensor model, delocalized random projections fall
exactly under the low-leverage central limit argument in Proposition 4.1:
their observations become jointly Gaussian with the witness and detect only
the linear regression vector `a`.

### 5.4 Tensorization

Tensorization controls the concentration constant of a product up to
universal factors, but it does not imply that a near-extremizing function
uses one factor.  The normalized sum (4.1) has the same variance and
first-moment scale as one block.  Any proposed "irreducibility" reduction
must therefore be proved for near-extremizers; it cannot be imposed by
discarding tensor sums.

---

## 6. What remains viable

The ridge route is a valid sufficient mechanism in the following precise
form.

> For some `eta<1/2` and `C_r<infinity`, every genuinely uncontrolled T3
> near-extremizer admits the raw approximation (1.5)--(1.6).

Theorem 1.1 would then give KLS immediately.  Proposition 4.1 shows that a
proof must use more than near extremality and second moments.  It must
establish at least one of these new statements.

1. **Linear-regression rigidity:** every hypothetical bad near-extremizer
   has `sigma=o(A)` in (3.2), followed by a proof that this produces an
   admissible 1-Lipschitz low-rank ridge.
2. **Tensor-sum exclusion:** a normalized sum over many nonlinear bad
   blocks is quantitatively separated from the true T3 optimum, despite
   (4.14)--(4.15).  Ordinary tensorization supplies only constant-factor
   comparison and cannot do this.
3. **Extremality-dependent convex/radial capture:** the nonlinear residual
   forces the witness into a convex or translated-radial class with a
   dimension-free first moment bound.

Absent one of these, the desired dimension reduction is not a consequence
of gradient covariance, conditional expectation, random projections, or
tensorization.  Conversely, proving any version with the centered error in
(1.13) and `eta<1`, or the raw error in (1.6) and `eta<1/2`, supplies the
explicit self-consistency contradiction (1.7).
