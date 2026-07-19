# Normalized self-convolution as a KLS renormalization

## Executive conclusion

Let `X,Y` be independent with the same centered isotropic log-concave law
`mu`, and put

\[
             S={X+Y\over\sqrt2},\qquad D={X-Y\over\sqrt2},
             \qquad \mathsf T\mu={\cal L}(S).                    \tag{0.1}
\]

The normalized convolution operator `T` has a clean but limited rigidity
theory.

* The Poincare constant is exactly nonincreasing:
  \[
                      C_P(\mathsf T\mu)\le C_P(\mu).             \tag{0.2}
  \]
  More strongly, every test function has an explicit nonnegative four-copy
  defect.  If equality of the two constants is attained, then the witness is
  affine, the common constant is one, and `mu` is Gaussian.
* No dimension-free quantitative contraction of the excess
  `C_P-1` is obtained.  Such a result requires a uniform stability theorem
  for an approximate additive Cauchy equation.  The exact defect isolates
  that missing statement; neither log-concavity nor the conditional symmetry
  below presently supplies it.
* The centered first-moment constant has no monotonicity in either direction.
  In one dimension its exact value is `E|X-EX|`.  A variance-one Laplace law
  increases from `1/sqrt(2)` to `3/4` after one normalized convolution,
  whereas a variance-one uniform law decreases from `sqrt(3)/2` to
  `sqrt(6)/3`.
* Fourth cumulants do contract exactly: the order-four cumulant tensor is
  halved, and
  \[
       \operatorname{Var}|S|^2-2n
       ={1\over2}\bigl(\operatorname{Var}|X|^2-2n\bigr).         \tag{0.3}
  \]
  This `M4` contraction does not control an arbitrary Poincare or Lipschitz
  witness.
* Conditionally on `S=s`, the law of `D` is symmetric and log-concave, but
  only its *average* covariance is fixed.  Trying to control the original
  witness on these fibers is precisely a symmetric conditional-KLS problem,
  not a free contraction argument.

Thus normalized convolution supplies an exact fixed-point theorem and a
useful stability certificate, but not a complete KLS renormalization.  It
also gives a decisive no-go for direct `M1` monotonicity.

## 1. Elementary structure of `(S,D)`

If `mu` has density `rho`, the orthogonal change of variables in (0.1) gives

\[
 p_{S,D}(s,d)=\rho\left({s+d\over\sqrt2}\right)
               \rho\left({s-d\over\sqrt2}\right).               \tag{1.1}
\]

Consequently the conditional law `kappa_s` of `D` given `S=s` has density

\[
 {d\kappa_s\over dd}(d)
 ={1\over Z(s)}
   \rho\left({s+d\over\sqrt2}\right)
   \rho\left({s-d\over\sqrt2}\right).                           \tag{1.2}
\]

It is symmetric under `d -> -d` and log-concave.  In particular,

\[
                 E[D\mid S]=0,
 \qquad E\,\operatorname{Cov}(D\mid S)=I.                       \tag{1.3}
\]

The second identity is only an average identity.  It gives no pointwise
upper bound for `Cov(D|S=s)`.

There is an exact fiber representation for any square-integrable `g`.  Set

\[
 A_s(d)={g((s+d)/\sqrt2)-g((s-d)/\sqrt2)\over\sqrt2}.            \tag{1.4}
\]

Then `A_s` is odd, so `E_{kappa_s}A_s=0`, and

\[
 \boxed{\quad
 \operatorname{Var}_\mu g
 =E_S\operatorname{Var}_{\kappa_S} A_S.
 \quad}                                                         \tag{1.5}
\]

If `g` is one-Lipschitz, then `A_s` is one-Lipschitz as a function of `d`.
For first moments,

\[
 {1\over\sqrt2}E|g-Eg|
 \le E_{S,D}|A_S(D)|
 \le\sqrt2 E|g-Eg|.                                            \tag{1.6}
\]

The two inequalities are the usual independent-copy symmetrization and the
triangle inequality.

Equations (1.3)--(1.6) explain both the appeal and the obstruction of the
conditional route.  A bad witness becomes an odd witness on symmetric
log-concave fibers, but those fibers are not isotropic and their covariance
is controlled only after averaging over `S`.  Applying a uniform Poincare or
concentration theorem to them would insert precisely the symmetric KLS input
which this route is supposed to prove.

## 2. Poincare nonexpansion and an exact four-copy defect

Let

\[
 C=C_P(\mu),\qquad \nu=\mathsf T\mu.
\]

For a locally Lipschitz `f` on the support of `nu`, write

\[
 F(x,y)=f\left({x+y\over\sqrt2}\right),\qquad m=E F,             \tag{2.1}
\]

and form its Hoeffding decomposition under `mu tensor mu`:

\[
 h(x)=E_YF(x,Y)-m,
 \qquad
 r(x,y)=F(x,y)-m-h(x)-h(y).                                    \tag{2.2}
\]

Then

\[
 E h=0,qquad E[r\mid X]=E[r\mid Y]=0,                          \tag{2.3}
\]

and the three summands are pairwise orthogonal.  Hence

\[
 \operatorname{Var}_\nu f=2\|h\|_{L^2(\mu)}^2
                              +\|r\|_{L^2(\mu^2)}^2.             \tag{2.4}
\]

On the other hand, conditional Poincare in the `Y` variable gives

\[
\begin{split}
 E\operatorname{Var}(F\mid X)
 &=\|h\|_2^2+\|r\|_2^2\\
 &\le {C\over2}E_\nu|\nabla f|^2.                               \tag{2.5}
\end{split}
\]

Combining (2.4)--(2.5) proves the following.

**Theorem 2.1 (convolution defect inequality).**  For every such `f`,

\[
 \boxed{\quad
 \operatorname{Var}_\nu f+\|r\|_2^2
 \le C_P(\mu)E_\nu|\nabla f|^2.
 \quad}                                                         \tag{2.6}
\]

In particular,

\[
                        C_P(\mathsf T\mu)\le C_P(\mu).          \tag{2.7}
\]

The residual has an exact four-copy expression.  For independent
`X,X',Y,Y'`, put

\[
\begin{split}
 \Delta_4 f={}&f((X+Y)/\sqrt2)+f((X'+Y')/\sqrt2)\\
 &-f((X+Y')/\sqrt2)-f((X'+Y)/\sqrt2).
                                                                    \tag{2.8}
\end{split}
\]

All one-variable Hoeffding terms cancel in `Delta_4`.  The four remaining
copies of `r` are pairwise orthogonal by (2.3), so

\[
                     E(\Delta_4f)^2=4\|r\|_2^2.                  \tag{2.9}
\]

Thus (2.6) is equivalently

\[
 \boxed{\quad
 \operatorname{Var}_\nu f+{1\over4}E(\Delta_4f)^2
 \le C_P(\mu)E_\nu|\nabla f|^2.
 \quad}                                                        \tag{2.10}
\]

This identity is dimension free and requires no conditional covariance
bound.  If `f` has Rayleigh quotient `R(f)`, it quantitatively says

\[
 {E(\Delta_4f)^2\over4E|\nabla f|^2}
 \le C_P(\mu)-R(f).                                            \tag{2.11}
\]

Accordingly, failure of strict contraction forces an approximate additive
Cauchy equation on four independent copies.

### 2.1 Exact equality and the Gaussian fixed point

Suppose `mu` is full-dimensional, equality

\[
 C_P(\nu)=C_P(\mu)<\infty                                      \tag{2.12}
\]

is attained by a nonconstant `f` for `nu`.  Equations (2.10)--(2.12) give
`Delta_4f=0` almost surely.  The support of a full-dimensional log-concave
law has convex interior.  The four-point equation on an open set implies,
in the sense of distributions, that every mixed second derivative of `f`
vanishes.  Hence `f` is affine on the support of `nu`.

Since `nu` is isotropic, an affine function has Rayleigh quotient one.
Thus the common constant in (2.12) is one.  Every coordinate function then
attains equality in the Poincare inequality for `mu`.  Polarizing the
first variation at a coordinate gives the Stein identities

\[
                    E[X_i\phi(X)]=E[\partial_i\phi(X)]           \tag{2.13}
\]

for all smooth compactly supported `phi`.  These characterize the standard
Gaussian.  We have proved:

**Corollary 2.2 (attained fixed-point rigidity).**  Attained equality in
(2.7) occurs only for the Gaussian law and affine witnesses.

This does not yet give a universal strict contraction.  If the relevant
Poincare constants lie at continuous spectral edges, there need not be an
attaining witness.  Equation (2.11) then yields a sequence of approximate
Cauchy equations.  A dimension-free theorem turning small `Delta_4` into
closeness to an affine function, with the correct gradient normalization,
is additional content.  No such estimate is being assumed here.

## 3. Exact `M4` contraction and distributional fixed points

Cumulants linearize normalized convolution.  For every order `k>=2`,

\[
             \operatorname{Cum}_k(S)
 =2^{1-k/2}\operatorname{Cum}_k(X).                             \tag{3.1}
\]

In particular, the third cumulant tensor is multiplied by `1/sqrt(2)` and
the fourth cumulant tensor by `1/2`.  For every unit vector `u`,

\[
 E\langle S,u\rangle^4-3
 ={1\over2}\left(E\langle X,u\rangle^4-3\right).                \tag{3.2}
\]

There is also a radial identity.  Since

\[
 |S|^2={|X|^2+|Y|^2+2\langle X,Y\rangle\over2}
\]

and `E<X,Y>^2=n`, all cross terms vanish and

\[
 \operatorname{Var}|S|^2
 ={1\over2}\operatorname{Var}|X|^2+n.                           \tag{3.3}
\]

Subtracting the Gaussian value `2n` gives (0.3).

At the level of laws, a finite-variance fixed point is also rigid.  If
`T mu=mu`, its characteristic function satisfies

\[
                         \widehat\mu(t)=\widehat\mu(t/\sqrt2)^2. \tag{3.4}
\]

Iterating (3.4) and using the centered covariance-one expansion at the
origin gives

\[
                         \widehat\mu(t)=e^{-|t|^2/2}.             \tag{3.5}
\]

Thus the only centered finite-covariance fixed point is Gaussian.

These are genuine `M4` and full-law fixed-point theorems.  Their limitation
is that fourth cumulants do not control a general nonlinear Rayleigh or
Lipschitz extremizer.  Upgrading (3.2)--(3.3) to such control would itself
require a new functional inequality.

## 4. `M1` has no convolution monotonicity

The one-dimensional centered first-moment problem is exactly solvable.

**Lemma 4.1 (one-dimensional `T3` extremizer).**  For every atomless
probability on the line with a finite first moment,

\[
 \boxed{\quad
 \sup_{\operatorname{Lip}(f)\le1}E|f-Ef|
 =E|X-EX|.
 \quad}                                                        \tag{4.1}
\]

The affine function `f(x)=x` is an extremizer.

**Proof.**  The cut--transport formula gives

\[
 {\cal D}(\mu)=2\sup_E p(1-p)
        W_1(\mu_E,\mu_{E^c}),\qquad p=\mu(E).                    \tag{4.2}
\]

For fixed `p`, the one-dimensional monotone-rearrangement inequality says
that the right side is maximized when `E` is a lower or upper halfline.
For completeness, approximate `mu` by an ordered finite measure and use the
monotone optimal matching formula for `W_1`.  Exchanging any inverted pair
of labels moves the selected point outward and the unselected point in the
opposite direction, and cannot decrease any matched absolute gap.  Iterating
leaves a threshold set; passage to atomless limits proves the claim.

After translating so `EX=0`, take `E=(-infinity,t]`.  The conditional laws
are ordered, and hence

\[
 2p(1-p)W_1(\mu_E,\mu_{E^c})
 =-2\int_{-\infty}^t x\,d\mu(x).                               \tag{4.3}
\]

The derivative in `t` is `-2t rho(t)` in the density case, so the maximum
occurs at `t=0`; approximation gives the general case.  The value is
`E|X|`.  QED.

Two elementary examples now give opposite behavior.

### 4.1 Laplace increases

Let `X` have the variance-one symmetric Laplace density

\[
                  p(x)={1\over\sqrt2}e^{-\sqrt2|x|}.             \tag{4.4}
\]

Then

\[
                         {\cal D}(\mu)=E|X|={1\over\sqrt2}.      \tag{4.5}
\]

The normalized sum has density

\[
                  q(s)={1\over2}(1+2|s|)e^{-2|s|},              \tag{4.6}
\]

and therefore

\[
                         {\cal D}(\mathsf T\mu)=E|S|={3\over4}. \tag{4.7}
\]

Thus `M1` increases by the factor `3sqrt(2)/4>1`.

### 4.2 Uniform decreases

Let `X` be uniform on `[-sqrt(3),sqrt(3)]`.  Then

\[
                  {\cal D}(\mu)={\sqrt3\over2}.                  \tag{4.8}
\]

For independent copies, `E|X+Y|=2sqrt(3)/3`, so

\[
                  {\cal D}(\mathsf T\mu)={\sqrt6\over3}
                  <{\sqrt3\over2}.                              \tag{4.9}
\]

Therefore neither nonincrease nor nondecrease is available for centered
`L1`, even in dimension one.  In particular, a direct `M1` renormalization
cannot underlie a general KLS proof.

## 5. The product-Laplace Poincare test

The same Laplace example behaves quite differently for the Poincare
constant.

### 5.1 Input constant

For the density `(a/2)e^{-a|x|}`, the exact Poincare constant is `4/a^2`.
This follows from the one-dimensional Hardy--Muckenhoupt criterion, with
equality approached by truncated exponential tail functions.  At
`a=sqrt(2)`,

\[
                              C_P(\mu)=2.                        \tag{5.1}
\]

### 5.2 A rigorous bracket after one convolution

For `q` in (4.6), symmetry reduces the Muckenhoupt quantity to

\[
\begin{split}
 B_*&=\sup_{x>0}\ q([x,\infty))\int_0^x{dt\over q(t)}\\
 &=\sup_{x>0}(x+1)e^{-2x}
        \int_0^x{e^{2t}\over1+2t}\,dt.                          \tag{5.2}
\end{split}
\]

The one-dimensional criterion gives

\[
                             B_*\le C_P(q)\le4B_*.               \tag{5.3}
\]

An elementary bound is `B_*<=49/100`.  Indeed, put

\[
 J(x)=\int_0^x{e^{2t}\over1+2t}dt,
 \qquad H(x)={49\over100}{e^{2x}\over x+1}-J(x).                \tag{5.4}
\]

The derivative of `H` changes sign only at `x=3/4`, from negative to
positive.  Convexity of `t ->(1+2t)^{-1}` gives

\[
 J(x)\le {e^{2x}-x-1\over2x+1}.                                 \tag{5.5}
\]

At `x=3/4`, (5.5) yields

\[
 H(3/4)\ge {7\over10}-{3\over25}e^{3/2}>0,                      \tag{5.6}
\]

where the last inequality follows, for example, by summing the exponential
series to fourth order and geometrically bounding its tail.  Hence `H>=0`
and the claimed bound follows.  Therefore

\[
                              C_P(q)\le {49\over25}=1.96<2.     \tag{5.7}
\]

There is also a useful explicit lower bound.  The absolute moments of `q`
are

\[
                         E|S|^k={(k+2)k!\over2^{k+1}}.           \tag{5.8}
\]

Minimizing the Rayleigh quotient over `span{x,x^3}` gives generalized
eigenvalue `(13-sqrt(15))/11`, and hence

\[
             C_P(q)\ge {11\over13-\sqrt{15}}>1.205.             \tag{5.9}
\]

For the `n`-fold product Laplace law, tensorization and a coordinate test
give the exact input constant two.  Normalized self-convolution acts
coordinatewise, so (5.7)--(5.9) apply in every dimension.  Product Laplace
therefore supports strict Poincare contraction; it is not a counterexample
to (2.7) or to a possible contraction of `C_P-1`.  It *is* a counterexample
to first-moment monotonicity by (4.5)--(4.7).

## 6. What remains for a renormalization proof

The strongest unconditional statement furnished by self-convolution is
(2.10).  A dimension-free strict contraction would follow from a stability
estimate of the schematic form

\[
 \operatorname{dist}_{L^2(\nu)}(f,\text{affine})^2
 \le C\,E(\Delta_4f)^2,                                       \tag{6.1}
\]

with the distance normalized compatibly with `E|grad f|^2`.  Exact vanishing
does imply affinity, but a universal quantitative constant in (6.1) has not
been proved.  Establishing it for all isotropic log-concave `mu` would be a
new additive-stability theorem; continuous spectral edges prevent replacing
it by an attainment argument.

Even such a contraction would have to be formulated relative to the
Gaussian fixed value:

\[
 C_P(\mathsf T\mu)-1\le q\,[C_P(\mu)-1],\qquad q<1,              \tag{6.2}
\]

because an absolute multiplicative contraction contradicts the Gaussian
fixed point.  This report neither proves nor disproves (6.2).  The exact
defect (2.10), the Gaussian equality case, and the product-Laplace check are
consistent with it.  On the other hand, fourth-moment contraction alone and
conditional symmetry do not imply it.

For `M1`, the Laplace and uniform calculations decisively rule out an
analogous monotone renormalization.  Any first-moment use of convolution
would need a different Lyapunov functional containing additional information
such as entropy, Fisher information, or the four-copy additive defect.  No
such functional is established here.
