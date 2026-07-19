# Angular stability for radial contractions and arbitrary conic sets

## 1. Statement

Let `G` be standard Gaussian in `R^n`, `n >= 2`.  Let

\[
 T(x)=a(|x|)x
\]

be a radial Brenier map such that, almost everywhere,

\[
 0\preceq \nabla T\preceq L I.                         \tag{1.1}
\]

Thus `a(r) >= 0`, and the tangential eigenvalue bound in (1.1) gives
`a(r) <= L` for almost every `r`.  Put `X=T(G)`.  Radial symmetry gives
`E X=0`.

Let `B` be an arbitrary conic Borel set: outside the origin,
membership in `B` depends only on `x/|x|`.  Put

\[
 h=\mathbf 1_B-g,\qquad g=\gamma_n(B),\qquad
 v=\mathbb E[hX],\qquad D=\mathbb E[hXX^T].          \tag{1.2}
\]

If `v ne 0`, set `u=v/|v|` and `P=I-uu^T`.

**Theorem 1 (radial-conic angular stability).**  For every
`delta in (0,1/2)` there is a finite `C_delta`, independent of `n`,
`a`, `L`, and `B`, such that, whenever

\[
 \delta\le g\le1-\delta,
\]

one has

\[
 \boxed{
  \|PD\|_{HS}^2
  \le C_\delta L^3\{L\mathcal I(g)-|v|\}.}           \tag{1.3}
\]

This is a genuine simultaneous nonlinear-map/non-halfspace case of the
proposed static contraction inequality.  The radial factor `a` may be
nonconstant, and the angular set may be any Borel subset of the sphere;
it need not be a hemisphere or convex.

For `n=1`, every nontrivial conic set is a half-line and `P=0`, so the
same conclusion is immediate.

## 2. Exact polar factorization

Write

\[
 G=R\Theta,
\]

where `R=|G|` and `Theta=G/|G|` are independent and `Theta` is uniform
on the unit sphere.  Define the angular moments

\[
 m=\mathbb E_\Theta[h(\Theta)\Theta],\qquad
 Q=\mathbb E_\Theta[h(\Theta)\Theta\Theta^T],        \tag{2.1}
\]

and the radial moments

\[
 r_1=\mathbb E R,\qquad
 b(R)=a(R)R,\qquad
 \beta_1=\mathbb E b(R),\qquad
 \beta_2=\mathbb E b(R)^2.                          \tag{2.2}
\]

Polar independence gives the four exact identities

\[
 \begin{aligned}
  \mathbb E[hG]&=r_1m,
  &\mathbb E[hGG^T]&=\mathbb E[R^2]Q=nQ,\\
  v&=\beta_1m,
  &D&=\beta_2Q.                                      \tag{2.3}
 \end{aligned}
\]

Because `a >= 0`, the direction of `v` is the direction of `m` whenever
`v ne 0`.  Hence the projection `P` in (1.2) is also the transverse
projection associated with the Gaussian centroid `r_1m`.

The pointwise bound `0 <= a <= L` yields

\[
 0\le\beta_1\le Lr_1,
 \qquad
 0\le\beta_2\le L^2\mathbb E R^2=L^2n.             \tag{2.4}
\]

No derivative of `a` is needed after (1.1); in particular, no smoothness
or strict monotonicity is hidden in the argument.

## 3. Reduction to the exact Gaussian angular theorem

Apply the Gaussian angular-stability theorem to the same source set `B`.
In the notation of (2.3), it says

\[
 \|P(nQ)\|_{HS}^2
 \le C_\delta\{\mathcal I(g)-r_1|m|\}.              \tag{3.1}
\]

Using (2.3), (2.4), and then `beta_1/L <= r_1`, we obtain

\[
 \begin{aligned}
 \|PD\|_{HS}^2
 &= {\beta_2^2\over n^2}\|P(nQ)\|_{HS}^2\\
 &\le C_\delta L^4
       \{\mathcal I(g)-r_1|m|\}\\
 &\le C_\delta L^4
       \left\{\mathcal I(g)-{\beta_1\over L}|m|\right\}\\
 &=C_\delta L^3\{L\mathcal I(g)-|v|\}.
                                                               \tag{3.2}
 \end{aligned}
\]

This proves (1.3).

The proof also identifies why radial nonlinearity cannot produce a
counterexample.  It changes the angular first and second correlated
moments only through the two scalars `beta_1` and `beta_2`.  The first is
at most `L` times its Gaussian value, while the square of the second is
at most `L^4` times its Gaussian value.  These are exactly the powers
required by the homogeneous inequality.

## 4. A dimension-free coarse bound for every contraction

For comparison, there is a useful estimate with no radial or conic
hypothesis.  Let `T` be any centered `1`-Lipschitz image of the Gaussian
(the gradient assumption is unnecessary here), and put `X=T(G)-ET(G)`.
For every symmetric `A` with `||A||_HS=1`, Gaussian Poincare gives

\[
 \begin{aligned}
 \operatorname {Var}(X^TAX)
 &\le \mathbb E|2(\nabla T)^TAX|^2\\
 &\le4\mathbb E|AX|^2
 =4\operatorname {tr}(A^2\operatorname {Cov}X)
 \le4.                                                   \tag{4.1}
 \end{aligned}
\]

The last inequality uses `Cov(X) preceq I`, itself an immediate scalar
Gaussian-Poincare consequence of `||nabla T||op <= 1`.  Since `Eh=0`,
Hilbert--Schmidt duality and Cauchy--Schwarz imply

\[
 \boxed{
  \|D\|_{HS}^2\le4g(1-g),\qquad
  \|PD\|_{HS}^2\le4g(1-g).}                         \tag{4.2}
\]

After scaling, the right side is `4L^4g(1-g)`.  Thus any failure of the
static angular-stability inequality must occur in the near-equality
regime `L I(g)-|v| -> 0`; it cannot be a raw high-dimensional quadratic
moment blowup.

## 5. Counterexample searches this theorem rules out

The exact factorization covers radial clipping, radial shrinkage with a
sharp transition, projections onto centered balls, and smooth radial
Moreau contractions, together with arbitrary angular branching of the
set.  In particular:

1. concentrating the change of radial scale on a thin Gaussian shell
   cannot amplify the Hilbert--Schmidt term;
2. allowing the angular set to have arbitrarily many components does not
   create a dimensional loss; and
3. radial maps cannot exploit high-dimensional shell concentration to
   make the first-moment deficit smaller than the angular second-moment
   defect.

The remaining intersection must therefore use a genuinely nonradial
nonlinear map and a pullback set whose Gaussian flux changes direction or
branches relative to that map.
