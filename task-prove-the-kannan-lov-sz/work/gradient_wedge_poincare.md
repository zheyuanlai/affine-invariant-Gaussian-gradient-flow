# A wedge-Poincare lower bound for high-rank gradient phases

## 1. Exact Hilbert-space inequality

Let `nu` be a probability measure with Poincare constant at most `C`, and
let `F:R^d->R^m` belong to the first-order Sobolev space. Put

\[
 m=E_\nu F,\qquad R=E_\nu[FF^T],\qquad H=\nabla F.
\]

Then

\[
\boxed{
 (trR)^2-tr(R^2)-\{|m|^2trR-m^TRm\}
 \le C E\,tr\{H^T(trR\,I-R)H\}.}                    \tag{1.1}
\]

In particular, if `||R||op<=kappa trR` with `kappa<1/2`, then

\[
\boxed{
 E||H||_{HS}^2
 \ge {1-2\kappa\over C}\,trR.}                      \tag{1.2}
\]

**Proof.** Let `X,Y` be independent with law `nu`. For fixed
`a=F(Y)`, apply Poincare componentwise to the exterior-product-valued map
`x mapsto a wedge F(x)`:

\[
 E_X|a\wedge(F(X)-m)|^2
 \le C E_X\sum_j|a\wedge H(X)e_j|^2.                \tag{1.3}
\]

Average in `Y`. The left side is

\[
 (trR)^2-tr(R^2)-\{|m|^2trR-m^TRm\}.                \tag{1.4}
\]

For every vector `b`,

\[
 E_Y|F(Y)\wedge b|^2=trR|b|^2-b^TRb,
\]

so the right side of (1.3) becomes the right side of (1.1). Since
`mm^T preceq R`, one has `|m|^2<=||R||op`. Therefore the left side of
(1.1) is at least

\[
 (trR)^2-||R||op\,trR-|m|^2trR
 \ge(1-2\kappa)(trR)^2.
\]

The right side is at most `C trR E||H||HS^2`, proving (1.2).

The proof uses neither symmetry of `H` nor that `F` is a gradient. The
gradient case is important only because the heat-flow identities provide a
separate upper budget for `H`.

## 2. Gaussian observation specialization

Let `X~mu`, where `mu` is isotropic and log-concave, let
`C_P(mu)=K`, let `Y=X+sqrt(s)G`, and let `q_s` be the law of `Y`.
Tensorization on `(X,G)` gives

\[
                         C_P(q_s)\le K+s.             \tag{2.1}
\]

For a binary label `B`, put

\[
 g(y)=P(B=1\mid Y=y),\qquad
 h(y)=2\arcsin\sqrt{g(y)}.
\]

The posterior differentiation formula gives

\[
 \nabla h={\nabla g\over\sqrt{g(1-g)}}
 ={Cov(B,X\mid Y)\over s\sqrt{g(1-g)}}.              \tag{2.2}
\]

Consequently, for `F=sqrt(s) grad h`, its second-moment matrix is exactly
the binary Fisher matrix in heat coordinates:

\[
 R=E_{q_s}[FF^T]
 ={1\over s}E\left[{vv^T\over g(1-g)}\right].        \tag{2.3}
\]

If `mu` is isotropic, conditional covariance and total covariance give

\[
                         R\preceq s^{-1}I.            \tag{2.4}
\]

Thus `trR>=r_0` at a scale `s asymp K` implies effective rank of order
`K`. More precisely, (1.2) applies once `sr_0>2`; for bounded `K` there is
nothing to prove in a KLS contradiction argument. Applying (1.2) in the
large-`K` regime,

\[
\boxed{
 E_{q_s}||\nabla^2h||_{HS}^2
 \ge {c r_0\over s(K+s)}\asymp {c r_0\over K^2}.}    \tag{2.5}


The equivalent natural-tilt statement follows from `C=s^{-1}Y` and
`C_P(Law(C))<=s^{-2}(K+s)`.

## 3. Exact relation to the heat eikonal variable

Writing `z=Phi^{-1}(g)`, define

\[
 a(z)={I(g)\over\sqrt{g(1-g)}}
 ={\varphi(z)\over\sqrt{\Phi(z)(1-\Phi(z))}}.
\]

Then

\[
 \nabla h=a(z)\nabla z,
\qquad
 \nabla^2h=a(z)\nabla^2z+a'(z)\nabla z\nabla z^T.    \tag{3.1}

The first term contains the angular Hessian controlled by the Bernstein
identity. The second is a purely longitudinal amplitude term. The scalar
functions `a` and `a'` are bounded, while log-concavity makes the posterior
`s^{-1}`-strongly log-concave and the sharp centroid bound gives
`s|grad z|^2<=1`; hence

\[
 sE||a'(z)\nabla z\nabla z^T||_{HS}^2\le {C\over s}. \tag{3.2}

At `s asymp K`, (3.2) has exactly the same `K^{-2}` order as the
multiplicity lower bound (2.5). Therefore (1.1) and the displayed upper
bound alone do not force angular energy: their orders do not exclude the
possibility that longitudinal amplitude wells pay the Poincare cost.

This is the differential form of the one-bit/phase-cell obstruction.
Gaussian halfspaces have a single direction and zero wedge term. Gaussian
parity and polyhedral phase cells change directions through regions where
the Fisher amplitude is small; the longitudinal term in (3.1) records that
transition. A closing theorem must either subtract this scalar amplitude
cost using extremality or show that many amplitude wells of a posterior of
one near-minimal cut have additional geometric cost.

## 4. Status

Equation (1.1) is an exact new rank-sensitive constraint and (2.5) has the
correct dimension-free normalization. It is not a KLS proof. Its lower
bound matches, rather than exceeds, the scalar Bernstein and mutual-
information budgets. A simultaneous near-equality inverse theorem for the
convolution Poincare proof could still be useful, but assuming that it
produces a Gaussian factor would be an unproved spectral-rigidity input.

All differentiations may first be made after Gaussian/convex smoothing and
with bounded truncations. If `Hess h` is not in `L^2`, the right side of
(1.1) is infinite and the inequality is automatic; otherwise closedness of
the weak gradient and lower semicontinuity pass the formulas to the limit.
