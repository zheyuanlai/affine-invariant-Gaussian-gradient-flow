# Independent audit of the posterior Hilbert--Schmidt third-moment theorem

## 0. Verdict

The Hilbert--Schmidt estimate in
`posterior_covariance_saturation_stability.md` is valid.  Precisely, if
\(Z\) has a centered \(1\)-strongly log-concave law \(\pi\) on its
intrinsic \(k\)-dimensional affine support, \(A=E[ZZ^T]\), and

\[
 \delta_\pi(u)=|u|^2-u^TAu,
\]

then

\[
 \left\|E[(u\cdot Z)ZZ^T]\right\|_{HS}
 \le4\sqrt{\delta_\pi(u)}.                           \tag{A.1}
\]

The constant \(4\) is justified.  No trace, rank, or dimension factor is
hidden in the proof.  The only approximation bookkeeping needed is to
center the smooth approximants and pass their moments using the uniform
subgaussian bounds supplied by Caffarelli contraction.

## 1. The quadratic-variance input

For every symmetric \(B\), the quadratic polynomial

\[
 q_B(z)=z^TBz-E[Z^TBZ]
\]

obeys

\[
 \boxed{
 E q_B(Z)^2
 \le4E|BZ|^2
 =4\operatorname {tr}(B^2A)
 \le4\|B\|_{HS}^2.}                                 \tag{A.2}
\]

This is an ordinary consequence of the Poincare constant one for a
\(1\)-strongly log-concave law: \(\nabla q_B=2Bz\).  It does not invoke
KLS.  For an extended convex potential or a hard convex support, the same
inequality follows directly by writing \(Z=T(G)\), where Caffarelli gives
a \(1\)-Lipschitz map from the standard Gaussian, and applying Gaussian
Poincare to \(q_B\circ T\).  Quadratic polynomials belong to the closed
form domain because this contraction gives all finite moments; truncation
then supplies a fully literal form-domain proof.

The numerical constant in (A.2) is \(4\), so
\(\|q_B\|_2\le2\|B\|_{HS}\).  This square-root constant \(2\) is the
first half of the final constant \(4\) in (A.1).

## 2. Caffarelli deficit identities

Let \(T=\nabla\varphi\) push \(G\sim N(0,I_k)\) to \(Z\), and write
\(J=DT\).  Caffarelli gives \(0\preceq J\preceq I\) almost everywhere.
Gaussian Poincare applied to \(u\cdot T\), followed by
\((I-J)^2\preceq I-J^2\), gives

\[
 E|(I-J)u|^2\le\delta_\pi(u).                       \tag{A.3}
\]

For

\[
 r=u\cdot(T(G)-G),\qquad e=(I-J)u,
\]

one has \(Er=0\), \(\nabla r=(J-I)u\), and a second Gaussian Poincare
inequality gives

\[
 Er^2\le\delta_\pi(u),\qquad E|e|^2\le\delta_\pi(u). \tag{A.4}
\]

The direction in the intermediate inequality is important and is correct:
\(\operatorname {Var}(u\cdot T)\le E|Ju|^2\) implies
\(E(|u|^2-|Ju|^2)\le |u|^2-u^TAu\).

## 3. Hilbert--Schmidt duality and matrix integration by parts

Since \(E[rZZ^T]\) is symmetric, Hilbert--Schmidt duality may be restricted
to symmetric \(B\).  Equations (A.2) and (A.4) give

\[
 \begin{aligned}
 \|E[rZZ^T]\|_{HS}
 &=\sup_{B=B^T,\,\|B\|_{HS}=1}|E[rq_B(Z)]|\\
 &\le2\sqrt{\delta_\pi(u)}.                         \tag{A.5}
 \end{aligned}
\]

For arbitrary \(B\) of Hilbert--Schmidt norm one,

\[
 |E[e^TBZ]|
 \le(E|e|^2)^{1/2}
      \{\operatorname {tr}(B^TBA)\}^{1/2}
 \le\sqrt{\delta_\pi(u)}.
\]

Therefore

\[
 \|E[eZ^T]\|_{HS}\le\sqrt{\delta_\pi(u)}.         \tag{A.6}
\]

The Sobolev Gaussian integration-by-parts identity is

\[
 E[(u\cdot G)ZZ^T]
 =E[(Ju)Z^T+Z(Ju)^T].                               \tag{A.7}
\]

It is valid because \(T\) is Lipschitz, \(J\) is bounded, and the two
sides are integrable; it follows first after a cutoff and then by Gaussian
dominated convergence.  Since \(Ju=u-e\), \(EZ=0\), and
\(u\cdot Z=u\cdot G+r\), (A.7) becomes

\[
 E[(u\cdot Z)ZZ^T]
 =E[rZZ^T]-E[eZ^T+Ze^T].                            \tag{A.8}
\]

Combining (A.5), (A.6), and invariance of the Hilbert--Schmidt norm under
transpose gives \(2\sqrt\delta+2\sqrt\delta=4\sqrt\delta\), proving
(A.1).

## 4. Posterior derivative and scaling checks

For the unit-Gaussian posterior, differentiation of the locally dominated
exponential family gives

\[
 D A_y[v]=E[Z_yZ_y^T(v\cdot Z_y)].                  \tag{A.9}
\]

The linear map \(v\mapsto(D A_y[v])u\) has matrix
\(M_{\pi_y}(u)=E[(u\cdot Z_y)Z_yZ_y^T]\).  Hence, for every orthonormal
basis \((e_j)\),

\[
 \sum_j\|(D A_y[e_j])u\|^2
 =\|M_{\pi_y}(u)\|_{HS}^2.                         \tag{A.10}
\]

Thus the Hilbert--Schmidt posterior derivative bound is exactly (A.1),
not merely a consequence of entrywise scalar estimates.  Since
\(D^2U=I-A_y\), the corresponding constants are \(4\) in the unscaled
third-derivative estimate, \(2|u|\) for the Lipschitz constant of
\(\sqrt{u^TD^2Uu}\), and \(8\) after the isotropic rescaling
\(s\mapsto\sqrt2s\).  These constants in the revised note are correctly
scaled.

## 5. Approximation checklist and scope

* For nonsmooth or hard-support targets, either use the generalized
  Caffarelli contraction directly or approximate the extended convex
  potential by smooth strongly convex potentials with the same curvature
  lower bound.
* Translate every approximant to its barycenter.  Translation preserves
  strong convexity and the contraction derivative bound.
* The contraction gives uniform one-dimensional subgaussian moments.
  In each fixed intrinsic dimension these imply uniform integrability of
  the covariance and third-moment entries, so (A.1) passes to the limit
  with the same constant.
* Proper affine supports must be treated intrinsically; ambient null
  directions are not part of the Hilbert--Schmidt norm.
* Posterior differentiation is legitimate for extended convex priors:
  on compact \(y\)-sets the Gaussian likelihood dominates every
  polynomial derivative by an integrable Gaussian-tail envelope.

The theorem is therefore formalization-ready as a local posterior
stability statement.  It removes the earlier trace loss in differentiating
an adaptively chosen near-saturation direction, but it does not by itself
prove the global one-form coercivity inequality; no such conclusion is
used in this audit.
