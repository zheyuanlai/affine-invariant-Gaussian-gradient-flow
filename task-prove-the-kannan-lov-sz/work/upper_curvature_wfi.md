# A pointwise WFI bound under an upper Hessian hypothesis

## 0. Statement and exact scope

This note records a valid positive subclass of the weighted-Fisher target.
The bounded-upper-Hessian subclass is now a valid quantitatively equivalent
target for KLS, but by a posterior argument distinct from the forward
normalized-convolution inequality.  More precisely, if

\[
 a=\lambda _1(\mu),\qquad
 b=\lambda _1\!\left(\mathcal L((X+G)/\sqrt2)\right),
\]

then posterior strong log-concavity gives

\[
 b\le\frac{2a}{1-a}\quad(a<1),\qquad
 a\ge\frac{b}{2+b};
 \tag{0.0}
\]

see `posterior_reverse_smoothing.md`.  The pointwise WFI estimate proved
below is not, by itself, a proof of a gap for that subclass: it supplies a
directional MMSE floor and closes only the near-linear gate.  A uniform
argument for the complementary nonlinear branch is still required.

Let \(X=(S,X_\perp)\) have a centered isotropic log-concave density
\(e^{-V}\) on \(\mathbb R^{1+d}\), and assume

\[
 0\preceq D^2V\preceq\kappa I
 \tag{0.1}
\]

in the smooth setting.  Let \(G_\perp\sim N(0,I_d)\) be independent, put
\(Z=X_\perp+G_\perp\), and write

\[
 r(s,z)=\operatorname {density}(S,Z)=e^{-U(s,z)},\qquad
 \rho(s)=\int r(s,z)\,dz,
\]

\[
 q_s(z)=r(s,z)/\rho(s),\qquad
 \ell_s(z)=\partial_s\log q_s(z).
\]

Then, for almost every \(s\),

\[
 \boxed{\quad
 I_\perp(s):=E_s\ell_s^2\le\kappa.
 \quad}
 \tag{0.2}
\]

Consequently, if \(\tau\) is the canonical Stein kernel of the centered
variance-one marginal \(S\),

\[
 \boxed{\quad
 \int\rho(s)\tau(s)^2I_\perp(s)\,ds
 \le400\kappa.
 \quad}
 \tag{0.3}
\]

The Stein--Fisher MMSE gate therefore gives, in every direction,

\[
 \boxed{\quad
 \mathbb E\operatorname {Var}(S\mid X+G)
 \ge\frac1{4(400+400\kappa)}
 =\frac1{1600(1+\kappa)}.
 \quad}
 \tag{0.4}
\]

In particular, a fixed unit-Gaussian regularization, rescaled to isotropic
coordinates, has \(0\preceq D^2V\preceq2I\), and hence satisfies (0.3)
with \(800\) and (0.4) with \(1/4800\).  By (0.0), proving a universal
spectral gap for all such regularized outputs would prove KLS for every
input.  The estimates (0.2)--(0.4) alone do not yet provide that spectral
gap.

## 1. The partial convolution preserves the directional upper bound

For fixed \((s,z)\), let \(\pi_{s,z}\) be the posterior law on \(x\) with
density proportional to

\[
 \exp\left[-V(s,x)-\frac12|z-x|^2\right].
\]

Differentiation under the integral gives

\[
 U_s=E_{\pi_{s,z}}V_s,
 \qquad
 U_{ss}=E_{\pi_{s,z}}V_{ss}
          -\operatorname {Var}_{\pi_{s,z}}(V_s).
 \tag{1.1}
\]

Prékopa gives \(D^2U\succeq0\), while (0.1) and (1.1) give

\[
 0\le U_{ss}\le E_{\pi_{s,z}}V_{ss}\le\kappa.
 \tag{1.2}
\]

Notice that no upper bound on the other blocks of \(D^2U\) is needed for
the next step.

## 2. Conditional Brascamp--Lieb and the Schur complement

Since

\[
 \ell_s=-U_s+\Phi'(s),\qquad \Phi=-\log\rho,
\]

we have

\[
 I_\perp(s)=\operatorname {Var}_{q_s}(U_s).
 \tag{2.1}
\]

Write the Hessian of \(U\) in the form

\[
 D^2U=
 \begin{pmatrix}
 U_{ss}&U_{sz}\\
 U_{zs}&H
 \end{pmatrix},
 \qquad H=D^2_{zz}U.
 \tag{2.2}
\]

First suppose \(H\succ0\).  Conditional Brascamp--Lieb applied to the
function \(z\mapsto U_s(s,z)\) gives

\[
 \operatorname {Var}_{q_s}(U_s)
 \le E_s\left[U_{sz}H^{-1}U_{zs}\right].
 \tag{2.3}
\]

The Schur complement of the positive-semidefinite matrix (2.2) gives,
pointwise,

\[
 U_{sz}H^{-1}U_{zs}\le U_{ss}\le\kappa.
 \tag{2.4}
\]

Equations (2.1)--(2.4) prove (0.2).

Here is a semidefinite passage which does not assume a uniform lower
curvature bound.  For \(\varepsilon>0\), use the conditional probability
with potential

\[
 U(s,z)+\frac\varepsilon2|z|^2.
\]

Apply Brascamp--Lieb to a bounded truncation \(T_R(U_s)\).  Its gradient
is \(T_R'(U_s)U_{zs}\), and block positivity implies

\[
 U_{sz}(H+\varepsilon I)^{-1}U_{zs}\le U_{ss}\le\kappa.
 \tag{2.5}
\]

Thus every truncated variance is at most \(\kappa\).  First send
\(\varepsilon\downarrow0\), using total-variation convergence of the
conditional probabilities for fixed \(R\).  Then use the pairwise formula

\[
 \operatorname {Var}(f)
 =\frac12\iint(f(z)-f(z'))^2q_s(z)q_s(z')\,dz\,dz'
\]

and Fatou as \(R\to\infty\).  This proves (0.2) without an invertibility
assumption on \(H\).  Standard convolution/confinement approximation
extends the argument to a distributional upper Hessian bound whenever
the latter is part of the hypotheses.

For completeness, here is the dimension-free Stein-kernel estimate used
in (0.3).  Let \(m\) be a median of the centered variance-one log-concave
density \(\rho\), and put \(h=\rho(m)\).  Cantelli's inequality gives
\(|m|\le1\).  On the half-line from \(m\) pointing away from a mode,
log-concavity gives \(\rho\le h\), and that half-line has mass \(1/2\).
Among nonnegative densities on a half-line bounded by \(h\) and having
mass \(1/2\), decreasing rearrangement minimizes the second moment about
the endpoint.  Therefore

\[
 \mathbb E(S-m)^2\ge \frac1{24h^2}.
\]

Since \(\mathbb E(S-m)^2=1+m^2\le2\), it follows that
\(h\ge1/\sqrt{48}>1/8\).  The right hazard
\(\rho/(1-F)\) is nondecreasing and the left reverse hazard \(\rho/F\)
is nonincreasing.  Both equal \(2h>1/4\) at the median.  Hence, for
\(s\ge m\), integration of the right-tail bound gives

\[
 \int_s^\infty t\rho(t)\,dt
 \le(1-F(s))(|s|+4),
 \qquad \frac{\rho(s)}{1-F(s)}\ge\frac14.
\]

For \(s\le m\), centering and the analogous left-tail calculation give

\[
 -\int_{-\infty}^s t\rho(t)\,dt
 \le F(s)(|s|+4),
 \qquad \frac{\rho(s)}{F(s)}\ge\frac14.
\]

The two numerators are equal to
\(N(s)=\rho(s)\tau(s)\).  Consequently

\[
 0\le\tau(s)\le4|s|+16,
 \qquad
 \mathbb E\tau(S)^2
 \le16\mathbb E(|S|+4)^2
 \le400,                                             \tag{2.6}
\]

where \(\mathbb E|S|\le1\) was used in the last step.  The argument also
covers finite support endpoints by taking the corresponding one-sided
hazard limits.  Combining (0.2) with (2.6) proves (0.3), and the already
audited Stein--Fisher MMSE inequality proves (0.4).

## 3. Stress tests

### 3.1 Products

If the slicing direction is a coordinate of a product law, then
\(q_s\) is independent of \(s\) and \(I_\perp(s)=0\).  For an oblique
direction, the conditional law generally changes with \(s\), but the
rotated Hessian still satisfies (0.1), so (2.2)--(2.5) give the same
dimension-free bound.  Tensor powers cannot accumulate a factor of the
dimension because the conclusion controls one scalar Hessian Schur
complement, not a trace.

### 3.2 The exponential cone

The transversely smoothed cone in gaussian_curvature_korn.md has a large
slice-local curvature/Korn ratio near its hard marginal endpoint.  It
does not satisfy a global bound on \(V_{ss}\), so it does not contradict
(0.2).  If the cone law is first convolved in **all** directions with a
unit Gaussian and rescaled by \(1/\sqrt2\), its potential satisfies
\(D^2V\preceq2I\).  A subsequent transverse channel therefore has
\(I_\perp(s)\le2\), including in the smoothed boundary layer.  Full
regularization removes exactly the unbounded \(s\)-curvature used by the
cone stress test.

### 3.3 High-Hermite Gaussian perturbations

Let

\[
 p_{\varepsilon,m}(x)\propto
 \exp\left[-\frac{x^2}{2}
 -\varepsilon\frac{H_m(x)}{\sqrt{m!}}\right],
 \qquad m\ge4\ \text{even},
\]

with \(\varepsilon>0\) small enough for convexity, and take products and
arbitrary rotated slicing directions.  The raw upper Hessian is not
uniform in \(x\), so (0.2) is not asserted before full regularization.
After convolution with a unit Gaussian and isotropic rescaling, the
potential has Hessian at most \(2I\), independently of \(m\), and (0.2)
gives \(I_\perp(s)\le2\).

Perturbatively this is consistent with the Hermite stress calculation:
Gaussian convolution attenuates a degree-\(m\) density perturbation by
\(2^{-m/2}\), while conditional Fisher information is quadratic in its
amplitude.  At leading order the non-Gaussian conditional Fisher
contribution is
\(O(m\varepsilon^2 2^{-m})+O_m(\varepsilon^3)\), not a quantity growing
like the inverse Hermite attenuation.  The earlier high-Hermite obstruction
concerns a reverse scalar Lyapunov budget and does not challenge the
pointwise Schur-complement bound.

## 4. What remains after the valid regularization reduction

For an arbitrary isotropic log-concave input, transverse convolution gives
\(0\preceq D^2_{zz}U\preceq I\), but it supplies no universal bound on
\(U_{ss}\).  Conditional Brascamp--Lieb then gives only

\[
 I_\perp(s)
 \le E_s\left[U_{sz}(D^2_{zz}U)^{-1}U_{zs}\right]
 \le E_sU_{ss},
 \tag{4.1}
\]

and the marginal identity is

\[
 \Phi''(s)=E_sU_{ss}-I_\perp(s).
 \tag{4.2}
\]

Thus (4.1) and (4.2) are tautological at the scale needed for general
WFI.  Replacing the missing \(U_{ss}\) bound by a dimension-free
integrated estimate remains the charge if one insists on proving the WFI
target directly for every unsmoothed input.

For spectral-gap purposes, full Gaussian regularization **can** be used as
a reverse reduction, but only through (0.0).  From the forward two-law
estimate

\[
 b(1-b)\le4(b-a)
\]

one obtains

\[
 a\le\frac{b(3+b)}4,
\]

not a lower bound for the original gap \(a\).  It remains correct but is not
the reverse mechanism.  The unresolved step is therefore sharper: prove a
universal gap for the analytic isotropic class
\(0\preceq D^2V\preceq2I\), for example by adding a dimension-free
nonlinear branch to (0.2)--(0.4).
