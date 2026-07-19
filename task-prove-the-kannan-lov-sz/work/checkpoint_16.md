# Checkpoint 16: reverse smoothing repaired and the regularized nonlinear gate

## 0. Status

This checkpoint does not prove KLS.  It corrects a sign error in the first
version of the two-law Gaussian reduction, supplies a different valid
posterior proof of the reduction, and isolates a single affine-orthogonal
interaction inequality which would close the remaining nonlinear branch.
That interaction inequality is not proved here.

The valid new general theorem is the posterior reverse-smoothing estimate

\[
 \lambda_1(\mathcal G\mu)
 \le \frac{2\lambda_1(\mu)}{1-\lambda_1(\mu)},
 \qquad
 \lambda_1(\mu)
 \ge\frac{\lambda_1(\mathcal G\mu)}
 {2+\lambda_1(\mathcal G\mu)}.                         \tag{RS}
\]

It has passed an independent audit covering the form domain, hard convex
supports, intrinsic lower-dimensional supports, and the \(1/\sqrt2\)
rescaling.

## 1. Mandatory sign correction

For isotropic input gaps with common lower bound \(a\), the unequal-law
forward convolution theorem remains

\[
 b(1-b)\le4(b-a),
 \qquad
 b\ge\frac{\sqrt{9+16a}-3}{2}.                         \tag{1.1}
\]

The correct rearrangement is

\[
 \boxed{a\le\frac{b(3+b)}4,}                           \tag{1.2}
\]

not the reverse inequality.  Thus (1.1) alone does not transfer a lower
bound on the smoothed gap back to the input.  The erroneous reverse claim
is retracted in
`fixed_gaussian_reduction_retraction.md`.

The forward Hoeffding proof, its continuous-spectral-edge passage, and the
lower amplification bound in (1.1) are unaffected.

## 2. The valid posterior reverse theorem

Let \(X\sim\mu\), \(G\sim N(0,I)\) be independent, \(Y=X+G\), and

\[
 F(y)=E[f(X)\mid Y=y]
\]

for a centered normalized input test with energy
\(q=E|\nabla f|^2<1\).  The posterior law is

\[
 d\mu_y(x)\propto
 \exp\left[-V(x)-\frac12|x-y|^2\right]dx.
\]

It is intrinsically \(1\)-strongly log-concave, so

\[
 A_y:=\operatorname{Cov}_{\mu_y}(X)\preceq I,
 \qquad
 \operatorname{Var}_{\mu_y}(f)
 \le E_{\mu_y}|\nabla f|^2.                            \tag{2.1}
\]

Total variance and disintegration give

\[
 \operatorname{Var}(F(Y))
 =1-E\operatorname{Var}(f(X)\mid Y)
 \ge1-q.                                               \tag{2.2}
\]

Posterior differentiation gives

\[
 \nabla F(y)=\operatorname{Cov}_{\mu_y}(X,f).
\]

For every unit \(u\), covariance Cauchy--Schwarz and (2.1) imply

\[
 |u\cdot\nabla F|^2
 \le\operatorname{Var}_{\mu_y}(u\cdot X)
       \operatorname{Var}_{\mu_y}(f)
 \le\operatorname{Var}_{\mu_y}(f).
\]

Taking the supremum in \(u\) and averaging yields

\[
 \int|\nabla F|^2d\mathcal L(Y)
 \le E\operatorname{Var}(f(X)\mid Y)
 \le q.                                                \tag{2.3}
\]

For \(\widetilde F(s)=F(\sqrt2s)\) on
\(S=(X+G)/\sqrt2\), the variance is at least \(1-q\) and the energy is
at most \(2q\).  Taking a form-domain minimizing sequence proves (RS)
without assuming an attained eigenfunction.

The distinction from the failed self-convolution projection is exact:
(2.3) uses the posterior unit-curvature Poincare inequality furnished by
the Gaussian likelihood.

## 3. Equivalent upper-curvature output class

If \(U=-\log(d\mathcal L(Y)/dy)\), the posterior formula gives

\[
 D^2U(y)=I-A_y.
\]

Thus the potential \(\widetilde U(s)=U(\sqrt2s)+\mathrm{const}\) of the
isotropic output satisfies

\[
 \boxed{0\preceq D^2\widetilde U\preceq2I.}           \tag{3.1}
\]

The density is positive and analytic.  By (RS), a universal gap for all
such fixed-Gaussian outputs implies a universal gap for every isotropic
log-concave input, with transfer \(c_0\mapsto c_0/(2+c_0)\).

Two dimension-free statements are now proved on this output class.

1. **Score rigidity.**  If \(0\preceq D^2V\preceq\beta I\), then with
   \(Y=\nabla V(X)\),
   \[
   E[(Y-X)(Y-X)^T]\preceq(\beta-1)I.
   \]
   Hence every affine-orthogonal \(g\) obeys
   \[
   |E\nabla g|\le\sqrt{\beta-1}\,\|g\|_2.             \tag{3.2}
   \]
   At \(\beta=2\), a bottom eigenfunction, or bottom spectral-window
   sequence, whose affine residual is at most \(1/2\) has gap at least
   \[
   1-1/\sqrt3.                                         \tag{3.3}
   \]

2. **Upper-curvature WFI.**  For any scalar/transverse splitting, the
   transversely smoothed conditional Fisher information obeys pointwise
   \[
   I_\perp(s)\le\beta.
   \]
   Therefore the Stein-weighted Fisher integral is at most
   \(400\beta\), and the audited Stein--Fisher estimate gives a
   directional unit-noise MMSE floor
   \(1/[1600(1+\beta)]\).  This is another closure of the near-affine
   posterior branch; it does not control genuinely nonlinear low modes.

## 4. The affine-orthogonal ANOVA target

Let \(S=(X+G)/\sqrt2\), let \(f\in L^2_0(\mathcal L(S))\), and put

\[
 R=f(S)-E[f(S)\mid X]-E[f(S)\mid G].                  \tag{4.1}
\]

The exact candidate is

\[
 \boxed{
 \operatorname{dist}_{L^2(\mathcal L(S))}
 (f,\mathrm{Aff})^2\le C_{\rm A}E R^2.}               \tag{AI}
\]

The operator audit proves

\[
 E R^2=\langle f,(I-Q_X-Q_G)f\rangle,                 \tag{4.2}
\]

and, for a full-dimensional log-concave input, the kernel of this positive
operator is exactly the affine space.  In the Gaussian case its spectrum
on Hermite degree \(m\) is \(1-2^{1-m}\), so \(C_{\rm A}=2\) is sharp.
For a one-sided exponential coordinate quadratic one needs
\(C_{\rm A}\ge3\).  No full-dimensional counterexample is known from the
canonical product, body, radial, or polynomial tests.

The lower-dimensional formulation must be intrinsic.  If ambient Gaussian
noise is retained in a null direction, a nonlinear function of that pure
Gaussian coordinate has \(R=0\), contradicting (AI).

The quadratic sector is already substantial.  For symmetric \(A\),
\(f_A(s)=s^TAs-\operatorname{tr}A\) satisfies

\[
 ER^2=\operatorname{tr}(A^2),                          \tag{4.3}
\]

while its affine-orthogonal variance contains
\(\frac14\operatorname{Var}(X^TAX)\).  Thus (AI) would imply a
dimension-free generalized quadratic-variance theorem.  It may not be
inserted as an elementary functional-equation lemma.

## 5. Why (AI) would close the regularized nonlinear branch

Let \(a=\lambda_1(\mu)\), \(b=\lambda_1(\mathcal G\mu)\), and choose a
centered normalized bottom spectral-window vector \(f\) for the output,
with energy \(q\le b+\varepsilon\).  Its unequal Hoeffding decomposition
has squared pieces \(A,B,r\), where \(r=ER^2\) and

\[
 A+B+r=1.
\]

Using the input gap \(a\), the Gaussian gap one, and the coordinatewise
Poincare inequalities for the residual gives

\[
 q\ge aA+B+(a+1)r
   =a+(1-a)B+r.                                       \tag{5.1}
\]

Therefore \(r\le q-a\).  The reverse theorem (RS) yields

\[
 a\ge\frac b{2+b},
 \qquad
 \limsup_{\varepsilon\downarrow0}r
 \le b-\frac b{2+b}
 =\frac{b(1+b)}{2+b}\le\frac23b.                     \tag{5.2}
\]

If (AI) held, the affine residual \(\delta\) of every bottom spectral
vector would satisfy

\[
 \delta^2\le\frac{2C_{\rm A}}3b+o(1).                \tag{5.3}
\]

For \(b<3/(8C_{\rm A})\), (5.3) puts the residual below \(1/2\), whereas
score rigidity forces \(b\ge1-1/\sqrt3\).  Hence

\[
 b\ge
 \min\left\{1-\frac1{\sqrt3},\frac3{8C_{\rm A}}\right\},
\]

and (RS) transfers this constant to the original input.  This proves that
(AI), with a universal constant and its full form-domain scope, is a
complete KLS target.

No proof of (AI) is supplied here.  It is the exact nonlinear gate left by
the new smoothing and score-rigidity reductions.

## 6. Other verified Checkpoint 16 progress

The clean-room audit of the scalar residual curvature--Hardy estimate
proves, in the smooth local setting,

\[
 B_s^0\le128\sqrt{\pi/2}\,\sigma_s C_s.
\]

If \(\sigma_s\le L(1+|s|)\), the integrated residual is below
\(6600L\).  A BV/Stieltjes proof also gives the nonsmooth marginal charge

\[
 \int(1+|s|)\tau(s)^2\rho(s)\,d(D^2\Phi)(s)\le41.
\]

This closes uniform-planar and scalar-Gaussian scale-growth subclasses.
For arbitrary slices, the scale-growth/direct weighted charge and the
nonsmooth convergence of conditional Poisson fields remain unproved.

## 7. Remaining gates

At this checkpoint a complete proof can proceed by either of two genuinely
new statements:

1. prove (AI), or a version restricted to bottom spectral windows, for
   every intrinsic fixed-Gaussian output; or
2. prove the general integrated nonlinear residual WFI estimate and a
   global mechanism transferring it from the near-affine branch to every
   low mode.

The forward convolution inequality, upper curvature, posterior MMSE,
thin shell, and the verified Gaussian/pure-translation slice subclasses
do not by themselves supply either statement.
