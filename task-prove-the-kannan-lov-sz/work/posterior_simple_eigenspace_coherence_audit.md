# Audit of simple posterior eigenspace coherence and the longitudinal Hardy gate

## 0. Verdict

This note independently audits
`posterior_simple_eigenspace_coherence.md`.  All of the displayed local
identities and constants (0.3), (2.1)--(3.3), and (5.1)--(5.6) are correct
for smooth fields for which the integrations by parts are justified.  In
particular, the coefficient (3) in (5.1), the constant
(64/(3\kappa)) in (5.4), and the constant (8) in (5.6) are correct.

There is, however, one important correction to the interpretation of
(5.4): it has **not** reduced the twice--Bochner obstruction to the signed
longitudinal term (5.5).  Absorbing the transverse part of the cubic tensor
leaves

\[
 -\frac{64}{3\kappa}E[W^THW].                         \tag{0.1}
\]

For a normalized exact bottom field,

\[
 E\|DW\|_{HS}^{2}+E[W^THW]=\lambda,
 \qquad \|\mathcal A_1W\|_2^2=\lambda^2.             \tag{0.2}
\]

Thus (0.1) is (O_\kappa(\lambda)), not
(O_\kappa(\lambda^2)).  It is an additional unclosed remainder at the
scale relevant to a gap proof.  Formula (5.4) does prove that there is no
ambient-dimensional trace loss, but by itself it does not prove
second-order transverse coercivity.

A one-dimensional disintegration along the low eigenline also does not
close the argument from the stated hypotheses.  The exact obstruction is
twofold:

1. the conditional densities on eigenline flow curves need not be
   log-concave; their effective curvature contains two terms not controlled
   by the simple eigengap or by the first-derivative posterior saturation
   estimate;
2. even a uniform Hardy inequality on every flow curve controls only the
   within-curve oscillation.  A separate quotient estimate is needed for
   the conditional amplitude means.

Sections 4--6 state the precise extra global conditions.  They are a
uniform Muckenhoupt bound for the flow-conditionals and a transverse
control of their means.  With these two inputs the simple-eigenline branch
closes already at the first Bochner level.  Neither input is proved by the
uniform simple eigengap.

## 1. Independent derivation of the twice--Bochner identity

Let

\[
 L=\Delta-\nabla U\cdot\nabla,
 \qquad \mathcal A_1=-L+H,
 \qquad H=D^2U,
\]

and let (W=\nabla f), (C=DW=D^2f).  For a scalar smooth function
(g), the integrated Bochner identity is

\[
 E(Lg)^2=E\|D^2g\|_{HS}^2+E[\nabla g^TH\nabla g].   \tag{1.1}
\]

Applying (1.1) to all components of (W) gives

\[
 E|-LW|^2=E\|DC\|_{HS}^2+E\operatorname{tr}(CHC),  \tag{1.2}
\]

where symmetry of (C) is used in the last term.  Integration by parts
in the cross term gives, with

\[
 (\mathcal T_W)_{aj}=\sum_b(\partial_jH_{ab})W_b,
\]

\[
\begin{aligned}
 E\langle-LW,HW\rangle
 &=E\sum_{a,j}(\partial_jW_a)\partial_j(HW)_a\\
 &=E\operatorname{tr}(CHC)+E\langle C,\mathcal T_W\rangle_{HS}.
                                                               \tag{1.3}
\end{aligned}
\]

The tensor (D^3U) is fully symmetric, so both (C) and
(\mathcal T_W) are symmetric.  Expanding
(\|-LW+HW\|_2^2), equations (1.2)--(1.3) give exactly

\[
 \boxed{
 \|\mathcal A_1W\|_2^2
 =E\|DC\|_{HS}^2+3E\operatorname{tr}(CHC)
  +2E\langle C,\mathcal T_W\rangle_{HS}+E|HW|^2.}   \tag{1.4}
\]

This verifies (5.1), including the coefficient (3).

For a unit-Gaussian posterior, the audited Hilbert--Schmidt saturation
estimate is

\[
 \|\mathcal T_W\|_{HS}^2
 =\sum_j|(D_jH)W|^2\le16W^THW.                       \tag{1.5}
\]

No derivative of the pointwise selected vector (W(y)) is taken in
(1.5); the assertion for an adaptive vector is therefore legitimate.

## 2. Block calculation and exact constants in (5.3)--(5.6)

At a fixed point choose an orthonormal eigenbasis with first vector (e),
write (P=e e^T), and write

\[
 H=\begin{pmatrix}h&0\\0&H_Q\end{pmatrix},
 \qquad H_Q\succeq(h+\kappa)I,
 \qquad
 C=\begin{pmatrix}\alpha&b^T\\b&D\end{pmatrix}.
\]

Then

\[
 C_0=C-\alpha P=\begin{pmatrix}0&b^T\\b&D\end{pmatrix},
 \qquad \|C_0\|_{HS}^2=2|b|^2+\|D\|_{HS}^2,
\]

and direct expansion yields

\[
\begin{aligned}
 \operatorname{tr}(CHC)
 &=h\alpha^2+h|b|^2
   +\operatorname{tr}\!\left(H_Q(bb^T+D^2)\right)\\
 &\ge h\alpha^2+\kappa|b|^2+\kappa\|D\|_{HS}^2\\
 &\ge h\alpha^2+\frac\kappa2\|C_0\|_{HS}^2.        \tag{2.1}
\end{aligned}
\]

This proves (5.3).  Put

\[
 t=\langle P,\mathcal T_W\rangle_{HS}
   =e^T\mathcal T_We.
\]

Then

\[
 2\langle C,\mathcal T_W\rangle
 =2\alpha t+2\langle C_0,\mathcal T_W\rangle.       \tag{2.2}
\]

Young's inequality with parameter (3\kappa/4) says

\[
 2\langle C_0,\mathcal T_W\rangle
 \ge-\frac{3\kappa}{4}\|C_0\|_{HS}^2
     -\frac4{3\kappa}\|\mathcal T_W\|_{HS}^2.       \tag{2.3}
\]

Combining (1.4), (1.5), (2.1), and (2.3) gives exactly

\[
\begin{aligned}
 \|\mathcal A_1W\|_2^2
 &\ge E\|DC\|_{HS}^2+3E[h\alpha^2]
      +\frac{3\kappa}{4}E\|C_0\|_{HS}^2\\
 &\quad+2E[\alpha t]+E|HW|^2
      -\frac{64}{3\kappa}E[W^THW].                  \tag{2.4}
\end{aligned}
\]

Thus (5.4) is correct.  Cauchy--Schwarz and (1.5) also give

\[
 2|E\langle C,\mathcal T_W\rangle|
 \le2(E\|C\|^2)^{1/2}(E\|\mathcal T_W\|^2)^{1/2}
 \le8\sqrt{E\|DW\|^2\,E[W^THW]},                  \tag{2.5}
\]

which verifies (5.6).

The scale warning (0.1)--(0.2) follows immediately from (2.4).  In
particular, replacing the transverse cubic by (2.3) is dimension-free but
not perturbative at a small spectral edge.

## 3. What the longitudinal scalar actually is

Full symmetry of (D^3U) and the eigenvalue differentiation formula give

\[
\begin{aligned}
 t
 &=D^3U[e,e,W]
   =e^T(D_WH)e
   =D_Wh=W\cdot\nabla h.                             \tag{3.1}
\end{aligned}
\]

If (a=e\cdot W), then

\[
 D_ea=e^T(DW)e+(D_ee)\cdot W
      =\alpha+(D_ee)\cdot QW.                       \tag{3.2}
\]

Hence (alpha) is the derivative of the scalar amplitude only when the
line is fixed or (QW=0).  The simple-eigenspace estimate gives

\[
 |D_ee|\le |De|_{HS}\le\frac4\kappa\sqrt h,         \tag{3.3}
\]

so (3.2) has a controlled error in a first-energy ledger, but it is not an
exact scalar identity in the quantitative problem.

When (e) is fixed and (W=a(s)e), (s=e\cdot y), equations
(3.1)--(3.2) reduce (5.5) to

\[
 E[a'(s)a(s)h'(s)].                                  \tag{3.4}
\]

This confirms the stated one-dimensional model, but it does not furnish a
sign.  The pointwise posterior estimate gives only

\[
 |h'|\le4\sqrt h,                                    \tag{3.5}
\]

and Cauchy--Schwarz applied to (3.4) again loses a multiple of
(E a^2).  Controlling that multiple is a global Hardy/Poincare question.

## 4. Exact flow-box disintegration and the missing curvature terms

Because the base space is contractible, a smooth rank-one projection can
be oriented locally (and globally when no regularity obstruction is
present); choose a unit vector field (e).  In a flow box let

\[
 \partial_s\Phi(s,z)=e(\Phi(s,z)),
 \qquad dx=J(s,z)\,ds\,dz.
\]

Liouville's formula gives

\[
 \partial_s\log J=(\operatorname{div}e)\circ\Phi.   \tag{4.1}
\]

The conditional density of (d\nu=e^{-U}dx/Z) on a flow curve is
therefore proportional to

\[
 q_z(s)=\exp[-\Psi_z(s)],
 \qquad
 \Psi_z(s)=U(\Phi(s,z))-\log J(s,z).                 \tag{4.2}
\]

Differentiating twice gives the exact formula

\[
 \boxed{
 \Psi_z''
 =h+\nabla U\cdot D_ee-D_e(\operatorname{div}e).}   \tag{4.3}
\]

The eigengap and posterior saturation control (D_ee) through (3.3), but
they do not control the product (\nabla U\cdot D_ee): the score is
unbounded on an unbounded support.  They also do not control
(D_e\operatorname{div}e), which contains second derivatives of the
spectral projection and hence second derivatives of (H).  Consequently
(h\ge0) does **not** imply (\Psi_z''\ge0).  Ordinary one-dimensional
log-concave Hardy estimates cannot simply be applied to these flow
conditionals.

### A convex high-frequency stress test

The failure of that inference is quantitative.  On (\mathbb R^2), let

\[
 U_m(x,y)=\frac18x^2+\frac38y^2
       +\frac{\varepsilon}{m^3}\sin(mx)\sin(my),
 \qquad \varepsilon=10^{-2},\quad m\ge1.             \tag{4.4}
\]

Its Hessian is

\[
 H_m=
 \begin{pmatrix}
 \frac14-\frac\varepsilon m\sin(mx)\sin(my)&
 \frac\varepsilon m\cos(mx)\cos(my)\\
 \frac\varepsilon m\cos(mx)\cos(my)&
 \frac34-\frac\varepsilon m\sin(mx)\sin(my)
 \end{pmatrix}.                                      \tag{4.5}
\]

The perturbation in (4.5) has operator norm at most (2\varepsilon/m).
Hence

\[
 0.23I\preceq H_m\preceq0.77I,                      \tag{4.6}
\]

and the two eigenvalues are separated by at least (1/2), uniformly in
(m).  Moreover each of the two matrices (D_jH_m) has operator norm at
most (2\varepsilon), so

\[
 \sum_{j=1}^2|(D_jH_m)u|^2
 \le8\varepsilon^2|u|^2
 \le16u^TH_mu.                                      \tag{4.7}
\]

Thus this family satisfies the same type of eigengap and first-derivative
Hilbert--Schmidt estimate used in the local argument, with much room in
the constants.

Let (e_m=(\cos\theta_m,\sin\theta_m)) be the low eigenvector.  If
(c=(\varepsilon/m)\cos(mx)\cos(my)), then near (c=0)

\[
 \frac{d\theta_m}{dc}=-2,                            \tag{4.8}
\]

because the diagonal gap is (1/2).  At

\[
 p_m=\left(\frac\pi{2m},-\frac\pi{2m}\right)
\]

one has (c=0), (e_m=(1,0)), (D_{e_m}e_m=0), and

\[
 \partial_{xy}c=-\varepsilon m,
 \qquad
 D_{e_m}(\operatorname{div}e_m)
 =\partial_{xy}\theta_m=2\varepsilon m.             \tag{4.9}
\]

Since the low eigenvalue at (p_m) is (1/4+\varepsilon/m), (4.3)
becomes

\[
 \Psi_z''(p_m)=\frac14+\frac\varepsilon m-2\varepsilon m,
                                                               \tag{4.10}
\]

which tends to (-\infty).  The measures (e^{-U_m}dx/Z_m) are uniformly
strongly log-concave and have uniformly nondegenerate, bounded covariance,
so this is not a tail or normalization degeneration.  This example is not
claimed to be a Gaussian-output counterexample; it is a clean obstruction
to deriving log-concavity of the eigenline flow conditionals from the
local hypotheses (eigengap plus first-derivative saturation).  Additional
Gaussian-output structure would have to control the new term in (4.3), not
merely repeat (0.1).

## 5. The exact one-dimensional condition: Muckenhoupt constants

Fix a reference point (s_0=s_0(z)) on a flow curve and normalize (q_z)
to be a probability density on its interval (I_z).  Define

\[
\begin{aligned}
 B_{z,+}&=\sup_{s>s_0}q_z([s,\sup I_z))
                   \int_{s_0}^s\frac{dt}{q_z(t)},\\
 B_{z,-}&=\sup_{s<s_0}q_z((\inf I_z,s])
                   \int_s^{s_0}\frac{dt}{q_z(t)},\\
 B_z&=\max(B_{z,+},B_{z,-}).                         \tag{5.1}
\end{aligned}
\]

The elementary one-dimensional Hardy theorem gives, for every locally
absolutely continuous (g) with (g(s_0)=0),

\[
 \int_{I_z}g^2q_z\,ds
 \le4B_z\int_{I_z}(g')^2q_z\,ds.                    \tag{5.2}
\]

The scale (B_z) is also necessary: testing a function which equals
(\int_{s_0}^{\min(s,r)}q_z^{-1}) on the right side shows that any
constant in (5.2) is at least (B_{z,+}), and reflection gives
(B_{z,-}).  Therefore a uniform conditional Hardy argument requires the
genuinely global estimate

\[
 \operatorname*{ess\,sup}_z B_z\le C.               \tag{5.3}
\]

Neither the pointwise eigengap nor (0.1) estimates the tail integrals in
(5.1).  Formula (4.3) explains why the usual log-concave shortcut to
(5.3) is unavailable.

Even (5.3) is not sufficient by itself.  If

\[
 \bar a(z)=\int_{I_z}a(s,z)q_z(s)\,ds,
\]

then (5.2), with the standard passage from an anchor to the conditional
mean, controls (a-\bar a(z)), not (\bar a(z)).  Two flow curves with
constant amplitudes (+1) and (-1) are the canonical stress test:
their global mean can vanish while every longitudinal derivative is zero.
A second, quotient-level estimate controlling the conditional means by
transverse deformation is indispensable.

One explicit sufficient form for the rank-one section (Z=a e) is

\[
 \int |\bar a(z)|^2\,d\eta(z)
 \le C_Q\left(
       E\|C_0\|_{HS}^2+E[h|W|^2]+E|QW|^2+|EW|^2
             \right),                                \tag{5.4}
\]

where (d\eta) is the quotient measure of the disintegration.  The exact
right side can be varied, but some estimate of this type cannot be omitted.
It is a global connectivity assertion, not a consequence of a
one-dimensional Hardy inequality.

## 6. What would close the branch

Suppose one proves, for every normalized low spectral-window exact field,
the longitudinal one-form Hardy estimate

\[
 E|PW|^2
 \le C_\kappa\left(
       E\|DW\|_{HS}^2+E[W^THW]+|EW|^2
                    \right).                        \tag{6.1}
\]

Conditions (5.3)--(5.4), together with (3.2)--(3.3), are a concrete route
to (6.1).  The already-audited transverse estimates give

\[
 E|QW|^2\le\kappa^{-1}E[W^THW].                     \tag{6.2}
\]

For an exact normalized bottom field, testing the eigenfunction equation against the coordinate functions and using covariance Cauchy--Schwarz give

\[
 |EW|^2\le 2\lambda E|W|^2, \tag{6.3a}
\]

because \(EW=\lambda E[Yf]\), \(\operatorname{Cov}(Y)=2I\), and \(E|W|^2=\lambda E f^2\). For a spectral window one retains the corresponding mean-error estimate. The first Bochner identity gives

\[
 E\|DW\|^2+E[W^THW]=\lambda E|W|^2.                 \tag{6.3}
\]

Equations (6.1)--(6.3) imply

\[
 E|W|^2\le(3C_\kappa+\kappa^{-1})\lambda E|W|^2,
\]

and hence

\[
 \lambda\ge(3C_\kappa+\kappa^{-1})^{-1}.            \tag{6.4}
\]

Thus a uniform flow Muckenhoupt bound plus quotient-mean connectivity
would close the simple-eigenline branch without needing the twice--Bochner
scalar term.  Conversely, presenting (5.5) as closed by “one-dimensional
Hardy” without proving both (5.3) and a quotient estimate such as (5.4)
would hide the global step.

## 7. Canonical stress tests

1. **Fixed Gaussian line.**  For a constant (e) and quadratic (U),
   (D e=0), (t=0), and (4.3) reduces to
   (\Psi''=h).  The audit reproduces equality in the harmless model.

2. **Fixed non-Gaussian product line.**  If
   (U(s,z)=u(s)+v(z)), the flow Jacobian is one and
   (\Psi''=u''=h\).  A one-dimensional log-concave Poincare estimate,
   with the actual variance of (s), closes the longitudinal amplitude.
   This is precisely the fixed-direction rank-one-defect mechanism.

3. **High-frequency rotating line.**  The family (4.4) has uniform
   convexity, uniform simple eigengap, uniformly controlled (DH), and
   uniformly controlled covariance, but its eigenline conditional
   curvature has no lower bound.  It invalidates the inference
   “posterior saturation plus eigengap makes the flow fibers log-concave.”

4. **Disconnected quotient amplitudes.**  Constant amplitudes of opposite
   sign on two flow components have zero within-fiber derivative and can
   have zero global mean.  This shows why a conditional Hardy theorem alone
   cannot replace the missing quotient connectivity estimate.

The verified conclusion is therefore narrower than the source note's
headline: simple posterior eigenspaces have dimension-free local rotation
control, and the twice--Bochner algebra has no trace loss, but the current
ledger retains both an (O_\kappa(\mathcal K)) saturation remainder and a
global longitudinal/quotient Hardy gate.
