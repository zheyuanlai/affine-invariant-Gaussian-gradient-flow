# Weighted Fisher information and reinforced Prékopa

## 0. Verdict

This note checks the weighted-Fisher route in the regular setting and
isolates its exact remaining lemma.  Let

\[
r(s,z)=e^{-U(s,z)},\qquad
\rho(s)=\int r(s,z)\,dz=e^{-\Phi(s)},\qquad
q_s(z)=r(s,z)/\rho(s).
\]

Write

\[
L_s=\Delta_z-\nabla_zU\cdot\nabla_z,\quad
\ell_s=\partial_s\log q_s=-U_s+\Phi',\quad
L_sg_s=\ell_s,\quad F_s=\nabla_zg_s,\quad H_s=D^2_{zz}U.
\]

Put

\[
C_s=\mathbb E_s\|D_zF_s\|_{\rm HS}^2,\qquad
B_s=\mathbb E_s\langle H_sF_s,F_s\rangle.
\]

The conditional Bochner identity gives

\[
I_\perp(s):=\mathbb E_s\ell_s^2=C_s+B_s. \tag{0.1}
\]

Define

\[
\alpha_s=\partial_sg_s-\frac12|F_s|^2,\qquad
\mathcal Q_s=D^2U[(1,-F_s),(1,-F_s)]\ge0.
\]

The two main exact identities are

\[
\boxed{L_s\alpha_s=\Phi''-\mathcal Q_s-\|D_zF_s\|_{\rm HS}^2}, \tag{0.2}
\]

\[
\boxed{\Phi''(s)=\mathbb E_s\mathcal Q_s+C_s}. \tag{0.3}
\]

If \(\tau\) is the canonical Stein kernel of a centered variance-one
marginal \(\rho\), then

\[
\int\rho\tau^2(C_s+\mathbb E_s\mathcal Q_s)\,ds
=\int\rho\tau^2\Phi''\,ds
=1-\mathbb E_\rho(\tau')^2\le1. \tag{0.4}
\]

Thus the deformation part \(C_s\) is bounded dimension-freely.  The exact
unresolved term is

\[
\mathcal B=\int\rho(s)\tau(s)^2B_s\,ds. \tag{0.5}
\]

The full-space Hodge identity gives exactly (0.4): the term \(B_s\)
cancels algebraically.  Consequently Hessian positivity and reinforced
Prékopa alone do not prove the weighted-Fisher bound.

## 1. Regularity and approximation

The literal calculations below hold, for example, when
\(U\in C^4(\mathbb R^{1+d})\), \(D^2U\succeq0\), \(e^{-U}\) has sufficient
decay, the displayed derivatives are integrable, and the centered Poisson
solution \(g_s\) and its \(s\)-derivative lie in the required operator
domains.  Cutoffs then remove compact-support assumptions.

No dimension-free Poincaré estimate is used to construct \(g_s\).  In each
fixed finite dimension, a nondegenerate log-concave probability has a
finite (possibly dimension-dependent) Poincaré constant, enough for the
centered weak Poisson solution when \(\ell_s\in L^2(q_s)\).

In the application, \(r\) is obtained by convolving the transverse
variable with a unit Gaussian.  Prékopa makes \(U\) convex, and Gaussian
differentiation gives

\[
D^2_{zz}U(s,z)
=I-\operatorname{Cov}(X_\perp\mid S=s,X_\perp+G=z),
\]

so

\[
0\preceq H_s\preceq I. \tag{1.1}
\]

For hard support or a nonsmooth potential, one may convolve the full law
with a small Gaussian, add vanishing quadratic confinement, prove cutoff
identities, and pass through the closed Dirichlet forms.  Equalities
(0.1)--(0.4) consist of nonnegative form budgets and survive this
regularization in the usual weak/lower-semicontinuous sense.  An upper
bound on (0.5) does not pass automatically: a proof of the missing
inequality must be uniform before the limit.  This is the precise
hard-support caveat.

## 2. Continuity and Bochner identities

Since \(L_sg_s=\ell_s=\partial_s\log q_s\),

\[
\partial_sq_s=\ell_sq_s=q_sL_sg_s
=\operatorname{div}_z(q_sF_s). \tag{2.1}
\]

Hence, for every regular test \(h(s,z)\),

\[
\frac d{ds}\mathbb E_sh
=\mathbb E_s(\partial_sh-F_s\cdot\nabla_zh). \tag{2.2}
\]

The integrated Bochner formula is

\[
\mathbb E_s(L_sg_s)^2
=\mathbb E_s\|D_z^2g_s\|_{\rm HS}^2
+ \mathbb E_s\langle H_s\nabla g_s,\nabla g_s\rangle,
\]

which proves (0.1).

## 3. Material acceleration

Differentiating \(L_sg_s=\ell_s\) gives

\[
L_s(\partial_sg_s)
=-U_{ss}+\Phi''+\langle\nabla_zU_s,F_s\rangle. \tag{3.1}
\]

The pointwise Bochner formula, together with
\(\nabla_z\ell_s=-\nabla_zU_s\), gives

\[
L_s\frac{|F_s|^2}{2}
=-\langle F_s,\nabla_zU_s\rangle
+ \|D_zF_s\|_{\rm HS}^2+\langle H_sF_s,F_s\rangle. \tag{3.2}
\]

Subtracting proves (0.2), because

\[
\mathcal Q_s
=U_{ss}-2\langle U_{sz},F_s\rangle+\langle H_sF_s,F_s\rangle.
\]

Integrating \(L_s\alpha_s\) against \(q_s\) proves (0.3).  Symmetry of
\(D_zF_s=D_z^2g_s\) also yields the material-acceleration identity

\[
\boxed{\partial_sF_s-(D_zF_s)F_s=\nabla_z\alpha_s}. \tag{3.3}
\]

## 4. Centroid identity and isotropy cancellation

Let \(m(s)=\mathbb E_sZ\).  Applying (2.2), first to \(Z\) and then to
\(F_s\), gives

\[
m'=-\mathbb E_sF_s,\qquad m''=-\mathbb E_s\nabla_z\alpha_s. \tag{4.1}
\]

Set

\[
R_s=\mathcal Q_s+\|D_zF_s\|_{\rm HS}^2.
\]

Self-adjointness of \(L_s\) gives, componentwise,

\[
\mathbb E_s[(Z-m)L_s\alpha_s]=-\mathbb E_s\nabla_z\alpha_s.
\]

Using \(L_s\alpha_s=\Phi''-R_s\), the constant term drops and therefore

\[
\boxed{m''(s)=-\mathbb E_s[(Z-m(s))R_s(Z)]}. \tag{4.2}
\]

Thus the nonnegative charge \(R_s\) has mass
\(\mathbb E_sR_s=\Phi''(s)\), while its transverse first moment is exactly
the acceleration of the conditional centroid.

If the original joint law is isotropic, adding independent transverse
Gaussian noise preserves \(\operatorname{Cov}(S,Z)=0\).  The scalar Stein
identity and (4.1) then imply

\[
0=\mathbb E[S\,m(S)]
=\mathbb E[\tau(S)m'(S)]
=-\int\rho(s)\tau(s)\mathbb E_sF_s\,ds. \tag{4.3}
\]

Equations (4.2)--(4.3) are the joint rigidity data unavailable to a
generic conditional Poincaré argument.

## 5. The one-dimensional Stein budget

For centered variance-one \(\rho=e^{-\Phi}\), the canonical kernel

\[
\tau(s)=\rho(s)^{-1}\int_s^\infty t\rho(t)\,dt
\]

satisfies

\[
(\rho\tau)'=-s\rho,\qquad
\tau'=\tau\Phi'-s,\qquad
\mathbb E_\rho\tau=1. \tag{5.1}
\]

Differentiating the middle identity, multiplying by \(\tau\), and
integrating by parts yields

\[
\mathbb E_\rho[\tau^2\Phi'']
=\mathbb E_\rho[\tau\tau''-\tau\tau'\Phi'+\tau]
=1-\mathbb E_\rho(\tau')^2. \tag{5.2}
\]

Together with (0.3), this proves (0.4).  Truncation proves the formula in
the finite regularized setting; lower semicontinuity gives the inequality
needed for the \(C_s\) budget in limits.

## 6. Full-space Hodge ledger

Define

\[
V(s,z)=(\tau(s),-\tau(s)F_s(z)).
\]

Equations (2.1) and (5.1) imply directly

\[
-\operatorname{div}(e^{-U}V)=s e^{-U}. \tag{6.1}
\]

For a regular vector field \(V\), weighted integration by parts gives

\[
\mathbb E(\operatorname{div}V-\nabla U\cdot V)^2
=\mathbb E\operatorname{Tr}[(DV)^2]+\mathbb ED^2U[V,V]. \tag{6.2}
\]

The first term is \(\operatorname{Tr}[(DV)^2]\), not
\(\|DV\|_{\rm HS}^2\).  Here

\[
DV=
\begin{pmatrix}
\tau'&0\\
-\tau'F-\tau\partial_sF&-\tau D_zF
\end{pmatrix},
\]

and hence

\[
\operatorname{Tr}[(DV)^2]
=(\tau')^2+\tau^2\|D_zF\|_{\rm HS}^2,\qquad
D^2U[V,V]=\tau^2\mathcal Q_s. \tag{6.3}
\]

Since \(\mathbb ES^2=1\), (6.1)--(6.3) give

\[
1=\mathbb E(\tau')^2+
\int\rho\tau^2(C_s+\mathbb E_s\mathcal Q_s)\,ds. \tag{6.4}
\]

This is exactly (0.4).  Although the expansion of \(\mathcal Q_s\)
contains \(B_s\), the horizontal and mixed Hessian terms cancel it.  The
Hodge identity therefore supplies no control of \(\mathcal B\).

## 7. Translation stress test

Let \(q_0=e^{-W}\) and

\[
q_s(z)=q_0(z-as),\qquad U(s,z)=\Phi(s)+W(z-as).
\]

Then

\[
\ell_s=a\cdot\nabla W(z-as),\qquad
g_s(z)=-a\cdot z+\text{constant},\qquad F_s=-a,
\]

so

\[
C_s=0,\qquad \mathcal Q_s=\Phi'',\qquad
B_s=\mathbb E_{q_0}\langle D^2W\,a,a\rangle. \tag{7.1}
\]

The pointwise ratio \(B_s/\Phi''\) can be arbitrarily large.  Convexity
alone cannot prove the desired estimate.  But
\(\operatorname{Cov}(S,Z)=a\operatorname{Var}(S)=a\), so isotropy forces
\(a=0\).  This shows why the global cancellation (4.3) is indispensable.

There is exact rigidity for a smoothly bending pure translate.  If

\[
U(s,z)=\Phi(s)+W(z-a(s)),
\]

then, along \((1,a'(s))\),

\[
D^2U[(1,a'),(1,a')]
=\Phi''-\nabla W(z-a(s))\cdot a''(s). \tag{7.2}
\]

If \(W\) is Legendre and \(\nabla W(\mathbb R^d)=\mathbb R^d\), joint
convexity forces \(a''=0\).  Full-support pure translations therefore
have constant velocity, which isotropy removes.  Bending remains possible
for hard or changing fibers; cones and wedges remain mandatory tests.

## 8. A direct Gaussian-channel reduction

There is a second exact formulation.  Let \(p_s(x)\) denote the
unsmoothed conditional slice and suppose it admits a continuity velocity
\(v_s\):

\[
\partial_sp_s+\operatorname{div}_x(p_sv_s)=0. \tag{8.1}
\]

For \(Z=X+G\), with \(G\) standard Gaussian independent of \(X\), Gaussian
integration by parts in (8.1) gives

\[
\partial_s\log q_s(Z)
=\mathbb E[v_s(X)\cdot G\mid Z,S=s]. \tag{8.2}
\]

Conditional Jensen and independence give the dimension-free contraction

\[
I_\perp(s)\le\mathbb E_{p_s}|v_s|^2. \tag{8.3}
\]

Consequently WFI would also follow from

\[
\int\rho(s)\tau(s)^2\mathbb E_{p_s}|v_s|^2\,ds\le C. \tag{8.4}
\]

Writing \(w=\tau v\), the field \(T=(\tau,w)\) satisfies

\[
\operatorname{div}_{s,x}(pT)=-sp,\qquad \mathbb Ew=0, \tag{8.5}
\]

where the last equality follows from isotropy.  Thus (8.4) asks for a
dimension-free \(L^2\) bound on the transverse part of one directional
Stein-kernel column.  Formula (8.3) is unconditional, but no universal
kinetic bound (8.4) is proved here; choosing a Poisson velocity and
bounding its energy by conditional Poincaré would again be circular.

## 9. Remaining formal lemma

A sufficient load-bearing statement is the following.

**Weighted slice-rigidity lemma (unproved).** There is a numerical
\(C<\infty\) such that, for every dimension and every regular jointly
log-concave density in the normalized setting above,

\[
\int\rho\tau^2\mathbb E_s\langle H_sF_s,F_s\rangle
\le C\int\rho\tau^2\left(
\mathbb E_sD^2U[(1,-F_s),(1,-F_s)]
+\mathbb E_s\|D_zF_s\|_{\rm HS}^2\right). \tag{WK}
\]

Its right side is at most one by (0.4), and (0.1) would then prove WFI.
Any proof must use global isotropy through (4.3) and joint slice geometry
through (4.2).  A generic one-form Poincaré estimate is KLS-strength and
cannot be inserted as an intermediate step.

## 10. Complete positive class: arbitrary pure translations

The kinetic formulation proves WFI for every translating conditional
family, including asymmetric exponential fibers.  The proof uses the
score-range inradius; a score-radius argument would be false.

Assume first that

\[
p_s(x)=Z_\varphi^{-1}e^{-\varphi(x-m(s))},\qquad
p(s,x)=e^{-\Phi(s)}p_s(x),
\]

is smooth and the joint law is isotropic.  If \(Y\) has density
proportional to \(e^{-\varphi}\), independence in this representation
gives

\[
\operatorname{Cov}(Y)\preceq I. \tag{10.1}
\]

Let

\[
K=\overline{\operatorname{conv}}\{\nabla\varphi(y):y\in\mathbb R^d\}.
\]

Then

\[
\boxed{K\supset c_0B_2^d,\qquad c_0=3/16.} \tag{10.2}
\]

Here is a self-contained proof.  Fix a unit vector \(u\), and let \(f_u\)
be the log-concave density of \(u\cdot Y\).  Its variance is at most one.
Chebyshev gives mass at least \(3/4\) in an interval of length four, hence

\[
M_u:=\|f_u\|_\infty\ge3/16. \tag{10.3}
\]

Write \(f_u=e^{-V_u}\), and let \(t_0\) be a mode.  If
\(a_u=\sup_tV_u'(t)<\infty\), convexity gives, for \(t\ge t_0\),

\[
f_u(t)\ge M_u e^{-a_u(t-t_0)}.
\]

Thus \(1\ge M_u/a_u\), so \(a_u\ge M_u\).  If \(a_u=\infty\), the same
conclusion is automatic.  Differentiating the marginal shows

\[
V_u'(t)=\mathbb E[u\cdot\nabla\varphi(Y)\mid u\cdot Y=t].
\]

Consequently \(h_K(u)\ge a_u\ge3/16\).  This holds for every unit \(u\);
the support-function characterization of convex containment proves
(10.2).  The same argument covers a bounded marginal endpoint by assigning
the corresponding score endpoint value \(+\infty\).

Test joint convexity on the vector \((1,m'(s))\).  The terms containing
\(D^2\varphi\) cancel exactly and give

\[
\Phi''(s)-\langle\nabla\varphi(y),m''(s)\rangle\ge0
\quad\text{for every }y. \tag{10.4}
\]

Equations (10.2)--(10.4) imply

\[
|m''(s)|\le c_0^{-1}\Phi''(s). \tag{10.5}
\]

Choose a mode \(s_0\) with \(\Phi'(s_0)=0\), and put \(b=m'(s_0)\).
Integration of (10.5) gives

\[
|m'(s)-b|\le c_0^{-1}|\Phi'(s)|. \tag{10.6}
\]

Cross-isotropy and the scalar Stein identity give

\[
\mathbb E[\tau(S)m'(S)]=\mathbb E[S\,m(S)]=0. \tag{10.7}
\]

Since \(\tau\Phi'=\tau'+S\), (5.2) implies

\[
\mathbb E[\tau^2(\Phi')^2]\le
2\mathbb E(\tau')^2+2\mathbb ES^2\le4. \tag{10.8}
\]

Using \(\mathbb E\tau=1\), (10.6)--(10.8) give
\(|b|\le2c_0^{-1}\).  The one-dimensional estimate
\(\mathbb E\tau^2\le400\) from the weighted-Fisher reduction now yields

\[
\boxed{\mathbb E[\tau^2|m'|^2]\le3208c_0^{-2}.} \tag{10.9}
\]

The translating slices satisfy (8.1) with \(v_s=m'(s)\).  Therefore
(8.3) and (10.9) prove WFI for this entire class.

The nonsmooth assertion has a direct finite-difference proof; no
regularization preserving the translating form is needed.  For \(h>0\),
apply joint convexity to the two points
\((s-h,m(s-h)+y)\) and \((s+h,m(s+h)+y)\).  With

\[
\delta_h=\frac{m(s-h)+m(s+h)}2-m(s),
\]

one obtains, for every \(y\),

\[
\varphi(y+\delta_h)-\varphi(y)
\le\frac{\Phi(s-h)+\Phi(s+h)-2\Phi(s)}2. \tag{10.10}
\]

Convexity of \(\varphi\) implies

\[
h_K(\delta_h)
\le\sup_y[\varphi(y+\delta_h)-\varphi(y)].
\]

Using (10.2) in (10.10) therefore gives the exact second-difference
domination

\[
c_0|m(s+h)-2m(s)+m(s-h)|
\le\Phi(s+h)-2\Phi(s)+\Phi(s-h). \tag{10.11}
\]

Passing to distributions shows that \(m'\) is locally of bounded
variation and

\[
|D m'|\le c_0^{-1}D\Phi' \tag{10.12}
\]

as vector and nonnegative Radon measures.  At a mode, choose the
representative \(b\) of the jump interval of \(m'\) corresponding to
\(0\in\partial\Phi\).  Integrating (10.12) on either side gives (10.6)
for almost every \(s\).  The one-dimensional Stein identity for locally
absolutely continuous \(m\), followed by truncation, gives
(10.7)--(10.9).  Thus the theorem covers asymmetric Laplace fibers and
hard one-dimensional endpoints, including curvature atoms at the mode.

The distinction from the unsmoothed Korn functional is essential.  For
the convex family with \(S\) Gaussian, asymmetric-Laplace score range
\([-\varepsilon,1]\), and
\(m(s)=-(s^2-1)/(2\varepsilon)\), the unsmoothed curvature term grows like
\(1/\varepsilon\).  After exact transverse whitening,
\(\sigma^2=1+3/(2\varepsilon^2)\), and unit Gaussian smoothing, the
translation velocity satisfies

\[
|F_s|^2=\frac{s^2}{\varepsilon^2\sigma^2}\le\frac23s^2.
\]

Since the smoothed transverse Hessian is at most one, its actual WFI is at
most \(2/3\).  The family refutes an unsmoothed shortcut but fully satisfies
the required post-noise estimate.

