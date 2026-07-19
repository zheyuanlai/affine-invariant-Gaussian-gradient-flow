# Two-law convolution amplification and fixed-Gaussian output structure

## 0. Result

Let \(\mu_0,\mu_1\) be centered isotropic probability measures on
\(\mathbb R^n\), each with finite Poincare constant, and let

\[
 X\sim\mu_0,\qquad Y\sim\mu_1
\]

be independent.  Put

\[
 S=\frac{X+Y}{\sqrt2},\qquad \nu=\mathcal L(S),
\]

and define

\[
 a=\min\{\lambda_1(\mu_0),\lambda_1(\mu_1)\},
 \qquad b=\lambda_1(\nu).
\]

Then, without assuming that any spectral edge is attained,

\[
 \boxed{\qquad b(1-b)\le4(b-a),\qquad
 b\ge\frac{\sqrt{9+16a}-3}{2}.\qquad}                 \tag{0.1}
\]

The proof is the unequal-law version of the normalized
self-convolution argument.  Equality of the two input laws is not used;
only their common lower spectral bound and isotropic covariance are used.

The fixed-Gaussian specialization supplies useful output structure.  The
forward theorem alone does not give a reverse reduction, but a separate
posterior theorem does.  Let \(\gamma=N(0,I_n)\) and

\[
 \mathcal G\mu=\mathcal L\!\left(\frac{X+G}{\sqrt2}\right),
 \qquad G\sim\gamma.
\]

For \(a=\lambda_1(\mu)\) and
\(b=\lambda_1(\mathcal G\mu)\), (0.1) rearranges only to

\[
 a\le\frac{b(3+b)}4.                                  \tag{0.2}
\]

This is a forward constraint and gives no lower bound on the input gap
from a lower bound on the regularized output gap.  Independently, posterior
unit strong log-concavity proves

\[
 \boxed{a\ge\frac{b}{2+b}.}                           \tag{0.2a}
\]

The proof, including form-domain spectral windows, hard supports, and
intrinsic affine supports, is in `posterior_reverse_smoothing.md`.

The regularized law has a positive analytic log-concave density.  If
\(\widetilde U\) is its potential in the isotropic \(S\)-coordinates, then

\[
 0\preceq D^2\widetilde U\preceq2I.                   \tag{0.3}
\]

Consequently the analytic isotropic class (0.3) is a quantitatively
equivalent KLS target: restriction gives one direction, while (0.2a)
transfers any universal gap for this class back to every isotropic
log-concave input.  This equivalence uses the separate posterior theorem,
not the forward inequality (0.1).

## 1. Unequal Hoeffding decomposition

Let \(f\) be centered and in the form domain of \(\nu\).  Define

\[
 F(x,y)=f\!\left(\frac{x+y}{\sqrt2}\right),\qquad
 h_0(x)=E_YF(x,Y),\qquad h_1(y)=E_XF(X,y),
\]

\[
 R(x,y)=F(x,y)-h_0(x)-h_1(y).
 \tag{1.1}
\]

Then

\[
 E[R\mid X]=E[R\mid Y]=0
\]

and the scalar Hoeffding decomposition is orthogonal:

\[
 \|f\|_{L^2(\nu)}^2
 =\|h_0\|_{L^2(\mu_0)}^2
  +\|h_1\|_{L^2(\mu_1)}^2
  +\|R\|_{L^2(\mu_0\otimes\mu_1)}^2.                 \tag{1.2}
\]

For a smooth \(f\), put

\[
 Z=\nabla f(S),\qquad
 U=E[Z\mid X],\qquad V=E[Z\mid Y].
\]

Differentiation under the conditional integrals gives

\[
 \nabla h_0=\frac{U}{\sqrt2},\qquad
 \nabla h_1=\frac{V}{\sqrt2},
\]

\[
 \nabla_xR=\frac{Z-U}{\sqrt2},\qquad
 \nabla_yR=\frac{Z-V}{\sqrt2}.                        \tag{1.3}
\]

These identities extend to the form domain by Sobolev approximation.
All cross terms vanish by conditional centering, so, with

\[
 q=\int|\nabla f|^2\,d\nu,\qquad
 E_R=\int|\nabla_{x,y}R|^2\,d(\mu_0\otimes\mu_1),
\]

one has

\[
 q=\int|\nabla h_0|^2\,d\mu_0
   +\int|\nabla h_1|^2\,d\mu_1+E_R.                  \tag{1.4}
\]

Since both input gaps are at least \(a\), applying the one-factor
Poincare inequality to \(R\) in each coordinate gives

\[
 E_R\ge2a\|R\|_2^2.                                   \tag{1.5}
\]

Applying the input inequalities also to \(h_0,h_1\), and using (1.2),
gives, for \(\|f\|_2=1\),

\[
 q\ge a(1-\|R\|_2^2)+E_R.
\]

Together with (1.5), this yields the same scalar defect estimate as in
the equal-law case:

\[
 \boxed{\qquad E_R\le2(q-a).\qquad}                   \tag{1.6}
\]

## 2. Unequal vector predictors

Let \(m=EZ\) and put

\[
 U_0=U-m,\qquad V_0=V-m,\qquad
 W=Z-U-V+m.
\]

The three vector fields \(U_0(X)\), \(V_0(Y)\), and \(W(X,Y)\) are
orthogonal in \(L^2\), although the first two need not have equal
variances.  Write

\[
 u=E|U_0|^2,\qquad v=E|V_0|^2,\qquad w=E|W|^2.
\]

Then

\[
 \operatorname {Var}(Z)=u+v+w.                       \tag{2.1}
\]

On the other hand, (1.3) gives

\[
 E_R=\frac12E|Z-U|^2+\frac12E|Z-V|^2
     =\frac{u+v}{2}+w.                                \tag{2.2}
\]

Therefore

\[
 \boxed{\qquad
 \operatorname {Var}_\nu(\nabla f)
 \le2E_R\le4(q-a).
 \qquad}                                              \tag{2.3}
\]

No equality of the predictor variances is needed.

## 3. Continuous bottom spectral edge

Fix \(\varepsilon>0\) and choose a centered unit vector \(f\) in the
spectral subspace
\(\mathbf1_{[b,b+\varepsilon]}(A_\nu)L^2_0(\nu)\).  Write

\[
 q=\langle A_\nu f,f\rangle=b+\alpha,\qquad
 z=(A_\nu-b)f,
\]

where

\[
 0\le\alpha\le\varepsilon,\qquad
 \|z\|_2^2\le\varepsilon\alpha\le\varepsilon^2.
\]

Because \(S\) is isotropic, \(\ell=E[Sf(S)]\) satisfies
\(|\ell|\le1\).  Testing the weak generator identity against coordinate
functions gives

\[
 E\nabla f=b\ell+E[Sz],
\qquad |E[Sz]|\le\varepsilon.
\]

Consequently

\[
 \operatorname {Var}_\nu(\nabla f)
\ge b(1-b)-2\varepsilon-\varepsilon^2.                \tag{3.1}
\]

Combining (2.3) and (3.1) and sending
\(\varepsilon\downarrow0\) gives

\[
 b(1-b)\le4(b-a).
\]

Solving \(b^2+3b-4a\ge0\) proves (0.1).

## 4. Fixed Gaussian regularization

Take \(\mu_0=\mu\) and \(\mu_1=\gamma\).  Since an isotropic law has
\(\lambda_1(\mu)\le1=\lambda_1(\gamma)\), the value of \(a\) in (0.1) is
exactly \(\lambda_1(\mu)\).  The first inequality in (0.1) gives

\[
 a\le b-\frac{b(1-b)}4=\frac{b(3+b)}4.                \tag{4.1}
\]

This is the forward direction already encoded in (0.1).  A lower bound on
\(b\) does not imply any lower bound on \(a\), so no fixed-noise reverse
transfer follows from (0.1) itself.  The separate posterior comparison in
`posterior_reverse_smoothing.md` gives

\[
 b\le\frac{2a}{1-a}\quad(a<1),
 \qquad a\ge\frac{b}{2+b}.                            \tag{4.2}
\]

Thus a lower bound on the fixed-Gaussian output gap does transfer to the
input, but through (4.2), not through the unequal-law forward argument.

It remains to verify the structural description (0.3).  Before the
\(1/\sqrt2\) rescaling, let \(Y=X+G\), let \(q\) be its density, and put
\(U=-\log q\).  The Gaussian posterior formula gives

\[
 D^2U(y)=I-\operatorname {Cov}(X\mid X+G=y).
\]

The posterior is \(1\)-strongly log-concave, including hard convex
supports by approximation, so

\[
 0\preceq\operatorname {Cov}(X\mid X+G=y)\preceq I.
\]

Hence \(0\preceq D^2U\preceq I\).  In the isotropic coordinate
\(s=y/\sqrt2\), the potential is
\(\widetilde U(s)=U(\sqrt2s)+\text{constant}\), and therefore

\[
 0\preceq D^2\widetilde U(s)\preceq2I.
\]

Gaussian convolution preserves log-concavity and makes the density
positive and analytic.  The same argument is performed intrinsically on
the affine hull for a lower-dimensional input.  A point mass remains
excluded by convention.

## 5. Exact scope

The unequal-law forward theorem by itself supplies no fixed-noise reverse
reduction: its correct implication is only \(a\le b(3+b)/4\).  The
posterior theorem (4.2) is an additional argument and makes the positive
analytic isotropic output class

\[
 0\preceq D^2\widetilde U\preceq2I
\]

quantitatively equivalent to the full KLS problem.  Upper curvature is not
a Bakry--Emery lower-curvature hypothesis, so establishing a universal gap
on this output class remains the substantive nonlinear task; once such a
gap \(b\ge c_0\) is proved, (4.2) gives
\(a\ge c_0/(2+c_0)\) for the unsmoothed input.
