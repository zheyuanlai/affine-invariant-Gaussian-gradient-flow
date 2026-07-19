# Stress audit of weighted slice rigidity

## 0. Corrected conclusion

For the post-noise weighted-Fisher setup, the general shape-changing-slice
estimate remains unproved.  Two conclusions are rigorous.

1. WFI holds dimension-freely for every pure translating conditional
   family.  The proof uses the inradius of the convex hull of the score
   range, not its radius.
2. A smooth asymmetric-Laplace construction makes the analogous
   unsmoothed curvature ratio diverge.  Exact whitening followed by the
   required unit transverse Gaussian noise restores a universal bound.
   Thus transverse smoothing is load-bearing.

The complete material-acceleration and Hodge identities are recorded in
work/weighted_fisher_prekopa.md.

## 1. Exact translated-slice formulas

Let

\[
p_s(x)=Z_\varphi^{-1}e^{-\varphi(x-m(s))},\qquad
p(s,x)=e^{-\Phi(s)}p_s(x), \tag{1.1}
\]

where \(U_0(s,x)=\Phi(s)+\varphi(x-m(s))\) is convex.  Let \(Y\) have
density proportional to \(e^{-\varphi}\).  The slices obey

\[
\partial_sp_s+\operatorname{div}_x(p_sm'(s))=0. \tag{1.2}
\]

After adding independent unit Gaussian noise \(G\), put
\(q_s=\mathcal L(m(s)+Y+G)\).  If \(\psi\) is the potential of \(Y+G\),
then

\[
\ell_s(z)=m'(s)\cdot\nabla\psi(z-m(s)),\qquad
g_s(z)=-m'(s)\cdot(z-m(s)),
\]

so

\[
F_s=-m'(s),\qquad D_zF_s=0,\qquad
I_\perp(s)=m'(s)^TJ(Y+G)m'(s)\le |m'(s)|^2. \tag{1.3}
\]

The last inequality is also the direct Gaussian-channel contraction from
(1.2).  Moreover,

\[
\mathbb E_s\mathcal Q_s=\Phi''(s),\qquad
0=\mathbb E[S\,m(S)]=\mathbb E[\tau(S)m'(S)]. \tag{1.4}
\]

## 2. Score-range inradius forced by isotropy

Assume the original joint law in (1.1) is isotropic.  Independence in
\(X=m(S)+Y\) gives

\[
\operatorname{Cov}(Y)\preceq I. \tag{2.1}
\]

In the smooth full-support case define

\[
K=\overline{\operatorname{conv}}\{\nabla\varphi(y):y\in\mathbb R^d\}.
\]

For extended-valued nonsmooth \(\varphi\), use the union of all
subdifferentials, including normal cones at hard boundaries.

**Lemma 2.1.** One has

\[
K\supset c_0B_2^d,\qquad c_0=3/16. \tag{2.2}
\]

**Proof.** Fix \(|u|=1\).  The marginal \(T=u\cdot Y\) has variance at
most one.  Chebyshev puts at least \(3/4\) of its mass in an interval of
length four, so its log-concave density \(f_u\) satisfies

\[
M_u:=\|f_u\|_\infty\ge3/16. \tag{2.3}
\]

Write \(f_u=e^{-V_u}\), and choose a mode \(t_0\).  If
\(a_u=\sup V_u'<\infty\), then for \(t\ge t_0\),

\[
f_u(t)\ge M_u e^{-a_u(t-t_0)}.
\]

Integration gives \(a_u\ge M_u\).  The conclusion is automatic when
\(a_u=\infty\), including a hard right endpoint.  Differentiating the
smooth marginal gives

\[
V_u'(t)=\mathbb E[u\cdot\nabla\varphi(Y)\mid u\cdot Y=t].
\]

Thus \(h_K(u)\ge a_u\ge3/16\).  This holds for every unit \(u\), proving
(2.2) by the support-function characterization of convex containment.
The subgradient version follows from the same marginal argument. \(\square\)

The inradius is essential.  A score norm upper bound alone gives no lower
bound in the opposite direction for an asymmetric Laplace law.

## 3. Curvature controls translation acceleration

In the smooth case, test \(D^2U_0\succeq0\) on \((1,m'(s))\).  The
\(D^2\varphi\) terms cancel, leaving

\[
\Phi''(s)-\langle\nabla\varphi(y),m''(s)\rangle\ge0
\quad\text{for all }y.
\]

Lemma 2.1 gives

\[
|m''(s)|\le c_0^{-1}\Phi''(s). \tag{3.1}
\]

The nonsmooth statement follows directly by finite differences.  Apply
joint convexity to \((s-h,m(s-h)+y)\) and
\((s+h,m(s+h)+y)\).  With

\[
\delta_h=\frac{m(s-h)+m(s+h)}2-m(s),
\]

one gets

\[
\varphi(y+\delta_h)-\varphi(y)
\le\frac{\Phi(s-h)+\Phi(s+h)-2\Phi(s)}2.
\]

The supremum of the left side dominates \(h_K(\delta_h)\), so

\[
c_0|m(s+h)-2m(s)+m(s-h)|
\le\Phi(s+h)-2\Phi(s)+\Phi(s-h). \tag{3.2}
\]

Passing to distributions gives the Radon-measure inequality

\[
|D m'|\le c_0^{-1}D\Phi'. \tag{3.3}
\]

At a mode, choose the representative \(b\) in the jump interval of \(m'\)
corresponding to \(0\in\partial\Phi\).  Integration gives

\[
|m'(s)-b|\le c_0^{-1}|\Phi'(s)| \quad\text{a.e.} \tag{3.4}
\]

## 4. Dimension-free kinetic bound

The canonical Stein kernel satisfies

\[
\tau\Phi'=\tau'+S,\qquad
\mathbb E(\tau')^2\le1,\qquad
\mathbb E\tau=1.
\]

Hence

\[
\mathbb E[\tau^2(\Phi')^2]\le4. \tag{4.1}
\]

Cross-isotropy in (1.4), followed by (3.4) and Cauchy--Schwarz, gives

\[
|b|\le2c_0^{-1}. \tag{4.2}
\]

Using the established one-dimensional estimate
\(\mathbb E\tau^2\le400\),

\[
\boxed{\mathbb E[\tau^2|m'|^2]\le3208c_0^{-2}.} \tag{4.3}
\]

Equations (1.3) and (4.3) prove WFI for every pure translating family,
uniformly in the transverse dimension.  The finite-difference proof
covers curvature atoms, asymmetric exponential tails, and hard endpoints.

## 5. Asymmetric-Laplace pre-noise obstruction

Fix \(0<\varepsilon<1\).  Let \(Y_\varepsilon\) have the centered
asymmetric-Laplace law whose potential slopes are
\(-\varepsilon\) left of its kink and \(1\) right of it.  Direct integration
gives

\[
\operatorname{Var}(Y_\varepsilon)=1+\varepsilon^{-2},
\qquad J(Y_\varepsilon)=\varepsilon.
\]

Let \(S\sim N(0,1)\) and

\[
m(s)=-\frac{s^2-1}{2\varepsilon},\qquad
U_0(s,z)=\frac{s^2}{2}+\varphi_\varepsilon(z-m(s)). \tag{5.1}
\]

The Schur expression is
\(1+\varphi_\varepsilon'/\varepsilon\ge0\), while the kink contributes a
positive rank-one curvature measure; hence \(U_0\) is convex.  Also
\(\mathbb Em(S)=\mathbb E[S\,m(S)]=0\).

The exact transverse variance is

\[
\sigma_\varepsilon^2=1+\frac{3}{2\varepsilon^2}. \tag{5.2}
\]

After replacing \(z\) by \(x=z/\sigma_\varepsilon\), the joint law is
isotropic.  Before Gaussian smoothing,

\[
F_s=-\frac{m'(s)}{\sigma_\varepsilon},\qquad
B_s=\frac{s^2}{\varepsilon},\qquad
\mathbb EB_s=\frac1\varepsilon,\qquad
\mathbb E\mathcal Q_s=1. \tag{5.3}
\]

Thus the unsmoothed weighted-Korn ratio has no universal constant.
Smooth convex soft-max approximations preserve the score interval and
Schur inequality, so this is not a kink artifact.

## 6. Unit smoothing repairs this family

After unit transverse Gaussian smoothing,
\(0\le\psi_\varepsilon''\le1\), and (5.2) gives

\[
B_s\le |F_s|^2
=\frac{s^2}{\varepsilon^2\sigma_\varepsilon^2}
\le\frac23s^2.
\]

The Gaussian marginal has \(\tau=1\), so

\[
\int\rho\tau^2B_s\,ds\le\frac23. \tag{6.1}
\]

This family refutes only the pre-noise shortcut.  It satisfies the actual
post-noise WFI bound.

## 7. Remaining scope

Pure translations are controlled.  The unresolved case is a genuinely
shape-changing family whose Poisson field has a substantial nonlinear
component.  Replacing that component by conditional Poincaré is circular.

