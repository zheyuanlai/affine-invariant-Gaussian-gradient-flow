# Shape-changing slices: centroid stability audit (clean version)

This note isolates the centroid part of the weighted-Fisher problem. Let
\(\rho=e^{-\Phi}\) be a centered, variance-one one-dimensional
log-concave marginal, let \(\tau\) be its canonical Stein kernel, and let
\(m(s)=E[Z\mid S=s]\). In the notation of
weighted_fisher_prekopa.md,
\[
 F_s=-m'(s)+F_s^0,\quad E_sF_s^0=0,\quad
 m''(s)=-E_s[y\,R_s],\quad
 R_s=\mathcal Q_s+\|D_zF_s\|_{\rm HS}^2 ,
\]
and
\[
 \int\rho(s)\tau(s)^2 E_s\|D_zF_s\|_{\rm HS}^2\,ds\le1,\qquad
 E[\tau(S)m'(S)]=0.                                      \tag{0.1}
\]
The last equality is the mixed-isotropy cancellation. The unresolved
general problem is to control the first moment \(E_s[yR_s]\) from the mass
budget in (0.1).

Two genuine shape-changing classes admit a dimension-free centroid bound.
The proofs below use no conditional Poincare/KLS inequality.

## 1. One-dimensional facts

We use standard elementary estimates for centered variance-one log-concave
laws:
\[
 |s_0|\le4\quad(s_0\ {\rm a\ mode}),\qquad E|S|^4\le M_4,\qquad
 E\tau^2\le400,                                      \tag{1.1}
\]
with a universal \(M_4\) (take \(M_4=10^5\)). Convexity of \(\Phi\) and
the tail definition of \(\tau\) give, at every differentiability point,
\[
                 \tau(s)|\Phi'(s)|\le1+|s|.             \tag{1.2}
\]
Indeed, in the direction where \(\Phi'\) has sign \(+\), write
\(\rho(s+u)\le\rho(s)e^{-u\Phi'(s)}\), bound
\(|s+u|\le|s|+u\), and use the right-tail formula for \(\tau\); the left
side is identical. Thus
\[
 E[\tau^2(1+|S|)^2(\Phi')^2]\le E(1+|S|)^4\le C_0.       \tag{1.3}
\]
We also use: if \(h\ge0\) is concave and \(Eh(S)\le\Lambda\), then
\[
 h(s)\le L_0\Lambda(1+|s|).                              \tag{1.4}
\]
This follows from the universal isotropic core \([-b,b]\) on which
\(\rho\ge c_b\), first bounding \(h(0)\) by \(2\Lambda/(bc_b)\), then using
the secant slope and nonnegativity.

## 2. Gaussian fibers with changing, noncommuting covariance

Consider
\[
 p(s,z)=Z^{-1}\exp\{-W(s)-\tfrac12(z-m(s))^TQ(s)(z-m(s))\},
 \quad R(s)=Q(s)^{-1}\succ0,                              \tag{2.1}
\]
on an interval times \(\mathbb R^d\). Assume joint log-concavity,
\(E[S\,m(S)]=0\), and \(ER(S)\preceq\Lambda I\). (For an isotropic law,
\(ER\preceq I\); after adding unit Gaussian transverse noise one may take
\(\Lambda=2\).) Then
\[
                    E[\tau(S)^2|m'(S)|^2]\le C_\Lambda. \tag{2.2}
\]

Write \(y=z-m(s)\). Direct differentiation and the Schur complement give
\[
 W''-m''{}^TQy+\tfrac12y^TQKQy\ge0,\qquad K:=-R''\succeq0. \tag{2.3}
\]
Equivalently, with \(x=Qy\),
\[
 W''-m''\!\cdot x+\tfrac12x^TKx\ge0\quad\forall x.        \tag{2.4}
\]
Hence \(m''\in{\rm Ran}\,K\) and \(m''{}^TK^\dagger m''\le2W''\). For
\(|u|=1\),
\[
 |u\!\cdot m''|^2\le2W''u^TKu
 \le2W''\lambda_{\max}(R)\operatorname{tr}(R^{-1}K).     \tag{2.5}
\]
The marginal potential is
\[
 \Phi=W-\tfrac12\log\det R+{\rm const},
\quad
 \Phi''=W''+\tfrac12\operatorname{tr}(R^{-1}K)
 +\tfrac12\|R^{-1/2}R'R^{-1/2}\|_{\rm HS}^2.             \tag{2.6}
\]
Therefore
\[
             |m''(s)|\le2\sqrt{\lambda_{\max}R(s)}\,\Phi''(s). \tag{2.7}
\]
Each \(u^TRu\) is nonnegative concave; (1.4) and \(ER\preceq\Lambda I\)
give \(\lambda_{\max}R(s)\le L_0\Lambda(1+|s|)\).

Choose a mode \(s_0\) with \(0\in\partial\Phi(s_0)\), set \(b=m'(s_0)\),
and integrate (2.7) from \(s_0\) to \(s\). Since \(|s_0|\le4\),
\[
 |m'(s)-b|\le C_1\sqrt{\Lambda}\sqrt{1+|s|}\,|\Phi'(s)|. \tag{2.8}
\]
By \(E[\tau m']=0\), \(E\tau=1\), (1.3), and Cauchy--Schwarz,
\[
 E[\tau^2|m'-b|^2]\le C_2\Lambda,\qquad |b|^2\le C_2\Lambda.
\]
Using \(E\tau^2\le400\) yields (2.2), with for instance
\(C_\Lambda=10^{12}\Lambda\). This remains valid in the
finite-difference/distributional limit.

For completeness, the centered Poisson field is affine:
\(F_s=-m'(s)+A_s y\), where \(A_s=A_s^T\) is the unique solution of
\[
                    QA_s+A_sQ=-QR'Q.                  \tag{2.9}
\]
Thus \(C_s=E_s\|D F_s\|_{\rm HS}^2=\|A_s\|_{\rm HS}^2\), and the budget in
(0.1) pays the shape component. No commutativity of \(R\) and \(R'\) is
used.

## 3. Hard wedges and cones (geometric stress test)

Let a planar convex body be
\[
 K=\{(s,z):s\in J,\ l(s)\le z\le u(s)\},\qquad
 m=(u+l)/2,\quad w=(u-l)/2.
\]
Then \(u\) is concave, \(l\) convex, and with
\(\mu_+=-D^2u\ge0,\ \mu_-=D^2l\ge0\),
\[
 |D m'|\le-D w'=\tfrac12(\mu_++\mu_-).                  \tag{3.1}
\]
The uniform marginal has \(\rho\propto w\), so
\[
 D\Phi'= (w'/w)^2\,ds+(-Dw')/w.                         \tag{3.2}
\]
If the transverse variance is normalized (\(E_\rho w^2\le3\)), (1.4)
gives \(w(s)\le L(1+|s|)\). Integrating (3.1) and using (1.2) gives
\[
 |m'(s)-b|\le C(1+|s|)|\Phi'(s)|,\qquad E[\tau m']=0,
\]
and (1.3) with the fourth moment bound yields
\(E[\tau^2|m'|^2]\le C\). Piecewise-affine \(u,l\) (triangles and
cones) are included; curvature atoms satisfy (3.1). The ideal uniform
fibers have centered continuity velocity \(v_s^0=(w'/w)y\), so its
deformation gradient has squared norm \((w'/w)^2\).

A sharper post-noise cone stress test is already explicit in
work/gaussian_curvature_korn.md, §§2--4. There
\[
 p(r,x)=\tfrac12e^{-r}{\bf1}_{\{r>0,\ |x|<r\}},
\quad S=(R-2)/\sqrt2,\quad Y=X/\sqrt2,
\]
is isotropic; \(Y\mid S=s\) is uniform on \([-a,a]\), \(a=s+\sqrt2\).
After adding \(G\sim N(0,1)\),
\(q_s={\rm Unif}[-a,a]*\gamma\), and the exact Poisson field is
\[
 F_s(z)=-\frac{E[Y\mid Y+G=z]}a,\qquad
 C_s=E_s|F_s'|^2\le a^{-2}.
\]
In the boundary layer \(z=a+t,\ 0\le t\le1/2\), the same note proves
\[
 B_s=E_s[U_s''F_s^2]\ge(256a)^{-1}.
\]
Hence the pointwise ratio \(B_s/(C_s+|E_sF_s|^2)\ge a/256\) diverges.
This disproves any slice-by-slice curvature--Korn estimate, even for the
exact Poisson field. Nevertheless \(\tau(s)=a/\sqrt2\), \(|F_s|\le1\),
\(0\le U_s''\le1\), and the full weighted integral obeys
\[
                         \int\rho\tau^2B_s\,ds\le3/2.   \tag{3.3}
\]
Thus tails and the scalar weighting are essential.

## 4. Rotating boxes / matrix shape

For Gaussian ellipsoids, \(R''\preceq0\). A volume-preserving pure rotation
would have constant eigenvalues; then \({\rm tr}\,R''=0\), so \(R''=0\), and
constant \({\rm tr}\,R^2\) forces \(R'=0\). Rotation must therefore pay by
changing eigenvalues. A concrete noncommuting path is
\[
 R(s)=I+sA-s^2B,\quad
 A=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
 B=\begin{pmatrix}1&r\\r&1\end{pmatrix},\quad0<|r|<1,
\]
on a sufficiently short interval. Then \(R''=-2B\preceq0\) but
\([A,B]\ne0\). Choosing \(m''\in{\rm Ran}\,B\) and
\(W''\ge\tfrac12m''{}^TB^\dagger m''\) satisfies (2.4); a convex barrier in
\(W\) makes a finite normalizable marginal. The Gaussian theorem applies
uniformly in the rotation angle.

## 5. Exact obstruction for general shape-changing fibers

The only universal identity currently available is
\[
 m''(s)=-E_s[y\,R_s],\qquad
 R_s=\mathcal Q_s+\|D_zF_s^0\|_{\rm HS}^2\ge0.            \tag{5.1}
\]
The budget controls
\[
 \int\rho\tau^2 E_sR_s\le1,                              \tag{5.2}
\]
but the needed quantity is the vector first moment
\[
 \int\rho\tau^2 E_s[yR_s]\,ds.                           \tag{5.3}
\]
There is no valid bound of (5.3) by (5.2) for an arbitrary nonnegative
charge: it can be concentrated at large \(|y|\). Conditional
Cauchy--Schwarz would introduce \(E_s|y|^2\) together with \(E_sR_s^2\),
and replacing this by \(E_sR_s\) is a conditional Poincare/KLS-strength
step. The cone calculation above shows that even the pointwise ratio
\(E_s\langle H_sF_s,F_s\rangle/(C_s+|E_sF_s|^2)\) can diverge. Any complete
extension must discover a joint convex-geometric analogue of the Gaussian
Schur inequality (2.4) or the interval endpoint inequality (3.1), using
the global isotropy cancellation and the material charge; no such lemma is
claimed here.

