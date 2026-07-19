# Displacement interpolation between two balanced halves

## 0. Scope and verdict

Let \(\mu\) be an isotropic log-concave probability on its affine
support \(F\simeq\mathbb R^d\), let \(E\) be Borel with
\(\mu(E)=1/2\), and set

\[
 \mu_+=2\,1_E\mu,\qquad \mu_-=2\,1_{E^c}\mu.             \tag{0.1}
\]

If \(T\) is the Brenier map from \(\mu_+\) to \(\mu_-\) and
\(F_t=(1-t)I+tT\), then the displacement interpolant
\(\nu_t=(F_t)_\#\mu_+\) satisfies the exact measure domination

\[
 \boxed{\nu_t\le 2\mu\quad(0\le t\le1).}                 \tag{0.2}
\]

Consequently \(\operatorname {Ent}_\mu(\nu_t)\le\log2\).  More is
true.  If \(A(x)=D_aT(x)\) denotes the absolutely continuous
Alexandrov derivative, then for \(0<t<1\)

\[
\begin{split}
 \log2-\operatorname {Ent}_\mu(\nu_t)
  =\int\bigg[&(1-t)V(x)+tV(Tx)-V(F_tx)\\
  &+\log\det((1-t)I+tA)-t\log\det A\bigg]\,d\mu_+(x).
\end{split}                                                   \tag{0.3}
\]

Both terms in brackets are nonnegative.  At the midpoint this becomes

\[
\boxed{
 \log2-\operatorname {Ent}_\mu(\nu_{1/2})
 =\mathbb E_{\mu_+}\left[
 {V(X)+V(TX)\over2}-V\left({X+TX\over2}\right)
 +\log\det {I+A\over2A^{1/2}}\right].}                   \tag{0.4}
\]

Thus the expected potential midpoint deficit plus the expected matrix
arithmetic--geometric-mean deficit is at most \(\log2\).

There is also an exact covariance identity.  With

\[
 D=T(X)-X,\qquad K=\mathbb E[D\otimes D],                 \tag{0.5}
\]

where \(K\) is the **uncentered** displacement second-moment matrix,

\[
 \boxed{\operatorname {Cov}(\nu_{1/2})=I-{1\over4}K.}    \tag{0.6}
\]

In particular \(K\preceq4I\), but this gives only
\(W_2^2(\mu_+,\mu_-)=\operatorname {tr}K\le4d\).

The entropy identity does not by itself upgrade this trace bound to a
universal constant.  It gives a dimension-free bound on a normalized
derivative strain, but it is blind to componentwise translations and to
the singular part of the Brenier Hessian.  The exact zero-gap objects are
locally translations when the Hessian is absolutely continuous, and more
generally may be separated by singular jumps and by target gaps on which
the displacement potential is not convex.
An explicit isotropic power-diagram example in Section 8 has zero entropy
gap and nonzero transport cost.

Define

\[
 \mathcal W_2(\mu)=
 \sup_{\mu(E)=1/2}W_2(2\,1_E\mu,2\,1_{E^c}\mu).           \tag{0.7}
\]

For log-concave probabilities, a universal bound on \(\mathcal W_2\) is
equivalent, up to universal constants, to KLS.  Precisely,

\[
 \mathcal W_2(\mu)\le4\sqrt{C_P(\mu)},
 \qquad
 D_1(\mu)\le {1\over2}\mathcal W_2(\mu),                \tag{0.8}
\]

and E. Milman's equivalence theorem gives
\(D_1(\mu)\asymp1/\psi_\mu\asymp\sqrt{C_P(\mu)}\) for
log-concave measures on their affine supports.  Therefore the missing
dimension-free displacement estimate is not a weaker consequence of
(0.2)--(0.6); it is another exact KLS target.  No counterexample to the
universal estimate is produced here.  Instead, the new usable output is
the strain/curvature versus singular-translation alternative made
quantitative in Section 5.

All statements below are intrinsic to the affine support.  A point mass
is excluded.  Thus lower-dimensional support causes no degeneracy: choose
an orthonormal affine identification with \(\mathbb R^d\), perform all
calculations there, and transport the conclusions back by the isometry.

## 1. Brenier and Jacobian facts with nonsmooth data

Write

\[
 d\mu(x)=e^{-V(x)}dx                                      \tag{1.1}
\]

on \(F\), with \(V:F\to(-\infty,+\infty]\) convex.  Both
\(\mu_+\) and \(\mu_-\) are absolutely continuous and have finite second
moment.  Brenier's theorem therefore supplies a convex function \(\phi\)
such that \(T=\nabla\phi\) is the unique optimal map
\(\mu_+\)-almost everywhere.

At \(\mu_+\)-almost every \(x\), all of the following hold:

1. \(T\) has an Alexandrov derivative
   \(A=D_aT=D_a^2\phi\succeq0\);
2. \(y=T(x)\) is a density point of the target;
3. the Monge--Ampere change-of-variables identity is

   \[
    2e^{-V(x)}=2e^{-V(y)}\det A.                          \tag{1.2}
   \]

Thus \(A\succ0\) and

\[
 \log\det A=V(y)-V(x).                                   \tag{1.3}
\]

For \(0<t<1\), monotonicity of \(T\) gives

\[
 \langle F_t(x)-F_t(x'),x-x'\rangle
 \ge(1-t)|x-x'|^2.                                       \tag{1.4}
\]

Hence \(F_t\) is injective on a full source-measure set and its inverse
on its range is \((1-t)^{-1}\)-Lipschitz.  In particular \(\nu_t\) is
absolutely continuous.  At the preceding regular points,

\[
 \rho_t(F_t x)\det((1-t)I+tA)=2e^{-V(x)},                \tag{1.5}
\]

where \(\rho_t\) is the Lebesgue density of \(\nu_t\).

These assertions do **not** assume that the distributional Hessian
\(D^2\phi\) is absolutely continuous.  To see where its singular part
goes, apply the area formula to the Lipschitz inverse of \(F_t\).  A set
on which \(F_t\) has a jump or singular expansion can have positive
Lebesgue image, but it has zero \(\nu_t\)-mass when its preimage has zero
\(\mu_+\)-mass.  Formula (1.5) holds at source-almost every regular point;
on the skipped image the density is zero.  Thus singular Hessian mass can
create holes in the interpolant, but it cannot increase its density and
it contributes no hidden positive term to (0.3).  Section 5 gives the
exact zero-gap interpretation of this phenomenon.

No smoothing of \(V\), \(E\), or the support is being used here.
Extended values \(V=+\infty\) are harmless: almost every endpoint has
finite \(V\), and convexity makes every interior point of its transport
chord finite.

## 2. Pointwise domination and the full entropy identity

Let \(z=F_t(x)\) and \(B_t=(1-t)I+tA\).  Scalar weighted AM--GM applied
to the eigenvalues of \(A\) gives

\[
 \det B_t\ge(\det A)^t.                                  \tag{2.1}
\]

Combining (1.3), (1.5), and convexity of \(V\),

\[
\begin{split}
 \log{\rho_t(z)\over e^{-V(z)}}
 &=\log2-V(x)-\log\det B_t+V(z)\\
 &\le\log2+V(z)-(1-t)V(x)-tV(Tx)\\
 &\le\log2.
\end{split}                                               \tag{2.2}
\]

At points outside the regular image, \(\rho_t=0\).  Hence (2.2) proves
the measure inequality \(\nu_t\le2\mu\).  The endpoints \(t=0,1\) are
immediate from (0.1).

In particular, if \(h_t=d\nu_t/d\mu\), then
\(0\le h_t\le2\), \(\int h_t\,d\mu=1\), and

\[
 0\le\operatorname {Ent}_\mu(\nu_t)
   =\int h_t\log h_t\,d\mu\le\log2.                     \tag{2.3}
\]

Changing variables through \(F_t\) in the entropy gives the sharper
identity

\[
\begin{split}
 \operatorname {Ent}_\mu(\nu_t)
 &=\log2+\mathbb E_{\mu_+}
     [V(F_tX)-V(X)-\log\det B_t]\\
 &=\log2-\mathbb E_{\mu_+}[\Delta_{V,t}(X,T X)
                            +\Delta_{A,t}(A)].
\end{split}                                               \tag{2.4}
\]

where

\[
\begin{split}
 \Delta_{V,t}(x,y)
   &=(1-t)V(x)+tV(y)-V((1-t)x+ty),\\
 \Delta_{A,t}(A)
   &=\log\det((1-t)I+tA)-t\log\det A.
\end{split}                                               \tag{2.5}
\]

Both are nonnegative.  This proves (0.3), including equality rather than
only an inequality.  At \(t=1/2\),

\[
 \Delta_{A,1/2}(A)
 =\log\det\left({I+A\over2A^{1/2}}\right),               \tag{2.6}
\]

which proves (0.4).  In particular,

\[
 \mathbb E\Delta_{V,1/2}\le\log2,
 \qquad
 \mathbb E\log\det\left({I+A\over2A^{1/2}}\right)
 \le\log2.                                               \tag{2.7}
\]

There is a useful nonsmooth curvature reading of the first term.  On a
transport chord put \(g(s)=V(x+s(y-x))\).  The distributional second
derivative \(g''\) is a nonnegative measure and

\[
 {g(0)+g(1)\over2}-g(1/2)
 ={1\over2}\int_0^1\min(s,1-s)\,dg''(s).                \tag{2.8}
\]

Thus kinks and singular curvature of \(V\) are included, not discarded.
If \(V\) is twice differentiable along the chord, the right side is

\[
 {1\over2}\int_0^1\min(s,1-s)
 \langle\nabla^2V(x+sD)D,D\rangle\,ds.                   \tag{2.9}
\]

In particular, if \(\nabla^2V\succeq\kappa I\) in the convex sense,

\[
 \Delta_{V,1/2}(x,y)\ge{\kappa\over8}|x-y|^2,
 \qquad
 W_2^2(\mu_+,\mu_-)\le{8\log2\over\kappa}.              \tag{2.10}
\]

This recovers a complete dimension-free answer for strongly log-concave
targets, but \(\kappa=0\) for the cube, simplex, uniform ball, and affine
exponential examples below.

## 3. The midpoint covariance identity, including barycenters

Let \(X\sim\mu_+\), \(Y=T(X)\sim\mu_-\), and \(Z=(X+Y)/2\).
Write

\[
 m_+=\mathbb E X,\qquad m_-=\mathbb E Y.                 \tag{3.1}
\]

Because \(\mu=(\mu_++\mu_-)/2\) and \(\mu\) is centered,

\[
 m_-=-m_+.                                                \tag{3.2}
\]

Moreover, isotropy gives the **raw** second-moment identity

\[
 {1\over2}\mathbb E[X\otimes X]
 +{1\over2}\mathbb E[Y\otimes Y]=I.                    \tag{3.3}
\]

For \(D=Y-X\) and \(K=\mathbb E[D\otimes D]\), direct expansion yields

\[
\begin{split}
 \mathbb E[Z\otimes Z]
 &={1\over4}\mathbb E[X\otimes X+Y\otimes Y
                       +X\otimes Y+Y\otimes X]\\
 &=I-{1\over4}K.
\end{split}                                               \tag{3.4}
\]

Equation (3.2) says \(\mathbb EZ=0\), proving (0.6).  If one insists on
the centered displacement covariance

\[
 K_c=\operatorname {Cov}(D),                              \tag{3.5}
\]

then \(\mathbb ED=-2m_+\) and the identity is

\[
 \boxed{
 \operatorname {Cov}(\nu_{1/2})
 =I-m_+\otimes m_+-{1\over4}K_c.}                        \tag{3.6}
\]

Also

\[
 I={1\over2}\operatorname {Cov}(\mu_+)
  +{1\over2}\operatorname {Cov}(\mu_-)
  +m_+\otimes m_+,                                       \tag{3.7}
\]

so \(|m_+|\le1\).  Positivity in (3.4) gives

\[
 K\preceq4I,\qquad
 4|m_+|^2\le W_2^2(\mu_+,\mu_-)=\operatorname {tr}K
 \le4d.                                                   \tag{3.8}
\]

Entropy and covariance alone cannot improve the trace conclusion.  For
example, if \(\gamma_d\) is standard Gaussian and \(\eta_d\) is
\(\gamma_d\) conditioned on its inner radial half, then
\(d\eta_d/d\gamma_d\le2\),
\(\operatorname {Ent}_{\gamma_d}(\eta_d)=\log2\), and

\[
 \operatorname {Cov}(\eta_d)
  =(1-c_d/\sqrt d)I,\qquad 0<c\le c_d\le C.              \tag{3.9}
\]

Thus the formal matrix \(4(I-\operatorname {Cov}\eta_d)\) has trace of
order \(\sqrt d\).  This is not asserted to be the midpoint of a bad
half-to-half Brenier geodesic; it isolates exactly why the facts
\(\nu_{1/2}\le2\mu\), entropy at most \(\log2\), and (0.6) alone cannot
close the argument.

## 4. A balanced-half \(W_2\) bound is KLS-equivalent

This section makes the remaining theorem-strength step precise.

### 4.1 Poincare implies a half-to-half transport bound

Suppose \(\mu\) satisfies

\[
 \operatorname {Var}_\mu(g)
 \le C_P\int|\nabla g|^2d\mu.                             \tag{4.1}
\]

For any probability \(\sigma=f\mu\) with
\(\chi^2(\sigma\mid\mu)=\int(f-1)^2d\mu<\infty\),

\[
 W_2(\sigma,\mu)
 \le2\sqrt{C_P\chi^2(\sigma\mid\mu)}.                  \tag{4.2}
\]

Here is a proof that also handles densities which vanish.  In the
mean-zero Sobolev closure, solve weakly

\[
 \int\langle\nabla u,\nabla h\rangle d\mu
 =\int(f-1)h\,d\mu.                                      \tag{4.3}
\]

Lax--Milgram and (4.1) give

\[
 \int|\nabla u|^2d\mu\le C_P\chi^2(\sigma\mid\mu).       \tag{4.4}
\]

Along the linear density path \(\rho_s=1+s(f-1)\), use flux
\(\rho_sv_s=\nabla u\).  This solves the continuity equation weakly.
Since \(\rho_s\ge1-s\), its metric speed is at most

\[
 \left(\int{|\nabla u|^2\over\rho_s}d\mu\right)^{1/2}
 \le{\|\nabla u\|_2\over\sqrt{1-s}}.                    \tag{4.5}
\]

Integrating the metric speed gives length at most \(2\|\nabla u\|_2\).
For complete rigor when \(f=0\) on a set, stop at \(s=1-\varepsilon\),
reparametrize by arc length to obtain finite Benamou--Brenier action, and
let \(\varepsilon\downarrow0\).  Finite second moments and lower
semicontinuity of \(W_2\) give (4.2).

For \(f=2\,1_E\), \(\chi^2=1\).  Applying (4.2) to both halves and using
the triangle inequality proves

\[
 \mathcal W_2(\mu)\le4\sqrt{C_P(\mu)}.                  \tag{4.6}
\]

### 4.2 A half-to-half transport bound implies T3

Let \(f:F\to\mathbb R\) be 1-Lipschitz and let \(m\) be a median.
Because a nonpoint log-concave probability is nonatomic on its affine
support, choose a Borel subset of \(\{f=m\}\), if necessary, so that a
set \(E\) of measure \(1/2\) satisfies

\[
 \{f<m\}\subset E\subset\{f\le m\}.                    \tag{4.7}
\]

Put \(a=\int f\,d\mu_+\) and \(b=\int f\,d\mu_-\).  Then

\[
 \int|f-m|d\mu
 ={1\over2}\big[(m-a)+(b-m)\big]
 ={b-a\over2}.                                           \tag{4.8}
\]

Kantorovich duality and \(W_1\le W_2\) imply

\[
 b-a\le W_1(\mu_+,\mu_-)
       \le W_2(\mu_+,\mu_-),                             \tag{4.9}
\]

so

\[
 D_1(\mu)\le{1\over2}\mathcal W_2(\mu).                  \tag{4.10}
\]

E. Milman's equivalence theorem applies to every log-concave probability
on a finite-dimensional Euclidean affine support: its Cheeger constant,
Poincare constant, exponential concentration scale, and median
first-moment Lipschitz scale are equivalent with numerical constants.
Applied intrinsically on \(F\), it gives

\[
 {c\over\psi_\mu}\le D_1(\mu)\le{C\over\psi_\mu},
 \qquad
 {c\over\psi_\mu}\le\sqrt{C_P(\mu)}
                    \le{C\over\psi_\mu}.                \tag{4.11}
\]

Together, (4.6), (4.10), and (4.11) show

\[
 \mathcal W_2(\mu)\asymp D_1(\mu)
 \asymp\sqrt{C_P(\mu)}\asymp{1\over\psi_\mu}           \tag{4.12}
\]

up to universal constants.  Therefore proving
\(\operatorname {tr}K\le C\) for every isotropic log-concave \(\mu\)
and every balanced \(E\) is equivalent to KLS.  It is not a consequence
that may be inserted after (0.4) without proof.

Notice that (4.12) concerns only bounded-density balanced halves.  It is
not a \(T_2\) inequality for arbitrary changes of measure and does not
imply a log-Sobolev inequality.  In particular, it is compatible with the
failure of \(T_2\) and log-Sobolev inequalities for exponential tails.

## 5. What the determinant gap actually controls

Let

\[
 H=(A-I)(A+I)^{-1}.                                      \tag{5.1}
\]

If \(\lambda_i\) are the eigenvalues of \(A\), the eigenvalues of \(H\)
are \(h_i=(\lambda_i-1)/(\lambda_i+1)\), and

\[
 \log{1+\lambda_i\over2\sqrt{\lambda_i}}
 =-{1\over2}\log(1-h_i^2)
 =\log\cosh\left({\log\lambda_i\over2}\right)
 \ge{1\over2}h_i^2.                                     \tag{5.2}
\]

Writing
\(G=\log2-\operatorname {Ent}_\mu(\nu_{1/2})\), (0.4)
therefore gives the exact dimension-free strain budget

\[
 \boxed{\mathbb E\|H\|_{HS}^2\le2G\le2\log2.}          \tag{5.3}
\]

For every \(R>1\),

\[
 \mathbb E\#\{i:\lambda_i\notin[R^{-1},R]\}
 \le {G\over\log\cosh((\log R)/2)},                    \tag{5.4}
\]

and, on the remaining eigenvalues,

\[
 \mathbb E\sum_{R^{-1}\le\lambda_i\le R}(\lambda_i-1)^2
 \le2(R+1)^2G.                                           \tag{5.5}
\]

Thus only a dimension-free expected number of derivative directions can
be substantially distorted.  Simultaneously,

\[
 \mathbb E\left[{V(X)+V(Y)\over2}-V(Z)\right]\le G.     \tag{5.6}
\]

Equations (5.3)--(5.6) are a genuine structural conclusion.  They do not
control \(\mathbb E|T(X)-X|^2\): derivatives do not determine integration
constants on disconnected cells, and (5.3) controls the Cayley strain,
not the unbounded quantity \(|A-I|\).

There is nevertheless a useful rank consequence when these estimates are
combined with (3.8).  If \(S>0\) and

\[
 S=W_2^2(\mu_+,\mu_-)=\operatorname {tr}K,
\]

then

\[
 {\operatorname {tr}K\over\|K\|_{\mathrm{op}}}\ge {S\over4}. \tag{5.6a}
\]

Thus a sequence with \(W_2\to\infty\) cannot be a single long
translation: its displacement second moment has effective rank tending
to infinity, while the **total** normalized absolutely continuous strain
in (5.3) remains bounded.  More explicitly,
\(K_c=K-4m_+\otimes m_+\preceq K\), so
\(\|K_c\|_{\mathrm{op}}\le4\) and

\[
 \operatorname {tr}K_c=S-4|m_+|^2\ge S-4.               \tag{5.6b}
\]

After removal of the rank-one barycenter part, a hypothetical large cost
therefore occupies at least \((S-4)/4\) effective displacement directions,
even though their total normalized local strain is bounded by (5.3).
This excludes a low-rank affine-deformation explanation.  Controlling the
integration constants and singular parts across those many low-strain
directions is precisely the missing global step.  The rank statement is
an exact counterexample constraint, not a bound on \(S\).

The exact zero-gap statement is especially clear.  If \(G=0\), then

\[
 A=I\quad\mu_+\text{-a.e.},
 \qquad
 V\text{ is affine on }[x,T(x)]\quad\mu_+\text{-a.e.}    \tag{5.7}
\]

The second assertion follows from equality at the midpoint for a convex
one-dimensional restriction.  On a connected open source region on which
the Brenier Hessian is absolutely continuous, the displacement \(T-I\)
is constant, so the map is a translation.  Globally this conclusion is
false without a singular-Hessian and source-gap audit.

Indeed, write distributionally

\[
 D^2\phi=A\,dx+(D^2\phi)^s.                              \tag{5.8}
\]

When \(A=I\) Lebesgue-almost everywhere on an entire open region (not
merely on its source half), one may write

\[
 \phi(x)={|x|^2\over2}+h(x),\qquad D^2h=(D^2\phi)^s\succeq0. \tag{5.9}
\]

Thus \(h\) is convex with purely singular Hessian in that stronger
subcase.  Its gradient can be piecewise constant with jumps on
power-diagram walls, or in one dimension it can have a singular continuous
Cantor part.  For a general zero-gap half-transport, however, \(A=I\) is
known only on \(E\).  On target gaps the convex Brenier potential may have
Hessian below \(I\), so \(h=\phi-|x|^2/2\) need not be convex.  These
features are all invisible to (5.3).  The correct alternative is therefore

* controlled absolutely continuous normalized strain and controlled
  potential bending; or
* large transport stored in cellwise translation constants and in the
  singular monotone part of the Brenier map.

To turn this alternative into KLS one would need a dimension-free global
connectivity/coercivity theorem controlling those constants and singular
parts from isotropy and log-concavity.  By Section 4, a theorem strong
enough to bound their total displacement for every balanced cut is itself
KLS-equivalent.

## 6. Uniform ball: inner half versus outer half

Let \(\mu\) be uniform on the isotropic ball \(B_R^d\), where
\(R=\sqrt{d+2}\), and put

\[
 E=B_{R2^{-1/d}}^d.                                      \tag{6.1}
\]

The Brenier map is radial.  If \(U\) is uniform on \((0,1)\), the coupled
radii are

\[
 r(U)=R(U/2)^{1/d},
 \qquad
 s(U)=R((1+U)/2)^{1/d}.                                  \tag{6.2}
\]

Hence

\[
 W_2^2(\mu_+,\mu_-)
 =R^2\int_0^1
 \left[((1+u)/2)^{1/d}-(u/2)^{1/d}\right]^2du.           \tag{6.3}
\]

Since

\[
 ((1+u)/2)^{1/d}-(u/2)^{1/d}
 \le {1\over d}\log(1+1/u),                             \tag{6.4}
\]

and \(\int_0^1\log^2(2/u)du<4\),

\[
 W_2^2(\mu_+,\mu_-)\le {12\over d}.                    \tag{6.5}
\]

For the radial derivative, the tangential eigenvalue is \(s/r\) with
multiplicity \(d-1\), the radial eigenvalue is \((r/s)^{d-1}\), and
\(\det A=1\).  The potential deficit vanishes because the density is
uniform and the ball is convex.  The full entropy gap is the integrable
matrix gap in (0.4), including the very anisotropic behavior near the
origin.  Thus this model passes the singular-endpoint and large-eigenvalue
tests without any dimension loss.

## 7. Gaussian halfspace and Gaussian radial cut

Let \(\mu=\gamma_d\).  Since \(V(x)=|x|^2/2+\text{constant}\), every
transport chord satisfies

\[
 {V(x)+V(y)\over2}-V((x+y)/2)={|x-y|^2\over8}.           \tag{7.1}
\]

Therefore (0.4) gives, for **every** balanced Gaussian cut,

\[
 W_2^2(\mu_+,\mu_-)\le8\log2.                           \tag{7.2}
\]

For the halfspace \(E=\{x_1\le0\}\), the other coordinates are fixed
and the first-coordinate monotone map is described by

\[
 x_1=\Phi^{-1}(u/2),\qquad
 T_1(x)=\Phi^{-1}((1+u)/2),\qquad 0<u<1.                 \tag{7.3}
\]

Thus

\[
 W_2^2=\int_0^1
 [\Phi^{-1}((1+u)/2)-\Phi^{-1}(u/2)]^2du,                \tag{7.4}
\]

and the barycenter bound gives \(W_2^2\ge8/\pi\).  Both
the potential term \(W_2^2/8\) and the one-dimensional determinant term
appear in (0.4).

For the radial cut \(E=\{|x|\le r_d\}\), where \(r_d\) is the median of
the \(\chi_d\) law with distribution function \(F_d\), the radial Brenier
map gives

\[
 W_2^2=\int_0^1
 [F_d^{-1}((1+u)/2)-F_d^{-1}(u/2)]^2du\le8\log2.          \tag{7.5}
\]

Thus the halfspace and radial models instantiate the same exact identity
through incompatible derivative geometries.

## 8. Cube halfspace and an exact power-diagram jump model

### 8.1 Cube halfspace

Let \(\mu\) be uniform on \([ -\sqrt3,\sqrt3]^d\) and
\(E=\{x_1\le0\}\).  The Brenier map is

\[
 T(x)=x+\sqrt3\,e_1.                                     \tag{8.1}
\]

Hence

\[
 W_2^2=3,\qquad A=I,\qquad \Delta_V=\Delta_A=0.           \tag{8.2}
\]

The midpoint law is exactly \(2\mu\) restricted to the central slab
\(\{|x_1|\le\sqrt3/2\}\), and its relative entropy is \(\log2\).
This is the simplest proof that the entropy gap does not see translation
amplitude.  Isotropy controls the single translation through the
barycenter, but gives no analogous estimate for many canceling cells.

### 8.2 A two-cell power diagram with singular Brenier Hessian

Start with the uniform probability on \([0,6]\) and let

\[
 E=[0,1]\cup[2,4],\qquad E^c=[1,2]\cup[4,6].             \tag{8.3}
\]

The increasing rearrangement is

\[
 T(x)=
 \begin{cases}
  x+1,&x\in[0,1],\\
  x+2,&x\in[2,4].
 \end{cases}                                             \tag{8.4}
\]

It is the restriction of the gradient of

\[
 \phi(x)={x^2\over2}+\max\{x+b_1,2x+b_2\},              \tag{8.5}
\]

with \(b_1-b_2=3/2\).  Thus the source cells are the
one-dimensional Laguerre, or power-diagram, cells cut by \(x=3/2\).
The absolutely continuous derivative is \(T'=1\), while the monotone
extension has a positive jump and its distributional derivative has a
Dirac mass.  That singular mass is not present in \(A=D_aT\) and must not
be silently discarded.

The midpoint images are \([1/2,3/2]\) and \([3,5]\), so
\(\nu_{1/2}=2\mu\) on their union.  Hence

\[
 \operatorname {Ent}_\mu(\nu_{1/2})=\log2,
 \qquad \Delta_V=\Delta_A=0.                             \tag{8.6}
\]

The two source cells have conditional weights \(1/3\) and \(2/3\), so

\[
 W_2^2={1\over3}\,1^2+{2\over3}\,2^2=3.                \tag{8.7}
\]

Centering at \(3\) and scaling by \(1/\sqrt3\) makes the uniform measure
isotropic and turns (8.7) into

\[
 W_2^2=1.                                                 \tag{8.8}
\]

In higher dimensions one may tensor this construction with any isotropic
convex factor and let the cut depend only on the displayed coordinate.
The model is an exact, log-concave, isotropic realization of the
piecewise-translation/singular-jump branch.  It is not a counterexample
to a universal \(W_2\) bound; it proves that any such bound needs a global
cell-connectivity argument absent from the entropy gap.

In fact this blindness is complete in one dimension.  For the uniform
law on any interval and any Borel half-partition, the increasing
rearrangement between the two restrictions satisfies \(T'=1\) at almost
every source density point, because the source and target Lebesgue
densities are equal.  Thus both midpoint deficits vanish for every such
partition; all dependence on the arrangement of the two halves is
carried by jumps and singular monotone increments.  Isotropic
one-dimensional geometry still bounds the resulting cost, but the
entropy gap records none of that bound.

This also prevents a common false inference.  On \(K=[0,4]\), take

\[
 E=[0,1]\cup[3,4],\qquad E^c=[1,3].
\]

The increasing Brenier map is \(T(x)=x+1\) on the first source interval
and \(T(x)=x-1\) on the second.  It has \(T'=1\) source-almost everywhere,
but \(T-I\) decreases from \(+1\) to \(-1\).  The convex Brenier potential
has the monotone extension \(T=2\) across the target gap \([1,3]\);
therefore, for \(h=\phi-|x|^2/2\), one has \(h''=-1\) on that gap.
Consequently \(DT=I\) on the source does **not** imply that the
displacement is the gradient of a globally convex function with purely
singular Hessian.  That power-diagram description is a useful stronger
subclass, not the whole zero-gap branch.

More generally, if

\[
 \phi(x)={|x|^2\over2}+\max_i\{a_i\cdot x+b_i\},         \tag{8.9}
\]

then \(T=x+a_i\) on the corresponding power cell.  Whenever source cells
\(P_i\) and translated cells \(P_i+a_i\) form the two halves of a uniform
convex body, the determinant and potential gaps vanish, while

\[
 W_2^2=\sum_i\mu_+(P_i)|a_i|^2.                          \tag{8.10}
\]

All information in (8.10) is stored in translation constants and jump
walls, not in the absolutely continuous Jacobian.

## 9. Isotropic simplex cap

Let \(v_0,\ldots,v_d\) be the vertices of a regular isotropic simplex,
normalized by

\[
 |v_i|^2=d(d+2),\qquad v_i\cdot v_j=-(d+2)\quad(i\ne j). \tag{9.1}
\]

Write \(X=\sum_{i=0}^dp_iv_i\), where
\((p_0,\ldots,p_d)\) is uniform Dirichlet.  Put

\[
 b=2^{-1/d},\qquad a=1-b,
 \qquad E=\{p_0\ge a\}.                                  \tag{9.2}
\]

Since \(p_0\sim\operatorname {Beta}(1,d)\), this cap has mass \(1/2\).
Its conditional barycenter is

\[
 m_+=(1-b)v_0,qquad m_-=-m_+.                            \tag{9.3}
\]

Consequently

\[
 W_2^2(\mu_+,\mu_-)
 \ge4d(d+2)(1-2^{-1/d})^2
 \longrightarrow4(\log2)^2.                             \tag{9.4}
\]

There is also a direct universal upper bound, independent of any KLS
input.  Conditional on \(p_0=s\), write

\[
 X=s v_0+(1-s)W,
\]

where \(W\) is uniform on the opposite facet and independent of \(s\).
Couple the upper and lower conditional laws of \(s\) by their common
quantile and use the same \(W\).  Then

\[
 \mathbb E|v_0-W|^2=(d+2)(d+3).                          \tag{9.5}
\]

The upper conditional variance of \(s\) is
\(b^2d/((d+1)^2(d+2))\le d^{-2}\); the lower conditional support has
length \(a\le(\log2)/d\); and the two conditional means differ by
\(2d(1-b)/(d+1)\le2(\log2)/d\).  Therefore the displayed coupling has
cost at most \(51\), and

\[
 4d(d+2)(1-2^{-1/d})^2\le W_2^2(\mu_+,\mu_-)\le51.       \tag{9.6}
\]

The simplex potential is constant on its convex support, so the potential
midpoint deficit is zero.  The Brenier map need not be the elementary
coupling just used, but its full normalized derivative strain still obeys
(5.3).  Thus neither the acute simplex geometry nor the cap limit causes
any hidden dimension dependence in the exact identities.

## 10. Product exponential cut by the maximum

Let \(Z_1,\ldots,Z_d\) be independent mean-one exponentials and put
\(X_i=Z_i-1\).  Their product law is isotropic and log-concave, with

\[
 V(x)=\sum_{i=1}^d(x_i+1)
 \quad\text{on }[-1,\infty)^d.                            \tag{10.1}
\]

Let \(M=\max_iX_i\) and choose its median

\[
 m_d=-1-\log(1-2^{-1/d}),
 \qquad E=\{M\le m_d\}.                                  \tag{10.2}
\]

The function \(M\) is 1-Lipschitz.  Hence

\[
 W_2(\mu_+,\mu_-)
 \ge W_1(\mu_+,\mu_-)
 \ge2\mathbb E|M-m_d|.                                  \tag{10.3}
\]

As \(d\to\infty\), \(M-(\log d-1)\) converges with uniform first
moment to a Gumbel variable \(G\), while
\(m_d-(\log d-1)\to-\log\log2\).  Thus the right side of (10.3)
converges to the finite positive constant

\[
 2\mathbb E|G+\log\log2|.                                \tag{10.4}
\]

For the upper bound, the one-dimensional exponential Poincare constant
is \(4\), and tensorization preserves it.  Equation (4.6) gives

\[
 W_2(\mu_+,\mu_-)\le8.                                   \tag{10.5}
\]

Most importantly for the entropy route, (10.1) is affine on every chord
inside the support.  Therefore

\[
 \Delta_{V,1/2}=0                                        \tag{10.6}
\]

for this highly nonsymmetric cut.  The whole entropy gap lies in the
matrix term.  This example confirms that one cannot covertly replace
log-concavity by positive curvature or a \(T_2\) hypothesis.

## 11. Final obstruction and reusable conclusions

The route supplies four exact, fully general facts:

1. every displacement interpolant between balanced halves obeys
   \(\nu_t\le2\mu\);
2. the entropy-convexity loss is exactly the sum of potential chord
   bending and matrix AM--GM strain in (0.3)--(0.4);
3. the midpoint covariance is exactly \(I-K/4\), including the
   barycenter correction (3.6);
4. the normalized absolutely continuous strain has total expected
   squared Hilbert--Schmidt norm at most \(2\log2\).

The remaining displacement can be stored in cellwise translations and
in the singular part of the Brenier Hessian.  Uniform bodies and affine
exponential potentials give no potential-curvature charge, and the cube
and power-diagram models show exact zero matrix charge for nontrivial
transport.  Isotropy bounds each displacement direction by (3.8), but
not the sum over directions.  If that sum is \(S\), its effective
displacement rank is at least \(S/4\) by (5.6a).

A dimension-free estimate of that sum for every balanced cut is exactly
the statement \(\mathcal W_2(\mu)\le C\), and Section 4 proves that this
is KLS-equivalent.  Thus the precise missing lemma is:

> **Singular-translation coercivity (KLS-equivalent).**  For every
> isotropic log-concave probability and every balanced Borel partition,
> the Brenier displacement satisfies
> \(\int|T(x)-x|^2d\mu_+(x)\le C\), including all
> componentwise translation constants and singular Hessian jumps.

It may not be assumed.  The entropy calculation reduces the possible
failure to a sharp global geometric branch, but does not exclude that
branch.  The mandatory ball, Gaussian, cube, simplex, product exponential,
and power-diagram tests all instantiate the formulas with universal
constants and expose no suppressed dimension parameter.
