# Scalar residual WFI from curvature--Hardy and conditional-scale growth

## 0. Result and scope

This note proves a genuinely global-in-\(s\) estimate for the residual
curvature term in transverse dimension one. It uses no conditional
Poincare inequality.

Let

\[
 r(s,z)=e^{-U(s,z)},\qquad
 \rho(s)=\int r(s,z)\,dz=e^{-\Phi(s)},\qquad
 q_s(z)=r(s,z)/\rho(s),
\]

where \(U\) is jointly convex and the \(z\)-variable has been convolved
with a unit Gaussian. Let \(S\sim\rho\) be centered with variance one,
and let \(\tau\) be its canonical Stein kernel. For the conditional
Poisson field, write

\[
 F_s=-m'(s)+F_s^0,\qquad E_sF_s^0=0,
\]

and put

\[
 H_s=U_{zz}(s,\cdot),\qquad
 C_s=E_s|(F_s^0)'|^2,\qquad
 B_s^0=E_s[H_s(F_s^0)^2].
\]

The Gaussian channel gives \(0\le H_s\le1\). Let
\(\sigma_s^2=\operatorname{Var}_{q_s}Z\).

> **Scalar residual theorem.** If, for some numerical \(L\ge1\),
> \[
>                         \sigma_s\le L(1+|s|)
>                         \quad\text{for a.e. }s,       \tag{0.1}
> \]
> then
> \[
> \boxed{\quad
>  \int\rho(s)\tau(s)^2B_s^0\,ds
>  \le 6600\,L .
> \quad}                                                \tag{0.2}
> \]

More intrinsically, without (0.1) one always has

\[
 \boxed{\quad
 \int\rho\tau^2B_s^0
 \le128\sqrt{\pi/2}
       \int\rho\tau^2\sigma_s C_s .
 \quad}                                                \tag{0.3}
\]

Thus the exact remaining scalar geometric charge is the
scale-weighted deformation energy. The cone in
work/gaussian_curvature_korn.md has
\(\sigma_s\asymp a\), \(C_s\lesssim a^{-2}\), and
\(B_s^0\gtrsim a^{-1}\); it shows that the factor \(\sigma_s\) in the
local estimate cannot be removed.

## 1. A curvature-tail lemma

Let \(q=e^{-V}\) be a positive \(C^2\) log-concave probability density on
\(\mathbb R\), and assume

\[
                         0\le V''\le1.                 \tag{1.1}
\]

Let \(x_0\) be a mode, \(M=q(x_0)=\|q\|_\infty\), and set
\(\kappa(dx)=V''(x)q(x)\,dx\). Then, with
\(c_G=\sqrt{\pi/2}\),

\[
 \kappa([x,\infty))\le c_Gq(x)\quad(x\ge x_0),\qquad
 \kappa((-\infty,x])\le c_Gq(x)\quad(x\le x_0).        \tag{1.2}
\]

To prove the right-hand inequality, put \(p=V'\) and
\(p_0=p(x)\ge0\). Since \(dp=V''dt\), monotonicity and (1.1) give

\[
 \begin{aligned}
 \frac{\kappa([x,\infty))}{q(x)}
 &=\int_x^\infty e^{-[V(t)-V(x)]}\,dp(t)\\
 &\le\int_{p_0}^\infty
       \exp\left[-\frac{u^2-p_0^2}{2}\right]du\\
 &\le\sqrt{\frac\pi2}.                                 \tag{1.3}
 \end{aligned}
\]

Indeed, while the score rises from \(p_0\) to \(u\), its slope is at most
one, so the accumulated potential is at least
\(\int_{p_0}^u v\,dv=(u^2-p_0^2)/2\); intervals on which the score is
constant add only nonnegative potential. The final Gaussian Mills ratio
is maximal at \(p_0=0\). Reflection proves the left-hand inequality.
Approximation covers non-strict curvature and a nonunique mode.

At the mode, (1.2) also gives

\[
                 \int V''q\le2c_GM.                   \tag{1.4}
\]

## 2. The scale-sensitive curvature--Hardy inequality

We need an elementary one-sided Hardy fact. If \(\mu(dx)=q(x)dx\),
\(\nu\) is a nonnegative measure on \([x_0,\infty)\), and \(g(x_0)=0\),
then

\[
 \int_{x_0}^\infty g^2\,d\nu
 \le4B_+\int_{x_0}^\infty(g')^2q\,dx,\qquad
 B_+=\sup_{x>x_0}\nu([x,\infty))
                    \int_{x_0}^x\frac{dt}{q(t)}.       \tag{2.1}
\]

For completeness, use the coordinate
\(y=\int_{x_0}^xq(t)^{-1}dt\). The Dirichlet energy becomes
\(\int|\partial_yg|^2dy\), while the pushed-forward tail is at most
\(B_+/y\). Integration by parts, Cauchy--Schwarz, and
\(\int_0^\infty g^2/y^2\le4\int_0^\infty(g')^2\)
give (2.1). The left-hand version is identical.

For \(x\ge x_0\), convexity gives

\[
 q(x)\int_{x_0}^x\frac{dt}{q(t)}\le\frac1M.            \tag{2.2}
\]

Indeed, normalize \(V(x_0)=0\), put \(a=V(x)\) and \(D=x-x_0\), and use
\(V(x_0+t)\le at/D\). Then

\[
 q(x)\int_{x_0}^xq^{-1}
 \le D\frac{1-e^{-a}}a
 \le\frac1M,
\]

because the mass of the same interval is at least
\(MD(1-e^{-a})/a\). The limiting case \(a=0\) is immediate.

Apply (2.1) with \(\nu=\kappa\). Equations (1.2) and (2.2) imply, for
every \(g\) with \(g(x_0)=0\),

\[
 \int g^2\,d\kappa
 \le\frac{4c_G}{M}\int(g')^2q.                         \tag{2.3}
\]

Here and below the two half-line estimates have been added.

We also require an anchored \(L^2(q)\) estimate. Log-concavity makes the
right hazard \(q(x)/q([x,\infty))\) increasing. Hence

\[
 \frac{q([x,\infty))}{q(x)}\le\frac1M\qquad(x\ge x_0).
\]

Together with (2.2), the Hardy criterion gives

\[
 \int g^2q\le\frac4{M^2}\int(g')^2q
 \quad\text{whenever }g(x_0)=0.                        \tag{2.4}
\]

The reverse hazard proves the other half-line statement.

Now let \(f\) be locally absolutely continuous and \(E_qf=0\), and set
\(g=f-f(x_0)\). Since \(f(x_0)=-E_qg\), (2.4) yields

\[
 |f(x_0)|^2\le\frac4{M^2}\int(f')^2q.                  \tag{2.5}
\]

Using \(f^2\le2g^2+2f(x_0)^2\), (1.4), (2.3), and (2.5), we obtain

\[
 \boxed{\quad
 \int V''f^2q
 \le\frac{24c_G}{M}\int(f')^2q .
 \quad}                                                \tag{2.6}
\]

If \(\sigma^2=\operatorname{Var}_qX\), Chebyshev puts mass at least
\(3/4\) in an interval of length \(4\sigma\), so
\(M\ge3/(16\sigma)\). Therefore

\[
 \boxed{\quad
 \int V''f^2q
 \le128\sqrt{\frac\pi2}\,
       \sigma\int(f')^2q .
 \quad}                                                \tag{2.7}
\]

The estimate is scale-sensitive rather than scale-invariant because the
unit Gaussian smoothing fixes the curvature cap in (1.1).

## 3. A Stein-weighted marginal-curvature moment

Let \(\rho=e^{-\Phi}\) be centered, variance one, and log-concave, and
let \(\tau\) be its canonical Stein kernel. In the smooth setting,

\[
 \tau'=\tau\Phi'-s,\qquad E\tau=1,\qquad
 E(\tau')^2\le1,\qquad E\tau^2\le400.                  \tag{3.1}
\]

For every nonnegative \(C^2\) function \(h\) for which the following
integrals are finite,

\[
 \boxed{\quad
 E[h\tau^2\Phi'']
 =E[h\tau]-E[h(\tau')^2]-E[h'\tau\tau'].
 \quad}                                                \tag{3.2}
\]

Indeed, differentiating (3.1) gives
\(\tau''=\tau'\Phi'+\tau\Phi''-1\). Multiply by \(h\tau\)
and integrate by parts against \(e^{-\Phi}\); the two
\(h\tau\tau'\Phi'\) terms cancel.

Approximate \(h(s)=1+|s|\) by smooth functions with
\(|h'|\le1\). Dropping the nonpositive middle term in (3.2) and using
Cauchy--Schwarz gives

\[
 \begin{aligned}
 E[(1+|S|)\tau^2\Phi'']
 &\le E[(1+|S|)\tau]+E|\tau\tau'|\\
 &\le 1+\sqrt{ES^2\,E\tau^2}
       +\sqrt{E\tau^2\,E(\tau')^2}\\
 &\le41.                                               \tag{3.3}
 \end{aligned}
\]

The same inequality passes to nonsmooth log-concave marginals by convex
regularization and lower semicontinuity of the nonnegative curvature
measure.

## 4. Proof of the residual theorem

Apply (2.7) conditionally with \(V=U(s,\cdot)+{\rm const}\) and
\(f=F_s^0\). Since \(E_sF_s^0=0\),

\[
 B_s^0\le128\sqrt{\frac\pi2}\,\sigma_sC_s.             \tag{4.1}
\]

This proves (0.3). Reinforced Prekopa gives pointwise

\[
 C_s+E_s\mathcal Q_s=\Phi''(s),
\qquad\text{hence}\qquad C_s\le\Phi''(s).              \tag{4.2}
\]

Under (0.1), equations (3.3), (4.1), and (4.2) yield

\[
 \boxed{\quad
 \int\rho\tau^2B_s^0
 \le5248\sqrt{\frac\pi2}\,L
 <6600L.
 \quad}                                                \tag{4.3}
\]

Condition (0.1) holds in the two principal scalar stress classes:

1. For uniform planar convex bodies, the slice half-width is concave and
   grows at most linearly after isotropic normalization; unit Gaussian
   convolution changes \(\sigma_s^2\) to \(\sigma_s^2+1\).
2. For conditional Gaussian fibers, the scalar covariance is concave
   under joint log-concavity, and its unconditional mean is bounded after
   isotropic normalization and unit smoothing.

It also holds for the isotropic exponential cone. There
\(\sigma_s^2=1+a^2/3\), \(a=s+\sqrt2\), and the lower bound
\(B_s^0\gtrsim a^{-1}\) together with \(C_s\lesssim a^{-2}\) shows that
(4.1) has the correct order.

## 5. Affine-shear audit and the remaining general issue

The exact affine-shear identity is

\[
 E_sD^2U[(1,m'(s)),(1,m'(s))]
 =\Phi''(s)+C_s+B_s^0.                                 \tag{5.1}
\]

It identifies \(B_s^0\) but does not control it: applying the full-space
Hodge formula to the affine-shear vector field reproduces (5.1)
tautologically. The curvature--Hardy estimate (4.1) is the additional
one-dimensional input.

One tempting scalar acceleration shortcut is invalid. If
\[
 A_s(z)=D^2U[(1,m'(s)),(1,m'(s))],
\]
the midpoint of
\((s-h,m(s-h)+y)\) and \((s+h,m(s+h)+y)\) is shifted by
\(\delta_h=\frac12m''(s)h^2+o(h^2)\). Taylor expansion gives
\[
 U+\tfrac12U_zm''h^2
 \le U+\tfrac12(A_s+U_zm'')h^2+o(h^2),
\]
so the \(U_zm''\) terms cancel and one recovers only \(A_s\ge0\).
There is no general inequality \(U_zm''\le A_s\).

The unresolved scalar question is therefore precise:

\[
 \text{Does every isotropic jointly log-concave scalar-slice family
 satisfy }\sigma_s\le C(1+|s|)?
 \tag{5.2}
\]

Alternatively it is enough to prove directly

\[
 \int\rho(s)\tau(s)^2\sigma_s\Phi''(s)\,ds\le C.       \tag{5.3}
\]

Both statements are purely scalar conditional-scale assertions. Neither
is replaced here by conditional Poincare, and neither is claimed without
proof.
