# Clean-room audit of the scalar curvature--Hardy residual estimate

## 0. Verdict

The smooth scalar estimate in scalar_residual_curvature_hardy.md is
correct.  More precisely, let

\[
 r(s,z)=e^{-U(s,z)},\qquad \rho(s)=\int r(s,z)\,dz,
 \qquad q_s(z)=r(s,z)/\rho(s),
\]

and suppose that \(U\) is jointly convex and that each \(z\)-section is a
unit-Gaussian convolution.  If \(F_s^0\) is the centered conditional
Poisson velocity and

\[
 C_s=E_s|(F_s^0)'|^2,\qquad
 B_s^0=E_s[U_{zz}(F_s^0)^2],\qquad
 \sigma_s^2=\operatorname {Var}_{q_s}Z,
\]

then

\[
 \boxed{
 B_s^0\le 128\sqrt{\frac\pi2}\,\sigma_s C_s.}
 \tag{0.1}
\]

If the marginal \(S\sim\rho\) is centered with variance one and

\[
 \sigma_s\le L(1+|s|)
 \tag{0.2}
\]

for almost every \(s\), then the reinforced Prekopa identity gives

\[
 \boxed{
 \int \rho(s)\tau(s)^2B_s^0\,ds
 \le 128\sqrt{\frac\pi2}\,41L<6600L.}
 \tag{0.3}
\]

Here \(\tau\) is the canonical Stein kernel of \(\rho\).  The constant
chain is dimension-free and the cone example has the same local scaling.

The audit does **not** establish (0.2) for an arbitrary jointly
log-concave family.  Thus (0.3) is a conditional scalar theorem, not yet
the general scalar instance of (WK).  The two subclass claims in the
original note are valid after inserting the concave-scale lemma in
Section 4 below.

## 1. Curvature-tail calculation

Let \(q=e^{-V}\) be a positive \(C^2\) log-concave probability density on
\(\mathbb R\), with \(0\le V''\le1\).  Choose a mode \(x_0\), put
\(M=q(x_0)\), and let \(p=V'\).  If \(x\ge x_0\), then \(p(x)\ge0\), and

\[
 \begin{aligned}
 &V(t)-V(x)-\frac12\bigl(p(t)^2-p(x)^2\bigr)\\
 &\hspace{30mm}{}'=p(t)(1-p'(t))\ge0
 \qquad(t\ge x).
 \end{aligned}
\]

Consequently

\[
 \begin{aligned}
 \int_x^\infty V''(t)q(t)\,dt
 &\le q(x)\int_x^\infty
 e^{-[p(t)^2-p(x)^2]/2}p'(t)\,dt\\
 &\le q(x)e^{p(x)^2/2}
       \int_{p(x)}^\infty e^{-u^2/2}\,du\\
 &\le \sqrt{\frac\pi2}\,q(x).
 \end{aligned}
 \tag{1.1}
\]

The change of variable is valid even when \(p'\) vanishes; equivalently
one uses the Lebesgue--Stieltjes measure induced by the monotone score.
Reflection gives the other half-line.  In particular,

\[
 \int V''q\le2\sqrt{\frac\pi2}\,M.
 \tag{1.2}
\]

This verifies the only nonstandard tail estimate in the argument without
requiring strict convexity or surjectivity of the score.

## 2. Hardy constants and centering

For a nonnegative measure \(\nu\) on \([x_0,\infty)\), the one-sided
Muckenhoupt--Hardy estimate is

\[
 \int g^2\,d\nu\le4B_+\int(g')^2q,
 \qquad
 B_+=\sup_{x>x_0}\nu([x,\infty))
                  \int_{x_0}^xq(t)^{-1}\,dt,
 \tag{2.1}
\]

for \(g(x_0)=0\).  The factor \(4\) follows after the change of variable
\(y=\int_{x_0}^xq^{-1}\) from the classical Dirichlet Hardy inequality
\(\int g^2/y^2\le4\int(g')^2\).

Writing \(q(x)=M e^{-a}\), convexity of \(V-V(x_0)\) and the mass bound
\(\int_{x_0}^xq\le1\) give

\[
 q(x)\int_{x_0}^xq^{-1}\le M^{-1}.
 \tag{2.2}
\]

Combining (1.1), (2.1), and (2.2), and adding the two half-lines, gives

\[
 \int V''g^2q\le\frac{4\sqrt{\pi/2}}M\int(g')^2q.
 \tag{2.3}
\]

The monotonicity of the hazard and reverse hazard similarly gives

\[
 \int g^2q\le\frac4{M^2}\int(g')^2q.
 \tag{2.4}
\]

For centered \(f\), put \(g=f-f(x_0)\).  Then

\[
 |f(x_0)|^2=|E_qg|^2\le\frac4{M^2}\int(f')^2q.
\]

Using \(f^2\le2g^2+2f(x_0)^2\), (1.2), and (2.3), one obtains

\[
 \int V''f^2q\le\frac{24\sqrt{\pi/2}}M\int(f')^2q.
 \tag{2.5}
\]

Chebyshev puts mass at least \(3/4\) in an interval of length \(4\sigma\),
so \(M\ge3/(16\sigma)\).  This proves (0.1), since unit-Gaussian
convolution gives \(0\le U_{zz}\le1\).

All factors have been counted after adding the two half-lines; there is no
missing factor of two.

## 3. Direct nonsmooth Stein-curvature passage

The marginal part does not require an informal smoothing assertion.  Let
\(\rho=e^{-\Phi}\) be any one-dimensional centered variance-one
log-concave density on the interior \(I\) of its support.  Let

\[
 N(s)=\int_s^{\sup I}t\rho(t)\,dt,
 \qquad \tau=N/\rho.
\]

Choose a right-continuous version \(p\) of the monotone derivative of
\(\Phi\), and let \(dp=D^2\Phi\) be its Radon curvature measure.  Locally
on \(I\),

\[
 \tau'=\tau p-s,
 \qquad
 d\tau'=p\tau'\,ds+\tau\,dp-ds.
 \tag{3.1}
\]

For every compactly supported locally Lipschitz \(h\), Stieltjes
integration by parts therefore gives the exact identity

\[
 \boxed{
 \int_I h\tau^2\rho\,dp
 =E[h\tau]-E[h(\tau')^2]-E[h'\tau\tau'].}
 \tag{3.2}
\]

Indeed, substitute the second identity in (3.1) and integrate
\(\int h\tau\rho\,d\tau'\); the two terms containing \(p\tau\tau'\)
cancel.  Compact support removes all endpoint terms.

Taking expanding cutoffs, the cutoff-derivative term tends to zero by
Cauchy--Schwarz.  With \(h=1\), (3.2) yields

\[
 E(\tau')^2\le E\tau=1.
 \tag{3.3}
\]

For \(h=1+|s|\), use cutoffs and smooth approximations to \(|s|\), drop
the nonpositive middle term, and invoke \(E\tau^2\le400\).  This gives

\[
 \begin{aligned}
 \int_I(1+|s|)\tau^2\rho\,dp
 &\le E[(1+|S|)\tau]+E|\tau\tau'|\\
 &\le1+\sqrt{ES^2E\tau^2}
       +\sqrt{E\tau^2E(\tau')^2}\\
 &\le41.
 \end{aligned}
 \tag{3.4}
\]

Thus atoms of \(D^2\Phi\), hard marginal endpoints, and kinks of the
marginal potential do not create a missing term.  What still requires
care in a fully nonsmooth WFI theorem is the convergence/definition of
the conditional Poisson energies \(B_s^0,C_s\).  An upper bound for
\(B_s^0\) cannot be obtained merely from lower semicontinuity.  In the
regular setting, or whenever reinforced Prekopa is available as the
measure inequality

\[
 \rho(s)C_s\,ds\le \rho(s)\,D^2\Phi(ds),
 \tag{3.5}
\]

(0.1)--(0.2) and (3.4) prove (0.3) directly.

## 4. The two advertised subclasses

We use the following elementary fixed-dimensional lemma.

> **Concave-scale lemma.**  Let \(S\) have a centered variance-one
> one-dimensional log-concave density, let \(a\ge0\) be concave on the
> support interval, and suppose \(E[a(S)^p]\le K^p\) for some \(p\ge1\).
> Then
> \[
> a(s)\le C K(1+|s|)
> \tag{4.1}
> \]
> throughout the support, where \(C\) is universal.

To prove it, use the elementary one-dimensional isotropic density bounds
to choose a universal \(c>0\) such that the support contains \([-c,c]\),
\(\rho\ge c\) there.  Concavity and nonnegativity show that
\(a(t)\ge a(0)/2\) on \([-c/2,c/2]\), hence \(a(0)\le C K\).  If \(s>0\),
concavity between \((-c,a(-c))\) and \((s,a(s))\), evaluated at zero,
gives

\[
 a(0)\ge\frac{c}{s+c}a(s).
\]

The case \(s<0\) is symmetric.  This proves (4.1).

For a planar convex body with interval fibers of width \(w(s)\), the
function \(w\) is nonnegative and concave.  Conditional uniformity and
isotropy give

\[
 E w(S)^2/12=E\operatorname {Var}(Z\mid S)\le1.
\]

Apply (4.1) with \(p=2\).  After unit-Gaussian convolution,
\(\sigma_s^2=1+w(s)^2/12\), so (0.2) holds universally.

For a jointly log-concave scalar conditional-Gaussian family, the scalar
conditional covariance \(R(s)\) is nonnegative and concave; this is the
\(z^2\)-coefficient in the Schur-complement convexity condition.  In the
isotropic-plus-unit-noise normalization,

\[
 ER(S)\le\operatorname {Var}Z=2.
\]

Apply (4.1) with \(p=1\) to \(R\).  Then
\(\sigma_s=\sqrt{R(s)}\le C\sqrt{1+|s|}\le C(1+|s|)\).

Neither argument extends automatically to a general log-concave slice:
the effective width

\[
 \frac{\int r(s,z)\,dz}{\sup_z r(s,z)}
\]

is comparable to \(\sigma_s\), but it need not be concave as a function
of \(s\).  Proving (0.2), or proving directly

\[
 \int\rho\tau^2\sigma_s\,D^2\Phi(ds)\le C,
\]

remains the precise scalar geometric charge.
