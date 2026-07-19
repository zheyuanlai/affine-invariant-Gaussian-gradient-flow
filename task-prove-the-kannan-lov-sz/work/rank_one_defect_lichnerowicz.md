# A rank-one-defect Lichnerowicz lemma

## 0. Statement and conclusion

Let \(u\in S^{n-1}\), let \(P=I-u\otimes u\), and let
\(\mu=e^{-V}dx\) be a log-concave probability on \(\mathbb R^n\).  Assume
in the sense of distributions that

\[
                         D^2V\succeq\kappa P               \tag{0.1}
\]

for some \(\kappa>0\).  Then

\[
\boxed{
 C_P(\mu)\le C\left(\kappa^{-1}
              +\operatorname {Var}_\mu\langle X,u\rangle\right),}   \tag{0.2}
\]

where \(C\) is numerical; the proof below gives \(C=96\).  In particular the assertion in the task is
true.  Centering and the extra hypothesis \(\operatorname {Cov}(\mu)\preceq
I\) are not needed for the lemma itself; under that hypothesis the second
term in (0.2) is at most one.

The proof has two ingredients.

1.  If \(\nu_t\) is the conditional law on \(u^\perp\), and
    \(e^{-U(t)}dt\) is the marginal law in the weak direction, there is a
    canonical continuity-equation velocity \(w_t\) satisfying

    \[
      \partial_t\nu_t+\operatorname {div}(\nu_tw_t)=0,
      \qquad
      \kappa\int|w_t|^2d\nu_t\le U''(t).                  \tag{0.3}
    \]

    The second inequality is a reinforced, infinitesimal
    Pr\'ekopa inequality.  It follows directly from Bochner's identity and
    the full block-matrix inequality (0.1).

2.  If a one-dimensional log-concave probability has potential \(U\),
    variance \(s^2\), and \(0\le W\le U''/\kappa\), then

    \[
      \operatorname {Var}(g)
      \le C(s^2+\kappa^{-1})
          \int {g'(t)^2\over1+W(t)}\,d\bar\mu(t).          \tag{0.4}
    \]

    This is an elementary weighted Hardy inequality.  The curvature part
    of its Muckenhoupt constant is universally bounded because

    \[
      \bar\mu([x,\infty))
      \int_m^x {U''(t)\over\bar\rho(t)}dt\le2.             \tag{0.5}
    \]

For \(g(t)=\int f(t,y)d\nu_t(y)\), the derivative term that usually
obstructs conditional/marginal decompositions is exactly

\[
 g'(t)=\int\bigl(\partial_tf+\nabla_yf\cdot w_t\bigr)d\nu_t. \tag{0.6}
\]

Consequently

\[
 {g'(t)^2\over1+\int|w_t|^2d\nu_t}
 \le\int|\nabla f|^2d\nu_t,                               \tag{0.7}
\]

which closes the two-scale estimate without discarding the derivative of
the conditional expectation.

## 1. A weighted one-dimensional lemma

We first prove the exact analytic input used for the marginal.

### Lemma 1.1 (curvature-weighted Poincare inequality)

Let \(I=(\alpha,\beta)\), where either endpoint may be infinite, and let

\[
 d\nu(t)=\rho(t)dt=Z^{-1}e^{-U(t)}1_I(t)dt                \tag{1.1}
\]

be a one-dimensional log-concave probability with \(U\in C^2(I)\) convex and
variance \(s^2\).  Let \(\kappa>0\), and let \(W:\mathbb R\to[0,\infty)\)
be measurable with

\[
                         \kappa W(t)\le U''(t).            \tag{1.2}
\]

Then every locally absolutely continuous \(g\in L^2(\nu)\) satisfies

\[
 \operatorname {Var}_\nu(g)
 \le 200(s^2+\kappa^{-1})
      \int {g'(t)^2\over1+W(t)}d\nu(t).                   \tag{1.3}
\]

#### Proof

Put

\[
 a_0(t)={\kappa\over\kappa+U''(t)}.                       \tag{1.4}
\]

By (1.2),

\[
                         {1\over1+W(t)}\ge a_0(t).         \tag{1.5}
\]

The measure has no atoms, so choose \(m\in I\) with
\(\nu((\alpha,m])=\nu([m,\beta))=1/2\).  The one-dimensional weighted
Muckenhoupt criterion states that the optimal constant in

\[
                         \operatorname {Var}_\nu g
 \le C_a\int a_0g'^2d\nu                                \tag{1.6}
\]

is at most \(4\max(B_+,B_-)\), where

\[
\begin{aligned}
 B_+&=\sup_{m<x<\beta}\nu([x,\beta))
              \int_m^x{dt\over a_0(t)\rho(t)},\\
 B_-&=\sup_{\alpha<x<m}\nu((\alpha,x])
              \int_x^m{dt\over a_0(t)\rho(t)}.           \tag{1.7}
\end{aligned}
\]

For completeness, the right-hand one-sided Hardy inequality is

\[
 \int_m^\beta (g-g(m))^2d\nu
 \le4B_+\int_m^\beta a_0g'^2d\nu,                       \tag{1.7a}
\]

and reflection gives its left-hand counterpart with \(B_-\).  Adding
them and using
\(\operatorname {Var}g=\inf_c\int(g-c)^2d\nu
\le\int(g-g(m))^2d\nu\) gives (1.6) with
\(C_a\le4\max(B_+,B_-)\).  This also records explicitly which version and
constant of the Muckenhoupt criterion is being used.

For the unweighted choice \(a_0=1\), the converse half of the same
criterion states
\[
 \max(B_+^{(0)},B_-^{(0)})\le C_P(\nu),                  \tag{1.7b}
\]
where \(B_\pm^{(0)}\) are (1.7) with \(a_0=1\).

Since

\[
 {1\over a_0\rho}={1\over\rho}+{1\over\kappa}{U''\over\rho}, \tag{1.8}
\]

we estimate the two terms separately.  For the first, the ordinary
Muckenhoupt criterion and the one-dimensional log-concave Poincare bound
give

\[
 \sup_{m<x<\beta}\nu([x,\beta))\int_m^x{dt\over\rho(t)}
 \le C_P(\nu)\le48s^2,                                   \tag{1.9}
\]

and the reflected estimate holds on the left.  The last inequality also
follows from the one-dimensional Cheeger bound
\(\psi_\nu\ge(2\sqrt3s)^{-1}\) and Cheeger's inequality.

It remains to prove the dimensionless curvature estimate.  Let

\[
                         Q(x)=\nu([x,\beta)),\qquad m<x<\beta. \tag{1.10}
\]

Integration by parts gives

\[
\begin{aligned}
 \int_m^x{U''(t)\over\rho(t)}dt
 &=Z\int_m^xU''(t)e^{U(t)}dt\\
 &\le Z\bigl[U'(x)e^{U(x)}-U'(m)e^{U(m)}\bigr]\\
 &={U'(x)\over\rho(x)}-{U'(m)\over\rho(m)},              \tag{1.11}
\end{aligned}
\]

because
\((U'e^U)'=(U''+(U')^2)e^U\ge U''e^U\).  If \(U'(x)>0\), convexity gives

\[
 Q(x)\le\int_x^\beta\rho(x)e^{-U'(x)(t-x)}dt
       \le\int_x^\infty\rho(x)e^{-U'(x)(t-x)}dt
       ={\rho(x)\over U'(x)},                             \tag{1.12}
\]

so \(Q(x)U'(x)/\rho(x)\le1\).  If \(U'(x)\le0\), this term is
nonpositive and may be discarded.

If \(U'(m)\ge0\), the second term on the right of (1.11) is also
nonpositive.  If \(U'(m)<0\), the tangent-line bound on the left gives

\[
 {1\over2}=\nu((\alpha,m])
 \le\int_{-\infty}^m\rho(m)e^{-|U'(m)|(m-t)}dt
 ={\rho(m)\over|U'(m)|},                                  \tag{1.13}
\]

and hence \(|U'(m)|/\rho(m)\le2\).  Since \(Q(x)\le1/2\), its
contribution is at most one.  Thus

\[
 \sup_{x>m}Q(x)\int_m^x{U''(t)\over\rho(t)}dt\le2.        \tag{1.14}
\]

Here is the reflected calculation explicitly.  Put
\(F(x)=\nu((\alpha,x])\) for \(\alpha<x<m\).  From the same integration by
parts,

\[
 \int_x^m{U''(t)\over\rho(t)}dt
 \le {U'(m)\over\rho(m)}-{U'(x)\over\rho(x)}.             \tag{1.14a}
\]

If \(U'(x)<0\), the left tangent bound gives
\(F(x)|U'(x)|/\rho(x)\le1\); if \(U'(x)\ge0\), the second
term in (1.14a) is nonpositive.  If \(U'(m)>0\), the right tangent
bound and \(\nu([m,\beta))=1/2\) give
\(U'(m)/\rho(m)\le2\), whose contribution is at most one because
\(F(x)\le1/2\); if \(U'(m)\le0\), the first term is nonpositive.  Hence

\[
 \sup_{\alpha<x<m}F(x)\int_x^m{U''(t)\over\rho(t)}dt\le2. \tag{1.14b}
\]

Equations
(1.8)--(1.14) yield

\[
                         B_\pm\le48s^2+{2\over\kappa}.    \tag{1.15}
\]

The weighted Muckenhoupt criterion and (1.5) prove (1.3). \(\square\)

Nothing in the proof assumes infinite support: the tangent estimates were
deliberately integrated first only to \(\alpha\) or \(\beta\), and then
bounded by the corresponding infinite exponential integral.  Thus finite
endpoints are included, with one-sided limits if \(x\) approaches an
endpoint.  A nondegenerate one-dimensional log-concave probability is
absolutely continuous on its support interval and therefore has no atoms.
The only atomic case is a point mass; then its variance and every
Poincare variance vanish, so it is handled separately.  For a nonsmooth
convex \(U\), \(U''\) can have atoms and (1.4) should not be read
pointwise.  We do not make that illegitimate identification: Section 5
first convolves the full measure with a Gaussian, applies the present
smooth lemma, and passes to the limit with constants independent of the
smoothing scale.

## 2. Conditional velocity and reinforced Prekopa

We now assume that \(V\in C^\infty(\mathbb R^n)\), that all derivatives
used below are integrable, and that a vanishing global quadratic
confinement has been added if necessary.  Section 5 removes these auxiliary
hypotheses quantitatively.

Use coordinates

\[
                         x=tu+y,\qquad y\in u^\perp.       \tag{2.1}
\]

Write the Hessian in blocks:

\[
 D^2V=
 \begin{pmatrix}
  a&b^T\\ b&C
 \end{pmatrix},qquad
 a=V_{tt},\quad b=\nabla_yV_t,\quad C=D_y^2V.             \tag{2.2}
\]

Hypothesis (0.1) says

\[
 \begin{pmatrix}a&b^T\\b&C-\kappa I\end{pmatrix}\succeq0,\qquad
                         C\succeq\kappa I.                \tag{2.3}
\]

Let

\[
 e^{-U(t)}=\int_{u^\perp}e^{-V(t,y)}dy                    \tag{2.4}
\]

after absorbing the global normalizing constant into \(U\), and let

\[
 d\nu_t(y)=e^{U(t)-V(t,y)}dy.                              \tag{2.5}
\]

Put

\[
                         s_t(y)=V_t(t,y)-\mathbb E_tV_t.  \tag{2.6}
\]

The conditional generator

\[
                         L_t=\Delta_y-\nabla_yV(t,\cdot)\cdot\nabla_y \tag{2.7}
\]

has spectral gap at least \(\kappa\).  Let \(\phi_t\) be the mean-zero
weak solution of

\[
                         L_t\phi_t=s_t,                   \tag{2.8}
\]

and set

\[
                         w_t=\nabla_y\phi_t,qquad
 W(t)=\mathbb E_t|w_t|^2.                                \tag{2.9}
\]

Existence in \(H^1(\nu_t)\) follows from Lax--Milgram and conditional
Poincare; smooth confinement gives the regularity required for the
calculation.  Since

\[
 \partial_t\nu_t=-s_t\nu_t,qquad
 \operatorname {div}_y(\nu_tw_t)=\nu_tL_t\phi_t=s_t\nu_t,\tag{2.10}
\]

the family satisfies the continuity equation

\[
                         \partial_t\nu_t+\operatorname {div}_y(\nu_tw_t)=0.
                                                                    \tag{2.11}
\]

### Lemma 2.1 (reinforced Prekopa inequality)

For every \(t\),

\[
                         \boxed{\kappa W(t)\le U''(t).}   \tag{2.12}
\]

#### Proof

The marginal differentiation formulas are

\[
 U'(t)=\mathbb E_tV_t,qquad
 U''(t)=\mathbb E_ta-\operatorname {Var}_t(V_t).          \tag{2.13}
\]

Bochner's identity for (2.8) gives

\[
 \mathbb E_ts_t^2
 =\mathbb E_t(L_t\phi_t)^2
 =\mathbb E_t\|D_y^2\phi_t\|_{HS}^2
   +\mathbb E_t\langle Cw_t,w_t\rangle.                  \tag{2.14}
\]

Also, integration by parts and \(\nabla_ys_t=b\) give

\[
 \mathbb E_t\langle b,w_t\rangle
 =\mathbb E_t\langle\nabla_ys_t,\nabla_y\phi_t\rangle
 =-\mathbb E_ts_tL_t\phi_t
 =-\mathbb E_ts_t^2.                                     \tag{2.15}
\]

Consequently

\[
 U''(t)=\mathbb E_t\bigl[a+\langle b,w_t\rangle\bigr].   \tag{2.16}
\]

Apply (2.3) pointwise to the vector \((1,w_t(y))\):

\[
 a+2\langle b,w_t\rangle+\langle Cw_t,w_t\rangle
 \ge\kappa|w_t|^2.                                      \tag{2.17}
\]

By (2.14)--(2.15),

\[
 \mathbb E_t\langle Cw_t,w_t\rangle
 \le\mathbb E_ts_t^2=-\mathbb E_t\langle b,w_t\rangle.  \tag{2.18}
\]

Thus the expectation of the left side of (2.17) is no larger than the
right side of (2.16).  Taking expectations in (2.17) proves (2.12).
\(\square\)

This proof uses the cross block \(b\).  Replacing (0.1) merely by
\(C\succeq\kappa I\) would lose (2.17) and is insufficient.

## 3. The derivative of conditional expectation

Let \(f\in C_c^\infty(\mathbb R^n)\), and put

\[
                         g(t)=\mathbb E_t f(t,\cdot).      \tag{3.1}
\]

Differentiating and using (2.10) gives

\[
\begin{aligned}
 g'(t)
 &=\mathbb E_t\partial_tf-\mathbb E_t(fs_t)\\
 &=\mathbb E_t\left(\partial_tf+\nabla_yf\cdot w_t\right).\tag{3.2}
\end{aligned}
\]

The second equality follows from

\[
 \mathbb E_t(fs_t)=\int f\operatorname {div}(\nu_tw_t)dy
                  =-\mathbb E_t(\nabla_yf\cdot w_t).      \tag{3.3}
\]

Cauchy--Schwarz on the product vectors
\((\partial_tf,\nabla_yf)\) and \((1,w_t)\) yields the exact estimate

\[
 \boxed{
 {g'(t)^2\over1+W(t)}
 \le\mathbb E_t\left((\partial_tf)^2+|\nabla_yf|^2\right)
 =\mathbb E_t|\nabla f|^2.}                               \tag{3.4}
\]

This is the term that a naive conditional/marginal proof misses.  Bounding
\(-\operatorname {Cov}_t(f,V_t)\) separately introduces uncontrolled cross
derivatives.  The continuity-equation representation packages the entire
term into the metric speed \(W(t)\), which is then paid for by marginal
curvature through Lemma 2.1.

## 4. Completion of the smooth proof

Let \(\bar\mu\) be the marginal law in (2.4), and put

\[
                         \sigma_u^2=\operatorname {Var}_\mu(t).     \tag{4.1}
\]

The law \(\bar\mu\) is log-concave by Pr\'ekopa.  The total-variance
identity gives

\[
 \operatorname {Var}_\mu f
 =\int\operatorname {Var}_{\nu_t}(f(t,\cdot))d\bar\mu(t)
   +\operatorname {Var}_{\bar\mu}g.                       \tag{4.2}
\]

By (2.3) and Lichnerowicz on each conditional law,

\[
 \int\operatorname {Var}_{\nu_t}(f(t,\cdot))d\bar\mu(t)
 \le{1\over\kappa}\int|\nabla_yf|^2d\mu.                \tag{4.3}
\]

Lemmas 1.1 and 2.1, followed by (3.4), give

\[
\begin{aligned}
 \operatorname {Var}_{\bar\mu}g
 &\le200(\sigma_u^2+\kappa^{-1})
       \int{g'(t)^2\over1+W(t)}d\bar\mu(t)\\
 &\le200(\sigma_u^2+\kappa^{-1})
       \int|\nabla f|^2d\mu.                             \tag{4.4}
\end{aligned}
\]

Combining (4.2)--(4.4),

\[
 \operatorname {Var}_\mu f
 \le201(\sigma_u^2+\kappa^{-1})\int|\nabla f|^2d\mu.     \tag{4.5}
\]

Thus one may take \(C=201\) in the smooth, confined setting.  No constant
depends on \(n\), on \(u\), or on any regularization parameter.

## 5. Nonsmooth convex potentials

Assume now only that \(V:\mathbb R^n\to(-\infty,+\infty]\) is proper,
lower semicontinuous and convex, that \(e^{-V}\) is integrable, and that
(0.1) holds distributionally.  Equivalently,

\[
                         V(x)-{\kappa\over2}|Px|^2        \tag{5.1}
\]

is convex.

Let \(\gamma_\varepsilon=N(0,\varepsilon I)\), and set

\[
                         \mu_\varepsilon=\mu*\gamma_\varepsilon,qquad
 d\mu_\varepsilon=e^{-V_\varepsilon}dx.                  \tag{5.2}
\]

The density is positive and smooth.  More importantly, the curvature loss
is explicit:

\[
                         D^2V_\varepsilon
 \succeq {\kappa\over1+\kappa\varepsilon}P.               \tag{5.3}
\]

To verify (5.3), condition \(X\sim\mu\) on
\(X+G=z\), where \(G\sim N(0,\varepsilon I)\).  The posterior potential is

\[
                         V(x)+{|z-x|^2\over2\varepsilon}, \tag{5.4}
\]

whose Hessian is bounded below by
\(\kappa P+\varepsilon^{-1}I\).  Brascamp--Lieb for linear functions,
valid for distributionally convex potentials by monotone approximation,
gives

\[
 \operatorname {Cov}(X\mid X+G=z)
 \preceq(\kappa P+\varepsilon^{-1}I)^{-1}.                \tag{5.5}
\]

The standard Gaussian-channel differentiation identity is

\[
 D^2V_\varepsilon(z)
 ={1\over\varepsilon}I
  -{1\over\varepsilon^2}
       \operatorname {Cov}(X\mid X+G=z).                  \tag{5.6}
\]

Equations (5.5)--(5.6) give zero in the \(u\)-direction and
\(\kappa/(1+\kappa\varepsilon)\) on \(u^\perp\), proving (5.3).

Also,

\[
 \operatorname {Var}_{\mu_\varepsilon}\langle X,u\rangle
 =\operatorname {Var}_\mu\langle X,u\rangle+\varepsilon. \tag{5.7}
\]

To match exactly the temporary confinement used in Section 2, first apply
the smooth proof to the probability with potential
\(V_\varepsilon+\delta|x|^2/2\).  Its curvature is still at least
\(\kappa_\varepsilon P\), where
\(\kappa_\varepsilon=\kappa/(1+\kappa\varepsilon)\).  As
\(\delta\downarrow0\), its density, weak-direction second moment, and all
compactly supported test-function integrals converge to those of
\(\mu_\varepsilon\) by dominated convergence; log-concavity guarantees
the required finite second moment.  Thus (4.5) applies to
\(\mu_\varepsilon\) itself:

\[
\begin{aligned}
 \operatorname {Var}_{\mu_\varepsilon}f
 &\le201\left(\kappa_\varepsilon^{-1}
       +\operatorname {Var}_{\mu_\varepsilon}\langle X,u\rangle\right)
       \int|\nabla f|^2d\mu_\varepsilon\\
 &\le201\left(\kappa^{-1}
       +\operatorname {Var}_\mu\langle X,u\rangle+2\varepsilon\right)
       \int|\nabla f|^2d\mu_\varepsilon.                 \tag{5.8}
\end{aligned}
\]

For \(f\in C_c^\infty\), convolution convergence gives convergence of the
variance and the Dirichlet integral in (5.8).  Letting
\(\varepsilon\downarrow0\) proves (4.5) for such \(f\) under the original
measure.

Here are the test-function limits explicitly.  For a locally Lipschitz
\(f\in L^2(\mu)\) with finite Dirichlet integral, let
\(T_M(r)=\max(-M,\min(r,M))\).  Then
\(T_M(f)\to f\) in \(L^2(\mu)\) and

\[
 |\nabla T_M(f)|=1_{\{|f|<M\}}|\nabla f|\quad\mu\text{-a.e.}             \tag{5.9}
\]

Choose a Lipschitz cutoff \(\chi_R\), equal to one on \(B_R\), zero
outside \(B_{2R}\), and with \(|\nabla\chi_R|\le2/R\).  For fixed \(M\),
\(\chi_RT_M(f)\to T_M(f)\) in \(L^2(\mu)\), while expansion of its gradient
and Cauchy--Schwarz show

\[
 \int|\nabla(\chi_RT_M(f))|^2d\mu
 \longrightarrow\int|\nabla T_M(f)|^2d\mu.              \tag{5.10}
\]

Indeed, the cutoff-gradient term has norm at most \(2M/R\), and the main
term converges by dominated convergence.  Each compactly supported
Lipschitz function \(h\) can next be convolved with a standard mollifier:
\(h_\delta\to h\) uniformly, \(|\nabla h_\delta|\le\operatorname{Lip}(h)\),
and \(\nabla h_\delta\to\nabla h\) Lebesgue-a.e.  Since a full-dimensional
log-concave measure is absolutely continuous, dominated convergence gives
convergence of both variance and Dirichlet integral.  Applying (4.5) in
the order \(\delta\downarrow0\), \(R\uparrow\infty\), and \(M\uparrow\infty\)
proves it for \(f\).  If the original Dirichlet integral is infinite the
inequality is automatic.  If the support is an affine subspace, this
entire argument, including Gaussian convolution and mollification, is
carried out intrinsically in that subspace, with \(P\) restricted to it.

This passage preserves the constant \(201\); there is no
regularization-dependent loss.

## 6. Adversarial models

### 6.1 Sheared products

Consider, first in one transverse dimension, and take
\(W(t)=ct^2/2\) with \(c>0\):

\[
                         V(t,y)=W(t)+{L\over2}(y-mt)^2,qquad L>\kappa.\tag{6.1}
\]

The Hessian blocks are

\[
 a=W''+Lm^2,qquad b=-Lm,qquad C=L.                       \tag{6.2}
\]

The condition \(D^2V\succeq\kappa P_{u^\perp}\) is equivalent to

\[
                         W''\ge {\kappa L\over L-\kappa}m^2.         \tag{6.3}
\]

The marginal potential is \(U(t)=ct^2/2+\mathrm{const}\), and
\(\nu_t=N(mt,L^{-1})\).  This model checks every sign in Sections 2--3.
Indeed,

\[
 s_t(y)=-Lm(y-mt),\qquad
 \phi_t(y)=m(y-mt),\qquad w_t=m,                         \tag{6.4}
\]

and hence

\[
 \partial_t\nu_t=Lm(y-mt)\nu_t=-s_t\nu_t,
 \qquad
 \partial_y(\nu_tw_t)=-Lm(y-mt)\nu_t=s_t\nu_t.          \tag{6.5}
\]

Thus the plus sign in the continuity equation and the plus sign in
(3.2) are forced.  Moreover,

\[
 \mathbb E_ts_t^2=Lm^2,qquad
 \mathbb E_t(bw_t)=-Lm^2,qquad
 \mathbb E_t(Cw_t^2)=Lm^2,                              \tag{6.6}
\]

so (2.14), (2.15), (2.16), and (2.18) are identities, with
\(U''=c\).  Evaluating the block quadratic form on \((1,w_t)=(1,m)\)
also gives exactly

\[
 a+2bw_t+Cw_t^2=c.                                      \tag{6.7}
\]

Consequently (0.3) reads \(\kappa m^2\le c\), a slightly weaker
consequence of (6.3).  Finally, for the affine test function
\(f(t,y)=\alpha t+\beta y\), formula (3.2) gives
\(g'=\alpha+\beta m\), and (3.4) becomes

\[
 { (\alpha+\beta m)^2\over1+m^2}\le\alpha^2+\beta^2,    \tag{6.8}
\]

with equality precisely when \((\alpha,\beta)\) is parallel to
\((1,m)\).  Thus the velocity normalization in (3.4) is sharp.  A large
shear is paid for by curvature of the weak marginal.  In the limiting
case \(L=\kappa\), the block inequality forces \(m=0\).

### 6.2 Sharply varying conditional means

Let

\[
                         m(t)=\int y\,d\nu_t(y).           \tag{6.9}
\]

The continuity equation gives

\[
                         m'(t)=\int w_t(y)d\nu_t(y).       \tag{6.10}
\]

Consequently Lemma 2.1 implies

\[
                         \kappa|m'(t)|^2
 \le\kappa\int|w_t|^2d\nu_t\le U''(t).                  \tag{6.11}
\]

Thus a rapidly moving conditional center cannot be inserted for free.  It
forces a proportionally large curvature spike in the marginal potential.
Lemma 1.1 shows exactly why such a spike costs at most \(C/\kappa\) in the
global Poincare constant.  This rules out both affine shears and sharply
varying conditional means as counterexamples.

## 7. Use after a rank-one localization stop

Suppose a localization or regularization argument produces a posterior
measure for which

\[
                         D^2V\succeq\kappa P_{u^\perp}    \tag{7.1}
\]

with one uncontrolled direction \(u\), while its covariance in that
direction is at most \(R^2\).  Then (0.2) gives immediately

\[
                         C_P(\mu)\le96(\kappa^{-1}+R^2).   \tag{7.2}
\]

Hence a localization stopped with one weak direction is sufficient if it
supplies a universal lower bound on the transverse curvature and a
universal upper bound on the survivor variance.  No control of the
derivative of the conditional mean is additionally required; it is already
contained in (0.3).

The constant \(96\), improving the elementary constant \(201\) obtained
in Sections 1--5, follows from the more general defect-subspace theorem
below.

## 8. A weak subspace of arbitrary dimension

The rank-one argument has a natural and fully provable extension.  It is
important that its marginal input is the *actual* Poincare constant of the
weak marginal; no dimension-free estimate for that constant is being
assumed.

### Theorem 8.1 (strong conditionals plus weak marginal)

Let \(\mathbb R^n=F\oplus E\) be an orthogonal decomposition, let
\(\pi_F\) be the orthogonal projection onto \(F\), and let
\(\bar\mu=(\pi_F)_\#\mu\).  Suppose that \(\mu=e^{-V}dx\) is log-concave
and, distributionally,

\[
                         D^2V\succeq\kappa P_E           \tag{8.1}
\]

for some \(\kappa>0\).  Then

\[
 \boxed{
 C_P(\mu)\le \kappa^{-1}
       +2\bigl(C_P(\bar\mu)+\kappa^{-1}\bigr)
 \le3\bigl(C_P(\bar\mu)+\kappa^{-1}\bigr).}            \tag{8.2}
\]

The assertion is intrinsic if the support is a proper affine subspace.  A
point-mass marginal is assigned \(C_P=0\), as all of its variances vanish.

#### Step 1: a matrix reinforced Prekopa inequality

First assume \(V\) is smooth and confined as in Section 2.  Write
\(x=z+y\), with \(z\in F\), \(y\in E\), and write

\[
 D^2V=\begin{pmatrix}A&B^*\\B&C\end{pmatrix}.            \tag{8.3}
\]

Thus (8.1) says

\[
 \begin{pmatrix}A&B^*\\B&C-\kappa I_E\end{pmatrix}\succeq0. \tag{8.4}
\]

Let \(e^{-U(z)}dz\) be the \(F\)-marginal and let \(\nu_z\) be the
conditional probability on \(E\).  For \(\xi\in F\), define the centered
conditional score

\[
 s_\xi(y)=\langle\nabla_zV(z,y),\xi\rangle
       -\mathbb E_z\langle\nabla_zV(z,\cdot),\xi\rangle, \tag{8.5}
\]

solve \(L_z\phi_\xi=s_\xi\) on \(E\), and set
\(\mathcal W_z(y)\xi=\nabla_y\phi_\xi(y)\).  The solution and hence
\(\mathcal W_z\) depend linearly on \(\xi\).  Define the positive
semidefinite operator \(G(z):F\to F\) by

\[
 \langle G(z)\xi,\xi\rangle
       =\mathbb E_z|\mathcal W_z\xi|^2.                  \tag{8.6}
\]

Applying exactly the scalar calculation (2.13)--(2.18) to the direction
\(\xi\), with \(a=\langle A\xi,\xi\rangle\),
\(b=B\xi\), and \(w=\mathcal W_z\xi\), gives

\[
 \langle D^2U(z)\xi,\xi\rangle
       \ge\kappa\mathbb E_z|\mathcal W_z\xi|^2.          \tag{8.7}
\]

Since this holds for every \(\xi\),

\[
                         \boxed{D^2U(z)\succeq\kappa G(z).} \tag{8.8}
\]

The conditional continuity equations in every \(F\)-direction also give,
for \(g(z)=\mathbb E_zf(z,\cdot)\),

\[
 \langle\nabla g(z),\xi\rangle
 =\mathbb E_z\left[
      \langle\nabla_zf,\xi\rangle
      +\langle\nabla_yf,\mathcal W_z\xi\rangle\right].  \tag{8.9}
\]

Put \(q=(\nabla_zf,\nabla_yf)\), and let
\(\mathcal B_z(y):F\to F\oplus E\) be
\(\mathcal B_z(y)\xi=(\xi,\mathcal W_z(y)\xi)\).  Then
\(\nabla g=\mathbb E_z\mathcal B_z^*q\) and
\(\mathbb E_z\mathcal B_z^*\mathcal B_z=I_F+G(z)\).  Matrix
Cauchy--Schwarz therefore yields

\[
 \boxed{
 \langle(I_F+G)^{-1}\nabla g,\nabla g\rangle
 \le\mathbb E_z|\nabla f|^2.}                           \tag{8.10}
\]

Indeed, with \(\alpha=(I_F+G)^{-1}\nabla g\),
\(\langle\alpha,\nabla g\rangle
=\mathbb E_z\langle\mathcal B_z\alpha,q\rangle\); ordinary
Cauchy--Schwarz and cancellation of
\(\langle\alpha,(I_F+G)\alpha\rangle\) prove (8.10).

#### Step 2: the needed matrix-weighted marginal inequality

We record and prove the exact improved Lichnerowicz inequality used here.

### Lemma 8.2 (improved Brascamp--Lieb from the actual gap)

Let \(d\eta=e^{-U}dz\) be a smooth log-concave probability on
\(\mathbb R^k\), let \(H=D^2U\succeq0\), assume
\(0<C_P(\eta)<\infty\), and put \(\lambda=C_P(\eta)^{-1}\).  Then every
smooth \(h\in L^2(\eta)\) satisfies

\[
 \operatorname {Var}_\eta h
 \le2\int\langle(H+\lambda I)^{-1}\nabla h,\nabla h\rangle d\eta. \tag{8.11}
\]

#### Proof

Let \(L=\Delta-\nabla U\cdot\nabla\), and let \(A_0=-L\) be the
nonnegative self-adjoint operator associated with the closed scalar
Dirichlet form, restricted to mean-zero functions.  Let \(\mathcal G\) be
the \(L^2(\eta;\mathbb R^k)\)-closure of gradients, and on it use the
closed one-form quadratic form associated with

\[
                         A_1=-L\otimes I+H               \tag{8.12}
\]

The gap makes \(A_0^{-1}\) bounded.  On a smooth core, and hence after
closure, the intertwining identity

\[
                         \nabla A_0=A_1\nabla             \tag{8.13}
\]

holds.  If \(r=\nabla a\), the scalar spectral gap gives

\[
\begin{aligned}
 \langle r,A_1r\rangle_{L^2(\eta)}
 &=\|A_0a\|_2^2
 \ge\lambda\langle a,A_0a\rangle
 =\lambda\|r\|_2^2.                                    \tag{8.14}
\end{aligned}
\]

On all vector fields, integration by parts also gives

\[
 \langle r,A_1r\rangle
 =\int\bigl(\|\nabla r\|_{HS}^2+\langle Hr,r\rangle\bigr)d\eta
 \ge\int\langle Hr,r\rangle d\eta.                     \tag{8.15}
\]

Thus, as quadratic forms on \(\mathcal G\),

\[
 A_1\succeq\tfrac12 P_{\mathcal G}(H+\lambda I)P_{\mathcal G}. \tag{8.16}
\]

The estimate (8.14) also shows that the operator associated with the
restriction of the \(A_1\)-form to \(\mathcal G\) has a bounded inverse.
For mean-zero \(h\), put \(a=A_0^{-1}h\).  Then
\(A_1\nabla a=\nabla h\) in the weak sense, and integration by parts gives
the Helffer--Sjostrand identity

\[
 \operatorname {Var}_\eta h
 =\langle\nabla h,\nabla a\rangle
 =\langle\nabla h,A_1^{-1}\nabla h\rangle_{\mathcal G}.   \tag{8.17}
\]

There is no hidden regularity assertion in the inverse comparison.  Its
variational formula, (8.16), and enlargement of the supremum from
\(\mathcal G\) to all vector fields give

\[
\begin{aligned}
 \langle\nabla h,A_1^{-1}\nabla h\rangle
 &=\sup_{r\in\mathcal G}
   \{2\langle\nabla h,r\rangle-\langle r,A_1r\rangle\}\\
 &\le\sup_{r\in\mathcal G}
   \{2\langle\nabla h,r\rangle-\tfrac12
     \langle r,(H+\lambda I)r\rangle\}\\
 &\le2\int\langle(H+\lambda I)^{-1}\nabla h,\nabla h\rangle d\eta.
                                                                    \tag{8.18}
\end{aligned}
\]

Approximation by the smooth core justifies the identities for every
function in the form domain.  This proves (8.11). \(\square\)

Return to the marginal \(\bar\mu\), and write
\(H=D^2U\), \(C_0=C_P(\bar\mu)\), and \(\lambda=C_0^{-1}\).  From
(8.8), inverse order gives

\[
 (I_F+G)^{-1}\succeq\kappa(\kappa I_F+H)^{-1}.           \tag{8.19}
\]

For every scalar \(r\ge0\),

\[
 {1\over r+\lambda}
 \le\left({1\over\lambda}+{1\over\kappa}\right)
       {\kappa\over\kappa+r};                            \tag{8.20}
\]

indeed, after clearing denominators the difference is
\(\kappa r+\lambda^2\ge0\).  Functional calculus for \(H\), followed by
(8.19), therefore gives

\[
 (H+\lambda I)^{-1}
 \preceq(C_0+\kappa^{-1})(I_F+G)^{-1}.                  \tag{8.21}
\]

Lemma 8.2 and (8.10) now imply

\[
 \operatorname {Var}_{\bar\mu}g
 \le2(C_0+\kappa^{-1})\int|\nabla f|^2d\mu.             \tag{8.22}
\]

Every conditional law on \(E\) is \(\kappa\)-strongly log-concave, so its
variance contribution is at most
\(\kappa^{-1}\int|\nabla_yf|^2d\mu\).  The total-variance identity proves
(8.2) in the smooth case.

#### Step 3: nonsmooth potentials and constants

For a distributionally convex potential, convolve \(\mu\) with
\(N(0,\varepsilon I)\).  Section 5 shows that (8.1) becomes

\[
 D^2V_\varepsilon\succeq
       {\kappa\over1+\kappa\varepsilon}P_E.              \tag{8.23}
\]

The weak marginal becomes
\(\bar\mu_\varepsilon=\bar\mu*N(0,\varepsilon I_F)\), and

\[
 C_P(\bar\mu_\varepsilon)\le C_P(\bar\mu)+\varepsilon. \tag{8.24}
\]

To see (8.24), apply Poincare first in the \(\bar\mu\) coordinate and then
in the independent Gaussian coordinate to \(h(Z+G)\); Jensen bounds both
gradient terms by the Dirichlet integral under \(Z+G\).  Apply the smooth
form of (8.2) and let \(\varepsilon\downarrow0\), using the explicit
test-function approximation of Section 5.  Since
\(\kappa_\varepsilon^{-1}=\kappa^{-1}+\varepsilon\), no constant is lost.
This proves (8.2) in the stated distributional generality.

For \(\dim F=1\), the standard one-dimensional log-concave estimate
\(C_P(\bar\mu)\le48\operatorname {Var}_{\bar\mu}(t)\) in (8.2) yields

\[
 C_P(\mu)\le3\kappa^{-1}
       +96\operatorname {Var}_\mu\langle X,u\rangle
 \le96\left(\kappa^{-1}
       +\operatorname {Var}_\mu\langle X,u\rangle\right). \tag{8.25}
\]

This proves (0.2) with \(C=96\).

Ordinary pointwise Brascamp--Lieb alone does **not** yield (8.21): its
weight \(H^{-1}\) blows up on flat marginal directions.  What makes the
matrix-weighted inequality true is the global marginal gap, inserted by
the exact-one-form estimate (8.14).  Thus the \(k\)-dimensional extension
does follow, but it follows from the proved improved-Lichnerowicz argument
in Lemma 8.2, not from a formal interpolation between ordinary Poincare
and ordinary Brascamp--Lieb.
