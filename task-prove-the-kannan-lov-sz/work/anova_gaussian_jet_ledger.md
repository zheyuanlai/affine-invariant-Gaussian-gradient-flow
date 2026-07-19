# The exact Gaussian-jet ledger for the affine-orthogonal ANOVA gate

## 0. Purpose and verdict

Let \(X\) be centered and isotropic in \(\mathbb R^n\), let
\(G\sim N(0,I_n)\) be independent, and put

\[
 S=\frac{X+G}{\sqrt2},\qquad \nu=\mathcal L(S).
\]

For centered \(f\in L^2(\nu)\), define

\[
 F(x,g)=f((x+g)/\sqrt2),\quad
 h(x)=E_G F(x,G),\quad v(g)=E_XF(X,g),
\]

\[
 R(X,G)=F(X,G)-h(X)-v(G).
\]

This note gives an exact diagonalization of the conditional-Gaussian
part of this decomposition.  It proves that the proposed ANOVA estimate

\[
 \operatorname {dist}_{L^2(\nu)}(f,\mathrm {Aff})^2
 \le C_A E R^2                                             \tag{AI}
\]

is equivalent to a dimension-free coercivity inequality for the complete
analytic jet of the heat transform \(h\).  No such coercivity estimate is
proved here.  Its quadratic restriction is exactly the generalized
quadratic-variance obstruction

\[
 \operatorname {Var}(X^TAX)\lesssim \operatorname {tr}(A^2),
\]

including every symmetric matrix \(A\), not only \(A=I\).  Thus Gaussian
chaos expansion does not make (AI) formal: it locates the conjecture-strength
step exactly in the base-law variation of the heat-transform jets.

## 1. Conditional Hermite coefficients

For a multi-index \(\alpha\in\mathbb N^n\), let
\(H_\alpha=\prod_i H_{\alpha_i}\) be the probabilists' Hermite polynomial
and put

\[
 \varphi_\alpha(g)=\frac{H_\alpha(g)}{\sqrt{\alpha!}}.
\]

The family \((\varphi_\alpha)_\alpha\) is an orthonormal basis of
\(L^2(\gamma_n)\).  First suppose \(f\) is bounded and smooth.  Gaussian
integration by parts gives, for every \(x\) and \(\alpha\),

\[
 \begin{aligned}
 c_\alpha(x)
 &:=E_G[F(x,G)\varphi_\alpha(G)]\\
 &=\frac1{\sqrt{\alpha!}}E_G[\partial_g^\alpha F(x,G)]
 =\frac{\partial^\alpha h(x)}{\sqrt{\alpha!}}.       \tag{1.1}
 \end{aligned}
\]

There is no missing power of two in the last equality: differentiation of
\(F\) with respect to \(x\) or to \(g\) gives the same factor
\(2^{-1/2}\) at each derivative.

For a general \(f\in L^2(\nu)\), the conditional Hermite coefficients are
defined for \(\mu\)-almost every \(x\), and (1.1) holds in the weak
heat-transform sense.  Indeed bounded smooth functions are dense in
\(L^2(\nu)\), while the maps

\[
 f\longmapsto c_\alpha(X)
 =E[F(X,G)\varphi_\alpha(G)\mid X]
\]

are contractions from \(L^2(\nu)\) to \(L^2(\mu)\).  All identities below
therefore pass by \(L^2\) closure; no pointwise differentiability of the
original \(f\) is required.

Conditional Parseval gives

\[
 F(x,g)=\sum_{\alpha}c_\alpha(x)\varphi_\alpha(g)
 \quad\hbox{in }L^2(\mu\otimes\gamma).                \tag{1.2}
\]

Since \(c_0=h\), and since the Hermite coefficient of \(v\) is

\[
 m_\alpha:=E_Xc_\alpha(X)
 =\frac{E_\mu\partial^\alpha h}{\sqrt{\alpha!}},      \tag{1.3}
\]

centering of \(f\) gives \(m_0=Eh=0\) and

\[
 v(g)=\sum_{|\alpha|\ge1}m_\alpha\varphi_\alpha(g). \tag{1.4}
\]

Subtracting (1.4) and the zeroth coefficient from (1.2) yields

\[
 R(x,g)=\sum_{|\alpha|\ge1}
       (c_\alpha(x)-m_\alpha)\varphi_\alpha(g).       \tag{1.5}
\]

Consequently

\[
 \boxed{
 ER^2=\sum_{|\alpha|\ge1}
 \frac{\operatorname {Var}_\mu(\partial^\alpha h)}{\alpha!}.}
                                                               \tag{1.6}
\]

This is an equality, not a Poincare or smoothness estimate.

## 2. The full norm and the affine projection

Write

\[
 V_\alpha=\frac{\operatorname {Var}_\mu(\partial^\alpha h)}{\alpha!},
 \qquad
 M_\alpha=\frac{|E_\mu\partial^\alpha h|^2}{\alpha!}.
\]

Parseval and (1.2) give

\[
 \boxed{
 \|f\|_{L^2(\nu)}^2
 =E_\mu h^2+
   \sum_{|\alpha|\ge1}(V_\alpha+M_\alpha).}          \tag{2.1}
\]

Let

\[
 p=E_\mu[Xh(X)],\qquad m_1=E_\mu\nabla h(X).
\]

Because \(E[XF]=E[Xh(X)]\), while the degree-one Gaussian Hermite
coefficient is \(E[GF]=E\nabla h\),

\[
 \ell:=E_\nu[Sf(S)]=\frac{p+m_1}{\sqrt2}.             \tag{2.2}
\]

The law \(\nu\) is isotropic, so orthogonal projection onto the affine
functions and (2.1)--(2.2) give the second exact identity

\[
 \boxed{
 \begin{aligned}
 \operatorname {dist}_{L^2(\nu)}(f,\mathrm {Aff})^2
 &=E h^2+\sum_{|\alpha|\ge1}(V_\alpha+M_\alpha)\\
 &\quad-\frac12|E[Xh]+E\nabla h|^2.
 \end{aligned}}                                      \tag{2.3}
\]

In particular, \(f\perp\mathrm {Aff}\) if and only if

\[
 Eh=0,qquad E[Xh]+E\nabla h=0.                       \tag{2.4}
\]

On that subspace, (AI) is exactly

\[
 \boxed{
 E h^2+\sum_{|\alpha|\ge1}M_\alpha
 \le (C_A-1)\sum_{|\alpha|\ge1}V_\alpha,}
                                                               \tag{AJ}
\]

for every \(h\) in the range of the Gaussian heat-transform operator
\(U_X:f\mapsto E[f((X+G)/\sqrt2)\mid X]\) satisfying (2.4).  Conversely,
(AJ) and (1.6)--(2.3) imply (AI).  Thus (AJ) is neither a consequence of
conditional Parseval nor a discarded regularity term; it is precisely the
missing coercive assertion.

## 3. Operator and maximal-correlation interpretation

Let \(T f=F\), and let \(P_X,P_G\) denote the centered conditional
expectation projections in \(L^2(\mu\otimes\gamma)\).  Then

\[
 Q_X=T^*P_XT,qquad Q_G=T^*P_GT,
\]

and

\[
 ER^2=\langle f,(I-Q_X-Q_G)f\rangle.                 \tag{3.1}
\]

Equations (1.5)--(1.6) identify (3.1) with the sum of the base-law
variances of every nonconstant conditional Gaussian-chaos coefficient.
A bound on either scalar maximal correlation alone does not establish
(AI): the two families of mean coefficients \((m_\alpha)\) and the base
term \(h\) are exactly the parts omitted by the right side of (1.6).
For Gaussian \(X\), orthogonal Hermite degree diagonalizes both channels
and gives defect \(1-2^{1-k}\) in degree \(k\), so \(C_A=2\).  For a
general log-concave input, (AJ) is the required replacement for that
Gaussian degree orthogonality.

If \(a=\lambda_1(\mu)\), ordinary Poincare applied to every derivative
only yields, with

\[
 V_k=\sum_{|\alpha|=k}V_\alpha,qquad
 M_k=\sum_{|\alpha|=k}M_\alpha,
\]

the one-sided hierarchy

\[
 E h^2\le\frac1a(V_1+M_1),
 \qquad
 V_k\le\frac{k+1}{a}(V_{k+1}+M_{k+1}).               \tag{3.2}
\]

The combinatorial factor follows from
\(\sum_i1/(\beta-e_i)!=(k+1)/\beta!\) for
\(|\beta|=k+1\).  The direction of (3.2) does not control the mean jets
\(M_k\) by the variance jets, and its explicit \(1/a\) is the original
spectral-gap obstruction.  Dropping either feature would assume the
content needed in (AJ).

### 3.1 Bottom-window transfer back to the input

The Dirichlet energy has an exact jet expansion.  If
\(f\in H^1(\nu)\), Gaussian Parseval applied to each weak derivative in
the \(g\)-variable gives

\[
 \boxed{
 \frac12\int |\nabla f|^2\,d\nu
 =\sum_{k\ge1}k(V_k+M_k).}                           \tag{3.3}
\]

Indeed \(\nabla_gF=2^{-1/2}\nabla f((x+g)/\sqrt2)\), while the Gaussian
Dirichlet form of the conditional Hermite series is
\(\sum_\alpha |\alpha|E|c_\alpha(X)|^2\).  The identity follows first
for bounded smooth \(f\) and then for every \(f\in H^1(\nu)\) by closure
of the Gaussian Dirichlet form and the coefficient maps from Section 1.

Now suppose

\[
 f\perp\operatorname {Aff},\qquad \|f\|_2=1,
 \qquad q:=\int|\nabla f|^2d\nu<1,
\]

and write \(m=E_\mu\nabla h\).  Equation (2.4) says
\(E[Xh]=-m\).  Therefore

\[
 h_0(x):=h(x)+m\cdot x                                  \tag{3.4}
\]

is centered and affine-orthogonal in \(L^2(\mu)\).  Isotropy, (2.1),
and (3.3) give the exact identities

\[
 \begin{aligned}
 \|h_0\|_2^2
 &=Eh^2-|m|^2\\
 &=1-\sum_{k\ge1}(V_k+M_k)-M_1,\\
 \int|\nabla h_0|^2d\mu&=V_1+4M_1.                  \tag{3.5}
 \end{aligned}
\]

Since

\[
 \sum_{k\ge1}(V_k+M_k)+M_1
 \le 2\sum_{k\ge1}k(V_k+M_k)=q,
\]

and \(V_1+4M_1\le4(V_1+M_1)\le2q\), it follows that

\[
 \boxed{
 \|h_0\|_2^2\ge1-q,
 \qquad
 \frac{\int|\nabla h_0|^2d\mu}{\|h_0\|_2^2}
 \le\frac{2q}{1-q}.}                                \tag{3.6}
\]

The same energy ledger also gives

\[
 \boxed{
 ER^2=\sum_{k\ge1}V_k\le\frac q2,
 \qquad
 \sum_{k\ge1}kM_k\le\frac q2.}                    \tag{3.7}
\]

Thus every mean-jet term is small on a bottom spectral window, whereas
the potentially order-one term in (AJ) is \(Eh^2\).  Formula (3.6)
identifies that term with an affine-orthogonal low-energy mode for the
input law.  Ordinary input Poincare controls it only with the factor
\(1/a\) in (3.2); eliminating that factor is the original KLS-type
obstruction, not a consequence of Gaussian smoothing.

### 3.2 What the weak output eigen-equation adds

The natural adjoint-channel test confirms, rather than removes, the
input-gap obstruction.  Let
\[
 ({\cal U}f)(x)=E[f(S)\mid X=x],\qquad
 ({\cal U}^*u)(s)=E[u(X)\mid S=s].
\]
The posterior law of \(X\) given \(S=s\) has potential
\[
 V(x)+\frac12|x-\sqrt2s|^2.
\]
It is \(1\)-strongly log-concave on its intrinsic support.  Posterior
differentiation, covariance Cauchy--Schwarz, and the posterior Poincare
inequality therefore give, pointwise in \(s\),
\[
 \begin{aligned}
 \nabla {\cal U}^*u(s)
 &=\sqrt2\operatorname {Cov}(X,u(X)\mid S=s),\\
 |\nabla {\cal U}^*u(s)|^2
 &\le2\operatorname {Var}(u(X)\mid S=s)
 \le2E[|\nabla u(X)|^2\mid S=s].
 \end{aligned}                                      \tag{3.8}
\]
Consequently
\[
 \int|\nabla {\cal U}^*u|^2d\nu
 \le2\int|\nabla u|^2d\mu.                       \tag{3.9}
\]
The statement holds first for smooth \(u\) and then on the closed form
domain; hard convex supports follow from the intrinsic strongly
log-concave posterior form.

Suppose in addition to the hypotheses of Section 3.1 that \(f\) is an
exact weak eigenfunction of the output generator:
\[
 \int\langle\nabla f,\nabla\phi\rangle d\nu
 =q\int f\phi\,d\nu
 \quad\text{for every form-domain }\phi.           \tag{3.10}
\]
Taking \(\phi={\cal U}^*h_0\) is legitimate by (3.9).  Moreover
\[
 \langle f,{\cal U}^*h_0\rangle_{L^2(\nu)}
 =E[h(X)h_0(X)]
 =Eh^2-M_1
 =\|h_0\|_2^2.                                     \tag{3.11}
\]
Cauchy--Schwarz in (3.10), followed by (3.9), yields
\[
 q\|h_0\|_2^2
 \le\sqrt q\left(2\int|\nabla h_0|^2d\mu\right)^{1/2}.
\]
Together with (3.6), this proves the two-sided comparison
\[
 \boxed{
 \frac{q(1-q)}2
 \le
 \frac{\int|\nabla h_0|^2d\mu}{\|h_0\|_2^2}
 \le
 \frac{2q}{1-q}.}                                   \tag{3.12}
\]
For \(q\le1/2\), the input Rayleigh quotient lies between \(q/4\) and
\(4q\).  Thus the weak eigen-equation shows that a nonlinear low output
mode transfers to a comparably low affine-orthogonal input mode.  It
does not supply an absolute lower bound or control \(Eh^2\) by
\(\sum_kV_k\); using an input spectral gap at this point is circular.

## 4. Exact quadratic specialization

Let \(A=A^T\) and

\[
 f_A(s)=s^TAs-\operatorname {tr}A.
\]

Then

\[
 h_A(x)=\frac12(x^TAx-\operatorname {tr}A),qquad
 \nabla h_A=Ax,qquad D^2h_A=A.                      \tag{4.1}
\]

Hence

\[
 \sum_{|\alpha|=1}V_\alpha=\operatorname {tr}(A^2),
 \qquad
 \sum_{|\alpha|=2}M_\alpha=\frac12\operatorname {tr}(A^2),
                                                               \tag{4.2}
\]

and all other nonzero jet terms vanish.  If

\[
 Q=X^TAX,qquad d_A=E[X(Q-\operatorname {tr}A)],
\]

then (1.6) and (2.3) become

\[
 ER^2=\operatorname {tr}(A^2),                       \tag{4.3}
\]

\[
 \operatorname {dist}(f_A,\mathrm {Aff})^2
 =\frac14\operatorname {Var}(Q)
   +\frac32\operatorname {tr}(A^2)-\frac18|d_A|^2.  \tag{4.4}
\]

For centrally symmetric \(X\), \(d_A=0\).  Therefore (AI) would imply

\[
 \operatorname {Var}(X^TAX)
 \le(4C_A-6)\operatorname {tr}(A^2)                  \tag{4.5}
\]

for every symmetric \(A\).  The dimension-free thin-shell theorem controls
the radial choice \(A=I\), but (4.5) is the full matrix-weighted statement;
replacing it by the radial theorem loses the directional information that
(AJ) requires.

## 5. Consequence for the proof search

The ANOVA operator has the correct exact kernel on full-dimensional inputs,
and the conditional Gaussian chaos expansion has no hidden dimension loss.
Nevertheless, a kernel statement plus Parseval does not provide a uniform
spectral angle.  A proof of (AI), even only on the affine-orthogonal
quadratic sector, must establish (4.5); a proof on bottom spectral windows
must establish the corresponding low-mode version of (AJ).  Either is a
genuinely new dimension-free coercivity mechanism.  Any argument which
replaces the right side of (AJ) by all jet energies, invokes an unproved
uniform Poincare inequality for \(\mu\), or bounds only the radial quadratic
is incomplete.

