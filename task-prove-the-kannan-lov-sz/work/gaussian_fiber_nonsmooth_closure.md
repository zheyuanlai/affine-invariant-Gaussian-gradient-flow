# Nonsmooth closure for Gaussian conditional fibers

## 0. Status and conclusion

This note closes the approximation gap for the Gaussian-conditional-fiber
subclass at a fixed positive transverse covariance floor.  It does not
remove the Gaussian-fiber hypothesis and it does not pass to zero transverse
noise.  Those remain separate, conjecture-strength problems.

The proof is dimension-free.  Its main point is that one must mollify the
expanded quadratic potential, rather than mollifying the covariance and
centroid separately.  This preserves joint convexity exactly.  The explicit
Gaussian Poisson formula then turns the target into a local first-derivative
quadratic form, for which strong local Sobolev convergence is available.

## 1. Closure theorem

Let \(J\subseteq\mathbb R\) be an interval with nonempty interior, and let
\(U:\mathbb R\times\mathbb R^d\to(-\infty,+\infty]\) be lower
semicontinuous and convex, finite on
\(\operatorname {int}J\times\mathbb R^d\), and \(+\infty\) off
\(J\times\mathbb R^d\).  Suppose that, for every
\(s\in\operatorname {int}J\),

\[
 U(s,z)=c(s)+\frac12z^TQ(s)z-r(s)^Tz,                 \tag{1.1}
\]

where \(Q(s)\in\mathbb S_{++}^d\).  Assume that
\(Z^{-1}e^{-U(s,z)}\,ds\,dz\) is a probability with finite second moment.
Write

\[
 R(s)=Q(s)^{-1},\qquad m(s)=R(s)r(s),                 \tag{1.2}
\]

and let \(\rho\) be the \(S\)-marginal.  Assume

\[
 \mathbb ES=0,\qquad \mathbb ES^2=1,\qquad
 \mathbb E[S\,m(S)]=0,                                \tag{1.3}
\]

where the vector moment is absolutely defined, and, for two finite
constants \(\kappa,\Lambda>0\),

\[
 R(s)\succeq\kappa I_d\quad(s\in\operatorname {int}J),
 \qquad
 \mathbb ER(S)\preceq\Lambda I_d.                     \tag{1.4}
\]

Let \(\tau\) be the canonical Stein kernel of \(\rho\).  Then \(Q,r,R,m\)
are locally Lipschitz on \(\operatorname {int}J\).  At their common
differentiability points define

\[
 J_R(s)=\frac12\operatorname {tr}
 \bigl(R^{-1}R'R^{-1}R'\bigr)                         \tag{1.5}
\]

and

\[
 \mathcal G(U)=\int_J\rho(s)\tau(s)^2
 \left\{J_R(s)+m'(s)^TR(s)^{-1}m'(s)\right\}\,ds.     \tag{1.6}
\]

There is a universal constant \(C_{\rm cent}\), independent of \(d\), such
that

\[
 \boxed{\qquad
 \mathcal G(U)\leq1+\frac{C_{\rm cent}\Lambda}{\kappa}.
 \qquad}                                               \tag{1.7}
\]

At almost every \(s\), the centered conditional Ornstein--Uhlenbeck
Poisson field is well defined by the explicit formula in Section 6 below,
and its charges satisfy

\[
 C_s+B_s=J_R(s)+m'(s)^TR(s)^{-1}m'(s).                \tag{1.8}
\]

Thus (1.7) is precisely the nonsmooth Gaussian-fiber weighted-Fisher
estimate, not merely a bound for a relaxed surrogate.

## 2. Coefficient regularity

Fix a compact interval \(K\Subset\operatorname {int}J\).  A finite convex
function is Lipschitz on every compact subset of the interior of its
domain.  In particular, for every fixed \(z\), the section
\(s\mapsto U(s,z)\) is Lipschitz on \(K\).  The coefficients in (1.1) are
recovered from finitely many such sections:

\[
 \begin{aligned}
 c(s)&=U(s,0),\\
 r_i(s)&=\frac12\{U(s,-e_i)-U(s,e_i)\},\\
 Q_{ii}(s)&=U(s,e_i)+U(s,-e_i)-2U(s,0),\\
 Q_{ij}(s)&=U(s,e_i+e_j)-U(s,e_i)-U(s,e_j)+U(s,0)
 \quad(i\ne j).
 \end{aligned}                                        \tag{2.1}
\]

Hence \(c,r,Q\) are locally Lipschitz.  Positivity and continuity of \(Q\)
give a positive lower eigenvalue on \(K\).  The inverse and product rules
for locally Lipschitz matrix functions then show that
\(R=Q^{-1}\) and \(m=Rr\) are locally Lipschitz as well.  All derivatives in
(1.5)--(1.6) therefore exist almost everywhere.

The covariance floor in (1.4) is equivalent to

\[
                         Q(s)\preceq\kappa^{-1}I_d.
 \tag{2.2}
\]

Continuity extends an almost-everywhere version of either inequality to
every interior point.

## 3. Convex quadratic-fiber approximation

Choose compact intervals \(K_k=[a_k,b_k]\) increasing to
\(\operatorname {int}J\).  Let \(\eta\) be a nonnegative, even,
\(C_c^\infty(-1,1)\) probability density.  Choose
\(\varepsilon_k\downarrow0\) so that
\([a_k-\varepsilon_k,b_k+\varepsilon_k]\Subset\operatorname {int}J\), and
put, on \(K_k\times\mathbb R^d\),

\[
 U_k^0(s,z)=\int\eta_{\varepsilon_k}(t)U(s-t,z)\,dt.   \tag{3.1}
\]

Set \(U_k^0=+\infty\) off \(K_k\times\mathbb R^d\).

For every fixed \(t\), the function \((s,z)\mapsto U(s-t,z)\) is convex.
Consequently (3.1), and then its sum with the indicator of the convex strip
\(K_k\times\mathbb R^d\), is convex.  Expanding (3.1) gives

\[
 U_k^0(s,z)=c_k(s)+\frac12z^TQ_k(s)z-r_k(s)^Tz,
 \quad
 (c_k,Q_k,r_k)=\eta_{\varepsilon_k}*(c,Q,r).           \tag{3.2}
\]

Thus the conditional fibers remain exactly Gaussian and all three
coefficients are smooth in the interior of \(K_k\).  Notice that neither
\(R\) nor \(m\) is mollified separately.

Equation (2.2) and positivity of the mollifier give

\[
 Q_k(s)\preceq\kappa^{-1}I_d,\qquad
 R_k(s):=Q_k(s)^{-1}\succeq\kappa I_d.                \tag{3.3}
\]

Because \(\eta\) is even, its barycenter is zero.  Convexity of the scalar
section \(s\mapsto U(s,z)\) and Jensen's inequality give the pointwise
domination

\[
 U_k^0(s,z)\geq U(s,z),\qquad
 {\bf1}_{K_k}(s)e^{-U_k^0(s,z)}\leq e^{-U(s,z)}.       \tag{3.4}
\]

On every compact subset of
\(\operatorname {int}J\times\mathbb R^d\),
\(U_k^0\to U\) locally uniformly.  Dominated convergence in (3.4) therefore
gives convergence of the normalizing constants and, for every function
\(g\) with \(|g|e^{-U}\) integrable,

\[
 \int g(s,z)\,p_k^0(s,z)\,ds\,dz
 \longrightarrow
 \int g(s,z)\,p(s,z)\,ds\,dz,                         \tag{3.5}
\]

provided the assertion is applied separately to the positive and negative
parts.  In particular, (3.5) holds with
\(g=1,s,z_i,s^2,s z_i,z_i z_j\).

The expected conditional covariances also converge.  To see this without
assuming a pointwise upper bound on \(R_k\), let
\(\bar\rho_k^0(s)=\int e^{-U_k^0(s,z)}\,dz\) denote the unnormalized
marginal.  Conditional Jensen gives

\[
 \bar\rho_k^0(s)|m_k(s)|^2
 \leq\int |z|^2e^{-U_k^0(s,z)}\,dz
 \leq\int |z|^2e^{-U(s,z)}\,dz.                       \tag{3.6}
\]

The conditional parameters converge pointwise on every fixed compact
subinterval.  Equations (3.4), (3.6), and dominated convergence yield

\[
 0\preceq R_k(s)
 \preceq\mathbb E_{p_k^0}[Z_kZ_k^T\mid S_k=s],
 \qquad
 \bar\rho_k^0(s)\operatorname {tr}R_k(s)
 \leq\int |z|^2e^{-U_k^0(s,z)}\,dz
 \leq\int |z|^2e^{-U(s,z)}\,dz.                     \tag{3.6a}
\]

Thus the covariance term itself is dominated in the positive-semidefinite
order by the same integrable conditional second-moment density; this is the
uniform-integrability input, rather than centroid Jensen alone.  Using also
the conditional identity
\(R_k=\mathbb E[Z_kZ_k^T\mid S_k]-m_km_k^T\), one obtains

\[
 \mathbb E_{p_k^0}R_k(S)
 =
 \mathbb E_{p_k^0}[ZZ^T]
 -\mathbb E_{p_k^0}[m_k(S)m_k(S)^T]
 \longrightarrow
 \mathbb E_pR(S).                                     \tag{3.7}
\]

All matrix convergence statements here may be read entrywise; in fixed
dimension this is equivalent to operator-norm convergence.

## 4. Exact normalization without changing conditional covariance

Let \((S_k,Z_k)\) have density \(p_k^0\), and define

\[
 \alpha_k=\mathbb ES_k,\qquad
 \sigma_k^2=\operatorname {Var}(S_k),\qquad
 \beta_k=\mathbb EZ_k,
 \tag{4.1}
\]

\[
 T_k=\frac{S_k-\alpha_k}{\sigma_k},\qquad
 \gamma_k=\mathbb E[T_k(Z_k-\beta_k)],\qquad
 Y_k=Z_k-\beta_k-\gamma_kT_k.                         \tag{4.2}
\]

For all sufficiently large \(k\), \(\sigma_k>0\).  Equations (3.5) and
(1.3) imply

\[
 \alpha_k\to0,\qquad \sigma_k\to1,\qquad
 \beta_k\to\mathbb EZ,\qquad \gamma_k\to0.            \tag{4.3}
\]

The map in (4.2) is an invertible affine map.  It preserves log-concavity
and the property of having Gaussian conditional fibers.  By construction,

\[
 \mathbb ET_k=0,\qquad \mathbb ET_k^2=1,\qquad
 \mathbb E[T_kY_k]=0.                                 \tag{4.4}
\]

Its conditional covariance is merely reparametrized:

\[
 \widehat R_k(t)=R_k(\alpha_k+\sigma_kt)\succeq\kappa I_d,
 \qquad
 \mathbb E\widehat R_k(T_k)=\mathbb ER_k(S_k).         \tag{4.5}
\]

The conditional centroid is

\[
 \widehat m_k(t)
 =m_k(\alpha_k+\sigma_kt)-\beta_k-\gamma_kt.          \tag{4.6}
\]

Put
\(\Lambda_k=\|\mathbb E\widehat R_k(T_k)\|_{\rm op}\).
By (3.7),

\[
                         \limsup_k\Lambda_k\leq\Lambda.
 \tag{4.7}
\]

To match explicitly the unit covariance-floor normalization of the smooth
Gaussian-fiber theorem, make the transverse affine change

\[
 \overline Y_k=\kappa^{-1/2}Y_k,\qquad
 \overline R_k=\kappa^{-1}\widehat R_k,\qquad
 \overline m_k=\kappa^{-1/2}\widehat m_k.             \tag{4.8a}
\]

The scalar variable, its marginal density, and its Stein kernel are
unchanged.  Moreover

\[
 \overline R_k\succeq I_d,\qquad
 \left\|\mathbb E\overline R_k(T_k)\right\|_{\rm op}
 =\frac{\Lambda_k}{\kappa},
 \qquad \mathbb E[T_k\overline m_k(T_k)]=0.          \tag{4.8b}
\]

Both quadratic energies are invariant under this transverse change:

\[
 \begin{aligned}
 \frac12\operatorname {tr}
  (\overline R_k^{-1}\overline R_k'
   \overline R_k^{-1}\overline R_k')
 &=\frac12\operatorname {tr}
  (\widehat R_k^{-1}\widehat R_k'
   \widehat R_k^{-1}\widehat R_k'),\\
 \overline m_k'{}^T\overline R_k^{-1}\overline m_k'
 &=\widehat m_k'{}^T\widehat R_k^{-1}\widehat m_k'.
 \end{aligned}                                       \tag{4.8c}
\]

Indeed, \(\overline R_k^{-1}=\kappa\widehat R_k^{-1}\),
\(\overline R_k'=\kappa^{-1}\widehat R_k'\), and
\(\overline m_k'=\kappa^{-1/2}\widehat m_k'\).  Thus applying the
unit-floor theorem to \((T_k,\overline Y_k)\) gives the stated factor
\(\Lambda_k/\kappa\), with no hidden change in \(\mathcal G_k\).

For completeness, its hard-endpoint integration by parts has no boundary
term.  If \([A_k,B_k]\) is the compact support interval of \(T_k\), then
for

\[
 N_k(s)=\widehat\rho_k(s)\widehat\tau_k(s)
       =\int_s^{B_k}t\widehat\rho_k(t)\,dt
\]

one has \(N_k(B_k)=0\), while centering gives
\(N_k(A_k)=\int_{A_k}^{B_k}t\widehat\rho_k(t)\,dt=0\).
Before the strip indicator is imposed, the mollified coefficients extend
smoothly to a neighborhood of both endpoints; their values and all
coefficient expressions occurring in the smooth proof are therefore
finite there.  Every Stein integration-by-parts boundary expression is a
finite coefficient expression multiplied by \(N_k\), and hence vanishes at
both endpoints.  (Equivalently, the positive marginal density has
\(\widehat\tau_k=N_k/\widehat\rho_k=0\) there.)

The smooth Gaussian-fiber theorem therefore applies on the compact interval
which is the image of \(K_k\) under (4.2), and gives

\[
 \mathcal G_k
 :=\int\widehat\rho_k\widehat\tau_k^2
 \left\{\frac12\operatorname {tr}
  (\widehat R_k^{-1}\widehat R_k'
   \widehat R_k^{-1}\widehat R_k')
 +\widehat m_k'{}^T\widehat R_k^{-1}\widehat m_k'\right\}
 \leq1+\frac{C_{\rm cent}\Lambda_k}{\kappa}.          \tag{4.8}
\]

## 5. Passage of the Stein weight and the quadratic energy

We record the convergence needed in (4.8), rather than invoking an
unspecified lower-semicontinuity principle.

Standard one-dimensional mollification and (2.1) give, for every compact
\(L\Subset\operatorname {int}J\) and every finite \(p\),

\[
 Q_k\to Q,\quad r_k\to r
 \quad\hbox{strongly in }W^{1,p}(L).                  \tag{5.1}
\]

The eigenvalues of \(Q_k,Q\) are bounded above and below on \(L\).
Smoothness of matrix inversion and multiplication, followed by (4.3),
therefore gives, after the affine reparametrization,

\[
 \widehat R_k\to R,\qquad
 \widehat m_k'\to m'
 \quad\hbox{strongly in }W^{1,p}(L)
 \ \hbox{and }L^p(L),\ \hbox{respectively}.           \tag{5.2}
\]

More explicitly, the first assertion in (5.2) is strong
\(W^{1,p}\) convergence; the second means
\(\widehat m_k\to m-\mathbb EZ\) strongly in \(W^{1,p}\), whose derivative
is \(m'\).

Let \(\widehat\rho_k\) be the marginal density of \(T_k\).  Equations
(3.4)--(4.3) imply

\[
 \int_{\mathbb R}(1+|s|)\,
       |\widehat\rho_k(s)-\rho(s)|\,ds\longrightarrow0. \tag{5.3}
\]

One way to verify the harmless affine change in (5.3) is first to truncate
to a compact interval, use continuity of translations and dilations in
\(L^1\), and then use the scalar second moments for the two tails.  In fact
the normalization gives exactly

\[
 \int s^2\widehat\rho_k(s)\,ds
 =\int s^2\rho(s)\,ds=1,
\]

and hence, uniformly in \(k\),

\[
 \int_{|s|>R}(1+|s|)
       (\widehat\rho_k(s)+\rho(s))\,ds
 \leq2R^{-2}+2R^{-1}\xrightarrow[R\to\infty]{}0.   \tag{5.3a}
\]

This supplies the weighted-\(L^1\) tail step explicitly.  The Gaussian
integral formula applied to (3.2) also gives
\(\widehat\rho_k\to\rho\) locally uniformly.

For a centered marginal write

\[
 N_k(s)=\widehat\rho_k(s)\widehat\tau_k(s)
       =\int_s^\infty t\widehat\rho_k(t)\,dt,
 \qquad
 N(s)=\rho(s)\tau(s)=\int_s^\infty t\rho(t)\,dt.      \tag{5.4}
\]

Equation (5.3) gives
\(\sup_s|N_k(s)-N(s)|\to0\).  Since \(\rho\) is positive on the interior
of its support, local uniform convergence of the marginal densities yields

\[
 \widehat\rho_k\widehat\tau_k^2
 =\frac{N_k^2}{\widehat\rho_k}
 \longrightarrow
 \frac{N^2}{\rho}
 =\rho\tau^2
 \quad\hbox{locally uniformly on }\operatorname {int}J. \tag{5.5}
\]

Take \(p=2\) in (5.2).  The coefficient matrices and the weights in (5.5)
are uniformly bounded on \(L\), while the derivatives converge strongly in
\(L^2(L)\).  Expanding the two nonnegative quadratic forms therefore gives

\[
 \begin{aligned}
 &\int_L\widehat\rho_k\widehat\tau_k^2
 \left\{\frac12\operatorname {tr}
  (\widehat R_k^{-1}\widehat R_k'
   \widehat R_k^{-1}\widehat R_k')
 +\widehat m_k'{}^T\widehat R_k^{-1}\widehat m_k'\right\}\\
 &\hspace{35mm}\longrightarrow
 \int_L\rho\tau^2
 \left\{J_R+m'^TR^{-1}m'\right\}.                     \tag{5.6}
 \end{aligned}
\]

Exhausting \(\operatorname {int}J\) by compact intervals and using
nonnegativity gives

\[
 \mathcal G(U)
 \leq\liminf_{k\to\infty}\mathcal G_k.                \tag{5.7}
\]

Combining (4.7), (4.8), and (5.7) proves (1.7).  No endpoint or tail energy
has been discarded: (5.7) uses compact exhaustion in the direction needed
for an upper bound on the limiting energy.

## 6. Identification with the conditional Poisson field

At almost every \(s\), all four coefficient functions are differentiable.
Put \(y=z-m(s)\) and let \(A=A^T\) be the unique solution of

\[
 QA+AQ=Q'=-QR'Q,
 \qquad\text{equivalently}\qquad AR+RA=-R'.           \tag{6.1}
\]

For the nonpositive conditional Ornstein--Uhlenbeck generator

\[
 L_s=\Delta_z-\langle Qy,\nabla_z\rangle,
\]

direct differentiation of the Gaussian conditional density shows that

\[
 g_s(z)=-m'(s)^Tz+\frac12y^TAy
\]

solves \(L_sg_s=\partial_s\log q_s\), up to an additive constant.  Hence

\[
 F_s=\nabla_zg_s=-m'(s)+Ay,\qquad
 C_s=\|A\|_{\rm HS}^2,                                \tag{6.2}
\]

\[
 B_s=m'(s)^TQm'(s)+B_s^0,\qquad
 B_s^0=\operatorname {tr}(QARA).                      \tag{6.3}
\]

Diagonalizing \(R\) at this fixed \(s\), without differentiating the
eigenbasis, gives the algebraic identity

\[
 C_s+B_s^0
 =\frac12\operatorname {tr}(R^{-1}R'R^{-1}R')
 =J_R(s).                                             \tag{6.4}
\]

Equations (6.2)--(6.4) prove (1.8) almost everywhere.  Since all terms are
nonnegative and (1.7) is finite, this a.e. field has exactly the required
weighted integrability.

## 7. Exact scope

1. The theorem covers nonsmooth dependence on the scalar slice variable,
   curvature atoms in the marginal, and hard scalar endpoints.
2. The positive covariance floor is essential.  The proof gives
   \(C_{\rm cent}\Lambda/\kappa\); it is not uniform as
   \(\kappa\downarrow0\).
3. The theorem does not approximate non-Gaussian conditional fibers by
   Gaussian ones.  A mixture created by generic smoothing in the transverse
   variables is not covered.
4. Nothing in this note supplies the nonlinear/near-linear spectral
   dichotomy needed to turn weighted-Fisher control into KLS.
