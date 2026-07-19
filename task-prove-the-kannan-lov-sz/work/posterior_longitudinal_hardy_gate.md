# Posterior longitudinal Hardy gate: the exact scalar closure and countermodels

## 0. Verdict

This note addresses the scalar term left open in
`posterior_simple_eigenspace_coherence.md`.  If the selected posterior
eigenline is already fixed, the longitudinal interaction has a complete
dimension-free closure.  More precisely, for a one-dimensional
log-concave probability $q=e^{-u}$ with variance at most two, put

\[
 L=\partial_s^2-u'\partial_s,\qquad h=u''.
\]

Then every smooth compactly supported amplitude $a$ satisfies

\[
 \boxed{
 \begin{aligned}
 &E_q(a'')^2+3E_q[h(a')^2]+2E_q[aa'h']+E_q[h^2a^2]\\
 &\hspace{24mm}\ge {1\over24}
       \{E_q(a')^2+E_q[ha^2]\}.
 \end{aligned}}                                      \tag{0.1}
\]

The constant $1/24$ uses the classical one-dimensional bound
$C_P(q)\le12\operatorname {Var}_q(s)$.  Formula (0.1) is exactly the
scalar part of the twice-Bochner identity, including the signed cubic
term.  Thus that cubic term cannot cancel a one-dimensional low mode.

Three tempting stronger statements are false, even for genuine normalized
Gaussian outputs of isotropic log-concave signals.

1. The longitudinal interaction has no sign for arbitrary exact fields.
2. It cannot be bounded below by a universal multiple of
   $-E[ha^2]$ without derivative terms.
3. The amplitude-biased law $|W|^2d\nu/E|W|^2$ has neither a universal
   Poincare constant nor a projective connectivity inequality controlled
   by its angular Dirichlet energy.  The latter failure persists for
   positive analytic exact gradients over the Gaussian output
   $N(0,2I_2)$, and is exponentially strong.

Consequently the remaining simple-eigenspace gate cannot be replaced by a
generic Hardy inequality for the amplitude tilt.  A valid continuation
must use simultaneously the bottom spectral equation, posterior
low-eigenspace alignment, and covariance normalization.  Sections 1--5
give the complete proofs and isolate a robust scalar inequality that would
be sufficient for the varying-line case.

## 1. Exact one-dimensional factorization

Let $q=e^{-u}$ be a positive $C^4$ log-concave probability density on
the line, and assume that its variance $v$ is finite.  On
$L^2(q)$, define

\[
 D={d\over ds},\qquad D^*=-{d\over ds}+u',
 \qquad \mathcal A_1=DD^*=-L+h,qquad h=u''.
                                                               \tag{1.1}
\]

For $a\in C_c^\infty(\mathbb R)$, put $g=D^*a$.  Adjointness gives

\[
 E_qg=\langle 1,D^*a\rangle=\langle D1,a\rangle=0,
 \qquad Dg=\mathcal A_1a.                         \tag{1.2}
\]

The Poincare inequality for $q$, applied to $g$, therefore gives the
factorized Hardy estimate

\[
 \boxed{
 \|\mathcal A_1a\|_{L^2(q)}^2
 \ge C_P(q)^{-1}\langle a,\mathcal A_1a\rangle
 =C_P(q)^{-1}\{E_q(a')^2+E_q[ha^2]\}.}            \tag{1.3}
\]

There is no spectral-attainment assumption in this argument.  Expanding
the square on the left and integrating once gives

\[
 \boxed{
 \begin{aligned}
 \|\mathcal A_1a\|_2^2
 &=E_q(a'')^2+3E_q[h(a')^2]\\
 &\quad+2E_q[aa'h']+E_q[h^2a^2].
 \end{aligned}}                                    \tag{1.4}
\]

Indeed, one-dimensional Bochner gives
$E_q(La)^2=E_q(a'')^2+E_q[h(a')^2]$, while

\[
 2\langle-L a,ha\rangle
 =2E_q[a'(ha)']=2E_q[h(a')^2]+2E_q[aa'h'].
\]

The classical one-dimensional log-concave Poincare theorem (the
one-dimensional variance bound in S. Bobkov, *Ann. Probab.* 27 (1999),
1903--1921) states, for every log-concave probability on the line with
finite variance and with no smoothness or strict-convexity requirement,
that

\[
 C_P(q)\le12\operatorname {Var}_q(s).              \tag{1.5}
\]

Applying (1.5) when $v\le2$ proves (0.1).  Approximation in the graph
norm of $\mathcal A_1$ extends (1.3)--(1.4) to its operator domain.  In
particular, if $\mathcal A_1a=\lambda a\ne0$, then

\[
 \lambda^2\|a\|_2^2
 \ge {1\over24}\lambda\|a\|_2^2,
 \qquad\text{hence}\qquad \lambda\ge {1\over24}. \tag{1.6}
\]

This is precisely the needed longitudinal absorption in the exact fixed
line case.  It is not circular: (1.5) is the elementary, dimension-one
log-concave theorem, and the variance two is fixed by
$Y=X+G$ when $\operatorname {Var}(X)=1$.

There is also a useful integration-by-parts form of the signed term:

\[
 \boxed{2E_q[aa'h']=-E_q[a^2Lh].}                 \tag{1.7}
\]

Thus a pointwise sign or a curvature-only estimate for $Lh$ would be a
stronger route to (0.1).  Section 3 shows that such estimates are false
for an actual Gaussian posterior.

## 2. What a fixed posterior eigenline forces

Let $d\nu=e^{-U(y)}dy$ on \(\mathbb R^k\), and suppose a fixed unit
vector $e$ is an eigenvector of $H=D^2U$ at every point.  Write
$y=se+z$, with $z\perp e$.  Since

\[
 D^2U(y)[e,v]=0\qquad(v\perp e),                  \tag{2.1}
\]

the function \(\partial_sU(s,z)\) is independent of $z$.  Hence, up
to an additive affine term which is absorbed into the two factors,

\[
 U(s,z)=u(s)+V(z),\qquad \nu=q\otimes\eta.         \tag{2.2}
\]

If an exact field is of the form $W=a(y)e=\nabla f(y)$, then the
transverse components of \(\nabla f\) vanish.  Therefore $f=f(s)$ and
$a=a(s)$.  If additionally \(\operatorname {Cov}(\nu)=2I\), then
$\operatorname {Var}_q(s)=2$, and (0.1) applies verbatim.  Thus the
longitudinal cubic in (5.5) of
`posterior_simple_eigenspace_coherence.md` is completely controlled once
the low line is fixed.  The unresolved issue is quantitative passage from
small amplitude-weighted rotation to a fixed line, not a scalar posterior
inequality after that passage.

## 3. A normalized posterior with no curvature sign

Let $E$ be a mean-one exponential random variable, set
$X=E-1$, and let $G\sim N(0,1)$ be independent.  Then $X$ is centered,
log-concave, and has variance one, while $Y=X+G$ has variance two.  Its
density is the explicit analytic Gaussian output

\[
 q(y)=e^{-y-1/2}\Phi(y),                           \tag{3.1}
\]

where \(\Phi\) and \(\phi\) are the standard Gaussian distribution and
density.  With

\[
 r(y)={\phi(y)\over\Phi(y)},\qquad
 U(y)=y+\tfrac12-\log\Phi(y),
\]

one has

\[
 U'=1-r,\qquad h=U''=-r'=r(y)(y+r(y)),
 \qquad 0<h<1.                                    \tag{3.2}
\]

This is also the posterior identity: conditional on $Y=y$, the variable
$E=X+1$ is $N(y,1)$ truncated to $[0,\infty)$, and its variance is
$1-h(y)$.

The right-tail asymptotics $r(y)=\phi(y)(1+o(1))$ give

\[
 h(y)=y\phi(y)(1+o(1)),\quad
 h'(y)=-y^2\phi(y)(1+o(1)),\quad
 h''(y)=y^3\phi(y)(1+o(1)).                       \tag{3.3}
\]

Consequently

\[
 Lh=h''-(1-r)h'=y^3\phi(y)(1+o(1))>0,
 \qquad {Lh\over h}=y^2(1+o(1))\longrightarrow\infty. \tag{3.4}
\]

At the other end, writing $A=-y\to\infty$, the inverse Mills expansion

\[
 r(-A)=A+A^{-1}-2A^{-3}+10A^{-5}+O(A^{-7})
\]

gives

\[
 \begin{aligned}
 h(-A)&=1-A^{-2}+6A^{-4}+O(A^{-6}),\\
 h'(-A)&=-2A^{-3}+24A^{-5}+O(A^{-7}),\\
 h''(-A)&=-6A^{-4}+O(A^{-6}),
 \end{aligned}                                    \tag{3.5}
\]

and therefore

\[
 Lh(-A)=-2A^{-2}+O(A^{-3})<0.                     \tag{3.6}
\]

Thus even the weighted Laplacian of the smallest posterior curvature
changes sign.  More quantitatively, for every $C>0$, (3.4) supplies an
interval on which $Lh>2Ch$.  Choosing a nonzero
$a\in C_c^\infty$ supported in that interval and using (1.7) yields

\[
 2E_q[aa'h']=-E_q[a^2Lh]<-2C E_q[ha^2].           \tag{3.7}
\]

Hence there is no universal curvature-only lower estimate

\[
 2E[aa'h']\ge-C E[ha^2]                           \tag{3.8}
\]

for Gaussian posteriors.  The $a''$ and $h(a')^2$ terms in (0.1) are
essential.

There is an even simpler sign test.  For any nonzero real $c$, let

\[
 a_c(y)=e^{c h(y)},\qquad
 f_c(y)=\int_0^y a_c(t)\,dt.
\]

Then $W=f_c'=a_c$ is a positive analytic exact field and, in the scalar
notation of the twice-Bochner formula,

\[
 \alpha=a_c'=c h'a_c,\qquad
 t=D^3U[1,1,W]=a_ch'.
\]

Therefore

\[
 \boxed{E_q[\alpha t]=cE_q[a_c^2(h')^2],}         \tag{3.9}
\]

which has either sign.  These fields are not bottom eigenfields; (3.9)
shows exactly why the eigen-equation must be used by any sign argument.

## 4. The amplitude tilt need not satisfy a Hardy or Poincare bound

The failure already occurs for the normalized Gaussian output

\[
 \nu=N(0,2),\qquad H=\tfrac12.
\]

This is the law of $X+G$ for $X,G\sim N(0,1)$.  Let
$a_R(x)=\cosh(Rx)$ and

\[
 d\sigma_R={a_R^2\,d\nu\over E_\nu a_R^2}.
\]

The law is symmetric, and Gaussian moment generating functions give

\[
 E_\nu a_R^2={e^{4R^2}+1\over2},
\]

\[
 \operatorname {Var}_{\sigma_R}(X)
 ={(2+16R^2)e^{4R^2}+2\over e^{4R^2}+1}.          \tag{4.1}
\]

Testing the Poincare inequality of \(\sigma_R\) with $x\mapsto x$
therefore yields

\[
 C_P(\sigma_R)\ge \operatorname {Var}_{\sigma_R}(X)
 \sim16R^2.                                       \tag{4.2}
\]

Thus no universal Poincare or two-sided Hardy constant is inherited by
the amplitude tilt, even when the base output is strongly log-concave and
the amplitude is positive and analytic.  In one dimension every such
amplitude is the derivative of a smooth function, so exactness alone does
not repair this failure.

## 5. Exponentially bad projective connectivity for an exact gradient

There is a sharper two-dimensional obstruction.  Let

\[
 \nu=N(0,2I_2),\qquad
 f_R(x,y)={e^{Rx}+e^{Ry}\over R},\qquad
 W_R=\nabla f_R=(a,b)=(e^{Rx},e^{Ry}).             \tag{5.1}
\]

This is again an actual normalized Gaussian output, now from
$X\sim N(0,I_2)$, and $W_R$ is positive, analytic, nonvanishing, and
exact.  Since

\[
 E_\nu|W_R|^2=2e^{4R^2},
\]

define

\[
 d\sigma_R={|W_R|^2d\nu\over2e^{4R^2}},\qquad
 N_R={W_RW_R^T\over|W_R|^2}.
\]

For any fixed rank-one projection $P_0=uu^T$,

\[
 E_{\sigma_R}\|N_R-P_0\|_{HS}^2
 =2-2u^T(E_{\sigma_R}N_R)u.                       \tag{5.2}
\]

Symmetry gives diagonal entries $1/2$ in $E_{\sigma_R}N_R$, while

\[
 (E_{\sigma_R}N_R)_{12}
 ={E_\nu[ab]\over2e^{4R^2}}
 ={1\over2}e^{-2R^2}.                              \tag{5.3}
\]

It follows exactly that

\[
 \boxed{
 \inf_{P_0:\,\operatorname {rank}P_0=1}
 E_{\sigma_R}\|N_R-P_0\|_{HS}^2
 =1-e^{-2R^2}.}                                    \tag{5.4}
\]

On the other hand, direct differentiation gives

\[
 \sum_{j=1}^2\|D_jN_R\|_{HS}^2
 ={4R^2a^2b^2\over(a^2+b^2)^2}.                  \tag{5.5}
\]

With $S=(x+y)/\sqrt2$ and $D=(x-y)/\sqrt2$, which are independent
$N(0,2)$ variables,

\[
 {a^2b^2\over a^2+b^2}
 ={e^{\sqrt2RS}\over2\cosh(\sqrt2RD)}.
\]

Consequently

\[
 \boxed{
 E_{\sigma_R}\sum_j\|D_jN_R\|_{HS}^2
 =R^2e^{-2R^2}E[\operatorname {sech}(\sqrt2RD)]
 \le R^2e^{-2R^2}.}                               \tag{5.6}
\]

Equations (5.4)--(5.6) disprove a universal amplitude-weighted
projective Poincare inequality

\[
 \inf_{P_0}E_\sigma\|N-P_0\|_{HS}^2
 \le C E_\sigma\|DN\|_{HS}^2.                   \tag{5.7}
\]

Adding a polynomial penalty in the local amplitude Fisher information
does not help.  Indeed, with $A_R=|W_R|$,

\[
 |\nabla\log A_R|^2
 =R^2{a^4+b^4\over(a^2+b^2)^2}\le R^2.           \tag{5.8}
\]

For every fixed $m<\infty$, the product of
$(1+E_{\sigma_R}|\nabla\log A_R|^2)^m$ and the right side of (5.6)
still tends to zero, while (5.4) tends to one.

The exact mixture interpretation makes the obstruction transparent:
\(\sigma_R\) is the equally weighted mixture of
$N((4R,0),2I_2)$ and $N((0,4R),2I_2)$.  The field points almost along
the first coordinate in the first packet and along the second coordinate
in the second packet; the transition has exponentially small
amplitude-biased mass.

This example does not satisfy the low-energy/simple-eigengap hypotheses:
$H=I/2$ has multiplicity two and
$E\|DW_R\|^2/E|W_R|^2=R^2$.  It therefore does not challenge the
simple-eigenspace program.  It does prove that neither exactness,
analyticity, Gaussian posterior structure, angular energy, nor finitely
many local amplitude moments can be the missing connectivity input by
itself.

## 6. The precise robust scalar target

Retain the notation of
`posterior_simple_eigenspace_coherence.md`: $P=ee^T$,
$W=ae+r$ with $r\perp e$, $C=DW$,

\[
 \alpha=e^TCe,\qquad t=D^3U[e,e,W]=Dh[W].
\]

The proved scalar inequality (0.1) says that when $P$ is fixed,
$r=0$, and exactness makes $a=a(e\cdot y)$,

\[
 E(a'')^2+3E[h\alpha^2]+2E[\alpha t]+E[h^2a^2]
 \ge {1\over24}E[\alpha^2+ha^2].                 \tag{6.1}
\]

A quantitatively stable version of (6.1) would close the longitudinal
branch.  One formal sufficient statement is the following: for each fixed
simple eigengap $\kappa>0$, prove dimension-free constants
$c_\kappa>0$, $C_\kappa<\infty$ such that

\[
 \begin{aligned}
 &E\|DC\|_{HS}^2+3E[h\alpha^2]+2E[\alpha t]
       +E[h^2a^2]\\
 &\quad\ge c_\kappa E[\alpha^2+ha^2]
       -C_\kappa\,\mathfrak T(W,P),               \tag{6.2}
 \end{aligned}
\]

where \(\mathfrak T(W,P)\) consists only of transverse terms already
coercive in (5.3)--(5.4), for example squared $C_0=C-\alpha P$, $r$,
and their posterior-weighted rotation terms, with coefficients small
enough to be absorbed by those displayed transverse positive terms.
Equation (6.1) proves (6.2) with \(\mathfrak T=0\) in the fixed-line
case.

The countermodels above impose two nonnegotiable requirements on a proof
of (6.2): it must use the bottom spectral equation (not arbitrary exact
fields), and it must use the low-curvature simple-eigenspace alignment
(not merely the amplitude tilt).  Proving (6.2), or an equivalent
amplitude-weighted nodal-capacity theorem under precisely those joint
hypotheses, remains the gate.  No full KLS claim is made here.
