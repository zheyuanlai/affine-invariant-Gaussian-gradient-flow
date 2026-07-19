# Additive localization windows and the non-rigid Gaussian baseline

## 1. Exact scalar reduction

Let \(g_t\) be the mass of a fixed balanced set under ordinary Gaussian
stochastic localization at precision time \(t\), and set

\[
 F(t)=\mathbb E I(g_t).
\]

Under the heat/localization correspondence, \(F(t)=J(1/t)\). If
\(r_t=|\operatorname {Cov}_{p_t}(\mathbf1_S,X)|\), put

\[
 \eta_t={\sqrt t\,r_t\over I(g_t)}\in[0,1].
\]

Itô's formula and \(I''=-1/I\) give

\[
 -F'(t)={1\over2}\mathbb E{r_t^2\over I(g_t)}
        ={1\over2t}\mathbb E[I(g_t)\eta_t^2].       \tag{1.1}
\]

Fix a reference scale \(K>0\), write \(\tau=Kt\), and normalize by
\(P_0=I(1/2)\):

\[
 b(\tau)={F(\tau/K)\over P_0},\qquad
 d\nu_\tau={I(g_{\tau/K})\over F(\tau/K)}d\mathbb P,
\qquad m(\tau)=\int\eta^2d\nu_\tau.                 \tag{1.2}
\]

Then

\[
 \boxed{m(\tau)=-{2\tau b'(\tau)\over b(\tau)}\in[0,1].}       \tag{1.3}
\]

The profile slack from the heat identity is

\[
 A=\mathbb E[I(g_t)(1-\eta_t^2)]
   =P_0\{b+2\tau b'\}=P_0b(1-m).                    \tag{1.4}
\]

Consequently the normalized endpoint functional in the Bernstein identity is

\[
             q_b(\tau)={b(\tau)(1-m(\tau))\over\sqrt\tau}.    \tag{1.5}
\]

Changing variables \(s=1/t=K/\tau\) in (4.5) of
heatflow_bernstein.md gives the exact decomposition

\[
 \boxed{
 -q_b'(\tau)
 =\mathcal G(\tau)
  +{b(\tau)\over2\tau^{3/2}}
       \int(1-\eta^2)^2d\nu_\tau,}                  \tag{1.6}
\]

where

\[
 \mathcal G(\tau)
 ={K^2\over P_0\tau^{7/2}}
   \int\rho_{K/\tau}
   \left(\|\nabla^2z\|_{HS}^2
       -2\langle\nabla^2\log q\,\nabla z,\nabla z\rangle\right)
 \ge0.                                               \tag{1.7}
\]

The first summand in (1.7) contains the angular/projected-Hessian energy;
the second is the nonnegative curvature energy of the observation density.
Jensen in (1.6) yields

\[
 -q_b'(\tau)
 \ge {b(\tau)(1-m(\tau))^2\over2\tau^{3/2}}.        \tag{1.8}
\]

Retain the exact nonnegative scalar remainder

\[
 \begin{aligned}
 \mathcal E_b(\tau)
 &:=-q_b'(\tau)
  -{b(\tau)(1-m(\tau))^2\over2\tau^{3/2}}\\
 &=\mathcal G(\tau)
   +{b(\tau)\over2\tau^{3/2}}
           \operatorname {Var}_{\nu_\tau}(\eta^2)\ge0.        \tag{1.9}
 \end{aligned}
\]

Because \(b'/b=-m/(2\tau)\), direct differentiation of (1.5) shows that
(1.8) is equivalent to the weak first-order condition

\[
 \boxed{m'(\tau)\ge-{m(\tau)(1-m(\tau))\over\tau}.}  \tag{1.10}
\]

Thus the scalar information in the full sum-of-squares identity is the pair
of constraints (1.3) and (1.10), together with a nonnegative remainder. It
does not distinguish curvature energy from angular energy.

## 2. Gaussian halfspace and the additive-window sequence

Take \(\mu=N(0,KI)\) and a centered halfspace. A direct one-dimensional
posterior computation gives

\[
 F(t)={P_0\over\sqrt{1+Kt}},\qquad
 b_G(\tau)=(1+\tau)^{-1/2},\qquad
 m_G(\tau)={\tau\over1+\tau}.                       \tag{2.1}
\]

At the additive windows \(t_j=j/K\), define

\[
 a_j={F(t_j)-F(t_{j+1})\over P_0}.
\]

Then exactly

\[
 \boxed{
 a_j={1\over\sqrt{1+j}}-{1\over\sqrt{2+j}},\qquad
 \sum_{j=0}^{\infty}a_j=1.}                        \tag{2.2}
\]

In particular

\[
 a_j={1+o(1)\over2j^{3/2}},\qquad
 j a_j={1+o(1)\over2\sqrt j}.                       \tag{2.3}
\]

This is the exact match between the harmonic relative defect \(j a_j\)
and a bridge cost of order \(j^{-1/2}\).

However, the Gaussian is not an equality case of the scalar inequality
(1.8). From (2.1),

\[
 q_G(\tau)={1\over\sqrt\tau(1+\tau)^{3/2}},
\]

and explicit subtraction gives

\[
 \boxed{
 \mathcal E_{b_G}(\tau)
 ={2\over\sqrt\tau(1+\tau)^{5/2}}>0.}               \tag{2.4}
\]

For the Gaussian halfspace, \(\nabla^2z=0\) and \(\eta^2\) is deterministic
under the boundary weight. Hence all of (2.4) is already the
one-dimensional curvature term in (1.7); none is angular phase energy. Any
proposed subtraction of a Gaussian baseline must therefore calibrate the
curvature term, not merely the scalar-square term.

## 3. Scalar variational problem and non-uniqueness

For a prescribed tail coefficient \(\lambda>0\), the natural scalar
admissible class is

\[
 \begin{aligned}
 \mathfrak A_\lambda=\{b\in C^2(0,\infty):{}&
 b(0+)=1,\quad b>0,\quad b'\le0,\\
 &\sqrt\tau b(\tau)\longrightarrow\lambda,\\
 &0\le m_b:=-2\tau b'/b\le1,\quad
 \mathcal E_b\ge0\}.                                \tag{3.1}
\end{aligned}
\]

For a finite-perimeter cut, the small-heat-time formula gives
\(\lambda=P_\mu(S)\sqrt K/P_0\), so this endpoint datum is precisely the
normalized original perimeter; it is not an extra regularity assumption.

For \(b\in\mathfrak A_\lambda\), put

\[
 a_j(b)=b(j)-b(j+1),\qquad \delta_j(b)=j a_j(b).     \tag{3.2}
\]

The hoped-for baseline principle would require the endpoint data in (3.1),
possibly together with finitely many seed values, to force \(b=b_G\) when
\(\lambda=1\), or at least to fix the integral of the non-angular part of
\(\mathcal E_b\). They do neither.

Choose a nonzero \(\phi\in C_c^\infty((2,4))\) and define

\[
                  b_\epsilon(\tau)=b_G(\tau)e^{\epsilon\phi(\tau)}.
                                                               \tag{3.3}
\]

For all sufficiently small \(|\epsilon|>0\), depending only on the
\(C^2\)-norm of \(\phi\),

\[
                  b_\epsilon\in\mathfrak A_1.        \tag{3.4}
\]

Indeed, outside \([2,4]\) it equals \(b_G\). On \([2,4]\), the strict
inequalities \(b_G'<0\), \(0<m_G<1\), and
\(\mathcal E_{b_G}>0\) have positive uniform margins. The maps

\[
 b\mapsto b',\qquad b\mapsto-2\tau b'/b,\qquad
 b\mapsto\mathcal E_b
\]

are continuous in the positive \(C^2\) topology on that compact interval,
so all defining inequalities persist for small \(\epsilon\). The initial
value, every seed value below time \(2\), and the asymptotic coefficient are
unchanged. If \(\phi(3)\ne0\), at least the adjacent increments \(a_2,a_3\)
change, by opposite amounts in their telescoping sum.

It follows that even seed data on an initial interval plus the sharp
asymptotic \(\sqrt\tau b(\tau)\to1\) do not select the Gaussian sequence
(2.2). The fixed total budget

\[
                  \sum_{j\ge0}a_j=b(0)-b(\infty)=1             \tag{3.5}
\]

can absorb a local positive angular contribution by a compensating change
of the scalar profile on nearby windows. In differential form, integration
of (1.9) from \(\tau_0\) to infinity gives

\[
 \int_{\tau_0}^{\infty}\mathcal E_b
 =q_b(\tau_0)
  -\int_{\tau_0}^{\infty}
       {b(1-m_b)^2\over2\tau^{3/2}}d\tau,            \tag{3.6}
\]

and neither term on the right is fixed by endpoint data.

## 4. Consequence for the power mismatch

The additive partition identifies the sharp obstruction but does not improve
it. For the Gaussian profile, \(a_j\asymp j^{-3/2}\), so the local quantity
\(j a_j\) is exactly of order \(j^{-1/2}\). The scalar admissible class
contains an open \(C^2\) neighborhood of this profile on every compact
interval. Therefore the Bernstein sum of squares, seed, and asymptotic data
alone do not prove that angular phase energy creates a strict excess over a
fixed Gaussian budget. No improvement over the existing
\(\sqrt\alpha\)-versus-\(\alpha\) mismatch follows without an additional
theorem that calibrates the curvature term or couples different windows more
rigidly.
