# Geometry and coarea of the terminal parallel-coupling alignment

## 0. Verdict

Let \(G_{t,B}(x)\) be the parallel-coupling flow from
`/tmp/parallel_src/thin-shell-arXiv-v2.tex`.  For a fixed Brownian path
\(B\), the terminal map

\[
 \Phi_B(x):=G_{1,B}(x)
 \tag{0.1}
\]

is a smooth orientation-preserving diffeomorphism of the *initial
exponential-tilt variable* \(x\).  Its forward Jacobian at the zero tilt is

\[
 \boxed{M_1=D_x\Phi_B(0).}
 \tag{0.2}
\]

It is not the inverse Jacobian, and it is not the derivative with respect
to the Brownian endpoint.

Define the time-one posterior scalar

\[
 F(\theta)=\mathbb E_{\mu_{1,\theta}}f.
 \tag{0.3}
\]

Then

\[
 \boxed{\nabla_\theta F(\theta)
 =\operatorname {Cov}_{\mu_{1,\theta}}(X,f).}
 \tag{0.4}
\]

Writing \(\Theta=\Phi_B(0)\), we have

\[
 c_1=\nabla F(\Theta).
\]

Consequently the proposed terminal alignment quantity is exactly

\[
 \boxed{
 {\langle M_1b,c_1\rangle^2\over|c_1|^2}
 =|\operatorname {proj}_{N(\Theta)}D\Phi_B(0)b|^2,\qquad
 N={\nabla F\over|\nabla F|}.}
 \tag{0.5}
\]

It is the squared normal stretch of the initial direction \(b\) through the
random tilt flow, measured against the terminal level set of \(F\).

There are exact change-of-variables, Stein, and coarea identities for this
geometry, but all have the wrong coercive direction for an upper bound:

* Change of variables controls the determinant divided by the tangential
  Jacobian, while (0.5) can concentrate in one normal singular direction.
* The exact transition-law integration by parts controls the conditional
  *mean* of \(M_1b\).  The alignment quotient is a conditional second
  moment.
* Coarea and the divergence theorem give a lower bound on the normal energy
  from its flux, not an upper bound.
* A second differentiation of the transition law produces the second
  derivative \(D_x^2\Phi_B(0)[b,b]\) with no sign.
* Wiener integration by parts does not see \(M_1\) at finite Cameron--Martin
  cost: changing the initial tilt is an impulse at time zero.

The obstruction already appears for expansive linear maps.  With

\[
 L_K=\operatorname {diag}(K,1,\ldots,1),\qquad
 F(y)=y_1,\qquad b=e_1,
\]

the alignment is \(K^2\).  Taking \(K=\sqrt n\) leaves
\(\operatorname {Tr}(L_K^TL_K)<2n\) and
\((\det L_K)^{1/n}\to1\), while the alignment grows like \(n\).
The level sets themselves are flat hyperplanes with unchanged tangential
area.  Thus no inequality based only on determinant, trace per dimension,
or level-set area can prove the desired bound.

The expectation in (0.5) is specifically over Wiener paths at the fixed
initial tilt \(x=0\).  Confusing it with Lebesgue integration over \(x\), or
with integration over the endpoint \(\theta\) without conditioning on the
path, is the invalid change-of-variables step.  The geometric reformulation
is exact, but it does not close the final-alignment lemma.

## 1. The forward tilt flow and its Jacobian

Recall

\[
 \Lambda_t(\theta)
 =\log\int\exp\left(\theta\cdot z-{t\over2}|z|^2\right)d\mu(z),
 \tag{1.1}
\]

and set

\[
 a(t,\theta)=\nabla\Lambda_t(\theta),
 \qquad
 A(t,\theta)=D^2\Lambda_t(\theta)
 =\operatorname {Cov}_{\mu_{t,\theta}}(X).
 \tag{1.2}
\]

For a fixed continuous path \(B\), \(G_{t,B}(x)\) is the solution of

\[
 G_{t,B}(x)
 =x+B_t+\int_0^t a(s,G_{s,B}(x))ds.
 \tag{1.3}
\]

For compactly supported \(\mu\), the vector field is smooth and globally
Lipschitz.  Thus \(x\mapsto G_{t,B}(x)\) is a global diffeomorphism; the
inverse is obtained by solving the corresponding ODE backward.

Let

\[
 M_t(x)=D_xG_{t,B}(x).
\]

Differentiating (1.3) gives

\[
 M_0(x)=I,
 \qquad
 {d\over dt}M_t(x)
 =A(t,G_{t,B}(x))M_t(x).
 \tag{1.4}
\]

The matrix in the low-mode alignment is

\[
 M_1=M_1(0)=D\Phi_B(0).
 \tag{1.5}
\]

Because \(A\succeq0\),

\[
 {d\over dt}|M_t(x)v|^2
 =2\langle M_t(x)v,A(t,G_{t,B}(x))M_t(x)v\rangle\ge0.
 \tag{1.6}
\]

Hence every singular value of \(M_t(x)\) is at least one.  Liouville's
formula gives the exact volume Jacobian

\[
 \boxed{
 J_B(x):=\det D\Phi_B(x)
 =\exp\left(
 \int_0^1\operatorname {Tr}A(t,G_{t,B}(x))dt
 \right)>0.}
 \tag{1.7}
\]

Although each \(A(t,G_t(x))\) is symmetric, the matrices at different
times need not commute.  Therefore \(M_t\) need not be symmetric and
\(\Phi_B\) is not known to be the gradient of a convex potential.  A
Monge--Ampere change of variables is unavailable; (1.7) is the applicable
Jacobian identity.

When changing variables \(y=\Phi_B(x)\), the inverse derivative is

\[
 D_y\Phi_B^{-1}(y)=D_x\Phi_B(x)^{-1},
 \qquad x=\Phi_B^{-1}(y),
 \tag{1.8}
\]

and

\[
 dx={dy\over J_B(\Phi_B^{-1}(y))}.
 \tag{1.9}
\]

Thus \(M^{-1}\) occurs only after passing to a preimage.  The alignment
itself uses the forward derivative \(M_1\) at \(x=0\).

## 2. The posterior scalar and the normal-stretch interpretation

At time one, define

\[
 F(\theta)
 ={\int f(z)e^{\theta\cdot z-|z|^2/2}d\mu(z)
   \over
   \int e^{\theta\cdot z-|z|^2/2}d\mu(z)}.
 \tag{2.1}
\]

Differentiation under the integral gives

\[
 \nabla F(\theta)
 =\mathbb E_{1,\theta}
 [(X-\mathbb E_{1,\theta}X)(f-\mathbb E_{1,\theta}f)].
 \tag{2.2}
\]

Let

\[
 \Theta=G_{1,B}(0).
\]

Then \(c_1=\nabla F(\Theta)\).  Put

\[
 H_B(x)=F(\Phi_B(x)).
 \tag{2.3}
\]

The chain rule yields

\[
 \nabla_xH_B(x)
 =D\Phi_B(x)^T\nabla F(\Phi_B(x)).
 \tag{2.4}
\]

In particular,

\[
 \partial_bH_B(0)
 =\langle M_1b,c_1\rangle.
 \tag{2.5}
\]

At a regular endpoint \(\theta\), let

\[
 N(\theta)={\nabla F(\theta)\over|\nabla F(\theta)|}
\]

be the unit normal to the level set of \(F\).  Then

\[
 {\langle M_1b,c_1\rangle^2\over|c_1|^2}
 =\langle D\Phi_B(0)b,N(\Theta)\rangle^2.
 \tag{2.6}
\]

Equivalently, the pullback conormal is

\[
 M_1^TN(\Theta),
\]

and (2.6) is its squared component in the initial direction \(b\).
Critical points of \(F\) contribute zero by convention.

## 3. Which probability measure is being used?

The expectation of interest is

\[
 \mathcal Q_b
 =\mathbb E_B\left[
 \langle M_1b,N(\Theta)\rangle^2
 \mathbf1_{\{|\nabla F(\Theta)|>0\}}
 \right],
 \tag{3.1}
\]

where:

1. \(B\) is a standard Brownian path;
2. the initial tilt is fixed at \(x=0\);
3. \(M_1=D_xG_{1,B}(0)\);
4. \(\Theta=G_{1,B}(0)\).

It is not an expectation over initial tilts \(x\), and it is not an
unconditional Lebesgue integral over terminal tilts.

The source proves the distributional identity

\[
 (G_{t,B}(x))_{t\ge0}
 \stackrel{d}{=}(x+B_t+tX_x)_{t\ge0},
 \tag{3.2}
\]

where \(X_x\) has the exponentially tilted law

\[
 d\mu_x(z)=e^{x\cdot z-\Lambda_0(x)}d\mu(z).
\]

At \(x=0,t=1\),

\[
 \Theta\stackrel d=X+G,
 \qquad X\sim\mu,quad G\sim N(0,I),quad X\perp G.
 \tag{3.3}
\]

Thus \(\Theta\) has density

\[
 q(\theta)
 =(2\pi)^{-n/2}e^{-|\theta|^2/2+\Lambda_1(\theta)}.
 \tag{3.4}
\]

The scalar \(F\) is precisely the posterior expectation

\[
 F(\theta)=\mathbb E[f(X)\mid X+G=\theta].
 \tag{3.5}
\]

However, \(M_1\) is generally not determined by \(\Theta\): it depends on
the entire Brownian path.  Define

\[
 \bar m_b(\theta)
 =\mathbb E[M_1b\mid\Theta=\theta]
 \tag{3.6}
\]

and

\[
 R_b(\theta)
 =\mathbb E[\langle M_1b,N(\theta)\rangle^2
             \mid\Theta=\theta].
 \tag{3.7}
\]

Then the exact endpoint representation is

\[
 \boxed{
 \mathcal Q_b=\int_{\mathbb R^n}q(\theta)R_b(\theta)d\theta.}
 \tag{3.8}
\]

The conditional decomposition

\[
 R_b(\theta)
 =\langle\bar m_b(\theta),N(\theta)\rangle^2
 +\operatorname {Var}(\langle M_1b,N(\theta)\rangle
                       \mid\Theta=\theta)
 \tag{3.9}
\]

shows the central issue: endpoint integration by parts identifies the mean
field \(\bar m_b\), but not the conditional variance in (3.9).

## 4. Exact transition density and Stein identity

The law of \(\Phi_B(x)=G_{1,B}(x)\) has a simple density.  From (3.2),

\[
 \Phi_B(x)\stackrel d=x+G+X_x.
\]

A direct completion of squares gives

\[
 \boxed{
 p_x(\theta)
 =q(\theta)
 \exp\left(
 x\cdot\theta-{|x|^2\over2}-\Lambda_0(x)
 \right).}
 \tag{4.1}
\]

Since \(\mu\) is centered and isotropic,

\[
 \nabla\Lambda_0(0)=0,
 \qquad D^2\Lambda_0(0)=I.
 \tag{4.2}
\]

Differentiate

\[
 \mathbb E_B\psi(\Phi_B(x))
 =\int\psi(\theta)p_x(\theta)d\theta
\]

at \(x=0\) in direction \(b\).  On the flow side the derivative is
\(\nabla\psi(\Theta)\cdot M_1b\); on the density side it is
\((b\cdot\theta)q(\theta)\).  Therefore

\[
 \boxed{
 \mathbb E_B[\nabla\psi(\Theta)\cdot M_1b]
 =\int(b\cdot\theta)\psi(\theta)q(\theta)d\theta.}
 \tag{4.3}
\]

Conditioning on \(\Theta\),

\[
 \int q\,\bar m_b\cdot\nabla\psi
 =\int q\,(b\cdot\theta)\psi.
\]

Equivalently,

\[
 \boxed{-\operatorname {div}(q\bar m_b)
 =(b\cdot\theta)q.}
 \tag{4.4}
\]

Thus \(\bar m_b\) is a Stein field for the convolved measure \(q\), not for
the original measure \(\mu\).

Taking \(\psi=F\) in (4.3) gives

\[
 \mathbb E\langle M_1b,c_1\rangle
 =\mathbb E[(b\cdot\Theta)F(\Theta)]
 =\mathbb E[(b\cdot(X+G))f(X)]
 =b\cdot\mathbb E[Xf].
 \tag{4.5}
\]

For \(b=a/|a|\), the last expression is \(|a|\), recovering the low-mode
martingale pairing.  Equation (4.3) is a first-moment identity.  By Jensen,

\[
 \langle\bar m_b,N\rangle^2\le R_b,
 \tag{4.6}
\]

so it supplies a lower bound on the desired second moment, not an upper
bound.

At the level of the endpoint identities this failure is sharp.  Replacing
the conditional vector by

\[
 M_1b\longmapsto M_1b+\zeta(\Theta)N(\Theta),
 \qquad\mathbb E[\zeta\mid\Theta]=0,
\]

preserves (4.3), (4.4), and every level-set flux, while increasing
\(R_b\) by \(\mathbb E[\zeta^2\mid\Theta]\).  This modification is not
claimed to arise from a parallel flow; it proves that first-order endpoint
integration by parts and flux data alone cannot give the required upper
bound.

## 5. The genuine change-of-variables identity

For a fixed path \(B\), ordinary change of variables gives

\[
 \int\psi(\Phi_B(x))e^{\Lambda_0(x)}dx
 =\int\psi(\theta)
 {e^{\Lambda_0(\Phi_B^{-1}(\theta))}
  \over J_B(\Phi_B^{-1}(\theta))}d\theta.
 \tag{5.1}
\]

There is also an exact averaged identity.  Multiply (4.1) by
\(e^{\Lambda_0(x)}\), integrate in \(x\), and use

\[
 \int_{\mathbb R^n}e^{x\cdot\theta-|x|^2/2}dx
 =(2\pi)^{n/2}e^{|\theta|^2/2}.
\]

Together with (3.4), this gives

\[
 \boxed{
 \mathbb E_B\int
 \psi(\Phi_B(x))e^{\Lambda_0(x)}dx
 =\int\psi(\theta)e^{\Lambda_1(\theta)}d\theta.}
 \tag{5.2}
\]

Combining (5.1)--(5.2),

\[
 \boxed{
 \mathbb E_B\left[
 {e^{\Lambda_0(\Phi_B^{-1}(\theta))}
  \over J_B(\Phi_B^{-1}(\theta))}
 \right]
 =e^{\Lambda_1(\theta)}}
 \tag{5.3}
\]

for almost every \(\theta\).  These are sigma-finite identities; compactly
supported test functions justify every interchange.

Equation (5.3) is the exact determinant statement associated with the
random diffeomorphism.  It does not evaluate the derivative at \(x=0\).
The terminal point \(\Theta=\Phi_B(0)\) is correlated with \(B\); inserting
this random \(\theta\) into (5.3) is not a valid conditioning operation.
Moreover, (5.3) controls an inverse determinant, whereas (0.5) is a forward
normal singular component.

For completeness, if one integrates the alignment over initial tilts, then
pathwise

\[
\begin{aligned}
 &\int e^{\Lambda_0(x)}
 {\langle D\Phi_B(x)b,\nabla F(\Phi_B(x))\rangle^2
  \over|\nabla F(\Phi_B(x))|^2}dx\\
 &\quad=\int
 {e^{\Lambda_0(\Phi_B^{-1}(\theta))}
  \over J_B(\Phi_B^{-1}(\theta))}
 \langle D\Phi_B(\Phi_B^{-1}(\theta))b,N(\theta)\rangle^2d\theta.
\end{aligned}
 \tag{5.4}
\]

Even after averaging, (5.3) cannot separate the inverse-Jacobian weight
from the correlated normal stretch in (5.4).  More importantly, the desired
quantity is the integrand at the single initial point \(x=0\), not the
left side of (5.4).

## 6. Coarea at the posterior level sets

Let

\[
 \Sigma_s=\{\theta:F(\theta)=s\}
\]

at regular values.  Applying coarea to (3.8) yields

\[
 \boxed{
 \mathcal Q_b
 =\int_{\mathbb R}
 \left[
 \int_{\Sigma_s}{q(\theta)R_b(\theta)\over|\nabla F(\theta)|}
 d\mathcal H^{n-1}(\theta)
 \right]ds.}
 \tag{6.1}
\]

Write

\[
 E_b(s)=\int_{\Sigma_s}{qR_b\over|\nabla F|}d\mathcal H^{n-1},
 \qquad
 P_F(s)=\int_{\Sigma_s}q|\nabla F|d\mathcal H^{n-1}.
 \tag{6.2}
\]

The known posterior-gradient energy is

\[
 \int_{\mathbb R}P_F(s)ds
 =\int q|\nabla F|^2d\theta
 =\mathbb E|c_1|^2\le\lambda.
 \tag{6.3}
\]

The Stein identity (4.4) gives an exact flux through each level.  With the
normal \(N=\nabla F/|\nabla F|\) pointing out of \(\{F<s\}\),

\[
 \boxed{
 \int_{\Sigma_s}q\,\bar m_b\cdot N,d\mathcal H^{n-1}
 =-\int_{\{F<s\}}(b\cdot\theta)q(\theta)d\theta.}
 \tag{6.4}
\]

Call this common left side \(J_b(s)\).  Conditional Jensen and
Cauchy--Schwarz on \(\Sigma_s\) imply

\[
\begin{aligned}
 |J_b(s)|^2
 &\le
 \left(\int_{\Sigma_s}{qR_b\over|\nabla F|}d\mathcal H^{n-1}\right)
 \left(\int_{\Sigma_s}q|\nabla F|d\mathcal H^{n-1}\right)\\
 &=E_b(s)P_F(s).
\end{aligned}
 \tag{6.5}
\]

Therefore coarea proves

\[
 \boxed{
 \mathcal Q_b=\int E_b(s)ds
 \ge\int {J_b(s)^2\over P_F(s)}ds.}
 \tag{6.6}
\]

The direction is lower, not upper.  Level-set area or flux can force normal
energy, but it cannot prevent energy from oscillating on a small portion of
a level set.  The conditional-variance term in (3.9) is invisible to
\(J_b(s)\).

## 7. Pullback coarea and the tangential Jacobian

Fix a path and a point \(x\), put

\[
 y=\Phi_B(x),\qquad L=D\Phi_B(x),\qquad
 n=N(y).
\]

The source normal to the pullback level set of \(H_B=F\circ\Phi_B\) is

\[
 \nu_x={L^Tn\over|L^Tn|}.
 \tag{7.1}
\]

The tangential Jacobian of \(\Phi_B\) on that hypersurface is

\[
 J_{\mathrm{tan}}\Phi_B(x)
 =|\det L|\,|L^{-T}\nu_x|
 ={J_B(x)\over|L^Tn|}.
 \tag{7.2}
\]

Thus

\[
 \boxed{
 \langle Lb,n\rangle^2
 =|L^Tn|^2\langle b,\nu_x\rangle^2
 =\left({J_B(x)\over J_{\mathrm{tan}}\Phi_B(x)}\right)^2
 \langle b,\nu_x\rangle^2.}
 \tag{7.3}
\]

Formula (7.3) is the exact coarea geometry behind the alignment quotient.
A determinant bound is useful only together with a lower bound on the
tangential Jacobian of the relevant pullback level set.  Neither the
parallel-coupling trace theorem nor the posterior-gradient estimate supplies
such a bound.  Coarea itself merely records the compensation between normal
stretch and tangential area.

## 8. Why Gaussian integration by parts does not close the estimate

### 8.1 Endpoint integration by parts is first order

The density \(q\) satisfies

\[
 \nabla\log q(\theta)
 =\nabla\Lambda_1(\theta)-\theta
 =\mathbb E[X\mid X+G=\theta]-\theta.
 \tag{8.1}
\]

Ordinary integration by parts under \(q\) gives identities for the mean
field \(\bar m_b\), equivalent to (4.3)--(4.4).  It cannot bound the
conditional second moment \(R_b\) without a conditional Poincare or
variance estimate for the Brownian bridge given \(\Theta\).  No such
dimension-free estimate appears in the source.

### 8.2 Second differentiation creates an uncontrolled flow Hessian

Let

\[
 N_{1,b}=D_x^2\Phi_B(0)[b,b].
\]

Differentiate the transition identity twice at \(x=0\).  From (4.1)--(4.2),

\[
 \partial_b^2p_x(\theta)|_{x=0}
 =\big((b\cdot\theta)^2-2|b|^2\big)q(\theta).
\]

The flow side gives the exact identity

\[
\boxed{
\begin{aligned}
 \mathbb E[ D^2\psi(\Theta)[M_1b,M_1b]
 +\nabla\psi(\Theta)\cdot N_{1,b}]
 =\int\psi(\theta)
 \big((b\cdot\theta)^2-2|b|^2\big)q(\theta)d\theta.
\end{aligned}}
 \tag{8.2}
\]

Taking \(\psi=\Psi\circ F\) exposes the desired square,

\[
 D^2\psi[M_1b,M_1b]
 =\Psi''(F)\langle c_1,M_1b\rangle^2
 +\Psi'(F)D^2F[M_1b,M_1b],
 \tag{8.3}
\]

but (8.2) also contains

\[
 \Psi'(F)c_1\cdot N_{1,b}
 \quad\text{and}\quad
 \Psi'(F)D^2F[M_1b,M_1b].
\]

They have no sign.  Moreover, obtaining the weight
\(|c_1|^{-2}\) would require \(\Psi''(F)=|\nabla F|^{-2}\), which is not a
function of \(F\) unless the gradient magnitude is constant on every level.
Second-order endpoint integration by parts therefore does not isolate the
alignment quotient.

### 8.3 The initial tilt is not a finite-energy Brownian shift

The derivative \(M_1b\) changes the initial condition in (1.3).  A
Cameron--Martin perturbation \(h\) of Brownian motion must satisfy
\(h(0)=0\) and has energy \(\int_0^1|h'(t)|^2dt\).  Approximating an initial
jump \(b\) by

\[
 h_\varepsilon(t)=
 \begin{cases}(t/\varepsilon)b,&0\le t\le\varepsilon,\\
 b,&t>\varepsilon,
 \end{cases}
\]

costs

\[
 \int_0^1|h_\varepsilon'(t)|^2dt
 ={|b|^2\over\varepsilon}\longrightarrow\infty.
 \tag{8.4}
\]

Thus Wiener Gaussian Poincare or integration by parts cannot control
\(M_1b\) at universal cost.  The Malliavin derivative with respect to a
noise increment at time \(s>0\) is instead

\[
 D_s\Theta=M_1M_s^{-1},
 \tag{8.5}
\]

not \(M_1\).  Letting \(s\downarrow0\) recovers the singular initial
perturbation and the divergent Cameron--Martin cost.

## 9. Linear-map falsification of determinant and area bounds

Before using any determinant or coarea inequality, consider

\[
 L_K=\operatorname {diag}(K,1,\ldots,1),
 \qquad K\ge1,
 \tag{9.1}
\]

and

\[
 F(y)=y_1,\qquad b=e_1,\qquad n=e_1.
\]

This is an expansive map: \(L_K^TL_K\succeq I\), and it is the time-one
product integral of the constant positive-semidefinite matrix
\(\operatorname {diag}(\log K,0,\ldots,0)\).  Its normal alignment is

\[
 \langle L_Kb,n\rangle^2=K^2.
 \tag{9.2}
\]

Yet

\[
 \det L_K=K,
 \qquad
 \operatorname {Tr}(L_K^TL_K)=K^2+n-1.
 \tag{9.3}
\]

Choose \(K=\sqrt n\).  Then

\[
 {1\over n}\operatorname {Tr}(L_K^TL_K)<2,
 \qquad
 (\det L_K)^{1/n}=n^{1/(2n)}\longrightarrow1,
 \tag{9.4}
\]

while (9.2) equals \(n\).  Therefore neither trace per dimension nor
dimension-normalized volume expansion controls a selected normal stretch.

The target level sets \(\{F=s\}\) are parallel hyperplanes.  Their
tangential Jacobian under \(L_K\) is exactly one, while

\[
 {|\det L_K|\over J_{\mathrm{tan}}}=K.
\]

Thus even perfectly flat level sets with unchanged tangential area allow
arbitrarily large normal amplification.  The ordinary change-of-variables
weight contributes only \(1/\det L_K=1/K\); multiplying by the alignment
still leaves \(K\), which diverges.  A squared inverse determinant would
cancel this example, but no such weight occurs in (5.1)--(5.4).

This example does not claim to satisfy the stochastic transition identity
(4.3).  Its precise role is to disprove every proposed deduction based only
on Jacobian, trace, expansion, or coarea geometry.  The transition identity
adds a first-moment constraint, but (3.9) shows that first moments still do
not control the conditional normal variance.

## 10. Bottom line

The terminal quotient has a clean geometric meaning:

\[
 {\langle M_1b,c_1\rangle^2\over|c_1|^2}
 =\text{squared normal stretch of }D\Phi_B(0)b
 \text{ across the level set of }F.
\]

The exact objects are:

\[
 \Phi_B(x)=G_{1,B}(x),
 \qquad M_1=D\Phi_B(0),
 \qquad F(\theta)=\mathbb E_{\mu_{1,\theta}}f,
 \qquad c_1=\nabla F(\Phi_B(0)).
\]

The expectation is over Brownian paths with initial tilt fixed at zero.
Pathwise change of variables integrates over all initial tilts and introduces
the inverse Jacobian at a random preimage; it does not transform this fixed
point expectation.  Endpoint Gaussian integration by parts determines only
the conditional mean Stein field, while coarea turns its flux into a lower
bound for normal energy.  Determinant identities cannot prevent rank-one
normal expansion, as the linear example shows.

Accordingly, determinant/Jacobian control, Gaussian integration by parts,
and posterior level-set area do not prove a universal upper bound for the
terminal alignment.  A successful continuation would need a genuinely
second-moment statement for the Brownian bridge conditional on
\(\Theta\), or a direct estimate on the conditional variance term in
(3.9).  That is new information beyond the full parallel-coupling trace
theorem.
