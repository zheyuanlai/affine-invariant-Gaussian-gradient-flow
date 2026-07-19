# Spectral extremizers for KLS: exact variations, tensor-square rigidity, and the normal-cone obstruction

## 1. Scope and conclusion

For an isotropic log-concave probability measure \(\mu\) on \(\mathbb R^n\), write
\[
 C_P(\mu)=\sup_{\substack{f\in H^1(\mu)\\ \int f\,d\mu=0}}
 \frac{\int f^2\,d\mu}{\int |\nabla f|^2\,d\mu},
 \qquad \lambda_1(\mu)=C_P(\mu)^{-1}.
\]
The KLS conjecture asks for a numerical \(c>0\), independent of \(n\), such that
\(\lambda_1(\mu)\ge c\).

This report develops the direct extremal route in which both the log-concave law and its first eigenfunction are varied. The main positive result is an exact tensor-square rigidity theorem:

> **A bad extremizer cannot be a two-sided interior point of the convex-potential cone.** More precisely, if \(\mu\otimes\mu\) is two-sided stationary for its first spectral gap under all moment-preserving potential perturbations, then a first eigenfunction of \(\mu\) is affine. On full space this gives a standard Gaussian factor and \(\lambda_1(\mu)=1\); on a bounded Neumann domain it is impossible.

The proof is dimension-free and uses the off-diagonal entry of the full eigenspace stress matrix. It also has a sharp Hilbert-space stability version.

The obstruction is equally exact. Convexity is an inequality constraint. At an extremizer its KKT multiplier is a positive-semidefinite matrix-valued measure \(M\), and the stationarity equation is
\[
 \bigl(G_\Theta+b\cdot x+A:(xx^T-I)\bigr)\mu=-\operatorname{div}\operatorname{div}M,
 \qquad M:D^2V=0.
\]
Here \(G_\Theta\) is a positive-semidefinite mixture of first-eigenspace stresses. Flat potentials, including all uniform measures on convex bodies, allow arbitrary nonzero \(M\). Facets have an analogous shape multiplier. The multiplier can carry precisely the off-diagonal tensor stress needed by the rigidity argument. A vanishing regularization coefficient does not make this defect vanish.

Thus the variational route identifies a concrete dimension-free mechanism, but it does **not** prove KLS. The remaining gap is a dimension-free estimate that controls or eliminates the convexity/shape normal-cone defects (at the scale of the normalized stress), together with a quantitative spectral stability statement. Sections 10--12 make this gap precise and test every tempting shortcut on Gaussian, cube, simplex, and radial examples.

## 2. Smooth weighted Neumann setting

Let \(\Omega\subset\mathbb R^n\) be a bounded, connected, convex \(C^3\) domain and let \(V\in C^3(\overline\Omega)\) be convex. Put
\[
 d\mu=\rho(x)\,dx,\qquad \rho=Z^{-1}e^{-V},\qquad
 Z=\int_\Omega e^{-V}\,dx.
\]
The positive weighted Neumann operator is
\[
 \mathcal A_V=-\Delta+\nabla V\cdot\nabla
             =-\rho^{-1}\operatorname{div}(\rho\nabla),
 \qquad \partial_n f|_{\partial\Omega}=0.
\]
Its spectrum is discrete,
\(0=\lambda_0<\lambda_1\le\lambda_2\le\cdots\), and
\[
 \lambda_1=\inf_{\int f\,d\mu=0}
 \frac{\int |\nabla f|^2\,d\mu}{\int f^2\,d\mu}.
\]
All formulas below remain valid on \(\mathbb R^n\) when the displayed functions and their derivatives have enough integrability to justify cutoff limits. Smooth strongly convex potentials are more than sufficient.

Assume throughout a displayed stationarity formula that \(\mu\) is isotropic:
\[
 \int x\,d\mu=0,\qquad \int xx^T\,d\mu=I_n.
\]
Let the multiplicity of \(\lambda=\lambda_1\) be \(m\), and choose an orthonormal eigenbasis
\[
 \mathcal A_V f_a=\lambda f_a,\qquad
 \int f_af_b\,d\mu=\delta_{ab},\qquad \int f_a\,d\mu=0,
 \quad 1\le a,b\le m.
\]
Define the stress matrix
\[
 g_{ab}(x)=\nabla f_a(x)\cdot\nabla f_b(x)-\lambda f_a(x)f_b(x).
\]
It is centered:
\[
 \int g_{ab}\,d\mu=0.
\]
For a simple eigenvalue we write \(f=f_1\) and
\(g=|\nabla f|^2-\lambda f^2\).

Two elementary identities will be used repeatedly. In the weak Neumann sense,
\[
 \mathcal A_V(f^2)=2\lambda f^2-2|\nabla f|^2=-2g.                 \tag{2.1}
\]
Also, testing the eigenvalue equation against \(f^{k+1}\), whenever integrable, gives
\[
 \int f^k|\nabla f|^2\,d\mu
   =\frac{\lambda}{k+1}\int f^{k+2}\,d\mu,
\]
and hence
\[
 \int g f^k\,d\mu
   =-\frac{k}{k+1}\lambda\int f^{k+2}\,d\mu.                  \tag{2.2}
\]
In particular,
\(\int gf^2\,d\mu=-(2\lambda/3)\int f^4\,d\mu\). These identities alone do not give a lower bound for \(\lambda\).

Finally, isotropy always gives the trial-function bound
\[
 \lambda_1(\mu)\le 1,                                        \tag{2.3}
\]
because every \(u\cdot x\), \(|u|=1\), has variance and Dirichlet energy both equal to one.

### 2.1 Euler equation in the function variable

For fixed \((\Omega,V)\), the joint constrained functional is
\[
 \mathscr L(f)=\int|\nabla f|^2\,d\mu
 -\lambda\left(\int f^2\,d\mu-1\right)-2\alpha\int f\,d\mu.
\]
Its first variation in an arbitrary \(\phi\in H^1(\mu)\) is
\[
 \int\nabla f\cdot\nabla\phi\,d\mu
 =\lambda\int f\phi\,d\mu+\alpha\int\phi\,d\mu.             \tag{2.4}
\]
Taking \(\phi=1\) gives \(\alpha=0\). Integration by parts then gives
\[
 \mathcal A_Vf=\lambda f\quad\text{in }\Omega,
 \qquad \partial_nf=0\quad\text{on }\partial\Omega.          \tag{2.5}
\]
Thus the eigenvalue equation and the Neumann condition are exactly the Euler equations in \(f\); they are not extra assumptions. At a multiple first eigenvalue, any normalized vector in the first eigenspace is a function-variable minimizer, while measure variation must use the branch matrix of Section 3.

The integrated weighted Bochner--Reilly identity provides a useful audit of signs. For a normalized eigenfunction,
\[
 \lambda^2
 =\int_\Omega\|D^2f\|_{\mathrm{HS}}^2\,d\mu
  +\int_\Omega D^2V(\nabla f,\nabla f)\,d\mu
  +\frac1Z\int_{\partial\Omega}e^{-V}
     \mathrm{II}(\nabla_\tau f,\nabla_\tau f)\,dS.           \tag{2.6}
\]
Here \(\mathrm{II}\succeq0\) for the outward-normal convention on a convex boundary. Formula (2.6) is dimension-free, but without a positive lower curvature bound it gives no lower bound on \(\lambda\); all three terms may degenerate along an extremizing sequence.

## 3. Exact potential variation

Let \(W\in C^2(\overline\Omega)\) and \(V_t=V+tW\), with the domain fixed. Then
\[
 \left.\frac d{dt}\right|_{0}d\mu_t
   =\bigl(\mathbb E_\mu W-W\bigr)d\mu.                       \tag{3.1}
\]

### Lemma 3.1 (eigenvalue-branch derivative)

The first derivatives of the \(m\) eigenvalue branches issuing from \(\lambda\) are the eigenvalues of the symmetric matrix
\[
 H(W)_{ab}=-\int W g_{ab}\,d\mu.                              \tag{3.2}
\]
Consequently, the right directional derivative of the bottom eigenvalue in the cluster is
\[
 D^+\lambda_1(V)[W]=\lambda_{\min}(H(W)).                     \tag{3.3}
\]

#### Proof

Use the generalized quadratic forms
\(a_t(u,v)=\int\nabla u\cdot\nabla v\,d\mu_t\) and
\(b_t(u,v)=\int uv\,d\mu_t\). On the eigenspace,
standard finite-dimensional degenerate perturbation gives the branch matrix
\(a'_0-\lambda b'_0\). Formula (3.1) yields
\[
 a'_0(f_a,f_b)-\lambda b'_0(f_a,f_b)
 =-\int W(\nabla f_a\cdot\nabla f_b-\lambda f_af_b)\,d\mu;
\]
the terms containing \(\mathbb EW\) cancel because both unperturbed forms equal
\(\lambda\delta_{ab}\) and \(\delta_{ab}\), respectively. The lowest branch has derivative the smallest eigenvalue of this matrix. \(\square\)

At an isotropic point, the first variations of the moment constraints are
\[
 \dot m=-\int W x\,d\mu,\qquad
 \dot\Sigma=-\int W(xx^T-I)\,d\mu.                            \tag{3.4}
\]
Thus the linear tangent space of moment-preserving potential perturbations is
\[
 \mathcal T=\left\{W:
   \int Wx\,d\mu=0,\quad
   \int W(xx^T-I)\,d\mu=0\right\}.                            \tag{3.5}
\]
Adding a constant to \(W\) changes neither the law nor these conditions.

Because \(\mu\) has positive density on an open set, the degree-one and centered degree-two statistics
\[
 x_i,\qquad x_ix_j-\delta_{ij}\quad (1\le i\le j\le n)
\]
are linearly independent in \(L^2(\mu)\). Their Gram matrix is therefore positive definite. It follows that (3.4) is onto the moment-constraint space and that every \(W\in\mathcal T\) is tangent to a genuine two-sided, exactly isotropic \(C^1\) path after adding a small degree-two correction, by the implicit function theorem.

### Proposition 3.2 (full matrix Euler identity at a two-sided interior minimum)

Suppose \(V\) is a local minimizer of \(\lambda_1\) among isotropic potentials and that every \(W\in\mathcal T\) is allowed with both signs; for example, it is enough that
\(D^2V\succeq\kappa I\) and \(W\) have bounded Hessian. Then for every \(a,b\),
\[
 g_{ab}(x)=b_{ab}\cdot x+A_{ab}:(xx^T-I)\quad\text{in }L^2(\mu), \tag{3.6}
\]
for some \(b_{ab}\in\mathbb R^n\) and symmetric \(A_{ab}\in\mathbb R^{n\times n}\).

#### Proof

Local minimality and (3.3) give
\(\lambda_{\min}(H(W))\ge0\) for every \(W\in\mathcal T\). Since \(-W\in\mathcal T\),
\[
 \lambda_{\min}(-H(W))=-\lambda_{\max}(H(W))\ge0.
\]
Hence \(H(W)=0\) for every tangent \(W\). Each \(g_{ab}\) therefore annihilates the common kernel (3.5). The annihilator of that kernel is exactly the span of the moment statistics, which gives (3.6). \(\square\)

The use of the *entire matrix* in this proposition is essential. A scalar nonsmooth KKT argument would generally supply only one mixture of stresses.

For a simple eigenvalue, (3.6) becomes
\[
 |\nabla f|^2-\lambda f^2=b\cdot x+A:(xx^T-I).                 \tag{3.7}
\]
It is tempting, but unjustified, to infer from (3.7) that \(f\) is affine or that \(A\) has rank one. Section 11 records exact counterexamples to the rank claim and shows where nonquadratic stresses occur.

## 4. Convexity KKT and the defect measure

The preceding proposition applies only at a two-sided interior point of the convex-potential cone. The correct identity at the boundary contains a matrix-valued multiplier.

Work in the Banach space \(C^2(\overline\Omega)\), and view convexity as the conic constraint
\[
 D^2V\in C(\overline\Omega;\mathbb S^n_+).
\]
The dual cone consists of finite positive-semidefinite matrix-valued Radon measures. The Hessian constraint by itself has the strict Slater direction \(W(x)=|x|^2\). For the Hessian and moment constraints together, assume the following combined Robinson condition: the derivative of the moment map is onto, and its kernel contains a direction \(W_0\) for which
\[
 D^2V+D^2W_0\succ0\quad\text{on }\overline\Omega.             \tag{4.0}
\]
This is the exact constraint qualification used below. In a finite-dimensional Galerkin restriction it is directly checkable; without it, one must allow additional abnormal multipliers. Weak-* limits of the resulting KKT multipliers give the continuum identity whenever their total variations stay bounded.

For an eigenvalue of multiplicity \(m\), let \(\Theta\in\mathbb S^m_+\) with
\(\operatorname{tr}\Theta=1\), and define
\[
 G_\Theta=\sum_{a,b=1}^m\Theta_{ab}g_{ab}.                    \tag{4.1}
\]

### Proposition 4.1 (potential KKT identity)

At a constrained local minimum of \(\lambda_1\), under the preceding constraint qualifications, there are
\[
 \Theta\succeq0,\quad\operatorname{tr}\Theta=1,\quad
 b\in\mathbb R^n,\quad A\in\mathbb S^n,\quad
 M\succeq0
\]
with \(M\) a finite matrix-valued measure on \(\overline\Omega\), such that for every \(W\in C^2(\overline\Omega)\),
\[
 \int W\bigl(G_\Theta+b\cdot x+A:(xx^T-I)\bigr)\,d\mu
   =-\int_{\overline\Omega}D^2W:M(dx),                        \tag{4.2}
\]
and
\[
 \int_{\overline\Omega}D^2V:M(dx)=0.                         \tag{4.3}
\]
In distribution notation, including the natural boundary action implicit in (4.2),
\[
 \bigl(G_\Theta+b\cdot x+A:(xx^T-I)\bigr)\mu
   =-\operatorname{div}\operatorname{div}M.                  \tag{4.4}
\]

#### Proof

The Clarke differential of the lowest eigenvalue cluster is the convex hull of the branch differentials (3.2), hence is represented by some density matrix \(\Theta\succeq0\), \(\operatorname{tr}\Theta=1\). Add Lagrange multipliers for (3.4). The derivative of the resulting scalar Lagrangian in direction \(W\) is
\[
 -\int W\bigl(G_\Theta+b\cdot x+A:(xx^T-I)\bigr)\,d\mu.
\]
The normal cone to \(D^2V\succeq0\) is the adjoint Hessian of a positive-semidefinite measure, with complementary slackness. Choosing the inequality Lagrangian with term
\(-\int D^2V:M\) gives (4.2)--(4.3). \(\square\)

If \(D^2V\succeq\kappa I\) on \(\overline\Omega\), (4.3) implies
\[
 0=\int D^2V:M\ge \kappa\int\operatorname{tr}M,
\]
so \(M=0\). Thus the scalar mixture in (4.1) is quadratic. At a genuine two-sided minimum Proposition 3.2 is stronger and makes every matrix entry quadratic.

If \(V\) is flat on an open set, (4.3) gives no control on \(M\) there. This is not a technical pathology. In one dimension the polar cone can be written explicitly. If a signed measure \(\sigma\) annihilates affine functions and
\(\int\phi\,d\sigma\le0\) for every convex \(\phi\), then
\[
 \sigma=-m''
\]
in distributions for a nonnegative function or measure \(m\). Indeed,
\(k(t)=\int(x-t)_+\,\sigma(dx)\le0\), \(k''=\sigma\), and one takes \(m=-k\). Complementarity is \(mV''=0\). On an interval where \(V''=0\), any nonnegative smooth compactly supported \(m\) produces a nonzero absolutely continuous defect \(-m''\). Hence flatness permits an infinite-dimensional family of stresses; it does not force a one-dimensional or rank-one stress.

## 5. Affine variations are a gauge

Let \(T_t(x)=x+t(a+Hx)\), and set \(\mu_t=(T_t)_\#\mu\). Pull a competitor on the target back to the source. The measure in the pulled-back Rayleigh quotient is fixed, while the inverse metric is
\((DT_t)^{-1}(DT_t)^{-T}\). Therefore the eigenvalue-branch matrix is
\[
 K(H)_{ab}=-2\int \nabla f_a^T(\operatorname{sym}H)\nabla f_b\,d\mu. \tag{5.1}
\]
The moment derivatives at isotropy are
\[
 \dot m=a,
 \qquad \dot\Sigma=H+H^T.                                    \tag{5.2}
\]
Thus a first-order affine variation tangent to the isotropic class has
\(a=0\) and \(H^T=-H\), and (5.1) vanishes identically.

In fact the degeneracy is exact, not merely infinitesimal. If \(X\) is isotropic and \(Y=TX+a\), then centering and isotropizing \(Y\) gives
\[
 (TT^T)^{-1/2}(Y-a)=QX,
 \qquad Q=(TT^T)^{-1/2}T,\qquad QQ^T=I.
\]
So every invertible affine deformation, followed by exact isotropic normalization, is just an orthogonal image of the original measure. It leaves \(\lambda_1\) unchanged. Affine variations can determine Lagrange multipliers, but cannot supply a coercive spectral identity.

## 6. Exact boundary-shape variation

Let \(V\) be defined in a neighborhood of \(\overline\Omega\). Move the boundary with normal velocity \(h\in C^2(\partial\Omega)\), keeping the ambient potential fixed. The weighted Neumann Hadamard formula gives the branch matrix
\[
 H^{\partial}(h)_{ab}
 =\frac1Z\int_{\partial\Omega}h e^{-V}
   \bigl(\nabla_\tau f_a\cdot\nabla_\tau f_b-\lambda f_af_b\bigr)\,dS. \tag{6.1}
\]
The tangential gradients appear because \(\partial_nf_a=0\). The normalization of the probability causes no additional term: the Rayleigh quotient is the quotient of two unnormalized integrals, and their common mass factor cancels.

The moment derivatives are
\[
 \dot m=\frac1Z\int_{\partial\Omega}he^{-V}x\,dS,
 \qquad
 \dot\Sigma=\frac1Z\int_{\partial\Omega}he^{-V}(xx^T-I)\,dS. \tag{6.2}
\]
Equations (6.1)--(6.2) follow directly from Reynolds transport and the envelope theorem for the Rayleigh quotient.

For a smooth strictly convex body containing the origin, write its support function as \(s\) on \(S^{n-1}\), with curvature-radius matrix
\[
 Q_s=\nabla^2_{S^{n-1}}s+sI\succeq0.                          \tag{6.3}
\]
A support perturbation \(s+t\eta\) has normal velocity \(\eta\) at the point with outer normal \(\theta\). If \(x_s(\theta)=\nabla_Ss+s\theta\), then the surface Jacobian is \(\det Q_s\).

The shape analogue of Proposition 4.1 is the following. Let
\[
 S_\partial(x)=G_\Theta^\partial(x)+b\cdot x+A:(xx^T-I),
\]
where
\[
 G_\Theta^\partial
 =\sum_{a,b}\Theta_{ab}
   (\nabla_\tau f_a\cdot\nabla_\tau f_b-\lambda f_af_b).
\]
Under the support-function Slater and moment-constraint qualifications, there is a positive-semidefinite tangential matrix measure \(N\) on \(S^{n-1}\) such that
\[
 \frac1Z\int_{S^{n-1}}\eta(\theta)e^{-V(x_s(\theta))}
 S_\partial(x_s(\theta))\det Q_s(\theta)\,d\theta
 =\int_{S^{n-1}}(\nabla_S^2\eta+\eta I):N(d\theta),           \tag{6.4}
\]
and
\[
 \int_{S^{n-1}}Q_s:N=0.                                     \tag{6.5}
\]
This is the conic KKT identity for (6.3). If \(Q_s\succeq\kappa_\partial I\), then \(N=0\), and (6.4) forces
\[
 G_\Theta^\partial+b\cdot x+A:(xx^T-I)=0
 \quad\text{pointwise on }\partial\Omega.                   \tag{6.6}
\]
At a two-sided local minimum, with the natural boundary moment map of full rank, the same \(H(h)\), \(H(-h)\) argument as in Proposition 3.2 yields (6.6) separately for every \(g_{ab}^\partial\), not only for one mixture.

For a polytope, arbitrary two-sided normal velocities are not available. On the relative interior of a flat facet, put the body locally in the form
\(\{(y,z):z\le0\}\). A perturbed boundary \(z=t h(y)\) bounds a convex hypograph for \(t>0\) only if \(h\) is concave. Requiring the same for \(-t\) forces \(h\) to be affine. Thus pointwise boundary stationarity cannot be inferred on facets; at most one obtains finitely many facet-average identities from parallel facet shifts.

## 7. A rigorous regularized compact extremizer

There are two separate roles for regularization: obtaining an actual smooth minimizer in a compact class, and approximating arbitrary log-concave laws. They must not be conflated.

### 7.1 Finite-dimensional compact regularization

Fix a finite-dimensional \(C^2\) family \(z\mapsto(\Omega_z,V_z)\) of smooth convex pairs. For example, take finitely many smooth basis functions for \(V\) and finitely many spherical harmonics for the support function, impose closed coefficient bounds, and impose
\[
 0\preceq D^2V_z\preceq KI,\qquad
 0\preceq Q_{s_z}\preceq K_\partial I,\qquad
 rB_2^n\subset\Omega_z\subset RB_2^n.                         \tag{7.1}
\]
Let \(K_N\) be the resulting compact parameter set, and assume the isotropy equations cut out a regular submanifold in its interior. Let \(B_N\) be any \(C^1\) proper barrier on that interior, tending to \(+\infty\) at the boundary. For \(\eta>0\), minimize
\[
 J_{N,\eta}(z)=\log\lambda_1(\Omega_z,V_z)+\eta B_N(z)         \tag{7.2}
\]
subject to isotropy.

Equivalently, this is the joint program
\[
 \min_{z,f}\left\{
 \log\!\int_{\Omega_z}|\nabla f|^2\,d\mu_z+\eta B_N(z):
 \int f\,d\mu_z=0,\ \int f^2\,d\mu_z=1,\ m_z=0,\ \Sigma_z=I
 \right\}.                                                    \tag{7.2a}
\]
After pulling functions to a fixed reference domain, compact Sobolev embedding gives a minimizing \(f\). Its Euler equation is exactly (2.4)--(2.5), so the measure, shape, and function variables are all varied in this regularized problem.

Because the weighted Neumann forms vary continuously under \(C^2\) domain diffeomorphisms and \(C^1\) convergence of \(V\), \(\lambda_1\) is continuous on this family. The barrier therefore gives an interior minimizer \(z_{N,\eta}\). At a simple eigenvalue the ordinary KKT equation is
\[
 D\log\lambda_1+\eta DB_N+b\cdot Dm+A:D(\Sigma-I)=0.         \tag{7.3}
\]
At a multiple eigenvalue, \(D\log\lambda_1\) is represented by a density matrix \(\Theta\) as in (4.1), divided by \(\lambda\). Equations (3.2), (5.1), and (6.1) give every term in (7.3) exactly.

Increasing the basis dimension makes such families dense in any fixed class with strict inequalities in (7.1). Thus (7.2) is a rigorous regularized compact smooth extremizer problem. It does not, however, justify deleting \(\eta DB_N\) as \(\eta\downarrow0\).

The elementary scalar model
\[
 \min_{0<s\le1}\{s-\eta\log s\}
\]
has minimizer \(s=\eta\), but
\(\eta(-1/s)=-1\) there. The barrier value \(\eta|\log\eta|\) tends to zero while its force remains of order one. A log-determinant curvature barrier has exactly this behavior: if an eigenvalue of \(D^2V_\eta\) is of order \(\eta\), then
\(\eta(D^2V_\eta)^{-1}\) has an order-one weak limit. That limit is the measure \(M\) in (4.2). The same phenomenon produces \(N\) when boundary curvature radii degenerate.

### 7.2 Compactness without a barrier

For completeness, fix \(k\ge3\), \(0<r<R\), positive upper regularity bounds, and positive lower and upper bounds on the Hessian of \(V\) and on \(Q_s\). The isotropic subclass of pairs obeying these closed bounds is compact in a weaker \(C^{k-1}\) topology by Arzela--Ascoli. Pulling all domains to a reference domain gives Mosco convergence of the weighted Neumann forms, so \(\lambda_1\) is continuous and attains its minimum. Such a minimizer can lie on any of the artificial bounds. The corresponding KKT multipliers must be retained. Letting the bounds recede again leads to (4.2) and (6.4), not automatically to the polynomial identity (3.6).

## 8. Approximation of arbitrary isotropic log-concave laws

The following reduction is sufficient for KLS and avoids any assertion that spectral gaps converge in the difficult direction.

### Lemma 8.1 (smooth, strongly convex, compact approximation)

Let \(\mu\) be an isotropic log-concave probability on \(\mathbb R^n\). There is a sequence \(\mu_j\) of isotropic probabilities with the following properties:

1. \(\mu_j\) has a \(C^\infty\), strictly positive density on a bounded \(C^\infty\), strictly convex domain \(\Omega_j\);
2. its interior potential is \(C^\infty\) and strongly convex, with a positive lower Hessian bound depending on \(j\);
3. \(\mu_j\Rightarrow\mu\) and all moments of every fixed order converge.

Consequently, if \(C_P(\mu_j)\le C\) with the same numerical \(C\) for all such smooth pairs, then \(C_P(\mu)\le C\).

#### Proof

Write the log-concave density of \(\mu\) as \(p=e^{-V}\), allowing \(V=+\infty\) off its convex support. Convolve with a centered Gaussian of covariance \(\delta I\):
\(p_\delta=p*\gamma_\delta\). Prekopa's theorem gives log-concavity, and Gaussian convolution gives a positive \(C^\infty\) density. Then set
\[
 \widetilde p_{\delta,\varepsilon}(x)
 =c_{\delta,\varepsilon}p_\delta(x)e^{-\varepsilon|x|^2/2}.
\]
Its potential has Hessian at least \(\varepsilon I\). Restrict it to a ball \(B_R\), normalize, and finally center and isotropize. The ball becomes an ellipsoid, which is smooth and strictly convex, and strong convexity is preserved under the invertible affine change.

Take first \(R\to\infty\), then \(\delta,\varepsilon\downarrow0\), and diagonalize. Log-concave exponential tails and the fixed second moments give uniform integrability of every fixed moment after isotropic normalization, so the claimed convergence follows.

For \(h\in C_c^\infty\), weak convergence gives
\[
 \operatorname{Var}_{\mu_j}h\to\operatorname{Var}_\mu h,
 \qquad
 \int|\nabla h|^2\,d\mu_j\to\int|\nabla h|^2\,d\mu.
\]
Passing the uniform Poincare inequality to the limit proves it for such \(h\). Cutoff and mollification are dense in the weighted Sobolev space on the interior of a convex support, so the inequality extends to every \(H^1(\mu)\) function. \(\square\)

This lemma shows that restricting a proof to smooth strongly convex compact pairs is legitimate only if every constant is uniform as the lower curvature and boundary-curvature bounds tend to zero and the diameter tends to infinity. Precisely those limits create the defect measures \(M\) and \(N\).

## 9. Tensor-square rigidity: the dimension-free mechanism

The strongest consequence of the full matrix identity is obtained by tensoring an extremizer with itself.

### Lemma 9.1 (degree-two tensor factorization)

Let \(\mu\) and \(\nu\) be centered, full-dimensional probabilities with finite fourth moments, and let
\(f\in L^2_0(\mu)\), \(h\in L^2_0(\nu)\) be nonzero. If
\[
 f(x)h(y)=P(x,y)\quad \mu\otimes\nu\text{-a.e.}              \tag{9.1}
\]
for a polynomial \(P\) of total degree at most two, then \(f\) and \(h\) are affine; centeredness makes them linear.

More quantitatively, suppose \(\mu,\nu\) are isotropic, \(\|f\|_2=\|h\|_2=1\), and
\[
 \operatorname{dist}_{L^2(\mu\otimes\nu)}
   (f\otimes h,\mathcal P_2)\le\varepsilon\le1,              \tag{9.2}
\]
where \(\mathcal P_2\) is the space of total-degree-at-most-two polynomials. If \(U_\mu=\operatorname{span}\{x_1,\ldots,x_n\}\) and similarly for \(U_\nu\), then
\[
 \|P_{U_\mu}f\|_2\,\|P_{U_\nu}h\|_2\ge\sqrt{1-\varepsilon^2}. \tag{9.3}
\]
In the symmetric case \(\mu=\nu\), \(f=h\),
\[
 \operatorname{dist}_{L^2(\mu)}(f,U_\mu)^2
 \le 1-\sqrt{1-\varepsilon^2}\le\varepsilon^2.               \tag{9.4}
\]

#### Proof

Apply double centering \((I-\mathbb E_x)(I-\mathbb E_y)\), an orthogonal contraction in \(L^2(\mu\otimes\nu)\). It fixes \(f\otimes h\). The double-centered part of a total-degree-two polynomial lies in
\(U_\mu\otimes U_\nu\), because only bilinear terms survive. Hence (9.2) implies
\[
 \operatorname{dist}(f\otimes h,U_\mu\otimes U_\nu)\le\varepsilon.
\]
The orthogonal projection is
\((P_{U_\mu}f)\otimes(P_{U_\nu}h)\), whose squared norm is the product of the two squared projection norms. This proves (9.3), and the symmetric specialization gives (9.4). With \(\varepsilon=0\), both functions lie in their linear-coordinate spaces, proving the exact assertion. \(\square\)

### Lemma 9.2 (an affine first eigenfunction gives a Gaussian factor)

Let \(d\mu\propto e^{-V}dx\) be isotropic on \(\mathbb R^n\), with \(V\in C^1\) convex. If a normalized first eigenfunction is
\(f(x)=u\cdot x\), \(|u|=1\), then
\[
 \lambda_1=1,
 \qquad
 \mu=\gamma_1\otimes\nu
\]
after rotating \(u\) to the first coordinate, where \(\gamma_1\) is the standard one-dimensional Gaussian and \(\nu\) is isotropic log-concave on \(u^\perp\).

On a bounded domain with Neumann boundary condition, no nonconstant affine function can be an eigenfunction.

#### Proof

The eigenvalue equation gives
\[
 u\cdot\nabla V=\lambda u\cdot x.
\]
Writing \(x=su+z\) and integrating in \(s\),
\[
 V(su+z)=\frac\lambda2s^2+W(z).
\]
Thus the \(u\)-marginal is Gaussian with variance \(1/\lambda\). Isotropy makes that variance one, so \(\lambda=1\) and the factorization follows.

On a bounded domain, Neumann data for \(u\cdot x\) would require \(u\cdot n=0\) almost everywhere on the boundary. But the divergence theorem applied to \((u\cdot x)u\) gives
\[
 |u|^2|\Omega|=\int_{\partial\Omega}(u\cdot x)(u\cdot n)\,dS=0,
\]
a contradiction. \(\square\)

### Theorem 9.3 (tensor-square interior rigidity)

Let \(\mu\) be an isotropic smooth log-concave probability, and let \(f\) be a normalized first eigenfunction with eigenvalue \(\lambda\). Suppose \(\mu\otimes\mu\) is a local minimum of \(\lambda_1\) among isotropic log-concave probabilities and its product potential is a two-sided interior point for all smooth moment-preserving potential perturbations. Then:

* on \(\mathbb R^n\), \(\lambda=1\) and \(\mu\) has a standard Gaussian factor;
* on a bounded Neumann domain, the hypotheses are inconsistent.

#### Proof

The product first eigenspace contains
\[
 F_1(x,y)=f(x),\qquad F_2(x,y)=f(y).
\]
Their gradients lie in orthogonal coordinate blocks, so their off-diagonal stress is
\[
 g_{12}^{\otimes}(x,y)
 =\nabla F_1\cdot\nabla F_2-\lambda F_1F_2
 =-\lambda f(x)f(y).                                         \tag{9.5}
\]
Proposition 3.2 says that every stress-matrix entry is a total-degree-two polynomial in \((x,y)\). Since \(\lambda>0\), Lemma 9.1 makes \(f\) affine. Lemma 9.2 finishes the proof. \(\square\)

This is a genuinely dimension-free extremal mechanism. It is stronger than trying to show directly that the scalar stress (3.7) has rank one. It exploits eigenvalue multiplicity created by tensorization and the mixed block, where the gradient term vanishes identically.

There is also a natural near-extremal scheme. Let
\[
 a_n=\inf\{\lambda_1(\mu):\mu\text{ isotropic log-concave on }\mathbb R^n\},
 \qquad a_*=\inf_n a_n.
\]
Tensorization gives
\(a_{n+m}\le\min(a_n,a_m)\). If \(a_n\le a_*+\delta\), then the tensor square of an almost minimizer in dimension \(n\) is a \(\delta\)-almost minimizer in dimension \(2n\). If one could obtain relative two-sided stationarity for
\(\log\lambda_1\), then (9.5) divided by \(\lambda\) is exactly
\(-f\otimes f\), with no small-gap loss. Lemma 9.1 would give dimension-free closeness of \(f\) to a linear function.

The phrase "relative two-sided stationarity" is essential. For \(\lambda\) itself, the mixed stress has size \(\lambda\), so an absolute error is useless when \(\lambda\) is small. For \(\log\lambda\), the normalized mixed stress has size one. The unresolved issue is that the normal-cone terms in Sections 4 and 6 are also of order one after this normalization.

There is an additional nonsmooth obstruction. Proposition 4.1 supplies only one density matrix \(\Theta\). At a tensor product it may choose a block-diagonal \(\Theta\), in which case the scalar KKT equation does not display \(g_{12}^{\otimes}\) at all. Recovering the mixed entry requires approximate *full-matrix* stationarity, as in Proposition 3.2, which in turn requires quantitative access to both signs of a perturbation.

## 10. What a successful closure would have to prove

The tensor argument reduces the variational route to two dimension-free statements. Neither is established here.

### Target A: normalized normal-cone control

For a tensor-square near-minimizer of \(\log\lambda_1\), one would first need to upgrade the one-sided scalar KKT condition to approximate two-sided full-matrix stationarity. In particular, the mixed-block normalized stress would have to satisfy
\[
 \operatorname{dist}_{L^2(\mu\otimes\mu)}
 \left(f\otimes f,\mathcal P_2\right)\le\varepsilon          \tag{10.1}
\]
with a numerical \(\varepsilon<1\), ideally \(o(1)\), after accounting for the potential multiplier \(M\) and the shape multiplier \(N\). This requires both (i) preventing the scalar spectral multiplier \(\Theta\) from discarding the mixed eigenspace entry and (ii) controlling the mixed block of
\(\operatorname{div}\operatorname{div}M\) in a norm strong enough to imply (10.1). Positivity of \(M\) alone only gives a block Cauchy--Schwarz inequality; complementarity gives no diagonal control where \(D^2V=0\). Uniform convex bodies have \(D^2V=0\) everywhere in their interiors.

### Target B: spectral stability near a linear eigenfunction

Lemma 9.1 turns (10.1) into \(L^2\)-closeness of \(f\) to \(u\cdot x\), with no dimension dependence. Exact linearity implies \(\lambda=1\), but \(L^2\)-closeness alone does not control the gradient of the residual without already using a Poincare estimate. A closure needs a statement of the form
\[
 \operatorname{dist}_{L^2}(f,U)\le\varepsilon_0
 \quad\text{and the normalized KKT equations}
 \quad\Longrightarrow\quad \lambda\ge c_0,                  \tag{10.2}
\]
with numerical \(\varepsilon_0,c_0>0\). The KKT equations, not mere \(L^2\) closeness, must supply the missing gradient control.

A proof of Targets A and B would yield KLS: choose dimensions approaching \(a_*\), tensor-square, regularize by a tensor-stable Ekeland procedure for \(\log\lambda_1\), apply (10.1)--(10.2), and then use Lemma 8.1. The current gap is Target A already at exact flat extremizers, and Target B at the quantitative level.

## 11. Countermodels and stress tests

### 11.1 Gaussian: quadratic stress does not select rank one

For \(\gamma_n=N(0,I_n)\), \(\lambda_1=1\) and the first eigenspace is
\(f_u(x)=u\cdot x\), \(|u|=1\). For an orthonormal coordinate basis,
\[
 g_{ij}=\delta_{ij}-x_ix_j.                                  \tag{11.1}
\]
Thus every entry satisfies the exact quadratic Euler identity. If
\(\Theta=I_n/n\), then
\[
 G_\Theta=1-\frac{|x|^2}{n},                                 \tag{11.2}
\]
whose quadratic coefficient has rank \(n\). Hence positivity of \(\Theta\), rotational symmetry, and exact stationarity do not force a rank-one multiplier or a distinguished one-dimensional direction. The full matrix identity, rather than a scalar mixture, is what makes the tensor-square argument work.

The Gaussian also shows that radial symmetry is compatible with a maximally high-rank scalar stress. It is not a counterexample to the affine conclusion of Theorem 9.3; it is the equality model.

### 11.2 Isotropic cube: product structure does not make the stress quadratic

Let
\[
 \Omega=[-\sqrt3,\sqrt3]^n,
 \qquad \mu=\text{uniform probability on }\Omega.
\]
Each coordinate has variance one. Put
\[
 a=\sqrt3,\qquad k=\frac\pi{2a},\qquad
 \lambda=k^2=\frac{\pi^2}{12},\qquad
 f_i(x)=\sqrt2\sin(kx_i).
\]
The product Neumann spectrum shows that \(\lambda_1=\pi^2/12\), with multiplicity \(n\), and the \(f_i\) are orthonormal first eigenfunctions. Direct calculation gives
\[
 g_{ii}(x)=2\lambda\cos(2kx_i)
           =\frac{\pi^2}{6}\cos\!\left(\frac{\pi x_i}{\sqrt3}\right), \tag{11.3}
\]
and, for \(i\ne j\),
\[
 g_{ij}(x)=-2\lambda\sin(kx_i)\sin(kx_j).                    \tag{11.4}
\]
Neither is quadratic. The equal mixture is
\[
 G_{I/n}=\frac{2\lambda}{n}\sum_{i=1}^n\cos(2kx_i),          \tag{11.5}
\]
also nonquadratic.

There is no contradiction with Proposition 3.2. The potential is \(V=0\). If both
\(V+tW\) and \(V-tW\) are convex for small \(t>0\), then
\(D^2W\succeq0\) and \(-D^2W\succeq0\), so \(W\) is affine. For an affine
\(W=c+b\cdot x\), the barycenter tangent equation in (3.4) gives \(b=0\); a constant does not change the law. Thus the two-sided moment-preserving potential tangent space is trivial. The multiplier \(M\) in (4.2) carries the nonquadratic stress.

The facets likewise have only affine two-sided normal graph variations. Product structure therefore supplies no polynomial stress identity at a flat-potential/faceted boundary point.

### 11.3 Isotropic regular simplex: the tangent equations can be vacuous

Let
\[
 H=\left\{z\in\mathbb R^{n+1}:\sum_{i=1}^{n+1}z_i=0\right\},
\]
let \(P=(P_1,\ldots,P_{n+1})\) be uniform on the probability simplex
\(\{p_i\ge0,\sum p_i=1\}\), and define
\[
 X=\sqrt{(n+1)(n+2)}
   \left(P-\frac1{n+1}{\bf1}\right)\in H.                    \tag{11.6}
\]
The Dirichlet moment formulas are
\[
 \operatorname{Var}(P_i)=\frac{n}{(n+1)^2(n+2)},
 \qquad
 \operatorname{Cov}(P_i,P_j)=-\frac1{(n+1)^2(n+2)}\quad(i\ne j).
\]
Therefore \(X\) is isotropic on the \(n\)-dimensional space \(H\). Its support is a regular simplex, and its interior potential is again \(V=0\).

Exactly as for the cube, a two-sided convex potential perturbation is affine and the isotropic tangent conditions kill its linear part. Hence bulk stationarity gives no information about the first eigenfunction.

There are only \(n+1\) parallel facet-offset parameters in the simplest shape family. If one insists that the raw shape variation itself be tangent to all barycenter and covariance constraints, those equations can kill every such direction. Alternatively, one may center and isotropize after each facet shift; this gives legitimate quotient paths, but it still supplies only \(n+1\) weighted facet-integral identities, with explicit affine correction terms from Section 5. It does not create pointwise normal degrees of freedom. General curved perturbations of a facet are only one-sided; requiring both signs forces the facet graph to be affine. Thus a claim that facet shifts recover a pointwise boundary Euler equation is false. Simplex symmetry can make all available facet integrals equal without imposing pointwise or rank-one structure.

Linear trial functions still give \(\lambda_1\le1\), but they are not Neumann eigenfunctions: their normal derivatives are nonzero on facets. Inferring affine eigenfunctions from the simplex's affine geometry is therefore also invalid.

### 11.4 Isotropic ball: radial symmetry does not make the first mode affine or its stress quadratic

Assume \(n\ge2\), and let \(\mu\) be uniform on the ball
\[
 B_R^n,\qquad R=\sqrt{n+2}.
\]
Since \(\mathbb E X_i^2=R^2/(n+2)\), this law is isotropic. The first nonzero Neumann modes have spherical-harmonic degree one. They can be written
\[
 f_i(x)=c\,h(r)\theta_i,\qquad r=|x|,\quad\theta=x/r,          \tag{11.7}
\]
where
\[
 h(r)=r^{-(n-2)/2}J_{n/2}(kr),\qquad k^2=\lambda_1,           \tag{11.8}
\]
and \(kR\) is the first positive root of
\[
 \left.\frac d{dr}\left(r^{-(n-2)/2}J_{n/2}(r)\right)\right|_{r=kR}=0.
\]
The constant \(c\) makes each \(f_i\) have \(L^2\)-norm one. Here is a short verification that this degree-one sector is the first nonzero one. A spherical-harmonic sector of degree \(\ell\) has radial Rayleigh quotient
\[
 \mathcal R_\ell[u]=
 \frac{\int_0^R\left(u'(r)^2r^{n-1}
  \ell(\ell+n-2)u(r)^2r^{n-3}\right)\,dr}
 {\int_0^R u(r)^2r^{n-1}\,dr}.                              \tag{11.8a}
\]
For \(\ell\ge1\), its minimum is increasing in \(\ell\), so the smallest nonradial eigenvalue occurs at \(\ell=1\). Let \(v\) be the first nonconstant radial Neumann eigenfunction. Differentiating its radial equation shows that \(w=v'\) solves the \(\ell=1\) radial equation with the same eigenvalue and with \(w(R)=0\). Thus \(w\) is an admissible trial profile in (11.8a), so the first \(\ell=1\) Neumann eigenvalue is no larger than the first nonconstant radial eigenvalue. Equality would make \(w\) satisfy both the Dirichlet condition \(w(R)=0\) and the natural Neumann condition \(w'(R)=0\), forcing \(w\equiv0\) by ODE uniqueness. Hence the inequality is strict. This proves that the first nonzero eigenvalue lies in the degree-one sector.

For the equal eigenspace mixture,
\[
 G_{I/n}(r)=\frac{c^2}{n}
 \left(h'(r)^2+(n-1)\frac{h(r)^2}{r^2}-k^2h(r)^2\right).      \tag{11.9}
\]
This is radial but not quadratic in \(x\). Indeed, up to a nonzero constant factor,
\[
 h(r)=r\left(1-\frac{k^2r^2}{2(n+2)}
 +\frac{k^4r^4}{8(n+2)(n+4)}+O(r^6)\right).
\]
Substitution into the bracket in (11.9) gives
\[
 n-2k^2r^2+
 \frac{3(n+3)}{2(n+2)^2}k^4r^4+O(r^6).                       \tag{11.10}
\]
The \(r^4\) coefficient is strictly positive. Thus rotational symmetry does not reduce the first Neumann mode to a linear function, and an averaged radial stress need not be a quadratic polynomial.

Again \(V=0\) is on the boundary of the convex-potential cone. The ball has a strictly convex shape and therefore admits two-sided shape perturbations, but it is not asserted to be a shape extremizer. One may use (6.1) to vary it; one may not use stationarity without an extremal hypothesis.

### 11.5 Vanishing regularization does not imply vanishing KKT force

The scalar example in Section 7 already disproves this inference. In potential notation, take a one-parameter curvature
\(D^2V_s=sI\) and a log barrier \(-\eta\log s\). Along a minimizing sequence with \(s=\eta\),
\[
 \eta(D^2V_s)^{-1}=I.
\]
Thus the weak limit of the regularizer's first variation is nonzero even though its coefficient tends to zero. The limiting object is precisely a complementary PSD multiplier supported on zero-curvature directions. Any proof that simply drops the regularizer term loses the principal KKT contribution.

## 12. Audited implications and remaining gap

The exact conclusions are:

1. Potential variation gives the branch matrix (3.2), not merely a scalar formula.
2. At a two-sided interior local minimum, every stress entry is a centered quadratic polynomial, (3.6).
3. At the convexity boundary, only the defect equation (4.2)--(4.3) is valid.
4. Affine deformations followed by isotropic normalization are rotations and yield no coercivity.
5. Shape variation gives the weighted Neumann Hadamard formula (6.1); facets admit only affine two-sided normal graphs.
6. Tensor-squaring converts the off-diagonal stress into \(-\lambda f\otimes f\). Exact two-sided stationarity forces a Gaussian factor and \(\lambda=1\).
7. The factorization lemma is quantitatively dimension-free, especially when the objective is \(\log\lambda\) and the stress is normalized by \(\lambda\).
8. Gaussian, cube, simplex, and ball models rule out rank-one, unrestricted-tangent, facet-pointwise, radial-affine, and vanishing-barrier shortcuts.

The unresolved statement needed for KLS is:

> **Remaining gap.** Recover dimension-free approximate two-sided, full-matrix stationarity for tensor-square near-extremizers; bound the normalized potential and shape normal-cone defects strongly enough to imply (10.1); and combine the normalized KKT equations with (10.1) to prove (10.2). Equivalently, show that a flat-potential/faceted extremizer can neither discard the mixed eigenspace stress through a scalar spectral multiplier nor carry that order-one stress through complementary PSD curvature multipliers, unless its spectral gap is already bounded below by a numerical constant.

No estimate in this report controls \(M\) or \(N\) on zero-curvature sets. Uniform measures on cubes, simplices, and balls show that those sets include the principal geometric models, not a removable exceptional class. Therefore the report does not establish \(C_P(\mu)\le C\); it isolates the exact variational rigidity mechanism and the precise normal-cone estimate still required to turn it into a dimension-free proof.
