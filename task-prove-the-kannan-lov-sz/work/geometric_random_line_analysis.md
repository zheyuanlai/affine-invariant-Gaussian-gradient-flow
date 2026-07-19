# Geometric random-line analysis for KLS

## 1. Random-line operator and its two Dirichlet forms

Let \(\mu\) be an isotropic log-concave probability measure on \(\mathbb R^n\). For
\(\theta\in S^{n-1}\), write
\[
 P_\theta f=\mathbb E[f(X)\mid \Pi_{\theta^\perp}X],
 \qquad P=\mathbb E_\theta P_\theta .
\]
Each \(P_\theta\) is an orthogonal projection in \(L^2(\mu)\), so \(P\) is a
self-adjoint Markov operator and
\[
 \mathcal D(f):=\langle f,(I-P)f\rangle
 =\mathbb E_\theta\mathbb E\operatorname{Var}
   (f(X)\mid\Pi_{\theta^\perp}X).
\]
Its spectral gap is
\[
 \gamma=\inf_{\operatorname{Var}_\mu f>0}
 \frac{\mathcal D(f)}{\operatorname{Var}_\mu f}.
\]

For a conditional one-dimensional law \(\nu(dt)=\rho(t)dt\), with mean \(m\),
define its Stein weight
\[
 \tau_\nu(t)=\frac{1}{\rho(t)}\int_t^\infty(s-m)\rho(s)\,ds.
\]
The one-dimensional weighted Poincare inequality gives
\[
 \operatorname{Var}_\nu g\leq \int \tau_\nu(t)g'(t)^2\,d\nu(t),
 \qquad \int\tau_\nu\,d\nu=\operatorname{Var}_\nu(t).
\]
Consequently
\[
 \mathcal D(f)\leq \mathcal W_\tau(f)
 :=\mathbb E_{\theta,X}\left[
 \tau_{\theta,\Pi X}(X)(\theta\cdot\nabla f(X))^2\right]
 =\int \nabla f(x)^TM_\tau(x)\nabla f(x)\,d\mu(x),
\]
where
\[
 M_\tau(x)=\mathbb E_\theta[
 \tau_{\theta,\Pi X}(x)\theta\theta^T].
\]
If \(S_\theta=\mathbb E_z\operatorname{Var}(t\mid z)\), isotropy and total
variance give \(S_\theta\leq\operatorname{Var}(\theta\cdot X)=1\). Hence
\[
 \int M_\tau\,d\mu=\mathbb E_\theta[S_\theta\theta\theta^T]
 \preceq \frac1n I,
 \qquad
 \int\operatorname{tr}M_\tau\,d\mu=\mathbb E_\theta S_\theta\leq1.
\]
This controls constant gradients, but not gradients correlated with the rare
locations at which \(M_\tau(x)\) is large.

For comparison, if \(c_{\theta,z}\) is the ordinary one-dimensional Poincare
constant, then
\[
 \mathcal D(f)\leq \mathcal W(f)
 :=\mathbb E[c_{\theta,z}(\theta\cdot\nabla f)^2].
\]

## 2. Exact logical obstruction

Define
\[
 \beta=\sup_{f\not\equiv\mathrm{const}}
 \frac{\mathcal D(f)}{\int|\nabla f|^2d\mu},
 \qquad
 \lambda=\inf_{f\not\equiv\mathrm{const}}
 \frac{\int|\nabla f|^2d\mu}{\operatorname{Var}_\mu f}.
\]
Testing \(\gamma\) on an approximate minimizer for \(\lambda\) gives
\[
 \gamma\leq\beta\lambda,
 \qquad\text{therefore}\qquad
 C_P(\mu)=\lambda^{-1}\leq\frac\beta\gamma.
\]
Thus any dimension-free cancellation \(\beta/\gamma=O(1)\) is already at least
as strong as KLS. Replacing \(\beta\) by the conditional-weight bound is not
possible: the pointwise conditional weight is much larger than its global
average, even for the ball and cube.

## 3. Isotropic cube

Let \(K=[-\sqrt3,\sqrt3]^n\). Fix a unit direction \(\theta\), put
\(s=\|\theta\|_1\), and let \(T_+\) be the forward distance from a uniform point
of \(K\) to the boundary along \(\theta\). Then
\[
 \Pr(T_+>t)=\prod_i\left(1-\frac{|\theta_i|t}{2\sqrt3}\right)_+.
\]
Using \(\prod(1-u_i)\leq e^{-\sum u_i}\), and using
\(\log(1-u)\geq-4u/3\) for \(0\leq u\leq1/4\), gives
\[
 \frac{e^{-1/3}(2\sqrt3)^2}{16s^2}
 \leq \mathbb ET_+^2
 \leq \frac{2(2\sqrt3)^2}{s^2}.
\]
If \(\ell=T_++T_-\) is the chord length, its conditional coordinate variance is
\(\ell^2/12\). From
\(T_+^2+T_-^2\leq\ell^2\leq2(T_+^2+T_-^2)\),
\[
 \frac{e^{-1/3}}{8s^2}\leq
 S_\theta:=\mathbb E_x\frac{\ell^2}{12}
 \leq\frac8{s^2}.
\]
For uniform \(\theta\), sphere concentration gives
\(\mathbb E\|\theta\|_1^{-2}=\Theta(1/n)\), so
\(S:=\mathbb E_\theta S_\theta=\Theta(1/n)\).

For \(f_j(x)=x_j\), summing over coordinates gives
\[
 \sum_j\mathcal D(f_j)=S.
\]
Hyperoctahedral symmetry therefore yields
\[
 \mathcal D(f_j)=\Theta(n^{-2}),
 \qquad \gamma\leq Cn^{-2}.
\]
The exact Poincare constant of a uniform interval of length \(\ell\) is
\(\ell^2/\pi^2\), hence
\[
 \mathcal W(f_j)=\frac{12}{\pi^2}\mathcal D(f_j)=\Theta(n^{-2}).
\]

At the cube center, the line interval is
\([-a,a]\), with \(a=\sqrt3/\|\theta\|_\infty\). Thus
\[
 c_\theta(0)=\frac{12}{\pi^2\|\theta\|_\infty^2},
 \quad
 M(0)=\frac{12}{\pi^2n}
 \mathbb E\|\theta\|_\infty^{-2}\,I
 =\Theta((\log n)^{-1})I.
\]
Here \(\mathbb E\|\theta\|_\infty^{-2}=\Theta(n/\log n)\). The Stein weight at
the center is \(a^2/2\), so likewise
\[
 M_\tau(0)=\frac{3}{2n}
 \mathbb E\|\theta\|_\infty^{-2}\,I
 =\Theta((\log n)^{-1})I.
\]
Consequently \(\gamma^{-1}M(0)\) and \(\gamma^{-1}M_\tau(0)\) are at least of
order \(n^2/\log n\). Localized high-frequency test functions show that a
quadratic-form bound \(\mathcal W\leq C\gamma\int|\nabla f|^2\) is impossible.

## 4. Isotropic Euclidean ball

Let \(K=B_{\sqrt{n+2}}\), so \(\operatorname{Cov}(X)=I\). For a line through
\(x\) in direction \(\theta\), set
\[
 a^2=R^2-|x|^2+(\theta\cdot x)^2.
\]
The conditional interval is \([-a,a]\), and therefore
\[
 \operatorname{Var}(t\mid z)=a^2/3,
 \qquad c_{\theta,z}=4a^2/\pi^2.
\]
Since \(R^2=n+2\), \(\mathbb E|X|^2=n\), and
\(\mathbb E(\theta\cdot X)^2=1\),
\[
 S_\theta=1
\]
for every \(\theta\). Hence, for a unit linear function,
\[
 \mathcal D(f)=1/n,
 \qquad \mathcal W(f)=12/(\pi^2n).
\]
Moreover \(Pf=(1-1/n)f\), so \(\gamma\leq1/n\).

The ordinary conditional-weight matrix is
\[
 M(x)=\frac4{\pi^2}\left[
 \frac{R^2-|x|^2}{n}I+
 \frac{|x|^2I+2xx^T}{n(n+2)}
 \right].
\]
In particular \(M(0)=4(n+2)/(\pi^2n)I=\Theta(I)\), already a factor \(n\)
larger than the linear/global scale.

For the Stein weight,
\[
 \tau=(a^2-t^2)/2=(R^2-|x|^2)/2,
 \qquad
 M_\tau(x)=\frac{R^2-|x|^2}{2n}I.
\]
Thus \(\int M_\tau d\mu=I/n\), but \(M_\tau(0)\sim I/2\). Even the optimal
one-dimensional weight cannot be bounded pointwise, or as a multiplication
quadratic form for arbitrary gradients, by its average scale.

## 5. Regular isotropic simplex

Choose regular-simplex vertices \(v_0,\dots,v_n\) satisfying
\[
 |v_i|^2=n(n+2),\qquad \langle v_i,v_j\rangle=-(n+2),\ i\ne j.
\]
Write \(\theta=\sum_i b_iv_i\), where \(\sum_i b_i=0\). Unit length means
\[
 \sum_i b_i^2=\frac1{(n+1)(n+2)}.
\]
For barycentric coordinates \(\lambda\sim\operatorname{Dirichlet}(1,\dots,1)\),
put
\[
 B=\sum_{b_i>0}b_i=\sum_{b_i<0}(-b_i)=\|b\|_1/2.
\]
The forward endpoint distance satisfies
\[
 T_+=\min_{b_i<0}\frac{\lambda_i}{-b_i},
 \qquad
 \Pr(T_+>t)=(1-Bt)^n,quad 0\leq t\leq B^{-1}.
\]
Therefore
\[
 \mathbb ET_+^2=\frac2{B^2(n+1)(n+2)}
\]
and
\[
 \frac1{3B^2(n+1)(n+2)}
 \leq S_\theta\leq
 \frac2{3B^2(n+1)(n+2)}.
\]
For random \(\theta\), the vector \(b\) is uniform on the sphere in the
zero-sum subspace and
\(\|b\|_1/\|b\|_2=\Theta(\sqrt n)\) in negative second moment. Hence
\[
 \mathbb E_\theta S_\theta=\Theta(1/n).
\]
Simplex symmetry and summation over an orthonormal basis give, for every unit
linear coordinate,
\[
 \mathcal D(f)=\Theta(n^{-2}),
 \qquad \mathcal W(f)=\Theta(n^{-2}),
 \qquad \gamma\leq Cn^{-2}.
\]

## 6. Convexification obstruction

For completeness, multiplicative convexification fails in dimension one. Let
\(\mu\) be uniform on \([-\sqrt3,\sqrt3]\), and
\[
 A=[-\sqrt3,-\sqrt3/2]\cup[\sqrt3/2,\sqrt3].
\]
Then \(\mu(A)=1/2\), \(\operatorname{conv}(A)=\operatorname{supp}\mu\), but
\[
 \int d(x,A)d\mu(x)=\frac{\sqrt3}{8},
 \qquad
 \int d(x,\operatorname{conv}A)d\mu(x)=0.
\]
An additive universal convexification estimate, combined with a universal
convex-set distance estimate, would itself yield the half-mass first-moment
form of KLS. Conversely KLS immediately implies the additive estimate since
\(0\leq d_A-d_{\operatorname{conv}A}\leq d_A\). It is therefore not an easier
reduction.
