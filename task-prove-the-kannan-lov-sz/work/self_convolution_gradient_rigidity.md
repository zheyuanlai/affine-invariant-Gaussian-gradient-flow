# Gradient rigidity under normalized self-convolution

## 0. Result and limitation

Let \(X,Y\) be independent with the same centered isotropic log-concave law
\(\mu\), and let

\[
 S=\frac{X+Y}{\sqrt2},\qquad \nu=\mathcal L(S)=\mathsf T\mu.
\]

Write

\[
 a=\lambda_1(\mu)=C_P(\mu)^{-1},\qquad
 b=\lambda_1(\nu)=C_P(\nu)^{-1}.
\]

The main unconditional conclusion of this note is

\[
 \boxed{\quad b(1-b)\leq4(b-a),\qquad
 b\geq\Phi(a):=\frac{\sqrt{9+16a}-3}{2}.\quad}        \tag{0.1}
\]

No spectral attainment is assumed.  In particular, if \(b\leq1/2\), then

\[
 b\geq\frac87a,
 \qquad C_P(\mathsf T\mu)\leq\frac78C_P(\mu).         \tag{0.2}
\]

Iteration forces the gaps of the dyadic normalized sums to converge to one.
This is a strict, dimension-free forward renormalization theorem.

It is **not** a proof of KLS.  The number of iterations needed to reach a
fixed gap is \(O(\log(1/a))\), and (0.1) gives no reverse bound on the
starting value \(a\).  A reverse/deconvolution principle or a universal
finite-step theorem is still missing.

The proof below actually needs only isotropy and a finite Poincare constant;
log-concavity supplies those hypotheses in each fixed finite dimension.

## 1. Hoeffding and energy decompositions

Take a centered \(f\) in the form domain of \(\nu\), and put

\[
 F(x,y)=f\!\left(\frac{x+y}{\sqrt2}\right),
 \qquad h(x)=\mathbb E_YF(x,Y),
\]

\[
 R(x,y)=F(x,y)-h(x)-h(y).                              \tag{1.1}
\]

Then

\[
 \mathbb E[R\mid X]=\mathbb E[R\mid Y]=0,
 \qquad
 \|f\|_{L^2(\nu)}^2=2\|h\|_{L^2(\mu)}^2+
                         \|R\|_{L^2(\mu^2)}^2.        \tag{1.2}
\]

For smooth \(f\), set \(Z=\nabla f(S)\) and

\[
 U=\mathbb E[Z\mid X],\qquad V=\mathbb E[Z\mid Y].    \tag{1.3}
\]

Differentiation under the conditional integral gives

\[
 \nabla h(X)=\frac{U}{\sqrt2},
 \quad \nabla_xR=\frac{Z-U}{\sqrt2},
 \quad \nabla_yR=\frac{Z-V}{\sqrt2}.                  \tag{1.4}
\]

Sobolev approximation and Fubini extend (1.4) to every form-domain
function.  Conditional centering makes all cross terms vanish, so if

\[
 q=\int|\nabla f|^2\,d\nu,
 \qquad E_R=\int|\nabla R|^2\,d\mu^2,
\]

then

\[
 \boxed{\quad q=2\int|\nabla h|^2\,d\mu+E_R.\quad}    \tag{1.5}
\]

The product Poincare inequality restricted to functions centered in each
variable gives

\[
 E_R\geq2a\|R\|_2^2.                                  \tag{1.6}
\]

Applying the one-factor inequality to \(h\), and using (1.2), yields

\[
 q\geq a(\|f\|_2^2-\|R\|_2^2)+E_R.                  \tag{1.7}
\]

For \(\|f\|_2=1\), (1.6)--(1.7) imply the sharp defect estimate

\[
 \boxed{\quad E_R\leq2(q-a).\quad}                    \tag{1.8}
\]

## 2. The two-predictor gradient estimate

The predictors \(U\) and \(V\) have the same mean
\(m=\mathbb EZ\); after centering, they are independent because one is a
function of \(X\) and the other a function of \(Y\).  Define the vector
Hoeffding residual

\[
 W=Z-U-V+m.
\]

It is orthogonal in \(L^2\) to both \(U-m\) and \(V-m\).  Therefore

\[
 \operatorname {Var}(Z)=2\operatorname {Var}(U)+\|W\|_2^2,
\]

whereas, by (1.4) and symmetry,

\[
 E_R=\mathbb E|Z-U|^2
     =\operatorname {Var}(U)+\|W\|_2^2.               \tag{2.1}
\]

It follows that

\[
 \boxed{
 \operatorname {Var}_\nu(\nabla f)
 \leq2E_R\leq4(q-a).}                                 \tag{2.2}
\]

The first factor two is algebraically sharp when \(W=0\), i.e. when the
gradient is additive in the two inputs.

## 3. Bottom spectral windows for \(\nu\)

Fix \(\varepsilon>0\), and choose a centered unit vector \(f\) in the
spectral subspace
\(\mathbf1_{[b,b+\varepsilon]}(A_\nu)L^2_0(\nu)\).  Write

\[
 q=\langle A_\nu f,f\rangle=b+\alpha,
 \qquad z=(A_\nu-b)f,
\]

where

\[
 0\leq\alpha\leq\varepsilon,
 \qquad \|z\|_2^2\leq\varepsilon\alpha\leq\varepsilon^2.
                                                               \tag{3.1}
\]

Let \(\ell=\mathbb E[Sf(S)]\).  Isotropy gives \(|\ell|\leq1\).  Testing
the weak generator equation against the coordinate functions gives

\[
 \mathbb E\nabla f=b\ell+\mathbb E[Sz],
 \qquad |\mathbb E[Sz]|\leq\|z\|_2\leq\varepsilon.    \tag{3.2}
\]

Consequently

\[
 \operatorname {Var}_\nu(\nabla f)
 =q-|\mathbb E\nabla f|^2
 \geq b(1-b)-2\varepsilon-\varepsilon^2.              \tag{3.3}
\]

Combining (2.2) and (3.3) gives

\[
 b(1-b)-2\varepsilon-\varepsilon^2
 \leq4(b+\varepsilon-a).
\]

Letting \(\varepsilon\downarrow0\) proves the first inequality in (0.1).
Solving the resulting quadratic

\[
 b^2+3b-4a\geq0
\]

proves the second.  This passage explicitly covers continuous spectral
edges.

## 4. Iteration and exact scope

Let \(\mu_0=\mu\), let \(\mu_{k+1}=\mathsf T\mu_k\), and put
\(a_k=\lambda_1(\mu_k)\).  Isotropy gives \(a_k\leq1\), while (0.1) gives

\[
 a_{k+1}\geq\Phi(a_k)>a_k\qquad(0<a_k<1).             \tag{4.1}
\]

Thus \(a_k\uparrow1\): a limit below one would contradict the strict
inequality \(\Phi(t)>t\) on \((0,1)\).  When \(a_{k+1}\leq1/2\), (0.2)
follows from

\[
 4a_k\leq a_{k+1}(3+a_{k+1})\leq\frac72a_{k+1}.
\]

This conclusion is consistent with exact one-dimensional models.  A
centered variance-one one-sided exponential has gap \(1/4\).  After one
normalized convolution it is a standardized \(\operatorname{Gamma}(2,1)\)
law, whose exact gap is \(4/9\).  More generally the standardized
\(\operatorname{Gamma}(2^k,1)\) iterate has

\[
 C_P=(1+2^{-k})^2.
\]

This Gamma family is distinct from the symmetric-Laplace convolution
density \(\tfrac12(1+2|s|)e^{-2|s|}\).

The forward convergence in (4.1) does not bound \(a_0\).  A sequence may
start arbitrarily close to zero and obey the recurrence for
\(O(\log(1/a_0))\) iterations.  Tensorization and the Hadamard rotation do
not reverse the estimate: the original low modes remain in the full product
even though every normalized-sum marginal improves.  This is the exact
remaining obstruction to using (0.1) as a complete KLS proof.
