# Integrated curvature and Langevin coupling: exact identities, no-go results, and a scale-resolved replacement

## 1. Conclusions

Let

\[
 d\mu(x)=Z^{-1}e^{-V(x)}\,dx,
 \qquad L=\Delta-\nabla V\cdot\nabla,
\]

with \(V\) smooth and convex.  The following conclusions are rigorous.

1. If \(X,Y\stackrel{\mathrm{iid}}\sim\mu\), then, without any isotropic assumption,
   \[
   \mathbb E\langle \nabla V(X)-\nabla V(Y),X-Y\rangle=2n.
   \]
   If \(\operatorname{Cov}(\mu)\preceq I\), this is at least the mean-square separation \(\mathbb E|X-Y|^2\).  There is in fact a whole weighted hierarchy of such identities.

2. This mean curvature is not pathwise curvature.  There are explicit, smooth, strictly convex, exactly isotropic product potentials for which the chord monotonicity quotient is arbitrarily close to zero at fixed interior pairs and on a set carrying almost all the mass, while \(\mathbb E\Delta V\) tends to infinity.

3. A quadratic \(T_2\) inequality is impossible for the class in the question.  A smooth, strictly convex, isotropic one-dimensional log-concave law with asymptotically linear potential fails every finite \(T_2\) inequality.  Thus a universal \(T_2\Rightarrow\) Poincare route cannot prove KLS.

4. No coupling can give uniform Euclidean \(W_1\) or \(W_2\) exponential contraction from every pair of points unless \(V\) is uniformly strongly convex.  More strongly, even a uniform time-integrated point-start \(W_2^2\) estimate fails for a fixed smooth isotropic smoothed-Laplace law: two starting points can be translated arbitrarily far into one linear tail and remain separated for an arbitrarily long time.

5. Starting the synchronous coupling from iid stationarity does give the desired curvature at time zero.  It is not stable.  If
   \[
   A(t)=\mathbb E\langle \nabla V(X_t)-\nabla V(Y_t),X_t-Y_t\rangle,
   \]
   then
   \[
   A(0)=2n,\qquad A'(0)=-4\mathbb E|\nabla V(X)|^2=-4\mathbb E\Delta V(X).
   \]
   In the explicit isotropic power family below, \(-A'(0)/A(0)\sim 2p/3\), with \(p\to\infty\).  Hence iid mean curvature can be destroyed on an arbitrarily short time scale.

6. The precise missing quantity is not entropy or mutual information.  A coupling of two stationary smoothed-Laplace marginals can use only one common sign bit (mutual information \(\log 2\)) and make the averaged chord curvature tend to zero.  The correct defect is a *radial pair score*.  For a pair density \(q=d\pi/d\mu^{\otimes2}\), an exact integration-by-parts identity isolates that defect at every distance scale.  This gives a concrete new functional to track under a hybrid or reflection coupling.

These facts rule out the proposed route in its naive forms.  They do not disprove KLS; they identify the additional correlation-control estimate a successful coupling proof would have to establish.

## 2. Basic integration-by-parts identities

Write \(g=\nabla V\), \(H=\nabla^2V\), and \(m=\mathbb EX\).  All identities below follow first for smooth potentials with sufficient decay.  Truncation gives them whenever the displayed quantities are integrable.

### Lemma 2.1 (score, chord, and mean-Hessian identities)

For \(X\sim\mu\),

\[
 \mathbb Eg(X)=0,
 \qquad
 \mathbb E\big[g(X)(X-m)^T\big]=I.
 \tag{2.1}
\]

For independent \(X,Y\sim\mu\),

\[
 \mathbb E\big[(g(X)-g(Y))(X-Y)^T\big]=2I,
 \tag{2.2}
\]

and hence, with

\[
 M(x,y):=\langle g(x)-g(y),x-y\rangle,
\]

one has

\[
 \boxed{\mathbb EM(X,Y)=2n.}
 \tag{2.3}
\]

Moreover,

\[
 \mathbb EH(X)=\mathbb E[g(X)g(X)^T]
 \succeq \operatorname{Cov}(X)^{-1}.
 \tag{2.4}
\]

In particular, \(\operatorname{Cov}(X)\preceq I\) implies \(\mathbb EH(X)\succeq I\).  It also implies

\[
 \frac{\mathbb EM(X,Y)}{\mathbb E|X-Y|^2}
 =\frac{n}{\operatorname{tr}\operatorname{Cov}(X)}\ge 1.
 \tag{2.5}
\]

#### Proof

Integration by parts gives

\[
 \int \partial_iV\,d\mu=0,
 \qquad
 \int (x_j-m_j)\partial_iV\,d\mu=\delta_{ij}.
\]

Independence and \(\mathbb Eg=0\) then give (2.2).  Taking the trace gives (2.3).  A second integration by parts gives \(\mathbb EH=\mathbb E[gg^T]\).  Finally, the covariance matrix of the random vector \((X-m,g(X))\) is

\[
 \begin{pmatrix}
 C&I\\ I&\mathbb E[gg^T]
 \end{pmatrix}\succeq0,
 \qquad C=\operatorname{Cov}(X).
\]

Its Schur complement is (2.4).  Also \(\mathbb E|X-Y|^2=2\operatorname{tr}C\), which proves (2.5).  \(\square\)

There is also a literal chord formula:

\[
 M(x,y)=\int_0^1 (x-y)^TH(y+s(x-y))(x-y)\,ds.
 \tag{2.6}
\]

Thus (2.3) says that curvature averaged simultaneously over random chords, chord directions, chord lengths, and chord locations has total mass \(2n\).  It says nothing about a prescribed chord or about the correlated pair law created by a coupling.

### Lemma 2.2 (the iid multiscale chord hierarchy)

Let \(R=|X-Y|\), and let \(w:[0,\infty)\to\mathbb R\) be \(C^1\), with enough integrability to justify the calculation.  Then

\[
 \boxed{
 \mathbb E\big[w(R)M(X,Y)\big]
 =2n\,\mathbb Ew(R)+2\mathbb E\big[Rw'(R)\big].
 }
 \tag{2.7}
\]

In particular, for every admissible \(a\ge0\),

\[
 \boxed{
 \mathbb E\big[R^aM(X,Y)\big]
 =2(n+a)\mathbb ER^a.
 }
 \tag{2.8}
\]

#### Proof

Apply integration by parts under \(\mu(dx)\mu(dy)\) to the vector field

\[
 (x-y)w(|x-y|).
\]

The two divergences contribute \(2nw(R)\).  Since

\[
 (\nabla_x-\nabla_y)R=2\frac{x-y}{R},
\]

the derivative of the weight contributes \(2Rw'(R)\).  This proves (2.7); taking \(w(r)=r^a\) gives (2.8).  \(\square\)

For a Gaussian, (2.8) reduces to the usual recursion for chi moments.  The point is that it holds for every smooth log-concave law.  It will be extended to correlated pair laws in Section 10.

## 3. A smooth isotropic flat-core family

The following family is useful because it requires no hard-support limiting argument.

Fix an even integer \(p\ge4\), and let \(W_p\) have density

\[
 c_p\exp(-|w|^p/p).
\]

Set

\[
 s_p^2=\mathbb EW_p^2
 =p^{2/p}\frac{\Gamma(3/p)}{\Gamma(1/p)},
 \qquad X_p=W_p/s_p.
\]

Then \(X_p\) has variance one and potential

\[
 V_p(x)=\frac{s_p^p|x|^p}{p}.
 \tag{3.1}
\]

This potential is \(C^\infty\), strictly convex, and has full support.  Products of (3.1) are exactly isotropic in every dimension.

As \(p\to\infty\), \(W_p\) converges to the uniform law on \([-1,1]\), so

\[
 s_p^2\longrightarrow\frac13.
 \tag{3.2}
\]

At the fixed pair \(x=0,y=1\),

\[
 \frac{M(0,1)}{|0-1|^2}=s_p^p\longrightarrow0
 \tag{3.3}
\]

exponentially fast.  Therefore no positive pointwise chord-curvature constant follows from isotropicity, even in dimension one.

The Hessian is

\[
 V_p''(x)=(p-1)s_p^p|x|^{p-2}
 =(p-1)s_p^2(s_p|x|)^{p-2}.
 \tag{3.4}
\]

For any fixed \(a<\sqrt3\), (3.2) gives constants \(C<\infty\) and \(\rho<1\) such that

\[
 \sup_{|x|\le a}V_p''(x)\le Cp\rho^p.
 \tag{3.5}
\]

On the other hand,

\[
 h_p:=\mathbb EV_p''(X_p)
 =p(p-1)
 \frac{\Gamma(3/p)\Gamma(1-1/p)}{\Gamma(1/p)^2}
 \sim \frac p3.
 \tag{3.6}
\]

For \(a=(1-\delta)\sqrt3\), the event in (3.5) has probability tending to \(1-\delta\), whereas the mean curvature (3.6) diverges.  Thus the mean Hessian is carried by a thinner and thinner boundary region.  In the \(n\)-fold product,

\[
 \mathbb E\Delta V_p(X)=nh_p\sim np/3,
\]

while every coordinate has the flat-core behavior (3.5).

This is the basic distinction between mean and pathwise curvature: \(\mathbb EH\succeq I\) can be true by an arbitrarily large margin while \(H\) is almost zero where a particular trajectory or test-field gradient is located.

## 4. Cube, product exponential, and simplex-like tests

### 4.1 Product exponential

The centered variance-one Laplace law has

\[
 V(x)=\sqrt2|x|,
 \qquad V'(x)=\sqrt2\,\operatorname{sgn}(x)
\]

away from the origin.  Hence

\[
 M(x,y)=0\quad\text{if }xy>0,
\]

whereas, if \(xy<0\),

\[
 M(x,y)=2\sqrt2(|x|+|y|).
\]

Direct integration gives \(\mathbb EM(X,Y)=2\).  All chord curvature is supplied by pairs that straddle the kink.  Products tensorize this example coordinatewise.

This nonsmooth example is only illustrative; the smoothed-Laplace family in Section 5 gives a smooth, strictly convex, exactly isotropic version of the same obstruction.

### 4.2 Cube

The family (3.1) converges to the isotropic interval \([-\sqrt3,\sqrt3]\), and its products converge to the isotropic cube.  Inside a fixed shrunken cube, the pointwise Hessian tends to zero exponentially in \(p\), while \(\mathbb E\Delta V_p\sim np/3\).  Consequently, passing to a hard cube loses the curvature as an ordinary bulk function: it becomes a boundary local-time or normal-cone term.

One must therefore not transfer (2.3) to the uniform cube by declaring \(\nabla V=0\) in the interior.  The missing \(2n\) is a boundary contribution.

### 4.3 Simplex-like potentials

The one-sided exponential \(E-1\), where \(E\sim\operatorname{Exp}(1)\), has mean zero, variance one, and potential

\[
 V(x)=x+1\quad (x>-1),\qquad V(x)=+\infty\quad (x\le-1).
\]

Its interior gradient is constant, so every interior chord has \(M=0\).  Again, the integration-by-parts mass lies entirely at the hard facet.  Products give a shifted orthant model.  Conditioning analogous one-sided exponentials by a sum constraint gives the same facet mechanism behind simplex models.

There is a standard smooth quantitative approximation for any isotropic convex body \(K\), including a regular simplex.  Let \(d_K\) be distance to \(K\), let \(\rho_\delta\) be a smooth mollifier, and set

\[
 U_{a,\varepsilon,\delta}(x)
 =a(\rho_\delta*d_K^2)(x)+\varepsilon|x|^2.
 \tag{4.1}
\]

This is smooth, convex, integrable, and strictly convex.  On the \(\delta\)-interior of \(K\), its Hessian is exactly \(2\varepsilon I\).  Taking \(a\to\infty\), then \(\delta,\varepsilon\to0\), and finally applying the covariance-whitening affine map gives smooth exactly isotropic measures converging to the uniform law on \(K\).  Fixed interior chord curvature tends to zero, while (2.3) remains \(2n\) at every smooth stage and is pushed into a narrowing facet layer.  Hence neither cubes nor simplices turn iid mean curvature into pathwise curvature.

## 5. Quadratic \(T_2\) is impossible

Here is a smooth, full-support, strictly convex, isotropic counterexample in dimension one.

For \(\varepsilon>0\), let \(W_\varepsilon\) have density proportional to

\[
 \exp\big(-\sqrt{w^2+\varepsilon^2}\big),
\]

let \(s_\varepsilon^2=\operatorname{Var}(W_\varepsilon)\), and put \(X_\varepsilon=W_\varepsilon/s_\varepsilon\).  Its potential is

\[
 V_\varepsilon(x)
 =\sqrt{s_\varepsilon^2x^2+\varepsilon^2},
 \tag{5.1}
\]

up to an additive constant.  It satisfies

\[
 |V_\varepsilon'(x)|\le s_\varepsilon,
 \qquad
 V_\varepsilon''(x)
 =\frac{s_\varepsilon^2\varepsilon^2}
 {(s_\varepsilon^2x^2+\varepsilon^2)^{3/2}}>0.
 \tag{5.2}
\]

Thus (5.1) is smooth, strictly convex, and exactly isotropic.

### Proposition 5.1

The law \(\mu_\varepsilon\) of \(X_\varepsilon\) satisfies no finite quadratic \(T_2\) inequality.

#### Proof

Let \(\nu_h\) be the translate of \(\mu_\varepsilon\) by \(h>0\).  The means differ by \(h\), so

\[
 W_2(\nu_h,\mu_\varepsilon)^2\ge h^2.
\]

On the other hand,

\[
 \operatorname{Ent}(\nu_h\mid\mu_\varepsilon)
 =\mathbb E\big[V_\varepsilon(X_\varepsilon+h)
                   -V_\varepsilon(X_\varepsilon)\big]
 \le s_\varepsilon h.
\]

If \(W_2^2\le2C\operatorname{Ent}\) held with finite \(C\), then

\[
 h^2\le2Cs_\varepsilon h
\]

for every \(h\), a contradiction.  \(\square\)

The same counterexample works in every dimension by taking products and translating one coordinate.  This is a tail obstruction, not a dimensional one: \(T_2\) forces Gaussian-type transport behavior, while general log-concave measures may have exponential tails.  Poincare remains possible for these laws, so \(T_2\) is strictly too strong a target.

A tail-compatible transport target would need a quadratic-linear cost, for example

\[
 \alpha_R(r)=
 \begin{cases}
 r^2,&r\le R,\\
 2Rr-R^2,&r>R.
 \end{cases}
 \tag{5.3}
\]

The local quadratic part can linearize to Poincare, while the linear part is consistent with the translation computation above.  This does not prove such an inequality, but it removes the categorical \(T_2\) obstruction.

## 6. Uniform point-start contraction is also impossible

### Lemma 6.1 (Euclidean contraction forces strong convexity)

Suppose that for some \(c>0\), all \(x,y\), and all sufficiently small \(t>0\),

\[
 W_2(P_t\delta_x,P_t\delta_y)
 \le e^{-ct}|x-y|.
 \tag{6.1}
\]

Then

\[
 \langle g(x)-g(y),x-y\rangle\ge c|x-y|^2
 \tag{6.2}
\]

for every \(x,y\).  The same conclusion follows from a uniform \(W_1\) contraction.

#### Proof

Let \(m_t(x)=\mathbb E_xX_t\).  As \(t\downarrow0\),

\[
 m_t(x)=x-tg(x)+o(t).
\]

Wasserstein distance is at least the distance between the means.  Squaring that lower bound and comparing it with (6.1) gives

\[
 |x-y|^2-2tM(x,y)+o(t)
 \le |x-y|^2-2ct|x-y|^2+o(t),
\]

which is (6.2).  For \(W_1\), project the difference of means onto \((x-y)/|x-y|\).  \(\square\)

The family (3.1), particularly (3.3), contradicts (6.2) for every universal \(c>0\).  Since the cost of any particular coupling dominates the optimal Wasserstein cost, neither synchronous nor reflection coupling can evade this obstruction for arbitrary starting pairs in the Euclidean metric.

There is also no uniform *integrated* point-start estimate.

### Proposition 6.2 (far-tail obstruction to integrated contraction)

Fix one potential (5.1).  There is no finite \(C\) such that

\[
 \int_0^\infty
 W_2(P_t\delta_x,P_t\delta_y)^2\,dt
 \le C|x-y|^2
 \tag{6.3}
\]

for every \(x,y\in\mathbb R\).

#### Proof

Let \(Z_t^z\) solve the Langevin equation from \(z\).  The derivative of the one-dimensional stochastic flow is

\[
 J_t^z:=\partial_zZ_t^z
 =\exp\left(-\int_0^tV_\varepsilon''(Z_u^z)\,du\right).
 \tag{6.4}
\]

Put \(T_R=\sqrt R\), and take \(z\in[R,R+1]\).  Since 

\[
 |V_\varepsilon'|\le s_\varepsilon,
 \qquad
 V_\varepsilon''(x)\le \frac{\varepsilon^2}{s_\varepsilon x^3}
 \quad(x>0),
\]

the Brownian reflection bound implies

\[
 \mathbb P\left(\inf_{0\le u\le T_R}Z_u^z<R/2\right)
 \le 2\exp(-cR^{3/2})
 \tag{6.5}
\]

for all sufficiently large \(R\), with \(c>0\) depending only on the fixed potential.  Indeed, the deterministic drift over this interval is \(O(\sqrt R)\), while a downward Brownian excursion of order \(R\) is required.

On the complementary event, the exponent in (6.4) is at most \(CR^{-5/2}\).  Hence, uniformly for \(0\le t\le T_R\) and \(z\in[R,R+1]\),

\[
 \partial_z\mathbb E Z_t^z=\mathbb EJ_t^z\ge\frac34
\]

once \(R\) is large.  Thus, for \(0<r\le1\),

\[
 |\mathbb EZ_t^{R+r}-\mathbb EZ_t^R|\ge\frac34r,
 \qquad 0\le t\le T_R.
\]

It follows that

\[
 \int_0^\infty W_2(P_t\delta_{R+r},P_t\delta_R)^2\,dt
 \ge \frac9{16}r^2\sqrt R.
\]

Letting \(R\to\infty\) contradicts (6.3).  \(\square\)

This proposition is relevant to both synchronous and reflection constructions.  Averaging the starting points against \(\mu\) could suppress the far-tail event, but an arbitrary-pair theorem cannot be the intermediate statement.

## 7. What synchronous coupling actually gives

Let \(X_t,Y_t\) solve

\[
 dX_t=\sqrt2\,dB_t-g(X_t)\,dt,
 \qquad
 dY_t=\sqrt2\,dB_t-g(Y_t)\,dt
\]

with the same Brownian motion.  Put

\[
 D(t)=\mathbb E|X_t-Y_t|^2,
 \qquad
 A(t)=\mathbb EM(X_t,Y_t).
\]

Convexity gives \(A(t)\ge0\), and exactly

\[
 \boxed{D'(t)=-2A(t).}
 \tag{7.1}
\]

Suppose \(X_0,Y_0\) are iid with law \(\mu\).  Each marginal remains \(\mu\), but the pair immediately becomes correlated.

### Lemma 7.1 (initial loss of iid curvature)

Under the preceding initialization,

\[
 D(0)=2\operatorname{tr}\operatorname{Cov}(\mu),
 \qquad A(0)=2n,
 \tag{7.2}
\]

and

\[
 \boxed{
 A'(0)=-4\mathbb E|g(X)|^2=-4\mathbb E\Delta V(X).
 }
 \tag{7.3}
\]

If \(\mu\) is isotropic and \(k(t)=A(t)/D(t)\), then

\[
 k(0)=1,
 \qquad
 k'(0)=2\left(1-\frac1n\mathbb E\Delta V(X)\right).
 \tag{7.4}
\]

#### Proof

Only (7.3) needs proof.  Relative to \(\mu^{\otimes2}\), the forward generator of the common-noise pair has the cross term

\[
 2\sum_i\partial_{x_i}\partial_{y_i}.
\]

If \(q_t\) is the pair density relative to \(\mu^{\otimes2}\), then at the independent initial law

\[
 \dot q_0(x,y)=2g(x)\cdot g(y).
\]

Consequently,

\[
 A'(0)=2\mathbb E\big[M(X,Y)g(X)\cdot g(Y)\big].
\]

Expanding \(M\), independence, \(\mathbb Eg=0\), and

\[
 \mathbb E[X_i g_j(X)]=\delta_{ij}
\]

leave exactly \(-4\mathbb E|g(X)|^2\).  The last equality in (7.3) is (2.4) before taking the trace.  Equation (7.4) follows from (7.1).  \(\square\)

For the \(n\)-fold product of (3.1), equations (3.6) and (7.3) give

\[
 A(0)=2n,
 \qquad
 A'(0)=-4nh_p,
 \qquad
 -\frac{A'(0)}{A(0)}=2h_p\sim\frac{2p}{3}.
 \tag{7.5}
\]

Therefore there is no universal \(C\) for which the iid curvature obeys even the local persistence estimate

\[
 A(t)\ge e^{-Ct}A(0)
\]

for every smooth isotropic log-concave potential.  The very large mean Hessian is not helpful: under common noise it makes the initially independent chord alignment disappear arbitrarily fast.

Equation (7.1) also explains why an unweighted time integral of \(A\) is not by itself a result:

\[
 2\int_0^T A(t)\,dt=D(0)-D(T).
\]

It is an exact accounting identity.  A Poincare-strength conclusion requires a lower bound on \(A(t)\) relative to the remaining distance, not merely its initial value or its total spent curvature.

## 8. Stationary marginals and small mutual information do not preserve curvature

Use the smooth isotropic law \(\mu_\varepsilon\) from (5.1).  Define a coupling \(\pi_\varepsilon\) by choosing a fair sign \(S\in\{-1,1\}\), then drawing \(X,Y\) independently from \(\mu_\varepsilon\) conditioned on having sign \(S\).  Both marginals are exactly \(\mu_\varepsilon\), and

\[
 \frac{d\pi_\varepsilon}{d\mu_\varepsilon^{\otimes2}}(x,y)
 =2\mathbf1_{\{xy>0\}}.
 \tag{8.1}
\]

Thus

\[
 \operatorname{Ent}(\pi_\varepsilon\mid
 \mu_\varepsilon^{\otimes2})=\log2.
 \tag{8.2}
\]

As \(\varepsilon\downarrow0\), the marginals converge to the variance-one Laplace law, \(V_\varepsilon'(x)\to\sqrt2\operatorname{sgn}(x)\) for \(x\ne0\), and dominated convergence gives

\[
 \mathbb E_{\pi_\varepsilon}M(X,Y)\longrightarrow0.
 \tag{8.3}
\]

Meanwhile, two conditionally positive variance-one Laplace variables are iid exponentials of rate \(\sqrt2\), so

\[
 \mathbb E_{\pi_\varepsilon}|X-Y|^2\longrightarrow1.
 \tag{8.4}
\]

Hence there is no universal \(c>0\) such that

\[
 \mathbb E_\pi M(X,Y)
 \ge c\mathbb E_\pi|X-Y|^2
\]

for all couplings with two stationary marginals, even if their mutual information is at most \(\log2\).  The iid identity can be destroyed by one common bit.  If a smooth positive pair density is desired, take an odd smooth \(h_\delta\), with \(|h_\delta|\le1\) and \(h_\delta\to\operatorname{sgn}\), and set

\[
 q_{\rho,\delta}(x,y)
 =1+\rho h_\delta(x)h_\delta(y),\qquad 0<\rho<1.
\]

Symmetry gives both marginals exactly \(\mu_\varepsilon\); the density is smooth and positive, and \(q_{\rho,\delta}\le2\) gives relative entropy at most \(\log2\).  Sending \(\rho\uparrow1\), \(\delta\downarrow0\), and then \(\varepsilon\downarrow0\) gives the same limit (8.3)--(8.4).

This is exactly the type of sign alignment that common noise encourages in a product-exponential potential.  Marginal stationarity and a coarse correlation budget therefore cannot close (7.1).

## 9. Averaged semigroup gradients: the exact endpoint and the failed shortcut

For \(G_f(t)=\int|\nabla P_tf|^2d\mu\), the Bochner identity is

\[
 G_f'(t)
 =-2\int\left(
 \|\nabla^2P_tf\|_{\mathrm{HS}}^2
 +\langle H\nabla P_tf,\nabla P_tf\rangle
 \right)d\mu.
 \tag{9.1}
\]

Although \(\mathbb EH\succeq I\), it is false that the curvature term in (9.1) controls \(G_f(t)\).  This already fails in dimension one for (3.1).  Fix

\[
 0<a<b<\sqrt3
\]

and choose a smooth \(0\le\eta\le1\) that is one on \([-a,a]\) and supported in \([-b,b]\).  Let \(f'=\eta\).  By (3.5),

\[
 \frac{\int V_p''(f')^2d\mu_p}
 {\int(f')^2d\mu_p}
 \le Cp\rho^p\longrightarrow0.
 \tag{9.2}
\]

The Hessian term \(\int(f'')^2d\mu_p\) has not vanished, so (9.2) is not a spectral-gap counterexample.  It proves precisely that replacing the weighted Hessian in (9.1) by its unweighted mean is invalid.  A successful proof would have to use the interaction of the two terms in (9.1), which is already a one-form spectral problem of KLS strength.

The exact semigroup accounting identity is

\[
 \operatorname{Var}_\mu(f)
 =2\int_0^\infty G_f(t)\,dt.
 \tag{9.3}
\]

Therefore an estimate

\[
 \int_0^\infty G_f(t)\,dt
 \le C\int|\nabla f|^2d\mu
\]

is exactly a Poincare inequality (up to the factor two).  It cannot be assumed as an intermediate averaged-gradient statement.  The iid curvature identity supplies no valid shortcut from \(\mathbb EH\) to the weighted quantity in (9.1).

## 10. The radial pair-score identity: a scale-resolved replacement functional

The loss of iid curvature under correlation admits an exact formula.

Let \(\pi\) have a smooth positive density \(q=d\pi/d\mu^{\otimes2}\).  Set

\[
 R=|x-y|,
 \qquad e=\frac{x-y}{|x-y|},
 \qquad
 S_q(x,y)=e\cdot(\nabla_x-\nabla_y)\log q(x,y).
 \tag{10.1}
\]

The scalar \(S_q\) is the score of the pair correlation in the relative radial direction.

### Lemma 10.1 (weighted correlated-chord identity)

For every admissible \(C^1\) radial weight \(w\),

\[
 \boxed{
 \begin{aligned}
 \mathbb E_\pi[w(R)M(X,Y)]
 &=2n\mathbb E_\pi w(R)
   +2\mathbb E_\pi[Rw'(R)]\\
 &\quad+\mathbb E_\pi[Rw(R)S_q(X,Y)].
 \end{aligned}
 }
 \tag{10.2}
\]

Equivalently, for \(a\ge0\),

\[
 \boxed{
 \mathbb E_\pi[R^aM]
 =2(n+a)\mathbb E_\pi R^a
  +\mathbb E_\pi[R^{a+1}S_q].
 }
 \tag{10.3}
\]

#### Proof

Apply the same pair integration by parts as in Lemma 2.2 to \(q(x,y)w(R)(x-y)\).  Differentiating \(q\) gives

\[
 (x-y)\cdot(\nabla_x-\nabla_y)\log q=RS_q;
\]

the other two terms are unchanged.  \(\square\)

Thus the only obstruction to the iid hierarchy is an explicit radial score term.  Define the scale-\(a\) curvature-retention coefficient

\[
 \mathfrak r_a(\pi)
 :=1+
 \frac{\mathbb E_\pi[R^{a+1}S_q]}
 {2(n+a)\mathbb E_\pi R^a}.
 \tag{10.4}
\]

Then (10.3) becomes

\[
 \mathbb E_\pi[R^aM]
 =2(n+a)\mathbb E_\pi R^a\,\mathfrak r_a(\pi).
 \tag{10.5}
\]

For iid pairs, \(q=1\) and \(\mathfrak r_a=1\) at every scale.  In the common-sign Laplace limit, \(\mathfrak r_0=0\): the distributional score at the sign boundary cancels all \(2n\) units of iid curvature.  For the synchronous law \(\pi_t\), (10.2) with \(w=1\) says

\[
 A(t)=2n+\mathbb E_{\pi_t}[R_tS_{q_t}],
\]

so (7.3) is equivalently the exact score-generation identity

\[
 \left.\frac d{dt}\right|_{t=0}
 \mathbb E_{\pi_t}[R_tS_{q_t}]
 =-4\mathbb E\Delta V.
\]

It measures how rapidly the radial score defect is generated from \(q_0=1\).

There is a useful quantitative sufficient condition.  Define the radial pair Fisher information

\[
 I_{\mathrm{rad}}(\pi)
 :=\mathbb E_\pi S_q^2.
 \tag{10.6}
\]

Cauchy-Schwarz in (10.2) with \(w=1\) gives

\[
 \mathbb E_\pi M
 \ge2n-
 \sqrt{\mathbb E_\pi R^2\,I_{\mathrm{rad}}(\pi)}.
 \tag{10.7}
\]

If both marginals equal a centered law with covariance at most \(I\), then

\[
 \mathbb E_\pi R^2\le4n.
\]

Consequently, if for some \(0\le\alpha<1\),

\[
 I_{\mathrm{rad}}(\pi)\le\alpha n,
 \tag{10.8}
\]

then

\[
 \mathbb E_\pi M
 \ge2(1-\sqrt\alpha)n
 \ge\frac{1-\sqrt\alpha}{2}\mathbb E_\pi R^2.
 \tag{10.9}
\]

Inserted in (7.1), (10.9) gives the genuine dimension-free contraction

\[
 D'(t)\le-(1-\sqrt\alpha)D(t).
 \tag{10.10}
\]

This is a concrete replacement for the false inference from marginal covariance to integrated curvature.  What must be controlled is the relative radial score generated by the coupling, not the mean Hessian and not ordinary mutual information.

The weighted version is

\[
 \left|\mathbb E_\pi[Rw(R)S_q]\right|
 \le
 \big(\mathbb E_\pi[R^2w(R)]\big)^{1/2}
 \big(\mathbb E_\pi[w(R)S_q^2]\big)^{1/2}
 \tag{10.11}
\]

for nonnegative \(w\).  Hence one can allow large score in thin distance bands and seek a dyadic budget for

\[
 I_w(\pi):=\mathbb E_\pi[w(R)S_q^2].
\]

This is the promised multiscale curvature functional.  A plausible coupling target is a time-integrated bound on the *negative* score terms in (10.2), separately on dyadic \(R\)-scales.  Such a bound is strictly more informative than pair entropy: Section 8 has bounded entropy but a score concentrating to infinity at the sign interface.

Synchronous coupling is poorly suited to controlling (10.6), because its relative coordinate receives no noise.  Reflection coupling does diffuse the radial relative coordinate and is therefore the natural mechanism for dissipating radial score.  The next section records the competing cost of doing so.

## 11. Reflection coupling: radial-score dissipation versus quadratic noise

Before coupling, reflection coupling has radial separation satisfying

\[
 dR_t=2\sqrt2\,d\beta_t-\frac{M(X_t,Y_t)}{R_t}\,dt.
 \tag{11.1}
\]

Thus, for smooth \(f\),

\[
 \frac d{dt}\mathbb Ef(R_t)
 =4\mathbb Ef''(R_t)
 -\mathbb E\left[f'(R_t)\frac{M(X_t,Y_t)}{R_t}\right].
 \tag{11.2}
\]

For a concave increasing \(f\), both displayed terms are nonpositive, but convexity alone gives no strict rate.  For the quadratic cost,

\[
 \frac d{dt}\mathbb ER_t^2=8-2\mathbb EM(X_t,Y_t).
 \tag{11.3}
\]

If the initial pair is iid, (2.3) gives the initial derivative \(8-4n\).  In dimension one it is positive.  Reflection noise therefore creates precisely the radial regularization missing from synchronous coupling, but it also injects a positive term into quadratic transport.

This suggests a genuinely hybrid target:

* use a linear or concave distance cost where reflection is active, so the \(4f''\) term is harmless;
* use synchronous evolution where a quadratic local cost is needed;
* track the dyadic radial score \(I_w\) from (10.11) to decide when sufficient chord curvature has survived for the synchronous phase.

The quadratic-linear cost (5.3) is compatible with this division: its tail is linear, so reflection creates no second-derivative penalty there, while its local quadratic part retains the second variation needed for Poincare.  The remaining hard theorem is a scale-resolved budget preventing the negative score term in (10.2) from cancelling the iid contribution at every scale.  Unlike a blanket integrated-gradient estimate, that theorem is not being assumed here; (10.2), (10.7), and the counterexamples above specify exactly what it must control and why entropy alone cannot control it.

## 12. Sharp obstruction summary

The proposed implications fail at distinct, noninterchangeable points:

\[
 \operatorname{Cov}(\mu)\preceq I
 \Longrightarrow
 \mathbb EH\succeq I
\]

is true, but

\[
 \mathbb EH\succeq I
 \centernot\Longrightarrow
 H(x)\succeq cI
\]

even on most of the mass (Section 3), and

\[
 \mathbb E_{\mu^{\otimes2}}M=2n
 \centernot\Longrightarrow
 \mathbb E_\pi M\gtrsim\mathbb E_\pi R^2
\]

for a correlated pair with the same marginals and only \(\log2\) mutual information (Section 8).  Synchronous coupling creates exactly such correlation and can erase the initial curvature at an arbitrarily large rate (Section 7).  Reflection coupling regularizes the missing radial score but adds quadratic separation noise (Section 11).  Finally, a global \(T_2\) target and an arbitrary-point integrated contraction target are both false for smooth isotropic measures (Sections 5 and 6).

The viable object left by this analysis is the scale-resolved retention coefficient \(\mathfrak r_a\), or equivalently the radial pair-score defect in (10.2).  Any coupling proof along this line must establish a dimension-free time/scale budget for that defect.  Mean curvature, marginal covariance, and pair entropy by themselves provably do not provide such a budget.
