# Random one-dimensional fibers do not give a dimension-descent recurrence

## 1. Scope and conclusion

Let \(\mu(dx)=Z^{-1}e^{-V(x)}dx\) be a smooth, full-dimensional,
isotropic log-concave probability measure on \(\mathbb R^n\).  For a unit
vector \(u\), write
\[
  X=Y+Tu,\qquad Y=P_{u^\perp}X,
\]
and denote the conditional law of \(T\), given \(Y=y\), by \(\nu_{u,y}\).
Put
\[
  m_u(y)=\mathbb E[T\mid Y=y],\qquad
  \sigma_u(y)^2=\operatorname{Var}(T\mid Y=y),
\]
and, for a locally Lipschitz \(f\),
\[
  g_u(y)=\mathbb E[f(X)\mid Y=y],
  \qquad
  Q_u(f)=\int \sigma_u(y)^2
       \mathbb E[(\partial_u f)^2\mid Y=y],d\bar\mu_u(y).
\]
Here \(\bar\mu_u=(P_{u^\perp})_*\mu\), which is itself isotropic on
\(u^\perp\).

The random-fiber proposal fails for two independent, explicit reasons.

1. On the isotropic cube, simplex, and shifted-exponential product,
   Haar-random fibers see only
   \[
      \mathbb E_u Q_u(f)=\Theta(n^{-2})\int |\nabla f|^2d\mu
   \]
   for linear modes.  On the cube the same upper scale holds for every
   vector in the genuine first Neumann eigenspace, including the normalized
   tensor-sum mode.  Thus a random fiber does not remove the required
   \(c/n\) fraction of a low mode.

2. Conditional expectation is not a contraction of transverse Dirichlet
   energy.  On the cube, for every linear \(f_a(x)=\langle a,x\rangle\),
   \[
     \mathbb E_u\int |\nabla g_u|^2d\bar\mu_u
       \ge \left(1-\frac1n+\frac{c}{\log n}\right)|a|^2.                 \tag{1.1}
   \]
   The excess is caused by the motion of the midpoint of a short chord.
   It occurs with \(\nabla^2f_a=0\), so the proposed Bochner Hessian budget
   does not control it.

Consequently these two terms cannot yield either
\[
 K_n\le (1-c/n)K_{n-1}+C/n
\]
or a fixed-factor dimension descent.  Pointwise inverse-chord reweighting
is nonintegrable on the cube and simplex.  Reweighting only the direction by
the inverse mean chord variance cancels upon normalization and leaves the
same \(n^{-2}\) scale.

All polytope calculations below can be transferred to smooth strictly
convex approximations by Hausdorff convergence, dominated chord integration,
and affine re-isotropization.  The counterexample therefore is not an
artifact of nonsmooth boundaries.

## 2. Exact disintegration identities

For every fixed \(u\), the law of total variance gives
\[
 \operatorname{Var}_\mu f
 =\int\operatorname{Var}_{\nu_{u,y}}(f(y+Tu)),d\bar\mu_u(y)
   +\operatorname{Var}_{\bar\mu_u}(g_u).                              \tag{2.1}
\]
The Euclidean energy splits as
\[
 \int|\nabla f|^2d\mu
 =\int(\partial_u f)^2d\mu+
   \int|P_{u^\perp}\nabla f|^2d\mu.                                  \tag{2.2}
\]
Therefore Haar averaging gives
\[
 \mathbb E_u\int(\partial_u f)^2d\mu={1\over n}\int|\nabla f|^2d\mu,
 \qquad
 \mathbb E_u\int|P_{u^\perp}\nabla f|^2d\mu
 ={n-1\over n}\int|\nabla f|^2d\mu.                                \tag{2.3}
\]

Every one-dimensional log-concave probability law \(\nu\) satisfies
\[
 \operatorname{Var}_\nu h
 \le C_1\operatorname{Var}_\nu(T)\int |h'|^2d\nu                    \tag{2.4}
\]
with a universal \(C_1\).  Applying (2.4) conditionally yields
\[
 \int\operatorname{Var}_{\nu_{u,y}}f,d\bar\mu_u(y)
 \le C_1 Q_u(f).                                                       \tag{2.5}
\]

The other term has an exact score correction.  If
\(q_u(y)=\int e^{-V(y+tu)}dt/Z\), differentiation under the integral gives,
for \(v\in u^\perp\),
\[
 \partial_v g_u(y)
 =\mathbb E[\partial_vf\mid Y=y]
  -\operatorname{Cov}_{\nu_{u,y}}(f,\partial_vV).                     \tag{2.6}
\]
Equivalently,
\[
 \nabla g_u
 =\mathbb E[P_{u^\perp}\nabla f\mid Y]
  -\operatorname{Cov}(f,P_{u^\perp}\nabla V\mid Y).                 \tag{2.7}
\]
Thus Jensen's inequality does not imply a transverse energy contraction.
Controlling the second term requires control of the conditional score, which
is absent for general convex \(V\) and becomes a boundary term for uniform
measures on convex bodies.

If one nevertheless had, uniformly in \(n,\mu,f\),
\[
 \mathbb E_u Q_u(f)\le {C\over n}\int|\nabla f|^2d\mu,
 \qquad
 \mathbb E_u\int|\nabla g_u|^2d\bar\mu_u
 \le(1-c/n)\int|\nabla f|^2d\mu,                                    \tag{2.8}
\]
then (2.1), (2.5), and the definition of \(K_{n-1}\) would give the desired
recurrence.  Sections 4 and 5 show that the mechanism required by (2.8)
fails: the first term can be too small to account for a \(c/n\) loss, and
the second term can be larger than the original energy.

## 3. What the Bochner budget does and does not say

For an exact eigenfunction of
\(L=\Delta-\langle\nabla V,\nabla\cdot\rangle\),
\[
  -Lf=\lambda f,
  \qquad
  \int|\nabla f|^2d\mu=\lambda\operatorname{Var}_\mu f,
\]
the integrated Bochner identity is
\[
 \int (Lf)^2d\mu
 =\int\|\nabla^2f\|_{HS}^2d\mu
  +\int\langle\nabla^2V\nabla f,\nabla f\rangle d\mu.
\]
Convexity hence gives
\[
 \int\|\nabla^2f\|_{HS}^2d\mu
 \le\lambda\int|\nabla f|^2d\mu.                                   \tag{3.1}
\]
But (3.1) controls variation of \(\nabla f\), not the conditional score in
(2.7).  The cube example below makes the latter arbitrarily large even for
a linear \(f\), for which the left side of (3.1) is zero.  Thus no estimate
that uses only the Rayleigh identity and (3.1) can control (2.7).

## 4. The isotropic cube: exact fiber scale

Let
\[
 K=[-\sqrt3,\sqrt3]^n,
 \qquad \mu=|K|^{-1}{\bf1}_Kdx.
\]
This measure is isotropic.  For a direction \(u\), let \(\ell_u(y)\) be the
length of the chord \((y+\mathbb Ru)\cap K\).  Conditional laws are uniform
on their chords, and hence
\[
 r(u):=\mathbb E\sigma_u(Y)^2
 ={1\over12|K|}\int_{u^\perp}\ell_u(y)^3dy.                           \tag{4.1}
\]
For every convex body,
\[
 \int_0^\infty t\,|K\cap(K-tu)|dt
 ={1\over6}\int_{u^\perp}\ell_u(y)^3dy.                              \tag{4.2}
\]
Indeed, on a chord of length \(\ell\), the integrand is
\(t(\ell-t)_+\), whose integral is \(\ell^3/6\).
Since
\[
 |K\cap(K-tu)|=\prod_{i=1}^n(2\sqrt3-t|u_i|)_+,
\]
(4.1)--(4.2), followed by \(t=2\sqrt3z\), give the exact formula
\[
 r(u)=6\int_0^{1/\|u\|_\infty}
 z\prod_{i=1}^n(1-z|u_i|)\,dz.                                      \tag{4.3}
\]
Put \(L=\|u\|_1\).  The inequalities
\[
 \prod_i(1-z|u_i|)\le e^{-zL}
\]
and, for \(0\le z\le(2L)^{-1}\),
\[
 \prod_i(1-z|u_i|)\ge e^{-2zL}
\]
imply universal constants \(c_0,C_0>0\) such that
\[
 {c_0\over\|u\|_1^2}\le r(u)\le {C_0\over\|u\|_1^2}.              \tag{4.4}
\]
For Haar \(u\), writing \(u=G/|G|\) with independent standard Gaussians,
standard Gaussian concentration gives
\[
 \mathbb E_u\|u\|_1^{-2}=\Theta(n^{-1}).                             \tag{4.5}
\]
The lower bound follows already from \(\|u\|_1^2\le n\).  For the upper
bound, with probability \(1-e^{-cn}\), at least a fixed fraction of the
\(|G_i|\)'s exceed \(1/2\) and \(|G|^2\le2n\); on that event
\(|G|^2/(\sum_i|G_i|)^2\le C/n\), and the complement is harmless.

For \(f_a(x)=\langle a,x\rangle\),
\[
 Q_u(f_a)=r(u)\langle a,u\rangle^2.
\]
Hyperoctahedral symmetry forces
\[
 \mathbb E_u[r(u)uu^T]={\mathbb E_ur(u)\over n}I,
\]
so (4.4)--(4.5) give
\[
 \boxed{\quad
 \mathbb E_uQ_u(f_a)=\Theta(n^{-2})|a|^2
 =\Theta(n^{-2})\operatorname{Var}_\mu(f_a).
 \quad}                                                               \tag{4.6}
\]
Here \(\nabla^2f_a=0\).

This persists for genuine first modes.  In one coordinate let
\[
 h(t)=\sqrt2\sin\!\left({\pi t\over2\sqrt3}\right).
\]
Then \(\operatorname{Var}h(X_1)=1\), its Neumann eigenvalue is
\(\lambda_1=\pi^2/12\), and
\[
 \int(h'')^2d\mu=\lambda_1\int(h')^2d\mu.
\]
Moreover, using \(\|h'\|_\infty^2=\pi^2/6\),
\[
 \mathbb E_uQ_u(h(x_1))
 \le {\pi^2\over6}\mathbb E_u[u_1^2r(u)]=O(n^{-2}).                 \tag{4.7}
\]
The quadratic form \(f\mapsto\mathbb E_uQ_u(f)\) is invariant under signed
coordinate permutations.  Its restriction to the irreducible first
eigenspace \(\operatorname{span}\{h(x_i):1\le i\le n\}\) is therefore a
scalar multiple of the identity.  Thus (4.7) also holds for the normalized
tensor sum
\[
 n^{-1/2}\sum_{i=1}^nh(x_i).                                         \tag{4.8}
\]
This is an exact low-Hessian first eigenmode, not merely a linear test.

## 5. The cube conditional-mean energy grows

The failure of transverse energy contraction can also be computed exactly.
Fix a generic \(u\), and let
\[
 m(y)=\mathbb E[T\mid Y=y]
\]
be the midpoint of the cube chord.  On every cell on which the lower and
upper endpoints lie on fixed facets \(i\) and \(j\), endpoint differentiation
in \(u^\perp\) gives
\[
 \nabla t_-=-{P_{u^\perp}e_i\over u_i},
 \qquad
 \nabla t_+=-{P_{u^\perp}e_j\over u_j}.                              \tag{5.1}
\]
If \(i\ne j\), then
\[
 |\nabla m|^2
 ={1\over4}\left({1\over u_i^2}+{1\over u_j^2}-4\right)
 \ge {1\over2\|u\|_\infty^2}-1.                                   \tag{5.2}
\]
If \(i=j\), the value is \(u_i^{-2}-1\), and the same lower bound holds.
Ties occur only on a null set.

For Haar \(u\), with probability at least \(1/2\),
\[
 \|u\|_\infty^2\le {8\log n\over n}.
\]
Consequently, for all sufficiently large \(n\),
\[
 \mathbb E_u\int|\nabla m|^2d\bar\mu_u
 \ge {c n\over\log n}.                                               \tag{5.3}
\]

Now let \(b_u(y)=\mathbb E[X\mid Y=y]=y+m(y)u\).  For
\(f_a(x)=\langle a,x\rangle\),
\[
 g_{u,a}(y)=\langle a,b_u(y)\rangle,
 \qquad
 \nabla g_{u,a}=P_{u^\perp}a+\langle a,u\rangle\nabla m.             \tag{5.4}
\]
The Hilbert--Schmidt identity
\[
 \sum_{k=1}^n|\nabla g_{u,e_k}|^2=n-1+|\nabla m|^2                    \tag{5.5}
\]
is pointwise: the differential of \(b_u\) maps
\(v\in u^\perp\) to \(v+u\langle\nabla m,v\rangle\), and the two summands
are orthogonal.  After averaging in \(y,u\), hyperoctahedral symmetry makes
the resulting quadratic form in \(a\) a scalar.  Taking its trace in
(5.5) and using (5.3) proves
\[
 \boxed{\quad
 \mathbb E_u\int|\nabla g_{u,a}|^2d\bar\mu_u
 \ge\left(1-{1\over n}+{c\over\log n}\right)|a|^2.
 \quad}                                                               \tag{5.6}
\]
Thus the conditional score in (2.7) does not merely cancel the formal
\(1/n\) transverse loss: it creates an energy excess much larger than
\(1/n\).  Since \(\nabla^2f_a=0\) and its Rayleigh quotient is \(1\), this
also shows that the Rayleigh identity plus the Bochner Hessian inequality
alone cannot imply the required contraction.  The linear quotient is within
the universal factor \(12/\pi^2\) of the true cube spectral quotient.

## 6. Simplex and shifted exponentials

### 6.1 Isotropic simplex

Let \(N=n+1\),
\[
 \Delta=\{z_i\ge0:\sum_{i=1}^Nz_i=1\},\qquad
 H=\{x:\sum_i x_i=0\},
\]
and
\[
 X=\sqrt{N(N+1)}\left(Z-N^{-1}{\bf1}\right),
 \qquad Z\sim\mathrm{Unif}(\Delta).
\]
Then \(X\) is isotropic on the \(n\)-dimensional space \(H\).  For
\(u\in H\), \(|u|=1\), put
\[
 A=\sum_{u_i>0}u_i=\sum_{u_i<0}|u_i|={1\over2}\|u\|_1.
\]
At a sampled point the unscaled distances to the two endpoints satisfy
\[
 \mathbb P(T_->s,T_+>t)=(1-A(s+t))_+^n.                               \tag{6.1}
\]
Integrating (6.1) gives
\[
 \mathbb ET_-^2=\mathbb ET_+^2={2\over A^2(n+1)(n+2)},
 \qquad
 \mathbb E(T_-T_+)={1\over A^2(n+1)(n+2)}.
\]
After the isotropic dilation, the squared chord length \(D^2\) has
\(\mathbb ED^2=6/A^2\).  Conditional fibers are uniform, hence
\[
 r(u)=\mathbb E\sigma_u(Y)^2={\mathbb ED^2\over12}
 ={1\over2A^2}={2\over\|u\|_1^2}.                                   \tag{6.2}
\]
Since \(\mathbb E_u\|u\|_1^{-2}=\Theta(n^{-1})\) on the sphere of \(H\),
irreducibility of the permutation action gives, for every linear mode,
\[
 \mathbb E_uQ_u(f_a)=\Theta(n^{-2})|a|^2.                             \tag{6.3}
\]
The midpoint is again piecewise affine with derivatives inverse to active
facet components, so the same moving-endpoint score obstruction is present.

### 6.2 Product of shifted exponentials

Let \(X_i=E_i-1\), with independent \(E_i\sim\mathrm{Exp}(1)\).  This law
is isotropic.  Put
\[
 A_+=\sum_{u_i>0}u_i,\qquad A_-=\sum_{u_i<0}|u_i|,
 \qquad a=\max(A_+,A_-).
\]
The two distances from a sampled point to the fiber endpoints are independent
exponentials of rates \(A_+\) and \(A_-\) (with an infinite endpoint allowed
when one rate is zero).  The conditional density along the interval is a
truncated exponential.  Direct integration gives universal constants
\[
 {c\over a^2}\le r(u)=\mathbb E\sigma_u(Y)^2\le {C\over a^2}.         \tag{6.4}
\]
Because \(\|u\|_1/2\le a\le\|u\|_1\),
\[
 r(u)=\Theta(\|u\|_1^{-2}),
 \qquad \mathbb E_ur(u)=\Theta(n^{-1}).                              \tag{6.5}
\]
Permutation symmetry makes the diagonal entries of
\(\mathbb E[r(u)uu^T]\) equal.  Hence each coordinate linear mode again has
\[
 \mathbb E_uQ_u(x_i)=\Theta(n^{-2}).                                 \tag{6.6}
\]
The normalized tensor-sum linear mode has the same order, by Gaussian
concentration of \(\sum_i u_i\) together with (6.5).

### 6.3 Radial control case

For an isotropic radial law, including an isotropically scaled density
proportional to \(e^{-a_n|x|^4}\), the conditional law of \(T\) given \(Y\)
is even.  Thus \(m_u(Y)=0\) and
\[
 \mathbb E\sigma_u(Y)^2=\operatorname{Var}\langle X,u\rangle=1.
\]
For linear \(f_a\),
\[
 \mathbb E_uQ_u(f_a)={|a|^2\over n},
 \qquad
 \mathbb E_u\int|\nabla g_u|^2d\bar\mu_u={n-1\over n}|a|^2.        \tag{6.7}
\]
The radial case is therefore ideal for the proposed descent.  Comparing
(6.7) with (4.6), (5.6), (6.3), and (6.6) shows that isotropy and a Hessian
budget do not distinguish the favorable radial geometry from short,
moving polytope chords.

## 7. Why inverse-chord reweighting does not repair the argument

There are two natural inverse-variance modifications, and both fail.

### 7.1 Pointwise inverse conditional variance is nonintegrable

For a uniform convex body, \(\sigma_u(y)^2=\ell_u(y)^2/12\), while the
marginal density of \(Y\) is \(\ell_u(y)/|K|\).  Therefore
\[
 \mathbb E_{\bar\mu_u}[\sigma_u(Y)^{-2}]
 ={12\over|K|}\int_{u^\perp}{dy\over\ell_u(y)}.                      \tag{7.1}
\]
For a generic direction in a polytope, the chord length vanishes linearly
with distance to an open piece of the boundary of the projection.  The
integral in (7.1) consequently contains \(\int_0^\varepsilon ds/s\) and is
infinite.  This applies to both the cube and the simplex.  For mixed-sign
directions in the shifted orthant, short truncated-exponential intervals
have variance comparable to \(\ell^2\), and the same logarithmic divergence
occurs.  Truncating the weight at \(\sigma^{-2}\le R\) introduces an explicit
\(\log R\) loss and no universal limiting estimate.

There is also a structural cancellation.  In the standardized fiber
coordinate \(z=(t-m_u(y))/\sigma_u(y)\), differentiation is
\(\partial_z=\sigma_u(y)\partial_t\).  Thus the one-dimensional Poincare
energy is still \(\sigma_u(y)^2(\partial_t f)^2\); a change of fiber units
does not remove the factor appearing in \(Q_u\).

### 7.2 Inverse mean variance cancels after normalization

Suppose directions are sampled with density proportional to \(r(u)^{-1}\).
On the cube, simplex, and exponential product,
\(r(u)=\Theta(\|u\|_1^{-2})\) and
\[
 \mathbb E_u r(u)^{-1}=\Theta(n).                                    \tag{7.2}
\]
For a coordinate linear mode, normalized inverse weighting gives
\[
 {\mathbb E_u[r(u)^{-1}Q_u(f_{e_i})]\over
   \mathbb E_u r(u)^{-1}}
 ={\mathbb E_u u_i^2\over\mathbb E_u r(u)^{-1}}
 =\Theta(n^{-2}).                                                     \tag{7.3}
\]
Thus the apparent factor \(r^{-1}\) is exactly lost in the normalization.
Any fixed polynomial weighting by \(r\) or \(r^{-1}\) remains concentrated
on the Haar-typical regime \(\|u\|_1\asymp\sqrt n\); reaching the sparse
coordinate directions of the cube requires an exponent growing with \(n\),
which is not a dimension-free mechanism.

Selecting a special basis of high-variance directions does work for the
cube (its coordinate axes), but no rotation-equivariant covariance datum
identifies that basis: the covariance is the identity.  A universal theorem
that simultaneously selects such directions and controls the score term
(2.7) would already provide the missing isoperimetric structure.  Treating
that selection as a harmless reweighting is therefore a hidden replacement
of KLS, not a proof of it.

## 8. Formal recurrence refutation

Combining (2.1), (2.5), and the \((n-1)\)-dimensional Poincare inequality
gives the only direct fiber recurrence:
\[
 \operatorname{Var}_\mu f
 \le C_1Q_u(f)+K_{n-1}\int|\nabla g_u|^2d\bar\mu_u.                   \tag{8.1}
\]
Averaging (8.1) does not improve it.  On the cube, for a unit linear mode,
the two right-hand geometric inputs obey
\[
 \mathbb E_uQ_u(f)=\Theta(n^{-2}),
 \qquad
 \mathbb E_u\int|\nabla g_u|^2d\bar\mu_u
 \ge1-{1\over n}+{c\over\log n}.                                   \tag{8.2}
\]
The Hessian budget is zero.  Hence no constants \(c,C>0\), independent of
\(n\), can turn the second term of (8.1) into
\((1-c/n)K_{n-1}\int|\nabla f|^2\), nor into a fixed-factor contraction.
The genuine cube first modes satisfy (3.1) with equality and have the same
\(O(n^{-2})\) fiber scale by (4.7)--(4.8), so imposing exact low-mode status
does not restore the missing fiber variance.

The obstruction is geometric and precise: Haar-random lines through a
faceted isotropic body have conditional variance of order \(1/n\), and the
projection of a fixed gradient onto such a line costs another factor
\(1/n\).  At the same time, the short chord's midpoint moves at speed
\(\Theta(\sqrt{n/\log n})\), which reappears as the conditional-score term.
Any viable dimension descent must use higher-dimensional fibers, an
extremality-dependent choice of directions, or a new control of the moving
conditional law; the low-Hessian condition by itself supplies none of these.
