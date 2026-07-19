# Random-line uncertainty: exact identities, sharp models, and strength

## 0. Verdict

For an isotropic log-concave probability \(\mu\) and a balanced Borel set
\(E\), define

\[
 \mathcal U_\mu(E)=
 \mathbb E_{\theta\sim{\rm Unif}(S^{n-1})}
 \int_{\theta^\perp}
 { \min(p_{\theta,y},1-p_{\theta,y})\over
   \sigma_{\theta,y}}\,
 d(\pi_{\theta^\perp}\mu)(y).                              \tag{0.1}
\]

Here \(p_{\theta,y}\) is the conditional mass of \(E\) on the line
\(y+\mathbb R\theta\), and \(\sigma_{\theta,y}\) is the conditional
standard deviation of its line coordinate.

The normalization in the proposed estimate

\[
                         \mathcal U_\mu(E)\ge {c\over\sqrt n} \tag{0.2}
\]

is exactly right. Slicing, one-dimensional log-concave isoperimetry, and
Crofton give

\[
 {1\over2\sqrt3}\mathcal U_\mu(E)
 \le a_n P_\mu(E),\qquad
 a_n=\mathbb E|\theta_1|
 ={ \Gamma(n/2)\over\sqrt\pi\,\Gamma((n+1)/2)}
 \asymp n^{-1/2}.                                         \tag{0.3}
\]

Thus (0.2) implies a dimension-free balanced Cheeger bound, hence KLS.
It is not a weaker reformulation obtained for free from Crofton: (0.3)
has the opposite direction. The random-line statement asserts additional
nonlocal uncertainty of every balanced set.

All mandatory models pass:

- a Gaussian halfspace has the exact value
  \[
  \mathcal U_\gamma(E)
  ={1\over\pi}\mathbb E_\theta\arcsin|\theta_1|
  \asymp n^{-1/2};
  \]
- a central coordinate cut of the isotropic cube satisfies
  \[
  {a_n\over48}\le\mathcal U(E)\le {a_n\over2};
  \]
- a balanced regular-simplex cap has
  \(\mathcal U(E)\asymp n^{-1/2}\);
- the median sphere of the isotropic radial exponential law has
  \(\mathcal U(E)\asymp n^{-1/2}\);
- the product-exponential maximum has
  \(\mathcal U(E)\asymp n^{-1/2}\), while a parity checkerboard has
  \(\mathcal U(E)\asymp\sqrt n\).

Halfspaces and smooth radial cuts are therefore sharp; many-facet or
oscillatory sets are not obstructions.

The exact Radon identity records a constant amount of cross-phase pair
mass, but with the kernel \(|x-z|^{n-1}\). Controlling that kernel by the
conditional variance uses \((n-1)\)-st moments and loses powers of \(n\);
it does not prove (0.2). No counterexample to (0.2) is produced here, but
the audit shows that it is a KLS-strength new theorem, not a completed
reduction. The missing assertion is a dimension-free lower bound for the
conductance of a variance-normalized hit-and-run form.

## 1. Conditional line laws and normalization

First assume that \(\mu\) is full-dimensional with density
\(\rho=e^{-V}\), \(V\) convex. Fix \(\theta\in S^{n-1}\). For
\(y\in\theta^\perp\), put

\[
\begin{aligned}
 Z_\theta(y)&=\int_{\mathbb R}\rho(y+t\theta)\,dt,\\
 f_{\theta,y}(t)&={\rho(y+t\theta)\over Z_\theta(y)},\\
 m_{\theta,y}&=\int t f_{\theta,y}(t)\,dt,\\
 \sigma_{\theta,y}^2
  &=\int(t-m_{\theta,y})^2f_{\theta,y}(t)\,dt.              \tag{1.1}
\end{aligned}
\]

The projection law has density \(Z_\theta(y)dy\), and
\(f_{\theta,y}\) is log-concave. For projection-almost every \(y\),
\(\sigma_{\theta,y}>0\). Define

\[
 p_{\theta,y}=\int_{\{t:y+t\theta\in E\}}
                    f_{\theta,y}(t)\,dt,\qquad
 q_{\theta,y}=\min(p_{\theta,y},1-p_{\theta,y}).            \tag{1.2}
\]

Then (0.1) is

\[
 \mathcal U_\mu(E)
 =\mathbb E_\theta\int_{\theta^\perp}
 {q_{\theta,y}\over\sigma_{\theta,y}}
 Z_\theta(y)\,dy.                                          \tag{1.3}
\]

Since \(\mu\) is isotropic, conditional variance decomposition gives, for
every fixed \(\theta\),

\[
 \int_{\theta^\perp}\sigma_{\theta,y}^2
 Z_\theta(y)\,dy
 \le\operatorname{Var}_\mu\langle X,\theta\rangle=1.        \tag{1.4}
\]

This explains why division by \(\sigma_{\theta,y}\) compensates for the
very short typical chords of bodies such as the cube.

### Same-line resampling

Given \((\theta,y)\), let \(T,T'\) be conditionally independent with
density \(f_{\theta,y}\). Then

\[
 \mathbb P\bigl(
 \mathbf1_E(y+T\theta)\ne\mathbf1_E(y+T'\theta)
 \mid\theta,y\bigr)=2p_{\theta,y}(1-p_{\theta,y}).          \tag{1.5}
\]

For \(0\le p\le1\),

\[
 {q\over2}\le p(1-p)\le q.                                 \tag{1.6}
\]

Consequently the rate-\(1/\sigma_{\theta,y}\) same-line flow

\[
 \mathcal F_\mu(E)=
 \mathbb E_\theta\int
 {2p_{\theta,y}(1-p_{\theta,y})\over\sigma_{\theta,y}}
 d(\pi_{\theta^\perp}\mu)(y)                               \tag{1.7}
\]

satisfies

\[
                         \mathcal U_\mu(E)
 \le\mathcal F_\mu(E)\le2\mathcal U_\mu(E).                 \tag{1.8}
\]

Thus (0.2) is precisely a balanced conductance bound for a
variance-normalized hit-and-run form.

## 2. Slicing, one-dimensional Cheeger, and Crofton

Let \(E\) have locally finite perimeter. Its weighted directional
variation is

\[
 P_{\mu,\theta}(E)
 =\int_{\partial^*E}\rho(x)|\nu_E(x)\cdot\theta|\,
                         d\mathcal H^{n-1}(x).              \tag{2.1}
\]

The BV slicing theorem and disintegration give the exact identity

\[
 P_{\mu,\theta}(E)=
 \int_{\theta^\perp}
 P_{\mu_{\theta,y}}\bigl(E\cap(y+\mathbb R\theta)\bigr)
 d(\pi_{\theta^\perp}\mu)(y).                              \tag{2.2}
\]

A one-dimensional log-concave probability of standard deviation
\(\sigma\) has Cheeger constant at least
\(1/(2\sqrt3\,\sigma)\). Equivalently, for every Borel subset \(A\) of
the line,

\[
 P_\nu(A)\ge
 {1\over2\sqrt3\,\sigma}
 \min(\nu(A),1-\nu(A)).                                    \tag{2.3}
\]

This follows from one-dimensional log-concave isoperimetry, whose
minimizers are half-lines, together with the isotropic quantile-density
bound \(f(F^{-1}(s))\ge\min(s,1-s)/(2\sqrt3)\).
Applying (2.3) in (2.2) yields

\[
                         P_{\mu,\theta}(E)
 \ge {1\over2\sqrt3}
 \int {q_{\theta,y}\over\sigma_{\theta,y}}
 d(\pi_{\theta^\perp}\mu)(y).                              \tag{2.4}
\]

Finally,

\[
\begin{aligned}
 \mathbb E_\theta P_{\mu,\theta}(E)
 &=\int_{\partial^*E}\rho(x)
       \mathbb E_\theta|\nu_E(x)\cdot\theta|\,
       d\mathcal H^{n-1}(x)\\
 &=a_nP_\mu(E).                                             \tag{2.5}
\end{aligned}
\]

Equations (2.4)--(2.5) prove (0.3).

For a log-concave measure the Minkowski isoperimetric profile is concave
on \((0,1)\) and symmetric under \(s\mapsto1-s\) (first for positive
smooth densities, and then for arbitrary log-concave probabilities by
restriction to the affine support and monotone regularization). Hence
\(I_\mu(s)/s\) is nonincreasing on \((0,1/2]\), and

\[
                         \psi_\mu=2I_\mu(1/2).              \tag{2.6}
\]

Therefore the balanced estimate supplied by (0.2)--(0.3) gives the full
Cheeger inequality.

Here the relaxed weighted BV perimeter is no larger than exterior
Minkowski content. Thus sets of infinite Minkowski content are harmless,
and approximation by finite-perimeter sets transfers (2.4) to every
Borel set for which it is nontrivial. Lower-dimensional log-concave
measures are first identified isometrically with their affine support;
isotropy and the spherical average are then taken in that support
dimension.

## 3. Exact Radon and pair identities

For a fixed line, use the unnormalized masses

\[
\begin{aligned}
 M_E&=\int_{E\cap(y+\mathbb R\theta)}\rho\,dt,\\
 M_{E^c}&=Z-M_E,\\
 S_2&=\int(t-m)^2\rho(y+t\theta)\,dt.                       \tag{3.1}
\end{aligned}
\]

Then

\[
 qZ=\min(M_E,M_{E^c}),\qquad
 Z^2\sigma^2=ZS_2
 ={1\over2}\int_{\mathbb R^2}
       (s-t)^2\rho(y+s\theta)\rho(y+t\theta)\,ds\,dt.       \tag{3.2}
\]

Also,

\[
 M_EM_{E^c}
 =\int_{\substack{y+s\theta\in E\\y+t\theta\notin E}}
       \rho(y+s\theta)\rho(y+t\theta)\,ds\,dt.              \tag{3.3}
\]

The affine Blaschke--Petkantschin formula, with \(\theta\) uniform on the
whole sphere, is

\[
\begin{aligned}
 &\mathbb E_\theta\int_{\theta^\perp}\int_{\mathbb R^2}
 H(y+s\theta,y+t\theta)|s-t|^{n-1}\,ds\,dt\,dy\\
 &\hspace{35mm}
 ={2\over|S^{n-1}|}
 \int_{\mathbb R^n}\int_{\mathbb R^n}H(x,z)\,dx\,dz.        \tag{3.4}
\end{aligned}
\]

For \(H(x,z)=\rho(x)\rho(z)\mathbf1_E(x)\mathbf1_{E^c}(z)\)
and \(\mu(E)=1/2\), this gives

\[
 \mathbb E_\theta\int_{\theta^\perp}
 \int_{E_y}\int_{E_y^c}
 \rho(y+s\theta)\rho(y+t\theta)|s-t|^{n-1}
 \,ds\,dt\,dy
 ={1\over2|S^{n-1}|}.                                     \tag{3.5}
\]

Equation (3.5) is the exact Radon form of balance.

It does not close (0.2). A one-dimensional log-concave law satisfies

\[
 \bigl(\mathbb E|T-T'|^{n-1}\bigr)^{1/(n-1)}
 \le Cn\,\sigma.                                           \tag{3.6}
\]

The direct substitution of (3.6) into (3.5) introduces
\((Cn)^{n-1}\), while \(\mathcal U\) contains only one inverse power of
\(\sigma\). The usual Hölder comparison therefore loses at least a
linear factor before the \((n-1)\)-st root is taken. The identity is
exact, but its natural kernel is of the wrong order for this direct
argument; a genuinely different use of (3.5) would be needed.

## 4. Gaussian halfspace: exact Hermite/noise calculation

Let \(\gamma_n\) be standard Gaussian measure and
\(E=\{x_1\le0\}\). Gaussian orthogonal coordinates are independent, so

\[
                         \sigma_{\theta,y}=1.               \tag{4.1}
\]

Put \(a=|\theta_1|\). The first coordinate of \(y\) is Gaussian with
variance \(1-a^2\), and, conditionally on \(y\),

\[
 q_{\theta,y}
 =\Phi\left(-{|y_1|\over a}\right).                         \tag{4.2}
\]

If \(Z,W\) are independent standard Gaussians, rotational invariance in
the \((Z,W)\)-plane gives

\[
\begin{aligned}
 \int q_{\theta,y}\,d(\pi_{\theta^\perp}\gamma_n)(y)
 &=\mathbb P\left(
 W\le-{\sqrt{1-a^2}\over a}|Z|\right)\\
 &={\arcsin a\over\pi}.                                    \tag{4.3}
\end{aligned}
\]

Hence

\[
 \boxed{\quad
 \mathcal U_{\gamma_n}(E)
 ={1\over\pi}\mathbb E_\theta\arcsin|\theta_1|.
 \quad}                                                     \tag{4.4}
\]

Since \(a\le\arcsin a\le(\pi/2)a\),

\[
 {a_n\over\pi}
 \le\mathcal U_{\gamma_n}(E)
 \le {a_n\over2}.                                          \tag{4.5}
\]

The same-line noise calculation gives the complementary identity. If
\(X,X'\) are conditionally independent on the line, then
\(\operatorname{Corr}(X_1,X_1')=1-a^2\). The Hermite expansion of
\(\operatorname{sign}\), equivalently Sheppard's noise-stability formula,
gives

\[
 \mathbb P(\operatorname{sign}X_1\ne\operatorname{sign}X_1')
 ={ \arccos(1-a^2)\over\pi}.                               \tag{4.6}
\]

Equations (4.3) and (4.6), together with (1.6), agree up to the required
universal factors. The Gaussian halfspace has exactly the conjectured
\(n^{-1/2}\) scale.

## 5. Uniform convex bodies: an exact halfspace formula

Let \(\mu\) be uniform on a convex body \(K\), let
\(E=K\cap\{x:u\cdot x\le h\}\), and put

\[
 \Sigma=K\cap\{x:u\cdot x=h\}.                              \tag{5.1}
\]

For \(z\in\Sigma\), define the two chord reaches

\[
 \tau_\pm(z,\theta)=
 \sup\{t\ge0:z\pm t\theta\in K\}.                           \tag{5.2}
\]

The conditional law on each chord is uniform, with

\[
 \sigma={\tau_++\tau_-\over\sqrt{12}},\qquad
 q={\min(\tau_+,\tau_-)\over\tau_++\tau_-}.                 \tag{5.3}
\]

Projection from the cut plane to \(\theta^\perp\) has Jacobian
\(|u\cdot\theta|\). Therefore

\[
\boxed{
 \begin{aligned}
 \mathcal U_{\mu,\theta}(E)
 &=\sqrt{12}\,|u\cdot\theta|\,P_\mu(E)\,
 \mathbb E_{z\sim{\rm Unif}(\Sigma)}
 { \min(\tau_+,\tau_-)\over\tau_++\tau_-}.
 \end{aligned}}                                             \tag{5.4}
\]

This identity isolates the only issue for a polyhedral halfspace:
whether a typical cut point is macroscopically centered in its random
chord.

## 6. The cube coordinate halfspace

Take \(K=[-\sqrt3,\sqrt3]^n\) and
\(E=K\cap\{x_1\le0\}\). Here

\[
                           P_\mu(E)={1\over2\sqrt3}.          \tag{6.1}
\]

For fixed \(\theta\), let \(z\) be uniform on the central slice. Put
\[
 R(z,\theta)={\min(\tau_+,\tau_-)\over\tau_++\tau_-}.        \tag{6.2}
\]

### Lemma 6.1 (a cube cut point is centered with fixed probability)

For every \(\theta\) with \(\theta_1\ne0\),

\[
                       {1\over48}\le\mathbb E_zR(z,\theta)
                       \le{1\over2}.                        \tag{6.3}
\]

#### Proof

Write \(a=\sqrt3\). The \(x_1\)-constraint gives the common cap
\(d_1=a/|\theta_1|\). For \(j\ge2\), after orienting by the sign of
\(\theta_j\), the positive and negative distances contributed by the
\(j\)-th coordinate are

\[
                       U_j,\qquad c_j-U_j,\qquad
 c_j={2a\over|\theta_j|},                                  \tag{6.4}
\]

where \(U_j\) is uniform on \([0,c_j]\), independently over \(j\).
Thus

\[
 \tau_+=\min(d_1,U_2,\ldots,U_n),\qquad
 \tau_-=\min(d_1,c_2-U_2,\ldots,c_n-U_n).                  \tag{6.5}
\]

Choose \(t_0\) as follows. If
\(\prod_{j\ge2}(1-d_1/c_j)_+\ge1/2\), put \(t_0=d_1\).
Otherwise choose the unique \(t_0<d_1\) satisfying

\[
                         \prod_{j\ge2}(1-t_0/c_j)_+={1\over2}.\tag{6.6}
\]

The joint event \(\{\tau_+\ge t_0/2,\tau_-\ge t_0/2\}\) has probability
at least \(1/2\), because its \(j\)-th factor is
\(1-t_0/c_j\). If \(t_0<d_1\), Bernoulli's inequality gives

\[
 \mathbb P(\tau_+>3t_0)
 \le\prod_j(1-3t_0/c_j)_+
 \le\left[\prod_j(1-t_0/c_j)_+\right]^3
 ={1\over8},                                               \tag{6.7}
\]

and the same holds for \(\tau_-\). If \(t_0=d_1\), both tails vanish.
Hence, with probability at least \(1/4\),

\[
                    {t_0\over2}\le\tau_\pm\le3t_0.          \tag{6.8}
\]

On this event \(R\ge1/12\), proving the lower bound \(1/48\).
The upper bound is pointwise. QED.

Equations (5.4), (6.1), and \(\sqrt{12}/(2\sqrt3)=1\) give

\[
                  {|\theta_1|\over48}
 \le\mathcal U_{\mu,\theta}(E)
 \le{|\theta_1|\over2}.                                    \tag{6.9}
\]

Averaging proves

\[
                         {a_n\over48}
 \le\mathcal U_\mu(E)\le {a_n\over2}.                       \tag{6.10}
\]

The inverse conditional standard deviation is essential here. Typical
stationary hit-and-run chords of the cube are short; the factor
\(1/\sigma\) restores the \(n^{-1/2}\) boundary scale.

## 7. A balanced cap of the isotropic regular simplex

Let \(K\) be the regular simplex in isotropic position, with vertices
\(v_0,\ldots,v_n\), and let \(\lambda_i\) be its barycentric coordinates.
The normalization may be chosen so that

\[
 |v_i|^2=n(n+2),\qquad v_i\cdot v_j=-(n+2)\quad(i\ne j).    \tag{7.1}
\]

The altitude from \(v_0\) to its opposite facet is

\[
                         h_0=(n+1)\sqrt{{n+2\over n}}.       \tag{7.2}
\]

Take

\[
 E=\{\lambda_0\ge t_n\},\qquad
                         (1-t_n)^n={1\over2}.               \tag{7.3}
\]

Since \(\lambda_0\sim{\rm Beta}(1,n)\), this cap is balanced. Its weighted
slice area is the density of the unit normal coordinate:

\[
 P_\mu(E)
 ={n(1-t_n)^{n-1}\over h_0}.                               \tag{7.4}
\]

For every \(n\ge2\),

\[
                              {1\over5}\le P_\mu(E)\le1.     \tag{7.5}
\]

### Lemma 7.1 (simplex chord balance)

There are universal \(c,C>0\) such that

\[
 ca_n\le
 \mathbb E_{\theta,z}\left[
 |u\cdot\theta|\,
 { \min(\tau_+(z,\theta),\tau_-(z,\theta))\over
   \tau_+(z,\theta)+\tau_-(z,\theta)}\right]
 \le {a_n\over2},                                          \tag{7.6}
\]

where \(\theta\) is uniform and \(z\) is uniform on the cap slice.

#### Proof

On the slice, write

\[
 (\lambda_1,\ldots,\lambda_n)
 =(1-t_n){(E_1,\ldots,E_n)\over S},\qquad
 S=E_1+\cdots+E_n,                                         \tag{7.7}
\]

with independent unit exponentials \(E_i\). If
\(b_i=D_\theta\lambda_i\), then \(\sum_{i=0}^n b_i=0\). Put

\[
 R_+=\sum_{\substack{1\le i\le n\\b_i>0}}b_i,\qquad
 R_-=\sum_{\substack{1\le i\le n\\b_i<0}}|b_i|.            \tag{7.8}
\]

Before imposing the deterministic coordinate \(\lambda _0=t_n\), the
two reaches are

\[
 {1-t_n\over S}W_+,\qquad {1-t_n\over S}W_-,
 \quad
 W_+\sim{\rm Exp}(R_-),\quad W_-\sim{\rm Exp}(R_+),         \tag{7.9}
\]

and \(W_+,W_-\) are independent.  The coordinate \(0\) truncates the
appropriate one of these clocks at

\[
 D_0={t_nS\over(1-t_n)|b_0|}                               \tag{7.10}
\]

in the \(W\)-scale.

The regular-simplex frame is tight. Consequently, after a scalar
normalization, \(b=(b_0,\ldots,b_n)\) is uniform on the unit sphere of
\(\{s:\sum_i s_i=0\}\). In particular its coordinates are exchangeable.
Let \(L=\sum_{i=1}^n|b_i|\). Since at most four coordinates can exceed
\(\frac15\sum_{i=0}^n|b_i|\), exchangeability gives, for \(n\ge8\),

\[
 \mathbb P_\theta\{|b_0|\le L/4\}
 \ge1-{4\over n+1}\ge {5\over9}.                           \tag{7.11}
\]

On this event, \(R_+-R_-=-b_0\), and therefore

\[
 {3L\over8}\le R_\pm\le {5L\over8}.                        \tag{7.12}
\]

Also

\[
 {1\over2}\le {nt_n\over1-t_n}\le1.                        \tag{7.13}
\]

If \(S\ge n/2\), (7.10)--(7.13) imply that the deterministic clock is at
least a universal multiple of both mean race times:

\[
                         D_0R_\pm\ge {3\over8}.             \tag{7.14}
\]

Conditional on a direction satisfying (7.11), there is a universal
event of positive probability on which

\[
 {1\over20R_-}\le W_+\le {1\over10R_-},\qquad
 {1\over20R_+}\le W_-\le {1\over10R_+}.                    \tag{7.15}
\]

On (7.14)--(7.15), the deterministic clock does not bind, the two actual
reaches differ by at most \(10/3\), and hence their unweighted ratio in
(7.6) is at least \(3/13\). In addition, for a uniform spherical
direction,

\[
 \mathbb P\left\{
 {c_0\over\sqrt n}\le|u\cdot\theta|
 \le {C_0\over\sqrt n}\right\}\ge p_0                     \tag{7.16}
\]

with universal positive constants. For all sufficiently large \(n\),
subtracting the \(4/(n+1)\) exceptional probability in (7.11) leaves a
fixed-probability intersection of (7.11) and (7.16). The lower gamma tail
\(\mathbb P(S<n/2)\le e^{-c n}\) may be subtracted from the fixed
probability of (7.15), without any independence assertion. This proves a
lower bound \(c/\sqrt n\asymp ca_n\) for all sufficiently large \(n\).
The finitely many dimensions \(2\le n<n_0\) have the same bound after
decreasing the constant: on an open set of directions with
\(|u\cdot\theta|>0\), both sign classes in (7.8) are nonempty, and on a
positive-measure box of exponential coordinates the two reaches are
comparable. The upper bound follows from the pointwise ratio \(1/2\) and
\(\mathbb E|u\cdot\theta|=a_n\).
QED.

Combining (5.4) and (7.5)--(7.6) gives

\[
                              \mathcal U_\mu(E)\asymp a_n
                              \asymp n^{-1/2}.               \tag{7.17}
\]

The simplex cap is therefore sharp but not a counterexample.

## 8. The radial exponential median sphere

Let

\[
 d\mu_n(x)=Z_n^{-1}e^{-c_n|x|}\,dx,\qquad c_n=\sqrt{n+1},   \tag{8.1}
\]

and let \(E_n=B(0,r_n)\), where \(c_nr_n\) is the median of
\({\rm Gamma}(n,1)\). This law is isotropic and the set is balanced.
Rotation invariance removes the \(\theta\)-average.

For \(u=|y|\), the conditional line density is

\[
 f_u(t)=Z(u)^{-1}e^{-c_n\sqrt{u^2+t^2}},                    \tag{8.2}
\]

and

\[
 p_u=
 \begin{cases}
 \displaystyle\int_{-\sqrt{r_n^2-u^2}}^{\sqrt{r_n^2-u^2}}
 f_u(t)\,dt,&u<r_n,\\[2mm]
 0,&u\ge r_n.
 \end{cases}                                                \tag{8.3}
\]

### Lemma 8.1 (tangent-window estimate)

There are universal \(c,C>0\) such that, for all sufficiently large \(n\),

\[
 c\le\sigma_u\le C,\qquad
 c\le h_{n-1}(u)\le C
 \quad\text{when }|u-r_n|\le1,                              \tag{8.4}
\]

where \(h_{n-1}\) is the radial density of the projection \(Y\).
Moreover, on

\[
 r_n-{2\over\sqrt n}\le u\le
 r_n-{1\over2\sqrt n},                                     \tag{8.5}
\]

\[
                              q_u\ge c.                     \tag{8.6}
\]

#### Proof

For \(u\asymp\sqrt n\), convexity and
\[
 \sqrt{u^2+t^2}-u
 \asymp {t^2\over u}
 \quad(|t|\le3)
\]
show that (8.2) is bounded above and below, on fixed \(t\)-windows, by
centered Gaussians of variance \(u/c_n\asymp1\); its tails are
sub-exponential. This proves the conditional variance bounds.

The projection radial density is proportional to
\[
 u^{n-2}e^{-c_nu}\int_{\mathbb R}
 e^{-c_n(\sqrt{u^2+t^2}-u)}\,dt.
\]
The last integral is between universal constants on the stated window.
More explicitly, division by the radial
\({\rm Gamma}(n,c_n)\) density gives
\[
 {h_{n-1}(u)\over
   c_n^nu^{n-1}e^{-c_nu}/\Gamma(n)}
 ={ |S^{n-2}|\over|S^{n-1}|}\,{1\over u}
   \int_{\mathbb R}
   e^{-c_n(\sqrt{u^2+t^2}-u)}\,dt\asymp1.
\]
Stirling's inequalities show that the gamma density is bounded above
and below at
\(u=r_n+O(1)=n/c_n+O(1)\). This proves the second part of (8.4).

For (8.5),
\[
 \sqrt{r_n^2-u^2}\in[c,C].
\]
The conditional Gaussian comparison then puts both \(p_u\) and \(1-p_u\)
above a universal constant, proving (8.6). QED.

Restricting (1.3) to (8.5) gives

\[
                         \mathcal U_{\mu_n}(E_n)\ge {c\over\sqrt n}. \tag{8.7}
\]

For the upper bound no global conditional-tail comparison is needed.
The weighted perimeter of the sphere is the density of
\({\rm Gamma}(n,c_n)\) at its median. Stirling's inequalities and
\(c_nr_n=n+O(1)\) give

\[
                         c\le P_{\mu_n}(E_n)\le C.          \tag{8.8}
\]

Equation (10.3), or directly (2.4)--(2.5), then gives
\[
                         \mathcal U_{\mu_n}(E_n)
                         \le Ca_nP_{\mu_n}(E_n)
                         \le {C\over\sqrt n}.               \tag{8.9}
\]

Thus a smooth high-rank radial interface is another sharp example.

## 9. Product exponential: maximum and checkerboard

Let \(X_i=Z_i-1\), where the \(Z_i\) are independent unit exponentials.
The product law is isotropic.

### 9.1 Median maximum

Put

\[
 E_{\max}=\{\max_iX_i\le L_n\},\qquad
 q_n=e^{-(L_n+1)},\qquad d_n=1-q_n,\qquad d_n^n={1\over2}.  \tag{9.1}
\]

Its weighted perimeter is

\[
 P_\mu(E_{\max})=nq_nd_n^{n-1}={nq_n\over2d_n}\in[1/4,1].  \tag{9.2}
\]

### Lemma 9.1 (facet-race estimate)

There are universal \(c,C>0\) such that

\[
                         {c\over\sqrt n}
 \le\mathcal U_\mu(E_{\max})
 \le {C\over\sqrt n}.                                      \tag{9.3}
\]

#### Proof

Write \(H_n=L_n+1=-\log q_n\), and sample a weighted point on the facet
\(Z_i=H_n\). The other \(Z_j\)'s are independent unit exponentials
conditioned by \(Z_j\le H_n\). Along a random direction, the two support
reaches are the minima of the disjoint races

\[
                  {Z_j\over|\theta_j|}
\quad\text{over the two sign classes of \(\theta_j\)}.      \tag{9.4}
\]

For a uniform spherical direction, with universal probability

\[
 c\sqrt n\le
 \sum_{\substack{j\ne i\\\theta_j>0}}|\theta_j|,
 \sum_{\substack{j\ne i\\\theta_j<0}}|\theta_j|
 \le C\sqrt n.                                             \tag{9.5}
\]

The same event may be intersected with
\(c/\sqrt n\le|\theta_i|\le C/\sqrt n\) while retaining universal
probability. In particular, if \(G_i\) denotes this direction event, then
\(\mathbb E[|\theta_i|\mathbf1_{G_i}]\ge ca_n\).

For \(t=O(n^{-1/2})\), the exact survival probability of either race is

\[
 \prod_j {e^{-|\theta_j|t}-q_n\over1-q_n}
 =\exp\left\{-t\sum_j|\theta_j|+O(t^2)\right\}.             \tag{9.6}
\]

Thus, conditional on (9.5), both support reaches lie in
\([c/\sqrt n,C/\sqrt n]\), with a fixed joint probability. A different
coordinate could hit the upper level \(H_n\) before a support endpoint
only if
\(H_n-Z_j\le C|\theta_j|/\sqrt n\). A union bound and
\(q_n\asymp1/n\) give

\[
 \mathbb P\{\text{a second maximum facet is hit}\mid\theta\}
 \le {Cq_n\over\sqrt n}\sum_j|\theta_j|
 \le {C\over n}.                                          \tag{9.7}
\]

Consequently there is a universal-probability event on which the facet
at \(t=0\) is the only \(E_{\max}\)-boundary crossing in the conditional
support chord and divides that chord in comparable proportions. Since
\(|\sum_j\theta_j|\le\sqrt n\), the log-density variation over a chord
of length \(O(n^{-1/2})\) is bounded by a universal constant. On this
event,

\[
 {q_{\theta,y}\over\sigma_{\theta,y}}
 \asymp f_{\theta,y}(t_{\rm facet}).                        \tag{9.8}
\]

Projection of the facet to \(\theta^\perp\) contributes the Jacobian
\(|\theta_i|\). On the event above a line is charged to only one facet,
and the preceding weighted direction estimate shows that summing (9.8)
over the \(n\) facets gives the lower bound

\[
 \mathcal U_\mu(E_{\max})
 \ge c\,\mathbb E_\theta\sum_i|\theta_i|\,a_i,
\]

where \(a_i\) is the weighted area of facet \(i\) and
\(\sum_i a_i=P_\mu(E_{\max})\). Hence
\[
 \mathbb E_\theta\sum_i|\theta_i|a_i
 =a_nP_\mu(E_{\max})\asymp n^{-1/2}.
\]
The upper bound follows immediately from (10.3) and (9.2). QED.

### 9.2 Parity checkerboard

Let \(m=\log2-1\), so \(\mathbb P(X_i\le m)=1/2\), and define

\[
E_{\rm par}=
 \left\{x:\sum_{i=1}^n\mathbf1_{\{x_i>m\}}
                    \text{ is even}\right\}.               \tag{9.9}
\]

The set is balanced. Its reduced boundary consists of the \(n\) threshold
hyperplanes, each of weighted area \(1/2\), so

\[
                              P_\mu(E_{\rm par})={n\over2}.  \tag{9.10}
\]

At a weighted point on the \(i\)-th threshold facet, the support
endpoints are again the two races
\(\min Z_j/|\theta_j|\), and hence are of order \(n^{-1/2}\) with fixed
probability. For \(j\ne i\), a further parity threshold lies within
distance \(t\) precisely when
\(|Z_j-\log2|\le|\theta_j|t\). Since the exponential density at
\(\log2\) is \(1/2\), for \(t=C/\sqrt n\)

\[
 \mathbb P\{\text{no further threshold before either endpoint}
             \mid\theta\}
 \ge \exp\left\{-{C\over\sqrt n}\sum_j|\theta_j|\right\}
 \ge c                                                     \tag{9.11}
\]

after restricting to the usual \(\ell_1\)-typical directions. Hence,
with universal probability, the facet at the origin is the only parity
flip in a chord of length \(\Theta(n^{-1/2})\) and divides it in
comparable proportions. Combining the product factors for the support
races and the threshold exclusions shows that these properties hold
jointly with universal probability, also after intersecting with
\(c/\sqrt n\le|\theta_i|\le C/\sqrt n\). The density changes by only a
universal factor. The projection calculation from Lemma 9.1 therefore
gives a contribution whose direction average is at least \(ca_n\) times
the weighted facet area. On the stated event the charged facet is unique,
so there is no multiplicity loss. Summing the \(n\) facet contributions
and using (9.10) gives the lower bound below; (10.3) gives the upper
bound:

\[
                              c\sqrt n
 \le\mathcal U_\mu(E_{\rm par})
 \le C\sqrt n.                                             \tag{9.12}
\]

The checkerboard strongly satisfies (0.2); oscillation increases
random-line uncertainty.

## 10. Strength of the proposed theorem

Let

\[
 \mathfrak U_n=
 \inf_{\substack{\mu\ {\rm isotropic\ log\!-\!concave}\\
                  \mu(E)=1/2}}
                    \mathcal U_\mu(E).                      \tag{10.1}
\]

Equations (0.3) and (2.6) prove the one-way implication

\[
                       \mathfrak U_n\ge {c\over\sqrt n}
                       \quad\Longrightarrow\quad
                       \psi_n\ge c'.                        \tag{10.2}
\]

The exact universal comparison available in the other direction is only

\[
                    \mathcal U_\mu(E)
                    \le2\sqrt3\,a_nP_\mu(E),                \tag{10.3}
\]

which follows by rearranging (2.4)--(2.5). A lower perimeter bound does
not reverse (10.3) for an individual set. Thus the uncertainty theorem is
not established as an equivalent reformulation of KLS; as stated for
every balanced Borel set, it contains additional hit-and-run conductance
content.

The model computations identify the only plausible extremal geometries:
one smooth interface, either nearly planar or nearly tangent to a radial
shell. Polyhedral multiplicity and checkerboard oscillation only increase
\(\mathcal U\). Any proof of (0.2) must supply a global uncertainty
principle preventing a balanced set from being almost line-measurable
after variance normalization in most directions. Neither covariance
decomposition (1.4) nor the high-order Radon identity (3.5) supplies that
principle without a new dimension-free argument.

## 11. Named inputs and scope

The only external inputs used above are the following established
statements.

1. One-dimensional log-concave isoperimetry: half-lines minimize
   perimeter, and an isotropic log-concave density has a universal
   Cheeger constant. A stronger normalization than (2.3),
   \(h^2\operatorname{Var}\ge1/3\), is recorded in Appendix B of
   Bobkov--Ledoux,
   [*One-dimensional empirical measures*](https://perso.math.univ-toulouse.fr/ledoux/files/2016/12/MEMO.pdf).
2. Concavity of the Minkowski isoperimetric profile for a log-concave
   probability, including non-smooth densities by approximation, is
   stated as Theorem 18 in Klartag's
   [*Isoperimetric inequalities in high-dimensional convex
   sets*](https://www.weizmann.ac.il/math/klartag/sites/math.klartag/files/uploads/lectures_IHP.pdf).
3. The BV slicing theorem and the affine
   Blaschke--Petkantschin formula are used only in the explicitly
   normalized forms (2.2) and (3.4); the factor in (3.4) can also be
   checked directly by the two-to-one parametrization
   \((\theta,y,s,t)\mapsto(y+s\theta,y+t\theta)\).
4. Sheppard's formula is used only for the exact Gaussian check (4.6).

None of these inputs contains the lower bound (0.2). In particular, the
report establishes the implication (10.2), the exact line/Radon
identities, and the model calculations, but not the proposed
random-line uncertainty theorem itself.
