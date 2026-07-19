# Random half-dimensional heat baths: exact algebra and a no-go for the naive KLS descent

## Executive verdict

Let \(X\sim\mu\) be isotropic and log-concave in \(\mathbb R^n\), let
\(\mathbb R^n=U\oplus V\), and let

\[
 P_Uf=\mathbb E[f(X)\mid X_U],\qquad
 P_Vf=\mathbb E[f(X)\mid X_V].
\]

These are exactly the two conditional-resampling (heat-bath) operators.  For
the random-scan kernel \(M_{U,V}=(P_U+P_V)/2\),

\[
 \langle f,(I-M_{U,V})f\rangle
 =\frac12\mathbb E\operatorname{Var}(f\mid X_U)
  +\frac12\mathbb E\operatorname{Var}(f\mid X_V).       \tag{0.1}
\]

For a fixed split its exact \(L^2\) gap is

\[
 \operatorname{gap}(I-M_{U,V})=\frac{1-\rho(X_U,X_V)}2, \tag{0.2}
\]

where \(\rho\) is maximal correlation.  This exact formula does **not** give
a dimension-free half-block estimate.  For the isotropic uniform measure on
the Euclidean ball, with \(k=\dim U\), \(l=\dim V\),

\[
 \rho(X_U,X_V)
 =\sqrt{\frac{kl}{(k+2)(l+2)}}.                         \tag{0.3}
\]

Thus for \(k=l=n/2\) the fixed-split gap equals
\(2/(n+4)\).  In particular, averaging fixed-split inequalities which use
\((1-\rho)^{-1}\) loses a factor of order \(n\), even on a body whose actual
Poincare constant is dimension-free.

One can instead average the Dirichlet forms *before* taking a gap.  This may
remove the rotating slow modes in the ball example.  However, the canonical
dimension-descent estimate for the conditional variances introduces the
conditional covariance matrices

\[
 C_{V\mid U}(y)=\operatorname{Cov}(X_V\mid X_U=y).
\]

The required averaged energy domination is false, even in dimension two.
For the isotropic log-concave product law

\[
 d\mu(x_1,x_2)=e^{-(x_1+x_2+2)}
   {\bf1}_{\{x_1,x_2\geq-1\}}\,dx_1dx_2,               \tag{0.4}
\]

let the split be Haar among pairs of orthogonal lines.  If
\(e=(1,-1)/\sqrt2\), \(x_t=(t,t)\), and

\[
 W(x)=\mathbb E_{U\oplus V}
 \big[\Pi_V C_{V\mid U}(\Pi_Ux)\Pi_V\big],             \tag{0.5}
\]

then

\[
 e^{\mathsf T}W(x_t)e\geq c(t+1).                      \tag{0.6}
\]

Consequently there is no universal \(B\) such that the Haar-averaged
covariance-weighted conditional energy is at most
\(B\int|\nabla f|^2d\mu\) for every smooth \(f\).  A localized bump near
\(x_t\) makes the ratio tend to infinity.  The same construction survives
truncation and affine re-isotropization, so it is not an artifact of an
infinite support.

The other natural descent, through the lower-dimensional marginal
\(P_Uf\), also fails at its formal commutation step.  For a smooth density
\(e^{-\Phi(y,z)}\),

\[
 \nabla_y\mathbb E[f(y,Z)\mid Y=y]
 =\mathbb E[\nabla_yf\mid y]
  -\operatorname{Cov}(f,\nabla_y\Phi\mid y),            \tag{0.7}
\]

with an additional boundary-flux term for convex bodies and one-sided
supports.  Neither term is controlled for free.  In (0.4), for the
sum/difference split, \(f(Y,Z)=Y\) satisfies

\[
 \mathbb E[Y\mid Z]=|Z|-1/\sqrt2,\qquad \partial_Zf=0, \tag{0.8}
\]

so the hoped-for gradient commutation fails in the strongest possible way.

**Conclusion.**  The exact heat-bath algebra is useful, but the canonical
random-half-block argument does not yield a dimension-free recurrence for
the KLS/Poincare constant.  The fixed-split maximal-correlation route loses
\(n\) on the ball, while the average-first conditional-covariance route has
an explicit Haar-averaged counterexample.  An average-first heat-bath gap by
itself is not refuted here; it is simply insufficient for dimension descent
without a new joint estimate that avoids pointwise conditional Poincare
weights.

---

## 1. Exact two-block identities

Write \(Y=X_U=\Pi_UX\), \(Z=X_V=\Pi_VX\), and work in
\(H=L^2(\mu)\).  Conditional expectation gives two self-adjoint orthogonal
projections

\[
 A=P_U=\mathbb E[\,\cdot\mid Y],\qquad
 B=P_V=\mathbb E[\,\cdot\mid Z].                       \tag{1.1}
\]

Operationally, \(A\) holds \(Y\) fixed and redraws \(Z\) from
\(\mu(dz\mid Y)\); \(B\) does the reverse.  Thus the random-scan block
heat bath is

\[
 M=\frac{A+B}{2}.                                      \tag{1.2}
\]

Since an orthogonal projection satisfies
\(\langle f,(I-A)f\rangle=\|f-Af\|_2^2\), one has the exact identities

\[
\begin{aligned}
 \mathcal D_{U,V}(f)
 &:=\langle f,(I-M)f\rangle \\
 &=\frac12\|f-Af\|_2^2+\frac12\|f-Bf\|_2^2 \\
 &=\frac12\mathbb E\operatorname{Var}(f\mid Y)
   +\frac12\mathbb E\operatorname{Var}(f\mid Z),      \tag{1.3}
\end{aligned}
\]

and, for centered \(f\),

\[
 \mathcal D_{U,V}(f)
 =\operatorname{Var}(f)
  -\frac12\operatorname{Var}(Af)
  -\frac12\operatorname{Var}(Bf).                     \tag{1.4}
\]

No log-concavity is needed for (1.3)--(1.4).

### 1.1 Maximal correlation is the exact fixed-split obstruction

Let

\[
 \rho(Y,Z)=\sup
 \left\{\mathbb E[g(Y)h(Z)]:
 \mathbb Eg=\mathbb Eh=0,\ \mathbb Eg^2=\mathbb Eh^2=1\right\}. \tag{1.5}
\]

Equivalently, \(\rho\) is the operator norm of the conditional-expectation
map from the centered \(Z\)-measurable subspace to the centered
\(Y\)-measurable subspace.  The two-projection theorem gives

\[
 \|A+B\|_{H_0}=1+\rho,                                 \tag{1.6}
\]

and therefore

\[
 \mathcal D_{U,V}(f)\geq\frac{1-\rho(Y,Z)}2
                         \operatorname{Var}(f).         \tag{1.7}
\]

The constant is optimal.  For completeness, on every principal two-plane
associated with a singular value \(s\) of \(AB\), the eigenvalues of
\((A+B)/2\) are \((1+s)/2\) and \((1-s)/2\).  On the orthogonal remainder
they are among \(0,1/2\).  The largest centered eigenvalue is consequently
\((1+\rho)/2\).

For a full-dimensional log-concave law, the common centered subspace of
\(A\) and \(B\) is zero.  Indeed, if
\(g(\Pi_Ux)=h(\Pi_Vx)\) on the interior of the convex support, local
variation first in a \(V\)-direction and then in a \(U\)-direction makes
both functions locally constant; connectedness makes them globally
constant.  Degenerate support should first be reduced to its affine hull.

### 1.2 Haar averaging

Let \(k+l=n\), and average over Haar orthogonal splits with these
dimensions.  Define

\[
 \overline M=\mathbb E_{U\oplus V}M_{U,V},\qquad
 \overline{\mathcal D}(f)=
 \mathbb E_{U\oplus V}\mathcal D_{U,V}(f).             \tag{1.8}
\]

Then

\[
 \overline{\mathcal D}(f)
 =\operatorname{Var}(f)-\frac12\mathbb E_{U\oplus V}
 \left(\operatorname{Var}(P_Uf)+\operatorname{Var}(P_Vf)\right). \tag{1.9}
\]

The averaged gap is exactly

\[
 \bar\gamma_{n;k,l}
 =\inf_{\mathbb Ef=0}\frac{\overline{\mathcal D}(f)}{\|f\|_2^2}. \tag{1.10}
\]

If \(k=l\), the laws of \(U\) and \(V\) agree and

\[
 \overline M=\mathbb E_U P_U,\qquad
 \overline{\mathcal D}(f)=
 \mathbb E_U\mathbb E\operatorname{Var}(f\mid X_U).   \tag{1.11}
\]

One always has the weaker consequence of the fixed-split gaps,

\[
 \bar\gamma_{n;k,l}\geq
 \frac12\mathbb E_{U\oplus V}[1-\rho(X_U,X_V)],        \tag{1.12}
\]

but (1.12) can be very far from the actual averaged gap: the optimizer in
(1.7) can rotate with the split.

---

## 2. The formal KLS recurrence and all of its load-bearing inputs

Let \(K_d\) be the supremum of \(C_P(\nu)\) over isotropic log-concave
probabilities in dimensions at most \(d\).  By affine normalization, every
log-concave \(\nu\) in dimension \(d\), with covariance \(C_\nu\), obeys

\[
 \operatorname{Var}_\nu h
 \leq K_d\int\langle C_\nu\nabla h,\nabla h\rangle d\nu. \tag{2.1}
\]

Apply (2.1) on the conditional fibers.  With

\[
 C_{V\mid U}(y)=\operatorname{Cov}(Z\mid Y=y),\qquad
 C_{U\mid V}(z)=\operatorname{Cov}(Y\mid Z=z),          \tag{2.2}
\]

one gets

\[
\begin{aligned}
 \mathbb E\operatorname{Var}(f\mid Y)
 &\leq K_l\,\mathbb E
 \langle C_{V\mid U}(Y)\nabla_Vf,\nabla_Vf\rangle,\\
 \mathbb E\operatorname{Var}(f\mid Z)
 &\leq K_k\,\mathbb E
 \langle C_{U\mid V}(Z)\nabla_Uf,\nabla_Uf\rangle.      \tag{2.3}
\end{aligned}
\]

Put \(m=\max(k,l)\) and

\[
\begin{aligned}
 \overline{\mathcal Q}(f):=\frac12\mathbb E_{U\oplus V,X}
 \big[&\langle C_{V\mid U}(X_U)\nabla_Vf(X),\nabla_Vf(X)\rangle\\
 &+\langle C_{U\mid V}(X_V)\nabla_Uf(X),\nabla_Uf(X)\rangle\big]. \tag{2.4}
\end{aligned}
\]

Equations (1.3) and (2.3) give the exact formal chain

\[
 \bar\gamma_{n;k,l}\operatorname{Var}(f)
 \leq\overline{\mathcal D}(f)
 \leq K_m\overline{\mathcal Q}(f).                    \tag{2.5}
\]

Thus the naive recurrence needs

\[
 \overline{\mathcal Q}(f)\leq\beta_{n;k,l}
 \int|\nabla f|^2d\mu,                                \tag{2.6}
\]

which would imply

\[
 C_P(\mu)\leq K_m\frac{\beta_{n;k,l}}{\bar\gamma_{n;k,l}}. \tag{2.7}
\]

There are three separate issues.

1. A fixed-split use of (1.7) requires a lower bound on
   \(1-\rho(X_U,X_V)\).  Section 3 shows that this is only \(O(1/n)\) for
   the ball.
2. Isotropy controls only the **mean** conditional covariance:

   \[
    \mathbb E C_{V\mid U}(Y)
    =I_V-\operatorname{Cov}(\mathbb E[Z\mid Y])\preceq I_V. \tag{2.8}
   \]

   It does not control (2.4), because \(\nabla f(X)\) can concentrate where
   the conditional covariance is large.  Section 4 shows that (2.6) fails
   even after Haar averaging.
3. Replacing (2.3) by a marginal induction for \(P_Uf\) requires a gradient
   estimate for conditional expectations.  Section 5 shows the missing
   score/boundary term and gives explicit failures.

There is also a constants issue.  Even a one-step bound
\(K_n\leq C K_{\lceil n/2\rceil}\) with a fixed \(C>1\) iterates to a
power of \(n\), not to KLS.  A closing recurrence would need a max-type or
additive estimate, or ratios whose product stays bounded.  The Gaussian
calibration has \(\beta=\bar\gamma=1/2\); any loss from separately estimating
the two sides of (2.5) is therefore consequential.

---

## 3. Maximal-correlation audit: the Euclidean ball

Let \(\mu\) be uniform on the ball \(B_R^n\), with
\(R=\sqrt{n+2}\), so that \(\mu\) is isotropic.  For any orthogonal split
of dimensions \(k,l\), write \(Y=X_U\), \(Z=X_V\).  The fourth moments are

\[
 \mathbb EX_i^4=\frac{3(n+2)}{n+4},\qquad
 \mathbb EX_i^2X_j^2=\frac{n+2}{n+4}\quad(i\ne j).     \tag{3.1}
\]

Consequently

\[
\begin{aligned}
 \operatorname{Cov}(|Y|^2,|Z|^2)&=-\frac{2kl}{n+4},\\
 \operatorname{Var}(|Y|^2)&=\frac{2k(l+2)}{n+4},\\
 \operatorname{Var}(|Z|^2)&=\frac{2l(k+2)}{n+4}.       \tag{3.2}
\end{aligned}
\]

Thus already the quadratic radial statistics give

\[
 \rho(Y,Z)\geq
 \sqrt{\frac{kl}{(k+2)(l+2)}}.                         \tag{3.3}
\]

In fact equality holds.  To see this, set

\[
 A=|Y|^2/R^2,\quad B=|Z|^2/R^2,\quad S=1-A-B.
\]

Then \((A,B,S)\) is Dirichlet with parameters
\((k/2,l/2,1)\).  Conditional on \((A,B)\), the two angular variables are
independent and uniform, so cross-correlation comes only from functions of
\(A\) and \(B\).  Orthogonalizing the monomials gives the Jacobi polynomial
blocks.  The squared canonical correlation in degree \(r\geq1\) is

\[
 \rho_r^2=
 \frac{(k/2)_r(l/2)_r}{(k/2+1)_r(l/2+1)_r}
 =\frac{kl}{(k+2r)(l+2r)}.                             \tag{3.4}
\]

This decreases with \(r\), and degree one gives (0.3).  Formula (3.4) can
also be checked directly from conditional beta moments: conditional
expectation preserves polynomial degree, is triangular on monomials, and
orthogonality removes all lower-degree entries.

For \(k=l=n/2\),

\[
 \rho=\frac{k}{k+2}=1-\frac{4}{n+4},\qquad
 \operatorname{gap}(I-M_{U,V})=\frac1{k+2}=\frac2{n+4}. \tag{3.5}
\]

This holds for every split, by rotational invariance.  Hence neither a
typical-split argument nor averaging the constants in (1.7) gives a
dimension-free gap.

The slow fixed-split statistic is essentially a normalized difference of
\(|Y|^2\) and \(|Z|^2\).  It rotates with \(U\).  This explains why (3.5)
does not by itself refute a gap for the average-first operator
\(\overline M\).  For example, if
\(f=|X|^2-n\), then

\[
 P_Uf=\frac{2}{l+2}|Y|^2+\text{constant},\qquad
 \frac{\operatorname{Var}(P_Uf)}{\operatorname{Var}(f)}
 =\frac{2k}{n(l+2)}.                                   \tag{3.6}
\]

For equal half-blocks this is \(2/(n+4)\), so the radial total norm is
almost completely killed rather than almost preserved.  Averaging can
repair rotating slow modes; it cannot repair the conditional-energy
counterexample in the next section.

Finally, the ball has particularly benign covariance weights:

\[
 C_{V\mid U}(y)=\frac{R^2-|y|^2}{l+2}I_V
 \preceq\frac{n+2}{l+2}I_V.                            \tag{3.7}
\]

For half-blocks the last constant is less than \(2\).  Thus (3.5) is a pure
maximal-correlation obstruction, not a bad-fiber-covariance obstruction.

---

## 4. Conditional covariance audit: an explicit Haar-averaged no-go

### 4.1 The law and one exceptional split

Let \(E_1,E_2\) be independent rate-one exponentials and
\(X_i=E_i-1\).  This gives (0.4), and each coordinate has mean zero and
variance one.  Hence \(\mu\) is isotropic and log-concave.

Use the sum/difference coordinates

\[
 Y=\frac{X_1+X_2}{\sqrt2},\qquad
 Z=\frac{X_1-X_2}{\sqrt2},\qquad
 T=E_1+E_2=\sqrt2Y+2.                                  \tag{4.1}
\]

Conditionally on \(T\), \(E_1/T\) is uniform on \([0,1]\).  Therefore

\[
 Z\mid Y\ \text{ is uniform on }[-T/\sqrt2,T/\sqrt2],
 \qquad
 \operatorname{Var}(Z\mid Y)=\frac{T^2}{6}.            \tag{4.2}
\]

The conditional covariance is unbounded along this split.  This alone is
not yet a Haar counterexample, because the exactly sum/difference split has
zero Haar measure.  A window of nearby splits supplies the missing mass.

### 4.2 Nearby directions and the \(t^2\times t^{-1}\) mechanism

Let \(v=(v_1,v_2)\) be a unit vector spanning \(V\), let \(U=V^\perp\), and
condition on the \(U\)-projection of \(x_t=(t,t)\).  Every point of the
fiber is \(x_t+sv\).  When \(v_1>0>v_2\), feasibility is the interval

\[
 I_v=\left[-\frac{t+1}{v_1},\frac{t+1}{-v_2}\right],   \tag{4.3}
\]

and the conditional density of \(s\) is proportional to

\[
 \exp[-(v_1+v_2)s]\,{\bf1}_{I_v}(s).                  \tag{4.4}
\]

Put \(e=(1,-1)/\sqrt2\).  If the angle between \(v\) and \(e\) is at most
\(c_0/(t+1)\), then

\[
 |I_v|\asymp t+1,\qquad
 |v_1+v_2|\,|I_v|\leq Cc_0.                            \tag{4.5}
\]

Choose \(c_0\) small.  The density ratio across the interval is then bounded
by a universal constant.  Comparing the left and right quarter-intervals
shows

\[
 \operatorname{Var}(s\mid\Pi_UX=\Pi_Ux_t)
 \geq c_1(t+1)^2.                                      \tag{4.6}
\]

For clarity, the elementary comparison used here is: if a density on an
interval of length \(L\) has maximum/minimum ratio at most \(R\), each end
quarter has probability at least \(1/(4R)\), and the independent-copy
formula for variance gives \(\operatorname{Var}\geq L^2/(64R^2)\).

The angular window in (4.5) has Haar measure at least \(c_2/(t+1)\), and
\(|e\cdot v|^2\geq1/2\) there.  It follows that the matrix in (0.5) satisfies

\[
 e^{\mathsf T}W(x_t)e
 \geq \frac{c_2}{t+1}\cdot\frac12\cdot c_1(t+1)^2
 \geq c_3(t+1).                                        \tag{4.7}
\]

This proves (0.6).

The estimate is uniform for \(x\) in a fixed small rectangle about \(x_t\):
replacing \(t+1\) in the two endpoints of (4.3) by
\(x_1+1,x_2+1\) changes (4.5)--(4.7) only by universal factors.  Here is a
fully local smooth test.  Put \(d=(1,1)/\sqrt2\), fix a small
\(\varepsilon>0\), and choose
\(\eta\in C_c^\infty((-2,2))\) with \(\eta=1\) on \([-1,1]\).  Define

\[
 f_t(x)=\big(e\cdot(x-x_t)\big)
 \eta\!\left(\frac{e\cdot(x-x_t)}{\varepsilon}\right)
 \eta\!\left(\frac{d\cdot(x-x_t)}{\varepsilon}\right). \tag{4.8a}
\]

On the inner square where both displayed coordinates have absolute value at
most \(\varepsilon\), \(\nabla f_t=e\).  On its support,
\(|\nabla f_t|\leq C_\eta\), uniformly in \(t\).  Moreover,
\(d\mu=e^{-2t-2-\sqrt2\,d\cdot(x-x_t)}dx\) there, so the density and the
inner-to-outer mass ratio vary by universal factors.  Positivity of every
conditional covariance and (4.7) on the inner square therefore yield

\[
 \frac{\mathbb E_{U\oplus V,X}
  \langle C_{V\mid U}(X_U)\nabla_Vf_t,\nabla_Vf_t\rangle}
 {\int|\nabla f_t|^2d\mu}
 \geq c_4(t+1).                                        \tag{4.8}
\]

Thus no finite universal \(\beta\) can satisfy (2.6).

This is robust.  Truncate each exponential to a long interval, center and
rescale the coordinates, and choose \(t\) well inside the truncation.  The
same lower bound holds up to universal factors, while the affine
normalization converges to the identity.  Hence compactly supported
log-concave approximations require arbitrarily large constants as well.

### 4.3 Why the mean covariance identity does not help

For every fixed split, total covariance gives

\[
 \mathbb E C_{V\mid U}(Y)\preceq I_V.                 \tag{4.9}
\]

In the example, the large weights occur far in a tail and in a shrinking
angular window, so (4.9) remains true.  But a test function is allowed to put
its gradient in precisely that tail and direction.  Matrix averaging before
multiplication by \(|\nabla f|^2\) is therefore invalid.  Notice also that
constant gradients are harmless:

\[
 \mathbb E\langle C_{V\mid U}(Y)\Pi_Va,\Pi_Va\rangle
 =\mathbb E\operatorname{Var}(a\cdot X\mid X_U)
 \leq|\Pi_Va|^2.                                       \tag{4.10}
\]

The obstruction is necessarily a localized, variable-gradient one.

---

## 5. Gradient of a conditional expectation

Another exact variance decomposition is

\[
 \operatorname{Var}(f)
 =\mathbb E\operatorname{Var}(f\mid Y)
  +\operatorname{Var}(P_Uf).                           \tag{5.1}
\]

The marginal law of \(Y\) is isotropic and log-concave, so it is tempting to
apply \(K_k\) to the second term.  The missing step is a bound on
\(\nabla P_Uf\).

Suppose first that
\(d\mu(y,z)=Z^{-1}e^{-\Phi(y,z)}dydz\), with a smooth full-support density
for which differentiation under the integral is valid.  If
\(m(y)=\mathbb E[f(y,Z)\mid Y=y]\), direct differentiation gives

\[
 \nabla m(y)
 =\mathbb E[\nabla_y f\mid y]
  -\operatorname{Cov}(f,\nabla_y\Phi\mid y).            \tag{5.2}
\]

Equivalently, if
\(s_y(z)=\nabla_y\log p(z\mid y)\) is the conditional score,

\[
 \nabla m(y)=\mathbb E[\nabla_yf\mid y]
              +\mathbb E[f\,s_y\mid y].               \tag{5.3}
\]

For a convex body or a one-sided density, the support of the conditional
law moves with \(y\).  Reynolds' transport formula adds a boundary flux.
Treating the convex constraint as \(\Phi=+\infty\), that flux is the
singular part of the score term.  It cannot be dropped.

The shifted-exponential example makes this completely explicit.  The
difference \(D=E_1-E_2\) has the Laplace density.  Given
\(D=d\), the sum satisfies

\[
 T\mid D=d=|d|+E,\qquad E\sim\operatorname{Exp}(1).     \tag{5.4}
\]

Consequently

\[
 \mathbb E[Y\mid Z=z]=|z|-\frac1{\sqrt2}.              \tag{5.5}
\]

For \(f(y,z)=y\), the right side of the hoped-for commutation estimate
\(|\nabla_zP_Vf|^2\leq\mathbb E[|\nabla_zf|^2\mid Z]\) is zero, whereas
the left side equals one almost everywhere.

The ball supplies a second transparent boundary example.  For
\(f(y,z)=|z|^2\), (3.7) gives

\[
 P_Uf(y)=\frac{l}{l+2}(R^2-|y|^2),\qquad
 \nabla_yP_Uf=-\frac{2l}{l+2}y,                        \tag{5.6}
\]

although \(\nabla_yf=0\).  Thus neither log-concavity nor rotational
symmetry makes conditional expectation commute with the retained-block
gradient.

Controlling the covariance/score term in (5.2) by conditional Poincare and
Fisher information merely creates a new curvature or conditional-score
hypothesis.  Such a hypothesis is not supplied by isotropy and is not a free
local-to-global principle.

---

## 6. Countertests on the requested model families

### 6.1 Cube

Let \(X_i\) be independent uniform variables on
\([-\sqrt3,\sqrt3]\).  For a coordinate split:

\[
 \rho(X_U,X_V)=0,\qquad C_{V\mid U}=I_V,               \tag{6.1}
\]

and conditional expectation commutes with derivatives in retained
coordinates.  This is the ideal product calibration.

It is not stable under rotation.  In dimension two set

\[
 Y=(X_1+X_2)/\sqrt2,\qquad Z=(X_1-X_2)/\sqrt2.
\]

Their support is the diamond
\(|Y|+|Z|\leq\sqrt6\).  Given \(Y=y\), \(Z\) is uniform on

\[
 [-L(y),L(y)],\qquad L(y)=\sqrt6-|y|,                  \tag{6.2}
\]

so

\[
 \operatorname{Var}(Z\mid Y=y)=L(y)^2/3\leq2.          \tag{6.3}
\]

The covariance weight stays bounded.  Dependence and boundary motion are
nevertheless visible.  Direct moments give

\[
 \operatorname{Corr}(Y^2,Z^2)=-\frac37,                \tag{6.4}
\]

so \(\rho\geq3/7\), and

\[
 P_U(Z^2)=L(Y)^2/3                                     \tag{6.5}
\]

has nonzero \(Y\)-gradient although \(Z^2\) has zero \(Y\)-gradient.

### 6.2 Regular simplex

Let \(P=(P_1,\ldots,P_m)\) be uniform on the probability simplex, so
\(P\sim\operatorname{Dirichlet}(1,\ldots,1)\), and put

\[
 X=\sqrt{m(m+1)}\left(P-\frac1m\mathbf1\right)
 \in\mathbf1^\perp.                                    \tag{6.6}
\]

This is the isotropic regular simplex in dimension \(n=m-1\).  Partition
the coordinates into \(A,B\), with sizes \(a,b\).  Let
\(V=\{v:\operatorname{supp}v\subset A,\ \sum_{i\in A}v_i=0\}\), and
\(U=V^\perp\cap\mathbf1^\perp\).  Then

\[
 \dim V=a-1,\qquad \dim U=b.                           \tag{6.7}
\]

For \(a=b=m/2\), this is a half-dimensional split up to one dimension.
The \(U\)-coordinate reveals
\(S=\sum_{i\in A}P_i\), and, conditionally on \(S=s\),
\(P_A/s\sim\operatorname{Dirichlet}(1^a)\).  Hence

\[
 C_{V\mid U}(s)
 =\frac{m(m+1)s^2}{a(a+1)}I_V.                         \tag{6.8}
\]

Since \(S\sim\operatorname{Beta}(a,b)\),

\[
 \mathbb EC_{V\mid U}=I_V,\qquad
 \sup_s\|C_{V\mid U}(s)\|_{\rm op}
 =\frac{m(m+1)}{a(a+1)}.                               \tag{6.9}
\]

For half-blocks the supremum tends to \(4\), so this natural simplex split
does not produce the covariance blow-up of Section 4.

It does exhibit both dependence and gradient creation.  Write
\(Q=\sum_{i\in A}(P_i/S-1/a)^2\), independent of \(S\).  Then

\[
 \|X_V\|^2=m(m+1)S^2Q,\qquad
 P_U\|X_V\|^2=\frac{m(m+1)(a-1)}{a(a+1)}S^2.          \tag{6.10}
\]

Thus a function with zero \(U\)-gradient acquires a nonzero conditional-
expectation gradient.  Also

\[
 \rho(X_U,X_V)\geq\operatorname{Corr}(S,S^2Q).          \tag{6.11}
\]

The correlation is explicit from

\[
\begin{gathered}
 \mathbb ES^r=\frac{(a)_r}{(m)_r},\qquad
 \mathbb EQ=\frac{a-1}{a(a+1)},\\
 \operatorname{Var}(Q)=
 \frac{4(a-1)}{(a+1)^2(a+2)(a+3)}.                    \tag{6.12}
\end{gathered}
\]

For \(a=b\to\infty\), (6.11) tends to \(1/\sqrt3\).  Thus the natural
simplex split has a constant, nonzero block dependence, but not the
near-deterministic \(1-O(1/n)\) dependence of the ball.

### 6.3 Euclidean ball

The complete audit is in Section 3:

\[
 \rho=\sqrt{\frac{kl}{(k+2)(l+2)}},\qquad
 C_{V\mid U}(y)=\frac{n+2-|y|^2}{l+2}I_V.              \tag{6.13}
\]

For half-blocks the covariance weight is bounded by two, but the fixed-split
gap is \(2/(n+4)\).  Conditional-expectation gradients have a boundary
velocity term, as (5.6) shows.

### 6.4 Product of one-sided exponentials

For a coordinate split all blocks are independent, so

\[
 \rho=0,\qquad C_{V\mid U}=I_V.                        \tag{6.14}
\]

For the two-dimensional sum/difference split, (4.2) gives an unbounded
conditional covariance.  The dependence can also be seen without solving
the full maximal-correlation problem.  Since \(T\sim\Gamma(2,1)\) and,
conditionally on \(T\),

\[
 |Z|=\frac{T A}{\sqrt2},\qquad A\sim\operatorname{Unif}[0,1]
 \text{ independent of }T,                            \tag{6.15}
\]

one obtains

\[
 \operatorname{Corr}(Y,|Z|)=\frac1{\sqrt2}.            \tag{6.16}
\]

Hence \(\rho\geq1/\sqrt2\).  More importantly, nearby orientations turn
the single bad split into the Haar-averaged divergence (4.7).

### 6.5 General radial law

Let \(X\) be spherically symmetric and isotropic, and put

\[
 q=\frac{\mathbb E|X|^4}{n(n+2)},\qquad
 \delta=\operatorname{Var}(|X|^2).
\]

Then

\[
 q-1=\frac{\delta-2n}{n(n+2)}.                         \tag{6.17}
\]

For every \(k+l=n\), spherical moment identities give

\[
\begin{aligned}
 \operatorname{Cov}(|Y|^2,|Z|^2)&=kl(q-1),\\
 \operatorname{Var}(|Y|^2)&=2kq+k^2(q-1),\\
 \operatorname{Var}(|Z|^2)&=2lq+l^2(q-1).              \tag{6.18}
\end{aligned}
\]

Therefore

\[
 \rho(Y,Z)\geq
 \frac{kl|q-1|}
 {\sqrt{(2kq+k^2(q-1))(2lq+l^2(q-1))}}.                \tag{6.19}
\]

For equal half-blocks and a thin radial law with \(\delta=o(n)\), the right
side tends to one.  If \(\delta=O(1)\), it is
\(1-O(1/n)\).  The uniform ball has
\(\delta=4n/(n+4)\) and attains the exact value in (3.5).  By contrast, the
Gaussian has \(\delta=2n\), \(q=1\), and independent orthogonal blocks.

This calculation identifies the geometry behind the maximal-correlation
failure: a thin total radius makes the two half-block radii almost
complementary.  It is not evidence for a large global Poincare constant.

---

## 7. Precise no-go statement and what remains viable

The preceding calculations prove the following method-level statement.

**Proposition (no canonical separated block descent).**  There are no
universal constants \(c,B>0\) for which the following program works for all
isotropic log-concave \(\mu\):

1. use the fixed-split maximal-correlation bound
   \(\mathcal D_{U,V}\geq c\operatorname{Var}\) on a typical Haar
   half-split; and
2. dominate the Haar-averaged affine-normalized conditional Poincare energy
   by \(B\int|\nabla f|^2d\mu\).

The ball disproves item 1, with the best \(c=2/(n+4)\).  The product of two
shifted exponentials disproves item 2, with the ratio in (4.8) unbounded.
Moreover, marginal induction cannot replace item 2 using the identity
\(\nabla P_Uf=\mathbb E[\nabla_Uf\mid X_U]\), because that identity is false
by (5.2), (5.5), and (5.6).

This does **not** prove that the average-first operator \(\overline M\) lacks
a universal spectral gap.  Establishing such a gap would be an interesting
standalone lemma.  But it would not close KLS through (2.3): Section 4 shows
that the natural lower-dimensional conditional Poincare upper bound has an
unbounded averaged coefficient.  A viable block argument must instead prove
a genuinely joint inequality, schematically

\[
 \overline{\mathcal D}(f)
 \leq C\int|\nabla f|^2d\mu                            \tag{7.1}
\]

by exploiting cancellation/localization inside the actual conditional
variance, not by first replacing each fiber with its worst affine Poincare
weight.  But (7.1) together with a lower gap is already a direct global
Poincare proof; calling it dimension descent supplies no additional
leverage.  The explicit tail fibers in Section 4 are the test that any such
joint lemma must pass.

The safe reusable outputs are therefore:

- the exact identities (1.3), (1.4), and (1.9);
- the exact fixed-split gap (1.7);
- the formal conditional recursion (2.5), with its covariance weights left
  intact;
- the mean identity (2.8), which may be used only against gradients
  independent of the conditioning variable;
- the conditional-gradient formula (5.2), including boundary flux in the
  nonsmooth case;
- the ball maximal-correlation formula (0.3); and
- the Haar-averaged covariance-weight counterexample (4.7)--(4.8).

None of the missing estimates may be replaced by KLS, a uniform conditional
covariance assertion, or an unproved local-to-global principle.
