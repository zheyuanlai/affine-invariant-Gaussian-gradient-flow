# Target T3: extremal Lipschitz witnesses

Let \(\mu\) be an isotropic log-concave probability measure on
\(\mathbb R^n\), and put
\[
  \mathcal O_1(\mu)=\sup_{\operatorname{Lip}(f)\leq 1}
       \int |f-\mu f|\,d\mu .
\]
All constants below are numerical and independent of \(n\).  The purpose of
this note is to derive necessary structure of a hypothetical witness with
large \(\mathcal O_1\), and then to isolate exactly where the proposed
exclusion mechanism stops.

## 0. Elementary structure of any bad witness

Fix a 1-Lipschitz \(f\), put \(g=f-\mu f\), and write
\[
 A=\mathbb E|g|,\qquad V=\mathbb Eg^2.
\]

### Lemma 0.1 (energy and dimension)

\[
 A^2\leq V\leq n,
 \qquad C_P(\mu)\geq V\geq A^2.                                    \tag{0.1}
\]

#### Proof

The first lower bound is Jensen.  If \(X'\) is an independent copy of
\(X\), then
\[
 V=\frac12\mathbb E(f(X)-f(X'))^2
 \leq\frac12\mathbb E|X-X'|^2=n.
\]
Rademacher's theorem gives \(\int|\nabla f|^2d\mu\leq1\), so the Rayleigh
quotient in the definition of \(C_P\) is at least \(V\).  \(\square\)

Thus a large T3 witness is already a low-energy spectral witness.  Using a
universal Poincare bound at this point would simply assume KLS.

### Lemma 0.2 (an infinitesimal bottleneck and a finite neck)

Let \(P_\mu(B)\) denote weighted BV perimeter.  Some sublevel set
\(S_t=\{f\leq t\}\) satisfies
\[
 \frac{P_\mu(S_t)}{\min(\mu(S_t),1-\mu(S_t))}\leq\frac2A.            \tag{0.2}
\]
There is also a number \(t\) for which, with
\[
 L=\{f\leq t\},\quad
 M=\{t<f\leq t+A/8\},\quad
 R=\{f>t+A/8\},
\]
both outer sets have positive mass and
\[
 \operatorname{dist}(L,R)\geq A/8,
 \qquad \mu(M)\leq\frac13\min(\mu(L),\mu(R)).                       \tag{0.3}
\]

#### Proof

Let \(m\) be a median of \(Y=f(X)\) and set
\(D=\mathbb E|Y-m|\).  A median minimizes absolute deviation, while
\(|m-\mathbb EY|\leq D\), and hence
\[
 A/2\leq D\leq A.                                                   \tag{0.4}
\]
If \(F(t)=\mathbb P(Y\leq t)\), layer cake and weighted coarea give
\[
 D=\int_{\mathbb R}\min(F(t),1-F(t))\,dt,
 \qquad
 \int_{\mathbb R}P_\mu(\{f\leq t\})\,dt
 =\int|\nabla f|\,d\mu\leq1.
\]
Averaging their ratio proves (0.2).

For the finite neck, fix \(0<r<D\) and define
\[
 N(t)=F(t+r)-F(t),\qquad
 H(t)=\min(F(t),1-F(t+r)).
\]
Fubini gives \(\int N(t)dt=r\).  Moreover,
\[
 H(t)\geq\min(F(t),1-F(t))-N(t),
\]
so \(\int H(t)dt\geq D-r\).  Some \(t\) therefore satisfies
\[
 N(t)\leq\frac r{D-r}H(t).
\]
Take \(r=A/8\).  By (0.4), \(r/(D-r)\leq1/3\).  Lipschitzness gives the
separation in (0.3).  \(\square\)

The smaller side in (0.2)--(0.3) need not have constant mass for an
arbitrary witness.

### Lemma 0.3 (half-mass concentration witness)

There is a closed set \(S\) with \(\mu(S)\geq1/2\) such that
\[
 \mathbb E\,d(X,S)\geq A/4.                                       \tag{0.5}
\]

#### Proof

At least one of \(\mathbb E(Y-m)_+\) and \(\mathbb E(m-Y)_+\) is at least
\(D/2\geq A/4\).  In the first case use \(S=\{f\leq m\}\) and the pointwise
inequality \((f-m)_+\leq d(\cdot,S)\); the other case is symmetric.  Replace
the level set by its closure if necessary.  \(\square\)

### Lemma 0.4 (no low-dimensional ridge explanation)

If \(P\) is an orthogonal projection of rank \(k\) and \(\phi\) is
1-Lipschitz, then
\[
 \mathbb E|f(X)-\phi(PX)|\geq\frac12(A-\sqrt k)_+.                  \tag{0.6}
\]

#### Proof

For \(h=\phi(PX)\) and an independent \(X'\),
\[
 \operatorname{Var}h
 =\frac12\mathbb E(h(X)-h(X'))^2
 \leq\frac12\mathbb E|P(X-X')|^2=k.
\]
If \(e=\mathbb E|f-h|\), then \(|\mathbb Ef-\mathbb Eh|\leq e\), and
\[
 A\leq2e+\mathbb E|h-\mathbb Eh|\leq2e+\sqrt k.
\]
\(\square\)

In particular, a witness with \(A\gg1\) is far in \(L^1\) from every
1-Lipschitz affine function, and one with \(A\gg\sqrt k\) is not explained
by a \(k\)-dimensional ridge.

## 1. Exact cut--transport duality

### Lemma 1 (cut--transport formula)

For every probability measure with a finite first moment,
\[
 \boxed{\quad
 \mathcal O_1(\mu)
 =2\sup_{0<\mu(E)<1}\mu(E)\mu(E^c)
   W_1(\mu_E,\mu_{E^c}).\quad}                                      \tag{1}
\]
Here \(\mu_E=\mu(\,\cdot\mid E)\).

#### Proof

For any integrable \(u\),
\[
 \int|u|\,d\mu=\sup_{\|h\|_\infty\leq1}\int hu\,d\mu.
\]
Consequently, after exchanging two suprema,
\[
 \mathcal O_1(\mu)
 =\sup_{\|h\|_\infty\leq1}\sup_{\operatorname{Lip}(f)\leq1}
   \int f(h-\mu h)\,d\mu .                                         \tag{2}
\]
The inner supremum is the Kantorovich--Rubinstein norm of the zero-mass
signed measure \((h-\mu h)\mu\).

It is enough to use sign functions.  Indeed, if \(U\) is uniform on
\([-1,1]\) and
\[
 s_U(x)=\begin{cases}1,&U\leq h(x),\\-1,&U>h(x),\end{cases}
\]
then \(\mathbb E_U s_U=h\).  Convexity of the KR norm gives
\[
 \|(h-\mu h)\mu\|_{KR}
 \leq \mathbb E_U\|(s_U-\mu s_U)\mu\|_{KR},
\]
so one of the signs does at least as well as \(h\).  Conversely signs are
admissible.  Write a sign as \(s=2\mathbf1_E-1\), and set
\(p=\mu(E)\), \(q=1-p\).  Then
\[
 (s-\mu s)\mu=2pq(\mu_E-\mu_{E^c}).
\]
Kantorovich--Rubinstein duality now turns (2) into (1).  \(\square\)

### Lemma 2 (existence and eikonal rigidity of an extremizer)

The supremum defining \(\mathcal O_1(\mu)\) is attained.  If \(f\) is an
extremizer, \(g=f-\mu f\), and \(\mathcal O_1(\mu)>0\), then
\[
 \mu(g=0)=0,
 \qquad |\nabla f|=1\quad\mu\text{-a.e.}                             \tag{3}
\]
Moreover, for \(E=\{g>0\}\), \(p=\mu(E)\), \(q=1-p\),
\[
 W_1(\mu_E,\mu_{E^c})=
 \frac{\mathcal O_1(\mu)}{2pq},                                    \tag{4}
\]
and \(g\) is an optimal Kantorovich potential.  Thus every optimal coupling
\(\pi\) is supported on pairs satisfying
\[
 g(x)-g(y)=|x-y|.                                                    \tag{5}
\]
If \(Z_0=\{g=0\}\), then in addition
\[
 \boxed{\quad |g(x)|=d(x,Z_0)\quad\text{for }\mu\text{-a.e. }x.\quad} \tag{5a}
\]
Thus an extremizer is not merely eikonal: it is a signed distance to its
zero separator on the support of the measure.

#### Proof

Subtract constants from a maximizing sequence so that \(f_j(0)=0\).  Then
\(|f_j(x)|\leq|x|\).  Arzela--Ascoli and a diagonal argument give a locally
uniformly convergent subsequence with a 1-Lipschitz limit \(f\).  Dominated
convergence, using the finite first moment, proves attainment.

Put \(A=\int|g|\,d\mu=\mathcal O_1(\mu)\), and initially take
\(E=\{g\geq0\}\).  Since
\[
 \int g_+\,d\mu=\int g_-\,d\mu=A/2,
\]
the admissible potential \(g\) gives
\[
 W_1(\mu_E,\mu_{E^c})\geq
 \mathbb E_{\mu_E}g-\mathbb E_{\mu_{E^c}}g
 =\frac{A}{2pq}.
\]
Formula (1) gives the reverse inequality, because \(A=\mathcal O_1(\mu)\).
Hence equality holds and every optimal coupling satisfies (5).

At every differentiability point \(x\), equality
\(g(x)-g(y)=|x-y|>0\) forces \(|\nabla g(x)|=1\).  To see this directly, set
\(v=(x-y)/|x-y|\).  The two Lipschitz inequalities on the two pieces of the
segment imply
\[
 g(x-tv)=g(x)-t\qquad(0\leq t\leq |x-y|),
\]
so \(\nabla g(x)\cdot v=1\), while \(|\nabla g(x)|\leq1\).
The same argument works at the other endpoint.

If \(\mu(g=0)>0\), the preceding construction assigns the zero level to
the source side.  Its points must be coupled to the strictly negative side,
so the equality argument gives \(|\nabla g|=1\) almost everywhere on
\(\{g=0\}\).  On the other hand, the gradient of a Lipschitz function is
zero almost everywhere on any one of its level sets.  This is a
contradiction.  Therefore the zero level has zero mass, and the endpoint
argument proves (3) on both sides.

Finally, equality in the Lipschitz inequality implies equality at every
point of the segment \([x,y]\).  That segment contains a point \(z\) with
\(g(z)=0\), and
\[
 |x-z|=g(x),\qquad |y-z|=-g(y).
\]
For any \(z'\in Z_0\), Lipschitzness gives
\(|g(x)|\leq|x-z'|\).  The calibrated point \(z\) gives the reverse
inequality.  Since the two marginals of \(\pi\) cover \(\mu\)-almost every
point, this proves (5a).  \(\square\)

The eikonal conclusion is genuinely restrictive, but a.e. regularity is
not enough: \(f(x)=|x|\) has \(|\nabla f|=1\) away from the origin.

### Lemma 3 (smooth eikonal functions are affine)

If \(f\in C^2(\mathbb R^n)\) and \(|\nabla f|=1\) everywhere, then \(f\) is
affine.

#### Proof

Let \(b=\nabla f\).  Differentiating \(|b|^2=1\) gives
\((\nabla b)b=0\).  Along the global flow of the bounded \(C^1\) vector
field \(b\), the vector \(b\) is therefore constant, so
\[
 \Phi_t(x)=x+t\nabla f(x).
\]
Every \(\Phi_t\) is a diffeomorphism, with inverse \(\Phi_{-t}\), and hence
\(I+t\nabla^2f(x)\) is invertible for every real \(t\).  Since the Hessian
is symmetric, a nonzero eigenvalue \(\lambda\) would make this derivative
singular at \(t=-1/\lambda\).  Thus \(\nabla^2f=0\).  \(\square\)

Accordingly, a nonlinear extremizer must use a singular/focal set from
which its nonbranching transport rays emanate.  Controlling that singular
set is one possible formulation of the unresolved step.

## 2. Constant-mass tails of a true extremizer

Let \(D=\mathcal O_1(\mu)\), let \(f\) be the extremizer from Lemma 2, set
\(g=f-\mu f\), and write
\[
 A=D,\qquad V=\int g^2\,d\mu,
 \qquad \kappa=V/A^2.
\]

Milman's convexity-based reverse concentration theorem gives, for every
fixed log-concave \(\mu\),
\[
 C_P(\mu)\leq C D^2.                                                 \tag{6}
\]
This is a quantitative equivalence theorem, not an assumed universal KLS
bound.  Since \(\int|\nabla f|^2d\mu=1\) by (3),
\[
 1\leq\kappa\leq C.                                                  \tag{7}
\]

### Lemma 4 (two macroscopic, separated tails)

There is a universal \(\delta>0\) such that
\[
 E_+=\{g\geq A/4\},\qquad E_-=\{g\leq-A/4\}
\]
satisfy
\[
 \mu(E_+),\mu(E_-)\geq\delta,
 \qquad \operatorname{dist}(E_+,E_-)\geq A/2.                       \tag{8}
\]

#### Proof

Paley--Zygmund applied separately to \(g_+\) and \(g_-\) gives
\[
 \mu(g\geq A/4)
 \geq\frac14\frac{(\mathbb Eg_+)^2}{\mathbb Eg_+^2}
 \geq\frac{A^2}{16V}=\frac1{16\kappa},
\]
and the same estimate holds for the negative side.  Use (7).  The distance
claim follows from Lipschitzness.  \(\square\)

Without the extremal comparison (6), constant tail masses are false.  In
one dimension take \(X=Z-1\) with \(Z\sim\operatorname{Exp}(1)\), and
\(f_s(x)=(x-s+1)_+\).  If \(\varepsilon=e^{-s}\), then
\[
 \mathbb Ef_s=\varepsilon,
 \quad \operatorname{Var}(f_s)=2\varepsilon-\varepsilon^2,
 \quad \mathbb E|f_s-\mathbb Ef_s|=2\varepsilon e^{-\varepsilon},
\]
so \(V/A^2\sim(2\varepsilon)^{-1}\to\infty\).

## 3. Barycentric, rank, and angular constraints

### Lemma 5 (conditional barycenters and covariance)

For any event \(B\) of mass \(s\in(0,1)\), if
\(b_B=\mathbb E[X\mid B]\) and \(\Sigma_B=\operatorname{Cov}(X\mid B)\),
then
\[
 |b_B|^2\leq\frac{1-s}{s},
 \qquad \Sigma_B\preceq \frac1s I.                                 \tag{9}
\]

#### Proof

Conditioning on the binary sigma-field generated by \(B\), total covariance
gives
\[
 I\succeq
 s b_Bb_B^T+(1-s)b_{B^c}b_{B^c}^T
 =\frac{s}{1-s}b_Bb_B^T,
\]
where \(b_{B^c}=-s b_B/(1-s)\).  This proves the first claim.  Also
\[
 I=s\,\mathbb E[XX^T\mid B]+(1-s)\mathbb E[XX^T\mid B^c]
 \succeq s\Sigma_B.
\]
\(\square\)

Apply this to the two tail events in (8).  Their barycenters are \(O(1)\),
even though their conditional laws are at \(W_1\)-distance at least \(A/2\).
The separation is therefore nonlinear rather than translational.

There is also a quantitative high-rank statement.  Let \(X_+\) and \(X_-\)
be independent samples from the two conditional tail laws, and put
\(Z=X_+-X_-\).  Then
\[
 |Z|\geq A/2\quad\text{a.s.},
 \qquad \|\operatorname{Cov}Z\|_{op}\leq2/\delta,
\]
whereas
\[
 \operatorname{tr}\operatorname{Cov}Z
 \geq A^2/4-4(1-\delta)/\delta.
\]
Consequently
\[
 \operatorname{srank}(\operatorname{Cov}Z)
 \geq \delta A^2/8-2(1-\delta).                                    \tag{10}
\]
Thus a bad witness requires separation spread over \(\Omega(A^2)\)
directions under independent tail sampling.

There is a complementary convexification constraint.  Since a conditional
barycenter belongs to the closed convex hull of its event, (9) gives
\[
 \operatorname{dist}\bigl(\overline{\operatorname{conv}}E_+,
                           \overline{\operatorname{conv}}E_-\bigr)
 \leq |b_{E_+}-b_{E_-}|\leq2\sqrt{(1-\delta)/\delta}=O(1),           \tag{10a}
\]
whereas the original sets are \(A/2\)-separated.  In particular, if both
tail sets were convex, then \(A=O(1)\).  A bad witness must therefore have
strongly nonconvex, angularly interlaced tails; convexification creates an
\(O(1)\) shortcut across an \(\Omega(A)\) metric gap.  This conclusion uses
only isotropy and is not an isoperimetric assertion.

### Lemma 6 (radial projection forces angular separation)

Define the thin-shell parameter
\[
 \tau^2=\mathbb E\bigl(|X|-\sqrt n\bigr)^2
\]
and let \(T(x)=\sqrt n\,x/|x|\) away from the origin.  For the tail laws in
Lemma 4,
\[
 W_1(T_\#\mu_{E_+},T_\#\mu_{E_-})
 \geq A/2-2\tau/\sqrt\delta.                                       \tag{11}
\]

#### Proof

For any event \(B\) of mass at least \(\delta\),
\[
 W_1(\mu_B,T_\#\mu_B)
 \leq\mathbb E[\,||X|-\sqrt n|\mid B]
 \leq\tau/\sqrt\delta.
\]
Use the triangle inequality and (8).  \(\square\)

Thus, once \(A\gg\tau\), the obstruction is necessarily angular.

### Lemma 6.1 (balanced transport needles)

Use the standard nonbranching transport-ray disintegration for the optimal
coupling in Lemma 2.  There are a quotient probability \(\lambda\) and
one-dimensional log-concave conditional measures \(\nu_\omega\), with
coordinate \(t=g\), such that
\[
 \mu=\int\nu_\omega\,d\lambda(\omega),
 \qquad
 \nu_\omega(t>0)=p,quad \nu_\omega(t<0)=q                         \tag{11a}
\]
for \(\lambda\)-almost every ray.  Furthermore, a universal positive
fraction of the rays satisfy
\[
 \int t_+\,d\nu_\omega\geq cA,qquad
 \int t_-\,d\nu_\omega\geq cA,                                   \tag{11b}
\]
and on each such ray both \(\{t\geq cA\}\) and \(\{t\leq-cA\}\)
have conditional mass at least \(c\).  All constants are dimension-free.

#### Proof

Outside the negligible branching set, let \(Q\) be the measurable quotient
map to maximal transport rays.  A calibrated pair has \(Q(x)=Q(y)\), so the
two marginals of the coupling give
\[
 Q_\#\mu_E=Q_\#\mu_F=\lambda.
\]
Disintegrating and using \(\mu=p\mu_E+q\mu_F\) proves (11a).  The fact that
the conditional densities are one-dimensional log-concave is the standard
localization/disintegration theorem for a log-concave measure along
nonbranching \(L^1\)-transport rays.  This theorem is used only here; it is
not a KLS estimate.

We record the elementary one-dimensional estimate needed next.  If a
log-concave density \(\varphi\) satisfies
\(\int_0^\infty\varphi=p\), \(\int_{-\infty}^0\varphi=q\), with
\(p,q\geq\delta\), and
\[
 m_+=\int_0^\infty t\varphi(t)dt,qquad
 m_-=\int_{-\infty}^0(-t)\varphi(t)dt,
\]
then
\[
 c_\delta m_+\leq m_-\leq C_\delta m_+,
 \qquad
 \int t_\pm^2\varphi(t)dt\leq C_\delta m_\pm^2.                    \tag{11c}
\]
Here is a proof.  Put \(h=\varphi(0)\) and \(H=\sup\varphi\).  If a mode
of height \(H\) lies at \(t_0>0\), concavity of \(\log\varphi\), with
\(a=\log(H/h)/t_0\), gives
\[
 p\geq(H-h)/a,qquad q\leq h/a;
\]
hence \(H/h\leq1+p/q=1/q\).  The case of a mode on the negative side is
symmetric, so \(H\leq h/\delta\).  A density bounded by \(H\) and having
mass \(p\) on the positive half-line satisfies
\(m_+\geq p^2/(4H)\): at least \(p/2\) of its mass lies beyond
\(p/(2H)\).  On the other hand, the survival function
\(S(t)=\int_t^\infty\varphi\) is log-concave, whence
\[
 S(t)\leq p\exp(-ht/p),
 \quad m_+\leq p^2/h,
 \quad \int_0^\infty t^2\varphi(t)dt
 =2\int_0^\infty tS(t)dt\leq2p^3/h^2.
\]
The same estimates hold on the negative side.  Since \(p,q\in[\delta,1]\),
they imply (11c).

Now set
\[
 s(\omega)=\int|t|\,d\nu_\omega.
\]
Then \(\int s\,d\lambda=A\), while conditional Jensen and (7) give
\[
 \int s^2d\lambda\leq\int g^2d\mu=V\leq C A^2.
\]
Paley--Zygmund shows that \(s(\omega)\geq A/2\) on a fixed positive
fraction of the rays.  Estimate (11c) makes the two one-sided first moments
comparable, proving (11b).  A final Paley--Zygmund application to \(t_+\)
and \(t_-\), using the second-moment part of (11c), proves the conditional
tail assertion.  \(\square\)

This rules out a loophole in which all long positive excursions live on
different rays from all long negative excursions.  It still does not bound
the common length: an abstract isotropic family can distribute long
two-sided tangent needles over many directions.

## 4. Exact transport rays are almost tangent

Let \(E=\{g>0\}\), \(F=\{g<0\}\), and let \(\pi\) be the optimal coupling
from Lemma 2, oriented from \(E\) to \(F\).  Equations (4)--(5), together
with the constant tail mass, imply that a fixed fraction of \(\pi\) consists
of rays of length between \(cA\) and \(C A\).  This can be made explicit.

Indeed, \(p,q\geq\delta\), and hence
\[
 \mathbb E_\pi|X-Y|=A/(2pq)\leq A/\delta.
\]
The source marginal gives mass at least \(\delta\) to \(\{g\geq A/4\}\),
and every pair starting there has length at least \(A/4\).  Markov's
inequality gives an upper bound \(8A/\delta^2\) outside a set of
\(\pi\)-mass at most \(\delta/8\).

Take
\[
 s=\sqrt8\,\tau/\delta.
\]
The two marginals together lose at most \(\delta/4\) upon restriction to
the annulus \(\{||x|-\sqrt n|\leq s\}\).  Thus a set of \(\pi\)-mass at
least \(\delta/2\) consists of pairs satisfying
\[
 A/4\leq|X-Y|\leq8A/\delta^2,
 \qquad ||X|-\sqrt n|,||Y|-\sqrt n|\leq s.                          \tag{12}
\]

For \(Z=X-Y\) and \(M=(X+Y)/2\),
\[
 \langle Z,M\rangle=(|X|^2-|Y|^2)/2.                               \tag{13}
\]
If additionally \(s\leq\sqrt n/4\) and
\(8A/\delta^2\leq\sqrt n\), then (12)--(13) give
\[
 \frac{|\langle Z,M\rangle|}{|Z|\,|M|}
 \leq C_\delta\frac{\tau}{A}.                                    \tag{14}
\]
Hence, if \(A\gg\tau\), a positive-mass family of exact, nonbranching
transport rays is almost tangent to the thin shell.  Formula (14) is the
most concrete geometric form of the remaining obstruction found here.

## 5. A midpoint exclusion mechanism, and its barrier

### Lemma 7 (thin-shell midpoint bound)

Let \(B,C\subset\mathbb R^n\) have \(\mu(B),\mu(C)\geq\delta\), and assume
\(\operatorname{dist}(B,C)\geq r\).  Then
\[
 r^2\leq64\left(
   \frac{\sqrt n\,\tau}{\sqrt\delta}
   +\frac{\tau^2}{\delta}
 \right).                                                          \tag{15}
\]

#### Proof

Set \(s=2\tau/\sqrt\delta\) and intersect both sets with
\[
 K_s=\{x:||x|-\sqrt n|\leq s\}.
\]
Chebyshev gives \(\mu(K_s^c)\leq\delta/4\), so the truncated sets
\(B',C'\) each have mass at least \(3\delta/4\).  Log-concavity of the
measure (the multiplicative Brunn--Minkowski inequality) gives
\[
 \mu\bigl((B'+C')/2\bigr)
 \geq\sqrt{\mu(B')\mu(C')}\geq3\delta/4.                            \tag{16}
\]
For \(z=(x+y)/2\) with \(x\in B'\), \(y\in C'\), the parallelogram
identity yields
\[
 |z|^2=\frac{|x|^2+|y|^2}{2}-\frac{|x-y|^2}{4}
 \leq(\sqrt n+s)^2-r^2/4.                                         \tag{17}
\]
Let the square root of the right side be \(\rho\).  From (16)--(17),
\(\mu(|X|\leq\rho)\geq3\delta/4\).  Chebyshev applied to the lower shell
tail therefore gives
\[
 \rho\geq \sqrt n-rac{2\tau}{\sqrt{3\delta}}.
\]
Substitution in (17), followed by expansion of the difference of squares,
proves (15) when the last lower bound is nonnegative.  If it is negative,
then \(\tau^2/\delta\gtrsim n\), and the trivial consequence
\(r\leq2(\sqrt n+s)\) gives the same bound after enlarging the numerical
constant to 64.  \(\square\)

Applying Lemma 7 to (8) gives the fully noncircular estimate
\[
 A\leq C\left(n^{1/4}\tau^{1/2}+\tau\right).                        \tag{18}
\]
This is a genuine exclusion mechanism using only log-concavity, isotropy,
Milman's fixed-measure equivalence, and the available thin-shell parameter.
It does not give a dimension-free bound.  The geometric reason is sharp:
a tangent chord of length \(A\) penetrates the radius only by order
\(A^2/\sqrt n\), so shell information cannot see lengths below the
\(n^{1/4}\) scale without extra angular input.

## 6. Why barycentric averaging does not repair the midpoint bound

One might average independent conditional samples.  Because the conditional
barycenters are \(O(1)\), such an average lies much farther inside the shell
than the deterministic midpoint estimate.  The missing inference would be
that the ambient measure gives appreciable mass near a typical independent
midpoint.  This inference is false even for the standard Gaussian, which has
optimal slicing, thin-shell, and KLS behavior.

Take
\[
 B=\{x_1\geq a\},\qquad C=\{x_1\leq-a\},
\]
where \(a>0\) is chosen so both have a fixed mass.  If
\(X_B,Y_C\) are independent conditional samples, then in coordinates
\(2,\dots,n\),
\[
 (X_B+Y_C)/2\sim N(0,I_{n-1}/2).
\]
Thus its radius in those coordinates is typically \(\sqrt{n/2}\), while a
standard Gaussian point has exponentially small probability of lying in a
fixed-width neighborhood of that radius.  At the same time the Minkowski
midpoint set \((B+C)/2\) is all of \(\mathbb R^n\).  Log-concavity controls
the mass of the midpoint *set*, not the distribution of independently
chosen midpoints.  No dimension-free density or entropy comparison between
that midpoint law and \(\mu\) can hold.

Likewise, a fixed-direction alignment claim is false.  For radial Gaussian
tail sets, both conditional barycenters vanish and the covariance of the
independent displacement is a scalar matrix, of stable rank \(n\).  Any
successful alignment statement must instead use local radial geometry and
the self-consistent, saturated transport-ray structure of Lemma 2.

## 7. The precise unresolved lemma and its KLS status

The preceding reductions leave the following candidate.

> **Tangential transport-ray exclusion.**  An isotropic log-concave measure
> cannot carry a fixed positive mass of nonbranching, saturated rays of a
> single 1-Lipschitz Kantorovich potential, with lengths in \([cA,CA]\),
> both endpoints in the thin shell, and angle \(O(\tau/A)\) from the local
> tangent space, unless \(A=O(1)\).

Lemmas 1--6 prove that this assertion would exclude every bad extremizer.
Therefore its unrestricted universal form is quantitatively equivalent to
T3 and hence, by Milman's reverse theorem, to KLS.  It must not be imported
as an auxiliary fact.  The only identified potentially noncircular route is
to exploit the singular/focal geometry forced by Lemma 3: a nonlinear
eikonal potential needs singularities, whereas long nonbranching tangent
rays force large reach and small curvature.  A new dimension-free estimate
relating the weighted size of that singular set to slicing or thin-shell
data would complete the argument.  No such estimate is proved here.

Other explicitly circular substitutions are:

1. a dimension-free upper bound for \(W_1(\mu_E,\mu_{E^c})\) for all
   fixed-mass cuts (this is exactly (1));
2. a dimension-free separation bound for arbitrary two fixed-mass sets
   (equivalent to one-scale concentration under log-concavity);
3. replacing the maximizing cut by a comparable halfspace;
4. a universal angular Poincare inequality for the conditional measures on
   typical spherical shells;
5. passing global isotropy to one-dimensional localization needles.
