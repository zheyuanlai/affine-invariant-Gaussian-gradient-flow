# Convex half-sets, mean distance, and the limits of convexification

## Executive conclusion

Let \(\mu\) be an isotropic log-concave probability on \(\mathbb R^n\),
let \(A\) be a closed convex set with \(\mu(A)=1/2\), and put
\[
 D=d(X,A),\qquad X\sim\mu,
 \qquad F(t)=\mu(A+tB_2^n)=\mathbb P(D\le t).             \tag{0.1}
\]
The first implication in the proposed route is true with an exact constant.
The function \(F\) is log-concave on \([0,\infty)\).  If
\[
 p=\mu^+(A)=F'_+(0),
\]
then
\[
 \boxed{
 \mathbb E D
 \ge \frac{2\log2-1}{4p}.
 }                                                        \tag{0.2}
\]
Consequently
\[
 \boxed{
 \mathbb E d(X,A)\le C
 \quad\Longrightarrow\quad
 \mu^+(A)\ge\frac{2\log2-1}{4C}.
 }                                                        \tag{0.3}
\]
The numerical constant is
\[
 \frac{2\log2-1}{4}=0.0965735\ldots .                   \tag{0.4}
\]
It is sharp if one assumes only log-concavity of \(F\); a one-dimensional
log-concave probability realizes equality.

As of July 2026, the universal distance bound is not a known consequence of
the thin-shell or slicing theorems.  Nor is there a known isotropic
log-concave family with \(\mathbb ED\to\infty\).  Such a family would be a
counterexample to KLS, because the distance function satisfies
\[
 \boxed{
 C_P(\mu)\ge(\mathbb ED)^2.
 }                                                        \tag{0.5}
\]
Indeed, \(D=0\) on half the mass and \(|\nabla D|\le1\) on the other
half.  Conversely KLS immediately gives the desired bound.  The current
dimension-dependent KLS theorem gives only
\[
 \mathbb Ed(X,A)\le C\sqrt{\log(n+1)}.                   \tag{0.6}
\]

The recently proved thin-shell theorem does settle the special case in
which \(A\) is a centered Euclidean ball, and half-spaces have an elementary
constant bound.  It does not control a general convex distance function.
The distinction is real: Section 4 constructs an explicit isotropic
absolutely continuous, but non-log-concave, family which simultaneously has

* uniformly bounded thin-shell width;
* uniformly bounded isotropic constant, the numerical conclusion of
  slicing; and
* a convex half-set \(A_n\) with
  \(\mathbb E d(X,A_n)\ge (17/40)\sqrt n\).

Thus the two deep theorems cannot be used as black boxes to prove the
distance claim.  A proof would need an additional consequence of
log-concavity which rules out the many-direction mixture in that model.

Finally, an arbitrary near-Cheeger witness cannot simply be convexified.
Even in one dimension:

* an arbitrarily near-optimal Gaussian witness can have convex hull equal
  to the entire line, changing its mass from \(1/2\) to \(1\);
* an exact Cheeger minimizer for the isotropic Laplace law can have convex
  hull equal to the entire line and stay at symmetric-difference distance
  at least \(1/2\) from every convex set of mass \(1/2\).  Smooth strictly
  log-concave approximations give arbitrarily near-optimal examples with the
  same separation.

A theorem replacing every balanced near-Cheeger witness by some balanced
convex set with comparable perimeter, when combined with (0.3), would
prove KLS.  It is therefore another load-bearing statement, not a harmless
geometric preprocessing step.

## 1. Log-concavity of parallel volume

For \(t\ge0\), write
\[
 A_t=A+tB_2^n.
\]
Convexity of \(A\) and of the Euclidean ball gives, for
\(s,t\ge0\) and \(0\le\theta\le1\),
\[
 (1-\theta)A_s+\theta A_t
 =A_{(1-\theta)s+\theta t}.                              \tag{1.1}
\]
The defining set inequality for a log-concave probability therefore yields
\[
 \begin{aligned}
 F((1-\theta)s+\theta t)
 &=\mu(A_{(1-\theta)s+\theta t})\\
 &\ge \mu(A_s)^{1-\theta}\mu(A_t)^\theta
 =F(s)^{1-\theta}F(t)^\theta.                            \tag{1.2}
 \end{aligned}
\]
Hence \(F\) is log-concave.  The argument extends from compact sets by
monotone approximation.  It also works intrinsically when \(\mu\) is
supported on a proper affine subspace.

Since \(F(0)=1/2\), \(F(t)\uparrow1\), and \(\log F\) is concave, its
right derivative exists.  If
\[
 p=F'_+(0)=\mu^+(A),                                     \tag{1.3}
\]
then
\[
 (\log F)'_+(0)=\frac{F'_+(0)}{F(0)}=2p.                 \tag{1.4}
\]
The tangent-line upper bound for a concave function gives
\[
 \log F(t)\le-\log2+2pt,
\]
or
\[
 F(t)\le\min\left\{1,\frac12e^{2pt}\right\}.            \tag{1.5}
\]
Here and below the displayed calculation is for \(0<p<\infty\).  The case
\(p=+\infty\) makes the desired perimeter lower bound automatic, while
\(p=0\) is impossible because (1.5) would contradict \(F(t)\uparrow1\).

## 2. The sharp mean-distance implication

The layer-cake identity gives
\[
 \mathbb ED=\int_0^\infty\mathbb P(D>t)dt
 =\int_0^\infty(1-F(t))dt.                               \tag{2.1}
\]
Let
\[
 t_0=\frac{\log2}{2p}.
\]
Using (1.5) until its exponential upper bound reaches one,
\[
 \begin{aligned}
 \mathbb ED
 &\ge\int_0^{t_0}\left(1-\frac12e^{2pt}\right)dt\\
 &=\frac{\log2}{2p}-\frac1{4p}
 =\frac{2\log2-1}{4p}.                                  \tag{2.2}
 \end{aligned}
\]
This proves (0.2)--(0.3).  Notice the direction: log-concavity gives an
upper tangent bound on \(F\), hence a lower bound on the mean distance for
a prescribed initial perimeter.

### 2.1 General starting mass

The same calculation is sometimes useful when \(q=\mu(A)\ne1/2\).  If
\(F(0)=q\) and \(p=F'_+(0)\), then
\[
 F(t)\le\min\{1,q e^{pt/q}\},
\]
and therefore
\[
 \boxed{
 \mathbb Ed(X,A)
 \ge\frac q p\left(\log\frac1q-1+q\right).
 }                                                        \tag{2.3}
\]
For \(q=1/2\), this is exactly (2.2).

### 2.2 Sharp one-dimensional model

Fix \(p>0\), put
\[
 L=\frac{\log2}{2p},
\]
and consider the log-concave density
\[
 \rho(x)=p e^{2px}\mathbf1_{(-\infty,L]}(x).             \tag{2.4}
\]
It integrates to one.  For \(A=(-\infty,0]\),
\[
 \mu(A)=\frac12,
 \qquad \mu^+(A)=\rho(0)=p,
\]
and
\[
 F(t)=
 \begin{cases}
  \frac12e^{2pt},&0\le t\le L,\\
  1,&t\ge L.
 \end{cases}                                             \tag{2.5}
\]
Thus equality holds in (2.2).  Centering and scaling this one-dimensional
law to isotropic position does not change the product
\(\mu^+(A)\mathbb Ed(X,A)\), so the constant remains sharp in isotropic
normalization.

## 3. Relation to Poincare and KLS

For a closed convex \(A\), the function \(D=d(\cdot,A)\) is convex and
one-Lipschitz.  It is differentiable almost everywhere, with
\[
 |\nabla D|=1\quad\text{almost everywhere on }A^c        \tag{3.1}
\]
for a full-dimensional absolutely continuous law.  The inequality
\(|\nabla D|\le1\) is enough below.

Put \(M=\mathbb ED\).  Since \(D=0\) on a set of mass \(1/2\),
\[
 \operatorname{Var}_\mu(D)
 =\int(D-M)^2d\mu
 \ge\int_A M^2d\mu
 =\frac{M^2}{2}.                                        \tag{3.2}
\]
Also
\[
 \int|\nabla D|^2d\mu\le\mu(A^c)=\frac12.              \tag{3.3}
\]
It follows directly from the definition of the Poincare constant that
\[
 C_P(\mu)\ge
 \frac{\operatorname{Var}_\mu(D)}
      {\int|\nabla D|^2d\mu}
 \ge M^2.                                                \tag{3.4}
\]

Thus:

1. KLS, namely \(C_P(\mu)\le C\), implies
   \(\mathbb Ed(X,A)\le\sqrt C\) for every half-mass set, convex or not.
2. An explicit isotropic log-concave convex-half-set family with
   \(\mathbb Ed\to\infty\) would give an explicit KLS counterexample.
3. The known theorem \(C_P(\mu)\le C\log(n+1)\) gives (0.6).

The desired distance statement is weaker than full KLS because it tests
only convex distance functions.  No converse from it to full Poincare is
known without a valid convexification theorem for isoperimetric witnesses.

There is also a useful functional interpretation.  If \(f\) is convex and
one-Lipschitz, \(m\) is a median, and \(A=\{f\le m\}\), then
\[
 (f(x)-m)_+\le d(x,A).                                  \tag{3.5}
\]
Indeed, for every \(y\in A\), convexity is not even needed for the estimate
\(f(x)-m\le f(x)-f(y)\le|x-y|\).  Conversely, distance to a convex set is
a convex one-Lipschitz function whose median is zero when the set has mass
one half.  Up to the harmless issue of a median level set having excess
mass (cut it inside by a half-space), the target is therefore the one-copy
\(L^1\) upper-deviation inequality for convex one-Lipschitz functions.  The
thin-shell theorem supplies this for the single function \(x\mapsto|x|\),
not for this entire class.

## 4. What thin-shell and slicing do and do not give

The thin-shell conjecture was proved by Klartag and Lehec in 2025: for an
isotropic log-concave \(X\),
\[
 \mathbb E(|X|-\sqrt n)^2\le C_{\rm TS}.                 \tag{4.1}
\]
Bourgain's slicing conjecture was proved by Klartag and Lehec, building on
Guan, in 2024.  In density language, it gives a universal bound on
the isotropic constant.  Neither theorem states concentration for every
convex one-Lipschitz function.

### 4.1 Centered balls

Let \(r\) be a median of \(|X|\), and let \(A=rB_2^n\).  From (4.1),
\[
 |r-\sqrt n|\le\sqrt{2C_{\rm TS}},                       \tag{4.2}
\]
because at least half the mass lies on the appropriate side of the median.
Therefore
\[
 \begin{aligned}
 \mathbb E d(X,rB_2^n)
 &=\mathbb E(|X|-r)_+\\
 &\le\mathbb E\big||X|-\sqrt n\big|+|r-\sqrt n|\\
 &\le(1+\sqrt2)\sqrt{C_{\rm TS}}.                       \tag{4.3}
 \end{aligned}
\]
Thus the new theorem proves the target for centered Euclidean balls.

### 4.2 Half-spaces

Let \(Z=\langle X,u\rangle\), \(|u|=1\), and let \(m\) be a median of
\(Z\).  Isotropy gives \(\mathbb EZ=0\) and \(\mathbb EZ^2=1\).  Since at
least half the mass lies beyond a median on either relevant side,
\[
 |m|\le\sqrt2.                                           \tag{4.4}
\]
For the half-space \(A=\{x:\langle x,u\rangle\le m\}\),
\[
 \mathbb Ed(X,A)
 =\mathbb E(Z-m)_+
 \le\sqrt{\mathbb E(Z-m)^2}
 \le\sqrt3.                                             \tag{4.5}
\]
No high-dimensional theorem is needed in this case.

### 4.3 A black-box obstruction to thin-shell plus slicing

The following family is not log-concave.  Its purpose is to prove that the
numerical conclusions of thin-shell and slicing alone cannot yield the
desired distance bound.

For \(n\ge2\), fix \(0<\varepsilon\le1/20\); one may take
\(\varepsilon=1/20\) throughout.
Let \(U\) be uniform on the Euclidean ball of
radius
\[
 r=\varepsilon\sqrt{n+2}.
\]
Let \(I\) be uniform on \(\{1,\ldots,n\}\), let \(S\) be an independent
uniform sign, and set
\[
 R=\sqrt{n(1-\varepsilon^2)},
 \qquad X=SR e_I+U.                                      \tag{4.6}
\]
Since
\[
 \operatorname{Cov}(U)=\varepsilon^2I,
 \qquad
 \mathbb E[R^2e_Ie_I^T]=(1-\varepsilon^2)I,
\]
the random vector \(X\) is centered and isotropic.

Its density is the equal mixture of \(2n\) translated uniform-ball
densities.  The components are disjoint for the chosen \(\varepsilon\),
and
\[
 \|p_X\|_\infty^{1/n}
 =\frac{(2n)^{-1/n}|B_2^n|^{-1/n}}
        {\varepsilon\sqrt{n+2}}
 \le\frac C\varepsilon.                                \tag{4.7}
\]
Thus its isotropic constant is bounded independently of \(n\), exactly the
kind of numerical bound supplied by slicing.

It also has a dimension-free thin shell.  The exact radial moments of a
uniform ball give
\[
 \operatorname{Var}(|U|^2)
 =\frac{4n}{n+4}\varepsilon^4\le4\varepsilon^4.          \tag{4.8}
\]
By symmetry the cross covariance vanishes, and hence
\[
 \begin{aligned}
 \operatorname{Var}(|X|^2)
 &=\operatorname{Var}(|U|^2)+4R^2\mathbb EU_I^2\\
 &\le4\varepsilon^4+4n\varepsilon^2.                    \tag{4.9}
 \end{aligned}
\]
Since \(\mathbb E|X|^2=n\),
\[
 \mathbb E(|X|-\sqrt n)^2
 \le\frac1n\operatorname{Var}(|X|^2)
 \le5\varepsilon^2.                                     \tag{4.10}
\]
For the concrete choice \(\varepsilon=1/20\), the last bound is \(1/80\).

Now define the closed convex set
\[
 A_n=\operatorname{conv}\{Re_1,\ldots,Re_n\}+rB_2^n.    \tag{4.11}
\]
It contains every positive component \(Re_i+rB_2^n\).  For \(n\ge2\), the Euclidean
distance from \(-Re_i\) to the positive simplex is
\[
 d\left(-Re_i,\operatorname{conv}\{Re_j\}_{j=1}^n\right)
 =R\sqrt{\frac n{n-1}}.                                 \tag{4.12}
\]
For \(\varepsilon\le1/20\), this is greater than \(2r\), so \(A_n\) is
disjoint from every negative component.  Consequently
\[
 \mathbb P(X\in A_n)=\frac12.                            \tag{4.13}
\]
For \(X=-Re_i+U\), the triangle inequality gives
\[
 d(X,A_n)
 \ge R\sqrt{\frac n{n-1}}-2r
 \ge\left(\sqrt{1-\varepsilon^2}-2\sqrt2\,\varepsilon\right)\sqrt n
 \ge\frac{17}{20}\sqrt n.                              \tag{4.14}
\]
It follows that
\[
 \boxed{
 \mathbb Ed(X,A_n)\ge\frac{17}{40}\sqrt n.
 }                                                        \tag{4.15}
\]

The density can be smoothed inside slightly enlarged component balls, with
an arbitrarily small covariance correction, without changing any estimate.
What cannot be done is make this separated mixture log-concave.  This is
precisely the extra structure a proof must exploit.

## 5. Canonical positive examples

The distance bound is automatic in several standard classes.

1. **Gaussian measures.**  The Gaussian Poincare constant is one, so
   (3.2)--(3.3) give \(\mathbb Ed(X,A)\le1\) for every half-mass set.
2. **Product log-concave measures.**  Tensorization of the one-dimensional
   Poincare inequality gives a numerical constant after isotropic
   normalization, hence a numerical distance bound.
3. **Strongly log-concave measures.**  If \(D^2V\succeq\kappa I\),
   Brascamp--Lieb gives \(\mathbb Ed(X,A)\le\kappa^{-1/2}\).
4. **Half-spaces and centered balls.**  These are covered by
   (4.3)--(4.5).

The cube, product exponential laws, Euclidean balls, and regular product
models therefore do not furnish counterexamples.  Any genuine divergent
family would have to exhibit the same unresolved global bottleneck as a
KLS counterexample.

## 6. Convexification of near-Cheeger witnesses

There are several inequivalent meanings of “convexify.”  Two natural ones
already fail in one dimension.

### 6.1 Convex hull can destroy balance for an arbitrarily near minimizer

Let \(\gamma\) be the standard Gaussian measure on \(\mathbb R\), with
density \(\varphi\) and distribution function \(\Phi\).  For
\(0<\delta<1/4\), choose
\[
 a_\delta=\Phi^{-1}(1/2-\delta),
 \qquad
 b_\delta=\Phi^{-1}(1-\delta),
\]
and put
\[
 E_\delta=(-\infty,a_\delta]\cup[b_\delta,\infty).       \tag{6.1}
\]
Then
\[
 \gamma(E_\delta)=\frac12,
 \qquad
 \gamma^+(E_\delta)=\varphi(a_\delta)+\varphi(b_\delta).
                                                                    \tag{6.2}
\]
Gaussian isoperimetry says the optimal balanced perimeter is
\(\varphi(0)\).  As \(\delta\downarrow0\),
\[
 \varphi(a_\delta)\longrightarrow\varphi(0),
 \qquad
 \varphi(b_\delta)\sim
 \delta\sqrt{2\log(1/\delta)}\longrightarrow0.          \tag{6.3}
\]
Hence
\[
 \frac{\gamma^+(E_\delta)}{\varphi(0)}\longrightarrow1. \tag{6.4}
\]
These are arbitrarily near-Cheeger balanced witnesses.  Nevertheless,
\[
 \operatorname{conv}(E_\delta)=\mathbb R,
 \qquad
 \gamma(\operatorname{conv}E_\delta)=1.                 \tag{6.5}
\]
Thus taking the convex hull can lose half the mass even in the most rigid
isoperimetric model.

There does exist a different half-line within symmetric-difference
\(2\delta\) of \(E_\delta\).  This distinction matters: (6.5) refutes the
literal convex-hull operation, not every possible stability-based
replacement theorem.

### 6.2 An exact Cheeger minimizer can be far from every convex half-set

Let \(\mu\) be the isotropic Laplace law on \(\mathbb R\), whose density is
\[
 f(x)=\frac1{\sqrt2}e^{-\sqrt2|x|},
 \qquad a=\frac{\log2}{\sqrt2},                          \tag{6.6}
\]
and set
\[
 E=(-\infty,-a]\cup[a,\infty).                           \tag{6.7}
\]
Each tail has mass \(1/4\), so \(\mu(E)=1/2\), and
\[
 \mu^+(E)=f(-a)+f(a)=\frac1{\sqrt2}.                    \tag{6.8}
\]
For this law
\(f(x)=\sqrt2\min\{F(x),1-F(x)\}\).  The standard one-dimensional
boundary formula therefore gives Cheeger constant \(\sqrt2\).  In
particular, every balanced set has perimeter at least \(1/\sqrt2\), so
\(E\) is an exact Cheeger minimizer, tied with either median half-line.
For completeness, on an interval component lying to one side of zero,
the sum of its endpoint densities is at least \(\sqrt2\) times its mass.
For a component crossing zero and having mass at most \(1/2\), that endpoint
sum equals \(\sqrt2\) times the mass of its complement and is again at least
\(\sqrt2\) times its own mass.  Summing components and approximating proves
the assertion for arbitrary finite-perimeter sets of mass at most \(1/2\).

On the other hand, \(\operatorname{conv}(E)=\mathbb R\).  Moreover, a
convex subset \(C\subseteq\mathbb R\) of mass \(1/2\) is an interval or ray.
If it met both tails in positive measure, it would contain the entire
central interval \((-a,a)\), which already has mass \(1/2\), a
contradiction.  Thus \(\mu(E\cap C)\le1/4\), and
\[
 \boxed{
 \mu(E\triangle C)\ge\frac12.
 }                                                        \tag{6.9}
\]
The same interval argument gives \(\mu(E\triangle C)\ge1/4\) whenever
\(\mu(C)\in[1/4,3/4]\).

This is not an artifact of the cusp of the Laplace density.  The smooth,
strictly log-concave densities proportional to
\(\exp(-\sqrt{x^2+\eta^2})\), put in isotropic position, converge to the
isotropic Laplace law as \(\eta\downarrow0\).  Taking the union of their two
outer quartile tails gives balanced sets whose perimeter divided by the
optimal balanced perimeter tends to one, while (6.9) remains valid exactly.
Thus even smooth strict log-concavity does not support a general
symmetric-difference convexification of near-Cheeger witnesses.

### 6.3 An abstract replacement theorem is KLS-strength

Suppose one had both of the following numerical statements.

1. Every convex \(A\) of mass \(1/2\) satisfies
   \(\mathbb Ed(X,A)\le C_0\).
2. Every balanced near-Cheeger witness \(E\) can be replaced by a convex
   balanced \(A\) with
   \[
    \mu^+(A)\le K\mu^+(E).                               \tag{6.10}
   \]

By (0.3),
\[
 \mu^+(E)
 \ge\frac{2\log2-1}{4KC_0}.                             \tag{6.11}
\]
The concavity of the isoperimetric profile permits the Cheeger constant to
be tested at mass \(1/2\), so (6.11) is a dimension-free Cheeger inequality,
and hence KLS.

Therefore a nonlocal replacement theorem of the form (6.10) is not ruled
out, but it is itself a principal missing theorem.  The elementary examples
above show that neither convex hull nor symmetric-difference approximation
provides it automatically.

## 7. Final status

The audited conclusions are:

1. The log-concavity argument is exact and yields the sharp constant
   \((2\log2-1)/4\).
2. As of July 2026, the universal convex-half-set distance bound remains
   unresolved by the cited results.  A
   divergent isotropic log-concave example would disprove KLS, so none can
   be honestly supplied from known examples.
3. Current Poincare/KLS technology gives
   \(\mathbb Ed\lesssim\sqrt{\log(n+1)}\).
4. The solved thin-shell theorem handles centered balls; isotropy handles
   half-spaces.  Neither the thin-shell nor slicing conclusion alone handles
   arbitrary convex \(A\).
5. The many-direction mixture (4.6) satisfies both black-box conclusions
   while having \(\mathbb Ed\asymp\sqrt n\); its failure of log-concavity is
   the exact point at which a new argument is needed.
6. Naive convexification of near-Cheeger witnesses is false, and a general
   perimeter-preserving balanced replacement theorem would be KLS-strength.

### Primary references for the current theorem status

* B. Klartag and J. Lehec, [*Thin-shell bounds via parallel
  coupling*](https://arxiv.org/abs/2507.15495), arXiv:2507.15495 (2025,
  revised 2026).
* B. Klartag and J. Lehec, [*Affirmative Resolution of Bourgain's Slicing
  Problem Using Guan's Bound*](https://arxiv.org/abs/2412.15044),
  arXiv:2412.15044 (2024).
* B. Klartag, [*Logarithmic bounds for isoperimetry and slices of convex
  sets*](https://arxiv.org/abs/2303.14938), arXiv:2303.14938 (2023), for the
  dimension-dependent KLS bound used in (0.6).
