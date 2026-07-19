# Sequential slab localization: exact transfer and the unavoidable needle

## Executive conclusion

There is a clean non-Gaussian localization mechanism, but it does not close
the dimension-free argument.

* Any sequential one-dimensional observation has an exact Bayesian
  disintegration.  If \(g_j=\mu_j(S)\), then \(g_j\) is a martingale, its
  expected binary-entropy loss is exactly the information learned about
  \(1_S\), and ordinary posterior perimeter averages back to the original
  perimeter.  For hard convex partitions the average posterior perimeter is
  at most the original one, with equality for generic cuts.  This is the
  useful direction.
* Additive log-concave noise, interval quantization of such noise, and hard
  interval or halfspace conditioning preserve log-concavity for every
  outcome.  Nonconvex quantization cells need not do so.
* A hard slab of width \(h\) gives the pathwise bound
  \(\operatorname{Var}(\langle u,X\rangle\mid Y)\le h^2/4\), but it adds no
  interior curvature and says nothing about transverse directions.  A merely
  log-concave noisy likelihood does not even give a pathwise variance bound:
  exponential prior plus matched Laplace noise has posterior variance of
  order \(y^2\) after a large outcome \(y\).
* A strongly log-concave one-dimensional noise likelihood does add rank-one
  curvature.  Cycling can make the total curvature positive definite only
  after a total directional precision of order \(d\).  In the weak-signal,
  many-round regime, every sufficiently regular noisy or quantized channel
  converges to the Gaussian likelihood experiment by local asymptotic
  normality.  Thus this regime is Gaussian stochastic localization in new
  notation, not a genuinely different endpoint.
* There is an especially attractive hard construction.  At every node choose
  a hyperplane which simultaneously bisects
  \(\mathbf 1_S\mu_j\) and \(\mathbf 1_{S^c}\mu_j\), and reveal a random side.
  Both posterior branches are log-concave, \(g_j\) stays exactly equal to
  \(g_0\), the observation reveals zero information about \(1_S\), and
  perimeter transfers exactly off codimension-two intersections.  However,
  this is precisely the classical two-constraint ham-sandwich localization.
  It can eliminate all but one affine direction, but in the last direction a
  common bisector need not exist.  The limiting objects are balanced
  one-dimensional needles with uncontrolled variance.

For an isotropic starting law, disintegration only gives

\[
 \mathbb E[\sigma_\omega^2\theta_\omega\theta_\omega^T]\preceq I,
 \qquad \mathbb E\sigma_\omega^2\le d.                 \tag{0.1}
\]

One-dimensional log-concave isoperimetry then recovers only the
\(d^{-1/2}\) bound.  Proving that the balanced needles have
\(\mathbb E\sigma_\omega^2=O(1)\), or uniformly bounded diameter, is already
the missing dimension-free isoperimetric assertion.  The memoryless upper
tail of a one-sided exponential shows concretely why no branchwise variance
contraction follows from bisection.

The cube, simplex, product exponential, radial exponential, Gaussian
coordinate-halfspace, and Gaussian parity tests all exhibit the same
phenomenon: scalar observations can remove or curve the directions they see,
but either spend exponential mass/information to control the full spectrum or
leave a quantitatively uncontrolled exceptional direction.  No KLS or
ball-walk expansion is used below.

## 1. Exact disintegration for an arbitrary sequential scalar channel

Let \(X\sim\mu\), put \(B=\mathbf 1_S(X)\), and let
\(\mathcal F_{j-1}\) be the previous observation history.  Predictably choose
a unit vector \(u_j\) and a Markov kernel

\[
 K_j(dy\mid x)=\ell_j(y,\langle u_j,x\rangle)\,\lambda_j(dy),
 \qquad
 \int \ell_j(y,z)\,\lambda_j(dy)=1.                  \tag{1.1}
\]

Writing \(\nu=\mu_{j-1}\), the predictive density and posterior are

\[
 q(y)=\int \ell_j(y,\langle u_j,x\rangle)d\nu(x),
 \qquad
 d\nu_y(x)=\frac{\ell_j(y,\langle u_j,x\rangle)}{q(y)}d\nu(x).
                                                               \tag{1.2}
\]

Consequently, for every integrable \(f\),

\[
 \mathbb E[\nu_Y(f)\mid\mathcal F_{j-1}]=\nu(f).       \tag{1.3}
\]

This is the exact martingale/disintegration; no diffusion limit is involved.
In particular \(g_j=\mu_j(S)\) is a bounded martingale.

There is also an exact binary experiment at each step.  Put \(g=\nu(S)\) and

\[
 r_b(y)=\mathbb E_\nu[\ell_j(y,\langle u_j,X\rangle)\mid B=b],
 \qquad b\in\{0,1\}.
\]

Then

\[
 q=gr_1+(1-g)r_0,
 \qquad
 g_y=\frac{g r_1(y)}{q(y)},
 \qquad
 \log\frac{g_y}{1-g_y}
 =\log\frac g{1-g}+\log\frac{r_1(y)}{r_0(y)}.        \tag{1.4}
\]

The conditional quadratic increment is

\[
 \mathbb E[(g_Y-g)^2\mid\mathcal F_{j-1}]
 =g^2(1-g)^2\int\frac{(r_1-r_0)^2}{q}\,d\lambda_j.   \tag{1.5}
\]

If \(h(s)=-s\log s-(1-s)\log(1-s)\), the exact information cost is

\[
 h(g)-\mathbb E[h(g_Y)\mid\mathcal F_{j-1}]
 = I(B;Y\mid\mathcal F_{j-1}).                       \tag{1.6}
\]

After any finite adaptive sequence, the chain rule gives

\[
 \mathbb E h(g_m)=h(g_0)-I(B;Y_1,\ldots,Y_m).         \tag{1.7}
\]

Thus information about \(X\) is not the relevant mass budget; only
information about the one bit \(B\) moves the posterior set mass.  If
\(g_0=1/2\), then for \(0<\delta<1/2\),

\[
 \mathbb P\{\delta\le g_m\le1-\delta\}
 \ge
 \left(1-\frac{I(B;Y_{1:m})}{\log2-h(\delta)}\right)_+ .       \tag{1.8}
\]

Indeed \(h(g_m)\le h(\delta)\) outside the displayed interval.  The Bayes
error \(e_m=\mathbb E\min(g_m,1-g_m)\) also obeys

\[
 e_m\ge h_{[0,1/2]}^{-1}(\log2-I(B;Y_{1:m})),         \tag{1.9}
\]

because \(\mathbb E h(g_m)\le h(e_m)\).  These bounds are sharp as
information statements.  They give no automatic bound on the information
that a full-spectrum localization learns about \(B\).

For a deterministic hard partition, a realized cell \(C\) has posterior
\(\mu(\cdot\mid C)\), path probability \(\mu(C)\), and realized surprise
\(-\log\mu(C)\).  The entropy of the random cell is the total information
learned about \(X\); its contribution concerning \(S\) is still exactly
(1.7).

## 2. Posterior perimeter transfers in the correct direction

Suppose first that \(\mu\) has density \(p\), \(S\) has locally finite
perimeter, and the likelihoods in (1.1) are continuous at
\(\mathcal H^{d-1}\)-almost every point of \(\partial^*S\).  Ordinary weighted
perimeter gives

\[
 P_{\mu_y}(S)
 =\frac1{q(y)}\int_{\partial^*S}
       \ell(y,\langle u,x\rangle)p(x)d\mathcal H^{d-1}(x).
\]

Tonelli and the normalization in (1.1) yield the exact identity

\[
 \boxed{\quad \int q(y)P_{\mu_y}(S)d\lambda(y)=P_\mu(S).\quad} \tag{2.1}
\]

This remains true step by step for predictable adaptive directions.

For hard conditioning, let \(\{C_a\}\) be a countable convex partition and
\(q_a=\mu(C_a)\).  Perimeter relative to the posterior support satisfies

\[
 \sum_a q_a P_{\mu(\cdot\mid C_a)}(S)
 =\int_{\partial^*S\setminus\bigcup_a\partial C_a}
       p\,d\mathcal H^{d-1}
 \le P_\mu(S).                                      \tag{2.2}
\]

Equality holds whenever the cell walls carry no \(S\)-perimeter.  A
continuous random shift of a slab grid, or a generic random cut, has equality
almost surely.  If a cell wall coincides with a face of \(S\), that face
disappears from the relative posterior perimeter, explaining the inequality.
The direction in (2.2) is useful:

\[
 P_\mu(S)\ge \mathbb E P_{\mu_Y}(S).
\]

Therefore, if every posterior had Cheeger scale \(a_Y\), then

\[
 P_\mu(S)\ge c\,\mathbb E[a_Y\min(g_Y,1-g_Y)].       \tag{2.3}
\]

The transfer step is not the obstruction.  The obstruction is obtaining a
dimension-free posterior scale while keeping the right side of (2.3)
nontrivial.

## 3. Which scalar outcomes preserve log-concavity?

The following closure facts are exact.

1. If \(k=e^{-W}\) is a log-concave noise density, then the additive-noise
   likelihood
   
   \[
   \ell_y(z)=k(y-z)
   \]
   
   is log-concave in \(z\) for every \(y\).  Thus \(X\mid X+Z=y\) is
   log-concave whenever \(X\) is.
2. If the observed outcome is that \(X+Z\) lies in an interval \(I\), then
   
   \[
   \ell_I(z)=\int_I k(y-z)dy
   \]
   
   is log-concave by Prékopa.  Half-lines are allowed, so both outcomes of
   a noisy threshold preserve log-concavity.
3. A noiseless interval or half-line has likelihood \(\mathbf1_I(z)\), which
   is log-concave in the extended-valued sense.

Convexity of the output cell is essential.  For example, a two-tail output
\(\{|X+Z|\ge a\}\) is a union of intervals; even for Gaussian \(Z\), its
likelihood is symmetric and increases away from zero and hence is not
log-concave.

Log-concavity alone does not imply variance contraction.  Here is an exact
one-dimensional counterexample.  Let \(T\sim\operatorname{Exp}(1)\) and let
\(Z\) have the rate-one Laplace density.  For an outcome \(Y=T+Z=y>0\), the
unnormalized posterior density is

\[
 e^{-t-|y-t|}\mathbf1_{t\ge0}
 =
 \begin{cases}
 e^{-y},&0\le t\le y,\\
 e^{y-2t},&t\ge y.
 \end{cases}                                        \tag{3.1}
\]

Its normalizer is \(e^{-y}(y+1/2)\), and

\[
 \\operatorname{Var}(T\mid Y=y)\sim y^2/12.           \tag{3.2}
\]

Thus a perfectly valid log-concave outcome can increase the observed
variance without bound.  Fixed-width quantization of large \(y\)'s retains
the same long nearly flat posterior segment.

There are only two robust ways around (3.2).

* A hard slab \(a\le\langle u,x\rangle\le a+h\) gives
  
  \[
  \\operatorname{Var}(\langle u,X\rangle\mid\text{slab})\le h^2/4,          \tag{3.3}
  \]
  
  but contributes no Hessian curvature in the slab interior.
* If \(W''\ge\kappa>0\), an additive-noise outcome adds at least
  \(\kappa uu^T\) to the posterior Hessian.  After observations in directions
  \(u_j\) with strengths \(\kappa_j\), the added curvature is
  
  \[
  Q=\sum_j\kappa_j u_ju_j^T.                         \tag{3.4}
  \]

  Hence \(Q\succeq cI\) forces
  
  \[
  \sum_j\kappa_j=\operatorname{tr}Q\ge cd             \tag{3.5}
  \]
  
  and requires spanning all \(d\) directions.  Adaptation and cycling do not
  evade this trace identity.

For exact zero information about \(S\), (1.4) requires \(r_1=r_0\) almost
everywhere.  In an additive channel whose noise characteristic function has
no zeros, this is equivalent to equality of the full one-dimensional laws

\[
 \mathcal L(\langle u,X\rangle\mid S)
 =\mathcal L(\langle u,X\rangle\mid S^c).             \tag{3.6}
\]

A common median supplies only one binary equality, not (3.6).  This is why a
hard common-bisecting halfspace can be information-free while a continuum of
noisy outcomes in the same direction generally is not.

## 4. Random slabs and the full-diameter cost

A canonical hard experiment illustrates both the benefit and its cost.
Choose a public shift \(\Theta\sim\operatorname{Unif}[0,h\)) and reveal

\[
 K=\left\lfloor\frac{\langle u,X\rangle-\Theta}{h}\right\rfloor .          \tag{4.1}
\]

Conditional on \((\Theta,K)\), the posterior is the original log-concave law
restricted to a slab of width \(h\).  Equations (1.3), (2.2), and (3.3) are
all exact.  Repeating in directions \(u_1,\ldots,u_m\) restricts the posterior
to an intersection of slabs.  If all widths equal \(h\), two points \(x,x'\)
in one cell satisfy

\[
 (x-x')^T\left(\sum_{j=1}^m u_ju_j^T\right)(x-x')\le mh^2.                \tag{4.2}
\]

For an orthonormal cycle this only gives diameter \(h\sqrt d\).  Random
polynomial-size frames have the same scale up to logarithms; hard slabs do
not secretly generate strong convexity.

There is a direction-free information obstruction to forcing very small
cells.  If a density \(f\) is bounded and a conditioning event \(A\) has
Euclidean diameter \(D\), then

\[
 \mu(A)\le \|f\|_\infty\operatorname{vol}B_2^d(D/2).                     \tag{4.3}
\]

For the isotropic cube, simplex, product exponential, radial exponential,
and Gaussian, Stirling gives, with a numerical \(C\),

\[
 \mu(A)\le\left(\frac{CD}{\sqrt d}\right)^d.           \tag{4.4}
\]

Thus a path ending in a cell of diameter \(O(1)\) has realized information
cost at least

\[
 -\log\mu(A)\ge \tfrac d2\log d-O(d).                 \tag{4.5}
\]

This does not by itself imply loss of \(S\)-balance: a specially chosen
partition can learn many bits about \(X\) and zero bits about \(B\).  It does
show that diameter localization is a high-information, many-cut operation,
not a constant-time rank-one effect.

The Gaussian gives a complementary covariance calculation.  If
\(q=\gamma_d(\cdot\mid A)\) and \(\operatorname{Cov}_q(X)\preceq\rho I\), the
Gaussian maximum-entropy inequality yields

\[
 \log\frac1{\gamma_d(A)}=D(q\|\gamma_d)
 \ge\frac d2(\rho-1-\log\rho).                       \tag{4.6}
\]

For a bounded filter \(0\le L\le1\), with
\(q=L\gamma_d/Z\), the same lower bound holds for \(-\log Z\), since
\(D(q\|\gamma_d)=\mathbb E_q\log L-\log Z\le-\log Z\).  Constant full-spectrum
covariance reduction therefore costs \(\exp[-\Theta(d)]\) survival on this
basic model.

## 5. Many gentle non-Gaussian observations return to the Gaussian channel

The diffusive weak-signal regime has a universal local limit.  Let
\(k=e^{-W}\) be a positive \(C^3\) density with finite Fisher information

\[
 J=\int ((\log k)'(y))^2k(y)dy,
\]

and enough domination for third-order Taylor remainders.  Observe independent
variables with density \(k(y-\varepsilon z)\), where \(z=\langle u,x\rangle\)
and \(N_\varepsilon\varepsilon^2\to t\).  Uniformly for \(z\) in compact
sets, the cumulative log likelihood relative to \(z=0\) is

\[
 \sum_{i=1}^{N_\varepsilon}
 \log\frac{k(Y_i-\varepsilon z)}{k(Y_i)}
 \ \Longrightarrow\ 
 \sqrt{tJ}\,Gz-\frac{tJ}{2}z^2,                     \tag{5.1}
\]

where \(G\sim N(0,1)\).  The proof is just Taylor expansion, the score CLT,
and
\(\mathbb E_k(\log k)''=-J\).  A finite interval quantization has the same
statement with \(J\) replaced by the Fisher information of the discrete
output probabilities, provided these probabilities are positive and
differentiable in quadratic mean.

With predictable changing directions, (5.1) gives the corresponding
rank-one Gaussian likelihood martingale and accumulated quadratic form.
Therefore splitting a finite-strength log-concave observation into many
weak noisy or quantized observations does not remove Gaussian posterior
tilts; it converges to them.  To remain genuinely non-Gaussian one must use
finite-strength or hard outcomes, where the information and exceptional-tail
obstructions above are present.

## 6. Exact common-bisector localization

Now take \(0<g=\nu(S)<1\).  Apply the ham-sandwich theorem to the two finite
measures

\[
 \nu_1=\mathbf1_S\nu,
 \qquad
 \nu_0=\mathbf1_{S^c}\nu.                            \tag{6.1}
\]

In dimension at least two there is a hyperplane with open sides \(H_+,H_-\)
such that

\[
 \nu_b(H_+)=\nu_b(H_-)=\tfrac12\nu_b(\mathbb R^d),
 \qquad b=0,1.                                      \tag{6.2}
\]

Consequently \(\nu(H_\pm)=1/2\), and both branch posteriors satisfy

\[
 \nu(S\mid H_+)=\nu(S\mid H_-)=g.                   \tag{6.3}
\]

The posterior is log-concave because a halfspace is convex.  If the branch
is sampled according to its probability, then

\[
 g_j\equiv g_0,
 \qquad I(B;Y_j\mid\mathcal F_{j-1})=0.              \tag{6.4}
\]

For a generic bisecting hyperplane, (2.2) gives

\[
 \frac12P_{\nu(\cdot\mid H_+)}(S)
 +\frac12P_{\nu(\cdot\mid H_-)}(S)=P_\nu(S).         \tag{6.5}
\]

Without genericity, the left side is no larger.  Thus this construction
solves both the mass-survival and the perimeter-transfer subproblems exactly.

It does not solve the geometry subproblem.  Given any fixed two-dimensional
subspace \(E\) of possible normals, project both measures in (6.1) to \(E\)
and apply the planar ham-sandwich theorem.  Hence a common-bisecting normal
can be found in every such \(E\).  This allows successive cuts transverse to
all but one affine direction.  Conversely, on a one-dimensional affine
support, a threshold bisects both measures only if their medians coincide,
which need not happen.

This is exactly the two-constraint localization-lemma recursion.  A compact
limiting construction cannot stop in affine dimension at least two, because
another common cut would be available; it can and generally does stop on a
one-dimensional needle.  No quantitative estimate from ham sandwich controls
the length or variance along that last line.

Suppose the resulting disintegration is

\[
 \mu=\int\mu_\omega\,d\pi(\omega),
 \qquad
 \mu_\omega(S)=g,
\]

where \(\mu_\omega\) is a one-dimensional log-concave law on direction
\(\theta_\omega\), with variance \(\sigma_\omega^2\).  The law of total
covariance and isotropy give exactly

\[
 \int \sigma_\omega^2\theta_\omega\theta_\omega^T\,d\pi(\omega)
 \preceq I,
 \qquad
 \int\sigma_\omega^2d\pi(\omega)\le d.              \tag{6.6}
\]

One-dimensional log-concave isoperimetry and (6.5) yield

\[
 P_\mu(S)
 \ge c\min(g,1-g)\int\frac{d\pi(\omega)}{\sigma_\omega}
 \ge \frac{c\min(g,1-g)}{\sqrt d}.                   \tag{6.7}
\]

The last step is Jensen and (6.6).  Improving (6.7) dimension-freely requires
a new estimate such as

\[
 \int\sigma_\omega^2d\pi(\omega)=O(1),               \tag{6.8}
\]

or a comparable weighted inverse-scale bound.  Such an estimate does not
follow from simultaneous bisection.  If every balanced needle had diameter
\(O(1)\), (6.5) and one-dimensional isoperimetry would immediately give the
desired dimension-free Cheeger inequality, so that assertion is the central
claim, not a preprocessing lemma.

The elementary memoryless countertest already rules out branchwise
contraction.  For \(T\sim\operatorname{Exp}(1)\) and \(a\ge0\),

\[
 T\mid\{T\ge a\}\ \stackrel d=\ a+\operatorname{Exp}(1),
 \qquad
 \\operatorname{Var}(T\mid T\ge a)=1.                 \tag{6.9}
\]

In particular the upper median half has mass \(1/2\) but retains the full
variance.  Repeated upper-half bisections only translate the law.  In a
product exponential, take \(S\) independent of the cut coordinate; these are
valid common bisectors of \(S\) and \(S^c\), yet the selected branch makes no
variance progress at all in the cut direction.

## 7. Mandatory model countertests

### 7.1 Isotropic cube

Let the coordinates be uniform on \([ -\sqrt3,\sqrt3]\).  Conditioning one
coordinate to an interval retaining fraction \(p\) gives

\[
 \\operatorname{Var}(X_i\mid A)=p^2,
 \qquad
 \\operatorname{Var}(X_j\mid A)=1\quad(j\ne i).       \tag{7.1}
\]

For a coordinate rectangle with retained fractions \(p_i\),

\[
 \mu(A)=\prod_i p_i,
 \qquad
 \operatorname{Cov}(X\mid A)=\operatorname{diag}(p_i^2),
 \qquad
 \operatorname{diam}(A)^2=12\sum_i p_i^2.            \tag{7.2}
\]

Hence full covariance contraction
\(\operatorname{Cov}(X\mid A)\preceq\rho I\) costs
\(\mu(A)\le\rho^{d/2}\), and diameter \(D\) costs

\[
 \mu(A)\le (D/\sqrt{12d})^d.                         \tag{7.3}
\]

After recentering and covariance normalization, every such rectangle is
again the original isotropic cube.  Raw coordinate shrinkage makes no
affine-normalized shape progress.

### 7.2 Isotropic simplex

For \(U\) uniform on the standard \(d\)-simplex, put \(T=U_1\) and
\(R=1-T\sim\operatorname{Beta}(d,1)\).  In a direction tangent to the opposite
facet, a coordinate likelihood changes variance by

\[
 \frac{\\operatorname{Var}(\langle v,U\rangle\mid A(T))}
      {\\operatorname{Var}\langle v,U\rangle}
 =\frac{\mathbb E[R^2\mid A]}{\mathbb E R^2}.         \tag{7.4}
\]

Among coordinate events of mass \(p\), rearrangement gives the sharp lower
bound \(p^{2/d}\).  Equality holds for the vertex cap

\[
 A_b=\{U_1\ge1-b\},
 \qquad p=b^d,
 \qquad \operatorname{Cov}(U\mid A_b)=b^2\operatorname{Cov}(U).           \tag{7.5}
\]

Thus constant full-spectrum shrinkage again costs \(\exp[-\Theta(d)]\), and
isotropizing the homothetic cap returns the same simplex.  Worse, the
opposite halfspace \(B_b=\{U_1\le1-b\}\) increases every opposite-facet
variance by

\[
 \frac{1-b^{d+2}}{1-b^d}>1.                          \tag{7.6}
\]

At half mass this is \(2-2^{-2/d}>1\).  An ideal central coordinate slice
changes each of the other \(d-1\) variances only by \(1-O(d^{-2})\).

### 7.3 Product one-sided exponentials

For \(X_i=Y_i-1\), \(Y_i\sim\operatorname{Exp}(1)\), the upper-tail identity
(6.9) leaves the entire covariance equal to \(I\); it only translates a
factor.  A soft exponential likelihood does shrink variance, but pays the
product cost exactly:

\[
 L=e^{-\lambda Y_i}:
 \quad Y_i\mid L\sim\operatorname{Exp}(1+\lambda),
 \quad \rho=(1+\lambda)^{-2},
 \quad Z=(1+\lambda)^{-1}=\sqrt\rho.                 \tag{7.7}
\]

Applying it in all coordinates gives \(Z=\rho^{d/2}\), and isotropization
again returns the original product-exponential family.  Equation (3.1) is
the complementary noisy-outcome failure: a log-concave Laplace observation
can make a single posterior variance arbitrarily large.

### 7.4 Radial exponential

Let

\[
 d\mu_d(x)\propto e^{-\sqrt{d+1}\|x\|}dx,
\]

which is isotropic.  An origin halfspace leaves all \(d-1\) transverse
variances exactly equal to one.  The ideal slice \(X_1=0\) leaves transverse
variance \(d/(d+1)\), and \(k\) orthogonal central slices leave

\[
 \frac{d-k+1}{d+1}I                                  \tag{7.8}
\]

on the surviving affine hull.  Thus even ideal observations need
\(\Theta(d)\) directions for constant full-spectrum change.

More generally, if a nonnegative likelihood depends only on a \(k\)-plane
\(E\), then every variance on \(E^\perp\) is at least
\((d-k+1)/(d+1)\).  To see this, write \(x=s+y\), \(s\in E\), and compare the
marginal of \(y\) with the exact slice density.  Their radial likelihood
ratio is

\[
 H(r)=\int_E L(s)
 e^{-\alpha(\sqrt{\|s\|^2+r^2}-r)}ds,
 \qquad \alpha=\sqrt{d+1},                           \tag{7.9}
\]

which is increasing in \(r\).  The marginal therefore radially dominates
the slice.  Cycling fewer than a linear number of independent directions
cannot give constant covariance contraction on this example.

### 7.5 Gaussian coordinate halfspaces and parity

If

\[
 q(x)\propto\gamma_d(x)\prod_j\ell_j(\langle u_j,x\rangle)
\]

and \(E=\operatorname{span}\{u_j\}\), then

\[
 q=q_E\otimes\gamma_{E^\perp}.                       \tag{7.10}
\]

Every unobserved direction keeps variance exactly one.  A coordinate sign
observation has

\[
 \\operatorname{Var}(Z\mid Z\ge0)=1-2/\pi,            \tag{7.11}
\]

but a full sign cycle has branch mass \(2^{-d}\).  For the balanced
coordinate halfspace \(S=\{X_1\ge0}\), all information-free observation
directions lie in \(e_1^\perp\); they can localize \(d-1\) directions but must
leave the \(S\)-normal direction exceptional.  Revealing its sign resolves
\(S\) completely.

Finally take the diagnostic set

\[
 S_{\rm par}=\left\{\prod_{i=1}^d\operatorname{sign}(X_i)=1\right\}.        \tag{7.12}
\]

Its Gaussian mass is \(1/2\), and its first-order correlation vector is zero.
After revealing any \(k<d\) coordinate signs, the posterior mass of
\(S_{\rm par}\) is still exactly \(1/2\); the final sign makes it zero or one.
Thus information and mass loss can occur as a cliff after \(d-1\) perfectly
balanced steps, invisible to covariance-based direction selection.  For
\(d\ge3\), the Gaussian conditioned on parity has covariance \(I\), so all
one- and two-coordinate marginal tests miss the global constraint.  The
parity-conditioned law is not log-concave; its role is to invalidate an
argument based only on scalar marginals or first-order correlations, not to
serve as a log-concave counterexample.

## 8. Final obstruction

The common-bisector construction is the strongest version of the proposed
idea:

\[
 \text{log-concavity preserved}
 \quad+\quad
 g_j\equiv g_0
 \quad+\quad
 \mathbb E P_{\mu_j}(S)\le P_\mu(S).
\]

What remains after those exact identities is a one-dimensional balanced
needle with uncontrolled scale.  Hard slabs can bound that scale only by
learning enough about the last coordinate to risk resolving \(S\); smooth
weak observations converge to Gaussian localization; strongly curved
observations require a linear-in-\(d\) precision budget and have no universal
binary-information bound.  The requested model tests show that neither
cycling nor affine renormalization supplies the missing contraction.

Accordingly, this route gives a rigorous non-Gaussian localization and a
correct perimeter transfer, but no dimension-free curvature or diameter
theorem.  Any lemma that bounded the surviving balanced-needle scale by a
universal constant would, through (6.5), immediately prove the desired
dimension-free isoperimetry and is exactly the unresolved central step.


