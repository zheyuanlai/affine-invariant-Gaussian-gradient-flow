# A balanced high-rank interface: what stability and the long-tube estimate do, and do not, prove

## 0. Verdict

Let \(\mu\) be isotropic and log-concave, let \(E\) be an attained
half-mass Cheeger minimizer, and write

\[
 \mu(E)=\frac12,\qquad p=P_\mu(E)=\frac{\psi}{2},\qquad
 Q=\frac1p\int_{\partial^*E}N\otimes N\,d\sigma_\mu .       \tag{0.1}
\]

Assume

\[
                         1-\operatorname {tr}Q^2\ge\rho,
                         \qquad 0<\rho\le1.                  \tag{0.2}
\]

The proposed combination does **not** close the high-rank branch.  The
precise conclusions are as follows.

1.  Balanced Cheeger minimality gives a global \(\psi\)-quasiminimality
    inequality, the multiplier bound \(|\lambda|\le\psi\), and the usual
    volume-constrained Jacobi stability.  Centered translation tests imply

    \[
     \mathbb E_\sigma\!\left[
       (|\mathrm {II}|^2+\nabla^2V(N,N))|N-m|^2\right]
     +\mathcal C_{\partial\Omega}(N-m)
     \le \mathbb E_\sigma|\mathrm {II}|^2,                  \tag{0.3}
    \]

    where \(d\sigma=p^{-1}d\sigma_\mu\),
    \(m=\mathbb E_\sigma N\), and the contact term is nonnegative.
    Every term in (0.3) is normalized by \(p\).  The affine Minkowski
    identities have the same homogeneity.  They contain no lower bound for
    \(p\).

2.  The impurity assumption forces the boundary normals away from one
    coherent direction, but it gives no upper bound on their effective
    rank:

    \[
      |m|\le(1-\rho)^{1/4},\qquad
      \inf_{|e|=1}\mathbb E_\sigma|N-e|^2
      \ge2\{1-(1-\rho)^{1/4}\}\ge\frac\rho2 .               \tag{0.4}
    \]

3.  The audited killed-normal-tube theorem can be applied with explicit
    parameters.  Put

    \[
                 h=\frac\rho{16},\qquad
                 \gamma=\frac{\rho^2}{1536}.                 \tag{0.5}
    \]

    There is a good part \(\Gamma\subset\partial^*E\) whose discarded
    surface fraction is at most \(\rho/16\), every ray from \(\Gamma\)
    survives a distance at least \(\log(1+\gamma)/\psi\), and

    \[
       I\succeq {A_\rho\over2\psi^2}Q_\Gamma,\qquad
       A_\rho\ge {e^{-1}\gamma^3\over96},                    \tag{0.6}
    \]

    where

    \[
       Q_\Gamma={1\over p}\int_\Gamma N\otimes N\,d\sigma_\mu,
       \qquad
       (\operatorname {tr}Q_\Gamma)^2-\operatorname {tr}Q_\Gamma^2
       \ge{7\rho\over8}.                                    \tag{0.7}
    \]

    Nevertheless (0.6) yields only

    \[
       \psi\ge {\gamma^{3/2}\over\sqrt{192e}}
                   \sqrt{1-\rho/16\over n}.                  \tag{0.8}
    \]

    High rank makes this operator bound weaker.  The factor \(n^{-1/2}\)
    is sharp for the available covariance/tube data.

4.  Equal-slope flat phases give the sharp obstruction.  On each regular
    phase the Jacobi form is only a Dirichlet energy, all centered
    translations are neutral, and finite mass-balanced translations cannot
    lower perimeter.  The balanced product-exponential maximum has
    \(Q=I_m/m\), zero regular Jacobi charge, and vanishing short-tube loss,
    yet its perimeter is bounded below only because of global
    one-dimensional slice normalization.  Its actual perimeter improvement
    comes from beveling codimension-two ridges.  A cyclically constrained
    version has the same local data but is affinely product-irreducible.

Thus the exact missing statement is a **global completion-or-saving
dichotomy**.  Persistent approximately planar packets must either complete
to enough full marginal slices, where one-dimensional isotropic
log-concavity gives \(p\ge c\), or their completion defect must produce a
bounded-reuse, exactly volume-preserving ridge/contact/medial competitor.
Neither (0.3), (0.6), nor (0.2) controls that defect.  A formal sufficient
version is stated in Section 6.  Proving it is the remaining
conjecture-strength step; it may not be replaced by the assertion that the
survivor is a product.

## 1. Exact consequences of balanced minimality

All perimeters in this note are relative to the affine support.  In the
smooth calculation let the support be a convex \(C^2\) domain \(\Omega\),
let \(V\in C^2(\Omega)\) be convex, and let
\(d\mu=Z^{-1}e^{-V}dx\).  The full-support case is obtained by taking
\(\Omega=\mathbb R^n\).  Write

\[
 d\sigma_\mu=Z^{-1}e^{-V}d\mathcal H^{n-1},\quad
 q=|\mathrm {II}|^2+\nabla^2V(N,N).                         \tag{1.1}
\]

### 1.1 A BV inequality which does not use smoothness

**Lemma 1.1 (balanced Cheeger quasiminimality).**  Suppose that the BV
perimeter and the exterior-Minkowski Cheeger constant have been identified
by the usual relative relaxation.  If (0.1) holds, then every finite-
perimeter set \(F\) satisfies

\[
             P_\mu(E)\le P_\mu(F)+\psi\,\mu(E\mathbin\triangle F).
                                                                    \tag{1.2}
\]

If \(E\mathbin\triangle F\Subset U\) and the exterior traces agree, the
same assertion holds with both perimeters restricted to \(U\).  If
\(\mu(F)=1/2\), then the stronger inequality

\[
                              P_\mu(E)\le P_\mu(F)             \tag{1.3}
\]

holds.

**Proof.**  Put \(v=\mu(F)\).  The definition of \(\psi\) gives

\[
 P_\mu(F)\ge\psi\min(v,1-v)
 =\psi\left({1\over2}-\left|v-{1\over2}\right|\right).
\]

Since \(|v-1/2|\le\mu(E\mathbin\triangle F)\) and
\(P_\mu(E)=\psi/2\), this is (1.2).  When \(v=1/2\) there is no error.
Locality of relative BV perimeter gives the localized statement. \(\square\)

In a smooth compact setting the fixed-volume Lagrange multiplier theorem
therefore applies.  The first variation gives

\[
 H_\mu:=\operatorname {tr}\mathrm {II}-\nabla V\cdot N
                  =\lambda.                                  \tag{1.4}
\]

The cusp of \(\min(v,1-v)\) at \(1/2\) is important.

**Lemma 1.2 (multiplier and stability).**  With the preceding smoothness
hypotheses,

\[
                              |\lambda|\le\psi,                \tag{1.5}
\]

and, for every smooth admissible normal speed satisfying
\(\int_\Sigma u\,d\sigma_\mu=0\),

\[
 \begin{split}
 Q_\Sigma(u):={}&\int_\Sigma
   \{ |\nabla_\Sigma u|^2-q u^2\}\,d\sigma_\mu\\
 &-\int_{\partial\Sigma}
      \mathrm {II}_{\partial\Omega}(N,N)u^2\,d\tau_\mu
 \ge0.
 \end{split}                                                 \tag{1.6}
\]

Here the last integral is absent in full support; the convention is
\(\mathrm {II}_{\partial\Omega}\ge0\) for a convex support.  In particular,
it is an additional favorable curvature charge after it is moved to the
right side.

**Proof.**  For a variation whose first volume derivative is \(a\), the
one-sided derivatives of
\(P/\min(v,1-v)\) at \(v=1/2\) are nonnegative only if
\(-\psi\le\lambda\le\psi\).  For \(a=0\), (1.3) and the standard second
variation of \(P-\lambda\mu\) give (1.6). \(\square\)

The claim sometimes made at this point that \(\lambda=\psi\), or even that
\(\lambda=0\), is false.  A nonsymmetric one-dimensional median can have a
nonzero multiplier, while a symmetric central halfspace has multiplier
zero.

### 1.2 Centered translations

Normalize the boundary law and put

\[
 d\sigma=p^{-1}d\sigma_\mu,\qquad m=\mathbb E_\sigma N.       \tag{1.7}
\]

For an ambient orthonormal basis \((e_i)\), the speeds
\(u_i=N_i-m_i\) have zero surface mean.  Since

\[
          \sum_i|\nabla_\Sigma N_i|^2=|\mathrm {II}|^2,
          \qquad \sum_i(N_i-m_i)^2=|N-m|^2,                   \tag{1.8}
\]

summing (1.6) proves the following exact statement.

**Lemma 1.3 (translation trace).**

\[
 \boxed{
 \mathbb E_\sigma[q|N-m|^2]
 +{1\over p}\int_{\partial\Sigma}
   \mathrm {II}_{\partial\Omega}(N,N)|N-m|^2\,d\tau_\mu
 \le\mathbb E_\sigma|\mathrm {II}|^2.}                       \tag{1.9}
\]

For a compactly supported cutoff \(\chi\), before centering one also has
the pointwise trace identity

\[
 \sum_iQ_\Sigma(\chi N_i)
 =\int_\Sigma\{ |\nabla_\Sigma\chi|^2
              -\chi^2\nabla^2V(N,N)\}\,d\sigma_\mu
 -\int_{\partial\Sigma}\chi^2
      \mathrm {II}_{\partial\Omega}(N,N)d\tau_\mu.          \tag{1.10}
\]

All shape-curvature terms cancel in (1.10).  Thus a localized translation
can save perimeter only when normal convexity/contact charge beats the
cutoff capacity.  Boundary normal rank controls neither side.

### 1.3 Affine Minkowski identities

The ambient field \(X(x)=a+Mx\) gives, at a CMC interface,

\[
                    \mathbb E_\sigma(\nabla V+\lambda N)=0,   \tag{1.11}
\]

\[
 \boxed{
 I=\mathbb E_\sigma\left[
       N N^T+(\nabla V+\lambda N)x^T\right].}                 \tag{1.12}
\]

Equation (1.12) is understood after pairing with arbitrary matrices \(M\);
only its symmetric part is needed for symmetric \(M\).  It follows by
equating the geometric first variation

\[
 \int_\Sigma\{\operatorname {tr}M-N^TMN
            -\nabla V\cdot(a+Mx)\}\,d\sigma_\mu
\]

with \(\lambda\int_\Sigma(a+Mx)\cdot N\,d\sigma_\mu\).
The dilation field is the case \(M=I\).

Equations (1.9), (1.11), and (1.12) are homogeneous in the unnormalized
surface measure.  After division by \(p\), only the boundary probability
law and \(\lambda\) remain.  The sole information about the desired scale
is \(|\lambda|\le2p\); the mixed term in (1.12) has no sign under convexity
of \(V\).  Isotropy is a bulk identity and does not turn it into a boundary
estimate.

### 1.4 What survives in BV

Lemma 1.1 is already a relative-BV statement and needs no regularity of
\(E\).  If the density is smooth and positive and the support boundary is
smooth, standard regularity for a perimeter minimizer gives (1.4)--(1.6)
on the regular stratum, for variations compactly supported away from the
singular stratum.  Passing from those compactly supported tests to the
global speeds \(N_i-m_i\) requires all of the following:

* cutoffs around the singular stratum whose surface \(W^{1,2}\) cost tends
  to zero;
* integrability of \(|\mathrm {II}|^2\), \(\nabla^2V(N,N)\), and the
  contact term against those cutoffs; and
* an exhaustion at infinity when the regular interface is noncompact.

Finite perimeter alone supplies none of these global estimates.  For a
nonsmooth convex support, a classical contact second fundamental form is
not even defined on corners.  Smooth approximation of the measure also
does not automatically select minimizers whose normal matrices converge to
the prescribed \(Q\).

Accordingly, (1.9)--(1.12) are used above only in the smooth case or under
these explicit cutoff hypotheses.  This limitation cannot rescue the
route: even granting all three identities globally leads only to Sections
2--3.  Conversely, full BV minimality is stronger than regular Jacobi
stability at a codimension-two corner.  It permits the finite bevel in
(4.9), although that corner has zero perimeter measure and is absent from
the regular second variation.  Any BV completion of the argument must
therefore retain finite ridge/contact/medial competitors, rather than claim
that (1.6) exhausts minimality.

## 2. What matrix impurity adds to stability

The matrix \(Q\) is positive semidefinite and \(\operatorname {tr}Q=1\).
Jensen and the operator bound give

\[
 |m|^2\le\|Q\|_{\mathrm {op}}\le
       \sqrt{\operatorname {tr}Q^2}\le\sqrt{1-\rho}.          \tag{2.1}
\]

If \(m\ne0\), choose \(e=m/|m|\).  Then

\[
 \inf_{|a|=1}\mathbb E_\sigma|N-a|^2
 =2(1-|m|)
 \ge2\{1-(1-\rho)^{1/4}\}\ge{\rho\over2}.                  \tag{2.2}
\]

The same lower bound is trivial if \(m=0\).  Also

\[
 |N-m|\ge1-|m|\ge1-(1-\rho)^{1/4}\ge{\rho\over4}.           \tag{2.3}
\]

Consequently (1.9) implies only the curvature comparison

\[
 {\rho^2\over16}\,\mathbb E_\sigma q
 \le\mathbb E_\sigma|\mathrm {II}|^2,                       \tag{2.4}
\]

with an additional nonnegative contact term omitted.  Since
\(q=|\mathrm {II}|^2+\nabla^2V(N,N)\), (2.4) is compatible
with \(q=0\) and gives no positive absolute scale.  In particular, it does
not put a universal upper bound on

\[
       r_{\mathrm {eff}}(Q)={\operatorname {tr}Q\over
                                  \|Q\|_{\mathrm {op}}}.
\]

Indeed \(Q=I_n/n\) satisfies (0.2) for every fixed \(\rho<1\) and has
effective rank \(n\).

## 3. Combining balance with the audited killed-tube theorem

This section assumes the regular killed-ray coverage theorem, including
stopping at first focal point, support contact, or loss of unique nearest
point.  This is automatic on the fully regular part in a smooth support;
the singular hard-wall extension remains a separate regularity input.

Along a surviving outward normal ray, the weighted Jacobian is

\[
 j_x(t)=\exp\{\lambda t-D_x(t)\},\qquad D_x(0)=0,\quad
 D_x(t)\ge0,\quad D_x'(t)\ge0.                               \tag{3.1}
\]

Let \(0<\gamma<1\), and stop when the outer parallel set has gained mass
\(\gamma/2\).  Since \(p=\psi/2\) and \(|\lambda|\le\psi\), the variable-
multiplier tube theorem gives

\[
 {\log(1+\gamma)\over\psi}\le T
 \le {\gamma\over(1-\gamma)\psi},                           \tag{3.2}
\]

and the normalized curvature/killing-free flux at \(T\) is at least

\[
                  a_\gamma=(1-\gamma)
                    \exp\{-\gamma/(1-\gamma)\}.              \tag{3.3}
\]

After deleting killed rays and surviving rays with \(D_x(T)>h\), the
discarded base fraction is at most

\[
                \delta\le{1-a_\gamma\over1-e^{-h}}.           \tag{3.4}
\]

For \(0<\gamma\le1/4\) and \(0<h\le1\),

\[
 1-a_\gamma\le\gamma+{\gamma\over1-\gamma}\le3\gamma,
 \qquad 1-e^{-h}\ge {h\over2},                               \tag{3.5}
\]

so \(\delta\le6\gamma/h\).  With the choice (0.5),

\[
                              \delta\le{\rho\over16}.         \tag{3.6}
\]

Let

\[
 D={1\over p}\int_{\partial^*E\setminus\Gamma}
             N\otimes N\,d\sigma_\mu,\qquad Q_\Gamma=Q-D.
\]

Then \(D\succeq0\), \(\operatorname {tr}D=\delta\), and

\[
\begin{split}
 &(\operatorname {tr}Q_\Gamma)^2-\operatorname {tr}Q_\Gamma^2\\
 &=(1-\delta)^2-\operatorname {tr}(Q-D)^2\\
 &\ge1-\operatorname {tr}Q^2-2\delta
 \ge {7\rho\over8}.
\end{split}                                                  \tag{3.7}
\]

The covariance part of the tube theorem is

\[
 I=\operatorname {Cov}(\mu)
 \succeq {A(\gamma,h)\over\psi^3}
                 \int_\Gamma N\otimes N\,d\sigma_\mu,
                                                                    \tag{3.8}
\]

where

\[
 A(\gamma,h)={\{\log(1+\gamma)\}^3\over12}
       \exp\{-h-\gamma/(1-\gamma)\}.                         \tag{3.9}
\]

Since \(p=\psi/2\), this is exactly (0.6).  Our parameters satisfy
\(h+\gamma/(1-\gamma)<1\), and
\(\log(1+\gamma)\ge\gamma/2\), whence

\[
                              A(\gamma,h)\ge{e^{-1}\gamma^3\over96}.
                                                                    \tag{3.10}
\]

Taking the largest eigenvalue in (3.8) proves

\[
       \psi^2\ge {A(\gamma,h)\over2}
                         \|Q_\Gamma\|_{\mathrm {op}}.        \tag{3.11}
\]

But \(\operatorname {tr}Q_\Gamma=1-\delta\) supplies only

\[
                    \|Q_\Gamma\|_{\mathrm {op}}
                       \ge{1-\delta\over n}.                 \tag{3.12}
\]

Equations (3.10)--(3.12) give (0.8).  The impurity lower bound (3.7)
cannot improve (3.12): matrices \((1-\delta)I_n/n\) have impurity tending
to \((1-\delta)^2\) and operator norm tending to zero.

There is also an exact algebraic saturation.  Put base mass \(p/n\) on
each of the rays \(\mathbb Re_i\), give every ray uniform density on
\([-T,T]\), and take \(T=h/p\).  The covariance contributed in direction
\(e_i\) is \(2h^3/(3np^2)\).  Choosing

\[
                         p=\sqrt{2/3}\,h^{3/2}n^{-1/2}        \tag{3.13}
\]

saturates unit covariance in all directions.  This current is not
log-concave and is not a KLS counterexample.  It proves that isotropy,
two-sided swept mass, tube length, and the normal matrix contain no stronger
moment inequality.  The missing input must use global log-concavity and
overlap, not another manipulation of covariance.

## 4. Sharp flat-phase obstruction to Jacobi and translation arguments

### 4.1 Exact finite identity

Consider a flat chart \((z,s)\) in which

\[
 E=\{s<0\},\qquad e^{-V(z,s)}=w(z)e^{-\kappa s}.              \tag{4.1}
\]

Replace the graph \(s=0\) by \(s=f(z)\), where \(f\) is compactly
supported.  The exact changes of mass and perimeter satisfy

\[
 \boxed{
 \Delta P+\kappa\Delta\mu
 =\int w(z)e^{-\kappa f(z)}
       \{\sqrt{1+|\nabla f(z)|^2}-1\}\,dz\ge0.}              \tag{4.2}
\]

Consequently, for any finite collection of disjoint equal-slope charts,

\[
                    \sum_j\Delta\mu_j=0
           \quad\Longrightarrow\quad
                    \sum_j\Delta P_j\ge0.                    \tag{4.3}
\]

The quadratic limit is

\[
                Q(f)=\int |\nabla f|^2\,d\sigma_\mu.         \tag{4.4}
\]

Thus regular second variation, centered translations, and their finite
graph versions are all exactly neutral on componentwise constants.  They
do not see the incidence between different flat phases.

### 4.2 Balanced independent exponentials

Let \(Y_1,\ldots,Y_m\) be independent unit exponentials and center by
\(X_i=Y_i-1\); this law is isotropic.  For

\[
 E_L=\{\max_iY_i\ge L\},\qquad
             (1-e^{-L})^m={1\over2},                         \tag{4.5}
\]

one has

\[
 p_m=P(E_L)={m\over2}(2^{1/m}-1)
       \in\left[{\log2\over2},{1\over2}\right],             \tag{4.6}
\]

\[
                         Q={I_m\over m},\qquad
                         1-\operatorname {tr}Q^2=1-{1\over m}.
                                                                    \tag{4.7}
\]

Every regular facet is flat, \(\nabla^2V=0\), the common weighted CMC
slope is one, and (4.4) is the complete regular Jacobi form.  If
\(T=h/p_m\), \(0<h\le1/10\), the two normalized tube losses are at most
\(3.5h\) and \(2.9h\).  Nevertheless the full coordinate slice area is
\(e^{-L}\), whereas the observed facet is truncated by all other
thresholds.  The total completion defect satisfies

\[
 {R_{\rm comp}\over p_m}=2(1-e^{-L})-1\longrightarrow1.      \tag{4.8}
\]

Thus killed loss tending to zero does not control global slice completion.

This set is a stress test, not an isoperimetric minimizer.  Its failure of
minimality occurs at the surface-measure-zero tie ridges.  Replacing the
largest-two-coordinate corner by the fixed bevel

\[
 C_m(q,.1)=\{Y_i\le q,\ Y_i+Y_j\le2q-.1\quad(i<j)\}
\]

and correcting \(q\) to retain mass \(1/2\) gives

\[
 \limsup_m P(C_m(q_m,.1))< {\log2\over2}-.0061.              \tag{4.9}
\]

This strict finite saving is invisible in (4.4).  It is a codimension-two
ridge-capacity effect.

### 4.3 A product-irreducible flat survivor

Let \(a=m^{-6}\) and condition the product exponential law on

\[
 \Omega_{m,a}=\{y_i\ge0,\ y_i+y_{i+1}\ge a\quad(i\bmod m)\}. \tag{4.10}
\]

At the median maximum level, the perimeter and normal matrix satisfy

\[
                  .33<p_m<.44,\qquad Q={I_m\over m}.          \tag{4.11}
\]

All regular shape, Hessian, and free-boundary contact-curvature charges
vanish.  For \(T=h/p_m\), the inward killed loss is at most
\(1.4T+4m^{-11}<4.3h+4m^{-11}\), and the outward loss is at most
\(T<3.04h\).  After centering and whitening, covariance lies between
\(.99I\) and \(1.01I\), so the same assertions hold with a universally
elliptic constant anisotropy.

The support is not affinely a nontrivial product.  Its irredundant facet
normals are

\[
                       e_i,\qquad e_i+e_{i+1}.                \tag{4.12}
\]

If they split between two dual direct summands, adjacency forces every
\(e_i\) into the same summand, which spans the whole dual space.  Hence a
high-rank, zero-Jacobi, low-killing survivor cannot simply be declared a
product.  Global support/Hessian coupling is a separate datum.

## 5. Gaussian, simplex, and \(\ell_1\)-ball checks

### 5.1 Gaussian halfspace

For standard Gaussian measure and \(E=\{x_1\le0\}\),

\[
 p=(2\pi)^{-1/2},\quad \psi=\sqrt{2/\pi},\quad
 \lambda=0,\quad Q=e_1e_1^T.                                \tag{5.1}
\]

The impurity is zero, so this is the coherent branch.  All centered
translation speeds vanish.  The normal tube bound sees its one eigenvalue
and gives a dimension-free estimate.  No high-rank argument should reject
this model.

### 5.2 Uniform regular simplex

For the standard simplex, one barycentric coordinate \(X_0\) has

\[
 \mathbb P(X_0\ge s)=(1-s)^n,\qquad
 \operatorname {Var}(X_0)={n\over(n+1)^2(n+2)}.              \tag{5.2}
\]

After whitening, its half-mass parallel-facet cut is rank one and has
Cheeger ratio

\[
 {P\over1/2}={n\sqrt n\over(n+1)\sqrt{n+2}}\,2^{1/n}
 \longrightarrow1.                                          \tag{5.3}
\]

Here \(V=0\) on the interior cut; the global scale is carried entirely by
the contact with the support facets.

For comparison, the centered homothetic inner simplex of mass \(1/2\) is
high rank.  In isotropic position the inradius is
\(r=\sqrt{(n+2)/n}\), its equal facet normals form a tight frame, and

\[
 Q={I_n\over n},\qquad
 p={n\over r}\,2^{-(n-1)/n}\asymp {n\over2}.                 \tag{5.4}
\]

It is already far above the desired scale and is unstable to ridge/global
competitors.  Local flat-facet second variation alone does not distinguish
it from the product-exponential phase pattern.

### 5.3 Uniform isotropic cross-polytope

For the unscaled \(\ell_1^n\) unit ball,

\[
 f_{X_1}(t)={n\over2}(1-|t|)^{n-1},\qquad
 \operatorname {Var}(X_1)={2\over(n+1)(n+2)}.                \tag{5.5}
\]

The isotropic central coordinate cut is rank one and has ratio

\[
 {P\over1/2}=n\sqrt{2\over(n+1)(n+2)}\longrightarrow\sqrt2. \tag{5.6}
\]

The centered homothetic inner cross-polytope has \(Q=I_n/n\).  Its
isotropic inradius is

\[
 r=\sqrt{(n+1)(n+2)\over2n},
\]

and its half-mass perimeter is

\[
 p={n\over r}\,2^{-(n-1)/n}\asymp\sqrt{n/2}.                 \tag{5.7}
\]

Again the high-rank candidate is harmless because of a global scale which
is absent from the normalized Jacobi identities; the competitive balanced
cut is the rank-one marginal cut.

## 6. The exact missing inequality

The preceding calculations isolate a concrete sufficient theorem.  It is
useful to state it in the planar form because every constant in its
completion branch is already known.

Let persistent good boundary packets \(\Gamma_i\) lie in hyperplanes with
unit normals \(u_i\).  Write \(a_i\) for their weighted areas, \(s_i\) for
the weighted areas of the corresponding **complete** hyperplane slices,
and \(q_i\le1/2\) for the masses of the smaller oriented halfspaces.  Put

\[
 p_G=\sum_i a_i,\qquad
 R_{\rm comp}=\sum_i(s_i-a_i).                               \tag{6.1}
\]

For an isotropic log-concave marginal, the one-dimensional estimate

\[
                              s_i\ge {q_i\over2\sqrt3}         \tag{6.2}
\]

holds.  Hence the following branch is already proved:

**Lemma 6.1 (full-slice closure).**  If

\[
        \sum_iq_i\ge\alpha,\qquad
        R_{\rm comp}\le\eta p_G,                              \tag{6.3}
\]

then

\[
                              p_G\ge{\alpha\over2\sqrt3(1+\eta)}.
                                                                    \tag{6.4}
\]

The swept tubes give a fixed \(\alpha\) once \(\rho\) is fixed.  What is
not proved is the complementary branch:

> **Global completion-or-saving theorem (missing).**  There are constants
> \(c_0,\eta_0>0\), depending only on the fixed impurity threshold \(\rho\),
> and a dimension-free packetization of the good rays such that either
> \(R_{\rm comp}\le\eta_0p_G\), or there is a finite-perimeter set \(F\)
> with
> \[
>       \mu(F)=\frac12,\qquad
>       P_\mu(F)\le P_\mu(E)-c_0R_{\rm comp}.                 \tag{6.5}
> \]
> Curved packets must be charged by the same statement through focal,
> contact, Hessian, or medial capacity, with bounded reuse.

For an actual minimizer the second alternative is impossible.  Equations
(6.3)--(6.4) would then give a dimension-free lower bound.  The
product-exponential maximum realizes the large-completion branch and its
simultaneous bevel verifies the saving alternative.  Equation (4.8) proves
that the premise of the saving alternative cannot be replaced by killed
flux.  The cyclic example proves that it cannot be replaced by an informal
product classification.

Equivalently, a coordinate-free high-rank closure strong enough for the
present task would be

\[
 (\operatorname {tr}Q_\Gamma)^2-\operatorname {tr}Q_\Gamma^2
       \le C\psi^2+C\delta_{\rm tube},                        \tag{6.6}
\]

where \(\delta_{\rm tube}\) includes, with a proved bounded-reuse rule,
all focal, contact, ridge, completion, and mixed-support/Hessian charge.
With (3.7), first take the tube parameters so that
\(C\delta_{\rm tube}\le\rho/2\); (6.6) gives
\(\psi\ge\sqrt{\rho/(2C)}\).  Current Jacobi and killed-ray estimates
control only the regular curvature and the flux actually visited by the
tube.  They do not control the global completion and coupling terms required
in \(\delta_{\rm tube}\).

This is the sharp stopping point of the balanced-interface
second-variation route.  No Poincare inequality, dimension-free boundary
capacity estimate, or KLS-equivalent assertion has been inserted.
