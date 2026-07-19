# Weighted boundary stability and the Cheeger-height obstruction

This note fixes the sign conventions and records the consequences of testing a
weighted isoperimetric hypersurface by translations, rotations, and linear
ambient fields.  Its main conclusion is negative but precise: at mass one half
all these identities are homogeneous in the boundary measure and do not contain
the Cheeger height.  Away from mass one half, a positive-Hessian regularization
precludes an attained smooth Cheeger minimizer altogether.

## 1. Variation formulas

Let
\[
 d\mu=\rho\,dx,\qquad \rho=Z^{-1}e^{-V},
\]
where \(V\in C^3(\mathbb R^n)\) is convex and \(\rho>0\).  Let \(A\) have a
compact \(C^3\) boundary \(\Sigma\), with outward unit normal \(N\).  Put
\[
 d\sigma_\mu=\rho\,d\mathcal H^{n-1},\quad
 P=\sigma_\mu(\Sigma),\quad
 \mathrm{II}(\tau)=D_\tau N,\quad H=\operatorname{tr}\mathrm{II},
\]
and define
\[
 H_\mu=H-\langle\nabla V,N\rangle,\qquad
 q=|\mathrm{II}|^2+\nabla^2V(N,N).
\]
For a compactly supported \(C^2\) flow with initial normal speed \(u\),
\[
 \delta\mu(A)=\int_\Sigma u\,d\sigma_\mu,
 \qquad
 \delta P=\int_\Sigma H_\mu u\,d\sigma_\mu.                 \tag{1}
\]
If \(H_\mu\equiv\lambda\), then the second variation of
\(P-\lambda\mu(A)\) is
\[
 Q(u)=\int_\Sigma\left(|\nabla_\Sigma u|^2-q u^2\right)d\sigma_\mu
     =-\int_\Sigma uJ_\Sigma u\,d\sigma_\mu,                \tag{2}
\]
where
\[
 J_\Sigma=\Delta_\Sigma-\langle\nabla_\Sigma V,\nabla_\Sigma\rangle+q.
\]
Indeed the material derivative of weighted mean curvature under a normal
variation is \((H_\mu)'=-J_\Sigma u\).  A fixed-volume local minimizer obeys
\[
 Q(u)\ge0\quad\hbox{whenever}\quad\int_\Sigma u\,d\sigma_\mu=0. \tag{3}
\]
The same formulas hold for a complete noncompact boundary after cutoff if all
terms in (1)--(2), including the cutoff-error terms, are integrable.  Without
these hypotheses the use of the constant, translation, or affine speeds below
is not justified.

For an actual weighted isoperimetric region, the regular part of the reduced
boundary is smooth, while in dimensions at least eight a singular set may be
present.  Formula (2) is initially available only for test functions compactly
supported in the regular part.  Passing to \(u=N_i\), rotations, or affine
speeds requires both (i) zero \(W^{1,2}\)-capacity cutoffs around the singular
set and (ii) global integrability of the displayed curvature and moment terms.
Neither follows merely from finite perimeter.  Accordingly every consequence
below is stated for a smooth compact boundary (or under the explicit cutoff
hypothesis); using it in a general KLS proof would require a separate regularity
and exhaustion lemma with dimension-free errors.

## 2. What Cheeger minimality gives

Suppose first that \(t=\mu(A)<1/2\) and that \(A\) is a smooth local minimizer
of \(P/\mu(A)\).  If \(h=P/t\), differentiating the quotient in arbitrary
directions gives
\[
 H_\mu\equiv h,\qquad Q(u)\ge0\quad\hbox{for every }u.       \tag{4}
\]
There is no mean-zero restriction in (4).

At \(t=1/2\), by contrast, the denominator \(\min(t,1-t)\) has a cusp.  If
\(A\) is stationary for fixed volume, write \(H_\mu\equiv\lambda\).  The two
one-sided first variations of the Cheeger quotient give only
\[
 |\lambda|\le 2P=h,                                        \tag{5}
\]
and second variation gives only (3).  In particular \(\lambda=h\) is false in
general; for every centrally symmetric example one has \(\lambda=0\) at the
central halfspace while \(h=2P>0\).

For an arbitrary set whose quotient is merely \((1+\varepsilon)\)-optimal,
none of (3)--(4) follows.  A separate quantitative variational-selection
argument is required.  This point cannot be hidden in a smooth approximation.

## 3. Translation obstruction and the exact half-mass dichotomy

When \(H_\mu\) is constant, translating \(\Sigma\) by a vector \(a\) gives
\[
 J_\Sigma\langle a,N\rangle=\nabla^2V(a,N).                 \tag{6}
\]
Consequently
\[
 \sum_{i=1}^n Q(N_i)=-\int_\Sigma\nabla^2V(N,N)\,d\sigma_\mu. \tag{7}
\]
Combining (4) and (7) proves:

**No attained sub-half minimizer under positive normal curvature.**  If
\(\nabla^2V(N,N)>0\) on a subset of positive \(\sigma_\mu\)-measure, no compact
smooth set of mass below \(1/2\) can locally minimize the Cheeger quotient.
In particular, after adding \(\varepsilon|x|^2/2\) to the potential, every
smooth attained Cheeger minimizer must occur at the half-mass cusp.  A limiting
sub-half witness can therefore disappear to infinity as \(\varepsilon\downarrow0\);
the exponential examples below do exactly this.

For a half-mass fixed-volume minimizer, normalize the surface law by
\(d\sigma=P^{-1}d\sigma_\mu\), and put
\[
 m=\mathbb E_\sigma N.
\]
The admissible translation speeds are \(u_i=N_i-m_i\).  Since
\(\sum_i|\nabla_\Sigma N_i|^2=|\mathrm{II}|^2\), summing (3) gives the exact
inequality
\[
 \boxed{\quad
 \mathbb E_\sigma\big[q|N-m|^2\big]\le
 \mathbb E_\sigma|\mathrm{II}|^2.\quad}                    \tag{8}
\]
Let \(r=|m|\), and, when \(r>0\), let \(\nu=m/r\).  Then
\[
 \mathbb E_\sigma|N-\nu|^2=2(1-r),\qquad
 \mathbb E_\sigma|N-m|^2=1-r^2.                            \tag{9}
\]
Moreover
\[
 \mathbb E_\sigma[q|N-m|^2]
 \ge (1-r)^2\mathbb E_\sigma q.                            \tag{10}
\]
Thus for every \(0<\eta<1\) there is the rigorous dichotomy
\[
 \mathbb E|\mathrm{II}|^2\ge\eta\mathbb E q,
 \quad\hbox{or}\quad
 \mathbb E|N-\nu|^2\le2\sqrt\eta.                          \tag{11}
\]
The first branch merely says that extrinsic curvature supplies a fixed fraction
of the Jacobi potential.  In the second branch the normal is aligned in boundary
\(L^2\), but converting this into closeness of \(A\) to a halfspace requires a
global weighted BV/Poincare estimate in the transverse directions.  No such
estimate follows from (8).

The limiting case does close: if \(N=\nu\) almost everywhere on a connected
boundary, the boundary lies in a hyperplane orthogonal to \(\nu\); a half-mass
side has perimeter equal to the density at a median of the one-dimensional
marginal \(\langle X,\nu\rangle\).  That marginal is isotropic and log-concave,
so its median density is bounded below explicitly as follows.  If \(f\) is a
one-dimensional log-concave density, \(M=\|f\|_\infty\), and \(m\) is a median,
then
\[
 f(m)\ge M/2.                                                \tag{11c}
\]
To prove this, suppose for definiteness that a mode \(x_0\) lies to the left of
\(m\), put \(r=f(m)/M\), and integrate the chord lower bound for \(\log f\) on
\([x_0,m]\) and the continuation upper bound to the right of \(m\).  They give
\[
 \int_{x_0}^m f\ge \frac{M(m-x_0)(1-r)}{-\log r},\qquad
 \int_m^\infty f\le \frac{Mr(m-x_0)}{-\log r}.
\]
The first integral is at most the left half-mass and the second is the right
half-mass, so \(1-r\le r\).  The other ordering is identical, and a nonattained
mode follows by approximation.  Also every density bounded by \(M\) obeys
\[
 \operatorname{Var}(X)\ge \frac1{12M^2};                    \tag{11d}
\]
indeed \(\mathbb P(|X-a|\le s)\le2Ms\), and integration of
\(2s(1-2Ms)\) over \(0\le s\le(2M)^{-1}\), with \(a=\mathbb EX\), gives
(11d).  Combining (11c)--(11d), an isotropic marginal satisfies
\[
 f(m)\ge\frac1{\sqrt{48}}=\frac1{4\sqrt3}.                  \tag{11e}
\]
Thus the exactly flat half-mass branch has \(P\ge1/(4\sqrt3)\) and
\(h=2P\ge1/(2\sqrt3)\).
For completeness, in one dimension the Cheeger constant itself is
\(2f(m)\): the distribution functions \(F\) and \(1-F\) are log-concave, so
\(f/F\) decreases up to the median and \(f/(1-F)\) increases after it; the
perimeter of a finite union of intervals is the sum of its endpoint densities.
Approximation extends the calculation to every Borel set of finite perimeter.

There is a concrete reason that an approximate version cannot be inserted for
free.  Let \(\gamma\) be standard Gaussian measure on \(\mathbb R\), let \(\nu\)
be any probability on \(\mathbb R^d\), and let \(B\) have
\(\nu(B)=1/2\) and perimeter \(p\).  In \(\gamma\otimes\nu\), define
\[
 A_L=\big((-\infty,L]\times B\big)
       \cup\big((-\infty,-L]\times B^c\big).
\]
Then \((\gamma\otimes\nu)(A_L)=1/2\), and, ignoring the codimension-two
intersections (or by direct BV approximation),
\[
 P(A_L)=\varphi(L)+p\,[\Phi(L)-\Phi(-L)].                    \tag{11a}
\]
The horizontal boundary has normal \(e_1\); the vertical boundary has a normal
orthogonal to \(e_1\).  Consequently
\[
 \frac1{P(A_L)}\int_{\partial^*A_L}|N-e_1|^2d\sigma
 =\frac{2p[\Phi(L)-\Phi(-L)]}{P(A_L)}.                       \tag{11b}
\]
If \(p\ll1\) and \(L\) is selected so that \(\varphi(L)=\sqrt p\), the
perimeter and normal error are both of order \(\sqrt p\).  Thus a purported
dimension-free implication from boundary-normal alignment to positive
perimeter would already have to exclude a lower-dimensional bottleneck.  In
analytic language, that exclusion is precisely the missing transverse
BV/Cheeger estimate.  Formula (11a) is a dimension-descent obstruction, not a
proof of such a bottleneck's existence in the isotropic log-concave class.

## 4. Affine first variations

For the ambient field \(X(x)=a+Mx\),
\[
 \delta P=\int_\Sigma
 \left(\operatorname{tr}M-N^TMN-\langle\nabla V,a+Mx\rangle\right)d\sigma_\mu,
\]
whereas \(\delta\mu(A)=\int_\Sigma\langle a+Mx,N\rangle d\sigma_\mu\).
At a weighted-CMC hypersurface these identities, divided by \(P\), become
\[
 \mathbb E_\sigma(\nabla V+\lambda N)=0,                    \tag{12}
\]
\[
 \boxed{\quad I=\mathbb E_\sigma\left[NN^T+
             (\nabla V+\lambda N)x^T\right].\quad}          \tag{13}
\]
They involve the normalized boundary law and \(\lambda\), but not \(P\).  At
half mass they consequently contain no Cheeger height \(h=2P\).

## 5. Rotation and full-linear stability sums

For \(i<j\), use the rotation speed
\[
 u_{ij}=x_jN_i-x_iN_j,\qquad \bar u_{ij}=\mathbb E_\sigma u_{ij}.
\]
Let \(S=\mathrm{II}\) be the shape operator and \(x_T=x-(x\cdot N)N\).
A direct differentiation in an orthonormal tangent frame gives
\[
 \sum_{i<j}u_{ij}^2=|x_T|^2,                                \tag{14}
\]
\[
 \sum_{i<j}|\nabla_\Sigma u_{ij}|^2
 =(n-1)+|x|^2|\mathrm{II}|^2-|Sx_T|^2-2(x\cdot N)H.         \tag{15}
\]
Therefore constrained stability implies
\[
 \mathbb E_\sigma\!\left[(n-1)+|x|^2|\mathrm{II}|^2
                  -|Sx_T|^2-2(x\cdot N)H\right]
 \ge \mathbb E_\sigma\!\left[q\sum_{i<j}(u_{ij}-\bar u_{ij})^2\right]. \tag{16}
\]

For the full linear basis put \(v_{ij}=N_i x_j\) and
\(B=\mathbb E_\sigma[Nx^T]\).  Pointwise,
\[
 \sum_{i,j}|\nabla_\Sigma v_{ij}|^2=|x|^2|\mathrm{II}|^2+n-1. \tag{17}
\]
Thus
\[
 \boxed{\quad
 \mathbb E_\sigma[|x|^2|\mathrm{II}|^2+n-1]
 \ge \mathbb E_\sigma[q|Nx^T-B|_F^2].\quad}                \tag{18}
\]
The explicit \(n-1\) in (16) and (18) is sharp: for a Gaussian central
halfspace, \(\mathrm{II}=0\), \(q=1\), and the right side equals \(n-1\).
Thus summing affine modes does not yield a dimension-free surplus.

## 6. Why the Cheeger height is invisible

Every identity (8), (12), (13), (16), and (18) is homogeneous of degree one in
the unnormalized surface measure.  After division by \(P\), it is independent
of \(P\).  At half mass \(h=2P\), whereas the CMC multiplier is the unrelated
quantity \(\lambda\) subject only to (5).  Hence these identities cannot by
themselves imply a lower bound for \(h\).  A global trace/normalization estimate
linking the normalized boundary law to the bulk probability normalization is
load-bearing.  For an exact hyperplane this missing estimate is precisely the
one-dimensional lower bound on the density of an isotropic log-concave marginal
at its median; for a curved surface it is already an isoperimetric statement.
In dimension one this is completely transparent: the boundary of a half-line
is a single point, so every mean-zero boundary variation is zero and constrained
stability is vacuous, while the Cheeger height is \(2f(m)\).  The lower bound
(11e) is global information about normalization and variance, not a consequence
of boundary stability.

## 7. Model checks

### Gaussian

For \(N(0,I_n)\) and \(A=\{x_1\le0\}\), one has
\[
 t=1/2,\quad P=(2\pi)^{-1/2},\quad h=\sqrt{2/\pi},\quad
 \lambda=0,\quad N=e_1,\quad \mathrm{II}=0,\quad q=1.
\]
All centered translation speeds vanish.  On the boundary, the nonzero affine
speeds are the tangent coordinates, and (18) is equality \(n-1=n-1\).  The
allowed modes are therefore completely neutral and do not recover \(h\).

### Cube

The isotropic cube is \([ -\sqrt3,\sqrt3]^n\).  Its central coordinate cut has
\[
 t=1/2,\qquad P=\frac1{2\sqrt3},\qquad \frac{P}{t}=\frac1{\sqrt3}.
\]
The uniform density and the contact edges violate the full-support smooth
hypotheses above.  A clean approximation is the isotropic product density with
one-dimensional potential proportional to \(|x|^p\), \(p>2\), followed by
\(p\to\infty\).  At the central coordinate plane, \(V_{11}(0)=0\), so
\(q=0\), and every stability inequality reduces to a nonnegative tangential
Dirichlet energy.  The perimeter value is supplied only by one-dimensional
normalization, not by stability.

### Simplex

Let \(X\) be uniform on the standard simplex and let \(X_0\) be one barycentric
coordinate.  Then
\[
 \mathbb P(X_0\ge s)=(1-s)^n,\qquad
 \operatorname{Var}(X_0)=\frac{n}{(n+1)^2(n+2)}=:\sigma_n^2.
\]
After whitening, the parallel-facet cut is a halfspace with unit normal, and
its boundary ratio is
\[
 \frac{P_s}{\mu(A_s)}=\frac{n\sigma_n}{1-s}.
\]
At half mass, \(1-s=2^{-1/n}\), hence
\[
 \frac{P_s}{1/2}=\frac{n\sqrt n}{(n+1)\sqrt{n+2}}\,2^{1/n},
\]
which tends to one.  Locally the cut is flat and \(V=0\); all information is in
the contact with the facets.  Any smooth approximation moves that information
into a thin region of large Hessian.  Dropping either the contact terms or that
Hessian region incorrectly makes the stability identities vacuous.

### The \(\ell_1^n\) ball

For the uniform law on the unit cross-polytope, the first coordinate has density
and variance
\[
 f_{X_1}(t)=\frac n2(1-|t|)^{n-1}1_{[-1,1]}(t),\qquad
 \operatorname{Var}(X_1)=\frac{2}{(n+1)(n+2)}.
\]
After isotropic scaling, the central coordinate cut therefore has Cheeger ratio
\[
 2f_{X_1/\sigma}(0)
 =n\sqrt{\frac{2}{(n+1)(n+2)}},
\]
which tends to \(\sqrt2\).  As for the simplex and cube, the flat interior
second-variation calculation has \(q=0\); the missing information is carried by
the contact with the nonsmooth support boundary.

### Product exponentials and a nonsymmetric example

For the variance-one symmetric exponential density
\(\rho(x)=2^{-1/2}e^{-\sqrt2|x|}\), the coordinate tail
\(A_a=\{x_1\ge a\}\), \(a\ge0\), satisfies
\[
 t=\tfrac12e^{-\sqrt2a},\qquad P=\tfrac1{\sqrt2}e^{-\sqrt2a},
 \qquad P/t=\sqrt2.
\]
Thus in one dimension exact Cheeger witnesses (and in every product, coordinate
witnesses with the same ratio) can have \(P\downarrow0\).  Along their boundary
\(\mathrm{II}=0\), \(q=0\), and the unconstrained stability form is just the
tangential Dirichlet energy.

For a strongly nonsymmetric variance-one example, take the shifted exponential
\(\rho(x)=e^{-(x+1)}1_{\{x\ge-1\}}\), which has mean zero and variance one.
For \(A_a=\{x_1\ge a\}\),
\[
 t=e^{-(a+1)},\qquad P=e^{-(a+1)},\qquad P/t=1.
\]
Again \(P\to0\) while the ratio stays fixed.  Smooth strictly convex
approximations replace these attained tail witnesses by minimizing sequences
escaping to infinity, exactly as predicted by (7).
