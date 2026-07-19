# Balanced-half displacement: entropy, covariance, and determinant audit

## 0. Verdict

Let \(\mu=e^{-V}dx\) be an isotropic log-concave probability on its
affine support, let \(\mu(E)=1/2\), and put

\[
 \mu_+=2\mathbf 1_E\mu,\qquad \mu_-=2\mathbf 1_{E^c}\mu.
\]

Let \(T=\nabla\phi\) be the Brenier map from \(\mu_+\) to \(\mu_-\),
write

\[
 X\sim\mu_+,\qquad Y=T(X),\qquad D=Y-X,
 \qquad K=\mathbb E[D D^{T}],
\]

and set \(F_t=(1-t)I+tT\), \(\nu_t=(F_t)_\#\mu_+\).  The conclusions
of this audit are as follows.

1. The whole interpolation satisfies \(\nu_t\le 2\mu\).  More
   precisely,

   \[
   \operatorname {Ent}_\mu(\nu_t)
   =\log 2-\mathbb E\big[\Delta_{V,t}(X,Y)+\Delta_{A,t}(D_aT(X))\big],
   \tag{0.1}
   \]

   where both deficits are nonnegative.  The second distributional
   derivative in \(t\) is exactly potential curvature along the transport
   chords plus the squared Eulerian velocity gradient.  The resulting
   one-bit budget is recorded with exact kernels in Section 2.

2. The often quoted formula

   \[
   \operatorname {Cov}(\nu_t)=I-t(1-t)K                 \tag{0.2}
   \]

   is **not true for general \(t\)**.  It is true at \(t=1/2\), and a
   symmetrized version is true for raw second moments.  The correct
   all-time formulas are

   \[
   \begin{aligned}
   \mathbb E[Z_tZ_t^T]
     &=(1-t)Q_++tQ_- -t(1-t)K,\\
   \operatorname {Cov}(Z_t)
     &=(1-t)C_++tC_- -t(1-t)K_c,
   \end{aligned}                                        \tag{0.3}
   \]

   where \(Z_t=X+tD\), \(Q_\pm\) are endpoint raw second moments,
   \(C_\pm\) are endpoint covariances, and
   \(K_c=\operatorname {Cov}(D)\).  At the midpoint,

   \[
   \boxed{\operatorname {Cov}(\nu_{1/2})=I-\tfrac14K.} \tag{0.4}
   \]

3. There is a universal spectral floor

   \[
   \frac1{4800}I\preceq \operatorname {Cov}(\nu_{1/2})\preceq I.
   \tag{0.5}
   \]

   Consequently a dimension-free midpoint determinant bound is neither a
   free entropy consequence nor a weaker substitute for the desired trace
   estimate: it is quantitatively equivalent to it.  If
   \(C=\operatorname {Cov}(\nu_{1/2})\), then

   \[
   \frac14\operatorname {tr}K
   \le -\log\det C
   \le 1200\operatorname {tr}K.                         \tag{0.6}
   \]

4. The natural covariance Brunn--Minkowski conjecture

   \[
   \log\det\operatorname {Cov}(\nu_{1/2})
   \ge \frac12\log\det C_+ +\frac12\log\det C_-       \tag{0.7}
   \]

   is false.  Section 5 gives a completely explicit one-dimensional
   log-concave counterexample in which \(C_+=C_-=1\), while

   \[
   \operatorname {Var}(\nu_{1/2})
     =\frac{9\sqrt5-5}{16}<1.                           \tag{0.8}
   \]

   In this example the potential deficit and matrix AM--GM deficit vanish
   for **every** \(t\), and
   \(\operatorname {Ent}_\mu(\nu_t)=\log2\) for every \(t\).  Thus even
   the entire entropy interpolation does not charge the covariance loss.

5. Radial Gaussian and radial exponential median balls, cube corridors,
   a simplex cap, the product-exponential maximum cut, and piecewise
   translations all have dimension-free transport cost by direct
   calculations which use no KLS input.  They expose complementary
   limitations of curvature, Fisher information, determinant, and
   Brunn--Minkowski arguments.

No dimension-free proof of \(\operatorname {tr}K\le C\) is obtained.
The exact remaining gap is the absolute midpoint determinant estimate

\[
 \boxed{-\log\det\operatorname {Cov}(\nu_{1/2})\le C,}  \tag{0.9}
\]

or, equivalently by (0.6), control of the cellwise translation constants
and singular part of the Brenier derivative which are invisible to
(0.1).  This is precisely the balanced-half transport form of the KLS
target, not an auxiliary lemma that can be assumed.

All arguments below are intrinsic to the affine support.  We therefore
write the proof in \(\mathbb R^d\) and omit the trivial affine-isometry
reduction.

## 1. Pointwise domination and the exact entropy identity

At almost every source point, the Alexandrov derivative

\[
 A=D_aT=D_a^2\phi
\]

exists and is positive definite.  Since the source and target densities
are respectively \(2e^{-V}\mathbf1_E\) and
\(2e^{-V}\mathbf1_{E^c}\), the Monge--Ampere identity gives

\[
 \det A=e^{V(Tx)-V(x)}.                                 \tag{1.1}
\]

For \(0<t<1\), put \(B_t=(1-t)I+tA\).  The map \(F_t\) is injective on
a full source-measure set: monotonicity of \(T\) gives

\[
 \langle F_t(x)-F_t(x'),x-x'\rangle\ge(1-t)|x-x'|^2.
\]

If \(\rho_t\) is the Lebesgue density of \(\nu_t\), the area formula
therefore gives

\[
 \rho_t(F_tx)\det B_t=2e^{-V(x)}.                       \tag{1.2}
\]

The formula is source-almost everywhere and is enough even when the
distributional Hessian of \(\phi\) has a singular part.  Images created
only by that singular part carry zero \(\nu_t\)-mass.

Define

\[
\begin{aligned}
 \Delta_{V,t}(x,y)
   &=(1-t)V(x)+tV(y)-V((1-t)x+ty),\\
 \Delta_{A,t}(A)
   &=\log\det((1-t)I+tA)-t\log\det A.
\end{aligned}                                           \tag{1.3}
\]

Convexity of \(V\) and scalar weighted AM--GM on the eigenvalues of
\(A\) show that both quantities are nonnegative.  Combining (1.1) and
(1.2), with \(z=F_tx\), yields

\[
 \log\frac{\rho_t(z)}{e^{-V(z)}}
 =\log2-\Delta_{V,t}(x,Tx)-\Delta_{A,t}(A)\le\log2.
\tag{1.4}
\]

Hence \(\nu_t\le2\mu\).  Changing variables through \(F_t\) in the
relative entropy proves the exact identity

\[
 H(t):=\operatorname {Ent}_\mu(\nu_t)
 =\log2-G(t),
 \qquad
 G(t)=\mathbb E[\Delta_{V,t}+\Delta_{A,t}].             \tag{1.5}
\]

Since relative entropy is nonnegative,

\[
 0\le G(t)\le\log2,qquad G(0)=G(1)=0.                 \tag{1.6}
\]

At the midpoint,

\[
 \Delta_{A,1/2}(A)
 =\log\det\frac{I+A}{2A^{1/2}}.                        \tag{1.7}
\]

If \(R=(A-I)(A+I)^{-1}\), then eigenvalue by eigenvalue

\[
 \log\frac{1+\lambda}{2\sqrt\lambda}
 =-\frac12\log\left(1-left(\frac{\lambda-1}{\lambda+1}\right)^2\right)
 \ge\frac12\left(\frac{\lambda-1}{\lambda+1}\right)^2.
\]

Thus

\[
 \boxed{\mathbb E\|R\|_{HS}^2\le2\log2.}             \tag{1.8}
\]

This controls normalized absolutely continuous strain, not displacement.

## 2. What the whole \(t\)-family and Fisher information add

The full family admits an exact second-variation representation.  On a
fixed chord let \(q_{V,x}\) be the nonnegative distributional second
derivative of

\[
 s\longmapsto V(X+sD).
\]

For the matrix term, direct differentiation gives

\[
 -\frac{d^2}{ds^2}\Delta_{A,s}(A)
 =\operatorname {tr}\left[((A-I)B_s^{-1})^2\right].    \tag{2.1}
\]

Define the nonnegative measure on \((0,1)\)

\[
 \mathcal Q(ds)
 =\mathbb E[q_{V,X}(ds)]
  +\mathbb E\operatorname {tr}\left[((A-I)B_s^{-1})^2\right]ds.
\tag{2.2}
\]

Then \(-G''=\mathcal Q\) distributionally.  The Green kernel for the
interval with zero boundary values gives

\[
 G(t)=\int_0^1 g_t(s)\,\mathcal Q(ds),                 \tag{2.3}
\]

where

\[
 g_t(s)=
 \begin{cases}
   s(1-t),&s\le t,\\
   t(1-s),&s\ge t.
 \end{cases}                                           \tag{2.4}
\]

In particular,

\[
 \frac12\int_0^1\min(s,1-s)\,\mathcal Q(ds)
 =G(1/2)\le\log2,                                     \tag{2.5}
\]

and integration over all \(t\) gives

\[
 \frac12\int_0^1s(1-s)\,\mathcal Q(ds)
 =\int_0^1G(t)dt\le\log2.                             \tag{2.6}
\]

Thus the entire interpolation gives a weighted curvature/deformation
budget.  It gives no unweighted control near the endpoints and, more
importantly, it vanishes identically for piecewise translations through
flat potential regions.

There is an equivalent Eulerian reading.  At a regular point
\(z=F_t(x)\), the velocity is

\[
 v_t(z)=D(x),\qquad
 \nabla v_t(z)=(A-I)B_t^{-1}.                           \tag{2.7}
\]

Consequently (2.2) is the exact Otto Hessian of relative entropy:

\[
 H''(t)=\int\left[
 \operatorname {tr}((\nabla v_t)^2)
 +\langle\nabla^2V,v_t,v_t\rangle
 \right]d\nu_t,                                       \tag{2.8}
\]

with the potential term interpreted distributionally when \(V\) is
nonsmooth.

Suppose the relative Fisher information

\[
 I_\mu(\nu_t)=\int\left|\nabla\log\frac{d\nu_t}{d\mu}\right|^2d\nu_t
\]

is finite.  The first variation and Cauchy--Schwarz give

\[
 |H'(t)|^2
 \le I_\mu(\nu_t)\int|v_t|^2d\nu_t
 =I_\mu(\nu_t)\operatorname {tr}K.                    \tag{2.9}
\]

The direction of (2.9) is unusable for an upper bound on
\(\operatorname {tr}K\): entropy gives no upper bound on Fisher
information.  The ordinary Cramer--Rao inequality has the same problem.
For a density \(\rho\) with covariance \(C\), integration by parts gives

\[
 \int (x-m)(\nabla\log\rho)^T\rho\,dx=-I,
\]

and hence its Fisher information matrix dominates \(C^{-1}\).  This is
again a lower bound on Fisher information.  For hard balanced cuts and
for the corridor examples below, the classical Fisher information is
infinite because the densities have jump boundaries, while \(H(t)\) is
constant.  Smoothing the boundaries makes the Fisher information diverge
and supplies no bounded quantity to combine with (2.9).

## 3. Exact covariance identities and the all-time correction

Let

\[
 m=\mathbb EX,qquad \mathbb EY=-m,
\]

the second identity following from
\(\mu=(\mu_++\mu_-)/2\) and centering of \(\mu\).  Put

\[
 Q_+=\mathbb E[XX^T],\quad Q_-=\mathbb E[YY^T],
 \quad C_+=Q_+-mm^T,\quad C_-=Q_--mm^T.
\]

Isotropy gives

\[
 \frac12(Q_++Q_-)=I.                                   \tag{3.1}
\]

For \(Z_t=(1-t)X+tY=X+tD\), elementary expansion proves

\[
 \mathbb E[Z_tZ_t^T]
 =(1-t)Q_++tQ_- -t(1-t)K.                              \tag{3.2}
\]

Since \(\mathbb EZ_t=(1-2t)m\) and
\(\mathbb ED=-2m\), if

\[
 K_c=\operatorname {Cov}(D)=K-4mm^T,
\]

then

\[
 \boxed{C_t:=\operatorname {Cov}(Z_t)
 =(1-t)C_++tC_- -t(1-t)K_c.}                           \tag{3.3}
\]

At \(t=1/2\), (3.1) and the zero midpoint mean give

\[
 C_{1/2}=I-\frac14K.                                   \tag{3.4}
\]

In particular,

\[
 0\preceq K\preceq4I.                                 \tag{3.5}
\]

Jensen applied to \(D\), together with \(\mathbb ED=-2m\), also gives
\(4mm^T\preceq K\).  Hence \(|m|\le1\).  If
\(S=\operatorname {tr}K\), then

\[
 \frac{\operatorname {tr}K}{\|K\|_{op}}\ge\frac S4,
 \qquad
 \operatorname {tr}K_c=S-4|m|^2\ge S-4.               \tag{3.5a}
\]

Thus a hypothetical large cost is necessarily high rank even after the
rank-one barycenter component is removed.  The entropy strain budget
(1.8), by contrast, remains dimension free.

For general \(t\), the exact symmetric raw-moment statement is

\[
 \frac12\left(\mathbb E[Z_tZ_t^T]
                 +\mathbb E[Z_{1-t}Z_{1-t}^T]\right)
 =I-t(1-t)K.                                           \tag{3.6}
\]

The symmetric covariance average has the additional barycenter term

\[
 \frac12(C_t+C_{1-t})
 =I-t(1-t)K-(1-2t)^2mm^T.                              \tag{3.7}
\]

Equations (3.6)--(3.7) are the only valid all-time readings of (0.2)
without extra endpoint symmetry.  A one-coordinate half-cube already
disproves (0.2): at \(t=0\), its conditional covariance is not \(I\).

Also note from (3.3) that

\[
 C_t''=2K_c\succeq0.                                   \tag{3.8}
\]

However,

\[
 \frac{d^2}{dt^2}\log\det C_t
 =\operatorname {tr}(C_t^{-1}C_t'')
  -\operatorname {tr}((C_t^{-1}C_t')^2),               \tag{3.9}
\]

which has no fixed sign.  Optimality of the Brenier coupling does not
repair this sign, as the explicit counterexample in Section 5 shows.

## 4. What a determinant argument would have to prove

### 4.1 A universal spectral floor

We first record two elementary one-dimensional facts with deliberately
generous constants.

**Lemma 4.1 (height of an isotropic log-concave density).**  If \(p\) is
a log-concave probability density on \(\mathbb R\) with variance one,
then \(\|p\|_\infty<10\).

**Proof.**  Translate a mode to zero and write \(M=p(0)\).  On the
positive side, let \(r\) be the first point at which \(p(r)=M/e\); if
the support ends first, use that endpoint and omit the tail below.  Since
\(p\ge M/e\) on \([0,r]\), one has \(r\le e/M\).  Concavity of
\(\log p\) implies, for \(x\ge r\),

\[
 p(x)\le M e^{-x/r}.
\]

Therefore

\[
 \int_0^\infty x^2p(x)dx
 \le Mr^3\left(\frac13+\int_1^\infty u^2e^{-u}du\right)
 =Mr^3\left(\frac13+\frac5e\right).
\]

The same argument on the negative side yields

\[
 \mathbb E(X-0)^2
 \le \frac{2e^3}{M^2}\left(\frac13+\frac5e\right)
 <\frac{88}{M^2}.
\]

Since \(\operatorname {Var}X\le\mathbb E(X-0)^2\) and the variance is
one, \(M<\sqrt{88}<10\). \(\square\)

**Lemma 4.2 (bounded density cannot have tiny variance).**  If a
probability density \(q\) on \(\mathbb R\) satisfies \(q\le L\), then

\[
 \operatorname {Var}(q)\ge\frac1{12L^2}.               \tag{4.1}
\]

**Proof.**  If \(a\) is the mean, then
\(\mathbb P(|X-a|\le r)\le2Lr\).  Hence

\[
 \mathbb E|X-a|^2
 =\int_0^\infty2r\,\mathbb P(|X-a|>r)dr
 \ge\int_0^{1/(2L)}2r(1-2Lr)dr
 =\frac1{12L^2}.\quad\square
\]

Apply these lemmas to a one-dimensional projection.  If \(u\) is a unit
vector, the law of \(u\cdot X_\mu\), \(X_\mu\sim\mu\), is isotropic and
log-concave.  The domination \(\nu_{1/2}\le2\mu\) passes to projections,
so the projected midpoint density is at most \(20\).  Lemma 4.2 gives

\[
 u^TC_{1/2}u\ge\frac1{4800}.                            \tag{4.2}
\]

Together with (3.4)--(3.5), this proves

\[
 c_0I\preceq C_{1/2}\preceq I,
 \qquad c_0=\frac1{4800}.                              \tag{4.3}
\]

### 4.2 Determinant and trace are equivalent here

Let \(a_1,\ldots,a_d\in[c_0,1]\) be the eigenvalues of \(C_{1/2}\).
By (3.4),

\[
 \operatorname {tr}K=4\sum_{i=1}^d(1-a_i).             \tag{4.4}
\]

For \(c_0\le a\le1\),

\[
 1-a\le-\log a\le\frac{1-a}{c_0}.                    \tag{4.5}
\]

Summing proves

\[
 \boxed{
 \frac14\operatorname {tr}K
 \le-\log\det C_{1/2}
 \le\frac1{4c_0}\operatorname {tr}K
 =1200\operatorname {tr}K.}                           \tag{4.6}
\]

Thus a universal lower bound on \(\det C_{1/2}\) would close the trace
problem, and conversely a trace bound would give a universal determinant
lower bound.  The determinant route has not weakened the target.

### 4.3 The true entropy--determinant inequality is dimensionful

Let \(h\) denote differential entropy and let \(\gamma_d\) be standard
Gaussian.  The endpoint entropies satisfy

\[
 \frac12h(\mu_+)+\frac12h(\mu_-)=h(\mu)-\log2.          \tag{4.7}
\]

The Jacobian identity gives at the midpoint

\[
 h(\nu_{1/2})
 =h(\mu)-\log2+\mathbb E\Delta_{A,1/2}.                \tag{4.8}
\]

The Gaussian maximum-entropy inequality at fixed covariance says

\[
 h(\nu_{1/2})
 \le h(\gamma_d)+\frac12\log\det C_{1/2}.
\]

Consequently the valid bound is

\[
 -\log\det C_{1/2}
 \le2D(\mu\|\gamma_d)+2\log2
       -2\mathbb E\Delta_{A,1/2}.                      \tag{4.9}
\]

The Gaussian entropy deficit \(D(\mu\|\gamma_d)\) can be of order \(d\)
for elementary product log-concave measures.  Thus (4.9) cannot give
(0.9).  Removing that term by a universal assertion would amount to new
geometric content, not to another manipulation of the entropy identity.

For a uniform measure on a convex body \(K\), (4.8) is precisely the
entropy form of Brunn--Minkowski:

\[
 h(\nu_t)=\log(|K|/2)+\mathbb E\Delta_{A,t}
 \ge\log(|K|/2).                                      \tag{4.10}
\]

Translations have \(A=I\), so equality in (4.10) is completely
independent of their length.  Brunn--Minkowski volume therefore has the
same translation blindness as (1.5).

## 5. A concrete failure of covariance-determinant concavity

This section gives the promised log-concave counterexample.  Start with
the uniform probability \(\mu_0\) on \([-1,1]\), and put

\[
 a=\frac{\sqrt5-1}{4},\qquad
 b=\frac{\sqrt5+1}{4}=a+\frac12,
 \qquad \ell=1-b=\frac{3-\sqrt5}{4}.                   \tag{5.1}
\]

Define the balanced set

\[
 E_0=[-1,-b]\cup[-a,a]\cup[b,1].                      \tag{5.2}
\]

Its Lebesgue length is one.  It is symmetric, and

\[
 a^3+1-b^3=\frac12.                                    \tag{5.3}
\]

Since the conditional density \(2\mathbf1_{E_0}d\mu_0\) is simply
Lebesgue density one on \(E_0\), (5.3) gives

\[
 \mathbb E_{\mu_{0,+}}X=0,qquad
 \mathbb E_{\mu_{0,+}}X^2
 =\frac23(a^3+1-b^3)=\frac13.                          \tag{5.4}
\]

The complement has the same mean and second moment, because their equal
mixture is \(\mu_0\).  After the dilation \(x\mapsto\sqrt3x\), the
ambient uniform law and **both** conditional endpoint laws have mean zero
and variance one.

The increasing, hence Brenier, map from \(E_0\) to \(E_0^c\) is the
following piecewise translation (endpoints are irrelevant):

\[
 T_0(x)=
 \begin{cases}
  x+\ell,&x\in[-1,-b],\\
  x-a,&x\in[-a,0],\\
  x+a,&x\in[0,a],\\
  x-\ell,&x\in[b,1].
 \end{cases}                                           \tag{5.5}
\]

Indeed these four pieces map, in increasing order, onto
\([-b,-a]\cup[a,b]\).  The derivative is one at every source density
point.  Its unscaled cost is

\[
 K_0=2\ell^3+2a^3=\frac{7-3\sqrt5}{4}.                 \tag{5.6}
\]

After isotropic dilation the cost is

\[
 \kappa=3K_0=\frac{3(7-3\sqrt5)}4>0.                  \tag{5.7}
\]

Since both endpoint covariances are one and both endpoint means are zero,
(3.3) gives for the entire interpolation

\[
 \operatorname {Var}(\nu_t)=1-t(1-t)\kappa.            \tag{5.8}
\]

In particular,

\[
 \operatorname {Var}(\nu_{1/2})
 =1-\frac\kappa4
 =\frac{9\sqrt5-5}{16}<1.                              \tag{5.9}
\]

This disproves

\[
 \det C_{1/2}\ge\sqrt{\det C_+\det C_-},              \tag{5.10}
\]

the logarithmic version (0.7), and the stronger Minkowski-style
determinant-root concavity.

There is more.  The ambient potential is constant on a convex interval,
and \(T_0'=1\) source-almost everywhere.  For every \(0<t<1\), the map
\((1-t)I+tT_0\) has derivative one on each source piece, is increasing,
and its image has length one.  Therefore

\[
 \nu_t=2\mu_0\big|_{F_t(E_0)},
 \qquad \operatorname {Ent}_{\mu_0}(\nu_t)=\log2,      \tag{5.11}
\]

and

\[
 \Delta_{V,t}=\Delta_{A,t}=G(t)=0                     \tag{5.12}
\]

for every \(t\).  Thus any proposed estimate of the form

\[
 -\log\det C_t
 \le C\,G(t)
\]

also fails, even with equal isotropic endpoint covariances.  All the
transport is stored in translation constants, flat target gaps, and the
singular jump of the monotone extension at the origin.

Tensoring this one-dimensional construction with any isotropic uniform
convex factor gives the same counterexample in every dimension.  It does
not produce a growing transport cost; it isolates the precise information
which entropy and determinant concavity fail to see.

## 6. Model tests

### 6.1 Radial Gaussian and radial exponential median balls

The following observation handles both models.  Let \(\mu\) be radial,
let \(R=|X_\mu|\), and let \(E\) be a median ball.  If \(R_0,R_1\) are
the radii conditioned respectively inside and outside the ball, then the
Brenier map preserves the direction and couples the radii monotonically.
Using instead independent \(R_0,R_1\) with the same common direction is
a valid, possibly nonoptimal, coupling.  For \(m=\mathbb ER\),

\[
\begin{aligned}
 W_2^2(\mu_+,\mu_-)
 &\le\mathbb E(R_1-R_0)^2\\
 &\le2\mathbb E(R_1-m)^2+2\mathbb E(R_0-m)^2
 =4\operatorname {Var}R.                              \tag{6.1}
\end{aligned}
\]

Rotational symmetry gives, with \(S=W_2^2(\mu_+,\mu_-)\),

\[
 K=\frac Sd I,qquad
 C_{1/2}=\left(1-\frac{S}{4d}\right)I,qquad
 \det C_{1/2}=\left(1-\frac{S}{4d}\right)^d.           \tag{6.2}
\]

For \(\mu=\gamma_d\), Gaussian Poincare applied to the 1-Lipschitz norm
gives \(\operatorname {Var}R\le1\).  Hence \(S\le4\).  For \(d\ge2\),

\[
 \det C_{1/2}\ge(1-1/d)^d\ge\frac14;                  \tag{6.3}
\]

for \(d=1\), \(\operatorname {Var}|G|=1-2/\pi\) gives an even larger
constant.  Independently, the quadratic potential gives the exact chord
term

\[
 \Delta_{V,t}=\frac12t(1-t)|D|^2,                      \tag{6.4}
\]

so the midpoint entropy identity alone yields \(S\le8\log2\).

For the isotropic radial exponential law

\[
 d\mu(x)=c_d e^{-\sqrt{d+1}|x|}dx,                     \tag{6.5}
\]

the radius is \(\operatorname {Gamma}(d,\sqrt{d+1})\).  Thus

\[
 \mathbb ER^2=d,qquad
 \operatorname {Var}R=\frac d{d+1}.                   \tag{6.6}
\]

Equations (6.1)--(6.2) imply

\[
 S\le\frac{4d}{d+1},
 \qquad
 \det C_{1/2}\ge\left(\frac d{d+1}\right)^d\ge e^{-1}.
\tag{6.7}
\]

Here the radial Brenier chords stay on rays and the potential
\(\sqrt{d+1}|x|\) is affine along each chord.  Therefore

\[
 \Delta_{V,t}=0\quad\text{for every }t,                \tag{6.8}
\]

and the whole one-bit entropy budget lies in the Jacobian AM--GM term.
The cost bound (6.7) comes instead from the one-dimensional radial
variance.  This cleanly demonstrates why potential curvature cannot be
the general mechanism.

For completeness, if \(F\) is the radial distribution function and
\(u\in(0,1)\), the radial coupling is

\[
 r(u)=F^{-1}(u/2),\qquad s(u)=F^{-1}((1+u)/2),          \tag{6.9}
\]

and the Brenier derivative has one radial eigenvalue \(s'(r)\) and
\(d-1\) tangential eigenvalues \(s(r)/r\).  Thus these examples test a
genuinely high-rank Jacobian rather than a hidden one-dimensional affine
map.

### 6.2 Cube corridors

Let \(\mu\) be uniform on the isotropic cube
\([ -\sqrt3,\sqrt3]^d\).  Divide the first coordinate interval into
\(2m\) consecutive slabs of common width

\[
 w=\frac{\sqrt3}{m},
\]

and let \(E_m\) be the union of the even-numbered slabs.  Translation by
\(we_1\) maps \(E_m\) exactly onto its complement, so the Brenier map is

\[
 T(x)=x+we_1.                                           \tag{6.10}
\]

Consequently

\[
 K=w^2e_1e_1^T,qquad \operatorname {tr}K=\frac3{m^2}.
\tag{6.11}
\]

The potential is constant, \(A=I\), and every interpolated image is a
union of slabs of half the cube volume.  Hence for every \(t\),

\[
 \nu_t=2\mu|_{F_t(E_m)},\qquad H(t)=\log2,qquad G(t)=0.
\tag{6.12}
\]

At \(m=1\) this is the half-cube translation of cost three.  Increasing
the number of corridors decreases the cost like \(m^{-2}\); entropy and
Brunn--Minkowski remain at exact equality and do not register either
scale.  Also, the common translation leaves covariance unchanged, and a
mixture calculation gives

\[
 C_t=I-\frac{w^2}{4}e_1e_1^T\quad\text{for every }t.    \tag{6.12a}
\]

This agrees with \(I-t(1-t)K\) only at \(t=1/2\), giving the simplest
explicit check that (0.2) is not an all-time covariance identity.

### 6.3 An isotropic simplex cap

Let \(v_0,\ldots,v_d\) be the vertices of a regular isotropic simplex,
normalized by

\[
 |v_i|^2=d(d+2),\qquad v_i\cdot v_j=-(d+2)\quad(i\ne j).
\tag{6.13}
\]

Write \(X=\sum_{i=0}^dp_iv_i\), with \(p\) uniform Dirichlet, and put

\[
 b=2^{-1/d},\qquad a=1-b,qquad E=\{p_0\ge a\}.         \tag{6.14}
\]

Since \(p_0\sim\operatorname {Beta}(1,d)\), this cap has mass one half.
Its conditional barycenter is \(m_+=av_0\), while the complementary
barycenter is \(-m_+\).  Therefore

\[
 W_2^2(\mu_+,\mu_-)
 \ge4|m_+|^2
 =4d(d+2)(1-2^{-1/d})^2
 \longrightarrow4(\log2)^2.                           \tag{6.15}
\]

There is a direct dimension-free upper coupling.  Conditional on
\(p_0=s\), write

\[
 X=sv_0+(1-s)W,
\]

where \(W\) is uniform on the opposite facet and independent of \(s\).
Use the same \(W\) for the cap and complement, and couple their scalar
\(s\)-laws in any way no worse than the independent coupling.  A direct
Dirichlet calculation gives

\[
 \mathbb E|v_0-W|^2=(d+2)(d+3).                        \tag{6.16}
\]

For the cap scalar \(s_+\),

\[
 \operatorname {Var}s_+
 =\frac{b^2d}{(d+1)^2(d+2)}\le\frac1{d^2}.             \tag{6.17}
\]

The complementary scalar is supported on \([0,a]\), so
\(\operatorname {Var}s_-\le a^2/4\).  Their means differ by
\(2da/(d+1)\), and \(a\le(\log2)/d\).  Hence

\[
 \mathbb E(s_+-s_-)^2
 \le\frac{1+(17/4)(\log2)^2}{d^2}.                     \tag{6.18}
\]

Since \((d+2)(d+3)/d^2\le12\), (6.16)--(6.18) give

\[
 W_2^2(\mu_+,\mu_-)<37.                                \tag{6.19}
\]

This uses no spectral or KLS estimate.  The simplex potential is constant
on its convex support, so \(\Delta_{V,t}=0\) for all \(t\); the entropy
budget again sees only normalized Brenier strain.

### 6.4 Product exponential cut by the maximum

Let \(Z_1,\ldots,Z_d\) be independent mean-one exponentials and set
\(X_i=Z_i-1\).  This product law is isotropic and has affine potential

\[
 V(x)=\sum_{i=1}^d(x_i+1)
 \quad\text{on }[-1,\infty)^d.                         \tag{6.20}
\]

Let \(M=\max_iX_i\) and choose its median

\[
 m_d=-1-\log(1-2^{-1/d}),
 \qquad E=\{M\le m_d\}.                               \tag{6.21}
\]

Because \(M\) is 1-Lipschitz, Kantorovich duality gives

\[
 W_2(\mu_+,\mu_-)
 \ge W_1(\mu_+,\mu_-)
 \ge2\mathbb E|M-m_d|.                                \tag{6.22}
\]

After centering by \(\log d-1\), the maximum converges with uniformly
integrable first moments to a Gumbel random variable.  The last expression
therefore converges to a finite positive constant.  Thus this is not a
vanishing-cost cut.

For a direct upper bound, we use two elementary facts.  First, the
mean-one exponential law has Poincare constant at most four.  Indeed, for
a mean-zero smooth \(g\), put

\[
 G(x)=\int_x^\infty g(s)e^{-s}ds.
\]

Then \(G(0)=0\), \(G'=-ge^{-x}\), and integration by parts gives

\[
 \int g^2e^{-x}dx=\int g'Gdx.
\]

The elementary weighted Hardy estimate

\[
 \int G^2e^x dx\le4\int(G')^2e^x dx                  \tag{6.23}
\]

follows by expanding
\(\int(G'+G/2)^2e^x dx\ge0\).  Cauchy--Schwarz now yields
\(\operatorname {Var}(g)\le4\int(g')^2e^{-x}dx\).
Tensorization preserves the constant four.

Second, if a probability \(\sigma=f\mu\) and \(\mu\) has Poincare
constant \(C_P\), then

\[
 W_2(\sigma,\mu)
 \le2\sqrt{C_P\int(f-1)^2d\mu}.                       \tag{6.24}
\]

To see this, solve weakly
\(\int\nabla u\cdot\nabla h\,d\mu=\int(f-1)h\,d\mu\).
Poincare gives
\(\int|\nabla u|^2d\mu\le C_P\int(f-1)^2d\mu\).
Along \(\rho_s=1+s(f-1)\), use flux \(\rho_sv_s=\nabla u\).  Since
\(\rho_s\ge1-s\), integration of the metric speed gives (6.24).

For either balanced half, \(f=2\mathbf1_E\) and the chi-square term is
one.  The triangle inequality and \(C_P\le4\) yield

\[
 W_2(\mu_+,\mu_-)\le8,
 \qquad \operatorname {tr}K\le64.                     \tag{6.25}
\]

This is an explicit product estimate, not an invocation of KLS.  Since
the potential (6.20) is affine on every chord in its convex support,

\[
 \Delta_{V,t}=0\quad\text{for all }t.                  \tag{6.26}
\]

The complicated maximum cut therefore spends its entire entropy budget
on the matrix term, while its cost is controlled by the independent
product Poincare calculation.

### 6.5 General one-dimensional piecewise translations

The phenomenon in Section 5 is not exceptional.  Let \(\mu\) be uniform
on an interval and let \(E\) be any measurable half-partition.  The
increasing rearrangement between the two conditional laws has

\[
 T'(x)=1
\]

at almost every source density point: both conditional Lebesgue densities
are equal.  The map can still have flat extensions across target gaps and
singular jumps between source components.  For every \(t\),
\(F_t'=1\) source-almost everywhere and all chords remain in the ambient
interval.  Hence

\[
 G(t)=0\quad(0\le t\le1),                              \tag{6.27}
\]

although the cost is generally positive.  One-dimensional isotropy and
(3.5) still give the crude sharp-form bound \(K\le4\), but none of it is
recorded by the entropy, curvature, Fisher, or Brunn--Minkowski gaps.

## 7. The exact remaining gap

The calculations leave two equivalent formulations of the missing
statement:

\[
 \sup_{\substack{\mu\ \mathrm{isotropic\ log\!\!-concave}\\
                  \mu(E)=1/2}}
 \mathbb E_{\mu_+}|T(X)-X|^2<\infty,                   \tag{7.1}
\]

or

\[
 \inf_{\substack{\mu\ \mathrm{isotropic\ log\!\!-concave}\\
                  \mu(E)=1/2}}
 \det\operatorname {Cov}(\nu_{1/2})>0.                \tag{7.2}
\]

Their quantitative equivalence is (4.6).  The full interpolation supplies
only

\[
 \mathbb E\left[
 \frac{V(X)+V(Y)}2-V\left(\frac{X+Y}2\right)
 +\log\det\frac{I+A}{2A^{1/2}}
 \right]\le\log2,                                    \tag{7.3}
\]

and its weighted second-variation refinements (2.5)--(2.6).  These terms
do not contain the integration constants of \(T-I\) on disconnected
source cells, the flat behavior across target gaps, or singular jumps of
the Brenier derivative.  Section 5 proves that those missing pieces can
change covariance while every quantity in the entire entropy family is
exactly zero.

For orientation, (7.1) is the balanced-half transport form of KLS.  A
Poincare inequality with constant \(C_P\), together with (6.24), gives

\[
 W_2(\mu_+,\mu_-)\le4\sqrt{C_P}.                       \tag{7.4}
\]

Conversely, if (7.1) is bounded by \(L^2\), take a 1-Lipschitz function
\(f\), choose a median half \(E\subset\{f\le m\}\), and use
Kantorovich duality to obtain

\[
 \mathbb E|f-m|
 \le\frac12W_1(\mu_+,\mu_-)
 \le\frac L2.                                         \tag{7.5}
\]

For log-concave measures this median \(L^1\) concentration scale is
equivalent, up to numerical constants, to the Cheeger/Poincare scale.
This final standard equivalence is mentioned only to identify the theorem
strength of the missing step; it is not used anywhere above.

Accordingly, the unresolved lemma cannot be stated as another local
entropy or covariance manipulation.  It must be a genuinely global
coercivity principle for the singular-translation branch of balanced
Brenier transport.  Proving that principle with a universal constant
would prove the desired KLS bound; assuming it would assume the target.
