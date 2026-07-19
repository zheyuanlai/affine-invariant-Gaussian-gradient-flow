# Ray-congruence rigidity: the endpoint-free theorem, an exact fan, and the incidence gap

## 0. Outcome

This note tests the mixed-BV ray-congruence target on the actual normal
cells of one signed-distance interface.  There is a complete theorem in the
literal no-contact branch.  If the normal exponential map has no cut,
focal, or support endpoint on either side of any basepoint, then the
interface is one affine hyperplane.  For an isotropic probability its
half-mass signed-distance moment is at most \(\sqrt2\).

On a long balanced cell of reciprocal median density \(s\), the normal
interval contains \((-s/4,s/4)\), every principal curvature is at most
\(4/s\), and the shape part of the completed curvature is at least

\[
                 {s\over36}\,\|S\|_{HS}^{2}.             \tag{0.1}
\]

Thus the fixed-mass core in
`singular_eikonal_curvature.md` satisfies, with explicit universal
constants,

\[
 \int_{\mathcal G}\|S_y\|_{HS}^{2}\,d\eta(y)
       \le {36A\over aD^{2}}                              \tag{0.2}
\]

whenever \(s_y\ge aD\) and \(\mathcal C_{\mathcal G}\le A/D\).
There is also an exact two-point reach inequality, (2.8) below.  These are
the strongest consequences of the smooth normal Jacobian which do not
already require a global incidence theorem.

The polyhedral/contact branch cannot be discarded.  For the uniform law on
the isotropic cube, the parity set has a cellwise balanced coordinate-fan
interface and

\[
 Q={I_n\over n},\qquad
 J={\sqrt3\over n+1},\qquad
 P={n\over2\sqrt3},\qquad
 (2P)J={n\over n+1}.                                     \tag{0.3}
\]

All its curvature is singular transition incidence.  This is an actual
log-concave normal congruence, not an abstract ray star, and it proves that
linear BV turning and its exact scale are unavoidable.

What remains is a dimension-free positive-incidence theorem.  A convenient
form is

\[
 \boxed{\quad
 \mathcal C_{\mathcal G}
       +\lambda_{\max}(Q_{\mathcal G})^{1/2}
       \ge c(A,\beta)>0 .\quad}                           \tag{CGI}
\]

Here the hypotheses are precisely the fixed quotient mass, long balanced
cells, basepoint covariance bound, and global log-concave normal-congruence
hypothesis of MBR.  Since
\(\lambda_{\max}(Q_{\mathcal G})\le\|Q_{\mathcal G}\|_{HS}
\le C/D\), (CGI) and \(\mathcal C_{\mathcal G}\le A/D\)
would give \(D=O_{A,\beta}(1)\).  The orthogonal abstract star violates
(CGI), while the hyperplane, radial, and parity models below satisfy it.
No proof of (CGI) is obtained here.  The obstruction is exact: stable rank
selects isolated directions but gives no positive-mass narrow direction
packets, and Prekopa--Leindler fills all cross-combinations rather than the
same-ray pairing.  Turning that filling into (CGI) is the remaining
conjecture-strength step.

## 1. Setting and conventions

Let \(\Sigma\subset\mathbb R^d\) be a closed \(C^2\) embedded hypersurface
without boundary, let \(N(y)\) be a chosen unit normal, and let

\[
                 F(y,t)=y+tN(y).                          \tag{1.1}
\]

The regular normal cell \(I_y=(a_y,b_y)\), with \(a_y<0<b_y\), is the
maximal interval on which \(y\) is the unique nearest point, the normal
Jacobian is positive, and the ray remains in the convex support of the
ambient log-concave density.  If \(S_y=D_\Sigma N(y)\), then

\[
 J_y(t)=\det(I+tS_y)>0\quad(t\in I_y).                    \tag{1.2}
\]

The normalized conditional density on the cell is

\[
 q_y(t)={1\over m(y)}e^{-V(y+tN(y))}J_y(t)1_{I_y}(t),     \tag{1.3}
\]

and is log-concave.  We assume it is balanced at zero and write

\[
 \int_{a_y}^{0}q_y=\int_0^{b_y}q_y={1\over2},\qquad
 s_y={1\over q_y(0)}.                                    \tag{1.4}
\]

Let \(W_y\) be the two-sided tail transport.  The completed cell curvature
is

\[
 \mathcal C_y=\int W_y\,dD^2(-\log q_y)+e_y^-+e_y^+,
 \qquad \mathcal C_y={2\over s_y}.                        \tag{1.5}
\]

In a smooth chart, the part of the first term arising from the normal
Jacobian is

\[
 \mathcal C_y^{\rm shape}
 =\int_{I_y}W_y(t)\sum_{j=1}^{d-1}
       {\kappa_j(y)^2\over(1+t\kappa_j(y))^2}\,dt.        \tag{1.6}
\]

Convexity of \(V\) makes the remaining interior Hessian measure
nonnegative.

## 2. Exact local consequences of a long balanced cell

### Lemma 2.1 (two-sided geometric length)

Under (1.3)--(1.4),

\[
                    -a_y\ge {s_y\over4},\qquad
                     b_y\ge {s_y\over4}.                 \tag{2.1}
\]

#### Proof

For a one-dimensional log-concave probability with median zero,
\(\|q_y\|_\infty\le2q_y(0)=2/s_y\).  Therefore

\[
 {1\over2}=\int_{a_y}^{0}q_y(t)dt
       \le {2\over s_y}(-a_y),
\]

and the right half is identical.  \(\square\)

### Corollary 2.2 (pointwise reach and curvature)

Every principal curvature obeys

\[
                         |\kappa_j(y)|\le {4\over s_y}.   \tag{2.2}
\]

#### Proof

By Lemma 2.1, \(I+tS_y\) is positive definite for
\(|t|<s_y/4\).  For an eigenvalue \(\kappa_j\), positivity of
\(1+t\kappa_j\) at both ends gives (2.2).  \(\square\)

### Lemma 2.3 (shape-energy lower bound)

The shape contribution (1.6) satisfies

\[
                 \mathcal C_y^{\rm shape}
                    \ge {s_y\over36}\|S_y\|_{HS}^{2}.    \tag{2.3}
\]

#### Proof

For \(0\le t\le s_y/8\),

\[
 W_y(t)={1\over2}-\int_0^tq_y(r)dr
       \ge {1\over2}-{2t\over s_y}\ge {1\over4};         \tag{2.4}
\]

the same estimate holds on \([-s_y/8,0]\).  By (2.2), on this interval

\[
                    0<{1\over2}\le1+t\kappa_j\le{3\over2}.
                                                                    \tag{2.5}
\]

Restricting (1.6) to this interval of length \(s_y/4\) gives, for each
\(j\),

\[
 \int_{-s_y/8}^{s_y/8}W_y(t)
       {\kappa_j^2\over(1+t\kappa_j)^2}dt
 \ge {s_y\over4}\,{1\over4}\,{4\over9}\kappa_j^2
 ={s_y\over36}\kappa_j^2.                               \tag{2.6}
\]

Summing proves (2.3).  \(\square\)

In particular, if \(s_y\ge aD\) on \(\mathcal G\), then

\[
 {aD\over36}\int_{\mathcal G}\|S_y\|_{HS}^2d\eta
 \le\int_{\mathcal G}\mathcal C_y d\eta
 \le {A\over D},                                        \tag{2.7}
\]

which is (0.2).

### Lemma 2.4 (two-point normal reach inequality)

Let \(y,z\in\Sigma\), and suppose the normal interval at \(y\) contains
\((-L_y,L_y)\).  Then

\[
       |\langle N(y),z-y\rangle|
                 \le {|z-y|^2\over2L_y}.                 \tag{2.8}
\]

If the analogous assertion holds at \(z\), then

\[
 |\langle N(y)-N(z),z-y\rangle|
 \le {|z-y|^2\over2}\left({1\over L_y}+{1\over L_z}\right).
                                                                    \tag{2.9}
\]

For balanced cells one may take \(L_y=s_y/4\).

#### Proof

Unique nearest projection of \(y+tN(y)\) to \(y\) implies, for every
\(|t|<L_y\),

\[
 |y+tN(y)-z|^2-t^2
   =|z-y|^2-2t\langle N(y),z-y\rangle\ge0.                \tag{2.10}
\]

Letting \(t\uparrow L_y\) and \(t\downarrow-L_y\) proves (2.8).
Apply it at both points and use the triangle inequality to obtain (2.9).
\(\square\)

The estimate controls only the component of the normal difference along
the chord.  Replacing its left side by \(|N(y)-N(z)|\) is false for skew
normal lines in dimension at least three.  Summing (2.8) directly costs
\(\mathbb E|Y-Y'|^2=2\operatorname{tr}\operatorname{Cov}Y\), which is
dimension dependent; the operator bound on the covariance does not remove
that trace.

## 3. Complete rigidity in the endpoint-free branch

### Theorem 3.1 (global normal exponential implies a hyperplane)

Assume that \(\Sigma\) is a closed \(C^2\) embedded hypersurface without
boundary and that, for every \(y\in\Sigma\), the normal cell is all of
\(\mathbb R\): \(I_y=\mathbb R\).  Equivalently, there is no focal value,
no cut value, and no support endpoint in either normal direction.  Then
\(\Sigma\) is one affine hyperplane and \(N\) is constant up to the global
choice of orientation.

#### Proof

For every eigenvalue \(\kappa_j(y)\), (1.2) says

\[
                         1+t\kappa_j(y)>0
                         \quad\hbox{for every }t\in\mathbb R.
\]

Hence \(\kappa_j(y)=0\).  Thus \(S_y=0\) everywhere.  On each connected
component the Gauss map is constant, and the component is an open subset
of an affine hyperplane.  Since the component is also closed and has no
boundary, it is the whole hyperplane.

Two different codimension-one affine hyperplanes are either nonparallel
and intersect, contradicting embedded disjoint components, or are
parallel.  In the parallel case the perpendicular bisector is a finite
cut locus for the inward normal rays from both planes.  This contradicts
\(I_y=\mathbb R\).  Hence there is exactly one component.  \(\square\)

### Corollary 3.2 (isotropic moment bound)

Let \(\mu\) be centered and isotropic, and suppose the interface in
Theorem 3.1 has half mass on either side.  Then

\[
                         \int d(x,\Sigma)d\mu(x)\le\sqrt2. \tag{3.1}
\]

#### Proof

Write \(\Sigma=\{x:\langle x,N\rangle=a\}\) and
\(Z=\langle X,N\rangle\).  Then \(\mathbb EZ=0\),
\(\operatorname{Var}Z=1\), and \(a\) is a median.  If \(a>0\), Cantelli's
inequality gives

\[
 {1\over2}\le\mathbb P\{Z\ge a\}\le {1\over1+a^2},
\]

so \(a\le1\); applying the same argument to \(-Z\) treats \(a<0\).
Consequently

\[
 \mathbb E|Z-a|\le\bigl(\mathbb E(Z-a)^2\bigr)^{1/2}
 =\sqrt{1+a^2}\le\sqrt2.                                 \tag{3.2}
\]

This is exactly the signed-distance moment.  \(\square\)

The radial sphere is not a counterexample to Theorem 3.1: its inward
normal interval has a focal endpoint at the center, even though the focal
endpoint has zero linear endpoint charge.

## 4. An exact realized high-rank polyhedral fan

Let \(n\ge2\), and let

\[
 K=[-a,a]^n,\qquad a=\sqrt3,
\]

and let \(\mu\) be normalized Lebesgue measure on \(K\).  Its coordinates
are independent, centered, and have variance \(a^2/3=1\), so \(\mu\) is
isotropic and log-concave.  Ignore the coordinate hyperplanes, which have
zero mass, and define the parity set

\[
 E=\left\{x\in K:\prod_{i=1}^n\operatorname{sign}(x_i)=+1\right\}.
                                                                    \tag{4.1}
\]

Independent uniform signs give \(\mu(E)=1/2\), and

\[
                         \partial E\cap\operatorname{int}K
                           =\bigcup_{i=1}^n\{x_i=0\}.     \tag{4.2}
\]

### Proposition 4.1 (normal cells and balance)

At a regular basepoint \(y\) of the facet \(y_i=0\), put

\[
                         b(y)=\min_{j\ne i}|y_j|.
\]

The normal Voronoi cell based at \(y\) is

\[
                         \{y+te_i:-b(y)<t<b(y)\},         \tag{4.3}
\]

up to reversal of \(e_i\) according to the orientation into \(E\).  Its
conditional density is uniform and even; hence every regular cell is
bisected at \(t=0\).

#### Proof

The distance from \(y+te_i\) to the facet \(x_i=0\) is \(|t|\).  Its
distance to \(x_j=0\) is \(|y_j|\).  Thus the nearest boundary point is
uniquely \(y\) precisely when \(|t|<\min_{j\ne i}|y_j|\).  Since
\(b(y)\le a\), support truncation occurs only on a null set where all
remaining coordinates lie on the cube boundary.  The ambient density is
constant, the facet is flat, and the interval is symmetric, proving
balance.  \(\square\)

### Proposition 4.2 (exact moment, perimeter, tensor, and charge)

For the parity fan,

\[
 \begin{aligned}
 J_\mu(E)&={a\over n+1},\\
 P_\mu(E)&={n\over2a},\\
 Q&={I_n\over n},\\
 \mathfrak I&=\mathfrak S=0,
 \qquad \mathfrak J=2P_\mu(E)={n\over a}.
 \end{aligned}                                           \tag{4.4}
\]

In particular \(\mathfrak J J_\mu(E)=n/(n+1)\).

#### Proof

For \(X\sim\mu\), the variables \(|X_i|\) are independent uniform on
\([0,a]\).  Equation (4.2) gives

\[
 d(X,\partial E)=\min_i|X_i|,
\]

and the tail formula yields

\[
 J_\mu(E)=\int_0^a\left(1-{t\over a}\right)^n dt
          ={a\over n+1}.                                 \tag{4.5}
\]

The exterior \(\varepsilon\)-tube adds only the complement half of the
two-sided tube.  Independence of signs and magnitudes gives

\[
 \mu(E_\varepsilon)-\mu(E)
 ={1\over2}\left[1-\left(1-{\varepsilon\over a}\right)^n\right]
 \quad(0<\varepsilon<a),                                 \tag{4.6}
\]

so \(P_\mu(E)=n/(2a)\).  Hyperoctahedral symmetry and
\(\operatorname{tr}Q=1\) force \(Q=I_n/n\).

Every regular facet is flat and \(V\) is constant in the interior, so the
regular interior charge vanishes.  As observed in Proposition 4.1, support
escape occurs only on a quotient-null set.  The completed Stieltjes
identity therefore puts the full charge \(2P\) on the medial transition
strata.  This proves (4.4).  \(\square\)

This example simultaneously rules out all of the following proposed
shortcuts:

* a claim that diffuse \(Q\) forces regular curvature;
* deletion of medial incidence in the polyhedral limit;
* replacement of the linear jump charge by an angle-squared charge;
* an argument which treats cellwise balance as a source of directional
  coherence.

## 5. Other mandatory stress tests

### 5.1 Cube halfspace

For the same cube and \(E=\{x_1>0\}\),

\[
 J={a\over2}={\sqrt3\over2},\qquad
 P={1\over2a}={1\over2\sqrt3},\qquad
 Q=e_1\otimes e_1.                                      \tag{5.1}
\]

There is no regular or medial turning.  The two support endpoints carry
\(\mathfrak S=2P=1/a\).  Thus a positive-incidence theorem must allow the
rank-one/contact alternative.

### 5.2 Isotropic radial exponential

Let \(d\mu=c_ne^{-\alpha|x|}dx\), where
\(\alpha=\sqrt{n+1}\), and let \(r_n\) be the median of \(R=|X|\).
Then \(R\) is Gamma with shape \(n\) and rate \(\alpha\), and
\(\mathbb ER^2=n\), so \(\mu\) is isotropic.  The median sphere has

\[
 Y=r_nN,qquad Q={I_n\over n},qquad
 \operatorname{Cov}(Y)={r_n^2\over n}I_n.                \tag{5.2}
\]

Moreover \(\operatorname{Var}R=n/(n+1)<1\) and the distance between the
mean and a median is at most one standard deviation.  Hence

\[
             J=\mathbb E|R-r_n|\le\sqrt2.                \tag{5.3}
\]

There is no codimension-one medial or support charge.  The inward endpoint
is focal, and the normal-Jacobian curvature supplies the charge.  This
shows why high-rank \(Q\), bounded basepoint covariance, and a smooth
connected interface do not imply near-affinity; the correct conclusion is
only that large cell scale is impossible.

### 5.3 Product maximum

For independent variance-one symmetric exponentials and the median box
\(\{\max_i|X_i|<r_n\}\), hyperoctahedral symmetry gives diffuse normalized
direction tensor, while the exact extreme-value calculation gives
\(J=O(1)\).  Its regular cells are cut when a competing coordinate becomes
active.  It is therefore a contact-incidence model, not evidence for a
smooth no-contact theorem.  It also warns that replacing the actual cut
endpoint by the support endpoint changes the cell aggregate by an
unbounded factor.

## 6. Why stable rank plus convexification does not yet prove incidence

On the long core, let

\[
 Q=\int N\otimes N\,d\nu,
 \qquad \operatorname{tr}Q=1,
 \qquad \|Q\|_{HS}\le {C\over D}.                        \tag{6.1}
\]

Linear-algebraic restricted-invertibility statements can select
\(k\gtrsim D^2\) individual directions, with weights, on a well-conditioned
subspace.  This is not a packet theorem.  A selected direction is a
\(\nu\)-null point, and (6.1) gives no lower bound on the mass of a narrow
cap around it.  The uniform law on a high-dimensional sphere makes this
failure quantitative: \(Q=I/n\), yet every cap of fixed angular radius
less than \(\pi/2\) has exponentially small mass.  This is an actual
normal-congruence direction law: it is the Gauss law of the isotropic
radial-exponential model in Section 5.2.  Consequently one cannot
obtain \(D^2\) almost-orthogonal packets, each of mass \(c/D^2\), from
stable rank alone.

There is nevertheless fixed conditional mass at macroscopic distance on
each long cell.  Indeed, \(q_y\le2/s_y\) implies

\[
 \int_{s_y/8}^{b_y}q_y(t)dt\ge {1\over4},\qquad
 \int_{a_y}^{-s_y/8}q_y(t)dt\ge {1\over4}.               \tag{6.2}
\]

Thus a quotient packet of mass \(p\) generates positive and negative
ambient sets of mass at least \(p/4\), separated along its normals at
scale \(D\).  If \(A_i^+,A_i^-\) are such sets, Prekopa--Leindler gives

\[
 \mu\left({A_i^++A_j^-\over2}\right)
       \ge\sqrt{\mu(A_i^+)\mu(A_j^-)}.                   \tag{6.3}
\]

The difficulty is that (6.3) concerns every cross-combination, not the
same-ray midpoint \(Y\).  Different midpoint sets can overlap completely,
and their transverse positions contain the uncontrolled trace of the
basepoint covariance.  The operator bound
\(\operatorname{Cov}_{op}Y\le A\) controls every fixed projection but not
\(\mathbb E|Y|^2\).  No valid summation of (6.3) was found which avoids
that trace.

This is also where the convexified orthogonal star changes character.  The
convex hull of \(\{\pm De_i:1\le i\le k\}\) is \(DB_1^k\).  Log-concavity
does force mass into its mixed-coordinate region, but the star hypotheses
provide neither positive-mass axial neighborhoods nor disjointness of the
many midpoint images.  Treating the selected atoms as packets silently
assumes precisely the missing incidence theorem.

## 7. The exact remaining geometric theorem

Here is a form which is strictly adapted to the estimates actually proved.

> **Core Gauss-incidence inequality (CGI).**  For every \(A<\infty\),
> \(\beta>0\), and fixed \(0<a<b<\infty\), there is
> \(c=c(A,\beta,a,b)>0\) with the following property.  Let \(\mu\) be an
> isotropic log-concave probability, and let \(\Sigma\) be one half-mass
> signed-distance interface whose regular normal cells are balanced.  Let
> \(\mathcal G\) have quotient mass at least \(\beta\), suppose
> \(aD\le s_y\le bD\) on \(\mathcal G\), and suppose
> \(\operatorname{Cov}_{op}(Y\mid\mathcal G)\le A\).  Normalize the
> transport-weighted Gauss law on \(\mathcal G\) to \(Q_{\mathcal G}\),
> and include regular Hessian curvature, medial incidence with
> multiplicity, and support escape in \(\mathcal C_{\mathcal G}\).  Then
> \[
> \mathcal C_{\mathcal G}
> +\lambda_{\max}(Q_{\mathcal G})^{1/2}\ge c.             \tag{7.1}
> \]

The endpoint-free theorem proves the geometric zero-contact branch of this
dichotomy: absence of every endpoint forces \(\lambda_{\max}(Q)=1\).  The
parity fan proves that the second term may be small only while transition
incidence is large.  The cube halfspace proves that the rank-one branch may
place all charge at support.  The radial exponential proves that regular
normal-Jacobian curvature is a genuine third possibility.

Under the long-core estimates,

\[
 \mathcal C_{\mathcal G}\le {A\over D},\qquad
 \lambda_{\max}(Q_{\mathcal G})^{1/2}
 \le\|Q_{\mathcal G}\|_{HS}^{1/2}\le {C\over\sqrt D}.    \tag{7.2}
\]

Therefore (7.1) gives a universal bound on \(D\).  Conversely, without the
requirement that the rays arise from one interface in one log-concave
ambient measure, the orthogonal star has
\(\mathcal C\asymp D^{-1}\), \(Q=I_{D^2}/D^2\), and violates (7.1).

The unproved content of (CGI) is thus exactly a quantitative statement
that log-concave convexification plus Frobenius integrability cannot hide
the transition strata connecting many long, almost orthogonal normal
packets.  Neither the normal Jacobian, the BV endpoint identity, stable
rank, nor basepoint covariance separately supplies it.

## 8. Approximation audit

There are two different approximation questions, and only one is benign.

### 8.1 Smoothing the ambient cube preserves the balanced fan exactly

Let \(U\) be uniform on \([-\sqrt3,\sqrt3]\), let \(G\) be standard
Gaussian, and put

\[
 Z_\varepsilon={U+\varepsilon G\over\sqrt{1+\varepsilon^2}}.
                                                                    \tag{8.1}
\]

Its density \(p_\varepsilon\) is positive, even, \(C^\infty\), and
log-concave, and \(Z_\varepsilon\) has variance one.  Hence
\(\mu_\varepsilon=p_\varepsilon^{\otimes n}dx\) is smooth, full-support,
isotropic, and log-concave.  Use the same parity set (4.1).  At a regular
basepoint of \(x_i=0\), the normal interval is still

\[
 \left(-\min_{j\ne i}|y_j|,\ \min_{j\ne i}|y_j|\right),  \tag{8.2}
\]

and its conditional density is proportional to the even function
\(p_\varepsilon(t)\).  Thus every cell remains exactly balanced; no limiting
Euler equation is being assumed.  Coordinate symmetry again gives
\(Q_\varepsilon=I_n/n\), while

\[
 J_\varepsilon=\mathbb E\min_i|Z_{\varepsilon,i}|,
 \qquad P_\varepsilon=n p_\varepsilon(0).                \tag{8.3}
\]

The perimeter formula follows by differentiating

\[
 {1\over2}\left[1-\mathbb P\{|Z_\varepsilon|>r\}^{,n}\right]
\]

at \(r=0\).  As \(\varepsilon\downarrow0\), bounded convergence for the
minimum and pointwise convergence of the one-dimensional densities give

\[
 J_\varepsilon\longrightarrow{\sqrt3\over n+1},
 \qquad P_\varepsilon\longrightarrow{n\over2\sqrt3}.     \tag{8.4}
\]

Thus neither hard support nor a nonsmooth potential is responsible for the
high-rank contact phenomenon.  With the standard Gaussian product itself,
the same construction gives the exact smooth-density formulas

\[
 P=n(2\pi)^{-1/2},\qquad
 J=\int_0^\infty[2(1-\Phi(t))]^n dt,qquad Q={I_n\over n}, \tag{8.5}
\]

and \(2PJ\to1\) as \(n\to\infty\).  The last limit follows after the
change of variables \(t=u/(2n\phi(0))\) and the local expansion
\(2(1-\Phi(t))=1-2\phi(0)t+O(t^3)\), with the tail dominated by an
exponential using log-concavity of the Gaussian survival function.

### 8.2 Smoothing the interface does not preserve the hypothesis

For each fixed dimension, the parity set can of course be approximated in
\(L^1\) and strictly in weighted BV by sets with smooth boundary.  The
distance functions then converge locally uniformly under a Hausdorff
rounding, so their first moments converge as well.  This does **not** give a
smooth balanced-cell approximation.  Resolving an intersection of two
coordinate facets chooses one of two pairings of the four incident
orthants; the reflected normal intervals cease to be symmetric, and the
cellwise Euler equation generally fails.  A subsequent volume correction
restores only global half mass, not balance on almost every normal cell.

Nor may the endpoint charge simply be identified with the limit of a
quadratic smooth curvature energy.  At the critical transport weight the
quadratic energy of a rounded Gauss jump depends on the rounding profile,
whereas the compactified-ray endpoint charge is the canonical linear BV
trace.  Therefore a legitimate approximation theorem for (CGI) must carry
the ray graph and its incidence measure through the limit.  Ordinary
smooth finite-perimeter approximation is insufficient.

The density smoothing in Section 8.1 is exact and harmless; the interface
smoothing is a genuine unresolved part of any proof which establishes
(CGI) only for smooth embedded interfaces.
