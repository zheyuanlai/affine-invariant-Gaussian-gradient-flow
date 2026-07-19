# What the equimeasurable BV minimizer actually satisfies

## 0. Outcome

Let \(G\) minimize

\[
 {\cal J}_\kappa(G)=|DG|_\mu+\kappa\|M(G)-M(F)\|_*
                                                        \tag{0.1}
\]

in the equimeasurability class of \(F\), where \(0<\kappa<1/3\).
There is a dimension-free BV first-variation theorem: one can choose a
**single compatible** nuclear subgradient \(H\), and the weighted
anisotropic first variation of the whole coarea varifold of \(G\) vanishes
under every compactly supported \(\mu\)-preserving flow.  This is proved in
Theorem 2.1 below.

That theorem is strictly weaker than

\[
 -\operatorname {div}_\mu D\Phi_H(DG)=\lambda(G)       \tag{0.2}
\]

and weaker than constant anisotropic mean curvature of almost every level.
Equation (0.2) is valid on a smooth noncritical band, by exact quantile
restoration (Theorem 3.1).  For a general BV minimizer, however, jump
interfaces are moved simultaneously for every level between their two
traces.  The admissible tangent cone does not contain independent
variations of those levels.  The exact finite-phase Euler equation is a
divided-difference law, (4.4), and it does not imply CMC for a threshold
cut.  This is an actual contact obstruction, not a missing regularity
citation.

In particular, neither the integrated profile deficit nor the matrix
retention estimate presently controls the jump/contact part of \(M(G)\).
The desired BV-to-killed-Wulff-tube step therefore does not follow from
(0.1).  A valid completion needs one additional theorem: either exclusion
of jump contacts on the central value band, a quantitative charge of all
incompatible contacts to the profile deficit, or a killed-tube theorem for
the resulting layered varifold.  All constants below are independent of
the ambient dimension.

## 1. Setup and the constant anisotropy

Let \(E\) be the affine hull of the measure and let
\(\Omega=\operatorname {ri}(\operatorname {dom}V)\subset E\).  In this
note the first-variation calculations are local in \(\Omega\).  Assume
there that

\[
 d\mu=\rho\,dx,\qquad \rho=e^{-V}>0,qquad
 \rho\in C^1_{\rm loc}(\Omega).                    \tag{1.1}
\]

The nonsmooth convex-potential case is obtained only after a separate
approximation; no differentiability of a general convex \(V\) is being
silently asserted here.

For \(u\in BV_{\rm loc}(\Omega)\), write

\[
 Du=\sigma_u|Du|,qquad
 |Du|_\mu=\rho|Du|,qquad
 M(u)=\int \sigma_u\sigma_u^T\,d|Du|_\mu.          \tag{1.2}
\]

Put \(X_0=M(G)-M(F)\) and

\[
 {\cal S}=\partial\|X_0\|_*
 =\{H=H^T:\|H\|_{op}\le1,
               \operatorname {tr}(HX_0)=\|X_0\|_*\}.          \tag{1.3}
\]

For \(H\in{\cal S}\), define

\[
 \Phi_H(\xi)=|\xi|+
       \kappa {\xi^TH\xi\over|\xi|},\qquad \Phi_H(0)=0.       \tag{1.4}
\]

For every unit \(n\) and \(h\perp n\),

\[
 (1-3\kappa)|h|^2
 \le D^2\Phi_H(n)[h,h]
 \le(1+3\kappa)|h|^2,                              \tag{1.5}
\]

and

\[
 (1-\kappa)|\xi|\le\Phi_H(\xi)
                    \le(1+\kappa)|\xi|.            \tag{1.6}
\]

These are the only ellipticity constants used below.  They are numerical
and dimension-free.

## 2. The exact BV inner-variation theorem

Let

\[
 {\mathscr X}_\mu=
 \{X\in C_c^1(\Omega;E):\operatorname {div}(\rho X)=0\}.
                                                               \tag{2.1}
\]

The flow \(\varphi_t\) of such an \(X\) preserves \(\mu\) exactly, so
\(G_t=G\circ\varphi_t^{-1}\) is equimeasurable with \(G\) for both signs
of sufficiently small \(t\).

**Theorem 2.1 (one compatible subgradient for all BV inner
variations).**  Suppose \(G\) is a minimizer of (0.1).  There is one
\(H\in{\cal S}\) such that, for every \(X\in{\mathscr X}_\mu\),

\[
 \boxed{
 \int_\Omega
 \left[
  \Phi_H(\sigma_G)(\operatorname {div}X-\nabla V\cdot X)
  -D\Phi_H(\sigma_G)\cdot(DX)^T\sigma_G
 \right]d|DG|_\mu=0.}                              \tag{2.2}
\]

Since \(X\in{\mathscr X}_\mu\), the first parenthesis in (2.2) is zero;
it is retained to display the full weighted first-variation formula.

**Proof.**  The BV change-of-variables formula under a \(C^1\)
diffeomorphism gives, for every fixed \(H\),

\[
 \begin{split}
 \int \rho\,\Phi_H(dDG_t)
 =\int &\rho(\varphi_t(x))\det D\varphi_t(x)\\
 &\times\Phi_H((D\varphi_t(x))^{-T}\sigma_G(x))\,d|DG|(x).
 \end{split}                                        \tag{2.3}
\]

All coefficients and their first derivatives are uniformly bounded on
the compact support of \(X\).  Dominated convergence therefore
differentiates (2.3) at zero.  The three derivatives are, respectively,
\(-\nabla V\cdot X\), \(\operatorname {div}X\), and
\(-D\Phi_H(\sigma_G)\cdot(DX)^T\sigma_G\).  This proves that the left
side of (2.2), denoted \(L_H(X)\), is the first variation of the
\(\Phi_H\)-variation.  It is linear and continuous in \(X\), and affine
and continuous in \(H\).

The one-sided derivative of the nuclear norm is

\[
 {d\over dt}\bigg|_{0+}\|X_0+tA\|_*
      =\max_{H\in{\cal S}}\operatorname {tr}(HA).   \tag{2.4}
\]

Minimality of \(G\) along \(G_t\), for both signs, hence says

\[
                         \max_{H\in{\cal S}}L_H(X)\ge0
 \quad\hbox{for every }X\in{\mathscr X}_\mu.       \tag{2.5}
\]

We use the following elementary separation lemma.  If \({\cal S}\) is a
compact convex subset of a finite-dimensional space and \(L_H\) is affine
in \(H\) and linear in \(X\) on a vector space \({\mathscr X}\), then

\[
 \sup_{H\in{\cal S}}L_H(X)\ge0\quad( X\in{\mathscr X})
 \quad\Longrightarrow\quad
 \exists H_*\in{\cal S}:L_{H_*}(X)=0\quad(X\in{\mathscr X}).  \tag{2.6}
\]

Indeed, if zero were not in the compact convex set
\(\{L_H:H\in{\cal S}\}\), viewed with the topology of pointwise
convergence on \({\mathscr X}\), strict separation would give one
\(X\) for which \(L_H(X)<0\) for every \(H\).  This contradicts (2.5).
Applying (2.6) proves (2.2).  QED.

By weighted coarea, (2.2) is equivalently the **aggregate** identity

\[
              \int_{\mathbb R}\delta P_{\Phi_H,\mu}
                  (\{G>r\})[X],dr=0
       \qquad(X\in{\mathscr X}_\mu).                \tag{2.7}
\]

It is important that \(X\) is the same at every value \(r\).  Equation
(2.7) does not say that its integrand vanishes for almost every \(r\).

Coarea does, of course, imply that \(\{G>r\}\) has locally finite
perimeter for almost every \(r\), and De Giorgi's theorem makes its reduced
boundary countably \((k-1)\)-rectifiable.  This supplies an almost-everywhere
normal but no levelwise generalized-mean-curvature bound.  Allard or
elliptic parametric-integrand regularity can only be invoked after such a
bound or a local quasiminimality inequality has been established; neither
is contained in (2.7).

There is also a pressure formulation of (2.2).  Locally set \(Y=\rho X\).
The functional in (2.2) annihilates every compactly supported
divergence-free \(Y\).  The distributional de Rham lemma supplies a scalar
distribution \(q\), unique up to a constant, for which the aggregate
first variation equals \(\langle q,\operatorname {div}Y\rangle\).
This is one pressure in physical space.  It is not a proof that
\(q=\lambda(G)\), and it gives no levelwise multiplier.

## 3. What is true on a smooth diffuse band

The formal quantile calculation becomes a theorem when the BV contact
issue is excluded explicitly.

**Theorem 3.1 (exact smooth-band Euler equation).**  In addition to
(1.1), assume \(\rho\in C^2\).  Let \(J\) be an open interval and suppose
that \(G\in C^2(\{G\in J\})\), \(|\nabla G|>0\) there, and the law of
\(G\) has a positive \(C^1\) density \(w\) on \(J\).  If \(G\) is a
local minimizer of (0.1) under equimeasurable perturbations supported in
that band, then there are one \(H\in{\cal S}\) and a measurable
\(\lambda:J\to\mathbb R\) such that

\[
 -\operatorname {div}_\mu D\Phi_H(\nabla G)
                         =\lambda(G)                \tag{3.1}
\]

distributionally on \(\{G\in J\}\).  Every regular level in this band
therefore has constant weighted \(\Phi_H\)-mean curvature
\(\lambda(r)\).

**Proof.**  For \(u\in C_c^2(\{G\in J\})\), let
\(K_t=G+tu\), let \(C_t\) be its distribution function, and put

\[
                         G_t=C_0^{-1}(C_t(K_t)).     \tag{3.2}
\]

The generalized-inverse convention makes \(G_t\) exactly equimeasurable
with \(G\).  Coarea differentiation gives

\[
 \partial_t C_t(s)|_{t=0}
     =-w(s)\,\mathbb E[u\mid G=s],
 \qquad
 \dot G_0=u-\mathbb E[u\mid G].                    \tag{3.3}
\]

Thus every smooth \(u\) with \(\mathbb E[u\mid G]=0\) is a two-sided
tangent direction.  Applying the separation argument (2.6) to this
linear tangent space gives one \(H\in{\cal S}\) for which

\[
 \int D\Phi_H(\nabla G)\cdot\nabla u\,d\mu=0
 \quad\hbox{whenever }\mathbb E[u\mid G]=0.         \tag{3.4}
\]

The orthogonal complement of
\(\{u:\mathbb E[u\mid G]=0\}\) consists of the functions measurable
with respect to \(G\).  Local truncation if necessary gives (3.1).
Zero-homogeneity of \(D\Phi_H\) turns (3.1) into the asserted level
equation.  QED.

The \(H\) in Theorems 2.1 and 3.1 can be chosen to be the same.  Indeed,
compose a \(\mu\)-preserving flow with a distribution-restored curve
(3.2).  The composition is exactly equimeasurable and its tangent is the
sum of the two tangents.  Apply the separation lemma to the vector space
of all such sums.  This observation is important: no level-dependent
anisotropy is introduced on the diffuse band.

If a regular level meets a \(C^2\) boundary of the support and all
admissible boundary-tangent flows are used, the corresponding natural
contact condition is

\[
                  D\Phi_H(n_{\{G>r\}})\cdot n_\Omega=0.        \tag{3.5}
\]

For a nonsmooth convex support, (3.5) is replaced by a normal-cone
variational inequality.  There is no dimension-free \(C^2\) boundary
regularity statement for an arbitrary convex support; a Wulff tube must be
killed at first support contact.

## 4. The BV contact obstruction in exact algebraic form

The obstruction can already be seen for an \(SBV\) function with finitely
many values.  Let

\[
 G=\sum_{i=0}^m a_i\mathbf1_{A_i},qquad
                       a_0<a_1<\cdots<a_m,          \tag{4.1}
\]

where the \(A_i\) form a Caccioppoli partition.  On a smooth interface
\(\Sigma_{ij}=\partial^*A_i\cap\partial^*A_j\), coarea gives

\[
 \int\rho\Phi_H(dDG)
   =\sum_{i<j}(a_j-a_i)
            P_{\Phi_H,\mu}(\Sigma_{ij}).            \tag{4.2}
\]

Equimeasurability fixes every phase volume \(\mu(A_i)\).  Consequently,
the ordinary multiphase first-variation calculation gives phase pressures
\(p_i\) such that

\[
 (a_j-a_i)\,{\mathcal H}_{ij}=p_j-p_i
                  \quad\hbox{on }\Sigma_{ij},       \tag{4.3}
\]

where \({\mathcal H}_{ij}\) is the signed weighted anisotropic mean
curvature, with one consistent orientation convention.  Thus

\[
 \boxed{{\mathcal H}_{ij}={p_j-p_i\over a_j-a_i}.}  \tag{4.4}
\]

For a threshold \(r\), the reduced boundary of \(\{G>r\}\) contains all
interfaces \(\Sigma_{ij}\) with \(a_i<r<a_j\).  It has one constant mean
curvature only if

\[
 {p_j-p_i\over a_j-a_i}
 \quad\hbox{is the same for every active interface }(i,j).     \tag{4.5}
\]

The phase-volume Euler equations do not impose (4.5).  Equivalently, they
do not force the pressures \(p_i\) to be affine functions of the labels
\(a_i\).  Formula (4.4) is the precise jump/contact replacement for the
smooth equation (3.1).

There is a useful multiplier interpretation.  Choose any scalar function
\(\Lambda\) on the finite set \(\{a_0,\ldots,a_m\}\) with
\(\Lambda(a_i)=p_i\).  Then (4.4) reads

\[
 {\mathcal H}_{ij}
   ={\Lambda(a_j)-\Lambda(a_i)\over a_j-a_i}.       \tag{4.6}
\]

If \(\Lambda\) is interpolated absolutely continuously and
\(\lambda=\Lambda'\), this is

\[
 {\mathcal H}_{ij}
   ={1\over a_j-a_i}\int_{a_i}^{a_j}\lambda(s)\,ds. \tag{4.7}
\]

Thus a jump sees the **average** of the value multiplier over its entire
trace interval.  A diffuse level sees the point value \(\lambda(r)\).
This is the correct SBV interpretation of the formal value-multiplier
calculus.  It also shows why writing
\(-\operatorname {div}_\mu D\Phi_H(DG)=\lambda(G)\) without specifying an
Anzellotti field and its jump traces is not a well-defined BV theorem.

Here is a concrete first-variation model.  In a large planar ball with
constant density, place a disk of radius \(1\) carrying value \(2\), a
disjoint disk of radius \(2\) carrying value \(1\), and put value \(0\)
elsewhere.  The interfaces are circles and the phase pressures can be
chosen so that (4.3) holds: the high--low interface has curvature \(1\)
and tension \(2\), whereas the middle--low interface has curvature
\(1/2\) and tension \(1\).  For every \(r\in(0,1)\), the cut
\(\{G>r\}\) is the union of the two disks and its two boundary components
have different mean curvatures.  Nevertheless the partition is stationary
under all separate phase-volume-preserving smooth deformations.  It is
also locally stable against transferring an area \(s\) of middle phase
from the radius-\(2\) disk to a coating of the radius-\(1\) disk, because

\[
 {d\over ds}\bigg|_{s=0}
 \left[2\pi\sqrt{1+s/\pi}
       +2\pi\sqrt{4-s/\pi}\right]
 =1-{1\over2}>0.                                   \tag{4.8}
\]

This configuration is not asserted to be the global rearrangement
minimizer: moving through the energy barrier and coating completely lowers
the perimeter.  Its purpose is sharper: it proves that compatible-
subgradient stationarity, ordinary partial regularity, and stability under
small smooth variations do **not** imply levelwise CMC.  Global
equimeasurable minimality must be used in an additional, genuinely
nonlocal no-contact argument.

For a general \(SBV\) function the same phenomenon is encoded by its jump
traces.  On the central interval \(J=(\alpha,\beta)\), define

\[
 \ell_J(x)=
  \big[\min\{G^+(x),\beta\}-\max\{G^-(x),\alpha\}\big]_+.
                                                               \tag{4.9}
\]

Then coarea gives the exact jump identities

\[
 \begin{split}
 \int_J P_\mu(\{G>r\};J_G)\,dr
   &=\int_{J_G}\ell_J\rho\,d{\mathcal H}^{k-1},\\
 \int_J M(\{G>r\};J_G)\,dr
   &=\int_{J_G}\ell_J\,n_Gn_G^T\rho\,d{\mathcal H}^{k-1}.
 \end{split}                                        \tag{4.10}
\]

One geometric interface is therefore counted for an entire interval of
levels.  A value-dependent deformation which moves it differently at two
values in that interval destroys the subgraph ordering and is not an
admissible two-sided BV variation.

## 5. Why profile deficit and matrix retention do not remove the issue

Write the BV decomposition

\[
 DG=\nabla G\,dx+D^cG+(G^+-G^-)n_G
                  {\mathcal H}^{k-1}\!\llcorner J_G.           \tag{5.1}
\]

The matrix retention estimate controls only the sum

\[
 M(G)=M^a(G)+M^c(G)+M^j(G).                        \tag{5.2}
\]

It gives no separate lower bound on \(M^a\), the part to which Theorem
3.1 applies.  Nor does the scalar integrated deficit give such a bound.
For example, if \(A\) is an isoperimetric set and
\(G=\mathbf1_A\), then every nontrivial level equals \(A\), the integrated
profile deficit is zero, and all of \(M(G)\) is jump matrix.  Taking
\(F=G\), this \(G\) is an actual global minimizer of (0.1), since every
equimeasurable competitor \(B\) satisfies
\[
 P_\mu(B)+\kappa\|M(B)-M(A)\|_*\ge P_\mu(B)\ge P_\mu(A).
\]
Thus this is not merely a stationary example.  It also shows that an
estimate of the form

\[
                     \operatorname {tr}M^j(G)
                  \le C\Delta_G                     \tag{5.3}
\]

is false for every finite universal \(C\).  A hypothesis on the value law
is indispensable.  Positivity of the value density by itself only proves
the smooth calculation after one already knows that the relevant part of
\(DG\) is diffuse; it does not eliminate local jump traces whose omitted
values are realized elsewhere.

On the jump part, (4.10) retains the normal projector perfectly, so angular
variance can in principle be retained there.  The smooth Wulff formula
cannot be applied level by level unless (4.5), a suitable generalized
version of it, or a quantitative contact-error estimate has first been
proved.

## 6. Exact conclusion for the KLS route

The following statements are now rigorous and dimension-free.

1. The minimizer has one compatible constant nuclear subgradient for all
   common \(\mu\)-preserving BV inner variations (Theorem 2.1).
2. On every smooth noncritical diffuse band with positive value density,
   exact distribution restoration gives
   \(-\operatorname {div}_\mu D\Phi_H(\nabla G)=\lambda(G)\)
   (Theorem 3.1).
3. A regular level from such a band has the constant-anisotropy CMC and
   Young-contact data required by the killed Wulff calculation, with
   ellipticity constants \(1\pm3\kappa\).

What is not a consequence of the minimization is equally precise.  The BV
theorem gives aggregate stationarity (2.7).  On jump contacts it gives the
divided-difference law (4.4), not CMC of every threshold.  Neither
\(\Delta_G\) nor the retained total matrix separates the diffuse and
contact charges.

Thus the load-bearing statement needed next must have one of the following
formal forms.

* **No-contact theorem:** on the central band, \(D^cG=0\) and the jump
  matrix in (4.10) is zero (or negligible compared with the retained angular
  variance).
* **Deficit-charge theorem:** the part of (4.10) on which (4.5) fails is at
  most \(C\Delta_G/\kappa\), with a numerical \(C\) independent of
  dimension.
* **Layered-tube theorem:** the divided-difference interfaces themselves
  admit a killed Wulff disintegration whose collision/focal charge is
  controlled by the same profile deficit.

Without one of these additions, replacing the formal smooth calculation
by the phrase “BV regularity” is invalid.  The note refutes that proposed
first-variation/regularity inference; it does not refute the possibility
of a stronger theorem that exploits global equimeasurable minimality to
exclude all incompatible contacts.
