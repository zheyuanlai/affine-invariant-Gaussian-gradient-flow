# Global flow-cell incidence: exact ray cells and a tangent-cell no-go

## 0. Verdict

There is an exact measurable smooth analogue of the *incidence* part of
the polyhedral central-cell argument on the normal basin: stopped
signed-normal rays partition that basin, outside a null cut locus, with
reuse exactly one. Along every ray the weighted Jacobian is log-concave.
If the normal basin has full measure, then for a balanced interface the
two families of cells each carry mass \(1/2\). A free-boundary contact
stratum can have a positive-measure shadow not covered by interior normal
rays; this is an additional global incidence obstruction, recorded
explicitly below.

This partition does **not** give the missing dimension-free completion
inequality. The canonical additive continuum target is to put \(q(x)\)
equal to the smaller mass of the tangent halfspaces at \(x\), put \(s(x)\)
equal to the full tangent-slice area, and find a positive measure \(\tau\)
on the interface such that

\[
 \int s\,d\tau\le C_0P_\mu(E),
 \qquad
 \int q\,d\tau\ge c_0.                                      \tag{0.1}
\]

The one-dimensional marginal inequality would imply
\(P_\mu(E)\ge c_0/(2\sqrt3C_0)\).

No universal \((C_0,c_0)\) statement of the form (0.1) follows from
balance, CMC, high normal rank, stability, or excellent normal tubes. For
the isotropic radial exponential law and its median sphere,

\[
                         {q(x)\over s(x)}
                  \le {2\over\sqrt{n+1}}\quad\text{for every }x. \tag{0.2}
\]

Thus the first inequality in (0.1) forces the second integral to be
\(O(C_0/\sqrt n)\). The sphere is balanced, smooth weighted CMC,
volume-constrained stable, has normalized normal matrix \(I/n\), and has
no killed rays on a fixed short tube. A chord followed by volume repair has
positive constrained second variation. Consequently

\[
 \text{failure of additive tangent completion}
 \ \Longrightarrow\ \text{macroscopic bevel saving}          \tag{0.3}
\]

is false if it uses only those local/tube/stability data.

The sphere is not asserted to be an exact global Cheeger minimizer. This
qualification is essential: exact global minimality is the only remaining
datum that could force a different, non-additive assignment. The report
also audits the centered-normal stability mechanism. It gives the sharp
estimate

\[
 \int_{\partial\Sigma}\mathrm{II}_{\partial K}(N,N)
 \le {2|m|\over(1-|m|)^2}\int_\Sigma|A|^2,\qquad
 |m|^2\le\|Q_N\|_{\rm op},                                  \tag{0.4}
\]

and, when combined with the tensor Minkowski identity, forces a large
reciprocal contact-curvature moment. It still does not create long support
cells: the smooth convex model
\(x_1^4+x_2^2+\cdots+x_n^2\le1\) has zero contact curvature in direction
\(e_1\) but no tangent segment in that direction. The report therefore
proves the exact ray partition and decisively rules out the canonical
additive tangent-cell and curvature-to-ray mechanisms, but does not prove
the global calibrated-cell theorem needed for KLS.

## 1. Exact stopped-normal flow cells

Let \(\Omega\subset\mathbb R^n\) be open and convex and let
\(d\mu=\rho dx\), where \(\rho=e^{-V}/Z\) is positive and smooth in
\(\Omega\), with \(V\) convex. Let \(E\subset\Omega\) have compact
\(C^2\) relative boundary

\[
                          \Sigma=\partial E\cap\Omega,        \tag{1.1}
\]

and let \(\nu(x)\) be its exterior unit normal. The hard-support boundary
is free and is not part of \(\Sigma\).

For \(x\in\Sigma\), let \(\tau_+(x)\) be the supremum of the \(t>0\) for
which

\[
                           F(x,t)=x+t\nu(x)                   \tag{1.2}
\]

lies in \(\Omega\), has \(x\) as its unique nearest point of \(\Sigma\),
and has no focal point on \([0,t]\). Define \(\tau_-(x)\) analogously for
\(x-t\nu(x)\). These stopping times are Borel: their defining conditions
can be tested on rational subintervals, and unique metric projection onto
a closed set is a Borel condition.

Put

\[
 d\sigma_\mu(x)=\rho(x)d\mathcal H^{n-1}(x),                 \tag{1.3}
\]

and, at a regular base point with principal curvatures
\(\kappa_1(x),\ldots,\kappa_{n-1}(x)\), set

\[
 j_x(t)=\prod_{r=1}^{n-1}(1+t\kappa_r(x))
       \exp[-V(x+t\nu(x))+V(x)].                             \tag{1.4}
\]

The factors in (1.4) are positive before the stopping time.

Let \(\mathcal B_+\) and \(\mathcal B_-\) be the points in the two phases
whose unique closest point in \(\overline\Sigma\) lies in the relative
interior \(\Sigma\), rather than in the contact stratum
\(\overline\Sigma\cap\partial\Omega\). Call these the two normal basins.

### Proposition 1.1 (reuse-one normal-basin disintegration)

Assume the cut locus inside the normal basins has \(\mu\)-measure zero.
Then for every nonnegative Borel function \(g\),

\[
\begin{aligned}
 \int_{\mathcal B_+}g\,d\mu
   &=\int_\Sigma\int_0^{\tau_+(x)}
       g(x+t\nu(x))j_x(t)\,dt\,d\sigma_\mu(x),\\
 \int_{\mathcal B_-}g\,d\mu
   &=\int_\Sigma\int_0^{\tau_-(x)}
       g(x-t\nu(x))j_x(-t)\,dt\,d\sigma_\mu(x).               \tag{1.5}
\end{aligned}
\]

Every point of the normal basins off the cut locus occurs exactly once.
If the contact shadows
\((E^c\cap\Omega)\setminus\mathcal B_+\) and
\((E\cap\Omega)\setminus\mathcal B_-\) are null, (1.5) holds with the
entire phases on the left.

#### Proof

At a point with unique closest boundary point \(x\), the minimizing
segment is perpendicular to \(T_x\Sigma\), hence lies on (1.2).
Conversely, every point strictly before a cut, focal, or support time has
unique projection \(x\). The normal maps are therefore injective and
cover the normal basins outside the cut locus. Their Euclidean Jacobians are
\(\prod_r(1\pm t\kappa_r)\). Multiplying by the density ratio gives
(1.4), and the area formula gives (1.5). QED.

For closed \(C^2\) hypersurfaces the contact shadows are empty and the cut
locus has Lebesgue measure zero. Several components are covered by
assigning a point to its unique nearest component; ties belong to the cut
locus. For a hypersurface ending on hard support, a positive-measure set
may have its closest point in the contact stratum. Such points require
separate support/contact cells and are not silently included in (1.5).

Define the one-sided cell masses

\[
 \ell_+(x)=\int_0^{\tau_+(x)}j_x(t)\,dt,
 \qquad
 \ell_-(x)=\int_0^{\tau_-(x)}j_x(-t)\,dt.                    \tag{1.6}
\]

For a balanced interface with full normal-basin coverage, (1.5) gives

\[
       \int_\Sigma\ell_+\,d\sigma_\mu
       =\int_\Sigma\ell_-\,d\sigma_\mu={1\over2}.             \tag{1.7}
\]

There is no overlap constant hidden in (1.7): reuse is one.

### Proposition 1.2 (each stopped ray is log-concave)

For every base point before its stopping times,
\(t\mapsto j_x(t)\) is log-concave. More explicitly,

\[
 {d^2\over dt^2}\log j_x(t)
 =-\sum_{r=1}^{n-1}{\kappa_r(x)^2\over(1+t\kappa_r(x))^2}
  -\nabla^2V(x+t\nu(x))[\nu(x),\nu(x)]\le0.                 \tag{1.8}
\]

This follows by differentiating (1.4). Restriction to the stopped interval
preserves log-concavity.

Propositions 1.1--1.2 are the strongest literal smooth replacement of a
finite polyhedral cell decomposition: measurable cells, log-concave
one-dimensional densities, exact mass accounting on the normal basin, and
optimal reuse. Contact shadows require an additional stratified
decomposition and are themselves part of the global support-incidence
problem.

## 2. Why reuse-one rays do not close the inverse

Put \(p=P_\mu(E)=\sigma_\mu(\Sigma)\). Equations (1.6)--(1.7) give

\[
 \int_\Sigma(\ell_++\ell_-)\,d\sigma_\mu=1,
 \qquad
 {1\over p}\int_\Sigma(\ell_++\ell_-)\,d\sigma_\mu={1\over p}. \tag{2.1}
\]

This section assumes full normal-basin coverage. If contact shadows carry
positive mass, (2.1) accounts only for the normal-basin portion and the
unassigned contact mass is an extra term, not a favorable estimate.

Small perimeter says exactly that the average flow-cell mass per unit
interface area is large. Log-concavity of the individual \(j_x\) gives no
upper bound on this scale.

The same fact appears through the one-Lipschitz signed distance:

\[
 \int|d_\Sigma|\,d\mu
 =\int_\Sigma\left[
   \int_0^{\tau_+}t j_x(t)\,dt+
   \int_0^{\tau_-}t j_x(-t)\,dt\right]d\sigma_\mu(x).         \tag{2.2}
\]

A fixed two-sided swept mass out to \(T=h/p\) makes (2.2) at least
\(c h^2/p\). A universal upper bound for (2.2) is the first-moment KLS
formulation, not an independent incidence estimate.

One-dimensional Cheeger on each normalized ray also retains its physical
scale: its density at zero is \((\ell_++\ell_-)^{-1}\), while its variance
may be large. Isotropy controls only a matrix-valued aggregate of centered
ray variances in the directions \(\nu(x)\); in the dispersed high-rank
regime this is the already-known operator-norm obstruction. No covariance
trace bound is claimed as a closure here.

## 3. The canonical additive tangent-cell target

For \(x\in\Sigma\), let

\[
 \Pi_x=\{y:(y-x)\cdot\nu(x)=0\},
 \qquad H_x^\pm=\{y:\pm(y-x)\cdot\nu(x)\ge0\}.                \tag{3.1}
\]

Define the Borel functions

\[
 q(x)=\min(\mu(H_x^+),\mu(H_x^-)),
 \qquad
 s(x)=\int_{\Pi_x}\rho(y)d\mathcal H^{n-1}(y).              \tag{3.2}
\]

Every one-dimensional marginal of an isotropic log-concave law is
isotropic and log-concave, so

\[
                              s(x)\ge {q(x)\over2\sqrt3}.     \tag{3.3}
\]

### Definition 3.1 (bounded-reuse additive tangent completion)

A finite positive Borel measure \(\tau\) on \(\Sigma\) is a
\((C_0,c_0)\)-completion measure if

\[
 \boxed{\int_\Sigma s\,d\tau\le C_0p,\qquad
        \int_\Sigma q\,d\tau\ge c_0.}                        \tag{3.4}
\]

The first inequality charges completed slices with bounded reuse. The
second is the continuum central-cell incidence bound.

### Lemma 3.2 (a completion measure closes the branch)

If \(C_0<\infty\) and \(c_0>0\) are universal, then

\[
                              p\ge {c_0\over2\sqrt3C_0}.      \tag{3.5}
\]

Indeed, integrate (3.3) against \(\tau\):

\[
 C_0p\ge\int s\,d\tau
 \ge {1\over2\sqrt3}\int q\,d\tau
 \ge {c_0\over2\sqrt3}.                                     \tag{3.6}
\]

For a polyhedral interface, \(\tau\) is counting measure on its distinct
facet planes. Completion gives the first inequality, and the central-cell
union bound gives \(c_0=1/2\).

## 4. Radial-exponential no-go

Let

\[
 d\mu_n(x)=Z_n^{-1}e^{-c_n|x|}dx,\qquad c_n=\sqrt{n+1},      \tag{4.1}
\]

and let \(E_n=B(0,r_n)\), where \(r_n=s_n/c_n\) and \(s_n\) is the median
of \(\operatorname{Gamma}(n,1)\). Then \(\mu_n\) is isotropic,
\(\mu_n(E_n)=1/2\), and

\[
                              n-{1\over3}\le s_n\le n.        \tag{4.2}
\]

As proved in "macroscopic_bevel_completion.md", the median sphere is
smooth weighted CMC, volume-constrained stable, has normalized normal
matrix \(I/n\), and its whole boundary survives fixed two-sided tubes with
defect \(O(T^2)\). Its perimeter \(p_n\) lies in \([1/10,1]\).

All tangent planes are rotationally equivalent. Let \(f_n\) be the
density of \(X_1\), and put

\[
 q_n=\mu_n\{X_1\ge r_n\},\qquad s_n^{\rm tan}=f_n(r_n).       \tag{4.3}
\]

### Lemma 4.1 (tangent-cap hazard is at least \(c_n/2\))

For every \(n\ge2\),

\[
                    {q_n\over s_n^{\rm tan}}
                    \le {2\over c_n}\le {2\over\sqrt n}.     \tag{4.4}
\]

#### Proof

Fix \(t>0\), write a point of the slice \(X_1=t\) as \((t,Y)\), and put
\(R=(t^2+|Y|^2)^{1/2}\). Differentiation under the integral gives

\[
 -{d\over dt}\log f_n(t)=c_n\,\mathbb E_t{t\over R},         \tag{4.5}
\]

where \(\mathbb E_t\) is expectation for the normalized slice density
proportional to \(e^{-c_nR}dY\). Integration by parts in the
\(d=n-1\) variables \(Y\) yields

\[
                         d=c_n\mathbb E_t{|Y|^2\over R}.      \tag{4.6}
\]

Consequently

\[
 \mathbb E_tR=t^2\mathbb E_t{1\over R}+{d\over c_n}
 \le t+{d\over c_n}.                                        \tag{4.7}
\]

Jensen's inequality and (4.7) give

\[
 \mathbb E_t{t\over R}\ge {t\over\mathbb E_tR}
 \ge {t\over t+d/c_n}.                                      \tag{4.8}
\]

At \(t=r_n=s_n/c_n\), (4.2) makes the last expression at least \(1/2\).
The marginal \(f_n\) is log-concave. Hence for \(u\ge r_n\),

\[
 f_n(u)\le f_n(r_n)e^{-(c_n/2)(u-r_n)}.                      \tag{4.9}
\]

Integrating proves (4.4). QED.

### Proposition 4.2 (no universal additive tangent completion)

For the median sphere, every finite positive \(\tau\) satisfying

\[
                         \int s\,d\tau\le C_0p_n             \tag{4.10}
\]

also satisfies

\[
                         \int q\,d\tau
                  \le {2C_0p_n\over\sqrt n}
                  \le {2C_0\over\sqrt n}.                    \tag{4.11}
\]

At every boundary point, Lemma 4.1 gives
\(q(x)\le2s(x)/\sqrt n\). Integrate and use \(p_n\le1\).

Thus no fixed \(c_0>0\) in Definition 3.1 is possible with universal
\(C_0\) on balanced, stable, high-rank CMC interfaces with excellent
tubes. The continuum of tiny tangent caps covers the exterior half of the
measure, but an additive weighting that sees constant mass costs at least
\(c\sqrt n\,p_n\) in completed tangent-slice area.

The sphere also falsifies the local form of (0.3): after volume repair its
chord replacement has positive leading constrained cost and its full
second variation is nonnegative. Exact global Cheeger minimality is not
among the data used here; any theorem for exact minimizers must encode it
explicitly rather than infer it from stationarity.

### 4.3 A globally smooth radial countermodel

The nonsmoothness of \(|x|\) at the origin is inessential. Choose a
nondecreasing \(C^\infty\) function
\(\eta:[0,\infty)\to[0,1]\) such that

\[
 \eta(r)=r\quad(0\le r\le1/4),\qquad
 \eta(r)=1\quad(r\ge1),
 \qquad \phi(r)=\int_0^r\eta(s)\,ds.                        \tag{4.12}
\]

Then \(x\mapsto\phi(|x|)\) is \(C^\infty\) and convex: it equals
\(|x|^2/2\) near the origin, while \(\phi(r)=r-a\) for \(r\ge1\), for a
fixed \(a\in(0,1)\). Let

\[
                    d\widetilde\mu_n(x)
       =\widetilde Z_n^{-1}e^{-\sqrt{n+1}\phi(|x|)}\,dx.     \tag{4.13}
\]

Outside the unit ball its radial density is, up to one normalizing
constant, exactly the gamma density from (4.1). The mass of the unit ball
is \(e^{-\Omega(n\log n)}\): indeed, its unnormalized mass is at most
\(|S^{n-1}|/n\), whereas the annulus
\([n/(2c_n),n/c_n]\), \(c_n=\sqrt{n+1}\), has unnormalized radial mass at
least

\[
 |S^{n-1}|\,{n\over2c_n}
 \left({n\over2c_n}\right)^{n-1}e^{-n+c_na}.                \tag{4.14}
\]

Consequently the radial median and second moment differ from those of
(4.1) by \(e^{-\Omega(n\log n)}\) relative errors. After the scalar
dilation making (4.13) isotropic, its median sphere is still entirely in
the linear region, as is every fixed-width two-sided tube around it.
On every tangent slice at that sphere the conditional density is exactly
proportional to \(e^{-\widetilde c_nR}\). Thus the integration-by-parts
proof of Lemma 4.1 applies verbatim, with fixed slack, and gives

\[
                    {q_n\over s_n^{\rm tan}}\le {3\over\sqrt n}
                                                                    \tag{4.15}
\]

for all sufficiently large \(n\); the finitely many smaller dimensions
are irrelevant to the no-go sequence. At the median sphere the radial
potential is affine, so the weighted CMC equation, the stability
quadratic form, the normal matrix \(I/n\), and the fixed-tube defect are
the same local formulas as for (4.1), after scaling. Its perimeter remains
bounded above and below by universal constants. Hence Proposition 4.2
has a positive \(C^\infty\) log-concave-density countermodel; no
regularization loophole remains.

## 5. Support-function cells and the continuum union problem

If \(E\) is convex, orient tangent planes so \(E\subset H_x^+\). Then

\[
                     E=\bigcap_{x\in\Sigma}H_x^+,\qquad
                     E^c=\bigcup_{x\in\Sigma}H_x^-.          \tag{5.1}
\]

For finitely many facets, the union bound applied to (5.1) is the
polyhedral central-cell argument. For a continuum, (5.1) has no canonical
counting measure of bounded total reuse. Proposition 4.2 quantifies the
failure on the sphere.

Normal-ray cells avoid overlap: on the sphere they are radial cones and
partition both phases with reuse one. But they are not halfspaces, and
their densities contain the polar Jacobian. Their one-dimensional scale
is universal here because the radial gamma law has universal median
hazard. For a general interface, proving a universal flow-cell hazard is
another form of the desired inverse.

## 6. Four stress tests

### 6.1 Cube flat cut

For the isotropic cube \([-\sqrt3,\sqrt3]^n\) and the central coordinate
cut \(E=\{x_1\le0\}\), there is one tangent plane. Its smaller halfspace
mass is \(q=1/2\), its full slice is the physical boundary \(s=p\), and
one atom gives

\[
                              C_0=1,\qquad c_0={1\over2}.      \tag{6.1}
\]

The stopped rays are the parallel coordinate chords of the cube, so the
tangent-cell and flow-cell descriptions coincide.

### 6.2 Radial exponential median sphere

The stopped rays are radial and partition the two phases. By symmetry,

\[
                    \ell_+(x)=\ell_-(x)={1\over2p_n}.         \tag{6.2}
\]

Flow-cell incidence is perfect, while tangent-halfspace incidence loses
\(\sqrt n\). A smooth theory must permit curved divergent cells and place
this example directly in a universal-perimeter branch.

### 6.3 Simplex-like caps

For a convex body \(K\) and a balanced cap \(E=K\cap H\), the interior
boundary is a complete slice, and one atom gives (6.1), regardless of
asymmetry. For a finite union of vertex caps, take one atom per cut plane.
The larger halfspaces have one convex central cell; balance and the finite
union bound give \(\sum_iq_i\ge1/2\). Asymmetry alone is no obstruction.

### 6.4 Product-exponential maximum interface

At the median maximum level for independent unit exponentials, put
\(q=e^{-L}\), \(d=1-q\), so \(d^m=1/2\). There are \(m\) coordinate
facet planes. Their full slice areas and smaller halfspace masses are
\(s_i=q_i=q\), while physical facet areas are \(a_i=qd^{m-1}\). Thus

\[
 \sum_iq_i=mq\in[1/2,\log2],\qquad
 \sum_is_i=mq=2d\,p\le2p.                                   \tag{6.3}
\]

Definition 3.1 holds with \(C_0=2,c_0=1/2\). The facets meet at genuine
right-angle ridges, and the explicit bevel in
"cheeger_facet_completion.md" shows the maximum set is not exact.
Smoothing only the ridges does not invalidate the macroscopic facet
grouping; replacing it by infinitesimal tangent planes would create an
artificial mesh loss. The cyclic constrained example has the same
conclusion up to \(O(m^{-11})\).

## 7. Centered-normal stability and the contact-tensor target

There is one useful exact consequence of the centered translation tests,
but it stops one derivative short of a support-cell theorem. The precise
statement is as follows.

Let \(K\) be a bounded \(C^2\) convex body, let \(\mu\) be normalized
Lebesgue measure on \(K\), and let \(\Sigma\) be a smooth two-sided
volume-constrained stable CMC interface meeting \(\partial K\)
orthogonally. Write \(\Gamma=\partial\Sigma\), let \(N\) be the unit normal
to \(\Sigma\), and normalize the area measures by \(|K|\):

\[
 d\sigma=|K|^{-1}d\mathcal H^{n-1}|_\Sigma,\qquad
 d\eta=|K|^{-1}d\mathcal H^{n-2}|_\Gamma,\qquad
 p=\sigma(\Sigma).                                         \tag{7.1}
\]

Put

\[
 m={1\over p}\int_\Sigma N\,d\sigma,\qquad
 Q_N={1\over p}\int_\Sigma N\otimes N\,d\sigma,\qquad
 \delta=|m|.                                                \tag{7.2}
\]

Here \(\operatorname{tr}Q_N=1\), so its effective rank is
\(\|Q_N\|_{\rm op}^{-1}\).

### Proposition 7.1 (the exact centered-normal estimate)

With the convention that the second fundamental form of the convex
support is nonnegative,

\[
 \int_\Sigma |A|^2\bigl(1-|N-m|^2\bigr)\,d\sigma
 \ \ge\
 \int_\Gamma \mathrm{II}_{\partial K}(N,N)|N-m|^2\,d\eta.  \tag{7.3}
\]

Consequently

\[
 \boxed{\;
 \int_\Gamma\mathrm{II}_{\partial K}(N,N)\,d\eta
 \le {2\delta\over(1-\delta)^2}
       \int_\Sigma|A|^2\,d\sigma ,
 \;}
                                                                    \tag{7.4}
\]

provided \(\delta<1\), and

\[
                         \delta^2\le\|Q_N\|_{\rm op}.        \tag{7.5}
\]

In particular, effective rank greater than \(235\) gives
\(\delta<1/\sqrt{235}<0.0653\), and the coefficient in (7.4) is less
than \(0.150\).

#### Proof

Free-boundary stability says, for every smooth \(u\) with
\(\int_\Sigma u\,d\sigma=0\),

\[
 \mathcal Q(u)=
 \int_\Sigma\bigl(|\nabla_\Sigma u|^2-|A|^2u^2\bigr)d\sigma
 -\int_\Gamma\mathrm{II}_{\partial K}(N,N)u^2d\eta\ge0.    \tag{7.6}
\]

Use \(u_i=N_i-m_i\) and sum over \(i\). Since
\(\sum_i|\nabla_\Sigma N_i|^2=|A|^2\), this is exactly (7.3).
Moreover,

\[
 1-|N-m|^2=2N\cdot m-\delta^2\le2\delta,\qquad
 |N-m|^2\ge(1-\delta)^2.                                   \tag{7.7}
\]

The support curvature is nonnegative, so (7.3)--(7.7) imply (7.4).
For \(v=m/\delta\), Jensen and (7.2) give

\[
 \delta^2=\left({1\over p}\int_\Sigma v\cdot N\,d\sigma\right)^2
 \le {1\over p}\int_\Sigma(v\cdot N)^2d\sigma
 \le\|Q_N\|_{\rm op},                                      \tag{7.8}
\]

with the case \(\delta=0\) immediate. QED.

For a smooth log-concave density, the stability form contains the
additional nonpositive potential term
\(-\int_\Sigma\nabla^2V(N,N)u^2d\sigma_\mu\). Convexity of \(V\) only
strengthens (7.3), provided all weighted first- and second-variation
hypotheses and the free-boundary term are available. The uniform case
already exposes the geometric obstruction.

### 7.2 What a genuinely global long tube would add

Suppose every normal ray survives to both \(\pm T\), and define its
two-sided geometric logarithmic defect by

\[
 \mathcal D_T(x)
   =-\sum_{j=1}^{n-1}\log(1-T^2\kappa_j(x)^2)
   \ge T^2|A(x)|^2.                                        \tag{7.9}
\]

If

\[
                       {1\over p}\int_\Sigma\mathcal D_T\,d\sigma
                       \le2\varepsilon,                    \tag{7.10}
\]

then (7.4) gives

\[
 \int_\Gamma\mathrm{II}_{\partial K}(N,N)\,d\eta
 \le {4\delta\varepsilon\over(1-\delta)^2T^2}\,p.           \tag{7.11}
\]

At the Cheeger scale \(T=\gamma/\psi\), the right side is

\[
              {4\delta\varepsilon\over
                (1-\delta)^2\gamma^2}\,\psi^2p.             \tag{7.12}
\]

This is dimension-free and is the strongest direct output of the proposed
centered-normal mechanism: total support curvature in the interface-normal
direction is subordinate to the global tube defect.

The word *global* in (7.10) is load-bearing. If a retained packet
\(G\subset\Sigma\) has relative area \(1-\alpha\), then

\[
 \|Q_N\|_{\rm op}\le
 (1-\alpha)\|Q_{N,G}\|_{\rm op}+\alpha,                     \tag{7.13}
\]

so normal-rank information transfers when \(\alpha\) is small. Curvature
energy does not: \(|A|^2\) can concentrate on \(\Sigma\setminus G\).
Applying stability to cutoff versions of \(N_i-m_i\) introduces
\(\int|\nabla\chi|^2\), for which the retained-ray argument supplies no
bound. Thus a long-ray theorem on a packet cannot be substituted into
(7.10).

### 7.3 Combining stability with the tensor Minkowski identity

The exact contact tensor gives a sharper target than scalar contact area.
For \(c\in\mathbb R^n\), put

\[
 X={1\over p}\int_\Sigma(x-c)\otimes N\,d\sigma,\qquad
 B={1\over p}\int_\Gamma(x-c)\otimes n_K\,d\eta.            \tag{7.14}
\]

Surface integration by parts gives

\[
                              I-Q_N=H X+B.                  \tag{7.15}
\]

Since
\(\|I-Q_N\|_F^2=n-2+\operatorname{tr}(Q_N^2)\), an
isotropic-scale bound

\[
             {1\over p}\int_\Sigma|x-c|^2d\sigma\le Cn     \tag{7.16}
\]

together with \(|H|\le2p\) implies the following dichotomy: either
\(p\ge c(C)>0\), or

\[
                              \|B\|_F\ge {1\over4}\sqrt n.   \tag{7.17}
\]

Thus a hypothetical small-perimeter interface must carry a
dimension-sized support-contact tensor.

Let \(k=\mathrm{II}_{\partial K}(N,N)\). Hilbert-space
Cauchy--Schwarz, with the convention \(1/0=+\infty\), yields the exact
reciprocal-curvature consequence

\[
 p^2\|B\|_F^2
 =\left\|\int_\Gamma(x-c)\otimes n_K\,d\eta\right\|_F^2
 \le\left(\int_\Gamma k\,d\eta\right)
       \left(\int_\Gamma{|x-c|^2\over k}\,d\eta\right).      \tag{7.18}
\]

Combining (7.11) and (7.17) gives

\[
 \int_\Gamma{|x-c|^2\over k}\,d\eta
 \ge {p\,n(1-\delta)^2T^2\over64\delta\varepsilon}.         \tag{7.19}
\]

For \(T=\gamma/\psi\) and a balanced Cheeger interface
\(\psi=2p\), this becomes

\[
 \int_\Gamma{|x-c|^2\over k}\,d\eta
 \ge {n(1-\delta)^2\gamma^2\over256\delta\varepsilon\,p}.   \tag{7.20}
\]

This is a rigorous contact-tensor/long-ray dichotomy. It says that the
large tensor flux is carried either where \(k=0\), or at contact points
having a very large reciprocal normal-curvature moment. It does **not**
yet produce support-to-support cells.

### 7.4 Why contact curvature does not create support rays

The missing implication is not a technical regularity detail. A
second-order vanishing statement on the contact stratum does not imply a
ruled support in the transverse direction. The bounded \(C^\infty\)
convex body

\[
             K=\left\{x\in\mathbb R^n:
                  x_1^4+x_2^2+\cdots+x_n^2\le1\right\}      \tag{7.21}
\]

has contact stratum

\[
 \Gamma=\{x_1=0,\ x_2^2+\cdots+x_n^2=1\}.                  \tag{7.22}
\]

At every point of \(\Gamma\), the tangent direction \(N=e_1\) satisfies

\[
                    \mathrm{II}_{\partial K}(N,N)=0.        \tag{7.23}
\]

Nevertheless \(x+tN\notin K\) for every \(t\ne0\), because its defining
function equals \(1+t^4\). There is no boundary ruling and not even a
nontrivial tangent support segment. This example remains so after an
affine isotropic normalization.

For the flat interface \(\Sigma=K\cap\{x_1=0\}\), completion still holds,
but for a different global reason: \(\Sigma\) itself is a full affine
slice. More generally, if \(|A|\equiv0\) on a connected compact interface
whose entire boundary lies on \(\partial K\), then \(\Sigma\) lies in an
affine hyperplane and is both relatively open and relatively closed in the
connected convex slice \(K\cap\Pi\); hence it is the full slice. The
finite collection of such components is handled by the polyhedral
central-cell argument.

Thus the exact zero-defect branch closes, but (7.23) alone does not.
Equations (7.18)--(7.20) control only the two-jet of the support at
\(\Gamma\). A contact-cell theorem needs an additional *longitudinal*
statement controlling how support normals turn after leaving
\(\Gamma\) in direction \(N\), or it must use global minimality to bypass
that motion altogether. Isotropic convex bodies have no universal lower
support-curvature bound; smoothed cubes make such a bound impossible.

There are three further audit points.

1. If \(m=0\), (7.3) forces \(k=0\) almost everywhere on \(\Gamma\), but
   gives no information on \(|A|\). With no hard support, the radial
   exponential sphere has \(m=0\), nonzero \(A\), and equality after
   summing its weighted stability tests.
2. A disconnected interface can have high global normal rank merely by
   mixing rank-one flat complete slices. This is the completion branch,
   not a curvature contradiction.
3. For anisotropic surface energy, the translation test functions are
   components of the Cahn--Hoffman field and the interior energy contains
   \(D^2\Phi\); formulas (7.3)--(7.9) do not transfer by replacing \(N\)
   with an anisotropic normal. Any anisotropic bevel approximation must
   redo this calculation.

The new load-bearing statement is therefore precise: convert the
dimension-sized tensor \(B\), together with either a longitudinal support
turning bound or exact global minimality, into complete macroscopic cells
whose smaller-side masses add to a universal constant with bounded reuse.
Neither stability nor the pointwise contact curvature supplies that
conversion.

## 8. Calibration flow lines

The direct Cheeger-deficit functional may give a field \(z\) with

\[
 |z|\le1,\qquad z\cdot\nu_E=1\quad|D\mathbf1_E|\text{-a.e.},
 \qquad |\operatorname{div}_\mu z|\le\psi.                  \tag{8.1}
\]

If \(z\) were smooth and generated a unique flow \(X_t\), its weighted
Jacobian would satisfy

\[
 {d\over dt}\log J_\mu(t,x)=\operatorname{div}_\mu z(X_t(x)),
 \qquad e^{-\psi|t|}\le J_\mu(t,x)\le e^{\psi|t|}.           \tag{8.2}
\]

This is adapted to the \(1/\psi\) scale and again gives reuse-one
longitudinal cells. But the BV calibration is merely bounded with
divergence in \(L^\infty\); a unique regular Lagrangian flow requires
additional Sobolev or BV control. Even for a smooth flow, reassigning a
collection of cells creates lateral perimeter not controlled by (8.2).
That lateral cost is the smooth completion defect.

## 9. Exact remaining target

A viable statement must distinguish the sphere from flat examples and use
global minimality before a local or mesh limit:

> **Calibrated flow-cell dichotomy.** There are universal \(c_0,C_0>0\)
> such that for every isotropic smooth log-concave law and every balanced
> exact Cheeger minimizer, either \(\psi_\mu\ge c_0\), or there is a
> countable family of macroscopic flow cells and full marginal slices
> satisfying (3.4), or a union/intersection reassignment decreases
> perimeter by \(c_0P_\mu(E)\) while changing volume by at most
> \(C_0P_\mu(E)^2\).

In the uniform convex-body branch, (7.15) makes the required theorem more
specific:

> **Contact-tensor flow-cell inverse.** Suppose the position moment
> satisfies (7.16), \(p<c_0\), and hence
> \(\|B\|_F\ge\sqrt n/4\). Decompose the contact flux
> \((x-c)\otimes n_K\,d\eta\), using exact global minimality, into a
> countable family of support-to-support cells with overlap at most
> \(C_0\). The family must either charge complete marginal slices by at
> most \(C_0p\) and have total smaller-side mass at least \(c_0\), or
> generate an admissible volume-repaired reassignment with strictly
> smaller perimeter.

Equations (7.18)--(7.20) identify where such cells have to start, but the
quartic-support example shows that their longitudinal extent cannot be
deduced from the contact second fundamental form. A proof must extract it
from convex support geometry over a nonzero distance, or directly from
global minimality.

The cube and simplex cap lie in the completion branch. The radial sphere
must lie in the first branch. The product maximum lies in the completion
branch and also exhibits the ridge competitor because it is not exact.

No proof of this dichotomy is obtained here. Proposition 4.2 proves that
one cannot replace macroscopic flow cells by an additive measure on
infinitesimal tangent halfspaces, and CMC cancellation proves that failure
of that additive measure does not itself yield a chord saving. The
unresolved datum is the lateral, globally calibrated incidence of flow
cells.
