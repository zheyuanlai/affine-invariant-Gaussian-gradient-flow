# Macroscopic smooth bevels: curvature extraction and a stable high-rank countermodel

## 0. Verdict

The polyhedral bevel mechanism does not extend to a smooth
constant-mean-curvature interface by replacing a curved patch with a
chord and repairing the volume.  The reason is an exact first-order
cancellation:

\[
                 \delta P=\lambda\,\delta V                 \tag{0.1}
\]

on a weighted CMC leaf.  A remote repair of the volume cancels the apparent
first-order chord saving, leaving precisely the constrained second
variation.  Spheres make this obstruction explicit and sharp.

There is nevertheless a useful positive consequence of a two-sided
killed-normal packet.  In the Euclidean/log-affine model, if both normal
rays survive to distance \(T\) with convexity defects at most \(\eta\),
then pointwise

\[
                           |II|^2\le {2\eta\over T^2}.         \tag{0.2}
\]

Thus the hypothetical \(T\asymp1/\psi\) packet is genuinely nearly flat
when \(\psi\) is small.  What (0.2) does not provide is tangential
clearance, slice completion, or bounded-reuse incidence.

To show that this is not a technical gap, Section 5 gives a rigorous
smooth countermodel satisfying all the **local** hypotheses simultaneously:
an isotropic radial exponential law and its median sphere.  The sphere is
balanced, weighted CMC with \(|\lambda|=O(n^{-1/2})\), volume-constrained
stable, has normalized normal matrix \(I/n\), and its entire boundary
supports two-sided tubes of fixed length with defect \(O(T^2)\).  Its
Cheeger scale is universal.  The model does not assert that the sphere is
a global Cheeger minimizer; rather, it proves that stationarity, stability,
high normal rank, and excellent tubes cannot replace exact global
minimality.  A successful smooth completion theorem must use global
support/contact/collision incidence or the common calibration.

## 1. First-order cancellation for every chord-and-repair scheme

Let \(d\mu=e^{-V}dx/Z\), and let \(\Sigma=\partial E\) be a smooth
weighted CMC hypersurface, away from hard support.  Fix the sign convention
in which an outward normal speed \(f\) satisfies

\[
 V'(0)=\int_\Sigma f\,d\sigma_\mu,
 \qquad
 P'(0)=\int_\Sigma H_\mu f\,d\sigma_\mu,                     \tag{1.1}
\]

and suppose \(H_\mu\equiv\lambda\).  Then

\[
                             P'(0)=\lambda V'(0).             \tag{1.2}
\]

Consider any local modification with small changes
\((\Delta P_{\rm loc},\Delta V_{\rm loc})\).  Repair its mass on a remote
regular patch.  By (1.2), a repair of volume
\(-\Delta V_{\rm loc}\) costs

\[
                         -\lambda\Delta V_{\rm loc}
                              +o(|\Delta V_{\rm loc}|).        \tag{1.3}
\]

Hence the leading fixed-volume cost is not the raw chord saving
\(\Delta P_{\rm loc}\), but

\[
                    \boxed{\Delta P_{\rm loc}
                              -\lambda\Delta V_{\rm loc}.}    \tag{1.4}
\]

For a polyhedral corner, \(\Delta P_{\rm loc}=-c\varepsilon+O(\varepsilon^2)\)
while \(\Delta V_{\rm loc}=O(\varepsilon^2)\), so (1.4) is negative.
For a smooth patch of curvature \(k\), both the arc--chord saving and
\(\lambda\) times the enclosed cap volume are of order
\(k^2\ell^{n+1}\).  They must be compared, and stability says that the
comparison need not be favorable.

### 1.1 Exact circle calculation

Take a Euclidean circle of radius \(R\), with locally constant density,
and replace the arc with central angles \([-\theta,\theta]\) by its chord.
The changes in perimeter and enclosed area are

\[
 \Delta P=2R(\sin\theta-\theta),
 \qquad
 \Delta V=-R^2(\theta-\sin\theta\cos\theta).                 \tag{1.5}
\]

The CMC multiplier is \(\lambda=1/R\).  Therefore

\[
\begin{aligned}
 \Delta P-\lambda\Delta V
  &=R\,[2\sin\theta-\theta-\sin\theta\cos\theta]\\
  &={R\theta^3\over3}+O(R\theta^5)>0                         \tag{1.6}
\end{aligned}
\]

for all sufficiently small nonzero \(\theta\).  The chord does save raw
length, but after its lost cap volume is restored the total perimeter
increases.  This is exactly the leading distinction between (2.1) in
`cheeger_facet_completion.md` and a macroscopic smooth bevel.

## 2. Distant patch pairing reduces to the stability form

Let \(f\) be a smooth normal speed satisfying

\[
                              \int_\Sigma f\,d\sigma_\mu=0.   \tag{2.1}
\]

After the standard second-order volume correction, the constrained
second variation is

\[
 Q_\Sigma(f)=\int_\Sigma
   \left(|\nabla_\Sigma f|^2-
       (|II|^2+\nabla^2V(\nu,\nu))f^2\right)d\sigma_\mu
       +Q_{\rm contact}(f).                                  \tag{2.2}
\]

For a convex hard support with the natural free-boundary sign,
\(Q_{\rm contact}\) is the usual boundary term; exact stability means
\(Q_\Sigma(f)\ge0\) for every (2.1).

Suppose two distant patches are moved with speeds \(f_1,f_2\), chosen so
that their volume changes cancel.  If their supports are disjoint, the
interior quadratic form (2.2) has no cross term:

\[
                              Q_\Sigma(f_1+f_2)
                                  =Q_\Sigma(f_1)+Q_\Sigma(f_2). \tag{2.3}
\]

Consequently, pairing remote normal patches does not by itself create a
negative interaction.  Any gain has to come from a genuinely global
operation that changes which sheets are connected (a union/intersection,
cut-and-paste, or flow-cell reassignment), rather than from balancing two
local chord replacements.

There is a second structural obstruction in a log-affine interior.  If
\(V(x)=a\cdot x+b\), translating a weighted CMC hypersurface by a constant
vector leaves its weighted mean curvature unchanged.  Linearizing gives

\[
          J_\Sigma(\theta\cdot\nu)=0,
 \qquad
          J_\Sigma=\Delta_\Sigma-
             \nabla_\Sigma V\cdot\nabla_\Sigma+|II|^2.       \tag{2.4}
\]

Thus the normal components of ambient translations are Jacobi fields.
After imposing the one volume constraint, as many as \(n-1\) independent
translation modes can remain neutral.  High normal rank is therefore not
in tension with local stability in a flat-potential region; hard-support
contacts and global collisions are the data that break this neutrality.

### 2.1 A genuine global gain: free-boundary anchoring

There is one elementary way in which exact global minimality goes beyond
the stability form.  Suppose \(\Omega\) is bounded and convex,
\(V(x)=a\cdot x+b\), and the closure of a smooth balanced candidate
\(E\) is compactly contained in \(\Omega\).  Choose a translation direction
\(\theta\) with \(a\cdot\theta\ge0\), and translate until first contact
with \(\partial\Omega\).  Such a sign can be chosen because boundedness
allows contact in both directions.  Before contact,

\[
 \mu(E+t\theta)=e^{-t a\cdot\theta}\mu(E),
 \qquad
 P_\mu(E+t\theta)=e^{-t a\cdot\theta}P_\mu(E).                \tag{2.5}
\]

The mass does not exceed one half, so the Cheeger ratio stays exactly
constant.  Just after a regular first contact, intersect the translate
with \(\Omega\).  The new cut lying in \(\partial\Omega\) is free.  At a
nondegenerate tangency of penetration depth \(s\), the removed interface
has weighted area \(\asymp s^{(n-1)/2}\), whereas the removed mass is
\(O(s^{(n+1)/2})\).  Hence

\[
 {P_\mu((E+t\theta)\cap\Omega)\over
   \mu((E+t\theta)\cap\Omega)}
 <{P_\mu(E)\over\mu(E)}                                     \tag{2.6}
\]

for small positive penetration.  Flat or higher-order contact gives the
same conclusion whenever the removed boundary-to-volume ratio diverges.
Thus a globally minimizing smaller phase cannot float strictly inside a
bounded log-affine support; it must be anchored to the free boundary.

The conclusion applies componentwise.  If a balanced Cheeger set is the
disjoint union of positive-mass components \(E_i\), then

\[
 P_\mu(E_i)\ge\psi\mu(E_i),
 \qquad
 \sum_iP_\mu(E_i)=\psi\sum_i\mu(E_i).                         \tag{2.7}
\]

Every inequality in (2.7) is therefore an equality.  Each component is
itself a smaller-side Cheeger minimizer and is subject to the anchoring
argument.

This anchoring move is genuinely global and is invisible to (2.2).  It
does not yet complete tangent slices, but it narrows the missing theorem:
in the compact log-affine case, dispersed nearly flat sheets must transmit
their incidence through support contacts or through collisions with other
sheets.

## 3. What two-sided killed rays really imply

The following pointwise identity is the strongest direct smooth
consequence of excellent two-sided normal tubes in the log-affine case.

### Lemma 3.1 (two-sided defect controls the full second fundamental form)

Let \(V(x)=a\cdot x+b\) near a smooth point of \(\Sigma\), and let
\(\kappa_1,\ldots,\kappa_{n-1}\) be the principal curvatures with the
orientation for which

\[
                         \lambda=\sum_j\kappa_j-a\cdot\nu.    \tag{3.1}
\]

Assume both normal rays are regular up to distance \(T\), so
\(|T\kappa_j|<1\).  Their weighted Jacobians have the exact forms

\[
\begin{aligned}
 j_+(T)&=e^{\lambda T-D_+(T)},
 &D_+(T)&=\sum_j[T\kappa_j-\log(1+T\kappa_j)],\\
 j_-(T)&=e^{-\lambda T-D_-(T)},
 &D_-(T)&=\sum_j[-T\kappa_j-\log(1-T\kappa_j)].              \tag{3.2}
\end{aligned}
\]

Consequently

\[
\begin{aligned}
 D_+(T)+D_-(T)
   &=-\sum_j\log(1-T^2\kappa_j^2)\\
   &\ge T^2\sum_j\kappa_j^2=T^2|II|^2.                      \tag{3.3}
\end{aligned}
\]

In particular, \(D_+(T),D_-(T)\le\eta\) implies (0.2).

#### Proof

The normal exponential map has Euclidean Jacobian
\(\prod_j(1\pm T\kappa_j)\), while the affine density ratio is
\(e^{\mp T a\cdot\nu}\).  Taking logarithms and using (3.1) proves
(3.2).  Adding the two expressions and applying
\(-\log(1-s)\ge s\) for \(0\le s<1\) proves (3.3).  QED.

For a packet with \(T=\gamma/\psi\), Lemma 3.1 gives

\[
                              |II|\le {\sqrt{2\eta}\over\gamma}\psi.
                                                                    \tag{3.4}
\]

Thus a hypothetical small-\(\psi\) packet is flat at the physical scale
\(1/\psi\).  However, (3.4) is pointwise along the surviving base set.  It
does not say that a good patch has tangential radius \(c/\psi\), that its
tangent hyperplane is filled, or that two patches with different normals
can be charged to disjoint competitors.  The next example shows that the
other local hypotheses do not supply those missing conclusions.

### 3.1 Exact zero defect does propagate globally

There is one exact smooth analogue of the polyhedral completion theorem.
Assume \(V\) is affine and a connected interior component of \(\Sigma\)
is a smooth embedded weighted-CMC hypersurface whose only boundary lies on
the hard support.  The local graph equation

\[
 \operatorname{div}{\nabla u\over\sqrt{1+|\nabla u|^2}}
   -a\cdot{(-\nabla u,1)\over\sqrt{1+|\nabla u|^2}}
                              =\lambda                         \tag{3.5}
\]

is analytic and uniformly elliptic on every bounded-gradient chart.
Hence the component is real analytic.  If the two-sided defect in Lemma
3.1 vanishes on a subset of positive surface measure, (3.3) gives
\(II=0\) on that subset.  The zero set of a nonzero real-analytic tensor
has surface measure zero, so \(II\equiv0\) on the connected component.
It is therefore contained in one affine hyperplane.  Embeddedness and the
assumption that it has no interior boundary make it the complete slice of
that hyperplane through \(\Omega\).

Thus **exact** zero defect plus positive-measure coverage recovers slice
completion and the central-cell proof of
`cheeger_facet_completion.md`.  What is missing in the actual packet is a
dimension-free quantitative unique-continuation theorem converting
\(D_++D_-\ll1\) into macroscopic hyperplane completion.  Ordinary analytic
continuation has no such stable conclusion: small curvature on one region
does not control remote incidence.

## 4. A radial-exponential family

For \(n\ge2\), put \(c_n=\sqrt{n+1}\), and let \(\mu_n\) have density

\[
                       d\mu_n(x)=Z_n^{-1}e^{-c_n|x|}\,dx.     \tag{4.1}
\]

If \(X=RU\), where \(U\) is uniform on \(S^{n-1}\), then \(R\) and
\(U\) are independent and

\[
                    c_nR\sim\operatorname{Gamma}(n,1),
 \qquad
                    \mathbb ER^2={n(n+1)\over c_n^2}=n.      \tag{4.2}
\]

Hence \(\mu_n\) is isotropic.  Let \(s_n\) be the median of
\(\operatorname{Gamma}(n,1)\), put \(r_n=s_n/c_n\), and set

\[
                             E_n=B(0,r_n).                    \tag{4.3}
\]

Then \(\mu_n(E_n)=1/2\).  The standard gamma-median bounds give

\[
                             n-\frac13\le s_n\le n.           \tag{4.4}
\]

The exterior Minkowski boundary measure of \(E_n\) is the radial density
at \(r_n\):

\[
 p_n=P_{\mu_n}(E_n)
     ={c_n s_n^{\,n-1}e^{-s_n}\over\Gamma(n)}.               \tag{4.5}
\]

Stirling's inequalities and (4.4) give universal constants, for example

\[
                              {1\over10}\le p_n\le1           \tag{4.6}
\]

for every \(n\ge2\).  No optimized constants are needed below.

The law itself has universal Cheeger scale.  This follows from Bobkov's
radial spectral-gap theorem: if a probability on \(\mathbb R^n\),
\(n\ge2\), is both log-concave and spherically symmetric, then its
Poincare constant satisfies

\[
 {\mathbb E|X|^2\over n}\le C_P(\mu)
      \le {13\mathbb E|X|^2\over n}.                         \tag{4.7}
\]

The hypotheses hold for (4.1), and \(\mathbb E|X|^2=n\) by (4.2).
Thus \(C_P(\mu_n)\le13\), and
the Buser--Ledoux inequality for log-concave measures gives

\[
                              c_0\le\psi_{\mu_n}\le2p_n\le2  \tag{4.8}
\]

with a universal \(c_0>0\).  This confirms that the countermodel is at
the correct universal scale; it is not a KLS counterexample.

## 5. The median sphere passes every local high-rank test

### 5.1 CMC multiplier

On \(\Sigma_n=\partial B(0,r_n)\), the density is constant and

\[
              \lambda_n={n-1\over r_n}-c_n
                  =c_n\left({n-1\over s_n}-1\right).          \tag{5.1}
\]

Equation (4.4) implies

\[
                              |\lambda_n|\le {3\over\sqrt n}. \tag{5.2}
\]

In view of (4.8), \(|\lambda_n|\le\psi_{\mu_n}\) for all sufficiently
large \(n\).

### 5.2 High normalized normal variance

Rotational invariance gives the exact matrix identity

\[
 \int_{\Sigma_n}\nu\otimes\nu\,d\sigma_{\mu_n}
                           ={p_n\over n}I_n.                  \tag{5.3}
\]

Thus the normalized normal matrix is \(Q_n=I_n/n\), of effective rank
\(n\).

### 5.3 Volume-constrained stability

Since \(V(x)=c_n|x|\), one has
\(\nabla^2V(\nu,\nu)=0\) on the sphere.  Moreover

\[
                              |II|^2={n-1\over r_n^2}.         \tag{5.4}
\]

For every smooth \(f\) with zero surface mean, the exact first nonzero
eigenvalue of the radius-\(r_n\) sphere is \((n-1)/r_n^2\).  Therefore

\[
 \int_{\Sigma_n}\left(|\nabla_\Sigma f|^2
          -{|II|^2}f^2\right)d\sigma_{\mu_n}\ge0.            \tag{5.5}
\]

The median sphere is thus volume-constrained stable.  The first spherical
harmonics are neutral, illustrating the translation-Jacobi obstruction in
(2.4).

### 5.4 Entire-boundary killed tubes

No outward ray is killed, and an inward ray is regular until time
\(r_n\).  For \(0\le t<r_n\), the exact normalized weighted Jacobians are

\[
\begin{aligned}
 j_+(t)&=\left(1+{t\over r_n}\right)^{n-1}e^{-c_nt}
              =e^{\lambda_nt-D_+(t)},\\
 D_+(t)&=(n-1)\left[{t\over r_n}
                     -\log\left(1+{t\over r_n}\right)\right],\\
 j_-(t)&=\left(1-{t\over r_n}\right)^{n-1}e^{c_nt}
              =e^{-\lambda_nt-D_-(t)},\\
 D_-(t)&=(n-1)\left[-{t\over r_n}
                     -\log\left(1-{t\over r_n}\right)\right]. \tag{5.6}
\end{aligned}
\]

For \(0\le T\le1\) and all sufficiently large \(n\), (4.4) gives
\(T/r_n\le1/2\) and

\[
                              D_+(T),D_-(T)\le2T^2.            \tag{5.7}
\]

Thus every boundary point survives a fixed two-sided tube with arbitrarily
small normalized defect once the fixed tube length is chosen small.  In
particular, if \(T=h/p_n\) and \(h\le10^{-3}\), then (4.6) gives
\(T\le10^{-2}\), the killed fraction is zero, and the convexity loss is
at most \(2\cdot10^{-4}\).  Together with (5.2)--(5.5), this verifies all
local/tube/stability hypotheses in the dispersed high-rank regime.

The density in (4.1) is nonsmooth only at the origin, which lies far from
the entire tube.  If global smoothness is desired, replace \(r\) near
zero by a smooth convex function \(\varphi(r)\) with
\(\varphi'(0)=0\) and \(\varphi(r)=r\) for \(r\ge1\), then recenter and
rescale isotropically.  In high dimension the changed mass is
superpolynomially small, the median sphere remains in the exactly linear
region, and (4.4)--(5.7) persist with fixed slack.

### 5.5 A macroscopically incompatible comparison has the same leading cost

The hyperplane through the origin also bisects \(\mu_n\).  Its boundary
measure is

\[
 h_n={c_n\over n-1}\,{\Gamma(n/2)\over
               \sqrt\pi\,\Gamma((n-1)/2)}.                   \tag{5.8}
\]

Stirling's formula, the interval (4.4), and the fact that the gamma
density changes by only \(O(1/n)\) between its mode \(n-1\) and that
interval give

\[
 p_n={1\over\sqrt{2\pi}}+O(n^{-1}),
 \qquad
 h_n={1\over\sqrt{2\pi}}+O(n^{-1}).                          \tag{5.9}
\]

Thus a high-rank spherical interface and a rank-one complete slice have
the same leading balanced perimeter, despite completely different normal
geometry.  Equation (5.9) does not claim either is the global minimizer;
it shows why a local normal-variance discriminator has no fixed gap even
in this explicit family.

## 6. What the countermodel rules out

The radial family proves that none of the following implications can be
the missing high-rank inverse:

\[
\begin{gathered}
 \text{high normal rank + bounded CMC + two-sided long rays}
       \Longrightarrow \text{one coherent normal direction},\\
 \text{the same hypotheses + volume-constrained stability}
       \Longrightarrow \text{a negative paired-patch variation},\\
 \text{small killed defect}
       \Longrightarrow \text{full tangent-hyperplane completion}.
                                                               \tag{6.1}
\end{gathered}
\]

The sphere has no complete tangent slice patch at all, yet it passes every
local test.  Its universal perimeter is detected globally by the radial
one-dimensional density (4.5), not by local beveling.

There is also a scale lesson.  Lemma 3.1 turns a tube of length
\(T\asymp1/\psi\) into \(|II|\lesssim\psi\) on the good packet.  Normal
directions separated by a fixed angle must then be connected, if at all
through the good set, by intrinsic paths of length \(\gtrsim1/\psi\).
Isotropy and the tube covariance bound only charge those paths direction
by direction and lose the effective rank.  A bounded-reuse theorem must
instead show that the *global phase incidence* of these long nearly flat
pieces either completes marginal slices or allows sheets to be
reconnected with a fixed fraction of perimeter saved.

## 7. The remaining macroscopic completion statement

The exact polyhedral proof suggests the following smooth replacement, but
the radial countermodel shows that its conclusion must include the
universal-perimeter alternative.

> **Macroscopic calibrated-cell alternative.**  Let \(E\) be an exact
> balanced Cheeger minimizer and let \(G\subset\partial^*E\) be a Wulff
> packet sweeping fixed mass to distance \(T=\gamma/\psi\), with two-sided
> defect at most \(\eta\).  Then at least one of the following holds:
>
> 1. \(\psi\ge c\);
> 2. a subpacket of anisotropic area \(cp\) admits tangent-plane
>    completion with total defect at most \(Cp\);
> 3. the common Cheeger calibration produces a union/intersection or
>    flow-cell reassignment competitor saving \(cp\), with every physical
>    interface charged at most \(C\) times.

The radial sphere belongs to branch 1; the exact polyhedral theorem proves
branches 2--3 when curvature is concentrated into genuine ridges.  A proof
for small nonzero curvature must use exact global minimality before taking
the smooth limit.  Local second variation cannot do it because of
(1.4), (2.3), and the neutral modes (2.4).

One concrete route is to construct Wulff flow cells from the **single
calibration shared by all direct-deficit level sets**.  The desired
measurable analogue of the polyhedral central-cell lemma would assert that
if all smaller completed halfspaces have total calibrated weight below a
fixed constant, their larger flow cells have intersection mass above one
half and carry a constant phase, contradicting balance.  Curvature error
must be integrated as a signed calibration flux, rather than bounded by
the absolute local chord deficit; otherwise the CMC cancellation in
Section 1 loses the required sign.

No such bounded-reuse calibrated-cell theorem is proved here.  The exact
advance is (0.2), and the exact obstruction is the radial family: the
remaining step is global, not a missing local curvature estimate.
