# Multi-cell chord cancellation: exact ledger and sharp obstruction

## 0. Verdict

Jointly canceling the volume changes of several disjoint replacements
removes the common CMC multiplier exactly, but it does not create a new
perimeter gain. For each cell \(i\), write

\[
 v_i=\Delta V_i,\qquad p_i=\Delta P_i,\qquad
 r_i=p_i-\lambda v_i,                                      \tag{0.1}
\]

where \(\lambda\) is the global volume multiplier. If the cells have
disjoint interiors, matching traces on their artificial boundaries, and

\[
                              \sum_i v_i=0,                 \tag{0.2}
\]

then the assembled competitor satisfies the exact identity

\[
                              \Delta P=\sum_i r_i.           \tag{0.3}
\]

Thus multi-cell cancellation is successful if and only if the selected
cells have negative total *constrained residual*. Contact-tensor mass and
small longitudinal support turning give capacity and small lateral errors,
but neither determines the sign of \(r_i\).

For a circular CMC arc, replacing an arc of central angle \(\theta\) by
its chord gives

\[
 r_i=\rho\left[
        2\sin{\theta\over2}-{\theta+\sin\theta\over2}\right]
 >0.                                                        \tag{0.4}
\]

All chord replacements of an oriented CMC interface also have the same
sign of volume change. Hence they cannot satisfy (0.2) without an inverse
operation, and after adding such operations stability makes the quadratic
form on the volume-zero subspace nonnegative. The same obstruction holds
in every dimension for a spherical patch: its exact residual is positive.
Consequently disk and spherical-cap minimizers never trigger the proposed
criterion.

There is one robust positive case. A genuine ridge of opening angle
\(\alpha<\pi\) has an \(O(\varepsilon)\) bevel saving and only an
\(O(\varepsilon^2)\) volume error:

\[
\begin{aligned}
 p_{\rm ridge}(\varepsilon)
   &=-2\rho_0\bigl(1-\sin(\alpha/2)\bigr)\varepsilon
       +O(\rho_0M\varepsilon^2),\\
 v_{\rm ridge}(\varepsilon)
   &=\pm{\rho_0\over2}\sin\alpha\,\varepsilon^2
       +O(\rho_0M\varepsilon^3).                            \tag{0.5}
\end{aligned}
\]

Its negative constrained residual survives exact joint volume
cancellation, because a smooth balancing cell needs amplitude
\(O(\varepsilon^2)\) and costs only \(O(\varepsilon^4)\) in constrained
perimeter. This recovers the product-exponential ridge exclusion.

The sharp algebraic obstruction is therefore:

> If every available smooth collar action has \(r_i\ge0\), no finite or
> infinitesimal joint choice with exact volume cancellation decreases
> perimeter. A dimension-free negative criterion must produce a universal
> negative residual-to-volume balance on two signed-volume packets, or a
> first-order ridge defect. Longitudinal turning and contact flux alone do
> not do so.

## 1. The exact finite-cell ledger

Let \(K\subset\mathbb R^n\) be convex and let
\(d\mu=\rho\,dx\), with \(\rho>0\) continuous. Let \(E\subset K\) have a
smooth weighted-CMC interior interface with volume multiplier
\(\lambda\). Let \(U_1,\ldots,U_N\) be pairwise disjoint open sets. In
each \(U_i\), choose a finite-perimeter replacement \(F_i\) whose trace
on \(\partial U_i\cap\operatorname{int}K\) agrees with that of \(E\).
Portions of \(\partial U_i\) lying on \(\partial K\) are free and are not
charged by relative perimeter.

Define

\[
\begin{aligned}
 v_i&=\mu(F_i\cap U_i)-\mu(E\cap U_i),\\
 p_i&=P_\mu(F_i;U_i)-P_\mu(E;U_i),\\
 r_i&=p_i-\lambda v_i.                                     \tag{1.1}
\end{aligned}
\]

Assemble

\[
 E'=\left(E\setminus\bigcup_iU_i\right)
       \cup\bigcup_i(F_i\cap U_i).                          \tag{1.2}
\]

### Proposition 1.1 (exact additivity)

The matching-trace and disjointness hypotheses imply

\[
 \mu(E')-\mu(E)=\sum_i v_i,\qquad
 P_\mu(E')-P_\mu(E)=\sum_i p_i.                             \tag{1.3}
\]

In particular, if \(\sum_i v_i=0\), then

\[
                   P_\mu(E')-P_\mu(E)=\sum_i r_i.           \tag{1.4}
\]

#### Proof

Volume is finitely additive on the disjoint \(U_i\). Relative perimeter
is local, and equality of traces removes every artificial jump term on
\(\partial U_i\cap\operatorname{int}K\). Hence (1.3) holds. Substituting
\(p_i=\lambda v_i+r_i\) and using \(\sum_i v_i=0\) gives
(1.4). QED.

No approximation or first-variation argument occurs in Proposition 1.1.
It is the exact ledger for finite chord and strip replacements.

## 2. Signed amplitudes and the stability obstruction

Suppose the \(i\)-th cell has a two-sided \(C^2\) variation with amplitude
\(t_i\), and write

\[
\begin{aligned}
 v_i(t_i)
  &=a_it_i+\frac12b_it_i^2+o(t_i^2),\\
 p_i(t_i)
  &=\lambda a_it_i+\frac12c_it_i^2+o(t_i^2).                \tag{2.1}
\end{aligned}
\]

The common linear coefficient is precisely weighted-CMC stationarity.
Choose

\[
                         t_i=\varepsilon z_i+\varepsilon^2w_i
                              +o(\varepsilon^2).             \tag{2.2}
\]

First-order volume cancellation requires

\[
                              \sum_i a_i z_i=0.              \tag{2.3}
\]

The second-order correction can be selected so that

\[
             \sum_i a_iw_i=-{1\over2}\sum_i b_i z_i^2.      \tag{2.4}
\]

Then the total volume is \(o(\varepsilon^2)\), and an implicit-function
correction makes it exactly zero whenever at least one \(a_i\ne0\).
Substitution into the perimeter expansion gives

\[
 \Delta P={\varepsilon^2\over2}
             \sum_i(c_i-\lambda b_i)z_i^2
             +o(\varepsilon^2).                             \tag{2.5}
\]

For disjoint supports the constrained second-variation form is additive,
and its restriction to (2.3) is exactly the quadratic form in (2.5).
Volume-constrained stability says

\[
                \sum_i(c_i-\lambda b_i)z_i^2\ge0
                \quad\text{whenever }\sum_i a_i z_i=0.      \tag{2.6}
\]

This is the infinitesimal sharp obstruction. Joint choice of signed
amplitudes eliminates the multiplier term \(\lambda\sum_i v_i\), but the
remaining quadratic form is the very form that stability makes
nonnegative.

There is an even simpler obstruction for one-sided chord actions. On an
oriented planar CMC interface with nonzero signed curvature, every
inscribed chord cuts volume on the same side. Thus all nontrivial
\(v_i\) have the same sign and cannot obey (0.2) at all. Opposite signed
volume requires an inverse bulging action or another kind of cell, after
which (2.6) applies.

## 3. Exact circular-arc calculation

Let an oriented CMC arc have radius \(\rho>0\), multiplier
\(\lambda=1/\rho\), and central angle \(\theta\in(0,\pi)\). Orient the
enclosed phase so replacing the arc by its chord removes the circular
segment. Put

\[
\begin{aligned}
 S(\theta)
  &=\rho\left(\theta-2\sin{\theta\over2}\right),\\
 A(\theta)
  &={\rho^2\over2}(\theta-\sin\theta).                      \tag{3.1}
\end{aligned}
\]

Here \(S>0\) is the raw length saving and \(A>0\) the lost area. Therefore

\[
                 p=-S(\theta),\qquad v=-A(\theta).          \tag{3.2}
\]

The constrained residual is

\[
\begin{aligned}
 r(\theta)
 &=p-\lambda v\\
 &=\rho\left[
       2\sin{\theta\over2}-{\theta+\sin\theta\over2}\right]
   =:\rho F(\theta).                                        \tag{3.3}
\end{aligned}
\]

Since

\[
 F(0)=0,\qquad
 F'(\theta)=\cos(\theta/2)\,[1-\cos(\theta/2)]>0,            \tag{3.4}
\]

we have

\[
 r(\theta)>0,\qquad
 r(\theta)={\rho\theta^3\over24}
              +O(\rho\theta^5).                             \tag{3.5}
\]

Equivalently, the chord saves strictly less perimeter than the CMC cost
of its volume error:

\[
                              S(\theta)<\lambda A(\theta).   \tag{3.6}
\]

For finitely many disjoint circular subarcs, every chord has \(v_i<0\)
and \(r_i>0\). Pure chords cannot cancel volume; if other cells add the
lost volume with nonnegative constrained residual, (1.4) is strictly
positive.

The exact relative isoperimetric caps in a Euclidean disk are circular
arcs meeting the disk orthogonally. Thus (3.3)--(3.6) audit a smooth
global minimizer, not only a formal stationary arc. The disk
classification is contained in A. Ros and E. Vergasta, *Stability for
hypersurfaces of constant mean curvature with free boundary*,
Geom. Dedicata **56** (1995), 19--33,
DOI 10.1007/BF01263611.

## 4. Exact spherical-patch calculation in every dimension

Let the interface dimension be \(m\ge1\). Consider a spherical patch of
radius \(\rho\) and polar half-angle
\(\theta\in(0,\pi/2]\), and replace it by its flat \(m\)-disk. Let
\(\kappa_m\) be the volume of the unit \(m\)-ball. The old area, new area,
and lost \((m+1)\)-volume are

\[
\begin{aligned}
 A_{\rm sph}
  &=m\kappa_m\rho^m\int_0^\theta\sin^{m-1}t\,dt,\\
 A_{\rm disk}
  &=\kappa_m\rho^m\sin^m\theta,\\
 V_{\rm seg}
  &=\kappa_m\rho^{m+1}\int_0^\theta\sin^{m+1}t\,dt.          \tag{4.1}
\end{aligned}
\]

The sum-mean-curvature multiplier is \(\lambda=m/\rho\).
With \(p=A_{\rm disk}-A_{\rm sph}\) and \(v=-V_{\rm seg}\),

\[
 r_m(\theta)=\kappa_m\rho^mG_m(\theta),                     \tag{4.2}
\]

where

\[
 G_m(\theta)
 =\sin^m\theta(1-\cos\theta)
       -\int_0^\theta\sin^{m+1}t\,dt.                       \tag{4.3}
\]

Indeed, integration of

\[
 {d\over dt}(\sin^m t\cos t)
 =m\sin^{m-1}t\cos^2t-\sin^{m+1}t
\]

reduces the direct expression for \(p-\lambda v\) to (4.3). Moreover,

\[
 G_m(0)=0,\qquad
 G_m'(\theta)
 =m\sin^{m-1}\theta\cos\theta(1-\cos\theta)>0               \tag{4.4}
\]

for \(0<\theta<\pi/2\). Hence

\[
 r_m(\theta)>0,\qquad
 G_m(\theta)
 ={m\over2(m+2)}\theta^{m+2}
       +O_m(\theta^{m+4}).                                  \tag{4.5}
\]

Thus the multi-cell chord criterion does not falsely fire on spherical
caps in any dimension. The positivity is exact and finite, rather than a
consequence of an untracked stability constant.

## 5. The sharp signed-packet criterion

For any admissible action \(a\), retain only its exact pair

\[
                         (v(a),r(a)).                       \tag{5.1}
\]

Suppose a positive-volume packet can realize total volume \(M>0\) with
minimal residual cost \(C_+(M)\), and a negative-volume packet can realize
total volume \(-M\) with minimal residual cost \(C_-(M)\). Proposition
1.1 gives the exact criterion

\[
 \boxed{\quad
 \text{a volume-preserving perimeter decrease exists}
 \iff
 C_+(M)+C_-(M)<0\ \text{for some }M>0.
 \quad}                                                     \tag{5.2}
\]

In a divisible infinitesimal packet, define its best residual cost per
unit absolute volume by

\[
 c_+=\inf_{v>0}{r\over v},\qquad
 c_-=\inf_{v<0}{r\over -v}.                                 \tag{5.3}
\]

Then the dimension-free local criterion is

\[
                              c_++c_-\le-\kappa             \tag{5.4}
\]

for some universal \(\kappa>0\), together with enough capacity on the two
packets to match a common \(M\). Conversely, if \(r(a)\ge0\) for every
available action, then \(C_\pm(M)\ge0\) and cancellation is impossible.
This converse is sharp by the circular and spherical computations.

The contact tensor can provide packet capacity. If \(U=B/\|B\|_F\) and

\[
 b_+(x)=\max\left(
   \langle U,(x-c)\otimes n_K(x)\rangle_F,0\right),
\]

then

\[
                              \int_\Gamma b_+\,d\eta\ge p\|B\|_F.  \tag{5.5}
\]

The longitudinal theorem in
“longitudinal_support_turning.md” can make the support-strip turning on a
fixed fraction of this flux at most
\(2\varepsilon+8\pi\zeta^2\). These facts may control the lateral-error
part of \(r(a)\) and the capacity \(M\). They do not supply the negative
term in (5.4):

### Proposition 5.1 (dimension-free packet criterion)

Suppose two disjoint contact packets, labeled \(+\) and \(-\), have
longitudinal turning at most \(\tau_+\) and \(\tau_-\). Assume that for
every \(0<M\le M_0\) they can realize actions with

\[
\begin{aligned}
 v_+(M)&=M,&
 r_+(M)&\le(C\tau_+-d_+)M,\\
 v_-(M)&=-M,&
 r_-(M)&\le(C\tau_--d_-)M.                                 \tag{5.6}
\end{aligned}
\]

Then their exact joint replacement preserves volume and satisfies

\[
 \Delta P\le
 \left[C(\tau_++\tau_-)-(d_++d_-)\right]M.                 \tag{5.7}
\]

It decreases perimeter by at least \(\kappa M\) whenever

\[
                  d_++d_-\ge C(\tau_++\tau_-)+\kappa.      \tag{5.8}
\]

All constants are dimension-free in isotropic units. If each packet
carries contact flux \(W_\pm\) and the geometric construction has capacity
\(M_0\ge a\min(W_+,W_-)\), (5.5) prevents the available \(M\) from
vanishing in the high-rank contact-tensor branch.

#### Proof

Add the two volume and residual identities in (5.6), and apply
Proposition 1.1. QED.

Proposition 5.1 is the requested turning/contact-packet formulation. Its
load-bearing quantities are \(d_\pm\), the interface straightening gains
per unit volume after subtracting the CMC multiplier. Longitudinal
turning controls only the displayed \(C\tau_\pm\) errors. The exact
cap computations show that \(d_\pm\) cannot be inferred from small
\(\tau_\pm\).

- a cube slice has large contact flux and zero turning, but \(r=0\);
- the quartic slice has large contact flux and
  \(O(n^{-1/2})\) isotropic turning, but again \(r=0\);
- a disk or spherical cap has smooth support incidence, but \(r>0\).

Therefore any proposed geometric version of (5.4) must include a measured
interface angle defect, nonstationarity, or another source of negative
constrained residual. A condition involving only \(B\) and longitudinal
support turning necessarily falsely fires on one of these three tests.

## 6. Genuine ridges survive exact volume cancellation

Let two planar interface facets meet at opening angle
\(\alpha\in(0,\pi)\), and let \(\rho\) be a positive \(C^1\) density near
the ridge point, with

\[
 \rho(0)=\rho_0,\qquad |\nabla\log\rho|\le M.                \tag{6.1}
\]

Cut points at distance \(\varepsilon\) along the two rays and replace the
two old segments by their chord. In the constant-density model,

\[
\begin{aligned}
 p_{\rm bevel}(\varepsilon)
  &=-2\rho_0\left(1-\sin{\alpha\over2}\right)\varepsilon,\\
 v_{\rm bevel}(\varepsilon)
  &=\pm{\rho_0\over2}\sin\alpha\,\varepsilon^2.              \tag{6.2}
\end{aligned}
\]

Put \(d_\alpha=1-\sin(\alpha/2)>0\). If
\(M\varepsilon\le d_\alpha/8\), the density bounds

\[
              e^{-M\varepsilon}\rho_0
              \le\rho\le e^{M\varepsilon}\rho_0             \tag{6.3}
\]

on the bevel triangle imply

\[
\begin{aligned}
 p_{\rm bevel}(\varepsilon)
   &\le-d_\alpha\rho_0\varepsilon,\\
 |v_{\rm bevel}(\varepsilon)|
   &\le e^{M\varepsilon}\rho_0\varepsilon^2.                \tag{6.4}
\end{aligned}
\]

Consequently, if also

\[
                         |\lambda|\varepsilon
                         \le {d_\alpha\over4},               \tag{6.5}
\]

then

\[
 r_{\rm bevel}(\varepsilon)
 =p_{\rm bevel}-\lambda v_{\rm bevel}
 \le-{d_\alpha\over2}\rho_0\varepsilon.                    \tag{6.6}
\]

The constants in (6.4) are deliberately weakened; they follow from
\(e^x\le1+2x\) and \(e^{-x}\ge1-x\) for \(0\le x\le1/2\).

### Proposition 6.1 (exact joint repair of a ridge)

Assume a disjoint smooth balancing patch has a two-sided \(C^2\) variation
\(G_s\) satisfying

\[
\begin{aligned}
 v_{\rm bal}(s)&=a s+O(A s^2),\qquad |a|\ge a_0>0,\\
 r_{\rm bal}(s)&=O(A s^2).                                  \tag{6.7}
\end{aligned}
\]

Then, for every sufficiently small ridge bevel, there is a unique
\(s(\varepsilon)=O(\rho_0\varepsilon^2/a_0)\) such that

\[
             v_{\rm bevel}(\varepsilon)
             +v_{\rm bal}(s(\varepsilon))=0.                \tag{6.8}
\]

For that exact volume-preserving two-cell replacement,

\[
 \Delta P
 \le-{d_\alpha\over2}\rho_0\varepsilon
      +O\left({A\rho_0^2\varepsilon^4\over a_0^2}\right)<0. \tag{6.9}
\]

#### Proof

The implicit-function theorem applied to (6.7) gives (6.8) and the stated
bound on \(s\). Proposition 1.1 and (6.6)--(6.7) then give (6.9). QED.

If two ridge packets have opposite volume signs, the balancing patch is
unnecessary: tune their bevel radii so their \(O(\varepsilon^2)\) volumes
agree. Both \(O(\varepsilon)\) residuals remain negative.

A ridge rounded only inside distance \(\delta\) changes the calculation
by \(O(\rho_0\delta)\). Choosing

\[
                 \delta\le {d_\alpha\over8C}\varepsilon     \tag{6.10}
\]

preserves a negative residual. At a right-angle ridge,

\[
                 d_{\pi/2}=1-2^{-1/2}.                     \tag{6.11}
\]

The product-exponential maximum interface has precisely these
right-angle facet transitions. In its log-affine chambers \(M\) and
\(|\lambda|\) are fixed-scale quantities. Taking
\(\delta\ll\varepsilon\ll1\) makes (6.6) hold. Hence the criterion fires
on the smoothed ridge model, as it must; exact global minimality excludes
that model before any smooth contact-tensor argument.

## 7. Stress tests

### 7.1 Exact disk and spherical caps

Equations (3.3) and (4.2) give \(r_i>0\) for every chord patch. Pure
chords have one-signed volume, and signed smooth combinations lie under
the nonnegative constrained Hessian. No false perimeter decrease occurs.

### 7.2 Cube flat slices

For a complete flat slice, chord replacement changes neither interface
nor volume:

\[
                              v_i=p_i=r_i=0.                 \tag{7.1}
\]

The contact tensor is large, but the example belongs to central-cell
completion. Criterion (5.4) does not fire.

### 7.3 Quartically flat support

For
\(K_n=\{x_1^4+|x'|^2\le1\}\) and
\(\Sigma=\{x_1=0\}\), the support has no exact \(e_1\)-ruling, while
\(B=I-e_1\otimes e_1\). Nevertheless the interface itself is flat, so
(7.1) still holds. The \(O(n^{-1/2})\) longitudinal turning obtained after
isotropic normalization controls only lateral approximation error; it
does not turn zero residual into a negative one.

### 7.4 Product-exponential smoothed ridges

At scale \(\varepsilon\) larger than the smoothing radius, (6.6) gives a
negative \(O(\varepsilon)\) residual and an \(O(\varepsilon^2)\) volume
change. Proposition 6.1 cancels volume exactly without changing the
leading sign. This is the one stress test on which the negative criterion
fires.

## 8. Final obstruction and revised target

The multi-cell proposal does not repair the smooth CMC chord mechanism:

\[
 \underbrace{\sum_i\lambda v_i}_{=0}
 \quad\text{is removed, but}\quad
 \underbrace{\sum_i r_i}_{\text{constrained defect}}
 \quad\text{remains}.                                      \tag{8.1}
\]

For exact smooth minimizers this residual is nonnegative infinitesimally
by stability and is positive for the canonical finite chord tests. The
only verified negative mechanism is a first-order angle or ridge defect.

A viable contact-tensor theorem must therefore do more than produce
oppositely signed volume cells. It must prove one of the following:

1. a universal negative cost gap (5.4) caused by a genuinely non-CMC
   macroscopic reassignment;
2. a nonsmooth angle defect, to which Proposition 6.1 applies; or
3. a global rearrangement not decomposable into cells with individually
   nonnegative constrained residual.

The third option is the only remaining smooth route. It must exploit
nonlocal changes in adjacency or reuse; merely summing disjoint chord and
strip replacements is algebraically incapable of defeating (8.1).
