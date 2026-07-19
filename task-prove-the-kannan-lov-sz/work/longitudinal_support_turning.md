# Longitudinal support turning in exact normal sections

## 0. Verdict

The exact two-dimensional calculation is favorable, but it identifies a
new incidence hypothesis rather than closing the contact-tensor inverse.

For a convex support graph \(r=g(s)\), tangent to \(r=0\) at the origin,
the one-sided reach \(R\) at inward depth \(h\) satisfies

\[
 h=g(R)=\int_0^R(R-s)g''(s)\,ds.                            \tag{0.1}
\]

Consequently, if the graph stays below depth \(h\) out to \(2T\), its
total normal turning on \([0,T]\) is at most \(h/T\). More generally,
turning can be controlled only on an interval separated by a fixed
fraction from the last known point: curvature may concentrate
arbitrarily close to the endpoint, where the kernel in (0.1) vanishes.

This couples cleanly to an exact planar free-boundary collar. If an
interface point at intrinsic distance \(u\le\varepsilon T\) from contact
has a two-sided normal segment of length \(4T\) inside the support, and
the interface normal has not backtracked by more than \(1/4\) radian
along the collar, then the support normal turns by at most

\[
                              2\varepsilon                 \tag{0.2}
\]

on each longitudinal interval of length \(T\). The estimate can be
integrated against the positive part of the contact-tensor flux. Thus a
dimension-sized contact tensor really does force a large
*flux-weighted, nearly flat support strip* if surviving rays can be paired
with the contact tensor through a short collar.

That pairing is not supplied by a retained normal-tube theorem. The
quartic body

\[
             K_n=\{x_1^4+x_2^2+\cdots+x_n^2\le1\}           \tag{0.3}
\]

quantifies the gap. Its central flat slice has
\(\|B\|_F=\sqrt{n-1}\) and zero contact curvature in direction \(e_1\),
but no nonzero tangent segment at contact. A normal segment of half-length
\(L\) first appears only at collar depth

\[
                         1-\sqrt{1-L^4}.                    \tag{0.4}
\]

After isotropic normalization, a fixed \(L\) survives on a fixed fraction
of the interface at depth \(\asymp L^4/\sqrt n\), while the support turns
by only \(O(L^3/\sqrt n)\). Hence the quartic example is compatible with
an *approximate* completion theorem, but falsifies every exact inference
from zero contact curvature to a ruled support.

Finally, global perimeter minimality does not turn a curved CMC strip into
a chord saving. Exact relative isoperimetric caps in the Euclidean disk
are circular arcs meeting the support orthogonally. Replacing a subarc of
angle \(\theta\) by its chord and repairing the lost area has constrained
cost

\[
 \rho\left[2\sin{\theta\over2}
           -{\theta+\sin\theta\over2}\right]
 ={\rho\theta^3\over24}+O(\rho\theta^5)>0.                  \tag{0.5}
\]

The missing theorem is therefore a global contact-to-collar incidence
statement followed by a multi-cell reassignment whose volume errors
cancel. A local chord plus scalar volume repair cannot supply the desired
first-order saving.

## 1. Exact reach and the curvature kernel

Let \(g:[0,R_*]\to[0,\infty)\) be \(C^2\) and convex, with

\[
                         g(0)=g'(0)=0.                      \tag{1.1}
\]

The epigraph \(D=\{(s,r):r\ge g(s)\}\) is the local convex support, with
\(r\) measuring inward depth and \(s\) the tangent direction at contact.
Define the one-sided tangent reach at depth \(h\) by

\[
 R_+(h)=\sup\{R\in[0,R_*]:g(s)\le h
                         \text{ for every }0\le s\le R\}.   \tag{1.2}
\]

Convexity and (1.1) make \(g\) nondecreasing, so this is equivalently the
largest \(R\) with \(g(R)\le h\).

### Proposition 1.1 (kernel identity and interior turning)

For every \(R\le R_*\),

\[
\begin{aligned}
 g'(R)&=\int_0^R g''(s)\,ds,\\
 g(R)&=\int_0^R(R-s)g''(s)\,ds.                            \tag{1.3}
\end{aligned}
\]

For every \(a\in(0,1)\),

\[
 g'(aR)\le {g(R)\over(1-a)R}.                              \tag{1.4}
\]

The total turning of the support normal on \([0,aR]\) is

\[
 \Theta_g(aR)=\int_0^{aR}
       {g''(s)\over1+g'(s)^2}\,ds
       =\arctan g'(aR)
       \le {g(R)\over(1-a)R}.                              \tag{1.5}
\]

In particular, if \(g(2T)\le h\), then

\[
                             \Theta_g(T)\le {h\over T}.      \tag{1.6}
\]

#### Proof

Twice integrating \(g''\) and using (1.1) gives (1.3). Since \(g''\ge0\),

\[
 g(R)\ge\int_0^{aR}(R-s)g''(s)\,ds
       \ge(1-a)R\int_0^{aR}g''(s)\,ds,
\]

which is (1.4). The tangent angle of the graph is
\(\arctan g'\), proving (1.5)--(1.6). QED.

If \(R_+(h)<R_*\) and \(g\) is not flat at the endpoint, continuity gives

\[
 h=g(R_+(h))
  =\int_0^{R_+(h)}(R_+(h)-s)g''(s)\,ds.                    \tag{1.7}
\]

Thus (0.1) is an equality, not an asymptotic expansion.

### Endpoint loss is unavoidable

There is no bound for \(g'(R)\) in terms of \(g(R)/R\) with a universal
constant. For any \(M>0\), choose a nonnegative smooth \(g''\) supported
in \([R-\delta,R]\), with mass \(M\), and then choose
\(\delta=h/M\) up to an inessential smoothing factor. The kernel
\(R-s\) makes \(g(R)\le h\) while \(g'(R)=M\). The fixed fractional
shrink in (1.4) is therefore necessary.

Everything above has a left-hand version. Put \(g_-(s)=g(-s)\) for
\(s\ge0\); its turning will be denoted \(\Theta_-(s)\), while
\(\Theta_+(s)\) denotes (1.5).

## 2. A surviving normal ray forces longitudinal flatness

We first work in an exact planar normal section. Let the support near the
origin be

\[
               D=\{(s,r):r\ge g(s)\},\qquad
               g(0)=g'(0)=0,\quad g''\ge0,                 \tag{2.1}
\]

where \(g\) is defined on both sides of zero. Let a \(C^2\) interface
collar be parametrized by arclength:

\[
 \gamma(v)=(a(v),h(v)),\quad
 \gamma(0)=(0,0),\quad \gamma'(0)=(0,1).                   \tag{2.2}
\]

Write

\[
 \gamma'(v)=(\sin\alpha(v),\cos\alpha(v)),\qquad
 N(v)=(\cos\alpha(v),-\sin\alpha(v)),\qquad \alpha(0)=0.    \tag{2.3}
\]

The interface meets the support orthogonally at \(v=0\). Define its total
collar turning

\[
                         K_u=\int_0^u|\alpha'(v)|\,dv.       \tag{2.4}
\]

### Proposition 2.1 (collar-ray turning lemma)

Fix \(L>0\) and \(0<u\le L/8\). Assume

\[
 K_u\le {1\over4}
 \quad\text{and}\quad
 \gamma(u)+tN(u)\in D\quad\text{for every }|t|\le L.        \tag{2.5}
\]

Then

\[
 \boxed{\qquad
            \Theta_+(L/4)\le {8u\over L},
            \qquad
            \Theta_-(L/4)\le {8u\over L}.
       \qquad}                                              \tag{2.6}
\]

Equivalently, if \(L=4T\) and \(u\le\varepsilon T\), then

\[
                         \Theta_\pm(T)\le2\varepsilon.      \tag{2.7}
\]

#### Proof

Equations (2.2)--(2.4) give

\[
 |a(u)|\le\int_0^u|\alpha(v)|\,dv\le uK_u\le {L\over32},
 \qquad 0\le h(u)\le u.                                    \tag{2.8}
\]

The two endpoints of the surviving segment are

\[
 (s_\pm,r_\pm)=
 \bigl(a(u)\pm L\cos\alpha(u),
       h(u)\mp L\sin\alpha(u)\bigr).                        \tag{2.9}
\]

Since both lie in \(D\) and \(g\ge0\),

\[
                 L|\sin\alpha(u)|\le h(u)\le {L\over8}.     \tag{2.10}
\]

Thus \(\cos\alpha(u)\ge\sqrt{63}/8\). Combining this with (2.8),

\[
              s_+>{L\over2},\qquad s_-<-{L\over2}.          \tag{2.11}
\]

Membership in \(D\) and (2.10) also give

\[
             g(s_\pm)\le r_\pm
             \le h(u)+L|\sin\alpha(u)|\le2h(u)\le2u.        \tag{2.12}
\]

Convexity makes \(g\) nondecreasing away from its minimum at zero.
Therefore

\[
                       g(\pm L/2)\le2u.                     \tag{2.13}
\]

Apply Proposition 1.1 with \(R=L/2\) and \(a=1/2\) on each side:

\[
 \Theta_\pm(L/4)\le
 {2u\over(1/2)(L/2)}={8u\over L}.
\]

This proves (2.6)--(2.7). QED.

### 2.2 CMC and curvature-energy versions

For a planar CMC interface, \(\alpha'(v)=\lambda\) is constant and

\[
                         K_u=|\lambda|u.                    \tag{2.14}
\]

For a balanced Cheeger interface, \(|\lambda|\le\psi\). If
\(T=\gamma/\psi\) and \(u\le\varepsilon T\), the turning hypothesis in
(2.5) holds whenever

\[
                              \varepsilon\gamma\le {1\over4}.\tag{2.15}
\]

For a non-CMC planar collar, Cauchy--Schwarz gives

\[
 K_u^2\le u\int_0^u|\alpha'(v)|^2\,dv.                     \tag{2.16}
\]

In a higher-dimensional interface, the corresponding estimate along an
inward conormal curve \(\gamma\subset\Sigma\) is

\[
 |N(\gamma(u))-N(\gamma(0))|
 \le\int_0^u|A(\gamma'(v))|\,dv
 \le\sqrt{u}
       \left(\int_0^u|A(\gamma'(v))|^2dv\right)^{1/2}.       \tag{2.17}
\]

This is the precise place where a long-ray curvature-energy estimate
could enter. But (2.17) alone does not put the collar and its normal ray
in the frozen affine plane
\(\operatorname{span}\{N(0),n_K(0)\}\). Proposition 2.1 transfers
literally only for an exact normal section, or after separately controlling
the transverse drift of both the interface and the support. Mean
curvature \(H=\operatorname{tr}A\) does not bound the individual
conormal curvature in dimension \(n>2\).

## 3. Integration against a large contact tensor

Let \(\Gamma_0\) be a contact patch admitting exact planar collars as in
Proposition 2.1. Fix \(c\in\mathbb R^n\), let

\[
 B={1\over p}\int_\Gamma(x-c)\otimes n_K(x)\,d\eta(x),
 \qquad U={B\over\|B\|_F},                                  \tag{3.1}
\]

and define the scalar flux

\[
 b(x)=\left\langle U,(x-c)\otimes n_K(x)\right\rangle_F,
 \qquad b_+(x)=\max(b(x),0).                                \tag{3.2}
\]

Here and below \(B\ne0\); the contact-tensor branch has
\(\|B\|_F\ge\sqrt n/4\).

Then

\[
              \int_\Gamma b_+\,d\eta
              \ge\int_\Gamma b\,d\eta=p\|B\|_F.             \tag{3.3}
\]

Thus, when \(\|B\|_F\ge\sqrt n/4\), the positive flux has size at least
\(p\sqrt n/4\).

For \(x\in\Gamma_0\), suppose a collar point at distance
\(u_x\le\varepsilon T\) has a two-sided surviving normal segment of
length \(4T\). Let \(K_x\) be its collar turning (2.4), and let

\[
                    \Theta_x={\Theta_{x,+}(T)+
                                     \Theta_{x,-}(T)\over2}. \tag{3.4}
\]

### Proposition 3.1 (flux-weighted longitudinal turning)

Assume that \(\Gamma_0\) carries a fraction \(\beta>0\) of the positive
contact flux:

\[
       \int_{\Gamma_0}b_+\,d\eta
       \ge\beta\int_\Gamma b_+\,d\eta.                      \tag{3.5}
\]

If \(K_x\le1/4\) throughout \(\Gamma_0\), then

\[
 { \int_{\Gamma_0}b_+(x)\Theta_x\,d\eta(x)
   \over \int_{\Gamma_0}b_+(x)\,d\eta(x)}
                              \le2\varepsilon.              \tag{3.6}
\]

More generally, define the probability measure

\[
 d\nu={b_+\,d\eta\over\int_{\Gamma_0}b_+\,d\eta}
       \quad\text{on }\Gamma_0.                             \tag{3.7}
\]

If

\[
                              \int K_x^2\,d\nu(x)\le\zeta^2, \tag{3.8}
\]

then

\[
                              \int\Theta_x\,d\nu(x)
                         \le2\varepsilon+8\pi\zeta^2.       \tag{3.9}
\]

#### Proof

Equation (3.6) is Proposition 2.1 integrated against \(b_+\).
For (3.9), Markov's inequality gives
\(\nu\{K_x>1/4\}\le16\zeta^2\). On the complement use (2.7);
on the exceptional set use the trivial convex-graph bound
\(\Theta_x\le\pi/2\). Therefore

\[
 \int\Theta_xd\nu
 \le2\varepsilon+{\pi\over2}(16\zeta^2),
\]

as claimed. QED.

This proves the requested integrated longitudinal estimate. It is
dimension-free and retains the full Frobenius-size contact tensor through
(3.3). Its unproved hypothesis is now explicit: a fixed portion of the
positive tensor flux must admit short collars ending at surviving
normal-ray points. A theorem saying merely that a large fraction of
\(\Sigma\) has surviving rays does not imply (3.5); it may discard the
entire collar of \(\Gamma\).

## 4. The quartic support, quantitatively

Consider (0.3) and its central flat interface

\[
                  \Sigma=\{x_1=0,\ |x'|\le1\},\qquad
                  N=e_1.                                   \tag{4.1}
\]

At a contact point \((0,\omega)\), \(|\omega|=1\), use the exact normal
section

\[
                        (s,r)\longmapsto (s,(1-r)\omega).    \tag{4.2}
\]

The lower support boundary is

\[
 g(s)=1-\sqrt{1-s^4},\qquad |s|\le1.                        \tag{4.3}
\]

It has

\[
\begin{aligned}
 g'(s)&={2s^3\over\sqrt{1-s^4}},\\
 g''(s)&={6s^2\over\sqrt{1-s^4}}
          +{4s^6\over(1-s^4)^{3/2}},\\
 g''(0)&=0.                                                 \tag{4.4}
\end{aligned}
\]

Nevertheless \(g(s)>0\) for every \(s\ne0\), so there is no nontrivial
tangent support segment. The exact reach at depth \(h\in[0,1]\) is

\[
                       R(h)=(2h-h^2)^{1/4},                 \tag{4.5}
\]

and

\[
 {h^{1/4}}\le R(h)\le(2h)^{1/4}.                            \tag{4.6}
\]

The interface point \((0,(1-u)\omega)\) has a surviving normal segment
\(\{(t,(1-u)\omega):|t|\le L\}\) exactly when

\[
 L^4+(1-u)^2\le1
 \quad\Longleftrightarrow\quad
 u\ge u_L:=1-\sqrt{1-L^4}.                                 \tag{4.7}
\]

For \(0\le L\le1\),

\[
                         {L^4\over2}\le u_L\le L^4.          \tag{4.8}
\]

The support turning through \(s\) is

\[
          \Theta(s)=\arctan{2s^3\over\sqrt{1-s^4}}
                    =2s^3+O(s^7).                          \tag{4.9}
\]

Thus the scale \(u_L/L\asymp L^3\) in Proposition 2.1 has the same cubic
order as the true support turning. Pointwise contact curvature
\(g''(0)=0\) sees none of this.

### 4.1 How much of the interface survives

The relative \((n-1)\)-area of the points of \(\Sigma\) whose normal
segment of half-length \(L\) survives is exactly

\[
                  \frac{|\{|x'|\le\sqrt{1-L^4}\}|}
                       {|B_2^{\,n-1}|}
                  =(1-L^4)^{(n-1)/2}.                       \tag{4.10}
\]

In the original coordinates a fixed \(L>0\) therefore kills almost the
entire interface as \(n\to\infty\). The relevant comparison, however, is
after isotropic normalization.

Let

\[
 \sigma_1^2=\mathbb E_{K_n}X_1^2,\qquad
 \sigma_\perp^2=\mathbb E_{K_n}X_j^2\quad(j\ge2).
\]

Direct beta integration gives

\[
\begin{aligned}
 \sigma_1^2
  &= {B(3/4,(n+1)/2)\over B(1/4,(n+1)/2)},\\
 \sigma_\perp^2&={2\over2n+3}.                              \tag{4.11}
\end{aligned}
\]

In particular,

\[
 \sigma_1^2\sim{\Gamma(3/4)\over\Gamma(1/4)}
                   \left({n+1\over2}\right)^{-1/2},
 \qquad \sigma_\perp\sim n^{-1/2}.                          \tag{4.12}
\]

Write \(L_{\rm iso}=L/\sigma_1\) and
\(u_{\rm iso}=u/\sigma_\perp\). For every fixed bounded
\(L_{\rm iso}\),

\[
 u_{{\rm iso},L}
 ={1-\sqrt{1-\sigma_1^4L_{\rm iso}^4}\over\sigma_\perp}
 \asymp {L_{\rm iso}^4\over\sqrt n},                        \tag{4.13}
\]

with universal comparison constants for all sufficiently large \(n\).
These constants can be made explicit. Put

\[
                         c_*={\Gamma(3/4)\over\Gamma(1/4)}.
\]

Wendel's gamma-ratio inequality applied with
\(x=(2n+3)/4\) and increment \(1/2\) gives, for \(n\ge3\),

\[
 {c_*^2\over\sqrt n}
 \le{\sigma_1^4\over\sigma_\perp}
 \le {4c_*^2\over\sqrt n}.                                 \tag{4.14}
\]

Since \(z/2\le1-\sqrt{1-z}\le z\) for \(0\le z\le1\), whenever
\(\sigma_1^4L_{\rm iso}^4\le1/2\),

\[
 {c_*^2L_{\rm iso}^4\over2\sqrt n}
 \le u_{{\rm iso},L}
 \le {4c_*^2L_{\rm iso}^4\over\sqrt n}.                    \tag{4.15}
\]

Moreover,

\[
 (1-\sigma_1^4L_{\rm iso}^4)^{(n-1)/2}
 \longrightarrow
 \exp\left[
 -\left({\Gamma(3/4)\over\Gamma(1/4)}\right)^2
       L_{\rm iso}^4\right].                                \tag{4.16}
\]

The isotropic support graph satisfies the explicit turning bound

\[
 \Theta_{\rm iso}(L_{\rm iso})
 \le {8\sqrt2\,c_*^2L_{\rm iso}^3\over\sqrt n}              \tag{4.17}
\]

under the same condition
\(\sigma_1^4L_{\rm iso}^4\le1/2\). Indeed, the derivative of the
isotropically rescaled graph is

\[
 {2\sigma_1^4L_{\rm iso}^3\over
   \sigma_\perp\sqrt{1-\sigma_1^4L_{\rm iso}^4}},
\]

and turning is at most slope. Thus a fixed isotropic tube survives on a
nonzero surface fraction in a
collar of depth \(O(n^{-1/2})\), and the associated support strip is
\(O(n^{-1/2})\)-flat, but it is never exactly ruled.

Finally, (4.1) is a full affine slice. The tensor Minkowski identity has
\(H=0\) and \(Q_N=e_1\otimes e_1\), hence

\[
 B=I-e_1\otimes e_1,\qquad
                       \|B\|_F=\sqrt{n-1}.                  \tag{4.18}
\]

This is the largest possible contact-tensor regime. The example shows
that a contact tensor and retained interior tubes can coexist while the
contact points themselves have zero reach. What rescues the example is
the independently visible full slice, not a support ruling inferred from
\(\mathrm{II}_{\partial K}(N,N)=0\).

## 5. Chord replacement on an exact smooth minimizer

The unit disk \(D=B_2^2\) gives an exact test. The classical relative
isoperimetric theorem in a Euclidean disk states that, for every prescribed
area \(v\in(0,|D|/2)\), a minimizer is the lens cut off by a circular arc
which meets \(\partial D\) orthogonally. At \(v=|D|/2\) the limiting arc
is a diameter. This follows from existence and planar regularity, the
constant-curvature and orthogonal-contact Euler equations, and circular
symmetrization; among the resulting candidates the boundary cap beats an
interior circle. Thus the nonflat caps are smooth exact global relative
perimeter minimizers, not merely stable critical points.

One precise source is A. Ros and E. Vergasta, *Stability for
hypersurfaces of constant mean curvature with free boundary*,
Geom. Dedicata **56** (1995), 19--33, DOI
10.1007/BF01263611: its ball partitioning theorem contains this planar
classification. The hypotheses used here are the unit Euclidean disk,
ordinary relative perimeter, fixed Lebesgue area, and embedded separating
curves.

Let a portion of such an interface lie on a circle of radius \(\rho\),
and let its central angle be \(\theta\in(0,\pi)\). Its length, chord
length, and the area between arc and chord are

\[
\begin{aligned}
 L_{\rm arc}&=\rho\theta,\\
 L_{\rm chord}&=2\rho\sin(\theta/2),\\
 A_{\rm seg}&={\rho^2\over2}(\theta-\sin\theta).             \tag{5.1}
\end{aligned}
\]

The CMC multiplier is \(\lambda=1/\rho\). Replacing the arc by its chord
loses \(A_{\rm seg}\) of volume. An optimal smooth volume repair on a
remote portion of the same CMC interface costs
\(\lambda A_{\rm seg}+O(A_{\rm seg}^2)\) in perimeter. The leading
constrained cost is therefore

\[
\begin{aligned}
 \Delta_\lambda(\theta)
 &=L_{\rm chord}-L_{\rm arc}+\lambda A_{\rm seg}\\
 &=\rho\left[
       2\sin{\theta\over2}-{\theta+\sin\theta\over2}\right].
                                                               \tag{5.2}
\end{aligned}
\]

Put

\[
 F(\theta)=2\sin(\theta/2)-{\theta+\sin\theta\over2}.        \tag{5.3}
\]

Then

\[
 F(0)=0,\qquad
 F'(\theta)=\cos(\theta/2)\,[1-\cos(\theta/2)]>0
                 \quad(0<\theta<\pi).                      \tag{5.4}
\]

Hence \(\Delta_\lambda(\theta)>0\) for every nontrivial subarc, and

\[
                    \Delta_\lambda(\theta)
                    ={\rho\theta^3\over24}+O(\rho\theta^5). \tag{5.5}
\]

This calculation also applies to a subarc beginning at the free-boundary
contact point: the chord is an admissible finite-perimeter competitor even
though it does not satisfy the stationary orthogonality condition.
For a short subarc, \(A_{\rm seg}=O(\rho^2\theta^3)\), so the
\(O(A_{\rm seg}^2)\) repair error is of strictly higher order than the
positive term in (5.5).

Therefore exact global minimality cannot be invoked to reverse the sign
of the CMC cancellation. A first-order saving is available at a genuine
ridge because the original surface is not stationary there. It is not
available on a smooth CMC contact strip, even when tangent completion
fails.

## 6. What remains

The completed implications are

\[
\begin{gathered}
 \text{contact flux paired to a depth-}\varepsilon T
 \text{ collar with surviving }4T\text{ rays}\\
 \Longrightarrow
 \text{flux-weighted support turning at scale }T
 \le2\varepsilon+8\pi\zeta^2,                              \tag{6.1}\\
 \text{zero global interface curvature}
 \Longrightarrow
 \text{full affine slices}
 \Longrightarrow
 \text{central-cell completion}.                           \tag{6.2}
\end{gathered}
\]

The invalid implications are

\[
\begin{gathered}
 \mathrm{II}_{\partial K}(N,N)=0\text{ on }\Gamma
 \ \not\Longrightarrow\
 \text{support ruled along }N,                              \tag{6.3}\\
 \text{failure of tangent completion on a smooth CMC strip}
 \ \not\Longrightarrow\
 \text{chord saving after volume repair}.                   \tag{6.4}
\end{gathered}
\]

A viable contact-tensor inverse now needs two genuinely global inputs:

1. a **contact-to-collar incidence theorem** showing that a universal
   fraction of the positive tensor flux in (3.3) reaches good normal rays
   at depth \(o(T)\), with controlled collar curvature energy and
   transverse drift; and
2. a **multi-cell reassignment** which cancels volume changes between
   cells before taking a local limit. Scalar repair of each chord
   separately has the wrong sign by (5.2).

The quartic calculation shows the correct form of the first statement:
it must permit approximately flat strips and a vanishing collar depth,
not demand exact support segments. No proof of either global input is
obtained here.
