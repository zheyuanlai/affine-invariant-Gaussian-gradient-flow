# Balanced swept tubes: an exact completion criterion and two stress tests

## 0. Verdict

Balance removes the rare-tail counterexample from
`high_rank_product_inverse.md`, but the aggregate balanced-tube data do not
by themselves give a dimension-free inverse.  There are three precise
conclusions.

1. For the isotropic product of one-sided exponentials, the balanced
   coordinate-maximum leaf has

   \[
       p_m={m\over2}\bigl(2^{1/m}-1\bigr)
       \in\left[{\log2\over2},{1\over2}\right].       \tag{0.1}
   \]

   Its normal matrix is exactly \((p_m/m)I_m\), every regular facet has
   zero Jacobi charge, and two-sided tubes of length \(h/p_m\) lose only
   \(O(h)\) of their flux.  Thus the balanced model has a universal
   perimeter; it is not a counterexample.

2. The cyclically constrained exponential example behaves the same way at
   its median maximum level.  With the parameter \(a=m^{-6}\), its balanced
   perimeter lies between `.33` and `.44`, its normal rank is \(m\), and its
   killed-flux loss before distance \(h/p\) is at most \(5h+O(m^{-11})\).
   After whitening this remains true for a universally elliptic constant
   anisotropy.  Hence balance genuinely excludes the vanishing-flux tail
   used in the earlier obstruction, even though it does not remove the
   cyclic support coupling.

3. There is an exact dimension-free closure if the planar patches can be
   completed to full hyperplane slices at cost \(O(p)\).  If \(a_i\) is the
   area of a planar patch, \(s_i\) is the area of the complete slice, and
   the corresponding smaller halfspaces have total mass at least
   \(\alpha\), then

   \[
      \sum_i(s_i-a_i)\le\eta p
      \quad\Longrightarrow\quad
      p\ge {\phi_-\alpha\over2\sqrt3(1+\eta)}.       \tag{0.2}
   \]

   Here \(\phi_-\) is the lower ellipticity constant of the anisotropy.
   Balanced swept mass supplies \(\alpha\simeq h\).  What it does **not**
   supply is the slice-completion estimate.  In the product example the
   finite killed loss tends to zero with \(h\), while the completion defect
   tends to \(p\).  Therefore no estimate of completion by the finite killed
   loss can be valid.

The missing theorem is consequently not a covariance estimate.  It is a
global completion/overlap theorem saying that persistent flat Wulff patches
either occupy a controlled part of their full marginal slices, or their
missing parts produce a bounded-reuse physical perimeter saving.  Without
that theorem, the scalar swept-mass statement merely constructs a bad
one-Lipschitz signed-distance function and is a reformulation of the desired
first-moment KLS bound.

## 1. What balanced swept mass says by itself

Let \(\mu\) be an isotropic log-concave probability on \(\mathbb R^n\).  Let
\(A\) be a finite-perimeter set of mass \(1/2\), and let \(d_A\) denote its
signed Euclidean distance, negative in \(A\) and positive in its complement.
The function \(d_A\) is one-Lipschitz.  Suppose that, for some \(0<h<1/2\)
and \(T>0\),

\[
 \mu\{-T<d_A<0\}\le h,
 \qquad
 \mu\{0<d_A<T\}\le h.                               \tag{1.1}
\]

The inequalities are written in the direction useful below; equality up to
\(O(h^2)\) is what the killed-tube calculation gives.  For every
\(0<t<T\),

\[
 \mu\{d_A\le-t\}\ge {1\over2}-h,
 \qquad
 \mu\{d_A\ge t\}\ge {1\over2}-h.                    \tag{1.2}
\]

Layer cake therefore gives the exact lower bound

\[
 \boxed{
 E_\mu|d_A|\ge(1-2h)T.}                              \tag{1.3}
\]

In the proposed saturated band, \(T=h/p\), so

\[
                    E_\mu|d_A|\ge{(1-2h)h\over p}.   \tag{1.4}
\]

Thus a hypothetical sequence \(p\to0\) already produces the bad
one-Lipschitz functions required by target (T3).  A dimension-free upper
bound for the left side of (1.4) is KLS itself.  Balance alone has not added
an independent inequality.

### 1.1 The strongest direct covariance estimate

The same point is visible at second order.  Let \(\Gamma\) carry base area
measure \(d\sigma\), let \(z(x)\) be the Wulff displacement, and suppose the
two-sided map

\[
                         F(x,t)=x+t z(x),\qquad |t|\le T,       \tag{1.5}
\]

is injective and its weighted Jacobian is at least \(\rho_0>0\).  Put

\[
                     M=\int_\Gamma z\otimes z\,d\sigma.        \tag{1.6}
\]

For every unit vector \(\theta\), isotropy and symmetry in \(t\) give

\[
\begin{aligned}
 1
 &=\int\langle x,\theta\rangle^2d\mu(x)\\
 &\ge\rho_0\int_\Gamma\int_{-T}^T
       \langle x+t z(x),\theta\rangle^2dt\,d\sigma(x)\\
 &\ge {2\rho_0T^3\over3}\,\theta^TM\theta.
                                                               \tag{1.7}
\end{aligned}
\]

Consequently

\[
                  {2\rho_0T^3\over3}M\preceq I.               \tag{1.8}
\]

If \(\operatorname {tr}M=p\), \(Q=M/p\), and \(T=h/p\), then

\[
 \boxed{
 p\ge\sqrt{2\rho_0\over3}\,h^{3/2}
                       \sqrt{\|Q\|_{\rm op}}.}                 \tag{1.9}
\]

For \(Q=I_k/k\), this is only \(p\ge c h^{3/2}/\sqrt{k}\).  High normal
rank weakens, rather than strengthens, the direct operator bound.

This loss is algebraically sharp.  Give each of the \(k\) rays
\(\mathbb Re_i\) base mass \(p/k\), put uniform tube density on
\([-T,T]\), and place the remaining mass at the origin.  The tube covariance
in direction \(e_i\) is

\[
                 {p\over k}\int_{-T}^Tt^2dt
                 ={2h^3\over3kp^2}.                            \tag{1.10}
\]

Choosing

\[
                         p=\sqrt{2\over3}{h^{3/2}\over\sqrt k} \tag{1.11}
\]

makes that covariance exactly one in every direction.  This ray-current
model is not log-concave and is not offered as a KLS counterexample.  It
proves that isotropy, swept mass, tube length, and normal rank contain no
stronger moment inequality.  Any improvement must use global log-concave
overlap or support/Hessian coupling.

## 2. A full-slice completion theorem

The exact one-dimensional input is recorded first.

**Lemma 2.1 (isotropic one-dimensional halfspaces).**  Let \(\nu\) be an
isotropic log-concave probability on \(\mathbb R\).  If \(H\) is a halfline
with \(\nu(H)\le1/2\), then

\[
                         \nu^+(H)\ge {1\over2\sqrt3}\nu(H).    \tag{2.1}
\]

**Proof.**  Write \(f\) for the density and
\(I(t)=f(F^{-1}(t))\).  One-dimensional log-concavity makes \(I\) concave
on \([0,1]\), with zero endpoint values.  If \(M=\|f\|_\infty\), concavity
gives \(I(1/2)\ge M/2\).  Among probability densities bounded by \(M\),
symmetric decreasing rearrangement and the bathtub principle show that the
variance is at least \(1/(12M^2)\).  Since the variance here is one,
\(M\ge1/\sqrt{12}\).  Concavity once more gives, for \(t\le1/2\),

\[
 {I(t)\over t}\ge {I(1/2)\over1/2}=2I(1/2)
       \ge M\ge {1\over\sqrt{12}}>{1\over2\sqrt3}.             \tag{2.2}
\]

The deliberately weaker last constant is (2.1).  QED.

The strict inequality in (2.2) is only typographical slack: in fact
\(1/\sqrt{12}=1/(2\sqrt3)\).

Let \(\Phi\) be a constant anisotropy satisfying

\[
                             \Phi(u)\ge\phi_->0
       \quad\hbox{for every unit }u.                           \tag{2.3}
\]

Suppose a flat boundary packet is partitioned into patches \(\Gamma_i\),
where \(\Gamma_i\) lies in

\[
                     \Pi_i=\{x:\langle x,u_i\rangle=b_i\}.    \tag{2.4}
\]

Let

\[
 a_i=\int_{\Gamma_i}\Phi(u_i)d\sigma_\mu,
 \qquad
 s_i=\Phi(u_i)P_\mu\{\langle x,u_i\rangle\ge b_i\},           \tag{2.5}
\]

where the orientation is chosen so that

\[
 q_i:=\mu\{\langle x,u_i\rangle\ge b_i\}\le {1\over2}.       \tag{2.6}
\]

Here the second expression in (2.5) means the anisotropic area of the full
hyperplane slice, not the perimeter of the halfspace multiplied twice.
Since \(\Gamma_i\subset\Pi_i\), one has \(0\le a_i\le s_i\).

**Proposition 2.2 (slice completion closes the flat branch).**  Assume

\[
 p=\sum_i a_i,
 \qquad \sum_iq_i\ge\alpha,
 \qquad R_{\rm comp}:=\sum_i(s_i-a_i)\le\eta p.                \tag{2.7}
\]

Then

\[
 \boxed{
                         p\ge {\phi_-\alpha
                                  \over2\sqrt3(1+\eta)}.}      \tag{2.8}
\]

**Proof.**  The marginal \(\langle X,u_i\rangle\) is isotropic and
log-concave.  Lemma 2.1 and (2.3) give

\[
                              s_i\ge{phi_-\over2\sqrt3}q_i.
                                                                    \tag{2.9}
\]

Summing and using (2.7),

\[
              (1+\eta)p\ge p+R_{\rm comp}
                =\sum_i s_i
                \ge{\phi_-\alpha\over2\sqrt3}.
\]

QED.

If a Wulff tube based on \(\Gamma_i\) sweeps mass \(w_i\) into the chosen
halfspace, then \(q_i\ge w_i\).  Thus a balanced band with total swept mass
\(\sum_iw_i\ge h(1-\varepsilon)\) supplies

\[
                             \alpha=h(1-\varepsilon).           \tag{2.10}
\]

For fixed \(h\), Proposition 2.2 is dimension free.  It uses neither normal
rank nor thin shell.  The load-bearing hypothesis is the global completion
bound in (2.7).

## 3. Exact balanced product-exponential calculation

Let \(Y_1,\ldots,Y_m\) be independent unit-rate exponentials and put
\(X_i=Y_i-1\).  The law of \(X\) is isotropic and log-concave.  For \(L>0\)
set

\[
                         A_L=\{\max_iY_i\ge L\}.                \tag{3.1}
\]

Write \(q=e^{-L}\) and \(d=1-q\).  Then

\[
 \mu(A_L)=1-d^m,
 \qquad
 P_\mu(A_L)=mqd^{m-1}.                              \tag{3.2}
\]

At the balanced level \(d^m=1/2\),

\[
 p_m=m(1-d)d^{m-1}
     ={m\over2}(d^{-1}-1)
     ={m\over2}(2^{1/m}-1).                         \tag{3.3}
\]

The lower bound in (0.1) follows from \(e^x-1\ge x\).  The upper bound
follows from \(2^{1/m}\le1+1/m\), which is equivalent to
\((1+1/m)^m\ge2\).

Up to null tie sets, the boundary is the disjoint union

\[
 \Gamma_i=\{Y_i=L,\ Y_j<L\ (j\ne i)\}.                         \tag{3.4}
\]

Every component has area \(p_m/m\), is flat, and has normal \(-e_i\).
The potential \(V(y)=\sum_i y_i\) has zero Hessian.  All regular Jacobi
terms therefore vanish, and

\[
                         M={p_m\over m}I_m.                     \tag{3.5}
\]

Moving from \(\Gamma_i\) into the central box by distance \(t\) survives
until no other coordinate reaches the new threshold \(L-t\).  Hence

\[
 {R_-(t)\over p_m}
   =\left({1-qe^t\over1-q}\right)^{m-1}.                       \tag{3.6}
\]

In the other direction there is no collision and the exponential density
gives

\[
                         {R_+(t)\over p_m}=e^{-t}.              \tag{3.7}
\]

Put

\[
                         b_m=(m-1)(2^{1/m}-1).                  \tag{3.8}
\]

For \(m\ge2\),

\[
 b_m\le (m-1){\log2\over m}2^{1/m}\le\log2;                  \tag{3.9}
\]

the last inequality follows from
\((1-1/m)2^{1/m}\le1\).  Bernoulli's inequality now yields, for
\(0\le t\le1\),

\[
 1-{R_-(t)\over p_m}
 \le b_m(e^t-1)
 \le(\log2)(e-1)t<1.2t,                                      \tag{3.10}
\]

while

\[
                         1-{R_+(t)\over p_m}\le t.             \tag{3.11}
\]

Take \(0<h\le1/10\) and \(T=h/p_m\).  By (0.1),

\[
                         T\le {2h\over\log2}<.29.              \tag{3.12}
\]

Therefore the two normalized flux losses before time \(T\) are at most
\(3.5h\) and \(2.9h\), respectively.  Moreover

\[
\begin{aligned}
 \int_0^TR_-(t)dt&\in[h(1-3.5h),h],\\
 \int_0^TR_+(t)dt&=p_m(1-e^{-T})
              \in[h(1-T/2),h].                                \tag{3.13}
\end{aligned}
\]

This is the requested two-sided global saturation with explicit constants.

The slice-completion quantities are equally explicit.  The full slice and
smaller halfspace mass associated with coordinate \(i\) are both \(q\),
whereas \(a_i=qd^{m-1}\).  Hence

\[
 \sum_iq_i=mq\in[1/2,\log2],
 \qquad
 {R_{\rm comp}\over p_m}
 ={mq-p_m\over p_m}=d^{-(m-1)}-1=2d-1\le1.                    \tag{3.14}
\]

Proposition 2.2 therefore gives a universal lower bound directly.  Notice
the important mismatch

\[
 {\hbox{finite killed loss}\over p_m}=O(h)
 \quad\hbox{but}\quad
 {R_{\rm comp}\over p_m}\longrightarrow1.                    \tag{3.15}
\]

as first \(m\to\infty\) and then \(h\downarrow0\).  Thus finite killed loss
cannot dominate completion defect.

The Bobkov--Houdre tensorization theorem also gives a universal lower bound
for the Cheeger constant of this product, while the balanced set above gives
\(\psi\le2p_m\).  Hence its perimeter is comparable, with universal
constants, to the balanced Cheeger scale.  No assertion that the maximum
set is an exact isoperimetric minimizer is needed here.

## 4. The cyclically constrained family at balance

Fix \(m\ge4\), put \(a=m^{-6}\), and define

\[
 \Omega_{m,a}=\{y\in[0,\infty)^m:y_i+y_{i+1}\ge a
                       \ (i\bmod m)\}.                          \tag{4.1}
\]

Let \(\mu_{m,a}\) have density proportional to
\(e^{-\sum_i y_i}\) on \(\Omega_{m,a}\).  Equivalently, it is the product
exponential law \(\pi_m\) conditioned on the event \(E=\Omega_{m,a}\).
Since

\[
 \pi_m\{Y_i+Y_{i+1}<a\}=1-e^{-a}(1+a)\le {a^2\over2},          \tag{4.2}
\]

the removed mass satisfies

\[
                         q_a:=\pi_m(E^c)\le {1\over2m^{11}}.    \tag{4.3}
\]

Consequently, for every Borel set \(B\),

\[
          |\mu_{m,a}(B)-\pi_m(B)|\le{q_a\over1-q_a}.           \tag{4.4}
\]

Let \(L_a\) be the unique median of \(\max_iY_i\) under \(\mu_{m,a}\), and
put \(d=1-e^{-L_a}\).  Equations (4.3)--(4.4) imply

\[
                         d^m\in[.49,.51].                       \tag{4.5}
\]

For \(L>a\), removing the coordinate fixed at \(L\) leaves a path of
\(m-1\) coordinates.  If

\[
 J_{m-1}(L)=\int_{[0,L)^{m-1}}
       e^{-\sum_jy_j}\mathbf1_{\{y_j+y_{j+1}\ge a\}}dy,       \tag{4.6}
\]

then

\[
 p_a(L)={m e^{-L}J_{m-1}(L)\over\pi_m(E)}.                     \tag{4.7}
\]

The union bound in the path gives

\[
 0\le d^{m-1}-J_{m-1}(L)\le{(m-2)a^2\over2}.                  \tag{4.8}
\]

For \(s=d^m\in[.49,.51]\), the unconstrained expression is

\[
 p_0(L)=m(1-d)d^{m-1}
       =ms(s^{-1/m}-1).                                        \tag{4.9}
\]

Using \(e^x-1\in[x,xe^x]\) and \(m\ge4\) gives

\[
                         .343<p_0(L)<.42.                       \tag{4.10}
\]

Since \(d^{m-1}\ge.49\), (4.3), (4.7), and (4.8) imply, with ample slack,

\[
 \boxed{
                         .33<p_a(L_a)<.44.}                     \tag{4.11}
\]

Cyclic symmetry again gives equal component areas and

\[
                         M(L_a)={p_a(L_a)\over m}I_m.           \tag{4.12}
\]

The regular components are planar, the potential is affine, and every
active support face is tangent to the corresponding normal translation.
Thus all smooth Jacobi and contact-curvature terms vanish exactly.

For the inward tube the exact survival ratio is

\[
 {R_{a,-}(t)\over p_a(L_a)}
             ={J_{m-1}(L_a-t)\over J_{m-1}(L_a)}.              \tag{4.13}
\]

For \(h\le1/10\), (4.11) gives \(T=h/p_a<.304\).  Throughout
\([0,T]\), the unconstrained truncated probability
\((1-e^{-(L_a-t)})^{m-1}\) is at least `.33`.  Applying (4.8) at both
endpoints and then the calculation (3.8)--(3.10), now with
\(d^m\in[.49,.51]\), yields

\[
 1-{R_{a,-}(t)\over p_a(L_a)}
       \le1.4t+4m^{-11}.                                      \tag{4.14}
\]

The outward ratio is still exactly \(e^{-t}\).  Hence at \(t=T\),

\[
 \boxed{
 \max\left\{1-{R_{a,-}(T)\over p_a},
              1-{R_{a,+}(T)\over p_a}\right\}
       \le5h+4m^{-11}.}                                       \tag{4.15}
\]

Both swept masses equal \(h(1-O(h)-O(m^{-11}))\).  This is the balanced
counterpart of the tail calculation in `high_rank_product_inverse.md`.

Finally, conditioning changes the covariance by less than one percent, as
proved there from the fourth moment of product exponentials:

\[
                         .99I\preceq\operatorname {Cov}(\mu_{m,a})
                                      \preceq1.01I.             \tag{4.16}
\]

After centering and whitening, the measure is isotropic and Euclidean
perimeter becomes the constant anisotropy

\[
                         \Phi(n)=|A^{-1/2}n|,                   \tag{4.17}
\]

whose ellipticity ratio is below `1.011`.  Wulff translation is exactly the
affine image of the rays used above, so (4.11)--(4.15) are unchanged in the
anisotropic normalization.  The support remains affinely product
irreducible.

The cyclic test therefore has a clear verdict: its rare tail leaves refute
a local high-rank-to-product claim, but its balanced leaf has universal
perimeter.  What has not been proved is that this particular balanced leaf
is an exact or asymptotically exact Euclidean Cheeger minimizer.  The test
does not use such an assertion.

## 5. Exact location of the remaining inverse

The balanced assumptions now split into three logically different levels.

* **Scalar saturation.**  Two-sided swept mass \(h(1-O(h))\) gives the
  signed-distance lower bound (1.4).  This is conjecture-strength, not a
  contradiction.
* **Moment saturation.**  Isotropy and a high-rank normal matrix give only
  (1.9), sharply witnessed by the ray-current model.
* **Global slice completion.**  If persistent flat patches cover their full
  marginal slices up to \(O(p)\), Proposition 2.2 gives a dimension-free
  lower bound immediately.

Finite killed flux distinguishes a balanced leaf from a vanishing tail,
but it does not control the third item.  Missing parts of a full slice may
already lie in another phase, may be hidden behind a tangent support face,
or may be coupled through a Hessian direction outside the swept band.  The
product calculation (3.15) rules out the simplest attempted estimate

\[
        R_{\rm comp}\le C\int_0^{h/p}(p-R(t))dt.                \tag{5.1}
\]

Indeed the right side is \(O(h^2)\), while the left side tends to the
nonzero constant \(p_m\).

A sufficient new lemma is therefore the following precise alternative.

> **Completion/overlap alternative.**  For the balanced near-isoperimetric
> Wulff foliation produced by joint equimeasurable replacement, either
> \(R_{\rm comp}\le C p\) on a packet sweeping fixed mass, or the missing
> full-slice pieces generate a physical competitor whose perimeter saving
> is at least \(c p\), with bounded reuse across levels.

The first branch closes by Proposition 2.2.  The second branch contradicts
the integrated profile deficit.  Neither scalar killed flux, covariance,
normal rank, nor one-dimensional marginal log-concavity proves this
alternative.  Establishing it requires the global incidence of phases with
all support and Hessian coupling directions.  This is exactly the datum
which the cyclic example makes invisible to a finite local tube.
