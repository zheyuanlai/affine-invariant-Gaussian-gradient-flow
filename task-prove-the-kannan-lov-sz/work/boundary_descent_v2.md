# Boundary descent, fiber repair, and the conditional-velocity obstruction

This note investigates the proposed induction from a half-mass isoperimetric
boundary in dimension `n` to a cut in the orthogonal marginal in dimension
`n-1`.  It records two estimates which are genuinely dimension free, followed
by the exact term which prevents them from giving a recurrence by themselves.
All perimeter statements below are first proved for smooth densities and sets
of finite perimeter and then pass to the usual BV approximation.  The
orthogonal marginal of an isotropic log-concave law is itself isotropic and
log-concave.

## 1. A one-dimensional oriented repair lemma

Let `eta` be a one-dimensional log-concave probability with density `f`,
variance `sigma^2`, and distribution function `F`.  If `E` has finite
`eta`-perimeter, call a boundary point *wrongly oriented* when the indicator
of `E` jumps from zero to one there.  Put

\[
 W_-(E)=\sum_{x:\,0\to1} f(x).
\]

Then there is a left ray `R=(-infinity,t]` (the empty and full limiting rays
are allowed) such that

\[
 \eta(E\mathbin\triangle R)\le 2\sqrt3\,\sigma W_-(E).       \tag{1}
\]

Here is a proof with constants.  Map `eta` to Lebesgue measure on `(0,1)` by
`z=F(x)`.  For a finite union `U` of intervals, there is a prefix `[0,a]`
such that

\[
 |U\mathbin\triangle[0,a]|
 \le \sum_{z:\,0\to1}\min(z,1-z).                           \tag{2}
\]

To see this, inspect the component or complementary component containing
`1/2`.  If it is a component of `U`, take `a` to be its right endpoint; if it
is a complementary component, take `a` to be its left endpoint.  Every
omitted gap to the left ends at a wrong jump `z<=1/2` and has length at most
`z`; every retained component to the right begins at a wrong jump `z>=1/2`
and has length at most `1-z`.  Summing proves (2).  Countable finite-perimeter
sets follow by monotone truncation of their endpoints.

The one-dimensional Cheeger formula and the elementary median-density bound
give

\[
 \psi_\eta=2f(m)\ge {1\over2\sqrt3\,\sigma}.
\]

Applying the Cheeger inequality to each half-line gives

\[
 f(x)\ge\psi_\eta\min(F(x),1-F(x)).
\]

Substitution in (2) proves (1).

## 2. Surface-normal alignment really gives a subgraph, not a halfspace

Let `mu` be isotropic and log-concave on `R x R^(n-1)`, write points as
`(s,y)`, and disintegrate

\[
 d\mu(s,y)=d\mu_y(s)\,d\nu(y).
\]

Let `sigma_y^2=Var_{mu_y}(s)`.  The law of total variance gives

\[
 \int\sigma_y^2d\nu(y)\le Var_\mu(s)=1.                    \tag{3}
\]

Let `A` be a finite-perimeter set, let `P=mu^+(A)`, and orient the first
coordinate by a unit vector `u`.  Define the negatively oriented horizontal
perimeter

\[
 P_-:=\int_{\partial^*A}(-N\cdot u)_+\,d\sigma_\mu.
\]

Slicing BV sets by the lines parallel to `u` and applying (1) on every fiber
gives a measurable threshold `g(y)` and the subgraph

\[
 G=\{(s,y):s\le g(y)\}
\]

such that, for every `R>0`,

\[
 \begin{split}
 \mu(A\mathbin\triangle G)
 &\le 2\sqrt3 R P_-+\nu\{\sigma_y>R\}\\
 &\le2\sqrt3 R P_-+R^{-2}.
 \end{split}                                                \tag{4}
\]

Optimizing at `R=(sqrt(3)P_-)^(-1/3)` yields

\[
 \boxed{\ \mu(A\mathbin\triangle G)
       \le3^{4/3}P_-^{2/3}.\ }                              \tag{5}
\]

The measurability of `g` can be obtained by minimizing the fiberwise
symmetric-difference functional over rational thresholds and taking the
least approximate minimizer.  The BV slicing theorem identifies the integral
of the wrong endpoint densities with `P_-`.

If

\[
 \delta={1\over P}\int_{\partial^*A}|N-u|^2d\sigma_\mu,
\]

then, since `(-a)_+ <= (1-a)/2` for `a in [-1,1]`,

\[
 P_-\le {1\over2}\int(1-N\cdot u)d\sigma_\mu={P\delta\over4}. \tag{6}
\]

Consequently

\[
 \boxed{\ \mu(A\mathbin\triangle G)
 \le3^{4/3}(P\delta/4)^{2/3}.\ }                            \tag{7}
\]

For a half-mass isoperimetric minimizer, `P` has a universal upper bound by
comparison with a coordinate halfspace and the standard one-dimensional
upper bound on the density of an isotropic log-concave law.  Thus (7) is
`O(delta^(2/3))`.  This is the strongest conclusion justified by normal
alignment alone: it gives closeness to a subgraph, not to a halfspace.  It is
consistent with the decorated product set `A_L` in the boundary-stability
note.

There is an important limitation.  The repair is fiberwise and need not
decrease, or even control, transverse perimeter.  Therefore (7) cannot be
inserted into an induction without a second argument.

## 3. Small horizontal perimeter produces two almost-pure phases

For each fiber put

\[
 a(y)=\mu_y(A_y),\qquad
 H_y=\operatorname{Per}_{\mu_y}(A_y),\qquad H=\int H_y d\nu.
\]

By one-dimensional Cheeger,

\[
 \min(a(y),1-a(y))\le2\sqrt3\,\sigma_yH_y.
\]

The same truncation used in (4), now using `H_y`, gives

\[
 U:=\int\min(a,1-a)d\nu
 \le2\sqrt3 R H+R^{-2}
 \le3^{4/3}H^{2/3}.                                       \tag{8}
\]

Here `H=int |N dot u| dsigma_mu<=P`.  If `mu(A)=1/2` and
`B={y:a(y)>=1/2}`, then

\[
 |\nu(B)-1/2|\le U.                                        \tag{9}
\]

Thus a small-perimeter half-mass cut has two nearly pure fiber phases, and
the phase classifier `B` is balanced.  This conclusion does not use normal
alignment.

If one could prove, in the aligned branch and for `P` below a fixed universal
constant,

\[
 \nu^+(B')\le C P                                           \tag{10}
\]

for some balanced modification `B'` of `B`, the inductive Cheeger bound in
dimension `n-1` would give

\[
 P\ge c\,h_{n-1}.
\]

Together with a universal lower bound in the curvature branch, this would
give the desired fixed-point recurrence

\[
 h_n\ge c\min(1,h_{n-1}).                                  \tag{11}
\]

The rest of the note explains why (10) is not a formal consequence of
slicing.

## 4. Exact derivative identity and the missing conditional velocity

Assume temporarily that `rho(s,y)` is positive and smooth.  Write

\[
 r(y)=\int\rho(s,y)ds,qquad
 f_y(s)=\rho(s,y)/r(y),qquad
 S_y(s)=\nabla_y\log f_y(s).
\]

For an arbitrary finite-perimeter set, the distributional derivative of the
conditional mass satisfies

\[
 \boxed{\quad
 r\,Da
 =r\,\mathbb E_{\mu_y}[(1_{A_y}-a)S_y],dy
   -\pi_\#(N_y\,d\sigma_\mu).\quad}                        \tag{12}
\]

This follows by differentiating
`r(y)a(y)=int 1_A(s,y)rho(s,y)ds`, using
`D_y1_A=-N_y H^(n-1)|partial*A`, and subtracting `a Dr`.

For a smooth subgraph `A={s<=g(y)}`, define

\[
 z=a(y)=F_y(g(y)),\qquad
 v_y(s)=-{\nabla_yF_y(s)\over f_y(s)}.
\]

Then

\[
 \nabla a=f_y(g)(\nabla g-v_y(g)),                          \tag{13}
\]

up to the harmless sign convention in the definition of `v`.  Equivalently,
the physical graph gradient is the product-coordinate gradient plus the
conditional-quantile velocity.  The graph perimeter is

\[
 \int r(y)\sqrt{f_y(g)^2+|\nabla a+f_y(g)v_y(g)|^2}\,dy,    \tag{14}
\]

again after choosing the matching sign convention.

There is a useful exact quantile-coordinate form of both the conditional
velocity and the CMC equation.  Put

\[
 Q(y,z)=F_y^{-1}(z),\qquad I(y,z)=f_y(Q(y,z)),\qquad
 v(y,z)=\nabla_yQ(y,z).
\]

For the boundary `g(y)=Q(y,a(y))`, set

\[
 w=\nabla a+Iv=I\nabla g,\qquad
 D=(I^2+|w|^2)^{1/2}=I(1+|\nabla g|^2)^{1/2}.
\]

Then the score term appearing below has the exact representation

\[
 \boxed{\quad J_u(A)=\int r(y)I(y,a(y))|v(y,a(y))|dy.\quad} \tag{14j}
\]

The two differential identities

\[
 I_z=-V_s(Q,y),
 \qquad \partial_z(Iv)=-S_y(Q)                              \tag{14k}
\]

follow by differentiating `F_y(Q(y,z))=z` and
`Iv=-nabla_y F_y(Q)`.  The perimeter and volume functionals are

\[
 \mathcal P(a)=\int rDdy,
 \qquad \mathcal V(a)=\int ra\,dy.
\]

Their Euler--Lagrange equation is the following exact form of weighted CMC:

\[
 \boxed{\quad
 -{1\over r}\operatorname{div}_y\!\left({rw\over D}\right)
 -{I V_s(Q,y)+w\cdot S_y(Q)\over D}=\lambda.
 \quad}                                                     \tag{14l}
\]

Thus CMC controls the quantile derivative `partial_z(Iv)=-S_y` and the
physical slope `w`, whereas the obstruction `J_u` is the boundary value of
the zeroth-order connection `Iv`.  In particular, on an exactly flat piece
`w=0`, (14l) reduces simply to `V_s(Q,y)=-lambda` and contains no `v` term at
all.  This explains algebraically why stationarity alone does not bound
`J_u`: the regular-simplex flat cut has `w=0`, `V_s=0` in the interior, and
`|v|` of order `n`.  In a smooth convex-body approximation the missing data
are carried by the support boundary layer/free-boundary contact terms.

If the conditional laws have full support and sufficient tail regularity so
that `Iv` tends to zero at `z=0,1`, (14k) also gives the exact tail
representation

\[
 |I v(y,a)|\le
 \min\left\{\int_0^a|S_y(Q(y,z))|dz,
             \int_a^1|S_y(Q(y,z))|dz\right\}.              \tag{14m}
\]

Consequently, with conditional Fisher information

\[
 \mathcal I_{\rm cond}
 =\int r(y)\int_0^1|S_y(Q(y,z))|^2dzdy,
\]

Cauchy--Schwarz and (8) give only

\[
 J_u(A)\le \sqrt{U\mathcal I_{\rm cond}}
 \le 3^{2/3}P^{1/3}\sqrt{\mathcal I_{\rm cond}}.           \tag{14n}
\]

There is no dimension-free upper bound on this Fisher information for
isotropic log-concave laws; it diverges under smooth approximation of a
uniform convex body.  Thus (14n) is not a permissible regularization route.
More importantly, the score in the CMC equation (14l) is multiplied by the
physical slope `w/D`.  In the alignment regime this coefficient tends to
zero, while `J_u` is the norm of `Iv` itself.  First-order stationarity is
therefore least coercive precisely in the regime where (14e) needs control.

The same cancellation persists at second order on a flat piece.  If
`g=c` and it is varied as `g_t=c+t h(y)`, then the induced quantile variation
is `dot a=Ih`, and

\[
 {d\over dt}_{|t=0}(\nabla a_t+I_t v_t)=I\nabla h.
\]

All zeroth-order occurrences of the potentially large velocity cancel.
Accordingly the Jacobi form is the usual physical second variation and does
not control `Iv`.  Any superlinear estimate in (14g) must use global
normalization/phase geometry beyond the CMC equation and local stability.

The score/velocity term in (12)--(14) is precisely the derivative term in a
Knothe conditional transport.  Neither isotropy nor log-concavity gives a
dimension-free pointwise bound on it.  Taking absolute values in (12) and
discarding this term is therefore invalid.

There is nevertheless an exact small-perimeter consequence which identifies
the load-bearing estimate.  Define

\[
 \begin{split}
 J_u(A)&=\int r(y)
   \left|\mathbb E_{\mu_y}[(1_{A_y}-a(y))S_y]\right|dy,\\
 T_u(A)&=\int_{\partial^*A}|N_y|d\sigma_\mu.
 \end{split}                                                \tag{14a}
\]

Formula (12) and the triangle inequality imply

\[
 \operatorname{TV}_\nu(a)\le J_u(A)+T_u(A).                \tag{14b}
\]

Assume `mu(A)=1/2` and let `h_(n-1)` denote the Cheeger lower bound for the
isotropic marginal `nu`.  If the uncertainty `U` in (8) is at most `1/20`,
then for every `t in [1/4,3/4]` the set `{a>t}` has `nu`-mass between `1/4`
and `3/4`.  Indeed, `B={a>=1/2}` has mass in
`[1/2-U,1/2+U]`, and every discrepancy between `B` and `{a>t}` in this range
has mass at most `4U`.  The BV coarea formula therefore gives

\[
 \operatorname{TV}_\nu(a)
 =\int_0^1\nu^+(\{a>t\})dt
 \ge {h_{n-1}\over8}.                                     \tag{14c}
\]

Since `U<=3^(4/3)P^(2/3)`, (14c) applies whenever `P` is below a fixed
numerical constant.  If `delta=P^(-1)int|N-u|^2 dsigma_mu`, then

\[
 T_u(A)\le P\sqrt\delta,                                   \tag{14d}
\]

because `|N_y|^2=1-(N dot u)^2<=|N-u|^2`.  Combining
(14b)--(14d) yields the exact descent inequality

\[
 \boxed{\quad {h_{n-1}\over8}
        \le J_u(A)+P\sqrt\delta.\quad}                     \tag{14e}
\]

In the aligned branch of the stability dichotomy in the boundary-stability
note, `delta<=2sqrt(eta)`, and hence

\[
 {h_{n-1}\over8}\le J_u(A)+\sqrt2 P\eta^{1/4}.             \tag{14f}
\]

This is stronger and more precise than asking directly for the perimeter of
the classifier `B`: all transverse geometry has been reduced to the single
conditional-score flux `J_u(A)`.

It also gives a sharp design criterion for a fixed-point recurrence.  Suppose
one could prove, with universal constants, the two estimates

\[
 \begin{array}{ll}
 \text{curvature branch:}&P\ge c_0\eta^\gamma,\\
 \text{aligned branch:}&J_u(A)\le C_0P^\alpha,
 \end{array}                                                \tag{14g}
\]

for some finite `gamma>0`.  Given a sufficiently small `P`, choose
`eta=(2P/c_0)^(1/gamma)`.  The curvature branch is then impossible, and
(14f) gives

\[
 h_{n-1}\le C\left(P^\alpha+P^{1+1/(4\gamma)}\right).
                                                                    \tag{14h}
\]

Writing `q=min(alpha,1+1/(4gamma))`, this would yield

\[
 P\ge c h_{n-1}^{1/q}.                                    \tag{14i}
\]

The recurrence has an attractive positive fixed point only when `q>1`.
Thus an estimate merely of the form `J=o(1)`, or `J<=CP^alpha` with
`alpha<1`, cannot close the induction.  One needs a genuinely superlinear
small-perimeter estimate `J=o(P)` together with a quantitative curvature
branch.  For `alpha=1`, (14i) is only a constant-loss linear recurrence and
does not prevent decay with the dimension unless its sharp coefficient is at
least one.

## 5. The isotropic regular simplex rules out projected BV contraction

The obstruction is already explicit for a flat halfspace.  Let `v_0,...,v_n`
be unit vectors with

\[
 \langle v_i,v_j\rangle=-1/n\quad(i\ne j),
\]

and consider the regular simplex

\[
 K=\{x:\langle v_i,x\rangle\le r,\ 0\le i\le n\}.
\]

Choose `r` so that the uniform law on `K` is isotropic; regular-simplex
symmetry makes its covariance a scalar before this rescaling.  Put `u=v_0`
and write `x=su+y`, `y perpendicular to u`.  The fiber over `y` is

\[
 [\ell(y),r],\qquad
 \ell(y)=\max_{1\le i\le n}n(\langle v_i,y\rangle-r).       \tag{15}
\]

On every cell where the `i`th affine function is active,

\[
 |\nabla\ell|=n|P_{u^\perp}v_i|=\sqrt{n^2-1}.              \tag{16}
\]

The conditional quantile is

\[
 Q_y(z)=\ell(y)+z(r-\ell(y)),
\]

and hence

\[
 |\nabla_yQ_y(z)|=(1-z)\sqrt{n^2-1}.                       \tag{17}
\]

Take the horizontal median cut `A_t={s<=t}`.  Its physical tangential
perimeter is zero and its total perimeter is the one-dimensional marginal
density at a median, hence is between two universal positive constants.
Nevertheless, on the set of fibers meeting the plane,

\[
 a(y)={t-\ell(y)\over r-\ell(y)},
 \qquad
 |\nabla a|={r-t\over(r-\ell(y))^2}\sqrt{n^2-1}.           \tag{18}
\]

The section `D_t={y:ell(y)<t}` is a regular `(n-1)`-simplex.  This assertion
can be quantified exactly.  Put

\[
 U={nr+t\over(n+1)r}.
\]

The marginal distribution of `S` has distribution function `U^n`, so at a
median `U=2^{-1/n}`.  If `alpha_1,...,alpha_n` are the barycentric
coordinates of a uniform point of `D_t`, then

\[
 t-\ell(y)=n(nr+t)\min_i\alpha_i.                           \tag{18a}
\]

For the uniform Dirichlet `(1,...,1)` law,

\[
 \mathbb P\{\min_i\alpha_i>\epsilon\}=(1-n\epsilon)^{n-1}.
\]

Taking `epsilon=(r-t)/(n(nr+t))` shows that the set
`{t-ell<=r-t}` has probability bounded below by a universal positive
constant (in fact it tends to `1/2`).  On this set one has

\[
 {r-t\over r-\ell(y)}\ge c.
\]

Consequently, using that `r(y)f_y(t)dy` is the surface law of the horizontal
cut and that its conditional distribution on `D_t` is uniform,

\[
 \int|\nabla a|d\nu
 =P(A_t)\,\mathbb E[|v_y(t)|\mid S=t]
 \ge c n P(A_t).                                           \tag{19}
\]

Thus even for an isotropic log-concave measure and an exactly flat median
boundary there is no estimate of the form

\[
 \int|Da|d\nu\le C\,\mu^+(A)                              \tag{20}
\]

with dimension-free `C`.  In (12), the conditional-score term cancels the
large derivative of `a` against the zero physical tangential flux.  This
example does not contradict the desired descent because its horizontal
perimeter is already universal.  It does prove that a valid induction must
contain an explicit dichotomy: either horizontal perimeter is universal, or
the conditional-velocity cancellation must be shown unable to create the
balanced phase transition in (9).

## 6. Product decoration gives the required consistency test

For the set `A_L` in `gamma tensor nu`, with `nu(B)=1/2` and
`Per_nu(B)=p`, put

\[
 w_L=\Phi(L)-\Phi(-L).
\]

Then

\[
 P_L=\varphi(L)+p w_L,
 \qquad
 \delta_L={2p w_L\over P_L}.                               \tag{21}
\]

If `P_L<varphi(1)`, then `L>1`, so `w_L>=w_1>0` and

\[
 P_L\ge w_1 p.                                             \tag{22}
\]

Thus the desired small-perimeter aligned-branch estimate `P>=c h_(n-1)` is
exactly consistent with the decorated bottleneck.  If instead
`varphi(L)=sqrt(p)` and `p` tends to zero, then

\[
 P_L\asymp\sqrt p,
 \qquad \delta_L\asymp\sqrt p,                             \tag{23}
\]

showing again that normal alignment cannot imply halfspace closeness.  Any
proof of (10) must recover the transverse wall `p w_L`, rather than erase it
by graph approximation.

## 7. What the curvature identities do and do not add

Let `A` be a smooth half-mass isoperimetric region with weighted CMC
`lambda`.  Along its outward parallel flow, if `F(t)=mu(A_t)` and `P(t)` is
the weighted perimeter, then at zero

\[
 F'(0)=P,
 \quad P'(0)=\lambda P,
 \quad P''(0)=P(\lambda^2-\mathbb E_\sigma q).              \tag{24}
\]

Global Cheeger minimality gives for small positive and negative `t`

\[
 P(t)\ge2P(1-F(t)),\qquad P(-t)\ge2P F(-t).                \tag{25}
\]

The first derivatives recover only `|lambda|<=2P`.  Formally, the second
derivative in (24) suggests that a large `E_sigma q` should force a perimeter
of order `sqrt(E_sigma q)`.  This inference is not rigorous without a
dimension-free control of the parallel-flow remainder up to time
`t comparable to P/E_sigma q`.  Neither (8), the affine identities, nor the
rotation sums control derivatives of `q`, focal times, or concentration of
curvature.  Pointwise Taylor expansion therefore cannot supply the curvature
branch of (11).

More fundamentally, the normalized stability and affine identities are
homogeneous in the surface measure.  A Gaussian central hyperplane satisfies
all of them with equality, and the same normalized boundary data continue to
satisfy them if the unnormalized surface weight `P` is replaced by an
arbitrary positive number.  Bulk normalization fixes the true Gaussian value,
but that information is absent from the identities.  Product cylinders over
a lower-dimensional stable half-mass boundary inherit the same identities and
the same perimeter as the lower-dimensional factor.  Therefore any curvature
recurrence must introduce a genuinely global normalization estimate; the
identities plus thin-shell information about the bulk radius do not by
themselves contain the Cheeger height.

## 8. Audit conclusion for this mechanism

The two dimension-free facts available without conjecture-strength input are

1. normal alignment repairs the cut to a subgraph with the quantitative error
   (7), and
2. small horizontal perimeter produces the balanced almost-pure phase
   classifier (8)--(9).

To obtain the fixed-point recurrence (11), one still needs one of the
following genuinely new statements:

* a horizontal-perimeter-versus-transverse-wall dichotomy controlling the
  conditional velocity in (12) only in the small-horizontal-perimeter regime;
  or
* a global normal-flow estimate converting the curvature branch into an
  absolute perimeter bound with a remainder uniform in dimension.

The exact quantile CMC equation (14l) shows that the first option cannot be
obtained from local stationarity: its conditional-score term is suppressed by
the physical slope and its zeroth-order velocity cancels from the Jacobi form.
The surviving problem is therefore global.  It must couple the phase change
of the conditional masses to either a global transport-ray/BV interface or to
bulk normalization along the entire conditional quantile foliation.  Any use
of conditional Fisher information is excluded by (14n) and the convex-body
approximation test.

The regular simplex rules out the stronger, tempting BV contraction (20), and
the product decoration rules out replacing a nearly graphical boundary by a
halfspace.  Neither missing statement follows from the current boundary
identities, so no valid recurrence has yet been obtained.
