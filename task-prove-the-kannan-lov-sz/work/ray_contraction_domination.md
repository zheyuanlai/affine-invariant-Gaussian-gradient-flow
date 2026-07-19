# Sharp domination under packetwise half-ray contractions

## Executive result

The proposed domination is true, with a sharp dimension-free constant.
Let `q` be a log-concave probability density on an interval, with median
zero, and let

\[
 S(t)=
 \begin{cases}
  a_+t,&t>0,\\
  a_-t,&t<0,
 \end{cases}
 \qquad 0<a_\pm\le1.                                    \tag{0.1}
\]

Then

\[
 \boxed{\quad
 {d(S_\#(q\,dt))\over q\,dt}\le
 \max_{\epsilon\in\{+,-\}}
 {2^{\,1-a_\epsilon}\over a_\epsilon}.
 \quad}                                                   \tag{0.2}
\]

For `a_\pm in [a_0,1]` the sharp uniform constant is

\[
 \boxed{\quad K(a_0)={2^{\,1-a_0}\over a_0}.\quad}       \tag{0.3}
\]

The factor `1/a_0` is attained by a uniform interval.  The additional
factor `2^{1-a_0}` is sharp for asymmetric two-sided exponentials and,
geometrically, for high-dimensional homothetic simplex fans.

More generally, if

\[
 p_-=\int_{-\infty}^0q,\qquad p_+=\int_0^\infty q,
\]

then the sharp sidewise constants are

\[
 K_+(a)={p_-^{-(1-a)}\over a},\qquad
 K_-(a)={p_+^{-(1-a)}\over a}.                          \tag{0.4}
\]

Thus `p_\pm>=delta` gives

\[
 K_\delta(a_0)={\delta^{-(1-a_0)}\over a_0}.            \tag{0.5}
\]

For a measurable calibrated-ray packet and measurable coefficients
`a_{y,\pm} in [a_0,1]`, the same conditional estimate integrates exactly:

\[
 \boxed{\quad T_\#\mu\le K_\delta(a_0)\mu.\quad}         \tag{0.6}
\]

There is no missing quotient or transverse Jacobian.  The ray quotient is
unchanged, and the one-dimensional factor `1/a_{y,\pm}` is the complete
coordinate Jacobian.  In a smooth normal chart the transverse Jacobian
ratio is already part of the conditional density `q_y`.  At nonsmooth
branching and focal sets, abstract disintegration proves (0.6) without
differentiation; those sets have zero conditional one-dimensional mass.

Domination gives useful packetwise constraints.  If `mu` is isotropic and
log-concave, every probability `nu<=Kmu` satisfies

\[
 |E_\nu X|\le C\log(2K),\qquad
 \sup_{|\theta|=1}
 \left(E_\nu|\langle\theta,X\rangle|^k\right)^{1/k}
 \le C\bigl(k+\log(2K)\bigr).                           \tag{0.7}
\]

For the ray half-moment vector measures

\[
 dM_+(y)=E[(T_y)_+\mid y]N_y\,d\eta(y),\qquad
 dM_-(y)=E[(-T_y)_+\mid y]N_y\,d\eta(y),                \tag{0.8}
\]

arbitrary packetwise choices imply

\[
 \|M_\pm\|_{\rm sv}
 :=\sup_{|g|\le1}\left|\int g\,dM_\pm\right|
 \le {C\log(2K)\over1-a_0}.                             \tag{0.9}
\]

There is an all-order version.  Put

\[
 m_{k,+}(y)=E[(T_y)_+^k\mid y],\qquad
 m_{k,-}(y)=E[(-T_y)_+^k\mid y].
\]

For every unit `theta`,

\[
 \boxed{\quad
 \int m_{k,\pm}(y)|\langle N_y,\theta\rangle|^k\,d\eta(y)
 \le
 \left[
 {C\bigl(k+\log(2K)\bigr)\over1-a_0}
 \right]^k.
 \quad}                                                   \tag{0.10}
\]

If the rays are `delta`-balanced and `sigma_y>=s` on `F`, then

\[
 \int_F|\langle N_y,\theta\rangle|^k\,d\eta(y)
 \le
 \left[
 {C_\delta\bigl(k+\log(2K)\bigr)\over(1-a_0)s}
 \right]^k.                                             \tag{0.11}
\]

Optimizing `k` recovers exponential smallness of every coherent long-ray
packet.  It also shows precisely why these contractions do not close the
inverse theorem.  A spherical direction law in rank `R>=s^4`, together
with bounded-support one-dimensional conditionals of scale `s`, satisfies
(0.9)--(0.11) and the new `s^5` projection law simultaneously.  The padded
incidence construction in `radial_endpoint_trichotomy.md` remains a concrete
non-log-concave survivor of every numerical constraint.

No explicit globally log-concave large-`s` survivor is produced: constructing
one would solve the unresolved global ray-realization problem.  Conversely,
(0.6)--(0.11) alone do not force `s=O(1)`.  They control conditional mass
and direction entropy but are insensitive to extremely remote endpoints of
exponentially small conditional density.  A global gluing or turning theorem
is still required.

## 1. The sharp one-dimensional theorem

### 1.1 Density at a quantile versus the mode

Let `q=e^\phi` be log-concave on its interval of support and let
`M=\|q\|_\infty`.  Put

\[
 p_-=\int_{-\infty}^0q(t)\,dt,\qquad p_+=1-p_-.
\]

**Lemma 1.1 (one-sided peak bound).**  If a mode lies to the right of zero,
then

\[
 q(0)\ge p_-M.                                          \tag{1.1}
\]

If a mode lies to the left, the reflected statement is
`q(0)>=p_+M`.  In particular, at a median,

\[
 q(0)\ge {M\over2}.                                     \tag{1.2}
\]

**Proof.**  Let `m>0` be a mode and set

\[
 \lambda={\phi(m)-\phi(0)\over m}.
\]

The assertion is trivial if `lambda=0`.  Concavity of `phi` gives

\[
 q(t)\le q(0)e^{\lambda t}\quad(t\le0),\qquad
 q(t)\ge q(0)e^{\lambda t}\quad(0\le t\le m).           \tag{1.3}
\]

Consequently

\[
 p_-\le {q(0)\over\lambda},\qquad
 p_+\ge {M-q(0)\over\lambda}.                           \tag{1.4}
\]

The first inequality gives `lambda<=q(0)/p_-`.  Insert this into the
second:

\[
 1-p_-\ge {p_-(M-q(0))\over q(0)}.
\]

Rearrangement yields (1.1).  Reflection proves the other side, and
`p_-=p_+=1/2` gives (1.2).  QED.

### 1.2 Exact contraction ratio

**Theorem 1.2 (sharp half-line contraction).**  For `a in (0,1]` and
`x>0`,

\[
 {q(x/a)\over a\,q(x)}
 \le {p_-^{-(1-a)}\over a}.                             \tag{1.5}
\]

For `x<0`,

\[
 {q(x/a)\over a\,q(x)}
 \le {p_+^{-(1-a)}\over a}.                             \tag{1.6}
\]

Both constants are sharp.

**Proof.**  Write `v=x/a`, so `x=av`.  If `q(v)<=q(av)`, the density
ratio before the factor `1/a` is at most one.  Otherwise a mode lies to the
right of zero.  Concavity between `0` and `v` gives

\[
 \phi(av)\ge(1-a)\phi(0)+a\phi(v).
\]

Therefore, using Lemma 1.1,

\[
 {q(v)\over q(av)}
 \le\left({q(v)\over q(0)}\right)^{1-a}
 \le\left({M\over q(0)}\right)^{1-a}
 \le p_-^{-(1-a)}.                                     \tag{1.7}
\]

This proves (1.5); reflection proves (1.6).

For sharpness, take the asymmetric Laplace density proportional to

\[
 q(t)=
 \begin{cases}
  e^{\lambda(t-m)},&t\le m,\\
  e^{-\mu(t-m)},&t\ge m,
 \end{cases}                                            \tag{1.8}
\]

with `m>0` chosen so that the left mass is `p_-`.  Direct integration gives

\[
 e^{-\lambda m}=p_-\left(1+{\lambda\over\mu}\right).
                                                               \tag{1.9}
\]

Letting `mu/lambda` tend to infinity gives
`exp(lambda m) to 1/p_-`.  At the image point `x=am`,

\[
 {q(m)\over a q(am)}
 ={e^{(1-a)\lambda m}\over a}
 \longrightarrow {p_-^{-(1-a)}\over a}.               \tag{1.10}
\]

Thus (1.5) is sharp, and reflection treats (1.6).  QED.

The function `a mapsto p^{-(1-a)}/a` is decreasing on `(0,1]`.  Equations
(0.2)--(0.5) follow immediately from the change-of-variables formula

\[
 {d(S_\#(q\,dt))\over dt}(x)
 ={1\over a_\pm}q(x/a_\pm)                              \tag{1.11}
\]

on the corresponding contracted half-interval.  Outside that image the
pushforward density is zero.

## 2. Passage through the ray quotient

Let `{\cal Y}` be the quotient of the nondegenerate transport rays.  Up to
a `mu`-null set there is a measurable one-to-one parametrization

\[
 F(y,t)=z_y+tN_y,\qquad
 d\mu(F(y,t))=q_y(t)\,dt\,d\eta(y).                     \tag{2.1}
\]

Let `P subset{\cal Y}` be measurable.  On `P` choose measurable
`a_{y,\pm} in[a_0,1]`, and put `a_{y,\pm}=1` off `P`.  Define

\[
 S_y(t)=a_{y,\operatorname {sgn}t}t,\qquad
 T(F(y,t))=F(y,S_y(t)).                                 \tag{2.2}
\]

**Theorem 2.1 (global domination).**  If every `q_y` has median zero, then

\[
 T_\#\mu\le {2^{1-a_0}\over a_0}\mu.                   \tag{2.3}
\]

If both sign masses are at least `delta`, then

\[
 T_\#\mu\le {\delta^{-(1-a_0)}\over a_0}\mu.            \tag{2.4}
\]

**Proof.**  For a Borel set `A`, let

\[
 A_y=\{t:F(y,t)\in A\}.
\]

The quotient label is unchanged by `T`.  Theorem 1.2 gives, for
`eta`-almost every `y`,

\[
 (S_y)_\#(q_y\,dt)(A_y)\le K(q_y,a_0)\int_{A_y}q_y(t)\,dt,
\]

with `K` bounded by (2.3) or (2.4).  Integrating in `y` proves the measure
inequality for arbitrary `A`.  QED.

### 2.1 Where the ambient Jacobian went

On a smooth ray chart, write the ambient volume element as

\[
 dx=J(y,t)\,dtdy
\]

and the original density as `rho(F(y,t))`.  Up to the quotient
normalization,

\[
 q_y(t)=\rho(F(y,t))J(y,t).                             \tag{2.5}
\]

The coordinate map is

\[
 (y,t)\longmapsto(y,a_y t).
\]

Even when `a_y` varies measurably or smoothly with `y`, its coordinate
Jacobian is triangular and has determinant `a_y`.  At a differentiability
point, the ambient Jacobian is

\[
 \det DT(F(y,t))
 =a_y{J(y,a_yt)\over J(y,t)}.                           \tag{2.6}
\]

Thus the ambient pushforward ratio at `F(y,s)` is exactly

\[
 {\rho(F(y,s/a_y))J(y,s/a_y)\over
   a_y\rho(F(y,s))J(y,s)}
 ={q_y(s/a_y)\over a_yq_y(s)}.                         \tag{2.7}
\]

There is no additional dimensional determinant.

### 2.2 Nonsmooth and branching rays

Normal charts can fail at medial, focal, fan, and branching sets.  The proof
of Theorem 2.1 never used a chart.  The Euclidean transport-ray equivalence
classes form a nonbranching measurable partition after deleting a null set.
Distinct ray interiors cannot cross without belonging to the same maximal
calibrated line.  Conditional densities have no atoms, so finite endpoints,
zero-level intersections, and branching endpoints carry zero conditional
mass.  Define `T` arbitrarily on that exceptional set.

The abstract disintegration proof then gives (2.3)--(2.4) verbatim.  In
particular, images from distinct ray interiors do not create an uncounted
multiplicity term.

## 3. Mandatory sharpness and geometry tests

### 3.1 Uniform interval

For `q=1/2` on `[-1,1]`, zero is a median.  Contracting either half by `a`
gives density `1/(2a)` on its image.  Hence

\[
 \left\|{dS_\#\mu\over d\mu}\right\|_\infty={1\over a}. \tag{3.1}
\]

This forces the Jacobian factor `1/a`.  The peak factor in (0.2) is absent
because the density is flat.

### 3.2 Asymmetric exponential

The family (1.8)--(1.10), with `p_-=1/2`, gives

\[
 \lim_{\mu/\lambda\to\infty}
 \left\|{dS_\#\mu\over d\mu}\right\|_\infty
 ={2^{1-a}\over a}.                                    \tag{3.2}
\]

Thus no smaller universal median constant is possible, already in one
dimension.

### 3.3 Homothetic simplex fan

Let a uniform `n`-simplex be parametrized from one vertex by a facet
coordinate `u` and a homothety coordinate `r in[0,1]`.  The volume element
is

\[
 dx=c(u)r^{n-1}\,dr\,du,
\]

so every quotient conditional is

\[
 q_n(r)=nr^{n-1}{\bf1}_{[0,1]}(r).                     \tag{3.3}
\]

Its median is

\[
 r_0=2^{-1/n}.                                         \tag{3.4}
\]

Contract about the median slice:

\[
 r'=r_0+a(r-r_0)
\]

on the outer half.  Direct ambient differentiation gives

\[
 \det DT=a\left({r'\over r}\right)^{n-1},
\qquad
 {dT_\#\mu\over d\mu}(r')
 ={1\over a}\left({r\over r'}\right)^{n-1}.             \tag{3.5}
\]

At the preimage `r=1`,

\[
 r'=r_0+a(1-r_0)
\]

and therefore

\[
 {1\over a}\left({1\over r'}\right)^{n-1}
 \longrightarrow {2^{1-a}\over a}.                    \tag{3.6}
\]

This fan realizes the sharp constant asymptotically in dimension.  It also
shows why treating the map as an isotropic contraction with determinant
`a^n` is wrong: the transverse factor is the ratio of the two ray
Jacobians and is already encoded in `q_n`.

The homothetic coordinate test does not require the fan to be the transport
fan of a particular signed-distance function; it tests the quotient and
Jacobian assertion itself.  For an actual polyhedral signed-distance fan,
each open cell has a smooth affine ray chart, the estimate holds cellwise,
and the interfaces and vertices are null.  Gluing the cells therefore
introduces no extra factor.

## 4. What `K`-domination gives without Poincare or Cheeger

The pushforward `T_\#mu` need not be log-concave.  Pointwise domination is
nevertheless enough for all one-dimensional moments.

**Lemma 4.1 (moments of a dominated perturbation).**  Let `mu` be isotropic
and log-concave and let `nu` be a probability with `nu<=Kmu`.  Then

\[
 |E_\nu X|\le C\log(2K),                                \tag{4.1}
\]

and, for every integer `k>=1`,

\[
 \sup_{|\theta|=1}
 \left(E_\nu|\langle\theta,X\rangle|^k\right)^{1/k}
 \le C\bigl(k+\log(2K)\bigr).                           \tag{4.2}
\]

In particular,

\[
 \|E_\nu XX^T\|_{op}
 \le C\bigl(1+\log K\bigr)^2.                           \tag{4.3}
\]

**Proof.**  For a unit `theta`, the isotropic log-concave marginal
`Z=<theta,X>` satisfies

\[
 P_\mu(|Z|\ge t)\le2e^{-ct}
\]

after changing universal constants.  Domination gives

\[
 P_\nu(|Z|\ge t)\le\min\{1,2Ke^{-ct}\}.                 \tag{4.4}
\]

Layer-cake integration, split at `C log(2K)`, proves (4.2).  Taking
`theta` in the direction of `E_\nu X` and using the one-sided version of
(4.4) proves (4.1).  The case `k=2` gives (4.3).  No functional inequality
is used.  QED.

## 5. Semivariation of the half-ray moment measures

For `k>=1` put

\[
 m_{k,+}(y)=\int_0^\infty t^kq_y(t)\,dt,\qquad
 m_{k,-}(y)=\int_{-\infty}^0(-t)^kq_y(t)\,dt.           \tag{5.1}
\]

### 5.1 First moment: vector semivariation

Let

\[
 M_\pm(A)=\int_A m_{1,\pm}(y)N_y\,d\eta(y).             \tag{5.2}
\]

Assume `a_0<1` and let `K=K_\delta(a_0)`.  Contract only the positive
half-rays in a packet `A` by `a_0`.  Since `E_\mu X=0`,

\[
 E_{T_\#\mu}X=-(1-a_0)M_+(A).                          \tag{5.3}
\]

Lemma 4.1 and Theorem 2.1 give

\[
 |M_+(A)|\le {C\log(2K)\over1-a_0}.                    \tag{5.4}
\]

Contracting the negative half gives the same estimate for `M_-(A)`.

More generally, for any measurable `g:{\cal Y}->[0,1]` choose

\[
 a_{y,+}=1-(1-a_0)g(y).
\]

The same calculation bounds `|\int g\,dM_+|`.  Decomposing a signed
`g in[-1,1]` into positive and negative parts proves

\[
 \boxed{\quad
 \|M_\pm\|_{\rm sv}
 \le {2C\log(2K)\over1-a_0}.
 \quad}                                                   \tag{5.5}
\]

Equivalently, for every unit `theta`,

\[
 \int m_{1,\pm}(y)|\langle N_y,\theta\rangle|\,d\eta(y)
 \le {2C\log(2K)\over1-a_0}.                            \tag{5.6}
\]

### 5.2 All half moments

Fix a packet `A`, a unit `theta`, and contract only its positive halves by
a common factor `a in[a_0,1]`.  Write

\[
 u_y=\langle\theta,z_y\rangle,\qquad
 v_y=\langle\theta,N_y\rangle.
\]

The difference of the `k`th raw moments is the polynomial

\[
 R_A(a)=
 \sum_{j=1}^k {k\choose j}(a^j-1)
 \int_Au_y^{k-j}v_y^jm_{j,+}(y)\,d\eta(y).              \tag{5.7}
\]

By (4.2),

\[
 \sup_{a\in[a_0,1]}|R_A(a)|
 \le2\left[C\bigl(k+\log(2K)\bigr)\right]^k.            \tag{5.8}
\]

The leading coefficient of (5.7) is

\[
 \int_Av_y^km_{k,+}(y)\,d\eta(y).                       \tag{5.9}
\]

If a degree-`k` polynomial is bounded by `H` on an interval of length
`1-a_0`, its leading coefficient is at most

\[
 \left({4\over1-a_0}\right)^kH;                         \tag{5.10}
\]

this is the extremal leading-coefficient bound for Chebyshev polynomials
after rescaling the interval.  Apply (5.10) to (5.8).  Since `A` was
arbitrary, the total variation of the scalar measure in (5.9) is at most
twice the supremum over packets.  Therefore

\[
 \boxed{\quad
 \int m_{k,+}(y)|\langle N_y,\theta\rangle|^k\,d\eta(y)
 \le
 \left[
 {C\bigl(k+\log(2K)\bigr)\over1-a_0}
 \right]^k.
 \quad}                                                   \tag{5.11}
\]

Reflection proves the negative-half estimate.  Bounds for every lower
coefficient in (5.7), and hence for mixed base/ray moment measures, follow
from the same polynomial argument.  The leading coefficient (5.11) is the
one free of the uncontrolled ray bases `z_y`.

For `k=2`, the same statement can be written as the covariance
semivariation bound

\[
 \sup_{|g|\le1}
 \left\|
 \int g(y)m_{2,\pm}(y)N_yN_y^T\,d\eta(y)
 \right\|_{op}
 \le
 \left[{C(2+\log(2K))\over1-a_0}\right]^2.              \tag{5.11a}
\]

Thus the covariance information from all packetwise contractions is already
contained in the `k=2` member of the moment hierarchy.

### 5.3 Balance converts half moments to direction moments

For a one-dimensional log-concave density of variance `sigma^2`,
`\|q\|_\infty<=C/sigma`.  If its positive mass is at least `delta`,
the density-cap rearrangement inequality gives

\[
 m_{k,+}\ge
 {\delta^{k+1}\sigma^k\over C^k(k+1)}
 \ge(c\delta^2\sigma)^k.                               \tag{5.12}
\]

Indeed, among nonnegative densities of mass `delta` bounded by
`C/sigma`, the smallest `k`th moment packs all mass immediately next to
zero.  The same estimate holds on the negative side.

Combining (5.11) and (5.12), on a family `F` with `sigma_y>=s`,
gives

\[
 \boxed{\quad
 \int_F|\langle N_y,\theta\rangle|^k\,d\eta(y)
 \le
 \left[
 {C_\delta\bigl(k+\log(2K)\bigr)\over(1-a_0)s}
 \right]^k.
 \quad}                                                   \tag{5.13}
\]

If `A subset F` satisfies `|<N_y,theta>|>=gamma`, then

\[
 \eta(A)\le
 \left[
 {C_\delta\bigl(k+\log(2K)\bigr)\over
        (1-a_0)s\gamma}
 \right]^k.                                             \tag{5.14}
\]

Choosing `k` to be a sufficiently small constant multiple of
`(1-a_0)s gamma` yields, whenever that quantity is larger than a constant,

\[
 \eta(A)\le
 C\exp[-c_{\delta,a_0}s\gamma].                         \tag{5.15}
\]

Thus arbitrary packetwise contractions recover the exponential
coherent-direction restriction by a moment method.  They also give the
stronger signed semivariation statement (5.5).

## 6. Interaction with the `s^5` projection law

The new radial projection theorem gives, for a boundary-weighted packet,

\[
 P_Fs^5E_\nu|PN|^4\le C\,\operatorname {rank}P.         \tag{6.1}
\]

For positive actual ray mass on one variance scale this forces normal
support rank at least order `s^4`.  Equations (5.13)--(5.15) force
exponential angular diffuseness.  These are compatible, not contradictory.

To see this explicitly, let `N` be uniform on `S^{R-1}`, take

\[
 R\ge s^4,
\]

and put the uniform conditional density of variance `s^2` on every abstract
ray:

\[
 q_y(t)={1\over2\sqrt3s}{\bf1}_{[-\sqrt3s,\sqrt3s]}(t). \tag{6.2}
\]

Then `P_F asymp1/s` and the full-rank case of (6.1) is saturated in scale.
For every unit `theta` and `k<=R`,

\[
 \left(E|\langle N,\theta\rangle|^k\right)^{1/k}
 \le C\sqrt{k\over R}.                                  \tag{6.3}
\]

Since `m_{k,\pm}^{1/k}<=sqrt3s`, the left side of (5.11), to the power
`1/k`, is at most

\[
 Cs\sqrt{k\over R}\le {C\sqrt k\over s},
\]

which is far below the allowed `Ck`.  For `k>R` the trivial bound also
passes.  Spherical cap probabilities are `exp(-cR)`, much smaller than
(5.15).

The padded projective-plane construction of
`radial_endpoint_trichotomy.md` adds near-maximal ideal endpoint gap and
large absolute non-Clifford defect.  Taking its active rank much larger than
`s^4` and its conditional density to be a truncated Gaussian of scale `s`
makes it satisfy (5.5), (5.11), and (6.1); its finite endpoint density
is

\[
 \exp[-\Theta(m/s^2)]
\]

times the central density.

These examples are abstract ray mixtures, not full-dimensional log-concave
signed-distance realizations.  They prove that the inequalities themselves
do not imply `s=O(1)`.  Producing an explicit globally log-concave survivor
would amount to solving the remaining global compatibility problem: the
ray mixture must be filled by log-concave bridge mass without destroying
calibration or the conditional laws.  No such survivor is supplied by the
interval, exponential, simplex, spherical, or padded-incidence tests.

## 7. What the contractions do and do not see

The contraction map moves conditional mass near the zero level.  Its
domination constant is independent of the distance to a finite endpoint.
Accordingly:

* it controls the vector semivariation of conditional half means;
* it controls every projected half moment and every coherent packet;
* together with thin shell it forces high active normal rank;
* it does not control a support endpoint carrying exponentially small
  conditional density.

The last point is exactly the far-endpoint escape in
`radial_endpoint_trichotomy.md`.  A truncated Gaussian can have standard
deviation `s` and support half-length `sqrt m` with endpoint density
`exp[-Theta(m/s^2)]/s`.  Every contraction considered here acts on the
central `O(s)` mass and is essentially insensitive to that support length.

Therefore the sharp domination theorem is a genuine new reusable tool, but
it does not by itself close the long-ray inverse step.  Closing it requires
an estimate that couples central conditional mass to endpoint/focal
geometry or charges the full BV turning needed to glue exponentially many
direction packets.
