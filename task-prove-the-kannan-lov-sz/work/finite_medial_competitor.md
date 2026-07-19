# Finite medial competitors: exact envelope inequalities and the endpoint obstruction

## Executive conclusion

There is an exact finite-amplitude replacement for the normal-height second
variation.  If `f` is a true `T3` extremizer, its balanced-ray quotient is
`Q`, and a one-Lipschitz competitor has the form

\[
                 g=f-h(Q)+r,
\]

then

\[
 \boxed{\qquad
   \int B_y(h(y)-\eta h)\,d\eta(y)\le 2\int|r|\,d\mu .
 \qquad}                                                     \tag{0.1}
\]

Here `B_y` is an explicit one-dimensional convex Bregman term.  On a
balanced log-concave ray of scale `s`, it is bounded below by a constant
times `min(h^2/s,|h|)`.  Thus any finite merge, offset, or rerouting which is
raywise constant up to a sufficiently small `L^1` defect strictly beats
`f`.

For a polyhedral min envelope, the defect in (0.1) is not implicit.  If

\[
 F=\min_i\ell_i,\qquad F_h=\min_i(\ell_i-h_i),
\]

then at every point

\[
 F_h=F-h_{I(x)}-e_h(x),\qquad
 e_h(x)=\max_j\{h_j-h_{I(x)}-(\ell_j-F)(x)\}_+.          \tag{0.2}
\]

Consequently extremality gives a nonlinear switching-layer inequality valid
at every amplitude.  Its small-amplitude limit is the medial graph energy,
but (0.2) also records its saturation at finite amplitude.

The finite inequality does **not**, by itself, exclude a complete or
expanding endpoint graph.  Random and multilevel offsets have no factor
depending on the number of packets: their total ideal gain is at most of
order the height range.  A complete envelope can turn every height vector
into one global translate and can saturate the inequality at every
amplitude.  The long-core cross-slack lemma proves that separated directions
cannot realize this collapse in the interior of long calibrated rays; all
such saturation must occur at medial or focal endpoints.  Bounding the mass
or the rematching capacity of that endpoint region remains the missing
Euclidean step.

There is a sharp finite rematching result under an additional packet
separation hypothesis.  A packet whose `(d+Delta)`-neighborhood fits in one
half of a flattened `d`-separated endpoint model gives a new half-cut and a
one-Lipschitz Kantorovich potential which improves the objective by exactly a
positive packet-mass multiple of `Delta`.  If same-sign packets are farther
apart than the calibrated cross length, extremality forces every cross
distance to be equal.  The resulting configuration is orthogonal-spherical
(the Clifford branch), and a concrete projected-radial competitor plus
translated thin shell bounds its scale universally.

Thus the finite program gives a rigorous dichotomy, but not a dimension-free
KLS proof:

1. a finite height has switching defect smaller than its explicit Bregman
   gain, and then it is a strict competitor;
2. finite rematching exposes a nontrivial cross-distance gap, and then a new
   cut is a strict competitor;
3. the endpoint relation is approximately complete, in which case one needs
   a robust orthogonal-spherical inverse theorem; or
4. endpoint/focal mass implements an all-amplitude expander saturation.

No argument presently derived from log-concavity, isotropy, and calibration
rules out item 4.  A full Euclidean/log-concave long balanced example in item
4 would itself be an unbounded `T3` example, so producing it is not a
legitimate local countertest; the bounded-scale Gaussian simplex is the
available sharp Euclidean countertest.

## 1. The exact finite calibration identity

For an integrable function `u`, write

\[
             \mathcal J(u)=\int|u-\mu u|\,d\mu .          \tag{1.1}
\]

Let `f` be an extremizer of `J` among one-Lipschitz functions.  Subtract a
constant so that `mu f=0`, put

\[
 \sigma=\operatorname {sgn}f,\qquad
 \alpha=\int\sigma\,d\mu,
\]

and recall that `mu(f=0)=0` for a nonconstant extremizer.

**Lemma 1.1 (exact finite Bregman identity).**  For every integrable `g`, put
`delta=g-f` and `v=delta-mu delta`.  Then

\[
 \boxed{\begin{aligned}
 \mathcal J(g)-\mathcal J(f)
   ={}&\int(\sigma-\alpha)\delta\,d\mu\\
     &+2\int |g-\mu g|\,
       {\bf1}_{\{f(g-\mu g)<0\}}\,d\mu .
 \end{aligned}}                                             \tag{1.2}
\]

If `Q` is the quotient map of the balanced calibrated rays, then

\[
 \mathbb E(\sigma\mid Q)=\alpha                           \tag{1.3}
\]

and hence

\[
 \int(\sigma-\alpha)\delta\,d\mu
 =\int(\sigma-\alpha)
       \bigl(\delta-\mathbb E(\delta\mid Q)\bigr)d\mu .  \tag{1.4}
\]

**Proof.**  For real `a!=0` and arbitrary `b`,

\[
 |a+b|-|a|
 =\operatorname {sgn}(a)b
  +2|a+b|{\bf1}_{\{a(a+b)<0\}}.                         \tag{1.5}
\]

Apply (1.5) with `a=f` and `b=v`.  Since

\[
 \int\sigma v\,d\mu
 =\int(\sigma-\alpha)\delta\,d\mu,
\]

this proves (1.2).  Exact raywise balance gives (1.3), and conditional
centering gives (1.4).  QED.

An important consequence is finite rather than infinitesimal.  If a
nonconstant ray height `h(Q)` could be added while preserving the Lipschitz
constraint exactly, the linear term in (1.2) would vanish and every positive
mass of sign crossings would strictly increase `J`.  Thus the whole issue is
the finite defect required to make a raywise height globally one-Lipschitz.

Disintegrate

\[
 d\mu=\int d\nu_y\,d\eta(y),\qquad T=f,
\]

and write

\[
 p=\nu_y(T>0),\qquad q=\nu_y(T<0),\qquad \alpha=p-q.      \tag{1.6}
\]

The two numbers do not depend on `y`.  Define

\[
 B_y(c)=\int|T-c|\,d\nu_y-\int|T|\,d\nu_y+\alpha c.     \tag{1.7}
\]

**Lemma 1.2 (exact raywise gain).**  If `h` is integrable and
`c(y)=h(y)-eta h`, then

\[
 \boxed{\quad
 \mathcal J(f-h\circ Q)-\mathcal J(f)
       =\int B_y(c(y))\,d\eta(y).\quad}                  \tag{1.8}
\]

Moreover

\[
 B_y(c)=
 \begin{cases}
  2\displaystyle\int_{0<T<c}(c-T)\,d\nu_y,&c\ge0,\\[6pt]
  2\displaystyle\int_{c<T<0}(T-c)\,d\nu_y,&c<0.
 \end{cases}                                             \tag{1.9}
\]

In particular `B_y>=0`.

**Proof.**  Centering `f-h(Q)` replaces `h` by `c`.  Add and subtract
`alpha c` in each conditional integral.  Its quotient integral is zero.
This proves (1.8).  For `c>=0`, differentiate the left side of (1.7) from
zero to `c` and use

\[
 {d\over dc}\int|T-c|d\nu_y=2\nu_y(T<c)-1.
\]

The derivative after adding `alpha c` is
`2 nu_y(0<T<c)`.  A second integration gives (1.9).  The negative case is
identical.  QED.

**Proposition 1.3 (finite height-defect inequality).**  Suppose `g` is
one-Lipschitz and, for some quotient height `h`,

\[
                       g=f-h(Q)+r.                       \tag{1.10}
\]

Then extremality of `f` implies

\[
 \boxed{\quad
       \int B_y(h(y)-\eta h)\,d\eta(y)
       \le2\int|r|\,d\mu .\quad}                        \tag{1.11}
\]

Consequently `g` is a strict competitor whenever the reverse strict
inequality holds.

**Proof.**  The centered `L^1` seminorm obeys

\[
 |\mathcal J(u+r)-\mathcal J(u)|
 \le\int|r-\mu r|d\mu\le2\int|r|d\mu.                  \tag{1.12}
\]

Apply this with `u=f-h(Q)`, use (1.8), and then use
`J(g)<=J(f)`.  QED.

This proposition applies to every finite construction below and also to any
future nonsmooth merge or rerouting.  It avoids differentiating a cut locus.

## 2. Quantitative one-dimensional gain with constants

The next estimates identify the exact scale of the left side of (1.11).

**Lemma 2.1 (finite Bregman lower bound).**  Let `nu` have a log-concave
density `varphi`, let

\[
 p=\nu(0,\infty),\quad q=\nu(-\infty,0),\quad b=\varphi(0),
\]

and define `B` by (1.7).  Then

\[
 B(c)\ge
 \begin{cases}
  \frac12\min\{bc^2,pc\},&c\ge0,\\[3pt]
  \frac12\min\{bc^2,q|c|\},&c<0.
 \end{cases}                                             \tag{2.1}
\]

Also, with

\[
 m_+=\int T_+d\nu,\qquad m_-=\int T_-d\nu,
\]

one has the large-amplitude bounds

\[
 B(c)\ge
 \begin{cases}
       2pc-2m_+,&c\ge0,\\
       2q|c|-2m_-,&c<0.
 \end{cases}                                             \tag{2.2}
\]

**Proof.**  For `c>=0`, let `S(t)=nu(T>=t)`.  The survival function of a
one-dimensional log-concave density is log-concave, `S(0)=p`, and
`S'_+(0)=-b`.  Therefore

\[
                         S(t)\le p e^{-bt/p}.             \tag{2.3}
\]

Formula (1.9) gives, with `x=bc/p`,

\[
 B(c)=2\int_0^c(p-S(t))dt
 \ge {2p^2\over b}\bigl(x-1+e^{-x}\bigr).               \tag{2.4}
\]

For `0<=x<=1`, Taylor's alternating bound gives
`x-1+e^{-x}>=x^2/3`.  For `x>=1`, the ratio
`(x-1+e^{-x})/x` is increasing and is at least `e^{-1}`.
These two estimates imply (2.1), with room in the constant.  Formula (1.9)
also gives

\[
 B(c)=2pc-2m_++2\int_{T>c}(T-c)d\nu,
\]

which proves (2.2).  Reflecting the line proves the negative statements.
QED.

**Lemma 2.2 (balanced scale normalization).**  Suppose, for constants
`delta,a_0,a_1>0`,

\[
 p,q\ge\delta,qquad
 a_0s\le m_+,m_-\le a_1s.                              \tag{2.5}
\]

Then

\[
 {\delta^3\over4a_1s}\le b\le {1\over a_0s}.            \tag{2.6}
\]

Consequently, whenever `|c|<=delta a_0s`,

\[
             B(c)\ge {\delta^3\over6a_1s}c^2.           \tag{2.7}
\]

**Proof.**  Put `H=sup varphi`.  A standard one-dimensional consequence of
log-concavity and `p,q>=delta` is

\[
                             H\le b/\delta.              \tag{2.8}
\]

For completeness, if a mode of height `H` is to the right of zero, compare
the secant slope of `log varphi` between zero and the mode on both half
lines.  This gives `p>=(H-b)/a` and `q<=b/a`, hence
`H/b<=1+p/q<=1/delta`; the other side is symmetric.

A density bounded by `H` and carrying positive-half-line mass `p` has

\[
                         m_+\ge {p^2\over4H},             \tag{2.9}
\]

while (2.3) gives `m_+<=p^2/b`.  Equations (2.8)--(2.9),
`p>=delta`, and (2.5) imply (2.6).  The same result follows from the
negative side.  If `|c|<= delta a_0s`, then `bc/p<=1`; use the quadratic
part of (2.4) and the lower bound in (2.6).  QED.

Thus, for heights in a fixed small fraction of the long-ray scale, (1.11)
is exactly a quotient `L^2` gain of order `Var(h)/s`.  For arbitrary
amplitude it crosses over to a linear gain, with no loss of exactness.

There is a useful localized large-height form.  If `eta c=0`, put

\[
                         L=\int c_+d\eta=\int c_-d\eta.
\]

If (2.5) holds raywise, then (2.2) gives

\[
 \int B_y(c_y)d\eta(y)
 \ge 2L-2a_1s\,\eta(c\ne0).                            \tag{2.10}
\]

Unlike a binary height centered by shifting its large complement, (2.10)
shows that a genuinely mean-zero multilevel height pays the first-moment
correction only on its nonzero packets.

## 3. Exact polyhedral min-envelope formula

Let `ell_1,...,ell_M` be one-Lipschitz functions.  They may be affine
polyhedral charts or distance-cone McShane charts.  Put

\[
 F=\min_i\ell_i,qquad F_h=\min_i(\ell_i-h_i).            \tag{3.1}
\]

Both functions are one-Lipschitz.  Choose the least active index

\[
 I(x)=\min\{i:\ell_i(x)=F(x)\},qquad
 r_j(x)=\ell_j(x)-F(x)\ge0.                              \tag{3.2}
\]

**Lemma 3.1 (exact finite switching defect).**  For every height vector
`h in R^M`,

\[
 \boxed{\begin{aligned}
 F_h(x)&=F(x)-h_{I(x)}-e_h(x),\\
 e_h(x)&=\max_j\{h_j-h_{I(x)}-r_j(x)\}_+.
 \end{aligned}}                                         \tag{3.3}
\]

In particular, writing `R=osc(h)`,

\[
 \boxed{\quad
 \int e_h d\mu
 =\int_0^R
   \mu\{x:\exists j,\ h_j-h_{I(x)}-r_j(x)>t\}\,dt .\quad} \tag{3.4}
\]

**Proof.**  Directly from (3.1),

\[
 F-F_h=\max_j\{h_j-r_j\}.
\]

The active term is `h_I`, so subtracting it gives (3.3).  Layer cake gives
(3.4).  QED.

Assume now that `F=f` almost everywhere and that the active label is a
function of the ray quotient: `I=i(Q)`.  Put

\[
 \eta_i=\eta\{I=i\},\qquad \bar h=\sum_i\eta_i h_i,
 \qquad c_i=h_i-\bar h.                                 \tag{3.5}
\]

**Proposition 3.2 (all-amplitude envelope inequality).**  If `f` is a true
extremizer, then every height vector satisfies

\[
 \boxed{\quad
   \sum_i\eta_i B_i(c_i)\le2\int e_h\,d\mu,
 \quad}                                                  \tag{3.6}
\]

where `B_i` is the average of `B_y` over the rays with active label `i`.

**Proof.**  Formula (3.3) writes `F_h=f-h_I-e_h`.  Apply Proposition 1.3.
QED.

For a binary merge, let `H subset {1,...,M}`, take `h_i=a 1_H(i)`, and put
`r=eta(H)`.  Then

\[
 e_H(x)={\bf1}_{\{I(x)\notin H\}}
       \left(a-\min_{j\in H}r_j(x)\right)_+,             \tag{3.7}
\]

and hence

\[
 \int e_Hd\mu
 =\int_0^a\mu\left\{I\notin H,
       \min_{j\in H}(\ell_j-F)<t\right\}dt.             \tag{3.8}
\]

If all the rays satisfy (2.5) and `a<=delta a_0s`, then

\[
 \boxed{\quad
 \int_0^a\mu\left\{I\notin H,
       \min_{j\in H}(\ell_j-F)<t\right\}dt
 \ge {\delta^3\over12a_1}
       {a^2\over s}r(1-r).
 \quad}                                                  \tag{3.9}
\]

Indeed

\[
 \sum_i\eta_i c_i^2=a^2r(1-r),                          \tag{3.10}
\]

so (3.9) follows from (2.7) and (3.6).  Formula (3.9) is the
finite-amplitude nonlinear conductance inequality.  Differentiating it at
`a=0` recovers the medial graph Dirichlet form, but (3.8) remains valid
after switching strips overlap and the quadratic graph approximation has
ceased to be accurate.

The max-envelope formula follows by applying Lemma 3.1 to `-ell_i`.  More
general finite signed polyhedral envelopes are built by a lattice expression
of min and max gates.  Every such expression remains one-Lipschitz after
leaf offsets.  Recursively applying (3.3) gives a defect supported on gates
whose unperturbed child gap is at most `osc(h)`; its absolute value is at
most `osc(h)`.  Proposition 1.3 then applies without any smoothness or cut
locus calculation.

For clarity, the last assertion has the following formal version.

**Lemma 3.3 (finite min/max lattice bound).**  Let `L` be any fixed binary
tree of min and max gates with one-Lipschitz leaves `ell_1,...,ell_M`.  Put

\[
                  F=\mathcal L(\ell_1,\ldots,\ell_M),
 \qquad F_h=\mathcal L(\ell_1-h_1,\ldots,\ell_M-h_M).     \tag{3.11}
\]

At each point, choose an active leaf `I(x)` by following the winning child
at every unperturbed gate.  If `R=osc(h)` and `S_R` is the set of points at
which at least one gate has unperturbed child-output gap at most `R`, then

\[
 |F_h(x)-(F(x)-h_{I(x)})|\le R{\bf1}_{S_R}(x).           \tag{3.12}
\]

If `I` is constant on ray labels and `f=F` almost everywhere, extremality
therefore forces

\[
                 \int B_y(h_{I(y)}-\eta h_I)d\eta(y)
                 \le2R\mu(S_R).                         \tag{3.13}
\]

**Proof.**  Min and max are translation-equivariant and one-Lipschitz in
the leaf values for the sup norm.  Hence the output perturbation lies between
`-max h` and `-min h`, proving the bound by `R`.  If every child gap exceeds
`R`, no gate can change its winning child under a perturbation whose leaf
range is `R`; the same leaf remains active and the defect is zero.  This
proves (3.12).  Apply Proposition 1.3 for (3.13).  QED.

## 4. Random and multilevel offsets do not gain an entropy factor

The hope behind using exponentially many packets is that a random height
might add one gain per packet while a min or max envelope pays only once at
each ambient point.  Quotient normalization prevents that amplification.

**Lemma 4.1 (range bound for every multilevel gain).**  Let `c` have mean
zero and range `R`.  Then

\[
 \int |c|d\eta\le {R\over2},
 \qquad
 \int c_+d\eta=\int c_-d\eta\le {R\over4}.              \tag{4.1}
\]

For every balanced ray,

\[
 B_y(c)\le
 \begin{cases}2pc,&c\ge0,\\2q|c|,&c<0,
 \end{cases}                                             \tag{4.2}
\]

and therefore

\[
                 \int B_y(c_y)d\eta(y)\le {R\over2}.     \tag{4.3}
\]

The right side contains no number of packets.

**Proof.**  The first inequality in (4.1) is the sharp mean absolute
deviation bound for a random variable in an interval of length `R`.  If the
interval is `[-a,b]`, mean zero gives
`int c_+=int c_-<=ab/(a+b)<=R/4`.  Formula (1.9) gives (4.2).  Since
`p+q=1`, (4.3) follows from the second part of (4.1).  QED.

This estimate remains true after averaging over random heights.  In
particular neither independent Bernoulli offsets nor a continuous multilevel
height distribution produces a factor `log M`, `sqrt(log M)`, or `M` in the
ideal gain.

The random binary calculation can be made exact.  Assume for simplicity
that the `M` active cells have equal quotient mass and choose `H` uniformly
among the `k`-subsets, with `rho=k/M`.  At a point `x`, put

\[
 N_t(x)=|\{j\ne I(x):r_j(x)<t\}|.                       \tag{4.3a}
\]

Then (3.8) and the hypergeometric law give

\[
 \boxed{\begin{aligned}
 \mathbb E_H\int e_Hd\mu
 ={}&(1-\rho)\int_0^a\int
 \left[1-{
   {M-1-N_t(x)\choose k}
  \over {M-1\choose k}}
 \right]d\mu(x)dt,
 \end{aligned}}                                         \tag{4.3b}
\]

with the binomial coefficient interpreted as zero when necessary.  If the
ray laws are identical, the ideal gain is the deterministic number

\[
       \rho B(a(1-\rho))+(1-\rho)B(-a\rho).             \tag{4.3c}
\]

For independent Bernoulli-`rho` labels, the bracket in (4.3b) is instead
`1-(1-rho)^{N_t(x)}`.  Thus randomization merely replaces the edge count by
a *saturated near-active multiplicity*.  Once `rho N_t` is large, the
bracket is already one, regardless of how many further packets are present.
This is the exact reason exponential packet entropy creates no additional
finite-amplitude gain.

There is an exact all-amplitude saturation model.  Suppose all charts agree
on the region charged by the measure, so `ell_i=F` there.  Then

\[
 F_h=F-\max_i h_i.                                      \tag{4.4}
\]

It is merely a global translate, so `J(F_h)=J(F)` for every height vector,
at every amplitude.  Formula (3.3) gives

\[
                     e_h=\max_i h_i-h_{I}.               \tag{4.5}
\]

Thus the finite defect, not only the quadratic junction energy, can absorb
every multilevel gain.  Declaring arbitrary complete-graph conductances
between duplicate labels changes no finite competitor.  In Euclidean space
duplicate charts have the same normal and belong to the harmless parallel
branch; attaching unrelated normals to the labels is precisely the abstract
non-Euclidean step in the expander countermodel.

There is also an exact graph criterion showing which noncomplete expanders
survive *all* multilevel amplitudes in the zero-slack endpoint model.  Let
`G` be an undirected graph on a probability space of labels, include a loop
at every vertex, and define

\[
 e_h(i)=\max_{j\sim i}h_j-h_i.                            \tag{4.5a}
\]

For `H` a label set, let its directed outer vertex boundary be

\[
 \partial^+H=\{i\notin H:\exists j\in H, i\sim j\}.     \tag{4.5b}
\]

**Lemma 4.2 (all-amplitude expander blocker).**  Suppose

\[
             \eta(\partial^+H)\ge\eta(H)(1-\eta(H))      \tag{4.5c}
\]

for every measurable `H`.  Then, for every height `h` and every family of
balanced ray laws,

\[
                  \int B_y(h_y-\eta h)d\eta(y)
                  \le2\int e_h d\eta.                   \tag{4.5d}
\]

Thus the associated zero-slack max/min envelope can absorb every random or
multilevel finite offset, not merely its second variation.

**Proof.**  Subtract `eta h`.  The layer-cake identity for the one-sided
local maximum is

\[
 \int e_hd\eta
 =\int_{-\infty}^{\infty}
       \eta(\partial^+\{h>t\})dt.                        \tag{4.5e}
\]

If `H,H'` are independent copies of the height, then

\[
 \int_{-\infty}^{\infty}
   \eta(h>t)(1-\eta(h>t))dt
 ={1\over2}\mathbb E|H-H'|.                             \tag{4.5f}
\]

Equations (4.5c)--(4.5f), Jensen, and Lemma 4.1 give

\[
 2\int e_hd\eta
 \ge\mathbb E|H-H'|
 \ge\mathbb E|H-\mathbb EH|
 \ge\int B_y(h_y-\eta h)d\eta.
\]

QED.

Condition (4.5c) is sharp at the level of binary large amplitudes.  For a
set of height `a` and quotient mass `rho`, (2.2) gives ideal gain

\[
                 2a\rho(1-\rho)-O(s),                  \tag{4.5g}
\]

while the zero-slack envelope pays exactly
`2a eta(partial^+H)`.  Letting `a/s` tend to infinity forces (4.5c).
Therefore a generic spectral expander need not survive finite offsets: it
fails if some directed vertex boundary is smaller than `rho(1-rho)`.  A
complete graph, and any sufficiently strong vertex expander satisfying
(4.5c), does survive at every amplitude.  The finite test upgrades the
quadratic spectral-gap requirement to this nonlinear vertex-expansion
profile; it does not rule out the profile.

A more geometric complete graph can have switching slacks distributed over
a finite interval.  For the symmetric uniform needle on `[-s,s]`,

\[
 B(c)=
 \begin{cases}
 c^2/(2s),&|c|\le s,\\
 |c|-s/2,&|c|\ge s.
 \end{cases}                                             \tag{4.6}
\]

If a switching gap `R_0` is uniform on `[0,s]`, then

\[
       \mathbb E(d-R_0)_+=B(d)\qquad(d\ge0).             \tag{4.7}
\]

Thus an abstract complete envelope can match the ray Bregman function both
in its quadratic regime and in its linear saturation regime.  Random or
multilevel amplitudes do not create a finite deficit against such a model.
This is an exact scalar countermodel, not a Euclidean realization with
dispersed long normals.

The globally log-concave Euclidean countertest is the Gaussian regular
simplex of `euclidean_expander_exclusion.md`, Theorem 6.1.  Its interior
signed distance is a min of `M` affine facet distances, every pair of facets
has a medial interface, and symmetry gives the complete graph with equal
conductances.  It shows that Euclidean incidence and global log-concavity do
not prevent complete medial switching geometry.  It fails exactly where a full
counterexample must fail: its conditional ray variances are at most one and
its rays are not exactly balanced.

For reference, the precise countertest is as follows.

**Proposition 4.3 (Gaussian simplex: Euclidean complete medial graph).**
In

\[
 H=\{x\in\mathbb R^M:\sum_i x_i=0\},
\]

let `u_1,...,u_M` be regular-simplex unit vectors and put

\[
 K_a=\{x:\langle u_i,x\rangle\le a\ \text{for all }i\}.
\]

For signed distance `F` to `partial K_a`, positive inside, and standard
Gaussian measure on `H`,

\[
                  F(x)=\min_i(a-\langle u_i,x\rangle)
                  \quad(x\in K_a).                       \tag{4.8}
\]

Every pair of branches has a relatively open medial interface.  Facet
heights therefore have exact junction form

\[
             \mathcal E_{med}(h)=w_{M,a}
                   \sum_{i<j}(h_i-h_j)^2,qquad w_{M,a}>0, \tag{4.9}
\]

so the normalized graph is `K_M`.  Nevertheless every open normal-ray
conditional is a one-sided truncation of a unit Gaussian and has variance at
most one.  Its positive mass varies with the truncation endpoint, so exact
raywise balance also fails.

**Proof.**  Formula (4.8) is the distance from an interior point to the
union of the complementary facet halfspaces.  For `i!=j`, the equality
region of branches `i,j` has a nonempty relatively open codimension-one
piece; permutation symmetry gives the same positive coarea weight on every
pair, proving (4.9).  Along the inward normal to facet `i`, the Gaussian
density is proportional to

\[
               e^{-(a-t)^2/2}{\bf1}_{\{t<r_i\}},
\]

which has variance at most one by one-dimensional integration by parts.
Its mass on `t>0` plainly depends on `r_i`.  QED.

Thus Proposition 4.3 is a fully Euclidean and globally log-concave
all-incidence countermodel, but not a long-balanced-extremal countermodel.
Constructing the latter would already construct the bad `T3` witness the
whole argument is meant to exclude.

## 5. Long-core cross slack localizes all saturation at endpoints

The preceding saturation cannot take place between separated long-ray
directions in the interior.  We record the precise way in which the
cross-slack lemma enters a finite envelope.

Let `B_i` and `C_i` be the positive and negative quantile bands on a
long-ray direction packet `i`, chosen as in Corollary 2.2 of
`euclidean_expander_exclusion.md`.  They have conditional masses at least
`delta/4`.  For `b in B_i`, `c in C_j`, put

\[
        s_f(b,c)=|b-c|-\bigl(f(b)-f(c)\bigr)\ge0.         \tag{5.1}
\]

**Lemma 5.1 (finite offsets cannot switch separated core packets).**  On a
common type family of rays of conditional scale at least `s`, one has

\[
 s_f(b,c)\ge {\delta s\over64}|N_i-N_j|^2               \tag{5.2}
\]

for all `b in B_i,c in C_j`.  Hence the shifted endpoint prescription

\[
                      u_h(x)=f(x)-h_i
 \quad(x\in B_i\cup C_i)                                \tag{5.3}
\]

satisfies every positive-to-negative Lipschitz inequality between packets
`i,j` whenever

\[
 |h_i-h_j|\le {\delta s\over64}|N_i-N_j|^2.              \tag{5.4}
\]

Equivalently, for distance-cone charts anchored on the opposite core band,
no chart from packet `j` can undercut packet `i` at amplitude below the
right side of (5.4).

**Proof.**  Equation (5.2) is precisely the two-ended interior calibration
bound: both chosen endpoints retain outward calibrated continuation of
length `delta s/4`.  Now

\[
 u_h(b)-u_h(c)
 =f(b)-f(c)-(h_i-h_j).
\]

The two choices of orientation reduce the Lipschitz condition to
`|h_i-h_j|<=s_f(b,c)`.  Use (5.2).  For a distance cone
`ell_j(b)=inf_{c in C_j}(f(c)+|b-c|)`, the chart gap is the infimum of
`s_f(b,c)`, giving the last statement.  QED.

If direction packets are chordally separated by `gamma`, every offset of
range less than `delta s gamma^2/64` is therefore raywise exact on the
opposite-sign long core.  Any defect in (1.11) must come from

* same-sign gluing constraints;
* the taper between the selected core and the rest of a ray; or
* medial/focal endpoint regions where outward continuation has been lost.

If those regions have total mass `beta` and a lattice envelope has no other
defect, then `|r|<=R` gives the conditional criterion

\[
 {\delta^3\over6a_1s}\int c^2d\eta>2\beta R            \tag{5.5}
\]

for a strict competitor in the quadratic amplitude range.  For a height of
size `R=theta s`, both sides of (5.5) are of order `s`; there is no packet
entropy factor.  Thus (5.5) closes only if the endpoint/taper mass has an
additional smallness property.  Neither global log-concavity nor isotropy
currently supplies such a bound.

This is the finite-amplitude version of the endpoint obstruction: interior
Euclidean slack is favorable, but an abstract endpoint expander can still
pay the whole Bregman gain.

## 6. A finite rematching competitor

There is a clean setting in which endpoint rerouting does close.  It also
shows exactly which hypothesis is absent from the general ray problem.

**Lemma 6.1 (Hall-rematching competitor).**  Let `V` be a set of `N` points,
`N` even, with the uniform probability, and suppose

\[
                         |x-y|\ge d                       \tag{6.1}
\]

for all distinct points.  Let `f_0` take the values `+d/2` and `-d/2` on
two sets of `N/2` points.  Then `f_0` is one-Lipschitz and
`J(f_0)=d/2`.

Fix `Delta>0` and let `G_Delta` join pairs at distance `<=d+Delta`.  If a
nonempty set `A` satisfies

\[
                    |A\cup N_{G_\Delta}(A)|\le N/2,      \tag{6.2}
\]

then there is an explicitly defined half-mass cut and an explicit
one-Lipschitz distance competitor `g` obeying

\[
 \boxed{\qquad
       \mathcal J(g)\ge {d\over2}+{|A|\over N}\Delta
       >\mathcal J(f_0).
 \qquad}                                                 \tag{6.3}
\]

**Proof.**  Choose a set `S` of `N/2` points containing
`A union N_{G_Delta}(A)`.  Every coupling from the uniform law on `S` to
the uniform law on `S^c` sends the conditional mass `2|A|/N` out of `A` a
distance at least `d+Delta`; all remaining mass moves at least `d`.
Therefore

\[
 W_1(\mu_S,\mu_{S^c})
 \ge d+{2|A|\over N}\Delta.                            \tag{6.4}
\]

In fact no transport optimizer has to be chosen: take

\[
                         g(x)=d(x,S^c).                  \tag{6.4a}
\]

This is one-Lipschitz, vanishes on `S^c`, is at least `d` on `S`, and is at
least `d+Delta` on `A`.  Testing its centered absolute moment against the
balanced sign of `S` gives

\[
 \mathcal J(g)\ge {1\over2}
    \left(d+{2|A|\over N}\Delta\right).                 \tag{6.5}
\]

which is (6.3).  QED.

At `Delta=0`, the same proof gives an exact matching condition.

**Corollary 6.2 (universal equality matching).**  If `f_0` is a global
extremizer in Lemma 6.1, then for every bisection `S`, the bipartite graph of
pairs at distance exactly `d` between `S` and `S^c` has a perfect matching.
In particular, for every nonempty `A`,

\[
                    |A\cup N_{G_0}(A)|>N/2.             \tag{6.6}
\]

**Proof.**  The cut formula and extremality give
`W_1(mu_S,mu_Sc)<=d`.  Assumption (6.1) gives the reverse inequality, so an
optimal coupling is supported on the distance-`d` relation.  Equal atomic
masses and the Birkhoff theorem give a perfect matching.  If (6.6) failed,
put `A union N(A)` inside one bisection; then `A` would have no equality edge
leaving that side.  QED.

The use of exactly balanced cuts in Lemma 6.1 loses information.  Allowing
the new cut to be slightly unbalanced gives the sharp quantitative expansion
alternative.

**Proposition 6.2a (all-cut neighborhood expansion).**  Under the hypotheses
of Lemma 6.1, assume `f_0` is a global extremizer.  For `Delta in (0,d]` and
nonempty `A subset V`, let

\[
       m=|A\cup N_{G_\Delta}(A)|,\qquad a=|A|.            \tag{6.6a}
\]

If `m<N`, then necessarily `m>N/2` and, with

\[
 t={m\over N}-{1\over2},\qquad
 \lambda={\Delta\over d},\qquad \theta={a\over N},
\]

one has

\[
 \boxed{\quad
 t^2\ge\lambda\theta\left({1\over2}-t\right).
 \quad}                                                   \tag{6.6b}
\]

Equivalently,

\[
 {m\over N}\ge {1\over2}
 +{\sqrt{\lambda^2\theta^2+2\lambda\theta}
          -\lambda\theta\over2}
 \ge {1\over2}+{1\over3}\sqrt{\lambda\theta}.          \tag{6.6c}
\]

The last inequality uses only `lambda theta<=1`.

**Proof.**  The case `m<=N/2` is excluded by Lemma 6.1.  If `m>N/2`, take
`S=A union N(A)`, so `|S|=m`.  The explicit function
`g_S(x)=d(x,S^c)` is at least `d+Delta` on `A`, at least `d` on the rest of
`S`, and zero on `S^c`.  Equivalently (or by the same transport lower bound),

\[
 \mathbb E_{\mu_S}g_S-\mathbb E_{\mu_{S^c}}g_S
 \ge d+{a\over m}\Delta.                                \tag{6.6d}
\]

The cut formula and global extremality give

\[
 2{m(N-m)\over N^2}
       \left(d+{a\over m}\Delta\right)\le {d\over2}.     \tag{6.6e}
\]

After cancellation this is

\[
              4a(N-m)\Delta\le d(N-2m)^2,               \tag{6.6f}
\]

which is (6.6b).  Solving the quadratic gives the first inequality in
(6.6c).  For `z=lambda theta<=1`,

\[
 {\sqrt{z^2+2z}-z\over2}\ge {\sqrt z\over3},
\]

which proves the second.  QED.

There is a useful colored interpretation.  Suppose `V=B disjoint union C`
with `|B|=|C|=N/2`.  For `A subset B`, even if every point of `B` is a near
neighbor, the same-sign part of its closed neighborhood has mass at most one
half.  Proposition 6.2a therefore forces

\[
 { |N_{G_\Delta}(A)\cap C|\over N}
 \ge {1\over3}\sqrt{ {\Delta\over d}{|A|\over N}}.       \tag{6.6g}
\]

Thus same-sign relays cannot absorb the entire finite rematching condition:
every packet must also have a square-root amount of opposite-sign near
mass.  For a singleton this is `Omega(sqrt(N))` opposite-sign neighbors
when `Delta` is comparable to `d`.

Here is the promised exact cross-distance conclusion.

**Corollary 6.3 (separated two-cloud model is Clifford).**  Let

\[
 B=\{b_1,\ldots,b_M\},\qquad C=\{c_1,\ldots,c_M\},
\]

and suppose, for some `Delta_0>0`,

\[
 \begin{aligned}
 |b_i-c_j|&\ge d &&\text{for all }i,j,\\
 |b_i-b_j|,|c_i-c_j|&\ge d+\Delta_0 &&\text{for }i\ne j.
 \end{aligned}                                           \tag{6.7}
\]

If the two-level function `f_0(B)=d/2`, `f_0(C)=-d/2` is a global
extremizer for the uniform law on `B union C`, then

\[
                         |b_i-c_j|=d                      \tag{6.8}
\]

for every `i,j`.  Hence `B,C` lie on spheres in mutually orthogonal affine
subspaces as in the cross-distance classification lemma.

**Proof.**  Apply (6.6) to the singleton `A={b_i}` in the exact-distance
graph `G_0`.  Its equality neighbors other than itself can only lie in `C`.
Since `|A union N(A)|>M`, all `M` points of `C` must be equality neighbors
of `b_i`.  Repeat for every `i`.  The Euclidean classification of two sets
with constant cross distance gives the final assertion.  QED.

Lemma 6.1 is a genuine finite one-Lipschitz competitor, not a second
variation.  It also explains the remaining difficulty.  Same-sign endpoint
packets may be closer than the calibrated positive-to-negative length and
can act as cheap relays in every new half-cut.  High medial conductance is
compatible with exactly such relay capacity.  The long-core slack lemma
controls cross-sign alternative chords, but it does not give the same-sign
separation in (6.7) at focal endpoints.

### 6.1 Combining expansion with long-core angular slack

Proposition 6.2a gives more than the qualitative phrase “endpoint mass.”
Here is the exact alternative obtained when the flattened packet model also
inherits the long-core slack.

**Hypothesis warning.**  Everything through (6.13) in this subsection is
conditional on the global packet separation (6.1).  The long-core lemma
supplies only the cross-sign excess (6.9); it does not supply (6.1) for
same-sign endpoints.  Thus the following is a sharp implication for the
flattened model, not a theorem about the unmodified endpoint bands of the
true extremizer.

Assume, conditionally, that the packet representatives satisfy the global
`d`-separation hypothesis (6.1), that the positive and negative endpoint
clouds have equal total packet mass, that their calibrated value gap `d`
lies in `[c_ds,C_ds]`, and that a cross pair whose normals differ by at least
`gamma` has distance at least

\[
                         d+\kappa s\gamma^2,             \tag{6.9}
\]

where one may take `kappa=delta/64` in the common-type core bands.  Put
`Delta=kappa s gamma^2`.  If a positive direction packet `A` has normalized
total endpoint mass `theta`, (6.6g) says that its opposite-sign angular
`gamma`-neighborhood must have mass at least

\[
 \boxed{\quad
       c_*\gamma\sqrt\theta,
 \qquad c_*={1\over3}\sqrt{\kappa/C_d}.
 \quad}                                                   \tag{6.10}
\]

Otherwise Lemma 6.1, with the slightly unbalanced refinement in Proposition
6.2a, supplies a strict finite competitor.

For a direction cap `P` of quotient mass `m(P)`, every cross-near neighbor
in (6.9) belongs to the chordal `gamma`-enlargement `P^gamma`.  Ignoring only
the fixed band-mass constants, the necessary cap-growth inequality is

\[
 \boxed{\qquad
        m(P^\gamma)\ge c_*\gamma\sqrt{m(P)}.
 \qquad}                                                  \tag{6.11}
\]

This is the strongest graph consequence of finite rematching plus the
interior cross-slack estimate.  Iterating nested enlargements gives, for
`P_k=P^{k gamma}` as long as the radii remain in the controlled range,

\[
 m(P_k)\ge(c_*\gamma)^{2(1-2^{-k})}
                 m(P)^{2^{-k}}.                         \tag{6.12}
\]

Thus a nontrivial cap is driven toward mass of order `gamma^2` after
`O(log log(1/m(P)))` steps.  In the exact finite graph, (6.6c) also says
that every threshold graph has minimum degree greater than `N/2`, and hence
diameter at most two.  Every missing near edge has a same-sign or cross-sign
two-step relay.  More quantitatively, applying (6.6c) to singletons gives

\[
 \deg_{G_\Delta}(v)
 \ge {N\over2}+{1\over3}\sqrt{N\Delta/d}-1,              \tag{6.12a}
\]

so any two vertices have at least
`(2/3)sqrt(N Delta/d)-2` common closed-neighborhood vertices.  A missing
cross edge must therefore possess that many near relays.  If same-sign
distances exceed the threshold, the relay is impossible and Corollary 6.3
upgrades the relation to complete cross equality.  If same-sign relays are
allowed, iteration alone yields (6.11), not complete cross equality.

The exponential coherent-packet estimate does not contradict (6.11) by
itself.  It is an upper bound

\[
                 m(\hbox{unit cap})\le C e^{-cs},        \tag{6.13}
\]

whereas (6.11) needs a lower bound on the starting mass of a smaller cap.
The direction law may split every unit cap into pieces far smaller than
`e^{-2cs}`; then the square-root lower bound is compatible with (6.13) for
all finitely many fixed-radius enlargements.  A covering lower bound for a
`gamma`-cap is of order `gamma^n`, which reintroduces dimension.  Choosing
`gamma` small and iterating (6.12) recovers at best a condition involving
`log n` or `log log n`; it is not dimension-free when `n` is arbitrarily
large relative to `s`.

Equivalently, for equal packet masses `1/N`, (6.10) demands
`Omega(gamma sqrt(N))` opposite-sign neighbors.  The cap estimate permits
as many as `N e^{-cs}` neighbors.  These bounds conflict only when an
independent upper bound on `N` (or a lower bound on packet mass) is supplied.
The entropy argument supplies the opposite inequality `N>=e^{cs}`.  Hence
packet entropy plus long-core slack does **not** ensure the neighborhood
hypothesis of Lemma 6.1.

The precise surviving expansion alternative is therefore (6.11): a
square-root angular neighborhood expansion at the endpoint law.  Proving
that a globally log-concave extremal endpoint law cannot satisfy (6.11)
without becoming approximately complete cross-distance/Clifford is new
content.  It is not a consequence of the current entropy bound.  Moreover,
the actual ray bands do not presently satisfy (6.1): same-sign endpoint
points can be arbitrarily close.  Quotienting those close relays without
losing endpoint mass is an additional prerequisite for applying (6.11) to
the true extremizer.

## 7. A sharp finite Euclidean defect model

Exact extremality alone cannot be upgraded from an almost-complete equality
relation to pointwise equality.  The following Euclidean atomic example is
sharp for that issue.

**Proposition 7.1 (one exceptional cross distance is invisible).**  Let
`N>=6` be even.  There are `N` points in `R^{N-1}` such that every pairwise
distance is one except one pair, whose distance is any prescribed

\[
              1<D\le 2\sqrt{ {N-1\over2(N-2)}}.          \tag{7.1}
\]

For the uniform law on these points,

\[
 \sup_{\operatorname {Lip}(g)\le1}\mathcal J(g)={1\over2}. \tag{7.2}
\]

In particular a balanced two-level extremizer can put the exceptional pair
on opposite sides even though its cross distance is strictly larger than
the calibrated value one.

**Proof.**  Take `N-2` vertices `q_3,...,q_N` of a regular simplex of side
one, centered in an `(N-3)`-dimensional space `W`.  Their common squared
radius is

\[
                         r_0^2={N-3\over2(N-2)}.          \tag{7.3}
\]

In a two-dimensional subspace orthogonal to `W`, choose `x_1,x_2` on the
sphere of squared radius

\[
                         R^2=1-r_0^2={N-1\over2(N-2)}     \tag{7.4}
\]

with chord length `D`.  Then `|x_i-q_j|=1`, all `q`-distances are one,
and only `|x_1-x_2|=D` differs.

The lower bound in (7.2) is any balanced assignment of the values
`+1/2,-1/2`.  For the upper bound use the exact cut formula.  If a cut has
`k` atoms with `2<=k<=N-2`, its two conditional laws admit a coupling which
avoids the single exceptional transport cell; the complete bipartite
transportation table with one cell removed is feasible because

\[
                         {1\over k}+{1\over N-k}\le1.     \tag{7.5}
\]

Thus its `W_1` distance is one and its cut objective is at most `1/2`.
For a singleton cut, the largest possible `W_1` is

\[
                         1+{D-1\over N-1},               \tag{7.6}
\]

and the cut objective is

\[
 {2(N-1)\over N^2}\left(1+{D-1\over N-1}\right)
 ={2(N+D-2)\over N^2}\le {1\over2}.                     \tag{7.7}
\]

The complementary cuts are identical.  This proves (7.2).  QED.

The exceptional incidence has product mass of order `N^{-2}`.  Therefore a
robust rematching theorem can at best conclude weighted approximate
cross-distance equality.  Proposition 7.1 is not log-concave; convexifying
its support destroys the atomic endpoint masses.  It is a sharp warning
about the strength of exact extremality, not a counterexample to the desired
Euclidean/log-concave theorem.

## 8. The Clifford endpoint pattern is harmless

Assume now that complete or approximately complete rematching has produced
the orthogonal-spherical alternative.  The exclusion is direct and uses no
quotient Poincare inequality.

**Lemma 8.0 (constant cross-distance classification).**  If nonempty
Euclidean sets `B,C` satisfy `|b-c|=d` for every `b in B,c in C`, then there
are mutually orthogonal linear subspaces `U,V,W`, points `b_0,c_0` with
`b_0-c_0 in W`, and radii `r_B,r_C>=0` such that

\[
 B\subset b_0+r_BS(U),\qquad C\subset c_0+r_CS(V),
 \qquad r_B^2+r_C^2+|b_0-c_0|^2=d^2.                    \tag{8.0}
\]

**Proof.**  Put `U=span(B-B)` and `V=span(C-C)`.  Subtracting the four
squared cross-distance identities for `b,b' in B` and `c,c' in C` gives
`<b-b',c-c'>=0`, so `U perpendicular V`.  Choose the closest affine fibers
`b_0+U` and `c_0+V`; their displacement is perpendicular to both.  Writing
`b=b_0+u`, `c=c_0+v`, the identity becomes

\[
 d^2=|b_0-c_0|^2+|u|^2+|v|^2.
\]

Fixing first `v` and then `u` shows that the two radii are constant.  QED.

Let `U,V,W` be mutually orthogonal, let `b_0-c_0=w in W`, and suppose the
positive and negative endpoint models are

\[
 B\subset b_0+r_BS(U),\qquad
 C\subset c_0+r_CS(V),
\]

with

\[
                         d^2=r_B^2+r_C^2+|w|^2.          \tag{8.1}
\]

When `w!=0`, put `e=w/|w|`; omit the third coordinate when `w=0`.  Define

\[
 T(x)=\left(
  |P_U(x-b_0)|,
  |P_V(x-c_0)|,
  \langle x-c_0,e\rangle
 \right).                                                \tag{8.2}
\]

The orthogonality of the underlying projections makes `T` one-Lipschitz.
On `B,C` its values are

\[
             t_B=(r_B,0,|w|),\qquad t_C=(0,r_C,0),
             \qquad |t_B-t_C|=d.                         \tag{8.3}
\]

**Lemma 8.1 (projected thin-shell competitor).**  Let `X` be isotropic and
log-concave.  Assume events `B,C` have mass at least `delta` and

\[
 \mathbb E[|T(X)-t_B|\mid B]
 +\mathbb E[|T(X)-t_C|\mid C]\le\varepsilon d           \tag{8.4}
\]

for some `epsilon<1`.  Then

\[
                    d\le
 {\sqrt{2C_{TS}+1}\over\delta(1-\varepsilon)},           \tag{8.5}
\]

where `C_TS` is the universal translated thin-shell constant.  In
particular an exact constant-mass Clifford pattern has universally bounded
ray length.

**Proof.**  Let `ell` be the unit linear functional in the direction
`t_B-t_C` and set

\[
                            g=\ell\circ T.               \tag{8.6}
\]

This is an explicit one-Lipschitz competitor.  The marginals `P_UX` and
`P_VX` are isotropic log-concave in their intrinsic subspaces.  Translated
thin shell gives

\[
 \operatorname {Var}|P_UX-P_Ub_0|,
 \operatorname {Var}|P_VX-P_Vc_0|\le C_{TS},             \tag{8.7}
\]

while isotropy gives variance at most one for the linear `W` coordinate.
By Cauchy--Schwarz in the three coefficients of `ell`,

\[
                         \operatorname {Var}g\le2C_{TS}+1. \tag{8.8}
\]

On the other hand, for `m=Eg`, conditional Jensen, (8.3), and (8.4) give

\[
 \begin{aligned}
 \mathbb E|g-m|
 &\ge\delta\bigl(|\ell(t_B)-m|+|\ell(t_C)-m|
          -\varepsilon d\bigr)\\
 &\ge\delta(1-\varepsilon)d.
 \end{aligned}                                           \tag{8.9}
\]

Combine (8.8)--(8.9).  QED.

This argument applies just as well to translated spheres and degenerate
subspaces.  It is the correct treatment of the sharp Clifford equality
pattern.  It also shows what an approximate inverse theorem must deliver:
feature concentration as in (8.4), not merely a set-theoretic reflection or
a small variance of an unweighted incidence count.

## 9. General signed distances and the precise remaining lemma

Nothing in Sections 1--2 assumes polyhedral regularity.  For an arbitrary
signed-distance extremizer and a finite partition of its ray quotient, choose
a mean-zero height `h` and form any one-Lipschitz minorant or majorant of the
ideal data `u=f-h(Q)`.  For example, on a prescribed core `K`, the McShane
minorant is

\[
               \mathcal M_Ku(x)=\inf_{z\in K}{u(z)+|x-z|\}. \tag{9.1}
\]

On `K` its exact failure to interpolate is

\[
 u(x)-\mathcal M_Ku(x)
 =\sup_{z\in K}{u(x)-u(z)-|x-z|\}_+.                   \tag{9.2}
\]

Thus (9.2) is the nonsmooth analogue of the polyhedral defect (3.3): height
differences are truncated by the actual calibration slack.  Sup envelopes
give the opposite orientation.  Once the construction is defined globally,
Proposition 1.3 converts its `L^1` interpolation defect into the exact
criterion for a strict competitor.

The long-core lemma supplies a linear-in-`s` lower bound for every
cross-packet slack in (9.2).  What is missing is a bound on the same-sign and
endpoint/focal terms.  A sufficient closing statement is the following.

> **Finite endpoint-defect lemma.**  On a fixed positive quotient-mass
> family of balanced scale-`s` rays of a true extremizer, with unit-separated
> direction packets, there is a mean-zero finite height `h` with
> `osc(h)<=c s` and a globally one-Lipschitz lattice/McShane competitor
> `g=f-h(Q)+r` such that
> \[
>             2\int|r|d\mu<\int B_y(h_y)d\eta(y),
> \]
> unless the endpoint laws satisfy the robust orthogonal-spherical condition
> (8.4), or the rays are approximately parallel or concurrent.

Proving this lemma would complete the finite competitor route.  The exact
formulas above show that neither an infinitesimal graph spectral gap nor the
mere fact that there are exponentially many packets proves it.  The needed
new input is a Euclidean, positive-density restriction on all-amplitude
endpoint rematching.  No quotient Poincare or KLS assertion was used in this
reduction.
