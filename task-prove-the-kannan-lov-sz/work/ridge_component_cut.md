# Ridge sweep cuts, physical componentization, and the interior-flow branch

## 0. Outcome

There is a complete abstract sweep theorem for a weighted ridge graph, and
there is an exact `BV` formula for the cost of turning a boundary-phase cut
into actual set components.  They do **not** combine into the proposed
two-way alternative.  The missing quantity is an interior filling capacity.

The expansion estimate used below is

\[
 {\lambda A\over 8}\left(1-{1\over R}\right).        \tag{0.1}
\]

More precisely:

1. If the admissible ridge graph has spectral gap `lambda`, (0.1) is the
   existing simultaneous-bevel saving.  If its weighted degrees obey
   `d_i<=D a_i`, a graph sweep gives a cut of conductance at most
   `sqrt(2D lambda)`.
2. For a Caccioppoli partition of a finite-perimeter set `E`, the exact
   perimeter-additivity error is twice the weighted area of the **interior**
   separating interface.  A codimension-two ridge cut times a bevel width
   does not control this filling area.  Long rectangles and smooth balls
   give sharp counterexamples.
3. The correct quantitative object is a terminal-capacitated physical
   min-cut.  It permits either discarding selected boundary terminals or
   inserting an interior separator.  Its exact dual is a divergence-free
   flow bounded by the physical density.  Small capacity gives the desired
   approximate components.  Large capacity gives, by the decomposition of
   acyclic normal currents, a positive measure of curves joining the two
   terminal phase groups.
4. The large-flow branch is genuine.  The curves need not be straight,
   normal-calibrated rays, or bounded-overlap chords.  The terminal flow may
   also select a low-rank submeasure of a high-rank phase law.  A further
   reach/cross-distance or heat-extremality theorem is required to turn it
   into an affine, radial, or orthogonal-radial configuration.
5. A selected-measure Poincare/Cheeger dichotomy does not bypass this third
   branch.  If `d nu` is proportional to
   `omega |nabla F| d mu`, its cut weight is
   `omega |nabla F| rho`, whereas physical componentization is paid by
   `rho`.  An explicit heat-smoothed Gaussian slab has zero selected
   Cheeger cost at its central critical plane and a strictly positive
   physical filling cost there.

Thus the graph algebra and the component-slope argument are both valid, but
the requested ridge-cut-to-component lemma is false for regular as well as
`BV` boundaries.  The exact corrected alternative has three branches:

\[
 \boxed{\text{bevel expansion}\quad\hbox{or}\quad
        \text{cheap physical min-cut}\quad\hbox{or}\quad
        \text{large density-bounded interior flow}.}           \tag{0.2}
\]

All perimeter statements below are on the affine support of the measure.
The density is denoted by `rho=e^{-V}` and is assumed positive and
continuous in the local statements.  Approximation gives the same `BV`
identities for a locally integrable log-concave density.

## 1. The exact weighted graph sweep theorem

Let `a_i>0`, `A=sum_i a_i`, and let `w_ij=w_ji>=0`.  Put

\[
 d_i=\sum_{j\ne i}w_{ij},\qquad
 {\cal E}(h)=\sum_{i<j}w_{ij}(h_i-h_j)^2,            \tag{1.1}
\]

and suppose

\[
                         d_i\le D a_i.               \tag{1.2}
\]

The generalized spectral gap and the one-sided conductance are

\[
 \lambda=\inf_{h\not\equiv const}
 {\mathcal E(h)\over\sum_i a_i(h_i-\bar h)^2},
 \qquad
 \Phi=\min_{0<a(S)\le A/2}{w(S,S^c)\over a(S)},      \tag{1.3}
\]

where `bar h=A^{-1}sum_i a_i h_i`.

**Lemma 1.1 (weighted sweep, with constants).**  Under (1.2),

\[
                  {\Phi^2\over2D}\le\lambda\le2\Phi.
                                                               \tag{1.4}
\]

In particular, some threshold sweep cut satisfies

\[
             w(S,S^c)\le\sqrt{2D\lambda}\,a(S),
             \qquad 0<a(S)\le A/2.                  \tag{1.5}
\]

**Proof.**  Testing (1.3) with `1_S` gives

\[
 \lambda\le {w(S,S^c)\over a(S)a(S^c)/A}
          \le {2w(S,S^c)\over a(S)},                 \tag{1.6}
\]

which proves the upper bound in (1.4).

For the other direction, choose a weighted median `m` of a real vector
`h`.  Both supports of

\[
 g_+=(h-m)_+,\qquad g_-=(m-h)_+                      \tag{1.7}
\]

have `a`-mass at most `A/2`.  For either nonnegative `g` with such support,
layer cake gives

\[
 \begin{aligned}
 \sum_{i<j}w_{ij}|g_i^2-g_j^2|
 &=\int_0^\infty
   2t\,w(\{g>t\},\{g\le t\})\,dt\\
 &\ge \Phi\sum_i a_i g_i^2.
 \end{aligned}                                       \tag{1.8}
\]

Cauchy--Schwarz and (1.2) give

\[
 \begin{aligned}
 \left(\sum_{i<j}w_{ij}|g_i^2-g_j^2|\right)^2
 &\le {\cal E}(g)
       \sum_{i<j}w_{ij}(g_i+g_j)^2\\
 &\le 2D\, {\cal E}(g)\sum_i a_i g_i^2.
 \end{aligned}                                       \tag{1.9}
\]

Thus `Phi^2 sum a_i g_i^2<=2D E(g)`.  Truncation is a contraction on every
edge, so

\[
 {\cal E}(g_+)+{\cal E}(g_-)\le {\cal E}(h),        \tag{1.10}
\]

whereas

\[
 \sum_i a_i(g_{+,i}^2+g_{-,i}^2)
 =\sum_i a_i(h_i-m)^2
 \ge\sum_i a_i(h_i-\bar h)^2.                       \tag{1.11}
\]

Apply (1.9) to both signs, add, and minimize over `h`.  This proves the
lower bound in (1.4).  The layer-cake proof also shows that one of the
threshold sets realizes (1.5).  QED.

### 1.1 Bevel expansion and the audited constants

Let boundary packets have selected areas `a_i`, representative unit normals
`n_i`, and

\[
 M=\sum_i a_i n_i n_i^T,\qquad
 R={A\over\|M\|_{op}}.                               \tag{1.12}
\]

Assume that `w_ij` are simultaneously admissible bevel conductances, in the
precise sense that an exactly mass-corrected finite splice saves at least

\[
              {1\over8}\sum_{i<j}w_{ij}|n_i-n_j|^2.  \tag{1.13}
\]

Applying the spectral inequality to all coordinate functions of `n_i`
and using

\[
 \left|A^{-1}\sum_i a_i n_i\right|^2\le {1\over R}  \tag{1.14}
\]

gives exactly

\[
 \operatorname{saving}\ge
 {\lambda A\over8}\left(1-{1\over R}\right).        \tag{1.15}
\]

At the conservative values `A=.004p`, `R=17`, and total permissible saving
`epsilon_*p=6.02*10^{-5}p`, absence of a contradiction implies only

\[
 \lambda\le
 {8\epsilon_*\over .004(1-1/17)}
 =0.127925.                                          \tag{1.16}
\]

Even with the optimal degree normalization `D=1`, Lemma 1.1 then supplies
only

\[
                         \Phi\le0.505817.            \tag{1.17}
\]

Using the sharper core values `A=.004489p` and `R=17.16` changes these to

\[
                         \lambda\le0.113924,\qquad
                         \Phi\le0.477333.            \tag{1.18}
\]

These are not rounding issues.  The synchronized-softmax argument permits
an approximate componentization error `e` only while

\[
 {2(\epsilon_*p+e)\over A}<{1\over17}.              \tag{1.19}
\]

For `A=.004489p`, this requires

\[
 \boxed{
 e< {A\over34}-\epsilon_*p
   =7.18294\,10^{-5}p,\qquad {e\over A}<0.016002.}   \tag{1.20}
\]

If `e` is twice an interior filling area, the allowed filling is only

\[
                         3.59147\,10^{-5}p.          \tag{1.21}
\]

Thus the graph sweep loss in (1.18) is almost thirty times too large even
before any geometric lifting loss.

There is also an abstract middle regime.  Take eighteen equal vertices,
`n_i=e_i`, and the complete graph with edge weight chosen so that
`lambda=.05`.  Then `R=18`, `D<.05`, and every balanced cut has
conductance at least `.025`, while (1.15) saves only

\[
 {0.05\over8}\left(1-{1\over18}\right)A
 < 2.66\,10^{-5}p.                                  \tag{1.22}
\]

Thus neither the audited bevel threshold nor the `.016002` component
threshold holds.  Any closing argument needs additional physical structure,
not a sharper invocation of graph Cheeger alone.

## 2. The exact `BV` componentization cost

Let `E` be a finite-perimeter set and let
`(E_1,...,E_m)` be a Caccioppoli partition of `E`: the sets are pairwise
disjoint modulo null sets, their union is `E`, and their perimeters are
locally summable.  The structure theorem for Caccioppoli partitions, applied
to the continuous weight `rho`, gives

\[
 \boxed{
 \sum_{i=1}^m P_\rho(E_i)
 =P_\rho(E)+2\sum_{i<j}
   \int_{\partial^*E_i\cap\partial^*E_j\cap E^{(1)}}
          \rho\,d\mathcal H^{n-1}.}                 \tag{2.1}
\]

For a binary partition `E=E_1 disjoint union E_2`, write

\[
 \operatorname{Fill}_\rho(E_1,E_2;E)
 =\int_{\partial^*E_1\cap\partial^*E_2\cap E^{(1)}}
       \rho\,d\mathcal H^{n-1}.                     \tag{2.2}
\]

Then (2.1) is

\[
                  P_\rho(E_1)+P_\rho(E_2)
                  =P_\rho(E)+2\operatorname{Fill}_\rho.        \tag{2.3}
\]

**Proof of (2.1).**  At `H^{n-1}`-almost every point of a Caccioppoli
partition interface, exactly two cells have density one on opposite sides
of one approximate tangent hyperplane.  An interface inside `E^{(1)}` is
therefore counted once in each of the two cell perimeters and not in
`P(E)`.  At `partial^*E`, exactly one cell has the trace of `E`; summing the
cell perimeter measures counts that point once, exactly as `P(E)` does.
The exceptional set is `H^{n-1}`-null.  Integrating the pointwise
multiplicities against `rho` proves (2.1).  QED.

Equation (2.3) identifies the dimensional mismatch.  A ridge conductance is
of the form

\[
       \ell\int_\Gamma\rho\,d\mathcal H^{n-2},       \tag{2.4}
\]

whereas componentization costs the area of an `(n-1)`-dimensional filling
of `Gamma` inside `E`.  There is no universal inequality from (2.4) to
(2.2).

### 2.1 A convex rectangular counterexample

Work first with constant density.  Let

\[
                     E_L=(-L,L)\times(-1,1),         \tag{2.5}
\]

and color the open top edge phase `+` and the other three open edges phase
`-`.  The phase interface consists of the two top corners.  A ridge proxy of
width `ell<=1` has weight `2ell`.

Suppose a binary Caccioppoli partition has trace one for `E_1` on the top
edge and trace zero on the bottom and side edges.  For almost every
`x in(-L,L)`, the one-dimensional slice
`y mapsto 1_{E_1}(x,y)` has top trace one and has to change value before it
can meet the zero exterior trace.  Its variation in `(-1,1)` is at least
one.  The slicing theorem for `BV` functions gives

\[
        \operatorname{Fill}(E_1,E_2;E_L)
        =P(E_1;E_L)\ge\int_{-L}^L1\,dx=2L.          \tag{2.6}
\]

Curves arbitrarily close to the top edge and joining the two top corners
show that the infimum is `2L`.  Consequently

\[
 {P(E_1)+P(E_2)-P(E_L)\over 2\ell}
 \ge {4L\over2\ell}={2L\over\ell}\longrightarrow\infty.       \tag{2.7}
\]

The example sits inside a log-concave probability by taking the uniform
law on a slightly larger rectangle.  The constant density multiplies both
sides of every formula and cancels in (2.7).

This is also the clean distinction between a ridge cut and an interior
cut.  The facet graph sees two zero-dimensional corner contacts.  The body
can be separated only by a horizontal-scale chord.

### 2.2 A smooth high-rank counterexample

Let `E=B_R` in `R^n`, again with constant density in a larger convex
support, and color the northern and southern hemispheres differently.  The
interface is the equator.  Its ridge proxy is

\[
 \ell\mathcal H^{n-2}(S_R^{n-2})
 =\ell(n-1)\kappa_{n-1}R^{n-2}.                     \tag{2.8}
\]

Any oriented interior separator with equatorial boundary has area at least
that of the equatorial disk,

\[
                         \kappa_{n-1}R^{n-1}.        \tag{2.9}
\]

Indeed the constant `(n-1)`-form
`dx_1 wedge ... wedge dx_{n-1}` has comass one.  Its integral over the
separator equals its integral over the disk because their difference is a
cycle in the ball.  This calibrates the disk and proves (2.9).  Hence

\[
 {\operatorname{Fill}\over\hbox{ridge proxy}}
 ={R\over(n-1)\ell}.                                \tag{2.10}
\]

For either hemisphere, the projective normal matrix is already full-rank:

\[
 \int_{S_R^{n-1}\cap\{x_n>0\}}nn^T,d\mathcal H^{n-1}
 ={\mathcal H^{n-1}(S_R^{n-1})\over2n}I.             \tag{2.11}
\]

This follows because `nn^T` is even under the antipodal map.  Thus smooth
regularity and effective normal rank do not repair the lifting lemma.  The
ball is the concurrent-radial branch: all its boundary normal lines meet at
the center.  A generic smooth strictly convex perturbation preserves the
large filling ratio and full normal rank but destroys exact concurrence;
geometric ridge/filling data alone do not classify that perturbation.

## 3. The corrected terminal-capacitated physical cut

The preceding failure has an exact repair, but the repair creates a third
branch.  Let `E` be a bounded Lipschitz set, let `rho` be positive and
continuous on its closure, and let `Gamma_+,Gamma_-` be disjoint Borel
parts of `partial E`.  Let

\[
 d\sigma_+=f_+\rho\,d\mathcal H^{n-1},\qquad
 d\sigma_-=f_-\rho\,d\mathcal H^{n-1},\qquad
 0\le f_\pm\le1,                                     \tag{3.1}
\]

be the selected terminal phase fluxes.

For a finite-perimeter `U subset E`, denote its interior perimeter by
`P_rho(U;E)` and its interior trace on `partial E` by `Tr_E 1_U`.  Define

\[
 \begin{aligned}
 \mathsf C_E(\sigma_+,\sigma_-)
 =\inf_{U\subset E}\bigg\{&P_\rho(U;E)
  +\int_{\Gamma_+}(1-\operatorname{Tr}_E1_U)\,d\sigma_+\\
  &+\int_{\Gamma_-}\operatorname{Tr}_E1_U\,d\sigma_-\bigg\}.
 \end{aligned}                                                   \tag{3.2}
\]

The first term inserts a physical separator.  The second and third terms
pay for selected phase terminals assigned to the wrong side.  The constant
choices `U=empty` and `U=E` show

\[
              \mathsf C_E\le\min\{\sigma_+(\Gamma_+),
                                      \sigma_-(\Gamma_-)\}.     \tag{3.3}
\]

Thus the terminal degree is automatically controlled by selected surface
flux, with constant one.

**Theorem 3.1 (physical max-flow/min-cut).**  The value in (3.2) equals

\[
 \begin{aligned}
 \sup_z\quad &\int_{\Gamma_+}[z,\nu_E]\,d\mathcal H^{n-1}\\
 \text{subject to}\quad
 &\operatorname{div}z=0\quad\hbox{in }E,\\
 &|z(x)|\le\rho(x)\quad\hbox{a.e.},\\
 &0\le[z,\nu_E]\le f_+\rho\quad\hbox{on }\Gamma_+,\\
 &-f_-\rho\le[z,\nu_E]\le0\quad\hbox{on }\Gamma_-,\\
 &[z,\nu_E]=0\quad\hbox{on the rest of }\partial E.
 \end{aligned}                                                   \tag{3.4}
\]

Here `[z,nu_E]` is the weak normal trace.  In the smooth case it is the
ordinary scalar product.  Both the primal and the dual may be relaxed and
their optima are attained.

**Proof.**  First relax the binary variable to `u in BV(E;[0,1])` and put

\[
 \begin{aligned}
 J(u)=&\int_E\rho\,d|Du|
 +\int_{\Gamma_+}(1-\operatorname{Tr}u)d\sigma_+
 +\int_{\Gamma_-}\operatorname{Tr}u d\sigma_- .
 \end{aligned}                                       \tag{3.5}
\]

The `BV` coarea formula and layer cake give

\[
 J(u)=\int_0^1\bigg[
 P_\rho(\{u>t\};E)
 +\sigma_+(\{\operatorname{Tr}u\le t\})
 +\sigma_-(\{\operatorname{Tr}u>t\})\bigg]dt.       \tag{3.6}
\]

Thus a threshold has binary cost no larger than `J(u)`, and the convex and
binary infima agree.

For every feasible `z`, Gauss--Green and `div z=0` give

\[
 \begin{aligned}
 \int_{\Gamma_+}[z,\nu_E]
 ={}&\int_{\Gamma_+}(1-\operatorname{Tr}u)[z,\nu_E]
    +\int_E z\cdot Du\\
   &+\int_{\Gamma_-}\operatorname{Tr}u(-[z,\nu_E])
 \le J(u).
 \end{aligned}                                       \tag{3.7}
\]

This proves weak duality.  Equality is the weighted `BV`
max-flow/min-cut duality.  For completeness, one obtains it by applying the
polar representation of weighted total variation to (3.5), then separating
the epigraph of the resulting convex lower-semicontinuous functional from a
strictly lower horizontal hyperplane.  The separating functional is a field
`z`; integration by parts gives `div z=0`, the polar constraint gives
`|z|<=rho`, and the two boundary support functions give exactly the normal
trace intervals in (3.4).  Conversely these constraints give (3.7), so the
separation value is the displayed supremum.  `BV` compactness, lower
semicontinuity of the total variation with boundary jump penalties, and
weak-* compactness of density-bounded divergence-measure fields give primal
and dual attainment.  Finally (3.6) produces an attaining binary threshold.
QED.

### 3.1 What small capacity gives

Let `U` have cost at most `mathsf C_E+delta`.  Put

\[
 c=P_\rho(U;E),\qquad
 d=\sigma_+(\operatorname{Tr}U=0)
   +\sigma_-(\operatorname{Tr}U=1).                  \tag{3.8}
\]

The partition `E_1=E cap U`, `E_2=E\setminus U` satisfies

\[
 P_\rho(E_1)+P_\rho(E_2)=P_\rho(E)+2c,              \tag{3.9}
\]

and all but terminal flux `d` is assigned to the intended component.  A
single conservative error parameter is therefore

\[
                    e=2(c+d)\le2(\mathsf C_E+\delta).            \tag{3.10}
\]

At the synchronized-softmax budget (1.20), the cheap-cut branch closes only
if

\[
 \boxed{\mathsf C_E<3.59147\,10^{-5}p}              \tag{3.11}
\]

after integration over all active levels and all reuse losses.  This is the
precise physical statement which a low ridge cut would have to imply.  The
rectangle and ball show that it need not.

### 3.2 Calibration and the large-flow branch

Let `z` and `U` be dual and primal optimizers.  Equality in (3.7) gives the
complementary-slackness identities

\[
 z\cdot D1_U=\rho|D1_U|\quad\hbox{on the interior cut},          \tag{3.12}
\]

with the sign fixed by orientation, and

\[
 \begin{array}{ll}
 [z,\nu_E]=f_+\rho
   &\text{where the `+` terminal is discarded},\\
 -[z,\nu_E]=f_-\rho
   &\text{where the `-` terminal is discarded}.
 \end{array}                                           \tag{3.13}
\]

Thus a large min-cut gives a large physical flow, not a hidden component.

Regard `z dx` as a normal one-current.  Delete its cyclic part; this does not
change its boundary flux and can only decrease its pointwise capacity use.
The decomposition theorem for acyclic normal one-currents then supplies a
measure `Pi` on simple rectifiable curves in `bar E` such that

\[
 \begin{aligned}
 z\,dx&=\int [\![\gamma]\!]\,d\Pi(\gamma),\\
 \int\operatorname{length}(\gamma)d\Pi(\gamma)
 &=\int_E|z|dx\le\int_E\rho dx=\mu(E),\\
 \Pi(\hbox{all curves})
 &=\mathsf C_E.
 \end{aligned}                                       \tag{3.14}
\]

The initial and terminal endpoint laws are precisely

\[
 d\beta=(-[z,\nu_E])d\mathcal H^{n-1}\le d\sigma_-,
 \qquad
 d\alpha=[z,\nu_E]d\mathcal H^{n-1}\le d\sigma_+,  \tag{3.15}
\]

and both have mass `mathsf C_E`.

This gives two quantitative facts.

1. The mean curve length is at most `mu(E)/mathsf C_E`.
2. If `M_term=int nn^T d sigma` has operator norm at most `A/R`, then the
   endpoint normal matrix of any dominated flux obeys

\[
 0\preceq M_\alpha\preceq M_{term},\qquad
 {\operatorname{tr}M_\alpha\over\|M_\alpha\|_{op}}
 \ge {\mathsf C_E R\over A}.                         \tag{3.16}
\]

   Thus a flow carrying a fixed fraction of a high-rank terminal law
   retains a corresponding fraction of its effective rank.

The theorem supplies no more.  In particular it does not imply that the
curves are straight, that their tangents agree with the boundary normals
away from their endpoints, or that their straight chords have bounded
overlap.  A divergence-free maximizing field can bend and merge.  If the
capacity is only `.008A`, (3.16) retains only `.008R`, which is not a useful
rank.  A graph decomposition would have to aggregate many terminal-disjoint
flows without reusing sink capacity before (3.16) becomes strong.

A further exact refinement is possible.  If terminal bands within distance
`ell` of the ridge carry at most `q_ridge`, and terminal pairs of opposite
color at Euclidean distance below `ell` outside those bands carry at most
`q_cross`, then at least

\[
              \mathsf C_E-2q_{ridge}-q_{cross}       \tag{3.17}
\]

of the curves have length at least `ell`.  This follows by removing curves
with an endpoint in a ridge band and then using
`length(gamma)>=|gamma(1)-gamma(0)|`.  The residual alternatives are exactly
ridge/end charge, short cross-contact, and long flow.  What is missing is a
theorem turning the last two terms into the already audited
reach/cross-distance rigidity with the required selected flux and bounded
reuse.

## 4. Why selected-measure Cheeger does not control physical filling

Let `F:R^n->[0,1]` be smooth and let

\[
 d\nu={\eta\over a}d\mu,\qquad
 \eta=\omega|\nabla F|,\qquad
 a=\int\eta d\mu.                                   \tag{4.1}
\]

For a smooth spatial cut `U`, its unnormalized selected-measure boundary
cost is

\[
                 \int_{\partial U}\omega|\nabla F|\rho\,d\mathcal H^{n-1}.
                                                               \tag{4.2}
\]

On the other hand, splitting every superlevel
`A_r={F>r}` by the same `U` creates, by (2.3), twice the interior cut

\[
 2\int_0^1\int_{\partial U\cap A_r}\rho\,d\mathcal H^{n-1}dr
 =2\int_{\partial U}F\rho\,d\mathcal H^{n-1}.        \tag{4.3}
\]

There is no pointwise upper bound on

\[
                         {F\over\omega|\nabla F|}.    \tag{4.4}
\]

The good-core condition `|m|/R>=1/2` gives an upper Gaussian bound on
`|nabla F|`; it gives no positive lower bound.  Replacing the Euclidean
metric in the selected Dirichlet form by the metric which makes (4.2) equal
to (4.3) inserts precisely the missing inverse-`|nabla F|` lapse.

### 4.1 An exact heat-smoothed Gaussian slab

The mismatch already occurs for an actual heat-smoothed indicator.  Let
`mu=gamma_n`, choose `L>0`, put

\[
 S=\{|x_1|\le L\},\qquad F=P_s1_S,                  \tag{4.5}
\]

and take `omega=1` (or any positive bounded analytic selector).  Then

\[
 F(x)=\Phi\left({L-x_1\over\sqrt s}\right)
      +\Phi\left({L+x_1\over\sqrt s}\right)-1.      \tag{4.6}
\]

It is even in `x_1`, strictly decreases with `|x_1|`, and

\[
 \partial_1F(x)={1\over\sqrt s}\left[
 \varphi\left({L+x_1\over\sqrt s}\right)
 -\varphi\left({L-x_1\over\sqrt s}\right)
 \right],\qquad \partial_1F(0)=0.                   \tag{4.7}
\]

The selected law `d nu` proportional to
`|partial_1F|d gamma_n` gives equal mass to the two halfspaces.  For
`U={x_1<0}`,

\[
 \nu^+(U)=0,                                        \tag{4.8}
\]

because `|partial_1F(x)|=O(|x_1|)` at the central plane, so the mass of its
one-sided `epsilon`-tube is `O(epsilon^2)`.  Thus the selected-measure
Cheeger constant is zero.  Equivalently, logarithmic transition functions
across the plane have arbitrarily small selected Dirichlet energy; the
selected Poincare constant is infinite.

For every `0<r<F(0)`, however,

\[
                         A_r=\{|x_1|<a_r\}            \tag{4.9}
\]

is a connected slab.  Splitting it by `x_1=0` inserts a Gaussian
hyperplane of perimeter `I_0=(2pi)^{-1/2}` in each component.  Integrated
over the active levels, the additivity error is

\[
                       2F(0)I_0>0.                  \tag{4.10}
\]

Equations (4.8)--(4.10) give an infinite separation between selected
Cheeger cost and physical component cost.  Choosing
`L=Phi^{-1}(3/4)` makes the original slab `S` a half-mass set.  The example
has bounded ambient Poincare constant and is not a KLS counterexample; its
purpose is exact: heat origin, analyticity, centrality, and a positive
selector do not control the lapse in (4.4).

## 5. Model audit

### 5.1 Cube

For an inner cube, the `2n` facet graph is the complete graph with the
opposite matching removed.  With equal packet areas and uniformly allocated
ridge strips, its generalized nonzero spectral gap is
`(2n-2)w/a`, the same scale as its normalized degree.  The projective normal
matrix has rank `n`.  Thus the full available ridge graph is an expander and
the simultaneous bevel branch applies.  Artificially choosing a tiny
subconductance is not legitimate when the remaining disjoint ridge strips
are available.

### 5.2 Simplex

Every pair of facets of a simplex meets in a ridge.  The facet graph is
`K_{n+1}`; with symmetric admissible widths its gap is
`(n+1)w/a`, again comparable to degree.  A regular simplex has a tight-frame
normal law.  Hence the simplex activates the bevel-expansion branch, not the
component branch.

### 5.3 Product exponentials

The upper facets of the exponential maximum box form `K_n`, and all
log-affine normal slopes equal one.  The exact rounded maximum construction
simultaneously bevels all higher-order corner overlaps.  The phase
proportions do not rotate because the slopes agree.  This model passes both
the ridge-expansion and synchronized-slope tests.

### 5.4 Long rectangle

The top-facet sweep has ridge conductance `2ell` and selected area `2L`, so
its ridge conductance tends to zero like `ell/L`.  Its terminal-capacitated
physical min-cut is instead `2L` (up to the common density).  The dual field
is the constant vertical flow, and its decomposition consists of coherent
straight vertical curves.  This is the affine/coherent instance of the
large-flow branch.

### 5.5 Thin neck

For two large lobes joined by a cylindrical neck of radius `r`, the
physical min-cut through the neck is of order
`kappa_{n-1}r^{n-1}`.  A ridge-band proxy at scale `ell` is of order
`ell(n-1)kappa_{n-1}r^{n-2}`.  Their ratio is

\[
                              {r\over(n-1)\ell}.      \tag{5.1}
\]

When the neck radius is at most the admissible surgery scale, the cheap-cut
branch is plausible.  When `r>>ell`, the flow branch carries the cross-section
capacity.  A smooth neck has no literal polyhedral ridge, so a `BV`/focal or
flow formulation is mandatory.

### 5.6 Smooth nonpolyhedral surfaces

The sphere calculation (2.8)--(2.11) has no ridge and has full projective
rank.  Its planar max flow is constant while its boundary normals are radial.
A smooth ellipsoid gives the affine-radial version.  Generic smooth strictly
convex perturbations and nonplanar terminal interfaces give legitimate
least-gradient cuts with nonconstant dual fields.  Max-flow/min-cut alone
does not classify them: the flow curves may bend and the straightened chords
may focus or overlap.  Excluding such arbitrary irreducible fillings from
the large-`K`, heat-extremal regime is a new inverse theorem, not a
consequence of ridge degree, normal rank, or log-concavity alone.

## 6. Exact implication for phase-slope synchronization

Suppose, after integrating over active levels and summing all reuse errors,
the terminal-capacitated cuts produce physical components with total error
`e` as in (3.10).  The exact component excess formula then changes only by
adding `e` to the scalar deficit.  With

\[
 A\ge.004489p,\qquad D_{co}\le6.02\,10^{-5}p,       \tag{6.1}
\]

the fraction of component flux with `kappa_i>=2psi` is at most

\[
                         q\le {2(D_{co}+e)\over A}.  \tag{6.2}
\]

If (1.20) holds, then `q<1/17`.  All remaining slopes lie in
`[psi,2psi)`, so their dimensionless softmax slopes differ by at most

\[
 \Delta c\le\sqrt{s}\,\psi\le2\sqrt\alpha=2\,10^{-5}.           \tag{6.3}
\]

On the central profile interval of length below `10.96`, the exact
likelihood-ratio estimate gives

\[
 \operatorname{Var}Q
 \le2\{e^{\sqrt2(10.96)\Delta c}-1\}^2+8q
 <{8\over17}.                                        \tag{6.4}
\]

Thus the cheap-physical-cut branch really would exclude the cross-level
softmax alternative, with no boundary Poincare inequality.  The failure is
solely the inference

\[
 \text{small ridge sweep}
 \quad\Longrightarrow\quad
 \mathsf C_E<3.59147\,10^{-5}p,                     \tag{6.5}
\]

which is disproved by Sections 2 and 4.

The formal replacement for (6.5) is the three-way alternative (0.2).  A
closing proof must charge the large-flow branch using the special
heat-extremal origin of the terminals.  At minimum it must show that a fixed
fraction of the selected high-rank endpoint flux travels either through
bounded-reuse short cross contacts (which can be finitely spliced) or along
nearly straight long calibrated curves (to which the existing
reach/cross-distance and affine/radial inverse theorems apply).  The physical
max-flow theorem supplies the curves and exact capacities, but not that
last geometric compatibility.
