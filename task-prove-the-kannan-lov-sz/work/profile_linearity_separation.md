# Separated phases, linear isoperimetric profiles, and exponential rigidity

## 0. Verdict

Perimeter additivity for spatially separated phase regions has a sharp scalar
consequence.  If their union is a near-isoperimetric set, then almost all of
their mass lies at volumes where

\[
                             q(v):={\mathcal I(v)\over v}            \tag{0.1}
\]

is close to its value at the total volume.  Exact additivity forces `q` to be
constant on a nontrivial interval.

There is a complete one-dimensional rigidity theorem.  On a log-concave
line, the logarithmic derivative of `q` is exactly the defect from an
exponential tail.  If `q` is constant on an interval, every active tail is
globally exponential up to the corresponding quantile.  If `q` is nearly
constant, then at most logarithmic-volume points the conditional tail is
not close in total variation to an exponential law.

What is not available is the required higher-dimensional equality/stability
theorem with hard convex support: linear profile on an interval would have
to turn separated isoperimetric components into affine exponential factors.
A smooth nested-foliation version follows from the second variation and is
proved below, but arbitrary log-concave measures need not have such a
foliation, and quantitative control of singular/contact sets is missing.
Consequently the low-ridge branch of the fixed physical-splice program is
reduced, not closed.

## 1. The exact component-additivity inequality

Let `mu` be log-concave on its affine support and let `mathcal I` be its
relaxed weighted-`BV` isoperimetric profile.  It is symmetric and concave,
with `mathcal I(0)=0`.  Hence `q(v)=mathcal I(v)/v` is nonincreasing on
`(0,1/2]`.

**Lemma 1.1 (separated-component deficit).**  Let `A_1,...,A_m` be
finite-perimeter sets with pairwise disjoint closures at positive distance,
or, more generally, suppose their canonical representatives satisfy

\[
 P_\mu\left(\bigcup_iA_i\right)=\sum_iP_\mu(A_i).   \tag{1.1}
\]

Put `v_i=mu(A_i)>0`, `v=sum_i v_i<=1/2`, and `A=union_iA_i`.  Then

\[
 \boxed{
 P_\mu(A)-\mathcal I(v)
 \ge\sum_i\{\mathcal I(v_i)-q(v)v_i\}
 =\sum_iv_i\{q(v_i)-q(v)\}.}                       \tag{1.2}
\]

If the left side is at most `epsilon mathcal I(v)`, then, for every
`kappa>0`,

\[
 \boxed{
 \sum_{\{i:q(v_i)>(1+\kappa)q(v)\}}v_i
 \le {\epsilon v\over\kappa}.}                     \tag{1.3}
\]

**Proof.**  Isoperimetry gives `P(A_i)>=mathcal I(v_i)`.  Equation (1.1),
followed by `sum v_i=v`, proves (1.2).  Every index in (1.3) contributes
more than `kappa q(v)v_i` to (1.2), whereas the assumed total deficit is
`epsilon q(v)v`.  This proves (1.3).

**Corollary 1.2 (exact and quantitative linearity).**  If the union is an
exact minimizer, `P_mu(A)=mathcal I(v)`, then every `A_i` is isoperimetric
and

\[
                         q(v_i)=q(v).                \tag{1.4}
\]

Consequently `mathcal I(t)=q(v)t` throughout `[0,v]`.  More quantitatively,
if the deficit is at most
`epsilon mathcal I(v)` and `v_i>=eta v`, then

\[
 q(v)\le q(t)\le(1+\epsilon/\eta)q(v)
 \quad(v_i\le t\le v).                              \tag{1.5}
\]

**Proof.**  Every summand in (1.2) and every individual isoperimetric
deficit is nonnegative.  Equality forces all of them to vanish.  The two
adjacent secant slopes

\[
 {\mathcal I(v_i)-\mathcal I(0)\over v_i}
 ={\mathcal I(v)-\mathcal I(v_i)\over v-v_i}=q(v)
\]

are equal.  Concavity therefore forces `mathcal I(t)=q(v)t` on all of
`[0,v]`.  In the approximate case, the single `i`-th summand is at most
`epsilon q(v)v`, so `q(v_i)-q(v)<=epsilon q(v)/eta`; monotonicity proves
(1.5).

There is also a useful extension below the component mass.  Suppose
`q(v_i)<=(1+delta)q(v)` and `v_i<=(1-eta_0)v`.  Then, for
`0<theta<=1` and `theta v_i<=t<=v_i`,

\[
 \boxed{q(t)\le q(v)
       \left(1+{\delta\over\theta\eta_0}\right).}    \tag{1.6}
\]

Indeed `J(t)=mathcal I(t)-q(v)t` is nonnegative and concave on `[0,v]`,
with `J(0)=J(v)=0` and `J(v_i)<=delta q(v)v_i`.  Concavity applied to the
chord from `(t,J(t))` to `(v,0)` gives

\[
 J(t)\le J(v_i){v-t\over v-v_i}
 \le {\delta q(v)v_i\over\eta_0},
\]

and division by `t>=theta v_i` proves (1.6).

The same proof works for countably many components by monotone convergence.
For phase regions separated only by a small ridge set, (1.1) acquires exactly
the weighted perimeter carried by that ridge.  That is the scalar form of
the bevel-versus-linearity alternative.

## 2. Exact one-dimensional exponential rigidity

Let `rho=e^phi` be a log-concave probability density on an interval, let
`F(x)=int_{-infinity}^x rho`, and consider a range of volumes on which the
left halfline is the active one-dimensional isoperimetric minimizer.  Then

\[
                         \mathcal I(v)=\rho(F^{-1}(v)).               \tag{2.1}
\]

At a differentiability point `x=F^{-1}(v)`, put

\[
 h(v)={\rho(x)\over F(x)}={\mathcal I(v)\over v}=q(v),
 \qquad d(v)=\phi'(x)=\mathcal I'(v).               \tag{2.2}
\]

Concavity of `phi` gives the tangent bound

\[
 \rho(x-u)\le\rho(x)e^{-d(v)u},\qquad u\ge0,         \tag{2.3}
\]

whenever `d(v)>0`.  Integration yields

\[
                         h(v)\ge d(v).               \tag{2.4}
\]

**Theorem 2.1 (linear profile if and only if exponential tail).**  If
`mathcal I(v)=cv` on a nonempty open volume interval of the active left-tail
branch, then for every quantile `x` corresponding to that interval,

\[
 \rho(t)=\rho(x)e^{c(t-x)}\quad\hbox{for almost every }t\le x.        \tag{2.5}
\]

In particular the density has a global one-sided exponential tail, not
merely an exponential slab.

**Proof.**  At almost every `v` in the interval, (2.2) gives `h(v)=d(v)=c`.
The right side of (2.3) integrates to `rho(x)/c=F(x)`, which is equality in
the integrated tangent bound.  The nonnegative pointwise gap in (2.3)
therefore vanishes almost everywhere, proving (2.5).  One such `x` proves
the assertion for its entire left tail; varying `x` gives the stated form.

The reflected assertion holds on an active right-tail branch.  A symmetric
linear profile is therefore the two-sided Laplace equality model, while a
one-sided support gives the ordinary exponential equality model.

### 2.1 A quantitative identity and stability statement

At every differentiability point with `d(v)>0`, define

\[
                         \varepsilon(v)=1-{d(v)\over h(v)}\in[0,1].   \tag{2.6}
\]

Since `q=h` and `q'=(d-h)/v`, one has the exact identity

\[
 \boxed{
 -{d\log q(v)\over d\log v}=\varepsilon(v).}        \tag{2.7}
\]

Consequently, for `0<a<b<=1/2` within one active branch,

\[
 \boxed{
 \int_a^b\varepsilon(v){dv\over v}
       =\log{q(a)\over q(b)}.}                      \tag{2.8}
\]

Thus if `q(a)<=e^eta q(b)`, then for every `gamma>0` the logarithmic measure
of volumes where `varepsilon>gamma` is at most `eta/gamma`.

There is also a distributional stability conclusion.  Let `U=x-X`
conditioned on `{X<=x}`, and let `k_x` be its density.  Equation (2.3) and
`d/h=1-epsilon` give pointwise

\[
 (1-\varepsilon)k_x(u)le d e^{-du},qquad u\ge0,   \tag{2.9}
\]

and the two sides have masses `1-epsilon` and `1`.  Hence

\[
 \boxed{
 \|k_x(u)du-\operatorname {Exp}(d)\|_{TV}
 \le\varepsilon(v).}                               \tag{2.10}
\]

Indeed the `L1` distance between the exponential density and
`(1-epsilon)k_x` is exactly `epsilon`; adding back the missing multiple of
`k_x` gives `L1` distance at most `2epsilon`.

Equations (2.8)--(2.10) are a quantitative rigidity theorem with no
unproved compactness step.

### 2.2 Isotropy forbids a slow exponential branch of fixed mass

The stability statement has a useful scale consequence.  Suppose the
one-dimensional law has variance one, the active tail has mass `v`, and
`varepsilon(v)<=1/50`.  Then

\[
                         \boxed{q(v)\ge {1\over5}\sqrt v.}           \tag{2.11}
\]

Indeed (2.10), after multiplying distance by `d=d(v)`, places the
conditional excess within total variation `1/50` of `Exp(1)`.  Under
`Exp(1)`, the intervals `[0,1/2]` and `[2,3]` have masses larger than `.39`
and `.085`, respectively.  The conditional law therefore gives them masses
larger than `.37` and `.065`.  For two independent conditional samples, the
identity `Var(Z)=E(Z-Z')^2/2` yields

\[
 Var(X\mid X\le x)>.05/d^2.
\]

Total variance gives `1=Var(X)>=v Var(X|X<=x)`, and hence
`d>sqrt(.05v)>sqrt(v)/5`.  Since `q=h>=d`, (2.11) follows.

Thus once a spatial rigidity theorem puts a fixed-mass phase packet onto
one affine marginal, even approximate exponential equality already forces a
universal profile slope.  The difficulty is obtaining that affine marginal,
not the final one-dimensional estimate.

## 3. What smooth `CD(0,infinity)` equality would give

The following proposition records the higher-dimensional calculation under
hypotheses strong enough to justify every variation.  It is not asserted for
an arbitrary nonsmooth log-concave measure.

**Proposition 3.1 (smooth nested-foliation rigidity).**  Let
`dmu=e^{-V}dx` on a smooth convex domain, with `V in C^2` convex.  Suppose
that on an open volume interval `(a,b)` there is a nested smooth family of
isoperimetric regions whose boundaries are generated by normal flow, with no
singular set and either no contact with the support boundary or vanishing
free-boundary error.  If `mathcal I` is affine on `(a,b)`, then every leaf
`Sigma_v` satisfies

\[
                         II_{\Sigma_v}=0,qquad
 \nabla^2V(n_v,n_v)=0.                               \tag{3.1}
\]

Each connected leaf is therefore contained in a hyperplane.  On every slab
swept out by parallel leaves, in coordinates `x=z+tn`,

\[
                         V(z+tn)=W(z)+ct.            \tag{3.2}
\]

**Proof.**  For a smooth isoperimetric leaf, weighted mean curvature is
constant.  Comparing with unit normal flow gives the standard second
variation support inequality

\[
 \mathcal I''(v)
 \le-{1\over\mathcal I(v)^2}
   \int_{\Sigma_v}\left\{|II|^2+\nabla^2V(n,n)\right\}
                       e^{-V}d\mathcal H^{n-1}.      \tag{3.3}
\]

Convexity makes the integrand nonnegative.  Affineness makes the left side
zero in the distributional sense, so the integral and both summands vanish
on every leaf.  A connected totally geodesic Euclidean hypersurface is a
piece of a hyperplane.  Finally a positive-semidefinite Hessian with
`n^T Hess(V)n=0` also has `Hess(V)n=0`.  Along the swept slab this says that
`partial_n V` is constant in both `t` and `z`, proving (3.2).

For a hard convex support, the omitted boundary term in (3.3) is
nonnegative.  Equality would additionally force its vanishing.  What is not
proved is that arbitrary equality components generate a common smooth
normal foliation, or that approximate equality controls singular and
free-boundary terms with dimension-free constants.

There is nevertheless an exact closure in the smooth full-support case.

**Proposition 3.2 (exact full-support closure).**  Assume `V` is smooth,
finite, and convex on all of `R^n`, standard weighted isoperimetric
minimizers exist, and their regular boundaries satisfy the usual codimension
at least eight singular-set theorem.  If an isotropic `mu` has

\[
                         \mathcal I(t)=ct\quad(0<t\le1/2),           \tag{3.4}
\]

then `c>=c_1` for a universal positive constant and hence
`C_P(mu)<=4/c_1^2`.

**Proof.**  At almost every regular volume in (3.4), the profile support
inequality (3.3) and `mathcal I''=0` show that every regular boundary piece
is contained in a hyperplane.  Two nonparallel hyperplanes meet in
codimension two; such a junction cannot be hidden in a singular set of
codimension at least eight.  Distinct complete hyperplanes also cannot be
disjoint unless they are parallel.  Thus, up to null sets, the minimizer is
the inverse image under one linear coordinate `<x,theta>` of a
one-dimensional finite-perimeter set `E`.

The marginal of `<X,theta>` is log-concave and has variance one.  Its
one-dimensional perimeter of `E`, divided by
`min(mu(E),1-mu(E))`, is at least its one-dimensional Cheeger constant,
which is bounded below by a universal `c_1` for every variance-one
log-concave law.  The ratio in (3.4) is `c`, so `c>=c_1`.  Finally
`psi_mu=c` on the small side and Cheeger's inequality gives the asserted
Poincare bound.

The full-support assumption is essential to this proof: in a hard convex
support, nonparallel flat sheets can end on different support faces without
creating an interior codimension-two junction.  The free-boundary equality
term must rule out or split that configuration; this is part of the missing
general theorem.

## 4. Consequence of a genuine affine branch

If the rigidity step produces an actual affine halfspace
`H={x:<x,theta><=t}` whose perimeter-to-small-side-mass ratio is
`(1+eta)psi_mu`, then no new conjecture is needed.  The one-dimensional
marginal `<X,theta>` is log-concave and has variance one under isotropy.  The
sharp one-dimensional log-concave Cheeger estimate gives

\[
 {P_\mu(H)\over\min(\mu(H),1-\mu(H))}\ge c_1        \tag{4.1}
\]

for a universal `c_1>0`.  Hence `psi_mu>=c_1/(1+eta)`, and Cheeger's
inequality gives `K<=4(1+eta)^2/c_1^2`.  Thus an affine/exponential equality
branch is automatically bounded-scale.

This observation also prevents a product of isotropic exponential factors
from being a large-`K` near-Cheeger obstruction.  If an additional factor
had large Poincare constant, its own balanced cut would have perimeter of
order `K^{-1/2}` by Cheeger--Buser equivalence and would beat the constant
exponential-coordinate cut.

## 5. Integration with the fixed physical rank seed

The physical transfer theorem in `fixed_scale_physical_splicing.md` gives a
coarea boundary submeasure of mass at least `.004p` and effective rank at
least `17`.  There are now two rigorous extreme cases.

1. If a fixed amount of this submeasure lies on nonparallel facets joined by
   positive ridge capacity, the finite bevel supplies an actual perimeter
   saving, charged by Lemma 2.1 of that report.
2. If a level set splits into genuinely separated components, Lemma 1.1
   charges the split unless `q(v)` is nearly constant at their component
   masses.  On a one-dimensional affine branch, Theorem 2.1 and (4.1) put
   the measure in the bounded-`K` case.

In the exact zero-deficit case this closes completely for the smooth
full-support model.  Indeed zero coarea deficit at a level of mass
`v<=1/2` gives `P(A)=psi_mu v`; since `P(A)>=mathcal I(v)>=psi_mu v`, one
has `q(v)=psi_mu`.  Exact component additivity makes the profile linear on
`[0,v]`, while monotonicity and `q(1/2)=psi_mu` make it linear on
`[v,1/2]`.  Proposition 3.2 then bounds `K`.  The task is to make this
argument quantitative and stable under hard-support contact.

The unresolved intermediate statement is formalizable as follows.

> **Missing separation rigidity.**  If the physical normal submeasure has
> effective rank `17`, the total ridge/bevel capacity is below the amount
> needed for a `10^{-4}p` saving, and the coarea deficit is at most
> `6.02*10^{-5}p`, then a fixed amount of boundary flux is carried by
> separated components to which a quantitative version of Proposition 3.1
> applies, producing affine exponential halfspaces.

Neither concavity of `mathcal I` nor (1.3) proves this statement: they contain
no information identifying a component boundary with a flat marginal
quantile.  Conversely, no genuine high-`K` log-concave countermodel is known;
constructing one would amount to constructing the counterexample sequence
excluded by KLS.  The exact gap is the passage from scalar near-linearity to
spatial affine splitting in the presence of singular/contact sets.

## 6. Model audit

1. **Symmetric Laplace.**  `mathcal I(v)=min(v,1-v)`, so (2.7) has zero
   defect on both active branches.  The two components of a two-tail set have
   exactly additive perimeter and one common projective normal: the affine
   equality branch.
2. **One-sided exponential.**  The upper-tail branch has
   `mathcal I(v)=v`; Theorem 2.1 is exact.  Products retain bounded Poincare
   constant by tensorization.  The overlapping max-box geometry is not a
   separated-component equality; its actual corner incidence is strictly
   improved by the rounded competitor in the physical-splice report.
3. **Gaussian.**  `mathcal I(v)=varphi(Phi^{-1}(v))` is strictly concave, and
   `q` strictly decreases.  Every nontrivial separated split has positive
   scalar deficit.  Halfspaces form the affine branch; balls form the
   concurrent radial test but are not Gaussian isoperimetric minimizers.
4. **Cube.**  A one-dimensional uniform marginal has constant profile, not
   a linear profile, so splitting a level into two components pays a fixed
   extra boundary.  Coordinate half-cubes are affine, while inner boxes have
   bevelled ridges.
5. **Simplex.**  Barycentric marginals have beta-type, nonexponential
   profiles, so exact linear equality is absent.  Hyperplane cuts are affine
   and inner homothetic simplices have nonzero ridge capacity.  The known
   simplex spectral estimate independently places the isotropic simplex in
   the bounded-scale branch.

Thus every required model lands in the intended alternative.  The model
audit does not supply the missing general quantitative splitting theorem.
