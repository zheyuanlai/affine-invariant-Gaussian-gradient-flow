# The high-rank flat branch does not imply a global product

## 0. Verdict

The proposed local-to-global inverse is false without an additional global
coverage or factor-coupling hypothesis.  There are full-dimensional
log-concave measures with a nested foliation having all of the following
properties:

1. every regular leaf is a union of flat components with one common
   weighted constant mean curvature;
2. the lapse is constant on every component, the Jacobi curvature is zero,
   and every codimension-one support contact is natural and has zero contact
   curvature;
3. the normalized normal matrix is exactly `I/m`, and hence has effective
   rank `m`;
4. the normalized killed/cut flux up to any prescribed fixed distance is as
   small as desired; but
5. the support, and therefore the measure, is not affinely isomorphic to any
   nontrivial Cartesian product.

The example is a cyclically constrained exponential measure.  It directly
models locally flat components which terminate on different support faces.
The support faces responsible for global irreducibility are tangent to the
observed leaf components, and are not detected by their Jacobi or finite
tube charge.

Thus high normal rank is not the missing invariant.  The correct additional
datum is the **global factor-coupling normal fan**: support normals and
curvature directions outside the swept tube must be assigned to the same
factor decomposition as the persistent phase normals.  Sections 6 and 7
give an exact polyhedral criterion and a quantitative candidate invariant.

The construction below does not refute a theorem which also assumes a
fixed lower bound on swept flux, a balanced volume interval, and quantitative
coverage of all globally coupling faces.  It proves that those assumptions
cannot be inferred from small Jacobi/contact/killed charge.

## 1. The cyclic exponential measure

Fix an integer `m>=4` and `a>0`, with indices read modulo `m`.  Define

\[
 \Omega_{m,a}=\{x\in\mathbb R^m:
       x_i\ge0,\quad x_i+x_{i+1}\ge a\quad(1\le i\le m)\}.       \tag{1.1}
\]

This is a full-dimensional, pointed, unbounded convex polyhedron.  Put

\[
 d\mu_{m,a}(x)=Z_{m,a}^{-1}e^{-\sum_i x_i}
          \mathbf1_{\Omega_{m,a}}(x)\,dx.                       \tag{1.2}
\]

The normalizing constant is finite because `Omega_(m,a)` is contained in
the positive orthant.  Hence `mu_(m,a)` is a log-concave probability with
affine potential

\[
                         V(x)=\sum_i x_i+\log Z_{m,a}.           \tag{1.3}
\]

For `L>a`, let

\[
 A_L=\{x\in\Omega_{m,a}:\max_i x_i\ge L\}.                     \tag{1.4}
\]

The sets `A_L` form a decreasing nested family.  Up to a set of
`H^(m-1)`-measure zero, their relative boundary is the disjoint union

\[
 \Gamma_i(L)=\{x\in\Omega_{m,a}:x_i=L, x_j<L\ (j\ne i)\}.
                                                                    \tag{1.5}
\]

The ties `x_i=x_j=L` are precisely the cut-locus ridges.  They have
codimension two and no perimeter mass.

## 2. Exact zero Jacobi and support-contact charge

Use Euclidean anisotropy `Phi(xi)=|xi|`.  On `Gamma_i(L)`, the unit normal
pointing out of `A_L` and into its complement is

\[
                              N_i=-e_i.                         \tag{2.1}
\]

The component is flat, so its shape operator is zero.  Since
`nabla V=(1,...,1)`, its weighted mean curvature, with the convention in
`anisotropic_foliation_inverse.md`, is

\[
 H_\mu=\operatorname {tr}S-\nabla V\cdot N_i=1.                 \tag{2.2}
\]

Thus every component of every leaf has the same constant weighted mean
curvature.  The potential is affine, so `nabla^2 V=0`.  Consequently the
Jacobi density is identically zero on every regular component:

\[
 |S|_{HS}^2+\nabla^2V(N_i,N_i)=0.                               \tag{2.3}
\]

The geometric speed when `L` changes is constant on all components.  More
explicitly, if `v(L)=mu_(m,a)(A_L)`, coarea gives `v'(L)=-P(L)`.  When the
family is oriented in the direction of increasing volume, its normal lapse
is therefore

\[
                              f={1\over P(L)}                       \tag{2.4a}
\]

on every component, and `int_(partial A_L) f dsigma_mu=1`.  Thus its
tangential logarithmic gradient is zero.  The CMC slope (2.2) is constant,
so its derivative is also zero.

It remains to check the hard wall, rather than silently discard it.  On
`Gamma_i(L)`, the two support faces

\[
 x_i+x_{i-1}=a,\qquad x_i+x_{i+1}=a
\]

are strictly inactive because `x_i=L>a` and all coordinates are
nonnegative.  Every active codimension-one support face is therefore one of

\[
 x_j=0\quad(j\ne i),
 \qquad x_j+x_{j+1}=a\quad(i\notin\{j,j+1\}).                   \tag{2.4}
\]

Its normal is orthogonal to `e_i`.  Hence `Gamma_i(L)` meets every active
support face in the natural orthogonal free-boundary angle.  All support
faces are flat, so the contact-curvature term is exactly zero.  Intersections
of several such faces cause no hidden loss: translation in the direction
`-e_i` is tangent to every face in (2.4), and, as long as the translated
`i`-th coordinate remains above `a`, it cannot reach either inactive face.

We have therefore proved the following exact statement.

**Lemma 2.1 (zero regular charge).**  On any band `L in [L_0,L_1]` with
`L_0>a`, the nested family (1.4) is a fixed-topology foliation away from its
codimension-two tie ridges.  On every regular spacetime component the lapse
gradient, anisotropic second fundamental form, potential-Hessian charge,
CMC slope drop, and support-contact curvature all vanish.  All
codimension-one support contacts satisfy the natural boundary condition.

## 3. High rank and an exact killed-flux estimate

Let

\[
 \mathcal P_k(a)=\{y\in[0,\infty)^k:
           y_j+y_{j+1}\ge a\ (1\le j<k)\}                     \tag{3.1}
\]

and let `nu_k^path` be the probability on this path polyhedron with density
proportional to `exp(-sum_j y_j)`.  Removing coordinate `i` from a cyclic
constraint leaves exactly `P_(m-1)(a)`.  Consequently, if

\[
 \vartheta_k(s)=\nu_k^{\rm path}\{\max_jY_j<s\},                \tag{3.2}
\]

then the weighted area of every component in (1.5) is

\[
 p_i(L)=Z_{m,a}^{-1}e^{-L}Z_{m-1,a}^{\rm path}
                            \vartheta_{m-1}(L).                 \tag{3.3}
\]

Cyclic symmetry makes these areas equal.  Thus, for the normal matrix and
total relative perimeter,

\[
 M(L)=\sum_i p_i(L)e_ie_i^T={P(L)\over m}I_m,
 \qquad P(L)=\sum_i p_i(L),                                   \tag{3.4}
\]

and in particular

\[
 Q(L):={M(L)\over\operatorname {tr}M(L)}={I_m\over m},
 \qquad {\operatorname {tr}Q(L)\over\|Q(L)\|_{op}}=m.         \tag{3.5}
\]

We now retain every support and cut event.  A point `x in Gamma_i(L)` sent
distance `t` along its outward normal becomes `x-te_i`.  If `0<=t<L-a`,
this point remains in `Omega_(m,a)`.  Because the support is upward closed,
the distance from a point `y` in the complement of `A_L` to `A_L` is

\[
                         L-\max_jy_j.                          \tag{3.6}
\]

The ray based on component `i` is therefore still the unique minimizing
normal ray at time `t` exactly when

\[
                         x_j<L-t\quad(j\ne i).                  \tag{3.7}
\]

It follows that the surviving flux `R_L(t)` satisfies the exact identity

\[
 {R_L(t)\over P(L)}={\vartheta_{m-1}(L-t)
                            \over\vartheta_{m-1}(L)}.           \tag{3.8}
\]

There is a dimension-free elementary tail bound for the path law.  Given
all coordinates except `Y_j`, its conditional distribution is a unit-rate
exponential translated by

\[
 \ell_j=\max\{0,a-Y_{j-1},a-Y_{j+1}\}\in[0,a],                 \tag{3.9}
\]

with the absent neighbor omitted at an endpoint.  Hence, for `s>=a`,

\[
 \nu_k^{\rm path}\{Y_j\ge s\mid(Y_r)_{r\ne j}\}
       =e^{-(s-\ell_j)}\le e^{a-s}.                            \tag{3.10}
\]

The union bound gives

\[
             1-\vartheta_k(s)\le k e^{a-s}.                    \tag{3.11}
\]

Combining (3.8) and (3.11), and writing
`delta=(m-1)e^(a-L)<1`, yields

\[
 \boxed{
  0\le1-{R_L(t)\over P(L)}
  \le {\delta e^t\over1-\delta}}
  \qquad(0\le t<L-a).                                        \tag{3.12}
\]

In particular, fix `T>0` and `0<epsilon<1/2`, and choose

\[
 L_0=a+\log{m-1\over\epsilon},
 \qquad \log{m-1\over\epsilon}>T.                            \tag{3.13}
\]

Then every leaf with `L>=L_0` loses at most

\[
                  {\epsilon e^T\over1-\epsilon}               \tag{3.14}
\]

of its base flux before distance `T`.  This tends to zero with `epsilon`,
while (3.5) retains effective rank `m` exactly.

The packet is not hiding an anomalous perimeter normalization.  Adding the
missing endpoint constraint to the path law produces the cyclic partition
function, and

\[
 e^{-a}\le
 \nu_{m-1}^{\rm path}\{Y_1+Y_{m-1}\ge a\}\le1.                \tag{3.15}
\]

The lower bound follows from the subevent `Y_1>=a` and (3.9).  Hence

\[
 1\le{Z_{m-1,a}^{\rm path}\over Z_{m,a}}\le e^a.               \tag{3.16}
\]

At `L=L_0`, (3.3), (3.11), and (3.16) imply

\[
 {m\over m-1}e^{-a}\epsilon(1-\epsilon)
       \le P(L_0)\le {m\over m-1}\epsilon.                    \tag{3.17}
\]

Thus the unnormalized retained trace is comparable to `epsilon`, with
constants independent of `m`.

There is also a near-log-affine tail profile.  Coarea gives
`mu(A_L)=int_L^infty P(s)ds`.  From (3.3) and (3.11),

\[
 1-\delta\le {P(L)\over\mu(A_L)}
       \le {1\over1-\delta/2}.                                \tag{3.18}
\]

Thus the counterexample retains not only local CMC, but also an
arbitrarily accurate one-sided exponential perimeter/volume law on the
chosen band.

## 4. The support is globally product-irreducible

Every inequality in (1.1) is irredundant.  The irredundant facet normals are

\[
                 e_i,qquad e_i+e_{i+1}\quad(1\le i\le m).     \tag{4.1}
\]

Suppose, for contradiction, that an affine image of `Omega_(m,a)` were a
nontrivial Cartesian product.  In the dual vector space there would then be
a nontrivial direct sum `U+W` such that every facet normal belongs either to
`U` or to `W`: facets of a product are a facet of one factor times the whole
other factor.

Each `e_i` must belong to one of the two summands.  If `e_i` and `e_(i+1)`
belonged to different summands, their sum `e_i+e_(i+1)` would belong to
neither summand, contradicting (4.1).  Hence every adjacent pair of
coordinate normals belongs to the same summand.  The cycle is connected,
so all `e_i` belong to one summand.  They span the full dual space, forcing
the other summand to be zero.  This contradiction proves:

**Lemma 4.1 (affine irreducibility).**  For every `m>=4` and `a>0`, the
polyhedron `Omega_(m,a)` is not affinely isomorphic to a nontrivial
Cartesian product.  In particular, `mu_(m,a)` is not an affine product of
one-dimensional exponential, slab, or any other factors.

This is an incidence obstruction, not an orthogonality argument.  It is
unchanged by arbitrary invertible affine maps.

Combining Lemmas 2.1 and 4.1 with (3.5) and (3.12) refutes the proposed
inverse: local zero Jacobi/contact charge and vanishing normalized killed
flux do not force a global affine product.  The mixed faces
`x_i+x_(i+1)=a` connect all would-be factors, but on `Gamma_i(L)` the two
mixed faces involving `i` are inactive and every other mixed face is tangent
to the normal trajectory.  A finite tube based at the tail leaf cannot see
the global coupling.

## 5. Isotropic and uniformly elliptic version

The obstruction is not an artifact of ill-conditioned affine coordinates.
Choose

\[
                              a_m=m^{-6}.                        \tag{5.1}
\]

Let `pi_m` be the product of `m` unit exponentials and let
`E={X in Omega_(m,a_m)}`.  Since

\[
 \mathbb P\{X_i+X_{i+1}<a_m\}
     =1-e^{-a_m}(1+a_m)\le {a_m^2\over2},                       \tag{5.2}
\]

the union bound gives

\[
 q_m:=\pi_m(E^c)\le {1\over2m^{11}}.                            \tag{5.3}
\]

For a unit vector `theta`, put
`Y=sum_i theta_i(X_i-1)`.  Independence gives

\[
                 \mathbb E_{\pi_m}Y^2=1,
 \qquad \mathbb E_{\pi_m}Y^4
       =3+6\sum_i\theta_i^4\le9.                               \tag{5.4}
\]

Cauchy--Schwarz and conditioning on `E` therefore give

\[
 \operatorname {Var}(Y\mid E)
 \ge {1-3\sqrt {q_m}\over1-q_m}-{q_m\over(1-q_m)^2},
 \qquad
 \operatorname {Var}(Y\mid E)\le {1\over1-q_m}.               \tag{5.5}
\]

For `m>=4`, the covariance `A_m` of `mu_(m,a_m)` consequently satisfies,
with ample slack,

\[
                            .99 I\preceq A_m\preceq1.01 I.      \tag{5.6}
\]

Center and whiten by

\[
                 y=A_m^{-1/2}(x-\mathbb E x).                   \tag{5.7}
\]

The resulting measure is isotropic.  Euclidean weighted perimeter in the
`x` variables becomes the constant ellipsoidal anisotropy

\[
                         \Phi_m(n)=|A_m^{-1/2}n|                 \tag{5.8}
\]

in the `y` variables.  Its ellipticity constants are bounded universally by
(5.6).  Affine covariance of anisotropic first and second variation
preserves all zero-charge statements in Lemma 2.1, the natural Young
contact condition, and the killed-flux ratio (3.12).  Indeed the Wulff ray
in `y` is exactly the inverse affine image of the Euclidean ray in `x`.

Cyclic symmetry makes the diagonal entries of `A_m` equal.  The transformed
unit normals are

\[
 n_i={A_m^{1/2}e_i\over|A_m^{1/2}e_i|}.
\]

Their equal-weight normalized matrix is

\[
 {1\over m}\sum_i n_in_i^T={A_m\over\operatorname {tr}A_m},    \tag{5.9}
\]

whose effective rank is at least

\[
 {\operatorname {tr}A_m\over\|A_m\|_{op}}
           \ge {.99\over1.01}m.                                \tag{5.10}
\]

The support remains affinely irreducible by Lemma 4.1.  Hence the
counterexample exists in isotropic position, with a universally elliptic
constant anisotropy and unbounded effective normal rank.

## 6. What exact global product rigidity would require

The obstruction suggests an exact replacement for informal phrases such as
"many flat phases are a product."  For a full-dimensional polyhedron

\[
                   \Omega=\bigcap_{r=1}^N\{a_r(x)\ge b_r\},     \tag{6.1}
\]

with irredundant facets, define its **facet-normal matroid** to be the linear
matroid of `(a_r)`.  The following elementary criterion does not assume
orthogonality.

**Proposition 6.1 (polyhedral product criterion).**  The polyhedron `Omega`
is affinely a product \(\Omega_1\times\cdots\times\Omega_k\) if and only if there
is a dual direct sum

\[
                       (\mathbb R^m)^*=U_1\oplus\cdots\oplus U_k \tag{6.2}
\]

such that every irredundant facet normal `a_r` lies in one `U_j`.  With this
decomposition, an affine potential `c(x)` splits automatically into the sum
of its restrictions to the primal dual factors.  Hence the corresponding
log-affine probability is a product.  One-dimensional factorization occurs
exactly when all `U_j` are one-dimensional.

**Proof.**  For a product, every facet is a facet of one factor times all
other factors, proving necessity.  Conversely, write
`x=x_1+...+x_k` in the primal decomposition dual to (6.2).  Every inequality
in (6.1) then involves exactly one `x_j`.  Grouping the inequalities by `j`
turns (6.1) into the Cartesian product of their solution sets.  A linear
functional is the sum of its restrictions to the direct summands, so its
exponential density factors as well.  QED.

Thus a global product theorem must control **all** irredundant support
normals, not only normals of the selected leaf.  For a nonsmooth convex
potential, directions in the matrix-valued Hessian measure play the same
role: mixed curvature can couple factors even when the support is a product.

## 7. Replacement invariant and the remaining viable theorem

After whitening, let `U_1+...+U_k` be a proposed orthogonal dual factor
decomposition and let `Pi_j` be its projections.  For a polyhedral support,
a quantitative support-coupling defect is

\[
 \mathfrak C_\Omega(U_1,\ldots,U_k)
 =\sum_r w_r\left(1-\max_j{|\Pi_j a_r|^2\over|a_r|^2}\right), \tag{7.1}
\]

where `w_r` is a specified global weighted surface or normal-trace weight.
For a general convex potential one must add the mixed Hessian defect

\[
 \mathfrak C_V(U_1,\ldots,U_k)
 =\sum_{i\ne j}\int
       \|\Pi_i\,d(D^2V)\,\Pi_j\|_* .                           \tag{7.2}
\]

The precise reference measure in (7.1)--(7.2) must be chosen by the global
argument; local leaf area is insufficient.  Vanishing of (7.1) is the
normal-fan condition in Proposition 6.1.  In the cyclic example, every
normal `e_i+e_(i+1)` couples two proposed one-dimensional phase factors, and
the resulting coupling graph is the full cycle.  Yet none of these edges is
charged by the finite tail tube.

A viable high-rank inverse must therefore have the following form.

* Either the swept Wulff tubes carry a fixed amount of flux to every mixed
  support/Hessian coupling, in which case that coupling is charged to
  contact, curvature, cut, or a perimeter competitor;
* or the unvisited coupling mass is itself included as a **coverage defect**;
* or the global coupling graph splits into blocks of universally bounded
  dimension.  Only in the last case may one apply finite-dimensional bounds
  to the blocks and then tensorize.

High rank of the retained normal matrix supplies none of these statements.
In particular, replacing (7.1) by pairwise orthogonality of the observed
normals is invalid: affine products need not be orthogonal before whitening,
and the cyclic example has coordinate phase normals which are already
orthogonal while its support remains irreducible.

## 8. Consequence for the joint-replacement route

The high-rank item 4 in Section 8 of
`anisotropic_foliation_inverse.md` cannot be closed by a theorem using only
the listed local Jacobi/contact/killed quantities.  Equations (1.1)--(5.10)
give a family for which those quantities tend to zero after normalization,
the normal rank is unbounded, and the global product conclusion is false.

The route must additionally extract, from balanced heat levels and joint
equimeasurable minimality, a fixed-flux **global saturation** statement.
Such a statement would rule out the tail construction because (3.17) shows
that its swept flux vanishes with the killed-flux parameter.  It would then
still have to bound the global coupling defects (7.1)--(7.2), or turn them
into a physical perimeter saving.  That is a strictly stronger and more
accurate remaining lemma than a rank-to-product inverse.
