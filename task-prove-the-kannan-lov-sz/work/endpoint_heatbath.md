# Endpoint incidence, heat-bath forms, and the failure of abstract gap rigidity

## Executive result

There is an exact ideal mechanism.  Weight the ray quotient by its boundary
density and let `B,C` be the positive and negative endpoint maps.  If the
junction form is the sum of the two conditional variances

\[
 \mathcal D_{HB}(h)
 =E\operatorname {Var}(h\mid B)+E\operatorname {Var}(h\mid C), \tag{0.1}
\]

then its spectral gap, normalized so that the largest possible gap is one,
is

\[
 \lambda_{HB}=1-\rho(B,C),                                \tag{0.2}
\]

where `rho(B,C)` is maximal correlation.  Exact gap one means that the two
endpoint maps are independent.  If every calibrated endpoint pair has the
same length `d`, independence upgrades same-ray equality to
`|b-c|=d` for the full product of the endpoint laws, and the endpoint sets
are spheres in orthogonal affine subspaces (the Clifford branch).

The geometric junction form is not generally (0.1).  On an endpoint fiber
with ray weights `pi_i` and interface conductances `w_ij`, it is a heat bath
if and only if

\[
 w_{ij}=c\,\pi_i\pi_j\quad(i\ne j).                     \tag{0.3}
\]

Generic medial coefficients are instead `rho/|N_i-N_j|`, with no reason for
(0.3).  Every such graph form has a reversible-kernel representation after
adding self-loops, but its time normalization is arbitrary.  The exact total
translation charge controls

\[
 {1\over2}\int |N-N'|^2\,dW(N,N'),                      \tag{0.4}
\]

not the scalar conductance `W(Omega^2)` that normalizes a Gibbs kernel.  A
two-sheet focal junction with normal gap `theta` has scalar coefficient
of order `1/theta` but translation charge of order `theta`.  Hence stability
plus total charge does not force (0.2) to be close to one.

Even granting the ideal heat bath, near-maximal gap alone does not give a
dimension-free approximate Clifford theorem.  A projective-plane incidence
model below has

\[
 \lambda_{HB}=1-\frac{\sqrt q}{q+1}\longrightarrow1,     \tag{0.5}
\]

exactly balanced long Euclidean ray data, mean-zero dispersed normals, and
constant same-ray cross distance.  Nevertheless the two endpoint marginals
are both isotropic on the *same* `(q^2+q)`-dimensional subspace and

\[
 \operatorname {Var}_{\mu_B\otimes\mu_C}|B-C|^2
 =4(q^2+q).                                               \tag{0.6}
\]

Thus it is neither aligned nor approximately Clifford in the squared
cross-distance/ANOVA sense.  After tensorization it also has effective normal
rank of order the square of the ray scale and exponentially small mass in
every coherent unit direction cap.  This strengthened model is still atomic,
and a separate radial-variance argument shows that its equal-radius long
segments cannot carry constant mass in any isotropic log-concave realization.
Thus it sharply refutes the abstract endpoint-incidence implication while
also identifying where global log-concavity genuinely rules out this
particular countergeometry.

## 1. Natural ray measure and endpoint maps

Let `eta` be the ray quotient, let

\[
 b_y=q_y(0),\qquad P=\int b_y\,d\eta(y),                 \tag{1.1}
\]

and define the boundary-weighted ray probability

\[
 d\pi(y)={b_y\over P}\,d\eta(y).                        \tag{1.2}
\]

This is the natural stationary measure for a junction spectral problem.  If
`bar h_eta=int h d eta`, then the left side of signed-distance stability is

\[
 \int b_y(h-\bar h_\eta)^2d\eta
 =P\left(\operatorname {Var}_\pi h
          +(E_\pi h-E_\eta h)^2\right)
 \geq P\operatorname {Var}_\pi h.                       \tag{1.3}
\]

Restrict for the moment to rays with two finite completed endpoints and let

\[
 B=e_+(Y),\qquad C=e_-(Y),\qquad Y\sim\pi.              \tag{1.4}
\]

The joint endpoint law is `gamma=(e_+,e_-)_#pi`.  Fibers of `e_+` are the
rays incident to one positive medial/focal endpoint; fibers of `e_-` are the
negative incidences.  The conditional expectation projections on
`L^2(pi)` are

\[
 P_+h=E[h\mid B],\qquad P_-h=E[h\mid C].                \tag{1.5}
\]

Their heat-bath forms are

\[
 \mathcal D_+(h)=\|h-P_+h\|_2^2=E\operatorname {Var}(h\mid B),
 \quad
 \mathcal D_-(h)=\|h-P_-h\|_2^2=E\operatorname {Var}(h\mid C). \tag{1.6}
\]

Rays with an infinite endpoint have no resampling block on that side; in a
Markov description they contribute a self-loop.  Ignoring them can only
overstate the two-block gap.

## 2. When is the medial form a conditional variance?

At a finite endpoint `x`, suppose the incident ray fiber is
`F_x={1,...,m}`, its conditional ray law is `pi^x=(pi_1,...,pi_m)`, and the
geometric junction form is

\[
 \mathcal J_x(h)=\sum_{i<j}w_{ij}(h_i-h_j)^2,
 \qquad w_{ij}=w_{ji}\geq0.                             \tag{2.1}
\]

For a generic medial interface, `w_ij` is the pushforward of

\[
 {\rho(z)\over|N_i(z)-N_j(z)|}\,d\mathcal H^{n-1}(z),   \tag{2.2}
\]

including the global sign weight `q` or `p` on the corresponding side.

**Lemma 2.1 (exact heat-bath factorization criterion).**  For a scalar
`c_x`, the identity

\[
 \mathcal J_x(h)=c_x\operatorname {Var}_{\pi^x}h
 \quad\hbox{for every }h                                \tag{2.3}
\]

holds if and only if

\[
 w_{ij}=c_x\pi_i\pi_j\qquad(i\ne j).                   \tag{2.4}
\]

For a two-ray fiber this always holds, but with the variable coefficient

\[
 c_x={w_{12}\over\pi_1\pi_2}.                           \tag{2.5}
\]

**Proof.**  Conditional variance has the complete-graph representation

\[
 \operatorname {Var}_{\pi^x}h
 ={1\over2}\sum_{i,j}\pi_i\pi_j(h_i-h_j)^2.             \tag{2.6}
\]

Comparing the off-diagonal entries of the two Laplacian matrices proves
(2.4).  Conversely, substitution proves (2.3).  QED.

Thus a multiway endpoint must have rank-one complete-graph conductances in
the ray weights.  Pairwise epi-derivatives such as (2.2) do not imply this.
In measure form, the identical criterion is that the off-diagonal
disintegration of the conductance over an endpoint `x` be
`c(x) pi_x tensor pi_x`.

There is always an exact weighted Markov analogue, but it is weaker and has
no intrinsic time normalization.

**Lemma 2.2 (Markovization with self-loops).**  Put
`d_i=sum_{j ne i}w_ij`.  For any

\[
 c_x\geq\max_i{d_i\over\pi_i},                           \tag{2.7}
\]

define

\[
 K_x(i,j)={w_{ij}\over c_x\pi_i}\quad(i\ne j),
 \qquad
 K_x(i,i)=1-{d_i\over c_x\pi_i}.                        \tag{2.8}
\]

Then `K_x` is a Markov kernel reversible for `pi^x` and

\[
 \mathcal J_x(h)
 =c_x\langle h,(I-K_x)h\rangle_{L^2(\pi^x)}.            \tag{2.9}
\]

**Proof.**  Nonnegativity and unit row sums follow from (2.7).  Detailed
balance is `pi_iK(i,j)=w_ij/c_x=pi_jK(j,i)`.  The standard reversible
Dirichlet identity gives

\[
 c_x\langle h,(I-K_x)h\rangle
 ={1\over2}\sum_{i,j}w_{ij}(h_i-h_j)^2=\mathcal J_x(h).
\]

QED.

Replacing `c_x` by `A c_x` replaces `K_x` by
`I-(I-K_x)/A` and leaves (2.9) unchanged, while dividing the Markov gap by
`A`.  Hence “the endpoint Markov gap” has no invariant meaning until the
coefficient is canonically fixed.  The heat bath is the special
normalization `K_x(i,j)=pi_j`, including the harmless self-loop `j=i`.

## 3. Exact spectrum of the ideal two-block heat bath

Let

\[
 \rho(B,C)=\sup
 \left\{E[f(B)g(C)]: Ef=Eg=0, Ef^2=Eg^2=1\right\}      \tag{3.1}
\]

be maximal correlation.  Assume the bipartite incidence relation is
connected, so the only functions measurable with respect to both endpoint
maps are constants.

**Theorem 3.1 (two-projection gap).**  The form

\[
 \mathcal D_{HB}=\mathcal D_++\mathcal D_-             \tag{3.2}
\]

has spectral gap

\[
 \boxed{\quad
 \inf_{Eh=0}{\mathcal D_{HB}(h)\over E h^2}
 =1-\rho(B,C).
 \quad}                                                  \tag{3.3}
\]

Its largest possible gap is one.  It equals one if and only if `B,C` are
independent.

**Proof.**  On `L^2(pi)`,

\[
 \mathcal D_{HB}(h)=\langle h,(2I-P_+-P_-)h\rangle.     \tag{3.4}
\]

The nonzero principal-angle cosines between the centered ranges of `P_+`
and `P_-` are the singular values of
`g(C) mapsto E[g(C)|B]`; the largest is (3.1).  On every corresponding
two-dimensional principal plane, `P_++P_-` has eigenvalues `1+rho_j` and
`1-rho_j`.  It is zero on the orthogonal complement of the two ranges.
After removing the constant eigenvector, its largest eigenvalue is
`1+rho(B,C)`.  Equation (3.3) follows from (3.4).

If `rho=0`, every centered product correlation vanishes, hence
`E[f(B)g(C)]=Ef Eg` for all bounded `f,g`; rectangles generate the product
sigma-algebra, so the endpoint variables are independent.  The converse is
immediate.  QED.

The random-scan Gibbs kernel `(P_++P_-)/2` has gap
`(1-rho)/2`; (3.2) is the normalization with maximum one.

Unequal sign weights change the normalization.  For `a,b>0`, the gap of

\[
 a\mathcal D_++b\mathcal D_-                            \tag{3.5}
\]

is

\[
 \lambda_{a,b}
 ={a+b-\sqrt{(a-b)^2+4ab\rho(B,C)^2}\over2}.            \tag{3.6}
\]

The maximal value is `min(a,b)`, attained at independence.  Formula (3.6)
follows from the same two-dimensional principal-plane calculation.  Thus
the factors `q,p` in the positive and negative medial charges must either be
retained or divided out; treating their sum as an unweighted gap changes the
claimed constant.

## 4. What stability would imply under the missing coefficient hypotheses

Write signed-distance stability, after (1.3), schematically as

\[
 P\operatorname {Var}_\pi h
 \leq\mathcal S(h)+\mathcal J(h),                         \tag{4.1}
\]

where `S` is the smooth normal-chart form and `J` is the completed
medial/focal form.

**Proposition 4.1 (conditional gap extraction).**  Suppose, for every
centered `h`,

\[
 \mathcal S(h)\leq\varepsilon_s P\operatorname {Var}_\pi h,
 \qquad
 \mathcal J(h)\leq(1+\varepsilon_c)P\mathcal D_{HB}(h). \tag{4.2}
\]

Then

\[
 \lambda_{HB}\geq{1-\varepsilon_s\over1+\varepsilon_c}
 \geq1-(\varepsilon_s+\varepsilon_c).                  \tag{4.3}
\]

**Proof.**  Substitute (4.2) in (4.1), cancel the smooth term, and take the
infimum over centered `h`.  QED.

This is the desired heat-bath mechanism, but neither inequality in (4.2) is
currently available.  The long-ray curvature estimate controls the smooth
form on special normal-coordinate tests; it is not an operator bound for all
endpoint heights.  Lemma 2.1 shows that the second comparison requires
fiberwise factorization and a coefficient bound.

The exact translation trace does not provide that coefficient bound.  With
`N_y` the ray normal, put

\[
 \mathcal S_N=\sum_k\mathcal S(\langle e_k,N\rangle),
 \qquad
 \mathcal J_N=\sum_k\mathcal J(\langle e_k,N\rangle).
\]

The full smooth-plus-focal identity is

\[
 \mathcal S_N+\mathcal J_N=P-{K\over2}\leq P,           \tag{4.4}
\]

where `K>=0` is the convexity charge.  If `W` is the symmetric junction
conductance measure, then

\[
 \mathcal J_N={1\over2}\int|N_y-N_{y'}|^2dW(y,y').     \tag{4.5}
\]

It does not control `W(Omega^2)`.  In particular, it pins only the product
of a Markov time coefficient and the Dirichlet energy of `N`; it cannot
separate that coefficient from the Markov gap.

At the abstract form level this already gives a counterexample.  Let an
ideal endpoint heat bath have any gap `lambda in (0,1)`, set

\[
 \mathcal J={P\over\lambda}\mathcal D_{HB},\qquad
 \mathcal S=0,                                             \tag{4.5a}
\]

and take a unit-variance bottom eigenfunction `h_0`.  Then (4.1) holds for
every height and `J(h_0)=P`, exactly saturating a prescribed trace charge,
while the unscaled Gibbs gap is the arbitrary number `lambda`.  Thus
stability plus one exact trace controls `(coefficient) times (gap)`, not the
gap.  The next example shows that the same coefficient pathology is
geometrically natural for unit normal jumps.

**Example 4.2 (two-sheet coefficient blow-up).**  Consider a two-chart
medial interface of weighted area `A` and let

\[
 \theta=|N_1-N_2|.
\]

The geometric form is

\[
 \mathcal J(h)={A\over\theta}(h_1-h_2)^2.               \tag{4.6}
\]

For conditional ray weights `(1/2,1/2)`,

\[
 \mathcal J(h)={4A\over\theta}\operatorname {Var}h,    \tag{4.7}
\]

while its translation charge is

\[
 \mathcal J_N=A\theta.                                  \tag{4.8}
\]

As `theta` tends to zero, the heat-bath coefficient diverges and the exact
total charge tends to zero.  This is precisely the focal regime in which
the denominator in the envelope epi-derivative vanishes.  It disproves any
attempt to read a heat-bath normalization from (4.4).

There is one valid alignment alternative.  If

\[
 \tau=E_\pi|N-E_\pi N|^2=1-|E_\pi N|^2,                \tag{4.9}
\]

then, when `tau<1`, for `v=E_pi N/|E_pi N|`,

\[
 E_\pi|N-v|^2=2(1-|E_\pi N|)\leq2\tau.                \tag{4.10}
\]

Thus small normal trace variance is exactly the aligned branch.  When
`tau` is bounded below, (4.4) merely says that stability nearly exhausts one
trace budget; without (4.2), it still says nothing about the bottom of the
full endpoint spectrum.

## 5. Exact and quantitative Clifford rigidity in the ideal model

**Theorem 5.1 (exact heat-bath/Clifford rigidity).**  Suppose
`lambda_HB=1` and the joint endpoint law is supported on

\[
 |B-C|=d                                                     \tag{5.1}
\]

for one constant `d`.  Then `gamma=mu_B tensor mu_C`, and

\[
 |b-c|=d\quad\text{for }\mu_B\otimes\mu_C\text{-a.e. }(b,c).
                                                               \tag{5.2}
\]

Consequently, after passing to essential supports, the two endpoint sets lie
on spheres in mutually orthogonal affine subspaces.

**Proof.**  Theorem 3.1 gives independence, hence (5.2).  Subtracting the
four squared identities for `(b,c),(b',c),(b,c'),(b',c')` gives

\[
 \langle b-b',c-c'\rangle=0.                            \tag{5.3}
\]

The two difference spans are orthogonal.  Choosing closest affine centers
and fixing one endpoint at a time in (5.2) shows that each endpoint law has
constant radius in its own span.  QED.

Near-maximal gap controls only an operator norm.  A quantitative product
conclusion requires an effective-rank hypothesis.

**Theorem 5.2 (quantitative rigidity with correlation rank).**  Suppose
`lambda_HB>=1-epsilon`, the joint law is absolutely continuous with respect
to `mu_B tensor mu_C`, and the centered conditional-expectation operator

\[
 T_0:g(C)\longmapsto E[g(C)\mid B]-Eg(C)                 \tag{5.4}
\]

has rank at most `r`.  Then

\[
 \chi^2(\gamma\|\mu_B\otimes\mu_C)\leq r\varepsilon^2,
 \qquad
 \|\gamma-\mu_B\otimes\mu_C\|_{TV}
 \leq {\sqrt r\,\varepsilon\over2}.                   \tag{5.5}
\]

If additionally (5.1) holds and

\[
 H=\operatorname*{ess\,sup}_{\mu_B\otimes\mu_C}
       \left||B-C|^2-d^2\right|<\infty,                 \tag{5.6}
\]

then

\[
 \operatorname {Var}_{\mu_B\otimes\mu_C}|B-C|^2
 \leq {H^2\sqrt r\,\varepsilon\over2}.                 \tag{5.7}
\]

Writing `b=EB`, `c=EC`, `U=B-b`, `V=C-c`, and `a=b-c`, the exact ANOVA
identity therefore gives

\[
\begin{split}
 &\operatorname {Var}(|U|^2+2\langle a,U\rangle)
 +\operatorname {Var}(|V|^2-2\langle a,V\rangle)\\
 &\hspace{35mm}+4\operatorname {Tr}(\Sigma_B\Sigma_C)
 \leq {H^2\sqrt r\,\varepsilon\over2}.                 \tag{5.8}
\end{split}
\]

Thus the endpoint laws are quantitatively radial and their covariance spans
are quantitatively orthogonal whenever `sqrt(r) epsilon` is small.

**Proof.**  By Theorem 3.1, `||T_0||_op<=epsilon`.  If
`k=d gamma/d(mu_B tensor mu_C)-1`, then `T_0` is the integral operator with
kernel `k`; hence

\[
 \chi^2(\gamma\|\mu_B\otimes\mu_C)
 =\|T_0\|_{HS}^2\leq r\|T_0\|_{op}^2.                  \tag{5.9}
\]

Cauchy--Schwarz gives the total-variation bound.  Since `gamma` gives full
mass to (5.1), the product law gives its complement mass at most the total
variation distance.  On that complement the squared defect in (5.6) is at
most `H^2`, proving (5.7).  Identity (5.8) is the orthogonal ANOVA
decomposition of the product cross-distance variance.  QED.

The same proof allows an approximate same-ray relation.  If
`G=|B-C|^2-d^2`, `E_gamma G^2<=zeta^2`, and `|G|<=H` under the product law,
then the right sides of (5.7)--(5.8) may be replaced by

\[
 \zeta^2+{H^2\sqrt r\,\varepsilon\over2}.               \tag{5.10}
\]

**Corollary 5.3 (conditional alignment/Clifford dichotomy).**  Under the
hypotheses of Proposition 4.1, put

\[
 \varepsilon_g
 =1-{1-\varepsilon_s\over1+\varepsilon_c}
 ={\varepsilon_s+\varepsilon_c\over1+\varepsilon_c}.    \tag{5.11}
\]

Assume also the rank, bounded-defect, and constant-length hypotheses of
Theorem 5.2.  For every `tau_0 in (0,1)`, one of the following holds:

\[
 \begin{array}{ll}
 \text{alignment:}&E_\pi|N-v|^2\leq2\tau_0
       \text{ for some unit }v,\\[1mm]
 \text{Clifford defect:}&
 \operatorname {Var}(|U|^2+2\langle a,U\rangle)
 +\operatorname {Var}(|V|^2-2\langle a,V\rangle)
 +4\operatorname {Tr}(\Sigma_B\Sigma_C)
 \leq {H^2\sqrt r\,\varepsilon_g\over2}.
 \end{array}                                             \tag{5.12}
\]

Indeed, use the first line when the normal variance `tau` in (4.9) is at
most `tau_0`; otherwise Proposition 4.1 and Theorem 5.2 give the second.
The second estimate in fact holds independently of `tau`; the split merely
records the harmless aligned branch.

The smooth operator bound, heat-bath coefficient comparison, bounded
endpoint length, and effective-rank hypothesis used in Corollary 5.3 are all
additional; none follows from total charge alone.

## 6. A Euclidean near-maximal-gap counterexample

The rank loss in Theorem 5.2 is necessary even for Euclidean endpoint data.

Let `q` be a prime power and take a projective plane of order `q`.  It has

\[
 M=q^2+q+1                                                   \tag{6.1}
\]

points and the same number of lines.  Every point lies on
`d=q+1` lines, and two distinct points lie on exactly one common line.  If
`A` is the point-line incidence matrix, then

\[
 AA^T=qI+J.                                                \tag{6.2}
\]

Let `gamma` be uniform on incidences.  Its endpoint marginals are uniform.
The nonconstant singular values of the normalized conditional-expectation
operator `A/d` are all `sqrt(q)/(q+1)`.  Therefore

\[
 \boxed{\quad
 \lambda_{HB}=1-{\sqrt q\over q+1}\longrightarrow1.
 \quad}                                                    \tag{6.3}
\]

Nevertheless

\[
 \|\gamma-\mu_B\otimes\mu_C\|_{TV}
 =1-{d\over M}\longrightarrow1,                         \tag{6.4}
\]

and

\[
 \chi^2(\gamma\|\mu_B\otimes\mu_C)
 ={M\over d}-1={q^2\over q+1}.                          \tag{6.5}
\]

This is exactly the accumulation of `M-1` correlation modes, each of size
`sqrt(q)/(q+1)`: here `sqrt(r) epsilon` is of order `sqrt(q)`, not small.

There is an explicit Euclidean realization with isotropic endpoint
marginals.  Work in

\[
 H=\{x\in\mathbb R^M:\langle x,{\bf1}\rangle=0\}.
\]

For a point `i` and a line `ell`, let `1_ell` be its incidence vector and
put

\[
 b_i=\sqrt M\left(e_i-{{\bf1}\over M}\right),
 \qquad
 c_\ell=\sqrt{M/q}\left(1_\ell-{d\over M}{\bf1}\right). \tag{6.6}
\]

Both uniform endpoint laws are centered and isotropic on `H`:

\[
 \operatorname {Cov}(B)=\operatorname {Cov}(C)=I_H.     \tag{6.7}
\]

Their norms are constant, `|b_i|^2=|c_ell|^2=M-1`.  Their inner products
are

\[
 \langle b_i,c_\ell\rangle=
 \begin{cases}
 q^{3/2},&i\in\ell,\\
 -(q+1)/\sqrt q,&i\notin\ell.
 \end{cases}                                             \tag{6.8}
\]

Hence every incidence edge has the same length

\[
 D^2=2(M-1)-2q^{3/2}=2q(q+1-\sqrt q),                  \tag{6.9}
\]

while every nonedge has squared length

\[
 D^2+{2M\over\sqrt q}>D^2.                              \tag{6.10}
\]

Assigning values `+D/2` to the `b_i` and `-D/2` to the `c_ell` is therefore
one-Lipschitz on the endpoint set, and a McShane extension calibrates every
incidence segment.

For the product endpoint law, incidence has probability `d/M`.  Equations
(6.9)--(6.10) give

\[
 \boxed{\quad
 \operatorname {Var}_{\mu_B\otimes\mu_C}|B-C|^2
 ={d\over M}\left(1-{d\over M}\right)
       \left({2M\over\sqrt q}\right)^2
 =4(M-1).
 \quad}                                                   \tag{6.11}
\]

Equivalently, the last term in the ANOVA identity is already
`4 Tr(I_H I_H)=4(M-1)`.  The endpoint covariance spans coincide rather than
being approximately orthogonal.

The calibrated ray normals

\[
 N_{i\ell}={b_i-c_\ell\over D},\qquad i\in\ell,          \tag{6.12}
\]

are mean zero.  Maximal correlation and (6.7) give

\[
 E NN^T\preceq {2(1+\sqrt q/(q+1))\over D^2}I_H
 \preceq {4\over M-1}I_H.                               \tag{6.13}
\]

Thus the normals are highly dispersed, not aligned.  If one puts the
uniform conditional on `[-D/2,D/2]` on every labelled edge, every ray is
exactly balanced, has standard deviation `D/sqrt(12) asymp sqrt(M)`, and

\[
 E[\sigma^2NN^T]\preceq {1\over3}I_H.                  \tag{6.14}
\]

So the model also satisfies the long-ray covariance scaling.  It is an
atomic labelled-ray construction, not the disintegration of one
full-dimensional log-concave law and not asserted to be a signed-distance
extremizer.  It proves that endpoint incidence, Euclidean calibration,
balance, long scale, covariance dispersion, and a heat-bath gap tending to
one still do not imply the proposed alignment/Clifford dichotomy abstractly.
The missing hypothesis is precisely global log-concave signed-distance
realizability (or an effective-rank substitute).

### 6.1 Tensorization meets the coherent-cap constraint

The single block has only polynomially many ray labels and therefore does
not pass the exponential coherent-direction test at its scale.  A product
of independent blocks does.

Take `k` orthogonal copies of `H`.  Put

\[
 \mathbf B=(B_1,\ldots,B_k),\qquad
 \mathbf C=(C_1,\ldots,C_k),                            \tag{6.15}
\]

with joint endpoint law `gamma^{tensor k}`.  A ray label is a `k`-tuple of
incidences.  Set

\[
 n=k(M-1)=kq(q+1),qquad
 D_k=\sqrt{k}\,D.                                       \tag{6.16}
\]

Then both endpoint marginals are isotropic on `H^{oplus k}`, every joint
edge has length `D_k`, and

\[
 \operatorname {Var}_{\mu_{\mathbf B}\otimes\mu_{\mathbf C}}
       |\mathbf B-\mathbf C|^2=4k(M-1)=4n.              \tag{6.17}
\]

The singular values of a tensor-product conditional expectation are all
products of single-block singular values.  A function of one block realizes
the largest nonconstant one, so

\[
 \lambda_{HB}^{(k)}=1-{\sqrt q\over q+1}.               \tag{6.18}
\]

Also

\[
 \|\gamma^{\otimes k}
       -\mu_{\mathbf B}\otimes\mu_{\mathbf C}\|_{TV}
 =1-\left({d\over M}\right)^k,                          \tag{6.19}
\]

because the joint law is supported on coordinatewise incidences, whose
product-law probability is `(d/M)^k`.

Let `L=Md` be the number of single-block ray labels.  The tensor normal is
the orthogonal concatenation

\[
 \mathbf N_{(e_1,\ldots,e_k)}
 ={1\over\sqrt k}(N_{e_1},\ldots,N_{e_k}).               \tag{6.20}
\]

The required coherent-cap estimate follows from a Bennett bound, including
for the chordal unit caps relevant to the small-event barycenter argument.
Fix a unit vector `v=(v_1,\ldots,v_k)` in `H^{oplus k}` and write

\[
 a_j=\sqrt k\,|v_j|,\qquad
 X_j=\left\langle N_{e_j},{v_j\over|v_j|}\right\rangle,
 \qquad \sum_j a_j^2=k,                                \tag{6.21}
\]

with `X_j=0` if `v_j=0`.  The variables `X_j` are independent, mean zero,
bounded above by one, and (6.13) gives

\[
 \operatorname {Var}X_j\leq {4\over M-1}\leq{4\over q^2}.
                                                               \tag{6.22}
\]

For any fixed `gamma in (0,1)`, membership in the angular cap
`<mathbf N,v>>=gamma` is equivalent to

\[
 \sum_j a_jX_j\geq\gamma k.                            \tag{6.23}
\]

Put `A=2/gamma` and discard the coordinates with `a_j>A`.  Their maximum
possible contribution is at most

\[
 \sum_{a_j>A}a_j
 \leq\sqrt{\#\{a_j>A\}\sum_j a_j^2}
 \leq {k\over A}={\gamma k\over2}.
\]

On the remaining coordinates, `Y_j=a_jX_j` has mean zero, is bounded above
by `A`, and its total variance is at most `V_0=4k/q^2`.  Bennett's
inequality, with
`h(u)=(1+u)log(1+u)-u`, therefore yields the uniform cap bound

\[
 \begin{split}
 P\{\langle\mathbf N,v\rangle\geq\gamma\}
 &\leq
 \exp\left\{-{\gamma^2k\over q^2}
                   h\left({q^2\over4}\right)\right\}\\
 &\leq
 \exp\left\{-{\gamma^2k\over8}
                  \log\left(1+{q^2\over4}\right)\right\}.
                                                               \tag{6.24}
 \end{split}
\]

Here the second line uses `h(u)>=(u/2)log(1+u)` for `u>=1`.
A chordal unit cap `{u:|u-v|<=1}` has
`<u,v>>=1/2`, so (6.24) applies with `gamma=1/2`.

Choose, for `q>=3`,

\[
 k=\left\lceil {q^2\over(\log q)^2}\right\rceil,
 \qquad
 s_k={D_k\over\sqrt{12}}.                               \tag{6.25}
\]

Then

\[
 s_k\asymp q\sqrt k\asymp {q^2\over\log q},
 \qquad n\asymp s_k^2,
 \qquad k\log L\asymp s_k.                             \tag{6.26}
\]

For every fixed `gamma in (0,1)`, (6.24) is therefore at most

\[
 \exp(-c_\gamma s_k).                                    \tag{6.27}
\]

Tensorization also preserves the ray covariance constraint.  Independence
of the blocks and the two bounds in (6.13) give

\[
 E\mathbf N\mathbf N^T\preceq {4\over n}I,
 \qquad
 s_k^2E\mathbf N\mathbf N^T
 \preceq {1+\sqrt q/(q+1)\over6}I
 \preceq {1\over3}I.                                  \tag{6.28}
\]

Thus the tensorized model simultaneously has ray scale `s_k`, effective
normal rank `Theta(s_k^2)`, exponentially small coherent direction caps,
and heat-bath gap tending to one.  The endpoint-incidence constraints still
do not force alignment or Clifford geometry.

### 6.2 Convexification creates abundant noncalibrated bridges

Each single-block endpoint set is a regular simplex: distinct `b_i` have
inner product `-1`, and the same is true of distinct `c_ell`.  The second
simplex is an orthogonal image of the first on `H`.  In the tensor product,

\[
 0\in\operatorname {conv}\{\mathbf b\}
   \cap\operatorname {conv}\{\mathbf c\}.              \tag{6.29}
\]

Hence the positive and negative endpoint convex hulls have distance zero,
although every calibrated matched chord has length `D_k asymp s_k`.  The
body

\[
 K=\operatorname {conv}
      (\{\mathbf b\}\cup\{\mathbf c\})                 \tag{6.30}
\]

is full-dimensional in `H^{oplus k}` and its uniform law is log-concave,
but the finite union of incidence segments has Lebesgue measure zero.  Thus
the most direct convex-hull realization creates all the bridges and gives
zero quotient mass to the proposed long rays.

The nonedge slack can be computed exactly.  For a product pair
`(mathbf b,mathbf c)`, let `h` be the number of coordinates which are not
incidences.  Then

\[
 |\mathbf b-\mathbf c|^2=D_k^2+h\Delta,
 \qquad \Delta={2M\over\sqrt q}.                        \tag{6.31}
\]

Under the product endpoint law,
`h` is binomial with parameter `1-d/M`.  For `q>=3`, Hoeffding gives

\[
 P\{h\geq k/2\}\geq1-e^{-k/20}.                        \tag{6.32}
\]

For every such pair,

\[
 |\mathbf b-\mathbf c|-D_k
 ={h\Delta\over
       \sqrt{D_k^2+h\Delta}+D_k}
 \geq {\sqrt{kq}\over6}.                               \tag{6.33}
\]

To check the last inequality, use `h>=k/2`,
`D_k<=2q sqrt(k)`, and
`sqrt(D_k^2+k Delta)+D_k<=3D_k`, valid for `q>=3`.
With (6.25), the bridge slack in (6.33) is of order
`q^(3/2)/log q`; it is smaller than the calibrated length by `q^(-1/2)`
but diverges absolutely.  Convexification therefore supplies overwhelmingly
many noncalibrated cross chords rather than alternative calibrated rays.

### 6.3 Radial variance excludes a log-concave long-ray realization

There is a stronger global obstruction than convex-hull volume: the matched
segments themselves violate translated thin shell as soon as they carry
constant ray mass with long log-concave conditionals.

**Lemma 6.1 (equal-radius chord forces radial variance).**  Let
`x_-,x_+` and `z` satisfy

\[
 |x_--z|=|x_+-z|=R,\qquad D=|x_+-x_-|,                 \tag{6.34}
\]

put `m=(x_++x_-)/2`, `N=(x_+-x_-)/D`, and let
`X=m+TN` be supported on the chord.  If the density of `T` is log-concave
with variance `sigma^2`, then

\[
 \boxed{\quad
 \operatorname {Var}|X-z|\geq{\sigma^4\over400R^2}.
 \quad}                                                   \tag{6.35}
\]

**Proof.**  Translate `z` to the origin, continuing to write `X,m` for the
translated variables.  Equal endpoint radii then give `<m,N>=0`, so

\[
 |X|^2=|m|^2+T^2.                                      \tag{6.36}
\]

For every one-dimensional log-concave variable,
`Var(T^2)>=sigma^4/100`.  One proof standardizes `T`, uses the universal
upper bound on its density, and observes that the sublevel set of
`|t^2+at-1|` has uniformly bounded length for every centering parameter
`a`; this is the elementary quadratic anti-concentration lemma.
Since the chord lies in `R B_2^n`,

\[
 \operatorname {Var}(|X|^2)
 ={1\over2}E[(|X|-|X'|)^2(|X|+|X'|)^2]
 \leq4R^2\operatorname {Var}|X|.                       \tag{6.37}
\]

Combine the two estimates.  QED.

For the tensorized incidence rays,

\[
 R^2=n,qquad
 D_k^2=2n\left(1-{\sqrt q\over q+1}\right)\geq n.      \tag{6.38}
\]

Consequently, if a conditional has

\[
 \sigma_y\geq\kappa D_k,                               \tag{6.39}
\]

then

\[
 \operatorname {Var}(|X|\mid Y=y)
 \geq {\kappa^4D_k^4\over400n}
 \geq {\kappa^4 n\over400}.                            \tag{6.40}
\]

If an isotropic log-concave law assigned quotient mass at least `alpha` to
such rays, total variance would give

\[
 \operatorname {Var}_\mu|X|
 \geq {\alpha\kappa^4 n\over400}.                       \tag{6.41}
\]

The translated thin-shell estimate

\[
 \sup_z\operatorname {Var}_\mu|X-z|\leq C_{TS}         \tag{6.42}
\]

contradicts (6.41) once
`n>400C_TS/(alpha kappa^4)`.  The center in (6.34) need not be the
barycenter of `mu`, because (6.42) is uniform in translations.

For the uniform conditional on the whole chord,
`kappa=1/sqrt(12)` and `sigma=s_k`.  Hence the tensorized countermodel passes
the covariance, coherent-cap, balance, endpoint-incidence, and heat-bath
tests but cannot be realized as a constant-mass ray family of one isotropic
log-concave law.  Uniform measure on the convex hull (6.30) is a genuine
full-dimensional log-concave realization of the endpoints and their
bridges, but it assigns zero mass to the ray skeleton.  Any attempted
thickening that preserves constant quotient mass and conditional scale
`Theta(D_k)` is ruled out by (6.41), independently of how the remaining
convex-hull volume is filled.

## 7. Smooth curvature and variable endpoint lengths

Three bookkeeping issues prevent direct use of the ideal theorem.

### 7.1 Smooth charge is a second form, not a small scalar error

On a regular normal chart the smooth second variation contains the full
base derivative of the height and the shape operator.  Long balanced rays
give

\[
 \sigma_y^2\|S_y\|_{HS}^2\leq C_\delta,                 \tag{7.1}
\]

but the resulting smooth translation charge has the same `1/sigma_y`
homogeneity as `b_y`.  It need not satisfy the operator estimate in (4.2)
for arbitrary endpoint heights.  The exact trace identity (4.4) combines
smooth and singular charge; it does not show that the smooth part is
negligible before taking the endpoint spectral gap.

### 7.2 Variable lengths change both the stationary law and the charge

The stationary ray weight is `b_y d eta`, and for a log-concave conditional
`b_y` is only comparable to `1/sigma_y` when the balance parameter is fixed.
Thus different ray lengths change the endpoint marginals in (1.4).  The
one-dimensional smooth deficit also contains the finite-endpoint terms

\[
 q\,q_y(b_y-)+p\,q_y(a_y+),                              \tag{7.2}
\]

as well as the density-convexity charge.  These terms remove budget from the
medial form but are not degrees of either endpoint heat bath.

Geometrically, if the endpoint half-lengths are functions `r(B)` and `s(C)`,
same-ray calibration gives

\[
 |B-C|=r(B)+s(C),                                        \tag{7.3}
\]

not a constant cross distance.  Even exact endpoint independence therefore
does not invoke the ordinary Clifford classification unless `r,s` are
constant (or quantitatively concentrated).  For example, take arbitrary
finite sets on two separated intervals of one line, with every `b` to the
right of every `c`, and put `r(b)=b-z_0`, `s(c)=z_0-c`.  Then (7.3) holds for
the full product law and the ideal heat-bath gap is one, but cross distances
vary.  This is the aligned branch: all normals are the same.

More generally, four instances of (7.3), after squaring, give

\[
 \langle b-b',c-c'\rangle
 =-(r(b)-r(b'))(s(c)-s(c')).                             \tag{7.4}
\]

Thus variable lengths produce orthogonality only after a one-dimensional
Minkowski-type lift; they do not give the Euclidean Clifford identity
`<b-b',c-c'>=0`.

### 7.3 Endpoint conductance is not endpoint probability

The pushforwards of `pi` under `e_+` and `e_-` record boundary-weighted ray
mass.  The interface measure (2.2) records density on medial strata divided
by angular jump.  Its fiber degrees need not equal either endpoint
pushforward.  Lemma 2.2 can force a reversible kernel by adding self-loops,
but the required coefficient is

\[
 c_x\geq\max_i{\sum_jw_{ij}\over\pi_i^x},               \tag{7.5}
\]

which is uncontrolled and diverges in Example 4.2.  This is the exact
coefficient mismatch hidden by the phrase “resample at an endpoint.”

## 8. Verdict

The endpoint heat-bath route is valid under the following additional
package:

1. the smooth form is a small operator perturbation of boundary variance;
2. each medial/focal fiber has product conductances in its conditional ray
   weights, with the correctly normalized coefficient;
3. finite endpoints carry essentially all relevant ray mass;
4. endpoint half-lengths are concentrated; and
5. the endpoint correlation operator has dimension-free effective rank (or
   another hypothesis converting operator mixing into product-law mixing).

Under these assumptions, Proposition 4.1 and Theorem 5.2 give the requested
quantitative dichotomy: small normal variance gives alignment, while the
nonaligned branch has a near-product constant-cross-distance endpoint law
and hence approximate Clifford geometry.

None of items 1, 2, 4, or 5 follows from signed-distance stability and the
exact total charge.  Example 4.2 disproves coefficient recovery, and the
tensorized projective-plane model disproves dimension-free near-gap rigidity
even with ideal heat-bath resampling, long balanced Euclidean rays, dispersed
normals, and exponentially small coherent caps.  On the other hand, Lemma
6.1 and translated thin shell rule out assigning constant log-concave mass
at comparable conditional scale to those equal-radius chords.  Therefore a
proof for a genuine extremizer must use global log-concavity and the common
signed-distance realization beyond the endpoint incidence spectrum; the
abstract heat-bath gap alone is insufficient, while this particular
high-rank countermodel cannot survive the global realization step.
