# Endpoint cells and the limits of a convex inradius inverse theorem

## Outcome

There is a clean dimension-free inverse theorem, but only in the
**covariance-saturated** regime.  If a family of long ray conditionals uses
almost all of the ambient covariance on a subspace `P`, then the projected
ray lines are almost concurrent.  Several saturated families have almost
orthogonal covariance subspaces.  Exact saturation gives the expected
orthogonal-radial/cylindrical geometry on every connected smooth chart.

Large two-sided reach, a John decomposition of the normals, and a large
convex inradius do **not** imply such saturation.  Three sharp obstructions
are proved below.

1. A centered isotropic normal measure defining a polar cell gives only the
   sharp containment
   `r B_2^n subset K subset n r B_2^n`; even antipodal symmetry improves
   `n` only to `sqrt(n)`.  Regular simplices and cubes attain the losses.
2. On a regular simplex, a fixed-reach core in every facet has normals in
   exact isotropic position, while the rms distance of its normal lines from
   every common center is `sqrt(n-1)` times the reach.
3. For `K=[-1,1]^n+R B_2^n`, uniform measure becomes isotropic after a scalar
   rescaling, every selected facet normal has two-sided reach `R`, and the
   normals are the exact tight frame `+/- e_i`.  Taking `R=n^(1/4)` gives
   isotropic reach at least `n^(1/4)/2`, but the radial distortion and the
   concurrence error divided by reach are both of order `n^(1/4)`.

The third example is not raywise balanced for the uniform law.  This is the
point: every proposed input coming only from convex endpoint cells, rolling
balls, polar bodies, John decompositions, isotropy, and dispersed normals is
already present, so any valid theorem for long balanced rays must use the
balance together with positive conditional mass in an essential way.

For comparison, an exact checkerboard signed-distance function under the
isotropic cube law has balanced rays and the same dispersed normal frame, but
its constant-mass ray reach is only `Theta(1/n)`.  Thus it also identifies the
endpoint-cell competition that the word "long" must overcome.

## 1. Setup: balance gives reach, and covariance gives capacity

Let an isotropic probability `mu` be disintegrated on affine calibrated
lines.  With `Y` the ray label, write

\[
 X=Z_Y+T U_Y=M_Y+S U_Y,
 \qquad |U_Y|=1,
\]

where

\[
 m_Y=E[T\mid Y],\qquad M_Y=Z_Y+m_YU_Y=E[X\mid Y],
 \qquad S=T-m_Y,qquad \sigma_Y^2=E[S^2\mid Y].
\]

The total-covariance identity is

\[
 I=\operatorname {Cov}(M_Y)+E[\sigma_Y^2U_YU_Y^T].       \tag{1.1}
\]

In particular, for a quotient event `Omega`,

\[
 A_\Omega:=E[1_\Omega\sigma_Y^2U_YU_Y^T]\preceq I.
                                                               \tag{1.2}
\]

This is the exact covariance capacity used below.

**Lemma 1.1 (balanced scale gives two-sided reach, with constants).**
Suppose the conditional density `q_Y` is log-concave, has standard deviation
`sigma_Y`, and

\[
 P(T\leq0\mid Y)\geq\delta,\qquad
 P(T\geq0\mid Y)\geq\delta.                               \tag{1.3}
\]

Then the support of `q_Y` extends at least

\[
 r_Y={\delta\sigma_Y\over10}                              \tag{1.4}
\]

on each side of zero.  Hence a calibrated ray has geometric two-sided reach
at least `r_Y`.

**Proof.**  A one-dimensional log-concave probability density of standard
deviation `sigma` obeys `||q||_infinity <= 10/sigma`.  (Translate a mode to
zero; log-concavity bounds each tail by the supporting exponential through
the first `1/e` height point, and the second moment about the mode is at most
`94/||q||_infinity^2`.)  Each half-support in (1.3) therefore has length at
least `delta/||q||_infinity >= delta sigma/10`.  The support is an interval,
which proves (1.4).  QED.

Thus balance supplies reach, but it does not by itself say that (1.2) is
close to equality on any subspace.

## 2. The covariance-saturation inverse theorem

The following is the precise positive result.

**Theorem 2.1 (saturation forces concurrent orthogonal blocks).**
Let `Omega_1,...,Omega_m` be disjoint quotient events and put

\[
 A_i=E[1_{\Omega_i}\sigma_Y^2U_YU_Y^T].                  \tag{2.1}
\]

Let `P_i` be orthogonal projections of ranks `k_i`.  Assume, for some
`0<=epsilon<1`,

\[
 A_i\succeq(1-\epsilon)P_i\qquad(1\leq i\leq m).         \tag{2.2}
\]

Then

\[
 E|P_iM_Y|^2\leq\epsilon k_i,                            \tag{2.3}
\]

and consequently

\[
 \int_{\Omega_i}\operatorname {dist}
       (0,P_iL_Y)^2\,d\eta(Y)\leq\epsilon k_i,           \tag{2.4}
\]

where `L_Y=Z_Y+span(U_Y)` is the affine ray line.  Moreover,

\[
 \sum_{j\ne i}\|P_iP_j\|_{HS}^2
 \leq {\epsilon\over1-\epsilon}k_i.                     \tag{2.5}
\]

In particular, exact saturation makes the `P_i` mutually orthogonal and
makes every projected line `P_iL_Y`, `Y in Omega_i`, pass through the origin
almost surely.

**Proof.**  Centering isotropic `mu` gives `EM_Y=EX=0`.  From (1.1),

\[
 \begin{split}
 E|P_iM_Y|^2
 &=\operatorname {Tr}(P_i\operatorname {Cov}M_Y)\\
 &=k_i-\operatorname {Tr}
       \left(P_iE[\sigma_Y^2U_YU_Y^T]\right)\\
 &\leq k_i-\operatorname {Tr}(P_iA_i)\leq\epsilon k_i.
 \end{split}                                             \tag{2.6}
\]

The point `P_iM_Y` lies on `P_iL_Y`, proving (2.4).  Also
`sum_j A_j <= E[sigma_Y^2U_YU_Y^T] <= I`, so

\[
 \sum_{j\ne i}\operatorname {Tr}(P_iA_j)
 \leq k_i-\operatorname {Tr}(P_iA_i)\leq\epsilon k_i.   \tag{2.7}
\]

Since `A_j-(1-epsilon)P_j` is positive semidefinite,

\[
 \operatorname {Tr}(P_iA_j)
 \geq(1-\epsilon)\operatorname {Tr}(P_iP_j)
 =(1-\epsilon)\|P_iP_j\|_{HS}^2.                        \tag{2.8}
\]

Summing proves (2.5).  QED.

Here is a directly checkable sufficient condition for (2.2).

**Corollary 2.2 (effective-dimension formulation).**
Let `eta(Omega)=alpha`, suppose `U_Y in ran(P)` and `sigma_Y>=s` on
`Omega`, and assume

\[
 E[U_YU_Y^T\mid\Omega]\succeq {1-\delta\over k}P,
 \qquad {\alpha s^2\over k}\geq1-\gamma.                 \tag{2.9}
\]

If `delta+gamma<1`, then Theorem 2.1 applies with

\[
 \epsilon=\delta+\gamma,                                \tag{2.10}
\]

and

\[
 \left(E[\operatorname {dist}(0,PL_Y)^2\mid\Omega]
 \right)^{1/2}
 \leq s\sqrt{{\delta+\gamma\over1-\gamma}}.             \tag{2.11}
\]

**Proof.**  Equations (2.9) give

\[
 A_\Omega\succeq {\alpha s^2(1-\delta)\over k}P
 \succeq(1-\delta)(1-\gamma)P
 \succeq(1-\delta-\gamma)P.                            \tag{2.12}
\]

Divide (2.4) by `alpha` and use
`k/alpha <= s^2/(1-gamma)`.  QED.

The ratio `alpha s^2/k` is decisive.  Isotropy always forces it to be at
most a constant; the inverse theorem requires it to be close to one.
Merely saying that the normals are dispersed gives a **lower** bound
`k >= alpha s^2`, and supplies no converse bound.  Thus dispersion alone
does not enter the saturation regime.

There is also an exact geometric interpretation.  Suppose the exact case
`delta=gamma=0` holds in Corollary 2.2 and a connected `C^1` zero-set chart
has normals in `ran(P)`.
Since `P M_Y=0`,

\[
 PZ_Y=-m_YU_Y.                                          \tag{2.13}
\]

Differentiate on the chart and take the scalar product with the normal
`U_Y`.  The normality of `U_Y` and `dU_Y perpendicular U_Y` give `dm_Y=0`.
Hence `m_Y` is constant and `PZ_Y` lies on one sphere in `ran(P)`.  The
unprojected directions are free only in `P^perp`; locally this is exactly a
radial cylinder.  With several saturated events, (2.5) makes these radial
blocks orthogonal.  This proves the exact radial/orthogonal-block conclusion
at saturation, without KLS.

## 3. Polar and John data have an unavoidable dimension loss

The next lemma is a sharp no-go for a polar-body proof based only on normal
covariance.

**Lemma 3.1 (sharp polar cell bound).**
Let `nu` be a probability on `S^{n-1}` such that

\[
 E_\nu U=0,
 \qquad E_\nu UU^T\succeq {\lambda\over n}I             \tag{3.1}
\]

for `0<lambda<=1`, and define

\[
 K_r=\{x:\langle u,x\rangle\leq r
                  \text{ for }\nu\text{-a.e. }u\}.     \tag{3.2}
\]

Then

\[
 rB_2^n\subset K_r\subset {nr\over\lambda}B_2^n.        \tag{3.3}
\]

If `nu` is antipodally symmetric, the second inclusion improves to

\[
 K_r\subset r\sqrt{n/\lambda}\,B_2^n.                  \tag{3.4}
\]

Both dimension dependences are sharp when `lambda=1`.

**Proof.**  Fix a unit vector `theta`, set
`A=<U,theta>`, and let `M=ess sup A`.  Since `-1<=A<=M`,

\[
 (M-A)(A+1)\geq0
 \quad\Longrightarrow\quad
 A^2\leq(M-1)A+M.                                       \tag{3.5}
\]

Taking expectations and using `EA=0` gives

\[
 M\geq EA^2\geq\lambda/n.                              \tag{3.6}
\]

If `x=|x|theta in K_r`, then `|x|M<=r`, proving (3.3).
The first inclusion is immediate from `|<u,x>|<=|x|`.
Under antipodal symmetry,
`ess sup A=ess sup |A| >= sqrt(EA^2)`, which proves (3.4).  QED.

For sharpness in (3.3), take the `n+1` regular-simplex normals
`u_0,...,u_n`, with

\[
 \langle u_i,u_j\rangle=-1/n\quad(i\ne j),
 \quad \sum_i u_i=0,
 \quad {1\over n+1}\sum_i u_iu_i^T={1\over n}I.         \tag{3.7}
\]

In direction `-u_0`, the largest scalar product with the support is `1/n`,
and `K_r` has the vertex `-nr u_0`.  For sharpness in (3.4), take `nu`
uniform on `+/-e_i`; then `K_r=[-r,r]^n`, whose circumradius is `sqrt(n)r`.
Thus neither a John decomposition nor passage to the polar body can produce
a dimension-free radial approximation from dispersed normals and inradius.

## 4. Exact endpoint-cell reach in a polytope

There is an exact formula showing what an endpoint cell does control.

**Lemma 4.1 (facet core formula).**
Let

\[
 K=\bigcap_j\{x:\langle u_j,x\rangle\leq a_j\},
 \qquad |u_j|=1,                                        \tag{4.1}
\]

and let `z` be in the relative interior of facet `i`, so
`<u_i,z>=a_i`.  Put

\[
 d_j(z)=a_j-\langle u_j,z\rangle,
 \qquad c_{ij}=\langle u_i,u_j\rangle.                  \tag{4.2}
\]

The inward normal segment `z-su_i`, `0<=s<=rho`, remains calibrated for the
signed distance to `partial K` if and only if

\[
 d_j(z)\geq\rho(1-c_{ij})\qquad(j\ne i).                \tag{4.3}
\]

Consequently its maximal inward reach is

\[
 \rho_i(z)=\min_{j:c_{ij}<1}{d_j(z)\over1-c_{ij}}.       \tag{4.4}
\]

If facets `i,j` meet at angle `theta_ij=arccos(c_ij)`, then (4.3) says that
`z` has tangential distance at least

\[
 \rho\tan(\theta_{ij}/2)                                \tag{4.5}
\]

from their ridge inside facet `i`.

**Proof.**  At `x=z-su_i`, the distances to the supporting hyperplanes are

\[
 s,\qquad d_j(z)+s c_{ij}\quad(j\ne i).                 \tag{4.6}
\]

The ball whose radius is the minimum of these quantities lies in every
defining halfspace, so this minimum is exactly `dist(x,partial K)`.
It equals `s` for every `s<=rho` exactly when (4.3) holds.  Within facet `i`,
the distance to the `ij` ridge is
`d_j(z)/sqrt(1-c_ij^2)`; substitute (4.3) and use
`(1-c)/sqrt(1-c^2)=tan(theta/2)`.  QED.

Formula (4.5) is the valid cell-inradius statement.  It still does not imply
concurrence, as the sharp simplex computation shows.

**Example 4.2 (simplex facet cores are `sqrt(n)`-nonconcurrent).**
Use the normals (3.7) and set

\[
 K=\{x:\langle u_i,x\rangle\leq R, 0\leq i\leq n\}.
                                                               \tag{4.7}
\]

Facet `F_i` has contact point `c_i=Ru_i` and vertices
`v_j=-nRu_j`, `j\ne i`.  For `0<=s<R`, Lemma 4.1 gives the exact reach core

\[
 F_i(s):=\{z\in F_i:\rho_i(z)\geq s\}
 =c_i+(1-s/R)(F_i-c_i).                                  \tag{4.8}
\]

Choose `I` uniformly from the facets, choose `Z` uniformly in `F_I(s)`, and
put `U=u_I`.  Then

\[
 EUU^T=I/n.                                             \tag{4.9}
\]

For every proposed common center `a`, the corresponding normal lines obey

\[
 E\operatorname {dist}(a,Z+\mathbb R U)^2
 =(1-s/R)^2(n-1)R^2+(1-1/n)|a|^2.                      \tag{4.10}
\]

In particular, at `s=R/2`, every line has two-sided reach at least `R/2`,
but the optimal rms concurrence error is

\[
 {R\over2}\sqrt{n-1}=\sqrt{n-1}\,s.                   \tag{4.11}
\]

**Proof.**  Since all competing inner products are `-1/n`, (4.3) turns the
facet inequalities into the homothety (4.8).  A uniform point in the
`(n-1)`-simplex `F_i` has tangential second moment about `c_i`

\[
 {1\over n(n+1)}\sum_{j\ne i}|v_j-c_i|^2
 =(n-1)R^2,                                             \tag{4.12}
\]

because `|v_j-c_i|^2=(n^2-1)R^2`.  Homothety gives the first term of
(4.10).  Averaging the center term uses (4.9):
`E|P_{U^perp}a|^2=(1-1/n)|a|^2`.  QED.

This is a scale-free counterexample to an `O(reach)` concurrence conclusion.
At `s=R/2`, however, (4.8) shows that the core occupies only
`2^{-(n-1)}` of each facet's relative volume.  Thus this construction also
shows exactly where a positive quotient-mass hypothesis could add content.
It is not an isotropically long-ray example: uniform measure on (4.7) has

\[
 \operatorname {Cov}(X)={nR^2\over n+2}I,               \tag{4.13}
\]

so isotropic normalization makes `R` a constant.  The next construction
keeps reach large after isotropic normalization.

## 5. An isotropic long-reach counterexample to the convex-geometric inputs

**Proposition 5.1 (rounded cube).**
Let

\[
 Q=[-1,1]^n,\qquad K_R=Q+R B_2^n,qquad 1\leq R\leq\sqrt n,
                                                               \tag{5.1}
\]

and let `mu_R` be uniform on `K_R`.  By signed-permutation symmetry,
`Cov(mu_R)=sigma_R^2I`.  Moreover

\[
 \sigma_R\leq{\sqrt n+R\over\sqrt n}\leq2.             \tag{5.2}
\]

Thus `K_R/sigma_R`, with its uniform law, is isotropic and log-concave.
Since `B_2^n subset Q`, its Euclidean inradius after this rescaling is at
least `(R+1)/sigma_R >= (R+1)/2`.

For `I` uniform on `{1,...,n}`, `epsilon` uniform on `{+/-1}`, and
`W_I=0`, `W_j` independent uniform on `[-1,1]` for `j\ne I`, put

\[
 U=\epsilon e_I,qquad Z=\epsilon(1+R)e_I+W.             \tag{5.3}
\]

Every line `Z+R U` is normal to a flat boundary facet and has two-sided
signed-distance reach at least `R`.  Its normal law is exactly isotropic:

\[
 EU=0,qquad EUU^T=I/n.                                 \tag{5.4}
\]

After isotropic rescaling, the reach is
`rho=R/sigma_R >= R/2`, while

\[
 \inf_a E\operatorname {dist}
   (a,Z/\sigma_R+\mathbb R U)^2
 ={n-1\over3\sigma_R^2}.                               \tag{5.5}
\]

Therefore

\[
 {\text{optimal rms concurrence error}\over\rho}
 ={\sqrt{(n-1)/3}\over R}.                              \tag{5.6}
\]

Taking `R=n^(1/4)` makes this ratio asymptotic to
`n^(1/4)/sqrt(3)` even though `rho>=n^(1/4)/2`.

The body itself is also far from radial.  Its widths in a coordinate and a
diagonal direction have ratio

\[
 {R+\sqrt n\over R+1};                                  \tag{5.7}
\]

for `R=n^(1/4)` this ratio is exactly `n^(1/4)`.  Translation does not change
widths, so no translate is within a dimension-free radial distortion.

Finally, the direction family cannot be approximated by a bounded number of
parallel packets.  For any unit vectors `v_1,...,v_k`,

\[
 E\min_{1\leq j\leq k}\sin^2\angle(U,v_j)
 \geq1-{k\over n}.                                      \tag{5.8}
\]

**Proof.**  The covariance is scalar by symmetry, and
`n sigma_R^2=E|X|^2 <= (sqrt(n)+R)^2`, proving (5.2).
Write `a=epsilon e_I+W in Q`; then `Z=a+RU` and
`B(a,R) subset K_R`.  For `0<=s<=R`,

\[
 B(Z-sU,s)\subset B(a,R)\subset K_R,                    \tag{5.9}
\]

while `Z` is a boundary point.  Thus the inward distance is exactly `s`.
The supporting hyperplane with normal `U` gives exact outward calibration
for every `s>=0`.

Projection perpendicular to `U` removes the first term of `Z`; averaging the
uniform `W` gives

\[
 E\operatorname {dist}(a,Z+\mathbb R U)^2
 ={n-1\over3}+(1-1/n)|a|^2,                             \tag{5.10}
\]

which proves (5.5) after rescaling.  The support function is
`h_K(theta)=||theta||_1+R`, giving (5.7).  Lastly,

\[
 {1\over n}\sum_i\max_j\langle e_i,v_j\rangle^2
 \leq {1\over n}\sum_{i,j}\langle e_i,v_j\rangle^2
 ={k\over n},                                           \tag{5.11}
\]

which is (5.8).  QED.

This proposition satisfies the law-level assumptions “isotropic” and
“log-concave” and all the proposed convex-geometric endpoint inputs, with
reach tending to infinity.  It deliberately does not satisfy raywise
balance: the chosen zero set is the boundary of the support.  Hence it is a
sharp refutation of any step that tries to obtain radial/orthogonal-block
geometry from inradius, rolling balls, John data, polar data, and normal
dispersion before using balance.

## 6. A genuine balanced signed-distance benchmark

The following exact example shows the endpoint-cell price paid when balance
and many directions coexist in a simple log-concave law.

**Proposition 6.1 (isotropic checkerboard).**
Let `a=sqrt(3)` and let `mu` be uniform on `[-a,a]^n`; it is isotropic.  Set

\[
 \mathcal Z=\bigcup_{i=1}^n\{x_i=0\},
 \qquad
 f(x)=\operatorname {sgn}\!\left(\prod_i x_i\right)
             \min_i|x_i|.                               \tag{6.1}
\]

Then `f` is globally one-Lipschitz and is the checkerboard signed distance to
`mathcal Z`.  Away from ties, if coordinate `i` is the unique minimum and

\[
 m=\min_{j\ne i}|x_j|,                                  \tag{6.2}
\]

the maximal calibrated ray is an interval `(-m,m)` parallel to `+/-e_i`.
Its conditional law is uniform, exactly balanced at zero, and has variance
`m^2/3`.  The quotient direction law is uniform on `+/-e_i`, so

\[
 EUU^T=I/n.                                             \tag{6.3}
\]

If `eta` denotes the ray quotient probability, then for `0<=r<=a`,

\[
 \eta\{m\geq r\}
 =\left(1-{r\over a}\right)^{n-1}
       \left(1+{(n-1)r\over a}\right).                  \tag{6.4}
\]

In particular,

\[
 \eta\{m\geq a/n\}
 =\left(1-{1\over n}\right)^{n-1}\left(2-{1\over n}\right)
 \geq {3\over2e}.                                       \tag{6.5}
\]

Thus a fixed positive mass of rays is balanced and has reach at least
`sqrt(3)/n` (and standard deviation at least `1/n`), but the constant-mass
reach scale is only `Theta(1/n)`.

**Proof.**  On same-sign checkerboard cells, (6.1) is the ordinary distance
to `mathcal Z` up to sign.  A segment joining opposite-sign cells crosses
`mathcal Z`; splitting the segment at a crossing proves the global
one-Lipschitz inequality.  With all other coordinates fixed, (6.1) equals
`+/-x_i` on `(-m,m)`, and uniform measure gives the asserted conditional.

For the quotient law, Fubini weights a base point by its ray length `2m`.
If `M` is the minimum of `n-1` independent uniforms on `[0,a]`, then

\[
 P(M\geq t)=(1-t/a)^{n-1},\qquad EM=a/n.                \tag{6.6}
\]

Therefore

\[
 \eta\{m\geq r\}
 ={E[M1_{M\geq r}]\over EM}
 ={r(1-r/a)^{n-1}+\int_r^a(1-t/a)^{n-1}dt\over a/n},   \tag{6.7}
\]

which simplifies to (6.4).  Equation (6.5) uses
`(1-1/n)^(n-1)>=e^(-1)` and `2-1/n>=3/2`.  QED.

## 7. Precise conclusion for the proposed inverse route

For a constant-mass family `Omega` with `sigma_Y>=s`, balance supplies
two-sided reach `c_delta s` by Lemma 1.1, while covariance supplies

\[
 A_\Omega\preceq I,
 \qquad \operatorname {Tr}A_\Omega\geq\alpha s^2.       \tag{7.1}
\]

This forces at least `alpha s^2` covariance directions.  It does not force
the effective dimension to be comparable to `alpha s^2`, which is exactly
the saturation needed in Corollary 2.2.  Lemma 3.1, Example 4.2, and
Proposition 5.1 show that convex endpoint cells, their inradii, polar bodies,
John decompositions, rolling reach, and dispersed normals cannot supply the
missing upper bound with dimension-free constants.
In particular, a full-dimensional John frame has `k=n`, so Theorem 2.1 sees
radial rigidity only when `alpha s^2` is already comparable to `n`; for
`s^2 much smaller than n` the covariance slack is genuine.

Accordingly, the strongest proved inverse statement is:

> **Saturated inverse.**  If the long balanced rays nearly saturate isotropic
> covariance on their effective normal span, then their projected lines are
> nearly concurrent; disjoint saturated spans are nearly orthogonal; exact
> smooth equality is orthogonal-radial/cylindrical.

Without a mechanism deriving saturation from balance and positive quotient
mass, the proposed endpoint-cell/inradius *geometric inference* is sharply
false.  The counterexamples above do not refute the full long-balanced
hypothesis.  A theorem covering all constant-mass long balanced ray families
would need a new argument using positive conditional mass and the global
signed-distance gluing to prove saturation (or to bypass it).  None of the
convex-geometric inputs tested here does so.
