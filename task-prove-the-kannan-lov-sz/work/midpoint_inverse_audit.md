# Audit of the midpoint inverse step

This note tests whether the midpoint multiplicity

\[
 K(x)=\sum_{i,j}{\bf 1}_{(B_i+C_j)/2}(x)                 \tag{1}
\]

can, by itself, support an inverse theorem for a family of calibrated
transport rays.  There are two independent obstructions.  First, (1) is not
invariant under null modifications of the bands.  Second, even an exact
reflection of the two endpoint sets need not respect the calibrated pairing;
the latter distinction is indispensable if the translated thin-shell theorem
is to be used.  An exact concurrence lemma is recorded at the end.

## 1. Minkowski multiplicity is not a measure-theoretic functional

**Lemma 1 (null-set instability).**  Let `mu` be a probability on
`R^n`, `n>=2`, which is absolutely continuous and gives positive mass to some
ball `B(z,R)`.  Let `B_i,C_j`, `1<=i,j<=m`, be arbitrary compact sets.  There
are compact sets `B_i'`, `C_j'` such that

\[
 \mu(B_i\mathbin\triangle B_i')=
 \mu(C_j\mathbin\triangle C_j')=0                         \tag{2}
\]

for all `i,j`, but the corresponding midpoint multiplicity satisfies

\[
 K'(x)=m^2\qquad(x\in B(z,R)).                            \tag{3}
\]

**Proof.**  Put `H=z+R S^{n-1}`.  Absolute continuity gives `mu(H)=0`.
Set

\[
 B_i'=B_i\cup H,\qquad C_j'=C_j\cup H.
\]

The sets remain compact and (2) holds.  Moreover

\[
 (H+H)/2=B(z,R).                                         \tag{4}
\]

Indeed, for `|v|<=R`, choose unit vectors `a,b` with
`(a+b)/2=v/R` (take equal components along `v` and opposite orthogonal
components).  Then `z+v=((z+Ra)+(z+Rb))/2`.  Thus every one of the `m^2`
midpoint sets contains `B(z,R)`, proving (3).  QED.

In dimension one the same conclusion follows by using an affine copy of the
middle-thirds Cantor set `H`, since `H+H` is an interval.  Thus the issue is
not dimensional.

The ray quotient and its conditional measures are defined only up to null
sets.  Consequently no conclusion about a positive quotient-mass family of
rays can follow from (1) unless a canonical representative is imposed and
the conclusion is proved to be independent of that representative.  In
particular, an inverse theorem stated only in terms of `K` and almost-everywhere
ray data is formally ill posed.

This is not repaired by the numerical lower bound in the midpoint lemma.  The
following incidence model realizes that lower bound without selecting one
distinguished center.  On the probability space `Z_m` with uniform measure,
put

\[
 M_{ij}=\{i+j\pmod m\}.
\]

Then every pair obeys `mu(M_ij)=1/m`, while at every `x`

\[
 \sum_{i,j}{\bf 1}_{M_{ij}}(x)=m.                         \tag{5}
\]

The graph at each `x` is a different perfect matching.  Hence the exact
information delivered by the midpoint lemma is compatible with a continuum
or a large finite family of reflection centers.  A geometric theorem must
exclude this many-center alternative; it is not excluded by the incidence
count itself.

## 2. Reflection of endpoint sets is weaker than concurrence of rays

The next model is finite, but it satisfies all of the local metric, covariance,
balance, and thin-shell constraints exactly.  It isolates why a set-level
reflection conclusion would still not permit the translated radial estimate.

Fix `m>=3`, let `pi(i)=i+1 mod m`, and work in `R^m`.  Put

\[
 p_i=\sqrt m\,e_i,\qquad c_i=-\sqrt m\,e_{\pi(i)},
 \qquad a=\sqrt{m/2}.                                    \tag{6}
\]

On the `2m` point set `S={p_i,c_i}`, prescribe

\[
 f(p_i)=a,\qquad f(c_i)=-a.                              \tag{7}
\]

For every `i,j`,

\[
 |p_i-c_j|=\sqrt m\,|e_i+e_{\pi(j)}|\ge\sqrt{2m}=2a.     \tag{8}
\]

Points of the same sign impose no further constraint.  Hence (7) is
1-Lipschitz on `S` and has a global 1-Lipschitz McShane extension.  For the
paired points `(p_i,c_i)`, equality holds in (8), so their joining segments
are calibrated rays of half-length `a`.

Let `mu_0` be uniform on `S`.  It is centered and

\[
 \operatorname{Cov}(\mu_0)=I_m,
 \qquad |X|^2=m\quad\mu_0\text{-a.s.}                    \tag{9}
\]

Thus its quadratic thin-shell variance is zero.  Every paired ray gives
conditional mass one half to each sign.  Its oriented direction and zero base
are

\[
 u_i={e_i+e_{\pi(i)}\over\sqrt2},\qquad
 z_i={\sqrt m\over2}(e_i-e_{\pi(i)}).                    \tag{10}
\]

The direction law has

\[
 {1\over m}\sum_i u_i u_i^T
 ={1\over 2m}(2I+P+P^T)\preceq {2\over m}I,              \tag{11}
\]

where `P` is the cyclic permutation matrix.  It therefore has stable rank at
least `m/2`, while the ray scale is `a=sqrt(m/2)`.

The rays are neither parallel nor concurrent.  Parallelism is excluded by
(11).  The line through `p_i,c_i` is contained in
`span(e_i,e_{pi(i)})` and satisfies

\[
 x_i-x_{\pi(i)}=\sqrt m.                                 \tag{12}
\]

The intersection of all these coordinate two-planes is `{0}`, whereas (12)
excludes `0`; hence no point belongs to every line.

Nevertheless the *unpaired endpoint sets* are exact reflections:

\[
 \{c_i:1\le i\le m\}=-\{p_i:1\le i\le m\}.              \tag{13}
\]

Thus reflection about one center, if it is allowed to permute the calibrated
rays, is not the radial/concurrent branch.  Indeed all points in (6) have the
same radius, while (7) takes both signs, so `f` is not a function of the radius
about the reflection center.

The sole missing hypothesis in this model is global log-concavity.  Convexifying
the support fills the bridges between the calibrated segments, and the
segments themselves then have zero volume.  Therefore a successful global
argument must use log-concavity to produce *positive representation mass on
the same-ray matching*, not merely membership of Minkowski midpoint sets or
set-level reflection.

## 3. Exact concurrence really does reduce to translated radial structure

The following elementary lemma gives the precise form of the branch on which
translated thin shell applies.

**Lemma 2 (normal concurrence rigidity).**  Let `Omega` be a connected
`C^1` submanifold of `S^{n-1}` and let `z:Omega->R^n` be a `C^1`
parametrization of a hypersurface whose oriented unit normal at `z(u)` is
`u`.  Suppose every normal line passes through a fixed point `z_0`.  Then
there is a constant `r` such that

\[
 z(u)=z_0+r u\qquad(u\in\Omega).                          \tag{14}
\]

Consequently, on any radial normal chart on which

\[
 F(u,t)=z(u)+tu,\qquad f(F(u,t))=t,                       \tag{15}
\]

and `r+t` has a fixed sign, `f` is a signed translate of the radial
function `|x-z_0|`.

**Proof.**  Concurrence gives `z(u)=z_0+r(u)u`.  For every
`v in T_u Omega`, normality and differentiation give

\[
 0=\langle Dz(u)v,u\rangle
  =D r(u)v+r(u)\langle v,u\rangle=D r(u)v.                \tag{16}
\]

Thus `r` is constant on the connected set `Omega`, proving (14).  Formula
(15) becomes `x=z_0+(r+t)u`, whence
`t=|x-z_0|-r` if `r+t>0` and the analogous signed formula on the other
branch.  QED.

For an isotropic log-concave random vector `X`, the translated thin-shell
estimate

\[
 \sup_{z_0}\operatorname{Var}|X-z_0|\le C               \tag{17}
\]

therefore controls every exactly concurrent calibrated chart.  Likewise, an
exactly constant normal gives an affine function and is controlled by
isotropy.  The word `reflection` cannot replace `concurrence` in this
deduction, as (6)--(13) show.

## 4. What a usable replacement would have to contain

A null-invariant candidate must count actual reflected representations, for
example through a convolution-type quantity

\[
 R_{ij}(x)=\int {\bf1}_{B_i}(x+h){\bf1}_{C_j}(x-h)
       \rho(x+h)\rho(x-h)\,dh.                            \tag{18}
\]

Positive `R_ij`-mass would survive null modifications and, after localization
in the `t` coordinate, could distinguish the calibrated matching from an
arbitrary permutation.  However log-concavity gives only

\[
 \rho(x)\ge\sqrt{\rho(x+h)\rho(x-h)},                    \tag{19}
\]

which is an upper, not a lower, comparison for the integrand in (18).
There is no universal lower bound transferring `mu((B_i+C_j)/2)` to
`R_ij`: even for Gaussian opposite halfspaces, the independent-midpoint law
lives at radius about `sqrt(n/2)` while `mu` lives at radius about `sqrt n`.

Accordingly the midpoint-overlap lemma, as presently formulated, cannot yield
the requested parallel-or-concurrent inverse theorem.  A closing lemma needs
a new, null-invariant weighted bridge estimate which (i) places positive mass
on reflected endpoint representations, (ii) respects the same-ray matching,
and (iii) rules out the many-center incidence pattern (5).  Proving any such
estimate with a universal constant is additional content, not a consequence
of the existing midpoint count.
