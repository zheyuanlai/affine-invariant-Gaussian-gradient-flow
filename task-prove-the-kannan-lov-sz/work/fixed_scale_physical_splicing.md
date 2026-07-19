# Fixed-scale phase data as physical boundary flux, and a finite corner splice

## 0. Outcome

This note supplies two proved pieces that were missing from the fixed-scale
report.

1. The good Fisher matrix at `s=10^{-10}K` transfers, with explicit
   constants, to a genuine submeasure of the physical coarea boundaries of
   `T_s1_S`.  For `K>=2*10^20`, that submeasure has mass at least `.004p`
   and effective rank at least `17`.  Thus the tilt-space rank is not merely
   a nonphysical certificate.
2. A finite, null-invariant bevel of all codimension-two corners of the
   one-sided-exponential max box gives a balanced convex competitor with a
   fixed perimeter saving.  In particular, if `mu_n` is the product of `n`
   standard one-sided exponentials, then

   \[
                    \limsup_{n\to\infty}\psi_{\mu_n}
                    \le \log2-0.012.                 \tag{0.1}
   \]

   This proves the alternative (7.8) in
   `fixed_scale_wedge_extremality.md` and removes that stress test.

The general incidence-to-saving statement is still unproved.  The exact
remaining datum is now clear: the physical normal submeasure controls facet
area, while a finite bevel is paid by actual codimension-two ridge capacity.
No dimension-free inequality from the former to the latter has been proved.

All perimeters below are relaxed weighted `BV` perimeters on the affine
support.  This is the perimeter in weighted coarea and has the same Cheeger
infimum as exterior Minkowski content.  For the explicit polytopes below the
two values agree directly.

## 1. Fisher rank transfers to physical coarea normals

Let the notation and fixed instantiation (3.17)--(3.20) of
`fixed_scale_wedge_extremality.md` be in force.  Thus

\[
 \alpha=10^{-10},\quad s=\alpha K,\quad \beta\le10^{-5},
 \quad \tau<0.098,                                   \tag{1.1}
\]

and, on the central good set `G`,

\[
 \operatorname {tr}R_G\ge {10^{-5}\over8\pi},\qquad
 {\operatorname {tr}R_G\over\|R_G\|_{op}}
       \ge {10^{-15}K\over8\pi}.                    \tag{1.2}
\]

Here, with `h=g(1-g)`,

\[
 R_G=\int_Gq_s(y){v(y)v(y)^T\over s h(y)}dy,
 \quad v=\operatorname {Cov}(1_S,X\mid Y=y).        \tag{1.3}
\]

Put

\[
 W(y)=\nabla g(y)={v(y)\over s},\quad
 u(y)={W(y)\over|W(y)|},\quad
 F_0(x)=T_s1_S(x),\quad m(x)=\nabla F_0(x)=E[W(Y)\mid X=x].                         \tag{1.4}
\]

The last identity follows by dominated differentiation of the Gaussian
average.  Define `theta=m/|m|` on `{m!=0}` and arbitrarily on `{m=0}`.

### 1.1 Fixed constants

Let

\[
 a_-:=\min_{\delta_0\le r\le1-\delta_0}
             {r(1-r)\over I(r)},\qquad
 a_+:=\max_{\delta_0\le r\le1-\delta_0}
             {r(1-r)\over I(r)}.                    \tag{1.5}
\]

For the fixed cutoff in (3.18),

\[
                         a_->0.13,qquad a_+<1.26.   \tag{1.6}
\]

Here is a certificate for these decimal bounds.  Symmetry reduces to
`r<=1/2`.  If `r<=1/4`, write `r=Phi(-x)`, `0<=x<=M_0`.  The Gaussian
Mills ratio is decreasing and

\[
 {\Phi(-M_0)\over\varphi(M_0)}\ge {M_0\over1+M_0^2}.
\]

Using `5.47<M_0<5.48` and `1-r>=3/4` gives a lower bound larger than
`.132`.  If `r in [1/4,1/2]`, then
`r(1-r)/I(r)>=(3/16)/I(1/2)>.46`.  Finally the elementary profile bound
`I(r)>=sqrt(2/pi)r(1-r)` gives `a_+<=sqrt(pi/2)<1.26`.

Define the good flux matrix

\[
                         B=E[1_G|W|uu^T],\qquad b=trB.             \tag{1.7}
\]

On `G`, `e=|v|^2/[sI(g)^2]` lies in `[1-tau,1]`, and pointwise

\[
 |W|uu^T={h\over\sqrt{s e}\,I(g)}
              {vv^T\over s h}.                       \tag{1.8}
\]

Consequently

\[
 {a_-\over\sqrt s}R_G\preceq B
 \preceq {a_+\over\sqrt{s(1-\tau)}}R_G.             \tag{1.9}
\]

Since `p(1-beta)<=1/sqrt K`, (1.2) and (1.9) give

\[
 b\ge {a_-\over8\pi\sqrt K}
   \ge b_0p,qquad
 b_0:={.13(1-10^{-5})\over8\pi}>.00517.             \tag{1.10}
\]

Moreover

\[
 {b\over\|B\|_{op}}
 \ge {a_-\sqrt{1-\tau}\over a_+}
       {10^{-15}K\over8\pi}
 \ge 3.8\cdot10^{-18}K.                              \tag{1.11}
\]

Thus the last effective rank is at least `500` when `K>=2*10^20`.

### 1.2 The exact alignment identity

Put

\[
 R(x)=E[|W(Y)|\mid X=x],\qquad d(x)=R(x)-|m(x)|.     \tag{1.12}
\]

The norm-Jensen deficit in the audited heat calculation gives

\[
 \Delta:=E d=E|W|-E|m|
 \le p\left\{\beta+{4\sqrt\alpha\over
                  c_G(1-\beta)}\right\}
 <\epsilon_*p,qquad \epsilon_*=6.02\cdot10^{-5}.   \tag{1.13}
\]

There is also the exact directional identity

\[
             E[|W||u-\theta(X)|^2]=2\Delta.          \tag{1.14}
\]

Indeed, conditionally on `X`,
`E[|W|u|X]=m=|m|theta`, and expansion of the square gives
`2(R-|m|)`.  The formula remains valid when `m=0`.

Let

\[
 r_G(x)=E[1_G|W(Y)|\mid X=x],qquad
 D=E[r_G(X)\theta(X)\theta(X)^T].                    \tag{1.15}
\]

Then `trD=b`.  For every unit vector `a`,

\[
 (a\cdot\theta)^2\le2(a\cdot u)^2+2|u-\theta|^2.
\]

Equations (1.7), (1.14), and (1.15) therefore imply

\[
                         \|D\|_{op}\le2\|B\|_{op}+4\Delta.       \tag{1.16}
\]

### 1.3 A genuine physical boundary submeasure

Fix `lambda=1/10` and let

\[
 H=\{x:d(x)\le\lambda r_G(x)\}.                     \tag{1.17}
\]

Since `int_{H^c}r_G<=Delta/lambda`,

\[
 D_H:=E[1_Hr_G\theta\theta^T],\qquad
 trD_H\ge b-\Delta/\lambda.                         \tag{1.18}
\]

On `H`, `r_G<=R=|m|+d<=|m|+lambda r_G`, and hence
`(1-lambda)r_G<=|m|`.  Define

\[
 \omega(x)=1_H(x){(1-\lambda)r_G(x)\over|m(x)|}     \tag{1.19}
\]

on `{m!=0}`, and set it equal to zero on `{m=0}`.  Then `0<=omega<=1` and

\[
 M_{phys}:=\int\omega|\nabla F_0|\theta\theta^Td\mu
          =(1-\lambda)D_H.                           \tag{1.20}
\]

Using (1.10), (1.13), (1.16), and the rank `500` in (1.11),

\[
 \boxed{
 trM_{phys}>.004p,qquad
 {trM_{phys}\over\|M_{phys}\|_{op}}>17.}            \tag{1.21}
\]

For completeness, after cancelling the factor `1-lambda`, the rank lower
bound used here is

\[
 {1-\epsilon_*/(\lambda b_0)
  \over 2/500+4\epsilon_*/b_0}>17,                  \tag{1.22}
\]

and the trace bound is
`.9(b_0-10epsilon_*)p>.004p`.

This matrix is physical in the literal coarea sense.  For
`A_r={F_0>r}`, weighted coarea, applied entrywise, gives

\[
 \boxed{
 M_{phys}=\int_0^1\int_{\partial^*A_r}
       \omega(x)n_r(x)n_r(x)^T e^{-V(x)}
       d\mathcal H^{n-1}(x)dr.}                     \tag{1.23}
\]

The sign difference between `n_r` and `theta` disappears in the projector.
The formula holds on the affine support and for extended-valued convex `V`
by weighted `BV` coarea.  No conditional law and no tilt-space separator
appears in (1.23).

### 1.4 A strictly positive analytic selector

The threshold in (1.17)--(1.19) is not needed.  There is a stronger transfer
whose selector has no discontinuity.  Since `0<=r_G<=R`, define

\[
                    \omega_{an}(x)={r_G(x)\over R(x)}.           \tag{1.24}
\]

The denominator is strictly positive.  Indeed `W` is not zero almost
everywhere because `b>0`, and convolution with the strictly positive
Gaussian kernel makes `R>0` everywhere.  The same argument and `b>0` give
`r_G>0` everywhere.  Both functions are real analytic Gaussian
convolutions of bounded nonnegative functions.  Thus `omega_an` is strictly
positive and real analytic, with `0<omega_an<=1`.

Put

\[
 M_{an}:=E\left[{r_G\over R}|m|\theta\theta^T\right]
        =\int\omega_{an}|\nabla F_0|\theta\theta^T d\mu.          \tag{1.25}
\]

Comparison with `D` in (1.15) is exact:

\[
 0\preceq D-M_{an},\qquad
 tr(D-M_{an})
 =E\left[{r_G\over R}(R-|m|)\right]\le E d=\Delta.              \tag{1.26}
\]

Consequently, using (1.10), (1.13), (1.16), and
`||B||op<=b/500`,

\[
 \boxed{
 trM_{an}>.0051p,\qquad
 {trM_{an}\over\|M_{an}\|_{op}}
 \ge {b-\Delta\over2b/500+4\Delta}>19.5.}          \tag{1.27}
\]

Weighted coarea gives the analogue of (1.23) with `omega_an`.  One can also
make the selector uniformly faithful without losing the rank.  The
heat-profile estimate gives `E|m|=int|nabla F_0|dmu<=p`.  Fix
`eta=10^{-5}` and set

\[
 \bar\omega=\eta+(1-\eta)\omega_{an},\qquad
 N=E[|m|\theta\theta^T],\qquad
 \bar M=\eta N+(1-\eta)M_{an}.                     \tag{1.28}
\]

Then `eta<=bar omega<=1` globally, `bar omega` is real analytic, and
`bar M=int bar omega|nabla F_0|theta theta^T dmu`.  Equations
(1.10), (1.13), (1.16), together with `trN<=p`, give

\[
 \begin{aligned}
 tr\bar M&\ge(1-\eta)(b-\Delta)>.0051p,\\
 \|\bar M\|_{op}
 &\le\eta p+2b/500+4\Delta,\\
 {tr\bar M\over\|\bar M\|_{op}}
 &>{ (1-10^{-5})(.00517-.0000602)
       \over .00001+.004(.00517)+4(.0000602)}>18.8. \tag{1.29}
 \end{aligned}
\]

Thus all rank-`17` arguments may use a globally positive analytic selector
which is pointwise comparable to the whole coarea boundary measure within
the fixed factor `10^5`.  This removes selector amplitude wells and
threshold switching.  It still gives no dimension-free bound on the
derivative of the selector or on the lapse of the level foliation.

## 2. The finite physical-splice functional

For a regular level `A_r`, call `B_r` a **finite physical splice** if

1. `mu(B_r)=mu(A_r)`;
2. `B_r triangle A_r` is contained in a finite union of compact balls in
   the affine support;
3. in every such ball, `B_r` is obtained by replacing the canonical
   finite-perimeter representative of `A_r` by another finite-perimeter
   set with the same trace on the boundary of the ball.

This definition includes chord bevels, cap exchanges, smooth normal pushes,
and a separate volume-correction patch.  It is invariant under changes of
`A_r` on a null set.

For a measurable family of such splices put

\[
 \mathfrak G(F_0;\{B_r\})=
 \int_0^1[P_\mu(A_r)-P_\mu(B_r)]_+dr,               \tag{2.1}
\]

and let `mathfrak G_fin(F_0)` be its supremum over finite splice families.
This is not a renaming of the coarea deficit: every term is the saving of an
actual same-mass physical set.

**Lemma 2.1 (exact splice charge).**

\[
 \boxed{\mathfrak G_{fin}(F_0)\le\mathcal D_{co}(F_0).}          \tag{2.2}
\]

**Proof.**  Put `a_r=mu(A_r)`.  Every same-mass competitor obeys

\[
 P_\mu(B_r)\ge\psi\min(a_r,1-a_r).
\]

Hence

\[
 [P_\mu(A_r)-P_\mu(B_r)]_+
 \le P_\mu(A_r)-\psi\min(a_r,1-a_r).
\]

Integrate and then take the supremum.  This proof also shows why exact mass
preservation, rather than a first-order volume constraint, is required.

At the fixed scale, (5.3) of the fixed-scale report gives

\[
                 \mathfrak G_{fin}(F_0)le6.02\cdot10^{-5}p.   \tag{2.3}
\]

Thus any geometric construction producing finite splices with total saving
`10^{-4}p` closes the fixed-scale route.

## 3. A local bevel lemma

The elementary operation behind the explicit example is the following.
Suppose a weighted finite-perimeter set agrees, in a product neighborhood
`Q times B_R^2`, with `Q times W_theta`, where `W_theta` is a planar convex
wedge of interior angle `theta in (0,pi)`, and the continuous density is
positive on `Q times {0}`.  At every point of `Q`, replace the two radial
segments of length `ell` by their chord.  The old cross-sectional boundary
length is `2ell`, the new length is `2ell sin(theta/2)`, and the removed
area is `ell^2 sin(theta)/2`.  Therefore

\[
 \begin{aligned}
 P(A)-P(A_{bevel})
 &=2\ell\{1-\sin(\theta/2)\}
       \int_Q e^{-V(q,0)}d\mathcal H^{n-2}(q)+o(\ell),\\
 \mu(A)-\mu(A_{bevel})&=O(\ell^2).                 \tag{3.1}
 \end{aligned}
\]

A normal displacement of size `O(ell^2)` on any separate regular patch
restores the mass and costs only `O(ell^2)` perimeter.  Hence every genuine
nonflat convex ridge admits a finite, exactly mass-preserving splice with
strict perimeter saving.  In terms of the angle `gamma=pi-theta` between
the two outward normals, the leading saving factor is

\[
                         1-\cos(\gamma/2).           \tag{3.2}
\]

The proof is just the area formula in the product chart plus uniform
continuity of the density.  Approximate wedge charts give the same statement
with an additional `o(ell)` error.  Concave ridges are treated by filling the
corresponding chord cap.

This lemma shows the correct prospective incidence functional: facet area
is not enough; one needs the weighted codimension-two measure of ridges,
multiplied by (3.2) and by an admissible bevel radius.

### 3.1 A proved graph-expansion-to-saving inequality

The algebraic part of the incidence argument can be completed.  Let facet
packets have weighted areas `a_i>0`, unit normals `n_i`, and

\[
 A=\sum_i a_i,qquad M=\sum_i a_in_in_i^T,qquad
 R={A\over\|M\|_{op}}.                               \tag{3.3}
\]

Let `w_ij>=0` be **admissible bevel conductances**: the corresponding ridge
tubes can be chosen disjointly and the simultaneous exactly-volume-corrected
bevel has saving at least

\[
 {1\over8}\sum_{i<j}w_{ij}|n_i-n_j|^2.              \tag{3.4}
\]

The factor `1/8` leaves room for density variation and volume correction;
the ideal wedge calculation gives at least `1/4`, because

\[
 2\{1-\cos(\gamma/2)\}
 \ge {1\over4}|n_i-n_j|^2.                          \tag{3.5}
\]

Suppose the weighted ridge graph has spectral gap `lambda`, in the exact
sense that for every scalar vector `(h_i)`,

\[
 \sum_{i<j}w_{ij}(h_i-h_j)^2
 \ge\lambda\sum_i a_i(h_i-\bar h)^2,qquad
 \bar h=A^{-1}\sum_i a_ih_i.                        \tag{3.6}
\]

**Lemma 3.2 (ridge expansion forces a finite saving).**  Under
(3.3)--(3.6), the simultaneous bevel saves at least

\[
 \boxed{{\lambda A\over8}\left(1-{1\over R}\right).}           \tag{3.7}
\]

**Proof.**  Apply (3.6) to `h_i=<n_i,e_k>` and sum over an orthonormal
basis.  This gives

\[
 \sum_{i<j}w_{ij}|n_i-n_j|^2
 \ge\lambda A(1-|\bar n|^2),qquad
 \bar n=A^{-1}\sum_i a_in_i.                        \tag{3.8}
\]

Cauchy--Schwarz in the discrete measure gives

\[
 |\bar n|^2\le {\|M\|_{op}\over A}={1\over R}.      \tag{3.9}
\]

Combine (3.4), (3.8), and (3.9).

For the physical seed `A>=.004p` and `R>=17`.  Thus a ridge graph with
`lambda>=1/4` produces saving at least

\[
 {1\over32}(.004p){16\over17}>1.17\cdot10^{-4}p,    \tag{3.10}
\]

which contradicts (2.3).  This is a proved fixed phase charge whenever the
actual ridge graph is a numerical expander.  If (3.6) has small gap, a graph
threshold decomposition produces groups of facet packets with small mutual
ridge conductance provided the usual weighted-degree normalization
`sum_jw_ij<=C a_i` is also available.  Without that normalization, even this
sweep conclusion carries the corresponding `sqrt(C)` loss.  Turning a
low-conductance packet group into a genuinely separated physical component
is the geometric step addressed, but not completed, in
`profile_linearity_separation.md`.

## 4. Exact simultaneous bevel of the exponential max box

Let

\[
 d\mu_n(x)=e^{-\sum_{i=1}^nx_i}1_{\{x_i\ge0\}}dx.              \tag{4.1}
\]

For `q>r>0`, define the convex rounded box

\[
 C_n(q,r)=\{x\in[0,\infty)^n:x_i\le q\ \forall i,
                 \ x_i+x_j\le2q-r\ \forall i<j\}.             \tag{4.2}
\]

Equivalently, the largest two coordinates have sum at most `2q-r`.
This simultaneously bevels every codimension-two upper corner and handles
all triple overlaps without inclusion-exclusion.

Write `F(t)=1-e^{-t}`.  Choosing the unique largest coordinate and then
integrating the remaining coordinates gives the exact mass

\[
 \boxed{
 M_n(q,r)=F(q-r/2)^n
 +n\int_{q-r/2}^{q}e^{-t}
        F(2q-r-t)^{n-1}dt.}                         \tag{4.3}
\]

The relative weighted boundary has two types of facets.  On `x_i=q`, every
other coordinate is at most `q-r`.  On
`x_i+x_j=2q-r`, all other coordinates are at most
`min(x_i,x_j)`, and the surface Jacobian is `sqrt2`.  Thus

\[
 \boxed{\begin{aligned}
 P_n(q,r)={}&ne^{-q}F(q-r)^{n-1}\\
 &+n(n-1)\sqrt2e^{-2q+r}
       \int_{q-r}^{q-r/2}F(t)^{n-2}dt .             \tag{4.4}
 \end{aligned}}
\]

There is no contribution from `x_i=0`: the set contains the support-side
neighborhood and the density is zero on the other side.  Equations
(4.3)--(4.4) therefore compute both relative `BV` perimeter and exterior
Minkowski content.

For each fixed `r=1/10` and all sufficiently large `n`, continuity and
strict nesting in `q` give a `q_n>r` with

\[
                         M_n(q_n,r)=1/2.             \tag{4.5}
\]

Put `lambda_n=ne^{-q_n}`.  Since

\[
 [0,q_n-r]^n\subset C_n(q_n,r)\subset[0,q_n]^n,
\]

(4.5) implies

\[
 e^{-r}n(1-2^{-1/n})\le\lambda_n
       \le n(1-2^{-1/n}).                           \tag{4.6}
\]

Hence every subsequential limit `lambda` belongs to
`[e^{-r}log2,log2]`.

If `ne^{-q}\to\lambda`, (4.3)--(4.4), after writing `t=q-y`, give

\[
 \begin{aligned}
 M(\lambda,r)
 &=e^{-\lambda e^{r/2}}
   +\lambda\int_0^{r/2}e^y
          e^{-\lambda e^{r-y}}dy,                  \tag{4.7}\\
 P(\lambda,r)
 &=\lambda e^{-\lambda e^r}
   +\sqrt2\lambda^2e^r
          \int_{r/2}^{r}e^{-\lambda e^y}dy.        \tag{4.8}
 \end{aligned}
\]

Dominated convergence is uniform for `lambda` in the compact interval
(4.6).  For `0<r<=1/10` and `lambda<=log2`, one has
`lambda e^r<1`.  Differentiating (4.8) in `lambda` shows that both terms are
increasing: the derivative of the second integrand is proportional to
`lambda(2-lambda e^y)>0`.  Therefore every subsequential perimeter limit is
at most `P(log2,r)`.

The derivative calculation at the unrounded box is exact.  With
`L=log2`,

\[
 M(L,0)=1/2,\quad \partial_\lambda M(L,0)=-1/2,
 \quad \partial_rM(L,0)=0,\quad
 P(L,0)=L/2,
\]

so the implicit half-mass parameter satisfies `lambda'(0)=0`.

and

\[
 \boxed{\partial_rP(L,0)
   =-L^2e^{-L}+{\sqrt2\over2}L^2e^{-L}
   ={L^2(\sqrt2-2)\over4}<0.}                      \tag{4.9}
\]

Thus the saving is first order while the mass correction is only second
order.

For the promised fixed numerical gap, monotonicity of the integrand in
(4.8) gives, at `r=1/10`,

\[
 \begin{aligned}
 P(L,1/10)
 &\le Le^{-Le^{1/10}}
 +{L^2e^{1/10}\over10\sqrt2}
       e^{-Le^{1/20}}\\
 &<0.3404.                                          \tag{4.10}
 \end{aligned}
\]

The last decimal is certified, for example, by the Taylor bounds with
remainders using
`0.6931<L<0.6932`, `1.1051<e^{.1}<1.1052`,
`1.0512<e^{.05}<1.0513`, and
`1.4142<sqrt2<1.4143`.  Since `L/2>0.34655`, (4.10) gives

\[
 \limsup_nP_{\mu_n}(C_n(q_n,.1))< {\log2\over2}-.0061.          \tag{4.11}
\]

The sets have mass one half, so

\[
 \boxed{\limsup_n\psi_{\mu_n}
 \le2\limsup_nP_{\mu_n}(C_n(q_n,.1))
 <\log2-.012.}                                      \tag{4.12}
\]

Translation by the barycenter makes `mu_n` isotropic because each
exponential coordinate has variance one; neither perimeter nor the Cheeger
constant changes.

## 5. The weakest fixed-scale geometric lemma actually needed

Sections 1--2 reduce the missing theorem to the following fixed numerical
form.  It is strictly weaker than the variable-parameter lemma in Section 6
of the fixed-scale report.

> **Fixed physical-splice lemma (unproved).**  There is a universal
> `C_geom` such that, whenever `K>=2*10^20` and the fixed-scale construction
> produces the physical coarea submeasure (1.23) with mass `.004p` and
> effective rank `17`, either a concurrent radial organization yields
> `K<=C_geom`, or there is a measurable family of finite physical splices
> satisfying
>
> \[
>                    \mathfrak G(F_0;\{B_r\})\ge10^{-4}p.       \tag{5.1}
> \]

The affine branch need not be stated separately: a submeasure with effective
rank larger than `17` cannot put as much as `1/17` of its trace on one line.
Lemma 2.1
and (2.3) show that (5.1) is impossible.  Thus this fixed statement alone
would bound `K` universally.

The finite bevel proves (5.1) whenever a fixed fraction of the submeasure is
carried by nonflat facet packets joined through ridges of fixed weighted
capacity.  What is not proved is a lower bound on that ridge capacity from
(1.21).  The matrix (1.23) is a vertex/facet-area statistic; the bevel is an
edge/ridge statistic.  Spatially separated planar patches can have arbitrary
normal rank and zero mutual ridge capacity.  Perimeter locality supplies no
cross-normal term between them.  Global near-Cheeger extremality and
log-concavity would have to rule out precisely that separated configuration,
or classify it as a bounded product/affine/radial branch.  No such theorem is
currently derived here.

## 6. Required model tests

1. **Symmetric Laplace halfline.**  There is only one projective normal, so
   (1.21) fails.  The affine branch correctly absorbs its positive
   longitudinal Fisher energy and zero coarea deficit.
2. **Radial Gaussian ball.**  The physical normal matrix has high rank, but
   all normal lines meet at the center and there are no polyhedral ridges.
   It is exactly the concurrent radial exception.  The ambient Gaussian has
   `K=1`.
3. **Cube.**  A coordinate half-cube is affine.  An inner box has orthogonal
   facets with positive ridge capacity; (3.1) bevels every corner with strict
   saving.  Independently, a coordinate cut already has smaller perimeter.
4. **Simplex.**  A barycentric halfspace is affine.  A homothetic inner
   simplex has nonzero dihedral ridge capacity and is bevelled by (3.1).
   The known simplex spectral estimate also places the isotropic simplex in
   the bounded-`K` branch.
5. **One-sided exponential product.**  Equations (4.2)--(4.12) are a global
   simultaneous bevel, including every higher-order corner overlap.  They
   prove a fixed gap below the limiting max-box value, so this model no
   longer survives the physical-splice test.

These tests show that the finite functional has the correct exceptional
branches.  They do not prove the incidence-to-saving assertion for an
arbitrary high-`K` near-Cheeger sequence; that assertion remains the sole
unproved step in (5.1).
