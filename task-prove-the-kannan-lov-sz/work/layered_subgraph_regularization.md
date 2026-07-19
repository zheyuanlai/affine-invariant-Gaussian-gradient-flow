# Layered regularization without a jump loss: the free-law calibration

## 0. Outcome

There is a way to remove the jump/contact obstruction which is stronger than
adding an elliptic term to the subgraph. Let

\[
 R(G):=\int_{\mathbb R}
       \min\{\mu(G>r),1-\mu(G>r)\}\,dr
      =\inf_{c\in\mathbb R}\int |G-c|\,d\mu .       \tag{0.1}
\]

The second expression says that \(R\) is exactly the norm on
\(L^1(\mu)/\mathbb R\). If \(p=\psi_\mu\), then coarea and the definition of
the Cheeger constant give

\[
             \operatorname {TV}_\mu(G)-pR(G)\ge0.  \tag{0.2}
\]

For the retained heat matrix \(M(F)\), define

\[
 \begin{split}
 K(G)&=\operatorname {TV}_\mu(G)
       +\kappa\|M(G)-M(F)\|_*,\\
 J(G)&=K(G)-pR(G),\qquad 0<\kappa<1/3 .             \tag{0.3}
 \end{split}
\]

The functional \(K\) is convex, although the map \(G\mapsto M(G)\) is not
linear. The reason is the exact representation of \(K\) as a supremum of
constant-anisotropy total variations. The functional \(J\) is nonnegative.
For the clipped heat comparator, \(J(F)\) is the already audited coarea
deficit

\[
 J(F)=\int
  [P_\mu(F>r)-p\min\{\mu(F>r),1-\mu(F>r)\}]\,dr.   \tag{0.4}
\]

Ekeland's principle, applied on \(L^1(\mu)/\mathbb R\), produces a \(G\)
with \(J(G)\le J(F)\) and an exact convex subgradient

\[
                   \xi\in\partial K(G),\qquad
 \int\xi\,d\mu=0,\qquad \|\xi\|_\infty\le p+\delta. \tag{0.5}
\]

One constant nuclear subgradient \(H\) and one Anzellotti calibration \(z\)
represent this subgradient. Weighted coarea then shows that the same \(z\)
calibrates almost every superlevel set \(A_r=\{G>r\}\). Precisely,

\[
 \boxed{
 P_{\Phi_H,\mu}(A_r)=\int_{A_r}\xi\,d\mu,
 \quad
 A_r\in\arg\min_B
   \left\{P_{\Phi_H,\mu}(B)-\int_B\xi\,d\mu\right\}.} \tag{0.6}
\]

Consequently every such leaf is a global bounded-forcing quasiminimizer,

\[
 P_{\Phi_H,\mu}(A_r)
 \le P_{\Phi_H,\mu}(B)+(p+\delta)\mu(A_r\triangle B), \tag{0.7}
\]

and

\[
 p\min(v_r,1-v_r)
 \le P_\mu(A_r)
 \le {p+\delta\over1-\kappa}\min(v_r,1-v_r).       \tag{0.8}
\]

This is the desired BV contact replacement. It neither excludes the jump
and Cantor parts nor pretends that they are smooth. Instead, the calibration
saturates the full measure \(DG\), and coarea transfers saturation to the
rectifiable boundaries of almost every level. On smooth pieces the weighted
anisotropic mean curvature has absolute value at most \(p+\delta\). A jump
interface has the exact divided-difference law described in Section 4, and
its curvature is again in this interval. Thus incompatible contacts cannot
carry unbounded curvature.

At the retuned constants in the heat construction, this operation preserves
the angular variance with room to spare: for \(\kappa=10^{-6}\), the loss in
normalized matrix purity is less than \(3.67\cdot10^{-5}\). No vertical
ellipticity constant, heat-gradient square, or value-law density occurs.

This note does **not** prove the remaining geometric inverse. It replaces
its invalid CMC hypothesis by the formalizable hypothesis (0.6)--(0.8): a
nested family of constant-anisotropy perimeter minimizers with one bounded
spatial forcing. Extending the killed Wulff inverse from constant curvature
to this bounded-forcing class remains load bearing.

## 1. Exact subgraph identities

Let \(E\) be the affine support and, initially, let
\(d\mu=\rho(x)dx\) on an open convex set \(\Omega\subset E\). For
\(G:\Omega\to[0,1]\), put

\[
                  U_G=\{(x,r):0<r<G(x)\}.          \tag{1.1}
\]

More generally, let \(U\subset\Omega\times(0,1)\) have locally finite
perimeter and write its measure-theoretic normal as \(N=(N_x,N_r)\).
Weighted slicing gives

\[
 \begin{split}
 P_x(U)&:=\int |N_x|\,\rho\,d|D\mathbf1_U|
       =\int_0^1P_\mu(U_r)\,dr,\\
 M_x(U)&:=\int {N_xN_x^T\over|N_x|}\,\rho\,d|D\mathbf1_U|
       =\int_0^1M(U_r)\,dr,                         \tag{1.2}
 \end{split}
\]

where the matrix integrand is zero when \(N_x=0\). For a subgraph,

\[
             P_x(U_G)=\operatorname {TV}_\mu(G),
 \qquad M_x(U_G)=M(G).                              \tag{1.3}
\]

This includes vertical walls over a jump of \(G\): a wall of label height
\(G^+-G^-\) is counted for every level in that trace interval.

The vertical variation is

\[
 P_r(U)=\int|N_r|\rho\,d|D\mathbf1_U|
       =\int_\Omega \operatorname {TV}_{(0,1)}
                    (\mathbf1_U(x,\cdot))\,d\mu(x). \tag{1.4}
\]

If the bottom and top traces are respectively one and zero, then
\(P_r(U)\ge1\). Equality holds if and only if almost every vertical section
is an initial interval. Thus \(P_r(U)-1\) is an exact local measure of the
failure to be a subgraph. It is nevertheless the wrong regularizer for the
present purpose; Section 7 gives an exact cylinder counterexample.

## 2. The free-law functional

For \(G\in BV(\mu)\), write

\[
 DG=\sigma_G|DG|,\qquad
 M(G)=\int\sigma_G\sigma_G^T\,d|DG|_\mu.           \tag{2.1}
\]

The trace of \(M(G)\) is \(\operatorname {TV}_\mu(G)\). Fix a target
\(M_F=M(F)\). Nuclear/operator duality gives

\[
 \begin{split}
 K(G)
 &=\sup_{H=H^T,\ \|H\|_{op}\le1}
   \left\{\int\Phi_H(dDG)-\kappa\operatorname {tr}(HM_F)\right\},\\
 \Phi_H(q)&=|q|+\kappa{q^THq\over|q|},\qquad \Phi_H(0)=0.       \tag{2.2}
 \end{split}
\]

For \(h\perp n\),

\[
 D^2\Phi_H(n)[h,h]
 =|h|^2+\kappa\{2h^THh-(n^THn)|h|^2\}
 \ge(1-3\kappa)|h|^2.                              \tag{2.3}
\]

Hence every \(\Phi_H\) is convex and one-homogeneous, and \(K\) is a proper
convex lower-semicontinuous functional on \(L^1(\mu)/\mathbb R\).

### 2.1 The quotient norm identity

Let \(Y\) be an integrable real random variable and let \(c\) be a median.
Layer cake on the two sides of \(c\) gives

\[
 \begin{split}
 \mathbb E|Y-c|
 &=\int_{-\infty}^{c}\mathbb P(Y\le r)\,dr
   +\int_c^\infty\mathbb P(Y>r)\,dr\\
 &=\int_{\mathbb R}
       \min\{\mathbb P(Y>r),1-\mathbb P(Y>r)\}\,dr. \tag{2.4}
 \end{split}
\]

The formula remains true when the median is not unique; on the interval of
medians both expressions are affine with zero slope. Therefore

\[
 R(G)=\inf_c\|G-c\|_{L^1(\mu)}                     \tag{2.5}
\]

is exactly the quotient norm. In particular,

\[
                       |R(G)-R(U)|
 \le \inf_c\|G-U-c\|_1.                            \tag{2.6}
\]

Weighted coarea and the Cheeger inequality for finite-perimeter sets imply

\[
 \operatorname {TV}_\mu(G)
 =\int P_\mu(G>r)\,dr
 \ge pR(G).                                        \tag{2.7}
\]

Thus \(J=K-pR\ge0\). Notice that (2.7) uses the value of the Cheeger
constant; it does not assume a dimension-free lower bound for it.

### 2.2 The comparator is exactly the audited deficit

Let \(F=h_C(F_0)\) be the clipped central heat function. Its nontrivial
superlevels are precisely the heat levels with \(r\in C\); the remaining
levels are empty or full and contribute zero. Since the matrix term
vanishes at \(F\),

\[
 \begin{split}
 J(F)&=\operatorname {TV}_\mu(F)-pR(F)\\
 &=\int_C[P_\mu(F_0>r)
       -p\min\{\mu(F_0>r),1-\mu(F_0>r)\}]\,dr
 =:\mathcal D_{co}(F_0;C).                         \tag{2.8}
 \end{split}
\]

This point is sign-sensitive. The source definition of
\(\mathcal D_{co}\) is the deficit above \(p\min(v,1-v)\), not merely the
deficit above the isoperimetric profile. Hence the fixed-scale audit gives
directly

\[
             0\le J(F)\le\varepsilon p,             \tag{2.9}
\]

with \(\varepsilon=6.02\cdot10^{-14}\) in the retuned hierarchy.

## 3. Ekeland produces a bounded calibration

The following theorem avoids every compactness issue caused by value spikes
or escape to the tail.

**Theorem 3.1 (free-law Ekeland calibration).** Let \(p=\psi_\mu>0\),
let \(F\in BV(\mu)\), let \(0<\kappa<1/3\), and let \(\delta>0\). There is
\(G\in BV(\mu)\), defined up to an additive constant, such that

\[
 J(G)\le J(F),\qquad
 \|M(G)-M(F)\|_*\le {J(F)\over\kappa},              \tag{3.1}
\]

and there is

\[
 \xi\in\partial K(G)\subset(L^1(\mu)/\mathbb R)^*,\qquad
 \int\xi\,d\mu=0,\qquad
 \|\xi\|_\infty\le p+\delta.                       \tag{3.2}
\]

More precisely, one can choose a quotient-norm subgradient
\(q\in\partial R(G)\) and \(e\in L^\infty(\mu)\), both of mean zero, so that

\[
                  \xi=pq+e,\qquad
                  \|q\|_\infty\le1,\qquad
                  \|e\|_\infty\le\delta.           \tag{3.2a}
\]

After choosing a representative for which zero is a median,
\(q(x)=\operatorname {sign}G(x)\) wherever \(G(x)\ne0\); on
\(\{G=0\}\), \(q\) takes values in \([-1,1]\) so that its mean is zero.

**Proof.** The quotient \(X=L^1(\mu)/\mathbb R\) is a Banach space. The
functional \(J\) is lower semicontinuous on \(X\), because \(K\) is lower
semicontinuous and \(R\) is continuous. It is bounded below by zero. Apply
Ekeland's variational principle to \(F\), with variational slope \(\delta\).
If \(J(F)=0\), use \(G=F\); otherwise use Ekeland with
\(\epsilon=J(F)\) and \(\lambda=J(F)/\delta\). It gives

\[
 J(G)\le J(F),\qquad
 J(G)\le J(U)+\delta\|U-G\|_X\quad(U\in X).        \tag{3.3}
\]

Using (2.6) in (3.3),

\[
 \begin{split}
 K(G)
 &\le K(U)+p[R(G)-R(U)]+\delta\|U-G\|_X\\
 &\le K(U)+(p+\delta)\|U-G\|_X.                  \tag{3.4}
 \end{split}
\]

Thus \(G\) minimizes the convex functional
\(U\mapsto K(U)+(p+\delta)\|U-G\|_X\). The convex subgradient sum rule,
applicable because the norm is continuous everywhere, gives a
\(\xi\in\partial K(G)\) with dual norm at most \(p+\delta\). The dual of
\(L^1/\mathbb R\) is the annihilator of the constants in \(L^\infty\),
which proves (3.2).

For the refinement, (3.3) says that \(G\) is a global minimizer of

\[
 U\longmapsto K(U)+\delta\|U-G\|_X-pR(U).           \tag{3.4a}
\]

If \(A\) is proper convex lower semicontinuous, \(B\) is continuous convex,
and \(G\) is a local minimizer of \(A-B\), then
\(\partial B(G)\subset\partial A(G)\). Indeed, for
\(q\in\partial B(G)\), local minimality gives
\(A(U)\ge A(G)+B(U)-B(G)\ge A(G)+\langle q,U-G\rangle\) near
\(G\). Convexity of \(A\), applied on the segment from \(G\) to an arbitrary
\(U\), extends this affine support inequality globally. Apply this with
\(A=K+\delta\|\cdot-G\|_X\) and \(B=pR\). The convex sum rule gives

\[
 p\,\partial R(G)
 \subset \partial K(G)+\delta B_{X^*}.             \tag{3.4b}
\]

Choose any \(q\in\partial R(G)\). Then \(pq=\xi-e\) for some
\(\xi\in\partial K(G)\) and \(\|e\|_{X^*}\le\delta\); changing the sign of
\(e\) gives (3.2a). The stated pointwise description of \(q\) is the
standard equality case in the duality between the quotient \(L^1\) norm
and its \(L^\infty\) annihilator.

Finally \(J(G)\le J(F)\), (2.7), and (0.3) imply

\[
 \kappa\|M(G)-M(F)\|_*
 \le J(G)\le J(F).                                 \tag{3.5}
\]

QED.

### 3.1 One constant anisotropy and one vector field

We record the precise BV representation of (3.2). First suppose that
\(\Omega\) is bounded and Lipschitz and that
\(0<c\le\rho\le C\), with \(\rho\in C^1(\overline\Omega)\). Then there are
one matrix

\[
 H=H^T,\qquad \|H\|_{op}\le1,\qquad
 \operatorname {tr}[H(M(G)-M_F)]=\|M(G)-M_F\|_*,   \tag{3.6}
\]

and one \(z\in L^\infty(\mu;E)\) with weighted divergence in
\(L^\infty(\mu)\) such that

\[
 \begin{split}
 z(x)\cdot q&\le\Phi_H(q) &&(q\in E,\ \text{a.e. }x),\\
 (z,DG)_\mu&=\Phi_H(DG) &&\text{as Radon measures},\\
 -\operatorname {div}_\mu z&=\xi,\qquad
 z\cdot n_\Omega=0 &&\text{in the weak normal-trace sense}.     \tag{3.7}
 \end{split}
\]

Here \((z,DG)_\mu\) is the weighted Anzellotti pairing and
\(\operatorname {div}_\mu z=\rho^{-1}\operatorname {div}(\rho z)\).

To justify the fact that \(H\) is single and constant, use (2.2). The
active matrices are exactly the compact convex nuclear subgradient set in
(3.6). The max rule for convex subgradients writes \(\xi\) as a weak-star
limit of convex combinations of subgradients of the active anisotropic
variations. Standard weighted-BV duality represents each of those by a
field \(z_j\) satisfying

\[
 z_j\cdot q\le\Phi_{H_j}(q),\qquad
 (z_j,DG)_\mu=\Phi_{H_j}(DG).                       \tag{3.8}
\]

Put \(H=\sum a_jH_j\) and \(z=\sum a_jz_j\), and then take the weak-star
limit. Since the active set is convex, \(H\) is still active. Since
\(\Phi_H\) is affine in \(H\),

\[
 z\cdot q\le\sum a_j\Phi_{H_j}(q)=\Phi_H(q),\qquad
 (z,DG)_\mu=\sum a_j\Phi_{H_j}(DG)=\Phi_H(DG).      \tag{3.9}
\]

This proves (3.6)--(3.7). In particular, no level-dependent anisotropy is
introduced by the nuclear norm.

For an unbounded convex support, apply the same duality on an exhaustion
and use the uniform polar bound in (3.7) to take a weak-star limit. The
weighted divergence equation passes distributionally. For a hard convex
support, relative weighted variation is the correct perimeter; its dual
condition is precisely the zero normal trace in (3.7). For an
extended-valued convex potential, work on
\(\operatorname {ri}(\operatorname {dom}V)\) and approximate \(\rho\) from
inside. No boundary perimeter is added in this passage. These operations
do not change the constants \(1-3\kappa\) and \(p+\delta\).
Lower-dimensional support is handled on its affine hull.

The exhaustion statement uses only the standard dual representation of
weighted relative \(BV\). It does not assert \(C^2\) regularity of a hard
support. At a smooth contact point (3.7) yields the anisotropic Young
condition; a Wulff tube must still be killed at its first support contact.

### 3.2 Exact boxed saddle: no Ekeland error

There is a complementary implementation which is preferable for second
variation. Minimize \(J=K-pR\) on

\[
              \mathcal C=\{G\in BV(\mu):0\le G\le1\}.           \tag{3.10}
\]

A minimizer exists. Indeed, on a minimizing sequence,
\[
 \operatorname {TV}_\mu(G)
 =[\,\operatorname {TV}_\mu(G)-pR(G)\,]+pR(G)
 \le J(F)+1+p/2.                                   \tag{3.11}
\]
Local weighted \(BV\) compactness, the bound \(0\le G\le1\), and tightness
of \(\mu\) give global \(L^1\) compactness. The box is closed, \(K\) is
lower semicontinuous, and \(R\) is continuous.

Let \(G\) be a minimizer and choose any \(q\in\partial R(G)\). Convex
interpolation \(G_t=(1-t)G+tU\) gives, for \(U\in\mathcal C\),

\[
 \begin{split}
 K(G_t)-K(G)
 &\ge p[R(G_t)-R(G)]
 \ge pt\int q(U-G)\,d\mu,\\
 K(G_t)-K(G)&\le t[K(U)-K(G)].
                                                               \tag{3.12}
 \end{split}
\]

Hence the same \(G\) minimizes

\[
                         K(U)-p\int qU\,d\mu        \tag{3.13}
\]

over \(\mathcal C\). To select one compatible \(H\), write (3.13) as
\(\min_U\sup_H L(U,H)\), using (2.2). Restrict \(U\) to a total-variation
sublevel large enough to contain \(G\) and a minimizer of
\(L(\cdot,H)\) for every \(\|H\|_{op}\le1\). Such a uniform sublevel follows
by comparison with \(U=0\) and
\(\Phi_H\ge(1-\kappa)|\cdot|\). It is compact and convex in \(L^1\).
The compact convex minimax theorem then gives a saddle
\((G,H_*)\). Thus \(H_*\) is active for the nuclear norm and \(G\) minimizes

\[
       \operatorname {TV}_{\Phi_{H_*},\mu}(U)
                         -p\int qU\,d\mu.           \tag{3.14}
\]

Define the set energy
\[
 {\cal E}_q(B)=P_{\Phi_{H_*},\mu}(B)-p\int_Bq\,d\mu,
 \qquad e_*=\inf_B{\cal E}_q(B)\le0.               \tag{3.15}
\]
Coarea and layer cake make (3.14) the integral of
\({\cal E}_q(\{G>r\})\) over \(0<r<1\). Indicators show that the infimum of
(3.14) is exactly \(e_*\). Therefore

\[
 \boxed{{\cal E}_q(\{G>r\})=e_*
              \quad\text{for almost every }0<r<1.} \tag{3.16}
\]

This exact setwise theorem does not require a global divergence
calibration. If one does introduce a calibration for the boxed problem,
normal-cone measures on \(\{G=0\}\) and \(\{G=1\}\) are indispensable; in
general they carry the negative constant \(e_*\). Thus (3.7) is the
obstacle-free statement for the quotient/Ekeland construction and must not
be copied verbatim to (3.10).

Every level in (3.16) satisfies

\[
 \begin{split}
 P_{\Phi_{H_*},\mu}(A)
 &\le P_{\Phi_{H_*},\mu}(B)+p\mu(A\triangle B),\\
 (1-\kappa)p\min(\mu A,1-\mu A)
 &\le P_{\Phi_{H_*},\mu}(A)
 \le p\min(\mu A,1-\mu A).                         \tag{3.17}
 \end{split}
\]

The first line follows by comparing the common minimum \(e_*\) and using
\(|q|\le1\). For the second, use \(e_*\le0\),
\(\int_Aq\le\min(\mu A,1-\mu A)\), anisotropic ellipticity, and the Cheeger
lower bound. In particular, every boxed leaf has generalized anisotropic
mean curvature bounded by \(p\), with no Ekeland error.

## 4. Coarea removes the jump and Cantor obstruction

Let \(A_r=\{G>r\}\). Vector-valued coarea and anisotropic coarea give

\[
 DG=\int_{\mathbb R}D\mathbf1_{A_r}\,dr,\qquad
 \Phi_H(DG)=\int_{\mathbb R}\Phi_H(D\mathbf1_{A_r})\,dr.        \tag{4.1}
\]

The polar inequality in (3.7) implies that

\[
 \Phi_H(D\mathbf1_{A_r})-(z,D\mathbf1_{A_r})_\mu
\]

is a nonnegative measure for almost every \(r\). Its integral in \(r\)
is zero by saturation on \(DG\). Hence, for almost every \(r\),

\[
 (z,D\mathbf1_{A_r})_\mu
             =\Phi_H(D\mathbf1_{A_r}).             \tag{4.2}
\]

Gauss--Green, (3.7), and the zero support trace now give

\[
 P_{\Phi_H,\mu}(A_r)
 =\int_{A_r}\xi\,d\mu.                             \tag{4.3}
\]

For every finite-perimeter \(B\), the polar inequality gives

\[
 P_{\Phi_H,\mu}(B)\ge\int_B\xi\,d\mu.             \tag{4.4}
\]

Equations (4.3)--(4.4) prove (0.6)--(0.7). Since
\(\int\xi d\mu=0\),

\[
 \int_{A_r}\xi\,d\mu
 \le\|\xi\|_\infty\min\{\mu(A_r),1-\mu(A_r)\}.     \tag{4.5}
\]

Combining (4.5), \(P_{\Phi_H}\ge(1-\kappa)P_\mu\), and the definition of
\(p\) proves (0.8).

This proof treats the Cantor part exactly. It is present in \(DG\), but
the coarea identity decomposes both its scalar and matrix charges into the
reduced boundaries of the sets \(A_r\). No claim of the form
\(|D^cG|=0\) is needed.

### 4.1 Exact finite-phase jump law

The preceding calibration is the robust statement. In a smooth
finite-phase model one can also see the multiplier explicitly. Let

\[
 G=\sum_{i=0}^m a_i\mathbf1_{E_i},\qquad
 a_0<a_1<\cdots<a_m,                               \tag{4.6}
\]

and put \(v(r)=\mu(G>r)\) and
\(m_0(v)=\min(v,1-v)\). Suppose first that \(G\) is an exact local
minimizer and that a compatible \(H\) has been selected. On a smooth
interface between phases \(a_i\) and \(a_j\), \(i<j\), the shape derivative
is

\[
 \boxed{
 {\mathcal H}_{ij}
 ={p\over a_j-a_i}\int_{a_i}^{a_j}\zeta(r)\,dr,
 \qquad
 \zeta(r)\in\partial m_0(v(r))\subset[-1,1].}       \tag{4.7}
\]

Indeed, moving the high phase across that interface changes \(v(r)\) by
the same swept volume for every \(r\in(a_i,a_j)\). The anisotropic
perimeter has tension \(a_j-a_i\), while the first variation of \(pR\) is
\(p\int_{a_i}^{a_j}\zeta(r)dr\). This proves (4.7), with the sign fixed by
the chosen normal orientation. In particular,

\[
                         |{\mathcal H}_{ij}|\le p. \tag{4.8}
\]

If the trace interval does not cross a median level, then \(\zeta\) is
identically \(1\) or \(-1\), and the curvature is exactly \(p\) or \(-p\).
If it crosses a median, (4.7) is the correct divided difference and remains
bounded by \(p\).

For the Ekeland point with slope \(\delta\), the same one-interface
variation has quotient-\(L^1\) size
\((a_j-a_i)\int|u|dP\) to first order. Hence (4.7) has an additive error
\(e_{ij}\) with \(|e_{ij}|\le\delta\), and

\[
                         |{\mathcal H}_{ij}|\le p+\delta.       \tag{4.9}
\]

This is the precise replacement for the uncontrolled divided-difference
law in the earlier BV stationarity audit.

### 4.2 The bounded-forcing killed-tube identity

Assume a regular patch of a calibrated leaf is \(C^2\), and assume \(V\)
is \(C^2\) there. Put

\[
 z_H=D\Phi_H(n),\qquad
 B_H=D^2\Phi_H(n)|_{n^\perp}.                       \tag{4.10}
\]

Along the anisotropic ray \(x+t z_H(x)\), before focal time, collision, or
support contact, the weighted flux Jacobian satisfies

\[
 {d^2\over dt^2}\log j_x(t)
 =-\operatorname {tr}
 [((I+tB_HS)^{-1}B_HS)^2]
  -\nabla^2V(x+t z_H)[z_H,z_H]\le0.                \tag{4.11}
\]

Its initial slope is the weighted anisotropic mean curvature \(h(x)\), and
the calibration gives \(|h(x)|\le p+\delta\). Therefore the exact
pointwise formula is

\[
                  j_x(t)=\exp\{h(x)t-D_x(t)\},
 \qquad D_x(t)\ge0.                                \tag{4.12}
\]

No common value of \(h\) is needed for (4.12). On a tube of length
\(|t|\le c/p\), its variable linear factor is bounded by

\[
                     e^{|h(x)t|}
 \le \exp\{c(1+\delta/p)\}.                        \tag{4.13}
\]

Thus the contact correction is dimension-free. What remains to be proved
for the full route is a singular killed-tube/inverse theorem using the
bounded-forcing variational formulation (0.6), rather than a smooth CMC
assertion.

### 4.3 Conditional flatness from unconstrained stability

There is a strong consequence on any leaf for which the median cusp is
inactive and exact unconstrained stability is available.

**Proposition 4.1 (stable one-sided leaves are flat).** Assume
\(\rho=e^{-V}\) is \(C^2\), \(\Phi\) is a constant \(C^3\), uniformly
elliptic one-homogeneous anisotropy, and a regular leaf \(\Sigma\) is a
two-sided stable critical point of

\[
                         P_{\Phi,\mu}(A)-s p\,\mu(A),
 \qquad s\in\{-1,1\},                              \tag{4.14}
\]

under all compactly supported variations. Allow a smooth convex hard wall
with the natural anisotropic Young condition. If
\(P_{\Phi,\mu}(\Sigma)<\infty\), then, almost everywhere on every regular
component,

\[
 S=0,\qquad
 \nabla^2V[D\Phi(n),D\Phi(n)]=0.                   \tag{4.15}
\]

At a smooth hard-wall contact, the corresponding nonnegative wall second
fundamental-form term also vanishes.

**Proof.** Put \(z=D\Phi(n)\) and
\(B=D^2\Phi(n)|_{n^\perp}\). The anisotropic Jacobi formula, written for a
Wulff-normal variation \(X=u z\), has the form

\[
 Q(u)=\int_\Sigma
 \left\{\langle C_\Phi(n)\nabla_\Sigma u,\nabla_\Sigma u\rangle
 -q_\Phi u^2\right\}\,d\sigma_{\Phi,\mu}
 -\int_{\partial\Sigma}b_{\partial\Omega,\Phi}u^2\,d\tau,       \tag{4.16}
\]

where

\[
 q_\Phi=\operatorname {tr}[(BS)^2]
       +\nabla^2V[z,z]\ge0,\qquad
 0\preceq C_\Phi(n)\preceq C(\kappa)I,\qquad
 b_{\partial\Omega,\Phi}\ge0.                      \tag{4.17}
\]

Formula (4.16) follows by differentiating the Wulff-normal graph
\(x\mapsto x+t u(x)z(x)\). For \(u\equiv1\), its potential is exactly
minus the second logarithmic derivative in (4.11); the terms containing
\(\nabla u\) form the positive tangent Hessian of the parametric
integrand. Convexity of \(V\), ellipticity of \(\Phi\), and convexity of the
hard wall give the three signs in (4.17). This is also the weighted
anisotropic free-boundary Jacobi formula; the last term is absent without
a hard wall.

Choose \(u_R(x)=\chi(|x|/R)\), where \(\chi=1\) on \([0,1]\),
\(\chi=0\) on \([2,\infty)\), and \(|\chi'|\le2\). Then
\[
 \int_\Sigma\langle C_\Phi\nabla u_R,\nabla u_R\rangle
      d\sigma_{\Phi,\mu}
 \le {4C(\kappa)\over R^2}P_{\Phi,\mu}(\Sigma)\longrightarrow0. \tag{4.18}
\]
The natural contact variation is admissible; the convex-wall term in
(4.16) only strengthens the inequality. Stability, Fatou's lemma, and
(4.18) imply \(\int_\Sigma q_\Phi\,d\sigma_{\Phi,\mu}=0\).
Both summands in (4.17) are nonnegative. Moreover \(BS\) is similar to the
symmetric matrix \(B^{1/2}SB^{1/2}\), so
\(\operatorname {tr}[(BS)^2]=0\) and ellipticity imply \(S=0\). This proves
(4.15). QED.

For the exact boxed construction of Section 3.2, Proposition 4.1 applies on
every neighborhood of a regular leaf on which \(q\) is the constant sign
\(s\), because (3.16) then makes the leaf a genuine local minimizer of
(4.14). This includes a diffuse regular level strictly above or below
every median value. It also includes a finite-phase interface when both
traces lie strictly on the same side of the chosen median. The more general
statement that a trace interval merely stays on one side of the
*half-mass profile* follows from the explicit phase law (4.7) only when the
median-plateau subgradient can be selected locally constant; this extra
selection must not be assumed at a contact.

There are two reasons not to apply Proposition 4.1 silently in every
branch.

1. If one uses the obstacle-free Ekeland construction rather than the exact
   boxed saddle, the Ekeland error is first order:
   \(J(G)\le J(G_t)+\delta\,O(|t|)\). Literal second variation at
   \(t\to0\) therefore has an error \(\delta/|t|\). One may take
   \(t\gg\delta\) along a sequence \(\delta\downarrow0\), but using this
   requires a uniform finite-difference Jacobi remainder and compactness of
   the varying leaves. Neither follows from (3.3).
2. If a jump trace interval contains a half-mass value, the function
   \(m_0(v)=\min(v,1-v)\) has a cusp. The interface only satisfies
   \(|\mathcal H|\le p\); it need not be an unconstrained stable critical
   point of either sign in (4.14).

The second residual can carry all of the matrix. Indeed, if \(A\) is a
balanced Cheeger minimizer, \(\mu(A)=1/2\) and \(P_\mu(A)=p/2\), then
\(G=\mathbf1_A\) satisfies

\[
 \operatorname {TV}(G)-pR(G)=0.                    \tag{4.19}
\]

If its normal matrix is the target, then \(J(G)=0\), and the entire matrix
lies on the one median-cusp jump. Thus the free-law functional by itself
does not prove flatness of the balanced extremal boundary; excluding or
classifying that branch is still conjecture-strength geometry.

It is useful to state this residual as a matrix decomposition. Fix a median
\(c\), and on the jump set put

\[
 M_{\rm cross}^j
 =\int_{J_G:\ G^-<c<G^+}
   (G^+-G^-)\,n_Gn_G^T\rho\,d\mathcal H^{k-1}.      \tag{4.20}
\]

The complementary jump matrix comes from interfaces whose two traces are
strictly on one side of \(c\), together with interfaces having one trace
equal to a median plateau. The latter are one-sided only when the
subgradient on that plateau can be selected locally constant; otherwise
they belong to the residual branch as well. The diffuse coarea charge at
levels \(r\ne c\) is one-sided. In the exact, regular, stable setting,
Proposition 4.1 makes the genuinely one-sided pieces flat/log-affine.
There is, however, no estimate

\[
                  \operatorname {tr}M_{\rm cross}^j
                 \le C\,J(F)                       \tag{4.21}
\]

with a finite universal \(C\): the balanced binary example has right-hand
side zero and can have a nonzero cross matrix. Any final inverse must
therefore accept the median-crossing matrix as a genuine branch, rather
than charge it to the scalar coarea deficit.

## 5. Matrix retention and the explicit constants

Put

\[
 M=M(F),\quad M'=M(G),\quad T=\operatorname {tr}M,
 \quad T'=\operatorname {tr}M',\quad Q=M/T,\quad Q'=M'/T'.       \tag{5.1}
\]

Let \(D=J(F)\). Theorem 3.1 gives

\[
 d:=\|M'-M\|_*\le D/\kappa,\qquad |T'-T|\le d.     \tag{5.2}
\]

If \(d<T\), then

\[
 \begin{split}
 \|Q'-Q\|_*
 &\le{\|M'-M\|_*\over T}
       +\left|{1\over T'}-{1\over T}\right|\|M'\|_*\\
 &\le {2d\over T}.                                 \tag{5.3}
 \end{split}
\]

Both \(Q,Q'\) are positive contractions with trace one. If \(E=Q'-Q\),
then \(\operatorname {tr}E=0\) and

\[
 \operatorname {tr}(Q'^2-Q^2)
 =\operatorname {tr}[E(Q'+Q-I)],\qquad
 \|Q'+Q-I\|_{op}\le1.                             \tag{5.4}
\]

Consequently

\[
 |\operatorname {tr}Q'^2-\operatorname {tr}Q^2|
 \le {2D\over\kappa T}.                            \tag{5.5}
\]

At the retuned fixed scale,

\[
 {D\over p}\le6.02\cdot10^{-14},\qquad
 {T\over p}>.0032827,\qquad
 1-\operatorname {tr}Q^2>.0032688.                 \tag{5.6}
\]

With \(\kappa=10^{-6}\),

\[
 {D\over\kappa T}<1.834\cdot10^{-5},\qquad
 {2D\over\kappa T}<3.668\cdot10^{-5}.             \tag{5.7}
\]

Thus

\[
                    1-\operatorname {tr}Q'^2>.0032321.          \tag{5.8}
\]

All of these numbers are independent of dimension. One may take, for
example, \(\delta=10^{-8}p\); the leafwise curvature and profile loss are
then \(1+O(10^{-6})\), dominated by the retained angular variance.

## 6. A finite layered incidence regularizer

The free-law construction retains actual nesting because its leaves are the
superlevels of one \(G\). For comparison, there is also a useful finite
layer regularizer which quantifies the price of abandoning nesting.

Let \(A_1^0\supset\cdots\supset A_N^0\) be reference heat levels with
weights \(w_i>0\), \(W=\sum_iw_i\le1\), and volumes \(v_i\). Put

\[
 M_0=\sum_iw_iM(A_i^0),\qquad
 D_N=\sum_iw_i[P(A_i^0)-I_\mu(v_i)],                \tag{6.1}
\]

and define the inversion charge

\[
             \mathcal C(A_1,\ldots,A_N)
 =\sum_{i<j}w_iw_j\mu(A_j\setminus A_i).           \tag{6.2}
\]

Minimize

\[
 \sum_iw_iP(A_i)
 +\kappa\left\|\sum_iw_iM(A_i)-M_0\right\|_*
 +\eta\mathcal C(A_1,\ldots,A_N)                  \tag{6.3}
\]

subject to \(\mu(A_i)=v_i\). In a bounded smooth weighted support a
minimizer exists by finite-product \(BV\) compactness. Comparison with the
nested reference and the isoperimetric lower bound give the exact budget

\[
 \boxed{
 \kappa\left\|\sum_iw_iM(A_i)-M_0\right\|_*
 +\eta\mathcal C(A_1,\ldots,A_N)\le D_N.}           \tag{6.4}
\]

The overlap term is Lipschitz under a change of one leaf:

\[
 |\mathcal C(\ldots,A_i,\ldots)
  -\mathcal C(\ldots,B_i,\ldots)|
 \le w_iW\mu(A_i\triangle B_i).                   \tag{6.5}
\]

The same compatible-subgradient separation as in Section 3 therefore gives
one constant \(H\) and, on every regular leaf, a volume multiplier
\(\lambda_i\) and a contact force \(h_i\) such that

\[
                  \mathcal H_{\Phi_H,\mu}(A_i)
                 =\lambda_i+h_i,\qquad
                  \|h_i\|_\infty\le\eta W.         \tag{6.6}
\]

At coincident interfaces the derivative of (6.2) is one-sided, but its
Clarke subgradient still has norm bounded by (6.5); this is exactly why the
bound in (6.6) survives contact.

For equal weights \(w_i=1/N\), sort each binary column
\((\mathbf1_{A_1}(x),\ldots,\mathbf1_{A_N}(x))\) into nonincreasing order,
and call the resulting nested tuple \(B_i\). If \(k(x)\) zeros occur before
the sorted cut, then exactly \(k(x)\) ones occur after it. All \(k(x)^2\)
cross pairs are inversions. Hence, pointwise,

\[
 {1\over N}\sum_i|\mathbf1_{A_i}-\mathbf1_{B_i}|
 ={2k\over N}\le2\sqrt{\mathcal C_x}.             \tag{6.7}
\]

After integration,

\[
 {1\over N}\sum_i\mu(A_i\triangle B_i)
 \le2\sqrt{\mathcal C}.                            \tag{6.8}
\]

If \(D_N\le\varepsilon p\) and one chooses
\(\eta=p\sqrt\varepsilon\), then (6.4), (6.6), and (6.8) give simultaneously

\[
 \mathcal C\le\sqrt\varepsilon,\qquad
 \|h_i\|_\infty\le p\sqrt\varepsilon,\qquad
 {1\over N}\sum_i\mu(A_i\triangle B_i)
 \le2\varepsilon^{1/4}.                            \tag{6.9}
\]

This finite theorem is a fallback if one wants explicit layer variables.
The free-law calibration is stronger: it keeps exact nesting and replaces
the contact force by one bounded spatial forcing without paying (6.9).

## 7. Why the obvious elliptic regularizers do not suffice

### 7.1 Full subgraph area does not remove jumps

Adding \(\epsilon\) times full subgraph perimeter gives, for a smooth graph,

\[
 \epsilon\int\sqrt{1+|\nabla G|^2}\,d\mu,           \tag{7.1}
\]

and its \(BV\) relaxation has recession function
\(\epsilon|D^sG|\). Thus a jump wall still has finite energy and is not
excluded. The heat comparator cost is bounded using only first variation,

\[
 \epsilon\int\sqrt{1+|\nabla F|^2}\,d\mu
 \le\epsilon[1+\operatorname {TV}_\mu(F)],          \tag{7.2}
\]

so fitting it into an \(o(p)\) matrix budget forces
\(\epsilon=o(p)\). Its ellipticity then degenerates with the very Cheeger
constant being estimated.

More seriously, small energy does not control its Euler error. On the unit
circle, let

\[
 \epsilon_m=m^{-1/2},\qquad
 u_m(x)=\tfrac12+m^{-3/2}\sin(mx).                  \tag{7.3}
\]

Then \(\|u_m'\|_\infty=m^{-1/2}\), and

\[
 \epsilon_m\int(\sqrt{1+|u_m'|^2}-1)dx=O(m^{-3/2})\to0,        \tag{7.4}
\]

whereas the Euler term

\[
 \epsilon_m{d\over dx}
       {u_m'\over\sqrt{1+|u_m'|^2}}                \tag{7.5}
\]

has amplitude tending to one. Thus no estimate of pointwise or \(L^1\)
curvature error follows from the comparator energy.

The same example defeats a quadratic regularizer. The cost
\(\epsilon_m\int|u_m'|^2=O(m^{-3/2})\), while
\(\epsilon_m u_m''\) has order-one amplitude. Superlinear growth removes
the literal jump only by introducing an uncontrolled second-derivative
term.

### 7.2 Vertical total variation has an unbounded contact multiplier

Consider in \(E\times(0,1)\) the anisotropy

\[
                     \Psi(N)=\Phi(N_x)+\epsilon|N_r|.           \tag{7.6}
\]

On a vertical cylinder \(\Sigma\times(a,b)\), \(N_r=0\). A Cahn--Hoffman
field may have vertical component any
\(\zeta(r)\in[-\epsilon,\epsilon]\). Its anisotropic mean curvature is

\[
                   \mathcal H_\Phi(\Sigma)+\zeta'(r).           \tag{7.7}
\]

Taking \(\zeta(r)=\epsilon\sin(Nr)\) makes the pressure oscillation have
amplitude \(\epsilon N\), with no change at all in the surface energy.
Thus a vertical-\(BV\) penalty gives only a negative-Sobolev control in the
label variable; it gives no pointwise-in-level pressure bound. This is the
same divided-difference obstruction in product notation.

The free-law functional avoids both failures. Its lower-order term is the
quotient \(L^1\) norm, whose entire subgradient lies in the unit ball of
\(L^\infty\). The forcing therefore has the exact pointwise bound (0.5),
with no derivative in either the physical or label variable.

## 8. Exact conclusion

The following pieces are now dimension-free and do not assume KLS.

1. The heat coarea deficit is the comparator value of the free-law
   functional, exactly and with the audited constant.
2. Ekeland supplies an obstacle-free almost minimizer without a compactness
   or value-range assumption; alternatively, the exact boxed saddle gives
   exact levelwise stability with explicit obstacle bookkeeping.
3. Convex duality supplies one constant matrix anisotropy; in the Ekeland
   version it also supplies one bounded obstacle-free calibration.
4. Coarea transfers that calibration to almost every level, including all
   jump and Cantor charge.
5. Every calibrated level is nested, globally quasiminimal, within
   \(1+O(\kappa+\delta/p)\) of the Cheeger profile, and has generalized
   anisotropic curvature bounded by \(p+\delta\).
6. The retained normal-projector variance remains above \(0.0032321\) at
   the retuned constants.

The remaining statement is geometric rather than variational: prove the
high-rank killed-Wulff inverse for a nested family satisfying (0.6)--(0.8),
with bounded measurable forcing in place of constant curvature. The
regularization itself no longer has a BV jump/contact gap.
