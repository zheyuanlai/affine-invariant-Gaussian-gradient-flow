# Clean-room audit of spatial phase coherence

## Verdict

Fix \(0<\delta\leq 1/2\), and write
\[
I(g)=\varphi(\Phi^{-1}(g))
\]
for the Gaussian isoperimetric profile.  Subject to the two conventions stated
below, **(L1)--(L4) are all true**.  In fact, (L1) follows from an exact
one-dimensional deficit decomposition, and the Jacobian part of (L3) follows
from (L2) with no loss.  The factor \(2\) in the product approximation is the
natural one when \(R\) is required to take values in \(u^\perp\); if that range
condition is omitted, a conditional-expectation choice even gives factor
\(1\).

The conventions are important.

1. Since \(u=v/|v|\), the statement is literally defined only when \(v\ne0\).
   If \(v=0\), one may choose any \(u\); then \(\varepsilon=1\), and the
   estimates become trivial after increasing \(C_\delta\).
2. By the **standardized Brenier contraction** I mean the following.  If
   \(T_0\) is the Brenier map from \(\gamma_n=N(0,I_n)\) to \(\pi\), then
   \[
   T(z)=\sqrt t\bigl(T_0(z)-\mathbb E_\pi X\bigr).
   \]
   Thus \(T_\#\gamma_n\) is the law of
   \(\sqrt t(X-\mathbb E X)\), and Caffarelli's theorem gives
   \(0\preceq H_T\preceq I\) almost everywhere.

Without this standardization, (L3) is false.  For example, let
\(\pi=N(0,t^{-1}I_n)\) and let \(S\) be a halfspace.  Then
\(\varepsilon=0\), while the unstandardized map has
\(H_{T_0}=t^{-1/2}I\), so
\(\mathbb E|H_{T_0}u-u|^2=(t^{-1/2}-1)^2\), which is nonzero if \(t\ne1\).

No counterexample exists for the four claims under the stated conventions.

This audit stops exactly at (L1)--(L4).  It does **not** audit any subsequent
claim of an exact full-measure product splitting, nor any equality or
compatibility conclusion obtained by applying these estimates in two distinct
tilt directions.  Those are separate rigidity questions.

---

## 1. Reduction to a scalar Gaussian quantile problem

Let \(m=\mathbb E_\pi X\), and set
\[
Y=\sqrt t\,\langle X-m,u\rangle .
\]
The law \(\eta\) of \(Y\) is one-strongly log-concave and has mean zero.
Indeed, after writing the density of \(\pi\) in the extended-valued form
\[
e^{-t|x|^2/2-\Psi(x)},\qquad \Psi:\mathbb R^n\to(-\infty,+\infty]
\quad\hbox{convex},
\]
Pr\'ekopa's theorem shows that the marginal density in direction \(u\), after
the factor \(\sqrt t\), is \(e^{-y^2/2-\psi(y)}\) with \(\psi\) convex.

Let \(Z\sim\gamma_1\), and let
\[
q(z)=F_\eta^{-1}(\Phi(z))
\]
be the increasing rearrangement from \(\gamma_1\) to \(\eta\).  The
one-dimensional Caffarelli contraction theorem gives
\[
q\text{ is absolutely continuous},\qquad 0\le q'(z)\le1
\quad\text{for a.e. }z,\qquad \mathbb E q(Z)=0.                \tag{1.1}
\]
This remains valid when the support is a proper interval; the approximation
needed for that assertion is given in Section 6.

Put
\[
c=\Phi^{-1}(1-g),\qquad h(z)=\mathbf 1_{\{z\ge c\}}.
\]
Then \(\mathbb Eh=g\) and
\[
\mathbb E[Zh(Z)]=\varphi(c)=I(g).                              \tag{1.2}
\]
Since a one-strongly log-concave marginal has a continuous density on the
interior of its interval support, we may realize
\(Y=q(Z)\) and disintegrate the event \(S\) over \(Z\).  Namely, define
\[
p(z)=\mathbb P(S\mid Z=z),\qquad 0\le p\le1.
\]
Then the exact bookkeeping is
\[
\mathbb Ep(Z)=g,                                                \tag{1.3}
\]
\[
\mathbb E[p(Z)q(Z)]
=\mathbb E[\mathbf1_S Y]
=\sqrt t\,\langle v,u\rangle
=\sqrt t\,|v|
=(1-\varepsilon)I(g),                                          \tag{1.4}
\]
and the halfspace \(H\) of mass \(g\) is, up to a null set,
\(\{Y\ge q(c)\}=\{Z\ge c\}\).  Consequently
\[
\pi(S\mathbin{\triangle} H)
=\mathbb E\,|p(Z)-h(Z)|.                                       \tag{1.5}
\]
The use of a conditional probability in (1.5) is essential: no assumption
that \(S\) is measurable with respect to the scalar projection is being made.

---

## 2. The exact deficit decomposition

Define
\[
\Delta=I(g)-\mathbb E[p(Z)q(Z)]=\varepsilon I(g),               \tag{2.1}
\]
and split it as
\[
\Delta=D_q+D_p,                                                 \tag{2.2}
\]
where
\[
D_q=I(g)-\mathbb E[h(Z)q(Z)],\qquad
D_p=\mathbb E[(h(Z)-p(Z))q(Z)].                                \tag{2.3}
\]
Both terms are nonnegative.  The second is nonnegative because \(q\) is
increasing and the upper tail \(h\) maximizes the integral of \(q\) among all
\([0,1]\)-valued functions of mean \(g\).

For the first term, let
\[
e(s)=1-q'(s)\in[0,1]
\]
and define the positive kernel
\[
K_g(s)=
\begin{cases}
g\Phi(s),&s<c,\\
(1-g)\overline\Phi(s),&s\ge c.
\end{cases}                                                     \tag{2.4}
\]
Since \(d(z)=z-q(z)\) has \(\mathbb Ed(Z)=0\) and \(d'=e\), the
one-dimensional covariance formula gives
\[
D_q=\mathbb E[h(Z)d(Z)]
=\int_{\mathbb R}e(s)K_g(s)\,ds.                               \tag{2.5}
\]
For completeness, the kernel in (2.5) is obtained from
\[
\operatorname{Cov}(\mathbf1_{\{Z\ge s\}},h(Z))
=
\begin{cases}
g\Phi(s),&s<c,\\
(1-g)\overline\Phi(s),&s\ge c.
\end{cases}
\]
Integrating this identity against \(d'(s)\,ds\) proves (2.5).  In particular,
\(D_q\ge0\), and (2.2) also proves the centroid bound
\(\sqrt t|v|\le I(g)\), hence \(0\le\varepsilon\le1\).

---

## 3. Proof of (L1): stability of the set

Let
\[
r=\mathbb E|p-h|=\pi(S\mathbin{\triangle} H).
\]
The mass moved across the threshold is the same on the two sides.  More
precisely,
\[
\alpha:=\mathbb E[p(Z)\mathbf1_{\{Z<c\}}]
=\mathbb E[(1-p(Z))\mathbf1_{\{Z\ge c\}}]
=\frac r2.                                                      \tag{3.1}
\]
If \(\alpha=0\), there is nothing to prove.  Otherwise let \(Z_-\) and
\(Z_+\) be independent variables with probability laws
\[
\frac{p(z)\mathbf1_{\{z<c\}}}{\alpha}\,d\gamma_1(z),
\qquad
\frac{(1-p(z))\mathbf1_{\{z\ge c\}}}{\alpha}\,d\gamma_1(z),   \tag{3.2}
\]
respectively.  Thus \(Z_-<c\le Z_+\), and
\[
D_p=\alpha\,\mathbb E[q(Z_+)-q(Z_-)].                           \tag{3.3}
\]
Introduce the corresponding Gaussian rearrangement cost
\[
D_G:=\mathbb E[(h-p)Z]
=\alpha\,\mathbb E[Z_+-Z_-].                                  \tag{3.4}
\]
Since \(z-q(z)=d(z)\) and \(d'=e\), Tonelli's theorem gives the exact identity
\[
D_G-D_p
=\alpha\,\mathbb E\int_{Z_-}^{Z_+}e(s)\,ds
=\int_{\mathbb R}e(s)W(s)\,ds,                                \tag{3.5}
\]
where
\[
W(s)=\alpha\,\mathbb P(Z_-\le s\le Z_+).
\]
The exchanged-mass laws in (3.2) give, with no hidden endpoint terms,
\[
W(s)=\int_{-\infty}^s p\,d\gamma_1\le\Phi(s)
\quad(s<c),                                                     \tag{3.6}
\]
and
\[
W(s)=\int_s^\infty(1-p)\,d\gamma_1\le\overline\Phi(s)
\quad(s\ge c).                                                 \tag{3.7}
\]
Because \(g,1-g\ge\delta\), (2.4), (3.6), and (3.7) imply
\[
W(s)\le \delta^{-1}K_g(s).
\]
Together with (2.5), this yields
\[
D_G\le D_p+\delta^{-1}D_q
\le\delta^{-1}(D_p+D_q)
=\delta^{-1}\Delta.                                           \tag{3.8}
\]

It remains only to quantify the Gaussian exchange cost.  Write
\(s_0=\Phi(c)=1-g\) and \(Q=\Phi^{-1}\).  The bathtub principle says that,
among all lower and upper exchanged pieces of Gaussian mass \(\alpha\), the
cost in (3.4) is minimized by taking the two intervals immediately adjacent
to \(c\).  Hence
\[
\begin{aligned}
D_G\ge{}&\int_{s_0-\alpha}^{s_0}[Q(s_0)-Q(s)]\,ds\\
&+\int_{s_0}^{s_0+\alpha}[Q(s)-Q(s_0)]\,ds.                    \tag{3.9}
\end{aligned}
\]
The endpoints are legitimate because
\(\alpha\le\min(g,1-g)\).  Moreover
\[
Q'(s)=\frac1{\varphi(Q(s))}\ge\sqrt{2\pi}.
\]
Therefore each integral in (3.9) is at least
\(\sqrt{2\pi}\,\alpha^2/2\), so
\[
D_G\ge\sqrt{2\pi}\,\alpha^2
=\frac{\sqrt{2\pi}}4r^2.                                     \tag{3.10}
\]
Combining (2.1), (3.8), and (3.10), and using
\(I(g)\le(2\pi)^{-1/2}\), gives the explicit estimate
\[
\boxed{\quad
\pi(S\mathbin{\triangle} H)
\le \sqrt{\frac{2}{\pi\delta}}\,\sqrt\varepsilon .
\quad}                                                         \tag{3.11}
\]
This proves (L1).

---

## 4. Proof of (L2): variance rigidity

Only the contraction part \(D_q\) is needed.  Let
\[
B=\mathbb E e(Z)=\int e(s)\varphi(s)\,ds.                       \tag{4.1}
\]
The elementary inverse Mills bounds imply, for a universal \(C\),
\[
\frac{\varphi(s)}{\Phi(s)}\le C(1+|s|),\qquad
\frac{\varphi(s)}{\overline\Phi(s)}\le C(1+|s|).              \tag{4.2}
\]
Using (2.4) and \(g,1-g\ge\delta\), for any \(R\ge1\),
\[
\begin{aligned}
B
&\le \int_{|s|\le R}e(s)\varphi(s)\,ds+\gamma_1(|Z|>R)\\
&\le \frac{C(1+R)}\delta D_q+\gamma_1(|Z|>R).                 \tag{4.3}
\end{aligned}
\]
If \(D_q=0\), then (2.5) and positivity of \(K_g\) imply \(e=0\) almost
everywhere, so \(B=0\).  If \(0<D_q\le1\), choose
\[
R=\sqrt{2\log(e/D_q)}.
\]
The Gaussian tail bound
\(\gamma_1(|Z|>R)\le 2\varphi(R)/R\) then gives
\[
B\le C_\delta D_q\sqrt{\log(e/D_q)}.                           \tag{4.4}
\]

Now \(d=\operatorname{id}-q\) is Lipschitz, \(\mathbb Ed=0\), and Gaussian
integration by parts gives
\[
\mathbb E[Zd(Z)]=\mathbb Ed'(Z)=B.
\]
It follows that
\[
1-\operatorname{Var}(q(Z))
=2B-\mathbb E d(Z)^2
\le2B.                                                         \tag{4.5}
\]
The function \(x\mapsto x\sqrt{\log(e/x)}\) is increasing on \([0,1]\).
Since
\[
D_q\le\Delta=\varepsilon I(g)\le\varepsilon\le1,
\]
(4.4)--(4.5) show
\[
1-\operatorname{Var}(Y)
\le C_\delta\varepsilon\sqrt{\log(e/\varepsilon)}.            \tag{4.6}
\]
Finally, \(\operatorname{Var}(Y)=t\operatorname{Var}\langle X,u\rangle\),
so (4.6) is exactly (L2).  At \(\varepsilon=0\), the right-hand side is
understood by continuity and equals zero.

---

## 5. Proof of (L3): Jacobian and product structure

Let \(G\sim\gamma_n\), let \(T\) be the standardized Brenier contraction, and
write
\[
f(G)=\langle T(G),u\rangle.
\]
This has the same law as \(Y\), so (L2) controls \(1-\operatorname{Var}f\).
Since \(T\) is a gradient and a contraction,
\[
0\preceq H_T\preceq I\qquad\text{a.e.}                         \tag{5.1}
\]
Gaussian Poincar\'e and (5.1) give
\[
\operatorname{Var}f
\le\mathbb E|\nabla f|^2
=\mathbb E|H_Tu|^2
\le\mathbb E\langle u,H_Tu\rangle.                            \tag{5.2}
\]
For every symmetric matrix \(0\preceq A\preceq I\),
\(|Au|^2\le\langle u,Au\rangle\).  Therefore
\[
\begin{aligned}
\mathbb E|H_Tu-u|^2
&=1-2\mathbb E\langle u,H_Tu\rangle
  +\mathbb E|H_Tu|^2\\
&\le1-\mathbb E\langle u,H_Tu\rangle\\
&\le1-\operatorname{Var}f.                                    \tag{5.3}
\end{aligned}
\]
Combining (5.3) with (L2) proves
\[
\boxed{\quad
\mathbb E|H_Tu-u|^2
\le C_\delta\varepsilon\sqrt{\log(e/\varepsilon)}.
\quad}                                                         \tag{5.4}
\]

For the asserted product approximation, let
\(P=I-u\otimes u\), decompose \(G=Zu+W\) with
\(Z\sim N(0,1)\), \(W=PG\), and define the measurable map
\[
R(w)=\mathbb E_Z\,[P T(w+Zu)]\in u^\perp.                      \tag{5.5}
\]
The two error components are orthogonal:
\[
T(G)-u\langle G,u\rangle-R(PG)
=u(f(G)-Z)+\bigl(PT(G)-R(W)\bigr).                             \tag{5.6}
\]
Because \(\mathbb E(f(G)-Z)=0\), the \(n\)-dimensional Gaussian Poincar\'e
inequality gives
\[
\mathbb E|f(G)-Z|^2
\le\mathbb E|\nabla f(G)-u|^2
=\mathbb E|H_T(G)u-u|^2.                                     \tag{5.7}
\]
For the second component, apply one-dimensional Gaussian Poincar\'e in the
\(Z\) variable, conditionally on \(W=w\), and sum its vector coordinates:
\[
\begin{aligned}
\mathbb E|PT(G)-R(W)|^2
&\le\mathbb E|P H_T(G)u|^2\\
&\le\mathbb E|H_T(G)u-u|^2,                                  \tag{5.8}
\end{aligned}
\]
where the last inequality uses \(Pu=0\).  Equations (5.6)--(5.8) yield
\[
\boxed{\quad
\mathbb E\big|T(G)-u\langle G,u\rangle-R(PG)\big|^2
\le2\,\mathbb E|H_T(G)u-u|^2.
\quad}                                                        \tag{5.9}
\]
Also \(\mathbb ER(PG)=0\).  If no requirement \(R(w)\in u^\perp\) is imposed,
the alternative choice
\(\widetilde R(w)=\mathbb E_Z[T(w+Zu)-Zu]\) and conditional vector-valued
Poincar\'e give the sharper factor \(1\).

---

## 6. Hard-support and nonsmooth approximation

Here are the approximation details used above.

### 6.1 Scalar marginal

An extended-valued one-strongly log-concave density can be written
\[
\eta(dy)=Z^{-1}e^{-y^2/2-\psi(y)}\,dy,
\qquad \psi:\mathbb R\to(-\infty,+\infty]\text{ convex}.       \tag{6.1}
\]
Its support is an interval, possibly bounded.  Let \(\psi_k\) be finite,
smooth convex approximants obtained, for example, by taking a Moreau envelope
with parameter tending to zero and then smoothing it by convolution.  Add an
irrelevant constant so that the normalizations are absorbed, and put
\[
\eta_k(dy)=Z_k^{-1}e^{-y^2/2-\psi_k(y)}\,dy.                    \tag{6.2}
\]
The approximants may be chosen to converge pointwise on the interior of the
support and to diverge off its closure.  Convex functions have an affine lower
support, so the Gaussian factor supplies an integrable dominating envelope
after a fixed affine recentering.  Consequently
\[
\eta_k\Rightarrow\eta,
\qquad F_{\eta_k}^{-1}(s)\to F_\eta^{-1}(s)\quad(0<s<1).       \tag{6.3}
\]
For every \(k\), the smooth one-dimensional Caffarelli argument gives that
\(q_k=F_{\eta_k}^{-1}\circ\Phi\) is 1-Lipschitz.  The pointwise limit in
(6.3) is therefore 1-Lipschitz as well.  This proves (1.1) for hard supports.
All identities in Sections 2--4 then apply directly to the limiting Lipschitz
map; no boundary mass occurs because (6.1) is absolutely continuous.

### 6.2 The multivariate Brenier map

Write
\[
\pi(dx)=Z^{-1}e^{-t|x|^2/2-\Psi(x)}\,dx,
\qquad \Psi:\mathbb R^n\to(-\infty,+\infty]\text{ convex}.     \tag{6.4}
\]
Approximate \(\Psi\) by finite smooth convex \(\Psi_k\) as above.  The target
laws converge weakly and in second moment, and stability of quadratic optimal
transport gives convergence of the Brenier maps in \(L^2(\gamma_n)\), after
passing to the unique limit map.  Each standardized map is a gradient
contraction.  Hence the limit \(T\) is also 1-Lipschitz and cyclically
monotone.  Rademacher's theorem and Alexandrov's theorem then give, almost
everywhere,
\[
H_T=H_T^{\mathsf T},\qquad 0\preceq H_T\preceq I.              \tag{6.5}
\]
Equations (5.2)--(5.9) use only weak derivatives of this Lipschitz map and
therefore remain valid.  This deals simultaneously with convex hard support
and nonsmooth potential; there is no unrecorded boundary term.

---

## 7. Proof of (L4): two-dimensional halfspace identifiability

We spell out what is meant by oriented normal distance.  Represent an oriented
halfspace as
\[
H(n,a)=\{x:\langle n,x\rangle\ge a\},\qquad |n|=1,
\]
and set
\[
d_{\mathrm{or}}((n,a),(n',a'))=|n-n'|+|a-a'|.                 \tag{7.1}
\]
Replacing the sum by the Euclidean product metric changes only an absolute
factor.

We prove the following fixed-dimensional lemma.

**Lemma 7.1 (uniform transversality).**  For \(d\in\{1,2\}\) and
\(0<\delta\le1/2\), there is \(c_\delta>0\) such that, for every isotropic
log-concave probability \(\nu\) on \(\mathbb R^d\) and every two oriented
halfspaces \(H(n_i,a_i)\) satisfying
\[
\nu(H(n_i,a_i))\in[\delta,1-\delta],
\]
one has
\[
\nu(H(n_1,a_1)\mathbin{\triangle} H(n_2,a_2))
\ge c_\delta\bigl(|n_1-n_2|+|a_1-a_2|\bigr).                  \tag{7.2}
\]

**Proof.**  We include the compactness and the infinitesimal cases, since the
latter is where linear, rather than merely qualitative, control enters.

First, every admissible threshold lies in a compact interval depending only
on \(\delta\).  Indeed, for \(Y=\langle X,n\rangle\), isotropy gives
\(\mathbb EY=0\), \(\mathbb EY^2=1\).  If \(a>0\), then
\(\delta\le\mathbb P(Y\ge a)\le a^{-2}\); if \(a<0\), apply the same argument
to \(-Y\) and the mass of the complementary side.  Thus
\[
|a|\le\delta^{-1/2}.                                           \tag{7.3}
\]

Suppose (7.2) were false.  There would be isotropic log-concave laws \(\nu_k\)
and admissible parameter pairs \(\theta_k=(n_k,a_k)\),
\(\theta'_k=(n'_k,a'_k)\), such that, with
\[
d_k=|n_k-n'_k|+|a_k-a'_k|,
\qquad
\rho_k=\nu_k(H(\theta_k)\mathbin{\triangle} H(\theta'_k)),
\]
we have \(\rho_k/d_k\to0\).  Isotropic log-concave laws in a fixed dimension
have a uniform exponential tail, hence form a tight family with uniformly
integrable second moments.  Passing to a subsequence,
\[
\nu_k\Rightarrow\nu,
\]
where \(\nu\) is again isotropic and log-concave.  In particular it is
full-dimensional.  By (7.3) and compactness of the unit sphere, the two
parameter sequences also converge.

If \(d_k\) does not tend to zero, the limiting oriented halfspaces are
distinct.  Their boundaries have \(\nu\)-measure zero, so convergence of the
measures and parameters gives
\[
\nu(H(\theta)\mathbin{\triangle} H(\theta'))=0.                \tag{7.4}
\]
This is impossible.  A full-dimensional log-concave law has a density that is
positive on the interior of its convex support.  Since each limiting
halfspace and its complement have mass at least \(\delta\), each boundary
crosses that interior.  Two distinct oriented affine halfspaces then disagree
on a nonempty open subset of the interior, which has positive measure.  Thus
necessarily \(d_k\to0\), and the two parameter limits agree; call the limit
\((n,a)\).

Normalize the increments.  After another subsequence,
\[
\frac{n'_k-n_k}{d_k}\to w,
\qquad
\frac{a'_k-a_k}{d_k}\to s,
\qquad |w|+|s|=1.                                              \tag{7.5}
\]
Because \(|n_k|=|n'_k|=1\), we have \(w\perp n\).  Let
\[
L=\{x:\langle n,x\rangle=a\}.
\]
The affine function
\[
M(x)=\langle w,x\rangle-s                                    \tag{7.6}
\]
is not identically zero on \(L\).  Indeed, an affine function vanishing on
\(L\) has linear part parallel to \(n\); together with \(w\perp n\), that
would force \(w=s=0\), contrary to (7.5).

The balanced-mass condition makes \(L\) cross the interior of the support of
\(\nu\).  Choose a compact line segment
\(J\subset L\cap\operatorname{int}(\operatorname{supp}\nu)\), away from the
at most one zero of \(M|_L\), on which
\[
|M|\ge\eta>0.                                                   \tag{7.7}
\]
The density \(f\) of \(\nu\) has a positive lower bound on a small rectangle
around \(J\).  A standard elementary convergence fact for full-dimensional
log-concave densities says that weak convergence \(\nu_k\Rightarrow\nu\)
implies locally uniform convergence of their densities on compact subsets of
the interior of the limiting support.  Hence, on a slightly smaller
rectangle, the densities \(f_k\) are bounded below by a common number
\(m>0\).

Use normal-tangential coordinates around the line
\(\langle n_k,x\rangle=a_k\).  On the tangential copy of \(J\), the two affine
defining functions differ by
\[
d_k(M(x)+o(1)).
\]
For each tangential coordinate, their signs therefore differ on a normal
interval of length at least \(\eta d_k/2\), for all sufficiently large
\(k\).  Integrating these intervals against the lower density bound \(m\)
gives
\[
\rho_k\ge \frac{m\eta\,\mathcal H^{d-1}(J)}2\,d_k,             \tag{7.8}
\]
contrary to \(\rho_k/d_k\to0\).  In dimension one the same argument has no
tangential coordinate: for equal orientations it is simply the positive
density at a balanced quantile; opposite orientations belong to the preceding
noninfinitesimal case.  This proves the lemma. \(\square\)

For reference, the two compactness facts used in the proof do include hard
supports.  Uniform exponential tails follow from the one-dimensional
log-concave tail lemma applied to finitely many directions in a fixed net.
Closure of log-concavity follows directly from the defining set inequality.
Uniform integrability preserves mean zero and covariance identity, so the
limit remains full-dimensional.  Finally, local uniform convergence of
densities follows by applying concavity to \(\log f_k\) on a simplex compactly
contained in the interior; pointwise bounds at the vertices give uniform
bounds and equicontinuity on smaller simplices.  Thus no smoothness or
full-support assumption is hidden in Lemma 7.1.

Taking \(C_\delta=c_\delta^{-1}\) in (7.2) proves (L4):
\[
\boxed{\quad
d_{\mathrm{or}}((n_1,a_1),(n_2,a_2))
\le C_\delta\,\nu(H_1\mathbin{\triangle} H_2).
\quad}                                                         \tag{7.9}
\]

---

## 8. Equality and endpoint checks

* If \(\varepsilon=0\), then \(D_q=D_p=0\).  Since \(K_g>0\) on every finite
  point, (2.5) gives \(q'=1\) almost everywhere.  The centering forces
  \(q(z)=z\), while (3.10) forces \(p=h\) almost everywhere.  Thus the scalar
  marginal is exactly Gaussian and \(S=H\) modulo null sets.  Equations
  (5.2)--(5.3) then force \(H_Tu=u\) almost everywhere.
* The restriction \(g\in[\delta,1-\delta]\) is used exactly once in the set
  estimate, through \(W\le K_g/\delta\), and once in the variance estimate,
  through the inverse Mills comparison with \(K_g\).  All quantile intervals
  in (3.9) remain inside \([0,1]\) because the exchanged mass obeys
  \(\alpha\le\min(g,1-g)\).
* For \(\varepsilon=1\), the right side in (L1) is a \(\delta\)-dependent
  constant and (L2)--(L3) are bounded by the contraction inequalities.  Thus
  no separate large-deficit case is missing.
* Convex hard support produces no atoms in a one-dimensional marginal and no
  mass on an affine hyperplane.  Hence the choices of \(\ge\) versus \(>\) in
  the definitions of the halfspaces do not affect any displayed identity.

This completes the clean-room proof of (L1)--(L4).
